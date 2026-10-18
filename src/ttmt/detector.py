"""Pick-free detector core: timing engine, peak reading, radial maps, MT-FW and BOMP-PK greedy loops with deflation.

Ported from ``code/explore_contrib/c10_lib.py`` (Eng, choose in mode "hard", _vectors, _deflate, mt_loop, bomp_energy,
_gram, bomp_loop), ``task2d_methods.peak_table``, ``c8_lib`` (stat_huber, mtfw_general with stat "huber_n") and
``c18c_oracle.compress``. Under the legacy conventions (Kernel 2 Hz, rounding "floor", radial map cap 6 or none) every
function returns the same numbers as its source; the conventions are arguments, never module state.

Shapes: MFn (ns, 3, NT) normalised matched-filter output; engine tables IDX, frac (N, T, 2, ns) for N nodes and T
origins; G (N, 6 ns, 6) radiation operators; Q (N, 6 ns, 6) their orthonormal bases.
"""
from dataclasses import dataclass

import numpy as np

from .conventions import LEGACY_ROUNDING, sample_index_clipped
from .kernel import LEGACY as LEGACY_KERNEL
from .station_terms import PeakDelay

SH = 4                                     # peak half-width in samples (0.25 s at 16 Hz): task1x_mtglrt.SH
EIG_FLOOR = 1e-10                          # task2c_methods.EIG_FLOOR
NT_LEGACY = 640                            # 40 s windows at 16 Hz


# ============================================================================================ radial maps
@dataclass(frozen=True)
class RadialMap:
    """map applied to each 3-component residual vector read at a station and phase, before projection.
    kind "none" (BOMP-PK, NOCAP), "cap" (registered MT-FW, c = 6: |v| limited to c), "log" (LOG6 of C18: |v|^2 ->
    c^2 log(1 + |v|^2 / c^2))."""
    kind: str = "none"
    c: float = 0.0

    def __post_init__(self):
        if self.kind not in ("none", "cap", "log"):
            raise ValueError(self.kind)

    def __call__(self, v):
        if self.kind == "cap":
            nv = np.linalg.norm(v, axis=-1, keepdims=True)
            return v * np.minimum(1.0, self.c / np.maximum(nv, 1e-12))
        if self.kind == "log":
            e = (v * v).sum(-1, keepdims=True)
            return v * np.sqrt(self.c * self.c * np.log1p(e / (self.c * self.c)) / np.maximum(e, 1e-24))
        return v

    @classmethod
    def from_capv(cls, capv):
        """legacy argument capv of c10_lib.mt_loop / c18d_detect: None = none, > 0 = cap, < 0 = log with c = -capv."""
        if capv is None:
            return cls("none")
        return cls("cap", float(capv)) if capv > 0 else cls("log", float(-capv))


NOCAP = RadialMap("none")
CAP6 = RadialMap("cap", 6.0)
LOG6 = RadialMap("log", 6.0)


def huber_rho(r, c):
    """c8_lib._rho: Huber loss of the norm r with threshold c."""
    return np.where(r <= c, 0.5 * r * r, c * r - 0.5 * c * c)


C_HUB = float(np.sqrt(7.814727903251178))      # chi2_3 0.95 quantile (c8_lib.C_HUB)
MED_CHI3 = float(np.sqrt(2.365973884375338))   # sqrt of the chi2_3 median


# ============================================================================================ engine
class Engine:
    """timing tables and operators of a node set (c10_lib.Eng, timing part and BOMP operators).

    tt (N, 2, ns): travel times of every node to every station (P, S); taus (T,): origin times in window time;
    terms: PeakDelay (energy-peak station terms, E54) added to the predicted times; rounding: sample-index convention
    (E52); G, Q: radiation operators and bases (N, 6 ns, 6)."""

    def __init__(self, XN, STA, tt, taus, terms, G, Q, dt=1.0 / 16, nt=NT_LEGACY, rounding=LEGACY_ROUNDING, bomp=False):
        if not isinstance(terms, PeakDelay):
            raise TypeError(f"an energy-peak engine needs PeakDelay terms, got {type(terms).__name__} (E54)")
        self.XN = np.asarray(XN, float)
        self.STA = np.asarray(STA, float)
        self.ns = len(self.STA)
        self.taus = np.asarray(taus, float)
        self.N, self.T = len(self.XN), len(self.taus)
        self.dt, self.nt, self.rounding = dt, nt, rounding
        self.terms = terms
        self.G, self.Q = G, Q
        self.st4 = np.arange(self.ns)[None, None, None, :]
        self.set_times(tt, terms)
        if bomp:
            Gb = self.G.reshape(self.N, 2, self.ns, 3, 6)
            self.GtG = np.einsum("nxk,nxl->nkl", self.G, self.G)
            self.Cps = np.einsum("nick,nicl->nikl", Gb[:, 0], Gb[:, 1])
            self.Kinv0 = pinv_sym(self.GtG)

    def set_times(self, tt, terms):
        """(re)build IDX and frac from travel times and terms (r1_run.engine and c18j_layered.layered_engine replace
        the tables of an existing engine the same way)."""
        if not isinstance(terms, PeakDelay):
            raise TypeError(f"an energy-peak engine needs PeakDelay terms, got {type(terms).__name__} (E54)")
        self.terms = terms
        x = (np.asarray(tt)[:, None] + self.taus[None, :, None, None] + terms.values.T[None, None]) / self.dt
        self.IDX, self.frac = sample_index_clipped(x, self.nt, self.rounding)
        return self


def pinv_sym(K):
    """task2d_methods.PKDict.pinv_sym: pseudo-inverse of symmetric matrices, eigenvalues below EIG_FLOOR x max cut."""
    ev, V = np.linalg.eigh(K)
    keep = ev > EIG_FLOOR * ev[..., -1:]
    inv = np.where(keep, 1.0 / np.where(keep, ev, 1.0), 0.0)
    return np.einsum("...kj,...j,...lj->...kl", V, inv, V)


# ============================================================================================ reading
def peak_table(R, sh=SH):
    """(ns, NT): for every station and centre sample t0, the sample of clip(t0 + o, 0, NT - 1), |o| <= sh, with the
    largest 3-component energy (first maximum). task2d_methods.peak_table (sh = 4) and xc.peak_table_w."""
    nt = R.shape[-1]
    offs = np.arange(-sh, sh + 1)
    ii = np.clip(np.arange(nt)[:, None] + offs[None, :], 0, nt - 1)
    Et = np.sum(np.ascontiguousarray(R.transpose(0, 2, 1)) ** 2, axis=-1)
    a = Et[:, ii].argmax(axis=-1)
    return ii[np.arange(nt)[None, :], a]


def read_peaks(E, R, sh=SH):
    """registered rule ("hard"): read every block at the energy peak within +-sh samples of its predicted sample."""
    return peak_table(R, sh)[E.st4, E.IDX]


def vectors(R, tpk, radial=NOCAP):
    """(..., 2, ns, 3) residual vectors at samples tpk (..., 2, ns), after the radial map (c10_lib._vectors)."""
    ns = R.shape[0]
    V = np.ascontiguousarray(R.transpose(0, 2, 1))
    v = V[np.arange(ns).reshape((1,) * (tpk.ndim - 1) + (ns,)), tpk]
    return radial(v) if radial.kind != "none" else v


def deflate(R, pred, tp, kernel=LEGACY_KERNEL):
    """subtract the fitted atom: at every station and phase, pred (3,) times the deflation response of the kernel
    centred on the read sample (c10_lib._deflate). In place."""
    ns, nt = R.shape[0], R.shape[-1]
    AC, H = kernel.ac, kernel.ac_half
    for b in range(2):
        for i in range(ns):
            t0 = int(tp[b, i])
            lo, hi = max(0, t0 - H), min(nt, t0 + H + 1)
            R[i, :, lo:hi] -= np.outer(pred[b, i], AC[H - (t0 - lo):H + (hi - t0)])
    return R


def score_map(E, MFn, radial=NOCAP):
    """first-iteration subspace statistic of every (node, origin) (c18c_oracle.score_map)."""
    tpk = read_peaks(E, MFn)
    vf = vectors(MFn, tpk, radial).reshape(E.N, E.T, -1)
    y = np.einsum("gtr,grk->gtk", vf, E.Q)
    return (y * y).sum(-1) / 2.0


# ============================================================================================ block statistics
class BlockStatistic:
    """statistic of one hypothesis from the peak vectors of its blocks (one block = one station and phase).

    score(vf, Q) -> (scores (N, T), aux): vf (N, T, 6 ns) raw peak vectors, Q (N, 6 ns, 6) bases.
    deflation(E, gi, tj, vf, aux) -> (2, ns, 3) fitted vectors of the selected hypothesis, removed by the loop.
    rho(r): contribution to the score of a block of norm r that the model explains exactly.
    influence(r): d rho / d r, how much one block of norm r can push the score."""
    name = "abstract"

    def score(self, vf, Q):
        raise NotImplementedError

    def deflation(self, E, gi, tj, vf, aux):
        raise NotImplementedError

    def rho(self, r):
        raise NotImplementedError

    def influence(self, r):
        raise NotImplementedError


class ProjectedEnergy(BlockStatistic):
    """||Q^T phi(v)||^2 / 2 with phi a radial map (none: MT-FW NOCAP of R1; cap c: registered MT-FW cap 6; log c:
    LOG6 of C18). Deflation: least-squares moment tensor of the mapped vector, G m (c10_lib.mt_loop)."""

    def __init__(self, radial=NOCAP):
        self.radial = radial
        self.name = radial.kind if radial.kind == "none" else f"{radial.kind}{radial.c:g}"

    def score(self, vf, Q):
        v = vf if self.radial.kind == "none" else self.radial(vf.reshape(vf.shape[:-1] + (-1, 3))).reshape(vf.shape)
        y = np.einsum("gtr,grk->gtk", v, Q)
        return (y * y).sum(-1) / 2.0, v

    def deflation(self, E, gi, tj, vf, aux):
        G = E.G[gi]
        m, *_ = np.linalg.lstsq(G, aux[gi, tj], rcond=None)
        return (G @ m).reshape(2, E.ns, 3)

    def rho(self, r):
        r = np.asarray(r, float)
        if self.radial.kind == "cap":
            return 0.5 * np.minimum(r, self.radial.c) ** 2
        if self.radial.kind == "log":
            c = self.radial.c
            return 0.5 * c * c * np.log1p(r * r / (c * c))
        return 0.5 * r * r

    def influence(self, r):
        r = np.asarray(r, float)
        if self.radial.kind == "cap":
            return np.where(r < self.radial.c, r, 0.0)
        if self.radial.kind == "log":
            return r / (1.0 + r * r / self.radial.c ** 2)
        return r


class Huber(BlockStatistic):
    """HUBER-N (c8_lib.stat_huber, scale_adapt False): sum_b rho_c(|v_b|) - min_a sum_b rho_c(|v_b - Q_b a|),
    c = sqrt(chi2_3 0.95). Deflation: Q a of the robust fit (c8_lib.mtfw_general)."""
    name = "huber"

    def __init__(self, c=None, scale_adapt=False):
        self.c, self.scale_adapt = c, scale_adapt

    def score(self, vf, Q):
        return huber_stat(vf, Q, self.scale_adapt, chub=self.c)

    def deflation(self, E, gi, tj, vf, aux):
        return (E.Q[gi] @ aux[gi, tj]).reshape(2, E.ns, 3)

    def rho(self, r):
        return huber_rho(np.asarray(r, float), C_HUB if self.c is None else self.c)

    def influence(self, r):
        return np.minimum(np.asarray(r, float), C_HUB if self.c is None else self.c)


# ============================================================================================ readers
class Reader:
    """per-station reading step: from residual traces R (ns, 3, NT) and the predicted samples idx (..., 2, ns) of
    every block (engine tables, rounding convention included), the samples read (same shape as idx)."""

    def read(self, R, idx):
        raise NotImplementedError


class EnergyArgmax(Reader):
    """registered rule: the 3-component energy peak within +-half_width samples of the predicted sample (first
    maximum; samples clipped to the window)."""

    def __init__(self, half_width=SH):
        self.half_width = int(half_width)

    def read(self, R, idx):
        ns = R.shape[0]
        st = np.arange(ns).reshape((1,) * (idx.ndim - 1) + (ns,))
        return peak_table(R, self.half_width)[st, idx]


ENERGY_ARGMAX = EnergyArgmax()


# ============================================================================================ MT-FW family
def greedy_loop(E, MFn, kmax, statistic, reader=ENERGY_ARGMAX, min_sep=5.0, kernel=LEGACY_KERNEL):
    """MT-FW family greedy loop (c10_lib.mt_loop in mode "hard" with W = SH and no refinement; c8_lib.mtfw_general
    for HUBER-N): read every block, score every (node, origin), take the best not excluded (lowest node, then first
    origin, among ties), exclude the nodes within min_sep km of its epicentre, deflate its fitted atom.
    Atoms: dicts x (3,), tau, s (statistic at its step), run (running minimum of s)."""
    ns = E.ns
    R = MFn.copy()
    excluded = np.zeros(E.N, bool)
    atoms, run = [], np.inf
    for _ in range(kmax):
        if excluded.all():
            break
        tpk = reader.read(R, E.IDX)
        vf = vectors(R, tpk, NOCAP).reshape(E.N, E.T, 6 * ns)
        sc, aux = statistic.score(vf, E.Q)
        sc[excluded] = -np.inf
        gi = int(np.argmax(sc.max(axis=1)))
        tj = int(np.argmax(sc[gi]))
        s = float(sc[gi, tj])
        if not np.isfinite(s):
            break
        src = np.asarray(E.XN[gi], float)
        run = min(run, s)
        atoms.append({"s": s, "run": float(run), "x": src, "tau": float(E.taus[tj])})
        excluded = excluded | (np.linalg.norm(E.XN[:, :2] - src[None, :2], axis=1) <= min_sep)
        deflate(R, statistic.deflation(E, gi, tj, vf, aux), tpk[gi, tj].copy(), kernel)
    return atoms


def mt_loop(E, MFn, kmax, capv=6.0, min_sep=5.0, kernel=LEGACY_KERNEL, reader=ENERGY_ARGMAX):
    """legacy signature of c10_lib.mt_loop (capv: None no cap, > 0 cap, < 0 log as in c18d_detect)."""
    return greedy_loop(E, MFn, kmax, ProjectedEnergy(RadialMap.from_capv(capv)), reader, min_sep, kernel)


# ============================================================================================ BOMP-PK
def bomp_energy(E, R, tpk, kernel=LEGACY_KERNEL):
    """e = h^T K(I)^+ h of every block at read samples tpk (c10_lib.bomp_energy, all nodes)."""
    N, T_ = tpk.shape[:2]
    v = vectors(R, tpk, NOCAP).reshape(N, T_, -1)
    h = np.matmul(v, E.G)
    lag = tpk[:, :, 1, :] - tpk[:, :, 0, :]
    wm, wp = kernel.ac_at(-lag), kernel.ac_at(lag)
    diff = np.any((wm != 0) | (wp != 0), axis=-1)
    Kinv = np.broadcast_to(E.Kinv0[:, None], (N, T_, 6, 6))
    if diff.any():
        Kinv = Kinv.copy()
        n_, t_ = np.nonzero(diff)
        K = E.GtG[n_] + np.einsum("mi,mikl->mkl", wm[n_, t_], E.Cps[n_]) + np.einsum("mi,milk->mkl", wp[n_, t_], E.Cps[n_])
        Kinv[n_, t_] = pinv_sym(0.5 * (K + np.swapaxes(K, -1, -2)))
    return np.einsum("ntk,ntk->nt", h, np.matmul(Kinv, h[..., None])[..., 0])


def _gram(Gb, I, kernel):
    GA, IA = np.stack(Gb), np.stack(I)
    Wt = kernel.ac_at(IA[:, None, :, None, :] - IA[None, :, None, :, :])
    n = len(Gb)
    K = np.einsum("apick,bqicl,abpqi->akbl", GA, GA, Wt, optimize=True).reshape(6 * n, 6 * n)
    return 0.5 * (K + K.T)


def bomp_loop(E, MFn, kmax, min_sep=5.0, kernel=LEGACY_KERNEL):
    """BOMP-PK greedy loop (c10_lib.bomp_loop, mode "hard", no refinement): block energy with the pseudo-inverse of
    the block Gram at the read samples, joint least-squares refit of all accepted atoms and re-deflation from MFn."""
    ns = E.ns
    stv = np.arange(ns)
    excluded = np.zeros(E.N, bool)
    Gs, Is, atoms = [], [], []
    R = MFn.copy()
    for _ in range(kmax):
        if excluded.all():
            break
        tpk = read_peaks(E, R)
        sc = bomp_energy(E, R, tpk, kernel) / 2.0
        e0 = sc.copy()
        sc[excluded] = -np.inf
        j = int(np.argmax(sc))
        gi, tj = divmod(j, E.T)
        s = float(sc[gi, tj])
        if not np.isfinite(s):
            break
        src, tau, tp, G = E.XN[gi], float(E.taus[tj]), tpk[gi, tj].copy(), E.G[gi]
        atoms.append({"s": s, "s0": float(e0[gi, tj]), "x": np.asarray(src, float), "tau": tau})
        Gs.append(G)
        Is.append(tp)
        excluded = excluded | (np.linalg.norm(E.XN[:, :2] - src[None, :2], axis=1) <= min_sep)
        Kg = _gram([g.reshape(2, ns, 3, 6) for g in Gs], Is, kernel)
        V = MFn[stv[None, None, :], :, np.stack(Is)]
        b = np.einsum("bx,bxk->bk", V.reshape(len(Is), -1), np.stack(Gs)).ravel()
        U, *_ = np.linalg.lstsq(Kg, b, rcond=None)
        U = U.reshape(len(Is), 6)
        R = MFn.copy()
        for g, I, u in zip(Gs, Is, U):
            deflate(R, (g @ u).reshape(2, ns, 3), I, kernel)
    return atoms


# ============================================================================================ HUBER-N
def huber_stat(vf, Q, scale_adapt=False, iters=4, chub=None):
    """c8_lib.stat_huber: robust likelihood ratio by iteratively reweighted least squares (4 iterations, ridge 1e-9).
    vf (N, T, 6 ns); Q (N, 6 ns, 6). Returns scores (N, T) and coefficients (N, T, 6) in the Q basis."""
    N, T_, R6 = vf.shape
    nb = R6 // 3
    vb = vf.reshape(N, T_, nb, 3)
    Qb = Q.reshape(N, nb, 3, 6)
    QQ = np.matmul(Qb.transpose(0, 1, 3, 2), Qb).reshape(N, nb, 36)
    PvT = np.matmul(Qb.transpose(0, 1, 3, 2), vb.transpose(0, 2, 3, 1))
    QT = Q.transpose(0, 2, 1)
    a = PvT.sum(1).transpose(0, 2, 1)
    nv = np.sqrt((vb * vb).sum(-1))

    def rnorm(a_):
        r_ = vf - np.matmul(a_, QT)
        return np.sqrt((r_.reshape(N, T_, nb, 3) ** 2).sum(-1))

    def wsolve(w, ridge=1e-9):
        A = np.matmul(w, QQ).reshape(N, T_, 6, 6) + ridge * np.eye(6)
        rhs = np.einsum("ntb,nbkt->ntk", w, PvT, optimize=True)
        return np.linalg.solve(A, rhs[..., None])[..., 0]
    r = rnorm(a)
    sig = np.maximum(1.0, np.median(r, axis=-1) / MED_CHI3) if scale_adapt else np.ones((N, T_))
    c = ((C_HUB if chub is None else chub) * sig)[..., None]
    for _ in range(iters):
        w = np.minimum(1.0, c / np.maximum(r, 1e-12))
        a = wsolve(w)
        r = rnorm(a)
    return huber_rho(nv, c).sum(-1) - huber_rho(r, c).sum(-1), a


MTFW_CAP6 = ProjectedEnergy(CAP6)
MTFW_NOCAP = ProjectedEnergy(NOCAP)
HUBER_N = Huber()


def huber_loop(E, MFn, kmax, min_sep=5.0, kernel=LEGACY_KERNEL, chub=None, reader=ENERGY_ARGMAX):
    """HUBER-N greedy loop (c8_lib.mtfw_general, stat "huber_n"); the registered R0 M4 score is the running minimum."""
    return greedy_loop(E, MFn, kmax, Huber(chub), reader, min_sep, kernel)
