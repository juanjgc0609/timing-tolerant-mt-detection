"""
c10_lib.py - C10: soft timing tolerance (Gaussian timing prior) inside BOMP-PK and MT-FW (cap 6).  SEED 999 ONLY.

Registered detectors read every station/phase at the energy peak within +-4 samples (+-0.25 s) of the predicted sample
of the hypothesis (node, origin); any shift inside the window is free and independent across stations.  Here the
chosen shifts delta_bi (peak sample time minus the CONTINUOUS predicted time d / V_b + tau + corr_bi) are charged a
Gaussian log-likelihood:

  hard   registered rule (window half-width W samples, no charge)
  indep  sample maximising  E_bi(k) / 2 - delta^2 / (2 s2_bi)   (E capped at 36 in the MT-FW family, as its vectors
         are); node statistic  s - sum_bi delta_bi^2 / (2 s2_bi).  s2_bi is a per (node, phase, station) variance.
  corr   prior covariance of the shift vector of a node  C = diag(sig_b^2) + J S J^T, J = d t / d (x, y, z, tau) of the
         straight-ray model, S = diag(h^2, h^2, hz^2, dtau^2) / 12 (uniform position of the event inside the lattice
         cell of the node and inside the origin step): the shifts are a location / origin-time term already searched
         plus independent residuals.  Peaks chosen with the marginal variances diag(C); theta = S J^T C^-1 delta;
         peaks re-chosen once among the same samples with (delta - J theta)^2 / (2 sig_b^2); node statistic
         s - delta^T C^-1 delta / 2.
  refine after an atom is selected by the HARD rule, it is relocated to the maximum of the indep objective on a 1 km
         lattice (same depths) within R_REF km (epicentral) and origins within T_REF s, outside the exclusion zones of
         the earlier atoms; ties and no gain keep the atom.  Exclusion and deflation use the relocated atom.

Timing level sig_b (phase b), fixed by rule from the catalogue of the network (sigma_catalogue): robust sd (1.4826
MAD) of the +SC peak delays of the gated catalogue arrivals minus the station correction, the event median removed
(common mode = origin-time error), per phase.  ORACLE: the same statistic at the TRUE catalogue positions and origin
times with the TRUE statics.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
for p in (HERE, WORK):
    if p not in sys.path:
        sys.path.insert(0, p)

import c9_lib as C9  # noqa: E402
import xc  # noqa: E402
from station_corrections import peak_delays  # noqa: E402
from task0_common import DT, NT, VP, VS  # noqa: E402

M, D2, X = C9.M, C9.D2, C9.X
SEED = 999
CAPV = 6.0
SH = D2.SH                      # 4 samples = 0.25 s
H_COARSE, HZ, DTAU = 3.0, 4.0, 0.1
H_FINE = 1.0
R_REF = 10.0                    # km
T_REF = 1.0                     # s
SIG_FLOOR = 0.05                # s
SIG_DEFAULT = 0.144             # s, used only if a catalogue has fewer than 10 usable residuals of a phase
VINV = np.array([1.0 / VP, 1.0 / VS])


# ============================================================================================ timing level
def em_sigma(r, half=1.0, iters=200):
    """sd of the Gaussian core of residuals r contaminated by wrong peaks: mixture pi N(0, sig^2) + (1 - pi) U(-half,
    half) fitted by EM (half = the +-1 s search of the +SC rule).  Returns sig, pi."""
    r = np.asarray(r, float)
    sig = float(np.clip(1.4826 * np.median(np.abs(r)), 0.05, 0.3))
    pi = 0.7
    u = 1.0 / (2.0 * half)
    for _ in range(iters):
        g = pi * np.exp(-0.5 * (r / sig) ** 2) / (sig * np.sqrt(2 * np.pi))
        w = g / (g + (1.0 - pi) * u)
        sig = float(np.sqrt(max(1e-12, (w * r * r).sum() / max(w.sum(), 1e-12))))
        pi = float(np.clip(w.mean(), 0.05, 0.999))
    return sig, pi


def catalogue_residuals(STA, cat, corr, gate=25.0):
    res = [[], []]
    for D, src, tau in cat:
        dl, en = peak_delays(STA, D, np.asarray(src, float), float(tau))
        ok = en >= gate
        m = int(ok.sum())
        if m < 3:
            continue
        r = dl - corr
        r = (r - np.median(r[ok])) * np.sqrt(m / (m - 1.0))
        for b in range(2):
            res[b].extend(r[ok[:, b], b].tolist())
    return res


def sigma_catalogue(STA, cat, corr, gate=25.0):
    """RULE (fixed on the catalogues of validation network 195 before any detection run): robust sd (1.4826 MAD about
    zero) of the +SC peak delays of the gated catalogue arrivals minus the station correction, the event median
    removed (common mode = catalogue origin-time error), POOLED over stations and phases; floor SIG_FLOOR; SIG_DEFAULT
    with fewer than 10 residuals.  cat: [(D, src, tau)].  Returns sig (2,) (the same value for P and S) and the
    counts per phase.  (The per-phase MAD and the Gaussian-core EM estimate em_sigma were tried on network 195 and
    discarded as unstable with 20-70 residuals per phase: 0.03-0.38 s.)"""
    res = catalogue_residuals(STA, cat, corr, gate)
    cnt = np.array([len(res[0]), len(res[1])])
    r = np.array(res[0] + res[1])
    sg = max(SIG_FLOOR, 1.4826 * float(np.median(np.abs(r)))) if len(r) >= 10 else SIG_DEFAULT
    return np.array([sg, sg]), cnt


def sigma_oracle(STA, cat_truth):
    """ORACLE timing level: sd per phase of the TRUE model error (travel time of the generator's gradient medium minus
    the detectors' homogeneous travel time) over the true catalogue sources, the station median (what a perfect station
    correction absorbs) and the event mean (origin time) removed.  Used with the TRUE statics as corrections."""
    import gen_realistic as GR
    me = []
    for src, _tau in cat_truth:
        d = np.linalg.norm(STA - src, axis=1)
        me.append([GR.travel_time_gradient(src, STA, ph, GR.make_cfg(None)) - d / V for ph, V in (("P", VP), ("S", VS))])
    me = np.array(me)                                                             # (K, 2, ns)
    me = me - np.median(me, axis=0, keepdims=True)
    me = me - me.mean(axis=(1, 2), keepdims=True)
    return np.maximum(0.03, me.std(axis=(0, 2)))


# ============================================================================================ engine
class Eng:
    """timing and amplitude operators of a node set."""

    def __init__(self, XN, STA, corr, taus, h, hz=HZ, dtau=DTAU, G=None, Q=None, bomp=False):
        self.XN = np.asarray(XN, float)
        self.STA = np.asarray(STA, float)
        self.ns = ns = len(self.STA)
        self.taus = np.asarray(taus, float)
        self.N, self.T = len(self.XN), len(self.taus)
        n, d = C9.rays_hom(self.XN, self.STA)
        self.n = n
        tt = np.stack([d / VP, d / VS], axis=1)                                   # (N, 2, ns)
        x = (tt[:, None] + self.taus[None, :, None, None] + np.asarray(corr, float).T[None, None]) / DT
        self.IDX = np.clip(x.astype(int), 0, NT - 1)
        self.frac = x - self.IDX
        self.G = C9.operators_hom(self.XN, self.STA) if G is None else G
        self.Q = C9.bases(self.G) if Q is None else Q
        # d t_bi / d (x, y, z, tau): -n_i / V_b, 1
        J = np.concatenate([-n[:, None] * VINV[None, :, None, None], np.ones((self.N, 2, ns, 1))], axis=-1)
        self.J = J.reshape(self.N, 2 * ns, 4)
        self.S = np.array([h * h, h * h, hz * hz, dtau * dtau]) / 12.0
        self.cellvar = np.einsum("nra,a->nr", self.J ** 2, self.S).reshape(self.N, 2, ns)
        self.st4 = np.arange(ns)[None, None, None, :]
        self.cache = {}
        if bomp:
            Gb = self.G.reshape(self.N, 2, ns, 3, 6)
            self.GtG = np.einsum("nxk,nxl->nkl", self.G, self.G)
            self.Cps = np.einsum("nick,nicl->nikl", Gb[:, 0], Gb[:, 1])
            self.Kinv0 = D2.PKDict.pinv_sym(self.GtG)

    def prior(self, sig):
        """sig (2,) per phase: marginal variances s2 (N, 2, ns), the residual-only variance sc2 (2, 1), C^-1
        (N, 2ns, 2ns) and the map theta = Tm delta (N, 4, 2ns)."""
        sc2 = np.repeat(np.asarray(sig, float) ** 2, self.ns)                     # (2ns,)
        Cm = np.einsum("nra,a,nsa->nrs", self.J, self.S, self.J) + np.diag(sc2)[None]
        Cinv = np.linalg.inv(Cm)
        Tm = np.einsum("a,nra,nrs->nas", self.S, self.J, Cinv)
        return {"s2": sc2.reshape(2, self.ns)[None] + self.cellvar, "sc2": sc2.reshape(2, self.ns), "Cinv": Cinv,
                "Tm": Tm}


def pick(Eh, IDX, frac, s2, W, shift=None):
    """Eh (ns, NT) half energies.  IDX, frac (..., 2, ns).  s2 broadcastable to IDX (None = no charge).
    shift (..., 2, ns) seconds: centre of the Gaussian relative to the predicted time.
    Returns the chosen samples and their shifts delta (s) relative to the continuous predicted time."""
    ns = Eh.shape[0]
    offs = np.arange(-W, W + 1)
    ii = np.clip(IDX[..., None] + offs, 0, NT - 1)
    val = Eh[np.arange(ns).reshape((1,) * (IDX.ndim - 1) + (ns, 1)), ii]
    if s2 is not None:
        dl = (offs - frac[..., None]) * DT
        if shift is not None:
            dl = dl - shift[..., None]
        val = val - dl * dl / (2.0 * np.asarray(s2)[..., None])
    a = val.argmax(axis=-1)
    tpk = np.take_along_axis(ii, a[..., None], axis=-1)[..., 0]
    return tpk, (a - W - frac) * DT


def _cache(E, pr, mode, W):
    """static arrays of the whole lattice of an engine: flat gather index, shift table dl, charge table."""
    key = ("flat", W)
    if key not in E.cache:
        offs = np.arange(-W, W + 1)
        ii = np.clip(E.IDX[..., None] + offs, 0, NT - 1)
        E.cache[key] = (ii + (np.arange(E.ns) * NT)[None, None, None, :, None],
                        ((offs - E.frac[..., None]) * DT).astype(np.float32))
    flat, dl = E.cache[key]
    ck = ("pen", mode, W)
    if ck not in pr:
        s2 = pr["s2"][:, None] if mode in ("indep", "corr") else np.broadcast_to(pr["sc2"][None, None], E.IDX.shape)
        pr[ck] = (dl * dl / (2.0 * s2[..., None])).astype(np.float32)
    return flat, dl, pr[ck]


def choose(E, Eh, mode, pr, W, sel=None, tsel=None, R=None):
    """peak samples and timing charge of the blocks of engine E (optionally nodes sel and origins tsel).
    mode: hard (registered rule: peak of the UNCAPPED energy of R) | indep (pr['s2']) | indepc (pr['sc2'], no cell
    term) | corr."""
    if mode == "hard":
        Kt = D2.peak_table(R) if W == SH else xc.peak_table_w(R, W)
        IDX = E.IDX if sel is None else E.IDX[sel]
        IDX = IDX if tsel is None else IDX[:, tsel]
        return Kt[E.st4, IDX], None
    if sel is not None or tsel is not None:
        assert mode == "indep"
        IDX, frac = E.IDX[sel][:, tsel], E.frac[sel][:, tsel]
        s2 = pr["s2"][sel][:, None]
        tpk, dl = pick(Eh, IDX, frac, s2, W)
        return tpk, (dl * dl / (2.0 * s2)).sum((-1, -2))
    flat, dl, pen = _cache(E, pr, mode, W)
    Ef = Eh.astype(np.float32).ravel()
    g = Ef[flat]
    a = (g - pen).argmax(axis=-1)
    if mode == "corr":
        N, T_ = E.N, E.T
        d1 = (a - W - E.frac) * DT
        th = np.einsum("nas,nts->nta", pr["Tm"], d1.reshape(N, T_, -1))
        sh = np.einsum("nra,nta->ntr", E.J, th).reshape(d1.shape).astype(np.float32)
        w2 = (0.5 / pr["sc2"]).astype(np.float32)[None, None, :, :, None]
        a = (g - (dl - sh[..., None]) ** 2 * w2).argmax(axis=-1)
        d = ((a - W - E.frac) * DT).reshape(N, T_, -1)
        tpk = np.clip(E.IDX + a - W, 0, NT - 1)
        return tpk, 0.5 * np.einsum("ntr,nrs,nts->nt", d, pr["Cinv"], d)
    tpk = np.clip(E.IDX + a - W, 0, NT - 1)
    return tpk, np.take_along_axis(pen, a[..., None], axis=-1)[..., 0].sum((-1, -2)).astype(float)


def _half_energy(R, capv):
    Et = np.sum(R * R, axis=1)
    return (Et if capv is None else np.minimum(Et, capv * capv)) / 2.0


def _vectors(R, tpk, capv):
    """(..., 2, ns, 3) residual vectors at samples tpk (..., 2, ns), capped."""
    ns = R.shape[0]
    V = np.ascontiguousarray(R.transpose(0, 2, 1))
    v = V[np.arange(ns).reshape((1,) * (tpk.ndim - 1) + (ns,)), tpk]
    if capv is not None:
        nv = np.linalg.norm(v, axis=-1, keepdims=True)
        v = v * np.minimum(1.0, capv / np.maximum(nv, 1e-12))
    return v


def _deflate(R, pred, tp):
    ns = R.shape[0]
    for b in range(2):
        for i in range(ns):
            t0 = int(tp[b, i])
            lo, hi = max(0, t0 - X.AC_HALF), min(NT, t0 + X.AC_HALF + 1)
            R[i, :, lo:hi] -= np.outer(pred[b, i], X.AC[X.AC_HALF - (t0 - lo):X.AC_HALF + (hi - t0)])


def _ref_sets(EF, src, tau, prev):
    """fine nodes within R_REF km of src and farther than 5 km from every earlier atom; origins within T_REF s."""
    dxy = np.linalg.norm(EF.XN[:, :2] - np.asarray(src)[None, :2], axis=1)
    ok = dxy <= R_REF
    for q in prev:
        ok &= np.linalg.norm(EF.XN[:, :2] - np.asarray(q)[None, :2], axis=1) > 5.0
    sel = np.where(ok)[0]
    tsel = np.where(np.abs(EF.taus - tau) <= T_REF + 1e-9)[0]
    return sel, tsel


# ============================================================================================ MT-FW family
def mt_loop(E, MFn, kmax, mode="hard", pr=None, W=SH, EF=None, prF=None, capv=CAPV, min_sep=5.0):
    """c6_lib.mtfw_dense / c9_lib.mtfw_c9 loop (identical for mode hard, W = SH, EF None) with the timing charge.
    EF: fine engine -> refine (selection by the hard rule on E, relocation by the indep objective on EF with prF).
    Atoms: dicts x (3,), tau, s (selection statistic, charge included), s0 (amplitude part), sref (refined objective)."""
    ns = E.ns
    R = MFn.copy()
    excluded = np.zeros(E.N, bool)
    atoms, prev = [], []
    for _ in range(kmax):
        if excluded.all():
            break
        Eh = _half_energy(R, capv)
        tpk, pen = choose(E, Eh, mode, pr, W, R=R)
        vf = _vectors(R, tpk, capv).reshape(E.N, E.T, 6 * ns)
        y = np.einsum("gtr,grk->gtk", vf, E.Q)
        s0 = (y * y).sum(-1) / 2.0
        sc = s0 if pen is None else s0 - pen
        sc[excluded] = -np.inf
        gi = int(np.argmax(sc.max(axis=1)))
        tj = int(np.argmax(sc[gi]))
        s = float(sc[gi, tj])
        if not np.isfinite(s):
            break
        src, tau, tp, vsel, G = E.XN[gi], float(E.taus[tj]), tpk[gi, tj].copy(), vf[gi, tj], E.G[gi]
        a = {"s": s, "s0": float(s0[gi, tj])}
        if EF is not None:
            sel, tsel = _ref_sets(EF, src, tau, prev)
            if len(sel):
                tpf, penf = choose(EF, Eh, "indep", prF, W, sel=sel, tsel=tsel)
                vff = _vectors(R, tpf, capv).reshape(len(sel), len(tsel), 6 * ns)
                yf = np.einsum("gtr,grk->gtk", vff, EF.Q[sel])
                of = (yf * yf).sum(-1) / 2.0 - penf
                k0 = np.where(np.all(np.abs(EF.XN[sel] - src[None]) < 1e-6, axis=1))[0]
                t0 = np.where(np.abs(EF.taus[tsel] - tau) < 1e-9)[0]
                o0 = float(of[k0[0], t0[0]]) if len(k0) and len(t0) else -np.inf
                kf = int(np.argmax(of.max(axis=1)))
                tf = int(np.argmax(of[kf]))
                a["sref"] = float(of[kf, tf])
                if of[kf, tf] > o0:
                    src, tau = EF.XN[sel[kf]], float(EF.taus[tsel[tf]])
                    tp, vsel, G = tpf[kf, tf].copy(), vff[kf, tf], EF.G[sel[kf]]
        a["x"], a["tau"] = np.asarray(src, float), tau
        atoms.append(a)
        prev.append(np.asarray(src, float))
        excluded = excluded | (np.linalg.norm(E.XN[:, :2] - src[None, :2], axis=1) <= min_sep)
        m, *_ = np.linalg.lstsq(G, vsel, rcond=None)
        _deflate(R, (G @ m).reshape(2, ns, 3), tp)
    return atoms


# ============================================================================================ BOMP-PK family
def bomp_energy(E, R, tpk, sel=None):
    """e = h^T K(I)^+ h of the blocks (task2d_methods.PKDict.energy with explicit operators)."""
    G = E.G if sel is None else E.G[sel]
    GtG = E.GtG if sel is None else E.GtG[sel]
    Cps = E.Cps if sel is None else E.Cps[sel]
    Kinv0 = E.Kinv0 if sel is None else E.Kinv0[sel]
    N, T_ = tpk.shape[:2]
    v = _vectors(R, tpk, None).reshape(N, T_, -1)
    h = np.matmul(v, G)
    lag = tpk[:, :, 1, :] - tpk[:, :, 0, :]
    wm, wp = M.ac_at(-lag), M.ac_at(lag)
    diff = np.any((wm != 0) | (wp != 0), axis=-1)
    Kinv = np.broadcast_to(Kinv0[:, None], (N, T_, 6, 6))
    if diff.any():
        Kinv = Kinv.copy()
        n_, t_ = np.nonzero(diff)
        K = GtG[n_] + np.einsum("mi,mikl->mkl", wm[n_, t_], Cps[n_]) + np.einsum("mi,milk->mkl", wp[n_, t_], Cps[n_])
        Kinv[n_, t_] = D2.PKDict.pinv_sym(0.5 * (K + np.swapaxes(K, -1, -2)))
    return np.einsum("ntk,ntk->nt", h, np.matmul(Kinv, h[..., None])[..., 0])


def _gram(Gb, I):
    GA, IA = np.stack(Gb), np.stack(I)
    Wt = M.ac_at(IA[:, None, :, None, :] - IA[None, :, None, :, :])
    n = len(Gb)
    K = np.einsum("apick,bqicl,abpqi->akbl", GA, GA, Wt, optimize=True).reshape(6 * n, 6 * n)
    return 0.5 * (K + K.T)


def bomp_loop(E, PD, MFn, kmax, mode="hard", pr=None, W=SH, EF=None, prF=None, min_sep=5.0):
    """task2d_methods.bomp_pk_record loop (identical for mode hard, W = SH, EF None) with the timing charge.
    PD (PKDict of the coarse grid) supplies the block energies when mode is hard and W = SH is not required: the
    energies are always computed by bomp_energy (validated against PD.energy)."""
    ns = E.ns
    stv = np.arange(ns)
    excluded = np.zeros(E.N, bool)
    Gs, Is, atoms, prev = [], [], [], []
    R = MFn.copy()
    for _ in range(kmax):
        if excluded.all():
            break
        Eh = _half_energy(R, None)
        tpk, pen = choose(E, Eh, mode, pr, W, R=R)
        e = bomp_energy(E, R, tpk) / 2.0
        sc = e if pen is None else e - pen
        sc[excluded] = -np.inf
        j = int(np.argmax(sc))
        gi, tj = divmod(j, E.T)
        s = float(sc[gi, tj])
        if not np.isfinite(s):
            break
        src, tau, tp, G = E.XN[gi], float(E.taus[tj]), tpk[gi, tj].copy(), E.G[gi]
        a = {"s": s, "s0": float(e[gi, tj])}
        if EF is not None:
            sel, tsel = _ref_sets(EF, src, tau, prev)
            if len(sel):
                tpf, penf = choose(EF, Eh, "indep", prF, W, sel=sel, tsel=tsel)
                of = bomp_energy(EF, R, tpf, sel=sel) / 2.0 - penf
                k0 = np.where(np.all(np.abs(EF.XN[sel] - src[None]) < 1e-6, axis=1))[0]
                t0 = np.where(np.abs(EF.taus[tsel] - tau) < 1e-9)[0]
                o0 = float(of[k0[0], t0[0]]) if len(k0) and len(t0) else -np.inf
                kf = int(np.argmax(of.max(axis=1)))
                tf = int(np.argmax(of[kf]))
                a["sref"] = float(of[kf, tf])
                if of[kf, tf] > o0:
                    src, tau = EF.XN[sel[kf]], float(EF.taus[tsel[tf]])
                    tp, G = tpf[kf, tf].copy(), EF.G[sel[kf]]
        a["x"], a["tau"] = np.asarray(src, float), tau
        atoms.append(a)
        prev.append(np.asarray(src, float))
        Gs.append(G)
        Is.append(tp)
        excluded = excluded | (np.linalg.norm(E.XN[:, :2] - src[None, :2], axis=1) <= min_sep)
        Kg = _gram([g.reshape(2, ns, 3, 6) for g in Gs], Is)
        V = MFn[stv[None, None, :], :, np.stack(Is)]
        b = np.einsum("bx,bxk->bk", V.reshape(len(Is), -1), np.stack(Gs)).ravel()
        U, *_ = np.linalg.lstsq(Kg, b, rcond=None)
        U = U.reshape(len(Is), 6)
        R = MFn.copy()
        for g, I, u in zip(Gs, Is, U):
            _deflate(R, (g @ u).reshape(2, ns, 3), I)
    return atoms


# ============================================================================================ hit rules
def label(atoms, truth, taus, tol_t=None, tol_km=5.0):
    """registered rule (tol_t None): nearest event within 5 km (task2c_methods.hit_of).  Strict (2I / verif_2f R5t05):
    nearest event with epicentre within 5 km AND |atom origin - true origin| <= tol_t."""
    out = []
    for a in atoms:
        h = -1
        if truth:
            dd = np.array([np.linalg.norm(a["x"][:2] - np.asarray(t)[:2]) for t in truth])
            if tol_t is None:
                j = int(np.argmin(dd))
                h = j if dd[j] <= tol_km else -1
            else:
                ok = (dd <= tol_km) & (np.abs(a["tau"] - np.asarray(taus)) <= tol_t + 1e-9)
                if ok.any():
                    h = int(np.where(ok, dd, np.inf).argmin())
        out.append(h)
    return out


def record(atoms, labels, key="s"):
    rec, run = [], np.inf
    for a, h in zip(atoms, labels):
        run = min(run, a[key])
        rec.append((float(run), int(h)))
    return rec
