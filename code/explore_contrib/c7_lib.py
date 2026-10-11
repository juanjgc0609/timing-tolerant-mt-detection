"""
c7_lib.py - C7: source-pulse-invariant moment-tensor subspace detectors (exploratory, SEED 999 ONLY).

Idea.  The registered detectors correlate the data with one fixed 2 Hz Ricker (matched filter, MF) and model each
block (node, origin) as 6 moment-tensor columns G[:, k] (x) Ricker(t - t_bi).  The generator uses per-event pulses
w = cos(phi) r_fc - sin(phi) H[r_fc], fc ~ U[1.4, 2.8] Hz, phi ~ U[-45, 45] deg.  A constant phase rotation lies
exactly in span{r_fc, H[r_fc]}; the centre-frequency family is close to a low-dimensional subspace.  With an
orthonormal pulse basis psi_1..psi_q the block becomes G (x) [psi_1 .. psi_q] (6q columns); the physical model has
coefficient matrix U = m c^T (6 x q, rank 1: one moment tensor, one source time function shared by all stations).

Whitened data domain (as task2c_methods): y = D sqrt(E_r) / sd, E_r = sum KER^2, sd the window's MF noise estimate.
Filter outputs F_j = y (*) psi_j (cross-correlation, unit-norm psi_j), so F_1 = MFn when psi_1 = KER / sqrt(E_r).
Atom column (k, j) read at sample I_bi: data waveform G_bi[:, k] psi_j(t - I_bi); its F_j' image is
G_bi[:, k] X[j', j](s - I_bi), X[j', j](tau) = sum_u psi_j(u) psi_j'(u - tau); inner products between columns are
G.G X(I_a - I_b).  With q = 1 and psi_1 = KER / sqrt(E_r), X[0, 0] = task1x_mtglrt.AC and everything reduces to
BOMP-PK (task2d_methods.bomp_pk_record).

Detector BOMP-PK-Q (per iteration, same grid / origins / +-4-sample peak search / 5 km exclusion / kmax / hit rule as
BOMP-PK):
  peak search  argmax over the +-4-sample window of the ENVELOPE energy sum_j ||R_j,i(:, k)||^2 (phase-invariant);
  statistic    V (6ns x q) = residual filter outputs at the peak samples; Y = Q^T V (Q = orthonormal basis of G,
               task1x_mtglrt.basis); R1 (rank-1, separable): e = lambda_max(Y^T Y) = max over unit c of
               ||Q^T V c||^2; FULL (relaxed 6q block): e = ||Y||_F^2.  No P/S cross term (as MT-FW; the C2 ablation
               showed this term immaterial).  Score e / 2, running minimum.
  cap          (optional, selection only, as BOMP-PK-CAP of c2_ablate) each station/phase 3q-vector clipped at norm
               CAPV = 6.
  refit        joint least squares of all active atoms with the exact Gram (P/S and inter-atom terms through X):
               FULL: unconstrained 6q coefficients; R1: bilinear u_a = m_a c_a^T by alternating least squares.
  residual     recomputed from F with the fitted atoms.
MT-FW-Q: c6_lib.mtfw_dense (capped leaf vectors, sequential deflation) with the envelope peak search, 3q-vector cap,
  R1 statistic and rank-1 deflation (top singular pair of the capped Y).
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

import gen_realistic as GR  # noqa: E402
import task1x_mtglrt as X  # noqa: E402
import task2c_methods as M  # noqa: E402
import task2d_methods as D2  # noqa: E402
from task0_common import DT, NT  # noqa: E402
from task1o_consistent import _KER  # noqa: E402

E_R = float(np.sum(_KER ** 2))
TB = np.arange(-2.0, 2.0 + 1e-9, DT)            # basis support, 65 samples, centred
FC_GRID = np.linspace(1.4, 2.8, 57)
PHI_GRID = np.radians(np.linspace(-45.0, 45.0, 37))


# ============================================================================================ pulse bases
def ricker_fc(t, fc):
    return GR.source_wavelet(t, fc, 0.0)


def hilbert_ricker(t, fc):
    return -GR.source_wavelet(t, fc, np.pi / 2)          # w(pi/2) = -H[r]


def family(t=TB, fcs=FC_GRID, phis=PHI_GRID):
    """unit-norm (discrete) members of the generator's pulse family, rows."""
    F = np.array([GR.source_wavelet(t, fc, ph) for fc in fcs for ph in phis])
    return F / np.linalg.norm(F, axis=1, keepdims=True)


def basis_q1():
    return (_KER / np.sqrt(E_R))[None, :]


def basis_ph():
    r = ricker_fc(TB, 2.0)
    h = hilbert_ricker(TB, 2.0)
    p1 = r / np.linalg.norm(r)
    h = h - (h @ p1) * p1
    return np.stack([p1, h / np.linalg.norm(h)])


def basis_phf(q):
    _U, s, Vt = np.linalg.svd(family(), full_matrices=False)
    return Vt[:q].copy(), s


def captured(B, F):
    """fraction of each unit-norm row of F captured by the span of the orthonormal rows of B (same support)."""
    L = B.shape[1]
    if F.shape[1] != L:
        raise ValueError
    return np.sum((F @ B.T) ** 2, axis=1)


def pad_to(B, L):
    if B.shape[1] == L:
        return B
    k = (L - B.shape[1]) // 2
    return np.pad(B, ((0, 0), (k, k)))


# ============================================================================================ filters and kernels
class Basis:
    def __init__(self, B, name):
        self.B = np.asarray(B, float)
        self.q, self.L = self.B.shape
        self.name = name
        self.is_q1 = (self.q == 1 and self.L == len(_KER) and np.allclose(self.B[0], _KER / np.sqrt(E_R), atol=0))
        if self.is_q1:
            self.XH = X.AC_HALF
            self.Xk = X.AC[None, None, :].copy()
        else:
            self.XH = self.L - 1
            Xk = np.zeros((self.q, self.q, 2 * self.XH + 1))
            for j in range(self.q):
                for jj in range(self.q):
                    # X[j, jj](tau) = sum_u psi_jj(u) psi_j(u - tau) = convolve(psi_jj, psi_j[::-1])[L-1+tau]
                    Xk[j, jj] = np.convolve(self.B[jj], self.B[j][::-1], mode="full")
            self.Xk = Xk

    def filt(self, D, sd, MFn=None):
        """(q, ns, 3, NT) filter outputs in whitened units."""
        if self.is_q1 and MFn is not None:
            return MFn[None].copy()
        ns = D.shape[0]
        F = np.empty((self.q, ns, 3, D.shape[-1]))
        for j in range(self.q):
            k = self.B[j][::-1]
            for i in range(ns):
                for c in range(3):
                    F[j, i, c] = np.convolve(D[i, c], k, mode="same")
        return F * (np.sqrt(E_R) / sd)

    def xat(self, lag):
        """X[:, :, lag] for integer lag array (...,) -> (q, q, ...), 0 outside the support."""
        lag = np.asarray(lag)
        ok = np.abs(lag) <= self.XH
        v = self.Xk[:, :, np.clip(lag + self.XH, 0, 2 * self.XH)]
        return np.where(ok, v, 0.0)


# ============================================================================================ dictionary
class QDict:
    """per network x grid x corrections x basis.  Uses the BlockDict of the registered detectors (indices, G)."""

    def __init__(self, BD, basis, Qst=None):
        self.BD, self.bs = BD, basis
        self.ns, self.N, self.T = BD.ns, BD.N, BD.T
        self.q = basis.q
        self.st4 = BD.st
        self.stv = np.arange(BD.ns)
        if Qst is None:
            from certified_fw import stack_Q
            Qst = stack_Q([X.basis(G) for G in BD.Gf])
        self.QstT = np.ascontiguousarray(np.swapaxes(Qst, 1, 2))[:, None]       # (N, 1, r, 6ns)
        self.Qst = Qst

    # ------------------------------------------------------------------ peak search on the envelope energy
    def peak_table(self, Rq):
        Et = np.sum(Rq ** 2, axis=(0, 2))                                          # (ns, NT)
        W = Et[:, D2.II_TAB]
        a = W.argmax(axis=-1)
        return D2.II_TAB[np.arange(NT)[None, :], a]

    def read(self, Rq, cap=False):
        """I (N, T, 2, ns) and sampled V (N, T, 6ns, q), optionally 3q-vector capped."""
        I = self.peak_table(Rq)[self.st4, self.BD.idx]
        RT = np.ascontiguousarray(Rq.transpose(1, 3, 2, 0))                        # (ns, NT, 3, q)
        V = RT[self.st4, I]                                                        # (N, T, 2, ns, 3, q)
        if cap:
            nv = np.sqrt(np.sum(V ** 2, axis=(-1, -2), keepdims=True))
            V = V * np.minimum(1.0, X.CAPV / np.maximum(nv, 1e-12))
        return I, V.reshape(self.N, self.T, 6 * self.ns, self.q)

    def stat(self, Rq, mode="R1", cap=False):
        I, V = self.read(Rq, cap)
        Y = np.matmul(self.QstT, V)                                                # (N, T, r, q)
        S = np.matmul(np.swapaxes(Y, -1, -2), Y)                                   # (N, T, q, q)
        if mode == "FULL" or self.q == 1:
            e = np.trace(S, axis1=-2, axis2=-1)
        elif self.q == 2:
            a, d, b = S[..., 0, 0], S[..., 1, 1], S[..., 0, 1]
            e = 0.5 * (a + d) + np.sqrt(0.25 * (a - d) ** 2 + b * b)
        else:
            e = np.linalg.eigvalsh(S)[..., -1]
        return e, I, S

    # ------------------------------------------------------------------ atoms (n, t, I)
    def cross_gram(self, A, B):
        """(na, 6, q, nb, 6, q) exact Gram between atom lists."""
        Gb = self.BD.Gb
        nA = np.array([a[0] for a in A], int)
        nB = np.array([b[0] for b in B], int)
        IA = np.stack([a[2] for a in A])
        IB = np.stack([b[2] for b in B])
        lag = IA[:, None, :, None, :] - IB[None, :, None, :, :]                   # (na, nb, p, p', ns)
        W = self.bs.xat(lag)                                                       # (q, q, na, nb, p, p', ns)
        # entry (a,k,j ; b,l,J) = sum_{p,p',i} G_a[p,i,:,k].G_b[p',i,:,l] X[j, J](I_a[p,i] - I_b[p',i])
        return np.einsum("apick,bqicl,jJabpqi->akjblJ", Gb[nA], Gb[nB], W, optimize=True)

    def corr_atoms(self, Fq, atoms):
        """b[a, k, j] = sum_{p,i} G_a[p,i,:,k] . F_j,i(:, I_a[p,i]) -> (na, 6, q)."""
        n = np.array([a[0] for a in atoms], int)
        I = np.stack([a[2] for a in atoms])                                        # (na, 2, ns)
        RT = Fq.transpose(1, 3, 2, 0)                                              # (ns, NT, 3, q)
        V = RT[self.stv[None, None, :], I]                                         # (na, 2, ns, 3, q)
        return np.einsum("apicj,apick->akj", V, self.BD.Gb[n])

    def synth_sub(self, Rq, n, I, U):
        """Rq -= contribution of atom (n, I) with coefficients U (6, q)."""
        pred = np.einsum("pick,kj->picj", self.BD.Gb[n], U)                         # (2, ns, 3, q)
        XH = self.bs.XH
        for b in range(2):
            for i in range(self.ns):
                t0 = int(I[b, i])
                lo, hi = max(0, t0 - XH), min(NT, t0 + XH + 1)
                Xw = self.bs.Xk[:, :, XH - (t0 - lo):XH + (hi - t0)]               # (q, q, len)
                Rq[:, i, :, lo:hi] -= np.einsum("cJ,jJs->jcs", pred[b, i], Xw)

    def residual(self, Fq, atoms, Us):
        R = Fq.copy()
        for a, U in zip(atoms, Us):
            self.synth_sub(R, a[0], a[2], U)
        return R


def fit_full(QD, Fq, atoms):
    K = QD.cross_gram(atoms, atoms)
    na, q = len(atoms), QD.q
    K = K.reshape(6 * q * na, 6 * q * na)
    K = 0.5 * (K + K.T)
    b = QD.corr_atoms(Fq, atoms).ravel()
    u, *_ = np.linalg.lstsq(K, b, rcond=None)
    return list(u.reshape(na, 6, q))


def fit_rank1(QD, Fq, atoms, C0, iters=60, tol=1e-10):
    """bilinear LS: u_a = m_a c_a^T (||c_a|| = 1), alternating least squares from c_a = C0[a]."""
    na, q = len(atoms), QD.q
    K = QD.cross_gram(atoms, atoms)                                                # (na, 6, q, na, 6, q)
    K = 0.5 * (K + K.transpose(3, 4, 5, 0, 1, 2))
    b = QD.corr_atoms(Fq, atoms)                                                   # (na, 6, q)
    C = np.array([c / np.linalg.norm(c) for c in C0])
    Mv = np.zeros((na, 6))
    prev = np.inf
    for _ in range(iters):
        Km = np.einsum("akjblJ,aj,bJ->akbl", K, C, C).reshape(6 * na, 6 * na)
        bm = np.einsum("akj,aj->ak", b, C).ravel()
        Mv = np.linalg.lstsq(Km, bm, rcond=None)[0].reshape(na, 6)
        Kc = np.einsum("akjblJ,ak,bl->ajbJ", K, Mv, Mv).reshape(q * na, q * na)
        bc = np.einsum("akj,ak->aj", b, Mv).ravel()
        C = np.linalg.lstsq(Kc, bc, rcond=None)[0].reshape(na, q)
        nrm = np.linalg.norm(C, axis=1, keepdims=True)
        nrm = np.where(nrm > 0, nrm, 1.0)
        Mv = Mv * nrm
        C = C / nrm
        U = np.einsum("ak,aj->akj", Mv, C)
        obj = float(np.einsum("akj,akjblJ,blJ->", U, K, U) - 2.0 * np.sum(U * b))
        if abs(prev - obj) <= tol * max(1.0, abs(obj)):
            break
        prev = obj
    return [np.outer(Mv[a], C[a]) for a in range(na)], C


def top_c(S):
    w, V = np.linalg.eigh(S)
    return V[:, -1]


# ============================================================================================ detectors
def bomp_q(QD, Fq, truth, kmax=5, mode="R1", cap=False, min_sep=5.0, diag=None):
    garr = QD.BD.garr
    excluded = np.zeros(QD.N, bool)
    S, Cs, rec, run, locs = [], [], [], np.inf, []
    R = Fq.copy()
    for _ in range(kmax):
        e, I, Sq = QD.stat(R, mode, cap)
        e[excluded] = -np.inf
        j = int(np.argmax(e))
        n, t = divmod(j, QD.T)
        s = float(e[n, t]) / 2.0
        if not np.isfinite(s):
            break
        run = min(run, s)
        src = garr[n]
        rec.append((float(run), int(M.hit_of(src, truth))))
        locs.append(src.tolist() + [float(QD.BD.taus[t])])
        S.append((n, t, I[n, t].copy()))
        Cs.append(top_c(Sq[n, t]))
        excluded = M.excl_mask(excluded, garr, src, min_sep)
        if mode == "FULL" or QD.q == 1:
            Us = fit_full(QD, Fq, S)
        else:
            Us, C = fit_rank1(QD, Fq, S, Cs)
            Cs = list(C)
        R = QD.residual(Fq, S, Us)
    if diag is not None:
        diag["atoms"], diag["U"], diag["C"] = S, Us if S else None, Cs
    return rec, locs


def mtfw_q(QD, IDX, Fq, truth, kmax=5, min_sep=5.0, cap=True):
    """dense MT-FW (c6_lib.mtfw_dense) with q pulse filters: envelope peak search, 3q cap, R1 statistic, rank-1
    sequential deflation.  IDX (N, T, 2, ns) the model's sample indices (= BD.idx for the homogeneous model)."""
    ns = QD.ns
    R = Fq.copy()
    garr = QD.BD.garr
    N, T = IDX.shape[:2]
    excluded = np.zeros(N, bool)
    rec, run, locs = [], np.inf, []
    for _ in range(kmax):
        if excluded.all():
            break
        e, I, Sq = QD.stat(R, "R1", cap)
        sc = e / 2.0
        sc[excluded] = -np.inf
        best = sc.max(axis=1)
        gi = int(np.argmax(best))
        tj = int(np.argmax(sc[gi]))
        s = float(sc[gi, tj])
        if not np.isfinite(s):
            break
        src = garr[gi]
        run = min(run, s)
        rec.append((float(run), int(M.hit_of(src, truth))))
        locs.append(src.tolist() + [float(QD.BD.taus[tj])])
        excluded = excluded | (np.linalg.norm(garr[:, :2] - src[:2], axis=1) <= min_sep)
        _I, V = QD.read(R, cap)
        Vb = V[gi, tj]                                                             # (6ns, q) capped
        Q = QD.Qst[gi]
        Y = Q.T @ Vb
        uu, ss, vt = np.linalg.svd(Y, full_matrices=False)
        P = np.outer(Q @ (uu[:, 0] * ss[0]), vt[0])                                # (6ns, q) rank-1 prediction
        pred = P.reshape(2, ns, 3, QD.q)
        tp = I[gi, tj]
        XH = QD.bs.XH
        for b in range(2):
            for i in range(ns):
                t0 = int(tp[b, i])
                lo, hi = max(0, t0 - XH), min(NT, t0 + XH + 1)
                Xw = QD.bs.Xk[:, :, XH - (t0 - lo):XH + (hi - t0)]
                R[:, i, :, lo:hi] -= np.einsum("cJ,jJs->jcs", pred[b, i], Xw)
    return rec, locs


# ============================================================================================ out-of-family pulses
def gabor(t, f0, th):
    """Gabor pulse, envelope sd-like width 0.6 / f0 (narrower band than the Ricker), carrier phase th."""
    s = 0.6 / f0
    t = np.asarray(t, float)
    return np.exp(-(t / s) ** 2) * np.cos(2 * np.pi * f0 * t + th)


def berlage(t, f0, th):
    """causal Berlage pulse t^2 exp(-2 pi f0 t) sin(2 pi f0 t + th), t >= 0 (onset at the arrival)."""
    t = np.asarray(t, float)
    tp = np.maximum(t, 0.0)
    a = 2 * np.pi * f0
    return np.where(t >= 0, tp * tp * np.exp(-a * tp) * np.sin(a * tp + th), 0.0)


_EN = {}


def _energy(fn, f0, th):
    key = (fn.__name__, round(f0, 12), round(th, 12))
    if key not in _EN:
        tt = np.arange(-6.0, 6.0, 1e-3)
        _EN[key] = float(np.sum(fn(tt, f0, th) ** 2) * 1e-3)
    return _EN[key]


class OOFPulse:
    """context manager: gen_realistic draws fc ~ U[1.4, 2.8] and phi ~ U[-phase_max, phase_max] exactly as before
    (same rng consumption) but the arrival waveform is replaced by an out-of-family pulse fn(t, fc, phi) scaled to
    the energy of the design Ricker of the same fc (so the SNR design is unchanged; coda energy unchanged).
    Use with cfg {'phase_max_deg': 90}."""

    def __init__(self, fn):
        self.fn = fn

    def __enter__(self):
        self._sw = GR.source_wavelet
        fn = self.fn

        def sw(t, fc, phi=0.0):
            return fn(t, fc, phi) * np.sqrt(GR.wavelet_energy(fc) / _energy(fn, fc, phi))
        GR.source_wavelet = sw
        return self

    def __exit__(self, *a):
        GR.source_wavelet = self._sw


OOF = {"GABOR90": (gabor, {"phase_max_deg": 90.0}), "BERLAGE90": (berlage, {"phase_max_deg": 90.0})}
