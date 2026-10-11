"""
rft_threshold.py - Component A: a threshold for the maximum of the moment-tensor certificate over
the hypothesis space, T(theta) = ||Q(theta)^T v(theta)||^2 / 2 (the WAVE-MT value of fields_x and the
Frank-Wolfe score of fw_record), with prescribed false-alarm probability alpha = P(M > u),
M = max over the 675 nodes x 31 origin times.

Predictors
  eec_threshold      expected Euler characteristic (Worsley 1994; Adler and Taylor 2007) of a chi2_6
                     random field in D = 4 dimensions (x, y, z, tau).  The non-chi-square marginal
                     of T is mapped by a quantile transform T* = F_chi2_6^{-1}(F_T(T)); the
                     Lipschitz-Killing curvatures (LKC) of the search box are computed from the
                     roughness Lambda = Var(grad Z) of the underlying Gaussian components, estimated
                     from the lag correlations of the transformed field (corr(T*(x), T*(x+h)) =
                     rho(h)^2 = exp(-h' Lambda h) for a chi2 field with Gaussian-shaped correlation).
                     Variant 'ecfit': the LKC are fitted to the empirical Euler characteristic curve
                     of the LATTICE excursion sets at moderate thresholds (Taylor and Worsley 2007).
  bonferroni         N (1 - F_T(u)) = alpha (reference).
  gev_threshold      GEV (and Gumbel) fit by maximum likelihood to maxima of calibration windows.
  conformal_threshold  split-conformal: the ceil((n+1)(1-alpha))-th smallest calibration maximum;
                     P(M_new > u) <= alpha under exchangeability of windows.

Functions: fast_field (vectorised WAVE-MT, checked against task1x_mtglrt.fields_x), null_max_samples,
ec_density_chi2, lkc_box, eec_threshold, bonferroni_threshold, gev_threshold, conformal_threshold,
selftest_eec (smooth chi2_6 field in 4-D: EEC formula against simulation), validate.

    python rft_threshold.py [out.txt]   # tables of rft_report.txt (about 6 min, one core)
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys
import time
from itertools import combinations
from math import comb, factorial

import numpy as np
from scipy import optimize, special, stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import task1x_mtglrt as X  # noqa: E402
import task1y_clutter as Y  # noqa: E402
from task0_common import AP, AS, DEFAULT_TAUS, DT, NT, VP, VS, default_grid, make_network, ricker  # noqa: E402
from task0_common import T as T_AX  # noqa: E402
from task1o_consistent import matched  # noqa: E402
from task1r_fair import SCEN  # noqa: E402
from task1w_mech import CP, CS  # noqa: E402

GRID = default_grid()
GARR = np.asarray(GRID)
AXES = [np.unique(GARR[:, k]) for k in range(3)]
LSHAPE = (len(AXES[0]), len(AXES[1]), len(AXES[2]), len(DEFAULT_TAUS))      # (15, 15, 3, 31)
# lattice index of each grid node (default_grid order: x outer, y, z inner)
_GIDX = np.full(LSHAPE[:3], -1, int)
_GIDX[tuple(np.searchsorted(AXES[k], GARR[:, k]) for k in range(3))] = np.arange(len(GARR))
STEP = np.array([AXES[0][1] - AXES[0][0], AXES[1][1] - AXES[1][0], AXES[2][1] - AXES[2][0],
                 DEFAULT_TAUS[1] - DEFAULT_TAUS[0]])                            # km, km, km, s
EXTENT = np.array([AXES[0][-1] - AXES[0][0], AXES[1][-1] - AXES[1][0], AXES[2][-1] - AXES[2][0],
                   DEFAULT_TAUS[-1] - DEFAULT_TAUS[0]])
NU = 6
ALPHAS = (0.01, 0.05, 0.1, 0.25)


# ----------------------------------------------------------------------------------------------
# fast WAVE-MT field
# ----------------------------------------------------------------------------------------------
def radiation_batch(STA, garr=GARR):
    """task1x_mtglrt.radiation_operator for all nodes at once: (n_nodes, 6 ns, 6)."""
    v = STA[None, :, :] - garr[:, None, :]                                     # (G, ns, 3)
    d = np.linalg.norm(v, axis=2)
    n = v / d[..., None]
    nx, ny, nz = n[..., 0], n[..., 1], n[..., 2]
    c = np.stack([nx * nx, ny * ny, nz * nz, 2 * nx * ny, 2 * nx * nz, 2 * ny * nz], -1)   # (G, ns, 6)
    gP = n[..., :, None] * c[..., None, :] * (AP / (d * np.sqrt(CP)))[..., None, None]     # (G, ns, 3, 6)
    z = np.zeros_like(nx)
    B = np.stack([np.stack([nx, z, z, ny, nz, z], -1), np.stack([z, ny, z, nx, z, nz], -1),
                  np.stack([z, z, nz, z, nx, ny], -1)], -2)                    # (G, ns, 3, 6)
    Pp = np.eye(3) - n[..., :, None] * n[..., None, :]
    gS = np.einsum("gsab,gsbk->gsak", Pp, B) * (AS / (d * np.sqrt(CS)))[..., None, None]
    return np.concatenate([gP.reshape(len(garr), -1, 6), gS.reshape(len(garr), -1, 6)], axis=1)


def basis_batch(G):
    U, s, _ = np.linalg.svd(G, full_matrices=False)
    keep = s > 1e-8 * s.max(axis=1, keepdims=True)
    return U * keep[:, None, :]                                                # zero columns if rank < 6


class Geometry:
    """Per-network precomputation: sample index of every (node, tau, phase, station) and bases."""

    def __init__(self, STA, taus=DEFAULT_TAUS):
        self.STA = STA
        self.ns = len(STA)
        d = np.linalg.norm(STA[None] - GARR[:, None], axis=2)                  # (G, ns)
        tt = np.stack([d / VP, d / VS], axis=1)[:, None] + np.asarray(taus)[None, :, None, None]
        self.idx = np.clip((tt / DT).astype(int), 0, NT - 1)                   # (G, nt, 2, ns)
        self.Q = basis_batch(radiation_batch(STA))                             # (G, 6ns, 6)


def fast_field(geo, MFn):
    """T(node, tau) = ||Q^T vt||^2 / 2 for all hypotheses (node_vectors semantics: energy argmax in
    the +-SH window, first occurrence, capped at CAPV).  Returns (n_nodes, n_taus)."""
    ns = geo.ns
    Et = np.sum(np.ascontiguousarray(MFn.transpose(0, 2, 1)) ** 2, axis=-1)   # (ns, NT)
    II = np.clip(np.arange(NT)[:, None] + X.OFFS[None, :], 0, NT - 1)
    K = II[np.arange(NT)[None, :], Et[:, II].argmax(axis=-1)]                  # (ns, NT)
    V = np.ascontiguousarray(MFn.transpose(0, 2, 1))
    nv = np.linalg.norm(V, axis=-1, keepdims=True)
    VC = V * np.minimum(1.0, X.CAPV / np.maximum(nv, 1e-12))
    st = np.arange(ns)[None, None, None, :]
    vf = VC[st, K[st, geo.idx]].reshape(geo.idx.shape[0], geo.idx.shape[1], 6 * ns)
    y = np.matmul(vf, geo.Q)
    return (y * y).sum(-1) / 2.0


def to_lattice(F):
    """(n_nodes, n_taus) -> (nx, ny, nz, nt)."""
    return F[_GIDX]


# ----------------------------------------------------------------------------------------------
# null windows
# ----------------------------------------------------------------------------------------------
def null_max_samples(scen, n_windows, seed, network="fixed", net_seed=None, n_clutter=0, clut_amp=0.35,
                     keep_fields=False, probe=None):
    """Maxima M of T over the hypothesis grid in windows WITHOUT sources.
    network 'fixed': one network make_network(kind, n, net_seed) for all windows (a real network);
            'per_window': make_network(kind, n, seed + k) per window (the task1x/1y protocol).
    n_clutter = 0: pure Gaussian noise (generate_y with no source and no clutter); otherwise
            generate_y with n_clutter pulses of amplitude clut_amp (pert condition; with no source the
            perturbations have no effect).
    Returns dict: M (n_windows,), fields (n_windows, 15, 15, 3, 31) float32 if keep_fields, and the
    values at the probe hypotheses (list of (node, tau_index))."""
    kind, nst, _ = SCEN[scen]
    rng = np.random.default_rng(seed + 41)
    rng_m = np.random.default_rng(seed + 233)
    rng_h = np.random.default_rng(seed + 61)
    geo = Geometry(make_network(kind, nst, seed=net_seed)) if network == "fixed" else None
    M = np.empty(n_windows)
    arg = np.empty((n_windows, 2), int)
    fields = np.empty((n_windows,) + LSHAPE, np.float32) if keep_fields else None
    pv = np.empty((n_windows, len(probe))) if probe else None
    for k in range(n_windows):
        g = geo if geo is not None else Geometry(make_network(kind, nst, seed=seed + k))
        D, _p, _t, _a = Y.generate_y(g.STA, [], [], n_clutter, rng, rng_m, rng_h, clut_amp)
        MF, sd = matched(D)
        F = fast_field(g, MF / sd)
        M[k] = F.max()
        arg[k] = np.unravel_index(int(F.argmax()), F.shape)
        if keep_fields:
            fields[k] = to_lattice(F)
        if probe:
            pv[k] = [F[a, b] for a, b in probe]
    return {"M": M, "argmax": arg, "fields": fields, "probe": pv}


# ----------------------------------------------------------------------------------------------
# Euler-characteristic densities of a chi2 field and LKC of a box
# ----------------------------------------------------------------------------------------------
def ec_density_chi2(j, u, nu=NU):
    """EC density rho_j(u) of a chi2_nu field built from unit-variance Gaussian fields with
    Var(grad Z) = I (Worsley 1994, Thm 3.5; Adler and Taylor 2007, eq. 15.10.x)."""
    u = np.asarray(u, float)
    if j == 0:
        return stats.chi2.sf(u, nu)
    pre = u ** ((nu - j) / 2.0) * np.exp(-u / 2.0) / ((2 * np.pi) ** (j / 2.0) * special.gamma(nu / 2.0)
                                                      * 2 ** (nu / 2.0 - 1))
    s = np.zeros_like(u)
    for l_ in range((j - 1) // 2 + 1):
        for m in range(j - 1 - 2 * l_ + 1):
            if nu >= j - m - 2 * l_:
                s = s + (comb(nu - 1, j - 1 - m - 2 * l_) * (-1) ** (j - 1 + m + l_) * factorial(j - 1)
                         / (factorial(m) * factorial(l_) * 2 ** l_)) * u ** (m + l_)
    return pre * s


def eec(u, L, nu=NU):
    return sum(L[j] * ec_density_chi2(j, u, nu) for j in range(len(L)))


def lkc_box(extent, Lam):
    """LKC L_0..L_D of a box with the given side lengths under the metric Var(grad Z) = Lam:
    L_j = sum_{|J| = j} prod_{a in J} extent_a * sqrt(det Lam_JJ) (a parallelotope in the induced
    Euclidean metric; the external angles of the 2^(D-j) parallel copies of a j-face sum to 1)."""
    D = len(extent)
    L = [1.0]
    for j in range(1, D + 1):
        s = 0.0
        for J in combinations(range(D), j):
            J = list(J)
            s += float(np.prod(np.asarray(extent)[J]) * np.sqrt(max(np.linalg.det(Lam[np.ix_(J, J)]), 0.0)))
        L.append(s)
    return np.array(L)


def lattice_ec(B):
    """Euler characteristic of the cubical complex of a boolean lattice set, for a stack of sets:
    B has shape (n_sets, d1, ..., dD); EC = sum_J (-1)^|J| #{|J|-cells along J with all vertices in B}."""
    D = B.ndim - 1
    ec = np.zeros(B.shape[0])
    for j in range(D + 1):
        for J in combinations(range(D), j):
            A = B
            for a in J:
                ax = a + 1
                s0 = [slice(None)] * A.ndim
                s1 = [slice(None)] * A.ndim
                s0[ax] = slice(0, -1)
                s1[ax] = slice(1, None)
                A = A[tuple(s0)] & A[tuple(s1)]
            ec += (-1) ** j * A.reshape(A.shape[0], -1).sum(1)
    return ec


# ----------------------------------------------------------------------------------------------
# marginal model and roughness
# ----------------------------------------------------------------------------------------------
class Marginal:
    """F_T: empirical CDF of the pooled calibration values for exceedance >= p_tail, and above it a
    scaled chi-square c*chi2_k fitted to the upper tail by matching the empirical 1 - p_tail and
    1 - p_tail/30 quantiles (continuous at the junction)."""

    def __init__(self, values, p_tail=1e-3):
        v = np.sort(np.asarray(values, float).ravel())
        self.v = v
        self.n = len(v)
        m, s2 = v.mean(), v.var()
        self.mom = (s2 / (2 * m), 2 * m * m / s2)                 # T ~ c chi2_k by moments: (c, k)
        self.p_tail = p_tail
        q1, q2 = np.quantile(v, [1 - p_tail, 1 - p_tail / 30.0])
        self.u_tail = q1

        def eqs(k):
            c = q1 / stats.chi2.isf(p_tail, k)
            return c * stats.chi2.isf(p_tail / 30.0, k) - q2
        try:
            k = optimize.brentq(eqs, 0.3, 200.0)
        except ValueError:
            k = self.mom[1]
        self.tail = (q1 / stats.chi2.isf(p_tail, k), k)

    def sf(self, u):
        u = np.asarray(u, float)
        emp = 1.0 - np.searchsorted(self.v, u, side="right") / self.n
        c, k = self.tail
        return np.where(u <= self.u_tail, emp, stats.chi2.sf(u / c, k))

    def isf(self, p):
        p = np.asarray(p, float)
        c, k = self.tail
        return np.where(p >= self.p_tail, np.quantile(self.v, np.clip(1 - p, 0, 1)), c * stats.chi2.isf(p, k))

    def to_chi2(self, T):
        """quantile transform to an equivalent chi2_6 value (ranks; ties broken by value order)."""
        r = (np.searchsorted(self.v, T, side="left") + np.searchsorted(self.v, T, side="right")) / 2.0
        p = np.clip(1.0 - r / self.n, 0.5 / self.n, 1 - 0.5 / self.n)
        return stats.chi2.isf(p, NU)


def lag_correlations(Fs, lags):
    """pooled Pearson correlation of the field at lattice lag h (integer 4-vector), over windows."""
    out = []
    Fs = np.asarray(Fs, float)
    for h in lags:
        s0, s1 = [slice(None)], [slice(None)]
        for a, k in enumerate(h):
            n = Fs.shape[a + 1]
            if k >= 0:
                s0.append(slice(0, n - k))
                s1.append(slice(k, n))
            else:
                s0.append(slice(-k, n))
                s1.append(slice(0, n + k))
        a_ = Fs[tuple(s0)].ravel()
        b_ = Fs[tuple(s1)].ravel()
        out.append(float(np.corrcoef(a_, b_)[0, 1]))
    return np.array(out)


def roughness(Fstar):
    """Lambda (4x4, units km^-2, km^-1 s^-1, s^-2) from the lag correlations of the chi2-equivalent
    field: -ln corr(h) = h' Lambda h (chi2 field, Gaussian-shaped component correlation); diagonal
    from unit lags, off-diagonal from (e_a + e_b) and (e_a - e_b)."""
    E = np.eye(4, dtype=int)
    diag = lag_correlations(Fstar, [E[a] for a in range(4)])
    q = -np.log(np.clip(diag, 1e-6, 1.0))
    Lam = np.diag(q / STEP ** 2)
    for a, b in combinations(range(4), 2):
        cp, cm = lag_correlations(Fstar, [E[a] + E[b], E[a] - E[b]])
        qp, qm = -np.log(np.clip(cp, 1e-6, 1)), -np.log(np.clip(cm, 1e-6, 1))
        Lam[a, b] = Lam[b, a] = (qp - qm) / (4 * STEP[a] * STEP[b])
    w, V = np.linalg.eigh(Lam)
    Lam = (V * np.maximum(w, 0)) @ V.T                                       # nearest PSD
    return Lam, diag


# ----------------------------------------------------------------------------------------------
# threshold predictors
# ----------------------------------------------------------------------------------------------
def _solve_eec(alpha, L):
    f = lambda u: np.log(max(eec(u, L), 1e-300)) - np.log(alpha)  # noqa: E731
    lo, hi = 1.0, 400.0
    grid = np.linspace(lo, hi, 4000)
    val = np.array([eec(u, L) for u in grid])
    # largest u where EEC crosses alpha from above
    idx = np.where(val >= alpha)[0]
    if len(idx) == 0:
        return lo
    i = idx[-1]
    return optimize.brentq(f, grid[i], grid[min(i + 1, len(grid) - 1)])


def eec_threshold(alpha, marginal, L):
    """u with EEC(u*) = alpha for the chi2_6-equivalent field, mapped back through the marginal."""
    us = _solve_eec(alpha, L)
    return float(marginal.isf(stats.chi2.sf(us, NU)))


def eec_time_threshold(alpha, marginal, lam_tau, n_nodes=None):
    """Lattice-aware variant: union bound over the grid nodes (adjacent nodes are nearly
    uncorrelated) of the 1-D EEC (Rice formula) of each node's chi2-equivalent process in tau:
    P(M > u) <= N_nodes [rho_0(u*) + L1_tau rho_1(u*)], L1_tau = extent_tau sqrt(Lambda_tautau)."""
    n_nodes = n_nodes or len(GARR)
    L = np.array([n_nodes, n_nodes * EXTENT[3] * np.sqrt(lam_tau)])
    us = _solve_eec(alpha, L)
    return float(marginal.isf(stats.chi2.sf(us, NU)))


def bonferroni_threshold(alpha, marginal, n_hyp):
    return float(marginal.isf(alpha / n_hyp))


def gev_threshold(M_cal, alpha, family="gev"):
    if family == "gumbel":
        loc, sc = stats.gumbel_r.fit(M_cal)
        return float(stats.gumbel_r.isf(alpha, loc, sc))
    c, loc, sc = stats.genextreme.fit(M_cal)
    return float(stats.genextreme.isf(alpha, c, loc, sc))


def conformal_threshold(M_cal, alpha):
    n = len(M_cal)
    k = int(np.ceil((n + 1) * (1 - alpha)))
    return float(np.sort(M_cal)[k - 1]) if k <= n else np.inf


def ecfit_lkc(Fstar, u_grid, min_mean=0.02):
    """LKC fitted (non-negative least squares) to the mean lattice EC curve of the chi2-equivalent
    fields at thresholds u_grid (moderate levels, where the curve is well estimated)."""
    ecs = np.array([lattice_ec(Fstar > u) for u in u_grid])               # (n_u, n_windows)
    m = ecs.mean(1)
    se = ecs.std(1) / np.sqrt(ecs.shape[1]) + 1e-3
    ok = m > min_mean
    A = np.stack([ec_density_chi2(j, u_grid) for j in range(5)], 1)
    L, _ = optimize.nnls(A[ok] / se[ok, None], m[ok] / se[ok])
    return L, m, se


def clopper_pearson(k, n, a=0.05):
    lo = stats.beta.ppf(a / 2, k, n - k + 1) if k > 0 else 0.0
    hi = stats.beta.ppf(1 - a / 2, k + 1, n - k) if k < n else 1.0
    return lo, hi


# ----------------------------------------------------------------------------------------------
# self-test of the EEC formula on a smooth chi2_6 field in 4-D
# ----------------------------------------------------------------------------------------------
def selftest_eec(n_rep=100, n=10, sigma=2.0, seed=7):
    """Gaussian-kernel smoothed white noise (sd sigma voxels): rho(h) = exp(-h^2 / (4 sigma^2)),
    Lambda = I / (2 sigma^2) per voxel.  Compare the EEC of the n^4 box with the observed mean
    lattice EC and with P(max > u)."""
    from scipy.ndimage import gaussian_filter
    rng = np.random.default_rng(seed)
    pad = int(np.ceil(4 * sigma))
    Mx, ecs = [], []
    ug = np.array([8.0, 12.0, 16.0, 20.0, 24.0, 28.0])
    for _ in range(n_rep):
        Z = []
        for _c in range(NU):
            w = rng.normal(size=(n + 2 * pad,) * 4)
            z = gaussian_filter(w, sigma, mode="wrap")[pad:pad + n, pad:pad + n, pad:pad + n, pad:pad + n]
            Z.append(z)
        Z = np.array(Z)
        Z /= np.sqrt(1.0 / ((2 * np.sqrt(np.pi) * sigma) ** 4))                # unit variance (theory)
        Tf = (Z ** 2).sum(0)
        Mx.append(Tf.max())
        ecs.append(lattice_ec(np.stack([Tf > u for u in ug])))
    Mx, ecs = np.array(Mx), np.array(ecs)
    L = lkc_box(np.full(4, n - 1.0), np.eye(4) / (2 * sigma ** 2))
    rows = []
    for k, u in enumerate(ug):
        rows.append((u, float(eec(u, L)), float(ecs[:, k].mean()), float(ecs[:, k].std() / np.sqrt(n_rep)),
                     float((Mx > u).mean())))
    return L, rows


# ----------------------------------------------------------------------------------------------
# validation
# ----------------------------------------------------------------------------------------------
def fit_eec_models(cal_fields):
    """marginal, chi2-equivalent fields, roughness LKC and EC-fit LKC from calibration fields."""
    marg = Marginal(cal_fields)
    Fstar = np.stack([marg.to_chi2(F) for F in cal_fields]).reshape(cal_fields.shape)
    Lam, diag = roughness(Fstar)
    L_r = lkc_box(EXTENT, Lam)
    L_e, mec, se = ecfit_lkc(Fstar, np.linspace(10.0, 30.0, 21))
    return marg, Lam, diag, L_r, L_e


def predict_all(M_cal, models, alphas=ALPHAS):
    marg, _Lam, _diag, L_r, L_e = models
    n_hyp = int(np.prod(LSHAPE))
    out = {}
    for a in alphas:
        out[("EEC-rough", a)] = eec_threshold(a, marg, L_r)
        out[("EEC-time", a)] = eec_time_threshold(a, marg, _Lam[3, 3])
        out[("EEC-ECfit", a)] = eec_threshold(a, marg, L_e)
        out[("Bonferroni", a)] = bonferroni_threshold(a, marg, n_hyp)
        out[("GEV", a)] = gev_threshold(M_cal, a)
        out[("Gumbel", a)] = gev_threshold(M_cal, a, "gumbel")
        out[("conformal", a)] = conformal_threshold(M_cal, a)
    return out


METHODS = ("EEC-rough", "EEC-time", "EEC-ECfit", "Bonferroni", "GEV", "Gumbel", "conformal")


def validate(M_test, thr, alphas=ALPHAS, methods=METHODS):
    """realized exceedance rate with Clopper-Pearson 95% interval for every (method, alpha)."""
    n = len(M_test)
    res = {}
    for m in methods:
        for a in alphas:
            u = thr[(m, a)]
            k = int((M_test > u).sum())
            lo, hi = clopper_pearson(k, n)
            res[(m, a)] = (u, k / n, lo, hi, lo <= a <= hi)
    return res


def fmt_table(res, title, alphas=ALPHAS, methods=None):
    methods = methods or [m for m in METHODS if any(k[0] == m for k in res)]
    L = [title, f"{'predictor':12s} " + " ".join(f"{'alpha=' + str(a):>31s}" for a in alphas)]
    L.append(f"{'':12s} " + " ".join(f"{'u   realized [95% CI]':>31s}" for _ in alphas))
    for m in methods:
        cells = []
        for a in alphas:
            if (m, a) not in res:
                cells.append(" " * 31)
                continue
            u, r, lo, hi, ok = res[(m, a)]
            cells.append(f"{u:6.2f} {r:5.3f} [{lo:5.3f},{hi:5.3f}]{'  ' if ok else ' *'}")
        L.append(f"{m:12s} " + " ".join(f"{c:>31s}" for c in cells))
    L.append("('*' = nominal alpha outside the 95% Clopper-Pearson interval of the realized rate)")
    return L


def check_fast_field(scen="ladder_6", n=3, seed=4242):
    """max |fast_field - fields_x WAVE-MT| over nodes (fields_x keeps the max over taus)."""
    kind, nst, _ = SCEN[scen]
    rng = np.random.default_rng(seed)
    worst = 0.0
    for k in range(n):
        S = make_network(kind, nst, seed=seed + k)
        D = rng.normal(scale=0.14, size=(nst, 3, NT))
        D[k % nst] += np.outer([0.6, -0.3, 0.5], ricker(T_AX - 9.3 - k))     # strong pulse: exercises the cap
        Qs = [X.basis(X.radiation_operator(g, S)) for g in GARR]
        f, MFn = X.fields_x(S, D, GRID, DEFAULT_TAUS, np.arange(len(GARR)), Qs)
        F = fast_field(Geometry(S), MFn)
        worst = max(worst, float(np.max(np.abs(F.max(1) - f["WAVE-MT"]))))
    return worst


def main():
    t0 = time.time()
    NCAL = int(os.environ.get("RFT_NCAL", "300"))
    NTEST = int(os.environ.get("RFT_NTEST", "300"))
    lines = []
    P = lambda s="": (print(s, flush=True), lines.append(s))  # noqa: E731
    P("RANDOM-FIELD THRESHOLD FOR THE MAXIMUM OF THE MOMENT-TENSOR CERTIFICATE (Component A) - rft_threshold.py")
    P("=" * 118)
    P(f"Statistic T(theta) = ||Q^T v||^2 / 2 (fields_x WAVE-MT, fw_record score); M = max over {len(GARR)} nodes x "
      f"{len(DEFAULT_TAUS)} taus = {int(np.prod(LSHAPE))} hypotheses.")
    P(f"fast_field vs task1x_mtglrt.fields_x (3 windows): max |difference| = {check_fast_field():.2e}")
    P("")
    # ---- self-test of the EC densities
    L_st, rows = selftest_eec()
    P("0. SELF-TEST of the chi2_6 EC densities and the LKC of a box: smooth chi2_6 field on a 10^4 lattice")
    P(f"   (Gaussian kernel sd 2.0 voxels, 100 realisations; LKC = {np.round(L_st, 1).tolist()})")
    P(f"   {'u':>5s} {'EEC':>8s} {'mean lattice EC':>18s} {'P(max>u)':>9s}")
    for u, e, m, se, p in rows:
        P(f"   {u:5.1f} {e:8.3f} {m:10.3f} +- {se:5.3f} {p:9.3f}")
    P("")
    scen_names = {"ladder_4": "random, 4 stations", "ladder_6": "random, 6 stations", "unilateral_c40": "one-sided, 8 stations"}
    NET_SEED = {"ladder_4": 999 + 25000, "ladder_6": 999 + 25000, "unilateral_c40": 999 + 25000}
    summary = []
    for scen in ("ladder_4", "ladder_6", "unilateral_c40"):
        kind, nst, _ = SCEN[scen]
        for network in ("fixed", "per_window"):
            P("-" * 118)
            P(f"SCENARIO {scen} ({scen_names[scen]}), network {network}"
              + (f" (make_network('{kind}', {nst}, seed={NET_SEED[scen]}))" if network == "fixed" else
                 " (a new network per window, as in task1x/1y)"))
            probe = [(int(_GIDX[7, 7, 1]), 15), (int(_GIDX[2, 12, 0]), 5), (int(_GIDX[12, 3, 2]), 25)]
            cal = null_max_samples(scen, NCAL, 100000, network, NET_SEED[scen], keep_fields=True, probe=probe)
            tst = null_max_samples(scen, NTEST, 200000, network, NET_SEED[scen], probe=probe)
            # (i) marginal at single hypotheses and pooled
            marg, Lam, diag, L_r, L_e = models = fit_eec_models(cal["fields"])
            pv = np.concatenate([cal["probe"], tst["probe"]])
            P("(i) marginal of 2T (chi2_6 if v were white Gaussian without selection):")
            P(f"    {'set':26s} {'mean':>6s} {'var':>7s} {'q.5':>6s} {'q.9':>6s} {'q.99':>6s} {'q.999':>6s} "
              f"{'KS vs chi2_6':>12s} {'moment fit c*chi2_k':>20s} {'KS vs fit':>9s}")
            chi_q = stats.chi2.ppf([0.5, 0.9, 0.99, 0.999], NU)
            P(f"    {'chi2_6 reference':26s} {6.0:6.2f} {12.0:7.2f} " + " ".join(f"{q:6.2f}" for q in chi_q))
            sets = [(f"hypothesis {j + 1} (n={len(pv)})", 2 * pv[:, j]) for j in range(pv.shape[1])]
            pool = 2 * cal["fields"].ravel().astype(float)
            sub = pool[np.random.default_rng(0).choice(len(pool), 200000, replace=False)]
            sets.append((f"pooled ({len(pool)} values)", sub))
            for name, x in sets:
                m, v = x.mean(), x.var()
                c, k = v / (2 * m), 2 * m * m / v
                ks0 = stats.kstest(x, "chi2", args=(NU,)).statistic
                ks1 = stats.kstest(x / c, "chi2", args=(k,)).statistic
                P(f"    {name:26s} {m:6.2f} {v:7.2f} " + " ".join(f"{q:6.2f}" for q in np.quantile(x, [.5, .9, .99, .999]))
                  + f" {ks0:12.3f} {'c=%.2f k=%.2f' % (c, k):>20s} {ks1:9.3f}")
            ct, kt = marg.tail
            P(f"    tail model of T above its 0.999 quantile ({marg.u_tail:.2f}): T ~ {ct:.3f} chi2_{kt:.2f}; "
              f"P(T > 20) pooled = {marg.sf(20.0):.2e}")
            # (ii) roughness and LKC
            P("(ii) EEC: lag-1 correlation of the chi2-equivalent field along x, y, z, tau: "
              + ", ".join(f"{c:.3f}" for c in diag))
            P("     Lambda (km^-2, km^-1 s^-1, s^-2) diagonal: " + ", ".join(f"{x:.3g}" for x in np.diag(Lam))
              + "; correlation of the off-diagonals: "
              + ", ".join(f"{Lam[a, b] / np.sqrt(Lam[a, a] * Lam[b, b] + 1e-300):+.2f}" for a, b in combinations(range(4), 2)))
            P("     LKC from roughness L0..L4: " + ", ".join(f"{x:.3g}" for x in L_r))
            P("     LKC fitted to the lattice EC curve: " + ", ".join(f"{x:.3g}" for x in L_e))
            thr = predict_all(cal["M"], models)
            res = validate(tst["M"], thr)
            uc = thr[("conformal", 0.05)]
            P(f"     effective number of independent hypotheses at alpha = 0.05: alpha / P(T > u_conformal) = "
              f"{0.05 / max(float(marg.sf(uc)), 1e-300):.0f} (lattice {int(np.prod(LSHAPE))}, nodes {len(GARR)}); "
              f"continuum EEC(u*)/P(T* > u*) at that u: {eec(stats.chi2.isf(float(marg.sf(uc)), NU), L_r) / float(marg.sf(uc)):.0f}")
            P(f"     calibration maxima: mean {cal['M'].mean():.2f}, sd {cal['M'].std():.2f}; held-out: mean "
              f"{tst['M'].mean():.2f}, sd {tst['M'].std():.2f} (n = {NCAL} / {NTEST})")
            P("")
            for s_ in fmt_table(res, f"(iii) Gaussian noise: predicted threshold u and realized P(M > u) on {NTEST} held-out windows"):
                P("   " + s_)
            summary.append((scen, network, "gauss", res))
            # clutter
            if network == "fixed":
                for amp in (0.35, 1.75):
                    ccal = null_max_samples(scen, NCAL, 300000 + int(amp * 100), network, NET_SEED[scen],
                                            n_clutter=35, clut_amp=amp, keep_fields=True)
                    ctst = null_max_samples(scen, NTEST, 400000 + int(amp * 100), network, NET_SEED[scen],
                                            n_clutter=35, clut_amp=amp)
                    res_g = validate(ctst["M"], thr)                       # Gaussian-calibrated thresholds
                    cmodels = fit_eec_models(ccal["fields"])
                    thr_c = predict_all(ccal["M"], cmodels)
                    res_c = validate(ctst["M"], thr_c)
                    P("")
                    P(f"   CLUTTER amplitude {amp} (35 pulses, no source): maxima mean {ctst['M'].mean():.2f} sd "
                      f"{ctst['M'].std():.2f} (Gaussian {tst['M'].mean():.2f} / {tst['M'].std():.2f}); "
                      f"clutter-cal LKC roughness {', '.join(f'{x:.3g}' for x in cmodels[3])}")
                    for s_ in fmt_table(res_g, "   thresholds calibrated on GAUSSIAN noise, applied to clutter windows:"):
                        P("   " + s_)
                    for s_ in fmt_table(res_c, "   thresholds re-calibrated on separate CLUTTER windows (same clutter law):"):
                        P("   " + s_)
                    summary.append((scen, network, f"clut{amp}-gausscal", res_g))
                    summary.append((scen, network, f"clut{amp}-recal", res_c))
            P("")
    # ---- summary: calibration score per method
    P("=" * 118)
    P("SUMMARY: number of (cell, alpha) combinations where the nominal alpha lies in the 95% interval of the realized")
    P("rate, and mean realized/nominal ratio (1 = calibrated; < 1 conservative; > 1 anti-conservative)")
    meths = METHODS
    for tag in ("gauss", "clut0.35-gausscal", "clut0.35-recal", "clut1.75-gausscal", "clut1.75-recal"):
        rr = [r for (_s, _n, t, r) in summary if t == tag]
        if not rr:
            continue
        cells = []
        for m in meths:
            ok = sum(r[(m, a)][4] for r in rr for a in ALPHAS)
            ratio = np.mean([r[(m, a)][1] / a for r in rr for a in ALPHAS])
            cells.append(f"{m} {ok}/{len(rr) * len(ALPHAS)} ({ratio:.2f})")
        P(f"{tag:18s} " + " | ".join(cells))
    P(f"(run time {time.time() - t0:.0f} s)")
    return lines


if __name__ == "__main__":
    out = main()
    if len(sys.argv) > 1:                       # optional: write the tables to the given path
        with open(sys.argv[1], "w") as fh:
            fh.write("\n".join(out) + "\n")
