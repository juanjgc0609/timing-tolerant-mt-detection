"""
c12_lib.py - C12: double-couple manifold statistic inside MT-FW (cap 6, registered hard +-0.25 s peak rule, +SC terms).
SEED 999 ONLY.

At a candidate (node g, origin t) the detector reads the capped vectors v (6 ns) at the energy peaks and uses
    FULL6   s6  = |Q_g^T v|^2 / 2                    (energy explained by the 6-D moment-tensor subspace; registered)
    DEV5    s5  = max over trace-free m of the explained energy = s6 - (n_g . y)^2 / 2,   y = Q_g^T v,
            n_g = unit vector of col(Rm_g) orthogonal to Rm_g {trace m = 0},  Rm_g = Q_g^T G_g   (G-metric projection)
    DC      sdc = max over unit double couples m of (m^T G^T v)^2 / (2 m^T G^T G m) = max (Rm m . y)^2 / (2 |Rm m|^2)
The DC manifold is the SO(3) orbit of M0 = diag(1, 0, -1) / sqrt(2) (strike, dip, rake = a rotation); the maximum is
found from 4 starts (the eigen-decomposition closed form of the LS estimate and the 3 best points of a fixed set of
NGRID Haar rotations) followed by NGN damped Gauss-Newton steps on the rotation vector.  sdc <= s5 <= s6 always, so
the lattice maximum of sdc is found EXACTLY (up to the local optimiser) by evaluating candidates in decreasing s5 until
the best sdc found exceeds the next s5.

Moment-tensor 6-vector convention of the operators (c9_lib.operators_hom): (Mxx, Myy, Mzz, Mxy, Mxz, Myz).

Tube formula (Hotelling-Weyl / Gaussian kinematic formula) for z ~ N(0, I_6) and the closed 3-manifold
U = {Rm m / |Rm m|} of S^5:  P(2 sdc > t^2) = P(max_U u.z > t) ~ L3 rho3(t) + L1 rho1(t),
rho3 = (2 pi)^-2 (t^2 - 1) exp(-t^2 / 2), rho1 = (2 pi)^-1 exp(-t^2 / 2), L3 = vol(U), L1 = (1 / 4 pi) int scalar
curvature (L2 = 0: no boundary; L0 = Euler characteristic = 0).  lk_curvatures computes L3 and L1 by Monte Carlo over
Haar rotations with the Gauss equation S = 6 + |H|^2 - |II|^2.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys

import numpy as np
from scipy import stats
from scipy.spatial.transform import Rotation

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
for p in (HERE, WORK):
    if p not in sys.path:
        sys.path.insert(0, p)

import c8_lib as C8  # noqa: E402
import c10_lib as C  # noqa: E402
import gen_realistic as GR  # noqa: E402
import xc  # noqa: E402

SEED = 999
NGRID = 3000
NGN = 8
NSTART = 3
CHUNK = 32
GRID = xc.make_grid()
GARR = np.asarray(GRID)
TAUS = np.asarray(xc.DEFAULT_TAUS, float)

# ============================================================================================ tensors
I6 = (np.array([0, 1, 2, 0, 0, 1]), np.array([0, 1, 2, 1, 2, 2]))
TRACE = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
M0 = np.diag([1.0, 0.0, -1.0]) / np.sqrt(2.0)
WGEN = np.zeros((3, 3, 3))
WGEN[0, 1, 2], WGEN[0, 2, 1] = -1.0, 1.0
WGEN[1, 0, 2], WGEN[1, 2, 0] = 1.0, -1.0
WGEN[2, 0, 1], WGEN[2, 1, 0] = -1.0, 1.0
ISO = np.diag([1.0, 1.0, 1.0, np.sqrt(2.0), np.sqrt(2.0), np.sqrt(2.0)])     # Frobenius metric in the 6-vector


def vec6(Mt):
    return Mt[..., I6[0], I6[1]]


def mat3(m):
    Mt = np.zeros(m.shape[:-1] + (3, 3))
    Mt[..., I6[0], I6[1]] = m
    Mt[..., I6[1], I6[0]] = m
    return Mt


def dc_of_rot(Rr):
    return np.einsum("...ia,ab,...jb->...ij", Rr, M0, Rr)


ROTG = Rotation.random(NGRID, random_state=20261005).as_matrix()
DCG = vec6(dc_of_rot(ROTG))                                                    # (NGRID, 6)


def _comm(Wk, Mt):
    return np.matmul(Wk, Mt) - np.matmul(Mt, Wk)


def _expm(d):
    return Rotation.from_rotvec(d).as_matrix()


# ============================================================================================ node tables
class Tab:
    """Rm = Q^T G (N, 6, 6), its pseudo-inverse, and the unit normal of the trace-free subspace in Q coordinates."""

    def __init__(self, E):
        self.Rm = np.einsum("nrk,nrl->nkl", E.Q, E.G)
        self.Ri = np.linalg.pinv(self.Rm)
        n = np.einsum("nlk,l->nk", self.Ri, TRACE)                             # Rm^-T trace
        self.nvec = n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-300)


def s5_of(y, nvec):
    return ((y * y).sum(-1) - (y * nvec).sum(-1) ** 2) / 2.0


def _fval(Rm, y, Rr):
    p = np.einsum("bkl,bl->bk", Rm, vec6(dc_of_rot(Rr)))
    pp = np.maximum((p * p).sum(-1), 1e-300)
    py = (p * y).sum(-1)
    return py * py / pp, p, py / pp


def dc_max(Rm, Ri, y, ngn=NGN, nstart=NSTART, full=False):
    """Rm, Ri (K, 6, 6), y (K, 6).  Returns sdc (K,) [and the rotation and the value of the closed-form start]."""
    K = len(y)
    # closed form: deviatoric part of the LS estimate, eigenvalues set to (+a, 0, -a)
    mh = np.einsum("bkl,bl->bk", Ri, y)
    Mh = mat3(mh)
    Mh = Mh - np.trace(Mh, axis1=-2, axis2=-1)[:, None, None] * np.eye(3) / 3.0
    _w, V = np.linalg.eigh(Mh)
    R0 = np.stack([V[:, :, 2], V[:, :, 1], V[:, :, 0]], axis=-1)
    R0[:, :, 1] *= np.sign(np.linalg.det(R0))[:, None]
    # grid starts
    Pg = np.einsum("bkl,gl->bkg", Rm, DCG)
    fg = np.einsum("bk,bkg->bg", y, Pg) ** 2 / np.maximum((Pg * Pg).sum(1), 1e-300)
    top = np.argpartition(-fg, nstart, axis=1)[:, :nstart]
    Rr = np.concatenate([R0[:, None], ROTG[top]], axis=1).reshape(-1, 3, 3)     # (K (1 + nstart), 3, 3)
    S = 1 + nstart
    RmB, yB = np.repeat(Rm, S, axis=0), np.repeat(y, S, axis=0)
    f, p, a = _fval(RmB, yB, Rr)
    f_closed = f.reshape(K, S)[:, 0].copy()
    lam = np.full(len(f), 1e-3)
    for _ in range(ngn):
        Mt = dc_of_rot(Rr)
        tk = np.stack([vec6(_comm(WGEN[k], Mt)) for k in range(3)], axis=-1)   # (B, 6, 3)
        J = np.concatenate([p[:, :, None], a[:, None, None] * np.matmul(RmB, tk)], axis=-1)   # (B, 6, 4)
        r = yB - a[:, None] * p
        JtJ = np.einsum("bka,bkc->bac", J, J)
        dg = np.einsum("baa->ba", JtJ)
        JtJ = JtJ + (lam[:, None] * np.maximum(dg, 1e-12))[:, :, None] * np.eye(4)[None]
        x = np.linalg.solve(JtJ, np.einsum("bka,bk->ba", J, r)[..., None])[..., 0]
        Rn = np.matmul(_expm(x[:, 1:]), Rr)
        fn, pn, an = _fval(RmB, yB, Rn)
        ok = fn >= f
        Rr = np.where(ok[:, None, None], Rn, Rr)
        f, p, a = np.where(ok, fn, f), np.where(ok[:, None], pn, p), np.where(ok, an, a)
        lam = np.where(ok, lam / 3.0, lam * 10.0)
    f = f.reshape(K, S)
    j = f.argmax(axis=1)
    out = f[np.arange(K), j] / 2.0
    if full:
        return out, Rr.reshape(K, S, 3, 3)[np.arange(K), j], f_closed / 2.0, fg.max(axis=1) / 2.0
    return out


def dc_brute(Rm, y, n=200000, seed=1):
    """reference maximum over a dense random set of rotations (validation only)."""
    Rr = Rotation.random(n, random_state=seed).as_matrix()
    P = vec6(dc_of_rot(Rr)) @ Rm.T
    return float(((P @ y) ** 2 / (P * P).sum(1)).max() / 2.0)


# ============================================================================================ detector loop
def loop(E, tab, MFn, kmax, rank="full6", capv=C.CAPV, min_sep=5.0, stat=None):
    """c10_lib.mt_loop (mode hard, W = SH, no refinement; identical atoms for rank = full6) with the node ranked by
    s6 | s5 | sdc.  Deflation: the full 6-D LS fit of the selected candidate in every variant (only the statistic
    changes).  Atoms carry s6, s5, sdc of the selected candidate."""
    ns = E.ns
    R = MFn.copy()
    excluded = np.zeros(E.N, bool)
    atoms = []
    for _ in range(kmax):
        if excluded.all():
            break
        tpk, _pen = C.choose(E, None, "hard", None, C.SH, R=R)
        vf = C._vectors(R, tpk, capv).reshape(E.N, E.T, 6 * ns)
        y = np.einsum("gtr,grk->gtk", vf, E.Q)
        s6 = (y * y).sum(-1) / 2.0
        s5 = s6 - np.einsum("gtk,gk->gt", y, tab.nvec) ** 2 / 2.0
        if rank == "full6":
            sc = s6.copy()
            sc[excluded] = -np.inf
            j = int(np.argmax(sc))
        elif rank == "dev5":
            sc = s5.copy()
            sc[excluded] = -np.inf
            j = int(np.argmax(sc))
        else:
            sc = s5.copy()
            sc[excluded] = -np.inf
            flat = sc.ravel()
            order = np.argsort(-flat, kind="stable")
            best, j, pos = -np.inf, -1, 0
            while pos < len(order) and np.isfinite(flat[order[pos]]) and flat[order[pos]] > best:
                idx = order[pos:pos + CHUNK]
                idx = idx[np.isfinite(flat[idx])]
                gi_, tj_ = np.divmod(idx, E.T)
                v = dc_max(tab.Rm[gi_], tab.Ri[gi_], y[gi_, tj_])
                k = int(np.argmax(v))
                if v[k] > best:
                    best, j = float(v[k]), int(idx[k])
                pos += CHUNK
                if stat is not None:
                    stat["ndc"] = stat.get("ndc", 0) + len(idx)
            if stat is not None:
                stat["nit"] = stat.get("nit", 0) + 1
        if j < 0:
            break
        gi, tj = divmod(j, E.T)
        if not np.isfinite(sc[gi, tj]):
            break
        sdc = float(dc_max(tab.Rm[[gi]], tab.Ri[[gi]], y[[gi], [tj]])[0])
        src, tau, tp, vsel, G = E.XN[gi], float(E.taus[tj]), tpk[gi, tj].copy(), vf[gi, tj], E.G[gi]
        atoms.append({"x": np.asarray(src, float), "tau": tau, "s6": float(s6[gi, tj]), "s5": float(s5[gi, tj]),
                      "sdc": sdc, "s": float(s6[gi, tj])})
        excluded = excluded | (np.linalg.norm(E.XN[:, :2] - src[None, :2], axis=1) <= min_sep)
        m, *_ = np.linalg.lstsq(G, vsel, rcond=None)
        C._deflate(R, (G @ m).reshape(2, ns, 3), tp)
    return atoms


def pack(atoms):
    return [[float(a["x"][0]), float(a["x"][1]), float(a["x"][2]), float(a["tau"]), a["s6"], a["s5"], a["sdc"]]
            for a in atoms]


# ============================================================================================ data
class Net:
    def __init__(self, cell, j, regime):
        self.nc = nc = xc.NetCfg(cell, j, regime)
        xops = xc.Ops(nc.STA, nc.corr, GRID)
        dops = C8.DOps(xops)
        self.E = C.Eng(GARR, nc.STA, nc.corr, TAUS, C.H_COARSE, G=dops.Gs, Q=dops.Qst)
        self.tab = Tab(self.E)


def draw_none(rng):
    return [], []


def stream(nc, tag, draw, nclut=None):
    """c9_run.stream (same generator call and rng streams) keeping the arrival amplitudes; nclut overrides the number
    of clutter arrivals (0 = noise only)."""
    b = [SEED, tag, nc.ci, nc.j, nc.ri]
    rng, rng_m, rng_h = (np.random.default_rng(b + [41]), np.random.default_rng(b + [233]),
                         np.random.default_rng(b + [61]))
    ncl = nc.nclut if nclut is None else nclut
    noise = GR.SNR * GR.make_cfg(None)["noise"]
    while True:
        srcs, taus = draw(rng)
        D, picks, site, amp = C8.generate_site(nc.STA, srcs, taus, ncl, rng, rng_m, rng_h, nc.clut_amp, cfg=None,
                                               statics=nc.statics)
        MF, sd = xc.matched(D)
        k = len(srcs)
        A = np.array([[amp[i][2 * e + b_] for i in range(len(nc.STA)) for b_ in range(2)] for e in range(k)]) / noise
        yield {"MFn": MF / sd, "truth": [np.asarray(s, float) for s in srcs], "taus": [float(t) for t in taus],
               "amp": np.sqrt((A * A).sum(axis=1)).tolist() if k else [],
               "amp_max": A.max(axis=1).tolist() if k else []}


# ============================================================================================ tube formula
def lk_curvatures(Rm, n=4000, seed=7):
    """L3 (volume) and L1 of U = {Rm m / |Rm m| : m unit double couple} in S^5; Monte Carlo over Haar rotations
    (the stabiliser of a double couple in SO(3) has order 4; vol SO(3) = 8 pi^2 in rotation-vector units)."""
    Rr = Rotation.random(n, random_state=seed).as_matrix()
    Mt = dc_of_rot(Rr)
    p = vec6(Mt) @ Rm.T
    rho = np.linalg.norm(p, axis=1)
    u = p / rho[:, None]
    c1 = [_comm(WGEN[k], Mt) for k in range(3)]
    dp = np.stack([vec6(c) @ Rm.T for c in c1], axis=1)                         # (n, 3, 6)
    du = (dp - np.einsum("nk,nc->nkc", np.einsum("nkc,nc->nk", dp, u), u)) / rho[:, None, None]
    g = np.einsum("nkc,nlc->nkl", du, du)
    sq = np.sqrt(np.maximum(np.linalg.det(g), 0.0))
    B = np.concatenate([u[:, None, :], du], axis=1).transpose(0, 2, 1)          # (n, 6, 4)
    Qb, _ = np.linalg.qr(B)
    II = np.zeros((n, 3, 3, 6))
    for j in range(3):
        for k in range(j, 3):
            h = 0.5 * (_comm(WGEN[j], c1[k]) + _comm(WGEN[k], c1[j]))
            d2 = vec6(h) @ Rm.T
            d2 = d2 - np.einsum("nca,na->nc", Qb, np.einsum("nca,nc->na", Qb, d2))
            II[:, j, k] = II[:, k, j] = d2 / rho[:, None]
    gi = np.linalg.inv(g)
    Hv = np.einsum("njk,njkc->nc", gi, II)
    H2 = (Hv * Hv).sum(-1)
    II2 = np.einsum("nja,nkb,njkc,nabc->n", gi, gi, II, II)
    S = 6.0 + H2 - II2
    L3 = 2.0 * np.pi ** 2 * sq.mean()
    L1 = 2.0 * np.pi ** 2 * (S * sq).mean() / (4.0 * np.pi)
    return float(L3), float(L1), float(2.0 * np.pi ** 2 * sq.std() / np.sqrt(n))


def tube_tail(x, L3, L1):
    """P(2 sdc > x), x = t^2 (two-term tube formula; capped at the chi-square_5 bound and at 1)."""
    x = np.asarray(x, float)
    e = np.exp(-x / 2.0)
    pr = L3 * (x - 1.0) * e / (4.0 * np.pi ** 2) + L1 * e / (2.0 * np.pi)
    return np.minimum(np.minimum(pr, stats.chi2.sf(x, 5)), 1.0)


def tube_quantile(alpha, L3, L1):
    """x with tube_tail(x) = alpha (upper tail)."""
    from scipy.optimize import brentq
    return brentq(lambda x: np.log(max(float(tube_tail(x, L3, L1)), 1e-300)) - np.log(alpha), 3.0, 400.0)


def eff_dof(x, pval):
    """k (real) with P(chi2_k > x) = pval."""
    from scipy.optimize import brentq
    return brentq(lambda k: stats.chi2.logsf(x, k) - np.log(pval), 0.2, 12.0)
