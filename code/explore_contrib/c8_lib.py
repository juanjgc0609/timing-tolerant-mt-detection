"""
c8_lib.py - C8: per-station/phase amplitude (site) terms inside the moment-tensor subspace certificate.  SEED 999 ONLY.

Generator worlds (gen_realistic.generate_real, copied as generate_site so that the site factors can be returned and,
optionally, overridden; the rng consumption is unchanged, so with mode "iid" the output is bit-identical to
generate_real; checked in c8_validate.py):
  mode "iid"   site factors exp(N(0, sigma)) per station and phase, REDRAWN IN EVERY WINDOW (this is what the default
               generator does: rng_h draws them at every call; the 20 catalogue events each get their own draw too).
  mode "pers"  site factors fixed per network (station, phase): exp(sigma z), z ~ N(0, 1) from rng
               [999, TAG_SITE, ci, j] (shared by the regimes, as the statics), used by the catalogue and by every
               window; the rng_h draw is still consumed.
  ps_same      the S factor equals the P factor of the same station (a per-station site term).
Catalogue (+SC timing corrections, the registered rule) is generated in the same world (same streams as
task2a_confirm.catalogue_real; bit-identical to xc.NetCfg for the default world).

Catalogue amplitude terms (SITE-CAT, SITE-CAT-ST):
  per catalogue event: peak times from c6_lib.pick_event (joint P/S energy peaks within +-1 s of the homogeneous
  prediction + the +SC correction), 3-component matched-filter vector at the peak sample, gate energy >= 25 (the +SC
  gate); radiation-averaged log-amplitude residuals with event terms removed, median per station/phase (see
  site_terms; estimator chosen on catalogue data of networks 60-63 before any detection run, c8_sitedev.py).
  -ST: one term per station (median over both phases).

Dense detector engine (mtfw_general): the c6_lib.mtfw_dense loop (energy peak search +-4 samples on the residual,
statistic over the whole (node, origin) lattice, lowest-index tie rule, running minimum, 5 km exclusion, sequential
deflation of the fitted prediction at the peak samples) with a pluggable statistic:
  "ls"     ||Q^T v||^2 / 2, Q = orthonormal basis of the (optionally row-scaled) operator; v optionally capped at capv
           (capv = 6, wrow = None reproduces c6_lib.mtfw_dense exactly)
  "huber"  robust likelihood ratio sum_b rho_c(|v_b|) - min_a sum_b rho_c(|v_b - Q_b a|), Huber rho on the norm of
           each station/phase 3-vector, IRLS (4 passes) from the LS fit; c = C_HUB * sigma, C_HUB = sqrt(chi2_3(0.95))
           = 2.7955 (the Gaussian block norm is left untouched with probability 0.95); sigma = 1 (noise units,
           HUBER-N) or sigma = max(1, median_b |r_b| / 1.5382) from the LS residual of the candidate (HUBER-S;
           1.5382 = sqrt(median chi2_3)).  Deflation with the IRLS fit.
  "gain"   profile likelihood with a free positive gain per station/phase block and a log-normal prior:
           T = |v|^2/2 - min_{m, u} [ sum_b |v_b - e^{u_b} G_b m|^2 / 2 + sum_b u_b^2 / (2 sp^2) ]
           (sp = 0.5 = the generator's nominal site sd: GAIN; sp = inf, u free, e^u >= 0: GAIN-FREE, i.e. only the
           polarisation directions and the sign of each block are modelled).  Alternating: u by a 13-point grid in
           [-3 sp, 3 sp] and 5 safeguarded Newton steps (GAIN) / closed form max(0, <v_b, g_b>) / |g_b|^2 (GAIN-FREE);
           m by weighted LS; 3 rounds from the LS fit.  Deflation with e^{u_b} G_b m.
Row scaling wrow (6 ns,) multiplies the operator rows (b, i, c) by the site factor of (b, i): the exact whitened
GLRT with known gains (ORACLE: true factors; SITE-CAT: catalogue estimates).
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
import c6_lib as L  # noqa: E402
import xc  # noqa: E402
import task2a_confirm as T  # noqa: E402
from station_corrections import estimate_station_corrections, peak_delays  # noqa: E402
from task1y_clutter import CLUT_AMP  # noqa: E402

X, M, D2 = L.X, L.M, L.D2
NT, DT = L.NT, L.DT
SEED = 999
TAG_SITE = 3809
C_HUB = float(np.sqrt(7.814727903251178))      # chi2_3 0.95 quantile
MED_CHI3 = float(np.sqrt(2.365973884375338))   # sqrt of the chi2_3 median


# ============================================================================================ generator
def generate_site(STA, sources, taus, n_clutter, rng, rng_m, rng_h, clut_amp, cfg=None, meta=None, statics=None,
                  site_override=None, ps_same=False):
    """gen_realistic.generate_real (verbatim logic) returning also the site factors used (ns, 2)."""
    c = GR.make_cfg(cfg)
    ns = len(STA)
    D = np.zeros((ns, 3, NT))
    picks = [[] for _ in range(ns)]
    amp = [[] for _ in range(ns)]
    st = rng_h.normal(0.0, c["static_sd"], size=(ns, 2))
    if statics is not None:
        st = np.asarray(statics, float).reshape(ns, 2)
    dn = rng_h.normal(0.0, float(np.tan(np.radians(c["dir_deg"]))), size=(ns, 2, 3))
    site = np.exp(rng_h.normal(0.0, c["site_sd"], size=(ns, 2)))
    if site_override is not None:
        site = np.asarray(site_override, float).reshape(ns, 2).copy()
    if ps_same:
        site[:, 1] = site[:, 0]
    rx = GR._child(rng_h)
    sos = GR._coda_filter(c["coda_band"])
    TT_ = GR.T
    for src, tau in zip(sources, taus):
        src = np.asarray(src, float)
        d = np.linalg.norm(STA - src, axis=1)
        Mt = GR.rand_dc(rng_m)
        fc = float(rx.uniform(*c["fc_range"]))
        phi = float(np.radians(rx.uniform(-c["phase_max_deg"], c["phase_max_deg"])))
        ew = GR.wavelet_energy(fc)
        if meta is not None:
            meta.append({"fc": fc, "phase_deg": float(np.degrees(phi))})
        tt = {b: GR.travel_time_gradient(src, STA, b, c) for b in ("P", "S")}
        for i in range(ns):
            for bi, (b, A) in enumerate((("P", GR.AP), ("S", GR.AS))):
                ta = float(tt[b][i]) + tau + float(st[i, bi])
                picks[i].append((ta, b))
                vec = GR.radiation_gradient(Mt, src, STA[i], b, c)
                a0 = np.linalg.norm(vec)
                if a0 > 0:
                    w = vec / a0 + dn[i, bi]
                    vec = a0 * site[i, bi] * w / np.linalg.norm(w)
                av = GR.SNR * A * float(np.linalg.norm(vec)) / d[i]
                amp[i].append(av)
                D[i] += GR.SNR * A * np.outer(vec, GR.source_wavelet(TT_ - ta, fc, phi)) / d[i]
                GR._add_coda(D[i], ta, av, ew, c, sos, rx)
    ref = np.asarray(sources[0], float) if len(sources) else np.zeros(3)
    scale = float(np.median(np.linalg.norm(STA - ref, axis=1)))
    for _ in range(n_clutter):
        i = int(rng.integers(ns))
        tc = float(rng.uniform(2, TT_[-1] - 4))
        q = rng.normal(size=3)
        q /= np.linalg.norm(q)
        D[i] += GR.SNR * clut_amp * np.outer(q, GR.ricker(TT_ - tc)) / scale
        picks[i].append((tc, "P" if rng.random() < 0.5 else "S"))
        amp[i].append(GR.SNR * clut_amp / scale)
    D += rng.normal(scale=GR.SNR * c["noise"], size=D.shape)
    for _ in range(int(c["spikes"])):
        i, k, j = int(rx.integers(ns)), int(rx.integers(3)), int(rx.integers(NT))
        D[i, k, j] += (1.0 if rx.random() < 0.5 else -1.0) * c["spike_amp"] * GR.SNR * c["noise"]
    return D, picks, site, amp


class World:
    def __init__(self, mode="iid", sigma=0.5, ps_same=False):
        assert mode in ("iid", "pers")
        self.mode, self.sigma, self.ps_same = mode, float(sigma), bool(ps_same)
        self.cfg = None if (sigma == 0.5) else {"site_sd": self.sigma}

    def name(self):
        return f"{self.mode}_s{self.sigma:g}" + ("_pssame" if self.ps_same else "")


class WNet:
    """network x regime x world: geometry, statics, (persistent) site, catalogue, +SC corrections, catalogue
    amplitude terms."""

    def __init__(self, cell, j, regime, world):
        self.cell, self.j, self.regime, self.world = cell, j, regime, world
        ci, _kind, nst, nclut = T.net_setup(SEED, cell)
        self.ci, self.nst, self.nclut = ci, nst, nclut
        self.ri = T.REGIMES.index(regime)
        STA, statics = T.network(SEED, cell, j)
        self.STA = STA
        self.statics = statics                                    # static_sd 0.3 (xc.NetCfg default)
        self.clut_amp = CLUT_AMP[regime]
        self.site_net = None
        if world.mode == "pers":
            z = np.random.default_rng([SEED, TAG_SITE, ci, j]).normal(size=(nst, 2))
            if world.ps_same:
                z[:, 1] = z[:, 0]
            self.site_net = np.exp(world.sigma * z)
        self.cat, self.cat_truth, self.cat_site = self._catalogue()
        corr, _n = estimate_station_corrections(STA, self.cat)
        self.corr = corr
        self.u_cat, self.u_info = site_terms(STA, self.cat, corr)
        self.u_cat_st = self.u_info["u_st"]

    def gen(self, srcs, taus, n_clutter, rng, rng_m, rng_h, meta=None):
        return generate_site(self.STA, srcs, taus, n_clutter, rng, rng_m, rng_h, self.clut_amp, cfg=self.world.cfg,
                             meta=meta, statics=self.statics, site_override=self.site_net,
                             ps_same=self.world.ps_same)

    def _catalogue(self, gate=25.0, min_arr=3):
        seeds = [SEED, T.TAG_CAT, self.ci, self.j, self.ri]
        rng = np.random.default_rng(list(seeds) + [0])
        rng_m = np.random.default_rng(list(seeds) + [1])
        rng_h = np.random.default_rng(list(seeds) + [2])
        rng_e = np.random.default_rng(list(seeds) + [3])
        out, truth, sites = [], [], []
        while len(out) < T.K_CAT:
            src = np.array([rng.uniform(-15, 15), rng.uniform(-15, 15), rng.uniform(5, 13)])
            tau = float(rng.uniform(5.5, 7.5))
            D, _p, site, _a = self.gen([src], [tau], self.nclut, rng, rng_m, rng_h)
            _dl, en = peak_delays(self.STA, D, src, tau)
            if (en >= gate).sum() < min_arr:
                continue
            src_c = src + rng_e.normal(0, T.CAT_ERR[0], 3)
            src_c[2] = float(np.clip(src_c[2], 1.0, 20.0))
            out.append((D, src_c, tau + float(rng_e.normal(0, T.CAT_ERR[1]))))
            truth.append((src, tau))
            sites.append(site)
        return out, truth, sites

    def stream(self, tag, draw):
        b = [SEED, tag, self.ci, self.j, self.ri]
        rng, rng_m, rng_h = (np.random.default_rng(b + [41]), np.random.default_rng(b + [233]),
                             np.random.default_rng(b + [61]))
        while True:
            srcs, taus = draw(rng)
            meta = []
            D, picks, site, _amp = self.gen(srcs, taus, self.nclut, rng, rng_m, rng_h, meta=meta)
            MF, sd = xc.matched(D)
            yield {"D": D, "MFn": MF / sd, "sd": sd, "truth": [np.asarray(s, float) for s in srcs],
                   "taus": [float(t) for t in taus], "picks": picks, "meta": meta, "site": site}


# ============================================================================================ catalogue amplitude terms
def _dc6(M_):
    return [M_[0, 0], M_[1, 1], M_[2, 2], M_[0, 1], M_[0, 2], M_[1, 2]]


def site_terms(STA, cat, corr, gate=25.0, rounds=3, min_obs=3):
    """radiation-averaged log-amplitude terms (estimator frozen after c8_sitedev.py, which compared it with the
    per-event moment-tensor ratio estimator on the catalogues of networks 60-63 only):
    l_kbi = log |v_kbi| - log A_kbi, A_kbi = RMS over 300 random double couples of |G_bi(src_k) m| (the expected
    radiation amplitude with 1/d spreading and the P/S constants of the model), minus the event term (median over the
    event's gated blocks of l - u), u_bi = median_k (>= min_obs gated picks, else 0), re-centred; 3 rounds.
    Returns u (2, ns) and info with the per-station variant u_st (median over both phases)."""
    from task1w_mech import rand_dc
    ns = len(STA)
    rng = np.random.default_rng(5)
    Ms = np.array([_dc6(rand_dc(rng)) for _ in range(300)])
    obs = []
    for D, src, tau in cat:
        MF, sd = xc.matched(D)
        MFn = MF / sd
        E = np.sum(MFn ** 2, axis=1)
        tobs, _tp, en = L.pick_event(E, STA, np.asarray(src, float), float(tau), L.Homog(), center_corr=corr)
        k = np.clip(np.rint(tobs / DT).astype(int), 0, NT - 1)
        v = np.array([[MFn[i, :, k[i, b]] for i in range(ns)] for b in range(2)])
        G = X.radiation_operator(np.asarray(src, float), STA).reshape(2, ns, 3, 6)
        rms = np.sqrt(np.mean(np.einsum("bick,mk->mbic", G, Ms) ** 2, axis=(0, 3)) * 3)
        obs.append((np.log(np.linalg.norm(v, axis=-1) / rms), (en >= gate).T))

    def collect(u):
        lr = [[[] for _ in range(ns)] for _ in range(2)]
        for l, ok in obs:
            if ok.sum() < 2:
                continue
            ev = np.median((l - u)[ok])
            for b in range(2):
                for i in range(ns):
                    if ok[b, i]:
                        lr[b][i].append(l[b, i] - ev)
        return lr
    u = np.zeros((2, ns))
    for _ in range(rounds):
        lr = collect(u)
        un = np.zeros((2, ns))
        est = np.zeros((2, ns), bool)
        for b in range(2):
            for i in range(ns):
                if len(lr[b][i]) >= min_obs:
                    un[b, i] = float(np.median(lr[b][i]))
                    est[b, i] = True
        if est.any():
            un[est] -= un[est].mean()
        u = un
    ust = np.zeros((2, ns))
    for _ in range(rounds):
        lr = collect(ust)
        un = np.zeros((2, ns))
        est_s = np.zeros(ns, bool)
        for i in range(ns):
            vals = lr[0][i] + lr[1][i]
            if len(vals) >= min_obs:
                un[:, i] = float(np.median(vals))
                est_s[i] = True
        if est_s.any():
            un[:, est_s] -= un[0, est_s].mean()
        ust = un
    return u, {"n_obs": [[len(lr[b][i]) for i in range(ns)] for b in range(2)], "u_st": ust}


# ============================================================================================ operators
class DOps:
    """dense-lattice operators: node positions, sample index table, radiation operators, bases (from xc.Ops)."""

    def __init__(self, xops):
        self.STA = np.asarray(xops.STA, float)
        self.ns = len(self.STA)
        self.grid = xops.grid
        self.garr = np.asarray(xops.grid, float)
        self.taus = np.asarray(xc.DEFAULT_TAUS, float)
        self.Gs = np.stack(xops.Gs)                                                  # (N, 6ns, 6)
        from certified_fw import stack_Q
        self.Qst = stack_Q(xops.Qs)
        self.IDX = L.idx_table(self.STA, xops.grid, L.Homog(), np.asarray(xops.corr, float), self.taus)
        self.Gp = np.linalg.pinv(self.Gs)                                            # (N, 6, 6ns)

    def scaled(self, wsite):
        """operators with rows (b, i, c) multiplied by wsite[b, i] (2, ns)."""
        w = np.repeat(np.asarray(wsite, float).reshape(-1), 3)                       # (6ns,)
        Gw = self.Gs * w[None, :, None]
        U, s, _ = np.linalg.svd(Gw, full_matrices=False)
        keep = s > 1e-8 * s.max(axis=1, keepdims=True)
        Q = U * keep[:, None, :]
        return Gw, Q


def _leaf(R, capv):
    Et = np.sum(np.ascontiguousarray(R.transpose(0, 2, 1)) ** 2, axis=-1)
    W = Et[:, D2.II_TAB]
    a = W.argmax(axis=-1)
    K = D2.II_TAB[np.arange(NT)[None, :], a]
    V = np.ascontiguousarray(R.transpose(0, 2, 1))
    if capv is not None:
        nv = np.linalg.norm(V, axis=-1, keepdims=True)
        V = V * np.minimum(1.0, capv / np.maximum(nv, 1e-12))
    return K, V


def _rho(r, c):
    return np.where(r <= c, 0.5 * r * r, c * r - 0.5 * c * c)


def _blocks(vf, Q):
    N, T_, R6 = vf.shape
    nb = R6 // 3
    vb = vf.reshape(N, T_, nb, 3)
    Qb = Q.reshape(N, nb, 3, 6)
    QQ = np.matmul(Qb.transpose(0, 1, 3, 2), Qb).reshape(N, nb, 36)                # (N, nb, 36)
    PvT = np.matmul(Qb.transpose(0, 1, 3, 2), vb.transpose(0, 2, 3, 1))              # (N, nb, 6, T)
    return N, T_, nb, vb, QQ, PvT


def _wsolve(w, QQ, PvT, ridge=1e-9):
    N, T_, _nb = w.shape
    A = np.matmul(w, QQ).reshape(N, T_, 6, 6) + ridge * np.eye(6)
    rhs = np.einsum("ntb,nbkt->ntk", w, PvT, optimize=True)
    return np.linalg.solve(A, rhs[..., None])[..., 0]


def stat_huber(vf, Q, scale_adapt, iters=4, chub=None):
    """vf (N, T, 6ns); Q (N, 6ns, 6).  Returns scores (N, T) and coefficients a (N, T, 6) in the Q basis."""
    N, T_, nb, vb, QQ, PvT = _blocks(vf, Q)
    QT = Q.transpose(0, 2, 1)
    a = PvT.sum(1).transpose(0, 2, 1)                                                # LS: Q^T v  (N, T, 6)
    nv = np.sqrt((vb * vb).sum(-1))

    def rnorm(a_):
        r_ = vf - np.matmul(a_, QT)
        return np.sqrt((r_.reshape(N, T_, nb, 3) ** 2).sum(-1))
    r = rnorm(a)
    if scale_adapt:
        sig = np.maximum(1.0, np.median(r, axis=-1) / MED_CHI3)
    else:
        sig = np.ones((N, T_))
    c = ((C_HUB if chub is None else chub) * sig)[..., None]
    for _ in range(iters):
        w = np.minimum(1.0, c / np.maximum(r, 1e-12))
        a = _wsolve(w, QQ, PvT)
        r = rnorm(a)
    sc = _rho(nv, c).sum(-1) - _rho(r, c).sum(-1)
    return sc, a


GRID_U = np.linspace(-3.0, 3.0, 13)


def _u_update(a_, c_, sp, newton=5):
    """argmin_u 0.5 a e^{2u} - c e^u + u^2 / (2 sp^2) on [-3 sp, 3 sp]: 13-point grid, then safeguarded Newton."""
    best = np.full(a_.shape, np.inf)
    u = np.zeros_like(a_)
    for uu in GRID_U * sp:
        e = np.exp(uu)
        f = 0.5 * a_ * e * e - c_ * e + uu * uu / (2 * sp * sp)
        upd = f < best
        best = np.where(upd, f, best)
        u = np.where(upd, uu, u)
    lo, hi = -3 * sp, 3 * sp
    for _ in range(newton):
        e = np.exp(u)
        g1 = a_ * e * e - c_ * e + u / (sp * sp)
        g2 = 2 * a_ * e * e - c_ * e + 1.0 / (sp * sp)
        step = np.where(g2 > 0, g1 / np.maximum(g2, 1e-12), 0.0)
        step = np.clip(step, -0.25 * sp, 0.25 * sp)
        un = np.clip(u - step, lo, hi)
        en = np.exp(un)
        fn = 0.5 * a_ * en * en - c_ * en + un * un / (2 * sp * sp)
        fo = 0.5 * a_ * e * e - c_ * e + u * u / (2 * sp * sp)
        u = np.where(fn <= fo, un, u)
    return u


def stat_gain(vf, Gs, Gp, sp, iters=3):
    """profile likelihood with per-block gains (log-normal prior sd sp; sp = inf: free non-negative gains).
    Returns scores (N, T), m (N, T, 6), gains (N, T, nb)."""
    N, T_, nb, vb, GG, PgT = _blocks(vf, Gs)
    GT = Gs.transpose(0, 2, 1)
    m = np.matmul(vf, Gp.transpose(0, 2, 1))                                         # LS start (N, T, 6)
    v2 = (vb * vb).sum(-1)

    def gains(m_):
        g = np.matmul(m_, GT).reshape(N, T_, nb, 3)
        a_ = (g * g).sum(-1)
        c_ = (vb * g).sum(-1)
        if np.isfinite(sp):
            u = _u_update(a_, c_, sp)
            return g, np.exp(u), (u * u).sum(-1) / (2 * sp * sp)
        return g, np.maximum(c_, 0.0) / np.maximum(a_, 1e-12), 0.0
    for _ in range(iters):
        _g, s, _pen = gains(m)
        m = _wsolve_gain(s, GG, PgT)
    g, s, pen = gains(m)
    res = ((vb - s[..., None] * g) ** 2).sum(-1).sum(-1)
    sc = 0.5 * v2.sum(-1) - (0.5 * res + pen)
    return sc, m, s


def _wsolve_gain(s, GG, PgT, ridge=1e-9):
    N, T_, _nb = s.shape
    A = np.matmul(s * s, GG).reshape(N, T_, 6, 6) + ridge * np.eye(6)
    rhs = np.einsum("ntb,nbkt->ntk", s, PgT, optimize=True)
    return np.linalg.solve(A, rhs[..., None])[..., 0]


def robust_argmax(vf, excluded, stat, Q, Gs, Gp, sp, chunk=1536):
    """exact argmax of a robust statistic over the non-excluded lattice by bound-ordered evaluation: the computed
    statistic never exceeds its bound (HUBER: sum_b rho_c(|v_b|) >= its LR because rho >= 0; GAIN: |v|^2/2 >=
    |v|^2/2 - objective because the objective is >= 0), so once the best value found exceeds the largest remaining
    bound the search stops.  Ties: lowest (node, origin) flat index, as the LS path."""
    N, T_, R6 = vf.shape
    nb = R6 // 3
    nv = np.sqrt((vf.reshape(N, T_, nb, 3) ** 2).sum(-1))
    if stat.startswith("huber"):
        if stat == "huber_s":
            bound = None                                                            # c depends on the candidate
        else:
            bound = _rho(nv, C_HUB).sum(-1)
    else:
        bound = 0.5 * (nv * nv).sum(-1)
    if bound is None:
        # scale-adaptive c >= C_HUB and rho_c(r) <= r^2/2, so the LS energy bound |v|^2/2 is valid
        bound = 0.5 * (nv * nv).sum(-1)
    bound = np.where(excluded[:, None], -np.inf, bound).ravel()
    order = np.argsort(-bound, kind="stable")
    best, bidx, bfit = -np.inf, -1, None
    for c0 in range(0, len(order), chunk):
        sel = order[c0:c0 + chunk]
        if not np.isfinite(bound[sel[0]]) or bound[sel[0]] < best:
            break
        n_, t_ = np.divmod(sel, T_)
        v_ = vf[n_, t_][:, None, :]
        if stat.startswith("huber"):
            sc, a_ = stat_huber(v_, Q[n_], stat == "huber_s")
            fits = (a_[:, 0],)
        else:
            sc, m_, g_ = stat_gain(v_, Gs[n_], Gp[n_], sp)
            fits = (m_[:, 0], g_[:, 0])
        sc = sc[:, 0]
        mx = float(sc.max())
        kr = np.where(sc == mx)[0]
        kk = int(kr[np.argmin(sel[kr])])
        idx = int(sel[kk])
        if mx > best or (mx == best and idx < bidx):
            best, bidx = mx, idx
            bfit = tuple(f[kk] for f in fits)
    if bidx < 0:
        return 0, 0, -np.inf, None
    gi, tj = divmod(bidx, T_)
    return int(gi), int(tj), best, bfit


def mtfw_general(dops, MFn, truth, kmax=5, stat="ls", capv=6.0, wsite=None, sp=0.5, min_sep=5.0, chub=None):
    ns = dops.ns
    R = MFn.copy()
    garr, IDX = dops.garr, dops.IDX
    N, T_ = IDX.shape[:2]
    st = np.arange(ns)[None, None, None, :]
    if wsite is not None:
        Gs, Q = dops.scaled(wsite)
        Gp = np.linalg.pinv(Gs) if stat == "gain" else None
    else:
        Gs, Q, Gp = dops.Gs, dops.Qst, dops.Gp
    excluded = np.zeros(N, bool)
    rec, run, locs = [], np.inf, []
    for _ in range(kmax):
        if excluded.all():
            break
        K, V = _leaf(R, capv if stat == "ls" else None)
        tpk = K[st, IDX]
        vf = V[st, tpk].reshape(N, T_, 6 * ns)
        if stat == "ls":
            y = np.einsum("gtr,grk->gtk", vf, Q)
            sc = (y * y).sum(-1) / 2.0
            sc[excluded] = -np.inf
            best = sc.max(axis=1)
            gi = int(np.argmax(best))                     # lowest node index among ties
            tj = int(np.argmax(sc[gi]))                   # first origin among ties
            s = float(sc[gi, tj])
        elif stat.startswith("huber"):                    # full lattice (the bound prunes little for Huber)
            sc, acoef = stat_huber(vf, Q, stat == "huber_s", chub=chub)
            sc[excluded] = -np.inf
            gi = int(np.argmax(sc.max(axis=1)))
            tj = int(np.argmax(sc[gi]))
            s = float(sc[gi, tj])
            fit = (acoef[gi, tj],)
        else:
            gi, tj, s, fit = robust_argmax(vf, excluded, stat, Q, Gs, Gp, sp)
        if not np.isfinite(s):
            break
        src = garr[gi]
        run = min(run, s)
        rec.append((float(run), int(M.hit_of(src, truth))))
        locs.append(src.tolist() + [float(dops.taus[tj])])
        excluded = excluded | (np.linalg.norm(garr[:, :2] - src[:2], axis=1) <= min_sep)
        if stat == "ls":
            m, *_ = np.linalg.lstsq(Gs[gi], vf[gi, tj], rcond=None)
            pred = (Gs[gi] @ m).reshape(2, ns, 3)
        elif stat.startswith("huber"):
            pred = (Q[gi] @ fit[0]).reshape(2, ns, 3)
        else:
            g = (Gs[gi] @ fit[0]).reshape(2 * ns, 3)
            pred = (g * fit[1][:, None]).reshape(2, ns, 3)
        tp = tpk[gi, tj]
        for b in range(2):
            for i in range(ns):
                t0 = int(tp[b, i])
                lo, hi = max(0, t0 - X.AC_HALF), min(NT, t0 + X.AC_HALF + 1)
                R[i, :, lo:hi] -= np.outer(pred[b, i], X.AC[X.AC_HALF - (t0 - lo):X.AC_HALF + (hi - t0)])
    return rec, locs


def norm_traces(MFn):
    """NORM-TR: each station trace divided by its peak 3-component amplitude over the window and multiplied by the
    median of those peaks (standard back-projection trace normalisation; noise no longer unit)."""
    A = np.linalg.norm(MFn, axis=1).max(axis=-1)                                     # (ns,)
    return MFn * (np.median(A) / np.maximum(A, 1e-12))[:, None, None]
