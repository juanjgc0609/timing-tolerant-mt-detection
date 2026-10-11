"""
xc.py - common harness for the exploratory contribution study (explore_contrib).  SEED 999 ONLY.

Networks: task2a_confirm.network(999, cell, j) (seed-999-derived; j may exceed the 5 calibration networks of 2A).
Catalogue and station corrections: task2a_confirm.catalogue_real / estimate_station_corrections, re-implemented here
only to pass a generator cfg (velocity error) and scaled statics; with cfg=None and static_sd=0.3 it reproduces
task2b_ablation.setup_net exactly.
Windows: gen_realistic.generate_real with new stream tags [999, TAG, ci, j, ri, 41|233|61] (TAG >= 3000, never used
before), optional custom source drawer (dense overlaps) and cfg.
Methods: BOMP-PK+SC with arbitrary peak half-width, TS-SFW+SC (f = 0.01, frozen in 2E), BOMP-PK-SLIDE+SC, MT-FW+SC, all
on an arbitrary lattice grid.  Oracle ceiling record.  Recall at caps and hierarchical bootstrap.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
if WORK not in sys.path:
    sys.path.insert(0, WORK)

import task2a_confirm as T  # noqa: E402
import task2c_methods as M  # noqa: E402
import task2d_methods as D2  # noqa: E402
import task2e_methods as E  # noqa: E402
import task1x_mtglrt as X  # noqa: E402
from certified_fw import certified_fw_record  # noqa: E402
from gen_realistic import generate_real  # noqa: E402
from station_corrections import estimate_station_corrections  # noqa: E402
from task0_common import DEFAULT_TAUS, DT, NT, VP, VS  # noqa: E402
from task1o_consistent import draw_sources, matched  # noqa: E402
from task1y_clutter import CLUT_AMP  # noqa: E402

SEED = 999
CAPS = (0.1, 0.25)
F_FROZEN = 0.01


# ============================================================================================ grids
def make_grid(step=3.0, depths=(5.0, 9.0, 13.0), half=21.0):
    ax = np.arange(-half, half + 1e-9, step) + 0.37
    return [np.array([x, y, z], float) for x in ax for y in ax for z in depths]


# ============================================================================================ networks
def catalogue_cfg(STA, K, statics, seeds, n_clutter, clut_amp, cfg=None, gate=25.0, min_arr=3, cat_err=T.CAT_ERR):
    """task2a_confirm.catalogue_real with a generator cfg (identical when cfg is None)."""
    from station_corrections import peak_delays
    rng = np.random.default_rng(list(seeds) + [0])
    rng_m = np.random.default_rng(list(seeds) + [1])
    rng_h = np.random.default_rng(list(seeds) + [2])
    rng_e = np.random.default_rng(list(seeds) + [3])
    out, truth, drawn = [], [], 0
    while len(out) < K:
        drawn += 1
        src = np.array([rng.uniform(-15, 15), rng.uniform(-15, 15), rng.uniform(5, 13)])
        tau = float(rng.uniform(5.5, 7.5))
        D, _p, _tr, _a = generate_real(STA, [src], [tau], n_clutter, rng, rng_m, rng_h, clut_amp, cfg=cfg,
                                       statics=statics)
        _dl, en = peak_delays(STA, D, src, tau)
        if (en >= gate).sum() < min_arr:
            continue
        src_c = src + rng_e.normal(0, cat_err[0], 3)
        src_c[2] = float(np.clip(src_c[2], 1.0, 20.0))
        out.append((D, src_c, tau + float(rng_e.normal(0, cat_err[1]))))
        truth.append((src, tau))
    return out, drawn, truth


class NetCfg:
    def __init__(self, cell, j, regime, static_sd=0.3, cfg=None, sc=True):
        self.cell, self.j, self.regime = cell, j, regime
        ci, _kind, nst, nclut = T.net_setup(SEED, cell)
        self.ci, self.nst, self.nclut = ci, nst, nclut
        self.ri = T.REGIMES.index(regime)
        STA, statics = T.network(SEED, cell, j)
        self.STA = STA
        self.statics = statics * (static_sd / T.STATIC_SD_2A)
        self.cfg = cfg
        self.clut_amp = CLUT_AMP[regime]
        self.cat, self.drawn, self.cat_truth = catalogue_cfg(STA, T.K_CAT, self.statics,
                                                             [SEED, T.TAG_CAT, ci, j, self.ri], nclut,
                                                             self.clut_amp, cfg=cfg)
        corr, _n = estimate_station_corrections(STA, self.cat)
        self.corr_cat = corr
        self.corr = corr if sc else np.zeros_like(corr)


def window_stream(nc, tag, draw=draw_sources):
    b = [SEED, tag, nc.ci, nc.j, nc.ri]
    rng, rng_m, rng_h = (np.random.default_rng(b + [41]), np.random.default_rng(b + [233]),
                         np.random.default_rng(b + [61]))
    while True:
        srcs, taus = draw(rng)
        D, picks, tr, _amp = generate_real(nc.STA, srcs, taus, nc.nclut, rng, rng_m, rng_h, nc.clut_amp, cfg=nc.cfg,
                                           statics=nc.statics)
        MF, sd = matched(D)
        yield {"D": D, "MFn": MF / sd, "sd": sd, "truth": [np.asarray(s, float) for s in srcs],
               "taus": [float(t) for t in taus], "picks": picks}


# ============================================================================================ operators
class Ops:
    """per network x grid x corrections: BlockDict, PKDict, TS-SFW Net, MT-FW operators."""

    def __init__(self, STA, corr, grid):
        t0 = time.time()
        self.STA, self.corr = STA, corr
        self.grid = grid
        garr = np.asarray(grid)
        self.Gs = [X.radiation_operator(g, STA) for g in garr]
        self.Qs = [X.basis(G) for G in self.Gs]
        self.BD = M.BlockDict(STA, grid, DEFAULT_TAUS, corr, Gs=self.Gs)
        self.PD = D2.PKDict(self.BD)
        self.net = E.Net(STA, corr, grid, DEFAULT_TAUS, self.Gs, self.Qs, self.BD)
        self.t_build = time.time() - t0


# ============================================================================================ peak search width
def peak_table_w(R, sh):
    offs = np.arange(-sh, sh + 1)
    ii = np.clip(np.arange(NT)[:, None] + offs[None, :], 0, NT - 1)
    Et = np.sum(np.ascontiguousarray(R.transpose(0, 2, 1)) ** 2, axis=-1)
    W = Et[:, ii]
    a = W.argmax(axis=-1)
    return ii[np.arange(NT)[None, :], a]


class PKDictW(D2.PKDict):
    def __init__(self, BD, sh):
        super().__init__(BD)
        self.sh = sh

    def read_idx(self, R, peak):
        if not peak:
            return self.BD.idx
        Kt = peak_table_w(R, self.sh)
        return Kt[self.st4, self.BD.idx]


def bomp_pk(PD, MFn, truth, kmax=5, diag=None):
    d = {} if diag is None else diag
    rec = D2.bomp_pk_record(PD, MFn, truth, kmax=kmax, diag=d)
    locs = [PD.BD.garr[a[0]].tolist() + [float(PD.BD.taus[a[1]])] for a in d["atoms"]]
    return [(float(s), int(h)) for s, h in rec], locs


def tssfw(ops, w, truth, f=None, lam=None, kmax=5):
    win = E.Win(w["D"], w["MFn"], w["sd"])
    dg = {}
    rec = E.tssfw_record(ops.net, win, truth, f=f, lam=lam, kmax=kmax, diag=dg)
    return [(float(s), int(h)) for s, h in rec], dg["P"][:, :4].tolist(), dg


def mtfw(ops, MFn, truth):
    rec = certified_fw_record(ops.STA, MFn, ops.grid, DEFAULT_TAUS, ops.Gs, ops.Qs, truth, tcorr=ops.corr)
    return [(float(s), int(h)) for s, h in rec]


# ============================================================================================ oracle ceiling
def true_arrivals(w, ns):
    k = len(w["truth"])
    A = np.zeros((k, 2, ns))
    for i in range(ns):
        for s in range(k):
            for b in range(2):
                A[s, b, i] = w["picks"][i][2 * s + b][0]
    return A


def oracle_scores(ops, w):
    """score of each true event with its TRUE position, origin and arrival times, after a joint least-squares fit of
    all true events (continuous atom model of 2E) and removal of the other events' fitted contributions:
    E_k / 2 with E_k the projection energy of the continuous residual at event k's true times (units of the record
    scores)."""
    net = ops.net
    k = len(w["truth"])
    if k == 0:
        return []
    win = E.Win(w["D"], w["MFn"], w["sd"])
    A = true_arrivals(w, net.ns)
    P = np.zeros((k, net.nv))
    for s in range(k):
        src = w["truth"][s]
        d = np.linalg.norm(net.STA - src, axis=1)
        T0 = d[None, :] * E.VINV[:, None] + w["taus"][s] + net.corrT
        P[s, :3] = src
        P[s, 3] = w["taus"][s]
        P[s, 4:] = (A[s] - T0).ravel()
    m, _ = E.convex_step(net, win, P, np.zeros((k, 6)), 0.0)
    out = []
    for s in range(k):
        oth = [q for q in range(k) if q != s]
        a_old, T_old = E.pred_old(net, win, P[oth], m[oth])
        v, _mu = E.refine_obj(net, win, P[s], a_old, T_old, grad=False)
        out.append(-float(v))            # = E/2
    return out


def ceiling_record(false_rec, osc):
    r = [(s, -1) for s, h in false_rec if h < 0] + [(s, k) for k, s in enumerate(osc)]
    return sorted(r, key=lambda x: -x[0])


# ============================================================================================ recall & bootstrap
def recall_caps(recs, nte, caps=CAPS):
    sc, hit, key = [], [], []
    for wi, r in enumerate(recs):
        for s_, h in r:
            sc.append(s_)
            hit.append(h)
            key.append(wi * 64 + h if h >= 0 else -1)
    if not sc or nte == 0:
        return [0.0 for _ in caps], [np.inf for _ in caps]
    sc, hit, key = np.array(sc), np.array(hit), np.array(key)
    o = np.argsort(-sc, kind="stable")
    sc, hit, key = sc[o], hit[o], key[o]
    _u, first = np.unique(key, return_index=True)
    newhit = np.zeros(len(sc), bool)
    newhit[first] = True
    newhit &= key >= 0
    H = np.cumsum(newhit)
    F = np.cumsum(hit < 0)
    last = np.r_[np.where(np.diff(sc) != 0)[0], len(sc) - 1]
    H, F, th = H[last], F[last], sc[last]
    out, thr = [], []
    for c in caps:
        ok = F / len(recs) <= c
        if ok.any():
            j = int(np.argmax(np.where(ok, H, -1)))
            out.append(float(H[j] / nte))
            thr.append(float(th[j]))
        else:
            out.append(0.0)
            thr.append(np.inf)
    return out, thr


def hier_boot(groups, methods, recs, ks, B=1000, seed=12345, pairs=(), caps=CAPS, strata=None):
    """groups: list of lists of window indices (one list per network).  recs[m]: list of per-window records.
    Returns point recalls, and for each pair (a, b) point diff and 95% percentile CI (networks then windows)."""
    rng = np.random.default_rng(seed)
    allw = [w for g in groups for w in g]
    nte = sum(ks[w] for w in allw)
    pt = {m: recall_caps([recs[m][w] for w in allw], nte, caps)[0] for m in methods}
    bs = {p: [] for p in pairs}
    bm = {m: [] for m in methods}
    strata = strata if strata is not None else [0] * len(groups)
    sidx = {}
    for g, s_ in enumerate(strata):
        sidx.setdefault(s_, []).append(g)
    for _ in range(B):
        gi = []
        for s_, gl in sidx.items():
            gi.extend(np.asarray(gl)[rng.integers(0, len(gl), len(gl))].tolist())
        ws = []
        for g in gi:
            gg = groups[g]
            ws.extend(np.asarray(gg)[rng.integers(0, len(gg), len(gg))].tolist())
        n_ = sum(ks[w] for w in ws)
        rc = {m: recall_caps([recs[m][w] for w in ws], n_, caps)[0] for m in methods}
        for m in methods:
            bm[m].append(rc[m])
        for a, b in pairs:
            bs[(a, b)].append(np.subtract(rc[a], rc[b]))
    ci = {m: np.percentile(np.array(bm[m]), [2.5, 97.5], axis=0).T.tolist() for m in methods}
    dif = {}
    for a, b in pairs:
        arr = np.array(bs[(a, b)])
        dif[(a, b)] = {"point": np.subtract(pt[a], pt[b]).tolist(),
                       "ci": np.percentile(arr, [2.5, 97.5], axis=0).T.tolist(),
                       "p_le0": (np.mean(arr <= 0, axis=0)).tolist()}
    return pt, ci, dif


def fmt_ci(p, ci):
    return f"{p:+.3f} [{ci[0]:+.3f}, {ci[1]:+.3f}]"


# ============================================================================================ wide-window TS-SFW
class _Wide:
    """temporarily widen the peak search of the TS-SFW oracle (task2d_methods.peak_table) and the shift bound."""

    def __init__(self, sh):
        self.sh = sh

    def __enter__(self):
        self._pt, self._sm = D2.peak_table, E.SHIFT_MAX
        sh = self.sh
        D2.peak_table = lambda R: peak_table_w(R, sh)
        E.SHIFT_MAX = sh * DT
        return self

    def __exit__(self, *a):
        D2.peak_table, E.SHIFT_MAX = self._pt, self._sm


def wide_net(ops, sh):
    import copy
    n2 = copy.copy(ops.net)
    ns = n2.ns
    n2.lo = n2.lo.copy()
    n2.hi = n2.hi.copy()
    n2.lo[4:] = -sh * DT
    n2.hi[4:] = sh * DT
    return n2


def tssfw_wide(ops, net_w, w, truth, sh, f=None, lam=None, kmax=5):
    win = E.Win(w["D"], w["MFn"], w["sd"])
    dg = {}
    with _Wide(sh):
        rec = E.tssfw_record(net_w, win, truth, f=f, lam=lam, kmax=kmax, diag=dg)
    return [(float(s), int(h)) for s, h in rec], dg["P"][:, :4].tolist(), dg


def bomp_hit_best(rec, k):
    best = [-np.inf] * k
    for s, h in rec:
        if h >= 0:
            best[h] = max(best[h], s)
    return best


def ceilings(bomp_rec, osc, truth, bomp_locs, excl_km=20.0):
    """CEIL-A: each true event scored max(oracle, best BOMP-PK hit score); false atoms = BOMP-PK's false atoms.
    CEIL-B: same, keeping only BOMP-PK false atoms farther than excl_km (epicentral) from every true event."""
    k = len(osc)
    bb = bomp_hit_best(bomp_rec, k)
    ts = [(max(osc[i], bb[i]), i) for i in range(k)]
    fa = [(s, -1) for (s, h) in bomp_rec if h < 0]
    fb = []
    for (s, h), loc in zip(bomp_rec, bomp_locs):
        if h >= 0:
            continue
        dmin = min([np.linalg.norm(np.asarray(loc[:2]) - t[:2]) for t in truth] or [np.inf])
        if dmin > excl_km:
            fb.append((s, -1))
    A = sorted(fa + ts, key=lambda x: -x[0])
    Bc = sorted(fb + ts, key=lambda x: -x[0])
    return A, Bc


# ============================================================================================ arrival masking (leakage)
def bomp_pk_mask(PD, MFn, truth, kmax=5, min_sep=5.0, mw=8, amp_gate=3.0, diag=None):
    """BOMP-PK+SC in which, after an atom is accepted, the samples within +-mw of each of its read arrival samples
    whose fitted arrival amplitude (norm of the 3-component prediction, MFn units) is >= amp_gate are removed from
    the search residual of later iterations (each arrival is explained by at most one event, as in pick association).
    The least-squares refit and the scores are otherwise those of bomp_pk_record."""
    garr = PD.BD.garr
    ns = PD.ns
    excluded = np.zeros(PD.N, bool)
    S, rec, run, locs = [], [], np.inf, []
    R = MFn.copy()
    mask = np.zeros((ns, MFn.shape[-1]), bool)
    for _ in range(kmax):
        Rs = R * (~mask)[:, None, :]
        e, I, _h = PD.energy(Rs, True)
        e[excluded] = -np.inf
        j = int(np.argmax(e))
        n, t = divmod(j, PD.T)
        s = float(e[n, t]) / 2.0
        if not np.isfinite(s):
            break
        run = min(run, s)
        src = garr[n]
        rec.append((run, M.hit_of(src, truth)))
        locs.append(src.tolist() + [float(PD.BD.taus[t])])
        a = (n, t, I[n, t].copy())
        S.append(a)
        excluded = M.excl_mask(excluded, garr, src, min_sep)
        K = PD.gram(S)
        b = PD.corr_atoms(MFn, S).ravel()
        U, *_ = np.linalg.lstsq(K, b, rcond=None)
        U = U.reshape(len(S), 6)
        R = PD.residual(MFn, S, U)
        pred = (PD.BD.Gf[n] @ U[-1]).reshape(2, ns, 3)
        amp = np.linalg.norm(pred, axis=-1)
        for bb in range(2):
            for i in range(ns):
                if amp[bb, i] >= amp_gate:
                    k0 = int(a[2][bb, i])
                    mask[i, max(0, k0 - mw):k0 + mw + 1] = True
    if diag is not None:
        diag["atoms"] = S
    return [(float(s), int(h)) for s, h in rec], locs


def ratio_filter(rec, locs, rho, radius=25.0):
    """leakage rule: drop atom k if its score < rho x the score of an earlier KEPT atom within radius km (epicentral).
    Running minimum re-applied to the kept atoms."""
    kept, out, run = [], [], np.inf
    for (s, h), l in zip(rec, locs):
        l = np.asarray(l[:2])
        if any(np.linalg.norm(l - np.asarray(lk[:2])) <= radius and s < rho * sk for sk, lk in kept):
            continue
        kept.append((s, l))
        run = min(run, s)
        out.append((run, h))
    return out


def mtfw_mask(ops, MFn, truth, kmax=5, min_sep=5.0, mw=6, amp_gate=3.0):
    """MT-FW+SC (certified_fw_record loop, enclosure oracle) with the arrival masking of bomp_pk_mask."""
    from certified_fw import Lattice, _exclusion_update, certified_oracle, stack_Q
    STA, ns = ops.STA, len(ops.STA)
    R = MFn.copy()
    lat = Lattice(ops.grid, DEFAULT_TAUS)
    garr = lat.garr
    excluded = np.zeros(len(garr), bool)
    Qst = stack_Q(ops.Qs)
    mask = np.zeros((ns, MFn.shape[-1]), bool)
    rec, run, locs = [], np.inf, []
    for _ in range(kmax):
        Rs = R * (~mask)[:, None, :]
        s, gi, vf, tpk = certified_oracle(STA, Rs, lat, DEFAULT_TAUS, ops.Qs, excluded, "enclosure", ops.corr, None, Qst)
        if gi is None:
            break
        src = garr[gi]
        run = min(run, s)
        rec.append((float(run), int(M.hit_of(src, truth))))
        locs.append(src.tolist())
        excluded = _exclusion_update(excluded, garr, src, min_sep)
        G = ops.Gs[gi]
        m, *_ = np.linalg.lstsq(G, vf, rcond=None)
        pred = (G @ m).reshape(2, ns, 3)
        for b in range(2):
            for i in range(ns):
                t0 = int(tpk[b, i])
                lo, hi = max(0, t0 - X.AC_HALF), min(NT, t0 + X.AC_HALF + 1)
                R[i, :, lo:hi] -= np.outer(pred[b, i], X.AC[X.AC_HALF - (t0 - lo):X.AC_HALF + (hi - t0)])
                if np.linalg.norm(pred[b, i]) >= amp_gate:
                    mask[i, max(0, t0 - mw):t0 + mw + 1] = True
    return rec, locs
