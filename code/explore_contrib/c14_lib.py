"""
c14_lib.py - C14: self-interference of the pick-free MT-FW detector (cap 6, registered +SC terms, hard +-0.25 s rule).
SEED 999 ONLY.

C13 showed that the false atoms of windows with events are created by the events themselves.  Here (1) every atom of the
unchanged detector is decomposed, with the truth, into the sources of its statistic, and (2) four truth-free remedies
are tested INSIDE the unchanged loop (c10_lib.mt_loop; mt_loop14 with an empty configuration reproduces it exactly).

Remedies (cfg keys of mt_loop14; all use only the residual and the atoms already accepted):
  r1 = (kind, Wc[, p])  arrival exclusivity.  After an atom is accepted, the samples within +-Wc of its read sample at
       every station (both phases) are CLAIMED.  kind 'hard': the residual read by later candidates is blanked there
       (weight 0): the statistic of a later candidate is the statistic of its unclaimed evidence alone.  kind 'soft':
       weight rho^p, rho = min(1, residual energy / original energy) of the claimed window (3 components) measured
       right after the deflation.  Weights multiply the residual used for peak reading and for the statistic; the
       deflation itself always works on the unweighted residual.
  r2 = (L, from0)  time-varying noise scale.  Per station and component, running robust scale of the CURRENT residual
       (sliding median of |R| over L samples / 0.6745) divided by the scale of the whole original trace of that station
       and component, floored at 1; the residual read by the candidates is divided by it.  from0 False: only after the
       first deflation (the first atom is the baseline atom).
  r3 = (kind, L)  adaptive deflation.  After the registered 6-parameter MT fit, every station-phase is subtracted with
       a small family around the template instead of the MT prediction:
         'amp'  predicted polarisation, free scalar amplitude, best shift within +-L samples
         'vec'  free 3-vector amplitude, best shift within +-L samples
         'fir'  free 3-vector on each of the 2L+1 shifted templates (short FIR shaping filter per component)
       The fit is least squares of the residual on the autocorrelation wavelet over the deflation window.
  r4 = (X, Y)  duplicate merge in origin time: a candidate within X km (epicentre) and within Y s (origin) of an
       accepted atom is the same event (excluded), in addition to the registered 5 km exclusion at all origins.

Diagnosis (truth used only here): the window is regenerated in parts (per event: direct waves and coda; clutter and
noise = the twin of C13), and the statistic of an atom, s0 = v . P v / 2 (v the capped residual vectors at the read
samples, P the projector of the node), is split LINEARLY over the parts of v:  direct waves of every event, coda of
every event, clutter + noise, and minus the deflations already applied.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import copy
import sys

import numpy as np
from scipy.ndimage import median_filter

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
for p in (HERE, WORK):
    if p not in sys.path:
        sys.path.insert(0, p)

import c8_lib as C8  # noqa: E402
import c10_lib as C  # noqa: E402
import gen_realistic as GR  # noqa: E402
import xc  # noqa: E402
from task0_common import DT, NT  # noqa: E402

X = C.X
SEED = 999
DEV_NETS = (1600, 1601, 1602)
EVAL_NETS = tuple(range(1610, 1622))
TAG_STD, TAG_DENSE = 1651, 1652
NEAR = 0.4                      # s, association of a read sample with a true arrival (diagnosis only)
AC, AH = X.AC, X.AC_HALF


# ============================================================================================ windows in parts
def clone(g, entropy):
    """copy of a generator including its seed sequence and spawn counter (as c13_twin.clone)."""
    ss = np.random.SeedSequence(entropy, n_children_spawned=g.bit_generator.seed_seq.n_children_spawned)
    bg = np.random.PCG64(ss)
    bg.state = copy.deepcopy(g.bit_generator.state)
    return np.random.Generator(bg)


def event_parts(STA, sources, taus, rng_m, rng_h, statics):
    """the event loop of c8_lib.generate_site (verbatim logic) keeping every event apart: direct waves Dd (K, ns, 3,
    NT), coda Dc (K, ns, 3, NT) and true arrival times arr (K, 2, ns)."""
    c = GR.make_cfg(None)
    ns = len(STA)
    K = len(sources)
    Dd, Dc, arr = np.zeros((K, ns, 3, NT)), np.zeros((K, ns, 3, NT)), np.zeros((K, 2, ns))
    rng_h.normal(0.0, c["static_sd"], size=(ns, 2))
    st = np.asarray(statics, float).reshape(ns, 2)
    dn = rng_h.normal(0.0, float(np.tan(np.radians(c["dir_deg"]))), size=(ns, 2, 3))
    site = np.exp(rng_h.normal(0.0, c["site_sd"], size=(ns, 2)))
    rx = GR._child(rng_h)
    sos = GR._coda_filter(c["coda_band"])
    TT_ = GR.T
    for e, (src, tau) in enumerate(zip(sources, taus)):
        src = np.asarray(src, float)
        d = np.linalg.norm(STA - src, axis=1)
        Mt = GR.rand_dc(rng_m)
        fc = float(rx.uniform(*c["fc_range"]))
        phi = float(np.radians(rx.uniform(-c["phase_max_deg"], c["phase_max_deg"])))
        ew = GR.wavelet_energy(fc)
        tt = {b: GR.travel_time_gradient(src, STA, b, c) for b in ("P", "S")}
        for i in range(ns):
            for bi, (b, A) in enumerate((("P", GR.AP), ("S", GR.AS))):
                ta = float(tt[b][i]) + tau + float(st[i, bi])
                arr[e, bi, i] = ta
                vec = GR.radiation_gradient(Mt, src, STA[i], b, c)
                a0 = np.linalg.norm(vec)
                if a0 > 0:
                    w = vec / a0 + dn[i, bi]
                    vec = a0 * site[i, bi] * w / np.linalg.norm(w)
                av = GR.SNR * A * float(np.linalg.norm(vec)) / d[i]
                Dd[e, i] += GR.SNR * A * np.outer(vec, GR.source_wavelet(TT_ - ta, fc, phi)) / d[i]
                GR._add_coda(Dc[e, i], ta, av, ew, c, sos, rx)
    return Dd, Dc, arr


def stream(nc, tag, draw, parts=False):
    """windows of c9_run.stream (same draws) with the event-free twin MFn0 (clutter and noise of the same window, as
    c13_twin) and, with parts, the matched-filtered direct waves and coda of every event and the true arrival times."""
    b = [SEED, tag, nc.ci, nc.j, nc.ri]
    ent = (41, 233, 61)
    rng, rng_m, rng_h = (np.random.default_rng(b + [e]) for e in ent)
    while True:
        srcs, taus = draw(rng)
        cp = [clone(g, b + [e]) for g, e in zip((rng, rng_m, rng_h), ent)]
        cq = [clone(g, b + [e]) for g, e in zip((rng, rng_m, rng_h), ent)]
        D, _p, _s, _a = C8.generate_site(nc.STA, srcs, taus, nc.nclut, rng, rng_m, rng_h, nc.clut_amp, cfg=None,
                                         statics=nc.statics)
        Ds, _p, _s, _a = C8.generate_site(nc.STA, srcs, taus, 0, cp[0], cp[1], cp[2], nc.clut_amp, cfg={"noise": 0.0},
                                          statics=nc.statics)
        MF, sd = xc.matched(D)
        MFs, _sd = xc.matched(Ds)
        w = {"MFn": MF / sd, "MFn0": (MF - MFs) / sd, "truth": [np.asarray(s, float) for s in srcs],
             "taus": [float(t) for t in taus]}
        if parts:
            Dd, Dc, arr = event_parts(nc.STA, srcs, taus, cq[1], cq[2], nc.statics)
            assert np.allclose(Dd.sum(0) + Dc.sum(0), Ds, atol=1e-9 * max(1.0, float(np.abs(Ds).max())))
            w["dir"] = np.array([xc.matched(x)[0] for x in Dd]).reshape(len(srcs), *MF.shape) / sd
            w["coda"] = np.array([xc.matched(x)[0] for x in Dc]).reshape(len(srcs), *MF.shape) / sd
            w["arr"] = arr
        yield w


# ============================================================================================ remedies
def _scale(R, S0, L):
    sc = median_filter(np.abs(R), size=(1, 1, L), mode="nearest") / 0.6745
    return np.maximum(1.0, sc / S0[:, :, None])


def _adaptive(R, pred, tp, kind, L):
    """subtract from R, per station-phase, the best fit of the family (see the module docstring)."""
    ns = R.shape[0]
    for b in range(2):
        for i in range(ns):
            t0 = int(tp[b, i])
            best = None
            u = pred[b, i] / max(np.linalg.norm(pred[b, i]), 1e-12)
            if kind == "fir":
                lo, hi = max(0, t0 - AH - L), min(NT, t0 + AH + L + 1)
                tt = np.arange(lo, hi)
                A = np.zeros((hi - lo, 2 * L + 1))
                for k, l in enumerate(range(-L, L + 1)):
                    j = tt - (t0 + l) + AH
                    ok = (j >= 0) & (j < len(AC))
                    A[ok, k] = AC[j[ok]]
                coef, *_ = np.linalg.lstsq(A, R[i, :, lo:hi].T, rcond=None)
                R[i, :, lo:hi] -= (A @ coef).T
                continue
            for l in range(-L, L + 1):
                tc = t0 + l
                lo, hi = max(0, tc - AH), min(NT, tc + AH + 1)
                a = AC[AH - (tc - lo):AH + (hi - tc)]
                v = R[i, :, lo:hi] @ a / (a @ a)
                if kind == "amp":
                    v = u * float(u @ v)
                gain = float(v @ v) * float(a @ a)
                if best is None or gain > best[0]:
                    best = (gain, lo, hi, a, v)
            _g, lo, hi, a, v = best
            R[i, :, lo:hi] -= np.outer(v, a)


# ============================================================================================ the loop
def mt_loop14(E, MFn, kmax, cfg=None, capv=C.CAPV, min_sep=5.0, diag=None, ora=None):
    """c10_lib.mt_loop (mode hard, W = SH, no refinement) with the remedies of cfg.  cfg None or {}: identical.
    diag: {'dir', 'coda' (K, ns, 3, NT), 'arr' (K, 2, ns), 'cn' (ns, 3, NT), 'truth', 'taus'} -> every atom gets the
          decomposition of its statistic (key 'dg').
    ora:  {'ev' (K, ns, 3, NT) what to remove of an event, 'truth'}: ORACLE deflation; an accepted atom within 5 km of
          an event not yet removed (the nearest) removes that event's 'ev' instead of the model fit (diagnosis only)."""
    cfg = cfg or {}
    r1, r2, r3, r4 = cfg.get("r1"), cfg.get("r2"), cfg.get("r3"), cfg.get("r4")
    ns = E.ns
    R = MFn.copy()
    excluded = np.zeros(E.N, bool)
    ex_nt = np.zeros((E.N, E.T), bool) if r4 else None
    Wt = np.ones((ns, NT)) if r1 else None
    S0 = np.maximum(np.median(np.abs(MFn), axis=2) / 0.6745, 1e-12) if r2 else None
    atoms, tps = [], []
    removed = set()
    for it in range(kmax):
        if excluded.all():
            break
        Rw = R
        if r1 and it > 0:
            Rw = Rw * Wt[:, None, :]
        if r2 and (it > 0 or r2[1]):
            Rw = Rw / _scale(R, S0, r2[0])
        Eh = C._half_energy(Rw, capv)
        tpk, _pen = C.choose(E, Eh, "hard", None, C.SH, R=Rw)
        vf = C._vectors(Rw, tpk, capv).reshape(E.N, E.T, 6 * ns)
        y = np.einsum("gtr,grk->gtk", vf, E.Q)
        sc = (y * y).sum(-1) / 2.0
        sc[excluded] = -np.inf
        if r4:
            sc[ex_nt] = -np.inf
        gi = int(np.argmax(sc.max(axis=1)))
        tj = int(np.argmax(sc[gi]))
        s = float(sc[gi, tj])
        if not np.isfinite(s):
            break
        src, tau, tp, G = E.XN[gi], float(E.taus[tj]), tpk[gi, tj].copy(), E.G[gi]
        a = {"s": s, "s0": s, "x": np.asarray(src, float), "tau": tau}
        if diag is not None:
            a["dg"] = _decompose(E, gi, vf[gi, tj], y[gi, tj], Rw, R, MFn, tp, tps, atoms, diag, capv)
        atoms.append(a)
        tps.append(tp)
        dxy = np.linalg.norm(E.XN[:, :2] - src[None, :2], axis=1)
        excluded = excluded | (dxy <= min_sep)
        if r4:
            ex_nt |= (dxy <= r4[0])[:, None] & (np.abs(E.taus - tau) <= r4[1] + 1e-9)[None, :]
        vsel = C._vectors(R, tp[None], capv).reshape(6 * ns)
        m, *_ = np.linalg.lstsq(G, vsel, rcond=None)
        pred = (G @ m).reshape(2, ns, 3)
        done = False
        if ora is not None and len(ora["truth"]):
            dd = np.array([np.linalg.norm(src[:2] - t[:2]) for t in ora["truth"]])
            j = int(dd.argmin())
            if dd[j] <= 5.0 and j not in removed:
                R -= ora["ev"][j]
                removed.add(j)
                done = True
        if not done:
            if r3:
                _adaptive(R, pred, tp, r3[0], r3[1])
            else:
                C._deflate(R, pred, tp)
        if r1:
            Wc = r1[1]
            for b in range(2):
                for i in range(ns):
                    lo, hi = max(0, int(tp[b, i]) - Wc), min(NT, int(tp[b, i]) + Wc + 1)
                    if r1[0] == "hard":
                        Wt[i, lo:hi] = 0.0
                    else:
                        rho = min(1.0, float((R[i, :, lo:hi] ** 2).sum()) / max(float((MFn[i, :, lo:hi] ** 2).sum()),
                                                                                1e-12))
                        Wt[i, lo:hi] = np.minimum(Wt[i, lo:hi], rho ** r1[2])
    return atoms


def _decompose(E, gi, v, y, Rw, R, MFn, tp, tps, atoms, dg, capv):
    """linear split of s0 = v . P v / 2 over the parts of the residual at the read samples (only meaningful for the
    baseline configuration, where Rw is R)."""
    ns = E.ns
    Pv = (E.Q[gi] @ y).reshape(2, ns, 3)
    st = np.arange(ns)[None, :]
    raw = Rw.transpose(0, 2, 1)[st, tp]                                            # (2, ns, 3)
    f = np.minimum(1.0, capv / np.maximum(np.linalg.norm(raw, axis=-1), 1e-12))    # cap factor (2, ns)
    c_bi = 0.5 * (Pv * v.reshape(2, ns, 3)).sum(-1)                                # contributions, sum = s0
    K = len(dg["truth"])
    arr = dg["arr"]
    t_read = tp * DT
    sh = lambda Z: 0.5 * f * (Pv * Z.transpose(0, 2, 1)[st, tp]).sum(-1)           # share (2, ns) of a part  # noqa
    s_cn = sh(dg["cn"])
    s_defl = -sh(MFn - R)
    live_same, live_cross, resid_e, coda_e = np.zeros(K), np.zeros(K), np.zeros(K), np.zeros(K)
    f_arr = f_same = f_cross = 0.0
    near_ev = -np.ones((2, ns), int)
    for e in range(K):
        sd_e, sc_e = sh(dg["dir"][e]), sh(dg["coda"][e])
        coda_e[e] = sc_e.sum()
        for b in range(2):
            for i in range(ns):
                bp = int(np.argmin(np.abs(t_read[b, i] - arr[e, :, i])))
                claimed = any(np.min(np.abs(q[:, i] * DT - arr[e, bp, i])) <= NEAR for q in tps)
                if claimed:
                    resid_e[e] += sd_e[b, i]
                elif bp == b:
                    live_same[e] += sd_e[b, i]
                else:
                    live_cross[e] += sd_e[b, i]
    # literal measure of the brief: share of the statistic read within NEAR of a true arrival of an event
    acc = set()
    for q in atoms:
        if K:
            dd = np.array([np.linalg.norm(q["x"][:2] - t[:2]) for t in dg["truth"]])
            if dd.min() <= 5.0:
                acc.add(int(dd.argmin()))
    f_acc = 0.0
    for b in range(2):
        for i in range(ns):
            if not K:
                continue
            dtm = np.abs(t_read[b, i] - arr[:, :, i])                              # (K, 2)
            e, bp = np.unravel_index(int(dtm.argmin()), dtm.shape)
            if dtm[e, bp] <= NEAR:
                f_arr += c_bi[b, i]
                near_ev[b, i] = e
                if bp == b:
                    f_same += c_bi[b, i]
                else:
                    f_cross += c_bi[b, i]
                if e in acc:
                    f_acc += c_bi[b, i]
    s0 = float(c_bi.sum())
    tot = float(s_cn.sum() + s_defl.sum() + live_same.sum() + live_cross.sum() + resid_e.sum() + coda_e.sum())
    assert abs(tot - s0) <= 1e-6 * max(1.0, abs(s0)), (tot, s0)
    return {"s0": s0, "cn": float(s_cn.sum()), "defl": float(s_defl.sum()), "live_same": live_same.tolist(),
            "live_cross": live_cross.tolist(), "resid": resid_e.tolist(), "coda": coda_e.tolist(),
            "f_arr": float(f_arr), "f_same": float(f_same), "f_cross": float(f_cross), "f_acc": float(f_acc),
            "acc": sorted(acc)}


# ============================================================================================ records
def pack(atoms):
    return [[float(a["x"][0]), float(a["x"][1]), float(a["x"][2]), float(a["tau"]), float(a["s"])] for a in atoms]


def unpack(lst):
    return [{"x": np.array(a[:3]), "tau": a[3], "s": a[4]} for a in lst]


def same(a, b):
    return len(a) == len(b) and all(np.array_equal(p["x"], q["x"]) and p["tau"] == q["tau"] and p["s"] == q["s"]
                                    for p, q in zip(a, b))
