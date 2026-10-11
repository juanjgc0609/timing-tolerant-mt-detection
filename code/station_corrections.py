"""
station_corrections.py - Component C: station time corrections from catalogued events.

MT-FW searches the peak of each station/phase within +-0.25 s of the predicted time; unmodelled
station delays beyond that window break it (Appendix Y.11: delays with sd 0.4 s dropped MT-FW in
the one-sided cell from 0.78 to 0.37 recall at 0.1 false events per window).  Here the delays are
a property of the network (fixed across windows, as real station terms are), estimated from K
catalogued events with known location and origin time, and applied as corrections to the
predicted P and S times of every hypothesis.

estimate_station_corrections(STA, cal, ...):  for every catalogue window k, station i and phase b,
    the 3-component matched-filter energy E_i(t) is searched within +-1 s of the predicted time
    d_ik / V_b + tau_k; the peak (parabolic sub-sample refinement) gives a delay; peaks below an
    energy gate (default 25, i.e. 5 noise sd) are discarded; the correction is the median of the
    accepted delays (0 when no event qualifies).

generate_c: generate_y (pert) with the station statics passed in (fixed per network) instead of
    drawn per window; direction and site perturbations and clutter remain per window.

    python station_corrections.py [out.txt]   # tables of station_corrections_report.txt
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import task1x_mtglrt as X  # noqa: E402
from certified_fw import certified_fw_record, dense_fw_record_c  # noqa: E402
from froc_exact import exact_HF  # noqa: E402
from task0_common import AP, AS, DEFAULT_TAUS, DT, NT, SNR, T, VP, VS, default_grid, make_network, ricker  # noqa: E402
from task1o_consistent import draw_sources, matched  # noqa: E402
from task1r_fair import SCEN, recall_at  # noqa: E402
from task1w_mech import CP, CS, rand_dc  # noqa: E402

STATIC_SD_NET = 0.4
CLUT_AMP = 0.35


def generate_c(STA, sources, taus, n_clutter, rng, rng_m, rng_h, clut_amp, statics):
    """task1y_clutter.generate_y with FIXED station statics (ns, 2) instead of per-window draws."""
    ns = len(STA)
    D = np.zeros((ns, 3, NT))
    dn = rng_h.normal(0.0, X.DIR_SD, size=(ns, 2, 3))
    site = np.exp(rng_h.normal(0.0, X.SITE_SD, size=(ns, 2)))
    for src, tau in zip(sources, taus):
        d = np.linalg.norm(STA - src, axis=1)
        M = rand_dc(rng_m)
        for i in range(ns):
            n = (STA[i] - src) / d[i]
            for bi, (b, A, V) in enumerate((("P", AP, VP), ("S", AS, VS))):
                ta = d[i] / V + tau + float(statics[i, bi])
                if b == "P":
                    vec = n * float(n @ M @ n) / np.sqrt(CP)
                else:
                    u = M @ n
                    vec = (u - float(u @ n) * n) / np.sqrt(CS)
                a0 = np.linalg.norm(vec)
                if a0 > 0:
                    w = vec / a0 + dn[i, bi]
                    vec = a0 * site[i, bi] * w / np.linalg.norm(w)
                D[i] += SNR * A * np.outer(vec, ricker(T - ta)) / d[i]
    ref = np.asarray(sources[0], float) if len(sources) else np.zeros(3)
    scale = float(np.median(np.linalg.norm(STA - ref, axis=1)))
    for _ in range(n_clutter):
        i = int(rng.integers(ns))
        tc = float(rng.uniform(2, T[-1] - 4))
        q = rng.normal(size=3)
        q /= np.linalg.norm(q)
        D[i] += SNR * clut_amp * np.outer(q, ricker(T - tc)) / scale
    D += rng.normal(scale=SNR * 0.02, size=D.shape)
    return D, [np.asarray(s, float) for s in sources]


MIN_PS = 6        # samples (0.375 s): the P peak must precede the S peak by at least this much


def peak_delays(STA, D, src, tau, search=1.0, joint=True):
    """per (station, phase): (delay of the matched-filter energy peak w.r.t. the predicted time,
    peak energy).  joint=True: the P and S peaks are chosen together, maximising E(kP) + E(kS) over
    the two +-search windows subject to kP + MIN_PS <= kS (P precedes S).  Searched independently,
    the P window of a station close to the event (S-P < 2 s) contains the stronger S arrival and
    the P delay is wrong (seen on a 4-station network: 7 of 13 P delays at one station were S).
    Parabolic refinement of each peak sample."""
    MF, sd = matched(D)
    E = np.sum((MF / sd) ** 2, axis=1)                                        # (ns, NT)
    d = np.linalg.norm(STA - src, axis=1)
    W = int(round(search / DT))
    out = np.zeros((len(STA), 2))
    en = np.zeros((len(STA), 2))
    for i in range(len(STA)):
        tp = [d[i] / VP + tau, d[i] / VS + tau]
        win = []
        for b in range(2):
            c = int(round(tp[b] / DT))
            win.append(np.arange(max(1, c - W), min(NT - 2, c + W) + 1))
        if joint:
            sc = E[i, win[0]][:, None] + E[i, win[1]][None, :]
            sc = np.where(win[1][None, :] >= win[0][:, None] + MIN_PS, sc, -np.inf)
            a, b_ = np.unravel_index(int(np.argmax(sc)), sc.shape)
            ks = (int(win[0][a]), int(win[1][b_]))
        else:
            ks = tuple(int(w[np.argmax(E[i, w])]) for w in win)
        for b in range(2):
            k = ks[b]
            e0, e1, e2 = E[i, k - 1], E[i, k], E[i, k + 1]
            den = e0 - 2 * e1 + e2
            dk = 0.5 * (e0 - e2) / den if den < 0 else 0.0
            out[i, b] = (k + float(np.clip(dk, -0.5, 0.5))) * DT - tp[b]
            en[i, b] = e1
    return out, en


def catalogue(STA, K, statics, rng_seed, n_clutter, gate=25.0, min_arr=3, cat_err=(0.0, 0.0)):
    """K catalogued events on this network: an event enters the catalogue if at least min_arr of
    its station/phase arrivals exceed the energy gate at their (true) peak (a detection criterion);
    its catalogue location and origin time are the true ones plus optional errors
    (cat_err = (location sd km, origin-time sd s))."""
    rng = np.random.default_rng(rng_seed)
    rng_m = np.random.default_rng(rng_seed + 1)
    rng_h = np.random.default_rng(rng_seed + 2)
    rng_e = np.random.default_rng(rng_seed + 3)
    out, drawn = [], 0
    while len(out) < K:
        drawn += 1
        src = np.array([rng.uniform(-15, 15), rng.uniform(-15, 15), rng.uniform(5, 13)])
        tau = float(rng.uniform(5.5, 7.5))
        D, _ = generate_c(STA, [src], [tau], n_clutter, rng, rng_m, rng_h, CLUT_AMP, statics)
        _dl, en = peak_delays(STA, D, src, tau)
        if (en >= gate).sum() < min_arr:
            continue
        src_c = src + rng_e.normal(0, cat_err[0], 3) * np.array([1, 1, 1])
        src_c[2] = float(np.clip(src_c[2], 1.0, 20.0))
        out.append((D, src_c, tau + float(rng_e.normal(0, cat_err[1]))))
    return out, drawn


def estimate_station_corrections(STA, cal, search=1.0, gate=25.0, joint=True):
    """Median over catalogue events of the per-(station, phase) delays whose peak energy passes the
    gate.  cal: list of (D, source xyz, origin time).  Returns corrections (ns, 2) and the number
    of accepted delays per (station, phase)."""
    ds, ens = [], []
    for D, src, tau in cal:
        dl, en = peak_delays(STA, D, np.asarray(src, float), float(tau), search, joint)
        ds.append(dl)
        ens.append(en)
    ds, ens = np.array(ds), np.array(ens)
    ok = ens >= gate
    corr = np.zeros(ds.shape[1:])
    nacc = ok.sum(0)
    for i in range(ds.shape[1]):
        for b in range(2):
            if nacc[i, b]:
                corr[i, b] = float(np.median(ds[ok[:, i, b], i, b]))
    return corr, nacc


def froc_recall(recs, nte, caps=(0.1, 0.25)):
    """recall at realized false events per window <= cap, identical to froc_exact.exact_HF followed
    by task1r_fair.recall_at (records are prefix-monotone per window), vectorised for bootstrap."""
    sc, hit, key = [], [], []
    for w, r in enumerate(recs):
        for s_, h in r:
            sc.append(s_)
            hit.append(h)
            key.append(w * 16 + h if h >= 0 else -1)
    if not sc:
        return [0.0 for _ in caps]
    sc, hit, key = np.array(sc), np.array(hit), np.array(key)
    o = np.argsort(-sc, kind="stable")
    sc, hit, key = sc[o], hit[o], key[o]
    newhit = np.zeros(len(sc), bool)
    seen = set()
    for j in range(len(sc)):
        if key[j] >= 0 and key[j] not in seen:
            seen.add(key[j])
            newhit[j] = True
    H = np.cumsum(newhit)
    F = np.cumsum(hit < 0)
    last = np.r_[np.where(np.diff(sc) != 0)[0], len(sc) - 1]      # end of each tie group
    H, F = H[last], F[last]
    out = []
    for c in caps:
        ok = F / len(recs) <= c
        out.append(float(H[ok].max() / max(1, nte)) if ok.any() else 0.0)
    return out


def boot_diff(recsA, recsB, ntrue_w, caps=(0.1, 0.25), B=400, seed=0):
    """paired bootstrap over windows of recall(A) - recall(B) at each cap: (point, lo, hi)."""
    rng = np.random.default_rng(seed)
    n = len(recsA)
    pt = np.subtract(froc_recall(recsA, sum(ntrue_w), caps), froc_recall(recsB, sum(ntrue_w), caps))
    bs = []
    for _ in range(B):
        ix = rng.integers(0, n, n)
        nte = int(sum(ntrue_w[i] for i in ix))
        a = froc_recall([recsA[i] for i in ix], nte, caps)
        b = froc_recall([recsB[i] for i in ix], nte, caps)
        bs.append(np.subtract(a, b))
    bs = np.array(bs)
    return [(float(pt[k]), float(np.quantile(bs[:, k], 0.025)), float(np.quantile(bs[:, k], 0.975))) for k in range(len(caps))]


def eval_windows(STA, statics, n, seed, n_clutter):
    rng = np.random.default_rng(seed + 41)
    rng_m = np.random.default_rng(seed + 233)
    rng_h = np.random.default_rng(seed + 61)
    out = []
    for _ in range(n):
        srcs, taus = draw_sources(rng)
        D, tr = generate_c(STA, srcs, taus, n_clutter, rng, rng_m, rng_h, CLUT_AMP, statics)
        out.append((D, tr))
    return out


def main():
    t0 = time.time()
    NEVAL = int(os.environ.get("SC_NEVAL", "40"))
    NETS = int(os.environ.get("SC_NETS", "4"))
    grid = default_grid()
    garr = np.asarray(grid)
    lines = []
    P = lambda s="": (print(s, flush=True), lines.append(s))  # noqa: E731
    P("STATION TIME CORRECTIONS FROM CATALOGUED EVENTS (Component C) - station_corrections.py")
    P("=" * 110)
    P(f"Statics N(0, {STATIC_SD_NET} s) per station and phase, FIXED per network; pert direction (~10 deg) and site")
    P(f"(lognormal sd 0.5) factors per window; weak clutter (amplitude {CLUT_AMP}, n_clutter of task1r_fair.SCEN).")
    P("Catalogue: events uniform in [-15,15]^2 x [5,13] km, origin 5.5-7.5 s, random double couple; an event is")
    P("catalogued if >= 3 station/phase arrivals exceed the energy gate (25 = 5 noise sd at the peak).")
    P("Estimate: 3-C matched-filter energy peak within +-1 s of the predicted time, parabolic refinement,")
    P("gate 25, median over events.  MT-FW = fw_record (via the certified oracle, identical records) with the")
    P(f"corrections added to the predicted P and S times.  {NEVAL} evaluation windows per network, fresh sources")
    P("(draw_sources: 0-2 events), same network and statics.  Recall at realized false events per window.")
    P("")
    variants = ["none", "K=5", "K=10", "K=20", "oracle"]
    pooled = {}
    resid_rows = []
    checked = False
    for scen in ("ladder_4", "ladder_6", "unilateral_c40"):
        kind, nst, ncl = SCEN[scen]
        for net in range(NETS):
            nseed = 77000 + 100 * net + {"ladder_4": 0, "ladder_6": 10, "unilateral_c40": 20}[scen]
            STA = make_network(kind, nst, seed=nseed)
            statics = np.random.default_rng(nseed + 5).normal(0.0, STATIC_SD_NET, size=(nst, 2))
            Gs = [X.radiation_operator(g, STA) for g in garr]
            Qs = [X.basis(G) for G in Gs]
            cat, drawn = catalogue(STA, 20, statics, nseed + 11, ncl)
            corr = {"none": None, "oracle": statics}
            for K in (5, 10, 20):
                c, nacc = estimate_station_corrections(STA, cat[:K])
                corr[f"K={K}"] = c
                res = c - statics
                resid_rows.append((scen, net, K, float(np.sqrt(np.mean(res ** 2))),
                                   float(np.sqrt(np.mean((res - res.mean()) ** 2))), float(np.mean(np.abs(res) > 0.25)),
                                   float((nacc == 0).mean())))
            # ablation: independent P and S searches (K = 20)
            ci, ni = estimate_station_corrections(STA, cat, joint=False)
            resi = ci - statics
            resid_rows.append((scen, net, "20 indep", float(np.sqrt(np.mean(resi ** 2))),
                               float(np.sqrt(np.mean((resi - resi.mean()) ** 2))), float(np.mean(np.abs(resi) > 0.25)),
                               float((ni == 0).mean())))
            # catalogue with location / origin-time errors (sensitivity, K = 20)
            cat_e, _ = catalogue(STA, 20, statics, nseed + 11, ncl, cat_err=(1.0, 0.1))
            ce, _n = estimate_station_corrections(STA, cat_e)
            rese = ce - statics
            resid_rows.append((scen, net, "20+err", float(np.sqrt(np.mean(rese ** 2))),
                               float(np.sqrt(np.mean((rese - rese.mean()) ** 2))), float(np.mean(np.abs(rese) > 0.25)),
                               float((_n == 0).mean())))
            ev = eval_windows(STA, statics, NEVAL, nseed + 50000, ncl)
            nte = sum(len(tr) for _D, tr in ev)
            recs = {v: [] for v in variants}
            for w, (D, tr) in enumerate(ev):
                MF, sd = matched(D)
                MFn = MF / sd
                for v in variants:
                    tc = None if corr[v] is None else np.asarray(corr[v], float)
                    r = certified_fw_record(STA, MFn, grid, DEFAULT_TAUS, Gs, Qs, tr, mode="enclosure", tcorr=tc)
                    recs[v].append(r)
                if not checked and w < 2:      # the corrected oracle against the dense corrected loop
                    for v in ("K=10", "oracle"):
                        rd = dense_fw_record_c(STA, MFn, grid, DEFAULT_TAUS, Gs, Qs, tr, tcorr=np.asarray(corr[v], float))
                        assert rd == recs[v][-1], "certified oracle with corrections differs from the dense loop"
                    rd = X.fw_record(STA, MFn, grid, DEFAULT_TAUS, Gs, Qs, tr)
                    assert rd == recs["none"][-1]
            checked = True
            row = []
            for v in variants:
                _t, H, F = exact_HF(recs[v])
                pooled.setdefault((scen, v), []).extend(recs[v])
                row.append((recall_at(H, F, nte, 0.1), recall_at(H, F, nte, 0.25)))
            pooled.setdefault((scen, "nte"), []).append(nte)
            pooled.setdefault((scen, "ntw"), []).extend(len(tr) for _D, tr in ev)
            for v in variants:                      # the vectorised readout equals exact_HF + recall_at
                _t, H, F = exact_HF(recs[v])
                assert np.allclose(froc_recall(recs[v], nte), [recall_at(H, F, nte, 0.1), recall_at(H, F, nte, 0.25)])
            P(f"{scen:15s} network {net} (seed {nseed}; |static| > 0.25 s in {np.mean(np.abs(statics) > 0.25):.0%} of "
              f"station-phases; catalogue drew {drawn} events for 20 catalogued; {nte} true events in evaluation): "
              + " | ".join(f"{v} {a:.2f}/{b:.2f}" for v, (a, b) in zip(variants, row)))
    P("")
    P("TABLE C1. Residual static after correction (estimated - true), per network; rms over the 2 ns station-phases.")
    P(f"{'cell':15s} {'net':>3s} {'K':>7s} {'rms (s)':>8s} {'rms, common':>12s} {'|res|>0.25 s':>13s} {'no estimate':>12s}")
    P(f"{'':15s} {'':>3s} {'':>7s} {'':>8s} {'offset out':>12s}")
    for scen, net, K, r, r0, f25, f0 in resid_rows:
        P(f"{scen:15s} {net:3d} {str(K):>7s} {r:8.3f} {r0:12.3f} {f25:13.2f} {f0:12.2f}")
    P(f"(uncorrected: rms = sd of the statics, {STATIC_SD_NET} s in expectation; '20+err' = K = 20 with catalogue location")
    P(" errors of 1 km per coordinate and origin-time errors of 0.1 s; '20 indep' = K = 20 with independent P and S")
    P(" searches instead of the joint P-before-S search, the ablation that motivated the joint search.)")
    P("")
    P("TABLE C2. MT-FW event recall at realized false events per window <= 0.1 / <= 0.25, pooled over the networks")
    P(f"of each cell ({NETS} x {NEVAL} windows).")
    P(f"{'cell':15s} {'events':>6s} " + " ".join(f"{v:>11s}" for v in variants))
    for scen in ("ladder_4", "ladder_6", "unilateral_c40"):
        nte = sum(pooled[(scen, "nte")])
        cells = []
        for v in variants:
            _t, H, F = exact_HF(pooled[(scen, v)])
            cells.append(f"{recall_at(H, F, nte, 0.1):.2f}/{recall_at(H, F, nte, 0.25):.2f}")
        P(f"{scen:15s} {nte:6d} " + " ".join(f"{c:>11s}" for c in cells))
    P("")
    P("TABLE C3. Paired differences of pooled recall at <= 0.1 / <= 0.25 false events per window, with 95% intervals")
    P("from a paired bootstrap over windows (400 resamples).")
    pairs = [("K=5", "none"), ("K=10", "none"), ("K=20", "none"), ("oracle", "none"), ("oracle", "K=20")]
    P(f"{'cell':15s} " + " ".join(f"{a + ' - ' + b:>34s}" for a, b in pairs))
    for scen in ("ladder_4", "ladder_6", "unilateral_c40"):
        cells = []
        for a, b in pairs:
            (p1, l1, h1), (p2, l2, h2) = boot_diff(pooled[(scen, a)], pooled[(scen, b)], pooled[(scen, "ntw")])
            cells.append(f"{p1:+.2f} [{l1:+.2f},{h1:+.2f}] {p2:+.2f} [{l2:+.2f},{h2:+.2f}]")
        P(f"{scen:15s} " + " ".join(f"{c:>34s}" for c in cells))
    P(f"(run time {time.time() - t0:.0f} s)")
    return lines


if __name__ == "__main__":
    out = main()
    if len(sys.argv) > 1:                       # optional: write the tables to the given path
        with open(sys.argv[1], "w") as fh:
            fh.write("\n".join(out) + "\n")
