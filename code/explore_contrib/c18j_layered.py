"""
c18j_layered.py - exploration C18j (Ridgecrest CAL only).  After C18i showed that, once the peak search is wide enough
(+-2 s), the layered model of Hadley and Kanamori predicts the held-out readings better than the registered
homogeneous medium, the detector is run with that medium: travel times of every coarse node to every station from the
layered model (c18_tt, first arrivals, VP/VS 1.73), catalogue station corrections re-estimated against the layered
times with a +-2 s search, and the cap 6 or LOG6 statistic, with the 2 Hz or 3 Hz matched filter.  Only the timing
table of the engine changes (sample index and fraction); radiation operators, subspaces and the loop are those of R0.
    python c18j_layered.py oracle FC SUBNET OUT.json
    python c18j_layered.py detect FC SUBNET VARIANTS OUT.json
"""
import json
import sys
import time

import numpy as np

sys.path.insert(0, "/home/claude/work/r0")
import r0_lib as R  # noqa: E402
import r0_run as RR  # noqa: E402

import c18_tt as TT  # noqa: E402

R_VPVS = 1.73


def tt_table(depths):
    tab = TT.table(TT.HK, R_VPVS, depths, xmax=120.0, dx=0.25)
    return tab


def tt_nodes(XN, STA, tab):
    """(N, 2, ns) layered first-arrival times, linear interpolation in epicentral distance at the node depth."""
    XN, STA = np.asarray(XN, float), np.asarray(STA, float)
    out = np.zeros((len(XN), 2, len(STA)))
    for z, (xs, tp, ts) in tab.items():
        m = np.abs(XN[:, 2] - z) < 1e-6
        if not m.any():
            continue
        X = np.hypot(XN[m, None, 0] - STA[None, :, 0], XN[m, None, 1] - STA[None, :, 1])
        out[m, 0] = np.interp(X, xs, tp)
        out[m, 1] = np.interp(X, xs, ts)
    assert np.all(out > 0), "node depth not in the table"
    return out


def corrections(net, tab, search=2.0, gate=25.0):
    """median over CAL events of (observed energy peak - layered time), peaks searched +-search s around the
    homogeneous prediction (station_corrections.peak_delays), energy gate 25 (as estimate_station_corrections)."""
    from station_corrections import peak_delays
    from task0_common import VP, VS
    ds, ens = [], []
    for D, src, tau in net.calwin:
        dl, en = peak_delays(net.STA, D, np.asarray(src, float), float(tau), search)
        d = np.linalg.norm(net.STA - np.asarray(src), axis=1)
        t_hom = np.stack([d / VP, d / VS], axis=1)
        t_lay = tt_nodes(np.asarray(src)[None], net.STA, {float(src[2]): TT.table(TT.HK, R_VPVS, [float(src[2])])[float(src[2])]})[0].T
        ds.append(dl + t_hom - t_lay)
        ens.append(en)
    ds, ens = np.array(ds), np.array(ens)
    ok = ens >= gate
    corr = np.zeros(ds.shape[1:])
    for i in range(ds.shape[1]):
        for b in range(2):
            if ok[:, i, b].any():
                corr[i, b] = float(np.median(ds[ok[:, i, b], i, b]))
    return corr


def layered_engine(net, tab):
    E = net.eng["SC"][0]
    corr = corrections(net, tab)
    tt = tt_nodes(E.XN, E.STA, tab)
    x = (tt[:, None] + E.taus[None, :, None, None] + corr.T[None, None]) / R.DT
    E.IDX = np.clip(x.astype(int), 0, R.NT - 1)
    E.frac = x - E.IDX
    E.cache = {}
    return E, corr


def main(mode, fc, sub, rest):
    if fc != 2.0:
        import c18h_kernel as K
        K.set_kernel(fc)
    import c10_lib as C
    import xc
    t0 = time.time()
    net = RR.Net(sub, 940.0)
    tab = tt_table(sorted({float(z) for z in net.eng["SC"][0].XN[:, 2]}))
    E, corr = layered_engine(net, tab)
    print("layered corrections (s), P / S:", np.round(corr[:, 0], 2).tolist(), np.round(corr[:, 1], 2).tolist(), flush=True)
    if mode == "oracle":
        import c18c_oracle as OC
        OC.VARIANTS = ("CAP6", "LOG6", "NOCAP")
        RR_Net = RR.Net
        RR.Net = lambda s, t: net                     # reuse the network with the layered engine
        try:
            OC.main(sub, rest[0])
        finally:
            RR.Net = RR_Net
        return
    import c18d_detect as D
    variants, outp = rest[0].split(","), rest[1]
    res = {"subnet": sub, "fc": fc, "corr": corr.tolist(), "atoms": {v: [] for v in variants}, "sec": {v: 0.0 for v in variants}}
    ws = RR.starts("cal")
    for wi, ts in enumerate(ws):
        MF, sd = xc.matched(R.window(net.d16, ts))
        MFn = MF / sd
        for v in variants:
            t1 = time.time()
            at = C.bomp_loop(E, None, MFn, R.KMAX) if v == "BOMP" else C.mt_loop(E, MFn, R.KMAX, capv=D.CAPS[v])
            for a in at:
                if R.OWN[0] - 1e-9 <= a["tau"] < R.OWN[1] - 1e-9:
                    x = a["x"]
                    res["atoms"][v].append([float(x[0]), float(x[1]), float(x[2]), float(ts + a["tau"]), float(a["s"]), wi])
            res["sec"][v] += time.time() - t1
        if wi % 50 == 0:
            print(time.strftime("%H:%M:%S"), sub, wi, len(ws), flush=True)
    json.dump(res, open(outp, "w"))
    print("DONE", sub, time.time() - t0, flush=True)


if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4:])
