"""c18j_loeo_bh.py - is the C18j gate reading (layered medium, +-2 s terms, 2 Hz, NOCAP; 5/16 on U8, 4/16 on N8 with
LOC+T at 30 false/h) helped by estimating the station terms on the SAME 16 CAL events that are then scored?
Re-run of c18j_layered.py detect with one change: for every 40 s window the layered station terms are re-estimated
from the CAL calibration events whose catalogue origin is NOT inside [ts - 1, ts + 41) s (so no event that the window
can contain, and in particular no event that the window can hit, contributes to the terms used to detect it).
Everything else is c18j_layered.main (same Net, table, engine, loop, owned range).
    python c18j_loeo_bh.py SUBNET VARIANTS OUT.json [narrow]
narrow: exclude only the calibration events that the window can HIT (origin in [ts + 4.5, ts + 8.5) s).
Per-event readings are computed once (same peak_delays / layered times as c18j_layered.corrections)."""
import copy
import json
import sys
import time

import numpy as np

sys.path[:0] = ["/home/claude/work/r0", "/home/claude/work/explore_contrib", "/home/claude/work"]
import r0_lib as R  # noqa: E402
import r0_run as RR  # noqa: E402
import c18j_layered as J  # noqa: E402


def main(sub, variants, outp):
    import c10_lib as C
    import c18d_detect as D
    import xc
    t0 = time.time()
    net = RR.Net(sub, 940.0)
    tab = J.tt_table(sorted({float(z) for z in net.eng["SC"][0].XN[:, 2]}))
    E0, corr_all = J.layered_engine(net, tab)
    tt = J.tt_nodes(E0.XN, E0.STA, tab)
    ev_t = np.array([e["t"] for e in net.calev if e["t"] - 6.5 >= 0])
    assert len(ev_t) == len(net.calwin)
    calwin_all = list(net.calwin)
    engines = {}
    narrow = len(sys.argv) > 4 and sys.argv[4] == "narrow"
    from station_corrections import peak_delays
    from task0_common import VP, VS
    import c18_tt as TT
    DS, EN = [], []
    for D_, src, tau in calwin_all:
        dl, en = peak_delays(net.STA, D_, np.asarray(src, float), float(tau), 2.0)
        d = np.linalg.norm(net.STA - np.asarray(src), axis=1)
        t_hom = np.stack([d / VP, d / VS], axis=1)
        z = float(src[2])
        t_lay = J.tt_nodes(np.asarray(src)[None], net.STA, {z: TT.table(TT.HK, J.R_VPVS, [z])[z]})[0].T
        DS.append(dl + t_hom - t_lay)
        EN.append(en)
    DS, EN = np.array(DS), np.array(EN)

    def corr_of(keep):
        ds, ok = DS[list(keep)], EN[list(keep)] >= 25.0
        c = np.zeros(ds.shape[1:])
        for i in range(ds.shape[1]):
            for b in range(2):
                if ok[:, i, b].any():
                    c[i, b] = float(np.median(ds[ok[:, i, b], i, b]))
        return c
    assert np.allclose(corr_of(range(len(calwin_all))), corr_all)

    def engine_for(ts):
        ex = (ev_t >= ts + 4.5) & (ev_t < ts + 8.5) if narrow else (ev_t >= ts - 1.0) & (ev_t < ts + 41.0)
        keep = tuple(int(k) for k in np.where(~ex)[0])
        if keep not in engines:
            corr = corr_of(keep)
            E = copy.copy(E0)
            x = (tt[:, None] + E.taus[None, :, None, None] + corr.T[None, None]) / R.DT
            E.IDX = np.clip(x.astype(int), 0, R.NT - 1)
            E.frac = x - E.IDX
            E.cache = {}
            engines[keep] = (E, corr)
        return engines[keep]

    res = {"subnet": sub, "corr_all": corr_all.tolist(), "atoms": {v: [] for v in variants}, "n_termsets": 0}
    ws = RR.starts("cal")
    for wi, ts in enumerate(ws):
        E, _corr = engine_for(ts)
        MF, sd = xc.matched(R.window(net.d16, ts))
        MFn = MF / sd
        for v in variants:
            at = C.bomp_loop(E, None, MFn, R.KMAX) if v == "BOMP" else C.mt_loop(E, MFn, R.KMAX, capv=D.CAPS[v])
            for a in at:
                if R.OWN[0] - 1e-9 <= a["tau"] < R.OWN[1] - 1e-9:
                    x = a["x"]
                    res["atoms"][v].append([float(x[0]), float(x[1]), float(x[2]), float(ts + a["tau"]), float(a["s"]), wi])
        if wi % 50 == 0:
            print(time.strftime("%H:%M:%S"), sub, wi, len(ws), "term sets", len(engines), flush=True)
    res["n_termsets"] = len(engines)
    res["corr_sets"] = {",".join(map(str, k)): v[1].tolist() for k, v in engines.items()}
    json.dump(res, open(outp, "w"))
    print("DONE", sub, round(time.time() - t0), "s; term sets", len(engines), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2].split(","), sys.argv[3])
