"""C4: false-alarm control in windows WITH events. Per network (seed 999): NE standard windows (TAG 3004, 0-2 events),
the network's 20 catalogue records (1 event each, catalogue location with 1 km error), NEMP empty windows (TAG 3005).
Methods with locations: BOMP-PK+SC, TS-SFW+SC (f 0.01), MT-FW+SC, QM-3C+SC (frozen 2A config, corrections in LUT).
Thresholds are computed in c4_analyze.py."""
import json, os, sys, tempfile, time
import numpy as np
import xc
import qm_wrapper as QM

CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
T = xc.T


def no_sources(rng):
    return [], []


def qm_rec_locs(events, truth, min_sep_km=5.0, max_events=5):
    acc, rec, locs = [], [], []
    for ev in sorted(events, key=lambda e: -e[0]):
        s, p = float(ev[0]), np.asarray(ev[1:3], float)
        if not np.isfinite(s):
            continue
        if any(np.linalg.norm(p - a) <= min_sep_km for a in acc):
            continue
        acc.append(p)
        hit = -1
        if truth:
            dd = [float(np.linalg.norm(p - np.asarray(t_, float)[:2])) for t_ in truth]
            jj = int(np.argmin(dd))
            if dd[jj] <= 5.0:
                hit = jj
        rec.append((s, hit))
        locs.append([float(x) for x in ev[1:5]])
        if len(acc) >= max_events:
            break
    return rec, locs


def methods(nc, ops, lut, qmc, D, MFn, sd, tr, tmpdir):
    rec, loc, sec = {}, {}, {}
    t = time.time(); r, l = xc.bomp_pk(ops.PD, MFn, tr); sec["BOMP-PK"] = time.time() - t
    rec["BOMP-PK"], loc["BOMP-PK"] = r, l
    t = time.time(); r, l = xc.bomp_pk_mask(ops.PD, MFn, tr, mw=6, amp_gate=3.0); sec["BOMP-PK-MASK"] = time.time() - t
    rec["BOMP-PK-MASK"], loc["BOMP-PK-MASK"] = r, l
    w = {"D": D, "MFn": MFn, "sd": sd}
    t = time.time(); r, l, _ = xc.tssfw(ops, w, tr, f=xc.F_FROZEN); sec["TS-SFW"] = time.time() - t
    rec["TS-SFW"], loc["TS-SFW"] = r, l
    t = time.time()
    r, l = xc.E.mtfw_sc_locs(nc.STA, MFn, ops.grid, xc.DEFAULT_TAUS, ops.Gs, ops.Qs, tr, nc.corr)
    sec["MT-FW"] = time.time() - t
    rec["MT-FW"], loc["MT-FW"] = [(float(s), int(h)) for s, h in r], l
    t = time.time(); r, l = xc.mtfw_mask(ops, MFn, tr); sec["MT-FW-MASK"] = time.time() - t
    rec["MT-FW-MASK"], loc["MT-FW-MASK"] = r, l
    t = time.time()
    with T.silent():
        evs = T.run_qm_variants(nc.STA, D, tmpdir, {"QM": (qmc["QM-3C"], nc.corr)}, lut=lut)
    r, l = qm_rec_locs(evs["QM"], tr)
    sec["QM-3C+SC"] = time.time() - t
    rec["QM-3C+SC"], loc["QM-3C+SC"] = r, l
    return rec, loc, sec


def run(regime, nets, ne, nemp, out):
    qmc = T.load_qm_cfgs()
    tmpdir = tempfile.mkdtemp(prefix="c4qm_")
    res = []
    for cell in CELLS:
        for j in nets:
            nc = xc.NetCfg(cell, j, regime)
            ops = xc.Ops(nc.STA, nc.corr, xc.make_grid())
            with T.silent():
                lut = QM.build_lut(nc.STA, xc.VP, xc.VS, {**QM.QM_CFG, **qmc["QM-3C"]})
            wins = []
            st = xc.window_stream(nc, 3004)
            for wi in range(ne):
                w = next(st)
                wins.append(("eval", wi, w["D"], w["MFn"], w["sd"], w["truth"]))
            for ci_, (D, src_c, tau_c) in enumerate(nc.cat):
                MF, sd = xc.matched(D)
                wins.append(("cat", ci_, D, MF / sd, sd, [np.asarray(src_c, float)]))
            st = xc.window_stream(nc, 3005, draw=no_sources)
            for wi in range(nemp):
                w = next(st)
                wins.append(("empty", wi, w["D"], w["MFn"], w["sd"], []))
            for typ, wi, D, MFn, sd, tr in wins:
                rec, loc, sec = methods(nc, ops, lut, qmc, D, MFn, sd, tr, tmpdir)
                res.append({"cell": cell, "net": j, "type": typ, "w": wi, "k": len(tr),
                            "truth": [np.asarray(t_).tolist() for t_ in tr], "rec": rec, "loc": loc, "sec": sec})
            print(cell, j, "done", time.strftime("%H:%M:%S"), flush=True)
            with open(out, "w") as fh:
                json.dump(res, fh)


if __name__ == "__main__":
    run(sys.argv[1], [int(x) for x in sys.argv[2].split(",")], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])
