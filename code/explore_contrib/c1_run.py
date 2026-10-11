"""C1: grid spacing (3/6/9 km) x depths (3 vs 2), same +-0.25 s tolerance (plus a wide-tolerance variant)."""
import json, sys, time
import numpy as np
import xc

CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
GRIDS = {"g3d3": (3.0, (5.0, 9.0, 13.0)), "g3d2": (3.0, (7.0, 11.0)), "g6d3": (6.0, (5.0, 9.0, 13.0)),
         "g6d2": (6.0, (7.0, 11.0)), "g9d3": (9.0, (5.0, 9.0, 13.0)), "g9d2": (9.0, (7.0, 11.0))}
SHW = {3.0: 4, 6.0: 8, 9.0: 12}
TAG = 3001


def run(regime, nets, nwin, out):
    res = []
    for cell in CELLS:
        for j in nets:
            nc = xc.NetCfg(cell, j, regime)
            st = xc.window_stream(nc, TAG)
            W = [next(st) for _ in range(nwin)]
            opsd = {}
            for g, (step, dep) in GRIDS.items():
                opsd[g] = xc.Ops(nc.STA, nc.corr, xc.make_grid(step, dep))
            base = opsd["g3d3"]
            for wi, w in enumerate(W):
                tr = w["truth"]
                osc = xc.oracle_scores(base, w)
                row = {"cell": cell, "net": j, "w": wi, "k": len(tr), "truth": [t.tolist() for t in tr],
                       "oracle": osc, "rec": {}, "loc": {}, "sec": {}}
                for g, (step, dep) in GRIDS.items():
                    ops = opsd[g]
                    t = time.time(); r, l = xc.bomp_pk(ops.PD, w["MFn"], tr); row["sec"][g + ":BOMP-PK"] = time.time() - t
                    row["rec"][g + ":BOMP-PK"], row["loc"][g + ":BOMP-PK"] = r, l
                    if g == "g3d3":
                        A, Bc = xc.ceilings(r, osc, tr, l)
                        row["rec"]["CEIL-A"], row["rec"]["CEIL-B"] = A, Bc
                    t = time.time(); r, l, _ = xc.tssfw(ops, w, tr, f=xc.F_FROZEN); row["sec"][g + ":TS-SFW"] = time.time() - t
                    row["rec"][g + ":TS-SFW"], row["loc"][g + ":TS-SFW"] = r, l
                    t = time.time(); r = xc.mtfw(ops, w["MFn"], tr); row["sec"][g + ":MT-FW"] = time.time() - t
                    row["rec"][g + ":MT-FW"] = r
                    sh = SHW[step]
                    if sh != 4:
                        PW = xc.PKDictW(ops.BD, sh)
                        t = time.time(); r, l = xc.bomp_pk(PW, w["MFn"], tr); row["sec"][g + ":BOMP-PK-W"] = time.time() - t
                        row["rec"][g + ":BOMP-PK-W"], row["loc"][g + ":BOMP-PK-W"] = r, l
                        nw = xc.wide_net(ops, sh)
                        t = time.time(); r, l, _ = xc.tssfw_wide(ops, nw, w, tr, sh, f=xc.F_FROZEN)
                        row["sec"][g + ":TS-SFW-W"] = time.time() - t
                        row["rec"][g + ":TS-SFW-W"], row["loc"][g + ":TS-SFW-W"] = r, l
                res.append(row)
            print(cell, j, "done", time.strftime("%H:%M:%S"), flush=True)
            with open(out, "w") as fh:
                json.dump(res, fh)


if __name__ == "__main__":
    regime = sys.argv[1]
    nets = [int(x) for x in sys.argv[2].split(",")]
    nwin = int(sys.argv[3])
    run(regime, nets, nwin, sys.argv[4])
