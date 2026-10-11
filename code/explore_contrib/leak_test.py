"""exploratory: leakage suppression variants on the C1 baseline windows (g3d3, same stream TAG 3001)."""
import json, sys, time
import numpy as np
import xc

CELLS = ("ladder_4", "ladder_6", "unilateral_c40")


def run(regime, nets, nwin, out, tag=3001, draw=None, kmax=5):
    res = []
    for cell in CELLS:
        for j in nets:
            nc = xc.NetCfg(cell, j, regime)
            st = xc.window_stream(nc, tag) if draw is None else xc.window_stream(nc, tag, draw=draw)
            ops = xc.Ops(nc.STA, nc.corr, xc.make_grid())
            for wi in range(nwin):
                w = next(st)
                tr = w["truth"]
                row = {"cell": cell, "net": j, "w": wi, "k": len(tr), "rec": {}, "loc": {}, "sec": {}}
                t = time.time(); r, l = xc.bomp_pk(ops.PD, w["MFn"], tr, kmax=kmax); row["sec"]["BOMP-PK"] = time.time() - t
                row["rec"]["BOMP-PK"], row["loc"]["BOMP-PK"] = r, l
                for mw in (6, 10):
                    for ag in (3.0, 6.0):
                        nm = f"MASK{mw}g{int(ag)}"
                        t = time.time(); r, l = xc.bomp_pk_mask(ops.PD, w["MFn"], tr, kmax=kmax, mw=mw, amp_gate=ag)
                        row["sec"][nm] = time.time() - t
                        row["rec"][nm], row["loc"][nm] = r, l
                res.append(row)
            print(cell, j, "done", time.strftime("%H:%M:%S"), flush=True)
            with open(out, "w") as fh:
                json.dump(res, fh)


if __name__ == "__main__":
    run(sys.argv[1], [int(x) for x in sys.argv[2].split(",")], int(sys.argv[3]), sys.argv[4])
