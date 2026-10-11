"""replication of arrival masking (MASK6g3 fixed after the leak_test pilot) on fresh seed-999 networks, TAG 3007."""
import json, sys, time
import numpy as np
import xc

CELLS = ("ladder_4", "ladder_6", "unilateral_c40")


def run(regime, nets, nwin, out):
    res = []
    for cell in CELLS:
        for j in nets:
            nc = xc.NetCfg(cell, j, regime)
            st = xc.window_stream(nc, 3007)
            ops = xc.Ops(nc.STA, nc.corr, xc.make_grid())
            for wi in range(nwin):
                w = next(st)
                tr = w["truth"]
                row = {"cell": cell, "net": j, "w": wi, "k": len(tr), "truth": [t.tolist() for t in tr],
                       "rec": {}, "loc": {}, "sec": {}}
                t = time.time(); r, l = xc.bomp_pk(ops.PD, w["MFn"], tr); row["sec"]["BOMP-PK"] = time.time() - t
                row["rec"]["BOMP-PK"], row["loc"]["BOMP-PK"] = r, l
                osc = xc.oracle_scores(ops, w)
                row["rec"]["CEIL-A"], row["rec"]["CEIL-B"] = xc.ceilings(r, osc, tr, l)
                t = time.time(); r, l = xc.bomp_pk_mask(ops.PD, w["MFn"], tr, mw=6, amp_gate=3.0)
                row["sec"]["BOMP-PK-MASK"] = time.time() - t
                row["rec"]["BOMP-PK-MASK"], row["loc"]["BOMP-PK-MASK"] = r, l
                t = time.time(); r, l = xc.E.mtfw_sc_locs(nc.STA, w["MFn"], ops.grid, xc.DEFAULT_TAUS, ops.Gs, ops.Qs, tr, nc.corr)
                row["sec"]["MT-FW"] = time.time() - t
                row["rec"]["MT-FW"], row["loc"]["MT-FW"] = [(float(s), int(h)) for s, h in r], l
                t = time.time(); r, l = xc.mtfw_mask(ops, w["MFn"], tr); row["sec"]["MT-FW-MASK"] = time.time() - t
                row["rec"]["MT-FW-MASK"], row["loc"]["MT-FW-MASK"] = r, l
                res.append(row)
            print(regime, cell, j, "done", time.strftime("%H:%M:%S"), flush=True)
            with open(out, "w") as fh:
                json.dump(res, fh)


if __name__ == "__main__":
    run(sys.argv[1], [int(x) for x in sys.argv[2].split(",")], int(sys.argv[3]), sys.argv[4])
