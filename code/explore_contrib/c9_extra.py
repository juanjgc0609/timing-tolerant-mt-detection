"""C9 extra (post hoc, after the first look at DENSE thin): polarization term inside the MT-FW selection statistic.
Regenerates the windows of a design from the same streams and runs MT-FW (identity check against the main run, enforced
in c9_analyze.load) plus MT-FW+SELPOL and MT-FW+PENPOL.
    python c9_extra.py DESIGN REGIME j0,j1,... NWIN out.json"""
import json
import os
import sys
import time

import numpy as np

import c2_run as C2
import c8_lib as C8
import c9_lib as C
import c9_run as R
import xc


def run(design, regime, nets, nwin, path):
    tag, kind, cfg, _ceil, _full = R.DESIGNS[design]
    draw = C2.draw_dense if kind == "dense" else xc.draw_sources
    kmax = 7 if kind == "dense" else 5
    res = json.load(open(path)) if os.path.exists(path) else {"rows": [], "nets": []}
    done = {(n["cell"], n["net"]) for n in res["nets"]}
    grid = xc.make_grid()
    garr = np.asarray(grid)
    for cell in R.CELLS:
        for j in nets:
            if (cell, j) in done:
                continue
            t0 = time.time()
            nc = xc.NetCfg(cell, j, regime)
            xops = xc.Ops(nc.STA, nc.corr, grid)
            dops = C8.DOps(xops)
            geo = C.Geo(garr, nc.STA)
            st = R.stream(nc, tag, draw, cfg)
            for wi in range(nwin):
                w = next(st)
                atoms = {"MT-FW": R.pack(C.mtfw_c9(dops, geo, w["MFn"], kmax)),
                         "MT-FW+SELPOL": R.pack(C.mtfw_c9(dops, geo, w["MFn"], kmax, sel="selpol")),
                         "MT-FW+PENPOL": R.pack(C.mtfw_c9(dops, geo, w["MFn"], kmax, sel="penpol"))}
                res["rows"].append({"cell": cell, "net": j, "w": wi, "truth": [t.tolist() for t in w["truth"]],
                                    "atoms": atoms})
            res["nets"].append({"cell": cell, "net": j})
            with open(path + ".tmp", "w") as fh:
                json.dump(res, fh)
            os.replace(path + ".tmp", path)
            print("EXTRA", design, regime, cell, j, "done", time.strftime("%H:%M:%S"), f"{time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], [int(x) for x in sys.argv[3].split(",")], int(sys.argv[4]), sys.argv[5])
