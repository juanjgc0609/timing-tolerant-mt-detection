"""C10 sensitivity sweep of the timing level (POST HOC in the sense of the task: not part of the rule; the design below
was fixed before the main results were analysed).  MT-FW family only, networks 160-169 (first half of the main sample),
the same windows (streams 4001 / 4002, regenerated; the MT-FW atoms are asserted identical to the main run).
Fixed sig (P = S) in SIGS for every network; variants T1c (indepc, +-0.25 s), T2w (corr, +-0.5 s), T3 (refine).
    python c10_sweep.py DESIGN REGIME j0,j1,... NWIN out.json"""
import json
import os
import sys
import time

import numpy as np

import c2_run as C2
import c9_run as R9
import c10_lib as C
import c10_run as R
import xc

SIGS = (0.05, 0.10, 0.20)
VARS = {"T1c": ("indepc", 4, False), "T2w": ("corr", 8, False), "T3": ("hard", 4, True)}


def run(design, regime, nets, nwin, path):
    tag, kind = R.DESIGNS[design]
    draw = C2.draw_dense if kind == "dense" else xc.draw_sources
    kmax = 7 if kind == "dense" else 5
    main = {(r["cell"], r["net"], r["w"]): r for r in json.load(open(f"out/c10_{design}_{regime}.json"))["rows"]}
    res = json.load(open(path)) if os.path.exists(path) else {"rows": [], "nets": []}
    done = {(n["cell"], n["net"]) for n in res["nets"]}
    for cell in R.CELLS:
        for j in nets:
            if (cell, j) in done:
                continue
            t0 = time.time()
            S = R.Setup(cell, j, regime, oracle=False)
            E, _pr, EF, _prF = S.eng["cat"]
            prs = {sg: (E.prior((sg, sg)), EF.prior((sg, sg))) for sg in SIGS}
            st = R9.stream(S.nc, tag, draw, None)
            for wi in range(nwin):
                w = next(st)
                tr = w["truth"]
                atoms = {"MT-FW": R.pack(C.mt_loop(E, w["MFn"], kmax))}
                m = main[(cell, j, wi)]
                assert np.allclose(atoms["MT-FW"], m["atoms"]["MT-FW"], equal_nan=True) and np.allclose(m["truth"], [t.tolist() for t in tr])
                atoms["BOMP-PK"] = m["atoms"]["BOMP-PK"]
                for sg in SIGS:
                    pr, prF = prs[sg]
                    for nm, (mode, W, ref) in VARS.items():
                        atoms[f"MT+{nm}@{sg:g}"] = R.pack(C.mt_loop(E, w["MFn"], kmax, mode=mode, pr=pr, W=W,
                                                                     EF=EF if ref else None, prF=prF))
                res["rows"].append({"cell": cell, "net": j, "w": wi, "k": len(tr), "truth": [t_.tolist() for t_ in tr],
                                    "taus": w["taus"], "atoms": atoms, "ceil": m["ceil"], "sec": {}})
            res["nets"].append({"cell": cell, "net": j, "regime": regime, "design": design, "t_total": time.time() - t0,
                                "t_setup": S.t_setup, "sig": S.sig.tolist(), "sig_ora": [0, 0],
                                "cell_sd": [0, 0], "cell_sd_fine": [0, 0]})
            with open(path + ".tmp", "w") as fh:
                json.dump(res, fh)
            os.replace(path + ".tmp", path)
            print("SWEEP", design, regime, cell, j, "done", time.strftime("%H:%M:%S"), f"{time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], [int(x) for x in sys.argv[3].split(",")], int(sys.argv[4]), sys.argv[5])
