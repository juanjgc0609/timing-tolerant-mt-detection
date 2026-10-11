"""C10 extra (POST HOC, after the DENSE analysis showed a gain for the oracle only): which half of the oracle matters.
MT-FW family with the TRUE statics as corrections but the timing level of the CATALOGUE rule (the level the detector
would really use): MT-ORA@cat+T2w (corr prior, +-0.5 s) and MT-ORA@cat+T3 (refinement).  The converse (catalogue
corrections, small fixed level) is in c10_sweep.py.  Networks 160-179, the main windows (regenerated; MT-FW atoms
asserted identical to the main run).
    python c10_extra.py DESIGN REGIME j0,... NWIN out.json"""
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
            S = R.Setup(cell, j, regime)
            E = S.eng["cat"][0]
            EO, _p, EFO, _pf = S.eng["ora"]
            pr, prF = EO.prior(S.sig), EFO.prior(S.sig)
            st = R9.stream(S.nc, tag, draw, None)
            for wi in range(nwin):
                w = next(st)
                tr = w["truth"]
                m = main[(cell, j, wi)]
                a0 = R.pack(C.mt_loop(E, w["MFn"], kmax))
                assert np.allclose(a0, m["atoms"]["MT-FW"], equal_nan=True)
                atoms = {k: m["atoms"][k] for k in ("MT-FW", "BOMP-PK", "MT-ORA", "MT-ORA+T2w", "MT-ORA+T3")}
                atoms["MT-ORA@cat+T2w"] = R.pack(C.mt_loop(EO, w["MFn"], kmax, mode="corr", pr=pr, W=8))
                atoms["MT-ORA@cat+T3"] = R.pack(C.mt_loop(EO, w["MFn"], kmax, mode="hard", pr=pr, W=4, EF=EFO, prF=prF))
                res["rows"].append({"cell": cell, "net": j, "w": wi, "k": len(tr), "truth": [t_.tolist() for t_ in tr],
                                    "taus": w["taus"], "atoms": atoms, "ceil": m["ceil"], "sec": {}})
            res["nets"].append({"cell": cell, "net": j, "regime": regime, "design": design, "t_total": time.time() - t0,
                                "t_setup": S.t_setup, "sig": S.sig.tolist(), "sig_ora": S.sig_ora.tolist(),
                                "cell_sd": [0, 0], "cell_sd_fine": [0, 0]})
            with open(path + ".tmp", "w") as fh:
                json.dump(res, fh)
            os.replace(path + ".tmp", path)
            print("EXTRA", design, regime, cell, j, "done", time.strftime("%H:%M:%S"), f"{time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], [int(x) for x in sys.argv[3].split(",")], int(sys.argv[4]), sys.argv[5])
