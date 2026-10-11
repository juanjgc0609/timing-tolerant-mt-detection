"""C10 post hoc: MT-FW family with the SC-R corrections (c10_screloc.sc_reloc; catalogue data only) instead of +SC.
MT@R (hard rule), MT@R+T2w (corr prior, +-0.5 s), MT@R+T3 (refinement); timing level by the catalogue rule recomputed
with the SC-R corrections.  Networks 160-179, the main windows (regenerated, MT-FW atoms asserted identical).
    python c10_screloc_run.py DESIGN REGIME j0,... NWIN out.json"""
import json
import os
import sys
import time

import numpy as np

import c2_run as C2
import c9_run as R9
import c10_lib as C
import c10_run as R
import c10_screloc as SR
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
            nc = xc.NetCfg(cell, j, regime)
            sig0 = C.sigma_catalogue(nc.STA, nc.cat, nc.corr)[0][0]
            corrR = SR.sc_reloc(nc.STA, nc.cat, nc.corr, sig0)
            S = R.Setup(cell, j, regime, oracle=False, corr=corrR)
            E, pr, EF, prF = S.eng["cat"]
            st = R9.stream(S.nc, tag, draw, None)
            for wi in range(nwin):
                w = next(st)
                tr = w["truth"]
                m = main[(cell, j, wi)]
                assert np.allclose(m["truth"], [t.tolist() for t in tr])
                atoms = {k: m["atoms"][k] for k in ("MT-FW", "BOMP-PK", "MT+T2w", "MT+T3", "MT-ORA")}
                atoms["MT@R"] = R.pack(C.mt_loop(E, w["MFn"], kmax))
                atoms["MT@R+T2w"] = R.pack(C.mt_loop(E, w["MFn"], kmax, mode="corr", pr=pr, W=8))
                atoms["MT@R+T3"] = R.pack(C.mt_loop(E, w["MFn"], kmax, mode="hard", pr=pr, W=4, EF=EF, prF=prF))
                res["rows"].append({"cell": cell, "net": j, "w": wi, "k": len(tr), "truth": [t_.tolist() for t_ in tr],
                                    "taus": w["taus"], "atoms": atoms, "ceil": m["ceil"], "sec": {}})
            res["nets"].append({"cell": cell, "net": j, "regime": regime, "design": design, "t_total": time.time() - t0,
                                "t_setup": S.t_setup, "sig": S.sig.tolist(), "sig_ora": [0, 0],
                                "err0": SR.err(nc.corr, nc.statics), "errR": SR.err(corrR, nc.statics),
                                "cell_sd": [0, 0], "cell_sd_fine": [0, 0]})
            with open(path + ".tmp", "w") as fh:
                json.dump(res, fh)
            os.replace(path + ".tmp", path)
            print("SCR", design, regime, cell, j, "done", time.strftime("%H:%M:%S"), f"{time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], [int(x) for x in sys.argv[3].split(",")], int(sys.argv[4]), sys.argv[5])
