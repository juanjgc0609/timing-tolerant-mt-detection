"""C6b: attribution of the BOMP-PK headroom to the generator/detector mismatches (exploratory, seed 999).
Networks 20-29 (as C6), windows TAG 3008 (new), 5 per network, both regimes.  The generator's new randomness (pulse,
coda) comes from a child of rng_h and is drawn identically in every variant, so the variants share sources,
mechanisms, statics, clutter and noise (paired).  Each variant removes ONE mismatch from the default generator:
  DEF     default gen_realistic (gradient medium, variable pulse, coda, site factors lognormal 0.5, 10 deg polarisation)
  MED     homogeneous medium (VS0 3.0, g 0): the detectors' travel times / ray directions are exact up to statics
  PULSE   fixed 2 Hz Ricker, no phase rotation: the matched filter is exact
  SITE    no site amplification factors, no polarisation deviation: the radiation pattern is exact
  CODA    no coda
Catalogue, station corrections (+SC) and the ORA model (the variant's true medium) are rebuilt per variant.
Methods: BOMP-PK:HOM, MT-FW:HOM, BOMP-PK:ORA, MT-FW:ORA, CEIL-GD (verif_2e definition, c6_ceil).
    python c6b_run.py thin 20,...,29 5 DEF,MED,PULSE,SITE,CODA out/c6b_thin.json"""
import json
import os
import sys
import time

import numpy as np

import c6_ceil as CZ
import c6_lib as L
import xc

VARIANTS = {"DEF": {}, "MED": {"VS0": 3.0, "g": 0.0}, "PULSE": {"fc_range": (2.0, 2.0), "phase_max_deg": 0.0},
            "SITE": {"site_sd": 0.0, "dir_deg": 0.0}, "CODA": {"coda_ratio": 0.0}}
TAG = 3008


def run(regime, nets, nwin, variants, path):
    res = json.load(open(path)) if os.path.exists(path) else {"rows": []}
    done = {(r["var"], r["cell"], r["net"]) for r in res["rows"]}
    grid = xc.make_grid()
    garr = np.asarray(grid)
    for var in variants:
        cfg = VARIANTS[var] or None
        for cell in ("ladder_4", "ladder_6", "unilateral_c40"):
            for j in nets:
                if (var, cell, j) in done:
                    continue
                t0 = time.time()
                nc = xc.NetCfg(cell, j, regime, cfg=cfg)
                ops0 = xc.Ops(nc.STA, nc.corr, grid)
                gc = L.GR.make_cfg(cfg)
                ORA = L.Gradient(gc["VS0"], gc["g"], gc["vpvs"])
                opsH = L.ModelOps(nc.STA, grid, L.Homog(), nc.corr)
                opsO = L.ModelOps(nc.STA, grid, ORA, L.corr_model(nc.STA, nc.cat, ORA))
                st = xc.window_stream(nc, TAG)
                for wi in range(nwin):
                    w = next(st)
                    tr = w["truth"]
                    row = {"var": var, "cell": cell, "net": j, "w": wi, "k": len(tr), "truth": [t.tolist() for t in tr],
                           "rec": {}, "loc": {}}
                    for name, ops in (("HOM", opsH), ("ORA", opsO)):
                        r, l = L.bomp_pk_m(ops, w["MFn"], tr)
                        row["rec"]["BOMP-PK:" + name], row["loc"]["BOMP-PK:" + name] = r, l
                        r, l = L.mtfw_dense(ops, w["MFn"], tr)
                        row["rec"]["MT-FW:" + name], row["loc"]["MT-FW:" + name] = r, l
                    R_or = CZ.oracle_residual_refined(ops0.net, w, tr, w["taus"])
                    row["rec"]["CEIL-GD"] = CZ.record(*CZ.ceil_scores(opsH.PD, garr, w["MFn"], R_or, tr))
                    res["rows"].append(row)
                with open(path + ".tmp", "w") as fh:
                    json.dump(res, fh)
                os.replace(path + ".tmp", path)
                print(regime, var, cell, j, "done", time.strftime("%H:%M:%S"), f"{time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    run(sys.argv[1], [int(x) for x in sys.argv[2].split(",")], int(sys.argv[3]), sys.argv[4].split(","), sys.argv[5])
