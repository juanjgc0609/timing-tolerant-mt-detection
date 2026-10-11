"""C7 run: source-pulse-invariant moment-tensor subspace detectors.  SEED 999 ONLY.

Networks task2a_confirm.network(999, cell, j), j = 40..59 (never used before: C1-C3 used 10-15, C6/C6b 20-29, leak2
30-35, 2A/2E calibration 0-4).  Catalogue / +SC corrections: xc.NetCfg (20 events, statics sd 0.3 s), rebuilt with
the out-of-family pulse in the OOF designs (the catalogue lives in the same world as the windows).
Designs (window stream seeds [999, TAG, ci, j, ri, 41|233|61], TAGs new):
  STD      TAG 3701, draw_sources (2A design, 0-2 events), kmax 5, default generator
  DENSE    TAG 3702, c2_run.draw_dense (3-5 events within 9 km, gaps >= 0.5 s), kmax 7, default generator
  GABOR90  TAG 3703, draw_sources, kmax 5, generator with the Gabor pulse (width 0.6/f0, phase U[-90, 90])
  BERLAGE90 TAG 3704, draw_sources, kmax 5, generator with the causal Berlage pulse (phase U[-90, 90])
Methods (grid 3 km x 3 depths, origins 5.0-8.0 s, +-4-sample peak search, 5 km exclusion, hit within 5 km):
  BOMP-PK            registered BOMP-PK+SC (task2d_methods.bomp_pk_record)
  BOMP-PK-NX         c7 code with q = 1 (2 Hz MF), no P/S cross term in the statistic (control)
  BOMP-PK-CAP        BOMP-PK-NX with the per-station/phase cap in the selection statistic
  BOMP-PK-PH         PRIMARY 1: q = 2 {Ricker 2 Hz, Hilbert}, rank-1 (separable) statistic and bilinear refit
  BOMP-PK-PHF        PRIMARY 2: q = 4 SVD basis of the design family (rule fixed before any detection run: smallest q
                     whose worst-case captured energy over the family is >= 0.95), rank-1
  BOMP-PK-PH-FULL, BOMP-PK-PHF-FULL   relaxed full 6q blocks (Frobenius statistic, unconstrained refit)
  BOMP-PK-PHF3       sensitivity: q = 3 SVD basis, rank-1
  BOMP-PK-PH-CAP, BOMP-PK-PHF-CAP     primaries with the cap (3q-vector norm <= 6) in the selection statistic
  MT-FW              dense MT-FW+SC (c6_lib.mtfw_dense, capped)
  MT-FW-PH, MT-FW-PHF   c7_lib.mtfw_q (capped, rank-1, sequential rank-1 deflation), q = 2 / 4
  CEIL-GD            verif_2e ceiling (c6_ceil: truth-start refined oracle deflation, HOM BOMP-PK statistic)
    python c7_run.py STD thin 40,...,59 5 out/c7_STD_thin.json
"""
import json
import os
import sys
import time

import numpy as np

import c2_run as C2
import c6_ceil as CZ
import c6_lib as L
import c7_lib as C
import xc

GR = C.GR
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
DESIGNS = {"STD": (3701, None, 5, None), "DENSE": (3702, C2.draw_dense, 7, None),
           "GABOR90": (3703, None, 5, "GABOR90"), "BERLAGE90": (3704, None, 5, "BERLAGE90")}
BASES = {"Q1": C.Basis(C.basis_q1(), "Q1"), "PH": C.Basis(C.basis_ph(), "PH"),
         "PHF": C.Basis(C.basis_phf(4)[0], "PHF4"), "PHF3": C.Basis(C.basis_phf(3)[0], "PHF3")}
BOMPS = {"BOMP-PK-NX": ("Q1", "R1", False), "BOMP-PK-CAP": ("Q1", "R1", True),
         "BOMP-PK-PH": ("PH", "R1", False), "BOMP-PK-PHF": ("PHF", "R1", False),
         "BOMP-PK-PH-FULL": ("PH", "FULL", False), "BOMP-PK-PHF-FULL": ("PHF", "FULL", False),
         "BOMP-PK-PHF3": ("PHF3", "R1", False),
         "BOMP-PK-PH-CAP": ("PH", "R1", True), "BOMP-PK-PHF-CAP": ("PHF", "R1", True)}
FWS = {"MT-FW-PH": "PH", "MT-FW-PHF": "PHF"}


class _Nullctx:
    def __enter__(self):
        return self

    def __exit__(self, *a):
        pass


def stream(nc, tag, draw):
    b = [xc.SEED, tag, nc.ci, nc.j, nc.ri]
    rng, rng_m, rng_h = (np.random.default_rng(b + [41]), np.random.default_rng(b + [233]),
                         np.random.default_rng(b + [61]))
    while True:
        srcs, taus = draw(rng)
        meta = []
        D, picks, _tr, _amp = GR.generate_real(nc.STA, srcs, taus, nc.nclut, rng, rng_m, rng_h, nc.clut_amp,
                                               cfg=nc.cfg, meta=meta, statics=nc.statics)
        MF, sd = xc.matched(D)
        yield {"D": D, "MFn": MF / sd, "sd": sd, "truth": [np.asarray(s, float) for s in srcs],
               "taus": [float(t) for t in taus], "picks": picks, "meta": meta}


def run(design, regime, nets, nwin, path):
    tag, draw, kmax, oof = DESIGNS[design]
    draw = draw or xc.draw_sources
    cfg = C.OOF[oof][1] if oof else None
    ctx = C.OOFPulse(C.OOF[oof][0]) if oof else _Nullctx()
    res = json.load(open(path)) if os.path.exists(path) else {"rows": [], "nets": []}
    done = {(n["cell"], n["net"]) for n in res["nets"]}
    grid = xc.make_grid()
    garr = np.asarray(grid)
    for cell in CELLS:
        for j in nets:
            if (cell, j) in done:
                continue
            t0 = time.time()
            with ctx:
                nc = xc.NetCfg(cell, j, regime, cfg=cfg)
                ops = xc.Ops(nc.STA, nc.corr, grid)
                mops = L.ModelOps(nc.STA, grid, L.Homog(), nc.corr)
                QDs = {k: C.QDict(ops.BD, b, mops.Qst) for k, b in BASES.items()}
                st = stream(nc, tag, draw)
                ws = [next(st) for _ in range(nwin)]
            rows = []
            for wi, w in enumerate(ws):
                tr = w["truth"]
                row = {"cell": cell, "net": j, "w": wi, "k": len(tr), "truth": [t.tolist() for t in tr],
                       "taus": w["taus"], "meta": w["meta"], "rec": {}, "loc": {}, "sec": {}}
                t = time.time()
                r, l = xc.bomp_pk(ops.PD, w["MFn"], tr, kmax)
                row["sec"]["BOMP-PK"] = time.time() - t
                row["rec"]["BOMP-PK"], row["loc"]["BOMP-PK"] = r, l
                t = time.time()
                F = {k: b.filt(w["D"], w["sd"], w["MFn"]) for k, b in BASES.items()}
                row["sec"]["filters"] = time.time() - t
                for name, (bk, mode, cap) in BOMPS.items():
                    t = time.time()
                    r, l = C.bomp_q(QDs[bk], F[bk], tr, kmax, mode, cap)
                    row["sec"][name] = time.time() - t
                    row["rec"][name], row["loc"][name] = r, l
                t = time.time()
                r, l = L.mtfw_dense(mops, w["MFn"], tr, kmax)
                row["sec"]["MT-FW"] = time.time() - t
                row["rec"]["MT-FW"], row["loc"]["MT-FW"] = r, l
                for name, bk in FWS.items():
                    t = time.time()
                    r, l = C.mtfw_q(QDs[bk], mops.IDX, F[bk], tr, kmax)
                    row["sec"][name] = time.time() - t
                    row["rec"][name], row["loc"][name] = r, l
                t = time.time()
                R_or = CZ.oracle_residual_refined(ops.net, w, tr, w["taus"])
                g1 = CZ.ceil_scores(ops.PD, garr, w["MFn"], R_or, tr)
                row["rec"]["CEIL-GD"] = CZ.record(*g1)
                row["sec"]["CEIL-GD"] = time.time() - t
                rows.append(row)
            res["rows"].extend(rows)
            res["nets"].append({"cell": cell, "net": j, "regime": regime, "design": design,
                                "t_total": time.time() - t0})
            with open(path + ".tmp", "w") as fh:
                json.dump(res, fh)
            os.replace(path + ".tmp", path)
            print(design, regime, cell, j, "done", time.strftime("%H:%M:%S"), f"{time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], [int(x) for x in sys.argv[3].split(",")], int(sys.argv[4]), sys.argv[5])
