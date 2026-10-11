"""C8 run: per-station/phase amplitude (site) terms in the moment-tensor subspace certificate.  SEED 999 ONLY.

Networks task2a_confirm.network(999, cell, j): main designs j = 60-79, sensitivity designs j = 80-89 (never used before:
2A/2E 0-4, C1-C3 10-15, C6/C6b 20-29, leak2 30-35, C7 40-59).  Window stream seeds [999, TAG, ci, j, ri, 41|233|61]:
  TAG 3801 DENSE-IID    world iid sigma 0.5 (= the default generator), c2_run.draw_dense, kmax 7
  TAG 3802 STD-IID      world iid sigma 0.5, draw_sources (0-2 events), kmax 5
  TAG 3803 DENSE-PERS   world pers sigma 0.5 (site factors fixed per network, rng [999, 3809, ci, j])
  TAG 3804 STD-PERS     world pers sigma 0.5
  TAG 3805 DENSE-IID-S25 / 3806 DENSE-IID-S80 (iid sigma 0.25 / 0.8), 3808 DENSE-PERS-S80 (pers sigma 0.8),
  TAG 3812 DENSE-PERS-PSSAME (pers sigma 0.5, S factor = P factor)   (sensitivity: networks 80-89, 4 windows each,
  without SITE-CAT-GAIN and CEIL-GD).  TAGs 3807 and 3811 reserved, not run.
  HUBER-N was first omitted from STD-IID, STD-PERS and the sensitivity designs (compute); it was put back after the
  DENSE-IID results (post hoc, before those designs ran).  DENSE-IID-HSWEEP (post hoc): the same DENSE-IID windows
  (TAG 3801), MT-FW and HUBER-N (identity check) plus the Huber constant at the chi2_3 0.90 / 0.99 quantiles.
  Window counts per network: DENSE-IID 5 (j 60-79), DENSE-PERS 4 (j 60-79), STD-IID 4 (j 60-79), STD-PERS 4 (j 60-69).
  TAG 3810 validation windows (c8_validate.py), TAG 3809 persistent site draws.
    python c8_run.py DESIGN REGIME j0,j1,... NWIN out.json
"""
import json
import os
import sys
import time

import numpy as np

import c2_run as C2
import c6_ceil as CZ
import c8_lib as C
import xc

HS = ("HUBER-N90", "HUBER-N99", "HUBER-C1.5", "HUBER-C5", "HUBER-C8")
DESIGNS = {                         # tag, kind, world, CEIL-GD, methods dropped, only these (None = all)
    "DENSE-IID": (3801, "dense", C.World("iid", 0.5), True, ("SITE-CAT-GAIN",) + HS, None),
    "STD-IID": (3802, "std", C.World("iid", 0.5), True, ("SITE-CAT-GAIN",) + HS, None),
    "DENSE-PERS": (3803, "dense", C.World("pers", 0.5), True, HS, None),
    "STD-PERS": (3804, "std", C.World("pers", 0.5), False, HS, None),
    "DENSE-IID-S25": (3805, "dense", C.World("iid", 0.25), False, ("SITE-CAT-GAIN",) + HS, None),
    "DENSE-IID-S80": (3806, "dense", C.World("iid", 0.8), False, ("SITE-CAT-GAIN",) + HS, None),
    "DENSE-PERS-S80": (3808, "dense", C.World("pers", 0.8), False, ("SITE-CAT-GAIN",) + HS, None),
    "DENSE-PERS-PSSAME": (3812, "dense", C.World("pers", 0.5, True), False, ("SITE-CAT-GAIN",) + HS, None),
    "DENSE-IID-HSWEEP": (3801, "dense", C.World("iid", 0.5), False, (), ("MT-FW", "HUBER-N", "HUBER-N90", "HUBER-N99")),
    # post hoc, wider constants, networks 60-69 only (same DENSE-IID windows)
    "DENSE-IID-HWIDE": (3801, "dense", C.World("iid", 0.5), False, (),
                        ("MT-FW", "HUBER-N", "HUBER-C1.5", "HUBER-C5", "HUBER-C8")),
}
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")


def methods(dops, net, w, kmax, drop=(), only=None):
    MFn, tr = w["MFn"], w["truth"]
    ecat, ecst = np.exp(net.u_cat), np.exp(net.u_cat_st)
    osite = np.asarray(w["site"]).T
    out = {
        "MT-FW": dict(stat="ls", capv=6.0),
        "MT-FW-NOCAP": dict(stat="ls", capv=None),
        "MT-FW-C3": dict(stat="ls", capv=3.0),
        "MT-FW-C8": dict(stat="ls", capv=8.0),
        "MT-FW-C12": dict(stat="ls", capv=12.0),
        "ORACLE": dict(stat="ls", capv=None, wsite=osite),
        "ORACLE-CAP": dict(stat="ls", capv=6.0, wsite=osite),
        "SITE-CAT": dict(stat="ls", capv=None, wsite=ecat),
        "SITE-CAT-CAP": dict(stat="ls", capv=6.0, wsite=ecat),
        "SITE-CAT-ST": dict(stat="ls", capv=None, wsite=ecst),
        "HUBER-N": dict(stat="huber_n"),
        "HUBER-S": dict(stat="huber_s"),
        "HUBER-N90": dict(stat="huber_n", chub=2.5003),          # sqrt(chi2_3 0.90)
        "HUBER-N99": dict(stat="huber_n", chub=3.3682),          # sqrt(chi2_3 0.99)
        "HUBER-C1.5": dict(stat="huber_n", chub=1.5),
        "HUBER-C5": dict(stat="huber_n", chub=5.0),
        "HUBER-C8": dict(stat="huber_n", chub=8.0),
        "GAIN": dict(stat="gain", sp=0.5),
        "GAIN-FREE": dict(stat="gain", sp=np.inf),
        "SITE-CAT-GAIN": dict(stat="gain", sp=0.5, wsite=ecat),
    }
    res = {}
    out = {k: v for k, v in out.items() if k not in drop and (only is None or k in only)}
    for name, kw in out.items():
        t = time.time()
        r, l = C.mtfw_general(dops, MFn, tr, kmax, **kw)
        res[name] = (r, l, time.time() - t)
    if only is None or "NORM-TR" in only:
        t = time.time()
        r, l = C.mtfw_general(dops, C.norm_traces(MFn), tr, kmax, stat="ls", capv=None)
        res["NORM-TR"] = (r, l, time.time() - t)
    return res


def run(design, regime, nets, nwin, path):
    tag, kind, world, ceil, drop, only = DESIGNS[design]
    draw = C2.draw_dense if kind == "dense" else xc.draw_sources
    kmax = 7 if kind == "dense" else 5
    res = json.load(open(path)) if os.path.exists(path) else {"rows": [], "nets": []}
    done = {(n["cell"], n["net"]) for n in res["nets"]}
    grid = xc.make_grid()
    garr = np.asarray(grid)
    for cell in CELLS:
        for j in nets:
            if (cell, j) in done:
                continue
            t0 = time.time()
            net = C.WNet(cell, j, regime, world)
            xops = xc.Ops(net.STA, net.corr, grid)
            dops = C.DOps(xops)
            st = net.stream(tag, draw)
            ws = [next(st) for _ in range(nwin)]
            rows = []
            for wi, w in enumerate(ws):
                tr = w["truth"]
                row = {"cell": cell, "net": j, "w": wi, "k": len(tr), "truth": [t.tolist() for t in tr],
                       "taus": w["taus"], "meta": w["meta"], "site": np.asarray(w["site"]).tolist(),
                       "rec": {}, "loc": {}, "sec": {}}
                if only is None or "BOMP-PK" in only:
                    t = time.time()
                    r, l = xc.bomp_pk(xops.PD, w["MFn"], tr, kmax)
                    row["sec"]["BOMP-PK"] = time.time() - t
                    row["rec"]["BOMP-PK"], row["loc"]["BOMP-PK"] = r, l
                for name, (r, l, s) in methods(dops, net, w, kmax, drop, only).items():
                    row["rec"][name], row["loc"][name], row["sec"][name] = r, l, s
                if ceil:
                    t = time.time()
                    R_or = CZ.oracle_residual_refined(xops.net, w, tr, w["taus"])
                    g1 = CZ.ceil_scores(xops.PD, garr, w["MFn"], R_or, tr)
                    row["rec"]["CEIL-GD"] = CZ.record(*g1)
                    row["sec"]["CEIL-GD"] = time.time() - t
                rows.append(row)
            res["rows"].extend(rows)
            res["nets"].append({"cell": cell, "net": j, "regime": regime, "design": design, "world": world.name(),
                                "u_cat": net.u_cat.tolist(), "u_cat_st": net.u_cat_st.tolist(),
                                "n_obs": net.u_info["n_obs"],
                                "site_net": None if net.site_net is None else net.site_net.tolist(),
                                "t_total": time.time() - t0})
            with open(path + ".tmp", "w") as fh:
                json.dump(res, fh)
            os.replace(path + ".tmp", path)
            print(design, regime, cell, j, "done", time.strftime("%H:%M:%S"), f"{time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], [int(x) for x in sys.argv[3].split(",")], int(sys.argv[4]), sys.argv[5])
