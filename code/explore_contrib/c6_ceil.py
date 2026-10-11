"""C6 ceilings, exactly the CEIL-GD of verif_2e/check4_ceiling_verif.py (oracle deflation by the joint LS fit of the
TRUTH-START refined continuous atoms, 500 L-BFGS-B iterations, shifts initialised at the energy peaks), on the C6
windows (seed 999, networks 20-29, TAG 3006, regenerated from the same streams).
  CEIL-GD    true events: max HOM BOMP-PK statistic e/2 within 5 km; false: greedy (5 km exclusion, <= 5) maxima of
             the HOM statistic on the oracle-deflated residual, > 5 km from every source
  CEIL-GD-O  the same with the ORA (gradient oracle model) statistic on the same residual
The first c6_run.py ceiling (deflation at the TRUE arrival times with the fixed 2 Hz atom) leaves the pulse-shape
mismatch in the residual and is NOT the verification's CEIL-GD; it is kept in the c6 json as CEIL-GD/-O(truetimes).
    python c6_ceil.py validate            (reproduces verif_2e check4 g_true / g_resid on 2A net 0, 4 windows)
    python c6_ceil.py thin 20,...,29 10 out/c6_ceil_thin.json"""
import json
import os
import sys
import time

import numpy as np

import c6_lib as L
import xc

sys.path.insert(0, os.path.join(xc.WORK, "verif_2e"))
E = xc.E


def near(garr, src):
    return np.where(np.linalg.norm(garr[:, :2] - np.asarray(src)[None, :2], axis=1) <= 5.0)[0]


def oracle_residual_refined(net, w, srcs, taus):
    import check4_ceiling_verif as C4
    win = E.Win(w["D"], w["MFn"], w["sd"])
    P_or = []
    for src, tau in zip(srcs, taus):
        _p0, _E0, p1, _E1 = C4.cont_refine(net, win, w["MFn"], src, float(tau))
        P_or.append(p1)
    if not P_or:
        return w["MFn"]
    Po = np.array(P_or)
    mo, _c = E.convex_step(net, win, Po, np.zeros((len(Po), 6)), 0.0)
    return w["MFn"] - E.predict_sampled(net, win, Po, mo)


def ceil_scores(PD, garr, MFn, R_or, srcs):
    e_full = PD.energy(MFn, True)[0] / 2.0
    g_true = [float(e_full[near(garr, s)].max()) for s in srcs]
    emr = (PD.energy(R_or, True)[0] / 2.0).max(axis=1)
    allowed = np.ones(len(garr), bool)
    for s in srcs:
        allowed[near(garr, s)] = False
    g_res = []
    for _ in range(5):
        cand = np.where(allowed, emr, -np.inf)
        n = int(np.argmax(cand))
        if not np.isfinite(cand[n]):
            break
        g_res.append(float(cand[n]))
        allowed[np.linalg.norm(garr[:, :2] - garr[n][None, :2], axis=1) <= 5.0] = False
    return g_true, g_res


def record(g_true, g_res):
    return sorted([(g, k) for k, g in enumerate(g_true)] + [(g, -1) for g in g_res], key=lambda x: -x[0])


def validate():
    import common_verif as V
    stored = json.load(open(os.path.join(xc.WORK, "verif_2e", "check4_parts_verif", "c_thin_ladder_4_net0.json")))
    S = V.setup("ladder_4", "thin", 0)
    err = 0.0
    for wi, w in enumerate(V.windows(S, "ladder_4", "thin", 0, 4)):
        srcs = [np.asarray(s, float) for s in w["srcs"]]
        R_or = oracle_residual_refined(S["net"], w, srcs, w["taus"])
        gt, gr = ceil_scores(S["PD"], S["net"].garr, w["MFn"], R_or, srcs)
        st = stored["windows"][wi]
        err = max(err, max([abs(a - b) for a, b in zip(gt, st["g_true"])] + [abs(a - b) for a, b in zip(gr, st["g_resid"])]
                           + [0.0]), float(len(gr) != len(st["g_resid"])))
        print(wi, len(srcs), np.round(gt, 3), np.round(gr, 3), "stored", np.round(st["g_true"], 3), np.round(st["g_resid"], 3))
    print("max abs diff vs verif_2e check4 parts:", err)


def run(regime, nets, nwin, path):
    res = json.load(open(path)) if os.path.exists(path) else {"rows": []}
    done = {(r["cell"], r["net"]) for r in res["rows"]}
    grid = xc.make_grid()
    garr = np.asarray(grid)
    for cell in ("ladder_4", "ladder_6", "unilateral_c40"):
        for j in nets:
            if (cell, j) in done:
                continue
            t0 = time.time()
            nc = xc.NetCfg(cell, j, regime)
            ops0 = xc.Ops(nc.STA, nc.corr, grid)
            ORA = L.Gradient(2.8, 0.04, 1.73)
            opsO = L.ModelOps(nc.STA, grid, ORA, L.corr_model(nc.STA, nc.cat, ORA))
            st = xc.window_stream(nc, 3006)
            for wi in range(nwin):
                w = next(st)
                R_or = oracle_residual_refined(ops0.net, w, w["truth"], w["taus"])
                g1 = ceil_scores(ops0.PD, garr, w["MFn"], R_or, w["truth"])
                g2 = ceil_scores(opsO.PD, garr, w["MFn"], R_or, w["truth"])
                res["rows"].append({"cell": cell, "net": j, "w": wi, "k": len(w["truth"]),
                                    "rec": {"CEIL-GD": record(*g1), "CEIL-GD-O": record(*g2)},
                                    "g": {"HOM": g1, "ORA": g2}})
            with open(path + ".tmp", "w") as fh:
                json.dump(res, fh)
            os.replace(path + ".tmp", path)
            print(regime, cell, j, "done", time.strftime("%H:%M:%S"), f"{time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "validate":
        validate()
    else:
        run(sys.argv[1], [int(x) for x in sys.argv[2].split(",")], int(sys.argv[3]), sys.argv[4])
