"""C6: physics candidate - 1D velocity models inside BOMP-PK and MT-FW.  SEED 999 ONLY.

Networks: task2a_confirm.network(999, cell, j) for j in 20..29 (not the 2A calibration networks 0-4, not the C1-C3
pilot networks 10-15).  Catalogue: xc.NetCfg (task2a catalogue_real: 20 events, errors 1 km / 0.1 s, statics sd 0.3 s).
Windows: xc.window_stream with TAG 3006 (never used before), draw_sources (the 2A window design), 10 per network.
Generator: gen_realistic default (1D gradient VS0 2.8, g 0.04, vpvs 1.73, variable pulses, coda, statics).

Models (each with its own station corrections; kmax 5, 5 km exclusion and hit tolerance as registered):
  HOM     homogeneous VS 3.0, vpvs 1.73, registered +SC corrections (= BOMP-PK+SC / MT-FW+SC)
  ORA     generator's gradient (2.8, 0.04, 1.73), curved rays and ray-direction polarisations; +SC re-referenced
  ORA-TS  ORA with the TRUE statics as corrections (ceiling of the travel-time correction)
  CAL-H   homogeneous (vs, vpvs) fitted by JHD on the 20 catalogue events; +SC re-referenced
  CAL-G   gradient (v0, g, vpvs) fitted by JHD; +SC re-referenced
  CAL-2L  two-layer (v1, v2, h, vpvs) fitted by JHD (mis-specified family); +SC re-referenced
Ceilings (genie, as verif_2e check 4 but with the oracle deflation at the TRUE arrival times, see REPORT):
  CEIL-GD   true events scored by the max HOM BOMP-PK block statistic within 5 km; false = greedy maxima (5 km
            exclusion, <= 5) of the HOM statistic on the residual after an oracle joint LS deflation of the true events
  CEIL-GD-O same with the ORA statistic
    python c6_run.py thin 20,21,...,29 10 out/c6_thin.json
"""
import json
import os
import sys
import time

import numpy as np

import c6_lib as L
import xc

E = xc.E
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
TAG = 3006
KMAX = 5


def near(garr, src):
    return np.where(np.linalg.norm(garr[:, :2] - np.asarray(src)[None, :2], axis=1) <= 5.0)[0]


def ceiling(PD, garr, MFn, R_or, truth):
    e_full = PD.energy(MFn, True)[0] / 2.0
    g_true = [float(e_full[near(garr, s)].max()) for s in truth]
    e_res = PD.energy(R_or, True)[0] / 2.0
    emr = e_res.max(axis=1)
    allowed = np.ones(len(garr), bool)
    for s in truth:
        allowed[near(garr, s)] = False
    g_res = []
    for _ in range(5):
        cand = np.where(allowed, emr, -np.inf)
        n = int(np.argmax(cand))
        if not np.isfinite(cand[n]):
            break
        g_res.append(float(cand[n]))
        allowed[np.linalg.norm(garr[:, :2] - garr[n][None, :2], axis=1) <= 5.0] = False
    rec = sorted([(g, k) for k, g in enumerate(g_true)] + [(g, -1) for g in g_res], key=lambda x: -x[0])
    return rec


def oracle_residual(ops0, w):
    net = ops0.net
    k = len(w["truth"])
    if k == 0:
        return w["MFn"]
    win = E.Win(w["D"], w["MFn"], w["sd"])
    A = xc.true_arrivals(w, net.ns)
    P = np.zeros((k, net.nv))
    for s in range(k):
        src = w["truth"][s]
        d = np.linalg.norm(net.STA - src, axis=1)
        T0 = d[None, :] * E.VINV[:, None] + w["taus"][s] + net.corrT
        P[s, :3] = src
        P[s, 3] = w["taus"][s]
        P[s, 4:] = (A[s] - T0).ravel()
    m, _ = E.convex_step(net, win, P, np.zeros((k, 6)), 0.0)
    return w["MFn"] - E.predict_sampled(net, win, P, m)


def build_models(nc):
    Ec = L.cat_energy(nc.cat)
    out, info = {}, {}
    ORA = L.Gradient(2.8, 0.04, 1.73)
    out["HOM"] = (L.Homog(), nc.corr)
    out["ORA"] = (ORA, L.corr_model(nc.STA, nc.cat, ORA, Ec))
    out["ORA-TS"] = (ORA, np.asarray(nc.statics, float))
    for name, fam in (("CAL-H", "homog"), ("CAL-G", "gradient"), ("CAL-2L", "twolayer")):
        t = time.time()
        mdl, inf = L.jhd_fit(nc.STA, nc.cat, fam, Ec, return_info=True)
        inf["sec"] = time.time() - t
        out[name] = (mdl, L.corr_model(nc.STA, nc.cat, mdl, Ec))
        info[name] = inf
    return out, info


def run(regime, nets, nwin, path):
    res = json.load(open(path)) if os.path.exists(path) else {"rows": [], "nets": []}
    done = {(n["cell"], n["net"]) for n in res["nets"]}
    grid = xc.make_grid()
    for cell in CELLS:
        for j in nets:
            if (cell, j) in done:
                continue
            t0 = time.time()
            nc = xc.NetCfg(cell, j, regime)
            models, info = build_models(nc)
            ops0 = xc.Ops(nc.STA, nc.corr, grid)
            mops = {name: L.ModelOps(nc.STA, grid, mdl, corr) for name, (mdl, corr) in models.items()}
            netinfo = {"cell": cell, "net": j, "regime": regime, "calib": info,
                       "params": {k: list(m.params) for k, (m, _c) in models.items()},
                       "corr_err": {k: float(np.median(np.abs(c - nc.statics))) for k, (_m, c) in models.items()},
                       "t_setup": time.time() - t0}
            st = xc.window_stream(nc, TAG)
            rows = []
            for wi in range(nwin):
                w = next(st)
                tr = w["truth"]
                row = {"cell": cell, "net": j, "w": wi, "k": len(tr), "truth": [t.tolist() for t in tr],
                       "taus": w["taus"], "rec": {}, "loc": {}, "sec": {}}
                for name, ops in mops.items():
                    t = time.time(); r, l = L.bomp_pk_m(ops, w["MFn"], tr, KMAX); row["sec"]["BOMP-PK:" + name] = time.time() - t
                    row["rec"]["BOMP-PK:" + name], row["loc"]["BOMP-PK:" + name] = r, l
                    t = time.time(); r, l = L.mtfw_dense(ops, w["MFn"], tr, KMAX); row["sec"]["MT-FW:" + name] = time.time() - t
                    row["rec"]["MT-FW:" + name], row["loc"]["MT-FW:" + name] = r, l
                t = time.time()
                R_or = oracle_residual(ops0, w)
                row["rec"]["CEIL-GD"] = ceiling(mops["HOM"].PD, grid_arr(grid), w["MFn"], R_or, tr)
                row["rec"]["CEIL-GD-O"] = ceiling(mops["ORA"].PD, grid_arr(grid), w["MFn"], R_or, tr)
                row["sec"]["CEIL"] = time.time() - t
                rows.append(row)
            netinfo["t_total"] = time.time() - t0
            res["rows"].extend(rows)
            res["nets"].append(netinfo)
            with open(path + ".tmp", "w") as fh:
                json.dump(res, fh)
            os.replace(path + ".tmp", path)
            print(regime, cell, j, "done", time.strftime("%H:%M:%S"), f"{netinfo['t_total']:.0f} s",
                  {k: [round(x, 3) for x in v] for k, v in netinfo["params"].items() if k.startswith("CAL")}, flush=True)


_GA = {}


def grid_arr(grid):
    key = id(grid)
    if key not in _GA:
        _GA[key] = np.asarray(grid, float)
    return _GA[key]


if __name__ == "__main__":
    run(sys.argv[1], [int(x) for x in sys.argv[2].split(",")], int(sys.argv[3]), sys.argv[4])
