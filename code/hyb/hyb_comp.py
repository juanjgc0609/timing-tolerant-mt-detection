"""
hyb_comp.py - the pick-based comparators given the same timing as the hybrid (Ridgecrest CAL, development only):
PhaseNet picks corrected by the pick-based station terms of hyb_dev2 (pick time minus term), PyOcto with a 1-D
velocity model (VelocityModel1D built from the Hadley-Kanamori layers, VS = VP / 1.73) and GaMMA with its eikonal 1-D
model, over the configuration grids of Task 2A.  So that a gain of the hybrid cannot come from a better velocity
model or station terms that the comparators did not get.
    python hyb_comp.py SUBNET OUT.json
"""
import json
import os
import sys
from datetime import timedelta

import numpy as np

sys.path[:0] = ["/home/claude/work/hyb/pydeps_pyrocko", "/home/claude/work/r0", "/home/claude/work/explore_contrib",
                "/home/claude/work", "/home/claude/work/gamma_src", "/home/claude/work/r0_verif/pydeps"]
import r0_lib as R  # noqa: E402
import r0_picks as RP  # noqa: E402
import r0_run as RR  # noqa: E402
import pandas as pd  # noqa: E402

import c18_tt as TT  # noqa: E402
import hyb_dev2 as H2  # noqa: E402
from task1r_baselines import GAMMA_GRID, PYOCTO_GRID  # noqa: E402

VPVS = 1.73
MODEL_PATH = "/home/claude/work/hyb/pyocto_hk.model"   # Task R1 sets its own path for JMA2001
POINTS = None   # Task R1: (depth, VP) points of a piecewise-linear model, given to PyOcto and GaMMA as such


def layers_df(model=None, zmax=40.0):
    import c18j_layered as LJ
    if model is None and POINTS is not None:
        z = [p[0] for p in POINTS if p[0] < zmax] + [zmax]
        vp = list(np.interp(z, [p[0] for p in POINTS], [p[1] for p in POINTS]))
        return pd.DataFrame({"depth": z, "vp": vp, "vs": np.array(vp) / VPVS})
    model = LJ.MODEL if model is None else model
    z, vp = [], []
    for k, (top, v) in enumerate(model):
        bot = model[k + 1][0] if k + 1 < len(model) else zmax
        z += [top, bot - 1e-3]
        vp += [v, v]
    return pd.DataFrame({"depth": z, "vp": vp, "vs": np.array(vp) / VPVS})


def pyocto_1d(STA, picks, n_picks, tol):
    import pyocto
    if not os.path.exists(MODEL_PATH):
        pyocto.VelocityModel1D.create_model(layers_df(), 0.25, 80.0, 40.0, MODEL_PATH)
    rows = [(f"S{i:02d}", b, RP.EPOCH.timestamp() + t) for i, pl in enumerate(picks) for t, b, _h in pl]
    if len(rows) < n_picks:
        return []
    pk = pd.DataFrame(rows, columns=["station", "phase", "time"]).sort_values("time").reset_index(drop=True)
    st = pd.DataFrame({"id": [f"S{i:02d}" for i in range(len(STA))], "x": STA[:, 0], "y": STA[:, 1], "z": STA[:, 2]})
    vm = pyocto.VelocityModel1D(MODEL_PATH, tolerance=float(tol))
    k = max(1, int(n_picks) // 3)
    assoc = pyocto.OctoAssociator(xlim=RP.XR, ylim=RP.XR, zlim=RP.ZR, velocity_model=vm, time_before=60.0,
                                  n_picks=int(n_picks), n_p_picks=k, n_s_picks=k, n_p_and_s_picks=k, n_threads=1)
    try:
        events, _a = assoc.associate(pk, st)
    except Exception as e:  # noqa: BLE001
        print("pyocto failed", n_picks, tol, repr(e)[:200], flush=True)
        return []
    if events is None or len(events) == 0:
        return []
    return [[float(r["x"]), float(r["y"]), float(r["z"]), float(r["time"]) - RP.EPOCH.timestamp(), float(r["picks"]), np.nan]
            for _, r in events.iterrows()]


def gamma_1d(STA, picks, min_picks, sigma):
    from gamma.utils import association
    rows = [(f"S{i:02d}", RP.EPOCH + timedelta(seconds=t), b.lower(), h) for i, pl in enumerate(picks) for t, b, h in pl]
    if len(rows) < min_picks:
        return []
    pk = pd.DataFrame(rows, columns=["id", "timestamp", "type", "prob"])
    pk["timestamp"] = pd.to_datetime(pk["timestamp"]).dt.tz_localize(None)
    pk["amp"] = 1.0
    pk = pk.sort_values("timestamp").reset_index(drop=True)
    st = pd.DataFrame({"id": [f"S{i:02d}" for i in range(len(STA))], "x(km)": STA[:, 0], "y(km)": STA[:, 1],
                       "z(km)": STA[:, 2]})
    df = layers_df()
    cfg = {"dims": ["x(km)", "y(km)", "z(km)"], "x(km)": list(RP.XR), "y(km)": list(RP.XR), "z(km)": list(RP.ZR),
           "vel": {"p": 6.0, "s": 6.0 / VPVS}, "method": "BGMM", "oversample_factor": 5, "use_amplitude": False,
           "use_dbscan": True, "dbscan_eps": 10.0, "dbscan_min_samples": 3, "min_picks_per_eq": int(min_picks),
           "max_sigma11": float(sigma), "max_sigma22": 1.0, "max_sigma12": 1.0,
           "bfgs_bounds": ((RP.XR[0] - 1, RP.XR[1] + 1), (RP.XR[0] - 1, RP.XR[1] + 1), RP.ZR, (None, None)),
           "covariance_prior": [5.0, 5.0], "initial_points": [1, 1, 1],
           "eikonal": {"vel": {"z": df["depth"].tolist(), "p": df["vp"].tolist(), "s": df["vs"].tolist()}, "h": 0.5,
                       "xlim": list(RP.XR), "ylim": list(RP.XR), "zlim": list(RP.ZR)}}
    np.random.seed(42)
    try:
        cat, _assign = association(pk, st, cfg, method="BGMM")
    except Exception as e:  # noqa: BLE001
        print("gamma failed", min_picks, sigma, repr(e)[:200], flush=True)
        return []
    out = []
    for e in (cat or []):
        t = (pd.Timestamp(e["time"]).tz_localize(None) - pd.Timestamp(RP.EPOCH).tz_localize(None)).total_seconds()
        out.append([float(e["x(km)"]), float(e["y(km)"]), float(e["z(km)"]), float(t), float(e.get("num_picks", 0)), np.nan])
    return out


def main(sub, outp):
    lo, hi = R.T_CAL
    os.chdir("/home/claude/work/hyb")
    net = RR.Net(sub, 940.0)
    pk = json.load(open(f"/home/claude/work/hyb/picks_cal_{sub}.json"))
    pa = H2.pick_arrays(pk["picks"])
    terms = H2.pick_terms(net, pa, lo, hi)
    picks_c = [sorted((t - terms[i, 0 if b == "P" else 1], b, h) for t, b, h in pl) for i, pl in enumerate(pk["picks"])]
    res = {"terms": terms.tolist(), "pyocto_1d": {}, "gamma_1d": {}}
    for cfg in PYOCTO_GRID:
        res["pyocto_1d"][json.dumps(cfg)] = pyocto_1d(net.STA, picks_c, **cfg)
    for cfg in GAMMA_GRID:
        res["gamma_1d"][json.dumps(cfg)] = gamma_1d(net.STA, picks_c, **cfg)
    json.dump(res, open(outp, "w"))
    import hyb_dev as H
    cat = net.cat
    for key in ("pyocto_1d", "gamma_1d"):
        best = max(((H.evaluate([e[:5] for e in v if lo <= e[3] < hi], cat, lo, hi), c) for c, v in res[key].items()),
                   key=lambda x: x[0][("LOC+T", 30)])
        print(sub, key, "best", best[1], {f"{u} {r}": round(x, 3) for (u, r), x in best[0].items()}, flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
