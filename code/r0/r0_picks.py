"""
r0_picks.py - Task R0 comparators: PhaseNet (original model of the repository, as real_picker.py) on the 100 Hz
velocity traces of the subnetwork, threshold 0.3 for P and S, minimum peak distance 0.5 s (PhaseNet defaults), no time
shift; GaMMA (BGMM, gamma_src of Task 2A) and PyOcto (0.2.0) on those picks, over the configuration grids of Task 2A
(task1r_baselines.GAMMA_GRID, PYOCTO_GRID), in the R0 frame and the R0 homogeneous medium (VP 6.0, VS 6.0/1.73).
Score of an associated event = number of picks associated to it.  Each associator runs on the raw picks and on
the picks corrected by the +SC station corrections of r0_run (pick time minus correction), as in Task 2A.
    python r0_picks.py SEGMENT SUBNET OUT.json
"""
import json
import os
import sys
import time
from datetime import datetime, timedelta, timezone

import numpy as np

import r0_lib as R

sys.path.insert(0, os.path.join(R.WORK, "gamma_src"))
import pandas as pd  # noqa: E402
from scipy.signal import find_peaks  # noqa: E402

import real_picker as RP  # noqa: E402
from task1r_baselines import GAMMA_GRID, PYOCTO_GRID  # noqa: E402

THR, MPD = 0.3, 0.5
CHUNK, PADS = 30000, 1500
XR = (-21.4, 21.4)
ZR = (0.0, 15.0)
EPOCH = datetime(2019, 7, 4, 19, 0, 0, tzinfo=timezone.utc)


def phasenet_hour(raw):
    """raw (ns, 3, n) m/s at 100 Hz, components E, N, Z (up).  Probabilities (ns, n, 3) by overlapping chunks."""
    sess, model, _v = RP.load_phasenet()
    ns, _, n = raw.shape
    out = np.zeros((ns, n, 3))
    for i in range(ns):
        X = raw[i].T[:, None, :]                                                    # (n, 1, 3)
        a = 0
        while a < n:
            lo, hi = max(0, a - PADS), min(n, a + CHUNK + PADS)
            seg = RP.normalize_long(X[lo:hi]).astype(np.float32)[None]
            P = sess.run(model.preds, feed_dict={model.X: seg, model.drop_rate: 0, model.is_training: False})
            P = P[0, :, 0, :]
            b = min(n, a + CHUNK)
            out[i, a:b] = P[a - lo:b - lo]
            a = b
    return out


def picks_from(P):
    mpd = int(round(MPD * 100))
    picks = []
    for pr in P:
        pl = []
        for k, b in ((1, "P"), (2, "S")):
            idx, prop = find_peaks(pr[:, k], height=THR, distance=mpd)
            pl += [(float(j / 100.0), b, float(h)) for j, h in zip(idx, prop["peak_heights"])]
        pl.sort()
        picks.append(pl)
    return picks


def gamma_run(STA, picks, min_picks, sigma):
    from gamma.utils import association
    rows = [(f"S{i:02d}", EPOCH + timedelta(seconds=t), b.lower(), h) for i, pl in enumerate(picks) for t, b, h in pl]
    if len(rows) < min_picks:
        return []
    pk = pd.DataFrame(rows, columns=["id", "timestamp", "type", "prob"])
    pk["timestamp"] = pd.to_datetime(pk["timestamp"]).dt.tz_localize(None)
    pk["amp"] = 1.0
    pk = pk.sort_values("timestamp").reset_index(drop=True)
    st = pd.DataFrame({"id": [f"S{i:02d}" for i in range(len(STA))], "x(km)": STA[:, 0], "y(km)": STA[:, 1],
                       "z(km)": STA[:, 2]})
    cfg = {"dims": ["x(km)", "y(km)", "z(km)"], "x(km)": list(XR), "y(km)": list(XR), "z(km)": list(ZR),
           "vel": {"p": R.T0.VP, "s": R.T0.VS}, "method": "BGMM", "oversample_factor": 5,
           "use_amplitude": False, "use_dbscan": True, "dbscan_eps": 10.0,
           "dbscan_min_samples": 3, "min_picks_per_eq": int(min_picks),
           "max_sigma11": float(sigma), "max_sigma22": 1.0, "max_sigma12": 1.0,
           "bfgs_bounds": ((XR[0] - 1, XR[1] + 1), (XR[0] - 1, XR[1] + 1), ZR, (None, None)),
           "covariance_prior": [5.0, 5.0], "initial_points": [1, 1, 1]}
    np.random.seed(42)
    try:
        cat, assign = association(pk, st, cfg, method="BGMM")
    except Exception as e:                                                          # noqa: BLE001
        print("gamma failed", min_picks, sigma, repr(e)[:200], flush=True)
        return []
    out = []
    for e in (cat or []):
        t = (pd.Timestamp(e["time"]).tz_localize(None) - pd.Timestamp(EPOCH).tz_localize(None)).total_seconds()
        out.append([float(e["x(km)"]), float(e["y(km)"]), float(e["z(km)"]), float(t), float(e.get("num_picks", 0)),
                    float(e.get("gamma_score", np.nan))])
    return out


def pyocto_run(STA, picks, n_picks, tol):
    import pyocto
    rows = [(f"S{i:02d}", b, EPOCH.timestamp() + t) for i, pl in enumerate(picks) for t, b, _h in pl]
    if len(rows) < n_picks:
        return []
    pk = pd.DataFrame(rows, columns=["station", "phase", "time"]).sort_values("time").reset_index(drop=True)
    st = pd.DataFrame({"id": [f"S{i:02d}" for i in range(len(STA))], "x": STA[:, 0], "y": STA[:, 1],
                       "z": STA[:, 2]})
    vm = pyocto.VelocityModel0D(p_velocity=R.T0.VP, s_velocity=R.T0.VS, tolerance=float(tol))
    k = max(1, int(n_picks) // 3)
    assoc = pyocto.OctoAssociator(xlim=XR, ylim=XR, zlim=ZR, velocity_model=vm, time_before=60.0,
                                  n_picks=int(n_picks), n_p_picks=k, n_s_picks=k, n_p_and_s_picks=k, n_threads=1)
    try:
        events, assignments = assoc.associate(pk, st)
    except Exception as e:                                                          # noqa: BLE001
        print("pyocto failed", n_picks, tol, repr(e)[:200], flush=True)
        return []
    out = []
    if events is None or len(events) == 0:
        return out
    for _, r in events.iterrows():
        out.append([float(r["x"]), float(r["y"]), float(r["z"]), float(r["time"]) - EPOCH.timestamp(),
                    float(r["picks"]), np.nan])
    return out


def main(segment, subnet, path):
    t0 = time.time()
    t_end = (R.T_CAL[1] if segment == "cal" else R.T_EVAL[1]) + 40.0
    names = R.SUBNETS[subnet]
    import r0_run as RR
    net = RR.Net(subnet, t_end)                      # same frame, stations, waveforms and +SC corrections as r0_run
    STA, raw, corr = net.STA, net.raw, net.corr_sc
    P = phasenet_hour(raw)
    picks = picks_from(P)
    t_pn = time.time() - t0
    res = {"segment": segment, "subnet": subnet, "names": names, "picks": picks, "t_phasenet": t_pn,
           "n_picks": {b: int(sum(1 for pl in picks for p in pl if p[1] == b)) for b in "PS"}}
    print("phasenet", round(t_pn), "s", res["n_picks"], flush=True)
    picks_sc = [sorted((t - corr[i, 0 if b == "P" else 1], b, h) for t, b, h in pl) for i, pl in enumerate(picks)]
    res["corr_sc"] = corr.tolist()
    for tag, pk in (("", picks), ("_sc", picks_sc)):
        res["gamma" + tag], res["pyocto" + tag] = {}, {}
        for cfg in GAMMA_GRID:
            res["gamma" + tag][json.dumps(cfg)] = gamma_run(STA, pk, **cfg)
            print("gamma" + tag, cfg, len(res["gamma" + tag][json.dumps(cfg)]), round(time.time() - t0), flush=True)
        for cfg in PYOCTO_GRID:
            res["pyocto" + tag][json.dumps(cfg)] = pyocto_run(STA, pk, **cfg)
            print("pyocto" + tag, cfg, len(res["pyocto" + tag][json.dumps(cfg)]), round(time.time() - t0), flush=True)
    res["t_total"] = time.time() - t0
    with open(path + ".tmp", "w") as fh:
        json.dump(res, fh)
    os.replace(path + ".tmp", path)
    print("DONE", round(res["t_total"]), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
