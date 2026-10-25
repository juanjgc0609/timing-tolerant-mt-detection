"""PyOcto 0.2.0 (Münchmeyer, 2024) on picks in the atlas frame.

Tuned (space.PYOCTO): as_min_picks -> n_picks; as_tolerance -> the velocity model's tolerance (s); the card's
"association cutoff in [0.5, 2] s" -> pick_match_tolerance (s), the only association cutoff of PyOcto in seconds
(association_cutoff_distance is in km); as_min_interevent_time -> min_interevent_time (s; ATLAS-0 amendment 1,
ADR-0009; the fixed value applies when a configuration does not carry it). n_p_picks = n_s_picks = n_p_and_s_picks = max(1, n_picks // 3), as R0 and
R1, unless fixed. Fixed (Frame.fixed["pyocto"] may override): one thread (n_threads = 1), the 1-D model of the frame as a PyOcto
table, time_before from the frame, PyOcto's defaults otherwise.
"""
import os

import numpy as np

from .base import EPOCH, detections, empty_result, quiet

FIXED = {"one_d": True, "table_delta": 0.25, "table_xdist": None, "association_cutoff_distance": None,
         "min_node_size": 10.0, "min_node_size_location": 1.5, "min_interevent_time": 3.0, "max_pick_overlap": 4,
         "refinement_iterations": 3, "time_slicing": 1200.0, "min_pick_fraction": 0.25, "n_p_picks": None,
         "n_s_picks": None, "n_p_and_s_picks": None}
PASS = ("min_node_size", "min_node_size_location", "min_interevent_time", "max_pick_overlap", "refinement_iterations",
        "time_slicing", "min_pick_fraction")


def velocity_model(frame, tolerance, fx):
    import pyocto
    if not fx["one_d"]:
        vp, vs = frame.velocity.mean_velocities(frame.zlim[1])
        return pyocto.VelocityModel0D(p_velocity=vp, s_velocity=vs, tolerance=float(tolerance),
                                      association_cutoff_distance=fx["association_cutoff_distance"])
    xdist = fx["table_xdist"] or frame.diagonal() + 20.0
    zdist = float(frame.zlim[1]) + 5.0
    key = abs(hash((frame.velocity.layers, frame.velocity.vpvs, round(xdist, 3), zdist, fx["table_delta"])))
    path = os.path.join(frame.workdir, f"pyocto_{key:x}.model")
    frame.velocity.pyocto_model(path, xdist, zdist, fx["table_delta"])
    return pyocto.VelocityModel1D(path, tolerance=float(tolerance),
                                  association_cutoff_distance=fx["association_cutoff_distance"])


def associate(frame, picks, cfg):
    import pandas as pd
    import pyocto
    n = int(cfg["as_min_picks"])
    if len(picks) < n:
        return empty_result(picks)
    fx = frame.fixed_for("pyocto", FIXED)
    codes = frame.codes
    pk = pd.DataFrame({"station": [codes[i] for i in picks.sta], "phase": ["P" if ph == 0 else "S" for ph in picks.phase],
                       "time": EPOCH.timestamp() + picks.t}, index=np.arange(len(picks)))
    st = pd.DataFrame({"id": codes, "x": frame.stations.xyz[:, 0], "y": frame.stations.xyz[:, 1],
                       "z": frame.stations.xyz[:, 2]})
    k = max(1, n // 3)
    assoc = pyocto.OctoAssociator(
        xlim=tuple(frame.xlim), ylim=tuple(frame.ylim), zlim=tuple(frame.zlim),
        velocity_model=velocity_model(frame, cfg["as_tolerance"], fx), time_before=float(frame.time_before),
        pick_match_tolerance=float(cfg["as_pick_match_tolerance"]), n_picks=n,
        n_p_picks=fx["n_p_picks"] or k, n_s_picks=fx["n_s_picks"] or k, n_p_and_s_picks=fx["n_p_and_s_picks"] or k,
        n_threads=1, **{key: fx[key] for key in PASS if key != "min_interevent_time"},
        min_interevent_time=float(cfg.get("as_min_interevent_time", fx["min_interevent_time"])))
    with quiet():
        events, assignments = assoc.associate(pk, st)
    if events is None or len(events) == 0:
        return empty_result(picks)
    row = {int(e): r for r, e in enumerate(events["idx"])}
    rows = [(r["x"], r["y"], r["z"], float(r["time"]) - EPOCH.timestamp()) for _, r in events.iterrows()]
    return detections(rows, picks, [(int(p), row[int(e)]) for p, e in zip(assignments["pick_idx"],
                                                                           assignments["event_idx"])])
