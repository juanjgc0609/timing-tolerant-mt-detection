"""GaMMA (Zhu et al., 2022; AI4EPS/GaMMA at 80394dd) on picks in the atlas frame.

Reproducibility (handoff of the R1 pick stage): without "ncpu" in its configuration GaMMA sets
ncpu = max(1, min(clusters // 4, min(32, cpu_count - 1))) and, with ncpu > 1, shuffles the DBSCAN clusters with the
global numpy generator and runs them in a process pool whose workers are reseeded with 42. The result then depends
on the machine's core count (on 16 cores the GaMMA-1D event counts of R1 changed and a worker hung). This wrapper
always passes ncpu = 1 (serial loop over the clusters in DBSCAN order, the same on every machine), seeds numpy with 42
before each call, and hands GaMMA a fresh copy of its configuration, because GaMMA replaces config["eikonal"] by the
initialised solver in place (a reused dictionary would be initialised twice). ncpu cannot be overridden.

Tuned (space.GAMMA): as_min_picks -> min_picks_per_eq, as_max_sigma11 -> max_sigma11 (time standard deviation, s),
as_method -> method. Only BGMM runs: in GaMMA 80394dd method "GMM" raises AttributeError before any fit, because
BaseMixture._initialize_parameters (gamma/_base.py:140, init_params "centers", the only one gamma.utils.association
uses) reads self.covariance_prior, which only BayesianGaussianMixture sets, and association() builds GaussianMixture
without a covariance prior (E63). The wrapper refuses GMM with that reason instead of calling GaMMA, and ATLAS-0
amendment 3 pins the tuned method to BGMM (space.GAMMA).
Fixed (Frame.fixed["gamma"] may override): the values of R0 and R1 (code/r0/r0_picks.py), so the untuned comparators
of the thesis are one point of this space's neighbourhood. A fixed value None removes the entry (GaMMA then uses its
own default, e.g. its data-driven covariance prior).
"""
import copy

import numpy as np

from .base import detections, empty_result, from_timestamp, quiet, to_datetime

FIXED = {"use_amplitude": False, "use_dbscan": True, "dbscan_eps": 10.0, "dbscan_min_samples": 3,
         "oversample_factor": 5, "max_sigma22": 1.0, "max_sigma12": 1.0, "covariance_prior": [5.0, 5.0],
         "initial_points": [1, 1, 1], "eikonal_h": 0.5, "one_d": True}
SEED = 42


def config_for(frame, cfg):
    fx = frame.fixed_for("gamma", FIXED)
    vp, vs = frame.velocity.mean_velocities(frame.zlim[1])
    c = {"dims": ["x(km)", "y(km)", "z(km)"], "x(km)": list(frame.xlim), "y(km)": list(frame.ylim),
         "z(km)": list(frame.zlim), "vel": {"p": vp, "s": vs}, "method": cfg["as_method"],
         "min_picks_per_eq": int(cfg["as_min_picks"]), "max_sigma11": float(cfg["as_max_sigma11"]),
         "bfgs_bounds": ((frame.xlim[0] - 1, frame.xlim[1] + 1), (frame.ylim[0] - 1, frame.ylim[1] + 1),
                         tuple(frame.zlim), (None, None)), "ncpu": 1}
    c.update({k: copy.deepcopy(v) for k, v in fx.items() if k not in ("eikonal_h", "one_d")})
    c = {k: v for k, v in c.items() if v is not None}           # a fixed value None removes GaMMA's entry
    if fx["one_d"]:
        c["eikonal"] = frame.velocity.gamma_eikonal(frame.xlim, frame.ylim, frame.zlim, fx["eikonal_h"])
    c["ncpu"] = 1
    return c


GMM_BROKEN = ("GaMMA 80394dd cannot run method GMM: gamma/_base.py:140 reads covariance_prior, which only "
              "BayesianGaussianMixture has (E63; ATLAS-0 amendment 3 pins BGMM)")


def associate(frame, picks, cfg):
    if cfg["as_method"] != "BGMM":
        raise ValueError(GMM_BROKEN if cfg["as_method"] == "GMM" else f"unknown GaMMA method {cfg['as_method']!r}")
    import pandas as pd
    from gamma.utils import association
    if len(picks) < int(cfg["as_min_picks"]):
        return empty_result(picks)
    codes = frame.codes
    pk = pd.DataFrame({"id": [codes[i] for i in picks.sta],
                       "timestamp": [to_datetime(t).replace(tzinfo=None) for t in picks.t],
                       "type": ["p" if ph == 0 else "s" for ph in picks.phase],
                       "prob": picks.p.astype(float), "amp": 1.0}, index=np.arange(len(picks)))
    pk["timestamp"] = pd.to_datetime(pk["timestamp"])
    pk = pk.sort_values("timestamp", kind="stable")              # index keeps the input pick number
    st = pd.DataFrame({"id": codes, "x(km)": frame.stations.xyz[:, 0], "y(km)": frame.stations.xyz[:, 1],
                       "z(km)": frame.stations.xyz[:, 2]})
    c = config_for(frame, cfg)
    np.random.seed(SEED)
    with quiet():
        cat, assign = association(pk, st, c, method=c["method"])
    cat = cat or []
    row = {e["event_index"]: k for k, e in enumerate(cat)}
    rows = [(e["x(km)"], e["y(km)"], e["z(km)"], from_timestamp(e["time"])) for e in cat]
    return detections(rows, picks, [(i, row[e]) for i, e, _p in assign if e in row])
