"""HARPA (Shi et al., Nat. Commun. 17:8140, 2026; DaDaCheng/phase_association at 7adf666) on picks in the atlas frame.

HARPA places a set of candidate events (locations and origin times), predicts their arrival times at every station
and moves the events by stochastic gradient Langevin dynamics (SGLD) to minimise a Wasserstein distance between the
predicted and the picked arrival times of each station; picks are then assigned to events by linear sum assignment
with a time-residual cutoff. It uses no amplitudes and no polarity.

Travel times (card: "the same 1-D model for every method that takes one"). HARPA's constant-speed model takes a
homogeneous medium; its neural-field model takes any travel-time function. With fixed ``one_d`` True (default) the
wrapper replaces the constant-speed arrival time by the frame's 1-D model, as a differentiable bilinear table of
:func:`ttmt.traveltime.first_arrival` (P and S) in epicentral distance and source depth (step ``table_h`` km), so
HARPA gets the model PyOcto, GaMMA and REAL get. The search, the loss, the assignment and every other step are
HARPA's own. With ``one_d`` False HARPA runs as published, with the slowness-averaged velocities of
:meth:`ttmt.associate.velocity.Velocity1D.mean_velocities` (also the velocities of its DBSCAN pre-clustering in both
modes, since that step only scales station distances).

Reproducibility. HARPA draws with the global ``random`` and ``numpy`` generators (station and phase of every SGLD
step) and seeds torch per slice. With more than one CPU it runs the DBSCAN slices in a process pool. The wrapper always
passes ``ncpu = 1`` (cannot be overridden), seeds ``random``, ``numpy`` and ``torch`` with SEED before each call, and
uses P and S picks (HARPA's default ``S_phase`` is False).

Tuned (space.HARPA): as_min_picks -> min_peak_pre_event (picks per event, P and S together), as_max_time_residual ->
max_time_residual (s; the analogue of PyOcto's tolerance, same range), as_dbscan_eps -> dbscan_eps (s, the time scale
of HARPA's DBSCAN slicing; same range as GaMMA's after ATLAS-0 amendment 4). Fixed: the "fast association" schedule of
HARPA's README (1000 SGLD steps before and 1000 after the decay), HARPA's defaults otherwise.
"""
import random

import numpy as np

from .base import detections, empty_result, from_timestamp, quiet, to_datetime

FIXED = {"one_d": True, "table_h": 0.5, "epoch_before_decay": 1000, "epoch_after_decay": 1000,
         "dbscan_min_samples": 1, "lr": 0.01, "noise": 1e-3, "lr_decay": 0.1, "noisy_pick": True, "LSA": True,
         "wasserstein_p": 2, "second_adjust": True, "n_event_max_rate": 1, "denoise_rate": 0}
TUNED = ("min_peak_pre_event", "max_time_residual", "dbscan_eps")
SEED = 42


class TravelTimeTable:
    """P and S first arrivals of a 1-D model on a (distance, depth) grid, read by differentiable bilinear
    interpolation in torch. Outside the grid the nearest edge value is used (HARPA's search box is padded by 5%)."""

    def __init__(self, velocity, rmax, zmax, h=0.5):
        from .. import traveltime as TT
        self.h = float(h)
        self.r = np.arange(0.0, float(rmax) + 2 * self.h, self.h)
        self.z = np.arange(0.0, float(zmax) + 2 * self.h, self.h)
        model = list(velocity.layers)
        X = np.maximum(self.r, 1e-3)
        from .velocity import table_depth          # no node on a layer top (E65)
        zq = [table_depth(model, z) for z in self.z]
        self.tab = np.stack([np.stack([TT.first_arrival(model, q, X, ratio) for q in zq])
                             for ratio in (1.0, velocity.vpvs)])          # (2, nz, nr)
        if not np.all(np.isfinite(self.tab)):
            raise ValueError("travel-time table has non-finite entries")
        self._t = None

    def __call__(self, r, z, phase):
        import torch
        if self._t is None:
            self._t = torch.tensor(self.tab, dtype=torch.float64)
        T = self._t[phase].to(r.dtype)
        fr = torch.clamp(r / self.h, 0, len(self.r) - 1 - 1e-6)
        fz = torch.clamp(z / self.h, 0, len(self.z) - 1 - 1e-6)
        i, k = fr.floor().long(), fz.floor().long()
        a, b = fr - i, fz - k
        return ((1 - a) * (1 - b) * T[k, i] + a * (1 - b) * T[k, i + 1] + (1 - a) * b * T[k + 1, i]
                + a * b * T[k + 1, i + 1])


def model_class(table):
    """HARPA's constant-speed AssignmentModel with the arrival time of the 1-D table."""
    import torch
    from harpa.model import AssignmentModel

    class AssignmentModel1D(AssignmentModel):
        def arrival_time(self, station_index, phase="P"):
            locs, t0 = self.loc_and_time()
            s = self.station_locs[station_index]
            r = torch.sqrt(((locs[:, :2] - s[:2]) ** 2).sum(dim=1) + 1e-12)
            return table(r, locs[:, 2], 0 if phase == "P" else 1).to(t0.dtype) + t0

    return AssignmentModel1D


def config_for(frame, cfg):
    fx = frame.fixed_for("harpa", FIXED)
    vp, vs = frame.velocity.mean_velocities(frame.zlim[1])
    c = {"x(km)": list(frame.xlim), "y(km)": list(frame.ylim), "z(km)": list(frame.zlim), "vel": {"P": vp, "S": vs},
         "P_phase": True, "S_phase": True, "DBSCAN": True, "time_before": float(frame.time_before),
         "min_peak_pre_event": int(cfg["as_min_picks"]), "max_time_residual": float(cfg["as_max_time_residual"]),
         "dbscan_eps": float(cfg["as_dbscan_eps"])}
    c.update({k: v for k, v in fx.items() if k not in ("one_d", "table_h")})
    c["ncpu"] = 1
    return c, fx


def _posix(t):
    """POSIX seconds of HARPA's origin time (an obspy UTCDateTime, whose .timestamp is a float)."""
    return float(t.timestamp)


def _seed():
    import torch
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)


def associate(frame, picks, cfg):
    import pandas as pd
    import harpa.harpa as H
    if len(picks) < int(cfg["as_min_picks"]):
        return empty_result(picks)
    c, fx = config_for(frame, cfg)
    codes = frame.codes
    pk = pd.DataFrame({"id": [codes[i] for i in picks.sta],
                       "timestamp": pd.to_datetime([to_datetime(t).replace(tzinfo=None) for t in picks.t]),
                       "type": ["p" if ph == 0 else "s" for ph in picks.phase], "prob": picks.p.astype(float)},
                      index=np.arange(len(picks)))
    st = pd.DataFrame({"id": codes, "x(km)": frame.stations.xyz[:, 0], "y(km)": frame.stations.xyz[:, 1],
                       "z(km)": frame.stations.xyz[:, 2]})
    original = H.AssignmentModel
    if fx["one_d"]:
        pad = 0.1 * max(frame.diagonal(), frame.zlim[1])
        H.AssignmentModel = model_class(TravelTimeTable(frame.velocity, frame.diagonal() + 2 * pad,
                                                        frame.zlim[1] + pad, fx["table_h"]))
    try:
        _seed()
        with quiet():
            out, cat = H.association(pk, st, c, verbose=0)
    finally:
        H.AssignmentModel = original
    if cat is None or len(cat) == 0:
        return empty_result(picks)
    row = {int(e): k for k, e in enumerate(cat["event_index"])}
    rows = [(r["x(km)"], r["y(km)"], r["z(km)"], from_timestamp(_posix(r["time"]))) for _, r in cat.iterrows()]
    return detections(rows, picks, [(int(i), row[int(e)]) for i, e in zip(out.index, out["event_index"])
                                    if e is not None and e == e and int(e) >= 0 and int(e) in row])
