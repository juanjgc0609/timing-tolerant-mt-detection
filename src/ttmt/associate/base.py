"""The frame every associator wrapper runs in, and the conversions between the atlas schema and the programs.

Frame: the station set (``ttmt.atlas.schema.Stations``, local km, z down), the search volume (km), the 1-D model,
the time padding of the program's sliding windows, a working folder for travel-time tables, and the fixed (untuned)
parameters of each program. Every wrapper takes ``Picks`` in seconds from the record start and returns
``(Detections, Picks)``: the detections in the same frame, scored by their number of associated picks (as R0 and
R1), and the input picks with ``det`` set to the index of the detection each pick was associated with (-1 if none).
"""
import os
import tempfile
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone

import numpy as np

from ..atlas.schema import Detections, Picks, Stations
from .velocity import Velocity1D

EPOCH = datetime(2000, 1, 1, tzinfo=timezone.utc)     # record time 0 for programs that need absolute times
PHASE_NAMES = ("P", "S")


@dataclass
class Frame:
    stations: Stations
    xlim: tuple
    ylim: tuple
    zlim: tuple
    velocity: Velocity1D
    time_before: float = 60.0
    workdir: str = None
    fixed: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.workdir is None:
            self.workdir = tempfile.mkdtemp(prefix="ttmt_assoc_")
        os.makedirs(self.workdir, exist_ok=True)

    @classmethod
    def around(cls, stations, velocity, pad=10.0, zmax=30.0, **kw):
        """the station box padded by pad km on every side, depth 0 to zmax."""
        xy = stations.xyz[:, :2]
        lo, hi = xy.min(axis=0) - pad, xy.max(axis=0) + pad
        return cls(stations, (float(lo[0]), float(hi[0])), (float(lo[1]), float(hi[1])), (0.0, float(zmax)),
                   velocity, **kw)

    @property
    def codes(self):
        return self.stations.code.tolist()

    def diagonal(self):
        return float(np.hypot(self.xlim[1] - self.xlim[0], self.ylim[1] - self.ylim[0]))

    def fixed_for(self, program, defaults):
        return {**defaults, **self.fixed.get(program, {})}


def empty_result(picks):
    return Detections(np.zeros(0), np.zeros((0, 3)), np.zeros(0)), with_det(picks, -np.ones(len(picks), int))


def with_det(picks, det):
    return Picks(picks.sta, picks.phase, picks.t, picks.p, np.asarray(det, int))


def detections(rows, picks, assign):
    """rows: [(x, y, z, t)] per event; assign: (pick index, event row) pairs. Score = associated picks per event."""
    det = -np.ones(len(picks), int)
    for i, e in assign:
        det[int(i)] = int(e)
    n = np.bincount(det[det >= 0], minlength=len(rows)) if len(rows) else np.zeros(0)
    a = np.asarray(rows, float).reshape(-1, 4)
    return Detections(a[:, 3], a[:, :3], n.astype(float)), with_det(picks, det)


def to_datetime(t):
    return EPOCH + timedelta(seconds=float(t))


def from_timestamp(x):
    """seconds from EPOCH of a pandas/py datetime, str, or POSIX seconds."""
    import pandas as pd
    if isinstance(x, (int, float, np.floating)):
        return float(x) - EPOCH.timestamp()
    ts = pd.Timestamp(x)
    if ts.tzinfo is not None:
        ts = ts.tz_convert(None)
    return (ts - pd.Timestamp(EPOCH).tz_convert(None)).total_seconds()


class quiet:
    """silence a program's prints (stdout and stderr, Python level) for the duration of a call."""

    def __enter__(self):
        import contextlib
        import io
        self._stack = contextlib.ExitStack()
        self._stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
        self._stack.enter_context(contextlib.redirect_stderr(io.StringIO()))
        return self

    def __exit__(self, *exc):
        self._stack.close()
        return False
