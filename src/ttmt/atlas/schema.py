"""Common frame of the atlas: every quantity in local kilometres (x east, y north, z down) and seconds from the start
of the record. Truth comes from the generator v2 thread (semi-synthetic or synthetic) or from a catalogue; picks and
detections come from the adapters."""
from dataclasses import dataclass, field

import numpy as np

PHASES = ("P", "S")


@dataclass
class Stations:
    """station codes (n,) and coordinates (n, 3) km."""
    code: np.ndarray
    xyz: np.ndarray

    def __post_init__(self):
        self.code = np.asarray(self.code, dtype=str)
        self.xyz = np.asarray(self.xyz, float).reshape(-1, 3)
        if len(self.code) != len(self.xyz) or len(set(self.code.tolist())) != len(self.code):
            raise ValueError("station codes must be unique and match the coordinates")

    def index(self, codes):
        lut = {c: i for i, c in enumerate(self.code.tolist())}
        return np.array([lut[c] for c in codes], int)


@dataclass
class Truth:
    """target and neutral events of one record.

    Per event (m,): id, t (origin, s), xyz (m, 3) km, mag, target (bool; False = catalogued event of the noise record,
    neutral in the scoring), pair (pair id or -1), kmec (K_mec/E of its pair, nan if unknown).
    Arrivals: arr_t (m, n, 2) true arrival times (s) of P and S at every station, nan where absent; arr_t_obs the
    same with the injected timing error (defaults to arr_t). snr (m, n) per-station signal-to-noise ratio (nan where
    unknown). Record span [t0, t1) and the noise hours, for the block bootstrap.
    """
    id: np.ndarray
    t: np.ndarray
    xyz: np.ndarray
    mag: np.ndarray
    arr_t: np.ndarray
    snr: np.ndarray
    t0: float
    t1: float
    target: np.ndarray = None
    pair: np.ndarray = None
    kmec: np.ndarray = None
    arr_t_obs: np.ndarray = None
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        m = len(self.t)
        self.id = np.asarray(self.id)
        self.t = np.asarray(self.t, float)
        self.xyz = np.asarray(self.xyz, float).reshape(m, 3)
        self.mag = np.asarray(self.mag, float)
        self.arr_t = np.asarray(self.arr_t, float)
        if self.arr_t.ndim != 3 or self.arr_t.shape[0] != m or self.arr_t.shape[2] != 2:
            raise ValueError("arr_t must be (events, stations, 2)")
        self.snr = np.asarray(self.snr, float).reshape(self.arr_t.shape[:2])
        self.target = np.ones(m, bool) if self.target is None else np.asarray(self.target, bool)
        self.pair = -np.ones(m, int) if self.pair is None else np.asarray(self.pair, int)
        self.kmec = np.full(m, np.nan) if self.kmec is None else np.asarray(self.kmec, float)
        self.arr_t_obs = self.arr_t.copy() if self.arr_t_obs is None else np.asarray(self.arr_t_obs, float)

    @property
    def n_events(self):
        return len(self.t)

    @property
    def hours(self):
        return (self.t1 - self.t0) / 3600.0

    def subset_stations(self, idx):
        """the same truth seen by the stations idx only."""
        idx = np.asarray(idx, int)
        return Truth(self.id, self.t, self.xyz, self.mag, self.arr_t[:, idx], self.snr[:, idx], self.t0, self.t1,
                     self.target, self.pair, self.kmec, self.arr_t_obs[:, idx], dict(self.meta))


@dataclass
class Picks:
    """picks of one method: station index (k,), phase (k,) 0 = P, 1 = S, time t (k,) s, probability p (k,), and the
    detection each pick was associated with (k,), -1 if none."""
    sta: np.ndarray
    phase: np.ndarray
    t: np.ndarray
    p: np.ndarray
    det: np.ndarray = None

    def __post_init__(self):
        self.sta = np.asarray(self.sta, int)
        self.phase = np.asarray(self.phase, int)
        self.t = np.asarray(self.t, float)
        self.p = np.asarray(self.p, float)
        self.det = -np.ones(len(self.t), int) if self.det is None else np.asarray(self.det, int)
        if not (len(self.sta) == len(self.phase) == len(self.t) == len(self.p) == len(self.det)):
            raise ValueError("pick arrays must have the same length")

    def __len__(self):
        return len(self.t)


@dataclass
class Detections:
    """detections of one method: origin t (d,) s, location xyz (d, 3) km, score s (d,), higher is more confident."""
    t: np.ndarray
    xyz: np.ndarray
    s: np.ndarray
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        self.t = np.asarray(self.t, float)
        self.xyz = np.asarray(self.xyz, float).reshape(len(self.t), 3)
        self.s = np.asarray(self.s, float)
        if len(self.s) != len(self.t):
            raise ValueError("scores and origins must have the same length")

    def __len__(self):
        return len(self.t)

    def within(self, t0, t1):
        k = (self.t >= t0) & (self.t < t1)
        return Detections(self.t[k], self.xyz[k], self.s[k], dict(self.meta))


@dataclass
class MethodRun:
    """one method on one record and one station subset: its detections, its picks if it has any, the configuration
    and the minimum number of picks its associator requires (for the picker/associator attribution of misses)."""
    method: str
    det: Detections
    picks: Picks = None
    config: dict = field(default_factory=dict)
    min_picks: int = None
