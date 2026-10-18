"""Data layers of the real-data tasks R0 (Ridgecrest 2019) and R1 (Hi-net swarm 2012), as configuration objects.

Ports of code/r0/r0_lib.py and code/r1/r1_lib.py. The legacy modules set the medium and the origin range in
task0_common at import time; here they are fields of the dataset object. Paths come from ``ttmt.config``.
"""
import os
from dataclasses import dataclass, field

import numpy as np

from . import config
from .partition import Partition

DT = 1.0 / 16
NT = 640


@dataclass(frozen=True)
class Dataset:
    name: str
    folder: str                              # key of config.DATASETS
    t0_utc: str
    cal: tuple
    eval: tuple
    subnets: dict
    m_target: float
    vp: float = 6.0
    vpvs: float = 1.73
    taus: np.ndarray = field(default_factory=lambda: np.round(np.arange(2.0, 11.95, 0.1), 10), compare=False)
    own: tuple = (5.0, 8.0)
    step: float = 3.0
    kmax: int = 6
    depths: tuple = (3.0, 7.0, 11.0)
    m_terms: float = None

    @property
    def vs(self):
        return self.vp / self.vpvs

    @property
    def data(self):
        return str(config.dataset_dir(self.folder))

    @property
    def partition(self):
        return Partition({"cal": self.cal, "eval": self.eval})

    def segment(self, name):
        return self.cal if name == "cal" else self.eval

    def starts(self, segment):
        """r0_run.starts: window starts of a segment (step 3 s); a window owns origins in [5, 8) s of window time."""
        lo, hi = self.segment(segment)
        out, t = [], max(0.0, lo - self.own[0])
        while t + self.own[0] < hi:
            out.append(round(t, 6))
            t += self.step
        return out

    # ---------------------------------------------------------------------------------------- catalogue and frame
    def read_catalogue(self):
        raise NotImplementedError

    def frame(self):
        c = [e for e in self.read_catalogue() if self.cal[0] <= e["t"] < self.cal[1]]
        return float(np.mean([e["lat"] for e in c])), float(np.mean([e["lon"] for e in c]))

    @staticmethod
    def to_xy(lat, lon, centre):
        lat0, lon0 = centre
        return ((np.asarray(lon) - lon0) * 111.19 * np.cos(np.radians(lat0)), (np.asarray(lat) - lat0) * 111.19)

    def catalogue_local(self, centre):
        out = []
        for e in self.read_catalogue():
            x, y = self.to_xy(e["lat"], e["lon"], centre)
            out.append(dict(e, x=float(x), y=float(y), z=float(e["dep"])))
        return out

    # ---------------------------------------------------------------------------------------- waveforms
    def _preprocess(self, raw):
        """demean, detrend, Butterworth 4 zero-phase 0.5-6 Hz, resample_poly 4/25 to 16 Hz, components (E, N, -Z)."""
        from scipy.signal import butter, resample_poly, sosfiltfilt
        raw -= raw.mean(-1, keepdims=True)
        tt = np.arange(raw.shape[-1])
        for i in range(raw.shape[0]):
            for k in range(3):
                p = np.polyfit(tt, raw[i, k], 1)
                raw[i, k] -= np.polyval(p, tt)
        sos = butter(4, [0.5, 6.0], btype="band", fs=100.0, output="sos")
        filt = sosfiltfilt(sos, raw, axis=-1)
        d16 = resample_poly(filt, 4, 25, axis=-1)
        d16[:, 2] *= -1.0
        return raw, d16

    def station_scale(self, d16):
        a, b = int(self.cal[0] / DT), int(self.cal[1] / DT)
        seg = d16[:, :, a:b]
        return np.median(np.abs(seg - np.median(seg, axis=-1, keepdims=True)), axis=(1, 2)) / 0.6745

    @staticmethod
    def window(d16n, t_start):
        """r0_lib.window: 40 s from the nearest sample of t_start, zero-padded past the end of the array."""
        k = int(round(t_start / DT))
        W = np.zeros((d16n.shape[0], 3, NT))
        seg = d16n[:, :, k:k + NT]
        W[:, :, :seg.shape[-1]] = seg
        return W


class Ridgecrest(Dataset):
    def read_catalogue(self):
        import csv
        from obspy import UTCDateTime
        t0 = UTCDateTime(self.t0_utc)
        out = []
        for r in csv.DictReader(open(os.path.join(self.data, "catalog.csv"))):
            out.append({"t": float(UTCDateTime(r["time"]) - t0), "M": float(r["magnitude"]),
                        "lat": float(r["latitude"]), "lon": float(r["longitude"]), "dep": float(r["depth_km"])})
        return out

    def stations(self, names, centre):
        import csv
        rows = {}
        for r in csv.DictReader(open(os.path.join(self.data, "stations.csv"))):
            if r["station"] in names and r["channel"].startswith("HH"):
                rows.setdefault(r["station"], {})[r["channel"][-1]] = r
        STA, sens = [], []
        for n in names:
            r = rows[n]
            x, y = self.to_xy(float(r["Z"]["latitude"]), float(r["Z"]["longitude"]), centre)
            STA.append([float(x), float(y), 0.0])
            sens.append([float(r[c]["response"]) for c in "ENZ"])
        return np.array(STA), np.array(sens)

    def load_waveforms(self, names, sens, t_end):
        from obspy import UTCDateTime, read
        t0 = UTCDateTime(self.t0_utc)
        st = read(os.path.join(self.data, "wave.mseed"))
        n100 = int(round(t_end * 100))
        raw = np.zeros((len(names), 3, n100))
        for i, n in enumerate(names):
            for k, c in enumerate("ENZ"):
                tr = st.select(station=n, channel="HH" + c)
                assert len(tr) == 1, (n, c, len(tr))
                tr = tr[0].copy()
                assert abs(tr.stats.sampling_rate - 100.0) < 1e-6
                i0 = int(round((t0 - tr.stats.starttime) * 100))
                x = tr.data.astype(float)
                seg = np.zeros(n100)
                a, b = max(0, i0), min(len(x), i0 + n100)
                seg[a - i0:b - i0] = x[a:b]
                raw[i, k] = seg / sens[i, k]
        return self._preprocess(raw)


class HiNet(Dataset):
    @staticmethod
    def _secs(hms):
        h, m, s = hms.split(":")
        return (int(h) - 3) * 3600 + int(m) * 60 + float(s) - 1200.0

    def read_jma(self):
        out = []
        for line in open(os.path.join(self.data, "catalog.dat")):
            f = line.split()
            if len(f) >= 7:
                out.append({"t": self._secs(f[1]), "lat": float(f[2]), "lon": float(f[3]), "dep": float(f[4]),
                            "M": float(f[6])})
        return out

    def read_catalogue(self):
        """truth list: Match&Locate detections; 'jma' True for those within 1.0 s and 5 km of a JMA event."""
        jma = [(e["t"], e["lat"], e["lon"]) for e in self.read_jma()]
        out = []
        for line in open(os.path.join(self.data, "MultipleTemplate", "DetectedFinal.dat")):
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            t, lat, lon, dep, mag = self._secs(f[2]), float(f[3]), float(f[4]), float(f[5]), float(f[6])
            isj = any(abs(t - tj) <= 1.0 and np.hypot((lat - a) * 111.19, (lon - b) * 111.19 * np.cos(np.radians(lat)))
                      <= 5.0 for tj, a, b in jma)
            out.append({"t": t, "M": mag, "lat": lat, "lon": lon, "dep": dep, "jma": bool(isj)})
        return out

    def frame(self):
        """centre = mean of the JMA epicentres with origin in CAL (registered)."""
        c = [e for e in self.read_jma() if self.cal[0] <= e["t"] < self.cal[1]]
        return float(np.mean([e["lat"] for e in c])), float(np.mean([e["lon"] for e in c]))

    def _sac(self, name, comp):
        from obspy import read
        return read(os.path.join(self.data, "Trace", "20120902", f"N.{name}.{comp}"))[0]

    def stations(self, names, centre):
        STA, sens = [], []
        for n in names:
            h = self._sac(n, "U").stats.sac
            x, y = self.to_xy(float(h.stla), float(h.stlo), centre)
            STA.append([float(x), float(y), 0.0])
            sens.append([1.0, 1.0, 1.0])
        return np.array(STA), np.array(sens)

    def load_waveforms(self, names, sens, t_end):
        from obspy import UTCDateTime
        t0 = UTCDateTime(self.t0_utc)
        n100 = int(round(t_end * 100))
        raw = np.zeros((len(names), 3, n100))
        for i, n in enumerate(names):
            for k, c in enumerate("ENU"):
                tr = self._sac(n, c)
                assert abs(tr.stats.sampling_rate - 100.0) < 1e-6
                i0 = int(round((t0 - tr.stats.starttime) * 100))
                x = tr.data.astype(float)
                seg = np.zeros(n100)
                a, b = max(0, i0), min(len(x), i0 + n100)
                seg[a - i0:b - i0] = x[a:b]
                raw[i, k] = seg / sens[i, k]
        return self._preprocess(raw)


RIDGECREST = Ridgecrest(
    name="R0", folder="ridgecrest_2019", t0_utc="2019-07-04T19:00:00", cal=(0.0, 900.0), eval=(900.0, 3560.0),
    subnets={"N8": ["CLC", "SRT", "CCC", "TOW2", "LRL", "WRC2", "WBM", "SLA"],
             "N6": ["CLC", "SRT", "CCC", "TOW2", "LRL", "WRC2"],
             "N4": ["CLC", "SRT", "CCC", "TOW2"],
             "U8": ["WBM", "SRT", "TOW2", "JRC2", "WCS2", "CLC", "WRC2", "MPM"]},
    m_target=2.0)

HINET = HiNet(
    name="R1", folder="hinet_2012", t0_utc="2012-09-02T03:20:00", cal=(0.0, 840.0), eval=(840.0, 1960.0),
    subnets={"H7": ["ATKH", "INWH", "NAZH", "ONIH", "THTH", "TSTH", "YNZH"]}, m_target=0.5, m_terms=1.0)
