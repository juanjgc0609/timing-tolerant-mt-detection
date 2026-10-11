"""
r0_lib.py - Task R0 (registered 10 Oct 2026 00:00:50 UTC): real-data pilot on one hour of the 2019 Ridgecrest
aftershock sequence (CI network, SCSN catalogue; benchflow-ai/skillsbench, Apache-2.0).

Import this module BEFORE any project module: it sets the registered homogeneous medium (VP = 6.0 km/s, VS = VP/1.73)
and the continuous origin range (0.0-14.9 s, step 0.1) in task0_common, so that every module that imports these
constants (xc, station_corrections, c10_lib, c11_lib, ...) builds its tables with them.

Time origin: 2019-07-04T19:00:00 UTC (t = 0).  CAL = origins in [0, 900) s; EVAL = origins in [900, 3560] s.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
EXPL = os.path.join(WORK, "explore_contrib")
for p in (EXPL, WORK):
    if p not in sys.path:
        sys.path.insert(0, p)
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import task0_common as T0  # noqa: E402

assert "xc" not in sys.modules and "c10_lib" not in sys.modules, "import r0_lib first"
VP_R0 = 6.0
T0.VP = VP_R0
T0.VS = VP_R0 / 1.73
T0.DEFAULT_TAUS = np.round(np.arange(2.0, 11.95, 0.1), 10)       # registered budget reduction (R0_DEV_LOG)
TAUS = T0.DEFAULT_TAUS
DT, NT = T0.DT, T0.NT                     # 1/16 s, 640 samples (40 s)

DATA = "/home/claude/realdata/skillsbench_eqassoc"
T_CAL = (0.0, 900.0)
T_EVAL = (900.0, 3560.0)
OWN = (5.0, 8.0)                          # origins owned by a window (s, window time)
STEP = 3.0                                # s
KMAX = 6                                  # registered budget reduction (R0_DEV_LOG)
DEPTHS = (3.0, 7.0, 11.0)
M_TARGET = 2.0
SUBNETS = {
    "N8": ["CLC", "SRT", "CCC", "TOW2", "LRL", "WRC2", "WBM", "SLA"],
    "N6": ["CLC", "SRT", "CCC", "TOW2", "LRL", "WRC2"],
    "N4": ["CLC", "SRT", "CCC", "TOW2"],
    "U8": ["WBM", "SRT", "TOW2", "JRC2", "WCS2", "CLC", "WRC2", "MPM"],
}
T0_UTC = "2019-07-04T19:00:00"
LAT0_DEG = None                           # frame centre: set by frame() from the CAL catalogue


# ============================================================================================ catalogue and frame
def read_catalogue():
    import csv
    from obspy import UTCDateTime
    t0 = UTCDateTime(T0_UTC)
    out = []
    for r in csv.DictReader(open(os.path.join(DATA, "catalog.csv"))):
        out.append({"t": float(UTCDateTime(r["time"]) - t0), "M": float(r["magnitude"]),
                    "lat": float(r["latitude"]), "lon": float(r["longitude"]), "dep": float(r["depth_km"])})
    return out


def frame(cat=None):
    """centre (lat, lon) = mean of the catalogued epicentres with origin in CAL (all magnitudes)."""
    cat = cat if cat is not None else read_catalogue()
    c = [e for e in cat if T_CAL[0] <= e["t"] < T_CAL[1]]
    return float(np.mean([e["lat"] for e in c])), float(np.mean([e["lon"] for e in c]))


def to_xy(lat, lon, centre):
    lat0, lon0 = centre
    return ((np.asarray(lon) - lon0) * 111.19 * np.cos(np.radians(lat0)), (np.asarray(lat) - lat0) * 111.19)


def catalogue_local(centre):
    out = []
    for e in read_catalogue():
        x, y = to_xy(e["lat"], e["lon"], centre)
        out.append(dict(e, x=float(x), y=float(y), z=float(e["dep"])))
    return out


def stations(names, centre):
    import csv
    rows = {}
    for r in csv.DictReader(open(os.path.join(DATA, "stations.csv"))):
        if r["station"] in names and r["channel"].startswith("HH"):
            rows.setdefault(r["station"], {})[r["channel"][-1]] = r
    STA, sens = [], []
    for n in names:
        r = rows[n]
        x, y = to_xy(float(r["Z"]["latitude"]), float(r["Z"]["longitude"]), centre)
        STA.append([float(x), float(y), 0.0])
        sens.append([float(r[c]["response"]) for c in "ENZ"])
    return np.array(STA), np.array(sens)


# ============================================================================================ waveforms
def load_waveforms(names, sens, t_end):
    """HH E, N, Z of the stations, counts -> m/s, demean, detrend; 100 Hz arrays (ns, 3, n100) from t = 0 to t_end and
    the detector version: Butterworth 4 zero-phase 0.5-6 Hz, resample_poly 4/25 to 16 Hz, components (E, N, -Z)."""
    from obspy import UTCDateTime, read
    from scipy.signal import butter, resample_poly, sosfiltfilt
    t0 = UTCDateTime(T0_UTC)
    st = read(os.path.join(DATA, "wave.mseed"))
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
    raw -= raw.mean(-1, keepdims=True)
    tt = np.arange(n100)
    for i in range(raw.shape[0]):
        for k in range(3):
            p = np.polyfit(tt, raw[i, k], 1)
            raw[i, k] -= np.polyval(p, tt)
    sos = butter(4, [0.5, 6.0], btype="band", fs=100.0, output="sos")
    filt = sosfiltfilt(sos, raw, axis=-1)
    d16 = resample_poly(filt, 4, 25, axis=-1)
    d16[:, 2] *= -1.0
    return raw, d16


def station_scale(d16):
    """per station: robust sd (MAD / 0.6745 over the three components) of the 16 Hz detector traces in the CAL
    segment; every station is divided by it (registered matched filter normalises by one sd per window)."""
    a, b = int(T_CAL[0] / DT), int(T_CAL[1] / DT)
    seg = d16[:, :, a:b]
    return np.median(np.abs(seg - np.median(seg, axis=-1, keepdims=True)), axis=(1, 2)) / 0.6745


def window(d16n, t_start):
    k = int(round(t_start / DT))
    W = np.zeros((d16n.shape[0], 3, NT))
    seg = d16n[:, :, k:k + NT]
    W[:, :, :seg.shape[-1]] = seg
    return W
