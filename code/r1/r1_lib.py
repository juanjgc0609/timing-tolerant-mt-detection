"""
r1_lib.py - Task R1 (registered 17 Oct 2026 18:49:11 UTC): data layer for the Hi-net swarm of Geolab-USTC/MatchLocate2
(7 stations, 2012-09-02 03:20:00-03:53:20 UTC).  Same interface as r0/r0_lib.py, so that the frozen R0 modules
(r0_run.Net, r0_picks, r0_analyze) run unchanged on these data: r1 scripts install this module as `r0_lib` in
sys.modules before importing them.

Import this module BEFORE any project module: it sets the homogeneous reference medium (VP 6.0 km/s, VS = VP / 1.73; it
only centres the +-2 s peak searches and the descriptive R0-style comparators, every registered method uses JMA2001)
and the origin range in task0_common, exactly as r0_lib.

Time origin: 2012-09-02T03:20:00 UTC (t = 0).  CAL = origins in [0, 840) s; EVAL = origins in [840, 1960) s.
Truth = Match&Locate detections (MultipleTemplate/DetectedFinal.dat); JMA events (catalog.dat) flagged.
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

assert "xc" not in sys.modules and "c10_lib" not in sys.modules, "import r1_lib first"
VP_R0 = 6.0
T0.VP = VP_R0
T0.VS = VP_R0 / 1.73
T0.DEFAULT_TAUS = np.round(np.arange(2.0, 11.95, 0.1), 10)
TAUS = T0.DEFAULT_TAUS
DT, NT = T0.DT, T0.NT

DATA = "/home/claude/realdata/MatchLocate2/example"
T_CAL = (0.0, 840.0)
T_EVAL = (840.0, 1960.0)
OWN = (5.0, 8.0)
STEP = 3.0
KMAX = 6
DEPTHS = (3.0, 7.0, 11.0)
M_TARGET = 0.5
M_TERMS = 1.0                             # truth events used for the station terms (CAL, magnitude >= 1.0)
SUBNETS = {"H7": ["ATKH", "INWH", "NAZH", "ONIH", "THTH", "TSTH", "YNZH"]}
T0_UTC = "2012-09-02T03:20:00"


def _secs(hms):
    h, m, s = hms.split(":")
    return (int(h) - 3) * 3600 + int(m) * 60 + float(s) - 1200.0


def read_catalogue():
    """truth list: Match&Locate detections; 'jma' True for the detections that coincide with a JMA event
    (within 1.0 s and 5 km)."""
    jma = []
    for line in open(os.path.join(DATA, "catalog.dat")):
        f = line.split()
        if len(f) >= 7:
            jma.append((_secs(f[1]), float(f[2]), float(f[3])))
    out = []
    for line in open(os.path.join(DATA, "MultipleTemplate", "DetectedFinal.dat")):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        t, lat, lon, dep, mag = _secs(f[2]), float(f[3]), float(f[4]), float(f[5]), float(f[6])
        isj = any(abs(t - tj) <= 1.0 and np.hypot((lat - a) * 111.19, (lon - b) * 111.19 * np.cos(np.radians(lat))) <= 5.0
                  for tj, a, b in jma)
        out.append({"t": t, "M": mag, "lat": lat, "lon": lon, "dep": dep, "jma": bool(isj)})
    return out


def read_jma():
    out = []
    for line in open(os.path.join(DATA, "catalog.dat")):
        f = line.split()
        if len(f) >= 7:
            out.append({"t": _secs(f[1]), "lat": float(f[2]), "lon": float(f[3]), "dep": float(f[4]), "M": float(f[6])})
    return out


def frame(cat=None):
    """centre = mean of the JMA epicentres with origin in CAL (registered)."""
    c = [e for e in read_jma() if T_CAL[0] <= e["t"] < T_CAL[1]]
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


def _sac(name, comp):
    from obspy import read
    return read(os.path.join(DATA, "Trace", "20120902", f"N.{name}.{comp}"))[0]


def stations(names, centre):
    """station coordinates from the SAC headers; 'sensitivity' 1 (counts: no response is available)."""
    STA, sens = [], []
    for n in names:
        h = _sac(n, "U").stats.sac
        x, y = to_xy(float(h.stla), float(h.stlo), centre)
        STA.append([float(x), float(y), 0.0])
        sens.append([1.0, 1.0, 1.0])
    return np.array(STA), np.array(sens)


def load_waveforms(names, sens, t_end):
    """E, N, U counts, demean, detrend; 100 Hz arrays (ns, 3, n100) from t = 0 to t_end, and the detector version:
    Butterworth 4 zero-phase 0.5-6 Hz, resample_poly 4/25 to 16 Hz, components (E, N, -U).  As r0_lib."""
    from obspy import UTCDateTime
    from scipy.signal import butter, resample_poly, sosfiltfilt
    t0 = UTCDateTime(T0_UTC)
    n100 = int(round(t_end * 100))
    raw = np.zeros((len(names), 3, n100))
    for i, n in enumerate(names):
        for k, c in enumerate("ENU"):
            tr = _sac(n, c)
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
    a, b = int(T_CAL[0] / DT), int(T_CAL[1] / DT)
    seg = d16[:, :, a:b]
    return np.median(np.abs(seg - np.median(seg, axis=-1, keepdims=True)), axis=(1, 2)) / 0.6745


def window(d16n, t_start):
    k = int(round(t_start / DT))
    W = np.zeros((d16n.shape[0], 3, NT))
    seg = d16n[:, :, k:k + NT]
    W[:, :, :seg.shape[-1]] = seg
    return W


def install():
    """make this module importable as r0_lib, so that the frozen R0 modules use it."""
    sys.modules["r0_lib"] = sys.modules[__name__]
    r0dir = os.path.join(WORK, "r0")
    if r0dir not in sys.path:
        sys.path.insert(0, r0dir)


def setup():
    """Task R1 configuration of the shared modules: r1_lib as r0_lib, JMA2001 for the layered engine, the hybrid and
    the 1-D associators (own PyOcto model file), the R1 epoch for the pick conversions."""
    install()
    for p in (os.path.join(WORK, "hyb"), os.path.join(WORK, "hyb", "pydeps_pyrocko"), os.path.join(WORK, "gamma_src"),
              os.path.join(WORK, "r0_verif", "pydeps")):
        if p not in sys.path:
            sys.path.append(p)
    from datetime import datetime, timezone
    import c18_tt as TT
    import c18j_layered as LJ
    import hyb_comp as HC
    import r0_picks as RP
    LJ.MODEL = TT.JMA2001
    HC.MODEL_PATH = os.path.join(HERE, "pyocto_jma2001.model")
    HC.POINTS = TT.JMA2001_POINTS
    RP.EPOCH = datetime(2012, 9, 2, 3, 20, 0, tzinfo=timezone.utc)


def t_end(segment):
    """waveforms read up to the end of the segment plus 40 s; CAL never reads EVAL samples."""
    return (T_CAL[1] if segment == "cal" else T_EVAL[1]) + 40.0


def make_net(segment):
    """r0_run.Net on the 7 stations; its calibration windows are the CAL truth events with magnitude >= M_TERMS
    (registered for both kinds of station terms)."""
    global M_TARGET
    import r0_run as RR
    keep = M_TARGET
    M_TARGET = M_TERMS
    try:
        net = RR.Net("H7", t_end(segment))
    finally:
        M_TARGET = keep
    return net
