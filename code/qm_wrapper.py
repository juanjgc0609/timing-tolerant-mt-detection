"""
qm_wrapper.py - QuakeMigrate (Winder et al., Seismica, https://seismica.library.mcgill.ca/article/view/1854;
github.com/QuakeMigrate/QuakeMigrate; PyPI `quakemigrate` 1.2.2) run on the synthetic three-component
windows of task1y, as a published comparator scored on the same exact FROC as the grid methods.

Public API
    run_quakemigrate(STA, D, tmpdir, vp, vs, cfg=None) -> [(coalescence_score, x_km, y_km, z_km, t0_s)]
    qm_records(events, truth, min_sep_km=5.0, max_events=5) -> [(score, hit_index or -1)] for froc_exact
    QM_CFG     "as published" channel practice (P onset on Z, S onset on N+E), parameters calibrated on
               seed 999
    QM_CFG_3C  same package and stages, P and S onsets from all three components (Z, N, E), parameters
               calibrated separately on the same seed-999 windows. RECOMMENDED comparator here: see below.

This is the real package (no re-implementation):
    lookup table   quakemigrate.lut.compute_traveltimes(method="homogeneous", vp, vs) on a Grid3D,
                   x, y in [-21, 22] km, z in [0, 20] km, 1 km nodes (44 x 44 x 21)
    onset          quakemigrate.signal.onsets.STALTAOnset, classic STA/LTA of the energy (the package
                   default for detect), per-channel onsets combined as RMS, clipped at 0.4
    detect         quakemigrate.QuakeScan.detect(): migration of the onsets and geometric-mean stacking
                   (the package's C library), maximum over the grid per time sample (.scanmseed)
    trigger        quakemigrate.Trigger.trigger(): static threshold on the maximum-coalescence trace,
                   marginal-window / minimum-event-interval merging, the package's own filtering of events
                   by trigger time window; writes TriggeredEvents.csv, which run_quakemigrate reads
    location       detect-stage location of the maximum coalescence at the triggered peak (COA_X/Y/Z of
                   TriggeredEvents.csv; the detect grid is NOT decimated). locate() is not run: checked on
                   5 windows (locate_check), it moves the epicentre by 0.2-4.5 km and changes no hit/miss.
    score          TRIG_COA = the maximum coalescence (geometric mean of onsets) of the triggered event.
    The calibration sweeps apply the same Trigger methods in memory (trigger_scan); on 192
    window x configuration pairs they gave events identical to the file-based Trigger.trigger().

Channel practice. The generator is a homogeneous full space with straight rays; at epicentral distance /
depth 2-6 the rays arrive 60-80 deg from vertical, so P energy is mostly on the horizontals and SV mostly on
the vertical. With the default practice the vertical-channel "P" onset peaks at the S arrival and the
horizontal "S" onset at the P arrival (checked on a clutter-free positive control): QM_CFG finds 3/8 events
even at 5x the nominal SNR, QM_CFG_3C finds 8/8 (median epicentral error 0.5 km) and 7/8 at nominal SNR.
The default practice is right for real Earth (velocity gradients steepen incidence), not for this generator.

Frames. make_network / generate_y use a local Cartesian frame (x, y, z) in km with stations at z=0 and
sources at z=+5..+13; the P polarisation (STA-src)/d points UP at the station for a source below it, so z is
depth (positive down). Written to MiniSEED: HHE = x (East), HHN = y (North), HHZ = -z (Up); the energy
STA/LTA is sign-blind, only the channel-to-axis assignment matters. The frame is placed at
(lon, lat) = (-76.5, 3.4) with the simple equirectangular conversion
    lon = lon0 + deg(x / (R cos lat0)),  lat = lat0 + deg(y / R),  R = 6371 km (sphere)
and the grid projection is pyproj 'eqc' (lat_ts = lat0, same sphere, units km), which is exactly that
conversion: grid-projection coordinates equal local km to 1e-12 km. Station elevations 0.

Time. Window sample k is at T_REF + k/16 s. QuakeMigrate needs data before and after the scanned origin
times (LTA + 3 STA + taper pre-pad, 11-12 s; max travel time (~24 s) + 2 LTA + taper post-pad, 40-41 s),
which a 40-s window cannot provide with LTA 6 s. The window is embedded in pad_s = 40 s of Gaussian noise on each side, with the
per-channel standard deviation estimated robustly (MAD) from the window itself (no signal information
added). Detect scans origin times scan_window = 3.0-11.5 s; the trigger keeps events whose coalescence peak
lies in trigger_window = (5.0, 8.0) + onset lag (classic STA/LTA peaks after the arrival; lag measured on the
calibration windows), i.e. a 3-s origin-time window like DEFAULT_TAUS (5.0-8.0 s) of the grid methods.
t0_s is the coalescence peak time and is therefore late by the onset lag (median 0.47 s on QM_CFG_3C hits).

Sampling. The 16 Hz synthetics are written at 16 Hz; QuakeMigrate resamples them itself to fs = 50 Hz
(Archive(resample=True, upfactor=25): linear interpolation to 400 Hz, then its decimation by 8).

Calibration (seed-999 stream of task1y_clutter.run_cell, thin clutter 0.35, pert, cells random-6 =
ladder_6 and onesided-8 = unilateral_c40; windows 30..59 of each cell, 47 sources; windows 0..29 are the
test of qm_report.txt): detect grid CAL_BANDPASS x CAL_STALTA (24 onsets), trigger grid floor quantile x
(marginal window, minimum event interval) x normalised coalescence x smoothing, two phases (reduced trigger
grid for all 24, full grid for the best 3); criterion recall at <= 0.25 realized false events per window
(froc_exact.exact_HF + task1r_fair.recall_at on qm_records), ties broken at 0.1 then 0.05.

    python qm_wrapper.py calibrate --out DIR --workers 2   # both variants (about 30 min on 2 cores)
    python qm_wrapper.py test --out DIR                      # the 30-window test, both variants + MT-FW
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MPLBACKEND", "Agg")
import contextlib
import json
import logging
import shutil
import sys
import tempfile
import time
import zlib

import numpy as np
import pandas as pd
from obspy import Stream, Trace, UTCDateTime
from pyproj import Proj, Transformer

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from froc_exact import exact_HF  # noqa: E402
from task1o_consistent import TOL_KM  # noqa: E402

LON0, LAT0 = -76.5, 3.4
R_EARTH_M = 6371000.0
T_REF = UTCDateTime("2024-01-01T12:00:00")
FS_NATIVE = 16                       # 1 / task0_common.DT
NET = "XX"

# ---------------------------------------------------------------------------------------------
# Configuration ("as published" channel practice). The detection-relevant values were chosen on the
# seed-999 calibration windows maximising recall at <= 0.25 realized false events per window
# (calibration recall 0.09 / 0.09 / 0.11 at 0.05 / 0.1 / 0.25; the package-default onset, 2-16 Hz and
# STA/LTA 0.2/1.0 s, gave 0.00). `threshold` is the static trigger FLOOR; the FROC is obtained by
# sweeping the coalescence score of the triggered events (qm_records + froc_exact).
# ---------------------------------------------------------------------------------------------
QM_CFG = {
    # lookup table (homogeneous, straight rays)
    "grid_x": (-21.0, 22.0), "grid_y": (-21.0, 22.0), "grid_z": (0.0, 20.0), "node_spacing": 1.0,
    # sampling
    "fs": 50, "upfactor": 25,                       # 16 Hz x 25 = 400 Hz -> decimate x8 -> 50 Hz
    # onset (STALTAOnset) - calibrated on seed 999 (detect grid CAL_BANDPASS x CAL_STALTA)
    "bp_p": (1.0, 4.0, 2), "bp_s": (1.0, 4.0, 2),  # [lowcut Hz, highcut Hz, corners] (two-pass)
    "stalta_p": (0.4, 6.0), "stalta_s": (0.4, 6.0),  # [STA s, LTA s]
    "position": "classic", "signal_transform": "energy", "min_onset_value": 0.4,
    "channel_maps": {"P": "*Z", "S": "*[N,E]"}, "channel_counts": {"P": 1, "S": 2},
    # detect / trigger - calibrated on seed 999 (trigger grid CAL_FLOOR_Q x CAL_MW_MEI x norm x smooth)
    "scan_window": (3.0, 11.5),                     # origin times computed by detect (s)
    "trigger_window": (5.14, 8.14),                 # TAU_PRIOR (5.0, 8.0) + onset lag 0.14 s
    "threshold": 1.641705,                          # static floor = median of the pooled COA trace
    "normalise_coalescence": False, "smooth_coa": True, "smooth_sigma": 0.2,
    "marginal_window": 0.5, "min_event_interval": 1.0,
    # data padding
    "pad_s": 40.0, "pad_seed": 12345,
    "threads": 1,
}


# QM-3C: the same package and stages, with the P and S onsets computed from ALL three components
# (channel_maps P "*[Z,N,E]", S "*[Z,N,E]", RMS combination of the per-channel STA/LTA, a documented
# STALTAOnset option). Needed here because the generator is a homogeneous full space with straight
# rays: at epicentral distance / depth 2-6 the rays arrive 60-80 deg from vertical, so P energy is
# mostly on the horizontals and SV mostly on the vertical, the opposite of the near-vertical incidence
# that motivates "P on Z, S on N/E" in real (velocity-gradient) Earth. Calibrated separately on the
# same seed-999 windows with the same grids (calibration recall 0.23 / 0.38 / 0.38 at 0.05 / 0.1 / 0.25;
# package-default onset 0.00 / 0.04 / 0.06).
CHANNELS_3C = {"channel_maps": {"P": "*[Z,N,E]", "S": "*[Z,N,E]"}, "channel_counts": {"P": 3, "S": 3}}
QM_CFG_3C = {**QM_CFG, **CHANNELS_3C,
             "bp_p": (1.0, 6.0, 2), "bp_s": (1.0, 6.0, 2), "stalta_p": (0.8, 6.0), "stalta_s": (0.8, 6.0),
             "trigger_window": (5.22, 8.22),        # TAU_PRIOR + onset lag 0.22 s
             "threshold": 1.4404625,                # 0.75 quantile of the pooled COA trace
             "normalise_coalescence": False, "smooth_coa": False,
             "marginal_window": 0.25, "min_event_interval": 0.5}

# ------------------------------------------------------------------------------------ frames
def _projections():
    coord_proj = Proj(proj="longlat", R=R_EARTH_M, no_defs=True)
    grid_proj = Proj(proj="eqc", lat_ts=LAT0, lat_0=LAT0, lon_0=LON0, R=R_EARTH_M, units="km",
                     no_defs=True)
    return coord_proj, grid_proj


_CP, _GP = _projections()
_KM2LL = Transformer.from_proj(_GP, _CP)
_LL2KM = Transformer.from_proj(_CP, _GP)


def km_to_lonlat(x, y):
    lon, lat = _KM2LL.transform(np.asarray(x, float), np.asarray(y, float))
    return lon, lat


def lonlat_to_km(lon, lat):
    x, y = _LL2KM.transform(np.asarray(lon, float), np.asarray(lat, float))
    return x, y


def station_names(n):
    return [f"S{i:02d}" for i in range(n)]


def stations_df(STA):
    lon, lat = km_to_lonlat(STA[:, 0], STA[:, 1])
    # QuakeMigrate convention (io.read_stations): Elevation stored as -elevation (positive down),
    # in the units of the grid projection (km). Stations are at elevation 0.
    return pd.DataFrame({"Latitude": np.atleast_1d(lat), "Longitude": np.atleast_1d(lon),
                         "Elevation": -np.asarray(STA[:, 2], float) * 0.0,
                         "Name": station_names(len(STA))})


# ------------------------------------------------------------------------------ lookup table
def build_lut(STA, vp, vs, cfg=None):
    from quakemigrate.lut import compute_traveltimes
    cfg = {**QM_CFG, **(cfg or {})}
    (x0, x1), (y0, y1), (z0, z1) = cfg["grid_x"], cfg["grid_y"], cfg["grid_z"]
    eps = 1e-6                     # km; ur corner nudged inwards so that ceil(43 + 1e-13) != 44
    lon0, lat0 = km_to_lonlat(x0, y0)
    lon1, lat1 = km_to_lonlat(x1 - eps, y1 - eps)
    grid_spec = {"ll_corner": [float(lon0), float(lat0), z0],
                 "ur_corner": [float(lon1), float(lat1), z1 - eps],
                 "node_spacing": [cfg["node_spacing"]] * 3,
                 "grid_proj": _GP, "coord_proj": _CP}
    with _quiet():
        lut = compute_traveltimes(grid_spec, stations_df(STA), method="homogeneous",
                                  phases=["P", "S"], vp=float(vp), vs=float(vs), log=False)
    want = np.round(np.array([x1 - x0, y1 - y0, z1 - z0]) / cfg["node_spacing"]).astype(int) + 1
    assert np.all(lut.node_count == want), (lut.node_count, want)
    gx = lut.grid_xyz
    assert abs(gx[0].min() - x0) < 1e-6 and abs(gx[1].min() - y0) < 1e-6 and abs(gx[2].min() - z0) < 1e-6
    return lut


# ----------------------------------------------------------------------------------- archive
def write_archive(D, STA, archive_dir, cfg=None):
    """Writes one MiniSEED file per station (HHE, HHN, HHZ at 16 Hz) in the QuakeMigrate
    'YEAR/JD/STATION' archive layout, with the window embedded in noise-only padding."""
    cfg = {**QM_CFG, **(cfg or {})}
    ns, nc, nt = D.shape
    npad = int(round(cfg["pad_s"] * FS_NATIVE))
    rng = np.random.default_rng([int(cfg["pad_seed"]), zlib.crc32(np.ascontiguousarray(D).tobytes())])
    t0 = T_REF - npad / FS_NATIVE
    day = os.path.join(archive_dir, f"{t0.year}", f"{t0.julday:03d}")
    os.makedirs(day, exist_ok=True)
    comp = (("E", D[:, 0]), ("N", D[:, 1]), ("Z", -D[:, 2]))    # x=E, y=N, Up=-z
    for i, name in enumerate(station_names(ns)):
        st = Stream()
        for c, arr in comp:
            x = np.asarray(arr[i], float)
            sd = float(np.median(np.abs(x - np.median(x))) / 0.6745)
            data = np.concatenate([rng.normal(0.0, sd, npad), x, rng.normal(0.0, sd, npad)])
            st += Trace(data=data.astype(np.float64),
                        header={"network": NET, "station": name, "location": "", "channel": f"HH{c}",
                                "sampling_rate": float(FS_NATIVE), "starttime": t0})
        st.write(os.path.join(day, f"{name}.mseed"), format="MSEED", encoding="FLOAT64")


@contextlib.contextmanager
def _quiet():
    """QuakeMigrate logs every step at INFO to stdout; silence INFO for the duration."""
    prev = logging.root.manager.disable
    logging.disable(logging.INFO)
    try:
        yield
    finally:
        logging.disable(prev)


def _onset(cfg):
    from quakemigrate.signal.onsets import STALTAOnset
    return STALTAOnset(position=cfg["position"], sampling_rate=int(cfg["fs"]),
                       signal_transform=cfg["signal_transform"],
                       min_onset_value=cfg["min_onset_value"],
                       phases=["P", "S"],
                       bandpass_filters={"P": list(cfg["bp_p"]), "S": list(cfg["bp_s"])},
                       sta_lta_windows={"P": list(cfg["stalta_p"]), "S": list(cfg["stalta_s"])},
                       channel_maps=dict(cfg["channel_maps"]),
                       channel_counts=dict(cfg["channel_counts"]))


def prepare(STA, D, workdir, vp, vs, cfg=None):
    """Writes the MiniSEED archive of one window under workdir and builds its lookup table."""
    cfg = {**QM_CFG, **(cfg or {})}
    os.makedirs(workdir, exist_ok=True)
    arch = os.path.join(workdir, "archive")
    if os.path.isdir(arch):
        shutil.rmtree(arch)
    write_archive(D, STA, arch, cfg)
    return build_lut(STA, vp, vs, cfg)


def detect(STA, workdir, lut, cfg=None):
    """QuakeScan.detect() over cfg['scan_window'] (a single timestep) on the archive in workdir."""
    from quakemigrate import QuakeScan
    from quakemigrate.io import Archive
    cfg = {**QM_CFG, **(cfg or {})}
    sdir = os.path.join(workdir, "qm", "detect")
    if os.path.isdir(sdir):
        shutil.rmtree(sdir)
    archive = Archive(archive_path=os.path.join(workdir, "archive"), stations=stations_df(STA),
                      archive_format="YEAR/JD/STATION", resample=(int(cfg["fs"]) != FS_NATIVE),
                      upfactor=int(cfg["upfactor"]))
    s0, s1 = cfg["scan_window"]
    fs = int(cfg["fs"])
    with _quiet():
        scan = QuakeScan(archive, lut, onset=_onset(cfg), run_path=workdir, run_name="qm",
                         log=False, threads=int(cfg["threads"]), timestep=float(s1 - s0))
        scan.detect(str(T_REF + s0), str(T_REF + s1 - 1.0 / fs))


def read_scan(workdir, lut, cfg=None, raw=False):
    """The continuous detect output (.scanmseed), as QuakeMigrate's own reader returns it (raw=True:
    columns DT, COA, COA_N, X=lon, Y=lat, Z) or with local-km coordinates and t in s."""
    from quakemigrate.io import Run, read_scanmseed
    cfg = {**QM_CFG, **(cfg or {})}
    s0, s1 = cfg["scan_window"]
    with _quiet():
        data, _st = read_scanmseed(Run(workdir, "qm", "", "detect"), T_REF + s0, T_REF + s1, 0.0,
                                   lut.unit_conversion_factor)
    if raw:
        return data
    x, y = lonlat_to_km(data["X"].to_numpy(float), data["Y"].to_numpy(float))
    return pd.DataFrame({"t": [float(t - T_REF) for t in data["DT"]], "COA": data["COA"].to_numpy(float),
                         "COA_N": data["COA_N"].to_numpy(float), "x": x, "y": y,
                         "z": data["Z"].to_numpy(float)})


def _trigger_kwargs(cfg):
    tw0, tw1 = cfg["trigger_window"]
    return dict(threshold_method="static", static_threshold=float(cfg["threshold"]),
                marginal_window=float(cfg["marginal_window"]),
                min_event_interval=float(cfg["min_event_interval"]),
                normalise_coalescence=bool(cfg["normalise_coalescence"]),
                smooth_coa=bool(cfg["smooth_coa"]), smoothing_kernel_sigma=float(cfg["smooth_sigma"]),
                plot_trigger_summary=False,
                pad=float(max(abs(tw0 - cfg["scan_window"][0]), abs(cfg["scan_window"][1] - tw1)) + 1.0))


def _events_from_df(ev):
    out = []
    if ev is None or len(ev) == 0:
        return out
    x, y = lonlat_to_km(ev["COA_X"].to_numpy(float), ev["COA_Y"].to_numpy(float))
    for k in range(len(ev)):
        out.append((float(ev["TRIG_COA"].iloc[k]), float(x[k]), float(y[k]),
                    float(ev["COA_Z"].iloc[k]), float(UTCDateTime(ev["CoaTime"].iloc[k]) - T_REF)))
    out.sort(key=lambda e: -e[0])
    return out


def run_trigger(workdir, lut, cfg=None, name="trig"):
    """QuakeMigrate Trigger.trigger() on the detect output (public API, reads the .scanmseed and
    writes TriggeredEvents.csv); returns [(score, x_km, y_km, z_km, t0_s)]."""
    from quakemigrate import Trigger
    cfg = {**QM_CFG, **(cfg or {})}
    tw0, tw1 = cfg["trigger_window"]
    evdir = os.path.join(workdir, "qm", "trigger", name)
    if os.path.isdir(evdir):
        shutil.rmtree(evdir)
    with _quiet():
        trig = Trigger(lut, run_path=workdir, run_name="qm", trigger_name=name, log=False,
                       **_trigger_kwargs(cfg))
        trig.trigger(str(T_REF + tw0), str(T_REF + tw1))
    evd = os.path.join(evdir, "events")
    if not os.path.isdir(evd):
        return []
    frames = [pd.read_csv(os.path.join(evd, f)) for f in sorted(os.listdir(evd))
              if f.endswith("_TriggeredEvents.csv")]
    frames = [f for f in frames if not f.empty]
    return _events_from_df(pd.concat(frames, ignore_index=True) if frames else None)


def trigger_scan(raw_scan, cfg, lut=None, _cache={}):
    """The same QuakeMigrate Trigger steps as Trigger._trigger_batch (smoothing, static threshold,
    _identify_candidates, _refine_candidates, _filter_events), applied in memory to a scan read by
    read_scan(raw=True), without plotting or file output. Used for the calibration sweeps; checked
    against run_trigger (file-based public API) in check_trigger_equivalence()."""
    from quakemigrate import Trigger
    kw = _trigger_kwargs(cfg)
    key = tuple(sorted((k, v) for k, v in kw.items()))
    if key not in _cache:
        tmp = tempfile.mkdtemp(prefix="qmtrig_")
        with _quiet():
            _cache[key] = Trigger(lut, run_path=tmp, run_name="cal", log=False, **kw)
        shutil.rmtree(tmp, ignore_errors=True)
    trig = _cache[key]
    tw0, tw1 = cfg["trigger_window"]
    data = raw_scan.copy()
    fs = 1.0 / float(data["DT"].iloc[1] - data["DT"].iloc[0])
    with _quiet():
        if trig.smooth_coa:
            data = trig._smooth_coa(data, fs)
        on = "COA_N" if trig.normalise_coalescence else "COA"
        thr = trig._get_threshold(data[on], fs)
        cand = trig._identify_candidates(data, on, thr)
        if cand.empty:
            return []
        ev = trig._filter_events(trig._refine_candidates(cand), T_REF + tw0, T_REF + tw1, None)
    return _events_from_df(ev)


def run_quakemigrate(STA, D, tmpdir, vp, vs, cfg=None):
    """QuakeMigrate detect + trigger on one window.

    STA (ns, 3) station coordinates in km (local frame, z=0); D (ns, 3, NT) window at 16 Hz with
    components (x, y, z) = (E, N, depth); tmpdir a scratch directory (a unique subdirectory is
    created and removed); vp, vs in km/s; cfg overrides QM_CFG.
    Returns [(coalescence_score, x_km, y_km, z_km, t0_s)] sorted by decreasing score, one entry per
    triggered event (score = the maximum coalescence value of the event, TRIG_COA; t0_s = time of
    the coalescence peak = origin time, s after the window start)."""
    cfg = {**QM_CFG, **(cfg or {})}
    os.makedirs(tmpdir, exist_ok=True)
    work = tempfile.mkdtemp(prefix="qm_", dir=tmpdir)
    try:
        lut = prepare(STA, D, work, vp, vs, cfg)
        detect(STA, work, lut, cfg)
        return run_trigger(work, lut, cfg)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def qm_records(events, truth, min_sep_km=5.0, max_events=5):
    """(score, hit) records for froc_exact.exact_HF, with the hit rule of score_events: the nearest
    true source within TOL_KM (5 km) epicentral distance is hit (its index), otherwise -1 (false).
    Events are taken by decreasing score with a greedy min_sep_km epicentral separation and at most
    max_events per window, as greedy_record does for the grid methods. Two accepted events hitting the
    same source are counted once by exact_HF (neither a second hit nor a false event)."""
    acc, rec = [], []
    for ev in sorted(events, key=lambda e: -e[0]):
        s, p = float(ev[0]), np.asarray(ev[1:3], float)
        if not np.isfinite(s):
            continue
        if any(np.linalg.norm(p - a) <= min_sep_km for a in acc):
            continue
        acc.append(p)
        hit = -1
        if truth:
            dd = [float(np.linalg.norm(p - np.asarray(t_, float)[:2])) for t_ in truth]
            j = int(np.argmin(dd))
            if dd[j] <= TOL_KM:
                hit = j
        rec.append((s, hit))
        if len(acc) >= max_events:
            break
    return rec


# ============================================================================================
# Calibration and test (seed-999 stream of task1y_clutter.run_cell; nothing here is used by
# run_quakemigrate / qm_records)
# ============================================================================================
CELLS = {"random-6": "ladder_6", "onesided-8": "unilateral_c40"}     # task1r_fair.SCEN keys
CAL_SEED, EVBASE, MOFF = 999, 25000, 233                            # run_cell defaults for the stream
N_TEST = 30                                                          # windows 0..29: test
CAL_OFFSET, N_CAL = 30, 30                                           # windows 30..59: calibration

# detect-stage grid (each needs its own detect run)
CAL_BANDPASS = [(1.0, 4.0), (1.0, 6.0), (2.0, 16.0)]                 # (2, 16): package default
CAL_STALTA = [(0.2, 1.0), (0.2, 3.0), (0.4, 1.5), (0.4, 3.0), (0.4, 6.0), (0.8, 3.0), (0.8, 6.0),
              (1.2, 6.0)]                                            # (0.2, 1.0): package default
# trigger-stage grid (in memory on the saved scans)
CAL_MW_MEI = [(0.25, 0.5), (0.5, 1.0), (1.0, 2.0), (0.5, 2.0), (2.0, 4.0)]   # (2, 4): package default
CAL_FLOOR_Q = (0.5, 0.75, 0.9, 0.95)          # static threshold = quantile of the pooled trigger trace
CAL_SCAN = (3.0, 11.5)
TAU_PRIOR = (5.0, 8.0)                        # origin-time prior of the grid methods (DEFAULT_TAUS)


def y_stream(scen, n, seed=CAL_SEED, clut="thin", evbase=EVBASE, moff=MOFF):
    """Windows 0..n-1 of task1y_clutter.run_cell's data stream for one cell: make_network(kind, nst,
    seed+evbase+k), draw_sources, generate_y (pert condition: statics, polarisation noise, site
    factors) with clutter amplitude CLUT_AMP[clut]. Returns [(STA, D, truth, taus)]."""
    from task0_common import make_network
    from task1o_consistent import draw_sources
    from task1r_fair import SCEN
    from task1y_clutter import CLUT_AMP, generate_y
    kind, nst, nclut = SCEN[scen]
    rng = np.random.default_rng(seed + 41 + evbase)
    rng_m = np.random.default_rng(seed + moff)
    rng_h = np.random.default_rng(seed + 61 + evbase)
    out = []
    for k in range(n):
        S = make_network(kind, nst, seed=seed + evbase + k)
        srcs, taus = draw_sources(rng)
        D, _p, tr, _amp = generate_y(S, srcs, taus, nclut, rng, rng_m, rng_h, CLUT_AMP[clut])
        out.append((S, D, tr, list(taus)))
    return out


def _detect_cfgs(extra=None):
    out = []
    for bp in CAL_BANDPASS:
        for w in CAL_STALTA:
            out.append({"bp_p": (bp[0], bp[1], 2), "bp_s": (bp[0], bp[1], 2), "stalta_p": w, "stalta_s": w,
                        "scan_window": CAL_SCAN, **(extra or {})})
    return out


def _scan_arrays(raw):
    return {"t": np.array([float(t - T_REF) for t in raw["DT"]]), "COA": raw["COA"].to_numpy(float),
            "COA_N": raw["COA_N"].to_numpy(float), "X": raw["X"].to_numpy(float),
            "Y": raw["Y"].to_numpy(float), "Z": raw["Z"].to_numpy(float)}


def _raw_from_arrays(a):
    return pd.DataFrame({"DT": [T_REF + float(t) for t in a["t"]], "COA": a["COA"], "COA_N": a["COA_N"],
                         "X": a["X"], "Y": a["Y"], "Z": a["Z"]})


def _cal_worker(args):
    scen, k, S, D, cfgs, root = args
    from task0_common import VP, VS
    work = os.path.join(root, f"{scen}_{k}")
    t0 = time.time()
    lut = prepare(S, D, work, VP, VS)
    scans = []
    for c in cfgs:
        detect(S, work, lut, c)
        scans.append(_scan_arrays(read_scan(work, lut, c, raw=True)))
    shutil.rmtree(work, ignore_errors=True)
    return scen, k, scans, (time.time() - t0) / len(cfgs)


def calibrate_detect(out_pkl, root, workers=2, n_cal=N_CAL, offset=CAL_OFFSET, extra=None):
    """Stage 1: detect for every detect-stage configuration on the calibration windows; the scans
    (maximum coalescence trace and its location, per time sample) are pickled for stage 2. `extra`
    fixes other detect settings for the whole grid (e.g. CHANNELS_3C)."""
    import pickle
    from concurrent.futures import ProcessPoolExecutor
    cfgs = _detect_cfgs(extra)
    jobs, meta = [], {}
    for cell, scen in CELLS.items():
        for k, (S, D, tr, taus) in enumerate(y_stream(scen, offset + n_cal)[offset:], start=offset):
            jobs.append((scen, k, S, D, cfgs, root))
            meta[(scen, k)] = (tr, taus)
    res = {}
    secs = []
    with ProcessPoolExecutor(workers) as ex:
        for scen, k, scans, sec in ex.map(_cal_worker, jobs):
            res[(scen, k)] = scans
            secs.append(sec)
            print(f"   detect {scen} {k}: {sec:.2f} s per config", flush=True)
    with open(out_pkl, "wb") as fh:
        pickle.dump({"cfgs": cfgs, "scans": res, "meta": meta, "sec_per_detect": float(np.mean(secs))}, fh)


def froc_summary(records, ntrue, caps=(0.05, 0.1, 0.25)):
    from task1r_fair import recall_at
    _t, H, F = exact_HF(records)
    return {str(c): recall_at(H, F, ntrue, c) for c in caps}


def onset_lag(scans, meta, keys, ci, tol_km=TOL_KM):
    """Median delay of the maximum-coalescence peak after the true origin time (classic STA/LTA onsets
    peak after the arrival), over calibration events whose peak in [tau-0.5, tau+2.5] s is located
    within tol_km of the source."""
    lags = []
    for key in keys:
        a = scans[key][ci]
        x, y = lonlat_to_km(a["X"], a["Y"])
        tr, taus = meta[key]
        for s, ta in zip(tr, taus):
            m = np.where((a["t"] >= ta - 0.5) & (a["t"] <= ta + 2.5))[0]
            j = m[np.argmax(a["COA"][m])]
            if np.hypot(x[j] - s[0], y[j] - s[1]) <= tol_km:
                lags.append(a["t"][j] - ta)
    return (float(np.median(lags)) if len(lags) >= 5 else None), len(lags)


def _trigger_rows(ci, dc, scans, meta, keys, ntrue, lut, norms, floor_q, mw_mei, smooths):
    fs = int(QM_CFG["fs"])
    lag, nl = onset_lag(scans, meta, keys, ci)
    if lag is None:
        lag = float(dc["stalta_p"][0])
    lag = round(lag * fs) / fs
    tw = (TAU_PRIOR[0] + lag, TAU_PRIOR[1] + lag)
    raws = {k: _raw_from_arrays(scans[k][ci]) for k in keys}
    rows = []
    for norm in norms:
        on = "COA_N" if norm else "COA"
        pooled = np.concatenate([scans[k][ci][on][(scans[k][ci]["t"] >= tw[0]) & (scans[k][ci]["t"] <= tw[1])]
                                 for k in keys])
        for q in floor_q:
            thr = float(np.quantile(pooled, q))
            for mw, mei in mw_mei:
                for sm in smooths:
                    c = {**QM_CFG, **dc, "trigger_window": tw, "threshold": thr, "normalise_coalescence": norm,
                         "marginal_window": mw, "min_event_interval": mei, "smooth_coa": sm}
                    recs = [qm_records(trigger_scan(raws[k], c, lut), meta[k][0]) for k in keys]
                    rows.append({"ci": ci, "bp": dc["bp_p"][:2], "stalta": dc["stalta_p"], "lag": lag,
                                 "n_lag": nl, "norm": norm, "floor_q": q, "threshold": thr, "mw": mw,
                                 "mei": mei, "smooth": sm, "recall_at": froc_summary(recs, ntrue),
                                 "n_events": int(sum(len(r) for r in recs))})
    return rows


def _rank(r):
    return (r["recall_at"]["0.25"], r["recall_at"]["0.1"], r["recall_at"]["0.05"])


def calibrate_trigger(pkl, out_json=None, lut=None, n_top=3):
    """Stage 2 (in memory, QuakeMigrate's Trigger methods): for each detect configuration the onset lag
    (trigger window = TAU_PRIOR + lag), then (a) a reduced trigger grid for every detect configuration
    and (b) the full trigger grid (floor quantile x marginal window / minimum event interval x
    normalised coalescence x smoothing) for the n_top detect configurations of (a). Criterion: recall
    at <= 0.25 realized false events per window over the pooled calibration windows (ties: 0.1, 0.05)."""
    import pickle
    with open(pkl, "rb") as fh:
        P = pickle.load(fh)
    cfgs, scans, meta = P["cfgs"], P["scans"], P["meta"]
    keys = sorted(scans)
    ntrue = sum(len(meta[k][0]) for k in keys)
    rows = []
    for ci, dc in enumerate(cfgs):
        rr = _trigger_rows(ci, dc, scans, meta, keys, ntrue, lut, (False, True), (0.5, 0.9),
                           [(0.25, 0.5), (0.5, 1.0)], (False, True))
        rows += [{**r, "phase": "a"} for r in rr]
        b = max(rr, key=_rank)
        print(f"   (a) cfg {ci} bp {dc['bp_p'][:2]} stalta {dc['stalta_p']} lag {b['lag']:.2f} ({b['n_lag']}) :: "
              f"{b['recall_at']}", flush=True)
    best_ci = sorted({r["ci"] for r in rows}, key=lambda ci: max(_rank(r) for r in rows if r["ci"] == ci),
                     reverse=True)[:n_top]
    for ci in best_ci:
        rr = _trigger_rows(ci, cfgs[ci], scans, meta, keys, ntrue, lut, (False, True), CAL_FLOOR_Q,
                           CAL_MW_MEI, (False, True))
        rows += [{**r, "phase": "b"} for r in rr]
        b = max(rr, key=_rank)
        print(f"   (b) cfg {ci} :: {b['recall_at']} norm {b['norm']} q {b['floor_q']} thr {b['threshold']:.4f} "
              f"mw/mei {b['mw']}/{b['mei']} smooth {b['smooth']}", flush=True)
    if out_json:
        with open(out_json, "w") as fh:
            json.dump({"ntrue": ntrue, "nwin": len(keys), "sec_per_detect": P.get("sec_per_detect"),
                       "rows": rows}, fh)
    return rows


def mtfw_record(S, D, truth):
    """MT-FW reference record exactly as task1y_clutter.run_cell calls it (default grid, DEFAULT_TAUS)."""
    import task1x_mtglrt as X
    from task0_common import DEFAULT_TAUS, default_grid
    from task1o_consistent import matched
    grid = default_grid()
    garr = np.asarray(grid)
    Gs = [X.radiation_operator(g, S) for g in garr]
    Qs = [X.basis(G) for G in Gs]
    MF, sd = matched(D)
    return X.fw_record(S, MF / sd, grid, DEFAULT_TAUS, Gs, Qs, truth)


def run_test(out_json, tmpdir, cfg=None, n=N_TEST, with_mtfw=True):
    """The qm_report.txt test: windows 0..n-1 of the seed-999 stream per cell, QuakeMigrate end to end
    (run_quakemigrate: archive, LUT, detect, public Trigger.trigger()) and MT-FW on the same windows."""
    from task0_common import VP, VS
    cfg = {**QM_CFG, **(cfg or {})}
    res = {"cfg": {k: (list(v) if isinstance(v, tuple) else v) for k, v in cfg.items()}, "cells": {}}
    for cell, scen in CELLS.items():
        wins = y_stream(scen, n)
        qrec, frec, evs, secs, fsecs = [], [], [], [], []
        for S, D, tr, taus in wins:
            t0 = time.time()
            ev = run_quakemigrate(S, D, tmpdir, VP, VS, cfg)
            secs.append(time.time() - t0)
            evs.append(ev)
            qrec.append(qm_records(ev, tr))
            if with_mtfw:
                t0 = time.time()
                frec.append(mtfw_record(S, D, tr))
                fsecs.append(time.time() - t0)
        ntrue = sum(len(w[2]) for w in wins)
        c = {"scen": scen, "nwin": len(wins), "ntrue": ntrue, "truth": [[list(map(float, s)) for s in w[2]] for w in wins],
             "taus": [w[3] for w in wins], "qm_events": evs, "qm_records": qrec, "qm_sec": secs,
             "qm_recall_at": froc_summary(qrec, ntrue)}
        if with_mtfw:
            c.update({"mtfw_records": [[(float(a), int(b)) for a, b in r] for r in frec], "mtfw_sec": fsecs,
                      "mtfw_recall_at": froc_summary(frec, ntrue)})
        res["cells"][cell] = c
        print(f"   {cell}: QM {c['qm_recall_at']} ({np.mean(secs):.2f} s/window)"
              + (f" | MT-FW {c['mtfw_recall_at']}" if with_mtfw else ""), flush=True)
    allq = [r for c in res["cells"].values() for r in c["qm_records"]]
    nt = sum(c["ntrue"] for c in res["cells"].values())
    res["pooled"] = {"ntrue": nt, "qm_recall_at": froc_summary(allq, nt)}
    if with_mtfw:
        allf = [r for c in res["cells"].values() for r in c["mtfw_records"]]
        res["pooled"]["mtfw_recall_at"] = froc_summary(allf, nt)
    with open(out_json, "w") as fh:
        json.dump(res, fh)
    return res


def locate_check(STA, D, tmpdir, vp, vs, cfg=None, position="centred", locate_marginal_window=1.0):
    """Runs QuakeMigrate's locate() on the triggered events of one window (by default centred STA/LTA
    onsets, as QuakeMigrate recommends for locate, same filters and windows otherwise; the marginal
    window of locate is widened to locate_marginal_window, otherwise the earlier centred-onset peak
    falls outside the calibrated 0.25-0.5 s window and locate discards the event) and returns
    {"detect": detect-stage events, "locate": [{spline/gau: (x, y, z), t0}]}."""
    from quakemigrate import QuakeScan
    from quakemigrate.io import Archive
    cfg = {**QM_CFG, **(cfg or {})}
    work = tempfile.mkdtemp(prefix="qmloc_", dir=tmpdir)
    try:
        lut = prepare(STA, D, work, vp, vs, cfg)
        detect(STA, work, lut, cfg)
        ev = run_trigger(work, lut, cfg, name="")
        if not ev:
            return []
        archive = Archive(archive_path=os.path.join(work, "archive"), stations=stations_df(STA),
                          archive_format="YEAR/JD/STATION", resample=(int(cfg["fs"]) != FS_NATIVE),
                          upfactor=int(cfg["upfactor"]))
        onset = _onset({**cfg, "position": position})
        tw0, tw1 = cfg["trigger_window"]
        with _quiet():
            scan = QuakeScan(archive, lut, onset=onset, run_path=work, run_name="qm", log=False,
                             threads=int(cfg["threads"]), marginal_window=float(locate_marginal_window),
                             plot_event_summary=False)
            scan.locate(starttime=str(T_REF + tw0 - 1.0), endtime=str(T_REF + tw1 + 1.0))
        out = []
        evd = os.path.join(work, "qm", "locate", "events")
        for f in sorted(os.listdir(evd)) if os.path.isdir(evd) else []:
            e = pd.read_csv(os.path.join(evd, f))
            locs = {}
            for tag, cx, cy, cz in (("spline", "X", "Y", "Z"), ("gau", "GAU_X", "GAU_Y", "GAU_Z"),
                                    ("cov", "COV_X", "COV_Y", "COV_Z")):
                if cx in e:
                    x, y = lonlat_to_km(e[cx].to_numpy(float), e[cy].to_numpy(float))
                    locs[tag] = (float(x[0]), float(y[0]), float(e[cz].iloc[0]))
            locs["t0"] = float(UTCDateTime(e["DT"].iloc[0]) - T_REF) if "DT" in e else None
            out.append(locs)
        return {"detect": ev, "locate": out}
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["calibrate", "test"])
    ap.add_argument("--out", default=os.path.join(tempfile.gettempdir(), "qm_out"))
    ap.add_argument("--workers", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.stage == "calibrate":
        for tag, extra in (("default", None), ("3c", CHANNELS_3C)):
            pkl = os.path.join(a.out, f"cal999_{tag}.pkl")
            calibrate_detect(pkl, os.path.join(a.out, "work"), workers=a.workers, extra=extra)
            calibrate_trigger(pkl, os.path.join(a.out, f"cal999_{tag}_trigger.json"))
    else:
        run_test(os.path.join(a.out, "test999_default.json"), os.path.join(a.out, "work"), cfg=QM_CFG)
        run_test(os.path.join(a.out, "test999_3c.json"), os.path.join(a.out, "work"), cfg=QM_CFG_3C,
                 with_mtfw=False)
