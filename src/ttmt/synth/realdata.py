"""Real continuous records at their native rate, and the empirical-Green's-function (EGF) library cut from them.

Semi-synthetic truth (ATLAS-0, truth design 1): recorded three-component waveforms of catalogued events are injected
into real continuous noise. This module reads the registered development datasets of ``ttmt.datasets`` (Ridgecrest
2019, Hi-net 2012) at 100 Hz, instrument-scaled, demeaned and detrended as ``Dataset._preprocess`` does, components
(E, N, Z up), and builds the EGF library:
- candidate events are catalogued events with no other catalogued event of magnitude >= M - ``iso_dm`` within
  [-``iso_pre``, +``iso_post``] s of their origin (so the cut window holds one event);
- arrival times are measured on the clean EGF itself: an AIC onset (Maeda 1985) on the 2-20 Hz band-passed vertical
  (P) and horizontal energy (S) in a window around the arrival predicted by a 1-D model from the catalogued
  hypocentre; a phase is kept when the onset is within ``max_dev`` s of the prediction and its pick SNR exceeds
  ``min_pick_snr``, otherwise its time is nan (absent in the truth);
- the time blocks of the record are split into an EGF role and a noise role by a seeded hash (disjoint by event and
  by time block, CR-06): EGFs come only from EGF blocks, injections go only into noise blocks.

The 1-D models are the published ones the project already uses for these networks (Hadley and Kanamori 1977 for
southern California, JMA2001 for Japan), re-declared here as data; arrival predictions only seed the AIC search.
"""
import hashlib
from dataclasses import dataclass, field

import numpy as np
from scipy.signal import butter, sosfiltfilt

from .media import LayeredModel, first_arrival

FS = 100.0
BLOCK_S = 300.0                                 # time block of the role split and of the tune/read split

_HK = LayeredModel(thickness=(5.5, 10.5, 16.0), vp=(5.5, 6.3, 6.7, 7.8), vs=tuple(v / 1.73 for v in (5.5, 6.3, 6.7, 7.8)),
                   rho=(2.4, 2.67, 2.8, 3.3), name="hadley_kanamori_1977")
_JMA_PTS = [(0.0, 4.800), (1.5, 5.220), (2.5, 5.450), (4.0, 5.720), (5.0, 5.820), (15.0, 6.220), (22.0, 6.544),
            (34.0, 7.188), (41.0, 7.494), (44.0, 7.598), (48.0, 7.688), (52.5, 7.750), (65.0, 7.845)]


def _jma2001(dz=0.5, zmax=60.0):
    z = np.arange(0.0, zmax, dz)
    vp = np.interp(z + dz / 2, [p[0] for p in _JMA_PTS], [p[1] for p in _JMA_PTS])
    return LayeredModel(thickness=tuple([dz] * (len(z) - 1)), vp=tuple(vp), vs=tuple(vp / 1.73),
                        rho=tuple([2.7] * len(z)), name="jma2001")


PRED_MODELS = {"ridgecrest_2019": _HK, "hinet_2012": None}


def pred_model(name):
    if name == "hinet_2012":
        if PRED_MODELS["hinet_2012"] is None:
            PRED_MODELS["hinet_2012"] = _jma2001()
        return PRED_MODELS["hinet_2012"]
    return PRED_MODELS[name]


# ------------------------------------------------------------------------------------------------ records
@dataclass
class Record:
    """a continuous record: data (ns, 3, n) at fs, components (E, N, Z up), sample 0 at record time t0 (s from the
    dataset's t0_utc); stations codes and local xyz (km, z down, 0 at the surface); catalogue in the same frame."""
    dataset: str
    codes: list
    xyz: np.ndarray
    data: np.ndarray
    fs: float
    t0: float
    catalogue: list
    centre: tuple
    depth_km: np.ndarray = None
    meta: dict = field(default_factory=dict)

    @property
    def t1(self):
        return self.t0 + self.data.shape[-1] / self.fs

    def cut(self, a, b):
        i, j = int(round((a - self.t0) * self.fs)), int(round((b - self.t0) * self.fs))
        if i < 0 or j > self.data.shape[-1]:
            raise ValueError("cut outside the record")
        return self.data[:, :, i:j]


def load_record(name):
    """the whole registered record of a dataset ('ridgecrest_2019' or 'hinet_2012') at 100 Hz."""
    from .. import datasets as DS
    ds = {"ridgecrest_2019": DS.RIDGECREST, "hinet_2012": DS.HINET}[name]
    centre = ds.frame()
    if name == "ridgecrest_2019":
        import csv
        import os
        names = sorted({r["station"] for r in csv.DictReader(open(os.path.join(ds.data, "stations.csv")))
                        if r["channel"].startswith("HH")})
        t_end = 3600.0
    else:
        names = ds.subnets["H7"]
        t_end = 2000.0
    STA, sens = ds.stations(names, centre)
    raw, _d16 = ds.load_waveforms(names, sens, t_end)
    cat = ds.catalogue_local(centre)
    depth = np.zeros(len(names))
    if name == "hinet_2012":                         # Hi-net sensors sit in boreholes (SAC header stdp or stel)
        for i, n in enumerate(names):
            h = ds._sac(n, "U").stats.sac
            depth[i] = float(getattr(h, "stdp", 0.0) or 0.0) / 1000.0
    return Record(name, list(names), np.asarray(STA, float), raw, FS, 0.0, cat, centre, depth,
                  {"t0_utc": ds.t0_utc, "units": "m/s (response-scaled)" if name == "ridgecrest_2019" else "counts",
                   "physical_units": name == "ridgecrest_2019"})


def blocks(rec, block_s=BLOCK_S):
    """(start, end) of the consecutive time blocks of a record (the last partial block dropped)."""
    n = int((rec.t1 - rec.t0) // block_s)
    return [(rec.t0 + i * block_s, rec.t0 + (i + 1) * block_s) for i in range(n)]


def block_hash(dataset, a, salt):
    return int.from_bytes(hashlib.sha256(f"{dataset}|{a:.3f}|{salt}".encode()).digest()[:8], "little") / 2 ** 64


def block_roles(rec, salt="egf-noise-v1", egf_frac=0.4, block_s=BLOCK_S):
    """role of each block, 'egf' or 'noise', by a seeded hash; at least one of each."""
    bl = blocks(rec, block_s)
    u = np.array([block_hash(rec.dataset, a, salt) for a, _b in bl])
    roles = np.where(u < egf_frac, "egf", "noise")
    if not (roles == "egf").any():
        roles[int(np.argmin(u))] = "egf"
    if not (roles == "noise").any():
        roles[int(np.argmax(u))] = "noise"
    return list(zip(bl, roles.tolist()))


def tune_read(dataset, block_start, salt="tune-read-v1", tune_frac=0.3):
    """'tune' or 'read' partition of a noise block (ATLAS-0, Tuning: 30/70 by a hash of the noise time)."""
    return "tune" if block_hash(dataset, block_start, salt) < tune_frac else "read"


def partition_blocks(dataset, starts, salt="tune-read-v1", tune_frac=0.3):
    """tune_read for a list of block starts, with at least one block in each partition when there are two or more
    (a record with few blocks can otherwise draw no tuning block: the block of lowest hash moves to 'tune', or of
    highest hash to 'read')."""
    u = np.array([block_hash(dataset, a, salt) for a in starts])
    part = np.where(u < tune_frac, "tune", "read")
    if len(starts) >= 2 and not (part == "tune").any():
        part[int(np.argmin(u))] = "tune"
    if len(starts) >= 2 and not (part == "read").any():
        part[int(np.argmax(u))] = "read"
    return dict(zip([float(a) for a in starts], part.tolist()))



def tune_read_split(rec, block_s=BLOCK_S):
    """ADR-0009 item 4 for one record: [((start, end), role)] for every block, role 'egf' (in neither partition),
    'tune' or 'read' (``partition_blocks`` on the noise-role blocks, with the fallback both ways)."""
    roles = block_roles(rec, block_s=block_s)
    part = partition_blocks(rec.dataset, [bl[0] for bl, r in roles if r == "noise"])
    return [(bl, r if r == "egf" else part[float(bl[0])]) for bl, r in roles]

# ------------------------------------------------------------------------------------------------ picking
_SOS = {}


def bandpass(x, fs, band=(2.0, 20.0)):
    key = (fs, band)
    if key not in _SOS:
        _SOS[key] = butter(4, band, btype="band", fs=fs, output="sos")
    return sosfiltfilt(_SOS[key], x, axis=-1)


def aic_onset(x):
    """index of the AIC minimum (Maeda 1985) of a 1-D trace."""
    n = len(x)
    k = np.arange(1, n - 1)
    c1 = np.cumsum(x * x)
    tot = c1[-1]
    v1 = c1[k - 1] / k
    v2 = (tot - c1[k - 1]) / (n - k)
    aic = k * np.log(np.maximum(v1, 1e-300)) + (n - k - 1) * np.log(np.maximum(v2, 1e-300))
    return int(k[np.argmin(aic)])


def pick_event(seg, fs, t_seg0, t_pred, half=(1.0, 1.5), max_dev=(0.5, 0.8), min_pick_snr=4.0):
    """AIC onsets of P (vertical) and S (horizontal energy) of one station segment (3, n) (E, N, Z up) starting at
    record time t_seg0; t_pred (2,) predicted times. Returns times (2,) (nan when rejected) and pick SNRs (2,)."""
    y = bandpass(seg, fs)
    out, snr = np.full(2, np.nan), np.full(2, np.nan)
    tr = (y[2], np.sqrt(y[0] ** 2 + y[1] ** 2))
    for b in range(2):
        a = int(round((t_pred[b] - half[b] - t_seg0) * fs))
        z = int(round((t_pred[b] + half[b] - t_seg0) * fs))
        if a < 1 or z > seg.shape[-1] - 1:
            continue
        x = tr[b][a:z]
        if b == 1:
            x = x - np.mean(x[: max(5, int(0.2 * fs))])
        k = aic_onset(x)
        t = t_seg0 + (a + k) / fs
        pre = tr[b][max(0, a + k - int(1.0 * fs)):a + k]
        post = tr[b][a + k:a + k + int(1.0 * fs)]
        s = float(np.max(np.abs(post)) / (np.std(pre) + 1e-30)) if len(pre) > 10 and len(post) else 0.0
        snr[b] = s
        if abs(t - t_pred[b]) <= max_dev[b] and s >= min_pick_snr:
            out[b] = t
    if np.isfinite(out).all() and out[1] <= out[0] + 0.2:          # S must follow P
        out[1] = np.nan
    return out, snr


# ------------------------------------------------------------------------------------------------ EGF library
@dataclass
class EGF:
    """one recorded event cut from the record: wave (ns, 3, n) from record time t_cut; origin t, xyz, mag; arr (ns, 2)
    measured arrival times (record time, nan absent); pick_snr (ns, 2); block start; id."""
    id: str
    t: float
    xyz: np.ndarray
    mag: float
    t_cut: float
    wave: np.ndarray
    arr: np.ndarray
    pick_snr: np.ndarray
    block: float
    noise_sd: np.ndarray


def egf_library(rec, pre=5.0, post=35.0, iso_pre=8.0, iso_post=30.0, iso_dm=2.0, min_mag=None, min_phases=6,
                salt="egf-noise-v1"):
    """EGF library of a record, from the EGF-role blocks only (see the module docstring)."""
    roles = block_roles(rec, salt)
    egf_blocks = [bl for bl, r in roles if r == "egf"]
    cat = rec.catalogue
    M = pred_model(rec.dataset)
    out = []
    for e in cat:
        blk = next((a for a, b in egf_blocks if a <= e["t"] - pre and e["t"] + post <= b), None)
        if blk is None or (min_mag is not None and e["M"] < min_mag):
            continue
        if any(o is not e and -iso_pre <= o["t"] - e["t"] <= iso_post and o["M"] >= e["M"] - iso_dm for o in cat):
            continue
        t_cut = e["t"] - pre
        try:
            seg = rec.cut(t_cut, e["t"] + post)
        except ValueError:
            continue
        ns = len(rec.codes)
        arr, ps = np.full((ns, 2), np.nan), np.full((ns, 2), np.nan)
        zs = max(float(e["z"]), 0.5)
        X = np.hypot(rec.xyz[:, 0] - e["x"], rec.xyz[:, 1] - e["y"])
        tp = first_arrival(M, zs, X, "P")[0]
        ts = first_arrival(M, zs, X, "S")[0]
        for i in range(ns):
            arr[i], ps[i] = pick_event(seg[i], rec.fs, t_cut, e["t"] + np.array([tp[i], ts[i]]))
        if np.isfinite(arr).sum() < min_phases:
            continue
        nsd = np.std(bandpass(seg[:, :, : int(0.8 * pre * rec.fs)], rec.fs), axis=-1)
        out.append(EGF(f"{rec.dataset}:{e['t']:.2f}", float(e["t"]), np.array([e["x"], e["y"], e["z"]], float),
                       float(e["M"]), t_cut, seg.copy(), arr, ps, blk, nsd))
    return out


def template_candidates(rec, pre=5.0, post=35.0, min_mag=None, salt="egf-noise-v1"):
    """default atom library for the CR-06 correlation: every catalogued event inside the noise-role blocks (the events a
    template or deflation method would take from the record the injections go into), cut [t - pre, t + post]."""
    roles = block_roles(rec, salt)
    nb = [bl for bl, r in roles if r == "noise"]
    out, ids = [], []
    for e in rec.catalogue:
        if min_mag is not None and e["M"] < min_mag:
            continue
        if not any(a <= e["t"] < b for a, b in nb):
            continue
        try:
            out.append(rec.cut(e["t"] - pre, e["t"] + post).copy())
            ids.append(f"{rec.dataset}:{e['t']:.2f}")
        except ValueError:
            continue
    return out, ids
