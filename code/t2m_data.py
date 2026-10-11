"""
t2m_data.py - Task 2M: real continuous noise of the Alto Tiberina array (IV network) prepared for injection, for ANY
day file set of avuan/PyMPA37 main.pympa.dir/24h/ (c16_data.py generalised; the processing is identical).
Reads untrusted data files; run with
    python -I t2m_data.py DAY STATIONS DATADIR INV ZMAP OUTDIR QC [EXTRA_CAT ...]
  DAY        YYMMDD (e.g. 101231)
  STATIONS   comma-separated station codes, in the order used for the station axis of the output
  DATADIR    directory with the files DAY.STA.CHA (three components per station, 25 Hz)
  INV        StationXML inventory (inv.IV); coordinates of the epoch active on DAY
  ZMAP       catalogue (templates.zmap); rows of DAY are excluded with +-60 s
  QC         "none" (C16: only padded samples are invalid) or "qc1" (step 1b below)
  EXTRA_CAT  optional PyMPA detection lists (e.g. 200.100723.cat): rows "<n> <UTC time> ..." of DAY are also excluded
No project module is imported (python -I).

Steps (identical to c16_data.py with QC = none; the day, the station set, the exclusion lists and QC are parameters):
  1. read the three-component day files (obspy), merge, trim/pad to [DAY 00:00:00, +86400 s) at 25 Hz (nearest
     sample; padded samples are invalid).  When the nearest-sample trim of a trace with a half-sample offset returns
     2160001 samples, the last one is dropped (first sample within 0.02 s of 00:00:00); C16's day never needed it.
  1b. QC "qc1" (data-quality rule; not in C16): on the raw 25 Hz counts of each channel, a sample is also INVALID if
     (i) it belongs to a run of >= 25 identical consecutive values (>= 1 s; flat-lined / sample-and-hold telemetry fill),
     or (ii) |x - median| > 1e4 * sd_rob, median and sd_rob = 1.4826 * MAD over the valid samples not flagged by (i)
     (digital spikes; also any real arrival of that size).  Invalid samples are handled exactly as padded ones (set to 0
     before the trend fit and the filters; every slot with an invalid sample of any station within 60 s is ineligible).
  2. per channel: mean and linear trend removed (fit on valid samples); zero-phase Butterworth high-pass 0.5 Hz (order 4)
     and anti-alias low-pass 6.5 Hz (order 8), sosfiltfilt at 25 Hz; resampling 25 -> 16 Hz by resample_poly(16, 25).
  3. components to the generator frame: x = East, y = North, z = DOWN (= -Z).
  4. coordinates: flat earth around 43.315 N, 12.388 E, then shifted by the centroid of the stations of THIS station set
     (for the 14 stations of 2010-07-23 this is exactly C16's frame).  Elevations ignored (z = 0).
  5. slots: 40 s segments starting at 0, 40, 80, ... s; ELIGIBLE if it does not intersect [t - 60 s, t + 60 s] of any
     catalogued origin of DAY (ZMAP rows of DAY and EXTRA_CAT rows of DAY), lies >= 60 s from both ends of the day and
     no sample of ANY station of the set is invalid within 60 s of the slot.
  6. amplitude rule per channel: s_ch = median over the eligible slots of the robust in-band level (median |MF| / 0.6745,
     MF = 2 Hz Ricker matched filter); channel multiplied by sd_gen / s_ch, sd_gen = 0.14 per sample.
Outputs: OUTDIR/noise.npy float32 (nsta, 3, 1382400); OUTDIR/meta.json (same keys as c16_meta.json plus day, t0,
valid_hours and per-station coverage).
"""
import json
import os
import sys

import numpy as np
from scipy.signal import butter, resample_poly, sosfiltfilt

LAT0, LON0 = 43.315, 12.388
FS_IN, FS_OUT = 25.0, 16.0
DAY = 86400
NSLOT = 640
SLOT_S = 40.0
EXCL = 60.0
R_EARTH = 6371.0
SNR, NOISE = 7.0, 0.02


def ricker_ker():
    dt = 1.0 / FS_OUT
    t = np.arange(-1.0, 1.0 + 1e-9, dt)
    x = (np.pi * 2.0 * t) ** 2
    return (1 - 2 * x) * np.exp(-x)


def qc1(d, bad, run_min=25, k_spike=1e4):
    """flat runs (>= run_min identical consecutive finite values) and spikes (|x - med| > k_spike * 1.4826 MAD)."""
    flat = np.zeros(len(d), bool)
    fin = np.isfinite(d)
    same = np.r_[False, (d[1:] == d[:-1]) & fin[1:] & fin[:-1]]
    edges = np.flatnonzero(np.diff(np.r_[0, same.astype(np.int8), 0]))
    for a, b in zip(edges[::2], edges[1::2]):
        if b - a + 1 >= run_min:
            flat[a - 1:b] = True
    ref = d[fin & ~flat]
    med = np.median(ref)
    sd = 1.4826 * np.median(np.abs(ref - med))
    spike = fin & ~flat & (np.abs(d - med) > k_spike * sd)
    return bad | flat | spike, int(flat.sum()), int(spike.sum())


def main(day, stations, ddir, inv_path, zmap_path, odir, qc, extra_cats):
    assert qc in ("none", "qc1"), qc
    import obspy
    from obspy import UTCDateTime
    os.makedirs(odir, exist_ok=True)
    yy, mm, dd_ = 2000 + int(day[:2]), int(day[2:4]), int(day[4:6])
    t0 = UTCDateTime(f"{yy:04d}-{mm:02d}-{dd_:02d}T00:00:00")
    inv = obspy.read_inventory(inv_path)
    coords = {}
    for net in inv:
        for sta in net:
            if sta.code in stations and sta.start_date <= t0 and (sta.end_date is None or sta.end_date > t0 + DAY):
                coords[sta.code] = (float(sta.latitude), float(sta.longitude), float(sta.elevation))
    assert set(coords) == set(stations), set(stations) - set(coords)
    lat = np.array([coords[s][0] for s in stations])
    lon = np.array([coords[s][1] for s in stations])
    x = R_EARTH * np.cos(np.radians(LAT0)) * np.radians(lon - LON0)
    y = R_EARTH * np.radians(lat - LAT0)
    cx, cy = float(x.mean()), float(y.mean())
    sos_hp = butter(4, 0.5, btype="highpass", fs=FS_IN, output="sos")
    sos_lp = butter(8, 6.5, btype="lowpass", fs=FS_IN, output="sos")
    n_in = int(DAY * FS_IN)
    n_out = int(DAY * FS_OUT)
    out = np.zeros((len(stations), 3, n_out), np.float32)
    valid = np.ones(n_out, bool)
    cover = {}
    info = {}
    for si, s in enumerate(stations):
        for ci, comp in enumerate("ENZ"):
            fn = [f for f in os.listdir(ddir) if f.startswith(f"{day}.{s}.") and f.endswith(comp) and len(f.split(".")) == 3]
            assert len(fn) == 1, (s, comp, fn)
            st = obspy.read(os.path.join(ddir, fn[0]))
            st.merge(method=1, fill_value=None)
            assert len(st) == 1
            tr = st[0]
            assert abs(tr.stats.sampling_rate - FS_IN) < 1e-6
            masked = np.ma.is_masked(tr.data)
            tr.trim(t0, t0 + DAY - 1.0 / FS_IN, pad=True, fill_value=None, nearest_sample=True)
            d = np.ma.filled(np.ma.masked_invalid(np.ma.asarray(tr.data, float)), np.nan) if np.ma.isMaskedArray(tr.data) \
                else np.asarray(tr.data, float)
            assert n_in <= len(d) <= n_in + 1, (s, comp, len(d))
            n_extra = len(d) - n_in
            d = d[:n_in]
            bad = ~np.isfinite(d)
            n_pad = int(bad.sum())
            n_flat = n_spike = 0
            if qc == "qc1":
                bad, n_flat, n_spike = qc1(d, bad)
            d = np.where(bad, 0.0, d)
            k = np.arange(n_in)
            ok = ~bad
            p = np.polyfit(k[ok], d[ok], 1)
            d = np.where(ok, d - np.polyval(p, k), 0.0)
            d = sosfiltfilt(sos_lp, sosfiltfilt(sos_hp, d))
            dd = resample_poly(d, 16, 25)[:n_out]
            if bad.any():
                tb = np.nonzero(bad)[0] / FS_IN
                kb = np.unique(np.clip(np.rint(tb * FS_OUT).astype(int), 0, n_out - 1))
                valid[kb] = False
            sign = -1.0 if comp == "Z" else 1.0
            out[si, ci] = (sign * dd).astype(np.float32)
            info[f"{s}.{tr.stats.channel}"] = {"start": str(tr.stats.starttime), "n_bad": int(bad.sum()),
                                               "had_gaps": bool(masked), "n_extra_dropped": int(n_extra),
                                               "n_pad": n_pad, "n_flat": int(n_flat), "n_spike": int(n_spike)}
            cover[f"{s}.{tr.stats.channel}"] = round(float(ok.mean()) * 24.0, 3)
    ev = []
    for line in open(zmap_path):
        f = line.split()
        if len(f) >= 10 and int(f[2]) == yy and int(f[3]) == mm and int(f[4]) == dd_:
            ev.append(int(f[7]) * 3600 + int(f[8]) * 60 + float(f[9]))
    n_tmpl = len(ev)
    for cat in extra_cats:
        for line in open(cat):
            f = line.split()
            if len(f) >= 2 and f[1].startswith(f"{yy:04d}-{mm:02d}-{dd_:02d}T"):
                ev.append(float(UTCDateTime(f[1]) - t0))
    ev = np.sort(np.array(ev))
    nsl = int(DAY // SLOT_S)
    elig = []
    for k in range(nsl):
        a, b = k * SLOT_S, (k + 1) * SLOT_S
        if a < EXCL or b > DAY - EXCL:
            continue
        if np.any((ev > a - EXCL) & (ev < b + EXCL)):
            continue
        i0, i1 = int(a * FS_OUT), int(b * FS_OUT)
        lo, hi = max(0, i0 - int(EXCL * FS_OUT)), min(n_out, i1 + int(EXCL * FS_OUT))
        if not valid[lo:hi].all():
            continue
        elig.append(k)
    elig = np.array(elig)
    ker = ricker_ker()
    sd_gen = SNR * NOISE * float(np.sqrt(np.sum(ker ** 2)))
    lev = np.zeros((len(stations), 3, len(elig)))
    for si in range(len(stations)):
        for ci in range(3):
            mf = np.convolve(out[si, ci].astype(float), ker[::-1], mode="same")
            seg = mf[: nsl * NSLOT].reshape(nsl, NSLOT)[elig]
            lev[si, ci] = np.median(np.abs(seg), axis=1) / 0.6745
    s_ch = np.median(lev, axis=2)
    scale = sd_gen / s_ch
    out *= scale[:, :, None].astype(np.float32)
    np.save(os.path.join(odir, "noise.npy"), out)
    rel = lev / s_ch[:, :, None]
    meta = {"stations": list(stations), "lat": lat.tolist(), "lon": lon.tolist(),
            "elev_m": [coords[s][2] for s in stations], "x_km": (x - cx).tolist(), "y_km": (y - cy).tolist(),
            "centroid_shift_km": [cx, cy], "slot_s": SLOT_S, "n_slots": nsl, "eligible": elig.tolist(),
            "events_s": ev.tolist(), "n_templates": n_tmpl, "sd_gen": sd_gen, "scale": scale.tolist(),
            "level_counts": s_ch.tolist(), "channels": info,
            "rel_level_q": {q: np.percentile(rel, q).tolist() for q in (5, 25, 50, 75, 95, 99)},
            "rel_level_per_slot_max": np.max(rel, axis=(0, 1)).tolist(),
            "day": day, "t0": str(t0), "extra_cats": [os.path.basename(c) for c in extra_cats], "qc": qc,
            "valid_hours": round(float(valid.mean()) * 24.0, 3), "coverage_hours": cover}
    with open(os.path.join(odir, "meta.json"), "w") as fh:
        json.dump(meta, fh)
    r = np.hypot(x - cx, y - cy)
    print("day", day, "stations", len(stations), "eligible slots", len(elig), "of", nsl, f"({len(elig) * SLOT_S / 3600:.2f} h)",
          "events", len(ev), "(templates", n_tmpl, ")", "valid hours (all stations)", meta["valid_hours"])
    print("x_km", np.round(x - cx, 1).tolist())
    print("y_km", np.round(y - cy, 1).tolist())
    print("distance to frame origin (km) min %.1f max %.1f" % (r.min(), r.max()))
    print("relative in-band level quantiles (5,25,50,75,95,99):",
          [round(float(np.percentile(rel, q)), 3) for q in (5, 25, 50, 75, 95, 99)])


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2].split(","), sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8:])
