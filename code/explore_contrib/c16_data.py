"""
c16_data.py - C16: real continuous noise of the Alto Tiberina array (IV network, 2010-07-23) prepared for injection.
SEED 999 ONLY.  Reads untrusted data files; run with  python -I c16_data.py DATADIR OUTDIR  (no project import).

Steps (all fixed before any detection run):
  1. read the 14 three-component day files (obspy), merge, trim/pad to [2010-07-23T00:00:00, +86400 s) at 25 Hz
     (nearest sample; the pad is reported and the padded samples are excluded through the slot rule).
  2. per channel: remove the mean and a linear trend; zero-phase Butterworth band-pass: high-pass 0.5 Hz (order 4) and
     anti-alias low-pass 6.5 Hz (order 8), both sosfiltfilt at 25 Hz; resampling 25 -> 16 Hz by scipy.signal.
     resample_poly(x, 16, 25) (polyphase FIR with its own Kaiser low-pass).  16 Hz = 1 / task0_common.DT.
  3. components to the generator frame: x = East, y = North, z = DOWN (= -Z).
  4. coordinates: station lat / lon of the epoch active on 2010-07-23 (inv.IV, StationXML); flat earth around the
     cluster 43.315 N, 12.388 E (x = R cos(lat0) dlon, y = R dlat, R = 6371 km), then shifted by the centroid of the 14
     stations so that the generator's source box [-15, 15]^2 km and the detection lattice (+-21 km) lie inside the real
     network.  Elevations are ignored (z = 0, as the generator's stations).
  5. slots: 40 s segments (640 samples, NT of the project) starting at 0, 40, 80, ... s; a slot is ELIGIBLE if it does
     not intersect [t - 60 s, t + 60 s] of any catalogued origin of that day (templates.zmap rows of 2010-07-23 and the
     detections listed in 200.100723.cat) and lies at least 60 s from both ends of the day and from any padded sample.
  6. amplitude rule (per channel, fixed): s_ch = median over the eligible slots of the robust in-band level of the slot
     (median |MF| / 0.6745, MF = the project's 2 Hz Ricker matched filter, xc.matched); the channel is multiplied by
     sd_gen / s_ch, sd_gen = SNR * noise * ||KER||_2 = the in-band level of the generator's white noise (0.14 per sample).
     The day-to-day level changes of a channel (busy hours, transients) are kept; the gain differences between
     channels (counts of different instruments, unknown without responses) are removed.
Outputs: OUTDIR/c16_noise.npy float32 (14, 3, 1382400); OUTDIR/c16_meta.json.
"""
import json
import os
import sys

import numpy as np
from scipy.signal import butter, resample_poly, sosfiltfilt

STATIONS = ["ATBU", "ATCA", "ATCC", "ATFO", "ATLO", "ATPC", "ATPI", "ATVO", "BADI", "FRON", "MURB", "PARC", "PIEI",
            "SSFR"]
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


def main(ddir, odir):
    import obspy
    from obspy import UTCDateTime
    os.makedirs(odir, exist_ok=True)
    t0 = UTCDateTime("2010-07-23T00:00:00")
    inv = obspy.read_inventory(os.path.join(ddir, "inv.IV"))
    coords = {}
    for net in inv:
        for sta in net:
            if sta.code in STATIONS and sta.start_date <= t0 and (sta.end_date is None or sta.end_date > t0 + DAY):
                coords[sta.code] = (float(sta.latitude), float(sta.longitude), float(sta.elevation))
    assert set(coords) == set(STATIONS), set(STATIONS) - set(coords)
    lat = np.array([coords[s][0] for s in STATIONS])
    lon = np.array([coords[s][1] for s in STATIONS])
    x = R_EARTH * np.cos(np.radians(LAT0)) * np.radians(lon - LON0)
    y = R_EARTH * np.radians(lat - LAT0)
    cx, cy = float(x.mean()), float(y.mean())
    sos_hp = butter(4, 0.5, btype="highpass", fs=FS_IN, output="sos")
    sos_lp = butter(8, 6.5, btype="lowpass", fs=FS_IN, output="sos")
    n_in = int(DAY * FS_IN)
    n_out = int(DAY * FS_OUT)
    out = np.zeros((len(STATIONS), 3, n_out), np.float32)
    valid = np.ones(n_out, bool)
    info = {}
    for si, s in enumerate(STATIONS):
        for ci, comp in enumerate("ENZ"):
            fn = [f for f in os.listdir(ddir) if f.startswith(f"100723.{s}.") and f.endswith(comp) and len(f.split(".")) == 3]
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
            assert len(d) == n_in, (s, comp, len(d))
            bad = ~np.isfinite(d)
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
                                               "had_gaps": bool(masked)}
    # catalogue origins of the day
    ev = []
    for line in open(os.path.join(ddir, "templates.zmap")):
        f = line.split()
        if len(f) >= 10 and int(f[2]) == 2010 and int(f[3]) == 7 and int(f[4]) == 23:
            ev.append(int(f[7]) * 3600 + int(f[8]) * 60 + float(f[9]))
    n_tmpl = len(ev)
    for line in open(os.path.join(ddir, "200.100723.cat")):
        f = line.split()
        if len(f) >= 2 and f[1].startswith("2010-07-23T"):
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
    lev = np.zeros((len(STATIONS), 3, len(elig)))
    for si in range(len(STATIONS)):
        for ci in range(3):
            mf = np.convolve(out[si, ci].astype(float), ker[::-1], mode="same")
            seg = mf[: nsl * NSLOT].reshape(nsl, NSLOT)[elig]
            lev[si, ci] = np.median(np.abs(seg), axis=1) / 0.6745
    s_ch = np.median(lev, axis=2)
    scale = sd_gen / s_ch
    out *= scale[:, :, None].astype(np.float32)
    np.save(os.path.join(odir, "c16_noise.npy"), out)
    rel = lev / s_ch[:, :, None]
    meta = {"stations": STATIONS, "lat": lat.tolist(), "lon": lon.tolist(),
            "elev_m": [coords[s][2] for s in STATIONS], "x_km": (x - cx).tolist(), "y_km": (y - cy).tolist(),
            "centroid_shift_km": [cx, cy], "slot_s": SLOT_S, "n_slots": nsl, "eligible": elig.tolist(),
            "events_s": ev.tolist(), "n_templates": n_tmpl, "sd_gen": sd_gen, "scale": scale.tolist(),
            "level_counts": s_ch.tolist(), "channels": info,
            "rel_level_q": {q: np.percentile(rel, q).tolist() for q in (5, 25, 50, 75, 95, 99)},
            "rel_level_per_slot_max": np.max(rel, axis=(0, 1)).tolist()}
    with open(os.path.join(odir, "c16_meta.json"), "w") as fh:
        json.dump(meta, fh)
    print("stations", len(STATIONS), "eligible slots", len(elig), "of", nsl, "events", len(ev), "(templates", n_tmpl, ")")
    print("x_km", np.round(x - cx, 1).tolist())
    print("y_km", np.round(y - cy, 1).tolist())
    print("relative in-band level quantiles (5,25,50,75,95,99):",
          [round(float(np.percentile(rel, q)), 3) for q in (5, 25, 50, 75, 95, 99)])


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
