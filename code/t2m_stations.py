"""
t2m_stations.py - Task 2M station-set rule for a noise day (data-quality only; no injection, no detector).
    python -I t2m_stations.py DAY STATIONS DATADIR ZMAP OUT.json [EXTRA_CAT ...]
For every station of STATIONS (the stations whose three day files exist) the invalid samples are those of t2m_data.py
with QC "qc1" (padding, flat runs >= 1 s, spikes > 1e4 robust sd).  E(S) = number of slots eligible for the station set
S (t2m_data.py slot rule: >= 60 s from the day ends, no catalogued origin within 60 s, no invalid sample of any station
of S within 60 s).  Rule (fixed before any detection run on these days):
    S = all stations; while |S| > 8: let s* maximise E(S minus s); if E(S minus s*) >= 1.2 E(S) (a station is dropped only
    when it costs more than a sixth of the usable slots), drop s*; otherwise stop.
    The day is USED only if |S| >= 8 (the unilateral cell needs 8 stations) and E(S) >= 300 slots.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t2m_data import DAY, EXCL, FS_IN, SLOT_S, qc1  # noqa: E402


def main(day, stations, ddir, zmap_path, outp, extra_cats):
    import obspy
    from obspy import UTCDateTime
    yy, mm, dd_ = 2000 + int(day[:2]), int(day[2:4]), int(day[4:6])
    t0 = UTCDateTime(f"{yy:04d}-{mm:02d}-{dd_:02d}T00:00:00")
    n_in = int(DAY * FS_IN)
    nsl = int(DAY // SLOT_S)
    ev = []
    for line in open(zmap_path):
        f = line.split()
        if len(f) >= 10 and int(f[2]) == yy and int(f[3]) == mm and int(f[4]) == dd_:
            ev.append(int(f[7]) * 3600 + int(f[8]) * 60 + float(f[9]))
    for cat in extra_cats:
        for line in open(cat):
            f = line.split()
            if len(f) >= 2 and f[1].startswith(f"{yy:04d}-{mm:02d}-{dd_:02d}T"):
                ev.append(float(UTCDateTime(f[1]) - t0))
    ev = np.array(ev)
    base = np.zeros(nsl, bool)
    for k in range(nsl):
        a, b = k * SLOT_S, (k + 1) * SLOT_S
        base[k] = a >= EXCL and b <= DAY - EXCL and not np.any((ev > a - EXCL) & (ev < b + EXCL))
    ok = {}
    info = {}
    for s in stations:
        inval = np.zeros(n_in, bool)
        nf = nsp = npad = 0
        for comp in "ENZ":
            fn = [f for f in os.listdir(ddir) if f.startswith(f"{day}.{s}.") and f.endswith(comp) and len(f.split(".")) == 3]
            assert len(fn) == 1, (s, comp, fn)
            st = obspy.read(os.path.join(ddir, fn[0]))
            st.merge(method=1, fill_value=None)
            tr = st[0]
            tr.trim(t0, t0 + DAY - 1.0 / FS_IN, pad=True, fill_value=None, nearest_sample=True)
            d = np.ma.filled(np.ma.masked_invalid(np.ma.asarray(tr.data, float)), np.nan) if np.ma.isMaskedArray(tr.data) \
                else np.asarray(tr.data, float)
            d = d[:n_in]
            bad = ~np.isfinite(d)
            npad += int(bad.sum())
            b2, f_, sp_ = qc1(d, bad)
            nf += f_
            nsp += sp_
            inval |= b2
        cs = np.r_[0, np.cumsum(inval)]
        o = np.zeros(nsl, bool)
        for k in range(nsl):
            lo = max(0, int((k * SLOT_S - EXCL) * FS_IN))
            hi = min(n_in, int(((k + 1) * SLOT_S + EXCL) * FS_IN))
            o[k] = base[k] and cs[hi] - cs[lo] == 0
        ok[s] = o
        info[s] = {"invalid_h": round(float(inval.mean()) * 24, 3), "pad_h_chan": round(npad / FS_IN / 3600, 3),
                   "flat_h_chan": round(nf / FS_IN / 3600, 3), "spike_samples": int(nsp), "eligible_alone": int(o.sum())}
        print(day, s, info[s], flush=True)

    def E(S):
        return int(np.all([ok[x] for x in S], axis=0).sum()) if S else 0
    S = list(stations)
    steps = [{"set": list(S), "E": E(S)}]
    while len(S) > 8:
        cand = sorted(((E([x for x in S if x != s]), -i, s) for i, s in enumerate(S)), reverse=True)
        e_new, _i, s_star = cand[0]
        if e_new >= 1.2 * E(S) and e_new > 0:
            S = [x for x in S if x != s_star]
            steps.append({"dropped": s_star, "set": list(S), "E": e_new})
        else:
            break
    used = len(S) >= 8 and E(S) >= 300
    out = {"day": day, "candidates": list(stations), "station_qc": info, "steps": steps, "selected": S, "E": E(S),
           "used": used}
    with open(outp, "w") as fh:
        json.dump(out, fh, indent=1)
    print(day, "SELECTED", ",".join(S), "eligible slots", E(S), "USED" if used else "NOT USED", flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2].split(","), sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6:])
