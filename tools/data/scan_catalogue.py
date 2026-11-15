"""Rank candidate sequences in a hypocentre catalogue, from catalogue fields only (no waveform is read).

Input: CSV files with the ComCat / NCEDC "ehpcsv" columns (time, latitude, longitude, depth, mag, magType, nst, gap,
rms, horizontalError, depthError, status, type, place). The catalogue is binned in cells of ``--cell`` degrees; for
every cell the busiest ``--days``-day window of the 3x3 neighbourhood is a candidate; candidates are kept greedily,
best first, so that no two overlap in space and time. For each one the script prints the counts at several
magnitudes, the completeness, location quality and the number of events followed by another within 0.5 / 1 / 2 s
and 10 km (overlapping events, the target regime).

Selecting evaluation sequences by these aggregates is allowed (protocol Section 2): it reads no waveform and no
per-event result of any method. The script output is archived with the data card of every sequence chosen.

  python tools/data/scan_catalogue.py --days 7 --top 25 ~/rawdata/ncedc_cat/20*.ehpcsv > candidates.txt
"""
import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from ttmt.qc import catalogue_qc  # noqa: E402


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def read(paths, types=("eq",)):
    ev = []
    for p in paths:
        for r in csv.DictReader(open(p, newline="", encoding="utf-8", errors="replace")):
            if types and r.get("type") not in types:
                continue
            M = f(r.get("mag"))
            if M is None:
                continue
            t = datetime.fromisoformat(r["time"].replace("Z", "+00:00")).timestamp()
            ev.append({"t": t, "lat": float(r["latitude"]), "lon": float(r["longitude"]), "dep": f(r["depth"]),
                       "M": M, "mag_type": r.get("magType"), "nst": f(r.get("nst")), "gap": f(r.get("gap")),
                       "rms": f(r.get("rms")), "erh": f(r.get("horizontalError")), "erz": f(r.get("depthError")),
                       "status": r.get("status"), "place": r.get("place", "")})
    ev.sort(key=lambda e: e["t"])
    return ev


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--cell", type=float, default=0.2)
    ap.add_argument("--days", type=float, default=7.0)
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--mmin", type=float, default=1.0, help="magnitude for the target count")
    ap.add_argument("--json", help="write the candidates with their full QC here")
    a = ap.parse_args()
    ev = read(a.files)
    t = np.array([e["t"] for e in ev])
    ci = np.array([int(np.floor(e["lat"] / a.cell)) for e in ev])
    cj = np.array([int(np.floor(e["lon"] / a.cell)) for e in ev])
    W = a.days * 86400.0
    cells = defaultdict(list)
    for k, key in enumerate(zip(ci, cj)):
        cells[key].append(k)
    cand = []
    for (i, j) in cells:
        idx = np.array(sorted(k for di in (-1, 0, 1) for dj in (-1, 0, 1) for k in cells.get((i + di, j + dj), [])))
        tt = t[idx]
        if tt.size < 200:
            continue
        hi = np.searchsorted(tt, tt + W)
        n = hi - np.arange(tt.size)
        b = int(np.argmax(n))
        cand.append((int(n[b]), (i, j), float(tt[b]), float(tt[b] + W)))
    cand.sort(reverse=True)
    kept = []
    for n, (i, j), t0, t1 in cand:
        if any(abs(i - k[1][0]) <= 2 and abs(j - k[1][1]) <= 2 and t0 < k[3] and k[2] < t1 for k in kept):
            continue
        kept.append((n, (i, j), t0, t1))
        if len(kept) >= a.top:
            break
    rows = []
    for n, (i, j), t0, t1 in kept:
        sel = [e for e in ev if t0 <= e["t"] < t1 and abs(int(np.floor(e["lat"] / a.cell)) - i) <= 1
               and abs(int(np.floor(e["lon"] / a.cell)) - j) <= 1]
        rel = [dict(e, t=e["t"] - t0) for e in sel]
        q = catalogue_qc(rel)
        mc = q["mc_maxc"]
        place = Counter(e["place"].split(",")[0] for e in sel).most_common(1)[0][0]
        row = {"start": datetime.fromtimestamp(t0, timezone.utc).isoformat()[:19],
               "days": a.days, "lat": round((i + 1.5) * a.cell - a.cell, 2), "lon": round((j + 1.5) * a.cell - a.cell, 2),
               "place": place, "n": len(sel),
               "n_M1": sum(e["M"] >= 1.0 for e in sel), "n_M1.5": sum(e["M"] >= 1.5 for e in sel),
               "n_Mc": sum(e["M"] >= mc for e in sel), "Mmax": q["M_max"], "mc": mc, "mc_gft90": q["mc_gft90"],
               "erh_p50": q["erh_km"]["p50"], "gap_p50": q["gap_deg"]["p50"],
               "pairs_1s_10km": q["next_within_s_10km"][1.0], "pairs_2s_10km": q["next_within_s_10km"][2.0],
               "qc": q}
        rows.append(row)
    hdr = ("start", "lat", "lon", "place", "n", "n_M1", "n_M1.5", "n_Mc", "Mmax", "mc", "mc_gft90", "erh_p50",
           "gap_p50", "pairs_1s_10km", "pairs_2s_10km")
    print("\t".join(hdr))
    for r in rows:
        print("\t".join(f"{r[h]:.2f}" if isinstance(r[h], float) else str(r[h]) for h in hdr))
    if a.json:
        Path(a.json).write_text(json.dumps(rows, indent=1, default=str))


if __name__ == "__main__":
    main()
