"""diag_subnets.py - re-derives the registered subnet rule from stations.csv: the 8 closest HH stations to the
pre-registration centroid (35.6773 N, 117.5418 W) and the narrowest azimuthal sector containing 8 HH stations."""
import csv, math
import numpy as np
rows = {r["station"]: r for r in csv.DictReader(open("/home/claude/realdata/skillsbench_eqassoc/stations.csv")) if r["channel"] == "HHZ"}
la0, lo0 = 35.6773, -117.5418
st = []
for n, r in rows.items():
    x = (float(r["longitude"]) - lo0) * 111.19 * math.cos(math.radians(la0)); y = (float(r["latitude"]) - la0) * 111.19
    st.append((math.hypot(x, y), math.degrees(math.atan2(x, y)) % 360, n))
st.sort()
print("HH stations by distance:", [(n, round(d, 1), round(a)) for d, a, n in st])
print("8 closest:", sorted(n for _, _, n in st[:8]))
az = sorted((a, n) for _, a, n in st)
best = None
for i in range(len(az)):
    sel = [az[(i + k) % len(az)] for k in range(8)]
    w = (sel[-1][0] - sel[0][0]) % 360
    if best is None or w < best[0]:
        best = (w, sel)
print("narrowest 8-sector width %.0f deg from %.0f to %.0f:" % (best[0], best[1][0][0], best[1][-1][0]), sorted(n for _, n in best[1]))
