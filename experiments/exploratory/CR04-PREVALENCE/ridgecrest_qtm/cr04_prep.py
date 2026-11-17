"""CR-04 prep: sanity check of the Ridgecrest QTM .cat, write qtm_cat.csv and sta_atlas.csv."""
import csv, hashlib, sys
from datetime import datetime, timedelta, timezone

cat, stfile, out_cat, out_sta = sys.argv[1:5]
CODES = "CLC SRT CCC TOW2 LRL WRC2 WBM SLA JRC2 WCS2 MPM".split()

rows = []
with open(cat, encoding="ascii") as f:
    for line in f:
        p = line.split()
        if not p:
            continue
        y, mo, d, h, mi = (int(x) for x in p[:5])
        sec = float(p[5])
        t = datetime(y, mo, d, h, mi, tzinfo=timezone.utc) + timedelta(seconds=sec)
        rows.append((t, float(p[7]), float(p[8]), float(p[9]), float(p[10])))
lat = [r[1] for r in rows]; lon = [r[2] for r in rows]; dep = [r[3] for r in rows]; c11 = [r[4] for r in rows]
print(f"rows {len(rows)}")
for n, v in (("lat", lat), ("lon", lon), ("depth", dep), ("col11", c11)):
    print(f"{n}: min {min(v)} max {max(v)}")
# Pre-run bounds. Col 11 was "roughly -1 to 7.1"; the observed -2.82..7.1 was accepted as a declared
# deviation (coordinator decision (a), 2026-11-15): the check only identifies the column, and magnitude
# is not used in the estimand. Bound below is the accepted observed range.
ok = (35 <= min(lat) and max(lat) <= 36.5 and -118.5 <= min(lon) and max(lon) <= -117
      and -2 <= min(dep) and max(dep) <= 35 and -2.82 <= min(c11) and max(c11) <= 7.1)
print("SANITY", "OK" if ok else "FAIL")
if not ok:
    sys.exit(2)
rows.sort(key=lambda r: r[0])
with open(out_cat, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["time", "lat", "lon", "depth_km"])
    for t, la, lo, de, _ in rows:
        w.writerow([t.strftime("%Y-%m-%dT%H:%M:%S.%fZ"), f"{la:.5f}", f"{lo:.5f}", f"{de:.3f}"])
h = hashlib.sha256(open(out_cat, "rb").read()).hexdigest()
print(f"qtm_cat.csv rows {len(rows)} sha256 {h} first {rows[0][0].isoformat()} last {rows[-1][0].isoformat()}")

with open(stfile, encoding="utf-8") as f:
    st = list(csv.DictReader(f))
print("stations.csv columns:", list(st[0].keys()))
found = {}
for r in st:
    code = r.get("station") or r.get("code")
    net = r.get("network", "")
    key = code.split(".")[1] if code.count(".") >= 1 else code
    if key in CODES and (not net or net == "CI") and (code.startswith("CI.") or "." not in code):
        v = (float(r["latitude"]), float(r["longitude"]), f"{net}.{code}")
        if key in found and found[key][:2] != v[:2]:
            print("COORD CONFLICT", key, found[key], v)
        found[key] = v
missing = [c for c in CODES if c not in found]
if missing:
    print("MISSING", missing); sys.exit(3)
with open(out_sta, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["code", "lat", "lon"])
    for c in CODES:
        w.writerow([c, found[c][0], found[c][1]])
        print(c, found[c][0], found[c][1], "src:", found[c][2])
