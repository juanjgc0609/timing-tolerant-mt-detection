"""Criterion 4 for the second-sequence candidates: stations with a usable 3-C set on the sampled day (day 8 of the
window). Metadata (FDSN station, text, level=channel) and MUSTANG percent_availability only; no waveforms.
Usage: python -I criterion4_stations.py"""
import csv, hashlib, io, json, math, os, sys, time, urllib.error, urllib.request
from collections import defaultdict
from datetime import datetime, timedelta, timezone

ROOT = r"C:\Users\juanj\ttmt\real_meta\candidates"
STATION = ["https://service.earthscope.org/fdsnws/station/1/query", "https://service.iris.edu/fdsnws/station/1/query"]
MUSTANG = "https://service.earthscope.org/mustang/measurements/1/query"
BANDS = ("HH", "BH", "EH", "HN", "SH")
R_KM, PA_MIN = 50.0, 90.0
CANDS = [  # rank, catalogue, window start (UTC), centre lat, lon
    (1, "hv", "2024-01-19T19:31:16", 19.30, -155.30),
    (2, "uu", "2020-03-18T13:09:31", 40.70, -112.10),
    (3, "hv", "2021-03-18T14:24:27", 19.30, -155.70),
    (4, "hv", "2024-06-27T10:00:23", 19.50, -155.30),
    (5, "uw", "2025-07-08T02:57:36", 46.90, -121.90),
    (6, "uu", "2021-06-21T14:33:47", 44.50, -110.50),
]


def fetch(url, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "ttmt-data-fetch/1.0"})
            with urllib.request.urlopen(req, timeout=300) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code in (204, 404):
                return b""
            print(f"  try {k + 1}: {url} -> HTTP {e.code}", file=sys.stderr)
        except Exception as e:  # noqa: BLE001
            print(f"  try {k + 1}: {url} -> {e}", file=sys.stderr)
        time.sleep(5 * (k + 1))
    raise RuntimeError(f"FAILED {url}")


def km(a, b, c, d):
    p1, p2 = math.radians(a), math.radians(c)
    h = math.sin((p2 - p1) / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(math.radians(d - b) / 2) ** 2
    return 2 * 6371.0 * math.asin(math.sqrt(h))


def ts(s):
    if not s or not s.strip():
        return None
    s = s.strip().rstrip("Z")
    for fmt in ("%Y-%m-%dT%H:%M:%S.%f", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d"):
        try:
            return datetime.strptime(s, fmt).replace(tzinfo=timezone.utc)
        except ValueError:
            pass
    raise ValueError(s)


def covers(epochs, t0, t1):
    iv = sorted((max(s, t0), min(e or t1, t1)) for s, e in epochs if s < t1 and (e is None or e > t0))
    cur = t0
    for s, e in iv:
        if s > cur:
            return False
        cur = max(cur, e)
    return cur >= t1


def save(path, data):
    with open(path, "wb") as f:
        f.write(data)
    return hashlib.sha256(data).hexdigest()


report = {}
for rank, cat, start, lat0, lon0 in CANDS:
    d0 = datetime.fromisoformat(start).replace(tzinfo=timezone.utc)
    day0 = datetime(d0.year, d0.month, d0.day, tzinfo=timezone.utc) + timedelta(days=7)
    day1 = day0 + timedelta(days=1)
    folder = os.path.join(ROOT, f"{rank}_{cat}_{start[:10]}")
    os.makedirs(folder, exist_ok=True)
    q = (f"?network=*&level=channel&latitude={lat0}&longitude={lon0}&maxradius=0.5"
         f"&starttime={day0:%Y-%m-%dT%H:%M:%S}&endtime={day1:%Y-%m-%dT%H:%M:%S}&format=text&nodata=404")
    body, used = None, None
    for s in STATION:
        try:
            body, used = fetch(s + q), s + q
            break
        except RuntimeError as e:
            print(e, file=sys.stderr)
    sha_st = save(os.path.join(folder, "stations_channels.txt"), body)
    files = {"stations_channels.txt": (sha_st, used)}
    stas, ch = {}, defaultdict(list)
    for line in body.decode("utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.split("|")
        net, sta, loc, cha = p[0], p[1], p[2], p[3]
        stas[(net, sta)] = (float(p[4]), float(p[5]))
        ch[(net, sta)].append((loc, cha, ts(p[15]), ts(p[16]) if len(p) > 16 else None))
    # candidate 3-C sets from metadata
    cand3c, fails = {}, {}
    for (net, sta), (la, lo) in stas.items():
        d = km(lat0, lon0, la, lo)
        if d > R_KM:
            fails[f"{net}.{sta}"] = (round(d, 2), "beyond 50 km")
            continue
        if net == "SY":
            fails[f"{net}.{sta}"] = (round(d, 2), "network SY (synthetics)")
            continue
        groups = defaultdict(lambda: defaultdict(list))
        for loc, cha, s, e in ch[(net, sta)]:
            if cha[:2] in BANDS:
                groups[(loc, cha[:2])][cha[2]].append((s, e))
        sets = []
        for (loc, band), comps in groups.items():
            okc = [c for c, ep in comps.items() if covers(ep, day0, day1)]
            hor = sorted(c for c in okc if c in "NE12")
            if "Z" in okc and len(hor) >= 2:
                sets.append((loc, band, ["Z"] + hor[:2]))
        if not sets:
            bands = sorted({f"{l}.{c}" for l, c, _, _ in ch[(net, sta)]})
            fails[f"{net}.{sta}"] = (round(d, 2), f"no 3-C set in bands {BANDS} covering the day; channels: {' '.join(bands)}")
            continue
        cand3c[(net, sta)] = (round(d, 2), sets)
    # MUSTANG percent_availability for the sampled day, per network
    pa = {}
    bynet = defaultdict(list)
    for net, sta in cand3c:
        bynet[net].append(sta)
    for net, sl in sorted(bynet.items()):
        url = (f"{MUSTANG}?metric=percent_availability&net={net}&sta={','.join(sorted(sl))}"
               f"&cha={','.join(b + '?' for b in BANDS)}&format=text&timewindow={day0:%Y-%m-%dT%H:%M:%S},{day1:%Y-%m-%dT%H:%M:%S}")
        raw = fetch(url)
        files[f"mustang_{net}.txt"] = (save(os.path.join(folder, f"mustang_{net}.txt"), raw), url)
        for line in raw.decode("utf-8").splitlines():
            if not line.startswith('"') or line.startswith('"value"') or line.count(",") < 3:
                continue
            v, target, st, en, _ = next(csv.reader(io.StringIO(line)))
            if st[:10].replace("/", "-") == f"{day0:%Y-%m-%d}":
                pa[target] = float(v)
    counted = []
    for (net, sta), (d, sets) in sorted(cand3c.items(), key=lambda x: x[1][0]):
        ok_sets, why = [], []
        for loc, band, comps in sets:
            vals = [pa.get(f"{net}.{sta}.{loc}.{band}{c}.M") for c in comps]
            if all(v is not None and v >= PA_MIN for v in vals):
                ok_sets.append((f"{loc}.{band}", vals))
            else:
                why.append(f"{loc}.{band} availability {vals}")
        if ok_sets:
            counted.append((f"{net}.{sta}", d, ok_sets))
        else:
            fails[f"{net}.{sta}"] = (d, "; ".join(why))
    report[f"{rank}_{cat}_{start[:10]}"] = {
        "centre": [lat0, lon0], "sampled_day": f"{day0:%Y-%m-%d}", "count": len(counted),
        "counted": [{"station": s, "dist_km": d, "sets": [{"set": k, "pa": v} for k, v in os_]} for s, d, os_ in counted],
        "fails": {k: {"dist_km": v[0], "reason": v[1]} for k, v in sorted(fails.items(), key=lambda x: x[1][0]) if v[0] <= R_KM},
        "n_stations_beyond_50km": sum(1 for v in fails.values() if v[0] > R_KM),
        "files": {k: {"sha256": a, "url": b} for k, (a, b) in files.items()},
    }
    with open(os.path.join(folder, "criterion4.json"), "w", encoding="utf-8", newline="\n") as f:
        json.dump(report[f"{rank}_{cat}_{start[:10]}"], f, indent=1)
    r = report[f"{rank}_{cat}_{start[:10]}"]
    print(f"== {rank} {cat} {start} centre {lat0},{lon0} day {r['sampled_day']}: COUNT {r['count']}")
    for c in r["counted"]:
        print(f"  + {c['station']} {c['dist_km']} km " + " ".join(f"{s['set']}{s['pa']}" for s in c["sets"]))
    for k, v in r["fails"].items():
        print(f"  - {k} {v['dist_km']} km: {v['reason']}")
    print(f"  ({r['n_stations_beyond_50km']} stations returned beyond 50 km not listed)")
    for k, v in r["files"].items():
        print(f"  file {k} sha256 {v['sha256']}")
