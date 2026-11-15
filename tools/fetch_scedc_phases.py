"""Fetch SCSN analyst phase picks of the Ridgecrest 2019 sequence from the SCEDC public bucket (s3://scedc-pds,
event_phases/YYYY/YYYY_DDD/<evid>.phase, anonymous HTTPS) and write a compact pick table for M2 (station terms in a
Gaussian process, plan T4.8 / spec M2).

    python tools/fetch_scedc_phases.py OUT_DIR [DAY0 DAY1]       (days of year 2019, inclusive; default 190 196)

Selection (fixed here, declared in experiments/exploratory/M2-GP/card.md):
  events: SCSN type "eq", epicentre in the box 35.45-36.15 N, 117.95-117.25 W, magnitude >= 1.5;
  picks:  every P and S line of those events (network, station, phase, weight, epicentral distance, travel time).
Days 190-196 (2019-07-09 to 07-15) are after the R0 hour (2019-07-04 19:00-20:00 UTC, day 185), which is spent
(CR-23) and is not touched. Citation: SCEDC (2013), Southern California Earthquake Data Center, Caltech,
doi:10.7909/C3WD3xH1.

Outputs: OUT_DIR/events.csv (evid, origin UTC, lat, lon, depth_km, mag), OUT_DIR/picks.csv (evid, net, sta, phase,
lat, lon, elev_m, weight, dist_km, tt_s), OUT_DIR/raw/*.phase (the source files), OUT_DIR/SHA256SUMS. Reads untrusted
text; nothing is executed from it.
"""
import concurrent.futures as cf
import csv
import hashlib
import os
import re
import sys
import time
import urllib.parse
import urllib.request

BUCKET = "https://scedc-pds.s3.amazonaws.com"
BOX = (35.45, 36.15, -117.95, -117.25)
MMIN = 1.5


def fetch(url, tries=5):
    """GET with retries (connection resets happen on long runs); the cache in raw/ makes reruns resume."""
    for k in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=60).read()
        except OSError:
            if k == tries - 1:
                raise
            time.sleep(2 ** (k + 1))


def list_day(day):
    keys, tok = [], None
    while True:
        u = f"{BUCKET}/?list-type=2&prefix=event_phases/2019/2019_{day:03d}/"
        if tok:
            u += "&continuation-token=" + urllib.parse.quote(tok)
        x = fetch(u).decode()
        keys += re.findall(r"<Key>([^<]+\.phase)</Key>", x)
        m = re.search(r"<NextContinuationToken>([^<]+)<", x)
        if not m:
            return keys
        tok = m.group(1)


def get(key, raw):
    p = os.path.join(raw, os.path.basename(key))
    if not os.path.exists(p):
        b = fetch(f"{BUCKET}/{key}")
        with open(p + ".tmp", "wb") as f:
            f.write(b)
        os.replace(p + ".tmp", p)
    return p


def parse(path):
    lines = open(path, encoding="ascii", errors="replace").read().splitlines()
    h = lines[0].split()
    ev = {"evid": h[0], "type": h[1], "origin": h[3].replace("/", "-").replace(",", "T"),
          "lat": float(h[4]), "lon": float(h[5]), "depth_km": float(h[6]), "mag": float(h[7])}
    picks = []
    for ln in lines[1:]:
        f = ln.split()
        if len(f) < 13 or f[7] not in ("P", "S"):
            continue
        picks.append({"evid": ev["evid"], "net": f[0], "sta": f[1], "phase": f[7], "lat": float(f[4]),
                      "lon": float(f[5]), "elev_m": float(f[6]), "weight": float(f[-3]), "dist_km": float(f[-2]),
                      "tt_s": float(f[-1])})
    return ev, picks


def main(out, d0=190, d1=196):
    raw = os.path.join(out, "raw")
    os.makedirs(raw, exist_ok=True)
    keys = [k for d in range(d0, d1 + 1) for k in list_day(d)]
    with cf.ThreadPoolExecutor(4) as ex:
        paths = sorted(ex.map(lambda k: get(k, raw), keys))
    evs, pks = [], []
    for p in paths:
        ev, picks = parse(p)
        if ev["type"] != "eq" or ev["mag"] < MMIN:
            continue
        if not (BOX[0] <= ev["lat"] <= BOX[1] and BOX[2] <= ev["lon"] <= BOX[3]):
            continue
        evs.append(ev)
        pks += picks
    evs.sort(key=lambda e: e["origin"])
    for name, rows, cols in (("events.csv", evs, ["evid", "origin", "lat", "lon", "depth_km", "mag"]),
                             ("picks.csv", pks, ["evid", "net", "sta", "phase", "lat", "lon", "elev_m", "weight",
                                                 "dist_km", "tt_s"])):
        with open(os.path.join(out, name), "w", newline="") as f:
            w = csv.DictWriter(f, cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
    with open(os.path.join(out, "SHA256SUMS"), "w") as f:
        for name in ("events.csv", "picks.csv"):
            f.write(f"{hashlib.sha256(open(os.path.join(out, name), 'rb').read()).hexdigest()}  {name}\n")
    print(f"{len(keys)} phase files, {len(evs)} events kept, {len(pks)} picks")


if __name__ == "__main__":
    a = sys.argv[1:]
    sys.exit(main(a[0], *(int(x) for x in a[1:3])))
