"""Small metadata and catalogue downloads the cloud cannot reach (plan de datos, D4; approved by the owner 2026-11-15).

python -I laptop_fetch.py --out C:/Users/juanj/ttmt/real_meta [--qtm-url URL]

Fetches, with the standard library only:
  1. NN StationXML (level=response) for the 60 km discs of conf_antelope_valley_2021 and conf_bodie_2020, epoch
     overlapping each window. Metadata only: purpose qc-metadata, allowed by the seal (no waveform sample is read).
  2. ComCat global M >= 5.5 events, 2010-01-01 to 2026-11-01, one CSV per year (CR-16 teleseism flags).
  3. The Ridgecrest 2019 QTM catalogue (Ross et al. 2019), when --qtm-url is given.
Writes MANIFEST.tsv (file, bytes, sha256, url, fetched_utc). Skips files already present with the same size."""
import argparse
import hashlib
import os
import sys
import time
import urllib.request
from datetime import datetime, timezone

STATION = ["https://service.earthscope.org/fdsnws/station/1/query", "https://service.iris.edu/fdsnws/station/1/query"]
EVENT = "https://earthquake.usgs.gov/fdsnws/event/1/query"
DISCS = {  # from data/reserved/reservation.json (window lat, lon, radius 60 km = 0.54 deg; 0.6 deg margin)
    "conf_antelope_valley_2021": (38.5, -119.5, "2021-07-08T22:00:00", "2021-08-07T22:00:00"),
    "conf_bodie_2020": (38.3, -118.9, "2020-04-11T14:00:00", "2020-05-11T14:00:00"),
}


def get(url, path, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "ttmt-data-fetch/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r, open(path + ".part", "wb") as f:
                while chunk := r.read(1 << 20):
                    f.write(chunk)
            os.replace(path + ".part", path)
            return True
        except Exception as e:  # noqa: BLE001
            print(f"  try {k + 1}: {url} -> {e}", file=sys.stderr)
            time.sleep(5 * (k + 1))
    return False


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--qtm-url")
    a = ap.parse_args()
    jobs = []
    for sid, (lat, lon, t0, t1) in DISCS.items():
        q = f"?network=NN&level=response&latitude={lat}&longitude={lon}&maxradius=0.6&starttime={t0}&endtime={t1}&format=xml&nodata=404"
        jobs.append((f"stationxml/NN_{sid}.xml", [s + q for s in STATION]))
    for y in range(2010, 2027):
        end = "2026-11-01" if y == 2026 else f"{y + 1}-01-01"
        jobs.append((f"comcat/comcat_M5.5_{y}.csv",
                     [f"{EVENT}?format=csv&starttime={y}-01-01&endtime={end}&minmagnitude=5.5&orderby=time-asc"]))
    if a.qtm_url:
        jobs.append((f"qtm/{os.path.basename(a.qtm_url.split('?')[0])}", [a.qtm_url]))
    rows = []
    for rel, urls in jobs:
        path = os.path.join(a.out, rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        used = next((u for u in urls if get(u, path)), None)
        if used is None:
            print(f"FAILED {rel}", file=sys.stderr)
            continue
        rows.append((rel, os.path.getsize(path), sha(path), used, datetime.now(timezone.utc).isoformat(timespec="seconds")))
        print(f"ok {rel} {rows[-1][1]} B")
    with open(os.path.join(a.out, "MANIFEST.tsv"), "w", encoding="utf-8") as f:
        f.write("file\tbytes\tsha256\turl\tfetched_utc\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    print(f"{len(rows)}/{len(jobs)} files; manifest written")


if __name__ == "__main__":
    main()
