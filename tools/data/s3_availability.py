"""Which stations of the NCEDC / SCEDC open buckets have continuous three-component data and StationXML near a
point on a given day. Reads S3 listings and StationXML only (no waveform), so it may be run for sealed candidates.

  python tools/data/s3_availability.py --day 2014-09-24 --lat 37.62 --lon -118.88 --radius 60 \
      --nets NC BK BG BP NP NN PB --cache ~/rawdata/stationxml
"""
import argparse
import json
import math
import re
import sys
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import date
from pathlib import Path

BUCKETS = {"ncedc": "https://ncedc-pds.s3.amazonaws.com", "scedc": "https://scedc-pds.s3.amazonaws.com"}


def list_keys(base, prefix):
    keys, token = [], None
    while True:
        url = f"{base}/?list-type=2&prefix={prefix}" + (f"&continuation-token={urllib.parse.quote(token)}" if token else "")
        body = urllib.request.urlopen(url, timeout=60).read().decode()
        keys += [(k, int(s)) for k, s in re.findall(r"<Key>([^<]+)</Key>.*?<Size>(\d+)</Size>", body, re.S)]
        m = re.search(r"<NextContinuationToken>([^<]+)</NextContinuationToken>", body)
        if not m:
            return keys
        token = m.group(1)


def day_files(bucket, net, day):
    """{(sta, loc, cha): size} of one network-day."""
    doy = day.timetuple().tm_yday
    out = {}
    if bucket == "ncedc":
        for k, s in list_keys(BUCKETS[bucket], f"continuous_waveforms/{net}/{day.year}/{day.year}.{doy:03d}/"):
            sta, n, cha, loc = Path(k).name.split(".")[:4]
            out[(sta, loc, cha)] = (s, k)
    else:  # scedc: continuous_waveforms/YYYY/YYYY_DOY/NNSSSSSCCCLL_YYYYDOY.ms
        for k, s in list_keys(BUCKETS[bucket], f"continuous_waveforms/{day.year}/{day.year}_{doy:03d}/{net}"):
            nm = Path(k).name
            sta, cha, loc = nm[2:7].rstrip("_"), nm[7:10], nm[10:12].rstrip("_")
            out[(sta, loc, cha)] = (s, k)
    return out


def stationxml(bucket, net, sta, cache):
    p = Path(cache) / bucket / f"{net}.{sta}.xml"
    if not p.exists():
        p.parent.mkdir(parents=True, exist_ok=True)
        url = (f"{BUCKETS[bucket]}/FDSNstationXML/{net}/{net}.{sta}.xml" if bucket == "ncedc"
               else f"{BUCKETS[bucket]}/FDSNstationXML/{net}/{net}_{sta}.xml")
        try:
            p.write_bytes(urllib.request.urlopen(url, timeout=60).read())
        except Exception:
            p.write_bytes(b"")
    return p if p.stat().st_size else None


def coords(xml):
    s = xml.read_text(errors="replace")
    lat = re.search(r"<Latitude[^>]*>([-\d.]+)<", s)
    lon = re.search(r"<Longitude[^>]*>([-\d.]+)<", s)
    return (float(lat.group(1)), float(lon.group(1))) if lat and lon else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bucket", default="ncedc", choices=BUCKETS)
    ap.add_argument("--day", required=True)
    ap.add_argument("--lat", type=float, required=True)
    ap.add_argument("--lon", type=float, required=True)
    ap.add_argument("--radius", type=float, default=60.0)
    ap.add_argument("--nets", nargs="+", default=["NC", "BK", "BG", "BP", "NP", "NN", "PB"])
    ap.add_argument("--cache", default=str(Path.home() / "rawdata" / "stationxml"))
    ap.add_argument("--json")
    a = ap.parse_args()
    day = date.fromisoformat(a.day)
    rows, no_xml = [], []
    for net in a.nets:
        files = day_files(a.bucket, net, day)
        by_sta = defaultdict(dict)
        for (sta, loc, cha), (size, key) in files.items():
            by_sta[sta][(loc, cha)] = size
        for sta, chans in sorted(by_sta.items()):
            xml = stationxml(a.bucket, net, sta, a.cache)
            c = coords(xml) if xml else None
            if c is None:
                no_xml.append(f"{net}.{sta}")
                continue
            d = 111.19 * math.hypot(c[0] - a.lat, (c[1] - a.lon) * math.cos(math.radians(a.lat)))
            if d > a.radius:
                continue
            groups = defaultdict(set)
            for loc, cha in chans:
                groups[(loc, cha[:2])].add(cha[2])
            three = sorted(f"{loc}.{bi}" for (loc, bi), cs in groups.items() if "Z" in cs and len(cs & set("NE12")) >= 2)
            rows.append({"nslc_station": f"{net}.{sta}", "dist_km": round(d, 1), "lat": c[0], "lon": c[1],
                         "three_component": three, "channels": sorted(f"{l}.{ch}" for l, ch in chans),
                         "mbytes": round(sum(chans.values()) / 1e6, 1)})
    rows.sort(key=lambda r: r["dist_km"])
    n3 = sum(bool(r["three_component"]) for r in rows)
    print(f"# {a.bucket} {a.day} within {a.radius} km of ({a.lat}, {a.lon}): {len(rows)} stations, {n3} with 3-C;"
          f" {len(no_xml)} stations of {a.nets} without StationXML in the bucket")
    for r in rows:
        print(f"{r['nslc_station']:10s} {r['dist_km']:6.1f} km  3C={','.join(r['three_component']) or '-':22s} "
              f"{r['mbytes']:7.1f} MB  {' '.join(r['channels'])}")
    if a.json:
        Path(a.json).write_text(json.dumps({"rows": rows, "no_stationxml": no_xml}, indent=1))


if __name__ == "__main__":
    sys.exit(main())
