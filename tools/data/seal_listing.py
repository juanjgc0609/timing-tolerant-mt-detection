"""Record the archive listing (key, size, ETag) of every continuous-waveform file of a reserved sequence: the
stations of the bucket inside the sealed disc, for every day of the window. Reads listings and StationXML only.

The listing pins the bytes without downloading them (S3 ETags are content hashes), so a later fetch can show it got
the files that existed when the sequence was sealed. Output: data/reserved/<id>.listing.tsv and its SHA-256.

  python tools/data/seal_listing.py conf_geysers_2024
"""
import hashlib
import math
import re
import sys
import urllib.parse
import urllib.request
from datetime import datetime, timedelta, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ttmt.registry import RESERVATION, check_read, load_reservation  # noqa: E402
from s3_availability import BUCKETS, coords, stationxml  # noqa: E402


def listing(base, prefix):
    out, token = [], None
    while True:
        url = f"{base}/?list-type=2&prefix={prefix}" + (f"&continuation-token={urllib.parse.quote(token)}" if token else "")
        body = urllib.request.urlopen(url, timeout=60).read().decode()
        for c in re.findall(r"<Contents>(.*?)</Contents>", body, re.S):
            k = re.search(r"<Key>([^<]+)", c).group(1)
            s = int(re.search(r"<Size>(\d+)", c).group(1))
            e = re.search(r"<ETag>([^<]+)", c).group(1).replace("&quot;", "").strip('"')
            out.append((k, s, e))
        m = re.search(r"<NextContinuationToken>([^<]+)</NextContinuationToken>", body)
        if not m:
            return out
        token = m.group(1)


def main():
    sid = sys.argv[1]
    doc, seqs = load_reservation()
    spec = next(s for s in doc["sequences"] if s["id"] == sid)
    seq = next(s for s in seqs if s.id == sid)
    cache = Path.home() / "rawdata" / "stationxml"
    rows, no_xml = [], set()
    day = datetime.fromtimestamp(seq.start, timezone.utc).date()
    while datetime(day.year, day.month, day.day, tzinfo=timezone.utc).timestamp() < seq.end:
        doy = day.timetuple().tm_yday
        for net in spec["networks"]:
            for k, s, e in listing(BUCKETS["ncedc"], f"continuous_waveforms/{net}/{day.year}/{day.year}.{doy:03d}/"):
                sta = Path(k).name.split(".")[0]
                xml = stationxml("ncedc", net, sta, cache)
                c = coords(xml) if xml else None
                if c is None:
                    no_xml.add(f"{net}.{sta}")
                    lat = lon = None
                else:
                    lat, lon = c
                d = None if c is None else \
                    111.19 * math.hypot(lat - seq.lat, (lon - seq.lon) * math.cos(math.radians(seq.lat)))
                if d is None or d <= seq.radius_km:
                    t0 = datetime(day.year, day.month, day.day, tzinfo=timezone.utc).timestamp()
                    check_read(net, lat, lon, t0, t0 + 86400, purpose="qc-metadata")
                    rows.append((k, s, e, "" if d is None else f"{d:.1f}"))
        day += timedelta(days=1)
    out = RESERVATION.parent / f"{sid}.listing.tsv"
    text = "key\tsize\tetag\tdist_km\n" + "".join(f"{k}\t{s}\t{e}\t{d}\n" for k, s, e, d in sorted(rows))
    out.write_text(text)
    print(f"{sid}: {len(rows)} files, {sum(r[1] for r in rows) / 1e9:.1f} GB, stations without StationXML in the "
          f"bucket (kept, distance unknown): {len(no_xml)}; sha256 {hashlib.sha256(text.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
