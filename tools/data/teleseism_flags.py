"""CR-16: windows with catalogued teleseismic arrivals in the development data used for time-slide nulls.

python -I tools/data/teleseism_flags.py OUT.tsv comcat_M5.5_*.csv      (ComCat CSV: time, latitude, longitude, depth, mag)

The rule is ttmt.qc.teleseism_windows (fixed 2026-11-15 before any flag was computed): M >= 5.5 anywhere, from the
iasp91 first P - 30 s to distance / 2.5 km/s + 600 s. Sites and spans are the development data of the atlas and of
the real-noise tests (2M noise days of Alto Tiberina, R0, R1); no reserved window is listed. Writes one row per
window: site, span, start and end UTC, duration s, number of events."""
import csv
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from ttmt.qc import teleseism_windows  # noqa: E402

SITES = {  # name: (lat, lon, start UTC, duration s)
    "altotiberina_101231": (43.315, 12.388, "2010-12-31T00:00:00", 86400.0),
    "altotiberina_100210": (43.315, 12.388, "2010-02-10T00:00:00", 86400.0),
    "altotiberina_100211": (43.315, 12.388, "2010-02-11T00:00:00", 86400.0),
    "altotiberina_100723": (43.315, 12.388, "2010-07-23T00:00:00", 86400.0),
    "ridgecrest_R0": (35.77, -117.60, "2019-07-04T19:00:00", 3560.0),
    "hinet_R1": (37.792, 139.99975, "2012-09-02T03:20:00", 1960.0),
}


def ep(s):
    return datetime.fromisoformat(s.replace("Z", "+00:00")).replace(tzinfo=timezone.utc).timestamp()


def iso(t):
    return datetime.fromtimestamp(t, timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def main():
    out, files = sys.argv[1], sys.argv[2:]
    ev = []
    for p in files:
        with open(p, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                try:
                    ev.append((ep(r["time"]), float(r["latitude"]), float(r["longitude"]),
                               float(r["depth"] or 0.0), float(r["mag"])))
                except (ValueError, KeyError):
                    continue
    rows = []
    for name, (lat, lon, start, dur) in SITES.items():
        t0 = ep(start)
        w = teleseism_windows(ev, lat, lon, t0, t0 + dur)
        flagged = sum(b - a for a, b, _ in w)
        print(f"{name}: {len(w)} windows, {flagged:.0f} s of {dur:.0f} s flagged ({100 * flagged / dur:.1f} %)")
        rows += [(name, start, iso(a), iso(b), f"{b - a:.0f}", n) for a, b, n in w]
    with open(out, "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f, delimiter="\t", lineterminator="\n")
        wr.writerow(["site", "span_start", "start", "end", "duration_s", "n_events"])
        wr.writerows(rows)
    print(f"{len(ev)} events read; {len(rows)} windows written to {out}")


if __name__ == "__main__":
    main()
