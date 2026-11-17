"""JMA arrival-time (phase-reading) files for 2012 from NIED Hi-net, for independent station terms of R1 (gap (a)).

Run on the laptop ONLY after the NIED registration-completion e-mail has arrived. Credentials come from the
environment (HINET_USER, HINET_PASS); never write them in a file, a commit or the chat.

python -I hinet_jma_fetch.py --out C:/Users/juanj/ttmt/hinet_jma [--year 2012] [--weeks N]   (--weeks 1: size probe)

Uses HinetPy (pip install HinetPy==0.9.1): Client(user, password).get_arrivaltime(startdate, span) with span <= 7 days
(checked in HinetPy 0.9.1, client.py: CatalogClient.get_arrivaltime -> data="measure"). Saves one raw file per week
under raw/, untouched, plus MANIFEST.tsv (file, bytes, sha256). NIED terms: raw files are not redistributed [to
verify]; only derived station terms and this script go to the repository.

Selection is done afterwards by a separate parser, written once the real format has been seen: events within 50 km of
the R1 centre (37.792 N, 139.99975 E), with readings at ATKH, INWH, NAZH, ONIH, THTH, TSTH or YNZH, and NOT on
2012-09-02 (the day of the evaluated 33-min record, whole day excluded). Criterion fixed before reading: at least 30
events per station; terms are the median P and S residual per station."""
import argparse
import hashlib
import os
import sys
import time
from datetime import date, timedelta


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--year", type=int, default=2012)
    ap.add_argument("--weeks", type=int, default=0, help="stop after N weeks (0: whole year)")
    a = ap.parse_args()
    user, pw = os.environ.get("HINET_USER"), os.environ.get("HINET_PASS")
    if not user or not pw:
        sys.exit("set HINET_USER and HINET_PASS in the environment")
    from HinetPy import Client
    client = Client(user, pw)
    raw = os.path.join(a.out, "raw")
    os.makedirs(raw, exist_ok=True)
    d, end = date(a.year, 1, 1), date(a.year + 1, 1, 1)
    rows = []
    while d < end and (not a.weeks or len(rows) < a.weeks):
        span = min(7, (end - d).days)
        name = os.path.join(raw, f"measure_{d:%Y%m%d}_{span}.txt")
        if not (os.path.exists(name) and os.path.getsize(name) > 0):
            for k in range(4):
                try:
                    client.get_arrivaltime(d, span, filename=name)
                    break
                except Exception as e:  # noqa: BLE001
                    print(f"  try {k + 1} {d}: {e}", file=sys.stderr)
                    time.sleep(10 * (k + 1))
        if os.path.exists(name):
            rows.append((os.path.basename(name), os.path.getsize(name), sha(name)))
            print(f"ok {rows[-1][0]} {rows[-1][1]} B")
        d += timedelta(days=span)
    with open(os.path.join(a.out, "MANIFEST.tsv"), "w", encoding="utf-8") as f:
        f.write("file\tbytes\tsha256\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    print(f"{len(rows)} files, {sum(r[1] for r in rows) / 1e6:.0f} MB")


if __name__ == "__main__":
    main()
