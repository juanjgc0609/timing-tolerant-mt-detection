"""CR-04 prevalence table of one sequence (docs/protocol/CR04_PREVALENCE.md). Runs on the author's machine only.

python run.py --catalogue cat.csv --stations sta.csv [--coverage cov.csv] --model hk|jma2001|homogeneous:VP \
    --name ridgecrest_2019 --out prevalence_ridgecrest_2019.json

cat.csv: time (ISO 8601 UTC or epoch s), lat, lon, depth_km. sta.csv: code, lat, lon. cov.csv: code, start, end
(gap-free data intervals; without it every station counts as available and the table is flagged "epochs only").
Writes the table and its SHA-256 next to it."""
import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone

import numpy as np

from ttmt import traveltime as TT
from ttmt.atlas import prevalence as P

KM_PER_DEG = 111.195


def _t(s):
    try:
        return float(s)
    except ValueError:
        return datetime.fromisoformat(s.replace("Z", "+00:00")).astimezone(timezone.utc).timestamp()


def _rows(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _local(lat, lon, lat0, lon0):
    return np.c_[(np.asarray(lon) - lon0) * KM_PER_DEG * np.cos(np.radians(lat0)), (np.asarray(lat) - lat0) * KM_PER_DEG]


def _tt(model):
    if model.startswith("homogeneous:"):
        vp = float(model.split(":")[1])
        return lambda src, sta: np.stack([np.linalg.norm(sta - src, axis=1) / vp,
                                          np.linalg.norm(sta - src, axis=1) / (vp / TT.R_VPVS)])
    return TT.SourceTimes({"hk": TT.HK, "jma2001": TT.JMA2001}[model])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalogue", required=True)
    ap.add_argument("--stations", required=True)
    ap.add_argument("--coverage")
    ap.add_argument("--model", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=20261101)
    ap.add_argument("--n-boot", type=int, default=2000)
    a = ap.parse_args()
    cat, sta = _rows(a.catalogue), _rows(a.stations)
    t = np.array([_t(r["time"]) for r in cat])
    lat, lon = np.array([float(r["lat"]) for r in cat]), np.array([float(r["lon"]) for r in cat])
    lat0, lon0 = float(lat.mean()), float(lon.mean())
    src = np.c_[_local(lat, lon, lat0, lon0), [float(r["depth_km"]) for r in cat]]
    codes = [r["code"] for r in sta]
    sxyz = np.c_[_local([float(r["lat"]) for r in sta], [float(r["lon"]) for r in sta], lat0, lon0), np.zeros(len(sta))]
    if a.coverage:
        cov = {c: [] for c in codes}
        for r in _rows(a.coverage):
            if r["code"] in cov:
                cov[r["code"]].append((_t(r["start"]), _t(r["end"])))
        from ttmt.atlas.cells import availability
        avail = availability(t, [cov[c] for c in codes])
    else:
        avail = np.ones((len(t), len(codes)), bool)
    out = P.table(t, src, sxyz, codes, avail, _tt(a.model), n_boot=a.n_boot, seed=a.seed)
    out.update(name=a.name, model=a.model, availability="data intervals" if a.coverage else "epochs only",
               catalogue_sha256=hashlib.sha256(open(a.catalogue, "rb").read()).hexdigest(),
               note="Y is a lower bound: reviewed catalogues under-list pairs below about 1 s (CR04_PREVALENCE.md, 5)")
    out["k"] = {str(k): v for k, v in out["k"].items()}
    txt = json.dumps(out, indent=1, default=float)
    with open(a.out, "w", encoding="utf-8") as f:
        f.write(txt)
    with open(a.out + ".sha256", "w") as f:
        f.write(hashlib.sha256(txt.encode()).hexdigest() + "\n")
    for k, v in out["k"].items():
        print(f"k={k}: n={v['n_counted']} Y={v['Y']} CI={v['Y_ci95']} Poisson={v['Y_poisson']}")


if __name__ == "__main__":
    main()
