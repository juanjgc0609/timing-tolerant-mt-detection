"""The run plan of ATLAS-0: which records every method runs on, with which stations (card, "Regimes and strata",
items 3 and CR-24; ADR-0009 item 4).

A run is one realisation of a semi-synthetic set seen by one controlled station subset: k in cells.STATION_COUNTS,
the k stations nearest the station centroid of the network (cells.station_subset; the injected records have no gaps,
so every station is available), ties broken by code. A network with fewer than k stations contributes nothing to that
k. The 4-station subset is also run on the contamination variants of the generator (CR-24), whose scaled stations
are indices of the full station list drawn from that same subset. The distance axis uses the centroid of the
network's development catalogue (card item 4), which the set manifest carries as "catalogue_centroid".

The plan only lists runs: it reads no waveform. Every method runs every row of its partition; the atlas reads
"read" rows only, tuning reads "tune" rows only (ttmt.associate.tune)."""
import csv
import json
from pathlib import Path

import numpy as np

from . import cells as C

PLAN_FIELDS = ("run", "set_id", "dataset", "label", "partition", "block_index", "r", "k", "variant", "stations",
               "station_idx", "scaled", "factor")


def subsets(stations_xyz, codes):
    """{k: station indices} of the controlled subsets that the network can fill (cells.station_subset, all stations
    available)."""
    out = {}
    for k in C.STATION_COUNTS:
        s = C.station_subset(stations_xyz, codes, np.ones(len(codes), bool), k)
        if s is not None:
            out[k] = [int(i) for i in s]
    return out


def run_name(label, k, variant=""):
    return f"{label}__k{k}" + (f"__{variant}" if variant else "")


def plan(manifest, partition=None):
    """rows of the run plan of one set manifest (ttmt.synth.sets), optionally of one partition only."""
    codes, xyz = list(manifest["stations"]), np.asarray(manifest["xyz"], float)
    subs = subsets(xyz, codes)
    rows = []
    for r in manifest["realisations"]:
        if partition is not None and r["partition"] != partition:
            continue
        for k, idx in subs.items():
            base = dict(set_id=manifest["set_id"], dataset=manifest.get("dataset", manifest["set_id"]),
                        label=r["label"], partition=r["partition"],
                        block_index=int(r.get("block_index", -1)), r=int(r.get("r", -1)), k=k,
                        stations=" ".join(codes[i] for i in idx), station_idx=" ".join(map(str, idx)))
            rows.append(dict(base, run=run_name(r["label"], k), variant="", scaled="", factor=""))
            if k != 4:
                continue
            for name, scaled, factor in r.get("variants", []):
                if not set(scaled) <= set(idx):
                    raise ValueError(f"{r['label']}: variant {name} scales stations {scaled} outside the 4-station "
                                     f"subset {idx}")
                rows.append(dict(base, run=run_name(r["label"], k, name), variant=name,
                                 scaled=" ".join(map(str, scaled)), factor=float(factor)))
    return rows


def catalogue_centroid(manifest):
    """the distance-axis centroid (x, y) in local km; the set manifest must carry it."""
    c = manifest.get("catalogue_centroid")
    if c is None:
        raise KeyError(f"manifest of {manifest.get('set_id')!r} has no catalogue_centroid (card item 4: distance to "
                       "the centroid of the network's development catalogue)")
    return tuple(map(float, c[:2]))


def write(rows, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, PLAN_FIELDS, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return path


def summary(rows):
    """counts of runs per (partition, k, contamination)."""
    out = {}
    for r in rows:
        key = (r["partition"], r["k"], "variant" if r["variant"] else "base")
        out[key] = out.get(key, 0) + 1
    return dict(sorted(out.items()))


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="write the ATLAS-0 run plan of set manifests")
    ap.add_argument("manifests", nargs="+")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parents[3] / "experiments" / "exploratory" /
                                         "ATLAS-0" / "run_plan.csv"))
    a = ap.parse_args()
    rows = [row for m in a.manifests for row in plan(json.loads(Path(m).read_text()))]
    print(write(rows, a.out), summary(rows))
