"""
verif_analysis.py - independent re-implementation of the R0 endpoint (does not import r0_analyze.py or r0_lib.py).

    python verif_analysis.py M2DIR CLOSEDIR PICKSDIR OUT.json [B]

Catalogue read with csv/datetime; frame centre = mean lat/lon of catalogued events with origin in [0, 900) s;
local frame x = east, y = north (equirectangular, 111.19 km/deg, cos(lat0)) -- the frame in which the detectors and
the associators were given station coordinates; projection error vs. a geodesic is reported separately.
Rules: LOC+T epicentre <= 5 km and |dt| <= 0.5 s; LOC <= 5 km and <= 3.0 s.  A detection matching any M >= 2 event
is a hit (on the nearest such event); matching only M < 2 events is neutral; otherwise false.
Pick-free atoms: de-duplicated (descending score; atom removed if within 5 km AND 1.0 s of an already kept atom).
Recall at <= cap false: cap = floor(rate * hours), hours = 2660/3600 (also reported with the 45-block convention of
r0_analyze, 0.75 h).  A threshold keeps detections with score strictly above it; the best threshold is the one
just above the (cap+1)-th largest false score (ties at that score excluded).  Associators: max over config grid.
"""
import csv
import glob
import json
import math
import os
import sys
from datetime import datetime

import numpy as np

DATA = "/home/claude/realdata/skillsbench_eqassoc"
T0 = datetime(2019, 7, 4, 19, 0, 0)
LO, HI = 900.0, 3560.0
SUBS = ("N8", "N6", "N4", "U8")
TOL = {"LOC+T": 0.5, "LOC": 3.0}
RATES = (10, 30, 60)


def catalogue():
    rows = []
    for r in csv.DictReader(open(os.path.join(DATA, "catalog.csv"))):
        t = (datetime.fromisoformat(r["time"]) - T0).total_seconds()
        rows.append((t, float(r["magnitude"]), float(r["latitude"]), float(r["longitude"]), float(r["depth_km"])))
    cal = [r for r in rows if 0.0 <= r[0] < 900.0]
    lat0 = sum(r[2] for r in cal) / len(cal)
    lon0 = sum(r[3] for r in cal) / len(cal)
    k = 111.19
    out = np.array([[(r[3] - lon0) * k * math.cos(math.radians(lat0)), (r[2] - lat0) * k, r[0], r[1], r[4]]
                    for r in rows])
    return out, (lat0, lon0)


def dedupe(atoms, dkm=5.0, dts=1.0):
    a = sorted(atoms, key=lambda r: -r[4])
    kept = []
    K = np.zeros((0, 3))
    for r in a:
        if len(K):
            close = (np.hypot(K[:, 0] - r[0], K[:, 1] - r[1]) <= dkm) & (np.abs(K[:, 2] - r[3]) <= dts)
            if close.any():
                continue
        kept.append(r)
        K = np.vstack([K, [r[0], r[1], r[3]]])
    return kept


def labels(det, cat, tol):
    """returns per detection: target index (>=0) for a hit on M>=2, -2 for neutral, -1 for false"""
    out = np.full(len(det), -1, int)
    for i, d in enumerate(det):
        dist = np.hypot(cat[:, 0] - d[0], cat[:, 1] - d[1])
        ok = (dist <= 5.0) & (np.abs(cat[:, 2] - d[3]) <= tol + 1e-9)
        big = ok & (cat[:, 3] >= 2.0)
        if big.any():
            out[i] = int(np.where(big, dist, np.inf).argmin())
        elif ok.any():
            out[i] = -2
    return out


def recall_at(scores, lab, targets, cap, w_det=None, w_tar=None):
    """exact recall at <= cap false (weighted for the bootstrap)."""
    if len(scores) == 0:
        return 0.0
    w_det = np.ones(len(scores)) if w_det is None else w_det
    fal = lab == -1
    fs, fw = scores[fal], w_det[fal]
    o = np.argsort(-fs, kind="stable")
    cum = np.cumsum(fw[o])
    k = int(np.searchsorted(cum, cap + 1e-9, side="right"))
    thr = fs[o][k] if k < len(fs) else -np.inf
    sel = (lab >= 0) & (scores > thr) & (w_det > 0)
    got = set(lab[sel].tolist())
    if w_tar is None:
        return len(got & set(targets.tolist())) / len(targets)
    tw = dict(zip(targets.tolist(), w_tar))
    tot = sum(w_tar)
    return sum(tw[g] for g in got if g in tw) / tot if tot > 0 else 0.0


def prep(det, cat):
    det = [d for d in det if LO <= d[3] < HI]
    s = np.array([d[4] for d in det], float)
    t = np.array([d[3] for d in det], float)
    return {"s": s, "t": t, "lab": {u: labels(det, cat, tol) for u, tol in TOL.items()}, "n": len(det),
            "det": det}


def load(m2dir, closedir, pdir, cat):
    sets = {}
    for s in SUBS:
        for m, f in (("M2", os.path.join(m2dir, f"det_eval_{s}.json")),
                     ("M1", os.path.join(closedir, f"det_eval_{s}.json")),
                     ("M3", os.path.join(closedir, f"detM3_eval_{s}.json"))):
            if os.path.exists(f):
                r = json.load(open(f))
                if m in r["atoms"]:
                    sets[(s, m)] = [prep(dedupe(r["atoms"][m]), cat)]
        f = os.path.join(pdir, f"picks_eval_{s}.json")
        if os.path.exists(f):
            r = json.load(open(f))
            for nm, key in (("GaMMA", "gamma"), ("PyOcto", "pyocto"), ("GaMMA+SC", "gamma_sc"),
                            ("PyOcto+SC", "pyocto_sc")):
                sets[(s, nm)] = [prep(v, cat) for v in r[key].values()]
    return sets


def point(sets, targets, hours):
    out = {}
    for (s, m), ss in sets.items():
        for u in TOL:
            for rate in RATES:
                cap = int(math.floor(rate * hours + 1e-9))
                out[(s, m, u, rate)] = max(recall_at(S["s"], S["lab"][u], targets, cap) for S in ss)
    for s in SUBS:
        for u in TOL:
            for rate in RATES:
                v = [out[k] for k in [(s, a, u, rate) for a in ("GaMMA", "PyOcto", "GaMMA+SC", "PyOcto+SC")]
                     if k in out]
                if v:
                    out[(s, "BESTPICK", u, rate)] = max(v)
    return out


def boot(sets, cat, targets, B, seed=7):
    """own block bootstrap (60 s blocks over [900, 3560), 45 blocks, last of 20 s), primary only:
    pooled M2 - BESTPICK, LOC+T, 30/h; cap from the resampled duration (sum of block lengths)."""
    nb = int(math.ceil((HI - LO) / 60.0))
    blen = np.array([min(60.0, HI - LO - 60.0 * b) for b in range(nb)])
    tblk = ((cat[targets, 2] - LO) // 60).astype(int)
    rng = np.random.default_rng(seed)
    res = []
    for _ in range(B):
        wb = np.bincount(rng.integers(0, nb, nb), minlength=nb).astype(float)
        cap = int(math.floor(30 * (wb * blen).sum() / 3600.0 + 1e-9))
        diff = []
        for s in SUBS:
            def rec(m):
                best = 0.0
                for S in sets[(s, m)]:
                    if S["n"] == 0:
                        continue
                    wd = wb[((S["t"] - LO) // 60).astype(int)]
                    best = max(best, recall_at(S["s"], S["lab"]["LOC+T"], targets, cap, wd, wb[tblk]))
                return best
            bp = max(rec(a) for a in ("GaMMA", "PyOcto", "GaMMA+SC", "PyOcto+SC"))
            diff.append(rec("M2") - bp)
        res.append(np.mean(diff))
    return np.percentile(res, [2.5, 97.5]).tolist()


def main(m2dir, closedir, pdir, outp, B=0):
    cat, centre = catalogue()
    targets = np.where((cat[:, 2] >= LO) & (cat[:, 2] < HI) & (cat[:, 3] >= 2.0))[0]
    sets = load(m2dir, closedir, pdir, cat)
    lines = [f"centre {centre}, targets {len(targets)}"]
    res = {}
    for hname, hours in (("2660s", (HI - LO) / 3600.0), ("45blk", 0.75)):
        pt = point(sets, targets, hours)
        res[hname] = {"|".join(map(str, k)): v for k, v in pt.items()}
        lines.append(f"\n== hours convention {hname} (caps " + ", ".join(
            f"{r}/h->{int(math.floor(r * hours + 1e-9))}" for r in RATES) + ")")
        meths = [m for m in ("M1", "M2", "M3", "GaMMA", "PyOcto", "GaMMA+SC", "PyOcto+SC", "BESTPICK")
                 if any(k[1] == m for k in pt)]
        for u in TOL:
            for rate in RATES:
                lines.append(f"  {u:5s} {rate:2d}/h " + "  ".join(
                    f"{m}:" + ",".join(f"{pt.get((s, m, u, rate), float('nan')):.3f}" for s in SUBS) for m in meths))
        for a in ("M2", "M1", "M3"):
            for u in TOL:
                for rate in RATES:
                    per = [pt.get((s, a, u, rate), np.nan) - pt.get((s, "BESTPICK", u, rate), np.nan) for s in SUBS]
                    lines.append(f"  {a}-BESTPICK {u:5s} {rate:2d}/h pooled {np.mean(per):+.3f}  per "
                                 + " ".join(f"{s}:{x:+.3f}" for s, x in zip(SUBS, per)))
        for u in TOL:
            for rate in RATES:
                per = [pt.get((s, "M1", u, rate), np.nan) - pt.get((s, "M3", u, rate), np.nan) for s in SUBS]
                lines.append(f"  M1-M3 {u:5s} {rate:2d}/h pooled {np.mean(per):+.3f}")
                per = [pt.get((s, "M2", u, rate), np.nan) - pt.get((s, "M1", u, rate), np.nan) for s in SUBS]
                lines.append(f"  M2-M1 {u:5s} {rate:2d}/h pooled {np.mean(per):+.3f}")
    lines.append("\nDescriptive (after de-dup, in segment): detections, any-score LOC+T and LOC hits (distinct "
                 "M>=2 targets, max over configs)")
    for (s, m), ss in sorted(sets.items()):
        n = max(S["n"] for S in ss)
        h1 = max(len(set(S["lab"]["LOC+T"][S["lab"]["LOC+T"] >= 0].tolist())) for S in ss)
        h1d = max(int((S["lab"]["LOC+T"] >= 0).sum()) for S in ss)
        h2 = max(len(set(S["lab"]["LOC"][S["lab"]["LOC"] >= 0].tolist())) for S in ss)
        lines.append(f"  {s} {m:9s} n {n:5d}  LOC+T hit-detections {h1d:3d} distinct targets {h1:3d}  "
                     f"LOC distinct {h2:3d}")
    if B:
        ci = boot(sets, cat, targets, B)
        res["boot_primary_ci"] = ci
        lines.append(f"\nOwn block bootstrap (B={B}, seed 7) pooled M2-BESTPICK LOC+T 30/h: [{ci[0]:+.3f}, {ci[1]:+.3f}]")
    txt = "\n".join(lines)
    print(txt)
    open(outp.replace(".json", ".txt"), "w").write(txt + "\n")
    json.dump(res, open(outp, "w"))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]) if len(sys.argv) > 5 else 0)
