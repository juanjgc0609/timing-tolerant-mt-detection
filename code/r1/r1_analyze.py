"""
r1_analyze.py - Task R1 analysis (registered endpoint and rules).
    python r1_analyze.py SEGMENT HYB_SEGMENT.json [HYB_CAL.json] OUT.txt
Targets = truth (Match&Locate) events in the segment with magnitude >= 0.5; a detection hitting only smaller events is
neutral, one hitting nothing is false (r0_analyze.label and froc, LOC+T 5 km / 0.5 s, LOC 5 km / 3.0 s).  Recall at
<= 10, 30 (primary), 60 false per hour by exact FROC; operating point chosen per method in the segment.
R1-a: HYB - BASE, LOC+T, 30/h: exact one-sided McNemar on the discordant targets (hit by one method at its operating
point and not by the other), met when p <= 0.025 and the difference is > 0.  R1-b: HYB - BESTPICK (the configuration
of the 1-D associators with the highest recall; first in order on ties), same test.  R1-c (descriptive): P1 - BASE,
JMA-only targets, 60 s block-bootstrap intervals (2000 replicates, operating point re-chosen), 10 and 60/h, LOC, the
deployable threshold (score threshold of the CAL operating point at 30/h, LOC+T, applied to the segment).
"""
import json
import os
import sys
from math import comb

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r1_lib as R  # noqa: E402

R.install()
import r0_analyze as A  # noqa: E402

RATES = (10, 30, 60)
RULES = ("LOC+T", "LOC")
B = 2000
SEED = 20261017


def operating(S, tmask, M, rate, rule, hours):
    """recall, set of target indices hit, score threshold (hits need s > thr) of one Set at full weights."""
    cap = int(np.floor(rate * hours + 1e-9))
    if not len(S.s):
        return 0.0, set(), np.inf
    lab = S.lab[rule]
    neutral = (lab >= 0) & (M[np.maximum(lab, 0)] < R.M_TARGET)
    fs = np.sort(S.s[lab < 0])[::-1]
    thr = fs[cap] if cap < len(fs) else -np.inf
    hit = (lab >= 0) & ~neutral & (S.s > thr)
    got = set(int(i) for i in lab[hit]) & set(int(i) for i in tmask[0])
    return len(got) / max(len(tmask[0]), 1), got, thr


def best_of(sets, tmask, M, rate, rule, hours):
    best = None
    for k, S in sets:
        r = operating(S, tmask, M, rate, rule, hours)
        if best is None or r[0] > best[1][0]:
            best = (k, r)
    return best


def mcnemar(b, c):
    n = b + c
    return 1.0 if n == 0 else sum(comb(n, k) for k in range(b, n + 1)) / 2.0 ** n


def load(d, cat, lo, hi):
    out = {m: [(m, A.Set(d[m], cat, lo, hi))] for m in ("BASE", "HYB", "P1")}
    out["BESTPICK"] = [(k, A.Set(v, cat, lo, hi)) for k, v in d["BESTPICK"].items()]
    for arm in ("gamma", "pyocto", "gamma_sc", "pyocto_sc"):
        out["R0-" + arm] = [(k, A.Set(v, cat, lo, hi)) for k, v in d["R0STYLE"].items() if k.split("|")[0] == arm]
    out["R0-BESTPICK"] = [x for arm in ("gamma", "pyocto", "gamma_sc", "pyocto_sc") for x in out["R0-" + arm]]
    return out


def main(segment, hp, calp, outp):
    lo, hi = R.T_CAL if segment == "cal" else R.T_EVAL
    hours = (hi - lo) / 3600.0
    cat = R.catalogue_local(R.frame())
    M = np.array([e["M"] for e in cat])
    tg = np.array([i for i, e in enumerate(cat) if lo <= e["t"] < hi and e["M"] >= R.M_TARGET])
    tj = np.array([i for i, e in enumerate(cat) if lo <= e["t"] < hi and e["jma"]])
    nb = int(np.ceil((hi - lo) / A.BLOCK))
    tmask = (tg, np.array([int((cat[i]["t"] - lo) // A.BLOCK) for i in tg]))
    jmask = (tj, np.array([int((cat[i]["t"] - lo) // A.BLOCK) for i in tj]))
    d = json.load(open(hp))
    sets = load(d, cat, lo, hi)
    L = [f"R1 analysis, segment {segment} [{lo}, {hi}) s ({hours:.4f} h), {len(tg)} targets (truth M >= {R.M_TARGET}), "
         f"{len(tj)} JMA events; BASE = {d['base'][0]} {d['base'][1]}",
         "detections: " + ", ".join(f"{m} {len(sets[m][0][1].s)}" for m in ("BASE", "HYB", "P1")), ""]
    res = {}
    L.append("recall (events) at <= rate false per hour, operating point chosen in the segment")
    L.append("   method          " + "  ".join(f"{u:>5s} {r:>2d}/h" for u in RULES for r in RATES))
    order = ("HYB", "BASE", "BESTPICK", "P1", "R0-BESTPICK", "R0-gamma", "R0-pyocto", "R0-gamma_sc", "R0-pyocto_sc")
    for m in order:
        row = []
        for u in RULES:
            for r in RATES:
                k, (rec, got, thr) = best_of(sets[m], tmask, M, r, u, hours)
                res[(m, u, r)] = (rec, got, k, thr)
                row.append(f"{rec:.3f} ({len(got):2d})")
        L.append(f"   {m:14s}  " + "  ".join(row))
    L.append(f"   BESTPICK configuration at LOC+T 30/h: {res[('BESTPICK', 'LOC+T', 30)][2]}")
    L.append("")
    L.append("JMA-only targets (descriptive)")
    for m in ("HYB", "BASE", "BESTPICK", "P1"):
        L.append(f"   {m:10s} " + "  ".join(f"LOC+T {r}/h {best_of(sets[m], jmask, M, r, 'LOC+T', hours)[1][0]:.3f}" for r in RATES))
    L.append("")
    tim = {i: (cat[i]["t"], cat[i]["M"]) for i in tg}
    for tag, a, b in (("R1-a", "HYB", "BASE"), ("R1-b", "HYB", "BESTPICK"), ("R1-c", "P1", "BASE")):
        ga, gb = res[(a, "LOC+T", 30)][1], res[(b, "LOC+T", 30)][1]
        nb_, nc_ = len(ga - gb), len(gb - ga)
        p = mcnemar(nb_, nc_)
        diff = (len(ga) - len(gb)) / len(tg)
        met = (p <= 0.025 and diff > 0)
        L.append(f"{tag}: {a} - {b}, LOC+T 30/h: {len(ga)} vs {len(gb)} of {len(tg)}, difference {diff:+.3f}; discordant "
                 f"{nb_} only {a}, {nc_} only {b}; exact one-sided McNemar p = {p:.4f}"
                 + ("" if tag == "R1-c" else f" -> {'MET' if met else 'NOT MET'}"))
        L.append(f"      only {a}: " + ", ".join(f"t {tim[i][0]:.1f} M {tim[i][1]:.2f}" for i in sorted(ga - gb, key=lambda i: tim[i][0])))
        L.append(f"      only {b}: " + ", ".join(f"t {tim[i][0]:.1f} M {tim[i][1]:.2f}" for i in sorted(gb - ga, key=lambda i: tim[i][0])))
    L.append("")
    # block bootstrap of the contrasts (operating point re-chosen in each replicate)
    rng = np.random.default_rng(SEED)
    con = (("HYB", "BASE"), ("HYB", "BESTPICK"), ("P1", "BASE"))
    vals = {c: [] for c in con}
    for _b in range(B):
        wb = np.bincount(rng.integers(0, nb, nb), minlength=nb).astype(float)
        f = {m: A.froc([S for _k, S in sets[m]], tmask, M, nb, wb, 30, "LOC+T") for m in ("HYB", "BASE", "BESTPICK", "P1")}
        for a, b in con:
            vals[(a, b)].append(f[a] - f[b])
    for a, b in con:
        v = np.array(vals[(a, b)])
        L.append(f"block bootstrap (60 s blocks, {B} replicates) {a} - {b}, LOC+T 30/h: median {np.median(v):+.3f}, "
                 f"95% percentile interval [{np.percentile(v, 2.5):+.3f}, {np.percentile(v, 97.5):+.3f}]")
    if calp:
        L.append("")
        L.append("deployable threshold (CAL operating point, LOC+T 30/h, applied to the segment)")
        dc = json.load(open(calp))
        clo, chi = R.T_CAL
        ctg = np.array([i for i, e in enumerate(cat) if clo <= e["t"] < chi and e["M"] >= R.M_TARGET])
        cmask = (ctg, np.zeros(len(ctg), int))
        csets = load(dc, cat, clo, chi)
        for m in ("HYB", "BASE", "P1"):
            _r, _g, thr = operating(csets[m][0][1], cmask, M, 30, "LOC+T", (chi - clo) / 3600.0)
            S = sets[m][0][1]
            lab = S.lab["LOC+T"]
            sel = S.s > thr
            got = set(int(i) for i in lab[sel & (lab >= 0)]) & set(int(i) for i in tg)
            nf = int((sel & (lab < 0)).sum())
            L.append(f"   {m:5s} threshold {thr:.3f}: recall {len(got) / len(tg):.3f} ({len(got)}), {nf} false "
                     f"({nf / hours:.1f} per hour)")
    open(outp, "w").write("\n".join(L) + "\n")
    print("\n".join(L), flush=True)


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], a[1], a[2] if len(a) == 4 else None, a[-1])
