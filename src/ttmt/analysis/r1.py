"""Task R1 analysis from the stored detection lists (port of code/r1/r1_analyze.py).

    python -m ttmt.analysis.r1 SEGMENT HYB_SEGMENT.json [HYB_CAL.json] OUT.txt [--match legacy|one_to_one]

Targets: truth (Match&Locate) events of the segment with magnitude >= 0.5; exact FROC at 10, 30 (primary) and 60 false
detections per hour; R1-a, R1-b (exact one-sided McNemar), R1-c, the 60 s block bootstrap and the deployable
threshold. Needs the catalogue files (raw data folder hinet_2012).
"""
import json
import sys
from math import comb

import numpy as np

from .. import scoring as S
from ..datasets import HINET as R

RATES = (10, 30, 60)
RULES = ("LOC+T", "LOC")
RULE_TOL = {"LOC+T": 0.5, "LOC": 3.0}
BLOCK = 60.0
B = 2000
SEED = 20261017


def _set(det, cat, lo, hi, match):
    return S.DetSet(det, cat, lo, hi, R.m_target, RULE_TOL, BLOCK, match)


def operating(St, tmask, M, rate, rule, hours):
    """recall, set of target indices hit, score threshold (hits need s > thr) of one DetSet at full weights."""
    cap = int(np.floor(rate * hours + 1e-9))
    if not len(St.s):
        return 0.0, set(), np.inf
    lab = St.lab[rule]
    neutral = (lab >= 0) & (M[np.maximum(lab, 0)] < R.m_target)
    fs = np.sort(St.s[lab < 0])[::-1]
    thr = fs[cap] if cap < len(fs) else -np.inf
    hit = (lab >= 0) & ~neutral & (St.s > thr)
    got = set(int(i) for i in lab[hit]) & set(int(i) for i in tmask[0])
    return len(got) / max(len(tmask[0]), 1), got, thr


def best_of(sets, tmask, M, rate, rule, hours):
    best = None
    for k, St in sets:
        r = operating(St, tmask, M, rate, rule, hours)
        if best is None or r[0] > best[1][0]:
            best = (k, r)
    return best


def mcnemar(b, c):
    n = b + c
    return 1.0 if n == 0 else sum(comb(n, k) for k in range(b, n + 1)) / 2.0 ** n


def load(d, cat, lo, hi, match):
    out = {m: [(m, _set(d[m], cat, lo, hi, match))] for m in ("BASE", "HYB", "P1")}
    out["BESTPICK"] = [(k, _set(v, cat, lo, hi, match)) for k, v in d["BESTPICK"].items()]
    for arm in ("gamma", "pyocto", "gamma_sc", "pyocto_sc"):
        out["R0-" + arm] = [(k, _set(v, cat, lo, hi, match)) for k, v in d["R0STYLE"].items() if k.split("|")[0] == arm]
    out["R0-BESTPICK"] = [x for arm in ("gamma", "pyocto", "gamma_sc", "pyocto_sc") for x in out["R0-" + arm]]
    return out


def analyze(segment, hyb, hyb_cal=None, match="legacy"):
    """hyb, hyb_cal: the loaded detection-list dicts (r1_hybrid output). Returns the report text."""
    lo, hi = R.segment(segment)
    hours = (hi - lo) / 3600.0
    cat = R.catalogue_local(R.frame())
    M = np.array([e["M"] for e in cat])
    tg = np.array([i for i, e in enumerate(cat) if lo <= e["t"] < hi and e["M"] >= R.m_target])
    tj = np.array([i for i, e in enumerate(cat) if lo <= e["t"] < hi and e["jma"]])
    nb = int(np.ceil((hi - lo) / BLOCK))
    tmask = (tg, np.array([int((cat[i]["t"] - lo) // BLOCK) for i in tg]))
    jmask = (tj, np.array([int((cat[i]["t"] - lo) // BLOCK) for i in tj]))
    d = hyb
    sets = load(d, cat, lo, hi, match)
    L = [f"R1 analysis, segment {segment} [{lo}, {hi}) s ({hours:.4f} h), {len(tg)} targets (truth M >= {R.m_target}), "
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
    rng = np.random.default_rng(SEED)
    con = (("HYB", "BASE"), ("HYB", "BESTPICK"), ("P1", "BASE"))
    vals = {c: [] for c in con}
    for _b in range(B):
        wb = np.bincount(rng.integers(0, nb, nb), minlength=nb).astype(float)
        f = {m: S.froc_blocks([St for _k, St in sets[m]], tmask, M, wb, 30, "LOC+T", R.m_target, BLOCK)
             for m in ("HYB", "BASE", "BESTPICK", "P1")}
        for a, b in con:
            vals[(a, b)].append(f[a] - f[b])
    for a, b in con:
        v = np.array(vals[(a, b)])
        L.append(f"block bootstrap (60 s blocks, {B} replicates) {a} - {b}, LOC+T 30/h: median {np.median(v):+.3f}, "
                 f"95% percentile interval [{np.percentile(v, 2.5):+.3f}, {np.percentile(v, 97.5):+.3f}]")
    if hyb_cal is not None:
        L.append("")
        L.append("deployable threshold (CAL operating point, LOC+T 30/h, applied to the segment)")
        clo, chi = R.cal
        ctg = np.array([i for i, e in enumerate(cat) if clo <= e["t"] < chi and e["M"] >= R.m_target])
        cmask = (ctg, np.zeros(len(ctg), int))
        csets = load(hyb_cal, cat, clo, chi, match)
        for m in ("HYB", "BASE", "P1"):
            _r, _g, thr = operating(csets[m][0][1], cmask, M, 30, "LOC+T", (chi - clo) / 3600.0)
            St = sets[m][0][1]
            lab = St.lab["LOC+T"]
            sel = St.s > thr
            got = set(int(i) for i in lab[sel & (lab >= 0)]) & set(int(i) for i in tg)
            nf = int((sel & (lab < 0)).sum())
            L.append(f"   {m:5s} threshold {thr:.3f}: recall {len(got) / len(tg):.3f} ({len(got)}), {nf} false "
                     f"({nf / hours:.1f} per hour)")
    return "\n".join(L) + "\n"


def main(argv=None):
    from ..runtime import pin_blas
    pin_blas(strict=False)
    a = list(sys.argv[1:] if argv is None else argv)
    match = "legacy"
    if "--match" in a:
        i = a.index("--match")
        match = a[i + 1]
        del a[i:i + 2]
    txt = analyze(a[0], json.load(open(a[1])), json.load(open(a[2])) if len(a) == 4 else None, match)
    with open(a[-1], "w", newline="\n") as fh:
        fh.write(txt)
    print(txt, end="")


if __name__ == "__main__":
    main()
