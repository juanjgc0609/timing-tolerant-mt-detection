"""Task R0 analysis from stored detections (port of code/r0/r0_analyze.py, main).

    python -m ttmt.analysis.r0 SEGMENT DIR [B] [--out OUT.txt] [--match legacy|one_to_one]

Reads DIR/det_, detM3_ and picks_<segment>_<subnet>.json; writes the text and the json next to OUT (default
DIR/analysis_<segment>.txt), as the legacy script. Needs the catalogue (raw data folder ridgecrest_2019).
"""
import json
import os
import sys

import numpy as np

from .. import scoring as S
from ..datasets import RIDGECREST as R

RULES = {"LOC+T": 0.5, "LOC": 3.0}
RATES = (10, 30, 60)
SUBS = ("N8", "N6", "N4", "U8")
DET = ("M1", "M2", "M3", "M5", "M6")
ASSOC = ("GaMMA", "PyOcto", "GaMMA+SC", "PyOcto+SC")
BLOCK = 60.0


def load(segment, d, match="legacy"):
    lo, hi = R.segment(segment)
    cat = R.catalogue_local(R.frame())
    M = np.array([e["M"] for e in cat])
    tg = np.array([i for i, e in enumerate(cat) if lo <= e["t"] < hi and e["M"] >= R.m_target])
    tblk = np.array([int((cat[i]["t"] - lo) // BLOCK) for i in tg])
    sets = {}
    for sub in SUBS:
        fp = os.path.join(d, f"picks_{segment}_{sub}.json")
        for pre in ("det", "detM3"):
            fd = os.path.join(d, f"{pre}_{segment}_{sub}.json")
            if os.path.exists(fd):
                r = json.load(open(fd))
                for m in DET:
                    if m in r["atoms"]:
                        sets[(sub, m)] = [S.DetSet(S.dedupe(r["atoms"][m]), cat, lo, hi, R.m_target, RULES, BLOCK, match)]
        if os.path.exists(fp):
            r = json.load(open(fp))
            for nm, key in (("GaMMA", "gamma"), ("PyOcto", "pyocto"), ("GaMMA+SC", "gamma_sc"), ("PyOcto+SC", "pyocto_sc")):
                if key in r:
                    sets[(sub, nm)] = [S.DetSet(v, cat, lo, hi, R.m_target, RULES, BLOCK, match) for v in r[key].values()]
    nb = int(np.ceil((hi - lo) / BLOCK))
    return cat, M, (tg, tblk), sets, nb, lo, hi


def endpoints(sets, tmask, M, wb):
    out = {}
    for (sub, m), ss in sets.items():
        for u in RULES:
            for r in RATES:
                out[(sub, m, u, r)] = S.froc_blocks(ss, tmask, M, wb, r, u, R.m_target, BLOCK)
    for sub in SUBS:
        for u in RULES:
            for r in RATES:
                v = [out[(sub, a, u, r)] for a in ASSOC if (sub, a, u, r) in out]
                if v:
                    out[(sub, "BESTPICK", u, r)] = max(v)
    return out


def analyze(segment, d, B=2000, match="legacy"):
    cat, M, tmask, sets, nb, lo, hi = load(segment, d, match)
    tg, _tblk = tmask
    lines = [f"R0 analysis, segment {segment} [{lo}, {hi}) s, {len(tg)} target events (M >= {R.m_target}), "
             f"{nb} blocks of {BLOCK:.0f} s, B = {B}"]
    pt = endpoints(sets, tmask, M, np.ones(nb))
    rng = np.random.default_rng(20261010)
    reps = []
    for _ in range(B):
        wb = np.bincount(rng.integers(0, nb, nb), minlength=nb).astype(float)
        reps.append(endpoints(sets, tmask, M, wb))
    keys = sorted(pt.keys(), key=str)

    def ci(vals):
        return np.percentile(vals, [2.5, 97.5])

    def pooled(rep, m, u, r):
        v = [rep[(s, m, u, r)] for s in SUBS if (s, m, u, r) in rep]
        return float(np.mean(v)) if len(v) == len(SUBS) else np.nan

    out = {}
    lines.append("\nRecall per subnet and method (point estimate), rule / rate per hour:")
    meths = [m for m in DET + ASSOC + ("BESTPICK",) if any(k[1] == m for k in keys)]
    for u in RULES:
        for r in RATES:
            lines.append(f"  {u:5s} {r:3d}/h  " + "  ".join(
                f"{m}:" + ",".join(f"{pt.get((s, m, u, r), np.nan):.3f}" for s in SUBS) for m in meths))
    contrasts = [("M2", "BESTPICK"), ("M2", "M1"), ("M1", "M3"), ("M2", "M5"), ("M1", "BESTPICK"), ("M3", "BESTPICK"),
                 ("M6", "BESTPICK"), ("M6", "M2")]
    lines.append("\nContrasts (pooled = mean of the 4 subnets), point [95% block-bootstrap interval]; per subnet points:")
    for a, b in contrasts:
        for u in RULES:
            for r in RATES:
                p = pooled(pt, a, u, r) - pooled(pt, b, u, r)
                v = np.array([pooled(rp, a, u, r) - pooled(rp, b, u, r) for rp in reps])
                if not np.isfinite(p):
                    continue
                c = ci(v)
                per = [pt.get((s, a, u, r), np.nan) - pt.get((s, b, u, r), np.nan) for s in SUBS]
                lines.append(f"  {a}-{b:8s} {u:5s} {r:3d}/h  {p:+.3f} [{c[0]:+.3f}, {c[1]:+.3f}]   per subnet "
                             + " ".join(f"{s}:{x:+.3f}" for s, x in zip(SUBS, per)))
                out[f"{a}-{b}|{u}|{r}"] = {"point": p, "ci": c.tolist(), "per": per}
    k = "M2-BESTPICK|LOC+T|30"
    if k in out:
        ra = out[k]["ci"][0] > -0.02
        rb = sum(x >= 0 for x in out[k]["per"]) >= 2
        sup = out[k]["ci"][0] > 0
        lines.append(f"\nREGISTERED DECISION: R0-a (lower bound > -0.02): {ra}; R0-b (>= 2 of 4 subnets): {rb}; "
                     f"superiority (lower bound > 0): {sup}")
    lines.append("\nDescriptive: detections in segment per subnet and method (after de-duplication / per config max):")
    for (sub, m), ss in sorted(sets.items()):
        n = max(len(S_.s) for S_ in ss)
        nh = max(int(((S_.lab['LOC+T'] >= 0) & (M[np.maximum(S_.lab['LOC+T'], 0)] >= R.m_target)).sum()) for S_ in ss)
        lines.append(f"  {sub} {m:7s} detections {n:5d}   LOC+T hits on M>=2 (any score) {nh}")
    return "\n".join(lines), {"point": {"|".join(map(str, k)): v for k, v in pt.items()}, "contrasts": out}


def main(argv=None):
    from ..runtime import pin_blas
    pin_blas(strict=False)
    a = list(sys.argv[1:] if argv is None else argv)
    opts = {"--out": None, "--match": "legacy"}
    for k in list(opts):
        if k in a:
            i = a.index(k)
            opts[k] = a[i + 1]
            del a[i:i + 2]
    segment, d = a[0], a[1]
    B = int(a[2]) if len(a) > 2 else 2000
    txt, js = analyze(segment, d, B, opts["--match"])
    outp = opts["--out"] or os.path.join(d, f"analysis_{segment}.txt")
    with open(outp, "w", newline="\n") as fh:
        fh.write(txt + "\n")
    json.dump(js, open(os.path.splitext(outp)[0] + ".json", "w"))
    print(txt)


if __name__ == "__main__":
    main()
