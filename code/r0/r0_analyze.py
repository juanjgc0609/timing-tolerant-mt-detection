"""
r0_analyze.py - Task R0 analysis (registered endpoint and contrasts).
    python r0_analyze.py SEGMENT DIR [B]      reads DIR/det_, detM3_ and picks_<segment>_<subnet>.json
                                              writes DIR/analysis_<segment>.txt (and .json)
Rules: LOC+T = epicentre <= 5 km and |origin error| <= 0.5 s of a catalogued event (any magnitude); LOC = <= 5 km and
<= 3.0 s.  A detection hitting an event of M >= 2.0 is a hit; one hitting only events of M < 2.0 is neutral; one hitting
nothing is false (conservative count).  Pick-free atoms of all windows are pooled and an atom within 5 km and 1.0 s of a
higher-scoring atom is removed.  Exact FROC at 10, 30 (primary) and 60 false detections per hour; for GaMMA and PyOcto
the envelope over their configuration grid (score = picks per event), with raw and +SC-corrected picks; BESTPICK = the
best of the four comparator sets in each subnet and replicate.  Block bootstrap of 60 s blocks of the segment
(same blocks for all subnets and methods), 2000 replicates, operating point re-chosen in each replicate.
"""
import json
import os
import sys

import numpy as np

import r0_lib as R

RULES = {"LOC+T": 0.5, "LOC": 3.0}
RATES = (10, 30, 60)
SUBS = ("N8", "N6", "N4", "U8")
DET = ("M1", "M2", "M3", "M5", "M6")
ASSOC = ("GaMMA", "PyOcto", "GaMMA+SC", "PyOcto+SC")
BLOCK = 60.0


def seg_range(segment):
    return R.T_CAL if segment == "cal" else R.T_EVAL


def dedupe(atoms, dkm=5.0, dts=1.0):
    a = sorted(atoms, key=lambda r: -r[4])
    keep = []
    for r in a:
        if all(np.hypot(r[0] - k[0], r[1] - k[1]) > dkm or abs(r[3] - k[3]) > dts for k in keep):
            keep.append(r)
    return keep


def label(det, cat, tol):
    """per detection: index of the hit event (nearest epicentre among those within 5 km and tol s, events of M >= 2.0
    first) or -1."""
    ex = np.array([[e["x"], e["y"], e["t"]] for e in cat])
    big = np.array([e["M"] >= R.M_TARGET for e in cat])
    out = []
    for d in det:
        dd = np.hypot(ex[:, 0] - d[0], ex[:, 1] - d[1])
        ok = (dd <= 5.0) & (np.abs(ex[:, 2] - d[3]) <= tol + 1e-9)
        if (ok & big).any():                                    # an event of M >= 2.0 is preferred
            ok = ok & big
        out.append(int(np.where(ok, dd, np.inf).argmin()) if ok.any() else -1)
    return np.array(out, int)


class Set:
    """one detection list of one method on one subnet: scores, labels per rule, block of each detection."""

    def __init__(self, det, cat, lo, hi):
        det = [d for d in det if lo <= d[3] < hi]
        self.s = np.array([d[4] for d in det], float)
        self.blk = np.array([int((d[3] - lo) // BLOCK) for d in det], int)
        self.lab = {u: label(det, cat, tol) if det else np.zeros(0, int) for u, tol in RULES.items()}
        self.det = det


def froc(sets, targets_mask, M, n_blocks, wb, rate, rule):
    """recall at <= rate false per hour for a method given as a list of Sets (envelope over configurations).
    wb: weight (multiplicity) of every block in this replicate.  targets: catalogue indices with M >= 2.0 in segment."""
    hours = wb.sum() * BLOCK / 3600.0
    cap = int(np.floor(rate * hours + 1e-9))
    tw = wb[targets_mask[1]]                                    # weight of every target (by its block)
    ntar = tw.sum()
    best = 0.0
    for S in sets:
        if not len(S.s):
            continue
        w = wb[S.blk]
        lab = S.lab[rule]
        neutral = (lab >= 0) & (M[np.maximum(lab, 0)] < R.M_TARGET)
        fal = lab < 0
        # false scores with multiplicity
        fs, fw = S.s[fal], w[fal]
        o = np.argsort(-fs)
        cum = np.cumsum(fw[o])
        k = np.searchsorted(cum, cap, side="right")             # number of false allowed = those with cum <= cap
        thr = fs[o][k] if k < len(fs) else -np.inf
        hit = (lab >= 0) & ~neutral & (S.s > thr)
        got = np.zeros(len(M), bool)
        got[lab[hit]] = True
        rec = float((got[targets_mask[0]] * tw).sum() / max(ntar, 1e-12))
        best = max(best, rec)
    return best


def load(segment, d):
    lo, hi = seg_range(segment)
    cat = R.catalogue_local(R.frame())
    M = np.array([e["M"] for e in cat])
    tg = np.array([i for i, e in enumerate(cat) if lo <= e["t"] < hi and e["M"] >= R.M_TARGET])
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
                        sets[(sub, m)] = [Set(dedupe(r["atoms"][m]), cat, lo, hi)]
        if os.path.exists(fp):
            r = json.load(open(fp))
            for nm, key in (("GaMMA", "gamma"), ("PyOcto", "pyocto"), ("GaMMA+SC", "gamma_sc"), ("PyOcto+SC", "pyocto_sc")):
                if key in r:
                    sets[(sub, nm)] = [Set(v, cat, lo, hi) for v in r[key].values()]
    nb = int(np.ceil((hi - lo) / BLOCK))
    return cat, M, (tg, tblk), sets, nb, lo, hi


def endpoints(sets, tmask, M, nb, wb):
    out = {}
    for (sub, m), ss in sets.items():
        for u in RULES:
            for r in RATES:
                out[(sub, m, u, r)] = froc(ss, tmask, M, nb, wb, r, u)
    for sub in SUBS:
        for u in RULES:
            for r in RATES:
                v = [out[(sub, a, u, r)] for a in ASSOC if (sub, a, u, r) in out]
                if v:
                    out[(sub, "BESTPICK", u, r)] = max(v)
    return out


def main(segment, d, B=2000):
    cat, M, tmask, sets, nb, lo, hi = load(segment, d)
    tg, tblk = tmask
    lines = [f"R0 analysis, segment {segment} [{lo}, {hi}) s, {len(tg)} target events (M >= {R.M_TARGET}), "
             f"{nb} blocks of {BLOCK:.0f} s, B = {B}"]
    pt = endpoints(sets, tmask, M, nb, np.ones(nb))
    rng = np.random.default_rng(20261010)
    reps = []
    for _ in range(B):
        wb = np.bincount(rng.integers(0, nb, nb), minlength=nb).astype(float)
        reps.append(endpoints(sets, tmask, M, nb, wb))
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
            row = f"  {u:5s} {r:3d}/h  " + "  ".join(
                f"{m}:" + ",".join(f"{pt.get((s, m, u, r), np.nan):.3f}" for s in SUBS) for m in meths)
            lines.append(row)
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
    # registered decision
    k = "M2-BESTPICK|LOC+T|30"
    if k in out:
        ra = out[k]["ci"][0] > -0.02
        rb = sum(x >= 0 for x in out[k]["per"]) >= 2
        sup = out[k]["ci"][0] > 0
        lines.append(f"\nREGISTERED DECISION: R0-a (lower bound > -0.02): {ra}; R0-b (>= 2 of 4 subnets): {rb}; "
                     f"superiority (lower bound > 0): {sup}")
    # descriptive: false per hour at the 30/h operating point is by construction <= 30; detections and hits counts
    lines.append("\nDescriptive: detections in segment per subnet and method (after de-duplication / per config max):")
    for (sub, m), ss in sorted(sets.items()):
        n = max(len(S.s) for S in ss)
        nh = max(int(((S.lab['LOC+T'] >= 0) & (M[np.maximum(S.lab['LOC+T'], 0)] >= R.M_TARGET)).sum()) for S in ss)
        lines.append(f"  {sub} {m:7s} detections {n:5d}   LOC+T hits on M>=2 (any score) {nh}")
    txt = "\n".join(lines)
    open(os.path.join(d, f"analysis_{segment}.txt"), "w").write(txt + "\n")
    json.dump({"point": {"|".join(map(str, k)): v for k, v in pt.items()}, "contrasts": out},
              open(os.path.join(d, f"analysis_{segment}.json"), "w"))
    print(txt)


def deploy(d, rate=30, rule="LOC+T"):
    """secondary: operating point fixed on CAL (threshold, and configuration for the associators, maximising CAL recall
    at <= rate false per hour), applied to EVAL; realised false per hour and recall on EVAL."""
    c_cat, c_M, c_t, c_sets, c_nb, c_lo, c_hi = load("cal", d)
    e_cat, e_M, e_t, e_sets, e_nb, e_lo, e_hi = load("eval", d)
    lines = [f"\nDeployable operating point (fixed on CAL at <= {rate}/h, {rule}), realised on EVAL:"]
    for key in sorted(c_sets):
        if key not in e_sets:
            continue
        best = None
        for ci_, S in enumerate(c_sets[key]):
            if not len(S.s):
                continue
            lab = S.lab[rule]
            fal = lab < 0
            cap = int(np.floor(rate * (c_hi - c_lo) / 3600.0 + 1e-9))
            fs = np.sort(S.s[fal])[::-1]
            thr = fs[cap] if len(fs) > cap else -np.inf
            hit = (lab >= 0) & (S.s > thr)
            got = set(lab[hit]) & set(c_t[0].tolist())
            rec = len(got) / max(1, len(c_t[0]))
            if best is None or rec > best[0]:
                best = (rec, ci_, thr)
        if best is None:
            continue
        S = e_sets[key][best[1]]
        lab = S.lab[rule]
        sel = S.s > best[2]
        fal = int(((lab < 0) & sel).sum())
        got = set(lab[(lab >= 0) & sel]) & set(e_t[0].tolist())
        hrs = (e_hi - e_lo) / 3600.0
        lines.append(f"  {key[0]} {key[1]:9s} CAL recall {best[0]:.3f}  EVAL recall {len(got) / max(1, len(e_t[0])):.3f}"
                     f"  EVAL false/h {fal / hrs:6.1f}")
    txt = "\n".join(lines)
    open(os.path.join(d, "analysis_deploy.txt"), "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    if sys.argv[1] == "deploy":
        deploy(sys.argv[2])
        sys.exit(0)
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 2000)
