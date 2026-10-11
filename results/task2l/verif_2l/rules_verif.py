"""rules_verif.py - own re-derivation of 2L-a, 2L-b, 2L-d (points, both seeds) and own hierarchical bootstrap of 2L-a/b."""
import sys
import numpy as np
from load_verif import *
W, N = load_all()
A, B_ = "MT@S2+T3", "MT@S0+T3"
BREP = int(sys.argv[1]) if len(sys.argv) > 1 else 600

def subset(seed, kind, reg=None, cell=None, fam_ok=None):
    return [w for w in W if w["seed"] == seed and w["kind"] == kind and (reg is None or w["reg"] == reg)
            and (cell is None or w["cell"] == cell) and (fam_ok is None or fam_ok in w["atoms"])]

def recs_of(ws, m, tol_t, running=True, key=4):
    return [window_rec(w["atoms"][m], label(w["atoms"][m], w["truth"], w["taus"], tol_t), running, key) for w in ws]

def point(ws, m, tol_t, **kw):
    ks = np.array([w["k"] for w in ws], float)
    return froc(recs_of(ws, m, tol_t, **kw), np.ones(len(ws)), ks.sum())

print("# Own re-derivation of the rule contrasts (points), recall at cap 0.1 / 0.25")
for seed in (7171, 9393):
    for kind in ("std", "dense"):
        ws = subset(seed, kind)
        for rule, tol in (("LOC+T", 0.5), ("LOC", None)):
            ra, rb, r0 = point(ws, A, tol)[0], point(ws, B_, tol)[0], point(ws, "MT@S0", tol)[0]
            print(f"{seed} {kind:5s} {rule:5s} nwin {len(ws)} nev {int(sum(w['k'] for w in ws))}: S0 {r0[0]:.3f}/{r0[1]:.3f}  S0+T3 {rb[0]:.3f}/{rb[1]:.3f}  "
                  f"S2+T3 {ra[0]:.3f}/{ra[1]:.3f}  diff(S2+T3 - S0+T3) {ra[0]-rb[0]:+.4f}/{ra[1]-rb[1]:+.4f}")
        for reg in REGIMES:
            ws2 = subset(seed, kind, reg=reg)
            ra, rb = point(ws2, A, 0.5)[0], point(ws2, B_, 0.5)[0]
            print(f"   2L-c {seed} {kind} {reg}: LOC+T diff at 0.25 {ra[1]-rb[1]:+.4f}")
    ws = subset(seed, "std", fam_ok="BOMP@S2+T3")
    ra, rb = point(ws, "BOMP@S2+T3", 0.5)[0], point(ws, "BOMP@S0+T3", 0.5)[0]
    print(f"{seed} BOMP std LOC+T nwin {len(ws)}: S0+T3 {rb[1]:.3f} S2+T3 {ra[1]:.3f} diff {ra[1]-rb[1]:+.4f}")
    # sensitivity: raw scores instead of running minimum; sref (refined objective) is NOT the record score
    ws = subset(seed, "std")
    ra, rb = point(ws, A, 0.5, running=False)[0], point(ws, B_, 0.5, running=False)[0]
    print(f"{seed} std LOC+T raw (no running-min) score: diff at 0.25 {ra[1]-rb[1]:+.4f}")

print(f"\n# Own hierarchical bootstrap ({BREP} replicates; strata regime x cell, networks then windows; operating point re-chosen)")
for seed in (7171, 9393):
    for kind, rule, tol in (("std", "LOC+T", 0.5), ("dense", "LOC+T", 0.5), ("std", "LOC", None), ("dense", "LOC", None)):
        ws = subset(seed, kind)
        ks = np.array([w["k"] for w in ws], float)
        ra_, rb_ = recs_of(ws, A, tol), recs_of(ws, B_, tol)
        strata = {}
        for i, w in enumerate(ws):
            strata.setdefault((w["reg"], w["cell"]), {}).setdefault(w["net"], []).append(i)
        rng = np.random.default_rng([20261005, seed, kind == "dense", tol is None])
        d = []
        for b in range(BREP):
            mult = np.zeros(len(ws))
            for st, nets in strata.items():
                keys = list(nets)
                for g in rng.integers(0, len(keys), len(keys)):
                    wl = np.asarray(nets[keys[g]])
                    np.add.at(mult, wl[rng.integers(0, len(wl), len(wl))], 1.0)
            nte = float((mult * ks).sum())
            d.append(froc(ra_, mult, nte)[0][1] - froc(rb_, mult, nte)[0][1])
        d = np.array(d)
        pt = froc(ra_, np.ones(len(ws)), ks.sum())[0][1] - froc(rb_, np.ones(len(ws)), ks.sum())[0][1]
        print(f"{seed} {kind:5s} {rule:5s}: point {pt:+.4f}  95% CI [{np.percentile(d,2.5):+.4f}, {np.percentile(d,97.5):+.4f}]  5th pct {np.percentile(d,5):+.4f}  "
              f"sd {d.std():.4f}  frac(d<=0) {np.mean(d<=0):.4f}", flush=True)
