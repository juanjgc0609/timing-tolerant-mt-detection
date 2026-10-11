"""attrib_analyze_verif.py - analysis of leak_attrib_verif.py parts: leakage flags and the terms-vs-T3-level attribution."""
import json, glob, os, sys
import numpy as np
from load_verif import *
W, N = load_all()
P = {}
for f in glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), "attrib_parts_verif", "*.json")):
    e = json.load(open(f)); P[(e["seed"], e["regime"], e["cell"], e["net"])] = e
print(f"# {len(P)} networks x regimes re-generated (fresh processes); nets {min(k[3] for k in P)}..{max(k[3] for k in P)}")
for fl in ("s0_reproduced", "s2_reproduced", "s2_blind_identical", "same_truth", "mt_s2t3_reproduced"):
    print(f"  {fl:20s}: {sum(e[fl] for e in P.values())}/{len(P)}")
dx = np.concatenate([np.array(e["cat_dx"]).ravel() for e in P.values()]); dt = np.concatenate([np.array(e["cat_dt"]) for e in P.values()])
print(f"  catalogue hypocentres seen by S0/S2 minus truth: position sd per coordinate {dx.std():.3f} km (mean {dx.mean():+.3f}; depth clipped), origin sd {dt.std():.3f} s (mean {dt.mean():+.4f}); "
      f"exactly equal to truth in {np.mean(np.abs(dx) < 1e-12):.4f} of coordinates")
ARMS = ("MT@S0+T3", "MT@S0+T3|sigS2", "MT@S2+T3|sigS0", "MT@S2+T3")
ws = []
for w in W:
    k = (w["seed"], w["reg"], w["cell"], w["net"])
    if w["kind"] == "std" and k in P:
        w2 = dict(w); w2["atoms"] = dict(w["atoms"])
        for a in ("MT@S0+T3|sigS2", "MT@S2+T3|sigS0"):
            w2["atoms"][a] = np.array(P[k]["rows"][w["w"]]["atoms"][a], float).reshape(-1, 7)
        ws.append(w2)
def recs(wl, m, tol):
    return [window_rec(w["atoms"][m], label(w["atoms"][m], w["truth"], w["taus"], tol)) for w in wl]
print("\n# Attribution on STANDARD windows of the re-generated subset: recall at cap 0.25 (0.1), LOC+T and LOC")
for seeds in ((7171,), (9393,), (7171, 9393)):
    for cells in (CELLS, ("unilateral_c40",)):
        wl = [w for w in ws if w["seed"] in seeds and w["cell"] in cells]
        ks = np.array([w["k"] for w in wl], float)
        for tol, rule in ((0.5, "LOC+T"), (None, "LOC")):
            r = {m: froc(recs(wl, m, tol), np.ones(len(wl)), ks.sum())[0] for m in ARMS}
            print(f"seeds {seeds} cells {'all' if len(cells) == 3 else cells[0]} {rule:5s} ({len(wl)} windows, {int(ks.sum())} events): " + "  ".join(f"{m} {r[m][1]:.3f} ({r[m][0]:.3f})" for m in ARMS))
wl = ws; ks = np.array([w["k"] for w in wl], float)
R_ = {m: recs(wl, m, 0.5) for m in ARMS}
strata = {}
for i, w in enumerate(wl):
    strata.setdefault((w["seed"], w["reg"], w["cell"]), {}).setdefault(w["net"], []).append(i)
rng = np.random.default_rng(777); B = 500; d = {m: [] for m in ARMS}
for _ in range(B):
    mult = np.zeros(len(wl))
    for nets in strata.values():
        keys = list(nets)
        for g in rng.integers(0, len(keys), len(keys)):
            wi = np.asarray(nets[keys[g]]); np.add.at(mult, wi[rng.integers(0, len(wi), len(wi))], 1.0)
    nte = float((mult * ks).sum())
    for m in ARMS:
        d[m].append(froc(R_[m], mult, nte)[0][1])
d = {m: np.array(v) for m, v in d.items()}
pt = {m: froc(R_[m], np.ones(len(wl)), ks.sum())[0][1] for m in ARMS}
print(f"\n# Both seeds pooled, LOC+T, cap 0.25, bootstrap {B} (strata seed x regime x cell, networks, windows)")
for a, b, nm in (("MT@S2+T3", "MT@S0+T3", "registered contrast (terms + level)"), ("MT@S2+T3|sigS0", "MT@S0+T3", "terms only (S0 level in both arms)"),
                 ("MT@S2+T3", "MT@S0+T3|sigS2", "terms only (S2 level in both arms)"), ("MT@S0+T3|sigS2", "MT@S0+T3", "level only (S0 terms)"),
                 ("MT@S2+T3", "MT@S2+T3|sigS0", "level only (S2 terms)")):
    x = d[a] - d[b]
    print(f"  {nm:38s} {a} - {b}: {pt[a]-pt[b]:+.4f} [{np.percentile(x,2.5):+.4f}, {np.percentile(x,97.5):+.4f}]")
