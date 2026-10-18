"""
c19_analyze.py - analysis of exploration C19 (c19_run.py parts): recall at 0.1 and 0.25 false events per window for
each mechanism family and method, the change relative to double couples (same windows except the mechanism), the
method contrasts of 2I within each family and their change relative to DC, and the subspace fractions.  Hierarchical
bootstrap (strata regime x cell, networks then windows), the same resamples for every family and method, operating
point re-chosen per family, method and replicate (as in 2I).

    python c19_analyze.py [--B 1000] > c19_out.txt
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
for p in (HERE, WORK):
    if p not in sys.path:
        sys.path.insert(0, p)

from xc import recall_caps  # noqa: E402

CAPS = (0.1, 0.25)
SHORT = {"BOMP-PK+SC": "BOMP", "MT-FW+SC": "MTFW", "MT-FW-NOCAP+SC": "NOCAP", "HUBER-N+SC": "HUBER"}
CONTRASTS = (("MT-FW+SC", "BOMP-PK+SC"), ("HUBER-N+SC", "BOMP-PK+SC"), ("HUBER-N+SC", "MT-FW+SC"),
             ("MT-FW+SC", "MT-FW-NOCAP+SC"))


def load(parts):
    P = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(parts, "c19_*.json")))]
    if not P:
        raise SystemExit("no parts")
    return P


def build(P, kind):
    fams, methods = P[0]["families"], P[0]["methods"]
    groups, strata, ks = [], [], []
    recs = {(f, m): [] for f in fams for m in methods}
    frac = {f: [] for f in fams}
    sk = {}
    for o in P:
        g = []
        for wi in range(o["n_win"]):
            g.append(len(ks))
            ks.append(o["windows"][fams[0]][kind][wi]["k"])
            for f in fams:
                w = o["windows"][f][kind][wi]
                assert w["k"] == ks[-1] and w["truth"] == o["windows"][fams[0]][kind][wi]["truth"]
                for m in methods:
                    recs[(f, m)].append(w["records"][m])
                frac[f].extend(w["frac"])
        groups.append(g)
        strata.append(sk.setdefault((o["regime"], o["cell"]), len(sk)))
    return fams, methods, groups, strata, ks, recs, {f: np.array(v) for f, v in frac.items()}


def boot(fams, methods, groups, strata, ks, recs, B, seed=12345):
    rng = np.random.default_rng(seed)
    keys = [(f, m) for f in fams for m in methods]

    def rec_of(ws):
        n_ = sum(ks[w] for w in ws)
        return {k: np.array(recall_caps([recs[k][w] for w in ws], n_, CAPS)[0]) for k in keys}
    pt = rec_of([w for g in groups for w in g])
    sidx = {}
    for g, s_ in enumerate(strata):
        sidx.setdefault(s_, []).append(g)
    bs = []
    for _ in range(B):
        ws = []
        for gl in sidx.values():
            for g in np.asarray(gl)[rng.integers(0, len(gl), len(gl))]:
                gg = groups[g]
                ws.extend(np.asarray(gg)[rng.integers(0, len(gg), len(gg))].tolist())
        bs.append(rec_of(ws))
    return pt, bs


def ci(pt_val, arr):
    lo, hi = np.percentile(arr, [2.5, 97.5], axis=0)
    return pt_val, lo, hi


def fmt(v, lo, hi, j):
    return f"{v[j]:+.3f} [{lo[j]:+.3f}, {hi[j]:+.3f}]"


def report(P, kind, B):
    fams, methods, groups, strata, ks, recs, frac = build(P, kind)
    pt, bs = boot(fams, methods, groups, strata, ks, recs, B)
    nev = sum(ks)
    print(f"\n=========================== {kind} windows: {len(groups)} networks, {len(ks)} windows, {nev} events")
    print("\nSubspace fraction of the noise-free arrival energy at the true source (mean, 10th percentile):")
    print(f"   {'family':<9} {'clean (gradient medium)':>26} {'with jitter + site':>22}")
    for f in fams:
        a = frac[f]
        print(f"   {f:<9} {a[:, 0].mean():>14.3f} ({np.percentile(a[:, 0], 10):.3f}) "
              f"{a[:, 1].mean():>13.3f} ({np.percentile(a[:, 1], 10):.3f})")
    for j, cap in enumerate(CAPS):
        print(f"\nRecall at {cap} false events per window (95% CI):")
        print(f"   {'family':<9}" + "".join(f"{SHORT[m]:>24}" for m in methods))
        for f in fams:
            row = ""
            for m in methods:
                v, lo, hi = ci(pt[(f, m)], [b[(f, m)] for b in bs])
                row += f"   {v[j]:.3f} [{lo[j]:.3f}, {hi[j]:.3f}]"
            print(f"   {f:<9}{row}")
        print(f"\nChange relative to DC, same method (family - DC) at {cap}:")
        print(f"   {'family':<9}" + "".join(f"{SHORT[m]:>27}" for m in methods))
        for f in fams[1:]:
            row = ""
            for m in methods:
                d = pt[(f, m)] - pt[("DC", m)]
                v, lo, hi = ci(d, [b[(f, m)] - b[("DC", m)] for b in bs])
                row += f"  {fmt(v, lo, hi, j):>25}"
            print(f"   {f:<9}{row}")
        print(f"\nMethod contrasts within each family at {cap} (and the change of the contrast relative to DC):")
        for a_, b_ in CONTRASTS:
            print(f"   {SHORT[a_]} - {SHORT[b_]}")
            for f in fams:
                d = pt[(f, a_)] - pt[(f, b_)]
                v, lo, hi = ci(d, [b[(f, a_)] - b[(f, b_)] for b in bs])
                line = f"      {f:<9} {fmt(v, lo, hi, j)}"
                if f != "DC":
                    dd = d - (pt[("DC", a_)] - pt[("DC", b_)])
                    v2, lo2, hi2 = ci(dd, [(b[(f, a_)] - b[(f, b_)]) - (b[("DC", a_)] - b[("DC", b_)]) for b in bs])
                    line += f"    vs DC {fmt(v2, lo2, hi2, j)}"
                print(line)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", default=os.path.join(HERE, "c19_parts"))
    ap.add_argument("--B", type=int, default=1000)
    a = ap.parse_args()
    P = load(a.parts)
    print(f"C19, seed 999, exploratory. {len(P)} parts (network x cell x regime); "
          f"cells {sorted({o['cell'] for o in P})}, regimes {sorted({o['regime'] for o in P})}, "
          f"nets {min(o['net'] for o in P)}-{max(o['net'] for o in P)}; numpy {P[0].get('numpy')}; B = {a.B}")
    for kind in ("dense", "std"):
        report(P, kind, a.B)
    t = np.array([o["t_total"] for o in P])
    print(f"\nTime per part: median {np.median(t):.0f} s")


if __name__ == "__main__":
    main()
