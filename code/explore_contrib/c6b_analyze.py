"""C6b analysis: per regime and generator variant, recall @0.1/@0.25 of BOMP-PK/MT-FW with the homogeneous (HOM) and
the variant's true-medium (ORA) model, CEIL-GD, headroom CEIL-GD - BOMP-PK:HOM, and paired differences between each
variant and DEF (the variants share sources, mechanisms, statics, clutter and noise).  Hierarchical bootstrap
(cells as strata, networks then windows, the same resampled (network, window) set for all variants), 1000 replicates.
    python c6b_analyze.py out/c6b_thin.json out/c6b_strong.json"""
import json
import os
import sys

import numpy as np

import xc

VARS = ("DEF", "MED", "PULSE", "SITE", "CODA")
METH = ("BOMP-PK:HOM", "MT-FW:HOM", "BOMP-PK:ORA", "MT-FW:ORA", "CEIL-GD")
B = int(os.environ.get("NBOOT", 1000))


def main():
    out = []
    for path in sys.argv[1:]:
        rows = json.load(open(path))["rows"]
        reg = os.path.basename(path).split("_")[1].split(".")[0]
        by = {}
        for r in rows:
            by.setdefault(r["var"], {})[(r["cell"], r["net"], r["w"])] = r
        vars_ = [v for v in VARS if v in by]
        keys = sorted(set.intersection(*[set(by[v]) for v in vars_]))
        for v in vars_:
            assert all(by[v][k]["k"] == by["DEF"][k]["k"] for k in keys)
        ks = [by["DEF"][k]["k"] for k in keys]
        recs = {(v, m): [[tuple(x) for x in by[v][k]["rec"][m]] for k in keys] for v in vars_ for m in METH}
        nets = sorted({(k[0], k[1]) for k in keys})
        idx = {n: [i for i, k in enumerate(keys) if (k[0], k[1]) == n] for n in nets}
        cells = sorted({n[0] for n in nets})

        def rc(ws):
            n = sum(ks[w] for w in ws)
            return {vm: np.array(xc.recall_caps([recs[vm][w] for w in ws], n)[0]) for vm in recs}
        pt = rc(list(range(len(keys))))
        rng = np.random.default_rng(21)
        reps = {vm: [] for vm in recs}
        for _ in range(B):
            ws = []
            for c in cells:
                nl = [n for n in nets if n[0] == c]
                for g in rng.integers(0, len(nl), len(nl)):
                    gg = idx[nl[g]]
                    ws.extend(np.asarray(gg)[rng.integers(0, len(gg), len(gg))].tolist())
            r_ = rc(ws)
            for vm in recs:
                reps[vm].append(r_[vm])
        reps = {vm: np.array(a) for vm, a in reps.items()}

        def fmt(p, a):
            lo, hi = np.percentile(a, [2.5, 97.5], axis=0)
            return f"{p[0]:+.3f} [{lo[0]:+.3f}, {hi[0]:+.3f}] / {p[1]:+.3f} [{lo[1]:+.3f}, {hi[1]:+.3f}]"
        out.append("=" * 118)
        out.append(f"C6b {reg}: {len(keys)} windows x {len(vars_)} variants, {sum(ks)} events, {len(nets)} networks")
        out.append("  recall @0.1 / @0.25 (point):")
        for v in vars_:
            out.append(f"    {v:6s} " + "  ".join(f"{m} {pt[(v, m)][0]:.3f}/{pt[(v, m)][1]:.3f}" for m in METH))
        out.append("  headroom CEIL-GD - BOMP-PK:HOM  and  CEIL-GD - MT-FW:HOM   (@0.1 / @0.25, 95% CI)")
        for v in vars_:
            out.append(f"    {v:6s} {fmt(pt[(v, 'CEIL-GD')] - pt[(v, 'BOMP-PK:HOM')], reps[(v, 'CEIL-GD')] - reps[(v, 'BOMP-PK:HOM')])}"
                       f"    {fmt(pt[(v, 'CEIL-GD')] - pt[(v, 'MT-FW:HOM')], reps[(v, 'CEIL-GD')] - reps[(v, 'MT-FW:HOM')])}")
        out.append("  variant - DEF (paired), @0.1 / @0.25:")
        for v in vars_[1:]:
            for m in ("BOMP-PK:HOM", "MT-FW:HOM", "CEIL-GD"):
                out.append(f"    {v:6s} {m:12s} {fmt(pt[(v, m)] - pt[('DEF', m)], reps[(v, m)] - reps[('DEF', m)])}")
        out.append("  ORA - HOM within variant (travel-time/ray model correction), @0.1 / @0.25:")
        for v in vars_:
            for d in ("BOMP-PK", "MT-FW"):
                a, b = (v, d + ":ORA"), (v, d + ":HOM")
                out.append(f"    {v:6s} {d:8s} {fmt(pt[a] - pt[b], reps[a] - reps[b])}")
        out.append("  atom coverage (event has an atom within 5 km, any score) BOMP-PK:HOM / MT-FW:HOM:")
        for v in vars_:
            s = []
            for m in ("BOMP-PK:HOM", "MT-FW:HOM"):
                cov = []
                for k in keys:
                    r = by[v][k]
                    Ls = [np.asarray(l[:2]) for l in r["loc"][m]]
                    for t in r["truth"]:
                        cov.append(min([np.linalg.norm(l - np.asarray(t[:2])) for l in Ls] or [np.inf]) <= 5.0)
                s.append(f"{np.mean(cov):.3f}")
            out.append(f"    {v:6s} " + " / ".join(s))
    print("\n".join(out))


if __name__ == "__main__":
    main()
