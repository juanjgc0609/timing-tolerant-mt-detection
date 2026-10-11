"""generic analysis: rows with 'cell', 'net', 'k', 'rec' {method: record}. Recall at 0.1/0.25 with hierarchical
bootstrap (networks resampled within cell, then windows), paired differences."""
import json, sys
import numpy as np
import xc


def load(path, filt=None):
    rows = json.load(open(path))
    if filt:
        rows = [r for r in rows if filt(r)]
    return rows


def table(rows, methods, pairs, B=int(__import__("os").environ.get("NBOOT", 1000)), seed=7, by_cell=True, title=""):
    out = []
    keys = sorted({(r["cell"], r["net"]) for r in rows})
    idx = {k: [] for k in keys}
    for i, r in enumerate(rows):
        idx[(r["cell"], r["net"])].append(i)
    groups = [idx[k] for k in keys]
    strata = [k[0] for k in keys]
    ks = [r["k"] for r in rows]
    recs = {m: [[tuple(x) for x in r["rec"][m]] for r in rows] for m in methods}
    pt, ci, dif = xc.hier_boot(groups, methods, recs, ks, B=B, seed=seed, pairs=pairs, strata=strata)
    out.append(f"{title}  windows {len(rows)}, events {sum(ks)}, networks {len(keys)}")
    for m in methods:
        out.append(f"  {m:28s} recall@0.1 {pt[m][0]:.3f} [{ci[m][0][0]:.3f},{ci[m][0][1]:.3f}]   "
                   f"@0.25 {pt[m][1]:.3f} [{ci[m][1][0]:.3f},{ci[m][1][1]:.3f}]")
    for a, b in pairs:
        d = dif[(a, b)]
        out.append(f"  {a} - {b}: @0.1 {xc.fmt_ci(d['point'][0], d['ci'][0])}   @0.25 {xc.fmt_ci(d['point'][1], d['ci'][1])}")
    if by_cell:
        for cell in sorted({r["cell"] for r in rows}):
            sub = [r for r in rows if r["cell"] == cell]
            n = sum(r["k"] for r in sub)
            s = "    " + cell + ": " + "  ".join(
                f"{m.split(':')[-1] if len(methods) < 0 else m} {'/'.join(f'{v:.2f}' for v in xc.recall_caps([[tuple(x) for x in r['rec'][m]] for r in sub], n)[0])}"
                for m in methods)
            out.append(s)
    return "\n".join(out), pt, ci, dif
