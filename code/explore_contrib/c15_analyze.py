"""C15 analysis.   python c15_analyze.py dev|eval [B]  -> out/c15_analysis_<mode>.txt
Uses the FROC / bootstrap machinery of c14_analyze (exact FROC, operating point re-chosen per replicate, hierarchical
bootstrap: strata regime x cell, networks then windows).  Paired differences against BASE (the registered report)."""
import glob
import json
import os
import sys

import numpy as np

import c14_analyze as A

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def load(mode):
    rows = []
    for f in sorted(glob.glob(os.path.join(OUT, f"c15_{mode}_parts", "*.json"))):
        p = json.load(open(f))
        for r in p["rows"]:
            r.update(regime=p["regime"], cell=p["cell"], net=p["net"])
            rows.append(r)
    return rows


def timing(rows, m):
    """origin error of atoms within 5 km of a true event (nearest), all scores."""
    e = []
    for r in rows:
        T = np.array(r["truth"])
        if not len(T):
            continue
        for a in r["atoms"][m]:
            d = np.linalg.norm(T[:, :2] - np.array(a[:2]), axis=1)
            i = int(d.argmin())
            if d[i] <= 5.0:
                e.append(a[3] - r["taus"][i])
    return np.array(e)


def main(mode, B):
    rows = load(mode)
    methods = list(rows[0]["atoms"].keys())
    fh = open(os.path.join(OUT, f"c15_analysis_{mode}.txt"), "w")
    nets = sorted({r["net"] for r in rows})
    fh.write(f"C15 analysis, mode {mode}: {len(rows)} windows, networks {nets}, methods {methods}, B = {B}\n")
    pl = np.array([p[0] for r in rows for p in r["plateau"]])
    fh.write(f"plateau size of the accepted atoms (origins with the maximal statistic): mean {pl.mean():.3f}, "
             f"fraction > 1: {(pl > 1).mean():.3f}\n")
    for design in ("STD", "DENSE"):
        rs = [r for r in rows if r["design"] == design]
        fh.write("\n" + "=" * 110 + f"\n{design}: {len(rs)} windows, {sum(r['k'] for r in rs)} events\n")
        for m in methods:
            e = timing(rs, m)
            fh.write(f"  origin error of atoms within 5 km ({m:5s}): n {len(e)}, median {np.median(e):+.3f} s, "
                     f"|e| <= 0.5: {np.mean(np.abs(e) <= 0.5):.3f}, |e| median {np.median(np.abs(e)):.3f}, "
                     f"early {np.mean(e < 0):.3f}\n")
        T = A.Tab(rs, methods)
        sel = np.ones(len(rs), bool)
        pt, reps = A.boot(T, sel, methods, B, 20261009)
        for u in A.RULES:
            for c in A.CAPS:
                fh.write(f"  {u:5s} <= {c}: ")
                for m in methods:
                    k = (m, u, c)
                    fh.write(f"{m} {pt[k]:.3f}  ")
                    if m != "BASE":
                        d = reps[k] - reps[("BASE", u, c)]
                        fh.write(f"(d {A.fmt(pt[k] - pt[('BASE', u, c)], A.ci(d))})  ")
                fh.write("\n")
        for m in methods:
            k = (m, "a1", 0.25)
            fh.write(f"  right place wrong time above the <= 0.25 LOC+T threshold, per window: {m} {pt[k]:.3f}; "
                     f"all scores {pt[(m, 'a1', 'all')]:.3f}\n")
        for m in methods:
            fh.write(f"  close pairs (10 km, 2 s) LOC+T <= 0.25: {m} {pt[(m, 'LOC+T', 0.25, 'close')]:.3f}"
                     f"   tight (7 km, 1 s): {pt[(m, 'LOC+T', 0.25, 'tight')]:.3f}\n")
        # per cell
        for cell in A.CELLS:
            s2 = np.array([r["cell"] == cell for r in rs])
            if not s2.any():
                continue
            pc, rc = A.boot(T, s2, methods, max(100, B // 2), 7)
            line = f"  cell {cell:15s} LOC+T <= 0.25: BASE {pc[('BASE', 'LOC+T', 0.25)]:.3f}"
            for m in methods[1:]:
                d = rc[(m, "LOC+T", 0.25)] - rc[("BASE", "LOC+T", 0.25)]
                line += f"  {m} {A.fmt(pc[(m, 'LOC+T', 0.25)] - pc[('BASE', 'LOC+T', 0.25)], A.ci(d))}"
            fh.write(line + "\n")
    fh.close()
    print(open(os.path.join(OUT, f"c15_analysis_{mode}.txt")).read())


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 400)
