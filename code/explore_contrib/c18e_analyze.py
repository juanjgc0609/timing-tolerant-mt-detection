"""c18e_analyze.py - exploration C18e: recall and paired contrasts (DENSE and STD windows, LOC and LOC+T, @0.1 and
@0.25 false events per window, cells pooled; hierarchical bootstrap of c10_analyze.Boot, B = 1000) for the statistics of
c18e_run.py on seed 999 with the qc1 noise of 2010-07-23.
    python c18e_analyze.py PARTSDIR OUT.txt"""
import os
import sys

import numpy as np

W = "/home/claude/work"
sys.path[:0] = [os.path.join(W, "explore_contrib"), W]
import c10_analyze as A  # noqa: E402
import c16_analyze as Z  # noqa: E402

PAIRS = [("MT@S0", "BOMP@S0"), ("MT-CAP24@S0", "BOMP@S0"), ("MT-LOG6@S0", "BOMP@S0"), ("HUBER-N@S0", "BOMP@S0"),
         ("MT-CAP24@S0", "MT@S0"), ("MT-LOG6@S0", "MT@S0"), ("MT-LOG6@S0", "HUBER-N@S0")]


def main(pdir, outp):
    nets, rows = Z.load(pdir, "REAL")
    names = list(rows[0]["atoms"])
    lines = [f"C18e: {len(nets)} REAL parts, seed 999, day 100723 qc1, networks {sorted({n['net'] for n in nets})}"]
    for des in ("DENSE", "STD"):
        for rule, tol in (("LOC", None), ("LOC+T", 0.5)):
            sub = [r for r in rows if r["design"] == des]
            bt = A.Boot(sub, Z.records(sub, names, tol), lambda r: True, 1000)
            lines.append(f"\n== {des} {rule}: {len(sub)} windows, {int(bt.ks.sum())} events; recall @0.1 / @0.25")
            for m in names:
                lines.append(f"   {m:12s} {bt.pt[m][0]:.3f}  {bt.pt[m][1]:.3f}")
            for a, b in PAIRS:
                d = bt.rep[a][:, 1] - bt.rep[b][:, 1]
                lines.append(f"   {a:12s} - {b:12s} @0.25 {bt.pt[a][1] - bt.pt[b][1]:+.3f} "
                             f"[{np.percentile(d, 2.5):+.3f}, {np.percentile(d, 97.5):+.3f}]")
    ku = np.array([r["kurt"] for r in rows if r["design"] == "DENSE"])
    lines.append(f"\nwindow excess kurtosis (dense) median {np.median(ku):.2f} p90 {np.percentile(ku, 90):.1f}")
    open(outp, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
