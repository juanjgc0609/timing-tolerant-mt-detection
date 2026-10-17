"""c18f_analyze.py - exploration C18f: the methods of C18e and C18f on the same windows (rows paired by cell, network,
design and window index; truths checked equal), recall and paired contrasts at @0.1 / @0.25 (c10_analyze.Boot,
B = 1000), DENSE and STD windows, LOC and LOC+T.
    python c18f_analyze.py C18E_PARTS C18F_PARTS OUT.txt"""
import os
import sys

import numpy as np

W = "/home/claude/work"
sys.path[:0] = [os.path.join(W, "explore_contrib"), W]
import c10_analyze as A  # noqa: E402
import c16_analyze as Z  # noqa: E402

PAIRS = [("BOMP-CAP6@S0", "BOMP@S0"), ("BOMP-STA@S0", "BOMP@S0"), ("MT@S0", "BOMP-CAP6@S0"), ("MT@S0", "BOMP-STA@S0"),
         ("HUBER-N@S0", "BOMP-CAP6@S0"), ("HUBER-N@S0", "BOMP-STA@S0"), ("MT-LOG6@S0", "BOMP-CAP6@S0")]


def main(pe, pf, outp):
    _ne, re_ = Z.load(pe, "REAL")
    _nf, rf = Z.load(pf, "REAL")
    key = {(r["cell"], r["net"], r["design"], r["w"]): r for r in rf}
    rows = []
    for r in re_:
        q = key.get((r["cell"], r["net"], r["design"], r["w"]))
        if q is None:
            continue
        assert r["k"] == q["k"] and np.allclose(np.array(r["truth"]).reshape(-1), np.array(q["truth"]).reshape(-1))
        r = dict(r)
        r["atoms"] = dict(r["atoms"], **q["atoms"])
        rows.append(r)
    names = list(rows[0]["atoms"])
    lines = [f"C18f: {len(rows)} paired rows (of {len(re_)} C18e rows), methods {names}"]
    for des in ("DENSE", "STD"):
        for rule, tol in (("LOC", None), ("LOC+T", 0.5)):
            sub = [r for r in rows if r["design"] == des]
            bt = A.Boot(sub, Z.records(sub, names, tol), lambda r: True, 1000)
            lines.append(f"\n== {des} {rule}: {len(sub)} windows, {int(bt.ks.sum())} events; recall @0.1 / @0.25")
            for m in names:
                lines.append(f"   {m:13s} {bt.pt[m][0]:.3f}  {bt.pt[m][1]:.3f}")
            for a, b in PAIRS:
                d = bt.rep[a][:, 1] - bt.rep[b][:, 1]
                lines.append(f"   {a:13s} - {b:13s} @0.25 {bt.pt[a][1] - bt.pt[b][1]:+.3f} "
                             f"[{np.percentile(d, 2.5):+.3f}, {np.percentile(d, 97.5):+.3f}]")
    open(outp, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
