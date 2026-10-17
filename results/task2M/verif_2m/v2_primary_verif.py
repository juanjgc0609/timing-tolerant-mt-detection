"""V2: independent recomputation of the 2M endpoints (point estimates and hierarchical bootstrap)."""
import sys
import time

import numpy as np

sys.path.insert(0, "/home/claude/work/t2m_verif")
import vlib_verif as V  # noqa: E402

B = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
OUT = "/home/claude/work/t2m_verif/v2_primary_verif.txt"
M = ["MT@S0", "HUBER-N@S0", "BOMP@S0", "MT@S0+T3", "MT@S2+T3"]
PAIRS = [("MT@S0", "BOMP@S0"), ("HUBER-N@S0", "BOMP@S0"), ("HUBER-N@S0", "MT@S0")]
lines = []


def P(s=""):
    print(s, flush=True)
    lines.append(s)
    with open(OUT, "w") as fh:
        fh.write("\n".join(lines) + "\n")


P(f"V2 independent recomputation; own bootstrap seed 20261016, B = {B}")
for seed in ("2718", "3141"):
    nets = V.load_parts(seed)
    for des, rule, tol in (("DENSE", "LOC", None), ("DENSE", "LOC+T", 0.5), ("STD", "LOC", None)):
        wins = V.windows(nets, des)
        prep = {m: V.prep(wins, m, 5.0, tol) for m in M}
        allw = np.arange(len(wins))
        pt = {m: V.recall_at(*prep[m], allw) for m in M}
        P("")
        P(f"seed {seed} {des} {rule}: {len(wins)} windows, {prep[M[0]][2].sum()} events")
        P("   " + "  ".join(f"{m} {pt[m][0]:.3f}/{pt[m][1]:.3f}" for m in M) + "   (@0.1/@0.25)")
        if not (des == "DENSE" and rule == "LOC") and seed == "3141":
            Bx = 500
        else:
            Bx = B if (des == "DENSE" and rule == "LOC") else 1000
        t0 = time.time()
        reps = V.boot_index(wins, Bx, 20261016)
        R = {m: np.array([V.recall_at(*prep[m], ws) for ws in reps]) for m in M[:3]}
        for a, b in PAIRS:
            d = R[a][:, 1] - R[b][:, 1]
            p = pt[a][1] - pt[b][1]
            P(f"   {a:11s} - {b:9s} @0.25 {p:+.4f} [{np.percentile(d, 2.5):+.4f}, {np.percentile(d, 97.5):+.4f}] "
              f"sd {d.std():.4f} p(d<=0) {(1 + np.sum(d <= 0)) / (Bx + 1):.4f}  p(d<=-0.03) "
              f"{(1 + np.sum(d <= -0.03)) / (Bx + 1):.4f}   @0.1 {pt[a][0] - pt[b][0]:+.4f}")
        P(f"   recall CIs @0.25: " + "  ".join(
            f"{m} [{np.percentile(R[m][:, 1], 2.5):.3f},{np.percentile(R[m][:, 1], 97.5):.3f}]" for m in M[:3]) +
          f"   (B {Bx}, {time.time() - t0:.0f} s)")
