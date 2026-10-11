"""Task 2M power planning from C16's evaluation parts (seed 999, REAL arm; exploratory data already reported).
    python t2m_power.py OUT.txt
Hierarchical-bootstrap sd of the G1-R contrasts (DENSE, LOC, @0.25) with (a) all 12 networks and 12 windows per network
and cell, (b) the first 6 windows of every network, (c) the first 6 networks of every cell; variance components
var = s_net^2 / N + s_win^2 / (N n) fitted to (a) and (b), checked against (c); projected sd for N networks per cell and n
DENSE windows per network; power of "2.5th percentile > 0 and point >= 0.03" per seed and in both seeds, with the
lower-tail inflation of C16's percentile intervals ((point - lower) / 1.96 / sd)."""
import os
import sys

import numpy as np
from scipy.stats import norm

WORK = os.path.dirname(os.path.abspath(__file__))
EXPL = os.path.join(WORK, "explore_contrib")
for p in (EXPL, WORK):
    if p not in sys.path:
        sys.path.insert(0, p)
import c10_analyze as A  # noqa: E402
import c16_analyze as Z  # noqa: E402

PAIRS = [("MT@S0", "BOMP@S0"), ("HUBER-N@S0", "BOMP@S0"), ("HUBER-N@S0", "MT@S0")]
B = 1000


def sds(rows, names):
    recs = Z.records(rows, names, None)
    bt = A.Boot(rows, recs, lambda r: True, B)
    out = {}
    for a, b in PAIRS:
        d = bt.rep[a][:, 1] - bt.rep[b][:, 1]
        p = bt.pt[a][1] - bt.pt[b][1]
        out[(a, b)] = (p, d.std(), np.percentile(d, 2.5))
    return out, len(bt.idx)


def main(outp):
    lines = []

    def P(s=""):
        print(s, flush=True)
        lines.append(s)
    _n, rows = Z.load(os.path.join(EXPL, "out", "c16_eval_parts"), "REAL")
    dense = [r for r in rows if r["design"] == "DENSE"]
    names = list(dense[0]["atoms"])
    full, nf = sds(dense, names)
    half_w, nw = sds([r for r in dense if r["w"] < 6], names)
    nets = sorted({(r["cell"], r["net"]) for r in dense})
    keep = {k for c in Z.CELLS for k in sorted(x for x in nets if x[0] == c)[:6]}
    half_n, nn = sds([r for r in dense if (r["cell"], r["net"]) in keep], names)
    P(f"C16 REAL DENSE LOC @0.25, B = {B}: (a) {nf} windows (12 nets x 12 windows per cell), (b) {nw} (12 x 6), "
      f"(c) {nn} (6 x 12)")
    comp = {}
    for pr in PAIRS:
        (p, s12, lo), (_p6, s6, _l6), (_pn, sn, _ln) = full[pr], half_w[pr], half_n[pr]
        N, n = 36, 12
        bw = max((s6 ** 2 - s12 ** 2) * N / (1 / 6 - 1 / 12), 0.0)
        an = max(N * (s12 ** 2 - bw / (N * n)), 0.0)
        pred_c = np.sqrt(an / 18 + bw / (18 * 12))
        infl = (p - lo) / 1.96 / s12
        comp[pr] = (an, bw, infl, p)
        P(f"  {pr[0]} - {pr[1]}: point {p:+.3f}; sd (a) {s12:.4f} (b) {s6:.4f} (c) {sn:.4f}; s_net^2 {an:.4f} "
          f"s_win^2 {bw:.4f} (share of var(a) from networks {an / N / s12 ** 2:.2f}); predicted (c) {pred_c:.4f}; "
          f"lower-tail inflation {infl:.2f}")
    P("")
    P("Projected sd and power per seed / both seeds (rule: 2.5th percentile > 0 and point >= 0.03), N networks per cell "
      "(3 cells), n DENSE windows per network; effects: C16 REAL point, half of it, C16 GEN point (0.069 / 0.113)")
    eff = {("MT@S0", "BOMP@S0"): (0.259, 0.130, 0.069), ("HUBER-N@S0", "BOMP@S0"): (0.294, 0.147, 0.113)}
    for N in (12, 16, 20, 24, 28, 32):
        for n in (8, 12, 16):
            s = []
            for pr in PAIRS[:2]:
                an, bw, infl, _p = comp[pr]
                sd = np.sqrt(an / (3 * N) + bw / (3 * N * n)) * infl
                pw = []
                for e in eff[pr]:
                    p1 = norm.cdf(e / sd - 1.96) * norm.cdf((e - 0.03) / (sd / infl))
                    pw.append((p1, p1 ** 2))
                s.append(f"{pr[0][:7]}: sd_eff {sd:.4f} power {' / '.join(f'{a:.2f}|{b:.2f}' for a, b in pw)}")
            P(f"  N {N:2d} n {n:2d}  " + "   ".join(s))
    with open(outp, "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main(sys.argv[1])
