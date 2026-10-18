"""power_mcnemar.py - number of targets a paired recall test needs (exploration P6-power).

Two methods score the same N catalogued targets; each target is hit or missed by each method. The decision rule is
the one of Task R1: an exact one-sided McNemar test on the discordant pairs at alpha = 0.05. Hits are drawn from a
Gaussian copula: target i is hit by method k when z_ik < Phi^-1(p_k), with corr(z_iA, z_iB) = rho. Targets are
independent here, so the numbers are a lower bound; with targets clustered in time blocks multiply N by the design
effect 1 + (m - 1) * ICC.

Usage: python3 power_mcnemar.py  (numpy only; about a minute on two cores). Output: power_mcnemar_out.txt next to it.
"""
from fractions import Fraction
from math import comb
from pathlib import Path
import numpy as np
from statistics import NormalDist

ALPHA = 0.05
P_BASE = 0.40          # recall of the comparator, as in R1 (BASE 0.407, BESTPICK 0.481)
DELTAS = (0.03, 0.05, 0.10, 0.15)
RHOS = (0.5, 0.7, 0.9)
NS = (54, 100, 150, 200, 300, 400, 600, 800, 1200, 1600, 2400)
SIMS = 4000
SEED = 20261018
OUT = Path(__file__).with_name("power_mcnemar_out.txt")


def critical(nmax):
    """crit[n] = smallest k with P(X >= k | X ~ Bin(n, 1/2)) <= ALPHA (n + 1 if none)."""
    crit = np.empty(nmax + 1, dtype=np.int64)
    for n in range(nmax + 1):
        a, tail, k = Fraction(ALPHA).limit_denominator(1000), 0, n + 1
        for j in range(n, -1, -1):
            tail += comb(n, j)                # integer arithmetic: 2**n overflows a float beyond n = 1023
            if tail > a * 2 ** n:
                break
            k = j
        crit[n] = k
    return crit


def power(n, p_a, p_b, rho, crit, rng):
    za = rng.standard_normal((SIMS, n))
    zb = rho * za + np.sqrt(1 - rho ** 2) * rng.standard_normal((SIMS, n))
    ha, hb = za < NormalDist().inv_cdf(p_a), zb < NormalDist().inv_cdf(p_b)
    b = (ha & ~hb).sum(1)               # A hits, B misses
    c = (~ha & hb).sum(1)
    return float(np.mean(b >= crit[b + c]))


def main():
    rng = np.random.default_rng(SEED)
    crit = critical(max(NS))
    lines = [f"exact one-sided McNemar, alpha {ALPHA}, comparator recall {P_BASE}, {SIMS} simulations, seed {SEED}", ""]
    lines.append("false-positive rate at delta 0 (should be <= alpha):")
    for rho in RHOS:
        lines.append(f"  rho {rho}: " + "  ".join(f"N={n}:{power(n, P_BASE, P_BASE, rho, crit, rng):.3f}" for n in (54, 300, 1200)))
    lines.append("")
    for rho in RHOS:
        lines.append(f"rho {rho} (correlation of the latent hit variables of the two methods)")
        lines.append("  delta  " + "".join(f"{n:>7d}" for n in NS) + "   N for 80 %")
        for d in DELTAS:
            pw = [power(n, P_BASE + d, P_BASE, rho, crit, rng) for n in NS]
            n80 = next((n for n, p in zip(NS, pw) if p >= 0.8), None)
            lines.append(f"  {d:+.2f}  " + "".join(f"{p:7.2f}" for p in pw) + f"   {n80 if n80 else '> ' + str(max(NS))}")
        lines.append("")
    OUT.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
