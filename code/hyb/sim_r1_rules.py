"""
sim_r1_rules.py - operating characteristics of candidate decision rules for Task R1 (registered real-data test on the
Hi-net swarm), by simulation, before any Hi-net waveform is read.

Model. N target events in EVAL. Event i has a latent difficulty z_i ~ N(0, 1); method m detects it when
z_i + e_mi > c_m, with e_mi ~ N(0, s^2) independent between methods (s sets how much the two methods disagree on the same
events: small s = highly correlated detections).  c_m is set so that the recall of the comparator is pB and that of the
candidate is pB + delta.  The operating point (threshold at 30 false detections per hour) is taken as given: the
simulation concerns the sampling of target events, which dominates with tens of events.

Rules evaluated (all one-sided at 0.025):
  BOOT0     percentile bootstrap over events of the paired recall difference, 2.5th percentile > 0
  MCN       exact McNemar (sign test on discordant events) p <= 0.025, candidate better
  MCN+PT    MCN and point difference >= +0.05
  NI(d)     bootstrap 2.5th percentile > -d (non-inferiority with margin d), d = 0.10, 0.20
Output: probability of meeting each rule for each (N, pB, delta, s): at delta = 0 this is the false-positive rate of the
superiority rules; at delta = -d the false-positive rate of NI(d); elsewhere the power.
"""
import itertools
import sys

import numpy as np
from scipy.stats import binom, norm

R = np.random.default_rng(20261017)
NSIM, B = 4000, 1000


def thresholds(pB, pH, s):
    """c such that P(z + e > c) = p with z ~ N(0,1), e ~ N(0, s^2)."""
    sd = np.sqrt(1 + s * s)
    return norm.ppf(1 - pB) * sd, norm.ppf(1 - pH) * sd


def sim(N, pB, delta, s):
    pH = min(max(pB + delta, 0.0), 1.0)
    cB, cH = thresholds(pB, pH, s)
    z = R.normal(size=(NSIM, N))
    hB = (z + R.normal(scale=s, size=(NSIM, N))) > cB
    hH = (z + R.normal(scale=s, size=(NSIM, N))) > cH
    d = hH.astype(float) - hB.astype(float)                       # per event paired difference
    pt = d.mean(1)
    idx = R.integers(0, N, size=(B, N))
    lo = np.empty(NSIM)
    for k in range(NSIM):
        lo[k] = np.percentile(d[k][idx].mean(1), 2.5)
    n10 = (hH & ~hB).sum(1)
    n01 = (~hH & hB).sum(1)
    nd = n10 + n01
    p_mcn = np.where(nd > 0, binom.sf(n10 - 1, np.maximum(nd, 1), 0.5), 1.0)
    return {"BOOT0": np.mean(lo > 0), "MCN": np.mean(p_mcn <= 0.025), "MCN+PT": np.mean((p_mcn <= 0.025) & (pt >= 0.05)),
            "NI(0.10)": np.mean(lo > -0.10), "NI(0.20)": np.mean(lo > -0.20), "mean_disc": nd.mean()}


def main(outp):
    lines = ["N  pB   s    delta  " + "  ".join(f"{k:>8s}" for k in ("BOOT0", "MCN", "MCN+PT", "NI(0.10)", "NI(0.20)")) + "  discordant"]
    for N, pB, s in itertools.product((9, 14, 40, 100), (0.3, 0.6), (0.5, 1.0)):
        for delta in (-0.20, -0.10, 0.0, 0.10, 0.20, 0.30):
            r = sim(N, pB, delta, s)
            lines.append(f"{N:3d} {pB:.1f} {s:.1f} {delta:+.2f}  " + "  ".join(f"{r[k]:8.3f}" for k in ("BOOT0", "MCN", "MCN+PT", "NI(0.10)", "NI(0.20)"))
                         + f"  {r['mean_disc']:6.1f}")
        print("\n".join(lines[-6:]), flush=True)
    open(outp, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main(sys.argv[1])
