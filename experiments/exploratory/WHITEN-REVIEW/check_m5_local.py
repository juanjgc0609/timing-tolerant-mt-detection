"""WHITEN-REVIEW part 1, check of the M5 --local verdict (see card.md).
M5-EVALUE's --local step ran on the container's fake i.i.d. t4 noise, one scale at every station. Per window,
run_c.mfn divides each station's matched-filter output by its own robust level. The mismatch with the unweighted
operator is then diag(1/lev_i); it is inert when max/min lev is close to 1. This measures that ratio on i.i.d. t4
windows through the legacy kernel (window length NT, 3 components), for 4 to 8 stations.
    python check_m5_local.py > out_m5_local.txt
"""
import os
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "code"))
sys.path.insert(0, os.path.join(ROOT, "code", "explore_contrib"))
sys.path.insert(0, os.path.join(ROOT, "code", "r0"))
import r0_lib as R  # noqa: E402,I001  (must come before c8_lib: it sets VP, VS)
from ttmt.kernel import LEGACY  # noqa: E402
import c8_lib as C8  # noqa: E402  (NT of the 2M / M5 blocks)

SEED = 20261126
NT = int(C8.NT)
rng = np.random.default_rng(SEED)
print(f"seed {SEED}; i.i.d. Student t (4 dof), equal scale; {NT} samples x 3 components per station; 2000 windows")
for ns in (4, 6, 8):
    r = []
    for _ in range(2000):
        D = rng.standard_t(4, size=(ns, 3, NT))
        MF, _sd = LEGACY.matched(D)
        lev = np.median(np.abs(MF), axis=(1, 2)) / 0.6745
        r.append(lev.max() / lev.min())
    r = np.array(r)
    print(f"ns {ns}: max/min per-station level, median {np.median(r):.3f}, q95 {np.quantile(r, .95):.3f}, "
          f"max {r.max():.3f}")
print("reading fixed in the card: median <= 1.10 -> inert (real noise in FE-FIX: 5.7-7.6)")

# Added after reading the ratios above (descriptive, not part of the declared reading): what such a ratio costs.
# Noise-free f of an event in the unweighted basis when stations are divided by i.i.d.-t4 window levels, on the
# Ridgecrest N8 geometry of R0 (any 8-station geometry would do; this one is in the repository).
import task1x_mtglrt as X  # noqa: E402
STA, _ = R.stations(R.SUBNETS["N8"], R.frame())
rng = np.random.default_rng(SEED + 1)
f = []
for _ in range(500):
    src = np.array([rng.uniform(-15, 15), rng.uniform(-15, 15), rng.uniform(3, 12)])
    G = X.radiation_operator(src, STA)
    Q = X.basis(G)
    MF, _sd = LEGACY.matched(rng.standard_t(4, size=(len(STA), 3, NT)))
    lev = np.median(np.abs(MF), axis=(1, 2)) / 0.6745
    w = np.tile(np.repeat(1.0 / lev, 3), 2)
    v = w * (G @ rng.normal(size=6))
    f.append(float(np.sum((Q.T @ v) ** 2) / np.sum(v ** 2)))
f = np.array(f)
print(f"descriptive: noise-free f in the unweighted basis at those levels (N8, 500 draws): median {np.median(f):.4f}, "
      f"min {f.min():.4f}")
