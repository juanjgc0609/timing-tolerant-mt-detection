"""M6-SPURIOUS reading of amendment 1 (card.md) on rows.json.    python analyse.py > analyse_out.txt"""
import json
from pathlib import Path

import numpy as np
from scipy.stats import beta, fisher_exact

from ttmt.stats import exact_mcnemar

HERE = Path(__file__).resolve().parent
rows = json.load(open(HERE / "rows.json"))
N = 1 + max(r["i"] for r in rows)
pos = {(r["kind"], r["i"], r["mis"]): (r["gain"] is not None and r["gain"] > 0) for r in rows}


def cp(k, n):
    return (beta.ppf(.025, k, n - k + 1) if k else 0.0, beta.ppf(.975, k + 1, n - k) if k < n else 1.0)


K = {(k, m): sum(pos[k, i, m] for i in range(N)) for k in ("C1", "C3") for m in (True, False)}
for (k, m), v in K.items():
    lo, hi = cp(v, N)
    print(f"{k} {'mismatched' if m else 'exact'}: {v}/{N} = {v / N:.3f} [{lo:.3f}, {hi:.3f}]")
b = sum(pos["C3", i, True] and not pos["C3", i, False] for i in range(N))
c = sum(pos["C3", i, False] and not pos["C3", i, True] for i in range(N))
p_mc = exact_mcnemar(b, c)
p_f = fisher_exact([[K["C1", True], N - K["C1", True]], [K["C3", True], N - K["C3", True]]], alternative="less")[1]
pos3 = [r for r in rows if r["kind"] == "C3" and r["mis"] and r["gain"] is not None and r["gain"] > 0]
near = []
for r in pos3:
    cost = [max(d / 5.0, abs(t) / 1.5) for d, t in zip(r["d"], r["dt"])]
    near.append(int(np.argmin(cost)))
share = np.bincount(near, minlength=3)
lo, hi = cp(int(share[1]), max(len(pos3), 1))
print(f"C3 paired McNemar: only mismatched {b}, only exact {c}, p {p_mc:.2e}")
print(f"Fisher C1 vs C3 (mismatched, one-sided): p {p_f:.2e}")
print(f"positive C3 mismatched, nearest target / neighbour / isolated: {share.tolist()} (neighbour share "
      f"{share[1] / max(len(pos3), 1):.2f} [{lo:.2f}, {hi:.2f}])")
pm, pe, pc1 = K["C3", True] / N, K["C3", False] / N, K["C1", True] / N
exact_removes = pm > 0 and pe <= 0.2 * pm and p_mc < 0.01
c1_small = pm > 0 and pc1 <= 0.2 * pm and p_f < 0.01
if pm == 0:
    v = "no tail reproduced"
elif exact_removes and c1_small:
    v = "H3 (mismatch x neighbour)"
elif exact_removes:
    v = "H1 (pulse mismatch)"
elif c1_small and len(pos3) and share[1] / len(pos3) >= 0.60:
    v = "H2 (neighbour)"
else:
    v = "no clear cause"
print(f"reading (amendment 1): {v}")
g = {m: np.array([r["gain"] for r in rows if r["kind"] == "C3" and r["mis"] == m and r["gain"] is not None])
     for m in (True, False)}
for m, a in g.items():
    if a.size:
        print(f"C3 {'mismatched' if m else 'exact'} gains: q95 {np.quantile(a, .95):.1f}, max {a.max():.1f}, n {a.size}")
