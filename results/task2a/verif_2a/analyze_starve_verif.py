import json, sys
import numpy as np
sys.path.insert(0, "/home/claude/work")
from froc_exact import exact_HF
from task1r_fair import recall_at
from task1r_baselines import GAMMA_GRID
R = json.load(open(sys.argv[1]))
print("consistency with stored run:", sum(R["consistency"]), "/", len(R["consistency"]))
pe = R["per_event"]
print("true events:", len(pe))
for key in ("matched_lab_0.3", "matched_anylabel_0.3", "matched_anylabel_1.0"):
    v = np.array([e[key] for e in pe]); n = np.array([e["n_arr"] for e in pe])
    print(f"{key}: mean per event {v.mean():.2f} of {n.mean():.1f} arrivals ({v.sum() / n.sum():.3f}); "
          + ", ".join(f">={m}: {np.mean(v >= m):.2f}" for m in (3, 4, 5, 6, 8)))
errs = np.concatenate([e["errs"] for e in pe if e["errs"]])
print(f"matched pick error (s): median {np.median(errs):+.3f}, mean {errs.mean():+.3f}, sd {errs.std():.3f}, n {len(errs)}")
for cell in sorted({e["cell"] for e in pe}):
    v = np.array([e["matched_lab_0.3"] for e in pe if e["cell"] == cell])
    print(f"   {cell}: matched picks per event mean {v.mean():.2f}, frac events >=3: {np.mean(v >= 3):.2f}, >=4: {np.mean(v >= 4):.2f}")
W = R["windows"]
k = [w["k"] for w in W]; K = sum(k); n = len(W)
def per(key):
    out = []
    for cap in (0.1, 0.25):
        best = 0.0
        for c in range(len(W[0][key])):
            H = sum(w[key][c][0] for w in W); F = sum(w[key][c][1] for w in W)
            if F / n <= cap + 1e-12: best = max(best, H / K)
        out.append(best)
    return out
def rec(key):
    _t, H, F = exact_HF([[tuple(x) for x in w[key]] for w in W])
    return [recall_at(H, F, K, c) for c in (0.1, 0.25)]
print(f"windows {n}, events {K}")
for key in ("G_true", "G_true_orSC", "G_stalta_hitsonly", "P_true", "G_DY_none", "G_DY_est", "G_DY_static", "G_DY_static+mis"):
    r = per(key)
    best_any = max(sum(w[key][c][0] for w in W) for c in range(len(W[0][key]))) / K
    print(f"{key:20s} recall@0.1 {r[0]:.2f}  @0.25 {r[1]:.2f}  (max recall any config {best_any:.2f})")
for key in ("MTFW_none", "MTFW_est", "MTFW_static", "MTFW_static+mis", "QM3C_none", "QM3C_est", "QM3C_static", "QM3C_static+mis"):
    r = rec(key)
    print(f"{key:20s} recall@0.1 {r[0]:.2f}  @0.25 {r[1]:.2f}")
