"""C6 diagnostic: size of the travel-time mismatch (gradient medium - homogeneous detector model) that a per
station/phase constant (a static correction) cannot absorb, over sources uniform in [-15, 15]^2 x [5, 13] km,
on the C6 networks (seed 999, j = 20..29)."""
import numpy as np
import c6_lib as L
import xc
rng = np.random.default_rng(0)
SRC = np.stack([rng.uniform(-15, 15, 4000), rng.uniform(-15, 15, 4000), rng.uniform(5, 13, 4000)], 1)
H, G = L.Homog(), L.Gradient()
for cell in ("ladder_4", "ladder_6", "unilateral_c40"):
    raw, res, p95 = [], [], []
    for j in range(20, 30):
        STA, _st = xc.T.network(999, cell, j)
        d = G.tt(SRC, STA) - H.tt(SRC, STA)                  # (N, 2, ns)
        raw.append(np.abs(d).mean())
        r = d - np.median(d, axis=0, keepdims=True)
        res.append(np.abs(r).mean())
        p95.append(np.percentile(np.abs(r), 95))
    print(f"{cell:15s} mean |tt_grad - tt_hom| {np.mean(raw):.3f} s; after a per-station/phase constant: mean {np.mean(res):.3f} s, P95 {np.mean(p95):.3f} s")
