"""Exploratory: network-level (cluster) bootstrap of MT-FW+SC minus the best of the 7 preregistered comparators
(re-selected per replicate), 2000 replicates, resampling the 5 networks of a cell with replacement."""
import json, sys
sys.path.insert(0, "/home/claude/work")
import numpy as np
import analyze_2a as A
B = 2000
for reg in A.REGIMES:
    for s in (2024, 7331):
        for ci, c in enumerate(json.load(open(f"task2a_seed{s}_{reg}.json"))):
            k = np.asarray(c["k_eval"], float); n = len(k)
            net = np.asarray(c["net_of_window"])
            M = {m: A.mats(c["methods"][m]) for m in ("MT-FW+SC",) + A.SOTA}
            rng = np.random.default_rng([s, ci, 77])
            C = np.zeros((B, n))
            for b in range(B):
                for j in rng.integers(0, 5, 5):
                    C[b, net == j] += 1
            kk = C @ k; nn = C.sum(1)
            row = []
            for cap in A.CAPS:
                def vc(HF):
                    H, F = HF
                    Hs, Fs = C @ H, C @ F
                    ok = Fs / nn[:, None] <= float(cap) + 1e-12
                    r = np.where(ok, Hs, -1.0).max(axis=1)
                    return np.where(r < 0, 0.0, r) / np.maximum(1.0, kk)
                d = vc(M["MT-FW+SC"]) - np.max(np.stack([vc(M[m]) for m in A.SOTA]), axis=0)
                p = A.point(M["MT-FW+SC"], k, float(cap)) - max(A.point(M[m], k, float(cap)) for m in A.SOTA)
                lo, hi = np.percentile(d, [2.5, 97.5])
                row.append(f"{p:+.2f} [{lo:+.2f}, {hi:+.2f}] P(d<0.05)={np.mean(d < 0.05):.2f}")
            print(f"{reg:6s} {s} {c['scenario']:15s} " + " | ".join(row))
