"""Own hierarchical bootstrap (networks, then windows) for two cells, TS-SFW+SC - BOMP-PK+SC, independent RNG."""
import json, os
import numpy as np
HERE = "/home/claude/work"

def HFmat(recs):
    thr = np.array(sorted({s for r in recs for s, h in r}, reverse=True))
    H = np.zeros((len(recs), len(thr))); F = np.zeros((len(recs), len(thr)))
    for w, r in enumerate(recs):
        for j, t in enumerate(thr):
            hs = {h for s, h in r if s >= t and h >= 0}
            H[w, j] = len(hs); F[w, j] = sum(1 for s, h in r if s >= t and h < 0)
    return H, F

def rec_cap(H, F, C, k, cap):
    Hs, Fs = C @ H, C @ F
    n = C.sum(1)[:, None]
    r = np.where(Fs / n <= cap + 1e-12, Hs, 0).max(1)
    return r / np.maximum(1, C @ k)

for seed, reg, cell in ((2024, "thin", "ladder_4"), (2024, "strong", "ladder_6")):
    recs = {"TS": [], "BP": []}; k = []; net = []
    for j in range(30):
        p = json.load(open(os.path.join(HERE, "task2e_parts", f"ev_seed{seed}_{reg}_{cell}_net{j}.json")))
        for w in p["windows"]:
            k.append(w["k"]); net.append(j)
            recs["TS"].append(w["records"]["TS-SFW+SC"]); recs["BP"].append(w["records"]["BOMP-PK+SC"])
    k = np.array(k, float); net = np.array(net)
    HF = {m: HFmat(recs[m]) for m in recs}
    rng = np.random.default_rng(424242)
    B = 2000
    C = np.zeros((B, len(k)))
    pos = [np.where(net == j)[0] for j in range(30)]
    for b in range(B):
        for jj in rng.integers(0, 30, 30):
            p_ = pos[jj]
            np.add.at(C[b], p_[rng.integers(0, len(p_), len(p_))], 1.0)
    one = np.ones((1, len(k)))
    for cap in (0.1, 0.25):
        pt = rec_cap(*HF["TS"], one, k, cap)[0] - rec_cap(*HF["BP"], one, k, cap)[0]
        d = rec_cap(*HF["TS"], C, k, cap) - rec_cap(*HF["BP"], C, k, cap)
        print(f"{reg} {cell} {seed} cap {cap}: point {pt:+.3f}  95% CI [{np.percentile(d, 2.5):+.3f}, {np.percentile(d, 97.5):+.3f}]  "
              f"p1 = {(1 + np.sum(d <= 0)) / (B + 1):.4f}")
