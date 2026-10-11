"""stats_verif.py - independent re-derivation of 2F point estimates and a hierarchical bootstrap (verification).

Own FROC code (per-window loop, no reuse of analyze_2f / froc_exact).  Registered hit rule: an atom with hit >= 0
counts as a hit of that true event the first time it is admitted; later atoms on the same event count neither as hit
nor as false (the 2A-2E rule).  Sensitivity: duplicates counted as false ("dupfalse").
Also threshold-free, within-window diagnostics: unique hits among the first k atoms (k = true count) and among all
atoms; per-window top-n operating point (rank threshold, no cross-window score comparison).
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
PARTS = os.path.join(WORK, "task2f_parts")
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
METHODS = ("MT-FW+SC", "BOMP-PK+SC", "BOMP-PK-CAP+SC", "MT-FW-NOCAP+SC")


def load(seed, arm, reg, nets):
    wins = []
    for cell in CELLS:
        for j in nets:
            e = json.load(open(os.path.join(PARTS, f"{arm}_seed{seed}_{reg}_{cell}_net{j}.json")))
            assert e["seed"] == seed and e["regime"] == reg and e["cell"] == cell and e["net"] == j
            for w in e["windows"]:
                wins.append({"cell": cell, "net": j, "k": w["k"], "rec": w["records"]})
    return wins


def froc_point(wins, m, caps=(0.1, 0.25), dupfalse=False, weights=None):
    """exact FROC: every distinct score is a threshold; recall = max over admissible thresholds."""
    wts = np.ones(len(wins)) if weights is None else np.asarray(weights, float)
    K = float(sum(wt * w["k"] for wt, w in zip(wts, wins)))
    nw = float(wts.sum())
    # per window: list of (score, delta_hit, delta_false) in record order
    thr = sorted({s for w in wins for s, _h in w["rec"][m]}, reverse=True)
    idx = {t: i for i, t in enumerate(thr)}
    dH = np.zeros(len(thr))
    dF = np.zeros(len(thr))
    for wt, w in zip(wts, wins):
        if wt == 0:
            continue
        seen = set()
        for s, h in w["rec"][m]:
            if h >= 0 and h not in seen:
                seen.add(h)
                dH[idx[s]] += wt
            elif h < 0 or dupfalse:
                dF[idx[s]] += wt
    H, F = np.cumsum(dH), np.cumsum(dF)
    out = {}
    for c in caps:
        ok = F / nw <= c + 1e-12
        out[c] = float(H[ok].max() / max(K, 1.0)) if ok.any() else 0.0
    return out


def within_window(wins, m):
    """mean unique hits among first k atoms / k, among all atoms / k (pooled ratio), and top-n recall at <=0.25 false."""
    K = sum(w["k"] for w in wins)
    hk = ha = 0
    for w in wins:
        r = w["rec"][m]
        hk += len({h for _s, h in r[:w["k"]] if h >= 0})
        ha += len({h for _s, h in r if h >= 0})
    topn = {}
    for n in range(1, 8):
        H = F = 0
        for w in wins:
            r = w["rec"][m][:n]
            H += len({h for _s, h in r if h >= 0})
            F += sum(1 for _s, h in r if h < 0)
        topn[n] = (H / K, F / len(wins))
    return hk / K, ha / K, topn


def hier_boot(wins, ma, mb, B=2000, seed=20260928, cap=0.25):
    """cells as strata; networks with replacement within cell; windows with replacement within network."""
    rng = np.random.default_rng(seed)
    cells = {}
    for i, w in enumerate(wins):
        cells.setdefault(w["cell"], {}).setdefault(w["net"], []).append(i)
    d = np.zeros(B)
    for b in range(B):
        wt = np.zeros(len(wins))
        for c, nets in cells.items():
            keys = list(nets)
            for kk in rng.integers(0, len(keys), len(keys)):
                ii = nets[keys[kk]]
                for q in rng.integers(0, len(ii), len(ii)):
                    wt[ii[q]] += 1
        d[b] = froc_point(wins, ma, (cap,), weights=wt)[cap] - froc_point(wins, mb, (cap,), weights=wt)[cap]
    return d


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "points":
        for arm, nets in (("dense", range(96)), ("std", range(96, 112))):
            for seed in (4242, 8888):
                for reg in ("thin", "strong"):
                    W = load(seed, arm, reg, nets)
                    pt = {m: froc_point(W, m) for m in METHODS}
                    pd = {m: froc_point(W, m, dupfalse=True) for m in METHODS}
                    print(f"{arm} {seed} {reg} nwin={len(W)} K={sum(w['k'] for w in W)}")
                    for m in METHODS:
                        print(f"   {m:16s} registered {pt[m][0.1]:.4f}/{pt[m][0.25]:.4f}   dup-as-false "
                              f"{pd[m][0.1]:.4f}/{pd[m][0.25]:.4f}")
                    for cap in (0.1, 0.25):
                        print(f"   diff MT-FW - BOMP-PK cap {cap}: registered {pt['MT-FW+SC'][cap] - pt['BOMP-PK+SC'][cap]:+.4f}"
                              f"  dup-as-false {pd['MT-FW+SC'][cap] - pd['BOMP-PK+SC'][cap]:+.4f}"
                              f"  MT-FW - NOCAP {pt['MT-FW+SC'][cap] - pt['MT-FW-NOCAP+SC'][cap]:+.4f}")
                    if arm == "dense":
                        for m in METHODS:
                            fk, fa, topn = within_window(W, m)
                            ndup = sum(len([h for _s, h in w["rec"][m] if h >= 0]) -
                                       len({h for _s, h in w["rec"][m] if h >= 0}) for w in W)
                            print(f"   {m:16s} hits@first-k {fk:.4f}  hits@all-7 {fa:.4f}  duplicate atoms {ndup}  "
                                  "top-n (recall, false/win): " +
                                  " ".join(f"{n}:({a:.3f},{b:.2f})" for n, (a, b) in topn.items()))
    elif mode == "boot":
        arm, seed, reg = sys.argv[2], int(sys.argv[3]), sys.argv[4]
        nets = range(96) if arm == "dense" else range(96, 112)
        W = load(seed, arm, reg, nets) if reg != "both" else load(seed, arm, "thin", nets) + load(seed, arm, "strong", nets)
        B = int(sys.argv[5]) if len(sys.argv) > 5 else 2000
        pt = froc_point(W, "MT-FW+SC")[0.25] - froc_point(W, "BOMP-PK+SC")[0.25]
        d = hier_boot(W, "MT-FW+SC", "BOMP-PK+SC", B=B)
        print(f"{arm} {seed} {reg}: point {pt:+.4f}  95% CI [{np.percentile(d, 2.5):+.4f}, {np.percentile(d, 97.5):+.4f}]"
              f"  90% CI [{np.percentile(d, 5):+.4f}, {np.percentile(d, 95):+.4f}]  sd {d.std():.4f}"
              f"  p(d<=0) {(1 + np.sum(d <= 0)) / (B + 1):.4f}  p(d<=-0.05) {(1 + np.sum(d <= -0.05)) / (B + 1):.4f}"
              f"  p(d>=0.05) {(1 + np.sum(d >= 0.05)) / (B + 1):.4f}")
