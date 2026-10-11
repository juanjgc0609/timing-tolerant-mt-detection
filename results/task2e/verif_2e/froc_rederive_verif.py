"""Independent re-derivation of Table 1 / 2E-a point estimates from task2e_parts (own FROC code, no project imports)."""
import json, os, glob
import numpy as np
HERE = "/home/claude/work"
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
M = ("TS-SFW+SC", "BOMP-PK-SLIDE+SC", "BOMP-PK+SC", "MT-FW+SC")

def recall_at(recs, k, caps=(0.1, 0.25)):
    """exact FROC: every distinct score is a threshold; at threshold t each window emits its atoms with score >= t;
    hits = distinct true sources matched, false = atoms with hit == -1. Recall at cap = max hits / events over
    thresholds with false / windows <= cap (0 if none)."""
    allsc = sorted({s for r in recs for s, h in r}, reverse=True)
    n = len(recs); K = sum(k)
    best = {c: 0 for c in caps}
    for t in allsc:
        H = F = 0
        for r in recs:
            hs = set()
            for s, h in r:
                if s >= t:
                    if h >= 0: hs.add(h)
                    else: F += 1
            H += len(hs)
        for c in caps:
            if F / n <= c + 1e-12 and H > best[c]:
                best[c] = H
    return {c: best[c] / K for c in caps}

out = []
for seed in (2024, 7331):
    for reg in ("thin", "strong"):
        for cell in CELLS:
            recs = {m: [] for m in M}; k = []
            for j in range(30):
                p = json.load(open(os.path.join(HERE, "task2e_parts", f"ev_seed{seed}_{reg}_{cell}_net{j}.json")))
                for w in p["windows"]:
                    k.append(w["k"])
                    for m in M:
                        recs[m].append([tuple(x) for x in w["records"][m]])
            r = {m: recall_at(recs[m], k) for m in M}
            d = {c: r["TS-SFW+SC"][c] - r["BOMP-PK+SC"][c] for c in (0.1, 0.25)}
            db = {c: r["TS-SFW+SC"][c] - r["BOMP-PK-SLIDE+SC"][c] for c in (0.1, 0.25)}
            out.append((reg, seed, cell, len(k), sum(k), r, d, db))
            print(f"{reg:6s} {seed} {cell:15s} win {len(k)} ev {sum(k)} | " +
                  " | ".join(f"{m} {r[m][0.1]:.3f}/{r[m][0.25]:.3f}" for m in M) +
                  f" || 2E-a d {d[0.1]:+.3f}/{d[0.25]:+.3f}  2E-b d {db[0.1]:+.3f}/{db[0.25]:+.3f}", flush=True)
json.dump([[a, b, c, n, K, {m: {str(x): y for x, y in v.items()} for m, v in r.items()}] for a, b, c, n, K, r, d, db in out],
          open(os.path.join(HERE, "verif_2e", "froc_rederive_verif.json"), "w"))
