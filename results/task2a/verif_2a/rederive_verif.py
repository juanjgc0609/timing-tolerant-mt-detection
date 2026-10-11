"""Independent re-derivation of recall at <= cap realized false events per window from the 2A JSON."""
import json, sys
import numpy as np

def rec_recall(records, k, cap):
    n = len(records); K = int(sum(k))
    ent = [(s, w, h) for w, r in enumerate(records) for s, h in r]
    thr = sorted({s for s, _w, _h in ent}, reverse=True)
    best = 0.0
    for t in thr:
        hits = {(w, h) for s, w, h in ent if s >= t and h >= 0}
        nf = sum(1 for s, w, h in ent if s >= t and h < 0)
        if nf / n <= cap + 1e-12:
            best = max(best, len(hits) / K)
    return best

def per_recall(per_w, k, cap):
    n = len(per_w); K = int(sum(k))
    best = 0.0
    for c in range(len(per_w[0])):
        H = sum(w[c][0] for w in per_w); F = sum(w[c][1] for w in per_w)
        if F / n <= cap + 1e-12:
            best = max(best, H / K)
    return best

for seed, reg, cell in [(2024, "thin", "ladder_6"), (7331, "strong", "unilateral_c40"), (7331, "thin", "ladder_4"), (2024, "strong", "ladder_4")]:
    C = {c["scenario"]: c for c in json.load(open(f"task2a_seed{seed}_{reg}.json"))}[cell]
    k = C["k_eval"]
    assert all(len(w["truth"]) == kk for w, kk in zip(C["windows_meta"], k))
    out = []
    for m in ("MT-FW+SC", "QM-3C+SC", "MT-FW"):
        out.append(f"{m} " + "/".join(f"{rec_recall(C['methods'][m]['records'], k, cap):.4f}" for cap in (0.1, 0.25)))
    for m in ("GaMMA-DY+SC", "GaMMA-STALTA+SC"):
        out.append(f"{m} " + "/".join(f"{per_recall(C['methods'][m]['per_window'], k, cap):.4f}" for cap in (0.1, 0.25)))
    print(seed, reg, cell, "n_win", len(k), "events", sum(k), "|", " | ".join(out))
