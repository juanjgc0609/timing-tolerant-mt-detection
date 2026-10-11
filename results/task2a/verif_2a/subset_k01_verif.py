"""Exploratory: recall at <= 0.1 / <= 0.25 false per window restricted to windows with 0 or 1 source (no
simultaneous events), and the fraction of 2-source windows whose origin times are within 0.75 s."""
import json
import numpy as np
from froc_exact import exact_HF
from task1r_fair import recall_at

def rec_at(recs, k, cap):
    _t, H, F = exact_HF([[tuple(x) for x in r] for r in recs])
    return recall_at(H, F, int(sum(k)), cap)

def per_at(pw, k, cap):
    n = len(pw); K = max(1, sum(k)); best = 0.0
    for c in range(len(pw[0])):
        if sum(w[c][1] for w in pw) / n <= cap + 1e-12:
            best = max(best, sum(w[c][0] for w in pw) / K)
    return best

M = ("MT-FW+SC", "MT-FW", "QM-3C+SC", "QM-3C")
G = ("GaMMA-DY+SC", "GaMMA-STALTA+SC")
close = tot2 = 0
for reg in ("thin", "strong"):
    for s in (2024, 7331):
        for c in json.load(open(f"task2a_seed{s}_{reg}.json")):
            k = c["k_eval"]
            idx = [i for i, kk in enumerate(k) if kk <= 1]
            for i, kk in enumerate(k):
                if kk == 2:
                    tot2 += 1
                    close += abs(c["windows_meta"][i]["taus"][0] - c["windows_meta"][i]["taus"][1]) < 0.75
            kk = [k[i] for i in idx]
            row = []
            for m in M:
                r = [c["methods"][m]["records"][i] for i in idx]
                row.append(f"{m} {rec_at(r, kk, 0.1):.2f}/{rec_at(r, kk, 0.25):.2f}")
            for m in G:
                pw = [c["methods"][m]["per_window"][i] for i in idx]
                row.append(f"{m} {per_at(pw, kk, 0.1):.2f}/{per_at(pw, kk, 0.25):.2f}")
            print(f"{reg:6s} {s} {c['scenario']:15s} win {len(idx)} ev {sum(kk)} | " + " | ".join(row))
print("2-source windows with |tau1 - tau2| < 0.75 s:", close, "/", tot2)
