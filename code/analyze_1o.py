"""Tablas del Apendice Y.5 a partir de task1o_seed*.json y task1o_baselines_seed*.json."""
import json
import os

import numpy as np
from scipy.stats import binomtest

M = ("TOL-TB", "KIN-LR", "PHYS-LR", "WAVE-E", "WAVE-P", "GAMMA", "PYOCTO")
LAB = {"unilateral_c40": "one-sided, 8 st.", "ladder_4": "4 stations", "ladder_6": "6 stations"}
CL = {"consistente": "consistent", "jit0.5": "control, 0.5 s"}


def load(seed):
    own = {(r["scenario"], r["cond"]): r for r in json.load(open(f"task1o_seed{seed}.json"))}
    fb = f"task1o_baselines_seed{seed}.json"
    if os.path.exists(fb):
        for r in json.load(open(fb)):
            key = (r["scenario"], r["cond"])
            if key in own:
                own[key]["methods"].update(r["methods"])
    return own


def pts(r, m):
    nt, n = sum(r["k_eval"]), len(r["k_eval"])
    d = r["methods"].get(m)
    if d is None:
        return None
    return [(sum(d[f"b{b}"]["false"]) / n, sum(d[f"b{b}"]["hits"]) / nt) for b in d["budgets"]]


def cap(p, x):
    return max([rc for fr, rc in p if fr <= x] + [0.0])


def mcn(a, b):
    a, b = np.asarray(a), np.asarray(b)
    n10, n01 = int(np.sum((a == 1) & (b == 0))), int(np.sum((a == 0) & (b == 1)))
    return n10, n01, (1.0 if n10 + n01 == 0 else binomtest(n10, n10 + n01, 0.5).pvalue)


D = {s: load(s) for s in (2024, 7331)}
print("TABLA A: recall por evento con tasa de falsos realizada <= 0.1 / <= 0.25 por ventana, semilla 2024 ; 7331")
print("| Generator | Cell | " + " | ".join(M) + " |")
print("|:--|:--|" + ":--|" * len(M))
for c in CL:
    for s_ in LAB:
        cells = []
        for m in M:
            vals = []
            for sd in (2024, 7331):
                r = D[sd].get((s_, c))
                p = pts(r, m) if r else None
                vals.append("-" if p is None else f"{cap(p, 0.1):.2f}/{cap(p, 0.25):.2f}")
            cells.append(" ; ".join(vals))
        print(f"| {CL[c]} | {LAB[s_]} | " + " | ".join(cells) + " |")
print("\nTABLA B: pares por evento al presupuesto 0.1 (detectado solo por A : solo por B, p)")
for c in CL:
    for s_ in LAB:
        row = []
        for a, b in (("KIN-LR", "PHYS-LR"), ("GAMMA", "PHYS-LR"), ("PYOCTO", "PHYS-LR"), ("PHYS-LR", "WAVE-P")):
            vv = []
            for sd in (2024, 7331):
                r = D[sd].get((s_, c))
                if r is None or a not in r["methods"] or b not in r["methods"]:
                    vv.append("-")
                    continue
                n10, n01, p = mcn(r["methods"][a]["b0.1"]["det"], r["methods"][b]["b0.1"]["det"])
                vv.append(f"{n10}:{n01}, p={p:.2g}")
            row.append(f"{a} vs {b}: " + " ; ".join(vv))
        print(f"{CL[c]} {LAB[s_]} :: " + " || ".join(row))
print("\nrealizado al presupuesto nominal 0.1 (recall, falsos/ventana):")
for c in CL:
    for s_ in LAB:
        out = []
        for m in M:
            vv = []
            for sd in (2024, 7331):
                r = D[sd].get((s_, c))
                p = pts(r, m) if r else None
                vv.append("-" if p is None else f"{p[1][1]:.2f},{p[1][0]:.2f}")
            out.append(f"{m} " + " ; ".join(vv))
        print(f"{CL[c]} {LAB[s_]} :: " + " | ".join(out))
