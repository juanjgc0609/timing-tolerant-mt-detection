"""analyze_1w.py - reglas W-a, W-b y W-c de la Tarea 1w y tabla Y.11. IC bootstrap pareado con el punto
de operacion re-elegido en cada replica; GaMMA con su configuracion re-elegida en cada replica."""
import json
import sys

import numpy as np

TAG = sys.argv[1] if len(sys.argv) > 1 else ""      # "_ev9000_m131_wave" para la replica

CAPS = ("0.1", "0.25")
LAB = {"unilateral_c40": "one-sided, 8 st.", "ladder_4": "4 stations", "ladder_6": "6 stations"}
B = 2000
SEEDS = (2024, 7331)


def mats(m):
    if "H" in m:
        return ("grid", np.asarray(m["H"], float), np.asarray(m["F"], float))
    per = m["per"]
    H = np.array([per[str(i)][0] for i in range(len(per))], float).T
    F = np.array([per[str(i)][1] for i in range(len(per))], float).T
    return ("cfg", H, F)


def val(x, i, k, cap):
    _t, H, F = x
    fr = F[i].sum(0) / len(i)
    ok = fr <= cap
    return float(H[i].sum(0)[ok].max()) / max(1, k[i].sum()) if ok.any() else 0.0


def boot(a, b, k, seed, cap=0.25):
    rng = np.random.default_rng(seed)
    n = len(k)
    allw = np.arange(n)
    d0 = val(a, allw, k, cap) - val(b, allw, k, cap)
    ds = []
    for _ in range(B):
        i = rng.integers(0, n, n)
        ds.append(val(a, i, k, cap) - val(b, i, k, cap))
    return d0, np.percentile(ds, 2.5), np.percentile(ds, 97.5)


R = {s: {r["scenario"]: r for r in json.load(open(f"task1w_seed{s}_r120{TAG}.json"))} for s in SEEDS}
WAVE_ONLY = "wave" in TAG
out = []
ok = {"W-a": {}, "W-b": {}, "W-c PHYS": {}, "W-c PT": {}}
tab = ["| Cell | detected | KIN-LR | PHYS-LR | PHYS-NODE | WAVE-E | WAVE-P | WAVE-PT | WAVE-PT-NODE | GaMMA |",
       "|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|"]
tab2 = ["| Cell | PHYS-LR - PHYS-NODE | PHYS-LR - KIN-LR | WAVE-PT - WAVE-PT-NODE | WAVE-PT - WAVE-E | PHYS-LR - GaMMA | WAVE-PT - GaMMA |",
        "|:--|:--|:--|:--|:--|:--|:--|"]
for s in SEEDS:
    na = nb = ncp = nct = 0
    excl = False
    for sc in LAB:
        r = R[s][sc]
        k = np.asarray(r["k_eval"], float)
        ra = {m: r["methods"][m]["recall_at"] for m in r["methods"]}
        out.append(f"{s} {LAB[sc]:17s} eventos {int(k.sum())} detectados {r['detected_fraction']:.2f} :: "
                   + " | ".join(f"{m} " + "/".join(f"{ra[m][c]:.2f}" for c in CAPS) for m in ra))
        g = lambda a, b: all(ra[a][c] - ra[b][c] >= 0.05 for c in CAPS) if a in ra and b in ra else False
        na += g("PHYS-LR", "PHYS-NODE")
        nb += g("WAVE-PT", "WAVE-PT-NODE") and g("WAVE-PT", "WAVE-E")
        ncp += g("PHYS-LR", "GAMMA")
        nct += g("WAVE-PT", "GAMMA")
        M = {m: mats(r["methods"][m]) for m in r["methods"]}
        diffs = {}
        pairs = (("PHYS-LR", "PHYS-NODE"), ("PHYS-LR", "KIN-LR"), ("WAVE-PT", "WAVE-PT-NODE"),
                 ("WAVE-PT", "WAVE-E"), ("PHYS-LR", "GAMMA"), ("WAVE-PT", "GAMMA"))
        for a, b in [p for p in pairs if p[0] in M and p[1] in M]:
            diffs[(a, b)] = boot(M[a], M[b], k, s)
            d, lo, hi = diffs[(a, b)]
            out.append(f"     {a} - {b} @0.25: {d:+.3f} [{lo:+.3f}, {hi:+.3f}]")
        if ("PHYS-LR", "PHYS-NODE") in diffs and diffs[("PHYS-LR", "PHYS-NODE")][1] > 0:
            excl = True
        tab.append(f"| {LAB[sc]} ({s}) | {r['detected_fraction']:.2f} | " + " | ".join(
            ("/".join(f"{ra[m][c]:.2f}" for c in CAPS) if m in ra else "-") for m in ("KIN-LR", "PHYS-LR", "PHYS-NODE", "WAVE-E",
                                                               "WAVE-P", "WAVE-PT", "WAVE-PT-NODE", "GAMMA")) + " |")
        tab2.append(f"| {LAB[sc]} ({s}) | " + " | ".join(
            f"{d:+.2f} [{lo:+.2f}, {hi:+.2f}]" for d, lo, hi in diffs.values()) + " |")
    ok["W-a"][s] = na >= 2 and excl
    ok["W-b"][s] = nb >= 2
    ok["W-c PHYS"][s] = ncp >= 2
    ok["W-c PT"][s] = nct >= 2
    out.append(f"   seed {s}: W-a celdas {na}/3 (IC excluye 0: {excl}); W-b celdas {nb}/3; "
               f"W-c PHYS-LR {ncp}/3; W-c WAVE-PT {nct}/3")
out.append(" ".join(f"{k_}: {all(v.values())}" for k_, v in ok.items()))
txt = "\n".join(out + [""] + tab + [""] + tab2)
print(txt)
open(f"analyze_1w{TAG}_out.txt", "w").write(txt + "\n")
