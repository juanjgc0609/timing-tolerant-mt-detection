"""analyze_1x.py - reglas X-a, X-b y X-c de la Tarea 1x, por condicion, y tablas. IC bootstrap pareado
con el punto de operacion (y la configuracion de GaMMA) re-elegido en cada replica."""
import json

import numpy as np

CAPS = ("0.1", "0.25")
LAB = {"unilateral_c40": "one-sided, 8 st.", "ladder_4": "4 stations", "ladder_6": "6 stations"}
B = 2000
SEEDS = (2024, 7331)
SHOW = ("WAVE-E", "WAVE-PT", "WAVE-MT", "WAVE-MT-NODE", "MT-FW", "GAMMA")


def mats(m):
    if "H" in m:
        return (np.asarray(m["H"], float), np.asarray(m["F"], float))
    per = m["per"]
    H = np.array([per[str(i)][0] for i in range(len(per))], float).T
    F = np.array([per[str(i)][1] for i in range(len(per))], float).T
    return (H, F)


def val(x, i, k, cap):
    H, F = x
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


out, tables = [], []
for cond in ("ideal", "pert"):
    R = {s: {r["scenario"]: r for r in json.load(open(f"task1x_seed{s}_{cond}_r100.json"))} for s in SEEDS}
    rules = {"X-a": {}, "X-b": {}, "X-b FW": {}, "X-c": {}}
    t1 = [f"TABLE recall ({cond})", "| Cell | detected | " + " | ".join(SHOW) + " |", "|:--|:--|" + ":--|" * len(SHOW)]
    t2 = [f"TABLE differences at <=0.25 ({cond})",
          "| Cell | WAVE-MT - WAVE-MT-NODE | WAVE-MT - WAVE-PT | WAVE-MT - GaMMA | MT-FW - WAVE-MT | MT-FW - GaMMA |",
          "|:--|:--|:--|:--|:--|:--|"]
    for s in SEEDS:
        na = nb = nbf = nc = 0
        for sc in LAB:
            r = R[s][sc]
            k = np.asarray(r["k_eval"], float)
            ra = {m: r["methods"][m]["recall_at"] for m in r["methods"]}
            g = lambda a, b: all(ra[a][c] - ra[b][c] >= 0.05 for c in CAPS)
            na += g("WAVE-MT", "WAVE-MT-NODE") and g("WAVE-MT", "WAVE-PT")
            nb += g("WAVE-MT", "GAMMA")
            nbf += g("MT-FW", "GAMMA")
            nc += ra["MT-FW"]["0.25"] - ra["WAVE-MT"]["0.25"] >= 0.05
            M = {m: mats(r["methods"][m]) for m in r["methods"]}
            dd = [boot(M[a], M[b], k, s) for a, b in (("WAVE-MT", "WAVE-MT-NODE"), ("WAVE-MT", "WAVE-PT"),
                                                       ("WAVE-MT", "GAMMA"), ("MT-FW", "WAVE-MT"), ("MT-FW", "GAMMA"))]
            out.append(f"{cond} {s} {LAB[sc]:17s} eventos {int(k.sum())} det {r['detected_fraction']:.2f} :: "
                       + " | ".join(f"{m} " + "/".join(f"{ra[m][c]:.2f}" for c in CAPS) for m in SHOW))
            t1.append(f"| {LAB[sc]} ({s}) | {r['detected_fraction']:.2f} | "
                      + " | ".join("/".join(f"{ra[m][c]:.2f}" for c in CAPS) for m in SHOW) + " |")
            t2.append(f"| {LAB[sc]} ({s}) | " + " | ".join(f"{d:+.2f} [{lo:+.2f}, {hi:+.2f}]" for d, lo, hi in dd) + " |")
        rules["X-a"][s], rules["X-b"][s], rules["X-b FW"][s], rules["X-c"][s] = na >= 2, nb >= 2, nbf >= 2, nc >= 1
        out.append(f"   {cond} seed {s}: X-a {na}/3; X-b (WAVE-MT) {nb}/3; X-b (MT-FW) {nbf}/3; X-c {nc}/3")
    out.append(f"{cond}: " + " ".join(f"{k_}: {all(v.values())}" for k_, v in rules.items()))
    tables += t1 + [""] + t2 + [""]
txt = "\n".join(out + [""] + tables)
print(txt)
open("analyze_1x_out.txt", "w").write(txt + "\n")
