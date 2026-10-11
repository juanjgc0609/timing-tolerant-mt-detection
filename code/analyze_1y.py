"""analyze_1y.py - reglas Y-a, Y-b y Y-c de la Tarea 1y (clutter consistente), por regimen, y tablas.
IC bootstrap pareado con el punto de operacion y la configuracion de GaMMA re-elegidos en cada replica."""
import json

import numpy as np

B = 2000


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
    ds = [val(a, i, k, cap) - val(b, i, k, cap) for i in (rng.integers(0, n, n) for _ in range(B))]
    return d0, np.percentile(ds, 2.5), np.percentile(ds, 97.5)

CAPS = ("0.1", "0.25")
LAB = {"unilateral_c40": "one-sided, 8 st.", "ladder_4": "4 stations", "ladder_6": "6 stations"}
SEEDS = (2024, 7331)
SHOW = ("WAVE-PT@20", "WAVE-MT@20", "MT-FW", "MT-FW-NODE", "GAMMA")
PAIRS = (("MT-FW", "GAMMA"), ("MT-FW", "MT-FW-NODE"), ("MT-FW", "WAVE-MT@20"), ("WAVE-MT@20", "GAMMA"),
         ("WAVE-PT@20", "GAMMA"))
out, tables = [], []
for reg in ("thin", "strong"):
    R = {s: {r["scenario"]: r for r in json.load(open(f"task1y_seed{s}_{reg}_r100.json"))} for s in SEEDS}
    rules = {"Y-a": {}, "Y-b": {}, "Y-c": {}}
    t1 = [f"TABLE recall ({reg})", "| Cell | picks clutter/true | " + " | ".join(SHOW) + " |",
          "|:--|:--|" + ":--|" * len(SHOW)]
    t2 = [f"TABLE differences at <=0.25 ({reg})", "| Cell | " + " | ".join(f"{a} - {b}" for a, b in PAIRS) + " |",
          "|:--|" + ":--|" * len(PAIRS)]
    for s in SEEDS:
        na = nb = nc = 0
        for sc in LAB:
            r = R[s][sc]
            k = np.asarray(r["k_eval"], float)
            ra = {m: r["methods"][m]["recall_at"] for m in r["methods"]}
            g = lambda a, b, caps=CAPS: all(ra[a][c] - ra[b][c] >= 0.05 for c in caps)
            na += g("MT-FW", "GAMMA")
            nb += g("MT-FW", "MT-FW-NODE")
            nc += ra["MT-FW"]["0.25"] - max(ra["WAVE-MT@20"]["0.25"], ra["WAVE-PT@20"]["0.25"]) >= 0.05
            M = {m: mats(r["methods"][m]) for m in r["methods"]}
            dd = [boot(M[a], M[b], k, s) for a, b in PAIRS]
            out.append(f"{reg} {s} {LAB[sc]:17s} eventos {int(k.sum())} picks {r['clutter_picks_per_window']:.1f}/"
                       f"{r['true_picks_per_window']:.1f} :: " + " | ".join(
                           f"{m} " + "/".join(f"{ra[m][c]:.2f}" for c in CAPS) for m in SHOW))
            t1.append(f"| {LAB[sc]} ({s}) | {r['clutter_picks_per_window']:.1f}/{r['true_picks_per_window']:.1f} | "
                      + " | ".join("/".join(f"{ra[m][c]:.2f}" for c in CAPS) for m in SHOW) + " |")
            t2.append(f"| {LAB[sc]} ({s}) | " + " | ".join(f"{d:+.2f} [{lo:+.2f}, {hi:+.2f}]" for d, lo, hi in dd) + " |")
        rules["Y-a"][s], rules["Y-b"][s], rules["Y-c"][s] = na >= 2, nb >= 2, nc >= 2
        out.append(f"   {reg} seed {s}: Y-a {na}/3; Y-b {nb}/3; Y-c {nc}/3")
    out.append(f"{reg}: " + " ".join(f"{k_}: {all(v.values())}" for k_, v in rules.items()))
    tables += t1 + [""] + t2 + [""]
txt = "\n".join(out + [""] + tables)
print(txt)
open("analyze_1y_out.txt", "w").write(txt + "\n")
