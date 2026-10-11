"""analyze_1u.py - reglas U-a y U-b de la Tarea 1u (medio no modelado). Usa la parte de forma de onda
con curva FROC exacta (task1u_seed*_r100_wave_exact.json) si existe, y GaMMA de task1u_seed*_r100.json
(mismas ventanas)."""
import json
import os

import numpy as np

CAPS = ("0.1", "0.25")
LAB = {"unilateral_c40": "one-sided, 8 st.", "ladder_4": "4 stations", "ladder_6": "6 stations"}
B = 2000


def op_wave(m, cap):
    H, F = np.asarray(m["H"]), np.asarray(m["F"])
    fr = F.sum(0) / H.shape[0]
    ok = np.where(fr <= cap)[0]
    if len(ok) == 0:
        return np.zeros(H.shape[0])
    return H[:, ok[np.argmax(H.sum(0)[ok])]].astype(float)


def op_gamma(m, cap, n):
    ok = [(rc, ci) for ci, (fr, rc) in enumerate(m["points"]) if fr <= cap]
    return np.asarray(m["per"][str(max(ok)[1])][0], float) if ok else np.zeros(n)


def boot(ha, hb, k, rng):
    n = len(k)
    ds = [(ha[i].sum() - hb[i].sum()) / max(1, k[i].sum()) for i in (rng.integers(0, n, n) for _ in range(B))]
    return (ha.sum() - hb.sum()) / max(1, k.sum()), np.percentile(ds, 2.5), np.percentile(ds, 97.5)


out, ua, ub = [], {}, {}
for seed in (2024, 7331):
    G = {r["scenario"]: r for r in json.load(open(f"task1u_seed{seed}_r100.json"))}
    fx = f"task1u_seed{seed}_r100_wave_exact.json"
    W = {r["scenario"]: r for r in json.load(open(fx))} if os.path.exists(fx) else G
    out.append(f"[{fx if os.path.exists(fx) else 'barrido registrado'}]")
    rng = np.random.default_rng(seed)
    na = nb = 0
    for sc in LAB:
        w, g = W[sc], G[sc]
        k = np.asarray(w["k_eval"], float)
        assert w["k_eval"] == g["k_eval"]
        ra = {m: w["methods"][m]["recall_at"] for m in ("WAVE-E", "WAVE-P", "WAVE-PSH")}
        ra["GAMMA"] = g["methods"]["GAMMA"]["recall_at"]
        out.append(f"{seed} {LAB[sc]:17s} eventos {int(k.sum()):4d} :: " + " | ".join(
            f"{m} " + "/".join(f"{ra[m][c]:.2f}" for c in CAPS) for m in ra))
        if "recall_at_coarse" in w["methods"]["WAVE-P"]:
            out.append("     barrido registrado: " + " | ".join(
                f"{m} " + "/".join(f"{w['methods'][m]['recall_at_coarse'][c]:.2f}" for c in CAPS)
                for m in ("WAVE-E", "WAVE-P", "WAVE-PSH")))
        na += all(ra["WAVE-P"][c] - ra["WAVE-E"][c] >= 0.05 for c in CAPS)
        nb += all(ra["WAVE-P"][c] - ra["GAMMA"][c] >= 0.05 for c in CAPS)
        for a, hb in (("WAVE-E", op_wave(w["methods"]["WAVE-E"], 0.25)),
                      ("WAVE-PSH", op_wave(w["methods"]["WAVE-PSH"], 0.25)),
                      ("GAMMA", op_gamma(g["methods"]["GAMMA"], 0.25, len(k)))):
            d, lo, hi = boot(op_wave(w["methods"]["WAVE-P"], 0.25), hb, k, rng)
            out.append(f"     WAVE-P - {a} @0.25: {d:+.3f} [{lo:+.3f}, {hi:+.3f}]")
    ua[seed], ub[seed] = na >= 2, nb >= 2
    out.append(f"   seed {seed}: U-a celdas {na}/3; U-b celdas {nb}/3")
out.append(f"U-a: {all(ua.values())}  U-b: {all(ub.values())}")
txt = "\n".join(out)
print(txt)
open("analyze_1u_out.txt", "w").write(txt + "\n")
