"""analyze_1r200.py - reglas R-a, R-b y R-d de la Tarea 1r con 200 ventanas de evaluacion.
    python analyze_1r200.py          (curva exacta si existe task1r_seed*_r200_exact.json)
    python analyze_1r200.py coarse   (barrido registrado)
Baselines: GaMMA con semilla de la Tarea 1s (mismas ventanas) y PyOcto de task1r_baselines --methods PYOCTO."""
import json
import os
import sys

import numpy as np

CAPS = ("0.1", "0.25")
LAB = {"unilateral_c40": "one-sided, 8 st.", "ladder_4": "4 stations", "ladder_6": "6 stations"}
B = 2000
COARSE = len(sys.argv) > 1 and sys.argv[1] == "coarse"


def op(m, cap):
    H, F = np.asarray(m["H"]), np.asarray(m["F"])
    fr = F.sum(0) / H.shape[0]
    ok = np.where(fr <= cap)[0]
    if len(ok) == 0:
        return np.zeros(H.shape[0])
    return H[:, ok[np.argmax(H.sum(0)[ok])]].astype(float)


def boot(ha, hb, k, rng):
    n = len(k)
    ds = [(ha[i].sum() - hb[i].sum()) / max(1, k[i].sum()) for i in (rng.integers(0, n, n) for _ in range(B))]
    return (ha.sum() - hb.sum()) / max(1, k.sum()), np.percentile(ds, 2.5), np.percentile(ds, 97.5)


out, Ra, Rb, Rd = [], {}, {}, {}
for seed in (2024, 7331):
    fx = f"task1r_seed{seed}_r200_exact.json"
    f = fx if (os.path.exists(fx) and not COARSE) else f"task1r_seed{seed}_r200.json"
    out.append(f"[{f}]")
    R = {r["scenario"]: r for r in json.load(open(f))}
    S = {r["scenario"]: r for r in json.load(open(f"task1s_seed{seed}_r200.json"))}
    fp = f"task1r_baselines_seed{seed}_r200_pyocto.json"
    P = {r["scenario"]: r for r in json.load(open(fp))} if os.path.exists(fp) else {}
    rng = np.random.default_rng(seed)
    na = nb = nd = 0
    excl = False
    for sc in LAB:
        r = R[sc]
        k = np.asarray(r["k_eval"], float)
        assert r["k_eval"] == S[sc]["k_eval"]
        ra = {m: r["methods"][m]["recall_at"] for m in ("KIN-LR", "PHYS-LR", "PHYS-SH", "WAVE-P", "TOL-TB")}
        ra["GAMMA"] = S[sc]["methods"]["GAMMA"]["recall_at"]
        if sc in P:
            ra["PYOCTO"] = P[sc]["methods"]["PYOCTO"]["recall_at"]
        out.append(f"{seed} {LAB[sc]:17s} eventos {int(k.sum()):4d} :: " + " | ".join(
            f"{m} " + "/".join(f"{ra[m][c]:.2f}" for c in CAPS) for m in ra))
        out.append("     variantes: " + ", ".join(f"{m} {r['methods'][m]['variant']}" for m in ("KIN-LR", "PHYS-LR", "PHYS-SH")))
        na += all(ra["PHYS-LR"][c] - ra["KIN-LR"][c] >= 0.05 for c in CAPS)
        nb += abs(ra["PHYS-SH"]["0.25"] - ra["KIN-LR"]["0.25"]) < 0.05
        best = {c: max(ra[m][c] for m in ("GAMMA", "PYOCTO") if m in ra) for c in CAPS}
        nd += all(ra["PHYS-LR"][c] - best[c] >= 0.05 for c in CAPS)
        hg = S[sc]["methods"]["GAMMA"]["chosen"]["0.25"]
        hg = np.asarray(hg[0], float) if hg else np.zeros(len(k))
        for a, b, ha, hb in (("PHYS-LR", "KIN-LR", op(r["methods"]["PHYS-LR"], 0.25), op(r["methods"]["KIN-LR"], 0.25)),
                             ("PHYS-LR", "PHYS-SH", op(r["methods"]["PHYS-LR"], 0.25), op(r["methods"]["PHYS-SH"], 0.25)),
                             ("PHYS-SH", "KIN-LR", op(r["methods"]["PHYS-SH"], 0.25), op(r["methods"]["KIN-LR"], 0.25)),
                             ("PHYS-LR", "GAMMA", op(r["methods"]["PHYS-LR"], 0.25), hg)):
            d, lo, hi = boot(ha, hb, k, rng)
            out.append(f"     {a} - {b} @0.25: {d:+.3f} [{lo:+.3f}, {hi:+.3f}]")
            if (a, b) == ("PHYS-LR", "KIN-LR") and lo > 0:
                excl = True
    Ra[seed], Rb[seed], Rd[seed] = na >= 2 and excl, nb >= 2, nd >= 2
    out.append(f"   seed {seed}: R-a celdas {na}/3 (IC excluye 0: {excl}); R-b celdas {nb}/3; R-d celdas {nd}/3")
out.append(f"R-a: {all(Ra.values())}  R-b: {all(Rb.values())}  R-d: {all(Rd.values())}")
txt = "\n".join(out)
print(txt)
open("analyze_1r200" + ("_coarse" if COARSE or "exact" not in out[0] else "_exact") + "_out.txt", "w").write(txt + "\n")
