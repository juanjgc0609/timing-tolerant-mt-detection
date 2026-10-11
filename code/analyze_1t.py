"""analyze_1t.py - reglas T-a, T-b y T-c de la Tarea 1t.  python analyze_1t.py 200"""
import json
import os
import sys

import numpy as np

REPS = sys.argv[1] if len(sys.argv) > 1 else "200"
CAPS = ("0.1", "0.25")
B = 2000
LAB = {"unilateral_c40": "one-sided, 8 st.", "ladder_4": "4 stations", "ladder_6": "6 stations"}


def op_wave(m, cap):
    """columna del barrido con mayor recall a falsos realizados <= cap (punto de operacion fijo)."""
    H, F = np.asarray(m["H"]), np.asarray(m["F"])
    fr = F.sum(0) / H.shape[0]
    ok = np.where(fr <= cap)[0]
    if len(ok) == 0:
        return np.zeros(H.shape[0])
    j = ok[np.argmax(H.sum(0)[ok])]
    return H[:, j].astype(float)


def load_base(seed):
    """baselines con 200 ventanas: GaMMA con semilla de la Tarea 1s (mismas ventanas, catalogo
    identico al de task1r_baselines con semilla) y PyOcto de task1r_baselines --methods PYOCTO."""
    base = {}
    f1 = f"task1s_seed{seed}_r{REPS}.json"
    if os.path.exists(f1):
        for r in json.load(open(f1)):
            g = r["methods"]["GAMMA"]
            base.setdefault(r["scenario"], {})["GAMMA"] = {
                "recall_at": g["recall_at"],
                "h": {c: (np.asarray(g["chosen"][c][0], float) if g["chosen"][c] else None) for c in g["chosen"]}}
    f2 = f"task1r_baselines_seed{seed}_r{REPS}_pyocto.json"
    if os.path.exists(f2):
        for r in json.load(open(f2)):
            m = r["methods"]["PYOCTO"]
            h = {}
            for c in ("0.05",) + CAPS:
                ok = [(rc, ci) for ci, (fr, rc) in enumerate(m["points"]) if fr <= float(c)]
                h[c] = np.asarray(m["per"][str(max(ok)[1])][0], float) if ok else None
            base.setdefault(r["scenario"], {})["PYOCTO"] = {"recall_at": m["recall_at"], "h": h}
    return base or None


def boot(ha, hb, k, rng):
    d0 = (ha.sum() - hb.sum()) / max(1, k.sum())
    n = len(k)
    ds = [(ha[i].sum() - hb[i].sum()) / max(1, k[i].sum()) for i in (rng.integers(0, n, n) for _ in range(B))]
    return d0, np.percentile(ds, 2.5), np.percentile(ds, 97.5)


out = []
ok_a, ok_b, ok_c = {}, {}, {}
for seed in (2024, 7331):
    fx = f"task1t_seed{seed}_r{REPS}_exact.json"
    R = json.load(open(fx if os.path.exists(fx) else f"task1t_seed{seed}_r{REPS}.json"))
    out.append(f"[{fx if os.path.exists(fx) else 'barrido registrado'}]")
    base = load_base(seed)
    rng = np.random.default_rng(seed)
    na = nc = 0
    excl = False
    b_ok = True
    for r in R:
        k = np.asarray(r["k_eval"], float)
        ra = {m: r["methods"][m]["recall_at"] for m in r["methods"]}
        line = f"{seed} {LAB[r['scenario']]:17s} eventos {int(k.sum()):4d} :: " + " | ".join(
            f"{m} " + "/".join(f"{ra[m][c]:.2f}" for c in CAPS) for m in ra)
        best = None
        if base is not None:
            bm = base[r["scenario"]]
            bb = {c: max(bm[x]["recall_at"][c] for x in bm) for c in CAPS}
            best = bb
            line += " || " + " ".join(f"{x} " + "/".join(f"{bm[x]['recall_at'][c]:.2f}" for c in CAPS) for x in bm)
        out.append(line)
        if "recall_at_coarse" in r["methods"]["WAVE-P"]:
            out.append("     barrido registrado: " + " | ".join(
                f"{m} " + "/".join(f"{r['methods'][m]['recall_at_coarse'][c]:.2f}" for c in CAPS) for m in ra))
        a_cell = all(ra["WAVE-P"][c] - ra["WAVE-E"][c] >= 0.05 for c in CAPS)
        na += a_cell
        if a_cell and ra["WAVE-P"]["0.25"] - ra["WAVE-PSH"]["0.25"] < 0.05:
            b_ok = False
        for a, b in (("WAVE-P", "WAVE-E"), ("WAVE-P", "WAVE-PSH"), ("WAVE-PSH", "WAVE-E")):
            d, lo, hi = boot(op_wave(r["methods"][a], 0.25), op_wave(r["methods"][b], 0.25), k, rng)
            out.append(f"     {a} - {b} @0.25: {d:+.3f} [{lo:+.3f}, {hi:+.3f}]")
            if (a, b) == ("WAVE-P", "WAVE-E") and lo > 0:
                excl = True
        if base is not None:
            nc += all(ra["WAVE-P"][c] - best[c] >= 0.05 for c in CAPS)
            bm = base[r["scenario"]]
            name = max(bm, key=lambda x: bm[x]["recall_at"]["0.25"])
            hb = bm[name]["h"]["0.25"]
            d, lo, hi = boot(op_wave(r["methods"]["WAVE-P"], 0.25), hb if hb is not None else np.zeros(len(k)), k, rng)
            out.append(f"     WAVE-P - {name} @0.25: {d:+.3f} [{lo:+.3f}, {hi:+.3f}]")
    ok_a[seed] = na >= 2 and excl
    ok_b[seed] = b_ok
    ok_c[seed] = nc >= 2 if base is not None else None
    out.append(f"   seed {seed}: T-a celdas {na}/3, IC excluye 0: {excl}; T-b: {b_ok}; T-c celdas {nc}/3")
out.append(f"T-a: {all(ok_a.values())}  T-b: {all(ok_b.values())}  T-c: {ok_c}")
txt = "\n".join(out)
print(txt)
open(f"analyze_1t_r{REPS}" + ("_exact" if "exact" in out[0] else "") + "_out.txt", "w").write(txt + "\n")
