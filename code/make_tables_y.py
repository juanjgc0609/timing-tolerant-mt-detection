"""make_tables_y.py - tablas Y.7 a Y.10 del Apendice Y desde los JSON (sin transcripcion manual)."""
import json
import os

import numpy as np

LAB = {"unilateral_c40": "one-sided, 8 st.", "ladder_4": "4 stations", "ladder_6": "6 stations"}
C = ("0.1", "0.25")
SEEDS = (2024, 7331)
B = 2000


def f2(ra):
    return "/".join(f"{ra[c]:.2f}" for c in C)


def cell(vals):
    return " ; ".join(vals)


def op(m, cap):
    H, F = np.asarray(m["H"]), np.asarray(m["F"])
    fr = F.sum(0) / H.shape[0]
    ok = np.where(fr <= cap)[0]
    if len(ok) == 0:
        return np.zeros(H.shape[0])
    return H[:, ok[np.argmax(H.sum(0)[ok])]].astype(float)


def _best(H, F, i, ntrue, cap):
    fr = F[i].sum(0) / len(i)
    ok = fr <= cap
    return float((H[i].sum(0)[ok]).max()) / max(1, ntrue) if ok.any() else 0.0


def boot(ma, mb, k, seed, cap=0.25):
    """IC bootstrap pareado por ventanas; en cada replica el punto de operacion de cada metodo de
    grilla se re-elige al tope (la verificacion mostro que fijarlo estrecha el intervalo). ma y mb son
    dicts de metodo con H y F, o un vector por ventana (GaMMA en su configuracion oraculo, fija)."""
    rng = np.random.default_rng(seed)
    n = len(k)
    mats = []
    for m in (ma, mb):
        if isinstance(m, dict):
            mats.append((np.asarray(m["H"], float), np.asarray(m["F"], float)))
        else:
            mats.append(np.asarray(m, float))

    def val(x, i):
        nt = k[i].sum()
        if isinstance(x, tuple):
            return _best(x[0], x[1], i, nt, cap)
        return x[i].sum() / max(1, nt)
    allw = np.arange(n)
    d0 = val(mats[0], allw) - val(mats[1], allw)
    ds = []
    for _ in range(B):
        i = rng.integers(0, n, n)
        ds.append(val(mats[0], i) - val(mats[1], i))
    return f"{d0:+.2f} [{np.percentile(ds, 2.5):+.2f}, {np.percentile(ds, 97.5):+.2f}]"


def load(pat):
    return {s: {r["scenario"]: r for r in json.load(open(pat.format(s=s)))} for s in SEEDS if os.path.exists(pat.format(s=s))}


R = load("task1r_seed{s}_r200_exact.json") or load("task1r_seed{s}_r200.json")
RC = load("task1r_seed{s}_r200.json")
S = load("task1s_seed{s}_r200.json")
P = load("task1r_baselines_seed{s}_r200_pyocto.json")
T = load("task1t_seed{s}_r200_exact.json") or load("task1t_seed{s}_r200.json")
UW = load("task1u_seed{s}_r100_wave_exact.json") or load("task1u_seed{s}_r100.json")
UG = load("task1u_seed{s}_r100.json")
L = []


def gam(s, sc):
    return S[s][sc]["methods"]["GAMMA"]["recall_at"]


def gam_h(s, sc):
    ch = S[s][sc]["methods"]["GAMMA"]["chosen"]["0.25"]
    return np.asarray(ch[0], float) if ch else np.zeros(len(S[s][sc]["k_eval"]))


L.append("TABLE Y.7 (1r, 200 windows): recall at realized false events per window <=0.1 / <=0.25; seed 2024 ; seed 7331")
L.append("| Cell | events | KIN-LR | PHYS-LR | PHYS-SH | GaMMA | PyOcto |")
L.append("|:--|:--|:--|:--|:--|:--|:--|")
for sc in LAB:
    ev = cell([str(sum(R[s][sc]["k_eval"])) for s in SEEDS])
    row = [LAB[sc], ev]
    for m in ("KIN-LR", "PHYS-LR", "PHYS-SH"):
        row.append(cell([f2(R[s][sc]["methods"][m]["recall_at"]) for s in SEEDS]))
    row.append(cell([f2(gam(s, sc)) for s in SEEDS]))
    row.append(cell([f2(P[s][sc]["methods"]["PYOCTO"]["recall_at"]) if s in P else "-" for s in SEEDS]))
    L.append("| " + " | ".join(row) + " |")
L.append("")
L.append("TABLE Y.7b differences at realized false <=0.25, paired bootstrap 95% CI over windows; seed 2024 ; seed 7331")
L.append("| Cell | PHYS-LR - KIN-LR | PHYS-LR - PHYS-SH | PHYS-SH - KIN-LR | PHYS-LR - GaMMA |")
L.append("|:--|:--|:--|:--|:--|")
for sc in LAB:
    row = [LAB[sc]]
    for a, b in (("PHYS-LR", "KIN-LR"), ("PHYS-LR", "PHYS-SH"), ("PHYS-SH", "KIN-LR")):
        row.append(cell([boot(R[s][sc]["methods"][a], R[s][sc]["methods"][b],
                              np.asarray(R[s][sc]["k_eval"], float), s) for s in SEEDS]))
    row.append(cell([boot(R[s][sc]["methods"]["PHYS-LR"], gam_h(s, sc),
                          np.asarray(R[s][sc]["k_eval"], float), s) for s in SEEDS]))
    L.append("| " + " | ".join(row) + " |")
L.append("")
L.append("TABLE Y.8 (1s): recall at realized false <=0.05 / <=0.1 / <=0.25; AUC of the polarization score between true and false GaMMA events; seed 2024 ; seed 7331")
L.append("| Cell | GaMMA | GaMMA+POL | GaMMA+SH | AUC POL | AUC SH | GaMMA ceiling |")
L.append("|:--|:--|:--|:--|:--|:--|:--|")
for sc in LAB:
    row = [LAB[sc]]
    for m in ("GAMMA", "GAMMA+POL", "GAMMA+SH"):
        row.append(cell(["/".join(f"{S[s][sc]['methods'][m]['recall_at'][c]:.2f}" for c in ("0.05",) + C) for s in SEEDS]))
    row.append(cell([f"{S[s][sc]['auc_pol']:.2f}" for s in SEEDS]))
    row.append(cell([f"{S[s][sc]['auc_sh']:.2f}" for s in SEEDS]))
    row.append(cell([f"{max(p[1] for p in S[s][sc]['methods']['GAMMA']['points']):.2f}" for s in SEEDS]))
    L.append("| " + " | ".join(row) + " |")
L.append("")
L.append("TABLE Y.9 (1t, 200 windows): recall at realized false <=0.1 / <=0.25; seed 2024 ; seed 7331")
L.append("| Cell | WAVE-E | WAVE-P | WAVE-PSH | GaMMA | PyOcto | WAVE-P - WAVE-E | WAVE-P - WAVE-PSH | WAVE-P - GaMMA |")
L.append("|:--|:--|:--|:--|:--|:--|:--|:--|:--|")
for sc in LAB:
    row = [LAB[sc]]
    for m in ("WAVE-E", "WAVE-P", "WAVE-PSH"):
        row.append(cell([f2(T[s][sc]["methods"][m]["recall_at"]) for s in SEEDS]))
    row.append(cell([f2(gam(s, sc)) for s in SEEDS]))
    row.append(cell([f2(P[s][sc]["methods"]["PYOCTO"]["recall_at"]) if s in P else "-" for s in SEEDS]))
    for b in ("WAVE-E", "WAVE-PSH"):
        row.append(cell([boot(T[s][sc]["methods"]["WAVE-P"], T[s][sc]["methods"][b],
                              np.asarray(T[s][sc]["k_eval"], float), s) for s in SEEDS]))
    row.append(cell([boot(T[s][sc]["methods"]["WAVE-P"], gam_h(s, sc),
                          np.asarray(T[s][sc]["k_eval"], float), s) for s in SEEDS]))
    L.append("| " + " | ".join(row) + " |")
L.append("")
L.append("TABLE Y.10 (1u, 100 windows, unmodelled statics 0.15 s and 10 deg): recall at realized false <=0.1 / <=0.25; seed 2024 ; seed 7331")
L.append("| Cell | WAVE-E | WAVE-P | WAVE-PSH | GaMMA | WAVE-P - GaMMA |")
L.append("|:--|:--|:--|:--|:--|:--|")


def gu_h(s, sc):
    m = UG[s][sc]["methods"]["GAMMA"]
    ok = [(rc, ci) for ci, (fr, rc) in enumerate(m["points"]) if fr <= 0.25]
    return np.asarray(m["per"][str(max(ok)[1])][0], float) if ok else np.zeros(len(UG[s][sc]["k_eval"]))


if UW:
    for sc in LAB:
        row = [LAB[sc]]
        for m in ("WAVE-E", "WAVE-P", "WAVE-PSH"):
            row.append(cell([f2(UW[s][sc]["methods"][m]["recall_at"]) for s in SEEDS]))
        row.append(cell([f2(UG[s][sc]["methods"]["GAMMA"]["recall_at"]) for s in SEEDS]))
        row.append(cell([boot(UW[s][sc]["methods"]["WAVE-P"], gu_h(s, sc),
                              np.asarray(UW[s][sc]["k_eval"], float), s) for s in SEEDS]))
        L.append("| " + " | ".join(row) + " |")
txt = "\n".join(L)
print(txt)
open("tables_y.md", "w").write(txt + "\n")
