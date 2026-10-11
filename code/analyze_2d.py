"""analyze_2d.py - reglas 2D-a, 2D-b y 2D-c de la Tarea 2D (PREREGISTRO_1j_1k.txt), por regimen de clutter, y tablas.
Escrito antes de ver cualquier resultado de 2D en las semillas 2024 y 7331.

Datos. Parte (i): MT-FW+SC y QM-3C+SC de las redes 0-4 de 2A (task2a_parts, 20 ventanas por red); BOMP+SC, MTSTACK+SC y
GL-FW+SC de 2C (task2c_parts/p1_*, referencia descriptiva); BOMP-PK+SC, MTSTACK-PK+SC, GL-FW-PK+SC, MT-FW-EX+SC y
MT-FW-NOCAP+SC de 2D (task2d_parts/p1_*); se comprueba ventana a ventana que verdad, origenes y k coinciden.
Parte (ii): maximos por ventana vacia de MT-FW+SC y QM-3C+SC en 200 ventanas vacias por celda y regimen de la semilla
999 (task2d_parts/p2_seed999_*).

Endpoint de la parte (i) como analyze_2a.py / analyze_2c.py: recall de eventos a eventos falsos realizados por ventana
<= 0.1 y <= 0.25 en la curva FROC exacta, punto de operacion elegido sobre la propia evaluacion; una diferencia cuenta
si es >= 0.05 - 1e-9.  Las estimaciones puntuales deciden 2D-a y 2D-b (la convencion de 2A-2C).
2D-a: MT-FW+SC - max(BOMP-PK+SC, MTSTACK-PK+SC) >= 0.05 a tope 0.25 en >= 2 de 3 celdas en ambas semillas.
2D-b: MT-FW+SC - MT-FW-EX+SC >= 0.05 a tope 0.25 en >= 2 de 3 celdas en ambas semillas.
Descriptivo: MT-FW-NOCAP+SC - MT-FW+SC y GL-FW-PK+SC - MT-FW+SC con IC; ganancia de la lectura en el pico de cada
comparador (X-PK+SC - X+SC).
2D-c: umbral split-conformal por celda, regimen y alfa (0.1, 0.25) sobre los n = 200 maximos por ventana vacia de la
semilla 999: q = el ceil((n+1)(1-alfa))-esimo menor (q = +inf si ese indice supera n; una ventana sin ningun evento
tiene maximo -inf); se declara un evento si su puntaje es > q (desigualdad estricta, la de la garantia
P(maximo de una ventana vacia nueva > q) <= alfa).  Umbral congelado y aplicado a las ventanas de evaluacion de 2A:
recall, falsos por ventana vacia (k = 0), fraccion de ventanas vacias con algun falso y falsos por ventana en todas las
ventanas.  "Tasa de falsos realizada ... en ventanas vacias" se lee como eventos falsos por ventana vacia (la definicion
de tasa de falsos de 2A; es >= la fraccion de ventanas vacias con algun falso, que es la cantidad que controla la
garantia conformal; ambas se reportan).  Regla, lectura primaria (por celda): una celda cuenta si la tasa de MT-FW+SC en
ventanas vacias es <= 2 alfa Y recall(MT-FW+SC) - recall(QM-3C+SC) >= 0.05; >= 2 de 3 celdas en ambas semillas, a
ambos alfa.  Lectura estricta (se reporta tambien): tasa <= 2 alfa en todas las celdas, semillas y alfa, y diferencia
>= 0.05 en >= 2 de 3 celdas en ambas semillas y alfa.
Estadistica: bootstrap jerarquico (2000 replicas; en cada una las 5 redes de la celda con reemplazo y, dentro de cada red
elegida, sus ventanas con reemplazo; default_rng([s, ci, ri, 2704])), las mismas remuestras para todos los metodos de la
celda; en la parte (i) se re-eligen el punto de operacion y el maximo sobre comparadores en cada replica; en 2D-c el
umbral queda fijo (la incertidumbre de la calibracion no entra en el IC).  Como lectura secundaria, bootstrap pareado
por ventanas (default_rng([s, ci, ri]), las remuestras de analyze_2a.py).  Holm: p unilateral p = (1 + #{d* <= 0}) /
(B + 1) con el bootstrap jerarquico, nivel 0.025, familia por regimen = 2D-a (diferencia con el maximo) y 2D-b a ambos
topes y 2D-c (diferencia de recall) a ambos alfa, en 3 celdas y 2 semillas: 36 comparaciones.
Exploratorio (declarado): fraccion de replicas en que se cumplirian las reglas; 2D-c con los umbrales congelados sobre
las 30 redes de 2C (2A + las 25 redes nuevas de la parte (ii) de 2C).

    python analyze_2d.py [--out analyze_2d_out.txt]
"""
import argparse
import json
import math
import os

import numpy as np

from analyze_2c import HF, counts_windows, holm, point, val

HERE = os.path.dirname(os.path.abspath(__file__))
B = 2000
CAPS = ("0.1", "0.25")
ALPHAS = (0.1, 0.25)
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
LAB = {"ladder_4": "4 st.", "ladder_6": "6 st.", "unilateral_c40": "one-sided 8"}
REGIMES = ("thin", "strong")
SEEDS = (2024, 7331)
EPS = 1e-9
THR = 0.05
ALPHA_HOLM = 0.025
M2D = ("BOMP-PK+SC", "MTSTACK-PK+SC", "GL-FW-PK+SC", "MT-FW-EX+SC", "MT-FW-NOCAP+SC")
M2C = ("BOMP+SC", "MTSTACK+SC", "GL-FW+SC")
M2A = ("MT-FW+SC", "QM-3C+SC")
M_TAB = ("MT-FW+SC",) + M2D + M2C + ("QM-3C+SC",)
M_CONF = ("MT-FW+SC", "QM-3C+SC")


def ge(d):
    return d >= THR - EPS


def fmt_ci(p, lo, hi):
    return f"{p:+.3f} [{lo:+.3f}, {hi:+.3f}]"


def load(path):
    with open(path) as fh:
        return json.load(fh)


def counts_hier(s, ci, ri, net, tag=2704):
    rng = np.random.default_rng([s, ci, ri, tag])
    ids = sorted(set(net))
    pos = {j: np.where(net == j)[0] for j in ids}
    C = np.zeros((B, len(net)))
    for b in range(B):
        for jj in rng.integers(0, len(ids), len(ids)):
            p = pos[ids[jj]]
            np.add.at(C[b], p[rng.integers(0, len(p), len(p))], 1.0)
    return C


# ------------------------------------------------------------------------------------------ conformal
def conformal_q(maxes, alpha):
    """umbral split-conformal: el ceil((n+1)(1-alpha))-esimo menor de los n maximos (+inf si el indice supera n)."""
    x = np.sort(np.asarray(maxes, float))
    n = len(x)
    kq = int(math.ceil((n + 1) * (1.0 - alpha) - 1e-12))
    return (float(x[kq - 1]) if kq <= n else float("inf")), kq, n


def apply_q(recs, q):
    """eventos con puntaje > q: aciertos (fuentes distintas) y falsos por ventana."""
    h, f = np.zeros(len(recs)), np.zeros(len(recs))
    for w, r in enumerate(recs):
        hits = set()
        for s, hh in r:
            if s > q:
                if hh >= 0:
                    hits.add(hh)
                else:
                    f[w] += 1
        h[w] = len(hits)
    return h, f


def conf_stats(h, f, k):
    empty = k == 0
    return {"recall": float(h.sum() / max(1.0, k.sum())), "fr_empty": float(f[empty].mean()) if empty.any() else 0.0,
            "frac_empty_any": float((f[empty] > 0).mean()) if empty.any() else 0.0, "fr_all": float(f.mean()),
            "n_empty": int(empty.sum())}


# ------------------------------------------------------------------------------------------ data
def cell_data(s, reg, cell):
    W = {"k": [], "net": [], "rec": {m: [] for m in M_TAB}, "gl": [], "shifts": [], "orth": []}
    ctl = {"windows": 0, "mismatch": 0, "controls": [], "kernels": set(), "kernel_required": []}
    for j in range(5):
        a = load(os.path.join(HERE, "task2a_parts", f"seed{s}_{reg}_{cell}_net{j}.json"))
        c = load(os.path.join(HERE, "task2c_parts", f"p1_seed{s}_{reg}_{cell}_net{j}.json"))
        d = load(os.path.join(HERE, "task2d_parts", f"p1_seed{s}_{reg}_{cell}_net{j}.json"))
        ctl["controls"].append(d["control"])
        ctl["kernels"].add(tuple(tuple(x) for x in d["kernel"]["blas"]))
        ctl["kernel_required"].append(d["kernel_required"])
        assert a["n_win"] == c["n_win"] == d["n_win"] == len(d["windows"])
        for wa, wc, wd in zip(a["windows"], c["windows"], d["windows"]):
            same = all(wa[x] == wc[x] == wd[x] for x in ("k", "truth", "taus"))
            ctl["windows"] += 1
            ctl["mismatch"] += not same
            W["k"].append(wa["k"])
            W["net"].append(j)
            for m in M2A:
                W["rec"][m].append(wa["records"][m])
            for m in M2C:
                W["rec"][m].append(wc["records"][m])
            for m in M2D:
                W["rec"][m].append(wd["records"][m])
            W["gl"].append(wd["gl"])
            W["shifts"].extend(wd["bomp_shifts"])
            W["orth"].extend(wd["bomp_orth"])
    W["k"] = np.asarray(W["k"], float)
    W["net"] = np.asarray(W["net"])
    return W, ctl


def cal_data(reg, cell):
    mx = {m: [] for m in M_CONF}
    ctl = []
    kern = set()
    fr_cal = {m: [] for m in M_CONF}
    for j in range(5):
        p = load(os.path.join(HERE, "task2d_parts", f"p2_seed999_{reg}_{cell}_net{j}.json"))
        ctl.append(p["control"])
        kern.add(tuple(tuple(x) for x in p["kernel"]["blas"]))
        for w in p["windows"]:
            for m in M_CONF:
                v = w["max"][m]
                mx[m].append(-np.inf if v is None else float(v))
                fr_cal[m].append(w["records"][m])
    return mx, ctl, kern, fr_cal


def c2_windows(s, reg, cell):
    """exploratorio: las 30 redes de 2C (2A redes 0-4 y 2C parte (ii) redes 5-29), MT-FW+SC y QM-3C+SC."""
    k, net, rec = [], [], {m: [] for m in M_CONF}
    for j in range(5):
        a = load(os.path.join(HERE, "task2a_parts", f"seed{s}_{reg}_{cell}_net{j}.json"))
        for w in a["windows"]:
            k.append(w["k"])
            net.append(j)
            for m in M_CONF:
                rec[m].append(w["records"][m])
    for j in range(5, 30):
        path = os.path.join(HERE, "task2c_parts", f"p2_seed{s}_{reg}_{cell}_net{j}.json")
        if not os.path.exists(path):
            continue
        c = load(path)
        for w in c["windows"]:
            k.append(w["k"])
            net.append(j)
            for m in M_CONF:
                rec[m].append(w["records"][m])
    return np.asarray(k, float), np.asarray(net), rec


# ------------------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "analyze_2d_out.txt"))
    a = ap.parse_args()
    L = []
    P = L.append
    verdicts = {}
    near = []
    P("# Tarea 2D - igualdad de tolerancia temporal, ablaciones de MT-FW y umbral desplegable (analyze_2d.py)")
    P("")
    P("Parte (i): recall de eventos a eventos falsos realizados por ventana <= 0.1 / <= 0.25 (curva FROC exacta, punto de "
      "operacion elegido sobre la propia evaluacion), 100 ventanas (5 redes x 20) por celda. Diferencias: punto [IC 95% H] "
      "(bootstrap jerarquico: redes y ventanas) y W[IC 95% bootstrap por ventanas]; 2000 replicas; punto de operacion y "
      "maximo re-elegidos en cada replica. Parte (ii): umbral split-conformal congelado en 200 ventanas vacias de la "
      "semilla 999 por celda y regimen.")
    data = {(reg, s, cell): cell_data(s, reg, cell) for reg in REGIMES for s in SEEDS for cell in CELLS}
    cal = {(reg, cell): cal_data(reg, cell) for reg in REGIMES for cell in CELLS}
    # ------------------------------------------------------------------------------------ controles
    P("")
    P("## Controles")
    P("")
    for reg in REGIMES:
        for s in SEEDS:
            for cell in CELLS:
                W, ctl = data[(reg, s, cell)]
                ok = all(all(v for k_, v in c.items() if k_ not in ("windows_truth_taus", "mtfw_sc_exact"))
                         and c["mtfw_sc_exact"] == 2 for c in ctl["controls"])
                P(f"- {reg} {s} {LAB[cell]}: 5 redes; red, catalogo, correcciones y MT-FW+SC exacto en las ventanas 0-1 "
                  f"{'correcto' if ok else 'FALLA'}; ventanas emparejadas con 2A y 2C {ctl['windows']} (distintas "
                  f"{ctl['mismatch']}); nucleo requerido {sorted(set(ctl['kernel_required']))}, cargado {sorted(ctl['kernels'])}")
                assert ok and ctl["mismatch"] == 0
    for reg in REGIMES:
        for cell in CELLS:
            mx, ctl, kern, _fr = cal[(reg, cell)]
            rep = [r for c in ctl for r in c["reproduces_2a"]]
            P(f"- calibracion 999 {reg} {LAB[cell]}: {len(mx['MT-FW+SC'])} ventanas vacias; red y retardos iguales a 2A "
              f"{all(c['STA'] and c['statics'] for c in ctl)}; correcciones identicas {sum(c['corr_exact'] for c in ctl)}/5 "
              f"(max |dif| {max(c['corr_maxabs'] for c in ctl):.1e} s); la tuberia reproduce 2A en las ventanas 0-1: "
              f"MT-FW+SC {sum(r['MT-FW+SC'] for r in rep)}/{len(rep)}, QM-3C+SC {sum(r['QM-3C+SC'] for r in rep)}/{len(rep)}; "
              f"nucleo {sorted(kern)}")
    # ------------------------------------------------------------------------------------ parte (i)
    tests_all = {reg: {} for reg in REGIMES}
    pts_all = {reg: {} for reg in REGIMES}
    for ri, reg in enumerate(REGIMES):
        P("")
        P(f"## Regimen {reg}: parte (i), redes 0-4 de 2A")
        P("")
        P(f"### Tabla 1 ({reg}). Recall a <= 0.1 / <= 0.25")
        P("")
        P("| Celda (semilla) | eventos | " + " | ".join(M_TAB) + " |")
        P("|:--|--:|" + ":-:|" * len(M_TAB))
        rowsA, rowsB, rowsD, rowsG = [], [], [], []
        cnt = {"2D-a": {s: 0 for s in SEEDS}, "2D-b": {s: 0 for s in SEEDS}}
        det = {"2D-a": {s: [] for s in SEEDS}, "2D-b": {s: [] for s in SEEDS}}
        rep = {"2D-a": {s: [] for s in SEEDS}, "2D-b": {s: [] for s in SEEDS}}
        for s in SEEDS:
            for ci, cell in enumerate(CELLS):
                W, _ctl = data[(reg, s, cell)]
                k, net = W["k"], W["net"]
                HFs = {m: HF(W["rec"][m]) for m in M_TAB}
                pt = {m: point(HFs[m], k) for m in M_TAB}
                P(f"| {LAB[cell]} ({s}) | {int(k.sum())} | " + " | ".join(
                    f"{pt[m]['0.1']:.3f}/{pt[m]['0.25']:.3f}" for m in M_TAB) + " |")
                Ch = counts_hier(s, ci, ri, net)
                Cw = counts_windows(s, ci, ri, len(k))
                bh = {m: val(HFs[m], Ch, k) for m in M_TAB}
                bw = {m: val(HFs[m], Cw, k) for m in M_TAB}
                ca, cb = [], []
                for cap in CAPS:
                    # 2D-a
                    p0 = pt["MT-FW+SC"][cap] - max(pt["BOMP-PK+SC"][cap], pt["MTSTACK-PK+SC"][cap])
                    dh = bh["MT-FW+SC"][cap] - np.maximum(bh["BOMP-PK+SC"][cap], bh["MTSTACK-PK+SC"][cap])
                    dw = bw["MT-FW+SC"][cap] - np.maximum(bw["BOMP-PK+SC"][cap], bw["MTSTACK-PK+SC"][cap])
                    best = "BOMP-PK+SC" if pt["BOMP-PK+SC"][cap] >= pt["MTSTACK-PK+SC"][cap] else "MTSTACK-PK+SC"
                    pv = (1.0 + np.sum(dh <= 0)) / (B + 1.0)
                    key = ("2D-a", s, cell, cap)
                    tests_all[reg][key] = pv
                    pts_all[reg][key] = (p0, *np.percentile(dh, [2.5, 97.5]))
                    ca.append(f"{best}: {fmt_ci(p0, *np.percentile(dh, [2.5, 97.5]))} W[{np.percentile(dw, 2.5):+.3f}, "
                              f"{np.percentile(dw, 97.5):+.3f}] p={pv:.4f}")
                    if cap == "0.25":
                        oka = ge(p0)
                        cnt["2D-a"][s] += oka
                        det["2D-a"][s].append((cell, oka, p0))
                        rep["2D-a"][s].append(dh >= THR - EPS)
                        if abs(p0 - THR) < 0.005:
                            near.append(f"2D-a {reg} {s} {cell}: {p0:+.4f}")
                    # 2D-b
                    p0 = pt["MT-FW+SC"][cap] - pt["MT-FW-EX+SC"][cap]
                    dh = bh["MT-FW+SC"][cap] - bh["MT-FW-EX+SC"][cap]
                    dw = bw["MT-FW+SC"][cap] - bw["MT-FW-EX+SC"][cap]
                    pv = (1.0 + np.sum(dh <= 0)) / (B + 1.0)
                    key = ("2D-b", s, cell, cap)
                    tests_all[reg][key] = pv
                    pts_all[reg][key] = (p0, *np.percentile(dh, [2.5, 97.5]))
                    cb.append(f"{fmt_ci(p0, *np.percentile(dh, [2.5, 97.5]))} W[{np.percentile(dw, 2.5):+.3f}, "
                              f"{np.percentile(dw, 97.5):+.3f}] p={pv:.4f}")
                    if cap == "0.25":
                        okb = ge(p0)
                        cnt["2D-b"][s] += okb
                        det["2D-b"][s].append((cell, okb, p0))
                        rep["2D-b"][s].append(dh >= THR - EPS)
                        if abs(p0 - THR) < 0.005:
                            near.append(f"2D-b {reg} {s} {cell}: {p0:+.4f}")
                sep = {m: " / ".join(fmt_ci(pt["MT-FW+SC"][cap] - pt[m][cap],
                                            *np.percentile(bh["MT-FW+SC"][cap] - bh[m][cap], [2.5, 97.5])) for cap in CAPS)
                       for m in ("BOMP-PK+SC", "MTSTACK-PK+SC")}
                rowsA.append(f"| {LAB[cell]} ({s}) | {ca[0]} | {ca[1]} | {'si' if ge(pts_all[reg][('2D-a', s, cell, '0.25')][0]) else 'no'} "
                             f"| {sep['BOMP-PK+SC']} | {sep['MTSTACK-PK+SC']} |")
                rowsB.append(f"| {LAB[cell]} ({s}) | {cb[0]} | {cb[1]} | {'si' if ge(pts_all[reg][('2D-b', s, cell, '0.25')][0]) else 'no'} |")
                dd = []
                for m in ("MT-FW-NOCAP+SC", "GL-FW-PK+SC"):
                    dd.append(" / ".join(fmt_ci(pt[m][cap] - pt["MT-FW+SC"][cap],
                                                *np.percentile(bh[m][cap] - bh["MT-FW+SC"][cap], [2.5, 97.5]))
                                         for cap in CAPS))
                rowsD.append(f"| {LAB[cell]} ({s}) | {dd[0]} | {dd[1]} |")
                gg = []
                for mp, me in (("BOMP-PK+SC", "BOMP+SC"), ("MTSTACK-PK+SC", "MTSTACK+SC"), ("GL-FW-PK+SC", "GL-FW+SC"),
                               ("MT-FW+SC", "MT-FW-EX+SC")):
                    gg.append(" / ".join(f"{pt[mp][cap] - pt[me][cap]:+.3f}" for cap in CAPS))
                gg.append(" / ".join(f"{pt['GL-FW-PK+SC'][cap] - pt['BOMP-PK+SC'][cap]:+.3f}" for cap in CAPS))
                gg.append(" / ".join(f"{max(pt['BOMP-PK+SC'][cap], pt['MTSTACK-PK+SC'][cap]) - pt['QM-3C+SC'][cap]:+.3f}"
                                     for cap in CAPS))
                rowsG.append(f"| {LAB[cell]} ({s}) | " + " | ".join(gg) + " |")
        P("")
        P(f"### Tabla 2 ({reg}). 2D-a: MT-FW+SC - max(BOMP-PK+SC, MTSTACK-PK+SC) (mejor comparador por celda y tope) "
          "[IC H] W[IC W] p unilateral H; y contra cada uno (H)")
        P("")
        P("| Celda (semilla) | a 0.1 | a 0.25 | 2D-a (>= 0.05 a 0.25) | - BOMP-PK+SC (0.1 / 0.25) | - MTSTACK-PK+SC (0.1 / 0.25) |")
        P("|:--|:-:|:-:|:-:|:-:|:-:|")
        L.extend(rowsA)
        P("")
        P(f"### Tabla 3 ({reg}). 2D-b: MT-FW+SC - MT-FW-EX+SC [IC H] W[IC W] p unilateral H")
        P("")
        P("| Celda (semilla) | a 0.1 | a 0.25 | 2D-b (>= 0.05 a 0.25) |")
        P("|:--|:-:|:-:|:-:|")
        L.extend(rowsB)
        P("")
        P(f"### Tabla 4 ({reg}, descriptiva). MT-FW-NOCAP+SC - MT-FW+SC y GL-FW-PK+SC - MT-FW+SC a 0.1 / 0.25 [IC H]")
        P("")
        P("| Celda (semilla) | NOCAP - MT-FW+SC | GL-FW-PK+SC - MT-FW+SC |")
        P("|:--|:-:|:-:|")
        L.extend(rowsD)
        P("")
        P(f"### Tabla 5 ({reg}, descriptiva). Efecto de la lectura en el pico (0.1 / 0.25, puntos)")
        P("")
        P("| Celda (semilla) | BOMP-PK - BOMP | MTSTACK-PK - MTSTACK | GL-FW-PK - GL-FW | MT-FW+SC - MT-FW-EX | "
          "GL-FW-PK - BOMP-PK | max(BOMP-PK, MTSTACK-PK) - QM-3C+SC |")
        P("|:--|:-:|:-:|:-:|:-:|:-:|:-:|")
        L.extend(rowsG)
        verdicts[reg] = {}
        for rule in ("2D-a", "2D-b"):
            okr = all(cnt[rule][s] >= 2 for s in SEEDS)
            joint = np.all([np.sum(rep[rule][s], axis=0) >= 2 for s in SEEDS], axis=0)
            P("")
            P(f"- Regla {rule} ({reg}): " + "; ".join(
                f"semilla {s}: {cnt[rule][s]}/3 celdas (" + ", ".join(
                    f"{LAB[c_]} {'si' if o else 'no'} {x:+.3f}" for c_, o, x in det[rule][s]) + ")" for s in SEEDS)
                + f" -> {'SE CUMPLE' if okr else 'NO SE CUMPLE'} (exploratorio: se cumple en {100 * joint.mean():.0f}% "
                f"de las replicas H)")
            verdicts[reg][rule] = okr
            verdicts[reg][rule + "_counts"] = cnt[rule]
            verdicts[reg][rule + "_points"] = {s: [(c_, x) for c_, _o, x in det[rule][s]] for s in SEEDS}
            verdicts[reg][rule + "_replicate_fraction"] = float(joint.mean())
        # diagnosticos
        gl = [g for s in SEEDS for cell in CELLS for g in data[(reg, s, cell)][0]["gl"]]
        pl = [p for g in gl for p in g["per_lambda"]]
        fin = [g["per_lambda"][-1] for g in gl]
        sh = np.array([x for s in SEEDS for cell in CELLS for at in data[(reg, s, cell)][0]["shifts"] for x in at], float)
        orth = [x for s in SEEDS for cell in CELLS for x in data[(reg, s, cell)][0]["orth"]]
        P("")
        P(f"Diagnosticos ({reg}). GL-FW-PK+SC, {len(gl)} ventanas: lambdas resueltos por ventana "
          f"{min(len(g['per_lambda']) for g in gl)}-{max(len(g['per_lambda']) for g in gl)} (mediana "
          f"{np.median([len(g['per_lambda']) for g in gl]):.0f}); parada por registro completo en "
          f"{sum(1 for g in gl if g['stopped'] and g['stopped'][0] == 'record_complete')}; pasos FW por ventana mediana "
          f"{np.median([sum(p[5] for p in g['per_lambda']) for g in gl]):.0f}, max {max(sum(p[5] for p in g['per_lambda']) for g in gl)}; "
          f"lambdas con tope de pasos FW {sum(p[9] for p in pl)}; brecha relativa max "
          f"{max(p[2] for p in pl if np.isfinite(p[2])):.1e} (ultimo lambda: mediana {np.median([p[2] for p in fin]):.1e}); "
          f"max c/lambda del oraculo {max(p[6] for p in pl):.9f}; max |c/lambda - 1| en el soporte {max(p[7] for p in pl):.1e}; "
          f"re-solves {sum(p[8] for p in pl)}; soporte final mediana {np.median([p[4] for p in fin]):.0f}, max "
          f"{max(p[4] for p in fin)}. BOMP-PK+SC: desplazamiento de las muestras de pico del atomo que entra respecto de "
          f"las predichas: {100 * np.mean(sh != 0):.0f}% de las lecturas estacion-fase desplazadas, |desplazamiento| medio "
          f"{np.mean(np.abs(sh)):.2f} muestras; ortogonalidad del reajuste max {max(orth):.1e}.")
    # ------------------------------------------------------------------------------------ parte (ii)
    thr_all = {}
    conf_res = {}
    for ri, reg in enumerate(REGIMES):
        P("")
        P(f"## Regimen {reg}: parte (ii), umbral split-conformal (200 ventanas vacias de la semilla 999 por celda)")
        P("")
        P(f"### Tabla 6 ({reg}). Calibracion: umbrales congelados q (indice k de n) y tasas en la propia calibracion")
        P("")
        P("| Celda | metodo | ventanas vacias con eventos | alfa 0.1: q (k/n), falsos/ventana, frac. con falso | "
          "alfa 0.25: q (k/n), falsos/ventana, frac. con falso |")
        P("|:--|:--|--:|:-:|:-:|")
        for cell in CELLS:
            mx, _ctl, _k, fr_cal = cal[(reg, cell)]
            for m in M_CONF:
                cc = []
                for al in ALPHAS:
                    q, kq, n = conformal_q(mx[m], al)
                    thr_all[(reg, cell, m, al)] = q
                    h, f = apply_q(fr_cal[m], q)
                    cc.append(f"{q:.4g} ({kq}/{n}), {f.mean():.3f}, {(f > 0).mean():.3f}")
                P(f"| {LAB[cell]} | {m} | {int(np.sum(np.isfinite(mx[m])))} | {cc[0]} | {cc[1]} |")
        P("")
        P(f"### Tabla 7 ({reg}). Evaluacion (ventanas de 2A) con el umbral congelado: recall; falsos por ventana vacia; "
          "fraccion de ventanas vacias con falso; falsos por ventana (todas)")
        P("")
        P("| Celda (semilla) | alfa | ventanas vacias | MT-FW+SC | QM-3C+SC | MT-FW+SC - QM-3C+SC recall [IC H] p | "
          "tasa MT-FW+SC vacias [IC H] | celda cuenta |")
        P("|:--|--:|--:|:-:|:-:|:-:|:-:|:-:|")
        cnt = {al: {s: 0 for s in SEEDS} for al in ALPHAS}
        cnt_rec = {al: {s: 0 for s in SEEDS} for al in ALPHAS}
        rate_ok_all = True
        det = {al: {s: [] for s in SEEDS} for al in ALPHAS}
        repc = {al: {s: [] for s in SEEDS} for al in ALPHAS}
        for s in SEEDS:
            for ci, cell in enumerate(CELLS):
                W, _ctl = data[(reg, s, cell)]
                k, net = W["k"], W["net"]
                empty = (k == 0).astype(float)
                Ch = counts_hier(s, ci, ri, net)
                for al in ALPHAS:
                    st = {}
                    hf = {}
                    for m in M_CONF:
                        h, f = apply_q(W["rec"][m], thr_all[(reg, cell, m, al)])
                        hf[m] = (h, f)
                        st[m] = conf_stats(h, f, k)
                    conf_res[(reg, s, cell, al)] = st
                    d0 = st["MT-FW+SC"]["recall"] - st["QM-3C+SC"]["recall"]
                    kk = Ch @ k
                    dstar = (Ch @ hf["MT-FW+SC"][0] - Ch @ hf["QM-3C+SC"][0]) / np.maximum(1.0, kk)
                    fr_star = (Ch @ (hf["MT-FW+SC"][1] * empty)) / np.maximum(1.0, Ch @ empty)
                    pv = (1.0 + np.sum(dstar <= 0)) / (B + 1.0)
                    key = ("2D-c", s, cell, str(al))
                    tests_all[reg][key] = pv
                    pts_all[reg][key] = (d0, *np.percentile(dstar, [2.5, 97.5]))
                    rate_ok = st["MT-FW+SC"]["fr_empty"] <= 2 * al + 1e-12
                    rate_ok_all &= rate_ok
                    okc = rate_ok and ge(d0)
                    cnt[al][s] += okc
                    cnt_rec[al][s] += ge(d0)
                    det[al][s].append((cell, okc, d0, st["MT-FW+SC"]["fr_empty"]))
                    repc[al][s].append((fr_star <= 2 * al + 1e-12) & (dstar >= THR - EPS))
                    if abs(d0 - THR) < 0.005:
                        near.append(f"2D-c {reg} {s} {cell} alfa {al}: {d0:+.4f}")
                    if abs(st["MT-FW+SC"]["fr_empty"] - 2 * al) < 0.005:
                        near.append(f"2D-c tasa {reg} {s} {cell} alfa {al}: {st['MT-FW+SC']['fr_empty']:.4f} vs {2 * al}")
                    fm = lambda x: (f"{x['recall']:.3f}; {x['fr_empty']:.3f}; {x['frac_empty_any']:.3f}; "  # noqa: E731
                                    f"{x['fr_all']:.3f}")
                    P(f"| {LAB[cell]} ({s}) | {al} | {st['MT-FW+SC']['n_empty']} | {fm(st['MT-FW+SC'])} | "
                      f"{fm(st['QM-3C+SC'])} | {fmt_ci(d0, *np.percentile(dstar, [2.5, 97.5]))} p={pv:.4f} | "
                      f"{st['MT-FW+SC']['fr_empty']:.3f} [{np.percentile(fr_star, 2.5):.3f}, {np.percentile(fr_star, 97.5):.3f}] "
                      f"| {'si' if okc else 'no'} |")
        okC = all(cnt[al][s] >= 2 for al in ALPHAS for s in SEEDS)
        okC_strict = rate_ok_all and all(cnt_rec[al][s] >= 2 for al in ALPHAS for s in SEEDS)
        joint = np.all([np.sum(repc[al][s], axis=0) >= 2 for al in ALPHAS for s in SEEDS], axis=0)
        P("")
        P(f"- Regla 2D-c ({reg}), lectura por celda (tasa en vacias <= 2 alfa y diferencia de recall >= 0.05 en la misma "
          "celda, >= 2 de 3 celdas): " + "; ".join(
              f"alfa {al} semilla {s}: {cnt[al][s]}/3 (" + ", ".join(
                  f"{LAB[c_]} {'si' if o else 'no'} {d:+.3f}, tasa {r:.3f}" for c_, o, d, r in det[al][s]) + ")"
              for al in ALPHAS for s in SEEDS) + f" -> {'SE CUMPLE' if okC else 'NO SE CUMPLE'}")
        P(f"- Regla 2D-c ({reg}), lectura estricta (tasa <= 2 alfa en todas las celdas, semillas y alfa; diferencia >= 0.05 "
          f"en >= 2 de 3): tasa en todas {'si' if rate_ok_all else 'no'}; diferencias " + ", ".join(
              f"alfa {al} semilla {s}: {cnt_rec[al][s]}/3" for al in ALPHAS for s in SEEDS)
          + f" -> {'SE CUMPLE' if okC_strict else 'NO SE CUMPLE'}")
        P(f"- Exploratorio: la lectura por celda se cumple en {100 * joint.mean():.0f}% de las replicas H (umbral fijo).")
        verdicts[reg]["2D-c"] = okC
        verdicts[reg]["2D-c_strict"] = okC_strict
        verdicts[reg]["2D-c_counts"] = {str(al): cnt[al] for al in ALPHAS}
        verdicts[reg]["2D-c_replicate_fraction"] = float(joint.mean())
        # Holm
        rej, adj = holm(tests_all[reg], ALPHA_HOLM)
        P("")
        P(f"### Tabla 8 ({reg}). Holm (familia de {len(tests_all[reg])} comparaciones del regimen, p unilateral H, nivel "
          f"{ALPHA_HOLM}): p ajustado, * = rechazo (MT-FW+SC mayor)")
        P("")
        P("| Comparacion | Celda (semilla) | tope/alfa 0.1 | tope/alfa 0.25 |")
        P("|:--|:--|:-:|:-:|")
        cntH = {r_: {s: 0 for s in SEEDS} for r_ in ("2D-a", "2D-b")}
        cntHc = {al: {s: 0 for s in SEEDS} for al in ALPHAS}
        for rule in ("2D-a", "2D-b", "2D-c"):
            for s in SEEDS:
                for cell in CELLS:
                    keys = [(rule, s, cell, c_) for c_ in (CAPS if rule != "2D-c" else ("0.1", "0.25"))]
                    P(f"| {rule} | {LAB[cell]} ({s}) | " + " | ".join(
                        f"{pts_all[reg][k_][0]:+.3f} pHolm={adj[k_]:.4f}{' *' if rej[k_] else ''}" for k_ in keys) + " |")
                    if rule != "2D-c":
                        k_ = (rule, s, cell, "0.25")
                        cntH[rule][s] += rej[k_] and ge(pts_all[reg][k_][0])
                    else:
                        for al in ALPHAS:
                            k_ = (rule, s, cell, str(al))
                            fr = conf_res[(reg, s, cell, al)]["MT-FW+SC"]["fr_empty"]
                            cntHc[al][s] += rej[k_] and ge(pts_all[reg][k_][0]) and fr <= 2 * al + 1e-12
        P("")
        for rule in ("2D-a", "2D-b"):
            okh = all(cntH[rule][s] >= 2 for s in SEEDS)
            P(f"- {rule} ({reg}) con Holm (rechazo y punto >= 0.05 a 0.25): " + "; ".join(
                f"semilla {s}: {cntH[rule][s]}/3" for s in SEEDS) + f" -> {'se cumple' if okh else 'no se cumple'}")
            verdicts[reg][rule + "_holm"] = okh
        okh = all(cntHc[al][s] >= 2 for al in ALPHAS for s in SEEDS)
        P(f"- 2D-c ({reg}) con Holm (lectura por celda, con rechazo de Holm de la diferencia de recall): " + "; ".join(
            f"alfa {al} semilla {s}: {cntHc[al][s]}/3" for al in ALPHAS for s in SEEDS)
          + f" -> {'se cumple' if okh else 'no se cumple'}")
        verdicts[reg]["2D-c_holm"] = okh
        verdicts[reg]["holm"] = {"|".join(map(str, k_)): [bool(rej[k_]), float(adj[k_]), float(tests_all[reg][k_])]
                                 for k_ in tests_all[reg]}
    # ------------------------------------------------------------------------------------ exploratorio: 30 redes de 2C
    P("")
    P("## Exploratorio: 2D-c con los mismos umbrales congelados sobre las 30 redes de 2C por celda (2A + 25 redes nuevas)")
    P("")
    P("| Regimen | Celda (semilla) | alfa | redes | ventanas vacias | recall MT-FW+SC / QM-3C+SC | dif [IC H] | "
      "tasa MT-FW+SC vacias [IC H] | frac. vacias con falso |")
    P("|:--|:--|--:|--:|--:|:-:|:-:|:-:|:-:|")
    for ri, reg in enumerate(REGIMES):
        for s in SEEDS:
            for ci, cell in enumerate(CELLS):
                k, net, rec = c2_windows(s, reg, cell)
                empty = (k == 0).astype(float)
                Ch = counts_hier(s, ci, ri, net, tag=2705)
                for al in ALPHAS:
                    hf = {m: apply_q(rec[m], thr_all[(reg, cell, m, al)]) for m in M_CONF}
                    st = {m: conf_stats(*hf[m], k) for m in M_CONF}
                    d0 = st["MT-FW+SC"]["recall"] - st["QM-3C+SC"]["recall"]
                    dstar = (Ch @ hf["MT-FW+SC"][0] - Ch @ hf["QM-3C+SC"][0]) / np.maximum(1.0, Ch @ k)
                    fr_star = (Ch @ (hf["MT-FW+SC"][1] * empty)) / np.maximum(1.0, Ch @ empty)
                    P(f"| {reg} | {LAB[cell]} ({s}) | {al} | {len(set(net))} | {int(empty.sum())} | "
                      f"{st['MT-FW+SC']['recall']:.3f} / {st['QM-3C+SC']['recall']:.3f} | "
                      f"{fmt_ci(d0, *np.percentile(dstar, [2.5, 97.5]))} | {st['MT-FW+SC']['fr_empty']:.3f} "
                      f"[{np.percentile(fr_star, 2.5):.3f}, {np.percentile(fr_star, 97.5):.3f}] | "
                      f"{st['MT-FW+SC']['frac_empty_any']:.3f} |")
    # ------------------------------------------------------------------------------------ veredictos
    P("")
    P("## Veredictos")
    P("")
    P("| Regla | " + " | ".join(REGIMES) + " |")
    P("|:--|" + ":-:|" * len(REGIMES))
    for r_ in ("2D-a", "2D-a_holm", "2D-b", "2D-b_holm", "2D-c", "2D-c_strict", "2D-c_holm"):
        P(f"| {r_} | " + " | ".join("se cumple" if verdicts[g][r_] else "no se cumple" for g in REGIMES) + " |")
    P("")
    P("Casos a menos de 0.005 del umbral: " + ("; ".join(near) if near else "ninguno"))
    txt = "\n".join(L)
    print(txt)
    with open(a.out, "w") as fh:
        fh.write(txt + "\n")
    with open(a.out.replace(".txt", ".json"), "w") as fh:
        json.dump({"verdicts": verdicts, "thresholds": {"|".join(map(str, k_)): v for k_, v in thr_all.items()}},
                  fh, indent=1, default=str)


if __name__ == "__main__":
    main()
