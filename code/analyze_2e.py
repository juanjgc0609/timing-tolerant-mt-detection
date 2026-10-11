"""analyze_2e.py - reglas 2E-a (con la lectura estricta de Holm) y 2E-b de la Tarea 2E (PREREGISTRO_1j_1k.txt), por
regimen de clutter; tablas descriptivas, errores de localizacion y diagnosticos de convergencia.  Escrito antes de ver
cualquier resultado de 2E en las semillas 2024 y 7331.

Datos por semilla, regimen y celda: 30 redes (redes 0-4 de 2A con 20 ventanas, redes 5-29 de 2C parte (ii) con 12),
task2e_parts/ev_*.  TS-SFW+SC y BOMP-PK-SLIDE+SC de 2E; BOMP-PK+SC: registros de 2D (task2d_parts/p1_*) en las redes 0-4
(el recalculo de 2E se compara con ellos) y recalculados en 2E (codigo de 2D sin cambios) en las redes 5-29; MT-FW+SC:
registros guardados de 2A (redes 0-4) y 2C (redes 5-29).  Se comprueba ventana a ventana que k, verdad y origenes
coinciden con los guardados.
Endpoint (como 2C y 2D): recall de eventos a eventos falsos realizados por ventana <= 0.1 y <= 0.25, curva FROC exacta
(froc_exact.exact_HF), punto de operacion elegido por metodo sobre la propia evaluacion; una diferencia cuenta si es
>= 0.05 - 1e-9.  IC 95%: bootstrap jerarquico (2000 replicas; en cada una las 30 redes de la celda con reemplazo y,
dentro de cada red elegida, sus ventanas con reemplazo; default_rng([s, ci, ri, 2801])), las mismas remuestras para
todos los metodos de la celda, punto de operacion re-elegido en cada replica (analyze_2c.val).
2E-a (primaria), por regimen: TS-SFW+SC - BOMP-PK+SC a tope 0.25 con estimacion puntual >= 0.05 y cota inferior del IC
95 jerarquico > 0 en >= 2 de 3 celdas en ambas semillas.  Lectura estricta adicional: Holm sobre las 6 pruebas celda x
semilla del regimen, p unilateral p = (1 + #{d* <= 0}) / (B + 1), nivel 0.025; una celda cuenta si Holm rechaza y la
estimacion puntual es >= 0.05; >= 2 de 3 celdas en ambas semillas.
2E-b (secundaria), por regimen: TS-SFW+SC - BOMP-PK-SLIDE+SC >= 0.05 a tope 0.25 (estimacion puntual) en >= 2 de 3
celdas en ambas semillas.
Descriptivo: TS-SFW+SC - MT-FW+SC y BOMP-PK-SLIDE+SC - BOMP-PK+SC con IC a ambos topes; error epicentral (y 3D en los
metodos continuos) mediano de los aciertos por metodo en el punto de operacion a tope 0.25 (un acierto por evento
verdadero: el primer atomo del registro que lo acierta); diagnosticos de convergencia (brecha relativa del paso convexo,
KKT, alternancias, iteraciones de L-BFGS-B, fraccion de desplazamientos en la cota, atomos con amplitud nula).
Exploratorio (declarado): fraccion de replicas en que se cumplirian las reglas.

    python analyze_2e.py [--out analyze_2e_out.txt]
"""
import argparse
import json
import os

import numpy as np

from analyze_2c import HF, holm, point, val
from froc_exact import exact_HF

HERE = os.path.dirname(os.path.abspath(__file__))
B = 2000
CAPS = ("0.1", "0.25")
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
LAB = {"ladder_4": "4 st.", "ladder_6": "6 st.", "unilateral_c40": "one-sided 8"}
REGIMES = ("thin", "strong")
SEEDS = (2024, 7331)
EPS = 1e-9
THR = 0.05
ALPHA_HOLM = 0.025
N_NETS_2A, N_NETS = 5, 30
METHODS = ("TS-SFW+SC", "BOMP-PK-SLIDE+SC", "BOMP-PK+SC", "MT-FW+SC")
CONT = ("TS-SFW+SC", "BOMP-PK-SLIDE+SC")
TAG = 2801


def ge(d):
    return d >= THR - EPS


def fmt_ci(p, lo, hi):
    return f"{p:+.3f} [{lo:+.3f}, {hi:+.3f}]"


def load(path):
    with open(path) as fh:
        return json.load(fh)


def counts_hier(s, ci, ri, net):
    rng = np.random.default_rng([s, ci, ri, TAG])
    ids = sorted(set(net))
    pos = {j: np.where(net == j)[0] for j in ids}
    C = np.zeros((B, len(net)))
    for b in range(B):
        for jj in rng.integers(0, len(ids), len(ids)):
            p = pos[ids[jj]]
            np.add.at(C[b], p[rng.integers(0, len(p), len(p))], 1.0)
    return C


def cell_data(s, reg, cell):
    W = {"k": [], "net": [], "rec": {m: [] for m in METHODS}, "loc": {m: [] for m in METHODS}, "truth": [],
         "diag": {m: [] for m in CONT}}
    ctl = {"nets": 0, "windows": 0, "mismatch": 0, "controls": [], "kernels": set(), "kernel_required": [],
           "bomp_eq": 0, "bomp_n": 0, "mtfw_loc_eq": 0, "fvals": set()}
    for j in range(N_NETS):
        e = load(os.path.join(HERE, "task2e_parts", f"ev_seed{s}_{reg}_{cell}_net{j}.json"))
        ctl["nets"] += 1
        ctl["controls"].append(e["control"])
        ctl["kernels"].add(tuple(tuple(x) for x in e["kernel"]["blas"]))
        ctl["kernel_required"].append(e["kernel_required"])
        ctl["fvals"].add(tuple(e["fvals"]))
        if j < N_NETS_2A:
            ref = load(os.path.join(HERE, "task2a_parts", f"seed{s}_{reg}_{cell}_net{j}.json"))
            d2 = load(os.path.join(HERE, "task2d_parts", f"p1_seed{s}_{reg}_{cell}_net{j}.json"))
        else:
            ref = load(os.path.join(HERE, "task2c_parts", f"p2_seed{s}_{reg}_{cell}_net{j}.json"))
            d2 = None
        assert e["n_win"] == ref["n_win"] == len(e["windows"])
        for w, (we, wr) in enumerate(zip(e["windows"], ref["windows"])):
            same = all(we[x] == wr[x] for x in ("k", "truth", "taus")) and \
                we["records"]["MT-FW+SC"] == wr["records"]["MT-FW+SC"]
            if d2 is not None:
                same &= all(we[x] == d2["windows"][w][x] for x in ("k", "truth", "taus"))
            ctl["windows"] += 1
            ctl["mismatch"] += not same
            W["k"].append(we["k"])
            W["net"].append(j)
            W["truth"].append(we["truth"])
            W["rec"]["TS-SFW+SC"].append(we["records"]["TS-SFW+SC"])
            W["rec"]["BOMP-PK-SLIDE+SC"].append(we["records"]["BOMP-PK-SLIDE+SC"])
            W["rec"]["MT-FW+SC"].append(wr["records"]["MT-FW+SC"])
            if d2 is not None:
                rb = d2["windows"][w]["records"]["BOMP-PK+SC"]
                eq = rb == we["records"]["BOMP-PK+SC"]
                ctl["bomp_eq"] += eq
                ctl["bomp_n"] += 1
            else:
                rb = we["records"]["BOMP-PK+SC"]
                eq = True
            W["rec"]["BOMP-PK+SC"].append(rb)
            for m in CONT:
                W["loc"][m].append(we["loc"][m])
                W["diag"][m].append(we["diag"][m])
            W["loc"]["BOMP-PK+SC"].append(we["loc"]["BOMP-PK+SC"] if eq else None)
            leq = we["mtfw_loc_eq_stored"]
            ctl["mtfw_loc_eq"] += leq
            W["loc"]["MT-FW+SC"].append(we["loc"]["MT-FW+SC(loc)"] if leq else None)
    W["k"] = np.asarray(W["k"], float)
    W["net"] = np.asarray(W["net"])
    return W, ctl


def op_threshold(recs, k, cap):
    """umbral del punto de operacion (el mayor umbral con el mayor numero de aciertos y falsos por ventana <= cap)."""
    thr, H, F = exact_HF([[tuple(x) for x in r] for r in recs])
    Hs, Fs = H.sum(0), F.sum(0)
    ok = Fs / len(recs) <= float(cap) + 1e-12
    if not ok.any():
        return np.inf
    best = Hs[ok].max()
    cand = np.where(ok & (Hs == best))[0]
    return float(thr[cand.min()])


def loc_errors(W, m, cap="0.25"):
    """errores de los aciertos (un acierto por evento: el primer atomo del registro que lo acierta) con puntaje >= el
    umbral del punto de operacion."""
    t = op_threshold(W["rec"][m], W["k"], cap)
    eh, e3, dz = [], [], []
    for recs, loc, tr in zip(W["rec"][m], W["loc"][m], W["truth"]):
        if loc is None:
            continue
        seen = set()
        for (sc, h), p in zip(recs, loc):
            if sc >= t and h >= 0 and h not in seen:
                seen.add(h)
                q = np.asarray(tr[h], float)
                p = np.asarray(p[:3], float)
                eh.append(float(np.linalg.norm(p[:2] - q[:2])))
                e3.append(float(np.linalg.norm(p - q)))
                dz.append(float(p[2] - q[2]))
    return np.array(eh), np.array(e3), np.array(dz)


def diag_summary(dl):
    it = [x for d in dl for x in d["iters"]]
    alt = np.array([x[3] for x in it])
    nit = np.array([n for x in it for n in x[4]])
    gaps = np.array([x[6] for x in it])
    kkt = np.array([x[8] for x in it])
    ref = np.array([x[10] for x in it if x[10] is not None])
    shb = np.array([d["frac_shift_at_bound"] for d in dl])
    zero = np.array([mn == 0.0 for d in dl for mn in d["m_norm"]])
    lam = np.array([d["lam"] for d in dl])
    s1 = np.array([d["iters"][0][0] for d in dl if d["iters"]])
    ratio = lam[:len(s1)] / np.sqrt(2 * s1) if len(s1) else np.zeros(0)
    return {"n_windows": len(dl), "n_iters": len(it), "alt_mean": float(alt.mean()), "alt_frac5": float(np.mean(alt >= 5)),
            "lbfgs_steps": int(len(nit)), "lbfgs_frac50": float(np.mean(nit >= 50)) if len(nit) else 0.0,
            "lbfgs_nit_med": float(np.median(nit)) if len(nit) else 0.0,
            "gap_max": float(gaps.max()), "gap_med": float(np.median(gaps)), "kkt_max": float(kkt.max()),
            "refine_frac50": float(np.mean(ref >= 50)) if len(ref) else 0.0,
            "shift_at_bound_mean": float(shb.mean()), "zero_amp_frac": float(zero.mean()) if len(zero) else 0.0,
            "lam_over_c1_med": float(np.median(ratio)) if len(ratio) else 0.0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "analyze_2e_out.txt"))
    a = ap.parse_args()
    L = []
    P = L.append
    verdicts = {}
    near = []
    P("# Tarea 2E - BLASSO de grupos con desplazamientos fuera de la rejilla (TS-SFW+SC) contra la busqueda por bloques "
      "con la misma tolerancia temporal (analyze_2e.py)")
    P("")
    P("Recall de eventos a eventos falsos realizados por ventana <= 0.1 / <= 0.25 (curva FROC exacta, punto de operacion "
      "elegido por metodo sobre la propia evaluacion), 400 ventanas (5 redes x 20 + 25 redes x 12) por celda, semilla y "
      "regimen. Diferencias: punto [IC 95% H] (bootstrap jerarquico: redes y ventanas; 2000 replicas; punto de operacion "
      "re-elegido en cada replica); p unilateral H = (1 + #{d* <= 0}) / 2001.")
    data = {(reg, s, cell): cell_data(s, reg, cell) for reg in REGIMES for s in SEEDS for cell in CELLS}
    P("")
    P("## Controles")
    P("")
    fall = set()
    for reg in REGIMES:
        for s in SEEDS:
            for cell in CELLS:
                W, ctl = data[(reg, s, cell)]
                fall |= ctl["fvals"]
                ok = all(all(v for k_, v in c.items() if k_ in ("STA", "statics", "cat_drawn", "cat_src", "cat_tau", "corr"))
                         and c["mtfw_sc_exact"] == 2 for c in ctl["controls"])
                P(f"- {reg} {s} {LAB[cell]}: {ctl['nets']} redes, {ctl['windows']} ventanas; red, catalogo, correcciones y "
                  f"MT-FW+SC exacto en las ventanas 0-1 {'correcto' if ok else 'FALLA'}; ventanas distintas de lo guardado "
                  f"{ctl['mismatch']}; BOMP-PK+SC recalculado = 2D en {ctl['bomp_eq']}/{ctl['bomp_n']} ventanas de las redes "
                  f"0-4; MT-FW+SC con posiciones = guardado en {ctl['mtfw_loc_eq']}/{ctl['windows']}; nucleo requerido "
                  f"{sorted(set(ctl['kernel_required']))}, cargado {sorted(ctl['kernels'])}")
                assert ok and ctl["mismatch"] == 0
    P(f"- f de TS-SFW+SC en todas las partes: {sorted(fall)}")
    tests = {reg: {} for reg in REGIMES}
    pts = {reg: {} for reg in REGIMES}
    locs = {}
    diags = {}
    for ri, reg in enumerate(REGIMES):
        P("")
        P(f"## Regimen {reg}")
        P("")
        P(f"### Tabla 1 ({reg}). Recall a <= 0.1 / <= 0.25")
        P("")
        P("| Celda (semilla) | eventos | " + " | ".join(METHODS) + " |")
        P("|:--|--:|" + ":-:|" * len(METHODS))
        rowsA, rowsB, rowsD = [], [], []
        cnt = {"2E-a": {s: 0 for s in SEEDS}, "2E-b": {s: 0 for s in SEEDS}}
        det = {"2E-a": {s: [] for s in SEEDS}, "2E-b": {s: [] for s in SEEDS}}
        rep = {"2E-a": {s: [] for s in SEEDS}, "2E-b": {s: [] for s in SEEDS}}
        for s in SEEDS:
            for ci, cell in enumerate(CELLS):
                W, _ctl = data[(reg, s, cell)]
                k, net = W["k"], W["net"]
                HFs = {m: HF(W["rec"][m]) for m in METHODS}
                pt = {m: point(HFs[m], k) for m in METHODS}
                P(f"| {LAB[cell]} ({s}) | {int(k.sum())} | " + " | ".join(
                    f"{pt[m]['0.1']:.3f}/{pt[m]['0.25']:.3f}" for m in METHODS) + " |")
                Ch = counts_hier(s, ci, ri, net)
                bh = {m: val(HFs[m], Ch, k) for m in METHODS}
                ca, cb = [], []
                for cap in CAPS:
                    p0 = pt["TS-SFW+SC"][cap] - pt["BOMP-PK+SC"][cap]
                    dh = bh["TS-SFW+SC"][cap] - bh["BOMP-PK+SC"][cap]
                    lo, hi = np.percentile(dh, [2.5, 97.5])
                    pv = (1.0 + np.sum(dh <= 0)) / (B + 1.0)
                    pts[reg][("2E-a", s, cell, cap)] = (p0, lo, hi, pv)
                    ca.append(f"{fmt_ci(p0, lo, hi)} p={pv:.4f}")
                    if cap == "0.25":
                        tests[reg][(s, cell)] = pv
                        oka = ge(p0) and lo > 0
                        cnt["2E-a"][s] += oka
                        det["2E-a"][s].append((cell, oka, p0, lo))
                        rep["2E-a"][s].append((dh >= THR - EPS) & (np.percentile(dh, 2.5) > 0))
                        if abs(p0 - THR) < 0.005 or abs(lo) < 0.005:
                            near.append(f"2E-a {reg} {s} {cell}: punto {p0:+.4f}, cota {lo:+.4f}")
                    p0 = pt["TS-SFW+SC"][cap] - pt["BOMP-PK-SLIDE+SC"][cap]
                    dh = bh["TS-SFW+SC"][cap] - bh["BOMP-PK-SLIDE+SC"][cap]
                    lo, hi = np.percentile(dh, [2.5, 97.5])
                    pv = (1.0 + np.sum(dh <= 0)) / (B + 1.0)
                    pts[reg][("2E-b", s, cell, cap)] = (p0, lo, hi, pv)
                    cb.append(f"{fmt_ci(p0, lo, hi)} p={pv:.4f}")
                    if cap == "0.25":
                        okb = ge(p0)
                        cnt["2E-b"][s] += okb
                        det["2E-b"][s].append((cell, okb, p0, lo))
                        rep["2E-b"][s].append(dh >= THR - EPS)
                        if abs(p0 - THR) < 0.005:
                            near.append(f"2E-b {reg} {s} {cell}: {p0:+.4f}")
                rowsA.append(f"| {LAB[cell]} ({s}) | {ca[0]} | {ca[1]} | "
                             f"{'si' if det['2E-a'][s][-1][1] else 'no'} |")
                rowsB.append(f"| {LAB[cell]} ({s}) | {cb[0]} | {cb[1]} | {'si' if det['2E-b'][s][-1][1] else 'no'} |")
                dd = []
                for m1, m2 in (("TS-SFW+SC", "MT-FW+SC"), ("BOMP-PK-SLIDE+SC", "BOMP-PK+SC"),
                               ("BOMP-PK-SLIDE+SC", "MT-FW+SC"), ("BOMP-PK+SC", "MT-FW+SC")):
                    dd.append(" / ".join(fmt_ci(pt[m1][cap] - pt[m2][cap],
                                                *np.percentile(bh[m1][cap] - bh[m2][cap], [2.5, 97.5])) for cap in CAPS))
                rowsD.append(f"| {LAB[cell]} ({s}) | " + " | ".join(dd) + " |")
                for m in METHODS:
                    locs[(reg, s, cell, m)] = loc_errors(W, m)
                for m in CONT:
                    diags[(reg, s, cell, m)] = diag_summary(W["diag"][m])
        P("")
        P(f"### Tabla 2 ({reg}). 2E-a: TS-SFW+SC - BOMP-PK+SC [IC 95% H] p unilateral H")
        P("")
        P("| Celda (semilla) | a 0.1 | a 0.25 | celda cuenta (punto >= 0.05 y cota > 0 a 0.25) |")
        P("|:--|:-:|:-:|:-:|")
        L.extend(rowsA)
        P("")
        P(f"### Tabla 3 ({reg}). 2E-b: TS-SFW+SC - BOMP-PK-SLIDE+SC [IC 95% H] p unilateral H")
        P("")
        P("| Celda (semilla) | a 0.1 | a 0.25 | celda cuenta (punto >= 0.05 a 0.25) |")
        P("|:--|:-:|:-:|:-:|")
        L.extend(rowsB)
        P("")
        P(f"### Tabla 4 ({reg}, descriptiva). Diferencias a 0.1 / 0.25 [IC 95% H]")
        P("")
        P("| Celda (semilla) | TS-SFW+SC - MT-FW+SC | BOMP-PK-SLIDE+SC - BOMP-PK+SC | BOMP-PK-SLIDE+SC - MT-FW+SC | "
          "BOMP-PK+SC - MT-FW+SC |")
        P("|:--|:-:|:-:|:-:|:-:|")
        L.extend(rowsD)
        verdicts[reg] = {}
        for rule in ("2E-a", "2E-b"):
            okr = all(cnt[rule][s] >= 2 for s in SEEDS)
            joint = np.all([np.sum(rep[rule][s], axis=0) >= 2 for s in SEEDS], axis=0)
            txt = "punto >= 0.05 y cota inferior > 0" if rule == "2E-a" else "punto >= 0.05"
            P("")
            P(f"- Regla {rule} ({reg}; {txt}, a 0.25, >= 2 de 3 celdas en ambas semillas): " + "; ".join(
                f"semilla {s}: {cnt[rule][s]}/3 (" + ", ".join(
                    f"{LAB[c_]} {'si' if o else 'no'} {x:+.3f} [cota {lo:+.3f}]" for c_, o, x, lo in det[rule][s]) + ")"
                for s in SEEDS) + f" -> {'SE CUMPLE' if okr else 'NO SE CUMPLE'}"
                + (f" (exploratorio: puntos >= 0.05 en >= 2 de 3 celdas por semilla en {100 * joint.mean():.0f}% de las "
                   f"replicas H)" if rule == "2E-b" else ""))
            verdicts[reg][rule] = okr
            verdicts[reg][rule + "_counts"] = cnt[rule]
            verdicts[reg][rule + "_points"] = {s: [(c_, x, lo) for c_, _o, x, lo in det[rule][s]] for s in SEEDS}
        rej, adj = holm(tests[reg], ALPHA_HOLM)
        cntH = {s: 0 for s in SEEDS}
        P("")
        P(f"### Tabla 5 ({reg}). 2E-a, lectura estricta: Holm sobre las 6 pruebas celda x semilla (unilateral, nivel "
          f"{ALPHA_HOLM})")
        P("")
        P("| Celda (semilla) | punto a 0.25 | p H | p Holm | rechazo | cuenta (rechazo y punto >= 0.05) |")
        P("|:--|:-:|:-:|:-:|:-:|:-:|")
        for s in SEEDS:
            for cell in CELLS:
                k_ = (s, cell)
                p0 = pts[reg][("2E-a", s, cell, "0.25")][0]
                c_ok = rej[k_] and ge(p0)
                cntH[s] += c_ok
                P(f"| {LAB[cell]} ({s}) | {p0:+.3f} | {tests[reg][k_]:.4f} | {adj[k_]:.4f} | {'si' if rej[k_] else 'no'} "
                  f"| {'si' if c_ok else 'no'} |")
        okH = all(cntH[s] >= 2 for s in SEEDS)
        P("")
        P(f"- 2E-a ({reg}) lectura estricta con Holm: " + "; ".join(f"semilla {s}: {cntH[s]}/3" for s in SEEDS)
          + f" -> {'SE CUMPLE' if okH else 'NO SE CUMPLE'}")
        verdicts[reg]["2E-a_holm"] = okH
        verdicts[reg]["holm"] = {f"{s}|{c_}": [bool(rej[(s, c_)]), float(adj[(s, c_)]), float(tests[reg][(s, c_)])]
                                 for s, c_ in tests[reg]}
        # ------------------------------------------------------------------ localizacion
        P("")
        P(f"### Tabla 6 ({reg}, descriptiva). Error de localizacion de los aciertos en el punto de operacion a 0.25: "
          "epicentral mediano (km) [n aciertos]; en los metodos continuos tambien 3D mediano y sesgo mediano de profundidad")
        P("")
        P("| Celda (semilla) | " + " | ".join(METHODS) + " |")
        P("|:--|" + ":-:|" * len(METHODS))
        for s in SEEDS:
            for cell in CELLS:
                cells = []
                for m in METHODS:
                    eh, e3, dz = locs[(reg, s, cell, m)]
                    if len(eh) == 0:
                        cells.append("-")
                    elif m in CONT:
                        cells.append(f"{np.median(eh):.2f} [{len(eh)}]; 3D {np.median(e3):.2f}; dz {np.median(dz):+.2f}")
                    else:
                        cells.append(f"{np.median(eh):.2f} [{len(eh)}]")
                P(f"| {LAB[cell]} ({s}) | " + " | ".join(cells) + " |")
        allm = {m: np.concatenate([locs[(reg, s, c_, m)][0] for s in SEEDS for c_ in CELLS]) for m in METHODS}
        P("")
        P(f"- Todas las celdas ({reg}): error epicentral mediano " + "; ".join(
            f"{m} {np.median(v):.2f} km (P90 {np.percentile(v, 90):.2f}, n {len(v)})" for m, v in allm.items() if len(v)))
        # ------------------------------------------------------------------ convergencia
        P("")
        P(f"### Tabla 7 ({reg}). Diagnosticos de convergencia de los metodos continuos")
        P("")
        P("| Metodo | Celda (semilla) | brecha rel. max / mediana | KKT max | alternancias media (frac. = 5) | L-BFGS-B "
          "deslizamiento: pasos, frac. con 50 it., mediana it. | refinamiento frac. con 50 it. | frac. desplazamientos en la "
          "cota | frac. atomos con m = 0 | lambda / sqrt(2 s_1) mediana |")
        P("|:--|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|")
        for m in CONT:
            for s in SEEDS:
                for cell in CELLS:
                    d = diags[(reg, s, cell, m)]
                    P(f"| {m} | {LAB[cell]} ({s}) | {d['gap_max']:.1e} / {d['gap_med']:.1e} | {d['kkt_max']:.1e} | "
                      f"{d['alt_mean']:.2f} ({d['alt_frac5']:.2f}) | {d['lbfgs_steps']}, {d['lbfgs_frac50']:.2f}, "
                      f"{d['lbfgs_nit_med']:.0f} | {d['refine_frac50']:.2f} | {d['shift_at_bound_mean']:.3f} | "
                      f"{d['zero_amp_frac']:.3f} | {d['lam_over_c1_med']:.3f} |")
    P("")
    P("## Veredictos")
    P("")
    P("| Regla | " + " | ".join(REGIMES) + " |")
    P("|:--|" + ":-:|" * len(REGIMES))
    for r_ in ("2E-a", "2E-a_holm", "2E-b"):
        P(f"| {r_} | " + " | ".join("se cumple" if verdicts[g][r_] else "no se cumple" for g in REGIMES) + " |")
    P("")
    P("Casos a menos de 0.005 del umbral o de 0: " + ("; ".join(near) if near else "ninguno"))
    a_any = any(verdicts[g]["2E-a"] for g in REGIMES)
    P("")
    P("Lectura preregistrada: " + (
        "2E-a falla en ambos regimenes: G1-M queda refutada en este diseno sintetico." if not a_any else
        "; ".join(f"{g}: " + ("2E-a y 2E-b se cumplen: G1-M apoyada en sintetico" if verdicts[g]["2E-a"] and verdicts[g]["2E-b"]
                              else "2E-a se cumple y 2E-b no: la ganancia se atribuye al deslizamiento fuera de la rejilla"
                              if verdicts[g]["2E-a"] else "2E-a no se cumple") for g in REGIMES)))
    txt = "\n".join(L)
    print(txt)
    with open(a.out, "w") as fh:
        fh.write(txt + "\n")
    with open(a.out.replace(".txt", ".json"), "w") as fh:
        json.dump({"verdicts": verdicts,
                   "points": {reg: {"|".join(map(str, k_)): list(v) for k_, v in pts[reg].items()} for reg in REGIMES},
                   "diagnostics": {"|".join(map(str, k_)): v for k_, v in diags.items()},
                   "loc_median_epi": {"|".join(map(str, k_)): (float(np.median(v[0])) if len(v[0]) else None)
                                      for k_, v in locs.items()}}, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
