"""analyze_2b.py - regla 2B de la Tarea 2B (PREREGISTRO_1j_1k.txt), por regimen de clutter, y tablas.

MT-FW+SC (registros de task2a_seed{s}_{regimen}.json) contra MT-FW-FREE+SC (task2b_seed{s}_{regimen}.json), en las
mismas ventanas (se comprueba ventana a ventana). Endpoint y estadistica como analyze_2a.py: recall a eventos falsos
realizados por ventana <= 0.1 y <= 0.25 sobre la curva FROC exacta, punto de operacion elegido sobre la propia
evaluacion; la regla se decide con estimaciones puntuales (diferencia >= 0.05 - 1e-9 a ambos topes en >= 2 de 3
celdas, en ambas semillas). IC 95%: bootstrap pareado por ventanas (2000 replicas, las mismas remuestras que
analyze_2a.py: default_rng([s, ci, ri])) y, como lectura adicional, bootstrap por red (2000 replicas, las 5 redes de
la celda con reemplazo, default_rng([s, ci, 77]) como verif_2a/cluster_boot_verif.py), re-eligiendo el punto de
operacion (y el maximo de los comparadores) en cada replica.
Exploratorio (no preregistrado): MT-FW-FREE+SC menos el mejor comparador de 2A (maximo por celda y tope de
GaMMA-STALTA(+SC), PyOcto-STALTA(+SC), QM, QM-3C, QM-3C+SC).

    python analyze_2b.py [--out analyze_2b_out.txt]
"""
import argparse
import json

import numpy as np

import analyze_2a as A

SEEDS = (2024, 7331)
FREE = "MT-FW-FREE+SC"
FULL = "MT-FW+SC"
CTX = ("MT-FW", "MT-FW-NODE+SC", "WAVE-PT@20")
NB = A.B


def boot_counts_windows(s, ci, ri, n):
    rng = np.random.default_rng([s, ci, ri])
    idx = rng.integers(0, n, size=(NB, n))
    C = np.zeros((NB, n))
    np.add.at(C, (np.repeat(np.arange(NB), n), idx.ravel()), 1.0)
    return C


def boot_counts_networks(s, ci, net):
    rng = np.random.default_rng([s, ci, 77])
    C = np.zeros((NB, len(net)))
    for b in range(NB):
        for j in rng.integers(0, 5, 5):
            C[b, net == j] += 1
    return C


def vc(HF, C, kk):
    """recall a cada tope para cada fila de C; la tasa de falsos se normaliza por el numero de ventanas de la replica
    (igual a n en el bootstrap por ventanas; variable en el bootstrap por red)."""
    H, F = HF
    Hs, Fs = C @ H, C @ F
    nn = C.sum(1)
    out = {}
    for cap in A.CAPS:
        ok = Fs / nn[:, None] <= float(cap) + 1e-12
        r = np.where(ok, Hs, -1.0).max(axis=1)
        out[cap] = np.where(r < 0, 0.0, r) / np.maximum(1.0, kk)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="analyze_2b_out.txt")
    a = ap.parse_args()
    L = []
    P = L.append
    verdicts, near, summary = {}, [], {}
    P("# Tarea 2B - ablacion del patron de radiacion (analyze_2b.py)")
    P("")
    P("Recall de eventos a eventos falsos realizados por ventana <= 0.1 / <= 0.25 (curva FROC exacta, punto de operacion "
      "elegido sobre la propia evaluacion). MT-FW+SC: registros de 2A. MT-FW-FREE+SC: misma llamada con G = Q = I. "
      f"IC 95%: bootstrap pareado por ventanas ({NB} replicas) y por red ({NB} replicas, lectura adicional). "
      "Diferencias: punto [IC 95%].")
    for ri, reg in enumerate(A.REGIMES):
        R2A = {s: {c["scenario"]: c for c in json.load(open(f"task2a_seed{s}_{reg}.json"))} for s in SEEDS}
        R2B = {s: {c["scenario"]: c for c in json.load(open(f"task2b_seed{s}_{reg}.json"))} for s in SEEDS}
        P("")
        P(f"## Regimen {reg} (clutter {R2A[SEEDS[0]][A.CELLS[0]]['clut_amp']})")
        P("")
        t1, t2, t3, t4 = [], [], [], []
        count = {s: 0 for s in SEEDS}
        detail = {s: [] for s in SEEDS}
        ind = {"win": {s: [] for s in SEEDS}, "net": {s: [] for s in SEEDS}}
        for s in SEEDS:
            for ci, cell in enumerate(A.CELLS):
                c2a, c2b = R2A[s][cell], R2B[s][cell]
                assert c2a["k_eval"] == c2b["k_eval"] and c2a["net_of_window"] == c2b["net_of_window"]
                assert [w["truth"] for w in c2a["windows_meta"]] == c2b["truth"]
                assert [w["taus"] for w in c2a["windows_meta"]] == c2b["taus"]
                assert c2a["methods"][FULL]["records"] == c2b["methods"][FULL]["records"]
                assert all(all(ct[key] for key in ("STA", "statics", "cat_drawn", "cat_src", "cat_tau", "corr"))
                           and ct["mtfw_sc_exact"] == 2 and ct["free_cert_eq_dense"] == 2 for ct in c2b["control"])
                k = np.asarray(c2a["k_eval"], float)
                n = len(k)
                net = np.asarray(c2a["net_of_window"])
                M = {m: A.mats(c2a["methods"][m]) for m in (FULL,) + CTX + A.SOTA}
                M[FREE] = A.mats(c2b["methods"][FREE])
                pt = {m: {cap: A.point(M[m], k, float(cap)) for cap in A.CAPS} for m in M}
                best = {cap: max(A.SOTA, key=lambda m: pt[m][cap]) for cap in A.CAPS}
                bpt = {cap: pt[best[cap]][cap] for cap in A.CAPS}
                d = {cap: pt[FULL][cap] - pt[FREE][cap] for cap in A.CAPS}
                ok = all(A.ge(d[cap]) for cap in A.CAPS)
                count[s] += ok
                detail[s].append((cell, ok, d))
                for cap in A.CAPS:
                    if abs(d[cap] - A.THR) < 0.005:
                        near.append(f"{reg} {s} {cell} tope {cap}: {d[cap]:+.4f}")
                summary[f"{reg}|{s}|{cell}"] = {"points": pt, "best_sota": best, "d_full_minus_free": d, "rule_2B_cell": ok}
                t1.append(f"| {A.LAB[cell]} ({s}) | {int(k.sum())} | "
                          + " | ".join(f"{pt[m]['0.1']:.2f}/{pt[m]['0.25']:.2f}" for m in (FULL, FREE) + CTX)
                          + f" | {bpt['0.1']:.2f}/{bpt['0.25']:.2f} ({best['0.1']} / {best['0.25']}) |")
                row2, row3, row4 = [], [], []
                for kind, C in (("win", boot_counts_windows(s, ci, ri, n)), ("net", boot_counts_networks(s, ci, net))):
                    kk = C @ k
                    bv = {m: vc(M[m], C, kk) for m in M}
                    ind[kind][s].append(np.all([A.ge(bv[FULL][cap] - bv[FREE][cap]) for cap in A.CAPS], axis=0))
                    for cap in A.CAPS:
                        dd = bv[FULL][cap] - bv[FREE][cap]
                        lo, hi = np.percentile(dd, [2.5, 97.5])
                        cell_txt = A.fmt_ci(d[cap], lo, hi)
                        if kind == "net":
                            cell_txt += f" P(<0.05)={np.mean(dd < A.THR - A.EPS):.2f}"
                        row2.append(cell_txt)
                        mx = np.max(np.stack([bv[m][cap] for m in A.SOTA]), axis=0)
                        de = bv[FREE][cap] - mx
                        lo, hi = np.percentile(de, [2.5, 97.5])
                        row3.append(A.fmt_ci(pt[FREE][cap] - bpt[cap], lo, hi))
                        if kind == "win":
                            dp = bv[FREE][cap] - bv["WAVE-PT@20"][cap]
                            lo, hi = np.percentile(dp, [2.5, 97.5])
                            row4.append(A.fmt_ci(pt[FREE][cap] - pt["WAVE-PT@20"][cap], lo, hi))
                            dn = bv[FULL][cap] - bv["MT-FW-NODE+SC"][cap]
                            lo, hi = np.percentile(dn, [2.5, 97.5])
                            row4.append(A.fmt_ci(pt[FULL][cap] - pt["MT-FW-NODE+SC"][cap], lo, hi))
                t2.append(f"| {A.LAB[cell]} ({s}) | " + " | ".join(row2) + f" | {'si' if ok else 'no'} |")
                t3.append(f"| {A.LAB[cell]} ({s}) | " + " | ".join(row3) + " |")
                t4.append(f"| {A.LAB[cell]} ({s}) | {row4[0]} | {row4[2]} | {row4[1]} | {row4[3]} |")
        P(f"### Tabla 1 ({reg}). Recall a <= 0.1 / <= 0.25")
        P("")
        P(f"| Celda (semilla) | eventos | {FULL} | {FREE} | " + " | ".join(CTX) + " | mejor comparador 2A (0.1 / 0.25) |")
        P("|:--|--:|" + ":-:|" * (2 + len(CTX)) + ":--|")
        L.extend(t1)
        P("")
        P(f"### Tabla 2 ({reg}). Regla 2B: {FULL} - {FREE}, punto [IC 95%]")
        P("")
        P("| Celda (semilla) | ventanas, tope 0.1 | ventanas, tope 0.25 | redes, tope 0.1 | redes, tope 0.25 | celda cumple |")
        P("|:--|:-:|:-:|:-:|:-:|:-:|")
        L.extend(t2)
        P("")
        P(f"### Tabla 3 ({reg}, exploratoria). {FREE} - mejor comparador de 2A (maximo re-elegido por replica)")
        P("")
        P("| Celda (semilla) | ventanas, tope 0.1 | ventanas, tope 0.25 | redes, tope 0.1 | redes, tope 0.25 |")
        P("|:--|:-:|:-:|:-:|:-:|")
        L.extend(t3)
        P("")
        P(f"### Tabla 4 ({reg}, exploratoria). {FREE} - WAVE-PT@20 y, como referencia, 2A-b ({FULL} - MT-FW-NODE+SC); "
          "bootstrap por ventanas")
        P("")
        P(f"| Celda (semilla) | FREE - WAVE-PT@20, 0.1 | FREE - WAVE-PT@20, 0.25 | {FULL} - NODE+SC, 0.1 | {FULL} - NODE+SC, 0.25 |")
        P("|:--|:-:|:-:|:-:|:-:|")
        L.extend(t4)
        ok = all(count[s] >= 2 for s in SEEDS)
        verdicts[reg] = ok
        P("")
        P(f"### Regla 2B ({reg})")
        P("")
        P("- 2B: " + "; ".join(f"semilla {s}: {count[s]}/3 celdas ("
                               + ", ".join(f"{A.LAB[c_]} {'si' if o else 'no'} " + "/".join(f"{x:+.2f}" for x in v.values())
                                           for c_, o, v in detail[s]) + ")" for s in SEEDS)
          + f" -> {'SE CUMPLE' if ok else 'NO SE CUMPLE'}")
        rr = []
        for kind, lab in (("win", "ventanas"), ("net", "redes")):
            hold = np.all([np.sum(ind[kind][s], axis=0) >= 2 for s in SEEDS], axis=0)
            rr.append(f"bootstrap por {lab} {np.mean(hold):.2f}")
        P("- Exploratorio: fraccion de replicas en que la regla 2B se cumple (remuestras independientes por celda): "
          + "; ".join(rr))
    P("")
    P("## Veredicto")
    P("")
    P("| Regla | " + " | ".join(A.REGIMES) + " |")
    P("|:--|" + ":-:|" * len(A.REGIMES))
    P("| 2B | " + " | ".join("se cumple" if verdicts[g] else "no se cumple" for g in A.REGIMES) + " |")
    P("")
    for g in A.REGIMES:
        if verdicts[g]:
            P(f"Interpretacion preregistrada ({g}): 2B se cumple: dentro del mismo solver, el patron de radiacion aporta "
              ">= 0.05 de recall a ambos topes sobre la energia libre de signo.")
        else:
            P(f"Interpretacion preregistrada ({g}): 2B falla: en este regimen la ganancia de 2A se atribuye a la migracion "
              "sin picks con deflacion y no al patron de radiacion, y la propuesta lo declara.")
    P("")
    P("Casos a menos de 0.005 del umbral 0.05: " + ("; ".join(near) if near else "ninguno"))
    txt = "\n".join(L)
    print(txt)
    with open(a.out, "w") as fh:
        fh.write(txt + "\n")
    with open(a.out.replace(".txt", ".json"), "w") as fh:
        json.dump({"verdicts": verdicts, "cells": summary}, fh, indent=1)


if __name__ == "__main__":
    main()
