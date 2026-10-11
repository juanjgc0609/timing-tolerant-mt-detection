"""
task0_curves.py — diagnóstico de dónde el canal de polarización ayuda y dónde estorba.

Los resultados v3 muestran un desdoblamiento: con el estimador (20c) el modelo restringido gana
a la ablación en el PUNTO CALIBRADO (6 de 8 celdas) y pierde en la MÉTRICA DE CURVA
falsos@recall=0.98 (5 de 8). Eso significa que el canal tiene un régimen de validez y que un
punto único lo oculta.

Este script no intenta mejorar nada. Traza la curva completa recall-vs-falsos de ABL y PHYS en
las dos celdas donde más discrepan y localiza el cruce, que es la cantidad que decide si el canal
tiene un régimen o simplemente no aporta.

    python task0_curves.py               # dos celdas, 30 realizaciones
    python task0_curves.py --reps 12     # rápido
"""
import argparse
import numpy as np

import task0_common as tc
from task0_common import (make_network, generate, field_own, operating_curve,
                          default_grid, DEFAULT_TAUS)

SRC = np.array([6.0, -8.0, 9.0])
TAU = [6.0]

CELLS = [
    # etiqueta, geometría, n_est, clutter, señuelo, soporte
    ("dispersa, clutter 40", "spread", 8, 40, None, 0),
    ("senuelo 12 km, 12/16", "spread", 8, 8, SRC + np.array([12.0, 0.0, 0.0]), 12),
]


def crossing(curve_a, curve_b, n=41):
    """
    Recall donde las dos curvas se cruzan en tasa de falsos. Interpola ambas sobre una rejilla
    común de recall y busca cambios de signo de (falsos_A - falsos_B).
    Devuelve lista de recalls de cruce y la diferencia media en cada tramo.
    """
    if not curve_a or not curve_b:
        return [], []
    ra = np.array([r for (_t, r, _f) in curve_a])
    fa = np.array([f for (_t, _r, f) in curve_a])
    rb = np.array([r for (_t, r, _f) in curve_b])
    fb = np.array([f for (_t, _r, f) in curve_b])
    lo = max(ra.min(), rb.min())
    hi = min(ra.max(), rb.max())
    if hi <= lo:
        return [], []
    gr = np.linspace(lo, hi, n)
    ia = np.interp(gr, ra, fa)
    ib = np.interp(gr, rb, fb)
    d = ia - ib                      # >0 : ABL peor que PHYS
    xs = []
    for k in range(len(gr) - 1):
        if d[k] == 0.0:
            xs.append(float(gr[k]))
        elif d[k] * d[k + 1] < 0:
            w = abs(d[k]) / (abs(d[k]) + abs(d[k + 1]))
            xs.append(float(gr[k] + w * (gr[k + 1] - gr[k])))
    return xs, (gr, d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--seed", type=int, default=777)
    args = ap.parse_args()

    grid = default_grid()
    print(f"curvas completas recall-falsos | {args.reps} realizaciones | seed {args.seed}\n")

    for label, kind, nst, nclut, dec, sup in CELLS:
        STA = make_network(kind, nst, seed=args.seed)
        rng = np.random.default_rng(args.seed)
        data = [generate(STA, [SRC], TAU, n_clutter=nclut, rng=rng,
                         coherent_decoy=dec, decoy_support=sup)
                for _ in range(args.reps)]
        truth = [t for (_D, _p, t) in data]

        curves = {}
        for name, use in (("ABL", set()), ("PHYS", {"CIII", "CIV"})):
            fields = [field_own(STA, D, p, use, grid, DEFAULT_TAUS) for (D, p, _t) in data]
            curves[name] = operating_curve(fields, truth, grid, n_points=25)

        print(f"  {label}")
        print("   El recall satura en 1.00 en casi todo el rango, asi que el eje informativo")
        print("   es la FRACCION del maximo, no el recall.\n")
        print(f"   {'fraccion':>9} | {'ABL rec':>8} {'ABL fal':>8} | {'PHYS rec':>9} {'PHYS fal':>9} | {'ABL-PHYS':>9}")
        print("   " + "-" * 64)
        cross = []
        prev = None
        for (fa_, ra_, va_), (fb_, rb_, vb_) in zip(
                sorted(curves["ABL"], key=lambda x: -x[0]),
                sorted(curves["PHYS"], key=lambda x: -x[0])):
            d = va_ - vb_
            if fa_ in (1.02,):
                continue
            if True:
                mark = "PHYS mejor" if d > 0 else ("ABL mejor" if d < 0 else "igual")
                print(f"   {fa_:9.2f} | {ra_:8.2f} {va_:8.3f} | {rb_:9.2f} {vb_:9.3f} | "
                      f"{d:+9.3f}  {mark}")
            if prev is not None and prev * d < 0:
                cross.append(fa_)
            if d != 0:
                prev = d

        allv = [(a[2] - b[2]) for a, b in zip(sorted(curves["ABL"], key=lambda x: -x[0]),
                                              sorted(curves["PHYS"], key=lambda x: -x[0]))]
        pos = sum(1 for v in allv if v > 0)
        neg = sum(1 for v in allv if v < 0)
        print(f"\n   puntos donde PHYS es mejor: {pos} | donde ABL es mejor: {neg} "
              f"| empates: {len(allv) - pos - neg}")
        if cross:
            print(f"   cambio(s) de signo en fraccion = " + ", ".join(f"{c:.2f}" for c in cross))
            print("   -> el canal tiene un REGIMEN: ayuda de un lado del umbral y no del otro")
        elif pos > neg:
            print("   -> PHYS domina en todo el rango donde hay diferencia")
        elif neg > pos:
            print("   -> ABL domina en todo el rango donde hay diferencia")
        else:
            print("   -> sin diferencia apreciable")
        print()

    print("Lectura: un cruce significa que el canal aporta por encima de cierto recall y resta")
    print("por debajo, lo cual es un resultado sobre su regimen de validez. La ausencia de cruce")
    print("significa que la conclusion es uniforme y no depende del punto de operacion.")


if __name__ == "__main__":
    main()
