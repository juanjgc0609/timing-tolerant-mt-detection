"""
task1_decide.py — Etapa 3 de la Tarea 1: la cantidad que decide.

Calcula la fracción de picks NO ASOCIADOS cuya polarización es coherente con un evento cercano, y
la compara contra la curva de respuesta del Apéndice T. El criterio está fijado aquí y no se
ajusta después de ver el número.

    python task1_decide.py
    python task1_decide.py --coh-deg 20 --rect-min 0.6

Salida: t1_data/decision.json y un veredicto en pantalla.

CRITERIO PRE-REGISTRADO (Apéndice T de la propuesta)
    fraccion coherente < 0.20  -> los apendices H, J, M mantienen su validez cuantitativa
    0.20 <= fraccion <= 0.50   -> reescalar los margenes absolutos con la curva del Apendice T
    fraccion > 0.50            -> re-medir la jerarquia de canales; A2 probablemente no sobrevive
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

OUT = "t1_data"
CSV = os.path.join(OUT, "picks_polarization.csv")

# curva de respuesta del Apendice T: tasa espuria de la polarizacion medida en sintetico
APPENDIX_T = {"clutter aleatorio": 0.00, "senuelo coherente": 0.40}


def bootstrap_ci(x, stat=np.mean, n=4000, alpha=0.05, seed=0):
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float)
    if x.size == 0:
        return (np.nan, np.nan)
    bs = [stat(rng.choice(x, x.size, replace=True)) for _ in range(n)]
    return (float(np.percentile(bs, 100 * alpha / 2)),
            float(np.percentile(bs, 100 * (1 - alpha / 2))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coh-deg", type=float, default=20.0,
                    help="desviacion angular por debajo de la cual un pick se considera coherente")
    ap.add_argument("--rect-min", type=float, default=0.60,
                    help="rectilinearidad minima para que el back-azimut sea fiable")
    ap.add_argument("--max-dist-km", type=float, default=200.0)
    args = ap.parse_args()

    if not os.path.exists(CSV):
        sys.exit("falta picks_polarization.csv; corre task1_analyze.py primero")
    df = pd.read_csv(CSV)

    tot = len(df)
    # bool("False") == True en Python: si 'associated' llegara como texto (dtype object) en vez
    # de bool nativo, .astype(bool) invertiria la logica y vaciaria 'na' por completo. En la
    # practica pandas infiere dtype=bool al leer un CSV con "True"/"False" (verificado), pero
    # se blinda igual porque no cuesta nada y elimina la fragilidad para siempre.
    if df["associated"].dtype == object:
        is_assoc = df["associated"].astype(str).str.strip().str.upper().isin(["TRUE", "1"])
        print("   AVISO: 'associated' llego como texto, no como bool nativo; se interpreto "
              "explicitamente en vez de usar .astype(bool).")
    else:
        is_assoc = df["associated"].astype(bool)
    na = df[~is_assoc].copy()
    n_na_raw = len(na)
    na_f1 = na[na["rectilinearity"].notna() & na["angle_dev_deg"].notna()]
    na_f2 = na_f1[na_f1["ev_dist_km"] <= args.max_dist_km]
    usable = na_f2[na_f2["rectilinearity"] >= args.rect_min]
    # Desglose explicito: sin esto, un 'no asociados 0' final no dice en CUAL de los tres
    # filtros se perdieron los picks, y se diagnostica a ciegas.
    if n_na_raw and len(usable) == 0:
        print("   DESGLOSE de por que no quedan picks utilizables:")
        print(f"     no asociados en el CSV              : {n_na_raw}")
        print(f"     tras exigir polarizacion y angulo   : {len(na_f1)}"
              f"  (perdidos: {n_na_raw - len(na_f1)})")
        print(f"     tras exigir ev_dist_km <= {args.max_dist_km:.0f}      : {len(na_f2)}"
              f"  (perdidos: {len(na_f1) - len(na_f2)})")
        print(f"     tras exigir rectilinearidad >= {args.rect_min:.2f}: {len(usable)}"
              f"  (perdidos: {len(na_f2) - len(usable)})")
        if len(na_f1) and len(na_f2) == 0:
            d = na_f1["ev_dist_km"]
            print(f"     -> el filtro de distancia los elimino todos. Rango real de "
                  f"ev_dist_km: {d.min():.0f}-{d.max():.0f} km. Sube --max-dist-km o "
                  f"amplia --box-margin-deg en task1_fetch.py para traer eventos mas cercanos.")
    na = na_f2

    print("=" * 68)
    print("TAREA 1 — COHERENCIA DE POLARIZACION DEL CLUTTER REAL")
    print("=" * 68)
    print(f"  picks totales                     {tot}")
    print(f"  no asociados                      {len(na)}")
    print(f"  con rectilinearidad >= {args.rect_min:.2f}      {len(usable)}")
    if len(usable) < 100:
        print("\n  AVISO: menos de 100 picks utilizables. El resultado sera ruidoso.")
        print("  Amplia el periodo antes de decidir nada.")

    if len(usable) == 0:
        sys.exit("sin picks utilizables")

    dev = usable["angle_dev_deg"].to_numpy(float)
    coherent = (dev <= args.coh_deg)
    frac = float(coherent.mean())
    lo, hi = bootstrap_ci(coherent.astype(float))

    print(f"\n  desviacion angular mediana        {np.median(dev):.1f} grados")
    print(f"  FRACCION COHERENTE (<{args.coh_deg:.0f} grados)   {frac:.3f}  IC95 [{lo:.3f}, {hi:.3f}]")
    print("\n  histograma de desviacion angular")
    h, edges = np.histogram(dev, bins=np.arange(0, 91, 10))
    for k in range(len(h)):
        bar = "#" * int(40 * h[k] / max(h.max(), 1))
        print(f"   {edges[k]:3.0f}-{edges[k+1]:3.0f} | {h[k]:5d} {bar}")

    # desglose por magnitud del vecino y por distancia
    print("\n  desglose por distancia al evento mas cercano")
    for lo_d, hi_d in ((0, 50), (50, 100), (100, 200)):
        m = (usable["ev_dist_km"] >= lo_d) & (usable["ev_dist_km"] < hi_d)
        if m.sum() >= 20:
            f = float((usable.loc[m, "angle_dev_deg"] <= args.coh_deg).mean())
            print(f"   {lo_d:3d}-{hi_d:3d} km | n={int(m.sum()):5d} | fraccion coherente {f:.3f}")

    # control: los picks ASOCIADOS deben ser mayoritariamente coherentes. Si no lo son, algo en
    # la orientacion de componentes o en la asociacion esta mal y el numero de arriba no vale.
    a = df[is_assoc]
    a = a[a["rectilinearity"] >= args.rect_min]
    ctrl = float((a["angle_dev_deg"] <= args.coh_deg).mean()) if len(a) >= 30 else np.nan
    print(f"\n  CONTROL con picks asociados: fraccion coherente {ctrl:.3f}  (n={len(a)})")
    ctrl_ok = (not np.isnan(ctrl)) and ctrl >= 0.60
    if not ctrl_ok:
        print("   FALLA EL CONTROL. Los picks de eventos conocidos deberian apuntar a su evento.")
        print("   Revisa: orientacion de componentes (rotacion a ZNE), banda, y ventana de")
        print("   covarianza. El numero de clutter NO es interpretable hasta que este control pase.")

    # veredicto pre-registrado
    print("\n" + "=" * 68)
    if not ctrl_ok:
        verdict = "INVALIDO: el control con picks asociados no pasa"
        action = "corregir orientacion/banda/ventana y repetir"
    elif hi < 0.20:
        verdict = "BAJA coherencia"
        action = ("los Apendices H, J y M mantienen su validez cuantitativa; el modelo de clutter "
                  "aleatorio es adecuado")
    elif lo > 0.50:
        verdict = "ALTA coherencia"
        action = ("re-medir la jerarquia de canales del Apendice M con clutter coherente; A2 "
                  "probablemente no sobrevive en su forma actual")
    else:
        verdict = "coherencia INTERMEDIA o indeterminada"
        action = ("reescalar los margenes absolutos con la curva del Apendice T y re-medir la "
                  "escalera con un modelo de clutter mixto")
    print(f"VEREDICTO: {verdict}")
    print(f"ACCION:    {action}")
    print("=" * 68)
    print("\ncurva de respuesta del Apendice T, para comparar:")
    for k, v in APPENDIX_T.items():
        print(f"   tasa espuria de la polarizacion con {k:20s}: {v:.2f}")

    res = dict(n_total=tot, n_unassociated=len(na), n_usable=len(usable),
               coh_deg=args.coh_deg, rect_min=args.rect_min,
               median_dev_deg=float(np.median(dev)),
               coherent_fraction=frac, ci95=[lo, hi],
               control_associated_fraction=None if np.isnan(ctrl) else ctrl,
               control_passed=bool(ctrl_ok), verdict=verdict, action=action)
    with open(os.path.join(OUT, "decision.json"), "w") as fh:
        json.dump(res, fh, indent=2)
    print(f"\nescrito {OUT}/decision.json")


if __name__ == "__main__":
    main()
