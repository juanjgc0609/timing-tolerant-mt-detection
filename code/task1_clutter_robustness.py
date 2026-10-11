"""
task1_clutter_robustness.py — el veredicto de clutter de la Tarea 1, por subpoblacion.

task1_decide.py da una sola fraccion coherente corregida por azar (eta) para todo el clutter de una
estacion. Si eta es bajo en conjunto pero alto en alguna subpoblacion, el veredicto depende de como
se muestreo el clutter y no es una propiedad del clutter real. Este script recalcula eta, con el
mismo estimador y el mismo control, por separado para:

    confianza del picker      prob < 0.5  /  prob >= 0.5
    retardo en la coda        < 30 min    /  >= 30 min despues de la P del evento de referencia
    distancia al evento       < 150 km    /  >= 150 km
    fase asignada por PhaseNet  P  /  S

    python task1_clutter_robustness.py --stations ZAR

Usa el CSV que dejo la ultima corrida de task1_analyze.py (debe ser la de la estacion y la
configuracion del veredicto). No cambia el veredicto registrado: lo somete a prueba.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

from task1_decide import clutter_share, strata

OUT = "t1_data"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stations", default="ZAR")
    ap.add_argument("--coh-deg", type=float, default=20.0)
    ap.add_argument("--rect-min", type=float, default=0.60)
    ap.add_argument("--max-dist-km", type=float, default=200.0)
    ap.add_argument("--min-n", type=int, default=30)
    args = ap.parse_args()

    csv = os.path.join(OUT, "picks_polarization.csv")
    if not os.path.exists(csv):
        sys.exit("falta picks_polarization.csv")
    df = pd.read_csv(csv)
    df = df[df["sta"].isin(args.stations.split(","))].copy()
    if df.empty:
        sys.exit(f"sin picks de {args.stations} en el CSV: corre antes task1_analyze.py para esa estacion")
    if df["associated"].dtype == object:
        assoc = df["associated"].astype(str).str.strip().str.upper().isin(["TRUE", "1"])
    else:
        assoc = df["associated"].astype(bool)
    ok = df["rectilinearity"] >= args.rect_min
    ctrl = df[assoc & ok & (df["phase"] == "P") & df["angle_dev_deg"].notna()]
    clut = df[~assoc & ok & (df["ev_dist_km"] <= args.max_dist_km) & df["angle_dev_deg"].notna()].copy()
    if "ref_lag_s" not in clut.columns:
        sys.exit("el CSV no tiene ref_lag_s: viene de una version anterior de task1_analyze.py")
    print(f"   control P {len(ctrl)} | clutter utilizable {len(clut)}")
    ctrl_m = strata(ctrl, args.coh_deg)

    subsets = [("todo el clutter", np.ones(len(clut), bool))]
    subsets += [("prob < 0.5", clut["prob"] < 0.5), ("prob >= 0.5", clut["prob"] >= 0.5),
                ("coda < 30 min", clut["ref_lag_s"] < 1800), ("coda >= 30 min", clut["ref_lag_s"] >= 1800),
                ("distancia < 150 km", clut["ev_dist_km"] < 150),
                ("distancia >= 150 km", clut["ev_dist_km"] >= 150),
                ("fase P", clut["phase"] == "P"), ("fase S", clut["phase"] == "S")]
    print(f"\n   {'subpoblacion':<22}{'n':>5}{'f':>8}{'azar f0':>9}{'eta':>8}   IC95 de eta          lectura")
    for name, m in subsets:
        sub = clut[np.asarray(m, bool)]
        if len(sub) < args.min_n:
            print(f"   {name:<22}{len(sub):>5}   (menos de {args.min_n}, no se estima)")
            continue
        s = clutter_share(strata(sub, args.coh_deg), ctrl_m)
        if s["hi"] < 0.20:
            lectura = "bajo 0.20"
        elif s["lo"] > 0.20:
            lectura = "SOBRE 0.20"
        else:
            lectura = "cruza 0.20"
        print(f"   {name:<22}{len(sub):>5}{s['f']:>8.3f}{s['f0']:>9.3f}{s['pi']:>8.3f}"
              f"   [{s['lo']:+.3f}, {s['hi']:+.3f}]   {lectura}")
    print("\n   Si todas las subpoblaciones con n suficiente quedan bajo 0.20 o cruzan con estimacion")
    print("   cercana a cero, el veredicto BAJA no depende de como se muestreo el clutter. Una")
    print("   subpoblacion con el intervalo entero SOBRE 0.20 lo contradiria y habria que reportarla.")
    print("   Son 8 comparaciones sin correccion: un intervalo que apenas roza el umbral en una sola")
    print("   de ellas es lo esperable por azar.")


if __name__ == "__main__":
    main()
