"""
diagnose_phase_split.py — separa el control de coherencia por fase (P vs S) usando el CSV que
ya existe, sin volver a correr el picker ni la asociacion.

Hipotesis a probar: el control 'picks asociados deberian apuntar a su evento' falla (32% en vez
de >=60%) no por un problema de orientacion/banda, sino porque se estan mezclando picks P y S.
El modelo de polarizacion validado en task1_verify.py es el de una onda P rectilinea, cuya
proyeccion horizontal apunta a lo largo del back-azimut. Una onda S (en particular SH) puede
estar polarizada TRANSVERSAL al back-azimut, no a lo largo de el; mezclarla con P en la misma
comparacion diluiria cualquier coherencia real de las P.

    python diagnose_phase_split.py
"""
import pandas as pd
import numpy as np

df = pd.read_csv("t1_data/picks_polarization.csv")
print(f"total picks: {len(df)}")
print(f"fases presentes: {df['phase'].value_counts().to_dict()}\n")

is_assoc = df["associated"].astype(bool) if df["associated"].dtype != object else \
    df["associated"].astype(str).str.strip().str.upper().isin(["TRUE", "1"])

a = df[is_assoc].copy()
a = a[a["rectilinearity"] >= 0.60]
print(f"picks ASOCIADOS con rectilinearidad >= 0.60: {len(a)}\n")

print("=" * 60)
print("CONTROL desglosado por fase (deberia ser >= 0.60 para 'pasar')")
print("=" * 60)
for ph in sorted(a["phase"].dropna().unique()):
    sub = a[a["phase"] == ph]
    if len(sub) == 0:
        continue
    coh = (sub["angle_dev_deg"] <= 20.0).mean()
    print(f"   fase {ph}: n={len(sub):4d}  fraccion coherente={coh:.3f}  "
          f"mediana_dev={sub['angle_dev_deg'].median():.1f} grados")

coh_all = (a["angle_dev_deg"] <= 20.0).mean()
print(f"\n   TODAS las fases juntas (como lo hace task1_decide.py hoy): "
      f"n={len(a)}  fraccion coherente={coh_all:.3f}")

print("\nSi 'fase P' pasa >=0.60 sola pero 'fase S' esta muy por debajo (o cerca de un pico en")
print("~90 grados de desviacion, indicando SH transversal), la hipotesis de mezcla de fases se")
print("confirma: la correccion es filtrar el control (y el clutter) a phase=='P' unicamente,")
print("no tocar orientacion/banda/ventana.")
