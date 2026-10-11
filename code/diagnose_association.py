"""
diagnose_association.py — revisa si el etiquetado associated=True tiene sentido temporal.

Para cada pick marcado como asociado, imprime el tiempo del pick, el tiempo del evento que se
le asigno (via ev_index), y la diferencia en DIAS. Si esa diferencia es de meses, confirma un
bug en la logica de comparacion, no un resultado real.

    python diagnose_association.py
"""
import pandas as pd
from obspy import read_events, UTCDateTime

df = pd.read_csv("t1_data/picks_polarization.csv")
cat = read_events("t1_data/catalog.xml")

print(f"dtype de 'associated': {df['associated'].dtype}")
print(f"valores unicos: {df['associated'].unique()}\n")

ev_times = {}
for k, e in enumerate(cat):
    o = e.preferred_origin() or e.origins[0]
    ev_times[k] = o.time

print(f"catalogo: {len(cat)} eventos, rango "
      f"{min(ev_times.values())} .. {max(ev_times.values())}\n")

assoc = df[df["associated"] == True]
print(f"picks marcados associated=True: {len(assoc)}\n")

if len(assoc):
    print(f"{'pick_t':>26} {'ev_index':>8} {'ev_t':>26} {'delta_dias':>12}")
    for _, r in assoc.head(15).iterrows():
        pt = UTCDateTime(r["t"])
        if pd.isna(r["ev_index"]):
            print(f"{str(pt):>26} {'NaN':>8}  (sin evento asignado; ev_index vacio en el CSV)")
            continue
        ek = int(r["ev_index"])
        et = ev_times.get(ek)
        delta_days = (pt - et) / 86400.0 if et else float("nan")
        print(f"{str(pt):>26} {ek:8d} {str(et):>26} {delta_days:12.1f}")
    all_deltas = []
    for _, r in assoc.iterrows():
        if pd.isna(r["ev_index"]):
            continue
        pt = UTCDateTime(r["t"]); et = ev_times.get(int(r["ev_index"]))
        if et:
            all_deltas.append(abs((pt - et) / 86400.0))
    if not all_deltas:
        print("\nsin deltas calculables (todos los ev_index eran NaN o -1)")
    else:
        print(f"\ndelta en dias: min {min(all_deltas):.1f}  max {max(all_deltas):.1f}  "
              f"mediana {sorted(all_deltas)[len(all_deltas)//2]:.1f}")
    if all_deltas and min(all_deltas) > 1.0:
        print("\nCONFIRMADO: los picks 'asociados' estan a DIAS de distancia de su evento.")
        print("Esto es imposible bajo la ventana de segundos que usa label_association.")
        print("Hay un bug en la comparacion de tiempos, no un resultado real.")
