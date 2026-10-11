"""
diagnose_pick_station.py — el veredicto anterior descarto la hipotesis simple: CHI SI tiene
forma de onda y esta a 47-73 km, la mas cercana de todas. Solo falta RUS (149 km).

Entonces la pregunta correcta es: en que ESTACION cayeron los 62 picks no asociados (clutter)
que dieron ev_dist_km = 365-648 km? Si ninguno cayo en CHI, eso explica todo: el picker no
declaro clutter en la unica estacion cercana, asi que todo el clutter que si se detecto viene
de estaciones lejanas cuyo evento mas proximo en tiempo esta, por definicion, lejos.

    python diagnose_pick_station.py
"""
import pandas as pd

df = pd.read_csv("t1_data/picks_polarization.csv")
print(f"picks totales: {len(df)}\n")

print("picks por estacion, y cuantos de cada una son clutter (no asociado):")
print(f"{'estacion':>10} {'total':>7} {'asociados':>10} {'clutter':>8}")
for s, g in df.groupby("sta"):
    assoc = int(g["associated"].astype(bool).sum())
    clut = len(g) - assoc
    print(f"{s:>10} {len(g):7d} {assoc:10d} {clut:8d}")

print("\nde los picks NO asociados (clutter), distribucion de ev_dist_km por estacion:")
na = df[~df["associated"].astype(bool)]
if len(na):
    print(na.groupby("sta")["ev_dist_km"].describe()[["count", "min", "50%", "max"]])
else:
    print("no hay picks no asociados")

print("\nveredicto:")
chi_clutter = na[na["sta"] == "CHI"] if len(na) else pd.DataFrame()
if len(chi_clutter) == 0:
    print("CHI (la estacion mas cercana, 47-73 km) NO aporto NINGUN pick de clutter.")
    print("Eso es fisicamente razonable: la coda de un sismo pequeno (prof ~10km, magnitud")
    print("baja segun USGS) decae rapido, y a corta distancia la senal es mas limpia, no mas")
    print("ruidosa. El clutter que SI aparecio viene de estaciones lejanas (HEL/ZAR/OCA/SJC,")
    print("300-650 km), cuyo evento mas cercano en tiempo esta, por definicion, lejos.")
    print()
    print("Esto NO es un bug de label_association. Es que la ventana de coda en la estacion")
    print("CERCANA simplemente no produjo detecciones espurias con este umbral de picker,")
    print("mientras que las estaciones LEJANAS si detectan clutter, pero ese clutter no tiene")
    print("un evento cercano con el que compararse.")
    print()
    print("Implicacion para la Tarea 1: el filtro --max-dist-km 200 excluye correctamente el")
    print("clutter de estaciones lejanas, porque para ESAS estaciones no hay ningun evento")
    print("real a menos de 200 km. El problema no es el filtro: es que la unica estacion util")
    print("(CHI) no genero suficiente clutter para construir una muestra.")
else:
    print(f"CHI SI aporto {len(chi_clutter)} picks de clutter, con ev_dist_km "
          f"{chi_clutter['ev_dist_km'].tolist()}")
