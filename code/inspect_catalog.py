"""
inspect_catalog.py — revisa el catalogo descargado ANTES de correr la Etapa 2 completa.

Verifica: cuantos eventos, su rango de magnitud, y sobre todo la distancia al centro de la red.
Si todos los eventos estan a cientos de km, la Etapa 3 tendra pocos picks verdaderos asociables
y el numero de coherencia sera ruidoso o indefinido, igual que ya paso una vez con ISC.

    python inspect_catalog.py
"""
import json
import numpy as np
from obspy import read_events
from obspy.geodetics import gps2dist_azimuth

OUT = "t1_data"

with open(f"{OUT}/fetch_report.json") as fh:
    rep = json.load(fh)

sta_lat = [r["lat"] for r in rep["channels"]]
sta_lon = [r["lon"] for r in rep["channels"]]
center_lat, center_lon = float(np.mean(sta_lat)), float(np.mean(sta_lon))
print(f"centro de la red: {center_lat:.3f}, {center_lon:.3f}\n")

cat = read_events(f"{OUT}/catalog.xml")
print(f"eventos en el catalogo: {len(cat)}\n")

rows = []
for e in cat:
    o = e.preferred_origin() or e.origins[0]
    m = e.preferred_magnitude() or (e.magnitudes[0] if e.magnitudes else None)
    d_m, _, _ = gps2dist_azimuth(float(o.latitude), float(o.longitude), center_lat, center_lon)
    rows.append((str(o.time)[:19], float(m.mag) if m else float("nan"),
                m.magnitude_type if m else "?", d_m / 1000.0, float(o.depth or 0) / 1000.0))

rows.sort(key=lambda r: r[3])
mags = [r[1] for r in rows if not np.isnan(r[1])]
dists = [r[3] for r in rows]

print(f"{'tiempo':>20} {'mag':>5} {'tipo':>5} {'dist_km':>9} {'prof_km':>8}")
for r in rows[:15]:
    print(f"{r[0]:>20} {r[1]:5.1f} {r[2]:>5} {r[3]:9.1f} {r[4]:8.1f}")
if len(rows) > 15:
    print(f"   ... y {len(rows)-15} mas")

print(f"\nmagnitud: min {min(mags):.1f}  max {max(mags):.1f}  mediana {np.median(mags):.1f}")
print(f"distancia al centro de la red: min {min(dists):.0f} km  max {max(dists):.0f} km  "
      f"mediana {np.median(dists):.0f} km")

near = sum(1 for d in dists if d <= 300)
print(f"\neventos a <=300 km del centro de la red: {near}/{len(rows)}")
if near < 10:
    print("\nAVISO: muy pocos eventos cercanos. USGS solo captura magnitud considerable")
    print("(tipicamente M>=4) para esta region, lo que sesga el catalogo hacia eventos")
    print("grandes y a menudo distantes. Con esto la Etapa 3 tendra pocos picks verdaderos")
    print("asociables. Considera: (a) bajar --minmag no ayuda con USGS porque no tiene mas")
    print("resolucion; (b) usar el catalogo del SGC directamente, que captura magnitud local")
    print("baja; (c) ampliar --cat-days para acumular mas eventos grandes con el tiempo.")
else:
    print("\nCobertura razonable de eventos cercanos. Puedes seguir a la Etapa 2.")
