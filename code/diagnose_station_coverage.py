"""
diagnose_station_coverage.py — por que 'load_streams' cargo solo 9 trazas y por que el
evento mas cercano a los picks esta a cientos de km.

Corre DESPUES de task1_fetch.py, ANTES de gastar mas tiempo en task1_analyze.py.

    python diagnose_station_coverage.py
"""
import glob
import json
import os
from collections import defaultdict

import numpy as np
from obspy import read_events, UTCDateTime
from obspy.geodetics import gps2dist_azimuth

OUT = "t1_data"

print("=" * 68)
print("PASO 1: archivos mseed por estacion")
print("=" * 68)
files = sorted(glob.glob(os.path.join(OUT, "waveforms", "*.mseed")))
by_sta = defaultdict(lambda: defaultdict(int))
for f in files:
    base = os.path.basename(f)
    parts = base.split(".")
    if len(parts) < 3:
        continue
    net, sta, loc = parts[0], parts[1], parts[2]
    cha = parts[3].split("__")[0] if len(parts) > 3 else "?"
    by_sta[sta][cha] += 1

print(f"total archivos: {len(files)}")
for sta in sorted(by_sta):
    chs = by_sta[sta]
    total = sum(chs.values())
    print(f"   {sta:6s} total={total:4d}  " + "  ".join(f"{c}={n}" for c, n in sorted(chs.items())))

print("\n" + "=" * 68)
print("PASO 2: rango temporal REAL cubierto por estacion (lee cabeceras, no todo el archivo)")
print("=" * 68)
from obspy.io.mseed.util import get_start_and_end_time
extent = defaultdict(lambda: [None, None])
bad = defaultdict(int)
for f in files:
    base = os.path.basename(f)
    sta = base.split(".")[1]
    try:
        fs, fe = get_start_and_end_time(f)
    except Exception as e:
        bad[sta] += 1
        continue
    lo, hi = extent[sta]
    extent[sta][0] = fs if lo is None or fs < lo else lo
    extent[sta][1] = fe if hi is None or fe > hi else hi

with open(os.path.join(OUT, "fetch_report.json")) as fh:
    rep = json.load(fh)
t0, t1 = UTCDateTime(rep["start"]), UTCDateTime(rep["end"])
print(f"ventana pedida en el fetch: {t0} .. {t1}\n")
for sta in sorted(extent):
    fs, fe = extent[sta]
    cov = "OK" if (fs is not None and fs <= t0 + 3600 and fe >= t1 - 3600) else "INCOMPLETA/AUSENTE"
    print(f"   {sta:6s} {str(fs):>26s} -> {str(fe):>26s}   [{cov}]"
          + (f"   ({bad[sta]} archivo(s) ilegible(s))" if bad[sta] else ""))
faltantes = set(rep["stations"].split(",")) - set(extent.keys())
if faltantes:
    print(f"\n   ESTACIONES SIN NINGUN ARCHIVO EN t1_data/waveforms: {sorted(faltantes)}")

print("\n" + "=" * 68)
print("PASO 3: distancia estacion -> cada evento del catalogo (para ver si alguna estacion")
print("        realmente tiene el swarm cerca, y cual no)")
print("=" * 68)
sta_coords = {r["sta"]: (r["lat"], r["lon"]) for r in rep["channels"]}
cat = read_events(os.path.join(OUT, "catalog.xml"))
ev = []
for k, e in enumerate(cat):
    o = e.preferred_origin() or e.origins[0]
    ev.append((k, o.time, float(o.latitude), float(o.longitude)))

for sta, (slat, slon) in sorted(sta_coords.items()):
    dists = []
    for k, ot, elat, elon in ev:
        d_m, _, _ = gps2dist_azimuth(elat, elon, slat, slon)
        dists.append(d_m / 1000.0)
    dists.sort()
    tag = " <-- tiene data" if sta in extent else " <-- SIN ARCHIVOS"
    print(f"   {sta:6s} dist a los {len(ev)} eventos (km), ordenadas: "
          f"{['%.0f' % d for d in dists[:5]]} ...{tag}")

print("\nSi una estacion listada como '<-- tiene data' muestra distancias de cientos de km")
print("a TODOS los eventos, esa estacion simplemente esta lejos del swarm (no es un bug).")
print("Si CHI (u otra estacion que sabes que esta a 47-73 km) aparece '<-- SIN ARCHIVOS' o")
print("'INCOMPLETA/AUSENTE' en el paso 2, ese es el problema: sus datos no llegaron esta vez.")
