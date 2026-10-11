"""
diagnose_distance_mismatch.py — por que find_good_window.py midio 47-73 km y
label_association reporta 365-648 km para los MISMOS eventos.

Hipotesis a confirmar: solo sobrevivieron 9 trazas de las 18 esperadas (6 estaciones x 3
componentes), es decir solo 3 estaciones tienen forma de onda real en esta ventana. Si esas 3
estaciones son justo las MAS LEJANAS de la secuencia del 8 de junio, mientras que las cercanas
(que find_good_window.py SI conto, porque solo mira el catalogo, no la forma de onda) no
tienen dato archivado, eso explica la contradiccion sin ningun bug de calculo.

    python diagnose_distance_mismatch.py
"""
import glob
import json
import numpy as np
from obspy import read, read_events
from obspy.geodetics import gps2dist_azimuth

print("=" * 70)
print("PASO 1: distancias evento-estacion con los datos reales de t1_data/")
print("=" * 70)

with open("t1_data/fetch_report.json") as f:
    rep = json.load(f)
sta = {r["sta"]: (r["lat"], r["lon"]) for r in rep["channels"]}
print(f"\nestaciones en fetch_report.json ({len(sta)}):")
for k, v in sta.items():
    print(f"   {k}: lat={v[0]:.3f} lon={v[1]:.3f}")

cat = read_events("t1_data/catalog.xml")
print(f"\neventos en catalog.xml: {len(cat)}\n")

all_dists = {s: [] for s in sta}
for e in cat:
    o = e.preferred_origin() or e.origins[0]
    ds = {}
    for s, (la, lo) in sta.items():
        d_m, _, _ = gps2dist_azimuth(float(o.latitude), float(o.longitude), la, lo)
        d_km = float(np.hypot(d_m / 1000.0, float(o.depth or 0) / 1000.0))
        ds[s] = d_km
        all_dists[s].append(d_km)
    ordered = sorted(ds.items(), key=lambda x: x[1])
    print(f"{str(o.time)[:19]}  lat={o.latitude:7.3f} lon={o.longitude:8.3f} "
          f"prof={float(o.depth or 0)/1000:5.1f}km")
    print("   " + "  ".join(f"{k}={v:.0f}km" for k, v in ordered))

print(f"\ndistancia MINIMA a cada estacion, sobre todos los eventos:")
for s, ds in sorted(all_dists.items(), key=lambda x: min(x[1])):
    print(f"   {s}: min={min(ds):.0f} km  mediana={np.median(ds):.0f} km")

print("\n" + "=" * 70)
print("PASO 2: que estaciones tienen forma de onda REAL en esta ventana")
print("=" * 70)

files = sorted(glob.glob("t1_data/waveforms/*.mseed"))
print(f"\narchivos mseed encontrados: {len(files)}")
stas_con_dato = {}
for f in files:
    try:
        st = read(f, headonly=True)
    except Exception:
        continue
    for tr in st:
        stas_con_dato.setdefault(tr.stats.station, []).append(
            (str(tr.stats.starttime)[:19], str(tr.stats.endtime)[:19], tr.stats.channel))

print(f"\nestaciones CON forma de onda en los archivos descargados: {sorted(stas_con_dato)}")
faltantes = set(sta) - set(stas_con_dato)
print(f"estaciones EN fetch_report.json pero SIN forma de onda: {sorted(faltantes) or 'ninguna'}")

print("\n" + "=" * 70)
print("PASO 3: veredicto")
print("=" * 70)

if faltantes:
    d_faltantes = {s: min(all_dists[s]) for s in faltantes if s in all_dists}
    d_presentes = {s: min(all_dists[s]) for s in stas_con_dato if s in all_dists}
    print(f"\nDistancia minima de las estaciones SIN dato: "
          f"{sorted(d_faltantes.items(), key=lambda x: x[1])}")
    print(f"Distancia minima de las estaciones CON dato: "
          f"{sorted(d_presentes.items(), key=lambda x: x[1])}")
    if d_faltantes and d_presentes:
        if min(d_faltantes.values()) < min(d_presentes.values()):
            print("\nCONFIRMADO: las estaciones MAS CERCANAS a la secuencia sismica no tienen")
            print("forma de onda archivada en esta ventana. find_good_window.py las conto")
            print("porque solo consulta el catalogo (no la forma de onda), asi que eligio un")
            print("periodo bueno para las estaciones EQUIVOCADAS.")
            print("\nSOLUCION: find_good_window.py debe filtrar, para cada periodo candidato,")
            print("solo las estaciones con disponibilidad de forma de onda CONFIRMADA (via")
            print("probe_waveforms.py o similar) antes de calcular distancias minimas.")
        else:
            print("\nNo se confirma la hipotesis simple. Revisar mas.")
else:
    print("\nTodas las estaciones de fetch_report.json tienen forma de onda. La causa de la")
    print("discrepancia de distancias debe buscarse en otro lugar (revisar si catalog.xml")
    print("cambio entre la corrida de find_good_window.py y la de task1_fetch.py).")
