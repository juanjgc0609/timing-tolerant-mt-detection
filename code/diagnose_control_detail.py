"""
diagnose_control_detail.py — mira el detalle de los picks ASOCIADOS con mayor desviacion
angular, para distinguir dos hipotesis distintas que ambas producen un control bajo:

  (a) el pick esta bien asociado al evento correcto, pero la polarizacion medida en dato real
      simplemente tiene mas dispersion que en el sintetico limpio de task1_verify.py (ruido de
      sitio, heterogeneidad cortical, ondas no planas -- limitacion FISICA real, no un bug), o

  (b) el pick quedo asociado al evento EQUIVOCADO dentro de la ventana de tolerancia (varios
      eventos del catalogo caen cerca en tiempo), y por eso se compara contra un azimut que no
      es el que realmente genero el pick.

Si los peores casos tienen ev_dist_km grande y/o quedan cerca en tiempo de OTRO evento del
catalogo, apunta a (b). Si tienen ev_dist_km chico y rectilinearidad alta y aun asi el azimut no
cuadra, apunta a (a) -- y ahi si hay que sospechar orientacion/banda/ventana, o aceptar el nivel
de dispersion como propio del dato real.

    python diagnose_control_detail.py
"""
import pandas as pd
import numpy as np
from obspy import read_events, UTCDateTime
from obspy.geodetics import gps2dist_azimuth

df = pd.read_csv("t1_data/picks_polarization.csv")
is_assoc = df["associated"].astype(bool) if df["associated"].dtype != object else \
    df["associated"].astype(str).str.strip().str.upper().isin(["TRUE", "1"])
a = df[is_assoc].copy()
a = a[a["rectilinearity"] >= 0.60].copy()

print(f"picks de control (asociados, rect>=0.60): {len(a)}\n")
print(f"{'sta':5s} {'ph':3s} {'rect':>5s} {'dev_deg':>8s} {'ev_dist_km':>11s} "
      f"{'best_dt_s':>10s} {'ev_index':>8s}")
a_sorted = a.sort_values("angle_dev_deg", ascending=False)
for _, r in a_sorted.iterrows():
    print(f"{r['sta']:5s} {r['phase']:3s} {r['rectilinearity']:5.2f} "
          f"{r['angle_dev_deg']:8.1f} {r['ev_dist_km']:11.0f} "
          f"{r['best_dt_s']:10.1f} {int(r['ev_index']):8d}")

print("\n" + "=" * 60)
print("por estacion")
print("=" * 60)
for sta in sorted(a["sta"].unique()):
    sub = a[a["sta"] == sta]
    coh = (sub["angle_dev_deg"] <= 20.0).mean()
    print(f"   {sta}: n={len(sub):4d}  fraccion coherente={coh:.3f}  "
          f"mediana_dev={sub['angle_dev_deg'].median():.1f}")

# revisa si hay OTRO evento del catalogo mas cerca en tiempo del que quedo asignado, para los
# 10 peores casos -- si lo hay, el desempate pudo elegir mal.
print("\n" + "=" * 60)
print("chequeo de posible evento alternativo para los 10 peores casos")
print("=" * 60)
cat = read_events("t1_data/catalog.xml")
ev = []
for k, e in enumerate(cat):
    o = e.preferred_origin() or e.origins[0]
    ev.append((k, o.time, float(o.latitude), float(o.longitude), float(o.depth or 0) / 1000.0))

import json
with open("t1_data/fetch_report.json") as fh:
    rep = json.load(fh)
sta_coords = {r["sta"]: (r["lat"], r["lon"]) for r in rep["channels"]}

for _, r in a_sorted.head(10).iterrows():
    pt = UTCDateTime(r["t"])
    slat, slon = sta_coords[r["sta"]]
    dists_t = []
    for k, ot, elat, elon, edep in ev:
        d_m, _, az_sta_ev = gps2dist_azimuth(elat, elon, slat, slon)
        d_km = float(np.hypot(d_m / 1000.0, edep))
        dists_t.append((abs(float(pt - ot)), k, d_km, az_sta_ev))
    dists_t.sort()
    print(f"\n  pick {r['sta']} {r['phase']} t={pt}  medido_az={r['backazimuth_meas']:.0f}  "
          f"ev_index_asignado={int(r['ev_index'])}  dev={r['angle_dev_deg']:.1f}")
    for dt, k, dkm, az in dists_t[:3]:
        marca = " <- asignado" if k == int(r["ev_index"]) else ""
        print(f"      evento {k:3d}  |t_pick-t_origen|={dt:8.1f}s  dist={dkm:6.0f}km  "
              f"az_verdadero={az:6.1f}{marca}")
