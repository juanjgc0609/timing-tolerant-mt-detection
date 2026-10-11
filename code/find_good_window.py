"""
find_good_window.py — elige un periodo con eventos CERCANOS antes de descargar nada.

El problema detectado: en 2025-10-01..10-31 solo hubo 6 eventos y todos lejanos (237-876 km),
asi que ninguna ventana de coda produce clutter con un vecino cercano y la Tarea 1 no puede
medir lo que pretende. Ampliar --box-margin-deg lo empeoro: trajo eventos aun mas lejanos.

Este script consulta SOLO el catalogo (rapido, sin descargar formas de onda) y busca la
ventana de N dias con mas eventos a menos de max_dist_km de alguna estacion.

    python find_good_window.py --days 30 --max-dist-km 150
"""
import argparse
import json
import numpy as np
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth

EVENT_PROVIDERS = ["USGS", "EMSC", "INGV", "ISC"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--days", type=float, default=30.0)
    ap.add_argument("--max-dist-km", type=float, default=150.0)
    ap.add_argument("--lookback-days", type=float, default=730.0,
                    help="cuanto hacia atras buscar")
    ap.add_argument("--minmag", type=float, default=2.5,
                    help="bajar esto trae mas eventos pero USGS tiene completitud ~M4 aqui")
    args = ap.parse_args()

    with open("t1_data/fetch_report.json") as fh:
        rep = json.load(fh)
    sta = {r["sta"]: (r["lat"], r["lon"]) for r in rep["channels"]}
    print(f"{len(sta)} estaciones: {sorted(sta)}\n")

    lats = [v[0] for v in sta.values()]; lons = [v[1] for v in sta.values()]
    # caja AJUSTADA: solo lo necesario para cubrir max_dist_km desde las estaciones.
    # Ampliarla de mas trae eventos lejanos que contaminan el analisis (lo aprendido).
    m = args.max_dist_km / 111.0
    box = (min(lats)-m, max(lats)+m, min(lons)-m, max(lons)+m)
    print(f"caja ajustada a {args.max_dist_km:.0f} km: "
          f"lat[{box[0]:.2f},{box[1]:.2f}] lon[{box[2]:.2f},{box[3]:.2f}]\n")

    now = UTCDateTime()
    t_end = now - 75 * 86400          # margen de rezago de catalogo
    t_start = t_end - args.lookback_days * 86400

    cat = None
    for prov in EVENT_PROVIDERS:
        try:
            cat = Client(prov, timeout=60).get_events(
                starttime=t_start, endtime=t_end, minmagnitude=args.minmag,
                minlatitude=box[0], maxlatitude=box[1],
                minlongitude=box[2], maxlongitude=box[3])
            if cat.count():
                print(f"catalogo: {cat.count()} eventos desde {prov}\n")
                break
        except Exception as e:
            print(f"   {prov} fallo: {type(e).__name__}")
    if not cat or not cat.count():
        raise SystemExit("ningun proveedor devolvio eventos en la caja ajustada")

    # distancia minima de cada evento a CUALQUIER estacion
    rows = []
    for e in cat:
        o = e.preferred_origin() or e.origins[0]
        dmin = min(
            np.hypot(gps2dist_azimuth(float(o.latitude), float(o.longitude), la, lo)[0] / 1000.0,
                     float(o.depth or 0) / 1000.0)
            for la, lo in sta.values())
        if dmin <= args.max_dist_km:
            rows.append((o.time, float(dmin)))
    rows.sort()
    print(f"eventos a <= {args.max_dist_km:.0f} km de alguna estacion: {len(rows)}\n")
    if not rows:
        raise SystemExit(f"ninguno. Sube --max-dist-km o baja --minmag "
                         f"(ojo: USGS tiene completitud ~M4 en esta region).")

    # ventana deslizante
    best = (0, None, [])
    times = [r[0] for r in rows]
    for i, t0 in enumerate(times):
        t1 = t0 + args.days * 86400
        grp = [r for r in rows if t0 <= r[0] <= t1]
        if len(grp) > best[0]:
            best = (len(grp), t0, grp)

    n, t0, grp = best
    print(f"MEJOR ventana de {args.days:.0f} dias: {n} eventos cercanos")
    print(f"   --start {str(t0)[:10]}\n")
    for t, d in grp:
        print(f"   {str(t)[:19]}   {d:6.0f} km")
    print(f"\nComando sugerido:")
    print(f"python task1_fetch.py --start {str(t0)[:10]} --days {args.days:.0f} \\")
    print(f"    --stations {','.join(sorted(sta))} \\")
    print(f"    --cat-start {str(t0)[:10]} --cat-days {args.days:.0f} \\")
    print(f"    --box-margin-deg {args.max_dist_km/111.0:.1f}")


if __name__ == "__main__":
    main()
