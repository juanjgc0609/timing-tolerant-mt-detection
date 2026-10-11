"""
task1_fetch.py — Etapa 1 de la Tarea 1: traer el dato.

Descarga registro CONTINUO (no ventanas de evento: el clutter es el objeto de estudio y las
ventanas de evento lo borran), el inventario con respuesta y orientación, y el catálogo del mismo
periodo desde el primer proveedor que responda.

    python task1_fetch.py --start 2024-03-01 --days 2 --stations OCA,RUS,BRR,PTB,ZAR,HEL,CHI,SJC

Salidas en ./t1_data/
    waveforms/*.mseed      trozos horarios
    stations/*.xml         StationXML con respuesta
    inventory.xml          inventario consolidado
    catalog.xml            catálogo QuakeML
    fetch_report.json      qué se obtuvo, de dónde, y los avisos
"""
import argparse
import json
import os
import sys
from datetime import datetime

import numpy as np
from obspy import UTCDateTime, read_inventory
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.mass_downloader import (RectangularDomain, Restrictions,
                                                MassDownloader)

OUT = "t1_data"
WAVE = os.path.join(OUT, "waveforms")
STAX = os.path.join(OUT, "stations")

# proveedores de ondas y de catálogo, en orden de intento.
# El SGC (Servicio Geologico Colombiano) opera la red CM directamente y tiene el continuo real;
# EarthScope/IRIS solo reflejan metadatos de CM sin el archivo de forma de onda completo
# (verificado: 12/12 combinaciones fallidas en EarthScope/IRIS para CM.OCA en dos fechas
# distintas, mientras que el SGC respondio con HHZ real en las 10/10 fechas probadas;
# BHZ no esta archivado alli, solo HHZ en loc=00).
SGC_URL = "http://sismo.sgc.gov.co:8080"
WAVE_PROVIDERS = [SGC_URL, "EARTHSCOPE", "IRIS"]
EVENT_PROVIDERS = ["USGS", "EMSC", "INGV", "ISC"]  # ISC al final: revisado, con meses de rezago


def ensure_dirs():
    for d in (OUT, WAVE, STAX):
        os.makedirs(d, exist_ok=True)


def make_client(prov, timeout=30):
    """
    SGC_URL es una URL completa (base_url); los demas son mnemonicos FDSN estandar.
    timeout explicito: sin esto, ObsPy usa 120s por defecto, y MassDownloader puede acumular
    docenas de reintentos silenciosos contra un servidor lento o caido (visto: mas de 3 horas
    sin progreso visible en una descarga de 30 dias x 8 estaciones). 30s es agresivo a
    proposito para fallar rapido y pasar al siguiente proveedor en vez de esperar.
    """
    if prov.startswith("http"):
        return Client(base_url=prov, timeout=timeout)
    return Client(prov, timeout=timeout)


def probe_availability(net, stations, days_back, probe_len_s=120):
    """
    El servicio fdsnws-availability de EarthScope esta permanentemente fuera de servicio
    (HTTP 410) desde su migracion a la nube (aviso oficial de EarthScope, 2026). No hay forma de
    consultar cobertura antes de pedir, asi que se tantea con get_waveforms en ventanas cortas
    (probe_len_s), retrocediendo desde HOY, hasta encontrar un dia con dato real para al menos
    una estacion. Devuelve (t0_encontrado, reporte_por_estacion) o (None, reporte) si nada sirvio.
    """
    from datetime import datetime, timedelta
    stas = stations.split(",")
    now = UTCDateTime(datetime.utcnow()) - 3 * 86400   # margen de latencia de archivado
    report = []
    for back in days_back:
        day0 = now - back * 86400
        day0 = UTCDateTime(day0.date)
        hit = []
        for prov in WAVE_PROVIDERS:
            c = make_client(prov)
            for sta in stas:
                try:
                    st = c.get_waveforms(net, sta, "*", "?H?", day0, day0 + probe_len_s,
                                         attach_response=False)
                    if len(st) and any(tr.stats.npts > 0 for tr in st):
                        hit.append((sta, prov, len(st)))
                except Exception:
                    continue
            if hit:
                break
        report.append(dict(days_back=back, date=str(day0.date), hits=hit, provider=prov))
        print(f"   sondeo {back:4d} dias atras ({day0.date}): "
              f"{len(hit)} estaciones con dato" + (f" via {prov}" if hit else ""))
        if hit:
            return day0, report
    return None, report


def fetch_inventory(net, stations, t0, t1):
    """Inventario a nivel de canal y respuesta. Aquí se decide si habrá que rotar."""
    last = None
    for prov in WAVE_PROVIDERS:
        try:
            c = make_client(prov)
            inv = c.get_stations(network=net, station=stations, starttime=t0, endtime=t1,
                                 level="response", channel="?H?")
            if len(inv.get_contents()["channels"]):
                return inv, prov
        except Exception as e:                                     # noqa: BLE001
            last = e
    raise RuntimeError(f"ningun proveedor devolvio inventario: {last}")


def orientation_report(inv):
    """
    Devuelve (filas, necesita_rotacion). Es la comprobación crítica de toda la tarea: si las
    horizontales no están en N y E, los back-azimuts salen rotados y la medición no vale.
    """
    rows, rotate = [], False
    for net in inv:
        for sta in net:
            for ch in sta:
                rows.append(dict(sta=sta.code, loc=ch.location_code or "--", cha=ch.code,
                                 sr=float(ch.sample_rate),
                                 az=float(ch.azimuth) if ch.azimuth is not None else None,
                                 dip=float(ch.dip) if ch.dip is not None else None,
                                 lat=float(sta.latitude), lon=float(sta.longitude),
                                 elev=float(sta.elevation)))
                code = ch.code[-1]
                if code in ("1", "2"):
                    rotate = True
                elif code in ("N", "E") and ch.azimuth is not None:
                    expect = 0.0 if code == "N" else 90.0
                    if abs(((ch.azimuth - expect + 180) % 360) - 180) > 1.0:
                        rotate = True
    return rows, rotate


def fetch_catalog(t0, t1, box, minmag):
    """
    ISC es un catalogo REVISADO con varios meses de rezago tipico: si t0 es reciente (el sondeo
    de disponibilidad puede fijar t0 en el dia de hoy), ISC devolvera 0 eventos aunque haya
    sismicidad real, y eso no significa que el periodo sea sismicamente tranquilo. USGS y EMSC
    son casi-tiempo-real y se prueban primero.
    """
    from obspy import UTCDateTime as _U
    age_days = (_U(datetime.utcnow()) - t0) / 86400.0
    if age_days < 60:
        print(f"   AVISO: la ventana empieza hace solo {age_days:.0f} dias. ISC (catalogo "
              f"revisado) tipicamente tarda meses en publicar; si los proveedores en tiempo "
              f"real tambien dan 0, NO asumas que no hubo sismicidad: es rezago de catalogo,")
        print("   no ausencia de eventos. Repite esta consulta mas adelante, o usa un --start")
        print("   de hace >=6 meses para que el catalogo revisado ya este disponible.")
    for prov in EVENT_PROVIDERS:
        try:
            c = Client(prov)
            cat = c.get_events(starttime=t0, endtime=t1,
                               minlatitude=box[0], maxlatitude=box[1],
                               minlongitude=box[2], maxlongitude=box[3],
                               minmagnitude=minmag)
            if cat.count():
                return cat, prov
        except Exception as e:                                     # noqa: BLE001
            print(f"   [catalogo] {prov} fallo: {e}", file=sys.stderr)
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", default="CM")
    ap.add_argument("--stations", default="OCA,RUS,BRR,PTB,ZAR,HEL,CHI,SJC")
    ap.add_argument("--start", default=None,
                    help="fecha fija YYYY-MM-DD; si se omite, se SONDEA automaticamente")
    ap.add_argument("--probe-back", default="7,30,90,180,365",
                    help="dias hacia atras desde hoy a sondear si --start no se da")
    ap.add_argument("--days", type=float, default=2.0)
    ap.add_argument("--skip-size-warning", action="store_true",
                    help="omite el aviso de tamano estimado de descarga para --days grande")
    ap.add_argument("--minmag", type=float, default=1.5)
    ap.add_argument("--cat-start", default=None,
                    help="fecha de INICIO del catalogo si difiere de --start (ej. si --start "
                         "es muy reciente y el catalogo revisado no ha llegado alli todavia)")
    ap.add_argument("--cat-days", type=float, default=None,
                    help="duracion de la ventana de catalogo si difiere de --days")
    ap.add_argument("--cat-margin-days", type=float, default=0.0,
                    help="dias de margen ANTES y DESPUES de la ventana de catalogo, SOLO para "
                         "inspeccionar contexto sismico regional (p.ej. con inspect_catalog.py). "
                         "NO ayuda a conseguir mas picks asociados: un pick solo puede asociarse "
                         "a un evento cuyo tiempo de llegada (origen + segundos de viaje, nunca "
                         "minutos) caiga DENTRO de la ventana de forma de onda descargada. Un "
                         "evento 10 dias antes o despues nunca puede producir un pick asociado, "
                         "sin importar el margen. Para conseguir mas picks asociados: (a) amplia "
                         "--days (la ventana de FORMA DE ONDA, no la de catalogo), o (b) usa un "
                         "catalogo con menor magnitud de completitud para esta region.")
    ap.add_argument("--box", default=None,
                    help="minlat,maxlat,minlon,maxlon; por defecto se calcula desde el "
                         "inventario con margen, para no descartar estaciones pedidas")
    ap.add_argument("--box-margin-deg", type=float, default=2.5,
                    help="margen en grados alrededor de las estaciones para buscar catalogo. "
                         "task1_decide.py mide coherencia hasta --max-dist-km 200 (por "
                         "defecto); 200 km equivale a ~1.8 grados de latitud, asi que el "
                         "margen debe ser mayor a eso o eventos reales entre 110 y 200 km "
                         "quedarian fuera de catalog.xml y se contarian como clutter "
                         "incoherente sin serlo. 2.5 grados (~275 km) da margen adicional.")
    args = ap.parse_args()

    ensure_dirs()
    if args.start is not None:
        t0 = UTCDateTime(args.start)
    else:
        print("0/3 sondeo de disponibilidad (el servicio fdsnws-availability de EarthScope")
        print("    esta fuera de servicio permanentemente; se tantea con get_waveforms)")
        back_list = [int(x) for x in args.probe_back.split(",")]
        t0, probe_report = probe_availability(args.net, args.stations, back_list)
        if t0 is None:
            sys.exit("ningun dia sondeado tuvo dato para ninguna estacion. Amplia "
                     "--probe-back (ej. 7,30,60,90,180,365,730) o revisa la red/estaciones.")
        print(f"   dia elegido: {t0.date}")
    t1 = t0 + args.days * 86400
    report = {"start": str(t0), "end": str(t1), "net": args.net,
              "stations": args.stations, "generated": datetime.utcnow().isoformat()}

    # ---------------- inventario y orientación ----------------
    print("1/3 inventario")
    inv, prov = fetch_inventory(args.net, args.stations, t0, t1)
    inv.write(os.path.join(OUT, "inventory.xml"), format="STATIONXML")
    rows, rotate = orientation_report(inv)
    report["inventory_provider"] = prov
    report["channels"] = rows
    report["needs_rotation"] = bool(rotate)

    if args.box is None:
        lats=[r["lat"] for r in rows]; lons=[r["lon"] for r in rows]
        box=[min(lats)-args.box_margin_deg, max(lats)+args.box_margin_deg,
             min(lons)-args.box_margin_deg, max(lons)+args.box_margin_deg]
        print(f"   caja calculada desde el inventario (+{args.box_margin_deg} deg): "
              f"lat[{box[0]:.2f},{box[1]:.2f}] lon[{box[2]:.2f},{box[3]:.2f}]")
    else:
        box=[float(x) for x in args.box.split(",")]
    print(f"   proveedor {prov}; {len(rows)} canales")
    print(f"   {'STA':6s} {'loc':4s} {'cha':5s} {'sps':>6s} {'az':>7s} {'dip':>7s}")
    for r in rows[:40]:
        az = "n/a" if r["az"] is None else f"{r['az']:.1f}"
        dp = "n/a" if r["dip"] is None else f"{r['dip']:.1f}"
        print(f"   {r['sta']:6s} {r['loc']:4s} {r['cha']:5s} {r['sr']:6.1f} {az:>7s} {dp:>7s}")
    if rotate:
        print("   AVISO: hay canales 1/2 o azimuts no estandar -> se rotara a ZNE en la etapa 2")
    else:
        print("   canales ya orientados en ZNE")
    pedidas = set(args.stations.split(","))
    obtenidas = {r["sta"] for r in rows}
    faltantes = pedidas - obtenidas
    if faltantes:
        print(f"   AVISO: estaciones pedidas SIN metadatos en {prov}: {sorted(faltantes)}")
        print("   -> revisa si existen en esta red/periodo antes de esperar sus formas de onda.")

    # ---------------- catálogo ----------------
    print("2/3 catalogo")
    if args.cat_start is not None:
        ct0 = UTCDateTime(args.cat_start)
    else:
        ct0 = t0
    ct1 = ct0 + (args.cat_days if args.cat_days is not None else args.days) * 86400
    if args.cat_margin_days > 0:
        ct0 = ct0 - args.cat_margin_days * 86400
        ct1 = ct1 + args.cat_margin_days * 86400
        print(f"   margen de catalogo +/-{args.cat_margin_days:.0f} dias aplicado: "
              f"{ct0} .. {ct1}")

    age_days = (UTCDateTime(datetime.utcnow()) - ct0) / 86400.0
    if age_days < 60 and args.cat_start is None:
        print(f"   AVISO: la ventana de catalogo empieza hace solo {age_days:.0f} dias.")
        print("   Los catalogos revisados (USGS/EMSC completos, ISC) tardan meses en publicar")
        print("   sismicidad de magnitud baja. Si esto da 0 eventos, NO uses --start reciente")
        print("   para el catalogo: pasa --cat-start con una fecha de hace >=6 meses, mientras")
        print("   las formas de onda (--start) pueden seguir siendo recientes si el SGC las")
        print("   tiene archivadas (verificado que sí, via task1_fetch.py con el SGC).")

    cat, cprov = fetch_catalog(ct0, ct1, box, args.minmag)
    report["catalog_start"] = str(ct0)
    report["catalog_end"] = str(ct1)
    if cat is None:
        print("   NINGUN proveedor devolvio eventos. Amplia el rango o baja --minmag.",
              file=sys.stderr)
        report["catalog_provider"] = None
        report["n_events"] = 0
    else:
        cat.write(os.path.join(OUT, "catalog.xml"), format="QUAKEML")
        report["catalog_provider"] = cprov
        report["n_events"] = int(cat.count())
        print(f"   {cat.count()} eventos desde {cprov}")
        if cat.count() < 20:
            print("   AVISO: menos de 20 eventos. La etapa 3 necesita vecinos para medir")
            print("   coherencia; considera otro periodo o una caja mayor.")

    # ---------------- ondas continuas ----------------
    if args.days > 5 and not args.skip_size_warning:
        n_sta = len(args.stations.split(","))
        est_mb = 0.174 * n_sta * args.days      # calibrado: 348 MB / 6 estaciones / 2 dias
        import math as _m
        rate = 62 / 365.0                        # tasa observada USGS M>=4 para esta region
        lam = rate * args.days
        p1 = 1 - _m.exp(-lam)
        print(f"   AVISO: --days {args.days:.0f} con {n_sta} estaciones -> "
              f"~{est_mb:.0f} MB estimados de descarga.")
        print(f"   Con la tasa observada de USGS (M>=4) para esta region, "
              f"P(>=1 evento en la ventana) = {p1:.0%}.")
        print("   Usa --skip-size-warning para omitir este aviso en corridas futuras.")
    print("3/3 formas de onda continuas")
    domain = RectangularDomain(minlatitude=box[0], maxlatitude=box[1],
                               minlongitude=box[2], maxlongitude=box[3])
    restrictions = Restrictions(
        starttime=t0, endtime=t1,
        chunklength_in_sec=3600,
        network=args.net, station=args.stations,
        # location="00" fijo asumia que TODO proveedor archiva en ese location; el SGC lo
        # confirma para HH/BH, pero el respaldo EarthScope/IRIS podria usar otro codigo para
        # alguna estacion, y con el valor fijo esa estacion se perderia en silencio.
        location_priorities=["00", "", "10"],
        # BH restaurado: el SGC no lo archiva (verificado), pero EarthScope/IRIS como respaldo
        # si podrian tenerlo para alguna estacion que no opere en HH en ese periodo.
        channel_priorities=["HH[ZNE12]", "BH[ZNE12]", "EH[ZNE12]"],
        reject_channels_with_gaps=False,
        minimum_length=0.0,                  # el clutter vive en los trozos sucios
        minimum_interstation_distance_in_m=0,
        sanitize=False)
    # MassDownloader necesita instancias de Client, no solo nombres/URLs
    mdl_clients = [make_client(p) for p in WAVE_PROVIDERS]
    mdl = MassDownloader(providers=mdl_clients)
    print(f"   descargando con timeout=30s por solicitud, threads_per_client=2 "
          f"(reducido para no saturar el SGC, que es un servidor publico modesto)")
    mdl.download(domain, restrictions, mseed_storage=WAVE, stationxml_storage=STAX,
                threads_per_client=2, download_chunk_size_in_mb=5)

    files = [f for f in os.listdir(WAVE) if f.endswith(".mseed")]
    report["n_mseed"] = len(files)
    print(f"   {len(files)} archivos en {WAVE}")
    if not files:
        print("   SIN DATOS. Causas frecuentes: fecha futura, location code, o canal.",
              file=sys.stderr)

    with open(os.path.join(OUT, "fetch_report.json"), "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"\nreporte en {OUT}/fetch_report.json")
    print("criterio para pasar a la etapa 2: >=6 estaciones con 3 componentes y >90% de cobertura")


if __name__ == "__main__":
    main()
