"""
probe_waveforms.py — prueba puntual de disponibilidad, SIN MassDownloader.

Objetivo: descartar si el problema es el downloader masivo o la disponibilidad real del dato en
el proveedor. Prueba doce combinaciones de proveedor/location/canal para una estación y una hora
fijas.

    python probe_waveforms.py
    python probe_waveforms.py --net CM --sta OCA --start 2024-06-01
"""
import argparse
from obspy import UTCDateTime
from obspy.clients.fdsn import Client


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", default="CM")
    ap.add_argument("--sta", default="OCA")
    ap.add_argument("--start", default="2024-06-01T00:00:00")
    ap.add_argument("--hours", type=float, default=1.0)
    args = ap.parse_args()

    t0 = UTCDateTime(args.start)
    t1 = t0 + args.hours * 3600

    print(f"probando {args.net}.{args.sta} desde {t0} hasta {t1}\n")
    print(f"{'proveedor':>10} {'loc':>6} {'cha':>5} | resultado")
    print("-" * 50)

    any_ok = False
    for prov in ["EARTHSCOPE", "IRIS"]:
        try:
            c = Client(prov)
        except Exception as e:
            print(f"{prov:>10} {'--':>6} {'--':>5} | no se pudo crear cliente: {e}")
            continue
        for loc in ["00", "", "10"]:
            for cha in ["HHZ", "BHZ"]:
                try:
                    st = c.get_waveforms(args.net, args.sta, loc, cha, t0, t1)
                    n = len(st)
                    npts = st[0].stats.npts if n else 0
                    print(f"{prov:>10} {loc!r:>6} {cha:>5} | OK, {n} trazas, npts={npts}")
                    any_ok = True
                except Exception as e:
                    print(f"{prov:>10} {loc!r:>6} {cha:>5} | FALLA -> {type(e).__name__}")

    print()
    if any_ok:
        print("Al menos una combinacion funciono: el downloader masivo tenia un problema de")
        print("configuracion (location/canal), no falta de dato. Corrige esos parametros en")
        print("task1_fetch.py.")
    else:
        print("TODAS fallaron: no hay continuo archivado para esta estacion/fecha en estos")
        print("proveedores. Prueba con otra fecha (--start) o pasa a descargar desde el SGC")
        print("directamente.")


if __name__ == "__main__":
    main()
