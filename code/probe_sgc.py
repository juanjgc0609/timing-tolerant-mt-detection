"""
probe_sgc.py — prueba el FDSN del SGC directamente (Servicio Geologico Colombiano).

El SGC opera la red CM y expone su propio servicio FDSN, independiente de EarthScope/IRIS:
    http://sismo.sgc.gov.co:8080/fdsnws/

Fuente: Poveda et al. 2022, Earth and Space Science, cita el mismo endpoint para descargar
miniSEED de la RSNC. No verificado en vivo desde aqui (sin acceso de red); este script lo hace
en tu maquina.

    python probe_sgc.py
    python probe_sgc.py --start 2024-06-01 --sta OCA
"""
import argparse
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

BASE = "http://sismo.sgc.gov.co:8080"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", default="CM")
    ap.add_argument("--sta", default="OCA")
    ap.add_argument("--start", default="2024-06-01T00:00:00")
    ap.add_argument("--hours", type=float, default=1.0)
    args = ap.parse_args()

    t0 = UTCDateTime(args.start)
    t1 = t0 + args.hours * 3600

    print(f"probando el FDSN del SGC en {BASE}\n")

    # paso 1: cliente
    try:
        c = Client(base_url=BASE)
        print("[OK] cliente creado")
    except Exception as e:
        print(f"[FALLA] no se pudo crear el cliente: {type(e).__name__}: {e}")
        print("El servidor puede estar caido, o el puerto 8080 bloqueado por tu red/firewall.")
        return

    # paso 2: metadatos de estacion (mas ligero, prueba conectividad basica)
    try:
        inv = c.get_stations(network=args.net, station=args.sta, level="channel", channel="?H?")
        n = len(inv.get_contents()["channels"])
        print(f"[OK] inventario: {n} canales para {args.net}.{args.sta}")
    except Exception as e:
        print(f"[FALLA] get_stations: {type(e).__name__}: {e}")
        print("Si esto falla, el servidor no responde o la red/estacion no existe alli.")
        return

    # paso 3: formas de onda, la prueba real
    print(f"\nprobando forma de onda {args.net}.{args.sta} {t0} -> {t1}")
    for loc in ["00", "", "10"]:
        for cha in ["HHZ", "BHZ"]:
            try:
                st = c.get_waveforms(args.net, args.sta, loc, cha, t0, t1)
                print(f"  loc={loc!r} cha={cha}: OK, {len(st)} trazas, npts={st[0].stats.npts}")
            except Exception as e:
                print(f"  loc={loc!r} cha={cha}: FALLA -> {type(e).__name__}")

    print("\nSi al menos una combinacion dio OK, el SGC SI tiene el continuo, y basta con")
    print("cambiar el proveedor en task1_fetch.py de EARTHSCOPE/IRIS a este base_url.")
    print("Si TODAS fallan aqui tambien, hay que escribir directamente al SGC")
    print("(radicacioncorrespondencia@sgc.gov.co) o usar el portal datos.sgc.gov.co.")


if __name__ == "__main__":
    main()
