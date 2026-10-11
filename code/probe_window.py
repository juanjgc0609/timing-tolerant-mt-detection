"""
probe_window.py — sondeo DIARIO de disponibilidad HH para una ventana especifica, antes de
gastar una descarga completa en ella.

Generaliza el chequeo puntual que ya se hizo a mano (probe_chi_rus.py) a cualquier ventana y
lista de estaciones, porque un sondeo MENSUAL (1 muestra por mes) puede pasar por alto un hueco
corto: visto en la practica, OCA cayo del 1 al 4 de junio de 2025 y esos 4 dias no aparecian en
un muestreo de "primer dia del mes".

    python probe_window.py --start 2025-02-23 --days 30 --stations CHI,RUS
"""
import argparse
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

SGC = Client(base_url="http://sismo.sgc.gov.co:8080", timeout=20)


def probe(sta, t0, dur=60):
    try:
        st = SGC.get_waveforms("CM", sta, "00", "HHZ", t0, t0 + dur, attach_response=False)
        return len(st) > 0 and any(tr.stats.npts > 0 for tr in st)
    except Exception:
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", required=True)
    ap.add_argument("--days", type=int, default=30)
    ap.add_argument("--stations", required=True, help="ej. CHI,RUS")
    args = ap.parse_args()

    stas = args.stations.split(",")
    t0 = UTCDateTime(args.start)
    print(f"{'fecha':12s}  " + "  ".join(f"{s:>5s}" for s in stas))
    huecos = {s: 0 for s in stas}
    for d in range(args.days):
        day = t0 + d * 86400
        row = [probe(s, day) for s in stas]
        for s, ok in zip(stas, row):
            if not ok:
                huecos[s] += 1
        marks = "  ".join(f"{'SI':>5s}" if ok else f"{'no':>5s}" for ok in row)
        print(f"{str(day.date):12s}  {marks}")

    print(f"\ndias sin dato en la ventana de {args.days}: "
          + ", ".join(f"{s}={huecos[s]}" for s in stas))
    if all(huecos[s] == 0 for s in stas):
        print("SIN HUECOS: la ventana esta limpia para todas las estaciones pedidas, se puede")
        print("descargar con confianza.")
    else:
        print("HAY HUECOS: antes de descargar, decide si acortas la ventana para evitarlos o")
        print("si aceptas la cobertura parcial (task1_fetch.py ya avisa por estacion al final).")


if __name__ == "__main__":
    main()
