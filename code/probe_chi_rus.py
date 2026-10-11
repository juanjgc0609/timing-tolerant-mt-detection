"""
probe_chi_rus.py — sondea si CHI y RUS tienen ALGUNA vez forma de onda HH archivada en el
SGC, en vez de asumir que el hueco es solo del 8-11 de junio de 2025.

Prueba ventanas cortas (60s) en varias fechas repartidas a lo largo de 2024-2025, y tambien
unos dias inmediatamente ANTES y DESPUES del 8 de junio de 2025 (por si el hueco es de solo
unos dias, no del mes completo).

    python probe_chi_rus.py
"""
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

SGC = Client(base_url="http://sismo.sgc.gov.co:8080", timeout=20)
STAS = ["CHI", "RUS", "OCA"]     # OCA de control: sabemos que SI tiene dato en este periodo

# barrido amplio: primer dia de cada mes en 2024 y 2025
dates_wide = [UTCDateTime(f"{y}-{m:02d}-01") for y in (2024, 2025) for m in range(1, 13)
              if UTCDateTime(f"{y}-{m:02d}-01") < UTCDateTime.now()]

# barrido fino: +/- 10 dias alrededor del 8 de junio de 2025
center = UTCDateTime("2025-06-08")
dates_fine = [center + d * 86400 for d in range(-10, 11)]

def probe(sta, t0, dur=60):
    try:
        st = SGC.get_waveforms("CM", sta, "00", "HHZ", t0, t0 + dur, attach_response=False)
        return len(st) > 0 and any(tr.stats.npts > 0 for tr in st)
    except Exception:
        return False

print("=" * 68)
print("BARRIDO FINO: +/- 10 dias alrededor de 2025-06-08 (HHZ, loc=00)")
print("=" * 68)
print(f"{'fecha':12s}  " + "  ".join(f"{s:>5s}" for s in STAS))
for t0 in dates_fine:
    row = [probe(s, t0) for s in STAS]
    marks = "  ".join(f"{'SI':>5s}" if ok else f"{'no':>5s}" for ok in row)
    print(f"{str(t0.date):12s}  {marks}")

print("\n" + "=" * 68)
print("BARRIDO AMPLIO: primer dia de cada mes, 2024-2025 (HHZ, loc=00)")
print("=" * 68)
print(f"{'fecha':12s}  " + "  ".join(f"{s:>5s}" for s in STAS))
for t0 in dates_wide:
    row = [probe(s, t0) for s in STAS]
    marks = "  ".join(f"{'SI':>5s}" if ok else f"{'no':>5s}" for ok in row)
    print(f"{str(t0.date):12s}  {marks}")

print("\nSi CHI/RUS dan 'no' en TODAS las fechas de ambos barridos: la estacion probablemente")
print("no tiene HH continuo archivado en el SGC en ningun periodo reciente (sensor caido,")
print("desmantelado, o solo transmite en tiempo real sin archivar). Si dan 'SI' en algunas")
print("fechas pero no en 2025-06-08/11: es un hueco puntual de esos dias especificos, y otro")
print("periodo con 'SI' para CHI y/o RUS serviria igual para la Tarea 1 (no hace falta que sea")
print("el mismo swarm; cualquier ventana con eventos catalogados cerca de una estacion con dato")
print("real sirve).")
