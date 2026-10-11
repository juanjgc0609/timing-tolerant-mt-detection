"""
check_chi_june.py — verifica si CHI tiene archivos de JUNIO 2025 en absoluto, o si el
downloader nunca los trajo para esta estacion en esta ventana especifica.

Confirmado por diagnose_chi_missing.py: los archivos leidos para CHI son de
2025-09-30/10-01, residuo de la corrida de 30 dias en octubre, NO de la ventana de junio que
se pidio despues. Esto se investiga aqui.

    python check_chi_june.py
"""
import glob
import re
from datetime import datetime

files = sorted(glob.glob("t1_data/waveforms/*CHI*.mseed"))
print(f"archivos CHI totales: {len(files)}\n")

# extraer fecha de inicio del nombre de archivo: patron __YYYYMMDDTHHMMSSZ__
dates = []
for f in files:
    m = re.search(r"__(\d{8})T\d{6}Z__", f)
    if m:
        dates.append(m.group(1))

from collections import Counter
by_month = Counter(d[:6] for d in dates)
print("archivos CHI por mes (YYYYMM):")
for ym, n in sorted(by_month.items()):
    print(f"   {ym}: {n} archivos")

junio = [f for f in files if "__202506" in f]
print(f"\narchivos CHI que empiezan en junio 2025: {len(junio)}")
if junio:
    print("primeros 5:", junio[:5])
else:
    print("NINGUNO. CHI no tiene ni un solo archivo de junio 2025.")
    print("La descarga de --start 2025-06-08 --days 3 nunca trajo CHI, pese a que")
    print("fetch_report.json y el inventario si lo declaran (verificado en diagnose_chi_missing).")
    print()
    print("Causa mas probable: task1_fetch.py NO limpia t1_data/waveforms/ entre corridas.")
    print("Si la descarga de junio para CHI fallo silenciosamente (por ejemplo timeout del")
    print("SGC para esa estacion en particular) mientras las otras 5 SI se descargaron, el")
    print("directorio quedo con una mezcla: 5 estaciones de junio + CHI residual de octubre.")
    print("Y como load_streams() filtra por VENTANA TEMPORAL pero NO verifica que cada")
    print("estacion tenga cobertura, dejo pasar los archivos de octubre de CHI creyendo que")
    print("son datos validos para analizar (no lo son: caen fuera de [2025-06-08,2025-06-11]).")
