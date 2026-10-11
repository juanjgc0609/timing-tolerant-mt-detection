"""
probe_sgc2.py — inspecciona el inventario COMPLETO del SGC para CM.OCA y barre location/canal
reales en varias fechas, incluyendo hoy. Corregido: bug de tipos en la lista de fechas.

    python probe_sgc2.py
"""
from datetime import datetime, timedelta
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

BASE = "http://sismo.sgc.gov.co:8080"
c = Client(base_url=BASE)

print("1) inventario completo, todos los canales\n")
inv = c.get_stations(network="CM", station="OCA", level="channel")
combos = set()
for net in inv:
    for sta in net:
        for ch in sta:
            loc = ch.location_code
            print(f"   loc={loc!r:6s} cha={ch.code:5s} sr={ch.sample_rate:7.2f} "
                  f"{ch.start_date.date} -> {ch.end_date.date if ch.end_date else 'presente'}")
            combos.add((loc, ch.code[:2]))

print(f"\n2) combinaciones (location, banda) encontradas: {sorted(combos)}\n")

print("3) barrido de forma de onda, combinaciones reales, varias fechas\n")
now_dt = datetime.utcnow()
offsets_days = [0, 1, 3, 7, 15, 30, 60, 90, 180, 365]
dates = [UTCDateTime((now_dt - timedelta(days=d)).date()) for d in offsets_days]

any_ok = False
for t0 in dates:
    t1 = t0 + 3600
    for loc, band in sorted(combos):
        if band not in ("HH", "BH"):     # los de banda ancha/corta son los utiles para 1-3 Hz
            continue
        cha = band + "Z"
        try:
            st = c.get_waveforms("CM", "OCA", loc, cha, t0, t1)
            print(f"   {t0.date}  loc={loc!r} cha={cha}: OK, {len(st)} trazas, "
                  f"npts={st[0].stats.npts}")
            any_ok = True
        except Exception as e:
            print(f"   {t0.date}  loc={loc!r} cha={cha}: FALLA -> {type(e).__name__}")

print()
if any_ok:
    print("Al menos una fecha/combinacion funciono. Usa esa fecha y esos codigos en task1_fetch.py.")
else:
    print("Nada funciono en ninguna fecha probada, pese a que el canal HH/00 esta 'operando'")
    print("segun el inventario desde 2022-12-08. Eso indica que el fdsnws-dataselect publico")
    print("no sirve el archivo completo, solo metadatos, o hay latencia/restriccion de acceso.")
    print("Siguiente paso: escribir a radicacioncorrespondencia@sgc.gov.co, o revisar si")
    print("datos.sgc.gov.co ofrece descarga directa de miniSEED por otra via.")
