"""
probe_sgc_events.py — verifica si el SGC expone fdsnws-event en el mismo servidor que
dataselect. Si lo hace, el catalogo tendria mucha menor magnitud de completitud que USGS
para esta region, lo que resolveria de raiz el problema de pocos eventos por ventana.

    python probe_sgc_events.py
"""
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

BASE = "http://sismo.sgc.gov.co:8080"

try:
    c = Client(base_url=BASE)
    print(f"[OK] cliente creado para {BASE}")
except Exception as e:
    print(f"[FALLA] {type(e).__name__}: {e}")
    raise SystemExit

t0 = UTCDateTime("2025-11-01")
t1 = t0 + 2 * 86400

try:
    cat = c.get_events(starttime=t0, endtime=t1, minmagnitude=1.0,
                       minlatitude=3.6, maxlatitude=10.9,
                       minlongitude=-76.6, maxlongitude=-72.0)
    print(f"[OK] fdsnws-event responde: {len(cat)} eventos M>=1.0 en 2 dias")
    for e in cat[:10]:
        o = e.preferred_origin() or e.origins[0]
        m = e.preferred_magnitude() or (e.magnitudes[0] if e.magnitudes else None)
        print(f"   {o.time}  M{m.mag if m else '?'}  {o.latitude:.2f},{o.longitude:.2f}  "
              f"prof={float(o.depth or 0)/1000:.1f}km")
except Exception as e:
    print(f"[FALLA] get_events: {type(e).__name__}: {e}")
    print("El SGC no expone fdsnws-event publico en este servidor, o el catalogo local no")
    print("esta disponible por esta via. Seguir usando USGS/EMSC/ISC como esta.")
