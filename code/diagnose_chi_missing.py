"""
diagnose_chi_missing.py — CHI no aparece NI UNA VEZ en picks_polarization.csv, pese a ser la
estacion mas cercana (47-73 km) a los 12 eventos del catalogo. Ni un solo pick, ni siquiera
asociado. Eso no es "poco clutter": es ausencia total, y hay que verificar cada eslabon de la
cadena para CHI en particular:

  1. tiene CHI archivos mseed reales en la ventana?
  2. inventario declarado para CHI
  3. las ventanas de evento+coda de CHI caen dentro del rango de datos disponible?
  4. si se le pasa CHI al picker de forma aislada, produce picks?

    python diagnose_chi_missing.py
"""
import glob
import json
import sys

sys.path.insert(0, ".")
import numpy as np
from obspy import read, UTCDateTime, read_events
from obspy.geodetics import gps2dist_azimuth

OUT = "t1_data"

print("=" * 70)
print("PASO 1: archivos mseed de CHI en el directorio de formas de onda")
print("=" * 70)
files = sorted(glob.glob(f"{OUT}/waveforms/*.mseed"))
chi_files = [f for f in files if "CHI" in f.upper()]
print(f"\narchivos totales: {len(files)}")
print(f"archivos con 'CHI' en el nombre: {len(chi_files)}")
if chi_files:
    print("primeros 5:", chi_files[:5])
    st_chi = read(chi_files[0])
    print("\nprimera traza de CHI:")
    print(st_chi)
else:
    print("\nNO HAY ARCHIVOS CON 'CHI' EN EL NOMBRE. La descarga nunca trajo CHI.")
    print("Revisa task1_fetch.py: la estacion pudo fallar silenciosamente en el downloader.")

print("\n" + "=" * 70)
print("PASO 2: inventario declarado para CHI")
print("=" * 70)
with open(f"{OUT}/fetch_report.json") as f:
    rep = json.load(f)
chi_channels = [r for r in rep["channels"] if r["sta"] == "CHI"]
print(f"\ncanales de CHI en fetch_report.json: {len(chi_channels)}")
for c in chi_channels:
    print(f"   {c}")

if not chi_files:
    print("\nDeteniendo aqui: sin archivos, los pasos 3 y 4 no tienen sentido.")
    sys.exit(0)

print("\n" + "=" * 70)
print("PASO 3: rango temporal REAL de los datos de CHI vs ventana de eventos")
print("=" * 70)
st_full = read(chi_files[0])
for f in chi_files[1:]:
    st_full += read(f)
print(f"\ntrazas leidas de CHI: {len(st_full)}")
for tr in st_full[:5]:
    print(f"   {tr.stats.channel}: {tr.stats.starttime} -> {tr.stats.endtime} "
          f"({tr.stats.npts} muestras)")

cat = read_events(f"{OUT}/catalog.xml")
chi_coord = (chi_channels[0]["lat"], chi_channels[0]["lon"]) if chi_channels else None
if chi_coord:
    print(f"\nCHI en ({chi_coord[0]:.3f}, {chi_coord[1]:.3f}). Llegadas previstas (vp=7.5):")
    for e in cat:
        o = e.preferred_origin() or e.origins[0]
        d_m, _, _ = gps2dist_azimuth(float(o.latitude), float(o.longitude), *chi_coord)
        d_km = float(np.hypot(d_m / 1000.0, float(o.depth or 0) / 1000.0))
        t_p = o.time + d_km / 7.5
        cubierto = any(tr.stats.starttime <= t_p <= tr.stats.endtime for tr in st_full)
        print(f"   {str(o.time)[:19]}  d={d_km:.0f}km  P_previsto={str(t_p)[:19]}  "
              f"cubierto_por_algun_archivo={cubierto}")

print("\n" + "=" * 70)
print("PASO 4: picker aislado sobre CHI solamente")
print("=" * 70)
try:
    import seisbench.models as sbm
    st_test = st_full.copy()
    st_test.merge(method=1)
    st_test = st_test.split()
    st_test.detrend("demean").detrend("linear")
    st_test.filter("bandpass", freqmin=1.0, freqmax=3.0, corners=4, zerophase=True)
    model = sbm.PhaseNet.from_pretrained("instance")
    out = model.classify(st_test, batch_size=256, P_threshold=0.05, S_threshold=0.05)
    print(f"\npicks de PhaseNet SOLO sobre CHI, umbral 0.05: {len(out.picks)}")
    for pk in out.picks[:10]:
        print(f"   {pk.trace_id}  {pk.phase}  {pk.peak_time}  prob={pk.peak_value:.3f}")
    if not out.picks:
        print("\nCERO PICKS incluso corriendo el picker aislado sobre el continuo COMPLETO")
        print("de CHI (no solo las ventanas recortadas). Esto descarta un problema del")
        print("recorte por eventos: el picker no encuentra nada en los datos de CHI, punto.")
        print("Revisar: nivel de ruido/ganancia de CHI, o si el canal esta realmente activo")
        print("en este periodo (podria estar en mantenimiento pese a tener metadatos).")
except Exception as ex:
    print(f"\nerror al picar CHI aislado: {type(ex).__name__}: {ex}")
