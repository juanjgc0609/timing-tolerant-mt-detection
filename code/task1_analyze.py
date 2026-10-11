"""
task1_analyze.py — Etapa 2 de la Tarea 1: picks, polarización y etiquetado.

Para cada pick declarado por el picker calcula la matriz de covarianza 3C, la rectilinearidad y
el back-azimut medido, y lo etiqueta como ASOCIADO a un evento del catálogo o NO ASOCIADO. Los no
asociados son el clutter cuya coherencia de polarización es el objeto de la Tarea 1.

    python task1_analyze.py --band 1.0 3.0 --win 0.6

Salida: t1_data/picks_polarization.csv con una fila por pick.

Detalles que deciden la validez de la medición y por eso están explícitos:
  - rotación a ZNE con el inventario si los canales no estan orientados
  - banda de trabajo igual a la del generador sintetico (2 Hz centrales)
  - ventana de covarianza igual a la del generador (0.6 s)
  - el back-azimut de polarizacion tiene ambiguedad de 180 grados porque el signo de un autovector
    es arbitrario; los angulos se pliegan a [0, 90] y eso se declara
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from obspy import read, read_inventory, read_events, UTCDateTime, Stream
from obspy.geodetics import gps2dist_azimuth

OUT = "t1_data"
VP_OVER_VS = 1.73
V_P_KMS = 6.0          # crustal, para sismicidad LOCAL SOMERA; ver --vp para otros regimenes
V_S_KMS = V_P_KMS / VP_OVER_VS


def load_streams(inv, band, rotate, t0=None, t1=None):
    """
    Carga SOLO los mseed cuya ventana temporal se solapa con [t0, t1]. Sin este filtro, cualquier
    archivo residual de una corrida anterior con fecha distinta hace que merge() intente
    construir un array continuo cubriendo la brecha completa entre ambas fechas: con una brecha
    de meses a 100 Hz eso es miles de millones de muestras y agota la memoria (visto: intento de
    reservar 6.9 GiB para 1.85e9 muestras, exactamente lo esperable para una brecha de ~7 meses).

    ROTA primero (rotate() si tolera pequenos huecos internamente es fragil; aqui se opta por
    partir en segmentos SIN huecos antes de cualquier operacion que no admita mascaras).
    merge(method=1) sin fill_value produce trazas ENMASCARADAS en los huecos, y
    demean/detrend/filter no soportan mascaras (NotImplementedError). split() separa la
    corriente en tramos continuos reales; se pierde continuidad exacta en los bordes de hueco,
    pero eso es preferible a rellenar con un valor artificial que contaminaria la polarizacion.
    """
    files = sorted(glob.glob(os.path.join(OUT, "waveforms", "*.mseed")))
    if not files:
        sys.exit("no hay mseed; corre task1_fetch.py primero")

    if t0 is not None and t1 is not None:
        from obspy.io.mseed.util import get_start_and_end_time
        kept, dropped = [], 0
        for f in files:
            try:
                fs, fe = get_start_and_end_time(f)
            except Exception:
                kept.append(f); continue          # si no se puede leer la cabecera, se intenta igual
            if fe < t0 or fs > t1:                 # sin solape alguno con la ventana pedida
                dropped += 1
                continue
            kept.append(f)
        if dropped:
            print(f"   {dropped} archivo(s) fuera de la ventana [{t0}, {t1}] descartados "
                  f"(residuo de una corrida anterior con otra fecha)")
        files = kept
        if not files:
            sys.exit(f"ningun mseed se solapa con la ventana [{t0}, {t1}]. Revisa --start/--days "
                     f"o borra archivos residuales de otras fechas en t1_data/waveforms.")

    st = read(files[0])
    for f in files[1:]:
        st += read(f)
    if t0 is not None and t1 is not None:
        st.trim(starttime=t0, endtime=t1)          # recorte duro antes de merge, por seguridad
    st.merge(method=1)             # sin fill_value: produce mascara donde hay huecos reales
    st = st.split()                # separa en tramos SIN mascara, listos para demean/filter
    n_gaps = len(st.get_gaps())
    if n_gaps:
        print(f"   AVISO: {n_gaps} huecos detectados en el continuo; se proceso por tramos")
    st.detrend("demean").detrend("linear")
    if rotate:
        print("   rotando a ZNE con el inventario")
        st.rotate("->ZNE", inventory=inv)
    st.filter("bandpass", freqmin=band[0], freqmax=band[1], corners=4, zerophase=True)
    st.merge(method=1)              # recompone para que polarization_at encuentre ventanas largas
    return st


def run_picker(st, p_thr, s_thr):
    """
    seisbench/PhaseNet no esta garantizado a soportar arrays enmascarados (huecos reales en el
    continuo). Se separa cualquier resto enmascarado en tramos limpios ANTES de picking, en vez
    de arriesgar un fallo silencioso o una excepcion a mitad de la etapa.
    """
    import seisbench.models as sbm
    if any(np.ma.is_masked(tr.data) for tr in st):
        st = st.split()
        print("   AVISO: quedaban tramos con huecos antes del picker; se separaron limpios")
    model = sbm.PhaseNet.from_pretrained("instance")
    out = model.classify(st, batch_size=256, P_threshold=p_thr, S_threshold=s_thr)
    rows = []
    for pk in out.picks:
        parts = pk.trace_id.split(".")
        if len(parts) < 2:
            print(f"   AVISO: trace_id inesperado '{pk.trace_id}', se descarta el pick")
            continue
        net, sta = parts[0], parts[1]
        rows.append(dict(net=net, sta=sta, phase=pk.phase,
                         t=UTCDateTime(pk.peak_time), prob=float(pk.peak_value)))
    return pd.DataFrame(rows)


def polarization_at(st, sta, t, win):
    """
    Covarianza 3C en una ventana centrada en t. Devuelve (rect, backaz_deg, inc_deg, n).

    Si hubo huecos reales en el continuo, load_streams() habra dejado, para un mismo canal,
    varios tramos disjuntos (Trace objects) en vez de uno solo. select(station=sta) devuelve
    TODOS esos tramos; hay que elegir, por canal, el tramo que realmente cubre [t-half, t+half],
    no el ultimo iterado (eso silenciosamente mediria la polarizacion de datos equivocados).

    Ademas se agrupa por PREFIJO de banda de instrumento (los dos primeros caracteres del
    codigo de canal, p.ej. HH, BH, EH), no solo por la componente final: si una estacion tuvo
    un cambio de instrumentacion dentro del periodo y coexisten trazas HH y EH, comp[k]=tr
    indexado solo por 'E'/'N'/'Z' mezclaria por ejemplo EHE con HHN, una covarianza sin
    sentido fisico. Solo se acepta un grupo cuyas tres componentes compartan banda Y tasa de
    muestreo.
    """
    half = win / 2.0
    sel = st.select(station=sta)
    bands = {}
    for tr in sel:
        if tr.stats.starttime <= t - half and tr.stats.endtime >= t + half:
            prefix = tr.stats.channel[:2]
            comp_code = tr.stats.channel[-1]
            bands.setdefault(prefix, {})[comp_code] = tr
    comp = None
    for prefix, group in bands.items():
        if all(k in group for k in ("Z", "N", "E")):
            sps = {group[k].stats.sampling_rate for k in ("Z", "N", "E")}
            if len(sps) == 1:
                comp = group
                break
    if comp is None:
        return None
    segs = []
    for k in ("E", "N", "Z"):                       # orden E,N,Z -> x,y,z
        tr = comp[k].slice(t - half, t + half)
        if tr.stats.npts < 8 or np.ma.is_masked(tr.data):
            return None
        # tr.data.astype(float) sobre un MaskedArray SIGUE siendo MaskedArray aunque no tenga
        # ningun valor enmascarado (is_masked() da False igual); X @ X.T con esa subclase
        # dispara la propagacion de mascaras de numpy.ma y falla al mezclar formas (3,n)/(n,3).
        # np.asarray() fuerza un ndarray plano de verdad, sin depender de is_masked().
        segs.append(np.asarray(tr.data, dtype=float))
    n = min(len(s) for s in segs)
    X = np.vstack([s[:n] for s in segs])
    X = np.asarray(X, dtype=float)      # doble seguro: X tambien plano antes del producto
    X = X - X.mean(axis=1, keepdims=True)
    C = X @ X.T
    w, V = np.linalg.eigh(C)
    lam = w[::-1]
    v1 = V[:, -1]
    if lam[0] <= 0:
        return None
    rect = float(1.0 - lam[1] / lam[0])
    h = v1[:2]
    nh = float(np.linalg.norm(h))
    if nh < 1e-12:
        return None
    # azimut de la direccion principal horizontal, medido desde el norte
    az = float((np.degrees(np.arctan2(h[0], h[1]))) % 180.0)      # plegado a [0,180)
    inc = float(np.degrees(np.arctan2(nh, abs(v1[2]))))
    return rect, az, inc, int(n)


def label_association(picks, cat, sta_coords, tol_s, max_reasonable_km=120.0, vp_kms=None):
    """
    max_reasonable_km: si la distancia al evento MAS CERCANO del catalogo excede esto para casi
    todos los picks, el catalogo no tiene cobertura temporal/espacial real de este periodo (viejo
    catalogo revisado con rezago, o red de eventos que no cubre la zona), y ningun pick se
    asociara nunca por bueno que sea el estimador de polarizacion. Eso se detecta y advierte
    aqui en lugar de dejar que produzca un control en n=0 sin explicacion.
    """
    from obspy import UTCDateTime as _U
    if not picks.empty and not isinstance(picks["t"].iloc[0], _U):
        picks = picks.copy()
        picks["t"] = picks["t"].apply(lambda x: _U(x))
    """
    Un pick es ASOCIADO si cae dentro de la ventana [t_orig + d/VP - tol, t_orig + d/VS + tol]
    de algun evento del catalogo en esa estacion. Los demas son clutter.
    Devuelve columnas: associated, ev_index, ev_az (azimut real estacion->evento), ev_dist_km.
    """
    ev = []
    for k, e in enumerate(cat):
        o = e.preferred_origin() or e.origins[0]
        m = (e.preferred_magnitude() or (e.magnitudes[0] if e.magnitudes else None))
        ev.append(dict(k=k, t=o.time, lat=float(o.latitude), lon=float(o.longitude),
                       dep=float(o.depth or 0) / 1000.0,
                       mag=float(m.mag) if m else np.nan))
    vp = vp_kms if vp_kms is not None else V_P_KMS
    vs = vp / VP_OVER_VS
    assoc, evi, evaz, evd, evdt = [], [], [], [], []
    for _, r in picks.iterrows():
        if r["sta"] not in sta_coords:
            assoc.append(False); evi.append(-1); evaz.append(np.nan); evd.append(np.nan); evdt.append(np.nan); continue
        slat, slon = sta_coords[r["sta"]]
        matches = []                          # TODOS los eventos cuya ventana cubre el pick
        closest = (False, -1, np.nan, np.nan, 1e18)
        for e in ev:
            d_m, az_ev_to_sta, az_sta_to_ev = gps2dist_azimuth(e["lat"], e["lon"], slat, slon)
            d_km = float(np.hypot(d_m / 1000.0, e["dep"]))
            t_p = e["t"] + d_km / vp
            t_s = e["t"] + d_km / vs
            lo = t_p - tol_s
            hi = t_s + tol_s
            dt = abs(float(r["t"] - e["t"]))
            if lo <= r["t"] <= hi:
                # residuo respecto de la llegada P prevista: criterio de desempate cuando hay
                # mas de un evento cuya ventana cubre el pick (secuencias de replicas densas).
                # Un break en el primer solapamiento asignaria arbitrariamente el indice mas
                # bajo del catalogo, sin importar si otro evento explica el pick con mucha
                # mayor precision.
                p_resid = abs(float(r["t"] - t_p))
                matches.append((p_resid, e["k"], az_sta_to_ev, d_km, dt))
            if dt < closest[4]:
                closest = (False, e["k"], az_sta_to_ev, d_km, dt)
        if matches:
            matches.sort(key=lambda m: m[0])
            _, ek, eaz, edk, edt = matches[0]
            best = (True, ek, eaz, edk, edt)
        else:
            best = closest
        assoc.append(best[0]); evi.append(best[1]); evaz.append(best[2]); evd.append(best[3]); evdt.append(best[4])
    picks = picks.copy()
    picks["associated"] = assoc
    picks["ev_index"] = evi
    picks["ev_azimuth"] = evaz
    picks["ev_dist_km"] = evd
    picks["best_dt_s"] = evdt
    near = picks["ev_dist_km"] <= max_reasonable_km
    if near.sum() == 0 and len(picks):
        import sys as _sys
        print(f"   AVISO CRITICO: NINGUN pick tiene un evento del catalogo a menos de "
              f"{max_reasonable_km:.0f} km. El catalogo descargado no cubre este periodo/zona "
              f"con la densidad necesaria (rezago de catalogo revisado o area mal elegida). "
              f"El control de la Etapa 3 fallara por esta razon, no por la polarizacion.",
              file=_sys.stderr)
    return picks


def slice_around_events(st, cat, sta_coords, vp_kms, pre_s=60.0, post_s=240.0,
                        coda_start_s=300.0, coda_win_s=300.0, coda_step_s=300.0,
                        noise_per_event=3):
    """
    Recorta el stream a ventanas alrededor de la llegada P prevista de cada evento en cada
    estacion, MAS ventanas de CODA del mismo evento para capturar clutter.

    Motivacion del recorte: PhaseNet sobre 30 dias continuos genera decenas de miles de
    ventanas de inferencia (visto: 11475 para solo 2 dias) y tarda horas en CPU.

    Motivacion de las ventanas de coda: picar SOLO alrededor de la llegada deja sin clutter
    (visto: 30 asociados, 0 no asociados). Y usar ruido desplazado una hora tampoco sirve: los
    picks resultantes quedaban a 237-876 km del evento catalogado mas cercano, asi que su
    coherencia se medía contra un sismo lejano cuyo azimut es casi identico desde toda la red
    (artefacto geometrico). La coda del MISMO evento es el clutter correcto: el vecino es el
    propio evento, a distancia conocida y pequena, que es justo el escenario del Apendice T.
    """
    vs = vp_kms / VP_OVER_VS
    out = Stream()
    n_ev, n_noise = 0, 0
    for e in cat:
        o = e.preferred_origin() or e.origins[0]
        for sta, (slat, slon) in sta_coords.items():
            d_m, _, _ = gps2dist_azimuth(float(o.latitude), float(o.longitude), slat, slon)
            d_km = float(np.hypot(d_m / 1000.0, float(o.depth or 0) / 1000.0))
            t_p = o.time + d_km / vp_kms
            t_s = o.time + d_km / vs
            w0 = t_p - pre_s
            w1 = max(t_s, t_p) + post_s
            seg = st.slice(starttime=w0, endtime=w1).select(station=sta)
            if len(seg):
                out += seg
                n_ev += 1
            # Ventanas de CODA del propio evento, no desplazadas una hora.
            #
            # Con ruido desplazado +1h, el evento catalogado mas cercano a esos picks estaba a
            # 237-876 km (medido), asi que angle_dev_deg se calculaba contra un sismo lejano.
            # Un evento a 800 km tiene casi el mismo azimut desde toda la red, de modo que la
            # "coherencia" resultante seria un artefacto geometrico y no el efecto que la
            # Tarea 1 busca.
            #
            # La coda del MISMO evento es el clutter correcto: el vecino es el propio evento,
            # a distancia conocida y pequena, que es exactamente el escenario del Apendice T
            # (clutter coherente con un evento cercano).
            for j in range(noise_per_event):
                c0 = t_s + coda_start_s + j * coda_step_s
                c1 = c0 + coda_win_s
                nseg = st.slice(starttime=c0, endtime=c1).select(station=sta)
                if len(nseg):
                    out += nseg
                    n_noise += 1
    print(f"   {n_ev} ventanas evento-estacion (-{pre_s:.0f}s / +{post_s:.0f}s) "
          f"+ {n_noise} ventanas de CODA (desde S+{coda_start_s:.0f}s, "
          f"{coda_win_s:.0f}s cada una)")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", nargs=2, type=float, default=[1.0, 3.0])
    ap.add_argument("--win", type=float, default=0.6)
    ap.add_argument("--p-thr", type=float, default=0.30)
    ap.add_argument("--s-thr", type=float, default=0.30)
    ap.add_argument("--assoc-tol", type=float, default=5.0,
                    help="holgura en segundos de la ventana de asociacion. Para sismicidad "
                         "SOMERA local, 5 s basta. Para sismicidad de PROFUNDIDAD INTERMEDIA "
                         "(p.ej. Nido de Bucaramanga, 75-160 km), la aproximacion de linea "
                         "recta con velocidad constante es cruda: usa --assoc-tol 15 o mas.")
    ap.add_argument("--vp", type=float, default=None,
                    help="velocidad P (km/s) para la ventana de asociacion. Por defecto 6.0 "
                         "(cortical, sismicidad somera). Para eventos de profundidad "
                         "intermedia (~75-160 km, ej. Nido de Bucaramanga) usa un valor mas "
                         "alto, ~7.5-8.0 km/s, que aproxima mejor una trayectoria Pn/Sn a "
                         "traves de corteza y manto superior.")
    ap.add_argument("--start", default=None,
                    help="ventana de forma de onda a analizar, si difiere de la ultima corrida "
                         "de task1_fetch.py. Por defecto se lee de fetch_report.json. Fijar "
                         "esto explicitamente evita que residuos .mseed de otras fechas en "
                         "t1_data/waveforms se mezclen y agoten la memoria al fusionar.")
    ap.add_argument("--days", type=float, default=None,
                    help="duracion de la ventana si --start se da explicitamente")
    ap.add_argument("--around-events", dest="around_events", action="store_true", default=True,
                    help="ACTIVO POR DEFECTO. Pica solo ventanas alrededor de las llegadas "
                         "previstas de los eventos del catalogo, en vez del continuo completo. "
                         "PhaseNet sobre 30 dias continuos genera decenas de miles de ventanas "
                         "de inferencia y tarda horas en CPU; este modo mide exactamente lo "
                         "mismo (control con picks de eventos + clutter de falsos positivos) "
                         "con dos ordenes de magnitud menos de computo.")
    ap.add_argument("--full-continuous", dest="around_events", action="store_false",
                    help="desactiva --around-events y pica el continuo COMPLETO. Solo si "
                         "necesitas clutter lejos de cualquier evento catalogado; cuenta con "
                         "horas de CPU.")
    ap.add_argument("--pre-s", type=float, default=60.0,
                    help="segundos antes de la llegada P prevista en cada ventana")
    ap.add_argument("--post-s", type=float, default=240.0,
                    help="segundos despues de la llegada S prevista, para cubrir coda")
    ap.add_argument("--coda-start-s", type=float, default=300.0,
                    help="segundos tras la S prevista donde empieza la primera ventana de "
                         "coda. Debe superar --assoc-tol con margen para que los picks de "
                         "coda se etiqueten como clutter y no como asociados.")
    ap.add_argument("--coda-win-s", type=float, default=300.0,
                    help="duracion de cada ventana de coda")
    ap.add_argument("--coda-step-s", type=float, default=300.0,
                    help="separacion entre ventanas de coda consecutivas")
    ap.add_argument("--noise-per-event", type=int, default=3,
                    help="cuantas ventanas de coda por evento-estacion. Es la palanca barata "
                         "para conseguir mas clutter sin picar el continuo completo.")
    args = ap.parse_args()

    import json
    with open(os.path.join(OUT, "fetch_report.json")) as fh:
        rep = json.load(fh)
    rotate = bool(rep.get("needs_rotation", False))
    # sta_coords se indexa solo por codigo de estacion: valido porque esta tarea trabaja con
    # UNA red a la vez (task1_fetch.py --net, por defecto CM). Si en el futuro se combinan
    # varias redes, esto colisionaria con nombres de estacion duplicados entre redes.
    sta_coords = {r["sta"]: (r["lat"], r["lon"]) for r in rep["channels"]}

    if args.start is not None:
        t0 = UTCDateTime(args.start)
        t1 = t0 + (args.days or 2.0) * 86400
    else:
        t0 = UTCDateTime(rep["start"])
        t1 = UTCDateTime(rep["end"])
    print(f"   ventana de analisis: {t0} .. {t1}  (de fetch_report.json salvo --start explicito)")

    inv = read_inventory(os.path.join(OUT, "inventory.xml"))
    print("1/4 cargando y filtrando formas de onda")
    st = load_streams(inv, args.band, rotate, t0=t0, t1=t1)
    print(f"   {len(st)} trazas, banda {args.band[0]}-{args.band[1]} Hz")

    # El catalogo se lee AQUI, antes del picker, porque en modo --around-events define que
    # tramos hay que picar. En la etapa 4 se reutiliza el mismo objeto.
    catf = os.path.join(OUT, "catalog.xml")
    if not os.path.exists(catf):
        sys.exit("falta catalog.xml; sin catalogo no hay etiquetado")
    cat = read_events(catf)
    print(f"   catalogo: {len(cat)} eventos")

    if args.around_events:
        vp_eff = args.vp if args.vp is not None else V_P_KMS
        st = slice_around_events(st, cat, sta_coords, vp_eff,
                                 pre_s=args.pre_s, post_s=args.post_s,
                                 coda_start_s=args.coda_start_s,
                                 coda_win_s=args.coda_win_s,
                                 coda_step_s=args.coda_step_s,
                                 noise_per_event=args.noise_per_event)
        if not len(st):
            sys.exit("el recorte dejo el stream vacio: revisa que el catalogo y la ventana de "
                     "forma de onda se solapen temporalmente.")
        print(f"   {len(st)} trazas tras el recorte")

    print("2/4 picker")
    picks = run_picker(st, args.p_thr, args.s_thr)
    print(f"   {len(picks)} picks (umbral bajo a proposito: queremos clutter)")
    if picks.empty:
        sys.exit("sin picks; baja los umbrales o revisa la banda")
    if picks["net"].nunique() > 1:
        # sta_coords se indexa solo por estacion (valido porque esta tarea trabaja con UNA red
        # a la vez, ver --net de task1_fetch.py); con mas de una red podria colisionar por
        # nombres de estacion duplicados entre redes.
        print(f"   AVISO: picks de mas de una red {sorted(picks['net'].unique())}; "
              f"sta_coords solo indexa por estacion y puede colisionar entre redes.")

    print("3/4 polarizacion")
    rec, az, inc, npts = [], [], [], []
    for _, r in picks.iterrows():
        out = polarization_at(st, r["sta"], r["t"], args.win)
        if out is None:
            rec.append(np.nan); az.append(np.nan); inc.append(np.nan); npts.append(0)
        else:
            rec.append(out[0]); az.append(out[1]); inc.append(out[2]); npts.append(out[3])
    picks["rectilinearity"] = rec
    picks["backazimuth_meas"] = az
    picks["incidence_meas"] = inc
    picks["n_samples"] = npts
    ok = picks["rectilinearity"].notna().sum()
    print(f"   {ok}/{len(picks)} picks con polarizacion calculable")

    print("4/4 asociacion contra el catalogo")
    picks = label_association(picks, cat, sta_coords, args.assoc_tol, vp_kms=args.vp)
    na = int((~picks["associated"]).sum())
    print(f"   asociados {int(picks['associated'].sum())} | NO asociados {na}")
    if na < 100:
        print("   AVISO: menos de 100 picks no asociados. La distribucion de coherencia sera")
        print("   ruidosa; amplia el periodo.")

    # desviacion angular plegada a [0,90]: el signo del autovector es arbitrario
    d = np.abs(picks["backazimuth_meas"] - (picks["ev_azimuth"] % 180.0))
    picks["angle_dev_deg"] = np.minimum(d, 180.0 - d)

    out = os.path.join(OUT, "picks_polarization.csv")
    picks.to_csv(out, index=False)
    print(f"\nescrito {out}")
    print("siguiente: python task1_decide.py")


if __name__ == "__main__":
    main()
