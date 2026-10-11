"""
task1c_network_survey.py — Tarea 1c: mapa de confiabilidad de la polarizacion de una estacion en
toda la red CM.

La Tarea 1 encontro que el control (la polarizacion P de una estacion apunta a eventos catalogados)
pasa en ZAR y falla en CHI y RUS, con OCA y SJC inconclusos por muestra. La pregunta siguiente es
de red: que fraccion de las estaciones CM tiene polarizacion utilizable, en que banda, y de que
depende. Para eso no hace falta el continuo: el control solo usa unos segundos alrededor de cada
arribo P catalogado. Este script descarga SOLO esas ventanas (unas 50 veces menos dato que el
continuo de la Tarea 1), pica P con PhaseNet dentro de cada una y mide la polarizacion en una
grilla de bandas y ventanas.

    python task1c_network_survey.py --start 2024-10-01 --days 90
    python task1c_network_survey.py --start 2024-10-01 --days 90 --stations ZAR,CHI,OCA --grid todas
    python task1c_network_survey.py --offline --grid propuesta        # re-analiza lo ya descargado

Salidas en t1c_survey/
    waveforms/<STA>.mseed      ventanas descargadas, reutilizadas en corridas siguientes
    inventory.xml, catalog.xml
    survey_picks.csv           un renglon por pick de control y configuracion
    survey_summary.csv         un renglon por estacion: estado del gate y diagnosticos

Estado por estacion (mismo criterio que diagnose_polarization_config.py)
    pasa        alguna configuracion llega a 0.60 con n >= 30, y tambien en las dos mitades
                disjuntas de eventos (par / impar)
    falla       n >= 30 en la configuracion de referencia y ninguna configuracion pasa
    inconcluso  menos de 30 picks de control
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from obspy import Stream, UTCDateTime, read, read_events, read_inventory
from obspy.geodetics import gps2dist_azimuth

from diagnose_polarization_config import (CTRL_PASS, GRIDS, REF_CONFIG, axial_R, fold_dev,
                                          frac_coh)
from task1_analyze import VP_OVER_VS, V_P_KMS, polarization_at, preprocess
from task1_fetch import SGC_URL, fetch_catalog, make_client

OUT = "t1c_survey"
WAVE = os.path.join(OUT, "waveforms")
MIN_N = 30
ASSOC_TOL_S = 5.0     # un pick P es de control si cae a <= 5 s de la P prevista
MAX_TRIES = 2         # intentos por ventana sin dato antes de darla por perdida


def grid_configs(name):
    if name == "todas":
        bands = GRIDS["base"][0] + GRIDS["propuesta"][0]
        wins = sorted(set(GRIDS["base"][1]) | set(GRIDS["propuesta"][1]))
    else:
        bands, wins = GRIDS[name]
    cfg = [(b, w, a) for b in bands for w in wins for a in ("centered", "onset")]
    if REF_CONFIG not in cfg:
        cfg.append(REF_CONFIG)
    return cfg


def cfg_key(c):
    (lo, hi), w, a = c
    return f"{lo:g}-{hi:g}|{w:g}|{a}"


# ------------------------------------------------------------------------------ inventario y catalogo
def get_inventory(client, t0, t1, stations):
    inv = client.get_stations(network="CM", station=stations or "*", channel="HH?",
                              starttime=t0, endtime=t1, level="channel")
    return inv


def station_table(inv):
    """Estaciones con las tres componentes HH, sus coordenadas y si hace falta rotar."""
    rows = {}
    for net in inv:
        for st in net:
            comps = {ch.code[-1]: ch for ch in st if ch.code.startswith("HH")}
            if not ({"Z", "N", "E"} <= set(comps) or {"Z", "1", "2"} <= set(comps)):
                continue
            rotate = "1" in comps or "2" in comps or any(
                ch.azimuth is not None and abs(((ch.azimuth - (0.0 if c == "N" else 90.0) + 180)
                                               % 360) - 180) > 1.0
                for c, ch in comps.items() if c in ("N", "E"))
            rows[st.code] = dict(sta=st.code, lat=float(st.latitude), lon=float(st.longitude),
                                 rotate=bool(rotate))
    return rows


def event_list(cat):
    ev = []
    for k, e in enumerate(cat):
        o = e.preferred_origin() or e.origins[0]
        m = e.preferred_magnitude() or (e.magnitudes[0] if e.magnitudes else None)
        ev.append(dict(k=k, t=o.time, lat=float(o.latitude), lon=float(o.longitude),
                       dep=float(o.depth or 0.0) / 1000.0, mag=float(m.mag) if m else np.nan))
    return ev


def station_events(ev, slat, slon, max_km, vp):
    vs = vp / VP_OVER_VS
    out = []
    for e in ev:
        d_m, _az, baz = gps2dist_azimuth(e["lat"], e["lon"], slat, slon)
        d = float(np.hypot(d_m / 1000.0, e["dep"]))
        if d > max_km:
            continue
        out.append(dict(e, dist=d, epi=d_m / 1000.0, baz=float(baz),
                        tP=e["t"] + d / vp, tS=e["t"] + d / vs))
    return sorted(out, key=lambda x: float(x["tP"]))


# ------------------------------------------------------------------------------ descarga
def download_windows(client, sta, sevs, pre, post, chunk=40, pause=0.5):
    """
    Descarga las ventanas [tP-pre, tP+post] que aun no estan en waveforms/<STA>.mseed. Pide por
    lotes con get_waveforms_bulk y, si el servidor no acepta POST, cae a una peticion por ventana.
    El SGC es un servidor publico modesto: lotes chicos y una pausa entre lotes.
    """
    path = os.path.join(WAVE, f"{sta}.mseed")
    donef = os.path.join(WAVE, f"{sta}.done.json")
    tries = {int(k): v for k, v in json.load(open(donef)).items()} if os.path.exists(donef) else {}
    # "hecho" = llego dato (-1) o ya se intento MAX_TRIES veces sin dato
    todo = [e for e in sevs if tries.get(e["k"], 0) != -1 and tries.get(e["k"], 0) < MAX_TRIES]
    if not todo:
        return 0
    got = 0
    for i in range(0, len(todo), chunk):
        part = todo[i:i + chunk]
        bulk = [("CM", sta, "*", "HH?", e["tP"] - pre, e["tP"] + post) for e in part]
        st = Stream()
        try:
            st = client.get_waveforms_bulk(bulk)
        except Exception:                                                   # noqa: BLE001
            for b in bulk:
                try:
                    st += client.get_waveforms(*b)
                except Exception:                                           # noqa: BLE001
                    pass
        st = st.split()          # un hueco dentro de la ventana llega enmascarado y no se puede escribir
        if len(st):
            with open(path, "ab") as fh:
                st.write(fh, format="MSEED")
            got += len(st)
        for e in part:
            covered = any(tr.stats.starttime <= e["tP"] <= tr.stats.endtime for tr in st)
            tries[e["k"]] = -1 if covered else tries.get(e["k"], 0) + 1
        with open(donef, "w") as fh:
            json.dump({str(k): v for k, v in tries.items()}, fh)
        time.sleep(pause)
    return got


# ------------------------------------------------------------------------------ medicion
def measure_station(sta, info, sevs, pre, post, configs, model, p_thr, inv):
    """Un renglon por (evento con pick P valido, configuracion)."""
    path = os.path.join(WAVE, f"{sta}.mseed")
    if not os.path.exists(path):
        return []
    # SIN merge global: unir ventanas separadas por dias estiraria arrays enmascarados sobre todo
    # el hueco entre eventos y agotaria la memoria (el mismo fallo de la Tarea 1 con el continuo).
    # Cada ventana se une y se procesa por separado en preprocess().
    allst = read(path)
    from task1_analyze import run_picker
    rows = []
    bands = sorted({c[0] for c in configs})
    for e in sevs:
        a, b = e["tP"] - ASSOC_TOL_S, e["tP"] + ASSOC_TOL_S
        raw = allst.slice(e["tP"] - pre, e["tP"] + post).copy()
        if len(raw) < 3:
            continue
        # el picker ve el tramo filtrado en la banda de referencia, como en task1_analyze.py
        try:
            ref = preprocess(raw.copy(), inv, REF_CONFIG[0], info["rotate"], taper_s=pre / 2.0)
        except Exception:                                                   # noqa: BLE001
            continue
        if model is None:
            t_pick, prob = e["tP"], np.nan
        else:
            pk = run_picker(model, ref, p_thr, p_thr)
            if pk.empty:
                continue
            pk = pk[(pk["phase"] == "P")]
            pk = pk[[(a <= t <= b) for t in pk["t"]]]
            if pk.empty:
                continue
            j = int(np.argmin([abs(float(t - e["tP"])) for t in pk["t"]]))
            t_pick, prob = pk["t"].iloc[j], float(pk["prob"].iloc[j])
        for band in bands:
            try:
                stb = preprocess(raw.copy(), inv, band, info["rotate"], taper_s=pre / 2.0)
            except Exception:                                               # noqa: BLE001
                continue
            for c in configs:
                if c[0] != band:
                    continue
                out = polarization_at(stb, sta, t_pick, c[1], c[2])
                if out is None:
                    continue
                rect, az, inc, _ = out
                rows.append(dict(sta=sta, ev_index=e["k"], t_pick=str(t_pick), prob=prob,
                                 resid_s=float(t_pick - e["tP"]), ev_baz=e["baz"],
                                 dist_km=e["dist"], epi_km=e["epi"], depth_km=e["dep"],
                                 mag=e["mag"], config=cfg_key(c), rect=rect, az=az, inc=inc))
    return rows


# ------------------------------------------------------------------------------ resumen
def chance_rate(meas, baz, coh_deg):
    """Valor esperado exacto de la fraccion coherente si el back-azimut se permutara entre picks."""
    d = np.abs(np.asarray(meas, float)[:, None] - (np.asarray(baz, float) % 180.0)[None, :])
    return float((np.minimum(d, 180.0 - d) <= coh_deg).mean())


def perm_p(meas, baz, coh_deg, n=4000, seed=0):
    rng = np.random.default_rng(seed)
    obs = frac_coh(fold_dev(meas, baz), coh_deg)
    null = np.array([frac_coh(fold_dev(meas, rng.permutation(baz)), coh_deg) for _ in range(n)])
    return float((np.sum(null >= obs) + 1) / (n + 1))


def best_fixed_transform(meas, baz, coh_deg, step=5):
    best = {"rotacion": (0, -1.0), "reflexion": (0, -1.0)}
    true = np.asarray(baz, float) % 180.0
    meas = np.asarray(meas, float)
    for off in range(0, 180, step):
        for kind, cand in (("rotacion", (meas + off) % 180.0), ("reflexion", (off - meas) % 180.0)):
            f = frac_coh(fold_dev(cand, true), coh_deg)
            if f > best[kind][1]:
                best[kind] = (off, f)
    return best


def summarize_station(df, sta, coh_deg, rect_min):
    s = df[(df["sta"] == sta) & (df["rect"] >= rect_min)]
    ref = s[s["config"] == cfg_key(REF_CONFIG)]
    row = dict(sta=sta, n_ref=len(ref))
    if len(ref):
        row["f_ref"] = frac_coh(fold_dev(ref["az"], ref["ev_baz"]), coh_deg)
        row["f0_ref"] = chance_rate(ref["az"], ref["ev_baz"], coh_deg)
        row["p_ref"] = perm_p(ref["az"].to_numpy(), ref["ev_baz"].to_numpy(), coh_deg)
        # exceso sobre el azar: 0 si la polarizacion no dice nada del evento, 1 si acierta siempre
        row["exceso_ref"] = (row["f_ref"] - row["f0_ref"]) / (1.0 - row["f0_ref"]) \
            if row["f0_ref"] < 1 else np.nan
        row["R_eventos"] = axial_R(ref["ev_baz"] % 180.0)
        row["R_medido"] = axial_R(ref["az"])
        bt = best_fixed_transform(ref["az"], ref["ev_baz"], coh_deg)
        row["mejor_rotacion"] = bt["rotacion"][1]
        row["offset_rotacion"] = bt["rotacion"][0]
        row["mejor_reflexion"] = bt["reflexion"][1]
        row["offset_reflexion"] = bt["reflexion"][0]
        if len(ref) >= 40:
            r2 = ref.sort_values("t_pick")
            edges = np.linspace(0, len(r2), 5).astype(int)
            q = [frac_coh(fold_dev(r2["az"].iloc[a:b], r2["ev_baz"].iloc[a:b]), coh_deg)
                 for a, b in zip(edges[:-1], edges[1:])]
            row["trimestres"] = "/".join(f"{x:.2f}" for x in q)
            row["salto_temporal"] = bool(max(q) - min(q) > 0.35)
        for col, lo, hi, name in (("mag", -9, 2.5, "f_mag_lt2.5"), ("mag", 3.5, 99, "f_mag_ge3.5"),
                                  ("depth_km", -1, 40, "f_prof_lt40"),
                                  ("depth_km", 100, 999, "f_prof_ge100")):
            m = (ref[col] >= lo) & (ref[col] < hi)
            row[name] = frac_coh(fold_dev(ref.loc[m, "az"], ref.loc[m, "ev_baz"]), coh_deg) \
                if m.sum() >= 10 else np.nan
            row["n_" + name[2:]] = int(m.sum())
    cfgrows = []
    for key, g in s.groupby("config"):
        dev = fold_dev(g["az"], g["ev_baz"])
        even = (g["ev_index"].to_numpy() % 2) == 0
        cfgrows.append(dict(config=key, n=len(g), f=frac_coh(dev, coh_deg),
                            f_even=frac_coh(dev[even], coh_deg),
                            f_odd=frac_coh(dev[~even], coh_deg)))
    ct = pd.DataFrame(cfgrows)
    if not ct.empty:
        ct = ct.sort_values("f", ascending=False)
        rob = ct[(ct["n"] >= MIN_N) & (ct["f"] >= CTRL_PASS) & (ct["f_even"] >= CTRL_PASS)
                 & (ct["f_odd"] >= CTRL_PASS)]
        top = rob.iloc[0] if not rob.empty else ct.iloc[0]
        row.update(mejor_config=top["config"], f_mejor=top["f"], n_mejor=int(top["n"]),
                   f_mejor_par=top["f_even"], f_mejor_impar=top["f_odd"])
        if not rob.empty:
            row["estado"] = "pasa"
        elif row["n_ref"] >= MIN_N:
            row["estado"] = "falla"
        else:
            row["estado"] = "inconcluso"
    else:
        row["estado"] = "sin datos"
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", default="2024-10-01")
    ap.add_argument("--days", type=float, default=90)
    ap.add_argument("--stations", default=None, help="lista separada por comas; por defecto toda la red CM con HH")
    ap.add_argument("--max-dist-km", type=float, default=200.0)
    ap.add_argument("--minmag", type=float, default=2.0)
    ap.add_argument("--max-events", type=int, default=400, help="tope de eventos por estacion")
    ap.add_argument("--pre", type=float, default=60.0,
                    help="segundos antes de la P prevista. 60 s dejan seis periodos a 0.1 Hz antes "
                         "del pick, para que el transitorio del filtro no llegue a la ventana")
    ap.add_argument("--post", type=float, default=60.0, help="segundos despues de la P prevista")
    ap.add_argument("--grid", choices=["base", "propuesta", "todas"], default="todas")
    ap.add_argument("--coh-deg", type=float, default=20.0)
    ap.add_argument("--rect-min", type=float, default=0.60)
    ap.add_argument("--p-thr", type=float, default=0.30)
    ap.add_argument("--vp", type=float, default=V_P_KMS)
    ap.add_argument("--no-picker", action="store_true",
                    help="medir en la P PREVISTA, sin PhaseNet (menos preciso; solo si no hay seisbench)")
    ap.add_argument("--offline", action="store_true",
                    help="no descargar nada: usar inventory.xml, catalog.xml y waveforms/ ya presentes")
    args = ap.parse_args()

    os.makedirs(WAVE, exist_ok=True)
    t0 = UTCDateTime(args.start)
    t1 = t0 + args.days * 86400
    invf, catf = os.path.join(OUT, "inventory.xml"), os.path.join(OUT, "catalog.xml")
    client = None if args.offline else make_client(SGC_URL, timeout=60)

    if args.offline:
        inv = read_inventory(invf)
    else:
        inv = get_inventory(client, t0, t1, args.stations)
        inv.write(invf, format="STATIONXML")
    stas = station_table(inv)
    if args.stations:
        stas = {k: v for k, v in stas.items() if k in set(args.stations.split(","))}
    print(f"   {len(stas)} estaciones CM con tres componentes HH: {', '.join(sorted(stas))}")
    if not stas:
        sys.exit("sin estaciones")

    if args.offline:
        cat = read_events(catf)
    else:
        lats = [v["lat"] for v in stas.values()]
        lons = [v["lon"] for v in stas.values()]
        m = args.max_dist_km / 111.0 + 0.2
        box = [min(lats) - m, max(lats) + m, min(lons) - m, max(lons) + m]
        cat, prov = fetch_catalog(t0, t1, box, args.minmag)
        if cat is None:
            sys.exit("ningun proveedor devolvio catalogo")
        cat.write(catf, format="QUAKEML")
        print(f"   catalogo {prov}: {len(cat)} eventos")
    ev = event_list(cat)

    configs = grid_configs(args.grid)
    model = None
    if not args.no_picker:
        from task1_analyze import load_picker
        model = load_picker()

    rows = []
    for sta in sorted(stas):
        info = stas[sta]
        sevs = station_events(ev, info["lat"], info["lon"], args.max_dist_km, args.vp)
        if len(sevs) > args.max_events:
            idx = np.linspace(0, len(sevs) - 1, args.max_events).astype(int)
            sevs = [sevs[i] for i in idx]
        if not args.offline:
            got = download_windows(client, sta, sevs, args.pre, args.post)
            print(f"   {sta}: {len(sevs)} eventos a <= {args.max_dist_km:.0f} km, {got} trazas nuevas")
        r = measure_station(sta, info, sevs, args.pre, args.post, configs, model, args.p_thr, inv)
        n_ev = len({x["ev_index"] for x in r})
        print(f"   {sta}: {n_ev} eventos con pick P valido")
        rows += r

    if not rows:
        sys.exit("ninguna medicion")
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "survey_picks.csv"), index=False)
    summ = pd.DataFrame([summarize_station(df, s, args.coh_deg, args.rect_min)
                         for s in sorted(df["sta"].unique())])
    summ.to_csv(os.path.join(OUT, "survey_summary.csv"), index=False)

    print("\n" + "=" * 100)
    print("MAPA DE CONFIABILIDAD DE LA POLARIZACION (control con picks P de eventos catalogados)")
    print("=" * 100)
    print(f"{'sta':<6}{'estado':<12}{'n ref':>6}{'f ref':>7}{'azar':>7}{'p':>8}"
          f"{'  mejor configuracion':<28}{'f':>6}{'par':>6}{'impar':>7}{'rot':>6}{'refl':>6}")
    for _, r in summ.iterrows():
        g = lambda k, fmt: (format(r[k], fmt) if k in r and pd.notna(r[k]) else "-")  # noqa: E731
        print(f"{r['sta']:<6}{r['estado']:<12}{int(r.get('n_ref', 0)):>6}{g('f_ref', '.3f'):>7}"
              f"{g('f0_ref', '.3f'):>7}{g('p_ref', '.4f'):>8}  {str(r.get('mejor_config', '-')):<26}"
              f"{g('f_mejor', '.3f'):>6}{g('f_mejor_par', '.2f'):>6}{g('f_mejor_impar', '.2f'):>7}"
              f"{g('mejor_rotacion', '.2f'):>6}{g('mejor_reflexion', '.2f'):>6}")
    counts = summ["estado"].value_counts().to_dict()
    print(f"\n   {counts}")
    print("   pasa = alguna configuracion >= 0.60 con n >= 30 y en ambas mitades de eventos;")
    print("   f ref / azar / p: configuracion de referencia (1-3 Hz, 0.6 s) contra el nulo por")
    print("   permutacion; rot / refl: mejor fraccion con una rotacion o reflexion fija (miswiring).")
    print(f"\n   escrito {OUT}/survey_summary.csv y {OUT}/survey_picks.csv")


if __name__ == "__main__":
    main()
