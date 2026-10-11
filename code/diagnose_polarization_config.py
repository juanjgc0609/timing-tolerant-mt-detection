"""
diagnose_polarization_config.py — por que falla el control de polarizacion y con que
configuracion del estimador (banda, ventana, alineacion) pasaria.

No re-corre PhaseNet: toma los picks de t1_data/picks_polarization.csv y recalcula la
polarizacion leyendo solo unos segundos de forma de onda alrededor de cada pick.

    python diagnose_polarization_config.py --stations CHI

Parte A (solo CSV)
  A1  control separado por fase (la S no se espera alineada con el back-azimut: SH es transversal)
  A2  control por residuo temporal del pick P (contaminacion de la clase "asociado")
  A3  control por incidencia (con incidencia casi vertical el azimut horizontal es ruido)
  A4  agrupamiento azimutal de los eventos y linea base por permutacion; con eventos agrupados
      en azimut la linea base no es 20/90 = 0.222

Parte B (formas de onda)
  Recalcula el control con picks P asociados para una grilla banda x ventana x alineacion con el
  mismo estimador de task1_analyze.py. La fila de referencia (banda 1-3, 0.6 s, centrada) debe
  reproducir el CSV; si no, la comparacion no vale. Elegir la mejor de muchas configuraciones sobre
  los mismos datos infla el resultado, asi que cada fila se evalua ademas en dos mitades disjuntas de
  eventos (ev_index par / impar). --grid propuesta prueba la banda primaria de la propuesta
  (0.1-2 Hz), que la Tarea 1 no cubrio.

Parte A tambien incluye
  A0  orientacion declarada en el inventario y epocas de canal
  A5  barrido de rotacion y de reflexion fijas (canal girado, invertido o intercambiado)
  A6  estabilidad temporal del control por trimestres
  A7  control por magnitud, profundidad y distancia epicentral (senal/ruido e incidencia)
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from obspy import read_inventory, UTCDateTime

from task1_analyze import polarization_at, WaveIndex, OUT, V_P_KMS

GRIDS = {
    # bandas de la Tarea 1 (1-10 Hz)
    "base": ([(1.0, 3.0), (2.0, 8.0), (1.0, 10.0)], [0.6, 1.0, 2.0, 3.0]),
    # banda primaria declarada en la Seccion 1.5 de la propuesta (0.1-2 Hz): ventanas mas largas
    # porque a 0.3 Hz un ciclo dura mas de 3 s
    "propuesta": ([(0.1, 2.0), (0.3, 1.0), (0.5, 2.0), (1.0, 2.0)], [1.0, 2.0, 4.0, 8.0]),
}
REF_CONFIG = ((1.0, 3.0), 0.6, "centered")      # la que produjo el CSV; se calcula siempre
ALIGNS = ["centered", "onset"]
CTRL_PASS = 0.60


def axial_R(deg):
    """Longitud resultante de datos axiales (periodo 180): 1 = todos iguales, 0 = uniformes."""
    a = np.radians(2.0 * np.asarray(deg, float))
    return float(np.abs(np.mean(np.exp(1j * a)))) if len(a) else np.nan


def axial_mean(deg):
    a = np.radians(2.0 * np.asarray(deg, float))
    return float((np.degrees(np.angle(np.mean(np.exp(1j * a)))) / 2.0) % 180.0)


def fold_dev(meas, ev_az):
    d = np.abs(np.asarray(meas, float) - (np.asarray(ev_az, float) % 180.0))
    return np.minimum(d, 180.0 - d)


def frac_coh(dev, coh_deg):
    dev = np.asarray(dev, float)
    return float((dev <= coh_deg).mean()) if dev.size else np.nan


def permutation_null(meas, ev_az, coh_deg, n=5000, seed=0):
    """Fraccion coherente si el azimut medido no tuviera relacion con el evento, conservando
    las distribuciones reales de ambos (se baraja el azimut del evento entre picks)."""
    rng = np.random.default_rng(seed)
    meas = np.asarray(meas, float)
    ev_az = np.asarray(ev_az, float)
    obs = frac_coh(fold_dev(meas, ev_az), coh_deg)
    null = np.array([frac_coh(fold_dev(meas, rng.permutation(ev_az)), coh_deg) for _ in range(n)])
    p = float((np.sum(null >= obs) + 1) / (n + 1))
    return obs, float(null.mean()), float(np.percentile(null, 95)), p


def line(label, sub, coh_deg, min_n=10):
    n = len(sub)
    if n < min_n:
        print(f"   {label:28s} n={n:4d}   (menos de {min_n}, no se reporta)")
        return
    f = frac_coh(sub["angle_dev_deg"], coh_deg)
    print(f"   {label:28s} n={n:4d}   fraccion coherente {f:.3f}   "
          f"desv. mediana {np.median(sub['angle_dev_deg']):5.1f}")


def metadata_check(sta, inv, rotate):
    """
    A0: lo que el inventario DECLARA sobre la orientacion de los canales, no lo que se mide.
    Un azimut/dip fuera de lo estandar aqui es un problema de METADATOS o instalacion (una
    causa corregible, distinta de un efecto de sitio real); varias epocas de un mismo canal que
    se solapan con el periodo analizado son la firma de un cambio de instrumento a mitad de
    periodo, que mezclaria dos orientaciones distintas bajo la misma etiqueta.
    """
    print("\nA0  metadatos de orientacion declarados en el inventario")
    print(f"   needs_rotation (fetch_report.json): {rotate}")
    if inv is None:
        print("   sin inventory.xml legible; no se puede verificar")
        return
    sel = inv.select(station=sta)
    chans = [(net.code, st.code, ch) for net in sel for st in net for ch in st]
    if not chans:
        print(f"   sin canales de {sta} en el inventario")
        return
    for net_code, st_code, ch in sorted(chans, key=lambda x: (x[2].code, str(x[2].start_date))):
        code = ch.code[-1]
        flag = ""
        if code == "Z" and abs(abs(ch.dip) - 90.0) > 5.0:
            flag = "  <-- dip no vertical para un canal Z"
        elif code == "N" and min(abs(ch.azimuth - 0.0), abs(ch.azimuth - 360.0)) > 5.0:
            flag = "  <-- azimut no 0 para un canal N (deberia estarlo si needs_rotation=False)"
        elif code == "E" and abs(ch.azimuth - 90.0) > 5.0:
            flag = "  <-- azimut no 90 para un canal E (deberia estarlo si needs_rotation=False)"
        print(f"   {ch.code:4s} loc={ch.location_code or '--':3s} azimut={ch.azimuth:6.1f} "
              f"dip={ch.dip:6.1f}  {str(ch.start_date)[:10]} - {str(ch.end_date)[:10] if ch.end_date else 'presente'}"
              f"{flag}")
    codes = {}
    for _, _, ch in chans:
        codes.setdefault(ch.code, []).append((ch.start_date, ch.end_date))
    for code, epochs in codes.items():
        if len(epochs) > 1:
            print(f"   AVISO: {code} tiene {len(epochs)} epocas en el inventario. Si dos se solapan "
                  f"con el periodo analizado, hay un cambio de instrumento a mitad de periodo que "
                  f"mezclaria dos orientaciones bajo la misma etiqueta de canal.")
    print("   si hay banderas arriba, el metadato mismo declara una orientacion no estandar: eso")
    print("   es un problema de metadatos/instalacion, no de fisica de propagacion en el sitio.")


def offset_reflection_sweep(cp, args, step=5):
    """
    A5: mismo barrido que diagnose_rus_offset.py, mas la familia de REFLEXIONES ademas de
    rotaciones puras. Una rotacion (proper, det=+1) desplaza el azimut medido pero conserva
    cuanto se agrupan los datos: meas -> (meas + off) % 180. Un canal N/E invertido, un canal con
    la polaridad cambiada, o E/N intercambiados son transformaciones IMPROPIAS (det=-1, reflexion)
    y no se pueden corregir con una rotacion: meas -> (off - meas) % 180. Si ninguna de las dos
    familias alcanza el umbral, el desajuste no es ni una rotacion fija ni una inversion fija de
    canal -- es dispersion real, no una transformacion lineal simple del sensor.
    """
    print("\nA5  barrido de rotacion Y reflexion fijas (miswiring corregible por una transformacion lineal)")
    if len(cp) < 10:
        print(f"   n={len(cp)} insuficiente para el barrido")
        return
    meas = cp["backazimuth_meas"].to_numpy(float)
    true = cp["ev_azimuth"].to_numpy(float) % 180.0
    best = {"rotacion": (0, -1.0), "reflexion": (0, -1.0)}
    for off in range(0, 180, step):
        for kind, cand in (("rotacion", (meas + off) % 180.0), ("reflexion", (off - meas) % 180.0)):
            d = np.abs(cand - true)
            d = np.minimum(d, 180.0 - d)
            coh = float((d <= args.coh_deg).mean())
            if coh > best[kind][1]:
                best[kind] = (off, coh)
    off0 = float((np.abs(meas - true) <= args.coh_deg).mean())
    print(f"   offset=0 (como esta ahora):      fraccion={off0:.3f}")
    for kind, (off, coh) in best.items():
        print(f"   mejor {kind:10s} (offset={off:3d}): fraccion={coh:.3f}")
    if off0 >= CTRL_PASS:
        print(f"   -> el offset actual (0) ya pasa {CTRL_PASS:.2f}: no hay error de orientacion "
              f"que corregir aqui (si el control global igual falla, la causa esta en otro lado).")
    else:
        kind = "rotacion" if best["rotacion"][1] >= best["reflexion"][1] else "reflexion"
        off, coh = best[kind]
        circ = min(off, 180 - off)
        if coh >= CTRL_PASS and circ > step:
            print(f"   -> {kind.upper()} FIJA DETECTADA (offset={off}): un desajuste de instalacion "
                  f"explicaria el problema. Corregible sin volver a instalar nada, ajustando la "
                  f"orientacion declarada en el inventario o en el post-proceso.")
        else:
            print("   -> ninguna transformacion lineal fija (rotacion o reflexion) alcanza "
                  f"{CTRL_PASS:.2f} con un offset real: no es un desajuste simple de orientacion.")


def temporal_stability(cp, args, n_bins=4):
    """
    A6: la fraccion coherente del control P por tramo cronologico. Un problema INSTRUMENTAL que
    aparece o desaparece durante el periodo (sensor reinstalado, cable reconectado) se ve como un
    salto entre tramos; un efecto de SITIO (scattering, incidencia) es estable en el tiempo.
    """
    print(f"\nA6  estabilidad temporal del control P ({n_bins} tramos cronologicos)")
    if len(cp) < 10 * n_bins:
        print(f"   n={len(cp)} insuficiente para {n_bins} tramos")
        return
    cp = cp.copy()
    cp["tt"] = cp["t"].apply(lambda x: float(pd.Timestamp(x).timestamp()) if not isinstance(x, (int, float)) else x)
    cp = cp.sort_values("tt")
    edges = np.linspace(0, len(cp), n_bins + 1).astype(int)
    fracs = []
    for k in range(n_bins):
        sub = cp.iloc[edges[k]:edges[k + 1]]
        f = frac_coh(sub["angle_dev_deg"], args.coh_deg)
        fracs.append(f)
        d0, d1 = str(sub["t"].iloc[0])[:10], str(sub["t"].iloc[-1])[:10]
        print(f"   tramo {k + 1}/{n_bins}  {d0} .. {d1}  n={len(sub):4d}  fraccion coherente {f:.3f}")
    valid = [f for f in fracs if np.isfinite(f)]
    if len(valid) >= 2 and (max(valid) - min(valid)) > 0.35:
        print("   AVISO: variacion grande entre tramos (>0.35). Podria ser un cambio a mitad de "
              "periodo (instrumental) en vez de un efecto de sitio estable; revisa A0 por si hay "
              "mas de una epoca de canal en ese rango de fechas.")
    else:
        print("   sin salto grande entre tramos: consistente con un efecto estable en el tiempo "
              "(sitio), no con un evento instrumental puntual.")


def part_a(df, is_assoc, args):
    rect_ok = df["rectilinearity"] >= args.rect_min
    near = df["ev_dist_km"] <= args.max_dist_km
    ctrl = df[is_assoc & rect_ok & df["angle_dev_deg"].notna()].copy()
    clut = df[~is_assoc & rect_ok & near & df["angle_dev_deg"].notna()].copy()

    print("\nA1  control por fase")
    line("todas las fases (como decide)", ctrl, args.coh_deg)
    for ph in ("P", "S"):
        line(f"fase {ph}", ctrl[ctrl["phase"] == ph], args.coh_deg)

    cp = ctrl[ctrl["phase"] == "P"].copy()
    print("\nA2  control P por residuo temporal |t_pick - t_P previsto|")
    if "best_dt_s" in cp.columns and len(cp):
        cp["resid"] = np.abs(cp["best_dt_s"] - cp["ev_dist_km"] / args.vp)
        for lo, hi in ((0, 1), (1, 3), (3, 1e9)):
            m = (cp["resid"] >= lo) & (cp["resid"] < hi)
            line(f"residuo {lo}-{'inf' if hi > 1e8 else hi} s", cp[m], args.coh_deg)
        print("   si solo los picks bien cronometrados son coherentes, la clase 'asociado' esta")
        print("   contaminada por picks falsos que caen dentro de la ventana P-S.")

    print("\nA3  control P por incidencia medida (grados desde la vertical)")
    for lo, hi in ((0, 20), (20, 40), (40, 91)):
        m = (cp["incidence_meas"] >= lo) & (cp["incidence_meas"] < hi)
        line(f"incidencia {lo}-{hi}", cp[m], args.coh_deg)

    if "prob" in cp.columns:
        print("\nA3b control P por probabilidad del picker")
        for lo, hi in ((0, 0.5), (0.5, 0.8), (0.8, 1.01)):
            m = (cp["prob"] >= lo) & (cp["prob"] < hi)
            line(f"prob {lo:.1f}-{min(hi, 1.0):.1f}", cp[m], args.coh_deg)

    print("\nA4  geometria y linea base por permutacion")
    print(f"   linea base ingenua (azimuts uniformes): {args.coh_deg / 90.0:.3f}")
    for name, sub in (("clutter", clut), ("control P", cp), ("control todas", ctrl)):
        if len(sub) < 10:
            continue
        R_ev = axial_R(sub["ev_azimuth"] % 180.0)
        R_me = axial_R(sub["backazimuth_meas"])
        obs, nm, n95, p = permutation_null(sub["backazimuth_meas"], sub["ev_azimuth"], args.coh_deg)
        print(f"   {name:14s} n={len(sub):4d} | R eventos {R_ev:.2f} (media {axial_mean(sub['ev_azimuth'] % 180):5.1f}) "
              f"| R medido {R_me:.2f} (media {axial_mean(sub['backazimuth_meas']):5.1f})")
        print(f"   {'':14s} observado {obs:.3f} | nulo por permutacion {nm:.3f} (p95 {n95:.3f}) | p={p:.4f}")
    print("   R cerca de 1 = azimuts concentrados. Si los eventos estan concentrados y el azimut")
    print("   medido tambien (por ejemplo ruido cultural de direccion fija), la coherencia del clutter")
    print("   puede superar 0.222 sin que el clutter 'apunte' a su evento: la comparacion correcta")
    print("   es contra el nulo por permutacion, no contra 0.222.")
    return cp


def event_covariates(cp, args):
    """
    A7: control P por magnitud, profundidad y distancia epicentral del evento, leidas de
    t1_data/catalog.xml por ev_index. Separa dos explicaciones que el resto no distingue: si solo
    fallan los eventos pequenos, la causa es la relacion senal/ruido y no el sitio; si fallan los
    profundos (incidencia casi vertical, por ejemplo el nido de Bucaramanga a ~150 km) la causa es
    la geometria del rayo. Un efecto de sitio falla parejo en todas las clases.
    """
    print("\nA7  control P por magnitud, profundidad y distancia epicentral del evento")
    catf = os.path.join(OUT, "catalog.xml")
    if not os.path.exists(catf) or cp.empty:
        print("   sin catalog.xml o sin control P; se omite")
        return
    from obspy import read_events
    cat = read_events(catf)
    mx = int(cp["ev_index"].max())
    if mx >= len(cat):
        print(f"   AVISO: ev_index llega a {mx} y el catalogo tiene {len(cat)} eventos: catalog.xml no")
        print("   es el que uso task1_analyze.py para este CSV (un fetch posterior lo sobrescribio).")
        return
    mag, dep = [], []
    for k in cp["ev_index"].astype(int):
        e = cat[k]
        o = e.preferred_origin() or e.origins[0]
        m = e.preferred_magnitude() or (e.magnitudes[0] if e.magnitudes else None)
        mag.append(float(m.mag) if m is not None else np.nan)
        dep.append(float(o.depth or 0.0) / 1000.0)
    cp = cp.copy()
    cp["mag"] = mag
    cp["depth_km"] = dep
    cp["epi_km"] = np.sqrt(np.clip(cp["ev_dist_km"] ** 2 - cp["depth_km"] ** 2, 0, None))
    for col, name, edges in (("mag", "magnitud", (-9, 2.5, 3.5, 99)),
                             ("depth_km", "profundidad km", (-1, 40, 100, 999)),
                             ("epi_km", "dist. epicentral km", (-1, 100, 150, 999))):
        for j, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            m = (cp[col] >= lo) & (cp[col] < hi)
            if j == 0:
                lab = f"< {hi:g}"
            elif j == len(edges) - 2:
                lab = f">= {lo:g}"
            else:
                lab = f"{lo:g}-{hi:g}"
            line(f"{name} {lab}", cp[m], args.coh_deg)
    print("   un sitio que degrada la polarizacion falla parejo en todas las clases; si la")
    print("   fraccion sube con la magnitud es senal/ruido, si baja con la profundidad es incidencia.")


def part_b(cp, args, sta, rotate, inv, bands, wins):
    if cp.empty:
        print("\nParte B: no hay picks P asociados")
        return None
    cp = cp.copy()
    cp["tt"] = cp["t"].apply(UTCDateTime)
    cp = cp.sort_values("tt")
    # margen de lectura: varios periodos de la frecuencia mas baja mas la ventana mas larga, para
    # que el transitorio del filtro no caiga sobre la ventana de covarianza
    PAD_S = max(30.0, 8.0 / min(b[0] for b in bands) + max(wins))
    idx = WaveIndex(cp["tt"].min() - 2 * PAD_S, cp["tt"].max() + 2 * PAD_S, stations={sta})
    print(f"\nParte B: {len(cp)} picks P asociados de {sta}, {idx.n_files(sta)} archivos mseed en el rango")
    if not idx.n_files(sta):
        print("   no hay mseed para esta estacion en el rango de los picks")
        return None

    configs = [(b, w, al) for b in bands for w in wins for al in ALIGNS]
    if REF_CONFIG not in configs:
        configs.append(REF_CONFIG)
    res = {c: [] for c in configs}
    bands_run = sorted({c[0] for c in configs})
    ref_cmp = []
    n_fail = 0
    for k, (_, r) in enumerate(cp.iterrows()):
        t = r["tt"]
        raw = idx.read_raw(sta, t - PAD_S, t + PAD_S)
        raw = raw.select(station=sta)
        if len(raw) < 3:
            n_fail += 1
            continue
        raw.merge(method=1)
        raw = raw.split()
        raw.detrend("demean").detrend("linear").taper(0.05)
        if rotate:
            try:
                raw.rotate("->ZNE", inventory=inv)
            except Exception:
                n_fail += 1
                continue
        for b in bands_run:
            st = raw.copy().filter("bandpass", freqmin=b[0], freqmax=b[1], corners=4, zerophase=True)
            for (bb, w, al) in configs:
                if bb != b:
                    continue
                out = polarization_at(st, sta, t, w, al)
                if out is None:
                    continue
                rect, az, inc, _ = out
                res[(b, w, al)].append((int(r["ev_index"]), rect, az, float(r["ev_azimuth"]), inc))
                if (b, w, al) == REF_CONFIG:
                    dd = abs(az - r["backazimuth_meas"])
                    ref_cmp.append(min(dd, 180.0 - dd))
        if (k + 1) % 25 == 0:
            print(f"   {k + 1}/{len(cp)}")
    if n_fail:
        print(f"   {n_fail} picks sin 3 componentes legibles alrededor del pick")

    print("\n   verificacion: fila de referencia (1-3 Hz, 0.6 s, centrada) contra el CSV")
    if ref_cmp:
        ref_cmp = np.asarray(ref_cmp)
        print(f"   diferencia de azimut con el CSV: mediana {np.median(ref_cmp):.1f} grados, "
              f"p90 {np.percentile(ref_cmp, 90):.1f} grados (n={len(ref_cmp)})")
        if np.median(ref_cmp) > 5.0:
            print("   AVISO: la referencia no reproduce el CSV (mediana > 5 grados). El filtrado por")
            print("   tramo corto difiere demasiado del continuo; no compares filas contra el CSV.")

    rows = []
    for (b, w, al), v in res.items():
        if not v:
            continue
        a = np.array(v, float)
        ok = a[:, 1] >= args.rect_min
        dev = fold_dev(a[ok, 2], a[ok, 3])
        even = (a[ok, 0] % 2) == 0
        rows.append(dict(band=f"{b[0]:g}-{b[1]:g}", win=w, align=al,
                         n_rect=int(ok.sum()), n_all=len(a),
                         frac=frac_coh(dev, args.coh_deg),
                         frac_even=frac_coh(dev[even], args.coh_deg),
                         frac_odd=frac_coh(dev[~even], args.coh_deg),
                         n_even=int(even.sum()), n_odd=int((~even).sum()),
                         med_dev=float(np.median(dev)) if dev.size else np.nan))
    tab = pd.DataFrame(rows).sort_values("frac", ascending=False)
    print(f"\n   control P (rect >= {args.rect_min:.2f}, coherente <= {args.coh_deg:.0f} grados), mejor primero")
    print(f"   {'banda':>7s} {'ventana':>7s} {'alineacion':>10s} {'n':>5s} {'fraccion':>8s} "
          f"{'mitad par':>9s} {'mitad impar':>11s} {'desv.med':>8s}")
    for _, r in tab.iterrows():
        mark = "  <- referencia" if (r["band"] == "1-3" and r["win"] == 0.6 and r["align"] == "centered") else ""
        print(f"   {r['band']:>7s} {r['win']:7.1f} {r['align']:>10s} {r['n_rect']:5d} {r['frac']:8.3f} "
              f"{r['frac_even']:9.3f} {r['frac_odd']:11.3f} {r['med_dev']:8.1f}{mark}")
    return tab


def verdict(tab, args):
    print("\n" + "=" * 72)
    if tab is None or tab.empty:
        print("sin resultado de la Parte B")
        return
    robust = tab[(tab["n_rect"] >= 30) & (tab["frac"] >= CTRL_PASS)
                 & (tab["frac_even"] >= CTRL_PASS) & (tab["frac_odd"] >= CTRL_PASS)]
    marginal = tab[(tab["n_rect"] >= 30) & (tab["frac"] >= CTRL_PASS)]
    if not robust.empty:
        r = robust.iloc[0]
        lo, hi = (float(x) for x in r["band"].split("-"))
        print(f"PASA en ambas mitades: banda {r['band']} Hz, ventana {r['win']} s, {r['align']} "
              f"(fraccion {r['frac']:.3f}, n={r['n_rect']})")
        print("El problema era la configuracion del estimador, no el sensor. Re-mide el clutter con:")
        print(f"   python task1_analyze.py --noise-per-event 10 --band {lo:g} {hi:g} "
              f"--win {r['win']:g} --win-align {r['align']}")
        print(f"   python task1_decide.py --stations {args.stations}")
    elif not marginal.empty:
        r = marginal.iloc[0]
        print(f"PASA en el total pero no en ambas mitades: banda {r['band']}, ventana {r['win']} s, "
              f"{r['align']} (total {r['frac']:.3f}, mitades {r['frac_even']:.3f}/{r['frac_odd']:.3f}).")
        print(f"No es robusto: puede ser efecto de elegir la mejor de {len(tab)} filas. No lo uses como")
        print("configuracion definitiva sin un periodo independiente que lo confirme.")
    else:
        best = tab.iloc[0]
        print(f"NINGUNA configuracion pasa {CTRL_PASS:.2f} (mejor: {best['frac']:.3f}, banda {best['band']}, "
              f"ventana {best['win']} s, {best['align']}).")
        print("El estimador de polarizacion de una estacion no recupera el back-azimut de eventos")
        print("conocidos en este sitio con esta precision. Mira A2-A4 para distinguir contaminacion")
        print("de la asociacion, incidencia casi vertical o un sitio con propagacion compleja.")
    print("=" * 72)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stations", default="CHI", help="una estacion (la Parte B se hace por estacion)")
    ap.add_argument("--coh-deg", type=float, default=20.0)
    ap.add_argument("--rect-min", type=float, default=0.60)
    ap.add_argument("--max-dist-km", type=float, default=200.0)
    ap.add_argument("--vp", type=float, default=V_P_KMS)
    ap.add_argument("--skip-waveforms", action="store_true", help="solo la Parte A")
    ap.add_argument("--grid", choices=["base", "propuesta", "todas"], default="base",
                    help="base: bandas 1-10 Hz de la Tarea 1; propuesta: banda primaria de la "
                         "propuesta, 0.1-2 Hz, con ventanas de 1 a 8 s; todas: la union")
    args = ap.parse_args()

    csv = os.path.join(OUT, "picks_polarization.csv")
    if not os.path.exists(csv):
        sys.exit("falta picks_polarization.csv; corre task1_analyze.py primero")
    df = pd.read_csv(csv)
    sta = args.stations.split(",")[0]
    df = df[df["sta"] == sta].copy()
    if df.empty:
        sys.exit(f"sin picks de {sta}")
    if df["associated"].dtype == object:
        is_assoc = df["associated"].astype(str).str.strip().str.upper().isin(["TRUE", "1"])
    else:
        is_assoc = df["associated"].astype(bool)
    print("=" * 72)
    print(f"DIAGNOSTICO DEL CONTROL DE POLARIZACION — {sta}")
    print("=" * 72)
    print(f"   picks {len(df)} | asociados {int(is_assoc.sum())} | no asociados {int((~is_assoc).sum())}")

    with open(os.path.join(OUT, "fetch_report.json")) as fh:
        rep = json.load(fh)
    rotate = bool(rep.get("needs_rotation", False))
    inv_path = os.path.join(OUT, "inventory.xml")
    inv = read_inventory(inv_path) if os.path.exists(inv_path) else None
    metadata_check(sta, inv, rotate)

    cp = part_a(df, is_assoc, args)
    offset_reflection_sweep(cp, args)
    temporal_stability(cp, args)
    event_covariates(cp, args)
    if args.skip_waveforms:
        return

    # Parte B usa TODOS los P asociados con polarizacion calculable, no solo rect >= rect_min del
    # CSV: la rectilinearidad cambia con la configuracion y se filtra por fila.
    cand = df[is_assoc & (df["phase"] == "P") & df["ev_azimuth"].notna()]
    if args.grid == "todas":
        bands = GRIDS["base"][0] + GRIDS["propuesta"][0]
        wins = sorted(set(GRIDS["base"][1]) | set(GRIDS["propuesta"][1]))
    else:
        bands, wins = GRIDS[args.grid]
    tab = part_b(cand, args, sta, rotate, inv, bands, wins)
    if tab is not None:
        out_csv = os.path.join(OUT, f"polarization_config_{sta}_{args.grid}.csv")
        tab.to_csv(out_csv, index=False)
        print(f"\n   escrito {out_csv}")
    verdict(tab, args)


if __name__ == "__main__":
    main()
