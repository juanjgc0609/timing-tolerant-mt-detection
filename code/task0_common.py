"""
task0_common.py — generador sintético y métricas para la Tarea 0. VERSIÓN CORREGIDA.

Correcciones respecto de la versión anterior:
  C1. `calibrate_threshold` ya no devuelve el umbral límite (el primero que roza el objetivo
      recorriendo de mayor a menor), que es el punto más frágil y sobreajusta el conjunto de
      calibración. Ahora elige un umbral con holgura y devuelve si el objetivo es alcanzable.
  C2. Se añade `operating_curve`: barre el umbral y devuelve la curva (recall, tasa de falsos)
      completa. Comparar curvas en vez de puntos elimina la discusión sobre igualación de
      puntos de operación.
  C3. `max_events` deja de estar fijo: forma parte del punto de operación y se expone.
  C4. `picks_to_dataframe(fmt='pyocto')` emite `time` NUMÉRICO POSIX, no datetime.
"""
import numpy as np
import pandas as pd
from datetime import datetime, timedelta, timezone

FC = 2.0
VS = 3.0
VP = 1.73 * VS
DT = 1 / (8 * FC)
T = np.arange(0, 40, DT)
NT = len(T)
SNR = 7.0
TOL = 0.6
AP, AS = 1.0, 1.5
VMIN = 2.0
EPOCH = datetime(2024, 1, 1, 0, 0, 0, tzinfo=timezone.utc)


def make_network(kind="spread", n=8, seed=None):
    rng = np.random.default_rng(seed)
    if kind == "spread":
        a = np.linspace(0, 2 * np.pi, n, endpoint=False) + 0.3
        r = np.linspace(9, 58, n)
    elif kind == "ring":
        a = np.linspace(0, 2 * np.pi, n, endpoint=False) + rng.uniform(0, 1)
        r = np.full(n, 28.0)
    elif kind == "onesided":
        a = rng.uniform(0, np.pi * 0.7, n)
        r = 45 * np.sqrt(rng.uniform(0.05, 1, n))
    elif kind == "random":
        a = rng.uniform(0, 2 * np.pi, n)
        r = 45 * np.sqrt(rng.uniform(0.05, 1, n))
    else:
        raise ValueError(kind)
    return np.stack([r * np.cos(a), r * np.sin(a), np.zeros(n)], axis=1)


def ricker(t):
    x = (np.pi * FC * t) ** 2
    return (1 - 2 * x) * np.exp(-x)


def polarization(src, sta, branch):
    v = np.asarray(sta, float) - np.asarray(src, float)
    n = v / np.linalg.norm(v)
    if branch == "P":
        return n
    ez = np.array([0.0, 0.0, 1.0])
    sv = ez - np.dot(ez, n) * n
    nn = np.linalg.norm(sv)
    return sv / nn if nn > 1e-9 else np.array([1.0, 0.0, 0.0])


def generate(STA, sources, taus, n_clutter=30, rng=None,
             coherent_decoy=None, decoy_support=0, decoy_tau=None, noise_scale=0.02,
             statics=None, dir_noise=None):
    """statics (ns, 2): retardo por estacion y fase que el modelo homogeneo no predice (medio
    heterogeneo). dir_noise (ns, 2, 3): perturbacion de la direccion de polarizacion por estacion y
    fase. Ambos None reproducen el generador original."""
    rng = rng or np.random.default_rng()
    ns = len(STA)
    D = np.zeros((ns, 3, NT))
    picks = [[] for _ in range(ns)]

    for src, tau in zip(sources, taus):
        d = np.linalg.norm(STA - src, axis=1)
        for i in range(ns):
            for bi, (b, A, V) in enumerate((("P", AP, VP), ("S", AS, VS))):
                ta = d[i] / V + tau + (0.0 if statics is None else float(statics[i, bi]))
                picks[i].append((ta, b))
                pv = polarization(src, STA[i], b)
                if dir_noise is not None:
                    pv = pv + dir_noise[i, bi]
                    pv = pv / np.linalg.norm(pv)
                D[i] += SNR * A * np.outer(pv, ricker(T - ta)) / d[i]

    if coherent_decoy is not None and decoy_support > 0:
        dd = np.linalg.norm(STA - coherent_decoy, axis=1)
        td = decoy_tau if decoy_tau is not None else taus[0] + 1.5
        slots = [(i, b) for i in range(ns) for b in ("P", "S")]
        rng.shuffle(slots)
        for i, b in slots[:decoy_support]:
            A, V = (AP, VP) if b == "P" else (AS, VS)
            ta = dd[i] / V + td
            picks[i].append((ta, b))
            D[i] += SNR * A * np.outer(polarization(coherent_decoy, STA[i], b),
                                       ricker(T - ta)) / dd[i]

    # La amplitud del clutter se normaliza por la distancia mediana ESTACION-FUENTE, no por la
    # distancia al origen de coordenadas: esta ultima no es invariante por traslacion y hacia
    # que el generador dependiera del sistema de referencia.
    ref = np.asarray(sources[0], float) if len(sources) else np.zeros(3)
    scale = float(np.median(np.linalg.norm(STA - ref, axis=1)))
    for _ in range(n_clutter):
        i = int(rng.integers(ns))
        tc = float(rng.uniform(2, T[-1] - 4))
        q = rng.normal(size=3)
        q /= np.linalg.norm(q)
        D[i] += SNR * 0.35 * np.outer(q, ricker(T - tc)) / scale
        picks[i].append((tc, "P" if rng.random() < 0.5 else "S"))

    if noise_scale > 0:
        D += rng.normal(scale=SNR * noise_scale, size=D.shape)
    return D, picks, [np.asarray(s, float) for s in sources]


def picks_to_dataframe(picks, fmt="pyocto"):
    rows = []
    for i, plist in enumerate(picks):
        for (t, ph) in plist:
            rows.append((f"S{i:02d}", ph, float(t)))
    df = pd.DataFrame(rows, columns=["station", "phase", "t_sec"])
    if fmt == "pyocto":
        out = pd.DataFrame({
            "station": df["station"],
            "phase": df["phase"],
            "time": EPOCH.timestamp() + df["t_sec"].to_numpy(float),
        })
        return out.sort_values("time").reset_index(drop=True)
    if fmt == "gamma":
        ts = [EPOCH + timedelta(seconds=float(s)) for s in df["t_sec"]]
        out = pd.DataFrame({
            "id": df["station"],
            "timestamp": pd.to_datetime(ts),
            "type": df["phase"].str.lower(),
            "prob": 1.0,
            "amp": 1.0,
        })
        return out.sort_values("timestamp").reset_index(drop=True)
    raise ValueError(fmt)


def stations_to_dataframe(STA, fmt="pyocto"):
    ids = [f"S{i:02d}" for i in range(len(STA))]
    if fmt == "pyocto":
        return pd.DataFrame({"id": ids, "x": STA[:, 0].astype(float),
                             "y": STA[:, 1].astype(float), "z": STA[:, 2].astype(float)})
    if fmt == "gamma":
        return pd.DataFrame({"id": ids, "x(km)": STA[:, 0].astype(float),
                             "y(km)": STA[:, 1].astype(float),
                             "z(km)": STA[:, 2].astype(float)})
    raise ValueError(fmt)


def evaluate_catalog(events_xyz, truth, tol_km=5.0):
    if len(truth) == 0:
        return False, len(events_xyz)
    rec, nf = False, 0
    for e in events_xyz:
        e = np.asarray(e, float)
        d = min(float(np.linalg.norm(e[:2] - np.asarray(s, float)[:2])) for s in truth)
        if d <= tol_km:
            rec = True
        else:
            nf += 1
    return rec, nf


def validate_positive_control(STA, picks, truth):
    if not truth:
        raise AssertionError("positive-control sin fuentes verdaderas")
    n_true = 2 * len(STA) * len(truth)
    n_tot = sum(len(p) for p in picks)
    if n_tot < n_true:
        raise AssertionError(f"faltan picks verdaderos: {n_tot} < {n_true}")
    return True


def _precompute(D, picks):
    ns = len(D)
    W = int(0.6 / DT)
    PD = np.zeros((ns, NT, 3))
    HIT = np.zeros((ns, NT), bool)
    NEAR = np.zeros((ns, NT))
    for i in range(ns):
        x = D[i]
        prod = np.einsum("at,bt->abt", x, x)
        cs = np.concatenate([np.zeros((3, 3, 1)), np.cumsum(prod, axis=2)], axis=2)
        lo = np.clip(np.arange(NT) - W // 2, 0, NT)
        hi = np.clip(np.arange(NT) + W // 2, 0, NT)
        C = (cs[:, :, hi] - cs[:, :, lo]).transpose(2, 0, 1)
        PD[i] = np.linalg.eigh(C)[1][:, :, -1]
        p = np.sort(np.array([t for (t, _b) in picks[i]], float))
        idx = np.clip(np.searchsorted(p, T), 0, len(p) - 1)
        cand = np.stack([p[np.clip(idx - 1, 0, len(p) - 1)], p[idx]])
        best = cand[np.argmin(np.abs(cand - T), axis=0), np.arange(NT)]
        NEAR[i] = best
        HIT[i] = np.abs(best - T) < TOL
    return PD, HIT, NEAR


def score_grid(STA, D, picks, use, grid, taus):
    ns = len(STA)
    PD, HIT, NEAR = _precompute(D, picks)
    DIJ = np.linalg.norm(STA[:, None, :] - STA[None, :, :], axis=2)
    st = np.arange(ns)
    out = np.zeros(len(grid))
    for gi, src in enumerate(grid):
        d = np.linalg.norm(STA - src, axis=1)
        best = 0.0
        for tau in taus:
            tt = np.stack([d / VP, d / VS]) + tau
            idx = np.clip((tt / DT).astype(int), 0, NT - 1)
            hit = HIT[st[None, :], idx]
            if hit.sum() == 0:
                continue
            if "CIII" in use:
                cl = NEAR[st[None, :], idx]
                bad = False
                for b in range(2):
                    ii = np.where(hit[b])[0]
                    if len(ii) > 1 and np.any(
                        np.abs(cl[b][ii][:, None] - cl[b][ii][None, :])
                        > DIJ[np.ix_(ii, ii)] / VMIN + 1e-9):
                        bad = True
                        break
                if bad:
                    continue
            s = hit.sum() / (2 * ns)
            if "CIV" in use:
                cs = [abs(float(np.dot(polarization(src, STA[i], bb), PD[i, idx[b, i]])))
                      for b, bb in enumerate(("P", "S")) for i in range(ns) if hit[b, i]]
                s *= float(np.mean(cs)) ** 2 if cs else 0.0
            best = max(best, s)
        out[gi] = best
    return out


def field_own(STA, D, picks, use, grid, taus):
    return score_grid(STA, D, picks, use, grid, taus)


def catalog_from_score_field(grid, field, thresh, min_sep_km=5.0, max_events=5,
                             relative=True):
    """
    C5. `thresh` se interpreta por defecto como FRACCION del maximo de ESTA realizacion, no como
    valor absoluto. Un umbral absoluto solo es valido si la escala del puntaje es la misma en
    todas las realizaciones, lo cual es cierto para la ablacion (su maximo es 1.000 por
    construccion) y FALSO para el modelo restringido (maximo 0.51-0.68 segun la realizacion).
    Umbralizar en absoluto favorecia estructuralmente a la ablacion.
    """
    field = np.asarray(field, float)
    thresh_abs = float(thresh) * float(field.max()) if relative else float(thresh)
    return _catalog_impl(grid, field, thresh_abs, min_sep_km, max_events)


def _catalog_impl(grid, field, thresh_abs, min_sep_km=5.0, max_events=5):
    grid = [np.asarray(g, float) for g in grid]
    field = np.asarray(field, float)
    if len(grid) != len(field):
        raise ValueError("grid y field deben tener la misma longitud.")
    if not np.isfinite(thresh_abs):
        raise ValueError("thresh_abs debe ser finito.")
    order = np.argsort(-field)
    events = []
    for k in order:
        if field[k] < thresh_abs:
            break
        p = grid[k]
        if all(np.linalg.norm(p[:2] - e[:2]) > min_sep_km for e in events):
            events.append(p.copy())
        if len(events) >= max_events:
            break
    return events


def _recall_false_at(fields, truths, grid, thr, min_sep_km, max_events, tol_km):
    rec = fal = 0
    for f, tr in zip(fields, truths):
        ev = catalog_from_score_field(grid, f, thr, min_sep_km, max_events)
        r, nf = evaluate_catalog(ev, tr, tol_km)
        rec += int(r)
        fal += int(nf > 0)
    n = len(fields)
    return rec / n, fal / n


def calibrate_threshold(fields, truths, grid, target_recall=0.99, min_sep_km=5.0,
                        max_events=5, tol_km=5.0):
    """
    C1. Los umbrales validos son todos los que estan POR DEBAJO del umbral limite, porque bajar
    el umbral solo puede aumentar el recall. Queremos el mas ALTO posible (menos eventos, menos
    falsos) pero con holgura respecto del limite, que es el punto fragil. De ahi el percentil
    ALTO (80 por defecto) del conjunto valido, no uno bajo.
    Devuelve (umbral, recall_calibracion, alcanzable).
    `alcanzable=False` no es un fallo del protocolo: es información sobre el método.
    """
    if not fields:
        raise ValueError("No hay campos para calibrar.")
    if len(fields) != len(truths):
        raise ValueError("fields y truths deben tener igual longitud.")
    cand = np.linspace(1.0, 0.40, 61)          # fracciones del maximo, de estricto a laxo

    valid, best = [], (float(cand[-1]), 0.0)
    for thr in cand:
        r, _f = _recall_false_at(fields, truths, grid, float(thr), min_sep_km,
                                 max_events, tol_km)
        if r > best[1]:
            best = (float(thr), float(r))
        if r >= target_recall:
            valid.append(float(thr))
    if not valid:
        return best[0], best[1], False
    # C1b. Entre los umbrales que alcanzan el recall objetivo, elegir el que MINIMIZA la tasa
    # de falsos. Elegir un percentil del conjunto valido optimiza una sola columna e ignora la
    # otra, que fue el defecto de la version anterior: PHYS calibraba a umbrales bajos con
    # recall 1.00 y tasa de falsos 1.00. Desempate por umbral mas alto, que es el mas estable.
    scored = []
    for thr in valid:
        r, f = _recall_false_at(fields, truths, grid, thr, min_sep_km, max_events, tol_km)
        scored.append((float(thr), float(r), float(f)))
    scored.sort(key=lambda s: (s[2], -s[0]))
    chosen, r, _f = scored[0]
    return float(chosen), float(r), True


def operating_curve(fields, truths, grid, n_points=25, min_sep_km=5.0,
                    max_events=5, tol_km=5.0):
    """C2. Curva característica: [(umbral, recall, tasa_falsos)], ordenada por recall."""
    obs = [np.asarray(f)[np.isfinite(f)] for f in fields]
    obs = [o for o in obs if o.size]
    if not obs:
        return []
    return _sweep_fractions(fields, truths, grid, n_points, min_sep_km, max_events, tol_km)


def _sweep_fractions(fields, truths, grid, n_points, min_sep_km, max_events, tol_km):
    """Barrido en FRACCION del maximo de cada realizacion (ver C5)."""
    fr = np.linspace(0.20, 1.02, n_points)
    out = []
    for f_ in fr[::-1]:
        r, fa = _recall_false_at(fields, truths, grid, float(f_), min_sep_km,
                                 max_events, tol_km)
        out.append((float(f_), float(r), float(fa)))
    out.sort(key=lambda t: t[1])
    return out


def _unused_abs_grid(obs, n_points):
    allv = np.concatenate(obs)
    nz = allv[allv > 1e-12]            # la grilla esta dominada por ceros; usarlos falsea
    base = nz if nz.size >= 8 else allv
    # cuantiles densos en la parte ALTA, que es donde vive el punto de operacion util
    qs = np.concatenate([np.linspace(0.0, 80.0, max(4, n_points // 3)),
                         np.linspace(82.0, 100.0, n_points - max(4, n_points // 3))])
    cand = np.unique(np.percentile(base, qs))
    out = []
    for thr in np.sort(cand)[::-1]:
        r, f = _recall_false_at(fields, truths, grid, float(thr), min_sep_km,
                                max_events, tol_km)
        out.append((float(thr), float(r), float(f)))
    out.sort(key=lambda t: t[1])
    return out


def false_rate_at_recall(curve, target_recall):
    """Tasa de falsos mínima entre los puntos que alcanzan el recall. None si no lo alcanza."""
    ok = [(r, f) for (_t, r, f) in curve if r >= target_recall]
    return None if not ok else min(f for _r, f in ok)


def default_grid(half=21.0, step=3.0, depths=(5.0, 9.0, 13.0)):
    ax = np.arange(-half, half + 1e-9, step) + 0.37   # verdad nunca sobre un nodo
    return [np.array([x, y, z], float) for x in ax for y in ax for z in depths]


DEFAULT_TAUS = np.arange(5.0, 8.01, 0.1)
