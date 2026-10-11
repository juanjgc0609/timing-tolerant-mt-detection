"""
froc_exact.py - curva FROC exacta para puntajes de grilla.

El catalogo que emite _catalog_impl a un umbral t es el prefijo, con puntaje >= t, de la lista voraz
completa (orden descendente con separacion minima y a lo sumo max_events): los eventos aceptados antes
de uno dado tienen mayor puntaje y no dependen de t. Por eso basta guardar, por ventana, esa lista con
su puntaje y el indice de la fuente verdadera que acierta (o -1 si es falso), y cada umbral posible es
uno de esos puntajes. El barrido por 80 cuantiles de sweep() muestreaba la curva con saltos de hasta
0.14 falsos por ventana en el extremo estricto (Tarea 1t), lo que hacia sensible la lectura a tope 0.1.
"""
import numpy as np

from task1o_consistent import TOL_KM


def greedy_record(grid, field, truth, min_sep_km=5.0, max_events=5):
    f = np.where(np.isfinite(field), field, -np.inf)
    order = np.argsort(-f)
    acc, rec = [], []
    for k in order:
        if not np.isfinite(f[k]):
            break
        p = np.asarray(grid[k], float)
        if all(np.linalg.norm(p[:2] - e[:2]) > min_sep_km for e in acc):
            acc.append(p)
            hit = -1
            if truth:
                dd = [float(np.linalg.norm(p[:2] - np.asarray(s)[:2])) for s in truth]
                j = int(np.argmin(dd))
                if dd[j] <= TOL_KM:
                    hit = j
            rec.append((float(f[k]), hit))
        if len(acc) >= max_events:
            break
    return rec


def exact_HF(records):
    """records: lista por ventana de [(puntaje, hit)]. Devuelve umbrales descendentes y matrices H, F
    (ventana x umbral) exactamente iguales a las de score_events sobre _catalog_impl a cada umbral."""
    thr = np.unique(np.concatenate([np.array([s for s, _h in r], float) for r in records if r]
                                   or [np.array([0.0])]))[::-1]
    H = np.zeros((len(records), len(thr)), int)
    F = np.zeros((len(records), len(thr)), int)
    for w, r in enumerate(records):
        hits, nf, pos = set(), 0, 0
        # r esta en orden descendente de puntaje
        for j, t in enumerate(thr):
            while pos < len(r) and r[pos][0] >= t:
                if r[pos][1] >= 0:
                    hits.add(r[pos][1])
                else:
                    nf += 1
                pos += 1
            H[w, j], F[w, j] = len(hits), nf
    return thr, H, F
