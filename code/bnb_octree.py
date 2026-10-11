"""
bnb_octree.py - busqueda por ramificacion y poda en un octree (x, y, profundidad, tiempo de origen)
con cotas CERTIFICADAS por la fisica, comparada con la grilla densa.

Cota. Para una celda con centro c y radio espacial h (distancia maxima del centro a cualquier punto
de la celda) y semiancho temporal dtau, el tiempo predicho de cualquier punto cambia a lo sumo
h / V_rama (Lipschitz del eikonal; en un medio heterogeneo, h / V_min) y el tiempo de origen a lo
sumo dtau. Como el nucleo de la verosimilitud cinematica decrece con el residuo, el puntaje de
cualquier punto de la celda esta acotado por el puntaje evaluado con residuos
max(0, r_i(c) - h/V - dtau). Se poda una celda solo si esa cota queda bajo el umbral de deteccion.

Consecuencia verificable: el conjunto de nodos de la grilla fina con puntaje sobre el umbral es
IDENTICO al de la grilla densa, y por tanto el catalogo emitido tambien: la poda no puede crear
falsos negativos (nada sobre el umbral se poda) ni falsos positivos (lo emitido se evalua exacto).
El script lo comprueba ventana por ventana y mide las evaluaciones ahorradas.
"""
import heapq
import time

import numpy as np

from task0_common import (DEFAULT_TAUS, DT, NT, VP, VS, _catalog_impl, _precompute, default_grid,
                          generate, make_network)
from task1o_consistent import EPS, degrade, draw_sources, score_events

SIG = 0.2
AX = np.arange(-21.0, 21.0 + 1e-9, 3.0) + 0.37
DEP = np.array([5.0, 9.0, 13.0])
TAUS = np.asarray(DEFAULT_TAUS)


def kin_terms(STA, NEAR_sorted, src, tau):
    d = np.linalg.norm(STA - src, axis=1)
    return np.stack([d / VP, d / VS]) + tau, d


def nearest_res(pk_sorted, tt):
    """residuo al pick mas cercano, estacion por estacion (picks ordenados)."""
    out = np.empty_like(tt)
    for i, p in enumerate(pk_sorted):
        for b in range(2):
            j = np.searchsorted(p, tt[b, i])
            c = [p[k] for k in (j - 1, j) if 0 <= k < len(p)]
            out[b, i] = min(abs(x - tt[b, i]) for x in c)
    return out


def score_from_res(res, u):
    nk = np.exp(-res ** 2 / (2 * SIG ** 2)) / (np.sqrt(2 * np.pi) * SIG) / u
    return float(np.log(EPS + (1 - EPS) * nk).sum())


def dense(STA, pk, u):
    grid = [np.array([x, y, z]) for x in AX for y in AX for z in DEP]
    f = np.full(len(grid), -np.inf)
    n = 0
    for gi, g in enumerate(grid):
        d = np.linalg.norm(STA - g, axis=1)
        for tau in TAUS:
            tt = np.stack([d / VP, d / VS]) + tau
            f[gi] = max(f[gi], score_from_res(nearest_res(pk, tt), u))
            n += 1
    return grid, f, n


def bnb(STA, pk, u, thr):
    """Devuelve el campo por nodo (exacto donde algun tau supera thr, -inf donde se podo)."""
    nx, ny, nz, nt = len(AX), len(AX), len(DEP), len(TAUS)
    field = np.full(nx * ny * nz, -np.inf)
    evals = 0
    stack = [((0, nx), (0, ny), (0, nz), (0, nt))]
    while stack:
        cx, cy, cz, ct = stack.pop()
        sizes = [cx[1] - cx[0], cy[1] - cy[0], cz[1] - cz[0], ct[1] - ct[0]]
        xs, ys, zs, ts = AX[cx[0]:cx[1]], AX[cy[0]:cy[1]], DEP[cz[0]:cz[1]], TAUS[ct[0]:ct[1]]
        c = np.array([(xs[0] + xs[-1]) / 2, (ys[0] + ys[-1]) / 2, (zs[0] + zs[-1]) / 2])
        h = float(np.sqrt(((xs[-1] - xs[0]) / 2) ** 2 + ((ys[-1] - ys[0]) / 2) ** 2 + ((zs[-1] - zs[0]) / 2) ** 2))
        tc, dt = (ts[0] + ts[-1]) / 2, (ts[-1] - ts[0]) / 2
        if max(sizes) == 1:
            ix = ((cx[0] * ny) + cy[0]) * nz + cz[0]
            d = np.linalg.norm(STA - c, axis=1)
            s = score_from_res(nearest_res(pk, np.stack([d / VP, d / VS]) + tc), u)
            evals += 1
            if s >= thr:
                field[ix] = max(field[ix], s)
            continue
        d = np.linalg.norm(STA - c, axis=1)
        tt = np.stack([d / VP, d / VS]) + tc
        res = nearest_res(pk, tt)
        slack = np.array([[h / VP + dt], [h / VS + dt]])
        ub = score_from_res(np.maximum(0.0, res - slack), u)
        evals += 1
        if ub < thr:
            continue
        # dividir la dimension mas larga en unidades de indice
        k = int(np.argmax(sizes))
        rng_ = [cx, cy, cz, ct]
        a, b = rng_[k]
        m = (a + b) // 2
        for part in ((a, m), (m, b)):
            r2 = list(rng_)
            r2[k] = part
            stack.append(tuple(r2))
    return field, evals


def main():
    rng = np.random.default_rng(5)
    rng_d = np.random.default_rng(6)
    for kind, nst, ncl in (("random", 4, 35), ("onesided", 8, 40)):
        data = []
        for k in range(24):
            S = make_network(kind, nst, seed=900 + k)
            srcs, taus = draw_sources(rng)
            D, p, tr = generate(S, srcs, taus, n_clutter=ncl, rng=rng)
            p2, lab = degrade(S, srcs, p, "consistent", rng_d)
            u = max(1e-3, sum(1 for li in lab for l in li if l == 0) / (nst * 40.0))
            pk = [np.sort([t for t, _b in pl]) for pl in p2]
            data.append((S, pk, tr, u))
        # umbral: calibrado a recall 0.9 en las 12 primeras ventanas con la grilla densa
        fields = []
        t0 = time.time()
        for S, pk, tr, u in data:
            grid, f, n_dense = dense(S, pk, u)
            fields.append(f)
        t_dense = (time.time() - t0) / len(data)
        cal = list(zip(fields[:12], data[:12]))
        ntr = sum(len(d[2]) for _f, d in cal)
        thr = None
        for q in np.unique(np.quantile(np.concatenate([f for f, _ in cal]), np.linspace(0.9, 1, 200)))[::-1]:
            h = sum(score_events(_catalog_impl(grid, f, float(q)), d[2])[0] for f, d in cal)
            if h / max(1, ntr) >= 0.9:
                thr = float(q)
                break
        if thr is None:
            thr = float(np.quantile(np.concatenate([f for f, _ in cal]), 0.9))
        same, ev_tot, t_b = 0, 0, 0.0
        for f, (S, pk, tr, u) in zip(fields[12:], data[12:]):
            t1 = time.time()
            fb, ev = bnb(S, pk, u, thr)
            t_b += time.time() - t1
            ev_tot += ev
            cat_d = _catalog_impl(grid, f, thr)
            cat_b = _catalog_impl(grid, fb, thr)
            eq = (len(cat_d) == len(cat_b) and all(np.allclose(a, b) for a, b in zip(cat_d, cat_b))
                  and np.array_equal(np.where(f >= thr)[0], np.where(fb >= thr)[0]))
            same += int(eq)
        print(f"{kind} {nst} est.: catalogos identicos {same}/12 | evaluaciones por ventana densa "
              f"{n_dense}, ramificacion y poda {ev_tot / 12:.0f} ({n_dense / (ev_tot / 12):.1f}x menos) | "
              f"tiempo por ventana densa {t_dense:.2f} s, poda {t_b / 12:.2f} s")


if __name__ == "__main__":
    main()
