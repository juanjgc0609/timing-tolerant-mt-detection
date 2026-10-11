"""
task1x_mtglrt.py - Tarea 1x: prueba piloto del metodo reformulado.

Idea. La forma de onda de tres componentes que una fuente puntual produce en todas las estaciones es
LINEAL en su tensor de momento m (6 componentes), para una hipotesis (x, t) dada: v = G(x) m + ruido.
El mecanismo desconocido deja de ser un estorbo: se maximiza analiticamente y la prueba de razon de
verosimilitud generalizada es la energia de los datos proyectada sobre el subespacio de dimension 6
que la fisica permite, conjuntamente en todas las estaciones y en P y S (detector de subespacio
adaptado, Scharf y Friedlander 1994). A diferencia de WAVE-PT, que suma energias por estacion, el
estadistico exige coherencia de signos y amplitudes relativas entre estaciones (patron de radiacion).
La misma cantidad es el certificado dual del problema convexo sobre medidas con amplitudes
tensoriales (BLASSO de grupo): su maximo global es el oraculo lineal de Frank-Wolfe, y restar la
contribucion ajustada y repetir es Frank-Wolfe con deflacion, que explica eventos superpuestos.

Metodos (todos sin picks, salvo GaMMA):
    WAVE-E        energia de tres componentes (ablacion sin direccion)
    WAVE-PT       P sobre el rayo, S en el plano transversal (Y.10)
    WAVE-MT       GLRT con el tensor de momento como estorbo lineal: ||Q(x)^T v||^2
    WAVE-MT-NODE  control: el subespacio de otro nodo a mas de 9 km (misma dimension, sin la geometria)
    MT-FW         Frank-Wolfe con deflacion sobre WAVE-MT (hasta 5 eventos por ventana)
    GAMMA         con semilla, rejilla de 1r, oraculo
Condiciones: "ideal" (fuentes de doble par, medio homogeneo conocido) y "pert" (ademas retardos por
estacion y fase N(0, 0.15 s), desviacion de polarizacion de ~10 grados y factor de sitio lognormal de
amplitud con sigma 0.5, desconocidos para todos los metodos: la debilidad conocida de los metodos
de tensor de momento a alta frecuencia es precisamente la amplitud).

    python task1x_mtglrt.py --seed 2024 --reps 100 --cond ideal
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import argparse
import json
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.stats import norm

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "gamma_src"))
warnings.filterwarnings("ignore")

from task0_common import AP, AS, DEFAULT_TAUS, DT, NT, SNR, T, VP, VS, default_grid, make_network, polarization, ricker
from task1o_consistent import _GD2, _KER, NOISE_SD, TOL_KM, draw_sources, matched, score_events
from task1r_baselines import GAMMA_GRID, run_gamma
from task1r_fair import CAPS, SCEN, recall_at
from task1w_mech import CP, CS, rand_dc
from froc_exact import exact_HF, greedy_record

STATIC_SD, DIR_SD, SITE_SD = 0.15, float(np.tan(np.radians(10.0))), 0.5
METH = ("WAVE-E", "WAVE-PT", "WAVE-MT", "WAVE-MT-NODE")
CAPV = 6.0                      # tope por estacion y fase de la norma del vector (energia 36), como en 1t
SH = int(round(0.25 / DT))
OFFS = np.arange(-SH, SH + 1)
# forma de la salida del filtro adaptado para un pulso: autocorrelacion normalizada a 1 en el pico
_ac = np.convolve(ricker(np.arange(-3.0, 3.0 + 1e-9, DT)), _KER[::-1], mode="same")
AC = _ac / _ac.max()
AC_HALF = len(AC) // 2


def generate_x(STA, sources, taus, n_clutter, rng, rng_m, pert, rng_h):
    """generate_dc de 1w, con perturbaciones opcionales; consume rng en el mismo orden que generate()."""
    ns = len(STA)
    D = np.zeros((ns, 3, NT))
    picks = [[] for _ in range(ns)]
    rad = [[] for _ in range(ns)]
    if pert:
        st = rng_h.normal(0.0, STATIC_SD, size=(ns, 2))
        dn = rng_h.normal(0.0, DIR_SD, size=(ns, 2, 3))
        site = np.exp(rng_h.normal(0.0, SITE_SD, size=(ns, 2)))
    for src, tau in zip(sources, taus):
        d = np.linalg.norm(STA - src, axis=1)
        M = rand_dc(rng_m)
        for i in range(ns):
            n = (STA[i] - src) / d[i]
            for bi, (b, A, V) in enumerate((("P", AP, VP), ("S", AS, VS))):
                ta = d[i] / V + tau + (float(st[i, bi]) if pert else 0.0)
                picks[i].append((ta, b))
                if b == "P":
                    vec = n * float(n @ M @ n) / np.sqrt(CP)
                else:
                    u = M @ n
                    vec = (u - float(u @ n) * n) / np.sqrt(CS)
                if pert:
                    amp = np.linalg.norm(vec)
                    if amp > 0:
                        w = vec / amp + dn[i, bi]
                        vec = amp * site[i, bi] * w / np.linalg.norm(w)
                rad[i].append(float(np.linalg.norm(vec)))
                D[i] += SNR * A * np.outer(vec, ricker(T - ta)) / d[i]
    ref = np.asarray(sources[0], float) if len(sources) else np.zeros(3)
    scale = float(np.median(np.linalg.norm(STA - ref, axis=1)))
    for _ in range(n_clutter):
        i = int(rng.integers(ns))
        tc = float(rng.uniform(2, T[-1] - 4))
        q = rng.normal(size=3)
        q /= np.linalg.norm(q)
        D[i] += SNR * 0.35 * np.outer(q, ricker(T - tc)) / scale
        picks[i].append((tc, "P" if rng.random() < 0.5 else "S"))
    D += rng.normal(scale=SNR * 0.02, size=D.shape)
    return D, picks, [np.asarray(s, float) for s in sources], rad


def degrade_x(STA, srcs, picks, rad, rng):
    ns, k = len(STA), len(srcs)
    out, lab = [], []
    for i in range(ns):
        pl, ll = [], []
        for j, (t, b) in enumerate(picks[i]):
            if j < 2 * k:
                d = float(np.linalg.norm(STA[i] - srcs[j // 2]))
                A = (AP if b == "P" else AS) * max(rad[i][j], 1e-6)
                a = SNR * A / d / NOISE_SD
                if rng.random() > norm.cdf((a - 2.0) / 0.5):
                    continue
                crb = NOISE_SD / (SNR * A / d * np.sqrt(_GD2))
                pl.append((t + float(rng.normal(0, max(0.02, 3.0 * crb))), b))
                ll.append(j // 2 + 1)
            else:
                pl.append((t, b))
                ll.append(0)
        if not pl:
            pl, ll = [(-1000.0, "P")], [0]
        out.append(pl)
        lab.append(ll)
    return out, lab


def radiation_operator(src, STA):
    """G (2*ns*3 x 6): vector de desplazamiento P y S (en unidades relativas) por componente de m,
    m = (Mxx, Myy, Mzz, Mxy, Mxz, Myz). Incluye 1/d y las amplitudes relativas de P y S."""
    rows = []
    for b in ("P", "S"):
        for i in range(len(STA)):
            v = STA[i] - src
            d = np.linalg.norm(v)
            n = v / d
            nx, ny, nz = n
            if b == "P":
                c = np.array([nx * nx, ny * ny, nz * nz, 2 * nx * ny, 2 * nx * nz, 2 * ny * nz])
                g = np.outer(n, c) * AP / (d * np.sqrt(CP))
            else:
                B = np.array([[nx, 0, 0, ny, nz, 0], [0, ny, 0, nx, 0, nz], [0, 0, nz, 0, nx, ny]])
                g = (np.eye(3) - np.outer(n, n)) @ B * AS / (d * np.sqrt(CS))
            rows.append(g)
    return np.vstack(rows)


def basis(G):
    U, s, _ = np.linalg.svd(G, full_matrices=False)
    return U[:, s > 1e-8 * s.max()]


def node_vectors(MFn, STA, src, taus):
    """para cada tau: vector firmado de tres componentes en el pico de energia de la ventana de cada
    estacion y fase (2, ns, 3), su indice de tiempo, y los segmentos para WAVE-E y WAVE-PT."""
    ns = len(STA)
    st = np.arange(ns)
    d = np.linalg.norm(STA - src, axis=1)
    tt = np.stack([d / VP, d / VS])[None] + taus[:, None, None]                  # (nt, 2, ns)
    idx = np.clip((tt / DT).astype(int), 0, NT - 1)
    ii = np.clip(idx[..., None] + OFFS, 0, NT - 1)                               # (nt, 2, ns, W)
    seg = MFn[st[None, None, :, None], :, ii]                                     # (nt, 2, ns, W, 3)
    ener = np.sum(seg ** 2, axis=-1)                                              # (nt, 2, ns, W)
    k = ener.argmax(axis=-1)                                                      # (nt, 2, ns)
    v = np.take_along_axis(seg, k[..., None, None], axis=3)[:, :, :, 0, :]        # (nt, 2, ns, 3)
    nv = np.linalg.norm(v, axis=-1, keepdims=True)
    v = v * np.minimum(1.0, CAPV / np.maximum(nv, 1e-12))
    tpk = np.take_along_axis(ii, k[..., None], axis=3)[..., 0]                    # (nt, 2, ns)
    return v, tpk, seg, ener


def fields_x(STA, D, grid, taus, nodemap, Qs):
    ns = len(STA)
    MF, sd = matched(D)
    MFn = MF / sd
    out = {m: np.full(len(grid), -np.inf) for m in METH}
    for gi, src in enumerate(grid):
        v, _tpk, seg, ener = node_vectors(MFn, STA, src, taus)
        ray = np.array([polarization(src, STA[i], "P") for i in range(ns)])          # (ns, 3)
        pP = np.einsum("ic,tiwc->tiw", ray, seg[:, 0]) ** 2
        pS = ener[:, 1] - np.einsum("ic,tiwc->tiw", ray, seg[:, 1]) ** 2
        e_ = (np.minimum(ener.max(axis=-1) - 2.0, 36.0).sum(axis=(1, 2))) / 2.0
        pt = (np.minimum(pP.max(axis=-1), 36.0).sum(axis=1) + np.minimum(pS.max(axis=-1) - 1.0, 36.0).sum(axis=1)) / 2.0
        vf = v.reshape(len(taus), -1)                                                  # orden (b, i, c)
        mt = np.sum((vf @ Qs[gi]) ** 2, axis=1) / 2.0
        mtn = np.sum((vf @ Qs[nodemap[gi]]) ** 2, axis=1) / 2.0
        out["WAVE-E"][gi], out["WAVE-PT"][gi] = e_.max(), pt.max()
        out["WAVE-MT"][gi], out["WAVE-MT-NODE"][gi] = mt.max(), mtn.max()
    return out, MFn


def fw_record(STA, MFn, grid, taus, Gs, Qs, truth, kmax=5, min_sep=5.0):
    """Frank-Wolfe con deflacion: el maximo global del estadistico es el oraculo lineal; se ajusta m por
    minimos cuadrados, se resta su prediccion en el dominio del filtro adaptado y se repite. Los
    puntajes se hacen monotonos (minimo acumulado) para que el catalogo a un umbral sea un prefijo."""
    ns = len(STA)
    st = np.arange(ns)
    R = MFn.copy()
    garr = np.asarray(grid)
    acc, rec, run = [], [], np.inf
    for _ in range(kmax):
        best = (-np.inf, None, None, None)
        for gi, src in enumerate(garr):
            if any(np.linalg.norm(src[:2] - a[:2]) <= min_sep for a in acc):
                continue
            v, tpk, _seg, _e = node_vectors(R, STA, src, taus)
            vf = v.reshape(len(taus), -1)
            sc = np.sum((vf @ Qs[gi]) ** 2, axis=1) / 2.0
            j = int(sc.argmax())
            if sc[j] > best[0]:
                best = (float(sc[j]), gi, vf[j], tpk[j])
        s, gi, vf, tpk = best
        if gi is None:
            break
        src = garr[gi]
        run = min(run, s)
        hit = -1
        if truth:
            dd = [float(np.linalg.norm(src[:2] - np.asarray(t_)[:2])) for t_ in truth]
            j = int(np.argmin(dd))
            if dd[j] <= TOL_KM:
                hit = j
        rec.append((run, hit))
        acc.append(src)
        G = Gs[gi]
        m, *_ = np.linalg.lstsq(G, vf, rcond=None)
        pred = (G @ m).reshape(2, ns, 3)
        for b in range(2):
            for i in range(ns):
                t0 = int(tpk[b, i])
                lo, hi = max(0, t0 - AC_HALF), min(NT, t0 + AC_HALF + 1)
                R[i, :, lo:hi] -= np.outer(pred[b, i], AC[AC_HALF - (t0 - lo):AC_HALF + (hi - t0)])
    return rec


def run_cell(args):
    scen, reps, seed, cond, evbase, moff = args
    kind, nst, nclut = SCEN[scen]
    grid = default_grid()
    garr = np.asarray(grid)
    rng = np.random.default_rng(seed + 41 + evbase)
    rng_d = np.random.default_rng(seed + 43 + evbase)
    rng_m = np.random.default_rng(seed + moff)
    rng_h = np.random.default_rng(seed + 61 + evbase)
    rng_n = np.random.default_rng(seed + 71 + evbase)
    t0 = time.time()
    ev = []
    for k in range(reps):
        S = make_network(kind, nst, seed=seed + evbase + k)
        srcs, taus = draw_sources(rng)
        D, p, tr, rad = generate_x(S, srcs, taus, nclut, rng, rng_m, cond == "pert", rng_h)
        p2, lab = degrade_x(S, srcs, p, rad, rng_d)
        ev.append((S, D, p2, tr, lab))

    def nodemap():
        m = rng_n.permutation(len(grid))
        for gi in range(len(grid)):
            tries = 0
            while np.linalg.norm(garr[m[gi], :2] - garr[gi, :2]) < 9.0 and tries < 50:
                m[gi] = rng_n.integers(len(grid))
                tries += 1
        return m

    fev, recs_fw = [], []
    for S, D, _p, tr, _l in ev:
        Gs = [radiation_operator(g, S) for g in garr]
        Qs = [basis(G) for G in Gs]
        f, MFn = fields_x(S, D, grid, DEFAULT_TAUS, nodemap(), Qs)
        fev.append(f)
        recs_fw.append(fw_record(S, MFn, grid, DEFAULT_TAUS, Gs, Qs, tr))
    nte = sum(len(e[3]) for e in ev)
    res = {"scenario": scen, "seed": seed, "cond": cond, "k_eval": [len(e[3]) for e in ev], "methods": {}}
    for m in METH:
        recs = [greedy_record(grid, x[m], e[3]) for x, e in zip(fev, ev)]
        _t, H, F = exact_HF(recs)
        res["methods"][m] = {"H": H.tolist(), "F": F.tolist(),
                             "recall_at": {str(c): recall_at(H, F, nte, c) for c in CAPS}}
    _t, H, F = exact_HF(recs_fw)
    res["methods"]["MT-FW"] = {"H": H.tolist(), "F": F.tolist(), "records": recs_fw,
                               "recall_at": {str(c): recall_at(H, F, nte, c) for c in CAPS}}
    pts, per = [], {}
    for ci, cfg in enumerate(GAMMA_GRID):
        Hs, Fs = [], []
        for S, _D, p, tr, _l in ev:
            a, b = score_events(run_gamma(S, [[(t, bb) for (t, bb) in pl if t > -100] for pl in p], **cfg), tr)
            Hs.append(a)
            Fs.append(b)
        pts.append((sum(Fs) / len(ev), sum(Hs) / max(1, nte)))
        per[str(ci)] = (Hs, Fs)
    ra = {}
    for c in CAPS:
        ok = [(rc, ci) for ci, (fr, rc) in enumerate(pts) if fr <= c]
        ra[str(c)] = max(ok)[0] if ok else 0.0
    res["methods"]["GAMMA"] = {"points": pts, "per": per, "recall_at": ra}
    ntrue_arr = sum(2 * len(e[3]) * nst for e in ev)
    res["detected_fraction"] = float(sum(l >= 1 for e in ev for lab_i in e[4] for l in lab_i) / max(1, ntrue_arr))
    res["seconds"] = time.time() - t0
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--cond", default="ideal", choices=["ideal", "pert"])
    ap.add_argument("--evbase", type=int, default=15000)
    ap.add_argument("--moff", type=int, default=173)
    ap.add_argument("--workers", type=int, default=1)
    args = ap.parse_args()
    cells = [(s, args.reps, args.seed, args.cond, args.evbase, args.moff) for s in SCEN]
    out = []
    with ProcessPoolExecutor(args.workers) as ex:
        for r in ex.map(run_cell, cells):
            out.append(r)
            print(f"   hecho {r['scenario']}/{r['cond']} en {r['seconds']:.0f} s :: " + " | ".join(
                f"{m} " + "/".join(f"{r['methods'][m]['recall_at'][str(c)]:.2f}" for c in CAPS)
                for m in r["methods"]), flush=True)
    with open(f"task1x_seed{args.seed}_{args.cond}_r{args.reps}.json", "w") as fh:
        json.dump(out, fh)


if __name__ == "__main__":
    main()
