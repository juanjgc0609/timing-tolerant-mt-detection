"""
task1w_mech.py - Tarea 1w: fuentes con mecanismo de doble par y controles que rompen la geometria.

La verificacion independiente de 1r a 1u encontro dos cosas. (1) El generador daba a toda onda S
polarizacion SV pura, sin mecanismo de fuente; WAVE-P proyecta S sobre SV y con fuentes de doble par
al azar su ventaja sobre WAVE-E casi desaparece. (2) PHYS-SH no es un control nulo: baraja las razones
entre las ranuras estacion-fase del mismo nodo y conserva parte de la informacion de polarizacion.

Generador: el de 1r (consistente) con un doble par al azar por fuente; P radia n (n.M.n) y S la parte
transversal de M n, normalizadas a amplitud cuadratica media 1 sobre la esfera; el factor de radiacion
entra en la forma de onda, en la probabilidad de deteccion y en el error de tiempo del pick.
Metodos:
    KIN-LR, PHYS-LR          como en 1r (f/g estimada en 40 ventanas propias de este generador)
    PHYS-NODE                PHYS-LR con la direccion predicha por OTRO nodo de la grilla (permutacion
                             fija por ventana): mismos grados de libertad, sin la geometria del nodo
    WAVE-E, WAVE-P           como en 1t (WAVE-P proyecta S sobre SV)
    WAVE-PT                  independiente del mecanismo: P proyectada sobre el rayo, S como energia en
                             el plano transversal al rayo
    WAVE-PT-NODE             WAVE-PT con el rayo de otro nodo (control)
    GAMMA                    con semilla, rejilla de 1r, oraculo
Endpoint: curva FROC exacta, recall a falsos realizados <= 0.1 y <= 0.25.

    python task1w_mech.py --seed 2024 --reps 120 --workers 1
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

from task0_common import (AP, AS, DEFAULT_TAUS, DT, NT, SNR, T, VP, VS, _precompute, default_grid,
                          make_network, polarization, ricker)
from task1o_consistent import _GD2, NB, NOISE_SD, draw_sources, matched, pol_hist, score_events
from task1r_baselines import GAMMA_GRID, run_gamma
from task1r_fair import CAPS, EPSS, SCEN, SIGMAS, recall_at
from froc_exact import exact_HF, greedy_record


def rand_dc(r):
    s = r.normal(size=3)
    s /= np.linalg.norm(s)
    d = r.normal(size=3)
    d -= d.dot(s) * s
    d /= np.linalg.norm(d)
    return np.outer(s, d) + np.outer(d, s)


_r = np.random.default_rng(0)
_M = [rand_dc(_r) for _ in range(300)]
_n = _r.normal(size=(3000, 3))
_n /= np.linalg.norm(_n, axis=1)[:, None]
CP = float(np.mean([np.mean(np.einsum("ki,ij,kj->k", _n, M, _n) ** 2) for M in _M]))
CS = float(np.mean([np.mean(np.sum((np.einsum("ij,kj->ki", M, _n)
                                    - np.einsum("ki,ij,kj->k", _n, M, _n)[:, None] * _n) ** 2, 1)) for M in _M]))


def generate_dc(STA, sources, taus, n_clutter, rng, rng_m):
    """generate() de task0_common con radiacion de doble par; consume rng en el mismo orden."""
    ns = len(STA)
    D = np.zeros((ns, 3, NT))
    picks = [[] for _ in range(ns)]
    rad = [[] for _ in range(ns)]
    for src, tau in zip(sources, taus):
        d = np.linalg.norm(STA - src, axis=1)
        M = rand_dc(rng_m)
        for i in range(ns):
            n = (STA[i] - src) / d[i]
            for b, A, V in (("P", AP, VP), ("S", AS, VS)):
                ta = d[i] / V + tau
                picks[i].append((ta, b))
                if b == "P":
                    vec = n * float(n @ M @ n) / np.sqrt(CP)
                else:
                    u = M @ n
                    vec = (u - float(u @ n) * n) / np.sqrt(CS)
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


def degrade_dc(STA, srcs, picks, rad, rng):
    """degrade() consistente de task1o con la amplitud multiplicada por el factor de radiacion."""
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


PFAM = ("KIN-LR", "PHYS-LR", "PHYS-NODE")
WMETH = ("WAVE-E", "WAVE-P", "WAVE-PT", "WAVE-PT-NODE")


def fields_w(STA, D, picks, grid, taus, ratio, u, nodemap):
    ns = len(STA)
    PD, _HIT, NEAR = _precompute(D, picks)
    MF, sd = matched(D)
    MFn = MF / sd
    st = np.arange(ns)
    sh = int(round(0.25 / DT))
    offs = np.arange(-sh, sh + 1)
    ratio = np.asarray(ratio)
    garr = np.asarray(grid)
    pred_all = np.array([[[polarization(g, STA[i], bb) for i in range(ns)] for bb in ("P", "S")] for g in garr])
    out = {f"{f}@{s}@{e}": np.full(len(grid), -np.inf) for f in PFAM for s in SIGMAS for e in EPSS}
    out.update({m: np.full(len(grid), -np.inf) for m in WMETH})
    for gi, src in enumerate(garr):
        d = np.linalg.norm(STA - src, axis=1)
        pred = pred_all[gi]
        pnode = pred_all[nodemap[gi]]
        bk = np.full((3, len(SIGMAS), len(EPSS)), -np.inf)
        bw = [-np.inf] * 4
        for tau in taus:
            tt = np.stack([d / VP, d / VS]) + tau
            idx = np.clip((tt / DT).astype(int), 0, NT - 1)
            cl = NEAR[st[None, :], idx]
            res = np.abs(cl - tt)
            kpk = np.clip((cl / DT).astype(int), 0, NT - 1)
            meas = PD[st[None, :], kpk]
            cosv = np.abs(np.einsum("bic,bic->bi", pred, meas))
            cosn = np.abs(np.einsum("bic,bic->bi", pnode, meas))
            rr = ratio[np.clip((cosv * NB).astype(int), 0, NB - 1)]
            rn = ratio[np.clip((cosn * NB).astype(int), 0, NB - 1)]
            nk = (np.exp(-res[None] ** 2 / (2 * SIGMAS[:, None, None] ** 2))
                  / (np.sqrt(2 * np.pi) * SIGMAS[:, None, None]) / u)
            e = EPSS[None, :, None, None]
            for fi, mult in enumerate((1.0, rr, rn)):
                v = np.log(e + (1 - e) * (nk * mult)[:, None]).sum(axis=(2, 3))
                bk[fi] = np.maximum(bk[fi], v)
            ii = np.clip(idx[:, :, None] + offs[None, None, :], 0, NT - 1)
            seg = MFn[st[None, :, None], :, ii]                                  # (2, ns, W, 3)
            ener = np.sum(seg ** 2, axis=-1)
            proj = np.einsum("bic,biwc->biw", pred, seg) ** 2
            ray, rayn = pred[0], pnode[0]
            pP = np.einsum("ic,iwc->iw", ray, seg[0]) ** 2
            pS = ener[1] - np.einsum("ic,iwc->iw", ray, seg[1]) ** 2
            nP = np.einsum("ic,iwc->iw", rayn, seg[0]) ** 2
            nS = ener[1] - np.einsum("ic,iwc->iw", rayn, seg[1]) ** 2
            vals = (float(np.minimum(ener.max(axis=2) - 2.0, 36.0).sum()) / 2.0,
                    float(np.minimum(proj.max(axis=2), 36.0).sum()) / 2.0,
                    float(np.minimum(pP.max(axis=1), 36.0).sum() + np.minimum(pS.max(axis=1) - 1.0, 36.0).sum()) / 2.0,
                    float(np.minimum(nP.max(axis=1), 36.0).sum() + np.minimum(nS.max(axis=1) - 1.0, 36.0).sum()) / 2.0)
            bw = [max(a, b) for a, b in zip(bw, vals)]
        for fi, f in enumerate(PFAM):
            for si, s in enumerate(SIGMAS):
                for ei, ep in enumerate(EPSS):
                    out[f"{f}@{s}@{ep}"][gi] = bk[fi, si, ei]
        for m, v in zip(WMETH, bw):
            out[m][gi] = v
    return out


def froc(fl, data, grid):
    return exact_HF([greedy_record(grid, f, d[3]) for f, d in zip(fl, data)])


EVBASE, MOFF, WAVE_ONLY = 5000, 97, False


def run_cell(args):
    scen, reps, ncal, seed = args[:4]
    if len(args) > 4:
        global EVBASE, MOFF, WAVE_ONLY
        EVBASE, MOFF, WAVE_ONLY = args[4:]
    kind, nst, nclut = SCEN[scen]
    grid = default_grid()
    off = 0 if EVBASE == 5000 else EVBASE          # replica con ventanas nuevas: flujos aleatorios distintos
    rng = np.random.default_rng(seed + 41 + off)
    rng_d = np.random.default_rng(seed + 43 + off)
    rng_m = np.random.default_rng(seed + MOFF)
    rng_n = np.random.default_rng(seed + 71)
    t0 = time.time()

    def dataset(n, base):
        out = []
        for k in range(n):
            S = make_network(kind, nst, seed=seed + base + k)
            srcs, taus = draw_sources(rng)
            D, p, tr, rad = generate_dc(S, srcs, taus, nclut, rng, rng_m)
            p2, lab = degrade_dc(S, srcs, p, rad, rng_d)
            out.append((S, D, p2, tr, lab))
        return out

    if WAVE_ONLY:
        cal, ev = [], dataset(reps, EVBASE)
        ratio, u = np.ones(NB), 1.0
    else:
        cal, ev = dataset(ncal, 1000), dataset(reps, EVBASE)
        pol = dataset(40, 3000)
        ratio, _nf, _ng = pol_hist(pol, grid)
        u = max(1e-3, sum(sum(1 for l in lab_i if l == 0) for c in pol for lab_i in c[4])
                / (40 * nst * float(T[-1])))

    def nodemap():
        m = rng_n.permutation(len(grid))
        # direcciones de un nodo distinto y a mas de 9 km, para que el control no conserve la geometria
        garr = np.asarray(grid)
        for gi in range(len(grid)):
            tries = 0
            while np.linalg.norm(garr[m[gi], :2] - garr[gi, :2]) < 9.0 and tries < 50:
                m[gi] = rng_n.integers(len(grid))
                tries += 1
        return m

    fcal = [fields_w(S, D, p, grid, DEFAULT_TAUS, ratio, u, nodemap()) for S, D, p, _t, _l in cal]
    fev = [fields_w(S, D, p, grid, DEFAULT_TAUS, ratio, u, nodemap()) for S, D, p, _t, _l in ev]
    ntc, nte = sum(len(c[3]) for c in cal), sum(len(e[3]) for e in ev)
    res = {"scenario": scen, "seed": seed, "k_eval": [len(e[3]) for e in ev], "pol_ratio": list(ratio),
           "methods": {}}
    for f in ([] if WAVE_ONLY else PFAM):
        best = None
        for s in SIGMAS:
            for ep in EPSS:
                v = f"{f}@{s}@{ep}"
                _t, H, F = froc([x[v] for x in fcal], cal, grid)
                key = float(np.mean([recall_at(H, F, ntc, c) for c in (0.1, 0.25)]))
                if best is None or key > best[0]:
                    best = (key, v)
        _t, H, F = froc([x[best[1]] for x in fev], ev, grid)
        res["methods"][f] = {"variant": best[1], "H": H.tolist(), "F": F.tolist(),
                             "recall_at": {str(c): recall_at(H, F, nte, c) for c in CAPS}}
    for m in WMETH:
        _t, H, F = froc([x[m] for x in fev], ev, grid)
        res["methods"][m] = {"H": H.tolist(), "F": F.tolist(),
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
    ap.add_argument("--reps", type=int, default=120)
    ap.add_argument("--ncal", type=int, default=60)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--evbase", type=int, default=5000)
    ap.add_argument("--moff", type=int, default=97)
    ap.add_argument("--wave-only", action="store_true")
    args = ap.parse_args()
    cells = [(s, args.reps, args.ncal, args.seed, args.evbase, args.moff, args.wave_only) for s in SCEN]
    out = []
    with ProcessPoolExecutor(args.workers) as ex:
        for r in ex.map(run_cell, cells):
            out.append(r)
            print(f"   hecho {r['scenario']} en {r['seconds']:.0f} s :: " + " | ".join(
                f"{m} " + "/".join(f"{r['methods'][m]['recall_at'][str(c)]:.2f}" for c in CAPS)
                for m in r["methods"]), flush=True)
    tag = "" if (args.evbase, args.moff, args.wave_only) == (5000, 97, False) else f"_ev{args.evbase}_m{args.moff}" + ("_wave" if args.wave_only else "")
    with open(f"task1w_seed{args.seed}_r{args.reps}{tag}.json", "w") as fh:
        json.dump(out, fh)


if __name__ == "__main__":
    main()
