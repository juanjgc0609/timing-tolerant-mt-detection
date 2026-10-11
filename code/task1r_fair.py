"""
task1r_fair.py - Tarea 1r: la comparacion de 1o repetida con las correcciones de la verificacion
independiente, para que un resultado positivo no dependa de ninguna ventaja de implementacion.

Cambios respecto de 1o:
  1. Endpoint sin sobrepaso de presupuesto: curva FROC en evaluacion (barrido de umbral), recall a
     tasa de falsos REALIZADA <= 0.05, 0.1 y 0.25 por ventana, e IC95 por bootstrap pareado de
     ventanas para la diferencia de recall. La variante de cada metodo se elige en calibracion.
  2. Control negativo PHYS-SH: la misma razon de polarizacion aplicada a cosenos permutados al azar
     entre arribos. Debe quedar igual a KIN-LR; si mejora, la ganancia no viene de la polarizacion.
  3. Condicion "exacto" (todos los arribos detectados, tiempo exacto): separa el efecto del
     generador del efecto de la forma de verosimilitud (PHYS-LR contra KIN-LR en el generador viejo).
  4. Esfuerzo simetrico: KIN-LR, PHYS-LR y PHYS-SH con la misma rejilla sigma x eps.
  5. f/g estimada en 40 ventanas propias, separadas de las de calibracion de umbral y variante.

    python task1r_fair.py --seed 2024 --workers 2
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from task0_common import (DEFAULT_TAUS, DT, NT, T, TOL, VMIN, VP, VS, _catalog_impl, _precompute,
                          default_grid, generate, make_network, polarization)
from task1o_consistent import NB, degrade, draw_sources, matched, pol_hist, score_events

SIGMAS = np.array([0.1, 0.2, 0.4])
EPSS = np.array([0.05, 0.1, 0.2])
SCEN = {"unilateral_c40": ("onesided", 8, 40), "ladder_4": ("random", 4, 35),
        "ladder_6": ("random", 6, 35)}
COND = {"consistente": "consistent", "exacto": "exact"}
FAM = ("KIN-LR", "PHYS-LR", "PHYS-SH")
VAR = [f"{f}@{s}@{e}" for f in FAM for s in SIGMAS for e in EPSS]
OWN = ["TOL-TB", "WAVE-P"] + VAR
CAPS = (0.05, 0.1, 0.25)


def degrade2(STA, srcs, picks, cond, rng):
    if cond == "exact":
        k = len(srcs)
        lab = [[(j // 2 + 1) if j < 2 * k else 0 for j in range(len(pl))] for pl in picks]
        return [pl if pl else [(-1000.0, "P")] for pl in picks], lab
    return degrade(STA, srcs, picks, cond, rng)


def fields(STA, D, picks, grid, taus, ratio, u, rng_sh):
    ns = len(STA)
    PD, HIT, NEAR = _precompute(D, picks)
    MF, sd = matched(D)
    MFn = MF / sd
    DIJ = np.linalg.norm(STA[:, None, :] - STA[None, :, :], axis=2)
    st = np.arange(ns)
    ratio = np.asarray(ratio)
    sh = int(round(0.25 / DT))
    offs = np.arange(-sh, sh + 1)
    out = {m: np.full(len(grid), -np.inf) for m in OWN}
    perm = rng_sh.permutation(2 * ns)
    for gi, src in enumerate(grid):
        d = np.linalg.norm(STA - src, axis=1)
        pred = np.array([[polarization(src, STA[i], bb) for i in range(ns)] for bb in ("P", "S")])
        bk = np.full((3, len(SIGMAS), len(EPSS)), -np.inf)
        bt = -np.inf
        bw = -np.inf
        for tau in taus:
            tt = np.stack([d / VP, d / VS]) + tau
            idx = np.clip((tt / DT).astype(int), 0, NT - 1)
            cl = NEAR[st[None, :], idx]
            res = np.abs(cl - tt)
            kpk = np.clip((cl / DT).astype(int), 0, NT - 1)
            cosv = np.abs(np.einsum("bic,bic->bi", pred, PD[st[None, :], kpk]))
            rr = ratio[np.clip((cosv * NB).astype(int), 0, NB - 1)]
            rsh = rr.reshape(-1)[perm].reshape(rr.shape)
            nk = (np.exp(-res[None] ** 2 / (2 * SIGMAS[:, None, None] ** 2))
                  / (np.sqrt(2 * np.pi) * SIGMAS[:, None, None]) / u)            # (S, 2, ns)
            e = EPSS[None, :, None, None]
            for fi, mult in enumerate((1.0, rr, rsh)):
                v = np.log(e + (1 - e) * (nk * mult)[:, None]).sum(axis=(2, 3))  # (S, E)
                bk[fi] = np.maximum(bk[fi], v)
            ii = np.clip(idx[:, :, None] + offs[None, None, :], 0, NT - 1)
            seg = MFn[st[None, :, None], :, ii]
            proj = np.einsum("bic,biwc->biw", pred, seg) ** 2
            bw = max(bw, float(np.minimum(proj.max(axis=2), 36.0).sum()) / 2.0)
            hit = HIT[st[None, :], idx]
            if hit.sum() == 0:
                continue
            s_hit = hit.sum() / (2 * ns)
            bt = max(bt, s_hit - 0.01 * float(np.sqrt(np.mean(res[hit] ** 2))) / TOL)
        out["TOL-TB"][gi] = bt
        out["WAVE-P"][gi] = bw
        for fi, f in enumerate(FAM):
            for si, s in enumerate(SIGMAS):
                for ei, ep in enumerate(EPSS):
                    out[f"{f}@{s}@{ep}"][gi] = bk[fi, si, ei]
    return out


def sweep(fl, data, grid, nthr=80):
    """matriz por ventana y umbral de (eventos verdaderos recuperados, falsos)."""
    vals = np.concatenate([f[np.isfinite(f)] for f in fl])
    thr = np.unique(np.quantile(vals, np.linspace(0.85, 1.0, nthr)))
    H = np.zeros((len(fl), len(thr)), int)
    F = np.zeros((len(fl), len(thr)), int)
    for w, (f, d) in enumerate(zip(fl, data)):
        ff = np.where(np.isfinite(f), f, -1e9)
        for j, t in enumerate(thr):
            a, b = score_events(_catalog_impl(grid, ff, float(t)), d[3])
            H[w, j], F[w, j] = a, b
    return thr, H, F


EXACT = False


def sweep_any(fl, data, grid, nthr=80):
    """barrido por cuantiles (registrado) o curva FROC exacta (froc_exact.py) segun EXACT."""
    if not EXACT:
        return sweep(fl, data, grid, nthr)
    from froc_exact import exact_HF, greedy_record
    return exact_HF([greedy_record(grid, f, d[3]) for f, d in zip(fl, data)])


def recall_at(H, F, ntrue, cap):
    rec = H.sum(0) / max(1, ntrue)
    fr = F.sum(0) / H.shape[0]
    ok = fr <= cap
    return float(rec[ok].max()) if ok.any() else 0.0


def run_cell(args):
    scen, cond, reps, ncal, seed = args
    kind, nst, nclut = SCEN[scen]
    grid = default_grid()
    rng = np.random.default_rng(seed + 41)
    rng_d = np.random.default_rng(seed + 43)
    rng_sh = np.random.default_rng(seed + 47)
    t0 = time.time()

    def dataset(n, base):
        out = []
        for k in range(n):
            S = make_network(kind, nst, seed=seed + base + k)
            srcs, taus = draw_sources(rng)
            D, p, tr = generate(S, srcs, taus, n_clutter=nclut, rng=rng)
            p2, lab = degrade2(S, srcs, p, COND[cond], rng_d)
            out.append((S, D, p2, tr, lab))
        return out

    cal, ev = dataset(ncal, 1000), dataset(reps, 5000)
    pol = dataset(40, 3000)
    ratio, nf_, ng_ = pol_hist(pol, grid)
    u = max(1e-3, sum(sum(1 for l in lab_i if l == 0) for c in pol for lab_i in c[4])
            / (40 * nst * float(T[-1])))
    fcal = [fields(S, D, p, grid, DEFAULT_TAUS, ratio, u, rng_sh) for S, D, p, _t, _l in cal]
    fev = [fields(S, D, p, grid, DEFAULT_TAUS, ratio, u, rng_sh) for S, D, p, _t, _l in ev]
    ntc, nte = sum(len(c[3]) for c in cal), sum(len(e[3]) for e in ev)
    res = {"scenario": scen, "cond": cond, "seed": seed, "k_eval": [len(e[3]) for e in ev],
           "pol_ratio": ratio, "u": u, "methods": {}}
    groups = {"TOL-TB": ["TOL-TB"], "WAVE-P": ["WAVE-P"]}
    for f in FAM:
        groups[f] = [v for v in VAR if v.startswith(f + "@")]
    for g, variants in groups.items():
        best = None
        for v in variants:
            _thr, H, F = sweep_any([x[v] for x in fcal], cal, grid, nthr=40)
            key = float(np.mean([recall_at(H, F, ntc, c) for c in (0.1, 0.25)]))
            if best is None or key > best[0]:
                best = (key, v)
        v = best[1]
        thr, H, F = sweep_any([x[v] for x in fev], ev, grid)
        res["methods"][g] = {"variant": v, "H": H.tolist(), "F": F.tolist(),
                             "recall_at": {str(c): recall_at(H, F, nte, c) for c in CAPS}}
    res["seconds"] = time.time() - t0
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--reps", type=int, default=60)
    ap.add_argument("--ncal", type=int, default=60)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--cond", default="all")
    ap.add_argument("--exact", action="store_true")
    args = ap.parse_args()
    global EXACT
    EXACT = args.exact
    conds = list(COND) if args.cond == "all" else [args.cond]
    cells = [(s, c, args.reps, args.ncal, args.seed) for c in conds for s in SCEN]
    out = []
    with ProcessPoolExecutor(args.workers) as ex:
        for r in ex.map(run_cell, cells):
            out.append(r)
            ra = {m: r["methods"][m]["recall_at"] for m in r["methods"]}
            print(f"   hecho {r['scenario']}/{r['cond']} en {r['seconds']:.0f} s :: "
                  + " | ".join(f"{m} " + "/".join(f"{ra[m][str(c)]:.2f}" for c in CAPS) for m in ra),
                  flush=True)
    with open(f"task1r_seed{args.seed}" + ("" if args.reps == 60 else f"_r{args.reps}")
              + ("_exact" if EXACT else "") + ".json", "w") as fh:
        json.dump(out, fh)


if __name__ == "__main__":
    main()
