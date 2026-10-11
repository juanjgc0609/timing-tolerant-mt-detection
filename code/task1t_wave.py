"""
task1t_wave.py - Tarea 1t: la fisica sin picks. El canal de polarizacion no puede recuperar los arribos
que el picker pierde, porque solo reordena lo que el picker entrega. Un puntaje que migra la forma de
onda de tres componentes sobre la hipotesis de fuente no depende del picker. Tres puntajes identicos
salvo en la polarizacion:
    WAVE-E    energia del filtro adaptado de las tres componentes en el tiempo predicho (familia de
              apilamiento por coherencia: QuakeMigrate, LOKI, MALMI), sin direccion
    WAVE-P    la misma energia proyectada sobre la direccion de polarizacion predicha por la hipotesis
    WAVE-PSH  control negativo: proyectada sobre la direccion predicha para OTRA estacion (permutacion
              fija por ventana), que conserva la distribucion de direcciones y rompe la geometria
Mismas ventanas, grilla y endpoint (FROC, recall a falsos realizados) que la Tarea 1r.

    python task1t_wave.py --seed 2024 --reps 200 --workers 1
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from task0_common import DEFAULT_TAUS, DT, NT, VP, VS, default_grid, generate, make_network, polarization
from task1o_consistent import draw_sources, matched
from task1r_fair import CAPS, COND, SCEN, degrade2, recall_at, sweep
from froc_exact import exact_HF, greedy_record

METH = ("WAVE-E", "WAVE-P", "WAVE-PSH")


def wave_fields(STA, D, grid, taus, rng_sh):
    ns = len(STA)
    MF, sd = matched(D)
    MFn = MF / sd
    st = np.arange(ns)
    sh = int(round(0.25 / DT))
    offs = np.arange(-sh, sh + 1)
    perm = rng_sh.permutation(ns)
    while ns > 1 and np.any(perm == st):
        perm = rng_sh.permutation(ns)
    out = {m: np.full(len(grid), -np.inf) for m in METH}
    for gi, src in enumerate(grid):
        d = np.linalg.norm(STA - src, axis=1)
        pred = np.array([[polarization(src, STA[i], bb) for i in range(ns)] for bb in ("P", "S")])
        psh = pred[:, perm, :]
        be = bp = bs = -np.inf
        for tau in taus:
            tt = np.stack([d / VP, d / VS]) + tau
            idx = np.clip((tt / DT).astype(int), 0, NT - 1)
            ii = np.clip(idx[:, :, None] + offs[None, None, :], 0, NT - 1)
            seg = MFn[st[None, :, None], :, ii]                              # (2, ns, W, 3)
            ener = np.sum(seg ** 2, axis=-1)
            proj = np.einsum("bic,biwc->biw", pred, seg) ** 2
            prsh = np.einsum("bic,biwc->biw", psh, seg) ** 2
            be = max(be, float(np.minimum(ener.max(axis=2) - 2.0, 36.0).sum()) / 2.0)
            bp = max(bp, float(np.minimum(proj.max(axis=2), 36.0).sum()) / 2.0)
            bs = max(bs, float(np.minimum(prsh.max(axis=2), 36.0).sum()) / 2.0)
        out["WAVE-E"][gi], out["WAVE-P"][gi], out["WAVE-PSH"][gi] = be, bp, bs
    return out


def run_cell(args):
    scen, cond, reps, ncal, seed = args
    kind, nst, nclut = SCEN[scen]
    grid = default_grid()
    rng = np.random.default_rng(seed + 41)
    rng_d = np.random.default_rng(seed + 43)
    rng_sh = np.random.default_rng(seed + 59)
    t0 = time.time()

    def dataset(n, base, keep):
        out = []
        for k in range(n):
            S = make_network(kind, nst, seed=seed + base + k)
            srcs, taus = draw_sources(rng)
            D, p, tr = generate(S, srcs, taus, n_clutter=nclut, rng=rng)
            p2, lab = degrade2(S, srcs, p, COND[cond], rng_d)
            if keep:
                out.append((S, D, p2, tr, lab))
        return out

    dataset(ncal, 1000, False)                 # consume el RNG como en 1r: mismas ventanas de evaluacion
    ev = dataset(reps, 5000, True)
    fev = [wave_fields(S, D, grid, DEFAULT_TAUS, rng_sh) for S, D, _p, _t, _l in ev]
    nte = sum(len(e[3]) for e in ev)
    res = {"scenario": scen, "cond": cond, "seed": seed, "k_eval": [len(e[3]) for e in ev], "methods": {}}
    for m in METH:
        _thr, Hc, Fc = sweep([x[m] for x in fev], ev, grid)
        recs = [greedy_record(grid, x[m], e[3]) for x, e in zip(fev, ev)]
        _t, H, F = exact_HF(recs)
        res["methods"][m] = {"H": H.tolist(), "F": F.tolist(), "records": recs,
                             "recall_at": {str(c): recall_at(H, F, nte, c) for c in CAPS},
                             "recall_at_coarse": {str(c): recall_at(Hc, Fc, nte, c) for c in CAPS}}
    res["seconds"] = time.time() - t0
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--ncal", type=int, default=60)
    ap.add_argument("--workers", type=int, default=1)
    args = ap.parse_args()
    cells = [(s, "consistente", args.reps, args.ncal, args.seed) for s in SCEN]
    out = []
    with ProcessPoolExecutor(args.workers) as ex:
        for r in ex.map(run_cell, cells):
            out.append(r)
            print(f"   hecho {r['scenario']} en {r['seconds']:.0f} s :: " + " | ".join(
                f"{m} " + "/".join(f"{r['methods'][m]['recall_at'][str(c)]:.2f}" for c in CAPS)
                for m in METH), flush=True)
    with open(f"task1t_seed{args.seed}_r{args.reps}_exact.json", "w") as fh:
        json.dump(out, fh)


if __name__ == "__main__":
    main()
