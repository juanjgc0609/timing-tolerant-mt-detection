"""
task1o_baselines.py - GaMMA y PyOcto sobre EXACTAMENTE las ventanas de la Tarea 1o (se regeneran con
las mismas semillas y el mismo orden de numeros aleatorios), con el mismo endpoint: umbral o
configuracion elegidos en calibracion para cada presupuesto de falsos por ventana, recall por evento
en evaluacion, y comparacion pareada por evento con PHYS-LR y KIN-LR al presupuesto 0.1.

    python task1o_baselines.py --seed 2024 --workers 2
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

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "gamma_src"))
warnings.filterwarnings("ignore")

from task0_common import generate, make_network
from task1jk_controls import GAMMA_GRID, PYOCTO_GRID, run_gamma, run_pyocto
from task1o_consistent import COND, SCEN, TOL_KM, degrade, draw_sources, score_events

BUDGETS = (0.05, 0.1, 0.25, 0.5)


def datasets(scen, cond, reps, ncal, seed):
    kind, nst, nclut = SCEN[scen]
    rng = np.random.default_rng(seed + 41)
    rng_d = np.random.default_rng(seed + 43)

    def dataset(n, base):
        out = []
        for k in range(n):
            S = make_network(kind, nst, seed=seed + base + k)
            srcs, taus = draw_sources(rng)
            D, p, tr = generate(S, srcs, taus, n_clutter=nclut, rng=rng)
            p2, _lab = degrade(S, srcs, p, COND[cond], rng_d)
            p2 = [[(t, b) for (t, b) in pl if t > -100] for pl in p2]
            out.append((S, p2, tr))
        return out

    return dataset(ncal, 1000), dataset(reps, 5000)


def run_cell(args):
    scen, cond, reps, ncal, seed = args
    t0 = time.time()
    cal, ev = datasets(scen, cond, reps, ncal, seed)
    ntc = sum(len(c[2]) for c in cal)
    res = {"scenario": scen, "cond": cond, "seed": seed, "k_eval": [len(e[2]) for e in ev], "methods": {}}
    for name, runner, grid in (("GAMMA", run_gamma, GAMMA_GRID), ("PYOCTO", run_pyocto, PYOCTO_GRID)):
        stats = []
        for cfg in grid:
            h = nf = 0
            for S, p, tr in cal:
                a, b = score_events(runner(S, p, **cfg), tr)
                h += a
                nf += b
            stats.append((cfg, h / max(1, ntc), nf / len(cal)))
        entry = {"budgets": list(BUDGETS)}
        for bgt in BUDGETS:
            ok = [s for s in stats if s[2] <= bgt]
            cfg = max(ok, key=lambda s: s[1])[0] if ok else None
            hits, fal, det = [], [], []
            for S, p, tr in ev:
                evs = runner(S, p, **cfg) if cfg is not None else []
                a, b = score_events(evs, tr)
                hits.append(a)
                fal.append(b)
                for s_ in tr:
                    det.append(int(any(np.linalg.norm(np.asarray(x)[:2] - s_[:2]) <= TOL_KM for x in evs)))
            entry[f"b{bgt}"] = {"config": cfg, "hits": hits, "false": fal, "det": det}
        res["methods"][name] = entry
    res["seconds"] = time.time() - t0
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--reps", type=int, default=60)
    ap.add_argument("--ncal", type=int, default=60)
    ap.add_argument("--workers", type=int, default=2)
    args = ap.parse_args()
    conds = os.environ.get("TASK1OB_CONDS", "consistente").split(",")
    cells = [(s, c, args.reps, args.ncal, args.seed) for c in conds for s in SCEN]
    out = []
    with ProcessPoolExecutor(args.workers) as ex:
        for r in ex.map(run_cell, cells):
            out.append(r)
            print(f"   hecho {r['scenario']}/{r['cond']} en {r['seconds']:.0f} s", flush=True)
    with open(f"task1o_baselines_seed{args.seed}.json", "w") as fh:
        json.dump(out, fh)


if __name__ == "__main__":
    main()
