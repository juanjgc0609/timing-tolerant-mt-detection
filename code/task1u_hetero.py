"""
task1u_hetero.py - Tarea 1u: el positivo de WAVE-P puede ser una propiedad de un medio homogeneo
conocido exactamente. Aqui el generador anade lo que un medio real anade y el modelo de inferencia no
predice: un retardo por estacion y fase N(0, 0.15 s) y una desviacion de la direccion de polarizacion
de unos 10 grados, fijos en cada ventana. Los picks se degradan como en 1r a partir de los tiempos
perturbados. WAVE-E, WAVE-P, WAVE-PSH y GaMMA (con semilla, rejilla de 1r, oraculo) sobre las mismas
ventanas; las fuentes y el clutter son los de las primeras ventanas de 1r (flujo aleatorio aparte).

    python task1u_hetero.py --seed 2024 --reps 100
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

from task0_common import DEFAULT_TAUS, default_grid, generate, make_network
from task1o_consistent import draw_sources, score_events
from task1r_baselines import GAMMA_GRID, run_gamma
from task1r_fair import CAPS, COND, SCEN, degrade2, recall_at, sweep
from task1t_wave import METH, wave_fields
from froc_exact import exact_HF, greedy_record

WAVE_ONLY = False
STATIC_SD = 0.15
DIR_SD = float(np.tan(np.radians(10.0)))


def run_cell(args):
    scen, cond, reps, ncal, seed = args
    kind, nst, nclut = SCEN[scen]
    grid = default_grid()
    rng = np.random.default_rng(seed + 41)
    rng_d = np.random.default_rng(seed + 43)
    rng_sh = np.random.default_rng(seed + 59)
    rng_h = np.random.default_rng(seed + 61)
    t0 = time.time()

    def dataset(n, base, keep):
        out = []
        for k in range(n):
            S = make_network(kind, nst, seed=seed + base + k)
            srcs, taus = draw_sources(rng)
            if keep:
                st = rng_h.normal(0.0, STATIC_SD, size=(nst, 2))
                dn = rng_h.normal(0.0, DIR_SD, size=(nst, 2, 3))
                D, p, tr = generate(S, srcs, taus, n_clutter=nclut, rng=rng, statics=st, dir_noise=dn)
            else:
                D, p, tr = generate(S, srcs, taus, n_clutter=nclut, rng=rng)
            p2, lab = degrade2(S, srcs, p, COND[cond], rng_d)
            if keep:
                out.append((S, D, p2, tr, lab))
        return out

    dataset(ncal, 1000, False)
    ev = dataset(reps, 5000, True)
    nte = sum(len(e[3]) for e in ev)
    res = {"scenario": scen, "cond": cond, "seed": seed, "k_eval": [len(e[3]) for e in ev], "methods": {},
           "static_sd": STATIC_SD, "dir_deg": 10.0}
    fev = [wave_fields(S, D, grid, DEFAULT_TAUS, rng_sh) for S, D, _p, _t, _l in ev]
    for m in METH:
        _thr, Hc, Fc = sweep([x[m] for x in fev], ev, grid)
        recs = [greedy_record(grid, x[m], e[3]) for x, e in zip(fev, ev)]
        _t, H, F = exact_HF(recs)
        res["methods"][m] = {"H": H.tolist(), "F": F.tolist(), "records": recs,
                             "recall_at": {str(c): recall_at(H, F, nte, c) for c in CAPS},
                             "recall_at_coarse": {str(c): recall_at(Hc, Fc, nte, c) for c in CAPS}}
    if WAVE_ONLY:
        res["seconds"] = time.time() - t0
        return res
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
    res["seconds"] = time.time() - t0
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--ncal", type=int, default=60)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--wave-only", action="store_true")
    args = ap.parse_args()
    global WAVE_ONLY
    WAVE_ONLY = args.wave_only
    cells = [(s, "consistente", args.reps, args.ncal, args.seed) for s in SCEN]
    out = []
    with ProcessPoolExecutor(args.workers) as ex:
        for r in ex.map(run_cell, cells):
            out.append(r)
            print(f"   hecho {r['scenario']} en {r['seconds']:.0f} s :: " + " | ".join(
                f"{m} " + "/".join(f"{r['methods'][m]['recall_at'][str(c)]:.2f}" for c in CAPS)
                for m in r["methods"]), flush=True)
    with open(f"task1u_seed{args.seed}_r{args.reps}" + ("_wave_exact" if WAVE_ONLY else "") + ".json", "w") as fh:
        json.dump(out, fh)


if __name__ == "__main__":
    main()
