"""Calibracion del tope por estacion CAPV del GLRT de tensor de momento, en ventanas propias
(evbase 20000, semilla 2024, condicion ideal y pert), antes de la corrida de evaluacion."""
import sys, numpy as np, json
from concurrent.futures import ProcessPoolExecutor
import task1x_mtglrt as X
from task0_common import default_grid, make_network, DEFAULT_TAUS
from task1o_consistent import matched, draw_sources
from task1r_fair import SCEN, recall_at
from froc_exact import exact_HF, greedy_record
CAPS_TRY = (6.0, 12.0, 24.0, 1e9)
def run(args):
    scen, cond = args
    kind, nst, nclut = SCEN[scen]; grid = default_grid(); garr = np.asarray(grid)
    seed, eb = 2024, 20000
    rng = np.random.default_rng(seed + 41 + eb); rm = np.random.default_rng(seed + 211); rh = np.random.default_rng(seed + 61 + eb)
    recs = {c: [] for c in CAPS_TRY}; ntrue = 0
    for k in range(40):
        S = make_network(kind, nst, seed=seed + eb + k); srcs, taus = draw_sources(rng)
        D, p, tr, rad = X.generate_x(S, srcs, taus, nclut, rng, rm, cond == "pert", rh); ntrue += len(tr)
        MF, sd = matched(D); MFn = MF / sd
        Qs = [X.basis(X.radiation_operator(g, S)) for g in garr]
        for c in CAPS_TRY:
            X.CAPV = c
            f = np.array([np.max(np.sum((X.node_vectors(MFn, S, g, DEFAULT_TAUS)[0].reshape(len(DEFAULT_TAUS), -1) @ Qs[gi]) ** 2, axis=1)) for gi, g in enumerate(garr)])
            recs[c].append(greedy_record(grid, f, tr))
    out = {}
    for c in CAPS_TRY:
        _t, H, F = exact_HF(recs[c]); out[c] = [recall_at(H, F, ntrue, cap) for cap in (0.1, 0.25)]
    return scen, cond, out
cells = [(s, c) for c in ("ideal", "pert") for s in SCEN]
with ProcessPoolExecutor(2) as ex:
    for scen, cond, out in ex.map(run, cells):
        print(scen, cond, {k: [round(x, 2) for x in v] for k, v in out.items()}, flush=True)
