"""Sensibilidad: la ventaja de MT-FW sobre WAVE-MT se debe solo a suprimir lobulos laterales? Se
compara con WAVE-MT y WAVE-PT con radio de supresion de 5, 12 y 20 km (12 km es la separacion minima
entre fuentes del generador: un conocimiento que el metodo no deberia tener). Mismas ventanas que 1x."""
import sys, numpy as np
import task1x_mtglrt as X
from task0_common import default_grid, make_network, DEFAULT_TAUS
from task1o_consistent import draw_sources
from task1r_fair import SCEN, recall_at
from froc_exact import exact_HF, greedy_record
scen, cond, seed, NW = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
kind, nst, nclut = SCEN[scen]; grid = default_grid(); garr = np.asarray(grid); eb = 15000
rng = np.random.default_rng(seed + 41 + eb); rm = np.random.default_rng(seed + 173); rh = np.random.default_rng(seed + 61 + eb)
rng_d = np.random.default_rng(seed + 43 + eb)
R = {k: [] for k in ("MT5", "MT12", "MT20", "PT5", "PT12", "PT20", "FW")}; nt = 0
for k in range(NW):
    S = make_network(kind, nst, seed=seed + eb + k); srcs, taus = draw_sources(rng)
    D, p, tr, rad = X.generate_x(S, srcs, taus, nclut, rng, rm, cond == "pert", rh); X.degrade_x(S, srcs, p, rad, rng_d)
    nt += len(tr)
    Gs = [X.radiation_operator(g, S) for g in garr]; Qs = [X.basis(G) for G in Gs]
    f, MFn = X.fields_x(S, D, grid, DEFAULT_TAUS, np.arange(len(grid)), Qs)
    for r_ in (5, 12, 20):
        R[f"MT{r_}"].append(greedy_record(grid, f["WAVE-MT"], tr, min_sep_km=float(r_)))
        R[f"PT{r_}"].append(greedy_record(grid, f["WAVE-PT"], tr, min_sep_km=float(r_)))
    R["FW"].append(X.fw_record(S, MFn, grid, DEFAULT_TAUS, Gs, Qs, tr))
for m, recs in R.items():
    _t, H, F = exact_HF(recs)
    print(scen, cond, seed, m, [round(recall_at(H, F, nt, c), 2) for c in (0.1, 0.25)], flush=True)
