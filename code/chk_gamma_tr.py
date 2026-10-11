"""GaMMA pierde eventos por el filtro de tiempo de origen [5, 8] s? Misma medicion con el filtro de
tiempo ampliado a [0, 15] s (x, y, z sin cambios), 60 ventanas de ladder_6 y unilateral, semilla 2024."""
import os, sys, warnings, numpy as np
os.chdir("/home/claude/work"); sys.path.insert(0, "gamma_src"); warnings.filterwarnings("ignore")
import task1r_baselines as B
from task0_common import generate, make_network
from task1o_consistent import draw_sources, score_events
from task1r_fair import SCEN, degrade2
seed, NW = 2024, 60
for scen in ("ladder_6", "unilateral_c40"):
    kind, nst, nclut = SCEN[scen]
    rng = np.random.default_rng(seed + 41); rng_d = np.random.default_rng(seed + 43)
    def gen(base, k):
        S = make_network(kind, nst, seed=seed + base + k); srcs, taus = draw_sources(rng)
        D, p, tr = generate(S, srcs, taus, n_clutter=nclut, rng=rng); p2, _ = degrade2(S, srcs, p, "consistent", rng_d)
        return S, [[(t, b) for (t, b) in pl if t > -100] for pl in p2], tr
    for k in range(60): gen(1000, k)
    ev = [gen(5000, k) for k in range(NW)]
    nt = sum(len(e[2]) for e in ev)
    for TR in ((5.0, 8.0), (0.0, 15.0)):
        B.TR = TR
        best = []
        for cfg in B.GAMMA_GRID:
            H = F = 0
            for S, p, tr in ev:
                a, b = score_events(B.run_gamma(S, p, **cfg), tr); H += a; F += b
            best.append((F / NW, H / nt))
        print(scen, "TR", TR, "max recall", round(max(r for _, r in best), 3), "puntos", [(round(f, 2), round(r, 2)) for f, r in best], flush=True)
