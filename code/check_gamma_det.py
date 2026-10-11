"""check_gamma_det.py - GaMMA depende del RNG global de numpy: se mide su variabilidad sin semilla y se
verifica que, con semilla fija, las dos rutas de codigo (task1r_baselines y task1s) dan el mismo catalogo."""
import os, sys, warnings
os.environ.setdefault("OMP_NUM_THREADS", "1")
warnings.filterwarnings("ignore")
import numpy as np
sys.path.insert(0, "gamma_src")
import task1r_baselines as B
import task1s_gamma_pol as S
from task0_common import generate, make_network
from task1o_consistent import draw_sources, score_events
from task1r_fair import COND, SCEN, degrade2

seed, scen, n = 2024, "unilateral_c40", 8
kind, nst, nclut = SCEN[scen]
rng = np.random.default_rng(seed + 41); rng_d = np.random.default_rng(seed + 43)
W = []
for k in range(60):                      # calibracion (se descarta), mismo orden que los scripts
    S_ = make_network(kind, nst, seed=seed + 1000 + k); s, t = draw_sources(rng)
    D, p, tr = generate(S_, s, t, n_clutter=nclut, rng=rng); degrade2(S_, s, p, COND["consistente"], rng_d)
for k in range(n):
    S_ = make_network(kind, nst, seed=seed + 5000 + k); s, t = draw_sources(rng)
    D, p, tr = generate(S_, s, t, n_clutter=nclut, rng=rng)
    p2, _ = degrade2(S_, s, p, COND["consistente"], rng_d)
    W.append((S_, [[(a, b) for (a, b) in pl if a > -100] for pl in p2], tr))

src = open("task1r_baselines.py").read()
unseeded = src.replace('        np.random.seed(42)  # GaMMA inicializa con el RNG global de numpy\n', '')
ns = {"__file__": os.path.abspath("task1r_baselines.py"), "__name__": "unseeded"}
exec(compile(unseeded, "unseeded", "exec"), ns)
diff_unseeded = diff_paths = diff_repeat = 0
tot = 0
for cfg in B.GAMMA_GRID:
    for S_, p, tr in W:
        tot += 1
        np.random.seed(1); a = ns["run_gamma"](S_, p, **cfg)
        np.random.seed(2); b = ns["run_gamma"](S_, p, **cfg)
        diff_unseeded += score_events(a, tr) != score_events(b, tr)
        c = B.run_gamma(S_, p, **cfg); c2 = B.run_gamma(S_, p, **cfg)
        d = [x for x, _m in S.gamma_with_assign(S_, p, **cfg)[0]]
        same = len(c) == len(d) and all(np.allclose(x, y) for x, y in zip(c, d))
        diff_paths += not same
        diff_repeat += not (len(c) == len(c2) and all(np.allclose(x, y) for x, y in zip(c, c2)))
print(f"llamadas {tot}: sin semilla, (aciertos, falsos) distintos entre dos estados del RNG: {diff_unseeded}; "
      f"con semilla, catalogo distinto entre rutas: {diff_paths}; entre repeticiones: {diff_repeat}")
