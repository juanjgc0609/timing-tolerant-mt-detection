"""Escalamiento de la poda certificada con la resolucion de la grilla (misma verificacion de exactitud)."""
import time
import numpy as np
import bnb_octree as B
from task0_common import _catalog_impl, generate, make_network
from task1o_consistent import degrade, draw_sources, score_events

for step, tstep in ((3.0, 0.1), (1.5, 0.05), (0.75, 0.05)):
    B.AX = np.arange(-21.0, 21.0 + 1e-9, step) + 0.37
    B.TAUS = np.arange(5.0, 8.0 + 1e-9, tstep)
    rng = np.random.default_rng(5); rng_d = np.random.default_rng(6)
    same = 0; ev = 0; nd = 0; td = tb = 0.0; nw = 0; thr = None
    # umbral fijo en unidades de log-verosimilitud, el mismo en todas las resoluciones
    thr = 2.0
    for k in range(6):
        S = make_network("random", 4, seed=900 + k)
        srcs, taus = draw_sources(rng)
        D, p, tr = generate(S, srcs, taus, n_clutter=35, rng=rng)
        p2, lab = degrade(S, srcs, p, "consistent", rng_d)
        u = max(1e-3, sum(1 for li in lab for l in li if l == 0) / (4 * 40.0))
        pk = [np.sort([t for t, _b in pl]) for pl in p2]
        t0 = time.time(); grid, f, n = B.dense(S, pk, u); td += time.time() - t0
        t0 = time.time(); fb, e = B.bnb(S, pk, u, thr); tb += time.time() - t0
        eq = np.array_equal(np.where(f >= thr)[0], np.where(fb >= thr)[0]) and \
             all(np.allclose(a, b) for a, b in zip(_catalog_impl(grid, f, thr), _catalog_impl(grid, fb, thr)))
        same += int(eq); ev += e; nd += n; nw += 1
    print(f"paso {step} km, {tstep} s: identicos {same}/{nw} | densa {nd // nw} eval., poda {ev // nw} "
          f"({nd / ev:.1f}x) | tiempo densa {td / nw:.1f} s, poda {tb / nw:.1f} s")
