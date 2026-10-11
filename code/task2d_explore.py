"""
task2d_explore.py - analisis EXPLORATORIO de la Tarea 2D (no preregistrado, escrito despues de ver analyze_2d_out.txt).

  E1  Donde caen los falsos de MT-FW+SC y QM-3C+SC con el umbral split-conformal congelado en las ventanas de 2A CON
      eventos (la calibracion en ventanas vacias no los ve): para cada falso sobre el umbral, distancia epicentral a la
      fuente verdadera mas cercana de su ventana.  MT-FW+SC se recalcula con el bucle denso (task2d_methods.fw_record_v,
      pico y tope) guardando el nodo de cada evento; control: (puntaje, acierto) identicos a los registros de 2A.
      QM-3C+SC: eventos guardados en task2a_parts con la seleccion de qm_wrapper.qm_records.
  E2  La validacion V6 de task2d_validate.py (error de tiempo inyectado de 0.2 s) con 90 ventanas por celda en lugar de
      30 (las 30 primeras son las mismas).

    OPENBLAS_CORETYPE=... python task2d_explore.py e1 --workers 2      (nucleo por red como la parte (i))
    OPENBLAS_CORETYPE=Prescott python task2d_explore.py e2 --workers 2
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import task2b_ablation as B  # noqa: E402
import task2c_run as C  # noqa: E402
import task2d_methods as D2  # noqa: E402
from task1o_consistent import TOL_KM  # noqa: E402

T = B.T
X = T.X
OUT_DIR = os.path.join(HERE, "task2d_explore_parts")


def fw_record_loc(STA, MFn, grid, taus, Gs, Qs, truth, tcorr, kmax=5, min_sep=5.0):
    """task2d_methods.fw_record_v (pico y tope) que devuelve ademas el nodo de cada evento."""
    ns = len(STA)
    R = MFn.copy()
    garr = np.asarray(grid)
    acc, rec, run = [], [], np.inf
    for _ in range(kmax):
        best = (-np.inf, None, None, None)
        for gi, src in enumerate(garr):
            if any(np.linalg.norm(src[:2] - a[:2]) <= min_sep for a in acc):
                continue
            v, tpk = D2.node_vectors_v(R, STA, src, taus, tcorr, True, True)
            vf = v.reshape(len(taus), -1)
            sc = np.sum((vf @ Qs[gi]) ** 2, axis=1) / 2.0
            j = int(sc.argmax())
            if sc[j] > best[0]:
                best = (float(sc[j]), gi, vf[j], tpk[j])
        s, gi, vf, tpk = best
        if gi is None:
            break
        src = garr[gi]
        run = min(run, s)
        hit = -1
        if truth:
            dd = [float(np.linalg.norm(src[:2] - np.asarray(t_)[:2])) for t_ in truth]
            j = int(np.argmin(dd))
            if dd[j] <= TOL_KM:
                hit = j
        rec.append((run, hit, [float(x) for x in src]))
        acc.append(src)
        G = Gs[gi]
        m, *_ = np.linalg.lstsq(G, vf, rcond=None)
        pred = (G @ m).reshape(2, ns, 3)
        for b in range(2):
            for i in range(ns):
                t0 = int(tpk[b, i])
                lo, hi = max(0, t0 - X.AC_HALF), min(D2.NT, t0 + X.AC_HALF + 1)
                R[i, :, lo:hi] -= np.outer(pred[b, i], X.AC[X.AC_HALF - (t0 - lo):X.AC_HALF + (hi - t0)])
    return rec


def e1_net(args):
    seed, regime, cell, j = args
    path = os.path.join(OUT_DIR, f"e1_seed{seed}_{regime}_{cell}_net{j}.json")
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    req = C.required_kernel(1, seed, regime, cell, j)
    assert C.numpy_blas_arch() == [req], (C.numpy_blas_arch(), req)
    with open(os.path.join(T.PARTS, f"seed{seed}_{regime}_{cell}_net{j}.json")) as fh:
        P2A = json.load(fh)
    _ci, _nst, nclut, STA, statics, _cat, _d, corr = B.setup_net(seed, regime, cell, j)
    assert corr.tolist() == P2A["corr"]
    grid = T.default_grid()
    garr = np.asarray(grid)
    Gs = [X.radiation_operator(g, STA) for g in garr]
    Qs = [X.basis(G) for G in Gs]
    rngs = T.window_rngs(seed, cell, j, regime)
    out = []
    n_same = 0
    for w in range(P2A["n_win"]):
        srcs, taus, D, _p, tr, _a = T.generate_window(STA, statics, nclut, T.CLUT_AMP[regime], rngs)
        st = P2A["windows"][w]
        if st["k"] == 0:
            out.append(None)
            continue
        MF, sd = T.matched(D)
        MFn = MF / sd
        r = fw_record_loc(STA, MFn, grid, T.DEFAULT_TAUS, Gs, Qs, tr, corr)
        same = [[float(a), int(b)] for a, b, _c in r] == st["records"]["MT-FW+SC"]
        n_same += same
        out.append({"rec": r, "same": same})
    res = {"seed": seed, "regime": regime, "cell": cell, "net": j, "windows": out, "n_same": n_same,
           "kernel": C.numpy_blas_arch()}
    with open(path, "w") as fh:
        json.dump(res, fh)
    return res


def e1(workers):
    os.makedirs(OUT_DIR, exist_ok=True)
    jobs = [(s, r, c, j) for s in (2024, 7331) for r in T.REGIMES for c in T.CELLS for j in range(T.N_NETS)
            if C.required_kernel(1, s, r, c, j) == os.environ.get("OPENBLAS_CORETYPE")]
    t0 = time.time()
    with ProcessPoolExecutor(workers) as ex:
        for r in ex.map(e1_net, jobs):
            print(f"   e1 {r['seed']} {r['regime']} {r['cell']} red {r['net']}: identicos {r['n_same']} "
                  f"({time.time() - t0:.0f} s)", flush=True)


def e2(workers):
    import task2d_validate as V
    V.N_V6 = 90
    jobs6 = [(cell, cond) for cell in T.CELLS for cond in ("sin error", "error")]
    with ProcessPoolExecutor(workers) as ex:
        res6 = list(ex.map(V.v6_cell, jobs6))
    L = [f"E2 (exploratorio): V6 con {V.N_V6} ventanas por celda; error de tiempo +-{V.ERR_S} s. Recall a 0.1 / 0.25:"]
    R6 = {(c, cd): (r, n, s) for c, cd, r, n, s in res6}
    for cell in T.CELLS:
        for cond in ("sin error", "error"):
            r, n, s = R6[(cell, cond)]
            L.append(f"   {cell:15s} {cond:9s} ({n} eventos, {s:.0f} s): " + " | ".join(
                f"{m} {r[m][0]:.2f}/{r[m][1]:.2f}" for m in V.V6_METHODS))
        r0, re_ = R6[(cell, "sin error")][0], R6[(cell, "error")][0]
        L.append(f"   {cell:15s} perdida por el error a 0.25: " + ", ".join(
            f"{m} {r0[m][1] - re_[m][1]:+.2f}" for m in V.V6_METHODS) + " ; PK - exacto con error a 0.25: " + ", ".join(
            f"{pk} {re_[pk][1] - re_[ex][1]:+.2f}" for ex, pk in (("MT-FW-EX", "MT-FW"), ("BOMP", "BOMP-PK"),
                                                                  ("MTSTACK", "MTSTACK-PK"), ("GL-FW", "GL-FW-PK"))))
    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(HERE, "task2d_explore_e2.txt"), "w") as fh:
        fh.write(txt + "\n")


if __name__ == "__main__":
    which = sys.argv[1]
    w = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 2
    os.makedirs(T.SCRATCH, exist_ok=True)
    {"e1": e1, "e2": e2}[which](w)
