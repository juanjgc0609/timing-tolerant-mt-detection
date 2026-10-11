"""
task2c_explore.py - analisis EXPLORATORIO de la Tarea 2C (no preregistrado, escrito despues de ver los resultados de la
parte (i)): por que los metodos en los tiempos predichos exactos (BOMP+SC, MTSTACK+SC, GL-FW+SC) quedan tan por debajo
de MT-FW+SC.

  E1  generador homogeneo IDEAL de 1x (task1x_mtglrt.generate_x, pert=False: el modelo de los detectores es exacto,
      sin retardos ni desajuste; clutter 0.35), semilla 999, 60 ventanas por celda (una red por ventana, como 1x):
      MT-FW, BOMP, MTSTACK y GL-FW sin correcciones.  Y la condicion "pert" de 1x (retardos N(0, 0.15 s),
      polarizacion y sitio) en las mismas fuentes, sin correcciones.
  E3  como E1 ideal, pero con cada fuente movida al nodo de la rejilla mas cercano y su origen redondeado a la
      rejilla de origenes (0.1 s): solo queda el redondeo de los tiempos a la muestra (<= 1 muestra).
  E2  las ventanas de 2A (semilla 2024, regimen thin, redes 0-4, 3 celdas) con correcciones "oraculo" expected_corr de
      task2a_parts (retardo verdadero + desajuste medio del modelo sobre las fuentes del catalogo) en lugar de las
      estimadas: MT-FW y BOMP.

    OPENBLAS_CORETYPE=Prescott python task2c_explore.py [all|e1|e2|e3]
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import task2b_ablation as B  # noqa: E402
import task2c_methods as M  # noqa: E402
from froc_exact import exact_HF  # noqa: E402
from task0_common import DEFAULT_TAUS, default_grid, make_network  # noqa: E402
from task1o_consistent import draw_sources, matched  # noqa: E402

T = B.T
X = T.X


def recall(recs, k):
    _t, H, F = exact_HF([[tuple(x) for x in r] for r in recs])
    H, F = np.asarray(H, float), np.asarray(F, float)
    out = []
    for cap in (0.1, 0.25):
        ok = F.sum(0) / len(recs) <= cap + 1e-12
        out.append(float(H.sum(0)[ok].max() / max(1, sum(k))) if ok.any() else 0.0)
    return out


def e1(n=60):
    grid = default_grid()
    garr = np.asarray(grid)
    L = []
    for cell in T.CELLS:
        kind, nst, nclut = T.SCEN[cell]
        for cond in ("ideal", "pert"):
            rng = np.random.default_rng([999, 2610, T.CELLS.index(cell)])
            rng_m = np.random.default_rng([999, 2611, T.CELLS.index(cell)])
            rng_h = np.random.default_rng([999, 2612, T.CELLS.index(cell)])
            recs = {m: [] for m in ("MT-FW", "BOMP", "MTSTACK", "GL-FW")}
            ks = []
            t0 = time.time()
            for w in range(n):
                S = make_network(kind, nst, seed=[999, 2613, T.CELLS.index(cell), w])
                srcs, taus = draw_sources(rng)
                D, _p, tr, _rad = X.generate_x(S, srcs, taus, nclut, rng, rng_m, cond == "pert", rng_h)
                MF, sd = matched(D)
                MFn = MF / sd
                Gs = [X.radiation_operator(g, S) for g in garr]
                Qs = [X.basis(G) for G in Gs]
                BD = M.BlockDict(S, grid, DEFAULT_TAUS, None, Gs=Gs)
                recs["MT-FW"].append(T.certified_fw_record(S, MFn, grid, DEFAULT_TAUS, Gs, Qs, tr))
                recs["BOMP"].append(M.bomp_record(BD, MFn, tr))
                recs["MTSTACK"].append(M.Stacker(BD).record(MFn, tr))
                recs["GL-FW"].append(M.glfw_record(BD, MFn, D, sd, tr))
                ks.append(len(srcs))
            line = (f"E1 {cond:5s} {cell:15s} {n} ventanas, {sum(ks)} eventos: " + " | ".join(
                f"{m} " + "/".join(f"{x:.3f}" for x in recall(r, ks)) for m, r in recs.items())
                + f"  ({time.time() - t0:.0f} s)")
            print(line, flush=True)
            L.append(line)
    return L


def e3(n=60):
    grid = default_grid()
    garr = np.asarray(grid)
    L = []
    for cell in T.CELLS:
        kind, nst, nclut = T.SCEN[cell]
        rng = np.random.default_rng([999, 2610, T.CELLS.index(cell)])
        rng_m = np.random.default_rng([999, 2611, T.CELLS.index(cell)])
        rng_h = np.random.default_rng([999, 2612, T.CELLS.index(cell)])
        recs = {m: [] for m in ("MT-FW", "BOMP", "MTSTACK", "GL-FW")}
        ks = []
        t0 = time.time()
        for w in range(n):
            S = make_network(kind, nst, seed=[999, 2613, T.CELLS.index(cell), w])
            srcs, taus = draw_sources(rng)
            srcs = [garr[int(np.argmin(np.linalg.norm(garr - s_, axis=1)))].copy() for s_ in srcs]
            taus = [float(DEFAULT_TAUS[int(np.argmin(np.abs(DEFAULT_TAUS - t_)))]) for t_ in taus]
            D, _p, tr, _rad = X.generate_x(S, srcs, taus, nclut, rng, rng_m, False, rng_h)
            MF, sd = matched(D)
            MFn = MF / sd
            Gs = [X.radiation_operator(g, S) for g in garr]
            Qs = [X.basis(G) for G in Gs]
            BD = M.BlockDict(S, grid, DEFAULT_TAUS, None, Gs=Gs)
            recs["MT-FW"].append(T.certified_fw_record(S, MFn, grid, DEFAULT_TAUS, Gs, Qs, tr))
            recs["BOMP"].append(M.bomp_record(BD, MFn, tr))
            recs["MTSTACK"].append(M.Stacker(BD).record(MFn, tr))
            recs["GL-FW"].append(M.glfw_record(BD, MFn, D, sd, tr))
            ks.append(len(srcs))
        line = (f"E3 ideal, fuentes en nodos {cell:15s} {n} ventanas, {sum(ks)} eventos: " + " | ".join(
            f"{m} " + "/".join(f"{x:.3f}" for x in recall(r, ks)) for m, r in recs.items())
            + f"  ({time.time() - t0:.0f} s)")
        print(line, flush=True)
        L.append(line)
    return L


def e2(seed=2024, regime="thin"):
    L = []
    for cell in T.CELLS:
        recs = {m: [] for m in ("MT-FW+SC", "BOMP+SC", "MT-FW+oracleSC", "BOMP+oracleSC")}
        ks = []
        t0 = time.time()
        for j in range(T.N_NETS):
            with open(os.path.join(T.PARTS, f"seed{seed}_{regime}_{cell}_net{j}.json")) as fh:
                P = json.load(fh)
            _ci, nst, nclut, STA, statics, _cat, _d, corr = B.setup_net(seed, regime, cell, j)
            assert corr.tolist() == P["corr"]
            eo = np.asarray(P["expected_corr"])
            grid = T.default_grid()
            garr = np.asarray(grid)
            Gs = [X.radiation_operator(g, STA) for g in garr]
            Qs = [X.basis(G) for G in Gs]
            BDo = M.BlockDict(STA, grid, T.DEFAULT_TAUS, eo, Gs=Gs)
            with open(os.path.join(HERE, "task2c_parts", f"p1_seed{seed}_{regime}_{cell}_net{j}.json")) as fh:
                C = json.load(fh)
            rngs = T.window_rngs(seed, cell, j, regime)
            for w in range(P["n_win"]):
                srcs, taus, D, _p, tr, _a = T.generate_window(STA, statics, nclut, T.CLUT_AMP[regime], rngs)
                assert [list(map(float, s)) for s in srcs] == P["windows"][w]["truth"]
                MF, sd = matched(D)
                MFn = MF / sd
                recs["MT-FW+SC"].append(P["windows"][w]["records"]["MT-FW+SC"])
                recs["BOMP+SC"].append(C["windows"][w]["records"]["BOMP+SC"])
                recs["MT-FW+oracleSC"].append(T.certified_fw_record(STA, MFn, grid, T.DEFAULT_TAUS, Gs, Qs, tr, tcorr=eo))
                recs["BOMP+oracleSC"].append(M.bomp_record(BDo, MFn, tr))
                ks.append(len(srcs))
        line = (f"E2 {seed} {regime} {cell:15s} {len(ks)} ventanas: " + " | ".join(
            f"{m} " + "/".join(f"{x:.3f}" for x in recall(r, ks)) for m, r in recs.items())
            + f"  ({time.time() - t0:.0f} s)")
        print(line, flush=True)
        L.append(line)
    return L


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    L = [f"task2c_explore.py (EXPLORATORIO, post hoc); nucleo {B.blas_arch()}"]
    if which in ("e1", "all"):
        L += e1()
    if which in ("e2", "all"):
        L += e2()
    if which in ("e3",):
        L += e3()
    with open(os.path.join(HERE, f"task2c_explore_{which}.txt"), "w") as fh:
        fh.write("\n".join(L) + "\n")
