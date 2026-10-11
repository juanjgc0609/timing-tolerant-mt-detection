"""task2b_kernel_check.py - sensibilidad (exploratoria) de MT-FW-FREE+SC y MT-FW+SC al nucleo de OpenBLAS.

Las 18 redes de la semilla 7331 que 2A produjo con el nucleo AVX-512 (SkylakeX) se recalculan aqui con el nucleo
por defecto de esta maquina (Prescott): se comparan los registros de MT-FW-FREE+SC y de MT-FW+SC (todas las ventanas)
con los de task2b_parts y se recalcula el recall de cada celda afectada a los topes 0.1 y 0.25.

    python task2b_kernel_check.py [--workers 2]            # sin OPENBLAS_CORETYPE
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import task2b_ablation as B  # noqa: E402

T = B.T
AFFECTED = ([("strong", "ladder_6", j) for j in (2, 3, 4)] + [("thin", "ladder_6", j) for j in range(5)]
            + [(r, "ladder_4", j) for r in ("strong", "thin") for j in range(5)])


def job(a):
    regime, cell, j = a
    seed = 7331
    with open(os.path.join(B.PARTS, f"seed{seed}_{regime}_{cell}_net{j}.json")) as fh:
        P = json.load(fh)
    _ci, nst, nclut, STA, statics, _cat, _d, corr = B.setup_net(seed, regime, cell, j)
    grid = T.default_grid()
    garr = np.asarray(grid)
    Gs = [T.X.radiation_operator(g, STA) for g in garr]
    Qs = [T.X.basis(G) for G in Gs]
    GsF, QsF = B.identity_ops(len(garr), nst)
    rngs = T.window_rngs(seed, cell, j, regime)
    out = {"corr_maxdiff": float(np.max(np.abs(corr - np.array(P["corr"])))), "free": [], "full": []}
    for w in range(P["n_win"]):
        _s, _t, D, _p, tr, _a = T.generate_window(STA, statics, nclut, T.CLUT_AMP[regime], rngs)
        MF, sd = T.matched(D)
        MFn = MF / sd
        out["free"].append(B.as_rec(T.certified_fw_record(STA, MFn, grid, T.DEFAULT_TAUS, GsF, QsF, tr, tcorr=corr)))
        out["full"].append(B.as_rec(T.certified_fw_record(STA, MFn, grid, T.DEFAULT_TAUS, Gs, Qs, tr, tcorr=corr)))
    return a, out, B.blas_arch()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    res = {}
    with ProcessPoolExecutor(a.workers) as ex:
        for key, out, arch in ex.map(job, AFFECTED):
            res[key] = out
    L = [f"nucleo OpenBLAS de esta corrida: {arch} (OPENBLAS_CORETYPE={os.environ.get('OPENBLAS_CORETYPE')})"]
    for regime in ("thin", "strong"):
        for cell in ("ladder_4", "ladder_6"):
            nets = list(range(5))
            W, alt = [], {"MT-FW+SC": [], "MT-FW-FREE+SC": []}
            n_rec_diff = {"MT-FW+SC": 0, "MT-FW-FREE+SC": 0}
            n_hit_diff = {"MT-FW+SC": 0, "MT-FW-FREE+SC": 0}
            maxrel = {"MT-FW+SC": 0.0, "MT-FW-FREE+SC": 0.0}
            cmax = 0.0
            for j in nets:
                with open(os.path.join(B.PARTS, f"seed7331_{regime}_{cell}_net{j}.json")) as fh:
                    P = json.load(fh)
                W.extend(P["windows"])
                r = res.get((regime, cell, j))
                for m, keyr in (("MT-FW+SC", "full"), ("MT-FW-FREE+SC", "free")):
                    for w, wd in enumerate(P["windows"]):
                        new = r[keyr][w] if r else wd["records"][m]
                        old = wd["records"][m]
                        alt[m].append(new)
                        if r:
                            n_rec_diff[m] += new != old
                            n_hit_diff[m] += [h for _s, h in new] != [h for _s, h in old]
                            if len(new) == len(old):
                                for (a_, _h), (b_, _g) in zip(new, old):
                                    maxrel[m] = max(maxrel[m], abs(a_ - b_) / abs(b_))
                            else:
                                maxrel[m] = float("inf")
                if r:
                    cmax = max(cmax, r["corr_maxdiff"])
            k = [w["k"] for w in W]
            orig = {m: T.summarize_records([w["records"][m] for w in W], k) for m in alt}
            newr = {m: T.summarize_records(alt[m], k) for m in alt}
            L.append(f"7331 {regime} {cell}: max|dSC| {cmax:.1e} s; ventanas con registros distintos "
                     + ", ".join(f"{m} {n_rec_diff[m]} (aciertos/orden distintos {n_hit_diff[m]}, max dif. rel. puntaje "
                                 f"{maxrel[m]:.1e})" for m in alt)
                     + "; recall 0.1/0.25 SkylakeX -> defecto: "
                     + ", ".join(f"{m} {orig[m]['0.1']:.3f}/{orig[m]['0.25']:.3f} -> {newr[m]['0.1']:.3f}/{newr[m]['0.25']:.3f}"
                                 for m in alt))
    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(HERE, "task2b_kernel_check.txt"), "w") as fh:
        fh.write(txt + "\n")


if __name__ == "__main__":
    main()
