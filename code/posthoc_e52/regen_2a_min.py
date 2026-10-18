"""
regen_2a_min.py - post hoc (E52): regenerate the per-network files of Task 2A that Tasks 2C and 2D read
(task2a_parts/seed{s}_{regime}_{cell}_net{j}.json), with only the fields they use.

The 2A parts are not stored in the repository. Tasks 2C and 2D read from them the network, the statics, the
catalogue draw, the station corrections and, per window, k, truth, taus and the records of MT-FW+SC and QM-3C+SC.
This script makes the same calls with the same seeds as task2a_confirm.run_job and skips everything else
(picker, degrade_y, GaMMA, PyOcto, MT-FW, MT-FW-NODE+SC, WAVE-PT, QM, QM-3C). Those steps draw from their own
generators (rng_d, rng_n, np.random.seed(42) inside GaMMA), never from the window generators, which is what
task2c_run.run_part1_net already relies on when it regenerates the 2A windows; the 2C and 2D controls (truth,
origins, MT-FW+SC recomputed and compared with the stored record) check it again.

Seed 999 is written with zero windows: Task 2D part (ii) only checks its network and corrections.

    OPENBLAS_CORETYPE=SkylakeX python posthoc_e52/regen_2a_min.py --kernel SkylakeX --workers 4
    OPENBLAS_CORETYPE=Prescott python posthoc_e52/regen_2a_min.py --kernel Prescott --workers 4
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse
import json
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if CODE not in sys.path:
    sys.path.insert(0, CODE)

import task2a_confirm as T  # noqa: E402
import task2c_run as C  # noqa: E402

SEEDS = (2024, 7331)


def kernel_of(seed, regime, cell, j):
    return C.required_kernel(1, seed, regime, cell, j)


def job(args):
    seed, regime, cell, j, n_win, out_path = args
    t_start = time.time()
    req = kernel_of(seed, regime, cell, j) if seed != T.CAL_SEED else "Prescott"
    if C.numpy_blas_arch() != [req] or os.environ.get("OPENBLAS_CORETYPE") != req:
        raise RuntimeError(f"{out_path}: kernel {C.kernel_info()} != {req}")
    ci, _kind, nst, nclut = T.net_setup(seed, cell)
    ri = T.REGIMES.index(regime)
    clut_amp = T.CLUT_AMP[regime]
    STA, statics = T.network(seed, cell, j)
    grid = T.default_grid()
    garr = np.asarray(grid)
    qmc = T.load_qm_cfgs()
    cat, drawn, cat_truth = T.catalogue_real(STA, T.K_CAT, statics, [seed, T.TAG_CAT, ci, j, ri], nclut, clut_amp)
    corr, nacc = T.estimate_station_corrections(STA, cat)
    W = []
    if n_win:
        Gs = [T.X.radiation_operator(g, STA) for g in garr]
        Qs = [T.X.basis(G) for G in Gs]
        lut = T.QM.build_lut(STA, T.VP, T.VS, {**T.QM.QM_CFG, **qmc["QM-3C"]})
        rngs = T.window_rngs(seed, cell, j, regime)
        tmpdir = os.path.join(T.SCRATCH, f"min_s{seed}_{regime}_{cell}_{j}")
        for _w in range(n_win):
            srcs, taus, D, _picks, tr, _amp = T.generate_window(STA, statics, nclut, clut_amp, rngs)
            MF, sd = T.matched(D)
            MFn = MF / sd
            rec = {"MT-FW+SC": T.certified_fw_record(STA, MFn, grid, T.DEFAULT_TAUS, Gs, Qs, tr, tcorr=corr)}
            evs = T.run_qm_variants(STA, D, tmpdir, {"QM-3C+SC": (qmc["QM-3C"], corr)}, lut=lut)
            rec["QM-3C+SC"] = T.QM.qm_records(evs["QM-3C+SC"], tr)
            W.append({"k": len(srcs), "truth": [list(map(float, s)) for s in srcs], "taus": [float(t) for t in taus],
                      "records": {m: [[float(s), int(h)] for s, h in r] for m, r in rec.items()},
                      "qm_events": {"QM-3C+SC": [[float(x) for x in e] for e in evs["QM-3C+SC"]]}})
        shutil.rmtree(tmpdir, ignore_errors=True)
    res = {"seed": seed, "regime": regime, "cell": cell, "net": j, "n_win": n_win, "minimal": True,
           "STA": STA.tolist(), "statics": statics.tolist(), "corr": corr.tolist(), "nacc": nacc.tolist(),
           "cat_drawn": drawn, "cat_src": [c[1].tolist() for c in cat], "cat_tau": [c[2] for c in cat],
           "resid_rms": float(np.sqrt(np.mean((corr - statics) ** 2))), "kernel": C.kernel_info(),
           "windows": W, "seconds": time.time() - t_start}
    tmp = out_path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(res, fh)
    os.replace(tmp, out_path)
    return out_path, res["seconds"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kernel", required=True, choices=("SkylakeX", "Prescott"))
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    os.makedirs(T.PARTS, exist_ok=True)
    os.makedirs(T.SCRATCH, exist_ok=True)
    jobs = []
    for seed in SEEDS + (T.CAL_SEED,):
        for regime in T.REGIMES:
            for cell in T.CELLS:
                for j in range(T.N_NETS):
                    if seed == T.CAL_SEED:
                        if a.kernel != "Prescott":
                            continue
                    elif kernel_of(seed, regime, cell, j) != a.kernel:
                        continue
                    path = os.path.join(T.PARTS, f"seed{seed}_{regime}_{cell}_net{j}.json")
                    if not os.path.exists(path):
                        jobs.append((seed, regime, cell, j, 0 if seed == T.CAL_SEED else T.N_WIN, path))
    jobs.sort(key=lambda x: (-T.SCEN[x[2]][1], x[0]))
    print(f"{len(jobs)} jobs ({a.kernel})", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(a.workers) as ex:
        for p, s in ex.map(job, jobs):
            print(f"   {os.path.basename(p)} {s:.0f} s (total {time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
