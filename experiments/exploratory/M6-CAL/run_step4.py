"""M6-CAL step 4 (declared in card.md before running, 14 November 2026): route A (fixed level calibrated with a
margin on set A, confirmed once on set B) against route B (the pair gain as a ranking score). Simulation only.

    python run_step4.py > run_step4_out.txt
"""
import importlib.util
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.argv = ["run.py", "200", "200", "20", "4", "op", "10"]          # near_km 10, the step-3 settings
spec = importlib.util.spec_from_file_location("m6cal", HERE / "run.py")
R = importlib.util.module_from_spec(spec)
sys.modules["m6cal"] = R
spec.loader.exec_module(R)

SEED_A, SEED_B, SEED_P = 20261115, 20261116, 20261117
GRID = (23.8, 34.0, 48.0, 68.0, 96.0)
TARGET, N, NPAIR, THR0 = 0.035, int(os.environ.get("M6CAL_N", 300)), int(os.environ.get("M6CAL_NPAIR", 40)), 23.8


def job_split(a):
    return R.job_split(a)[0] > 0


def job_gain(a):
    return R.job_gain(a)


def job_pair(args):
    key, dtau, sep, thrs = args
    srcs, Y, q = R._scene(key, "pair", True, dtau, sep)
    out = {"greedy": len(R.m6.hits_loctv(R.v2_loop(R.E, Y, q, R.KMAX, R.LAG), srcs)) == 2}
    for th in thrs:
        a = R.m6.adcg_loop(R.E, Y, q, R.KMAX, R.LAG, thr=th, free="joint", near_km=R.NEAR_KM)
        out[th] = len(R.m6.hits_loctv(a, srcs)) == 2
    a = R.m6.adcg_loop(R.E, Y, q, R.KMAX, R.LAG, thr=-np.inf, free="joint", near_km=R.NEAR_KM)
    sp = {id(x) for x in R.MB.spurious(a, srcs)}
    out["true_gains"] = [x["pair_gain"] for x in a if x["pair_gain"] is not None and id(x) not in sp]
    return out


def main():
    t0 = time.time()
    print(f"M6-CAL step 4  near_km {R.NEAR_KM}  kmax {R.KMAX}  sets A {SEED_A}, B {SEED_B}, pairs {SEED_P}  numpy {np.__version__}")
    with mp.Pool(4) as pool:
        thr_star = None
        for th in GRID:
            r = np.mean(pool.map(job_split, [((SEED_A, 2, 2, 1, c), "C3", True, th) for c in range(N)], chunksize=1))
            print(f"set A, with neighbours, thr {th:5.1f}: scenes split {r:.4f} ({int(round(r * N))}/{N})", flush=True)
            if r <= TARGET:
                thr_star = th
                break
        print(f"thr* = {thr_star}\n", flush=True)
        g_sp = [x for xs in pool.map(job_gain, [((SEED_A, 2, 2, 1, c), "C3", True) for c in range(N)], chunksize=1)
                for x in xs]
        okA = thr_star is not None
        if okA:
            for kind, lab in (("C3", "with neighbours"), ("C1", "one source")):
                k = int(sum(pool.map(job_split, [((SEED_B, 2, 2, 1, c), kind, True, thr_star) for c in range(N)],
                                     chunksize=1)))
                lo, hi = R.cp(k, N)
                print(f"set B, {lab:15s} thr* {thr_star}: scenes split {k}/{N}  95 % CI [{lo:.3f}, {hi:.3f}]  "
                      f"{'ok' if hi <= 0.10 else 'FAILS the bound'}", flush=True)
                okA &= hi <= 0.10
        thrs = (THR0,) + ((thr_star,) if thr_star and thr_star != THR0 else ())
        g_true, near = [], {k: [] for k in ("greedy",) + thrs}
        print("\npair recall (LOC+T-V), mismatched pulses, 40 pairs per cell")
        print(f"{'stratum':7s} {'dtau':>5s} {'greedy':>8s} " + " ".join(f"{'m6@' + str(t):>9s}" for t in thrs))
        for sname, sep in R.J.STRATA.items():
            for dtau in R.DTAUS:
                out = pool.map(job_pair, [((SEED_P, int(100 * dtau), len(sname), c), dtau, sep, thrs)
                                          for c in range(NPAIR)], chunksize=1)
                print(f"{sname:7s} {dtau:5.2f} {np.mean([o['greedy'] for o in out]):8.3f} "
                      + " ".join(f"{np.mean([o[t] for o in out]):9.3f}" for t in thrs), flush=True)
                if sname == "near":
                    for k in near:
                        near[k] += [o[k] for o in out]
                    g_true += [x for o in out for x in o["true_gains"]]
    rec = {k: float(np.mean(v)) for k, v in near.items()}
    print("\nnear pairs pooled (n %d): " % len(near["greedy"]) + ", ".join(f"{k} {v:.3f}" for k, v in rec.items()))
    if okA:
        loss = rec[THR0] - rec[thr_star]
        print(f"recall loss of thr* against {THR0}: {loss:+.3f} ({'ok' if loss <= 0.05 else 'FAILS'})")
        okA &= loss <= 0.05
    g_true, g_sp = np.array(g_true, float), np.array(g_sp, float)
    auc = float(np.mean(g_true[:, None] > g_sp[None, :]) + 0.5 * np.mean(g_true[:, None] == g_sp[None, :])) \
        if g_true.size and g_sp.size else float("nan")
    print(f"\nroute B: gains of tested true atoms in near pairs n {g_true.size}, median {np.median(g_true):.1f}, "
          f"q10 {np.quantile(g_true, .1):.1f}; largest spurious gain per control scene n {g_sp.size} of {N}, "
          f"median {np.median(g_sp):.1f}, q90 {np.quantile(g_sp, .9):.1f}; AUC {auc:.3f}")
    print(f"\nDECISION (rule of card.md): route A {'ADOPTED, thr* = ' + str(thr_star) if okA else 'fails'}; "
          f"route B {'supported' if auc >= 0.90 else 'not supported'} (AUC {auc:.3f})")
    print(f"{time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
