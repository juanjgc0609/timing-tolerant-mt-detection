"""M6-CAL step 5 (declared in card.md before running, 15 November 2026): one confirmation of route B (the pair gain
as a ranking score) on the unread control set B, and its scope at a lower SNR on new seeds. Simulation only.

    python run_step5.py > run_step5_out.txt
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

SEED_B, SEED_P = 20261116, 20261117            # set B (unread until now) and the pairs of step 4, SNR 12
SEED_C8, SEED_P8 = 20261118, 20261119          # new control and pair seeds for SNR 8
N, NPAIR = int(os.environ.get("M6CAL_N", 300)), int(os.environ.get("M6CAL_NPAIR", 40))
PROCS = int(os.environ.get("M6CAL_PROCS", 4))


def job_gain(args):
    snr, a = args
    R.J.SNR = snr
    return R.job_gain(a)


def job_true(args):
    snr, key, dtau, sep = args
    R.J.SNR = snr
    srcs, Y, q = R._scene(key, "pair", True, dtau, sep)
    a = R.m6.adcg_loop(R.E, Y, q, R.KMAX, R.LAG, thr=-np.inf, free="joint", near_km=R.NEAR_KM)
    sp = {id(x) for x in R.MB.spurious(a, srcs)}
    return [x["pair_gain"] for x in a if x["pair_gain"] is not None and id(x) not in sp]


def auc(t, s):
    return float(np.mean(t[:, None] > s[None, :]) + 0.5 * np.mean(t[:, None] == s[None, :]))


def main():
    t0 = time.time()
    print(f"M6-CAL step 5  near_km {R.NEAR_KM}  kmax {R.KMAX}  numpy {np.__version__}")
    sep = R.J.STRATA["near"]
    res = {}
    with mp.Pool(PROCS) as pool:
        for lab, snr, sc, sp in (("SNR 12, set B", 12.0, SEED_B, SEED_P), ("SNR 8, new seeds", 8.0, SEED_C8, SEED_P8)):
            g_sp = np.array([x for xs in pool.map(job_gain, [(snr, ((sc, 2, 2, 1, c), "C3", True)) for c in range(N)],
                                                  chunksize=1) for x in xs], float)
            g_t = np.array([x for xs in pool.map(job_true, [(snr, (sp, int(100 * d), len("near"), c), d, sep)
                                                            for d in R.DTAUS for c in range(NPAIR)], chunksize=1)
                            for x in xs], float)
            a = auc(g_t, g_sp)
            q95 = float(np.quantile(g_sp, .95))
            res[lab] = a
            print(f"{lab}: true gains of near pairs n {g_t.size} (of {len(R.DTAUS) * NPAIR} pairs), median "
                  f"{np.median(g_t):.1f}, q10 {np.quantile(g_t, .1):.1f}; largest spurious gain per control scene "
                  f"n {g_sp.size} of {N}, median {np.median(g_sp):.1f}, q90 {np.quantile(g_sp, .9):.1f}, q95 {q95:.1f}, "
                  f"max {g_sp.max():.1f}; AUC {a:.3f}; true gains above the spurious q95: {np.mean(g_t > q95):.3f}",
                  flush=True)
    a12, a8 = res["SNR 12, set B"], res["SNR 8, new seeds"]
    print(f"\nDECISION (rule of card.md, step 5): route B {'CONFIRMED' if a12 >= 0.90 else 'NOT confirmed'} on set B "
          f"(AUC {a12:.3f}); scope at SNR 8: {'holds' if a8 >= 0.90 else 'does not hold'} (AUC {a8:.3f})")
    print(f"{time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
