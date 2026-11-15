"""M6-CAL: calibration of the resolvability threshold and the M6 pair recall under the project's hit rule, after the
theory review T5.2 (9 November 2026) and the adversarial pass 2 (CR-29, CR-30, CR-37). Simulation (cloud), the step
that comes before any reading of G-M6a.

Changes against M6-JOINT, each from a finding:
  scoring    LOC+T-V (``ttmt.v2.m6.hits_loctv`` = ``ttmt.stats.label_voronoi``, acceptance order), with the
             duplicator negative control (ADR-0008 item 3) on the greedy loop (CR-29);
  kmax       4 everywhere (CR-37);
  one-event  the joint alternative (``free_one_event_joint``: P and S together, +-1 s, neighbours refitted; T5.2 b1-b3);
  pulses     mismatched by default: per arrival, centre frequency x U[0.7, 1.4] and phase rotation U[-45, 45] deg
             (the ranges of code/gen_realistic.py; the kernel does not change), the case where the criticism found a
             false split (CR-30); exact pulses reported as a reference;
  threshold  thr = the larger of the 95th percentiles of the spurious-split gain on NCAL >= 200 one-source scenes and
             on NCAL >= 200 scenes WITH NEIGHBOURS (target + overlapping neighbour + isolated source, M6-BIAS),
             mismatched pulses; no floor (T5.2 a). Both, because the gain's scale differs between them (smoke run).
             The level is per WINDOW (T5.2 d): the statistic is the largest spurious gain of each control scene, since
             with kmax 4 a scene makes up to three tests (smoke run: per-test calibration split 7 of 12 atoms).
Read-outs: false splits per scene on fresh controls (one source; with neighbours; exact and mismatched) with a 95 %
Clopper-Pearson interval; pair recall (LOC+T-V) of near / far pairs at dtau 0.1, 0.25, 0.5, 1 s for greedy (5 km
exclusion), greedy0 (no exclusion), m6, and the duplicator of greedy.

    python run.py [NCAL] [NTEST] [NPAIR] [PROCS] [op] > run_out.txt        (op: step 2 of card.md, the threshold
    calibrated on the operating loop by bisection: the smallest thr with <= 5 % of the calibration scenes split in
    both controls)
"""
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import importlib.util  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.signal import hilbert  # noqa: E402
from scipy.stats import beta  # noqa: E402

HERE = Path(__file__).resolve().parent


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_argv, sys.argv = sys.argv, sys.argv[:1]
J = _load("m6joint_run", HERE.parent / "M6-JOINT" / "run.py")
MB = _load("m6bias_run", HERE.parent / "M6-BIAS" / "run.py")
sys.argv = _argv
S = J.S

from ttmt.radiation import node_operators  # noqa: E402
from ttmt.station_terms import OnsetDelay  # noqa: E402
from ttmt.traveltime import homogeneous  # noqa: E402
from ttmt.v2 import m6  # noqa: E402
from ttmt.v2.loop import v2_loop  # noqa: E402
from ttmt.v2.m1 import OnsetEngine  # noqa: E402

SEED = 20261110
NCAL = int(sys.argv[1]) if len(sys.argv) > 1 else 200
NTEST = int(sys.argv[2]) if len(sys.argv) > 2 else 200
NPAIR = int(sys.argv[3]) if len(sys.argv) > 3 else 20
PROCS = int(sys.argv[4]) if len(sys.argv) > 4 else 4
OP = len(sys.argv) > 5 and sys.argv[5] == "op"
NEAR_KM = float(sys.argv[6]) if len(sys.argv) > 6 else 5.0      # radius of the resolvability test (step 3: 10)
KMAX = 4
DTAUS = (0.1, 0.25, 0.5, 1.0)
_RIC0 = S.ricker
MODE = {"mismatch": True, "rng": None}


def _ricker(t, fc):
    """the data pulse: the kernel's Ricker, or (mismatch) a stretched and phase-rotated one, drawn per call (= per
    arrival and phase in M6C-SUBSAMPLE.field)."""
    if not MODE["mismatch"]:
        return _RIC0(t, fc)
    r = MODE["rng"]
    w = _RIC0(t, fc * r.uniform(0.7, 1.4))
    ph = np.deg2rad(r.uniform(-45.0, 45.0))
    return np.cos(ph) * w + np.sin(ph) * np.imag(hilbert(w))


S.ricker = _ricker
STA, XN = S.network(), S.nodes()
TAUS = np.round(np.arange(4.0, 12.0 + 1e-9, 0.1), 10)
_tt, _ = homogeneous(XN, STA, S.VP, S.VS)
_Gs, _Qs = node_operators(XN, STA)
E = OnsetEngine(XN, STA, _tt, TAUS, OnsetDelay(np.full((8, 2), -S.LAG)), _Gs, _Qs, rounding="nearest")
LAG = np.full((8, 2), S.LAG / S.DT)
ZSD = float(np.median(np.abs(J.noise_mf(np.random.default_rng(0), 8))) / 0.6745)


def _scene(key, kind, mismatch, dtau=None, sep=None):
    r = np.random.default_rng(key)
    MODE["mismatch"], MODE["rng"] = mismatch, np.random.default_rng(list(key) + [99])
    if kind == "pair":
        s1 = J.source(r)
        srcs = [s1, J.source(r, s1[0], sep, 6.0 + dtau)]
    else:
        srcs = MB.null_scene(r, kind)
    Y, q, _F = J.scene(r, STA, srcs, ZSD)
    return srcs, Y, q


def job_gain(args):
    key, kind, mismatch = args
    srcs, Y, q = _scene(key, kind, mismatch)
    a = m6.adcg_loop(E, Y, q, KMAX, LAG, thr=-np.inf, free="joint", near_km=NEAR_KM)
    g = [x["pair_gain"] for x in MB.spurious(a, srcs) if x["pair_gain"] is not None]
    return [max(g)] if g else []              # level per window (T5.2 d): the largest spurious gain of the scene


def job_split(args):
    key, kind, mismatch, thr = args
    srcs, Y, q = _scene(key, kind, mismatch)
    return (len(MB.spurious(m6.adcg_loop(E, Y, q, KMAX, LAG, thr=thr, free="joint", near_km=NEAR_KM), srcs)),
            len(MB.spurious(v2_loop(E, Y, q, KMAX, LAG, min_sep=0.0), srcs)))


def job_pair(args):
    key, dtau, sep, thr = args
    srcs, Y, q = _scene(key, "pair", True, dtau, sep)
    g = v2_loop(E, Y, q, KMAX, LAG)
    res = {"greedy": g, "greedy0": v2_loop(E, Y, q, KMAX, LAG, min_sep=0.0),
           "m6": m6.adcg_loop(E, Y, q, KMAX, LAG, thr=thr, free="joint", near_km=NEAR_KM),
           "dup(greedy)": m6.duplicate_atoms(g, dtau / 2.0)}
    return {k: len(m6.hits_loctv(a, srcs)) == 2 for k, a in res.items()}


def cp(k, n):
    lo = beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return lo, hi


def main():
    t0 = time.time()
    print(f"near_km {NEAR_KM}")
    print(f"M6-CAL  seed {SEED}  NCAL {NCAL}  NTEST {NTEST}  NPAIR {NPAIR}  kmax {KMAX}  numpy {np.__version__}")
    print("8-station ring, homogeneous medium, synthetic picker, Gaussian noise, SNR 12; LOC+T-V scoring\n")
    with mp.get_context("fork").Pool(PROCS) as pool:
        if OP:
            def rate(th):
                r = []
                for kind in ("C1", "C3"):
                    o = pool.map(job_split, [((SEED, 1, len(kind), c), kind, True, th) for c in range(NCAL)],
                                 chunksize=1)
                    r.append(np.mean([x[0] > 0 for x in o]))
                print(f"  thr {th:8.1f}: scenes split C1 {r[0]:.3f}  C3 {r[1]:.3f}", flush=True)
                return max(r)
            lo, hi = 5.3, 400.0
            if rate(hi) > 0.05:
                raise SystemExit("even thr 400 splits more than 5 %")
            for _ in range(6):
                mid = 0.5 * (lo + hi)
                lo, hi = (lo, mid) if rate(mid) <= 0.05 else (mid, hi)
            thr = hi
            print(f"thr (operating-loop calibration, step 2) = {thr:.1f}\n")
        else:
            thr, parts = -np.inf, []
            for kind in ("C1", "C3"):
                g = np.array([x for xs in pool.map(job_gain, [((SEED, 1, len(kind), c), kind, True)
                                                              for c in range(NCAL)], chunksize=1) for x in xs])
                p95 = float(np.quantile(g, 0.95)) if g.size else -np.inf
                thr = max(thr, p95)
                print(f"calibration {kind}: {NCAL} scenes, mismatched pulses: {g.size} scenes with a spurious near "
                      f"atom, max gain median {np.median(g) if g.size else np.nan:.1f} p95 {p95:.1f} max "
                      f"{g.max() if g.size else np.nan:.1f}")
            print(f"thr = max of the two p95 = {thr:.1f}\n")
        print(f"{'fresh control':34s} {'scenes':>6s} {'m6 splits':>9s} {'scenes with >=1':>15s} {'95 % CI':>15s} "
              f"{'greedy0 splits':>14s}")
        for kind, mis, n in (("C1", True, NTEST), ("C3", True, NTEST), ("C1", False, NTEST // 2),
                             ("C3", False, NTEST // 2)):
            out = pool.map(job_split, [((SEED, 2, len(kind), int(mis), c), kind, mis, thr) for c in range(n)],
                           chunksize=1)
            sp = np.array([o[0] for o in out])
            k = int((sp > 0).sum())
            lo, hi = cp(k, n)
            lab = f"{'one source' if kind == 'C1' else 'with neighbours'}, {'mismatched' if mis else 'exact'}"
            print(f"{lab:34s} {n:6d} {sp.sum():9d} {k:15d} {f'[{lo:.3f}, {hi:.3f}]':>15s} "
                  f"{sum(o[1] for o in out):14d}")
        names = ("greedy", "greedy0", "m6", "dup(greedy)")
        print(f"\npair recall (LOC+T-V), mismatched pulses, {NPAIR} pairs per cell")
        print(f"{'stratum':7s} {'dtau':>5s} " + " ".join(f"{k:>11s}" for k in names))
        for sname, sep in J.STRATA.items():
            for dtau in DTAUS:
                out = pool.map(job_pair, [((SEED, 3, int(100 * dtau), len(sname), rep), dtau, sep, thr)
                                          for rep in range(NPAIR)], chunksize=1)
                print(f"{sname:7s} {dtau:5.2f} " + " ".join(f"{np.mean([o[k] for o in out]):11.2f}" for k in names))
    print(f"\n{time.time() - t0:.0f} s")


if __name__ == "__main__":
    sys.exit(main())
