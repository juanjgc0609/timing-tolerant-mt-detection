"""M6-BIAS: biases of the one-event alternative of the resolvability test, and the threshold recalibrated with
neighbouring events (review T5.2 of 9 November 2026, teoria-fase5/revision-A1-A3a-M6.md, items b1-b3 and section 3).

Simulation as M6-JOINT (8-station ring, homogeneous medium, Ricker P/S, synthetic picker, Gaussian noise, SNR 12).
Null scenes (no true pair, so every accepted near atom is a false split):
  C1   one source (the control of M6-JOINT);
  C3   three sources: the target, a neighbour 6-14 km away and 0.5-1.5 s later (its arrivals overlap the target's),
       and an isolated one 6-14 km away and 3.0 s later.
For each null scene adcg_loop runs with every candidate kept (thr = -inf, kmax = number of sources + 1) and records
the RSS gain of each spurious near atom, with the legacy one-event fit (bounded +-8 samples, P then S, others fixed)
and with the joint one (``free_one_event_joint``: +-16 samples, P and S together, others refitted). The bias is the
difference legacy - joint on the same candidates.
Threshold: thr = 95th percentile of the joint gain on C3 (no floor: z_e^2 is not a test level, review item a).
Then, on fresh seeds: false splits of C1 and C3 at thr, and pair recall (LOC+T) of the M6-JOINT near/far pairs at
dtau 0.25 and 0.5 s with thr (kmax 2), against the greedy loop with the 5 km exclusion.

    python run.py [NCAL] [NTEST] [NPAIR] > run_out.txt
"""
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
_argv, sys.argv = sys.argv, sys.argv[:1]
_spec = importlib.util.spec_from_file_location("m6joint_run", HERE.parent / "M6-JOINT" / "run.py")
J = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(J)
sys.argv = _argv
S = J.S

from ttmt.radiation import node_operators  # noqa: E402
from ttmt.station_terms import OnsetDelay  # noqa: E402
from ttmt.traveltime import homogeneous  # noqa: E402
from ttmt.v2 import m6  # noqa: E402
from ttmt.v2.loop import v2_loop  # noqa: E402
from ttmt.v2.m1 import OnsetEngine  # noqa: E402

SEED = 20261109
NCAL = int(sys.argv[1]) if len(sys.argv) > 1 else 30
NTEST = int(sys.argv[2]) if len(sys.argv) > 2 else 30
NPAIR = int(sys.argv[3]) if len(sys.argv) > 3 else 20


def null_scene(r, kind):
    s1 = J.source(r)
    if kind == "C1":
        return [s1]
    s2 = J.source(r, s1[0], (6.0, 14.0), 6.0 + r.uniform(0.5, 1.5))
    s3 = J.source(r, s1[0], (6.0, 14.0), 9.0)
    return [s1, s2, s3]


def spurious(atoms, srcs):
    """atoms within 5 km and 1.5 s of a true source that are not the one atom matched to it."""
    used = set()
    for s in srcs:
        for j, a in enumerate(atoms):
            if j not in used and np.linalg.norm(a["x"][:2] - s[0][:2]) <= 5.0 and abs(a["tau"] - s[1]) <= 0.5:
                used.add(j)
                break
    return [a for j, a in enumerate(atoms) if j not in used and any(
        np.linalg.norm(a["x"][:2] - s[0][:2]) <= 5.0 and abs(a["tau"] - s[1]) <= 1.5 for s in srcs)]


def main():
    t0 = time.time()
    STA, XN = S.network(), S.nodes()
    taus = np.round(np.arange(4.0, 12.0 + 1e-9, 0.1), 10)
    tt, _ = homogeneous(XN, STA, S.VP, S.VS)
    Gs, Qs = node_operators(XN, STA)
    E = OnsetEngine(XN, STA, tt, taus, OnsetDelay(np.full((8, 2), -S.LAG)), Gs, Qs, rounding="nearest")
    lag = np.full((8, 2), S.LAG / S.DT)
    zsd = float(np.median(np.abs(J.noise_mf(np.random.default_rng(0), 8))) / 0.6745)
    print(f"M6-BIAS  seed {SEED}  NCAL {NCAL}  NTEST {NTEST}  NPAIR {NPAIR}  numpy {np.__version__}\n")
    gains = {}
    for kind in ("C1", "C3"):
        for free in ("legacy", "joint"):
            g = []
            for c in range(NCAL):
                r = np.random.default_rng([SEED, 1, len(kind), c])
                srcs = null_scene(r, kind)
                Y, q, _F = J.scene(r, STA, srcs, zsd)
                a = m6.adcg_loop(E, Y, q, len(srcs) + 1, lag, thr=-np.inf, free=free)
                g += [x["pair_gain"] for x in spurious(a, srcs) if x["pair_gain"] is not None]
            gains[kind, free] = np.array(g)
            print(f"{kind} {free:6s}: {len(g):3d} spurious near atoms, gain median {np.median(g) if g else np.nan:9.1f} "
                  f"p95 {np.quantile(g, .95) if g else np.nan:9.1f} max {max(g) if g else np.nan:9.1f}")
    print("\n(the legacy and joint runs share the candidate stream only until the first rejection; the bias is read on "
          "the distributions)")
    for kind in ("C1", "C3"):
        L, Jg = gains[kind, "legacy"], gains[kind, "joint"]
        if len(L) and len(Jg):
            print(f"bias {kind}: median legacy - joint {np.median(L) - np.median(Jg):+.1f}, "
                  f"p95 {np.quantile(L, .95) - np.quantile(Jg, .95):+.1f}")
    thr = float(np.quantile(gains["C3", "joint"], .95)) if len(gains["C3", "joint"]) else 0.0
    thr_old = J.Z_E2 + max(0.0, float(np.quantile(gains["C1", "legacy"], .95)) if len(gains["C1", "legacy"]) else 0.0)
    print(f"\nthr (joint, C3 p95) {thr:.1f}; the first version's thr (z_e^2 + C1 legacy p95) {thr_old:.1f}\n")
    for kind in ("C1", "C3"):
        n = {"legacy@old": 0, "joint@new": 0}
        for c in range(NTEST):
            r = np.random.default_rng([SEED, 2, len(kind), c])
            srcs = null_scene(r, kind)
            Y, q, _F = J.scene(r, STA, srcs, zsd)
            for key, free, th in (("legacy@old", "legacy", thr_old), ("joint@new", "joint", thr)):
                n[key] += len(spurious(m6.adcg_loop(E, Y, q, len(srcs) + 1, lag, thr=th, free=free), srcs))
        print(f"fresh {kind}: false splits in {NTEST} scenes: legacy at old thr {n['legacy@old']}, "
              f"joint at new thr {n['joint@new']}")
    print(f"\n{'stratum':7s} {'dtau':>5s}  pair recall: greedy  m6 legacy@old  m6 joint@new")
    for sname, sep in J.STRATA.items():
        for dtau in (0.25, 0.5):
            pr = np.zeros(3)
            for rep in range(NPAIR):
                r = np.random.default_rng([SEED, 3, int(100 * dtau), len(sname), rep])
                s1 = J.source(r)
                s2 = J.source(r, s1[0], sep, 6.0 + dtau)
                srcs = [s1, s2]
                Y, q, _F = J.scene(r, STA, srcs, zsd)
                for k, a in enumerate((v2_loop(E, Y, q, 2, lag),
                                       m6.adcg_loop(E, Y, q, 2, lag, thr=thr_old, free="legacy"),
                                       m6.adcg_loop(E, Y, q, 2, lag, thr=thr, free="joint"))):
                    pr[k] += len(J.hits(a, srcs)) == 2
            print(f"{sname:7s} {dtau:5.2f}  " + "  ".join(f"{x / NPAIR:6.2f}" for x in pr))
    print(f"\n{time.time() - t0:.0f} s")


if __name__ == "__main__":
    sys.exit(main())
