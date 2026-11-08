"""M6-JOINT: does the joint refit with the resolvability test (ttmt.v2.m6.adcg_loop) recover overlapping pairs that
the greedy v2 loop (one atom at a time, 5 km exclusion) loses? Simulation proxy of gate G-M6a (plan T4.19).

Scenes: two double-couple sources, origin separation dtau in {0.1, 0.25, 0.5, 1.0} s, epicentral separation in two
strata: "far" (6-14 km, outside the 5 km exclusion) and "near" (1.5-5 km, inside it). Raw-trace Gaussian noise
through the matched filter, record normalised by its robust sd (as Kernel.matched); source scale such that the
median strong-block peak is about 12 noise sd (SNR 12.4 of Appendix AB). Same network, lattice, kernel and
synthetic picker as M6C-SUBSAMPLE (onset bumps of 0.04 s, onset error N(0, 0.02 s)).

Methods (kmax 2 each):
  greedy   ttmt.v2.loop.v2_loop (onset prior reading, ACAP, waveform-refined deflation, 5 km exclusion);
  m6-none  adcg_loop with no resolvability test (thr None: every candidate kept);
  m6       adcg_loop with the test, thr = z_e^2 (3.29^2, Corollary A.1) + max(0, 95th percentile of the RSS gain of
           a spurious second atom near the first on a noise-only control: one source + noise, separate seeds).
A source is hit by an atom within 5 km and 0.5 s (LOC+T), one atom per source. Pair recall = both sources hit.
Also reported: the control's false-split rate at thr on fresh seeds.

    python experiments/exploratory/M6-JOINT/run.py [NREP] [NCTRL] > experiments/exploratory/M6-JOINT/run_out.txt
"""
import sys
import time
from pathlib import Path

import numpy as np

import importlib.util  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "m6c_scenes", Path(__file__).resolve().parents[1] / "M6C-SUBSAMPLE" / "run.py")
S = importlib.util.module_from_spec(_spec)          # scene helpers of M6C-SUBSAMPLE, loaded by path (both are run.py)
_spec.loader.exec_module(S)

from ttmt.kernel import LEGACY as K  # noqa: E402
from ttmt.radiation import _rand_dc, node_operators  # noqa: E402
from ttmt.station_terms import OnsetDelay  # noqa: E402
from ttmt.traveltime import homogeneous  # noqa: E402
from ttmt.v2.loop import v2_loop  # noqa: E402
from ttmt.v2.m1 import OnsetEngine  # noqa: E402
from ttmt.v2.m6 import adcg_loop  # noqa: E402

SEED = 20261108
NREP = int(sys.argv[1]) if len(sys.argv) > 1 else 20
NCTRL = int(sys.argv[2]) if len(sys.argv) > 2 else 40
DTAUS = (0.1, 0.25, 0.5, 1.0)
STRATA = {"far": (6.0, 14.0), "near": (1.5, 5.0)}
SNR = 12.0
Z_E2 = 3.29 ** 2           # Corollary A.1, total error 0.1: the pair must beat the free one-event fit by z_e^2 ...


def noise_mf(r, ns):
    D = r.normal(size=(ns, 3, S.NT))
    return np.array([[np.convolve(D[i, c], K.taps[::-1], mode="same") for c in range(3)] for i in range(ns)])


def source(r, centre=None, sep=None, tau=6.0):
    if centre is None:
        x = np.r_[r.uniform(-9, 9, 2), 8.0]
    else:
        a = r.uniform(0, 2 * np.pi)
        d = r.uniform(*sep)
        x = np.r_[centre[:2] + d * np.array([np.cos(a), np.sin(a)]), 8.0]
    m = _rand_dc(r).ravel()[[0, 4, 8, 1, 2, 5]]
    return [x, tau + r.uniform(0, S.DT), m, r.uniform(-0.15, 0.15, (2, 8))]


def scene(r, STA, srcs, zsd):
    """record (normalised) of the sources + noise; each source scaled so its median strong-block peak is SNR zsd."""
    F = []
    for s in srcs:
        MF, g, tarr = S.field(STA, *s)
        pk = np.linalg.norm(MF, axis=1).max(axis=-1)
        sc = SNR * zsd / np.median(np.sort(pk)[-6:])
        s[2] = s[2] * sc
        F.append(S.field(STA, *s))
    Y = sum(f[0] for f in F) + noise_mf(r, len(STA))
    sd = float(np.median(np.abs(Y)) / 0.6745)
    q = S.synth_q([f[2] for f in F], len(STA), r)
    return Y / sd, q, F


def hits(atoms, srcs):
    used, got = set(), set()
    for k, s in enumerate(srcs):
        for j, a in enumerate(atoms):
            if j not in used and np.linalg.norm(a["x"][:2] - s[0][:2]) <= 5.0 and abs(a["tau"] - s[1]) <= 0.5:
                used.add(j)
                got.add(k)
                break
    return got


def main():
    t0 = time.time()
    STA, XN = S.network(), S.nodes()
    taus = np.round(np.arange(4.0, 12.0 + 1e-9, 0.1), 10)
    tt, _ = homogeneous(XN, STA, S.VP, S.VS)
    Gs, Qs = node_operators(XN, STA)
    E = OnsetEngine(XN, STA, tt, taus, OnsetDelay(np.full((8, 2), -S.LAG)), Gs, Qs, rounding="nearest")
    lag = np.full((8, 2), S.LAG / S.DT)
    zsd = float(np.median(np.abs(noise_mf(np.random.default_rng(0), 8))) / 0.6745)
    print(f"M6-JOINT  seed {SEED}  {NREP} scenes per cell, {NCTRL} controls  numpy {np.__version__}")
    print(f"SNR {SNR} (median strong-block peak / noise sd), kmax 2, LOC+T hit (5 km, 0.5 s)\n")

    # calibration of thr on a noise-only control (one source), then a fresh control for the false-split rate
    gains = []
    for c in range(NCTRL):
        r = np.random.default_rng([SEED, 1, c])
        srcs = [source(r)]
        Y, q, _F = scene(r, STA, srcs, zsd)
        a = adcg_loop(E, Y, q, 2, lag, thr=-np.inf)          # test computed, every candidate kept
        g = [x["pair_gain"] for x in a if x["pair_gain"] is not None]
        gains += g
    thr = Z_E2 + max(0.0, float(np.quantile(gains, 0.95)) if gains else 0.0)
    print(f"control: {NCTRL} single sources, {len(gains)} spurious near atoms, RSS gain median / thr / max "
          f"{np.median(gains) if gains else float('nan'):.1f} / {thr:.1f} / {max(gains) if gains else float('nan'):.1f}")
    splits = 0
    for c in range(NCTRL):
        r = np.random.default_rng([SEED, 2, c])
        srcs = [source(r)]
        Y, q, _F = scene(r, STA, srcs, zsd)
        a = adcg_loop(E, Y, q, 2, lag, thr=thr)
        splits += sum(1 for x in a[1:] if np.linalg.norm(x["x"][:2] - a[0]["x"][:2]) <= 5.0
                      and abs(x["tau"] - a[0]["tau"]) <= 1.5)
    print(f"fresh control: false splits at thr {thr:.1f}: {splits} of {NCTRL}\n")

    print(f"{'stratum':7s} {'dtau':>5s} | {'pair recall':^26s} | {'source recall':^26s}")
    print(f"{'':7s} {'':>5s} | {'greedy':>8s} {'m6-none':>8s} {'m6':>8s} | {'greedy':>8s} {'m6-none':>8s} {'m6':>8s}")
    for sname, sep in STRATA.items():
        for dtau in DTAUS:
            pr = {k: 0 for k in ("greedy", "m6-none", "m6")}
            sr = dict(pr)
            for rep in range(NREP):
                r = np.random.default_rng([SEED, 3, int(100 * dtau), len(sname), rep])
                s1 = source(r)
                s2 = source(r, s1[0], sep, 6.0 + dtau)
                srcs = [s1, s2]
                Y, q, _F = scene(r, STA, srcs, zsd)
                res = {"greedy": v2_loop(E, Y, q, 2, lag),
                       "m6-none": adcg_loop(E, Y, q, 2, lag, thr=None),
                       "m6": adcg_loop(E, Y, q, 2, lag, thr=thr)}
                for k, a in res.items():
                    h = hits(a, srcs)
                    pr[k] += len(h) == 2
                    sr[k] += len(h)
            print(f"{sname:7s} {dtau:5.2f} | " + " ".join(f"{pr[k] / NREP:8.2f}" for k in pr) + " | "
                  + " ".join(f"{sr[k] / (2 * NREP):8.2f}" for k in sr))
    print(f"\n{time.time() - t0:.0f} s")


if __name__ == "__main__":
    sys.exit(main())
