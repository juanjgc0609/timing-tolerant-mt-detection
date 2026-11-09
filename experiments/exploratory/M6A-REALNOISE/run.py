"""M6A-REALNOISE: gate G-M6a of the v2 specification (M6) with REAL noise, first step (plan T4.23, owner's approval of
8 November 2026). Does the M6 loop (joint refit + resolvability test, ``ttmt.v2.m6.adcg_loop``) beat the greedy v2
loop (5 km exclusion) by >= 0.10 of pair recall on pairs at 0.1-1 s whose K_mec/E is in the upper tercile, when the
noise is real continuous noise and the threshold of the test is calibrated on that noise?

One factor changes with respect to M6-JOINT (simulation, PR #39): the noise. Two arms on the SAME sources:
  REAL  raw 3C noise of eligible 40 s slots of a Task 2M noise day (Alto Tiberina, 16 Hz, written by
        code/t2m_data.py: band 0.5-6.5 Hz, per-channel level rule; no catalogued event within 60 s), through the
        same matched filter;
  GEN   Gaussian raw noise as M6-JOINT (the control that isolates the effect of the real noise).
Geometry: the 8 stations of the day nearest to its centroid (real positions, z = 0); lattice, homogeneous medium,
Ricker P/S pulses, block statics, kernel and synthetic picker (onset bumps of 0.04 s, onset error N(0, 0.02 s)) are
those of M6-JOINT / M6C-SUBSAMPLE. Each source is scaled so that its median strong-block peak is SNR = 12 times the
robust sd of THIS scene's noise through the filter. Sources come from one random stream, noise from another, so the
two arms see identical sources.

Methods (kmax 2): greedy (v2_loop, 5 km exclusion), greedy0 (no exclusion: the control column of the review), m6-none
(adcg_loop, no test), m6 (adcg_loop, thr = z_e^2 + max(0, 95th percentile of the spurious-split gain on a
single-source control of the same arm, separate slots)). LOC+T hit (5 km, 0.5 s), one atom per source.

Gate (declared before the run, card.md): pairs at dtau in {0.1, 0.25, 0.5, 1.0} s, both separation strata pooled,
K_mec/E (ttmt.synth.inject.kmec_over_e on the scaled block vectors) in the upper tercile of all pairs of the run;
REAL arm: pair recall(m6) - pair recall(greedy) >= 0.10 and the 2.5th percentile of its paired bootstrap > 0
(B = 2000, resampling pairs). Secondary: false splits of a fresh single-source control (m6 must stay <= 2 of NCTRL),
the same numbers on GEN, per stratum, per tercile.

    python run.py NOISE_DIR [NREP] [NCTRL] > run_out.txt
NOISE_DIR: a directory written by code/t2m_data.py (noise.npy, meta.json with "eligible"), e.g. day 101231.
"""
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_argv, sys.argv = sys.argv, sys.argv[:1]     # the loaded scripts read sys.argv at import: give them their defaults
J = _load("m6joint_run", HERE.parent / "M6-JOINT" / "run.py")   # loaded by path (several experiments use run.py)
sys.argv = _argv
S = J.S                                                          # M6C-SUBSAMPLE scene helpers

from ttmt.kernel import LEGACY as K  # noqa: E402
from ttmt.radiation import node_operators  # noqa: E402
from ttmt.station_terms import OnsetDelay  # noqa: E402
from ttmt.synth.inject import kmec_over_e  # noqa: E402
from ttmt.traveltime import homogeneous  # noqa: E402
from ttmt.v2.loop import v2_loop  # noqa: E402
from ttmt.v2.m1 import OnsetEngine  # noqa: E402
from ttmt.v2.m6 import adcg_loop  # noqa: E402

SEED = 20261109
NOISE_DIR = sys.argv[1] if len(sys.argv) > 1 else None
NREP = int(sys.argv[2]) if len(sys.argv) > 2 else 30
NCTRL = int(sys.argv[3]) if len(sys.argv) > 3 else 40
NSTA = 8
DTAUS = J.DTAUS
STRATA = J.STRATA
SNR = J.SNR
Z_E2 = J.Z_E2
NAMES = ("greedy", "greedy0", "m6-none", "m6")
B = 2000


def mf(D):
    return np.array([[np.convolve(D[i, c], K.taps[::-1], mode="same") for c in range(3)] for i in range(len(D))])


class Noise:
    """real noise slots of one Task 2M day, station subset of the NSTA stations nearest to the centroid."""

    def __init__(self, d):
        self.arr = np.load(os.path.join(d, "noise.npy"), mmap_mode="r")
        self.meta = json.load(open(os.path.join(d, "meta.json")))
        x, y = np.array(self.meta["x_km"]), np.array(self.meta["y_km"])
        self.sub = np.sort(np.argsort(np.hypot(x - x.mean(), y - y.mean()), kind="stable")[:NSTA])
        self.STA = np.stack([x[self.sub], y[self.sub], np.zeros(NSTA)], 1)
        self.stations = [self.meta["stations"][i] for i in self.sub]
        el = np.array(self.meta["eligible"], int)
        self.slots = el[np.random.default_rng([SEED, 9]).permutation(len(el))]   # fixed order, used once each
        self.k = 0
        assert self.arr.shape[2] >= (el.max() + 1) * S.NT, "slot length does not match NT"

    def take(self):
        s = int(self.slots[self.k])
        self.k += 1
        return s, np.asarray(self.arr[self.sub, :, s * S.NT:(s + 1) * S.NT], float)


def scene(srcs, nmf):
    """normalised record (sources scaled to SNR x robust sd of this noise), onset evidence, scaled fields."""
    zsd = float(np.median(np.abs(nmf)) / 0.6745)
    F = []
    for s in srcs:
        MF, _g, _t = S.field(STA, *s)
        pk = np.linalg.norm(MF, axis=1).max(axis=-1)
        s[2] = s[2] * SNR * zsd / np.median(np.sort(pk)[-6:])
        F.append(S.field(STA, *s))
    Y = sum(f[0] for f in F) + nmf
    sd = float(np.median(np.abs(Y)) / 0.6745)
    return Y / sd, F


def picker(F, r):
    return S.synth_q([f[2] for f in F], NSTA, r)


def noise_for(arm, rn, real):
    if arm == "GEN":
        return None, J.noise_mf(rn, NSTA)
    s, D = real.take()
    return s, mf(D)


def run_arm(arm, real, E, lag):
    out = {"arm": arm}
    gains = []
    for c in range(NCTRL):
        rs, rn = np.random.default_rng([SEED, 1, c]), np.random.default_rng([SEED, 11, c])
        srcs = [J.source(rs)]
        _s, nmf = noise_for(arm, rn, real)
        Y, F = scene(srcs, nmf)
        a = adcg_loop(E, Y, picker(F, rs), 2, lag, thr=-np.inf)
        gains += [x["pair_gain"] for x in a if x["pair_gain"] is not None]
    thr = Z_E2 + max(0.0, float(np.quantile(gains, 0.95)) if gains else 0.0)
    out.update(thr=thr, n_gain=len(gains), gain_med=float(np.median(gains)) if gains else np.nan,
               gain_max=float(max(gains)) if gains else np.nan)
    splits = {"greedy0": 0, "m6": 0}
    for c in range(NCTRL):
        rs, rn = np.random.default_rng([SEED, 2, c]), np.random.default_rng([SEED, 12, c])
        srcs = [J.source(rs)]
        _s, nmf = noise_for(arm, rn, real)
        Y, F = scene(srcs, nmf)
        q = picker(F, rs)
        for k, a in (("greedy0", v2_loop(E, Y, q, 2, lag, min_sep=0.0)), ("m6", adcg_loop(E, Y, q, 2, lag, thr=thr))):
            splits[k] += sum(1 for x in a[1:] if np.linalg.norm(x["x"][:2] - a[0]["x"][:2]) <= 5.0
                             and abs(x["tau"] - a[0]["tau"]) <= 1.5)
    out["splits"] = splits
    rows = []
    for sname, sep in STRATA.items():
        for dtau in DTAUS:
            for rep in range(NREP):
                key = [SEED, 3, int(100 * dtau), len(sname), rep]
                rs, rn = np.random.default_rng(key), np.random.default_rng(key + [1])
                s1 = J.source(rs)
                s2 = J.source(rs, s1[0], sep, 6.0 + dtau)
                srcs = [s1, s2]
                slot, nmf = noise_for(arm, rn, real)
                Y, F = scene(srcs, nmf)
                q = picker(F, rs)
                ke = kmec_over_e(F[0][1], F[1][1])
                res = {"greedy": v2_loop(E, Y, q, 2, lag), "greedy0": v2_loop(E, Y, q, 2, lag, min_sep=0.0),
                       "m6-none": adcg_loop(E, Y, q, 2, lag, thr=None), "m6": adcg_loop(E, Y, q, 2, lag, thr=thr)}
                row = {"stratum": sname, "dtau": dtau, "rep": rep, "slot": slot, "kmec": ke}
                for k, a in res.items():
                    row[k] = len(J.hits(a, srcs))
                rows.append(row)
    out["rows"] = rows
    return out


def report(o, cut):
    rows = o["rows"]
    print(f"== arm {o['arm']}: control {NCTRL} single sources, {o['n_gain']} spurious near atoms, RSS gain median / "
          f"thr / max {o['gain_med']:.1f} / {o['thr']:.1f} / {o['gain_max']:.1f}")
    print(f"   fresh control: false splits of {NCTRL} single events: greedy0 {o['splits']['greedy0']}, "
          f"m6 {o['splits']['m6']}")
    print(f"   {'stratum':7s} {'dtau':>5s} | pair recall " + " ".join(f"{k:>8s}" for k in NAMES)
          + " | upper tercile (n) " + " ".join(f"{k:>8s}" for k in NAMES))
    for sname in STRATA:
        for dtau in DTAUS:
            rr = [r for r in rows if r["stratum"] == sname and r["dtau"] == dtau]
            up = [r for r in rr if r["kmec"] >= cut]
            f = lambda R, k: np.mean([r[k] == 2 for r in R]) if R else np.nan  # noqa: E731
            print(f"   {sname:7s} {dtau:5.2f} |             " + " ".join(f"{f(rr, k):8.2f}" for k in NAMES)
                  + f" | {len(up):17d} " + " ".join(f"{f(up, k):8.2f}" for k in NAMES))
    for lab, R in (("all pairs", rows), ("upper tercile", [r for r in rows if r["kmec"] >= cut]),
                   ("lower two terciles", [r for r in rows if r["kmec"] < cut])):
        d = np.array([(r["m6"] == 2) - (r["greedy"] == 2) for r in R], float)
        d0 = np.array([(r["m6"] == 2) - (r["greedy0"] == 2) for r in R], float)
        rb = np.random.default_rng([SEED, 7, len(lab)])
        bs = d[rb.integers(0, len(d), (B, len(d)))].mean(1)
        print(f"   {lab:19s} n {len(R):4d}  m6 - greedy {d.mean():+.3f} [{np.quantile(bs, .025):+.3f}, "
              f"{np.quantile(bs, .975):+.3f}]  m6 - greedy0 {d0.mean():+.3f}")
    return rows


def main():
    global STA
    t0 = time.time()
    real = Noise(NOISE_DIR) if NOISE_DIR else None
    STA = real.STA if real else None
    if STA is None:
        raise SystemExit("usage: run.py NOISE_DIR [NREP] [NCTRL]")
    XN = S.nodes()
    taus = np.round(np.arange(4.0, 12.0 + 1e-9, 0.1), 10)
    tt, _ = homogeneous(XN, STA, S.VP, S.VS)
    Gs, Qs = node_operators(XN, STA)
    E = OnsetEngine(XN, STA, tt, taus, OnsetDelay(np.full((NSTA, 2), -S.LAG)), Gs, Qs, rounding="nearest")
    lag = np.full((NSTA, 2), S.LAG / S.DT)
    print(f"M6A-REALNOISE  seed {SEED}  {NREP} pairs per cell, {NCTRL} controls  numpy {np.__version__}")
    print(f"noise day {real.meta.get('day')} ({NOISE_DIR}), stations {','.join(real.stations)}, "
          f"{len(real.slots)} eligible slots")
    print(f"SNR {SNR} (median strong-block peak / robust sd of the scene's filtered noise), kmax 2, LOC+T (5 km, 0.5 s)")
    res = {arm: run_arm(arm, real, E, lag) for arm in ("REAL", "GEN")}
    ke = np.array([r["kmec"] for r in res["REAL"]["rows"]])
    cut = float(np.nanquantile(ke, 2 / 3))
    print(f"K_mec/E of the pairs: terciles {np.nanquantile(ke, 1 / 3):.3f} / {cut:.3f}; same sources in both arms\n")
    for arm in ("REAL", "GEN"):
        report(res[arm], cut)
        print()
    used = sorted(r["slot"] for r in res["REAL"]["rows"])
    print(f"REAL slots used by pairs: {len(used)} (+ {2 * NCTRL} control), unique {len(set(used))}")
    gate = [r for r in res["REAL"]["rows"] if r["kmec"] >= cut]
    d = np.array([(r["m6"] == 2) - (r["greedy"] == 2) for r in gate], float)
    lo = np.quantile(d[np.random.default_rng([SEED, 7, len("upper tercile")]).integers(0, len(d), (B, len(d)))]
                     .mean(1), .025)
    ok = d.mean() >= 0.10 and lo > 0 and res["REAL"]["splits"]["m6"] <= 2
    print(f"G-M6a (REAL, upper tercile): m6 - greedy {d.mean():+.3f}, 2.5th pct {lo:+.3f}, m6 false splits "
          f"{res['REAL']['splits']['m6']} -> {'PASS' if ok else 'FAIL'}")
    print(f"\n{time.time() - t0:.0f} s")


STA = None
if __name__ == "__main__":
    sys.exit(main())
