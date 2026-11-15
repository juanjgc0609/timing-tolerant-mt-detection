"""M6A-REALNOISE: gate G-M6a of the v2 specification (M6) with REAL noise, first step (plan T4.23, owner's approval of
8 November 2026). Does the M6 loop (joint refit + resolvability test, ``ttmt.v2.m6.adcg_loop``) beat the greedy v2
loop (5 km exclusion) by >= 0.10 of pair recall at the registered false-detection budget on pairs at 0.1-1 s whose
K_mec/E is in the upper tercile, when the noise is real continuous noise?

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

Methods, settings of M6-CAL: kmax 4; greedy (v2_loop, 5 km exclusion), greedy0 (no exclusion, the review's control
column), m6-none (adcg_loop, no test), m6 (adcg_loop, one-event alternative "joint", near_km 10, every candidate
kept, thr = -inf), m6-gain (secondary, same atoms as m6, other score). Hit rule: LOC+T-V (ttmt.stats.label_voronoi,
ADR-0008 item 2).

Main criterion (card.md, "Declared change of criteria", approved by the owner on 15 November 2026; route B of
M6-CAL steps 4-5): no fixed threshold. Each method's detections carry a score (rule in `scored`), all scenes of an
arm (NFRESH fresh control scenes with neighbours + all pair scenes) are pooled, and the threshold is the exact-FROC
operating point at BUDGET = 30 false detections per hour (ttmt.stats.operating_point). Hours scanned: SCAN_S = 8.1 s
of origin times per scene (81 bins of 0.1 s, 4-12 s). Gate: REAL arm, upper K_mec/E tercile of the pairs, both
strata and the four separations pooled: pair recall(m6) - pair recall(greedy) >= 0.10 and the 2.5th percentile of
its paired bootstrap > 0 (B = 2000, resampling pairs, operating points held fixed). False splits of m6 at the
reference threshold 23.8 (scenes, of 200; reference 12) are descriptive, not in the gate.
The gate is NOT read until the per-station S split and the calibration on real isolated events (CR-30) are closed
(card.md); until then a run is a diagnostic.

    python run.py NOISE_DIR [NREP] [NFRESH] > run_out.txt
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
MB = _load("m6bias_run", HERE.parent / "M6-BIAS" / "run.py")      # null scenes' false-split rule
sys.argv = _argv
S = J.S                                                          # M6C-SUBSAMPLE scene helpers

from ttmt.kernel import LEGACY as K  # noqa: E402
from ttmt.radiation import node_operators  # noqa: E402
from ttmt.station_terms import OnsetDelay  # noqa: E402
from ttmt.synth.inject import kmec_over_e  # noqa: E402
from ttmt.traveltime import homogeneous  # noqa: E402
from ttmt.v2.loop import v2_loop  # noqa: E402
from ttmt.v2.m1 import OnsetEngine  # noqa: E402
from ttmt.stats import FALSE, false_count, hits_at, label_voronoi, operating_point  # noqa: E402
from ttmt.v2.m6 import adcg_loop  # noqa: E402

SEED = 20261109
NOISE_DIR = sys.argv[1] if len(sys.argv) > 1 else None
NREP = int(sys.argv[2]) if len(sys.argv) > 2 else 30
NFRESH = int(sys.argv[3]) if len(sys.argv) > 3 else 200   # fresh control scenes with neighbours
SPLIT_MAX = 12           # reference: scenes with a false split among 200 at THR_REF (descriptive, criterion 2)
THR_REF = 23.8           # reference threshold of the descriptive false-split count
BUDGET = 30.0            # registered false-detection budget per hour (ADR-0008 item 4)
SCAN_S = 8.1             # seconds of origin time scanned per scene (taus 4-12 s, step 0.1 s)
NSTA = 8
KMAX = 4                 # M6-CAL (CR-37)
NEAR_KM = 10.0           # radius of the resolvability test, twice the LOC+T radius (M6-CAL step 3)
DTAUS = J.DTAUS
STRATA = J.STRATA
SNR = J.SNR
NAMES = ("greedy", "greedy0", "m6-none", "m6", "m6-gain")
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
        if self.k >= len(self.slots):
            raise SystemExit(f"not enough eligible noise slots ({len(self.slots)}); lower NFRESH and declare it")
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


def null_scene(rs):
    """fresh false-detection control with neighbours (review T5.2): the target, a neighbour 6-14 km away and
    0.5-1.5 s later (overlapping arrivals) and an isolated source 6-14 km away and 3 s later."""
    s1 = J.source(rs)
    return [s1, J.source(rs, s1[0], (6.0, 14.0), 6.0 + rs.uniform(0.5, 1.5)), J.source(rs, s1[0], (6.0, 14.0), 9.0)]


def scored(name, atoms):
    """(scores, xy, tau) of the detections of one method in one scene. Score rule (card.md, criterion 2, declared
    before any run): greedy, greedy0, m6-none: the ranking statistic s of each atom; m6: s for an atom that was not
    tested, min(s, pair gain) for a tested one (it must pass the threshold on both); m6-gain (secondary, the literal
    reading of route B): the pair gain alone for a tested atom."""
    sc = []
    for a in atoms:
        g = a.get("pair_gain")
        if name == "m6" and g is not None:
            sc.append(min(float(a["s"]), float(g)))
        elif name == "m6-gain" and g is not None:
            sc.append(float(g))
        else:
            sc.append(float(a["s"]))
    xy = np.array([np.asarray(a["x"], float)[:2] for a in atoms]).reshape(-1, 2)
    return np.array(sc, float), xy, np.array([float(a["tau"]) for a in atoms], float)


def run_methods(E, Y, q, lag):
    a6 = adcg_loop(E, Y, q, KMAX, lag, thr=-np.inf, near_km=NEAR_KM)     # every candidate kept: route B
    return {"greedy": v2_loop(E, Y, q, KMAX, lag), "greedy0": v2_loop(E, Y, q, KMAX, lag, min_sep=0.0),
            "m6-none": adcg_loop(E, Y, q, KMAX, lag, thr=None, near_km=NEAR_KM), "m6": a6, "m6-gain": a6}


class Pool:
    """detections of all scenes of one arm, labelled per scene under LOC+T-V (prefix-stable, so the labels of the
    full list hold at every threshold); targets numbered globally."""

    def __init__(self):
        self.s = {k: [] for k in NAMES}
        self.lab = {k: [] for k in NAMES}
        self.nt = 0

    def add(self, res, srcs):
        exy = np.array([np.asarray(s[0], float)[:2] for s in srcs])
        et = np.array([float(s[1]) for s in srcs])
        ids = np.arange(self.nt, self.nt + len(srcs))
        for k in NAMES:
            sc, xy, tau = scored(k, res[k])
            lab = label_voronoi(sc, xy, tau, exy, et) if len(sc) else np.zeros(0, int)
            self.s[k].append(sc)
            self.lab[k].append(np.where(lab >= 0, ids[np.maximum(lab, 0)], FALSE))
        self.nt += len(srcs)
        return ids

    def hit(self, k, hours):
        s, lab = np.concatenate(self.s[k]), np.concatenate(self.lab[k])
        thr = operating_point(s, lab, hours, BUDGET)
        return thr, false_count(s, lab, thr), hits_at(s, lab, thr, self.nt)


def run_arm(arm, real, E, lag):
    out = {"arm": arm}
    pool, gains = Pool(), []
    splits = {"greedy0": 0, "m6@23.8": 0}     # SCENES with at least one false split (descriptive, criterion 2)
    for c in range(NFRESH):
        rs, rn = np.random.default_rng([SEED, 2, c]), np.random.default_rng([SEED, 12, c])
        srcs = null_scene(rs)
        _s, nmf = noise_for(arm, rn, real)
        Y, F = scene(srcs, nmf)
        q = picker(F, rs)
        res = run_methods(E, Y, q, lag)
        pool.add(res, srcs)
        g = [x["pair_gain"] for x in MB.spurious(res["m6"], srcs) if x["pair_gain"] is not None]
        gains += [max(g)] if g else []
        splits["greedy0"] += len(MB.spurious(res["greedy0"], srcs)) > 0
        splits["m6@23.8"] += len(MB.spurious(adcg_loop(E, Y, q, KMAX, lag, thr=THR_REF, near_km=NEAR_KM),
                                             srcs)) > 0
    out.update(splits=splits, n_gain=len(gains), gain_med=float(np.median(gains)) if gains else np.nan,
               gain_q95=float(np.quantile(gains, .95)) if gains else np.nan,
               gain_max=float(max(gains)) if gains else np.nan)
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
                ids = pool.add(run_methods(E, Y, q, lag), srcs)
                rows.append({"stratum": sname, "dtau": dtau, "rep": rep, "slot": slot, "ids": ids,
                             "kmec": kmec_over_e(F[0][1], F[1][1])})
    hours = (NFRESH + len(rows)) * SCAN_S / 3600.0
    out["hours"] = hours
    out["op"] = {}
    for k in NAMES:
        thr, nf, hit = pool.hit(k, hours)
        out["op"][k] = (thr, nf)
        for r in rows:
            r[k] = int(hit[r["ids"]].sum())
    out["rows"] = rows
    return out


def report(o, cut):
    rows = o["rows"]
    print(f"== arm {o['arm']}: {NFRESH} fresh control scenes with neighbours + {len(rows)} pair scenes, "
          f"{o['hours']:.3f} h scanned ({SCAN_S} s of origin times per scene), budget {BUDGET:g}/h "
          f"(at most {int(np.floor(BUDGET * o['hours'] + 1e-9))} false detections)")
    print("   operating point (threshold, false detections): "
          + ", ".join(f"{k} {t:.2f} ({n})" for k, (t, n) in o["op"].items()))
    print(f"   control: {o['n_gain']} scenes with a spurious near atom; its largest pair gain median / q95 / max "
          f"{o['gain_med']:.1f} / {o['gain_q95']:.1f} / {o['gain_max']:.1f}")
    print(f"   descriptive: scenes with a false split, of {NFRESH}: greedy0 {o['splits']['greedy0']}, "
          f"m6 at thr {THR_REF} {o['splits']['m6@23.8']} (reference {SPLIT_MAX} of 200, not in the gate)")
    print(f"   pair recall at the operating point")
    print(f"   {'stratum':7s} {'dtau':>5s} | all pairs " + " ".join(f"{k:>8s}" for k in NAMES)
          + " | upper tercile (n) " + " ".join(f"{k:>8s}" for k in NAMES))
    for sname in STRATA:
        for dtau in DTAUS:
            rr = [r for r in rows if r["stratum"] == sname and r["dtau"] == dtau]
            up = [r for r in rr if r["kmec"] >= cut]
            f = lambda R, k: np.mean([r[k] == 2 for r in R]) if R else np.nan  # noqa: E731
            print(f"   {sname:7s} {dtau:5.2f} |           " + " ".join(f"{f(rr, k):8.2f}" for k in NAMES)
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
        raise SystemExit("usage: run.py NOISE_DIR [NREP] [NFRESH]")
    XN = S.nodes()
    taus = np.round(np.arange(4.0, 12.0 + 1e-9, 0.1), 10)
    assert abs(len(taus) * 0.1 - SCAN_S) < 1e-9
    tt, _ = homogeneous(XN, STA, S.VP, S.VS)
    Gs, Qs = node_operators(XN, STA)
    E = OnsetEngine(XN, STA, tt, taus, OnsetDelay(np.full((NSTA, 2), -S.LAG)), Gs, Qs, rounding="nearest")
    lag = np.full((NSTA, 2), S.LAG / S.DT)
    print(f"M6A-REALNOISE  seed {SEED}  {NREP} pairs per cell, {NFRESH} fresh controls  numpy {np.__version__}")
    print(f"noise day {real.meta.get('day')} ({NOISE_DIR}), stations {','.join(real.stations)}, "
          f"{len(real.slots)} eligible slots")
    print(f"SNR {SNR} (median strong-block peak / robust sd of the scene's filtered noise), kmax {KMAX}, "
          f"near_km {NEAR_KM}, LOC+T-V, budget {BUDGET:g} false/h on the exact FROC")
    res = {arm: run_arm(arm, real, E, lag) for arm in ("REAL", "GEN")}
    ke = np.array([r["kmec"] for r in res["REAL"]["rows"]])
    cut = float(np.nanquantile(ke, 2 / 3))
    print(f"K_mec/E of the pairs: terciles {np.nanquantile(ke, 1 / 3):.3f} / {cut:.3f}; same sources in both arms\n")
    for arm in ("REAL", "GEN"):
        report(res[arm], cut)
        print()
    used = [r["slot"] for r in res["REAL"]["rows"]]
    print(f"REAL slots used: {len(used)} by pairs + {NFRESH} by controls, {real.k} in all, each once")
    gate = [r for r in res["REAL"]["rows"] if r["kmec"] >= cut]
    d = np.array([(r["m6"] == 2) - (r["greedy"] == 2) for r in gate], float)
    lo = np.quantile(d[np.random.default_rng([SEED, 7, len("upper tercile")]).integers(0, len(d), (B, len(d)))]
                     .mean(1), .025)
    ok = d.mean() >= 0.10 and lo > 0
    print(f"G-M6a (REAL, upper tercile, n {len(d)}, at {BUDGET:g} false/h): m6 - greedy {d.mean():+.3f}, "
          f"2.5th pct {lo:+.3f} -> {'PASS' if ok else 'FAIL'} by the formula; NOT a gate reading until CR-30 and the "
          f"per-station S split are closed (card.md). Descriptive: m6 scenes with a false split at thr {THR_REF}: "
          f"{res['REAL']['splits']['m6@23.8']} of {NFRESH} (reference {SPLIT_MAX} of 200)")
    print(f"\n{time.time() - t0:.0f} s")


STA = None
if __name__ == "__main__":
    sys.exit(main())
