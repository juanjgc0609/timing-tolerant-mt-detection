"""M7-TIMESLIDE, step 0 of M7 (v2 specification M7.1; theory D.3): does the time-slide null hold on real noise?

Question. Theorem D.3 makes the time-slide p-value exact under "no coherent source": every station's noise in a
window is circularly shift-invariant and the stations are independent. Day 101231 (seed 2718) has transients that
coincide across stations (M5-EVALUE: tail blocks 1.7x the independent rate). Independent per-station slides
scramble that coherence. If the coherence is in the noise windows, the real (unslid) maxima are larger than the slid
ones and the time-slide null is anti-conservative there. Measured here before M7 is built on it.

Design (fixed before the run; card.md):
  networks   the Task 2M REAL-arm networks of results/task2M/cfg_SEED.json (code/t2m_lib.py unchanged).
  slots      every slot of the network's noise block except its catalogue slots (no fit: slides need none).
  detector   the registered MT@S0 first atom: c10_lib.mt_loop(E, MFn, 1, mode="hard", pr, W=4), E the coarse
             lattice with the catalogue station corrections, pr = E.prior(sigma_catalogue). T = its score s.
  slides     B per slot: each station's MFn (all three channels) rolled by an independent uniform offset in
             [0, NT), seeded by (seed, tag, cell, net, slot, b), independent of the data. sd of MFn is global, so a
             roll leaves it unchanged. p = (1 + #{T_b >= T}) / (B + 1).
  positive   the catalogue windows of the network (synthetic events in real noise): p must be small (power check).
Reported: P(p <= a) for a in {0.01, 0.05, 0.1, 0.2} (exact null: <= a), with 95 % bounds by a bootstrap of
networks (the three cells share noise: clusters are noise blocks); the pooled tail: the rate of real T above the
pooled slid quantiles q(1 - a), a in {1e-2, 1e-3}; and per cell.

    python run_ts.py fake OUTDIR [nsta] [hours] [seed]                 (container smoke: i.i.d. t4 noise day)
    python run_ts.py run SEED NOISE_DIR[,..] OUTDIR B [nets|all] [cells] [--procs N]
    python run_ts.py analyze OUTDIR[,OUTDIR] [nboot]
"""
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
for p in (os.path.join(ROOT, "code"), os.path.join(ROOT, "code", "explore_contrib"), os.path.join(ROOT, "src")):
    if p not in sys.path:
        sys.path.insert(0, p)
import numpy as np  # noqa: E402

import t2m_lib as L2  # noqa: E402
import c10_lib as C  # noqa: E402
import c10_run as R  # noqa: E402
from task1o_consistent import matched  # noqa: E402
from ttmt.v2.m7 import slide_stats  # noqa: E402

T2M = os.path.join(ROOT, "results", "task2M")
TAG_SLIDE = 7701
ALPHAS = (0.01, 0.05, 0.1, 0.2)
TAILS = (1e-2, 1e-3)


def configure(seed, ndirs):
    fn = str(seed) if str(seed).endswith(".json") else os.path.join(T2M, f"cfg_{seed}.json")
    with open(fn) as fh:
        c = json.load(fh)
    L2.configure(c["seed"], ndirs, c["j0"], c["nblocks"], c["nwin"], ["MT@S0"])
    return c


def detector(net):
    E = C.Eng(R.GARR, net.STA, net.corr, R.TAUS, C.H_COARSE)
    sig, _cnt = C.sigma_catalogue(net.STA, net.cat, net.corr)
    pr = E.prior(sig)

    def stat(MFn, atom=False):
        a = C.mt_loop(E, MFn, 1, mode="hard", pr=pr, W=4)
        s = float(a[0]["s"]) if a else -np.inf
        if atom:
            return s, (np.r_[a[0]["x"], a[0]["tau"]] if a else np.full(4, np.nan))
        return s
    return stat


def energy_track(MFn):
    """per station: sample and value of the maximum of the three-channel energy (coherence diagnostic)."""
    e = (MFn * MFn).sum(1)
    return e.argmax(-1), e.max(-1)


def window_test(stat, MFn, B, rng):
    return slide_stats(MFn, stat, B, rng)


def run_part(args):
    seed, ndirs, cell, j, B, odir = args
    configure(seed, ndirs)
    fn = os.path.join(odir, f"{cell}_{j}.npz")
    if os.path.exists(fn):
        return fn
    t0 = time.time()
    net = L2.Net2M(cell, j, "REAL")
    stat = detector(net)
    _di, blk = L2.block(j)
    used = set(int(s) for s in net.cat_slots)
    slots = [int(s) for s in blk if int(s) not in used]
    T, TB, X, KT, KE = [], [], [], [], []
    for s in slots:
        MF, sd = matched(L2.real_segment(net.di, net.sub, s))
        rng = np.random.default_rng([int(L2.CFG["seed"]), TAG_SLIDE, net.ci, j, s])
        t, tb = window_test(stat, MF / sd, B, rng)
        T.append(t)
        TB.append(tb)
        X.append(stat(MF / sd, atom=True)[1])
        kt, ke = energy_track(MF / sd)
        KT.append(kt)
        KE.append(ke)
    PT, PTB = [], []                                       # positive control: catalogue windows (events in noise)
    for k, (D, _src, _tau) in enumerate(net.cat):
        MF, sd = matched(D)
        rng = np.random.default_rng([int(L2.CFG["seed"]), TAG_SLIDE + 1, net.ci, j, k])
        t, tb = window_test(stat, MF / sd, B, rng)
        PT.append(t)
        PTB.append(tb)
    np.savez(fn, T=np.array(T), TB=np.array(TB), slots=np.array(slots), X=np.array(X), KT=np.array(KT),
             KE=np.array(KE), PT=np.array(PT), PTB=np.array(PTB),
             cell=cell, j=j, ns=net.ns, day=str(net.day), secs=time.time() - t0)
    print(f"{cell} {j}: {len(slots)} slots, {len(PT)} cat, {time.time() - t0:.0f} s", flush=True)
    return fn


def run(seed, ndirs, odir, B, nets, cells, procs):
    c = configure(seed, ndirs)
    os.makedirs(odir, exist_ok=True)
    nets = L2.CFG["nets"] if nets in (None, "all") else [int(x) for x in nets.split(",")]
    cells = list(L2.CELLS) if not cells else cells.split(",")
    jobs = [(seed, ndirs, cell, j, B, odir) for j in nets for cell in cells]
    with open(os.path.join(odir, "config.json"), "w") as fh:
        json.dump({"cfg": c, "noise_dirs": ndirs, "B": B, "nets": nets, "cells": cells}, fh, indent=1)
    if procs > 1:
        from multiprocessing import Pool
        with Pool(procs) as pool:
            list(pool.imap_unordered(run_part, jobs))
    else:
        for jb in jobs:
            run_part(jb)


def pvals(T, TB):
    return (1.0 + (TB >= T[:, None]).sum(1)) / (TB.shape[1] + 1.0)


def load(odirs):
    parts = []
    for d in odirs:
        for f in sorted(os.listdir(d)):
            if f.endswith(".npz"):
                z = np.load(os.path.join(d, f))
                parts.append({k: z[k] for k in z.files} | {"src": os.path.basename(d.rstrip("/"))})
    return parts


def summary(parts, nboot, rng):
    """rates of p <= a and pooled tail rates; bootstrap over noise blocks (src, j)."""
    keys = sorted({(str(p["src"]), int(p["j"])) for p in parts})
    by = {k: [p for p in parts if (str(p["src"]), int(p["j"])) == k] for k in keys}

    def stats(sel):
        T = np.concatenate([p["T"] for p in sel])
        TB = np.concatenate([p["TB"] for p in sel])
        pv = pvals(T, TB)
        out = [float(np.mean(pv <= a)) for a in ALPHAS]
        pooled = np.sort(TB.ravel())
        for a in TAILS:
            q = pooled[int(np.ceil((1 - a) * len(pooled))) - 1]
            out.append(float(np.mean(T > q)))
        return out, len(T)
    est, n = stats(parts)
    bs = []
    for _ in range(nboot):
        pick = rng.integers(len(keys), size=len(keys))
        bs.append(stats([p for i in pick for p in by[keys[i]]])[0])
    bs = np.array(bs)
    return est, n, np.quantile(bs, 0.025, axis=0), np.quantile(bs, 0.975, axis=0)


def analyze(odirs, nboot=500):
    parts = load(odirs)
    rng = np.random.default_rng(7702)
    hdr = ["p<=%.2g" % a for a in ALPHAS] + ["tail %.0e" % a for a in TAILS]
    print("columns (nominal):", ", ".join(f"{h} ({a:g})" for h, a in zip(hdr, ALPHAS + TAILS)))
    groups = [("ALL", parts)] + [(c, [p for p in parts if str(p["cell"]) == c]) for c in L2.CELLS]
    groups += [(f"src {s}", [p for p in parts if str(p["src"]) == s]) for s in sorted({str(p["src"]) for p in parts})]
    for name, sel in groups:
        if not sel:
            continue
        est, n, lo, hi = summary(sel, nboot, rng)
        print(f"{name:16s} n={n:6d} " + "  ".join(f"{e:.4f} [{a:.4f},{b:.4f}]" for e, a, b in zip(est, lo, hi)))
    coherence(parts)
    PT = np.concatenate([p["PT"] for p in parts])
    PTB = np.concatenate([p["PTB"] for p in parts])
    if len(PT):
        pv = pvals(PT, PTB)
        print(f"positive control (catalogue windows) n={len(PT)}: P(p<=0.01)={np.mean(pv <= 0.01):.3f}, "
              f"median p={np.median(pv):.4f}")
    secs = sum(float(p["secs"]) for p in parts)
    print(f"parts {len(parts)}, CPU time {secs / 60:.1f} min")


def coherence(parts, a=0.01, tol_s=1.0, k=3.0):
    """diagnostic of the windows with p <= a against all windows: share of stations whose energy maximum exceeds k x
    that station's median over the part ('loud'), and whether the loud stations peak within +-tol_s of each other
    (simultaneous: a zero-moveout transient rather than an earthquake, whose moveout is seconds)."""
    if not all("KT" in p for p in parts):
        return
    rows = []
    for p in parts:
        pv = pvals(p["T"], p["TB"])
        med = np.median(p["KE"], axis=0)
        for w in range(len(pv)):
            loud = p["KE"][w] > k * med
            kt = p["KT"][w][loud]
            sim = loud.sum() >= 2 and (kt.max() - kt.min()) * L2.DT <= 2 * tol_s
            rows.append((pv[w] <= a, loud.mean(), loud.sum() >= 2, sim))
    r = np.array(rows, float)
    for name, sel in (("all windows", np.ones(len(r), bool)), (f"p <= {a:g}", r[:, 0] > 0)):
        if sel.sum():
            print(f"coherence {name:12s} n={int(sel.sum()):6d}: loud share {r[sel, 1].mean():.3f}, "
                  f">=2 loud {r[sel, 2].mean():.3f}, simultaneous (+-{tol_s:g} s) {r[sel, 3].mean():.3f}")


def fake(odir, nsta=10, hours=4.0, seed=1):
    """an i.i.d. Student t4 noise day in the t2m_data.py format (container smoke; a null with no coherence)."""
    os.makedirs(odir, exist_ok=True)
    rng = np.random.default_rng(seed)
    NT = L2.NT
    nsl = int(hours * 3600 / (NT * L2.DT))
    x = rng.uniform(-20, 20, nsta)
    y = rng.uniform(-20, 20, nsta)
    out = rng.standard_t(4, size=(nsta, 3, nsl * NT))
    _MF, lev = matched(out[:, :, :50 * NT])                  # robust MF level -> 0.14, as t2m_data.py (sd_gen)
    out = (out * (0.14 / lev)).astype(np.float32)
    np.save(os.path.join(odir, "noise.npy"), out)
    meta = {"stations": [f"F{i:02d}" for i in range(nsta)], "x_km": x.tolist(), "y_km": y.tolist(),
            "n_slots": nsl, "eligible": list(range(nsl)), "day": "fake", "slot_s": NT * L2.DT}
    with open(os.path.join(odir, "meta.json"), "w") as fh:
        json.dump(meta, fh)
    print(odir, nsl, "slots")


if __name__ == "__main__":
    a = sys.argv[1:]
    procs = 1
    if "--procs" in a:
        i = a.index("--procs")
        procs = int(a[i + 1])
        del a[i:i + 2]
    if a[0] == "fake":
        fake(a[1], *(int(a[2]),) if len(a) > 2 else (), *(float(a[3]),) if len(a) > 3 else ())
    elif a[0] == "run":
        run(a[1], a[2].split(","), a[3], int(a[4]), a[5] if len(a) > 5 else None, a[6] if len(a) > 6 else None,
            procs)
    elif a[0] == "analyze":
        analyze(a[1].split(","), int(a[2]) if len(a) > 2 else 500)
