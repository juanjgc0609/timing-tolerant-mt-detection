"""V3: symmetry / leakage / BOMP-collapse anatomy.
 a) identical windows across methods (structural) and finiteness / ranges of all atom scores;
 b) slot mapping: noise statistics stored in rows recomputed from noise.npy for a sample of rows;
 c) which stations carry the transients of day 101231 (per station and hour, from noise.npy);
 d) BOMP false atoms: score vs window max|MF|/sd and kurtosis; location of the top false atoms relative to stations;
 e) exact FROC curves (recall vs realised false per window) per method, 2718 DENSE LOC.
"""
import json
import sys

import numpy as np

sys.path.insert(0, "/home/claude/work/t2m_verif")
sys.path.insert(0, "/home/claude/work")
import vlib_verif as V  # noqa: E402
from task1o_consistent import _KER  # noqa: E402  (the matched-filter kernel only, to reproduce stored noise stats)

OUT = "/home/claude/work/t2m_verif/v3_symmetry_verif.txt"
lines = []


def P(s=""):
    print(s, flush=True)
    lines.append(s)
    with open(OUT, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def mf(D):
    return np.array([[np.convolve(D[i, c], _KER[::-1], mode="same") for c in range(3)] for i in range(D.shape[0])])


NT = 640
M = ["MT@S0", "HUBER-N@S0", "BOMP@S0", "MT@S0+T3", "MT@S2+T3"]
for seed in ("2718", "3141"):
    nets = V.load_parts(seed)
    nonfin = {m: 0 for m in M}
    natoms = {m: 0 for m in M}
    smax = {m: 0.0 for m in M}
    nonmono = {m: 0 for m in M}
    for n in nets:
        for r in n["rows"]:
            assert set(r["atoms"]) == set(M)
            for m in M:
                a = np.array(r["atoms"][m], float).reshape(-1, 7)
                natoms[m] += len(a)
                nonfin[m] += int((~np.isfinite(a[:, :6])).sum())
                if len(a):
                    smax[m] = max(smax[m], float(np.nanmax(a[:, 4])))
                    nonmono[m] += int(np.sum(np.diff(a[:, 4]) > 1e-9))
    P(f"seed {seed}: atoms per method {natoms}; non-finite entries {nonfin}; max score {({m: round(v, 1) for m, v in smax.items()})}")
    P(f"   atoms whose score exceeds the previous atom's (non-monotone emission; record() takes running min): {nonmono}")

# b) slot mapping check
P("")
P("b) stored noise statistics (kurtosis, max|MF|/sd of the REAL noise part) recomputed from noise.npy[sub, :, slot]")
rng = np.random.default_rng(5)
for seed, ndirs in (("2718", ["101231"]), ("3141", ["100210", "100211"])):
    nets = V.load_parts(seed)
    arrs = {d: np.load(f"/home/claude/work/t2m_out/noise/{d}/noise.npy", mmap_mode="r") for d in ndirs}
    nchk, nbad, worst = 0, 0, 0.0
    for n in [nets[i] for i in rng.choice(len(nets), 12, replace=False)]:
        for r in [n["rows"][i] for i in rng.choice(len(n["rows"]), 4, replace=False)]:
            N = np.asarray(arrs[n["day"]][n["sub"], :, r["slot"] * NT:(r["slot"] + 1) * NT], float)
            x = mf(N).ravel()
            x = x / (np.median(np.abs(x)) / 0.6745)
            ku = float(np.mean(x ** 4) / np.mean(x ** 2) ** 2 - 3.0)
            mx = float(np.max(np.abs(x)))
            nchk += 1
            err = max(abs(ku - r["kurt"]), abs(mx - r["maxabs"]))
            worst = max(worst, err)
            nbad += err > 0.01
    P(f"   seed {seed}: {nchk} rows checked, mismatches > 0.01: {nbad}, worst abs diff {worst:.4f}")

# c) station anatomy of transients, day 101231
P("")
P("c) day 101231: per station, share of eligible slots whose station max|MF| exceeds 10x the station robust sd "
  "(robust sd = median over eligible slots of the per-slot robust sd), by 8-hour block")
meta = json.load(open("/home/claude/work/t2m_out/noise/101231/meta.json"))
arr = np.load("/home/claude/work/t2m_out/noise/101231/noise.npy", mmap_mode="r")
el = np.array(meta["eligible"])
ns = len(meta["stations"])
mx = np.zeros((ns, len(el)))
rs = np.zeros((ns, len(el)))
for k, s in enumerate(el):
    seg = np.asarray(arr[:, :, s * NT:(s + 1) * NT], float)
    m_ = mf(seg)
    a = np.linalg.norm(m_, axis=1)
    mx[:, k] = a.max(1)
    rs[:, k] = np.median(np.abs(m_), axis=(1, 2)) / 0.6745
ref = np.median(rs, axis=1)
hr = el * 40 / 3600
np.savez("/home/claude/work/t2m_verif/v3_station_slotstats_verif.npz", mx=mx, rs=rs, el=el)
for lo, hi in ((0, 8), (8, 16), (16, 24)):
    sel = (hr >= lo) & (hr < hi)
    P(f"   hours {lo:2d}-{hi:2d} ({sel.sum()} slots): " + "  ".join(
        f"{meta['stations'][i]} {np.mean(mx[i, sel] / ref[i] > 10):.2f}" for i in range(ns)))
P("   share of slots where the network max is on each station (hours 8-16): " + "  ".join(
    f"{meta['stations'][i]} {np.mean(np.argmax(mx[:, (hr >= 8) & (hr < 16)] / ref[:, None], 0) == i):.2f}" for i in range(ns)))

# d) BOMP false atoms vs transients, 2718 DENSE
P("")
P("d) 2718 DENSE LOC: first atom of each window, by window max|MF|/sd class: share false, median score per method")
nets = V.load_parts("2718")
wins = V.windows(nets, "DENSE")
mxw = np.array([w["maxabs"] for w in wins])
cls = [("max<=6", mxw <= 6), ("6<max<=10", (mxw > 6) & (mxw <= 10)), ("10<max<=30", (mxw > 10) & (mxw <= 30)),
       ("max>30", mxw > 30)]
for m in ("MT@S0", "HUBER-N@S0", "BOMP@S0"):
    S, Lb, K = V.prep(wins, m)
    s1 = np.array([s[0] if len(s) else np.nan for s in S])
    f1 = np.array([lb[0] < 0 if len(lb) else False for lb in Lb])
    P(f"   {m:11s} " + "  ".join(f"{c} n {sel.sum()}: false {f1[sel].mean():.2f} med s {np.nanmedian(s1[sel]):.0f}"
                                  for c, sel in cls))
# top-scored false BOMP atoms: distance to nearest station of the network, station holding the window maximum
S, Lb, K = V.prep(wins, "BOMP@S0")
recs = []
for wi, w in enumerate(wins):
    for ai, (s, lb) in enumerate(zip(S[wi], Lb[wi])):
        recs.append((s, lb, wi, ai))
recs.sort(key=lambda t: -t[0])
top = recs[:288]
nf = sum(1 for t in top if t[1] < 0)
P(f"   BOMP@S0: of the top 288 (= 0.25 x 1152) scored atoms, {nf} are false; windows of those false atoms: median "
  f"max|MF|/sd {np.median([mxw[t[2]] for t in top if t[1] < 0]):.1f}; share with max > 10: "
  f"{np.mean([mxw[t[2]] > 10 for t in top if t[1] < 0]):.2f} (all windows: {np.mean(mxw > 10):.2f})")
netd = {(n["cell"], n["net"]): np.array(n["STA"]) for n in nets}
dmin = []
for t in top:
    if t[1] < 0:
        w = wins[t[2]]
        a = np.array(w["atoms"]["BOMP@S0"][t[3]][:3])
        dmin.append(np.min(np.linalg.norm(netd[(w["cell"], w["net"])][:, :2] - a[None, :2], axis=1)))
P(f"   epicentral distance of those false atoms to the nearest network station: quartiles "
  f"{np.percentile(dmin, [25, 50, 75]).round(1).tolist()} km")
for m in ("MT@S0", "HUBER-N@S0"):
    S2, Lb2, _ = V.prep(wins, m)
    r2 = sorted([(s, lb, wi) for wi in range(len(wins)) for s, lb in zip(S2[wi], Lb2[wi])], key=lambda t: -t[0])[:288]
    P(f"   {m}: of its top 288 scored atoms {sum(1 for t in r2 if t[1] < 0)} are false; their windows share max>10: "
      f"{np.mean([mxw[t[2]] > 10 for t in r2 if t[1] < 0]) if any(t[1] < 0 for t in r2) else float('nan'):.2f}")

# e) FROC curves
P("")
P("e) exact FROC, 2718 / 3141 DENSE LOC: recall at realised false events per window <= cap")
caps = (0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0)
for seed in ("2718", "3141"):
    wins_ = V.windows(V.load_parts(seed), "DENSE")
    P(f"   seed {seed} caps " + " ".join(f"{c:5.2f}" for c in caps))
    for m in M:
        S, Lb, K = V.prep(wins_, m)
        r = V.recall_at(S, Lb, K, np.arange(len(wins_)), caps)
        P(f"   {m:11s}      " + " ".join(f"{x:5.3f}" for x in r))
