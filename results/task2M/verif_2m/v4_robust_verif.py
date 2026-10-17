"""V4: robustness of the 2718 primary contrasts (DENSE, @0.25): hit tolerance, cells, hour blocks, leave-one-network-out,
leave-one-hour-block-out, greedy removal of the most favourable networks, networks without a given station, quiet windows,
atom truncation.  Own FROC and bootstrap (vlib_verif)."""
import sys
import time

import numpy as np

sys.path.insert(0, "/home/claude/work/t2m_verif")
import vlib_verif as V  # noqa: E402

OUT = "/home/claude/work/t2m_verif/v4_robust_verif.txt"
lines = []


def P(s=""):
    print(s, flush=True)
    lines.append(s)
    with open(OUT, "w") as fh:
        fh.write("\n".join(lines) + "\n")


A, H, Bm = "MT@S0", "HUBER-N@S0", "BOMP@S0"
nets = V.load_parts("2718")
wins = V.windows(nets, "DENSE")
stas = {(n["cell"], n["net"]): n["stations"] for n in nets}


def contrasts(prep, ws, B=0, sub=None, seed=99):
    r = {m: V.recall_at(*prep[m], ws)[1] for m in (A, H, Bm)}
    out = [r[A] - r[Bm], r[H] - r[Bm], r[A], r[H], r[Bm]]
    if B:
        ww = [wins[i] for i in ws]
        reps = V.boot_index(ww, B, seed)
        d1, d2 = [], []
        for rep in reps:
            g = ws[rep]
            rr = {m: V.recall_at(*prep[m], g)[1] for m in (A, H, Bm)}
            d1.append(rr[A] - rr[Bm])
            d2.append(rr[H] - rr[Bm])
        out += [np.percentile(d1, 2.5), np.percentile(d2, 2.5)]
    return out


def fmt(c):
    s = f"MT-BOMP {c[0]:+.3f}  HUB-BOMP {c[1]:+.3f}  (MT {c[2]:.3f} HUB {c[3]:.3f} BOMP {c[4]:.3f})"
    if len(c) > 5:
        s += f"  lower 2.5%: {c[5]:+.3f} / {c[6]:+.3f}"
    return s


t0 = time.time()
allw = np.arange(len(wins))
P("V4 robustness, seed 2718, REAL, DENSE, @0.25 (bootstrap lower 2.5 percentiles with own hierarchical bootstrap)")
P("1) hit tolerance (LOC rule, epicentre within tol)")
for tol in (2.0, 3.0, 5.0, 7.0, 10.0):
    prep = {m: V.prep(wins, m, tol, None) for m in (A, H, Bm)}
    P(f"   tol {tol:4.1f} km: " + fmt(contrasts(prep, allw, 1000)))
prep = {m: V.prep(wins, m, 5.0, None) for m in (A, H, Bm)}
P("2) per cell (LOC 5 km)")
for c in V.CELLS:
    ws = np.array([i for i, w in enumerate(wins) if w["cell"] == c])
    P(f"   {c:15s}: " + fmt(contrasts(prep, ws, 1000)))
P("3) per block of networks (hours of the day; cells share hours)")
for lo, hi in ((2800, 2808), (2808, 2816), (2816, 2824)):
    ws = np.array([i for i, w in enumerate(wins) if lo <= w["net"] < hi])
    hrs = [w["slot"] * 40 / 3600 for i, w in enumerate(wins) if lo <= w["net"] < hi]
    P(f"   nets {lo}-{hi - 1} (hours {min(hrs):.1f}-{max(hrs):.1f}): " + fmt(contrasts(prep, ws, 1000)))
P(f"   ({time.time() - t0:.0f} s)")
# 4) leave-one-network-out and leave-one-hour-block-out
units = sorted({(w["cell"], w["net"]) for w in wins})
P("4) leave-one-network-out (72 (cell, net) units): range of the point contrasts")
lono = []
for u in units:
    ws = np.array([i for i, w in enumerate(wins) if (w["cell"], w["net"]) != u])
    lono.append(contrasts(prep, ws))
lono = np.array(lono)
for k, nm in ((0, "MT-BOMP"), (1, "HUB-BOMP")):
    P(f"   {nm}: min {lono[:, k].min():+.3f} ({units[int(lono[:, k].argmin())]}), max {lono[:, k].max():+.3f}")
P("   leave-one-hour-block-out (drop network index j in all 3 cells; 24 options):")
loho = []
for j in range(2800, 2824):
    ws = np.array([i for i, w in enumerate(wins) if w["net"] != j])
    loho.append(contrasts(prep, ws))
loho = np.array(loho)
for k, nm in ((0, "MT-BOMP"), (1, "HUB-BOMP")):
    P(f"   {nm}: min {loho[:, k].min():+.3f} (j {2800 + int(loho[:, k].argmin())}), max {loho[:, k].max():+.3f}")
# 5) greedy removal of the most favourable networks
P("5) greedy removal of the (cell, net) unit whose removal most reduces the contrast, until the point < +0.03")
for k, nm in ((0, "MT-BOMP"), (1, "HUB-BOMP")):
    keep = list(units)
    path = []
    cur = contrasts(prep, allw)[k]
    while cur >= 0.03 and len(keep) > 3:
        best = None
        for u in keep:
            ks = set(keep) - {u}
            ws = np.array([i for i, w in enumerate(wins) if (w["cell"], w["net"]) in ks])
            v = contrasts(prep, ws)[k]
            if best is None or v < best[0]:
                best = (v, u)
        keep.remove(best[1])
        cur = best[0]
        path.append((best[1], round(cur, 3)))
        if len(path) in (1, 3, 5, 10, 15, 20, 30):
            P(f"   {nm}: after removing {len(path):2d} units point {cur:+.3f}")
    P(f"   {nm}: units removed to bring the point below +0.03: {len(path)} of 72 (remaining {len(keep)}); "
      f"first five {path[:5]}")
    hours = [u[1] for u, _ in path]
    P(f"   {nm}: removed units by hour block: 2800-07 {sum(h < 2808 for h in hours)}, 2808-15 "
      f"{sum(2808 <= h < 2816 for h in hours)}, 2816-23 {sum(h >= 2816 for h in hours)}")
P(f"   ({time.time() - t0:.0f} s)")
# 6) networks without a given station
P("6) only the networks whose station subset does NOT contain station X (bootstrap B = 500)")
allst = sorted({s for v in stas.values() for s in v})
for s in allst:
    ws = np.array([i for i, w in enumerate(wins) if s not in stas[(w["cell"], w["net"])]])
    nu = len({(wins[i]["cell"], wins[i]["net"]) for i in ws})
    if nu < 6:
        P(f"   without {s}: only {nu} units")
        continue
    P(f"   without {s} ({nu:2d} units): " + fmt(contrasts(prep, ws, 500)))
# 7) quiet windows only
P("7) windows by noise character (window max|MF|/sd of the noise part)")
mx = np.array([w["maxabs"] for w in wins])
for lab, sel in (("max <= 10", mx <= 10), ("max <= 6", mx <= 6), ("max > 10", mx > 10)):
    ws = np.nonzero(sel)[0]
    P(f"   {lab:9s} ({len(ws)} windows): " + fmt(contrasts(prep, ws, 500)))
# 8) atom truncation
P("8) only the first k atoms of each window")
for kc in (3, 5):
    pk = {m: V.prep(wins, m, 5.0, None, kc) for m in (A, H, Bm)}
    P(f"   k <= {kc}: " + fmt(contrasts(pk, allw)))
P(f"done ({time.time() - t0:.0f} s)")
