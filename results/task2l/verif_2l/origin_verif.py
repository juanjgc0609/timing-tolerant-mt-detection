"""origin_verif.py - origin-time checks: signed origin errors of epicentral hits, hit windows 0.3/0.5/0.75/1.0 s,
median-debiased origin error (global per method, per cell x regime per method, per network per method), offsets."""
import numpy as np
from load_verif import *
W, N = load_all()
A, B_, H0 = "MT@S2+T3", "MT@S0+T3", "MT@S0"

def sel(seed, kind, cell=None, reg=None):
    return [w for w in W if w["seed"] == seed and w["kind"] == kind and (cell is None or w["cell"] == cell) and (reg is None or w["reg"] == reg)]

def signed(ws, m):
    """signed origin error (atom - true) of every atom whose nearest epicentre is within 5 km (LOC hits)."""
    out = []
    for w in ws:
        at = w["atoms"][m]
        lab = label(at, w["truth"], w["taus"], None)
        out.append([(at[i, 3] - w["taus"][h]) for i, h in enumerate(lab) if h >= 0])
    return out

def rec_point(ws, m, tol, shifts=None, cap_i=1):
    ks = np.array([w["k"] for w in ws], float)
    recs = [window_rec(w["atoms"][m], label(w["atoms"][m], w["truth"], w["taus"], tol, shift=0.0 if shifts is None else shifts[i]))
            for i, w in enumerate(ws)]
    return froc(recs, np.ones(len(ws)), ks.sum())[0][cap_i]

print("# 1. Signed origin-time error (atom origin - true origin, s) of LOC hits (all atoms, epicentre within 5 km)")
for seed in (7171, 9393):
    for kind in ("std", "dense"):
        for cell in (None,) + CELLS:
            ws = sel(seed, kind, cell)
            row = []
            for m in (H0, B_, A) + (("MT@ORA+T3",) if cell is None else ()):
                wsm = [w for w in ws if m in w["atoms"]]
                e = np.concatenate([np.asarray(x, float) for x in signed(wsm, m)])
                row.append(f"{m}: n {len(e)} mean {e.mean():+.3f} median {np.median(e):+.3f} q10 {np.quantile(e,.1):+.2f} q90 {np.quantile(e,.9):+.2f} "
                           f"|e|<=0.3 {np.mean(np.abs(e)<=0.3+1e-9):.3f} <=0.5 {np.mean(np.abs(e)<=0.5+1e-9):.3f} <=0.75 {np.mean(np.abs(e)<=0.75+1e-9):.3f} <=1.0 {np.mean(np.abs(e)<=1.0+1e-9):.3f}")
            print(f"{seed} {kind} {cell or 'all cells'}:\n    " + "\n    ".join(row))

print("\n# 2. Contrast MT@S2+T3 - MT@S0+T3 at cap 0.25 (and 0.1) under hit windows of 0.3 / 0.5 / 0.75 / 1.0 s and LOC only")
for seed in (7171, 9393):
    for kind in ("std", "dense"):
        ws = sel(seed, kind)
        row = []
        for tol in (0.2, 0.3, 0.5, 0.75, 1.0, None):
            a, b = rec_point(ws, A, tol), rec_point(ws, B_, tol)
            a1, b1 = rec_point(ws, A, tol, cap_i=0), rec_point(ws, B_, tol, cap_i=0)
            row.append(f"tol {tol}: {b:.3f} -> {a:.3f} ({a-b:+.4f}; cap 0.1 {a1-b1:+.4f})")
        print(f"{seed} {kind}: " + " | ".join(row))

print("\n# 3. Median-debiased origin error, hit window 0.5 s: the method's median signed origin error of its LOC hits is removed")
print("#    (a) one median per method (pooled), (b) per method x regime x cell, (c) per method x network (truth-aided, DIAGNOSTIC ONLY)")
for seed in (7171, 9393):
    for kind in ("std", "dense"):
        ws = sel(seed, kind)
        res = {}
        for m in (A, B_):
            sg = signed(ws, m)
            allv = np.concatenate([np.asarray(x, float) for x in sg])
            g = float(np.median(allv))
            byc, byn = {}, {}
            for w, x in zip(ws, sg):
                byc.setdefault((w["reg"], w["cell"]), []).extend(x)
                byn.setdefault((w["reg"], w["cell"], w["net"]), []).extend(x)
            sc = [float(np.median(byc[(w["reg"], w["cell"])])) for w in ws]
            sn = [float(np.median(byn[(w["reg"], w["cell"], w["net"])])) if len(byn[(w["reg"], w["cell"], w["net"])]) >= 3 else g for w in ws]
            res[m] = (rec_point(ws, m, 0.5), rec_point(ws, m, 0.5, [g] * len(ws)), rec_point(ws, m, 0.5, sc), rec_point(ws, m, 0.5, sn), g)
        print(f"{seed} {kind}: raw {res[B_][0]:.3f}->{res[A][0]:.3f} ({res[A][0]-res[B_][0]:+.4f}) | (a) medians S0+T3 {res[B_][4]:+.3f}, S2+T3 {res[A][4]:+.3f}: "
              f"{res[B_][1]:.3f}->{res[A][1]:.3f} ({res[A][1]-res[B_][1]:+.4f}) | (b) {res[B_][2]:.3f}->{res[A][2]:.3f} ({res[A][2]-res[B_][2]:+.4f}) | "
              f"(c) {res[B_][3]:.3f}->{res[A][3]:.3f} ({res[A][3]-res[B_][3]:+.4f})")

print("\n# 4. Common offset of the terms, median(c - statics) per network (s), and its spread; rms error (median removed)")
for seed in (7171, 9393):
    for cell in CELLS:
        for reg in REGIMES:
            o0 = np.array([N[(seed, reg, cell, j)]["off"]["S0"] for j in range(160)])
            o2 = np.array([N[(seed, reg, cell, j)]["off"]["S2"] for j in range(160)])
            e0 = np.array([N[(seed, reg, cell, j)]["err"]["S0"] for j in range(160)])
            e2 = np.array([N[(seed, reg, cell, j)]["err"]["S2"] for j in range(160)])
            print(f"{seed} {cell:15s} {reg:6s}: offset S0 mean {o0.mean():+.3f} sd {o0.std():.3f} | S2 mean {o2.mean():+.3f} sd {o2.std():.3f} | "
                  f"rms err S0 median {np.median(e0):.3f} S2 median {np.median(e2):.3f} | S2 worse than S0 in {np.mean(e2 > e0):.3f} of networks")
