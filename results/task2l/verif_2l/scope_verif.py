"""scope_verif.py - exploratory scope from the stored parts: per cell (own bootstrap), term error vs gain per network,
networks where S2 is worse than S0."""
import numpy as np
from load_verif import *
W, N = load_all()
A, B_ = "MT@S2+T3", "MT@S0+T3"

def recs_of(ws, m, tol):
    return [window_rec(w["atoms"][m], label(w["atoms"][m], w["truth"], w["taus"], tol)) for w in ws]

def boot(ws, tol, B, seed):
    ks = np.array([w["k"] for w in ws], float)
    ra, rb = recs_of(ws, A, tol), recs_of(ws, B_, tol)
    strata = {}
    for i, w in enumerate(ws):
        strata.setdefault((w["reg"], w["cell"]), {}).setdefault(w["net"], []).append(i)
    rng = np.random.default_rng(seed)
    one = np.ones(len(ws))
    pa, pb = froc(ra, one, ks.sum())[0][1], froc(rb, one, ks.sum())[0][1]
    d = []
    for _ in range(B):
        mult = np.zeros(len(ws))
        for nets in strata.values():
            keys = list(nets)
            for g in rng.integers(0, len(keys), len(keys)):
                wl = np.asarray(nets[keys[g]])
                np.add.at(mult, wl[rng.integers(0, len(wl), len(wl))], 1.0)
        nte = float((mult * ks).sum())
        d.append(froc(ra, mult, nte)[0][1] - froc(rb, mult, nte)[0][1])
    return pb, pa, np.percentile(d, [2.5, 97.5])

print("# 1. Per cell (and per cell x regime): MT@S2+T3 - MT@S0+T3, cap 0.25, operating point per subset; own bootstrap 400 reps")
for kind in ("std", "dense"):
    for tol, rule in ((0.5, "LOC+T"), (None, "LOC")):
        for cell in CELLS:
            row = []
            for seed in (7171, 9393):
                ws = [w for w in W if w["seed"] == seed and w["kind"] == kind and w["cell"] == cell]
                pb, pa, ci = boot(ws, tol, 400, [5, seed, CELLS.index(cell), kind == "std", tol is None])
                row.append(f"{seed}: {pb:.3f}->{pa:.3f} {pa-pb:+.3f} [{ci[0]:+.3f}, {ci[1]:+.3f}]")
            print(f"{kind:5s} {rule:5s} {cell:15s} " + " | ".join(row), flush=True)
for cell in CELLS:
    for reg in REGIMES:
        row = []
        for seed in (7171, 9393):
            ws = [w for w in W if w["seed"] == seed and w["kind"] == "std" and w["cell"] == cell and w["reg"] == reg]
            pb, pa, ci = boot(ws, 0.5, 400, [6, seed, CELLS.index(cell), reg == "thin"])
            row.append(f"{seed}: {pb:.3f}->{pa:.3f} {pa-pb:+.3f} [{ci[0]:+.3f}, {ci[1]:+.3f}]")
        print(f"std   LOC+T {cell:15s} {reg:6s} " + " | ".join(row), flush=True)

print("\n# 2. Per-network decomposition at the POOLED operating points (thresholds of the pooled FROC at cap 0.25, per seed and window type, LOC+T)")
print("#    hits = events detected above the method's pooled threshold; d = (hits S2+T3 - hits S0+T3); networks split by the sign of the")
print("#    change of the term rms error (err S2 - err S0, median removed, vs the true statics: truth used only for this split)")
for seed in (7171, 9393):
    for kind in ("std", "dense"):
        ws = [w for w in W if w["seed"] == seed and w["kind"] == kind]
        ks = np.array([w["k"] for w in ws], float)
        hits, fal = {}, {}
        for m in (A, B_):
            rc = recs_of(ws, m, 0.5)
            thr = froc(rc, np.ones(len(ws)), ks.sum())[1][1]
            hits[m] = np.array([float((r[2] & (r[0] >= thr)).sum()) for r in rc])
            fal[m] = np.array([float((r[1] & (r[0] >= thr)).sum()) for r in rc])
        key = [(w["reg"], w["cell"], w["net"]) for w in ws]
        nets = sorted(set(key))
        idx = {k: i for i, k in enumerate(nets)}
        gi = np.array([idx[k] for k in key])
        dH = np.bincount(gi, hits[A] - hits[B_], len(nets))
        nE = np.bincount(gi, ks, len(nets))
        de = np.array([N[(seed,) + k]["err"]["S2"] - N[(seed,) + k]["err"]["S0"] for k in nets])
        cellv = np.array([k[1] for k in nets])
        print(f"{seed} {kind}: total hits S0+T3 {hits[B_].sum():.0f} S2+T3 {hits[A].sum():.0f} (diff {dH.sum():+.0f} of {nE.sum():.0f} events = {dH.sum()/nE.sum():+.4f}); "
              f"false/window S0+T3 {fal[B_].mean():.3f} S2+T3 {fal[A].mean():.3f}")
        worse = de > 0
        for nm, msk in (("S2 terms WORSE than S0", worse), ("S2 terms better", ~worse)):
            print(f"    {nm:24s}: {msk.mean():.3f} of networks ({msk.sum()}), events {nE[msk].sum():.0f}, recall change {dH[msk].sum()/max(1,nE[msk].sum()):+.4f} "
                  f"(networks with net loss {np.mean(dH[msk] < 0):.2f}, net gain {np.mean(dH[msk] > 0):.2f}, tie {np.mean(dH[msk] == 0):.2f})")
        for cell in CELLS:
            mc = cellv == cell
            print(f"    {cell:15s}: S2 worse in {worse[mc].mean():.3f}; recall change where worse {dH[mc & worse].sum()/max(1,nE[mc & worse].sum()):+.4f}, where better {dH[mc & ~worse].sum()/max(1,nE[mc & ~worse].sum()):+.4f}; cell total {dH[mc].sum()/nE[mc].sum():+.4f}")
        # dose-response: quartiles of error reduction
        q = np.quantile(-de, [0.25, 0.5, 0.75])
        b = np.digitize(-de, q)
        print("    quartiles of error reduction (err S0 - err S2): " + " | ".join(
            f"Q{i+1} mean red. {(-de)[b==i].mean():+.3f} s: recall change {dH[b==i].sum()/nE[b==i].sum():+.4f}" for i in range(4)))
        r = dH / np.maximum(nE, 1)
        rk = lambda v: np.argsort(np.argsort(v)).astype(float)
        print(f"    Spearman(error reduction, per-network recall change) = {np.corrcoef(rk(-de), rk(r))[0,1]:+.3f} over {len(nets)} networks")
