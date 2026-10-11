"""exploratory robustness of the 2I contrasts from the stored parts only (not preregistered)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats_verif import *
P = print
B2 = int(os.environ.get("VB", 400))
for seed in (5151, 6262):
    wins, _, _ = load(seed)
    for kind in ("dense", "std"):
        ww = [w for w in wins if w["kind"] == kind]
        rng = np.random.default_rng([777, seed, len(kind)])
        C = boot_counts(ww, rng, B2)
        P(f"=== seed {seed} {kind}")
        # ---- (A) strict hit rules
        for tol_t in (0.5, 1.0):
            lab = {m: [relabel(w, m, 5.0, tol_t) for w in ww] for m in METH}
            for a, b in ((H, BO), (H, MT)):
                c = contrast(ww, a, b, C, labels=lab)
                P(f"  [hit: <=5 km AND |dt|<={tol_t}s] {a} - {b} @0.25: {fmt(c)} | @0.1 {c[0.1]['pt']:+.4f} [{c[0.1]['lo']:+.4f}, {c[0.1]['hi']:+.4f}]")
        dts = {m: [] for m in METH}
        for w in ww:
            ta = np.asarray(w["taus"])
            for m in METH:
                for d in w["det"][m]:
                    if d[1] >= 0: dts[m].append(abs(d[5] - ta[d[1]]))
        P("  among location hits: fraction with |origin error| <= 0.5 s:", {m: round(float(np.mean(np.array(v) <= 0.5 + 1e-12)), 3) for m, v in dts.items()},
          " median |dt|:", {m: round(float(np.median(v)), 2) for m, v in dts.items()})
        # ---- (B) per regime, per cell (own threshold per subset)
        for reg in REGS:
            idx = [i for i, w in enumerate(ww) if w["reg"] == reg]; sub = [ww[i] for i in idx]
            for a, b in ((H, BO), (H, MT)):
                P(f"  [regime {reg}] {a} - {b} @0.25: {fmt(contrast(sub, a, b, C[:, idx]))}")
        if kind == "dense":
            for reg in REGS:
                for cell in CELLS:
                    idx = [i for i, w in enumerate(ww) if w["reg"] == reg and w["cell"] == cell]; sub = [ww[i] for i in idx]
                    c1 = contrast(sub, H, BO, C[:, idx]); c2 = contrast(sub, H, MT, C[:, idx])
                    P(f"  [cell {reg} {cell}] H-BOMP @0.25 {c1[0.25]['pt']:+.4f} [{c1[0.25]['lo']:+.4f}, {c1[0.25]['hi']:+.4f}]   H-MTFW {c2[0.25]['pt']:+.4f} [{c2[0.25]['lo']:+.4f}, {c2[0.25]['hi']:+.4f}]")
            # macro: mean of per-cell contrasts (per-cell operating point) instead of one pooled threshold
            dd = []; pp = []
            for reg in REGS:
                for cell in CELLS:
                    idx = [i for i, w in enumerate(ww) if w["reg"] == reg and w["cell"] == cell]; sub = [ww[i] for i in idx]
                    fa, fb = Froc(sub, H), Froc(sub, BO)
                    dd.append(fa.recall(C[:, idx])[0.25] - fb.recall(C[:, idx])[0.25])
                    pp.append(fa.recall(np.ones((1, len(idx))))[0.25][0] - fb.recall(np.ones((1, len(idx))))[0.25][0])
            d = np.mean(dd, 0)
            P(f"  [per-cell operating point, mean over 6 cells] H-BOMP @0.25 {np.mean(pp):+.4f} [{np.percentile(d, 2.5):+.4f}, {np.percentile(d, 97.5):+.4f}]")
        # ---- (C) mechanism: within-window separation vs cross-window score reliability
        kk = sum(w["k"] for w in ww)
        for n in (1, 2, 3, 5):
            row = []
            for m in (BO, NC, MT, H):
                hits = sum(len({d[1] for d in w["det"][m][:n] if d[1] >= 0}) for w in ww)
                fal = sum(sum(1 for d in w["det"][m][:n] if d[1] < 0) for w in ww)
                row.append(f"{m} recall {hits / kk:.3f} false/win {fal / len(ww):.2f}")
            P(f"  [threshold-free, first {n} atoms of every window] " + "; ".join(row))
        row = []
        for m in (BO, NC, MT, H):
            hits = sum(len({d[1] for d in w["det"][m] if d[1] >= 0}) for w in ww)
            row.append(f"{m} {hits / kk:.3f}")
        P("  [ceiling: recall with all atoms accepted] " + "; ".join(row))
        # first-atom discrimination across windows (AUC of first-atom score: hit vs false)
        row = []
        for m in (BO, NC, MT, H):
            s = np.array([w["det"][m][0][0] for w in ww]); h = np.array([w["det"][m][0][1] >= 0 for w in ww])
            r = np.argsort(np.argsort(s)) + 1.0
            auc = (r[h].sum() - h.sum() * (h.sum() + 1) / 2) / (h.sum() * (~h).sum())
            row.append(f"{m} P(first atom hit) {h.mean():.3f} AUC {auc:.3f}")
        P("  [first atom: hit fraction (within-window) and cross-window AUC of its score] " + "; ".join(row))
        # score-swap by rank: labels of X with scores of Y (both records are non-increasing in rank, so within-window order is kept)
        def swapped(src):
            return lambda w, m, a: w["det"][src][a][0] if a < len(w["det"][src]) else -np.inf
        one = np.ones((1, len(ww)))
        for labm in (BO, H, MT):
            row = []
            for scm in (BO, H, MT):
                f = Froc(ww, labm, score=swapped(scm)); r = f.recall(one)
                row.append(f"scores of {scm}: {r[0.1][0]:.3f}/{r[0.25][0]:.3f}")
            P(f"  [rank-matched score swap] atoms+labels of {labm:11s} recall 0.1/0.25 with " + "; ".join(row))
        fa = Froc(ww, BO, score=swapped(H)); fb = Froc(ww, BO); d = fa.recall(C)[0.25] - fb.recall(C)[0.25]
        P(f"  [swap] (BOMP atoms, HUBER scores) - BOMP @0.25: {fa.recall(one)[0.25][0] - fb.recall(one)[0.25][0]:+.4f} [{np.percentile(d, 2.5):+.4f}, {np.percentile(d, 97.5):+.4f}]")
        fa = Froc(ww, H); fb = Froc(ww, H, score=swapped(BO)); d = fa.recall(C)[0.25] - fb.recall(C)[0.25]
        P(f"  [swap] HUBER - (HUBER atoms, BOMP scores) @0.25: {fa.recall(one)[0.25][0] - fb.recall(one)[0.25][0]:+.4f} [{np.percentile(d, 2.5):+.4f}, {np.percentile(d, 97.5):+.4f}]")
        fa = Froc(ww, H, score=swapped(BO)); fb = Froc(ww, BO); d = fa.recall(C)[0.25] - fb.recall(C)[0.25]
        P(f"  [swap] (HUBER atoms, BOMP scores) - BOMP @0.25 (pure atom-selection effect): {fa.recall(one)[0.25][0] - fb.recall(one)[0.25][0]:+.4f} [{np.percentile(d, 2.5):+.4f}, {np.percentile(d, 97.5):+.4f}]")
        same = np.mean([w["det"][H][0][2:5] == w["det"][BO][0][2:5] for w in ww])
        P(f"  fraction of windows whose first atom is at the same node for HUBER-N and BOMP-PK: {same:.3f}")
