import sys, os, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats_verif import *
P = print
for seed in (5151, 6262):
    wins, miss, meta = load(seed)
    P(f"=== seed {seed}: windows {len(wins)}, missing parts {len(miss)}, meta variants {len(meta)}: {sorted(meta)}")
    cnt = collections.Counter((w["kind"], w["reg"], w["cell"]) for w in wins)
    nets = collections.defaultdict(set)
    for w in wins: nets[(w["kind"], w["reg"], w["cell"])].add(w["net"])
    P("  nets per (kind, regime, cell):", sorted(set(len(v) for v in nets.values())), "windows per cell:", sorted(set(cnt.values())),
      "n cells:", len(cnt))
    for kind in ("dense", "std"):
        ww = [w for w in wins if w["kind"] == kind]
        P(f"  {kind}: {len(ww)} windows, {sum(w['k'] for w in ww)} events, k set {sorted(set(w['k'] for w in ww))}; mean detections/window",
          {m: round(np.mean([len(w['det'][m]) for w in ww]), 2) for m in METH})
        # stored labels vs my recomputation of the registered rule
        nd = sum(int(h != d[1]) for w in ww for m in METH for h, d in zip(relabel(w, m), w["det"][m]))
        mono = all(all(w["det"][m][i][0] >= w["det"][m][i + 1][0] for i in range(len(w["det"][m]) - 1)) for w in ww for m in METH)
        sep = min([np.hypot(a[2] - b[2], a[3] - b[3]) for w in ww for m in METH for i, a in enumerate(w["det"][m]) for b in w["det"][m][:i]] or [9e9])
        P(f"    stored hit labels differing from own recomputation (<=5 km epicentral, nearest): {nd}; scores non-increasing in every record: {mono}; min epicentral separation between atoms of a record {sep:.2f} km")
        P("    mean seconds per window (the 'sec' field):", {m: round(np.mean([w['sec'][m] for w in ww]), 3) for m in METH})
        rng = np.random.default_rng([20261002, seed, 0 if kind == "dense" else 1])
        C = boot_counts(ww, rng, B)
        for m in METH:
            f = Froc(ww, m); r = f.recall(np.ones((1, len(ww))))
            P(f"    {m:15s} recall cap0.1 {r[0.1][0]:.4f} (brute {f.brute(0.1):.4f})  cap0.25 {r[0.25][0]:.4f} (brute {f.brute(0.25):.4f})")
        for a, b in ((H, BO), (H, MT), (MT, BO)):
            c = contrast(ww, a, b, C)
            P(f"    {kind} pooled {a} - {b} @0.25: {fmt(c)}   | @0.1: {c[0.1]['pt']:+.4f} [{c[0.1]['lo']:+.4f}, {c[0.1]['hi']:+.4f}]  (B={B})")
