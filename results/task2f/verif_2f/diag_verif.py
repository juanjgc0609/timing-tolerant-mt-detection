"""diag_verif.py - exploratory diagnostics on the confirmatory 2F parts (verification; not registered).
(a) 2F-b point estimates pooled over regimes; (b) kmax truncation (records are sequential, so the kmax = 5 / 3 record is
the prefix of the kmax = 7 record); (c) MT-FW - BOMP-PK by true event count k; (d) per-window score scale."""
import numpy as np

from stats_verif import load, froc_point

M1, M0, NC = "MT-FW+SC", "BOMP-PK+SC", "MT-FW-NOCAP+SC"


def trunc(W, n):
    return [{**w, "rec": {m: r[:n] for m, r in w["rec"].items()}} for w in W]


for seed in (4242, 8888):
    Wt = load(seed, "dense", "thin", range(96))
    Ws = load(seed, "dense", "strong", range(96))
    Wb = Wt + Ws
    p = {m: froc_point(Wb, m) for m in (M1, M0, NC)}
    print(f"seed {seed} pooled regimes: MT-FW - NOCAP {p[M1][0.25] - p[NC][0.25]:+.4f}; NOCAP - BOMP-PK "
          f"{p[NC][0.25] - p[M0][0.25]:+.4f}; MT-FW - BOMP-PK {p[M1][0.25] - p[M0][0.25]:+.4f} (cap 0.25)")
    for reg, W in (("thin", Wt), ("strong", Ws)):
        for n in (7, 5, 3, 2, 1):
            T_ = trunc(W, n)
            a, b = froc_point(T_, M1), froc_point(T_, M0)
            print(f"   {reg:6s} kmax {n}: MT-FW {a[0.1]:.4f}/{a[0.25]:.4f}  BOMP-PK {b[0.1]:.4f}/{b[0.25]:.4f}  diff "
                  f"{a[0.1] - b[0.1]:+.4f}/{a[0.25] - b[0.25]:+.4f}")
        for k in (3, 4, 5):
            Wk = [w for w in W if w["k"] == k]
            a, b = froc_point(Wk, M1), froc_point(Wk, M0)
            print(f"   {reg:6s} k={k} ({len(Wk)} windows): diff {a[0.1] - b[0.1]:+.4f}/{a[0.25] - b[0.25]:+.4f}")
        for m in (M1, M0):
            s1 = np.array([w["rec"][m][0][0] for w in W])
            fh = np.array([w["rec"][m][0][1] >= 0 for w in W])
            print(f"   {reg:6s} {m:12s} first-atom score: median {np.median(s1):.1f}, p90/p10 ratio "
                  f"{np.percentile(s1, 90) / np.percentile(s1, 10):.2f}; first-atom hit rate {fh.mean():.3f}; "
                  f"first-atom hit rate among top-20% scores {fh[s1 >= np.percentile(s1, 80)].mean():.3f}")
