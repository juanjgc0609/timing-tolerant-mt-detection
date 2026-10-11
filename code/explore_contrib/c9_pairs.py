"""C9 P4 supplementary pairs (cross-base differences and per-cell values for the polarization-weighted variants).
    python c9_pairs.py > out/c9_pairs.txt"""
import os

import numpy as np

import c9_analyze as A

A.B = 1000
KEY = ["BOMP-PK", "MT-FW", "BOMP-PK/RR-pol", "BOMP-PK/RR-pol^2", "MT-FW/RR-pol", "MT-FW/RR-pol^2", "MT-FW+SELPOL",
       "MT-FW+PENPOL", "MT-FW+PEN", "CEIL-GD"]
PAIRS = [("MT-FW+SELPOL", "MT-FW"), ("MT-FW+SELPOL", "BOMP-PK"), ("MT-FW+PENPOL", "MT-FW"), ("MT-FW+PENPOL", "BOMP-PK"),
         ("MT-FW/RR-pol", "MT-FW"), ("MT-FW/RR-pol", "BOMP-PK"), ("BOMP-PK/RR-pol", "BOMP-PK"),
         ("BOMP-PK/RR-pol^2", "BOMP-PK"), ("MT-FW", "BOMP-PK")]


def rep(rows, title):
    built = [A.build(r) for r in rows]
    recs = [b[0] for b in built]
    meths = [m for m in KEY if all(m in rc for rc in recs)]
    pt, reps = A.boot(rows, recs, meths)
    print(f"\n{title}: windows {len(rows)}, events {sum(r['k'] for r in rows)}")
    for a, b in PAIRS:
        if a not in meths or b not in meths:
            continue
        d = reps[a] - reps[b]
        c = np.percentile(d, [2.5, 97.5], axis=0)
        p = np.subtract(pt[a], pt[b])
        line = (f"  {a:17s} - {b:8s} @0.1 {A.f3(p[0], c[0][0], c[1][0])}  @0.25 {A.f3(p[1], c[0][1], c[1][1])}  "
                f"P(<=0) {np.mean(d[:, 1] <= 0):.3f}  SD {d[:, 1].std():.4f}")
        if "CEIL-GD" in meths and b == "BOMP-PK":
            H = reps["CEIL-GD"][:, 1] - reps["BOMP-PK"][:, 1]
            fr = d[:, 1] / np.where(np.abs(H) > 1e-9, H, np.nan)
            fr = fr[np.isfinite(fr)]
            line += (f"  headroom {p[1] / (pt['CEIL-GD'][1] - pt['BOMP-PK'][1]):+.2f} [{np.percentile(fr, 2.5):+.2f}, "
                     f"{np.percentile(fr, 97.5):+.2f}]")
        print(line)


for design in ("DENSE", "STD", "DENSE-POL45"):
    rows = A.load(design)
    rep(rows, f"{design} pooled")
    if design != "DENSE-POL45":
        for cell in ("ladder_4", "ladder_6", "unilateral_c40"):
            rep([r for r in rows if r["cell"] == cell], f"{design} {cell} (both regimes)")
