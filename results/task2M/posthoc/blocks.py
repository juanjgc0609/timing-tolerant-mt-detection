"""Task 2M, POST HOC exploratory breakdown (written after reading analyze_2m_out.txt): the REAL contrasts of 2M-a/2M-b
(DENSE, LOC, @0.25, cells pooled) by groups of networks, i.e. by blocks of hours of the noise day, with the noise
statistics of their windows.  Same loader, records and bootstrap as t2m_analyze.py (c16_analyze / c10_analyze.Boot)."""
import os, sys
import numpy as np
W = "/home/claude/work"
sys.path[:0] = [os.path.join(W, "explore_contrib"), W]
import c10_analyze as A
import c16_analyze as Z
B = 1000
for lab, pdir, groups in (("2718", "t2m_out/parts_2718", [(2800, 2807), (2808, 2815), (2816, 2823)]),
                          ("3141", "t2m_out/parts_3141", [(2900, 2905), (2906, 2909), (2910, 2916)])):
    nets, rows = Z.load(os.path.join(W, pdir), "REAL")
    first_slot = {n["net"]: min(min(v) for v in n["slots"].values() if v) for n in nets}
    sub_all = [r for r in rows if r["design"] == "DENSE"]
    names = list(sub_all[0]["atoms"])
    for g0, g1 in groups:
        sub = [r for r in sub_all if g0 <= r["net"] <= g1]
        recs = Z.records(sub, names, None)
        bt = A.Boot(sub, recs, lambda r: True, B)
        hrs = sorted({first_slot[j] * 40 / 3600 for j in range(g0, g1 + 1) if j in first_slot})
        ku = np.array([r["kurt"] for r in sub]); mx = np.array([r["maxabs"] for r in sub])
        out = []
        for a, b in (("MT@S0", "BOMP@S0"), ("HUBER-N@S0", "BOMP@S0")):
            d = bt.rep[a][:, 1] - bt.rep[b][:, 1]
            out.append(f"{a}-{b} {bt.pt[a][1] - bt.pt[b][1]:+.3f} [{np.percentile(d, 2.5):+.3f},{np.percentile(d, 97.5):+.3f}]")
        print(f"{lab} nets {g0}-{g1} (first slot hour {hrs[0]:.1f}-{hrs[-1]:.1f}): {len(sub)} windows; "
              f"BOMP {bt.pt['BOMP@S0'][1]:.3f} MT {bt.pt['MT@S0'][1]:.3f} HUBER {bt.pt['HUBER-N@S0'][1]:.3f}; " + "; ".join(out)
              + f"; window kurtosis median {np.median(ku):.2f} p90 {np.percentile(ku, 90):.1f}, max|MF|/sd median {np.median(mx):.1f} "
              f"frac>10 {np.mean(mx > 10):.2f}", flush=True)
