"""diag_corrections.py - CAL-only diagnostic of the +SC station corrections (post hoc): for a subnet, per station and
phase, the registered corr_sc (search +-1.0 s) and nacc, the same estimator with a +-2.0 s search, and the
fraction of accepted CAL delays that sit at the edge of the +-1 s search (|delay| >= 0.95 s)."""
import sys
import numpy as np
sys.path.insert(0, "/home/claude/work/r0")
import os
os.chdir("/home/claude/work/r0")
import r0_lib as R
import r0_run as RR
from station_corrections import estimate_station_corrections, peak_delays
s = sys.argv[1]
net = RR.Net(s, 940.0)
c1, n1 = net.corr_sc, net.nacc
c2, n2 = estimate_station_corrections(net.STA, net.calwin, search=2.0)
edge = np.zeros_like(c1)
for D, src, tau in net.calwin:
    dl, en = peak_delays(net.STA, D, src, tau)
    edge += (en >= 25.0) & (np.abs(dl) >= 0.95)
print(f"{s}: {len(net.calwin)} CAL windows; sigma_SC {net.eng['SC'][4][0]:.3f} s, sigma_S2 {net.eng['S2'][4][0]:.3f} s")
for i, nm in enumerate(net.names):
    print(f"  {nm:5s} P corr {c1[i,0]:+.2f} (n {n1[i,0]:2d}, edge {int(edge[i,0]):2d}) wide {c2[i,0]:+.2f} (n {n2[i,0]:2d})"
          f" | S corr {c1[i,1]:+.2f} (n {n1[i,1]:2d}, edge {int(edge[i,1]):2d}) wide {c2[i,1]:+.2f} (n {n2[i,1]:2d})"
          f" | S2 {net.corr_s2[i,0]:+.2f} {net.corr_s2[i,1]:+.2f}")
