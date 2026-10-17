"""
diag_alignment.py - adversarial diagnostics (post hoc, descriptive):
  (1) origin-time and epicentre error of the detections nearest to each EVAL M>=2 event (detections within +-3 s;
      nearest epicentre), per subnet and method: a systematic origin offset would indicate a time-alignment bug.
  (2) PhaseNet pick residuals against the catalogue (homogeneous VP 6.0, VS 6.0/1.73, straight rays, station at z=0):
      median P and S residual per station on CAL and EVAL events (absolute time alignment of waveforms vs catalogue,
      and whether the +SC corrections estimated on CAL match the EVAL residuals).
    python diag_alignment.py M2DIR CLOSEDIR PICKSDIR
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import verif_analysis as V  # noqa: E402

SUBS = V.SUBS


def nearest_stats(sets, cat, targets):
    print("(1) nearest detection (within +-3 s) to each M>=2 EVAL event: median epicentral distance, "
          "median/mean origin error (det - cat) of those within 5 km, n events with a det within 5 km & 3 s")
    for (s, m), ss in sorted(sets.items()):
        best = None
        for S in ss:
            det = np.array([[d[0], d[1], d[3]] for d in S["det"]]) if S["n"] else np.zeros((0, 3))
            dist, dts = [], []
            for i in targets:
                if not len(det):
                    continue
                ok = np.abs(det[:, 2] - cat[i, 2]) <= 3.0
                if not ok.any():
                    dist.append(np.inf)
                    continue
                dd = np.hypot(det[ok, 0] - cat[i, 0], det[ok, 1] - cat[i, 1])
                j = int(dd.argmin())
                dist.append(dd[j])
                if dd[j] <= 5.0:
                    dts.append(det[ok][j, 2] - cat[i, 2])
            n5 = len(dts)
            if best is None or n5 > best[0]:
                best = (n5, np.median(dist) if dist else np.nan, np.median(dts) if dts else np.nan,
                        np.mean(dts) if dts else np.nan, np.mean(np.abs(dts) <= 0.5) if dts else np.nan)
        print(f"  {s} {m:9s} med dist {best[1]:6.2f} km  n(<=5km) {best[0]:2d}  med dt {best[2]:+.2f}  "
              f"mean dt {best[3]:+.2f}  frac |dt|<=0.5 {best[4]:.2f}")


def pick_residuals(pdir, cat):
    import r0_lib as R  # stations in the frame used by the run
    print("(2) PhaseNet pick residuals vs catalogue (pick - (t0 + d/V)), matched within +-1.5 s of prediction, "
          "per station: median P / S residual on CAL events | EVAL events, and the CAL +SC correction (P, S)")
    centre = R.frame()
    for s in SUBS:
        f = os.path.join(pdir, f"picks_eval_{s}.json")
        if not os.path.exists(f):
            continue
        r = json.load(open(f))
        STA, _ = R.stations(r["names"], centre)
        corr = np.array(r["corr_sc"])
        print(f"  {s}")
        for i, nm in enumerate(r["names"]):
            pk = r["picks"][i]
            out = []
            for seg in ((0, 900), (900, 3560)):
                res = {"P": [], "S": []}
                for e in cat:
                    if not (seg[0] <= e[2] < seg[1]) or e[3] < 2.0:
                        continue
                    d = math.sqrt((STA[i, 0] - e[0]) ** 2 + (STA[i, 1] - e[1]) ** 2 + e[4] ** 2)
                    for b, v in (("P", 6.0), ("S", 6.0 / 1.73)):
                        tp = e[2] + d / v
                        c = [p[0] - tp for p in pk if p[1] == b and abs(p[0] - tp) <= 1.5]
                        if c:
                            res[b].append(min(c, key=abs))
                out.append(res)
            f2 = lambda x: f"{np.median(x):+.2f}({len(x):2d})" if x else "   nan( 0)"  # noqa: E731
            print(f"    {nm:5s} CAL P {f2(out[0]['P'])} S {f2(out[0]['S'])} | EVAL P {f2(out[1]['P'])} "
                  f"S {f2(out[1]['S'])} | corr_sc P {corr[i, 0]:+.2f} S {corr[i, 1]:+.2f}")


if __name__ == "__main__":
    cat, centre = V.catalogue()
    targets = np.where((cat[:, 2] >= V.LO) & (cat[:, 2] < V.HI) & (cat[:, 3] >= 2.0))[0]
    sets = V.load(sys.argv[1], sys.argv[2], sys.argv[3], cat)
    nearest_stats(sets, cat, targets)
    sys.path.insert(0, "/home/claude/work/r0")
    pick_residuals(sys.argv[3], cat)
