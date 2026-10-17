"""
diag_ceiling.py - robustness of the primary to ranking and de-duplication (post hoc, descriptive):
  (a) ranking ceiling: fraction of the 45 targets with ANY LOC+T (and LOC) hit at any score, per subnet, for M2 with
      and without de-duplication, vs. the BESTPICK recall at 30/h; pooled bound on M2 - BESTPICK under perfect ranking.
  (b) M2 recall at 30/h without de-duplication.
  (c) primary contrast with LOC+T time tolerance widened to 1.0 and 1.5 s and epicentral radius 10 km (sensitivity).
    python diag_ceiling.py M2DIR CLOSEDIR PICKSDIR
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import verif_analysis as V  # noqa: E402

cat, _ = V.catalogue()
tg = np.where((cat[:, 2] >= V.LO) & (cat[:, 2] < V.HI) & (cat[:, 3] >= 2.0))[0]
hours = (V.HI - V.LO) / 3600.0
cap = int(math.floor(30 * hours + 1e-9))
m2dir, cdir, pdir = sys.argv[1:4]


def lab_custom(det, tol, rkm):
    out = np.full(len(det), -1, int)
    for i, d in enumerate(det):
        dist = np.hypot(cat[:, 0] - d[0], cat[:, 1] - d[1])
        ok = (dist <= rkm) & (np.abs(cat[:, 2] - d[3]) <= tol + 1e-9)
        big = ok & (cat[:, 3] >= 2.0)
        if big.any():
            out[i] = int(np.where(big, dist, np.inf).argmin())
        elif ok.any():
            out[i] = -2
    return out


def rec(det, tol, rkm):
    det = [d for d in det if V.LO <= d[3] < V.HI]
    if not det:
        return 0.0, 0.0
    s = np.array([d[4] for d in det])
    lab = lab_custom(det, tol, rkm)
    anyh = len(set(lab[lab >= 0].tolist()) & set(tg.tolist())) / len(tg)
    return V.recall_at(s, lab, tg, cap), anyh


rows = {}
for sub in V.SUBS:
    m2 = json.load(open(os.path.join(m2dir, f"det_eval_{sub}.json")))["atoms"]["M2"]
    pk = json.load(open(os.path.join(pdir, f"picks_eval_{sub}.json")))
    assoc = [v for key in ("gamma", "pyocto", "gamma_sc", "pyocto_sc") for v in pk[key].values()]
    r = {}
    for name, tol, rkm in (("LOC+T", 0.5, 5.0), ("LOC", 3.0, 5.0), ("T1.0", 1.0, 5.0), ("T1.5", 1.5, 5.0),
                           ("R10T0.5", 0.5, 10.0), ("R10T1.5", 1.5, 10.0)):
        a_dd = rec(V.dedupe(m2), tol, rkm)
        a_raw = rec(m2, tol, rkm)
        bp = max(rec(v, tol, rkm)[0] for v in assoc)
        bp_any = max(rec(v, tol, rkm)[1] for v in assoc)
        r[name] = dict(m2_30=a_dd[0], m2_any=a_dd[1], m2_nodedup_30=a_raw[0], m2_nodedup_any=a_raw[1], bp_30=bp,
                       bp_any=bp_any)
    rows[sub] = r
print("per subnet: M2 recall@30/h | M2 any-score ceiling | M2 no-dedup @30/h, any | BESTPICK @30/h (any-score)")
for name in rows["N8"]:
    print(f" rule {name}")
    for sub in V.SUBS:
        x = rows[sub][name]
        print(f"   {sub} M2 {x['m2_30']:.3f} | ceil {x['m2_any']:.3f} | nodedup {x['m2_nodedup_30']:.3f}, "
              f"{x['m2_nodedup_any']:.3f} | BESTPICK {x['bp_30']:.3f} ({x['bp_any']:.3f})")
    p = np.mean([rows[s][name]["m2_30"] - rows[s][name]["bp_30"] for s in V.SUBS])
    c = np.mean([rows[s][name]["m2_any"] - rows[s][name]["bp_30"] for s in V.SUBS])
    c2 = np.mean([max(rows[s][name]["m2_any"], rows[s][name]["m2_nodedup_any"]) - rows[s][name]["bp_30"]
                  for s in V.SUBS])
    print(f"   pooled M2 - BESTPICK @30/h {p:+.3f}; perfect-ranking bound (M2 any-score ceiling - BESTPICK@30/h) "
          f"{c:+.3f}; same without de-dup {c2:+.3f}")
