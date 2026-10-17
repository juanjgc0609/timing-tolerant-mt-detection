"""C18a supplement: are the large timing residuals of CAL readings (HOM6 + leave-one-out station-phase constant)
explained by other catalogued events arriving near the same station within +-2 s (overlap contamination)?"""
import json, sys
import numpy as np
sys.path.insert(0, "/home/claude/work/r0")
import r0_lib as R
from task0_common import VP, VS
j = json.load(open(sys.argv[1]))
rows = j["rows"]
centre = R.frame(); cat = R.catalogue_local(centre)
STA, _ = R.stations(R.SUBNETS[j["subnet"]], centre)
evs = [e for e in cat if R.T_CAL[0] <= e["t"] < R.T_CAL[1] and e["M"] >= j["minmag"] and e["t"] >= 6.5]
def pred(e, i, b):
    d = np.linalg.norm(STA[i] - np.array([e["x"], e["y"], e["z"]])); return e["t"] + d / (VP if b == 0 else VS)
err, ovl = [], []
for r in rows:
    tr = [q for q in rows if q["ev"] != r["ev"] and q["i"] == r["i"] and q["b"] == r["b"]]
    if len(tr) < 3: continue
    e = evs[r["ev"]]
    res = lambda q: q["tobs"] - (np.linalg.norm(STA[q["i"]] - np.array(q["src"])) / (VP if q["b"] == 0 else VS))
    err.append(res(r) - np.median([res(q) for q in tr]))
    tabs = e["t"] - 6.5 + r["tobs"]          # absolute time of the reading (window starts at origin - 6.5 s)
    others = [abs(pred(o, r["i"], bb) - tabs) for o in cat if o is not e for bb in (0, 1)]
    ovl.append(min(others))
err, ovl = np.abs(np.array(err)), np.array(ovl)
for lab, m in (("another arrival within 2 s", ovl <= 2.0), ("none within 2 s", ovl > 2.0)):
    print(f"{j['subnet']} M>={j['minmag']}: {lab}: n {m.sum()}, fraction |e| <= 0.25: {np.mean(err[m] <= 0.25):.2f}, median |e| {np.median(err[m]):.3f}")
# breakdown by phase and station, and residual sign
B = np.array([r["b"] for r in rows if len([q for q in rows if q["ev"] != r["ev"] and q["i"] == r["i"] and q["b"] == r["b"]]) >= 3])
I = np.array([r["i"] for r in rows if len([q for q in rows if q["ev"] != r["ev"] and q["i"] == r["i"] and q["b"] == r["b"]]) >= 3])
for b in (0, 1):
    print(f"  phase {'PS'[b]}: n {np.sum(B == b)}, |e| <= 0.25 {np.mean(err[B == b] <= 0.25):.2f}, |e| > 0.5 {np.mean(err[B == b] > 0.5):.2f}")
print("  per station |e|<=0.25: " + " ".join(f"{R.SUBNETS[j['subnet']][i]}:{np.mean(err[I == i] <= 0.25):.2f}({np.sum(I == i)})" for i in np.unique(I)))
