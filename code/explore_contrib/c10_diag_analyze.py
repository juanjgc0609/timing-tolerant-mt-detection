"""C10 first-iteration separability diagnostic -> out/c10_diag.txt"""
import json

import numpy as np

RULES = ("hard", "T2w-cat", "T2w-ora", "fine-ora")
lines = []


def P(s=""):
    print(s)
    lines.append(s)


P("C10 DIAGNOSTIC: first-iteration separability (no deflation).  Per true event: best node statistic within 5 km versus")
P("the best statistic of the window among nodes farther than 5 km from every true event.  [95%: bootstrap over networks]")
rng = np.random.default_rng(99)
for design in ("STD",):                 # the DENSE version was withdrawn before it ran
    rows = []
    for reg in ("thin", "strong"):
        for r in json.load(open(f"out/c10d_{design}_{reg}.json")):
            r["regime"] = reg
            rows.append(r)
    g = {}
    for i, r in enumerate(rows):
        g.setdefault((r["regime"], r["cell"], r["net"]), []).append(i)
    gl = list(g.values())
    draws = [rng.integers(0, len(gl), len(gl)) for _ in range(400)]
    P(f"== {design}: {len(rows)} windows with events, {sum(r['k'] for r in rows)} events, {len(gl)} network x regime")

    def stat(nm, ws, sub=None):
        m = np.array([t - rows[w]["ev"][nm]["false"] for w in ws for t in rows[w]["ev"][nm]["true"]
                      if sub is None or rows[w]["regime"] == sub])
        return np.mean(m > 0), np.median(m)
    allw = list(range(len(rows)))
    for sub in (None, "thin", "strong"):
        P(f"   regimes {sub or 'both'}: share of events whose best value beats the best false alternative; median margin")
        base = stat("hard", allw, sub)
        for nm in RULES:
            pt = stat(nm, allw, sub)
            bs = np.array([np.r_[stat(nm, [w for q in d for w in gl[q]], sub),
                                 stat(nm, [w for q in d for w in gl[q]], sub)[0]
                                 - stat("hard", [w for q in d for w in gl[q]], sub)[0]] for d in draws])
            P(f"      {nm:9s} share {pt[0]:.3f} [{np.percentile(bs[:, 0], 2.5):.3f},{np.percentile(bs[:, 0], 97.5):.3f}]  "
              f"- hard {pt[0] - base[0]:+.3f} [{np.percentile(bs[:, 2], 2.5):+.3f},{np.percentile(bs[:, 2], 97.5):+.3f}]  "
              f"margin {pt[1]:+.2f} [{np.percentile(bs[:, 1], 2.5):+.2f},{np.percentile(bs[:, 1], 97.5):+.2f}]")
    e = np.array([x for r in rows for x in r["ev"]["time-ora"]])
    P(f"   maximum of the oracle joint objective (1 km lattice) within 10 km of each true event: median distance to the "
      f"event {np.median(e):.2f} km; within 1 / 2 / 5 km: {np.mean(e <= 1):.2f} / {np.mean(e <= 2):.2f} / "
      f"{np.mean(e <= 5):.2f}")
    # the same split by whether the event wins under the hard rule
    win = np.array([t > r["ev"]["hard"]["false"] for r in rows for t in r["ev"]["hard"]["true"]])
    P(f"      among events that beat the best false alternative under the hard rule ({win.mean():.2f}): median "
      f"{np.median(e[win]) if win.any() else np.nan:.2f} km, within 2 km {np.mean(e[win] <= 2) if win.any() else np.nan:.2f}; "
      f"among the others: median {np.median(e[~win]):.2f} km, within 2 km {np.mean(e[~win] <= 2):.2f}")
open("out/c10_diag.txt", "w").write("\n".join(lines) + "\n")
