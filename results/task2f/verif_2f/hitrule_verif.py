"""hitrule_verif.py - exploratory: does the 2F gain survive stricter hit rules?  (verification; not registered)
The registered hit rule (2A-2E) is epicentre within 5 km of the nearest true event, origin time ignored.  In dense
windows all 3-5 events lie within 9 km of a centre, so location-only matching is loose.  Here hits are recomputed from
the stored atom locations and origins (MT-FW+SC, MT-FW-NOCAP+SC, BOMP-PK+SC store them; BOMP-PK-CAP+SC does not):
  R5      registered rule recomputed (must equal the stored labels)
  R5t1    within 5 km AND |origin - true origin| <= 1.0 s (nearest such event)
  R5t05   within 5 km AND |origin error| <= 0.5 s
  R3t05   within 3 km AND |origin error| <= 0.5 s
Scores (running minimum) unchanged; exact pooled FROC from stats_verif.froc_point.
"""
import numpy as np

from stats_verif import load, froc_point

METH = ("MT-FW+SC", "MT-FW-NOCAP+SC", "BOMP-PK+SC")
RULES = {"R5": (5.0, None), "R5t1": (5.0, 1.0), "R5t05": (5.0, 0.5), "R3t05": (3.0, 0.5)}


def relabel(W, raw):
    out, mism = {r: [] for r in RULES}, 0
    for w, rw in zip(W, raw):
        tr = np.array(rw["truth"], float)
        taus = np.array(rw["taus"], float)
        for rname, (dk, dt) in RULES.items():
            rec = {}
            for m in METH:
                lab = []
                for (s, h0), loc in zip(w["rec"][m], rw["loc"][m]):
                    d = np.linalg.norm(tr[:, :2] - np.array(loc[:2]), axis=1)
                    ok = d <= dk
                    if dt is not None:
                        ok &= np.abs(taus - loc[3]) <= dt
                    h = int(np.argmin(np.where(ok, d, np.inf))) if ok.any() else -1
                    if rname == "R5" and h != h0:
                        mism += 1
                    lab.append((s, h))
                rec[m] = lab
            out[rname].append({"k": w["k"], "rec": rec})
    return out, mism


if __name__ == "__main__":
    import json
    import os
    from stats_verif import PARTS, CELLS
    for seed in (4242, 8888):
        for reg in ("thin", "strong"):
            W = load(seed, "dense", reg, range(96))
            raw = []
            for cell in CELLS:
                for j in range(96):
                    e = json.load(open(os.path.join(PARTS, f"dense_seed{seed}_{reg}_{cell}_net{j}.json")))
                    raw.extend(e["windows"])
            rel, mism = relabel(W, raw)
            print(f"seed {seed} {reg}: R5 relabel mismatches vs stored labels: {mism}")
            for rname in RULES:
                p = {m: froc_point(rel[rname], m) for m in METH}
                print(f"   {rname:6s} " + "  ".join(f"{m} {p[m][0.1]:.4f}/{p[m][0.25]:.4f}" for m in METH) +
                      f"  | MT-FW - BOMP-PK {p[METH[0]][0.1] - p[METH[2]][0.1]:+.4f}/{p[METH[0]][0.25] - p[METH[2]][0.25]:+.4f}"
                      f"  MT-FW - NOCAP {p[METH[0]][0.25] - p[METH[1]][0.25]:+.4f}")
