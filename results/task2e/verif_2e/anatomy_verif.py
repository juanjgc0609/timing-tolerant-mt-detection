"""Where the detectable-but-missed events go (seed 999 calibration windows): distance from each true event to the
nearest atom of BOMP-PK+SC / TS-SFW+SC (final positions), for events never hit although their best near-truth grid
score beats every clutter+noise-only candidate of the window."""
import glob, json, os
import numpy as np
out = []
C = {}
for f in sorted(glob.glob("/home/claude/work/verif_2e/check4_parts_verif/c_*.json")):
    p = json.load(open(f))
    st = json.load(open(f"/home/claude/work/task2e_parts/cal_seed999_{p['regime']}_{p['cell']}_net{p['net']}.json"))
    for w, ws in zip(p["windows"], st["windows"]):
        nmax = max(w["g_noise"]) if w["g_noise"] else 0.0
        for m, key in (("BOMP-PK+SC", "BOMP-PK+SC"), ("TS-SFW+SC", "TS-SFW+SC f=0.01")):
            hits = {h for s, h in ws["records"][key] if h >= 0}
            locs = np.array(ws["loc"][key])
            recs = ws["records"][key]
            for i, g in enumerate(w["g_true"]):
                if g <= nmax or i in hits:
                    continue
                src = np.array(ws["truth"][i])
                d = np.linalg.norm(locs[:, :2] - src[None, :2], axis=1)
                j = int(np.argmin(d))
                # score of the nearest atom vs the near-truth grid score
                C.setdefault(m, []).append((float(d[j]), float(recs[j][0]) >= g, recs[j][1], p["regime"]))
for m, v in C.items():
    d = np.array([x[0] for x in v])
    s = (f"{m}: {len(v)} detectable-but-never-hit events; nearest atom epicentral distance median {np.median(d):.1f} km; "
         f"in (5, 10] km: {np.mean((d > 5) & (d <= 10)):.2f}; in (10, 15] km: {np.mean((d > 10) & (d <= 15)):.2f}; "
         f"> 15 km: {np.mean(d > 15):.2f}; nearest atom is a false atom whose score >= the event's near-truth grid "
         f"score: {np.mean([x[1] and x[2] < 0 for x in v]):.2f}; nearest atom duplicates another event: "
         f"{np.mean([x[2] >= 0 for x in v]):.2f}")
    print(s); out.append(s)
open("/home/claude/work/verif_2e/anatomy_verif.txt", "w").write("\n".join(out) + "\n")
