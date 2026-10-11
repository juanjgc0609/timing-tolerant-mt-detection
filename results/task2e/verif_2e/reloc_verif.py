"""Decomposition of TS-SFW+SC - BOMP-PK+SC from stored parts (no recomputation): TS-SFW+SC with hits evaluated at the
ENTRY grid node of each atom (oracle node) instead of its final slid position. Difference TS(final) - TS(entry) is the
contribution of the off-grid relocation; TS(entry) - BOMP-PK is the contribution of the different residual/deflation."""
import json, os, glob, sys
import numpy as np
sys.path.insert(0, "/home/claude/work")
from task0_common import default_grid
import analyze_verif as A
garr = np.asarray(default_grid())
TOL = 5.0

def hit_of(src, truth):
    if not truth:
        return -1
    dd = [float(np.linalg.norm(np.asarray(src)[:2] - np.asarray(t)[:2])) for t in truth]
    j = int(np.argmin(dd))
    return j if dd[j] <= TOL else -1

out = []
for tag, pat, key in (("seed 999 (calibration, 600 windows)", "cal_seed999_{r}_*.json", "TS-SFW+SC f=0.01"),
                      ("seeds 2024+7331 (evaluation, 2400 windows per regime)", "ev_seed*_{r}_*.json", "TS-SFW+SC")):
    for reg in ("thin", "strong"):
        files = sorted(glob.glob(os.path.join("/home/claude/work/task2e_parts", pat.format(r=reg))))
        per_cell = {}
        for f in files:
            p = json.load(open(f))
            c = (p["cell"], p.get("seed"))
            for w in p["windows"]:
                d = per_cell.setdefault(c, {"k": [], "fin": [], "ent": [], "bp": []})
                it = w["diag"][key]["iters"]
                rec = w["records"][key]
                assert len(it) == len(rec)
                d["k"].append(w["k"]); d["fin"].append(rec); d["bp"].append(w["records"]["BOMP-PK+SC"])
                d["ent"].append([(s, hit_of(garr[x[1]], w["truth"])) for (s, h), x in zip(rec, it)])
        rows = []
        for c, d in sorted(per_cell.items()):
            k = np.array(d["k"], float); one = np.ones((1, len(k)))
            r = {m: A.recall(A.HFmat(d[m]), one, k) for m in ("fin", "ent", "bp")}
            rows.append((c, {m: {cap: float(r[m][cap][0]) for cap in A.CAPS} for m in r}))
        fin = np.mean([[x[1]["fin"][c] - x[1]["ent"][c] for c in A.CAPS] for x in rows], 0)
        ent = np.mean([[x[1]["ent"][c] - x[1]["bp"][c] for c in A.CAPS] for x in rows], 0)
        tot = np.mean([[x[1]["fin"][c] - x[1]["bp"][c] for c in A.CAPS] for x in rows], 0)
        s = (f"{tag}, {reg}: mean over {len(rows)} cell(-seed) pools of recall differences at 0.1 / 0.25: "
             f"TS-SFW - BOMP-PK {tot[0]:+.3f} / {tot[1]:+.3f} = relocation (final vs entry node) {fin[0]:+.3f} / "
             f"{fin[1]:+.3f} + residual/ordering (TS-SFW at entry nodes - BOMP-PK) {ent[0]:+.3f} / {ent[1]:+.3f}")
        print(s); out.append(s)
open("/home/claude/work/verif_2e/reloc_verif.txt", "w").write("\n".join(out) + "\n")
