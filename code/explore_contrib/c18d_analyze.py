"""c18d_analyze.py - recall on CAL at 10 / 30 / 60 false detections per hour (r0_analyze.py functions: de-duplication,
labels, exact FROC; 60 s block bootstrap, B = 1000) for the C18d variants.
    python c18d_analyze.py OUT.txt FILE.json [FILE.json ...]"""
import json
import sys

import numpy as np

sys.path.insert(0, "/home/claude/work/r0")
import r0_analyze as A  # noqa: E402
import r0_lib as R  # noqa: E402


def main(outp, files):
    lo, hi = R.T_CAL
    cat = R.catalogue_local(R.frame())
    M = np.array([e["M"] for e in cat])
    tg = np.array([i for i, e in enumerate(cat) if lo <= e["t"] < hi and e["M"] >= R.M_TARGET])
    tblk = np.array([int((cat[i]["t"] - lo) // A.BLOCK) for i in tg])
    nb = int(np.ceil((hi - lo) / A.BLOCK))
    lines = [f"C18d: CAL [{lo}, {hi}) s, {len(tg)} target events (M >= 2.0), {nb} blocks"]
    rng = np.random.default_rng(18)
    W = [np.bincount(rng.integers(0, nb, nb), minlength=nb).astype(float) for _ in range(1000)]
    for f in files:
        r = json.load(open(f))
        sets = {v: [A.Set(A.dedupe(a), cat, lo, hi)] for v, a in r["atoms"].items()}
        for rule in ("LOC+T", "LOC"):
            for rate in (10, 30, 60):
                pt = {v: A.froc(s, (tg, tblk), M, nb, np.ones(nb), rate, rule) for v, s in sets.items()}
                row = f"{r['subnet']} {rule:5s} {rate:2d}/h  " + "  ".join(f"{v} {pt[v]:.3f}" for v in sets)
                if rate == 30 and "CAP6" in sets:
                    reps = [{v: A.froc(s, (tg, tblk), M, nb, wb, rate, rule) for v, s in sets.items()} for wb in W]
                    row += "   vs CAP6: " + "  ".join(
                        f"{v} {pt[v] - pt['CAP6']:+.3f} [{np.percentile([x[v] - x['CAP6'] for x in reps], 2.5):+.3f},"
                        f"{np.percentile([x[v] - x['CAP6'] for x in reps], 97.5):+.3f}]" for v in sets if v != "CAP6")
                lines.append(row)
        nh = {v: int(len(s[0].s)) for v, s in sets.items()}
        lines.append(f"{r['subnet']} detections after de-duplication: {nh}; seconds {r['sec']}")
    open(outp, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
