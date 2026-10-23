"""c18j_eval_analyze.py - declared post hoc (E59): recall of the C18j configuration on the R0 evaluation segment, with
the functions of c18d_analyze (r0_analyze de-duplication, labels and exact FROC) and the same target definition as R0
(M >= 2.0). Events recalled are given as counts, with the C18j gate (5 of 16 CAL events, 0.31) and the registered R0
readings of M1 (homogeneous MT-FW) and BESTPICK on the same segment for reference; 95 per cent intervals of the 60 s
block bootstrap of R0 (B = 1000) are given for orientation only (see E61 on what this resampler can and cannot show).
    python c18j_eval_analyze.py OUT.txt FILE.json [FILE.json ...]"""
import json
import sys

import numpy as np

sys.path.insert(0, "/home/claude/work/r0")
import r0_analyze as A  # noqa: E402
import r0_lib as R  # noqa: E402

R0_EVAL = {  # results/taskR0/analysis_eval.txt, LOC+T, columns N8 and U8
    "N8": {10: (0.022, 0.400), 30: (0.044, 0.600), 60: (0.044, 0.644)},
    "U8": {10: (0.000, 0.000), 30: (0.000, 0.289), 60: (0.000, 0.333)},
}


def main(outp, files):
    lo, hi = R.T_EVAL
    cat = R.catalogue_local(R.frame())
    M = np.array([e["M"] for e in cat])
    tg = np.array([i for i, e in enumerate(cat) if lo <= e["t"] < hi and e["M"] >= R.M_TARGET])
    tblk = np.array([int((cat[i]["t"] - lo) // A.BLOCK) for i in tg])
    nb = int(np.ceil((hi - lo) / A.BLOCK))
    L = [f"C18j out of sample: EVAL [{lo}, {hi}) s, {len(tg)} target events (M >= 2.0), {nb} blocks; station terms "
         "estimated on the 16 CAL events only", ""]
    rng = np.random.default_rng(18)
    W = [np.bincount(rng.integers(0, nb, nb), minlength=nb).astype(float) for _ in range(1000)]
    for f in files:
        r = json.load(open(f))
        sub = r["subnet"]
        sets = {v: [A.Set(A.dedupe(a), cat, lo, hi)] for v, a in r["atoms"].items()}
        for rule in ("LOC+T", "LOC"):
            for rate in (10, 30, 60):
                cells = []
                for v, s in sets.items():
                    pt = A.froc(s, (tg, tblk), M, nb, np.ones(nb), rate, rule)
                    bs = [A.froc(s, (tg, tblk), M, nb, wb, rate, rule) for wb in W] if rate == 30 else None
                    c = f"{v} {pt:.3f} ({int(round(pt * len(tg)))}/{len(tg)})"
                    if bs is not None:
                        c += f" [{np.percentile(bs, 2.5):.3f}, {np.percentile(bs, 97.5):.3f}]"
                    cells.append(c)
                row = f"{sub} {rule:5s} {rate:2d}/h  " + "  ".join(cells)
                if rule == "LOC+T" and sub in R0_EVAL:
                    m1, bp = R0_EVAL[sub][rate]
                    row += f"   | R0 registered: M1 {m1:.3f}, BESTPICK {bp:.3f}"
                L.append(row)
        nh = {v: int(len(s[0].s)) for v, s in sets.items()}
        L.append(f"{sub} detections after de-duplication: {nh}; corrections P {np.round(np.asarray(r['corr'])[:, 0], 2).tolist()}"
                 f" S {np.round(np.asarray(r['corr'])[:, 1], 2).tolist()}")
        L.append("")
    L.append("C18j gate on CAL (in sample): 5/16 = 0.313 on U8 and 4/16 = 0.250 on N8 (NOCAP, LOC+T, 30/h).")
    open(outp, "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
