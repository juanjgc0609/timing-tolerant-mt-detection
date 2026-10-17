"""
r1_hybrid.py - Task R1 detection lists of one segment: BASE (the configuration of PyOcto-1D or GaMMA-1D with the highest
recall on CAL, LOC+T at 30 false per hour; ties: fewer CAL detections, then PyOcto, then the first of the grid), HYB
(hyb_dev2 PICKSCORE-UNION on BASE and P1), P1, BESTPICK (every 1-D configuration) and the R0-style comparators.
    python r1_hybrid.py SEGMENT PICKS.json P1.json OUT.json
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r1_lib as R  # noqa: E402

R.setup()
import hyb_dev as H  # noqa: E402
import hyb_dev2 as H2  # noqa: E402
import r0_analyze as A  # noqa: E402
import r1_run  # noqa: E402
from task1r_baselines import GAMMA_GRID, PYOCTO_GRID  # noqa: E402

TOL = 0.5          # s, pick support
W_WAVE = 0.5       # weight of the waveform-score percentile
DKM, DTS = 5.0, 1.5


def in_seg(evs, lo, hi):
    return [list(e[:5]) for e in evs if lo <= e[3] < hi]


def select_base(pk, cat):
    rows = []
    for fam, (key, grid) in enumerate((("pyocto_1d", PYOCTO_GRID), ("gamma_1d", GAMMA_GRID))):
        for gi, cfg in enumerate(grid):
            c = json.dumps(cfg)
            ev = in_seg(pk[key][c], *R.T_CAL)
            rec = H.evaluate(ev, cat, *R.T_CAL)[("LOC+T", 30)]
            rows.append((-rec, len(ev), fam, gi, key, c))
    rows.sort()
    return rows[0][4], rows[0][5], [(k, c, -r, n) for r, n, _f, _g, k, c in rows]


def hybrid(base, p1, wb, wp, STA, pa, terms):
    cands = [e[:4] + [e[4] + 100.0, wb[i]] for i, e in enumerate(base)] + [a[:4] + [wp[i], wp[i]] for i, a in enumerate(p1)]
    sc = H2.score_candidates(cands, STA, pa, terms, TOL)
    dets = sorted([[c[0], c[1], c[2], c[3], c[4] + W_WAVE * c[5]] for c in sc], key=lambda r: -r[4])
    keep = []
    for r in dets:
        if all(np.hypot(r[0] - k[0], r[1] - k[1]) > DKM or abs(r[3] - k[3]) > DTS for k in keep):
            keep.append(r)
    return keep


def main(segment, pkp, p1p, outp):
    lo, hi = R.T_CAL if segment == "cal" else R.T_EVAL
    net = R.make_net(segment)
    pk = json.load(open(pkp))
    E = r1_run.engine(net, pk["pick_terms"])[0]
    cat = net.cat
    bkey, bcfg, table = select_base(pk, cat)
    base = in_seg(pk[bkey][bcfg], lo, hi)
    p1 = A.dedupe(in_seg(json.load(open(p1p))["atoms"]["NOCAP"], lo, hi))
    wb = H.pct(H.waveform_scores(net, E, base, "NOCAP")).tolist()
    wp = H.pct([a[4] for a in p1]).tolist()
    pa = H2.pick_arrays([[tuple(p) for p in pl] for pl in pk["picks"]])
    terms = np.array(pk["pick_terms"])
    hyb = hybrid(base, p1, wb, wp, net.STA, pa, terms)
    out = {"segment": segment, "base": [bkey, bcfg], "selection": table, "BASE": base, "HYB": hyb, "P1": p1,
           "BESTPICK": {f"{k}|{c}": in_seg(v, lo, hi) for k in ("pyocto_1d", "gamma_1d") for c, v in pk[k].items()},
           "R0STYLE": {f"{k}|{c}": in_seg(v, lo, hi) for k in ("gamma", "pyocto", "gamma_sc", "pyocto_sc")
                       for c, v in pk[k].items()}}
    json.dump(out, open(outp + ".tmp", "w"))
    os.replace(outp + ".tmp", outp)
    print(segment, "BASE", bkey, bcfg, "CAL recall", table[0][2], "n", table[0][3], "| BASE", len(base), "P1", len(p1),
          "HYB", len(hyb), flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:5])
