"""
r1_picks.py - Task R1 pick-based methods: PhaseNet (r0_picks: original model, threshold 0.3, 0.5 s minimum peak
distance) on the 100 Hz counts of the 7 stations; pick-based station terms (hyb_dev2.pick_terms over the CAL truth
events with magnitude >= 1.0, +-2 s around the JMA2001 times); PyOcto-1D and GaMMA-1D (hyb_comp: JMA2001 layers,
VS = VP / 1.73) on the corrected picks over the Task 2A grids (BASE and BESTPICK); R0-style comparators (homogeneous
VP 6.0, raw and +SC picks; descriptive).  Each stage is written to OUT.json as soon as it ends (resumable).
    python r1_picks.py SEGMENT OUT.json
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r1_lib as R  # noqa: E402

R.setup()
import hyb_comp as HC  # noqa: E402
import hyb_dev2 as H2  # noqa: E402
import r0_picks as RP  # noqa: E402
from task1r_baselines import GAMMA_GRID, PYOCTO_GRID  # noqa: E402


def save(res, path):
    json.dump(res, open(path + ".tmp", "w"))
    os.replace(path + ".tmp", path)


def main(segment, path):
    t0 = time.time()
    res = json.load(open(path)) if os.path.exists(path) else {"segment": segment, "names": R.SUBNETS["H7"]}
    net = R.make_net(segment)
    if "picks" not in res:
        P = RP.phasenet_hour(net.raw)
        res["picks"] = RP.picks_from(P)
        res["n_picks"] = {b: int(sum(1 for pl in res["picks"] for p in pl if p[1] == b)) for b in "PS"}
        save(res, path)
        print("phasenet", round(time.time() - t0), "s", res["n_picks"], flush=True)
    picks = [[tuple(p) for p in pl] for pl in res["picks"]]
    pa = H2.pick_arrays(picks)
    terms = H2.pick_terms(net, pa, R.T_CAL[0], R.T_CAL[1], mmin=R.M_TERMS)
    res["pick_terms"] = terms.tolist()
    picks_c = [sorted((t - terms[i, 0 if b == "P" else 1], b, h) for t, b, h in pl) for i, pl in enumerate(picks)]
    picks_sc = [sorted((t - net.corr_sc[i, 0 if b == "P" else 1], b, h) for t, b, h in pl) for i, pl in enumerate(picks)]
    stages = [("pyocto_1d", PYOCTO_GRID, lambda c: HC.pyocto_1d(net.STA, picks_c, **c)),
              ("gamma_1d", GAMMA_GRID, lambda c: HC.gamma_1d(net.STA, picks_c, **c)),
              ("pyocto", PYOCTO_GRID, lambda c: RP.pyocto_run(net.STA, picks, **c)),
              ("pyocto_sc", PYOCTO_GRID, lambda c: RP.pyocto_run(net.STA, picks_sc, **c)),
              ("gamma", GAMMA_GRID, lambda c: RP.gamma_run(net.STA, picks, **c)),
              ("gamma_sc", GAMMA_GRID, lambda c: RP.gamma_run(net.STA, picks_sc, **c))]
    for key, grid, fn in stages:
        res.setdefault(key, {})
        for cfg in grid:
            k = json.dumps(cfg)
            if k in res[key]:
                continue
            res[key][k] = fn(cfg)
            save(res, path)
            print(key, cfg, len(res[key][k]), round(time.time() - t0), "s", flush=True)
    res["corr_sc"] = net.corr_sc.tolist()
    res["done"] = True
    save(res, path)
    print("DONE", round(time.time() - t0), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
