"""Task 2M runner: one part file per (arm, cell, network); resumable (finished parts are skipped, writes are atomic).
    python t2m_run.py CONFIG.json ARM PARTSDIR [j0,j1,...|all] [cell,cell,...]
CONFIG.json: {"seed", "noise_dirs", "j0", "nblocks", "nwin": {"STD", "DENSE", "NOISE"}, "methods"} (t2m_lib.configure).
Windows per network: nwin["STD"] STD (kmax 5), nwin["DENSE"] DENSE (kmax 7), nwin["NOISE"] noise-only (kmax 7).
Part content = C16's part (c16_run.py) plus seed, day and station names; the rows are C16's rows.
"""
import json
import os
import sys
import time

import numpy as np

import t2m_lib as L2

KMAX = {"STD": 5, "DENSE": 7, "NOISE": 7}


def run_net(arm, cell, j, pdir):
    fn = os.path.join(pdir, f"{arm}_{cell}_{j}.json")
    if os.path.exists(fn):
        return
    t0 = time.time()
    net = L2.Net2M(cell, j, arm)
    t_cat = time.time() - t0
    net.setup()
    t_setup = time.time() - t0 - t_cat
    rows = []
    for design in ("STD", "DENSE", "NOISE"):
        nwin = L2.CFG["nwin"][design]
        if nwin <= 0:
            continue
        st = net.noise_windows() if design == "NOISE" else net.stream(design)
        for wi in range(nwin):
            w = next(st)
            atoms, sec = {}, {}
            for nm in L2.CFG["methods"]:
                t = time.time()
                atoms[nm] = L2.pack(net.run(nm, w["MFn"], KMAX[design]))
                sec[nm] = round(time.time() - t, 4)
            if design == "NOISE":
                w["Ds"] = np.zeros_like(w["D"])
                w["picks"] = [[] for _ in range(net.ns)]
                w["N"] = w["D"]
            snr, lev, kurt, mx = L2.snr_stats(w, net.ns)
            rows.append({"design": design, "w": wi, "k": len(w["truth"]), "truth": [t_.tolist() for t_ in w["truth"]],
                         "taus": w["taus"], "atoms": atoms, "sec": sec, "slot": w["slot"], "sd": w["sd"], "snr": snr,
                         "lev": lev, "kurt": kurt, "maxabs": mx})
    out = {"arm": arm, "cell": cell, "net": j, "sub": net.sub.tolist(), "STA": net.STA.tolist(),
           "statics": net.statics.tolist(), "terms": {k: np.asarray(v).tolist() for k, v in net.terms.items()},
           "err": {k: L2.L.err(np.asarray(v), net.statics) for k, v in net.terms.items()},
           "sig": net.sig, "s2info": net.s2info, "cat_drawn": net.drawn, "cat_slots": net.cat_slots,
           "slots": net.sl, "t_cat": t_cat, "t_setup": t_setup, "t_total": time.time() - t0, "rows": rows,
           "seed": L2.CFG["seed"], "day": net.day, "stations": net.stations}
    with open(fn + ".tmp", "w") as fh:
        json.dump(out, fh)
    os.replace(fn + ".tmp", fn)
    print(arm, cell, j, "day", net.day, "done", time.strftime("%H:%M:%S"), f"{time.time() - t0:.0f} s (cat {t_cat:.1f}, "
          f"setup {t_setup:.1f}); err", {k: round(v, 3) for k, v in out["err"].items()}, "drawn", net.drawn, flush=True)


if __name__ == "__main__":
    cfg, arm, pdir = sys.argv[1], sys.argv[2], sys.argv[3]
    assert arm in ("REAL", "GEN")
    L2.configure_file(cfg)
    nets = L2.CFG["nets"] if len(sys.argv) <= 4 or sys.argv[4] == "all" else [int(x) for x in sys.argv[4].split(",")]
    cells = sys.argv[5].split(",") if len(sys.argv) > 5 else L2.CELLS
    os.makedirs(pdir, exist_ok=True)
    for j in nets:
        for cell in cells:
            run_net(arm, cell, j, pdir)
