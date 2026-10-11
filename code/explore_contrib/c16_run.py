"""C16 runner: one part file per (arm, cell, network).  SEED 999 ONLY.
    python c16_run.py ARM j0,j1,... NWIN PARTSDIR [cell,cell]
ARM in REAL | GEN.  Windows per network: NWIN STD (kmax 5), NWIN DENSE (kmax 7), NWIN noise-only (kmax 7).
Methods c16_lib.METHODS (unchanged code paths: c10_lib.mt_loop / bomp_loop, c8_lib.mtfw_general huber_n, c11_lib S2).
"""
import json
import os
import sys
import time

import numpy as np

import c16_lib as L16

KMAX = {"STD": 5, "DENSE": 7, "NOISE": 7}


def run_net(arm, cell, j, nwin, pdir):
    fn = os.path.join(pdir, f"{arm}_{cell}_{j}.json")
    if os.path.exists(fn):
        return
    t0 = time.time()
    net = L16.Net16(cell, j, arm)
    t_cat = time.time() - t0
    net.setup()
    t_setup = time.time() - t0 - t_cat
    rows = []
    for design in ("STD", "DENSE", "NOISE"):
        st = net.noise_windows() if design == "NOISE" else net.stream(design)
        for wi in range(nwin):
            w = next(st)
            atoms, sec = {}, {}
            for nm in L16.METHODS:
                t = time.time()
                atoms[nm] = L16.pack(net.run(nm, w["MFn"], KMAX[design]))
                sec[nm] = round(time.time() - t, 4)
            if design == "NOISE":
                w["Ds"] = np.zeros_like(w["D"])
                w["picks"] = [[] for _ in range(net.ns)]
                w["N"] = w["D"]
            snr, lev, kurt, mx = L16.snr_stats(w, net.ns)
            rows.append({"design": design, "w": wi, "k": len(w["truth"]), "truth": [t_.tolist() for t_ in w["truth"]],
                         "taus": w["taus"], "atoms": atoms, "sec": sec, "slot": w["slot"], "sd": w["sd"], "snr": snr,
                         "lev": lev, "kurt": kurt, "maxabs": mx})
    out = {"arm": arm, "cell": cell, "net": j, "sub": net.sub.tolist(), "STA": net.STA.tolist(),
           "statics": net.statics.tolist(), "terms": {k: np.asarray(v).tolist() for k, v in net.terms.items()},
           "err": {k: L16.L.err(np.asarray(v), net.statics) for k, v in net.terms.items()},
           "sig": net.sig, "s2info": net.s2info, "cat_drawn": net.drawn, "cat_slots": net.cat_slots,
           "slots": net.sl, "t_cat": t_cat, "t_setup": t_setup, "t_total": time.time() - t0, "rows": rows}
    with open(fn + ".tmp", "w") as fh:
        json.dump(out, fh)
    os.replace(fn + ".tmp", fn)
    print(arm, cell, j, "done", time.strftime("%H:%M:%S"), f"{time.time() - t0:.0f} s (cat {t_cat:.1f}, setup "
          f"{t_setup:.1f}); err", {k: round(v, 3) for k, v in out["err"].items()}, "drawn", net.drawn, flush=True)


if __name__ == "__main__":
    arm, nets, nwin, pdir = sys.argv[1], [int(x) for x in sys.argv[2].split(",")], int(sys.argv[3]), sys.argv[4]
    cells = sys.argv[5].split(",") if len(sys.argv) > 5 else L16.CELLS
    os.makedirs(pdir, exist_ok=True)
    for j in nets:
        for cell in cells:
            run_net(arm, cell, j, nwin, pdir)
