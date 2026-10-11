"""C13 correction run (declared post hoc in REPORT_C13.txt): decoy D3c = D3 with lags restricted so that the samples read
by the hypotheses stay inside the part of the trace where the generator places clutter.

Reason: the generator draws clutter arrivals in 2 s .. 36 s of the 40 s trace; the hypotheses read samples 95 .. 451.
The circular lags of D3 (and the offsets of D1o) are uniform over the whole trace and bring the clutter-free ends into
the samples read, so those decoys see fewer clutter arrivals than the targets (visible in the strong-clutter windows
without events).  D3c: every station trace is shifted (np.roll) by a lag between -7.5 s and +3.5 s; the lags of a
window are pairwise at least 1 s apart, uniformly over the configurations that satisfy this.
The windows are regenerated from the same streams as c13_run.py (checked against the stored targets on the first DENSE
window of every part).  Decoy draws: stream [999, 1405, ci, j, ri, design index, window, 7].
    python c13_run_d3c.py   -> out/c13_d3c_parts/REGIME_CELL_J.json
"""
import glob
import json
import os
import time

import numpy as np

import c9_run as R9
import c10_lib as C
import c10_run as R
import c13_lib as L
import c13_run as RR
from task0_common import DT

PARTS = os.path.join(L.HERE, "out", "c13_d3c_parts")
LAG_LO, LAG_HI, GAP = int(round(-7.5 / DT)), int(round(3.5 / DT)), int(round(1.0 / DT))


def lags_c(rng, n):
    """n integer lags in [LAG_LO, LAG_HI], pairwise at least GAP apart (uniform over such configurations)."""
    free = (LAG_HI - LAG_LO) - (n - 1) * GAP
    assert free >= 0
    base = np.sort(rng.integers(0, free + 1, n)) + GAP * np.arange(n) + LAG_LO
    return rng.permutation(base)


def main():
    os.makedirs(PARTS, exist_ok=True)
    for f in sorted(glob.glob(os.path.join(RR.PARTS, "*.json"))):
        out = os.path.join(PARTS, os.path.basename(f))
        if os.path.exists(out):
            continue
        t0 = time.time()
        p = json.load(open(f))
        S = R.Setup(p["cell"], p["net"], p["regime"], oracle=False)
        nc, E = S.nc, S.eng["cat"][0]
        rows = []
        for di, (design, (tag, draw, kmax, clut)) in enumerate(RR.DESIGNS.items()):
            st = R9.stream(nc if clut else RR.NoClutter(nc), tag, draw, None)
            stored = [r for r in p["rows"] if r["design"] == design]
            for wi, r in enumerate(stored):
                w = next(st)
                assert len(w["truth"]) == r["k"] and np.allclose(w["taus"], r["taus"])
                if design == "DENSE" and wi == 0:
                    assert np.allclose(np.array(L.pack(C.mt_loop(E, w["MFn"], kmax))), np.array(r["atoms"]["T"]))
                rng = np.random.default_rng([L.SEED, L.TAG_DECOY, nc.ci, p["net"], nc.ri, di, wi, 7])
                lag = lags_c(rng, E.ns)
                t = time.time()
                MFd = np.stack([np.roll(w["MFn"][i], int(lag[i]), axis=-1) for i in range(E.ns)])
                rows.append({"design": design, "w": wi, "D3c": L.pack(C.mt_loop(E, MFd, kmax)), "sec": time.time() - t})
        with open(out + ".tmp", "w") as fh:
            json.dump({"regime": p["regime"], "cell": p["cell"], "net": p["net"], "rows": rows}, fh)
        os.replace(out + ".tmp", out)
        print(time.strftime("%H:%M:%S"), os.path.basename(f), f"{time.time() - t0:.0f} s", flush=True)
    print("ALL DONE", time.strftime("%H:%M:%S"), flush=True)


if __name__ == "__main__":
    main()
