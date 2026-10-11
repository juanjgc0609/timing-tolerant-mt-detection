"""C13: target and decoy runs of MT-FW (cap 6, registered +SC terms, hard peak rule).  SEED 999 ONLY.
Networks task2a_confirm.network(999, cell, j), j = 1400-1411 (never used before).  Streams [999, TAG, ci, j, ri, ...]:
  1401 STD windows (xc.draw_sources, 0-2 events, kmax 5)      1402 DENSE windows (c2_run.draw_dense, 3-5 events, kmax 7)
  1403 CLUT windows (no event, clutter of the regime, kmax 5)  1404 NOISE windows (no event, no clutter, kmax 5)
  1405 decoy draws ([999, 1405, ci, j, ri, design index, window])          1406 smoke test
One part file per (regime, cell, network): out/c13_parts/REGIME_CELL_J.json; finished parts are skipped.
    python c13_run.py j0,j1,... [NSTD NDENSE NCLUT NNOISE]
"""
import json
import os
import sys
import time

import numpy as np

import c2_run as C2
import c9_run as R9
import c10_lib as C
import c10_run as R
import c13_lib as L
import xc

CELLS = R.CELLS
REGIMES = ("strong", "thin")
PARTS = os.path.join(L.HERE, "out", "c13_parts")


def no_events(_rng):
    return [], []


# name: (tag, draw, kmax, clutter on)
DESIGNS = {"STD": (L.TAG_STD, xc.draw_sources, 5, True), "DENSE": (L.TAG_DENSE, C2.draw_dense, 7, True),
           "CLUT": (L.TAG_CLUT, no_events, 5, True), "NOISE": (L.TAG_NOISE, no_events, 5, False)}


class NoClutter:
    """view of a network configuration without clutter arrivals (noise-only windows)."""

    def __init__(self, nc):
        self.__dict__.update(nc.__dict__)
        self.nclut = 0


def part(regime, cell, j, nwin):
    t0 = time.time()
    S = R.Setup(cell, j, regime, oracle=False)
    nc = S.nc
    E = S.eng["cat"][0]
    dec = L.Decoys(E, nc.corr)
    t_setup = time.time() - t0
    rows = []
    for di, (design, (tag, draw, kmax, clut)) in enumerate(DESIGNS.items()):
        st = R9.stream(nc if clut else NoClutter(nc), tag, draw, None)
        for wi in range(nwin[design]):
            w = next(st)
            t = time.time()
            tgt = C.mt_loop(E, w["MFn"], kmax)
            t_t = time.time() - t
            rng = np.random.default_rng([L.SEED, L.TAG_DECOY, nc.ci, j, nc.ri, di, wi])
            t = time.time()
            d = dec.run_all(rng, w["MFn"], kmax)
            t_d = time.time() - t
            atoms = {"T": L.pack(tgt)}
            atoms.update({k: L.pack(v) for k, v in d.items()})
            rows.append({"design": design, "w": wi, "k": len(w["truth"]), "truth": [x.tolist() for x in w["truth"]],
                         "taus": w["taus"], "atoms": atoms, "sec_t": t_t, "sec_d": t_d})
    return {"regime": regime, "cell": cell, "net": j, "ns": int(E.ns), "t_setup": t_setup, "t_total": time.time() - t0,
            "rows": rows}


def main(nets, nwin):
    os.makedirs(PARTS, exist_ok=True)
    for j in nets:
        for regime in REGIMES:
            for cell in CELLS:
                path = os.path.join(PARTS, f"{regime}_{cell}_{j}.json")
                if os.path.exists(path):
                    continue
                res = part(regime, cell, j, nwin)
                with open(path + ".tmp", "w") as fh:
                    json.dump(res, fh)
                os.replace(path + ".tmp", path)
                print(time.strftime("%H:%M:%S"), regime, cell, j, f"{res['t_total']:.0f} s", flush=True)
    print("ALL DONE", time.strftime("%H:%M:%S"), flush=True)


if __name__ == "__main__":
    nets = [int(x) for x in sys.argv[1].split(",")]
    n = [int(x) for x in sys.argv[2:6]] if len(sys.argv) >= 6 else [20, 20, 10, 4]
    main(nets, dict(zip(DESIGNS, n)))
