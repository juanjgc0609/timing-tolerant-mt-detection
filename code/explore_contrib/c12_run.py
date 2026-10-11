"""C12 runs.  SEED 999 ONLY.  One process; one part file per (kind, regime, cell, network) in out/c12_parts (resumable).
Networks task2a_confirm.network(999, cell, j): EVALUATION j = 1200-1215; NULL with clutter j = 1220-1225; NULL without
clutter (noise only) j = 1226-1229; validation j = 1290.  Stream TAGs [999, TAG, ci, j, ri, 41|233|61], all new:
4301 STD evaluation windows (xc.draw_sources, kmax 5), 4302 DENSE evaluation windows (c2_run.draw_dense, kmax 7),
4303 event-free windows with the clutter of the regime (kmax 5), 4304 event-free windows without clutter (kmax 5),
4305 validation.
    python c12_run.py eval j0,j1,... NWIN      |   python c12_run.py nullc j0,... NWIN   |   python c12_run.py nulln j0,... NWIN
"""
import json
import os
import sys
import time

import numpy as np

import c2_run as C2
import c10_lib as C
import c12_lib as L
import xc

CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
REGIMES = ("thin", "strong")
DESIGNS = {"STD": (4301, xc.draw_sources, 5), "DENSE": (4302, C2.draw_dense, 7)}
TAG_NULLC, TAG_NULLN = 4303, 4304
PARTS = os.path.join(L.HERE, "out", "c12_parts")
NCAND = 150


def detect(net, MFn, kmax):
    atoms, sec, stat = {}, {}, {}
    for key, rk in (("F6", "full6"), ("D5", "dev5"), ("DC", "dc")):
        t = time.time()
        atoms[key] = L.pack(L.loop(net.E, net.tab, MFn, kmax, rk, stat=stat))
        sec[key] = time.time() - t
    return atoms, sec, stat


def candidates(net, MFn, rng):
    """single-candidate statistics at NCAND random (node, origin) of the raw window: with the registered peak read
    (capped vectors) and at the predicted sample (no peak search, no cap).  Rows [s6, s5, sdc]."""
    E, tab = net.E, net.tab
    g, t = rng.integers(0, E.N, NCAND), rng.integers(0, E.T, NCAND)
    out = {}
    tpk, _ = C.choose(E, None, "hard", None, C.SH, R=MFn, sel=g)
    for name, tp, cap in (("peak", tpk[np.arange(NCAND), t], C.CAPV), ("fixed", E.IDX[g, t], None)):
        v = C._vectors(MFn, tp, cap).reshape(NCAND, 6 * E.ns)
        y = np.einsum("br,brk->bk", v, E.Q[g])
        s6 = (y * y).sum(-1) / 2.0
        out[name] = np.stack([s6, L.s5_of(y, tab.nvec[g]), L.dc_max(tab.Rm[g], tab.Ri[g], y)], axis=1).tolist()
    return out


def save(path, obj):
    with open(path + ".tmp", "w") as fh:
        json.dump(obj, fh)
    os.replace(path + ".tmp", path)


def run_eval(nets, nwin):
    for j in nets:
        for regime in REGIMES:
            for cell in CELLS:
                path = os.path.join(PARTS, f"eval_{regime}_{cell}_{j}.json")
                if os.path.exists(path):
                    continue
                t0 = time.time()
                net = L.Net(cell, j, regime)
                t_setup = time.time() - t0
                rows = []
                for design, (tag, draw, kmax) in DESIGNS.items():
                    st = L.stream(net.nc, tag, draw)
                    for wi in range(nwin):
                        w = next(st)
                        atoms, sec, stat = detect(net, w["MFn"], kmax)
                        rows.append({"design": design, "cell": cell, "net": j, "regime": regime, "w": wi,
                                     "k": len(w["truth"]), "truth": [t_.tolist() for t_ in w["truth"]], "taus": w["taus"],
                                     "amp": w["amp"], "amp_max": w["amp_max"], "atoms": atoms, "sec": sec,
                                     "ndc": stat.get("ndc", 0), "nit": stat.get("nit", 0)})
                save(path, {"rows": rows, "t_setup": t_setup, "t_total": time.time() - t0})
                print("eval", regime, cell, j, time.strftime("%H:%M:%S"), f"{time.time() - t0:.0f} s", flush=True)


def run_null(kind, nets, nwin):
    tag = TAG_NULLC if kind == "nullc" else TAG_NULLN
    for j in nets:
        for regime in (REGIMES if kind == "nullc" else ("thin",)):
            for cell in CELLS:
                path = os.path.join(PARTS, f"{kind}_{regime}_{cell}_{j}.json")
                if os.path.exists(path):
                    continue
                t0 = time.time()
                net = L.Net(cell, j, regime)
                st = L.stream(net.nc, tag, L.draw_none, nclut=None if kind == "nullc" else 0)
                rng = np.random.default_rng([L.SEED, tag, net.nc.ci, j, net.nc.ri, 977])
                rows = []
                for wi in range(nwin):
                    w = next(st)
                    atoms, sec, stat = detect(net, w["MFn"], 5)
                    rows.append({"cell": cell, "net": j, "regime": regime, "w": wi, "atoms": atoms, "sec": sec,
                                 "cand": candidates(net, w["MFn"], rng)})
                save(path, {"rows": rows, "t_total": time.time() - t0})
                print(kind, regime, cell, j, time.strftime("%H:%M:%S"), f"{time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    os.makedirs(PARTS, exist_ok=True)
    nets, nwin = [int(x) for x in sys.argv[2].split(",")], int(sys.argv[3])
    if sys.argv[1] == "eval":
        run_eval(nets, nwin)
    else:
        run_null(sys.argv[1], nets, nwin)
