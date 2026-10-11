"""Check 3: budget handicap of TS-SFW+SC (f = 0.01) on seed-999 calibration windows (2A networks 0-4, windows 0..NW-1
of each, 3 cells x 2 regimes).  Variants: REG (registered: 50 it., 5 alternations), X10 (500 it., 20 alternations),
MS (9-start refinement, registered caps), X10MS (both).  The oracle score stays the registered one.
    OPENBLAS_CORETYPE=Prescott python3 check3_budget_verif.py NW"""
import json
import os
import sys
import time

import numpy as np

import common_verif as V

E, D2, T = V.E, V.D2, V.T
OUT = os.path.join(V.HERE, "verif_2e", "check3_parts_verif")
F = 0.01
_ORIG_ORACLE = E.oracle
_ORIG_REFINE = E.refine_atom
_STATE = {"R": None, "net": None}


def oracle_stash(net, R, excluded, stat="mtfw", PD=None):
    _STATE["R"] = R
    return _ORIG_ORACLE(net, R, excluded, stat, PD)


def refine_multistart(net, win, p0, a_old, T_old, maxiter=E.MAXIT_LOCAL):
    """9 starts: chosen node at t0 and t0 +- 0.1 s, and its 6 axis neighbours (+-3 km x, y; +-4 km z) at t0; the shifts
    of every start re-initialised from the energy peaks of the current sampled residual (as init_atom)."""
    R = _STATE["R"]
    Kt = D2.peak_table(R)
    x0 = p0[:3]
    starts = []
    for dx, dt in [((0, 0, 0), 0.0), ((0, 0, 0), -0.1), ((0, 0, 0), 0.1), ((3, 0, 0), 0.0), ((-3, 0, 0), 0.0),
                   ((0, 3, 0), 0.0), ((0, -3, 0), 0.0), ((0, 0, 4), 0.0), ((0, 0, -4), 0.0)]:
        x = np.clip(x0 + np.array(dx, float), net.lo[:3], net.hi[:3])
        t0 = float(np.clip(p0[3] + dt, net.lo[3], net.hi[3]))
        d = np.linalg.norm(net.STA - x, axis=1)
        T0 = d[None, :] * E.VINV[:, None] + t0 + net.corrT
        idx = np.clip((T0 / E.DT).astype(int), 0, E.NT - 1)
        I = Kt[np.arange(net.ns)[None, :], idx]
        dl = np.clip(I * E.DT - T0, -E.SHIFT_MAX, E.SHIFT_MAX)
        starts.append(np.concatenate([x, [t0], dl.ravel()]))
    starts[0] = p0.copy()            # the registered start itself
    best = None
    for s in starts:
        p, inf = _ORIG_REFINE(net, win, s, a_old, T_old, maxiter=maxiter)
        if best is None or inf["E1"] > best[1]["E1"]:
            best = (p, inf)
    best[1]["n_starts"] = len(starts)
    return best


def set_variant(v):
    big = v in ("X10", "X10MS")
    E.MAX_ALT = 20 if big else 5
    it = 500 if big else 50
    E.slide_step.__defaults__ = (it,)
    E.oracle = oracle_stash
    if v in ("MS", "X10MS"):
        E.refine_atom = lambda net, win, p0, a_old, T_old, maxiter=it: refine_multistart(net, win, p0, a_old, T_old,
                                                                                         maxiter)
    else:
        E.refine_atom = lambda net, win, p0, a_old, T_old, maxiter=it: _ORIG_REFINE(net, win, p0, a_old, T_old,
                                                                                    maxiter)


def run(nw, variants=("REG", "X10", "MS", "X10MS"), w0=0, suffix=""):
    os.makedirs(OUT, exist_ok=True)
    for cell in T.CELLS:
        for regime in T.REGIMES:
            for j in range(T.N_NETS):
                path = os.path.join(OUT, f"b{suffix}_{regime}_{cell}_net{j}.json")
                if os.path.exists(path):
                    continue
                t_net = time.time()
                S = V.setup(cell, regime, j)
                net = S["net"]
                stored = json.load(open(os.path.join(V.HERE, "task2e_parts", f"cal_seed999_{regime}_{cell}_net{j}.json")))
                W = []
                for w_i, w in enumerate(V.windows(S, cell, regime, j, nw)):
                    if w_i < w0:
                        continue
                    st = stored["windows"][w_i]
                    assert [list(map(float, s_)) for s_ in w["srcs"]] == st["truth"]
                    win = E.Win(w["D"], w["MFn"], w["sd"])
                    rec, dia = {}, {}
                    for v in variants:
                        set_variant(v)
                        dg = {}
                        t0 = time.time()
                        r = E.tssfw_record(net, win, w["truth"], f=F, diag=dg)
                        sec = time.time() - t0
                        rec[v] = [[float(s), int(h)] for s, h in r]
                        its = dg["iters"]
                        dia[v] = {"sec": sec, "J": [x["J"] for x in its], "pos": dg["P"][:, :4].tolist(),
                                  "alt": [x["alternations"] for x in its],
                                  "loc_nit": [x["local_nit"] for x in its],
                                  "ref_nit": [x["refine"]["nit"] for x in its],
                                  "ref_E1": [x["refine"]["E1"] for x in its],
                                  "frac_bound": float(np.mean(np.abs(dg["P"][:, 4:]) >= E.SHIFT_MAX - 1e-6)),
                                  "node": [x["node"] for x in its]}
                    rec["stored_REG"] = st["records"][f"TS-SFW+SC f={F}"]
                    rec["BOMP-PK+SC"] = st["records"]["BOMP-PK+SC"]
                    rec["MT-FW+SC"] = st["records"]["MT-FW+SC"]
                    rec["BOMP-PK-SLIDE+SC"] = st["records"]["BOMP-PK-SLIDE+SC"]
                    W.append({"w": w_i, "k": len(w["srcs"]), "truth": st["truth"], "rec": rec, "diag": dia})
                    print(f"{regime} {cell} net{j} w{w_i} k={len(w['srcs'])} " +
                          " ".join(f"{v}:{dia[v]['sec']:.1f}s" for v in variants) +
                          f" REG==stored {rec['REG'] == rec['stored_REG']}", flush=True)
                with open(path + ".tmp", "w") as fh:
                    json.dump({"cell": cell, "regime": regime, "net": j, "windows": W,
                               "seconds": time.time() - t_net}, fh)
                os.replace(path + ".tmp", path)


if __name__ == "__main__":
    if len(sys.argv) > 2:
        run(int(sys.argv[1]), tuple(sys.argv[3].split(",")), w0=int(sys.argv[2]), suffix="x")
    else:
        run(int(sys.argv[1]))
