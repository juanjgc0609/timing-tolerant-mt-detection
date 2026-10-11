"""C14: self-interference remedies inside MT-FW (cap 6, registered +SC terms, hard peak rule).  SEED 999 ONLY.
Networks task2a_confirm.network(999, cell, j): DEVELOPMENT j = 1600-1602, HELD-OUT j = 1610-1621 (never used before).
Streams [999, TAG, ci, j, ri, 41|233|61]: 1651 STD windows (xc.draw_sources, 0-2 events, kmax 5), 1652 DENSE windows
(c2_run.draw_dense, 3-5 events, kmax 7).  One part file per (regime, cell, network); finished parts are skipped.
    python c14_run.py dev  j0,j1,... NSTD NDENSE          -> out/c14_dev_parts   (grid DEV + diagnosis + ceilings)
    python c14_run.py dev2 j0,j1,... NSTD NDENSE          -> out/c14_dev2_parts  (combinations, same dev windows)
    python c14_run.py eval j0,j1,... NSTD NDENSE          -> out/c14_eval_parts  (frozen set FROZEN + twins)
"""
import json
import os
import sys
import time

import numpy as np

import c2_run as C2
import c10_lib as C
import c10_run as R
import c14_lib as L
import xc

CELLS = R.CELLS
NTWIN = 8                       # eval: event-free twins of the first NTWIN windows of every design (compute budget)
REGIMES = ("strong", "thin")
DESIGNS = {"STD": (L.TAG_STD, xc.draw_sources, 5), "DENSE": (L.TAG_DENSE, C2.draw_dense, 7)}

# development grid (round 1): one family at a time
DEV = {
    "R1h4": {"r1": ("hard", 4)}, "R1h8": {"r1": ("hard", 8)},
    "R1s8p1": {"r1": ("soft", 8, 1.0)}, "R1s8p05": {"r1": ("soft", 8, 0.5)},
    "R2a32": {"r2": (32, False)}, "R2a64": {"r2": (64, False)}, "R2f64": {"r2": (64, True)},
    "R3amp2": {"r3": ("amp", 2)}, "R3vec2": {"r3": ("vec", 2)}, "R3fir2": {"r3": ("fir", 2)},
    "R4x10t05": {"r4": (10.0, 0.5)}, "R4x10t10": {"r4": (10.0, 1.0)},
}
# round 2 (combinations and the constants added after round 1) and the frozen set: filled in after round 1
DEV2 = json.load(open(os.path.join(L.HERE, "c14_dev2.json"))) if os.path.exists(os.path.join(L.HERE, "c14_dev2.json")) \
    else {}
FROZEN = json.load(open(os.path.join(L.HERE, "c14_frozen.json"))) if os.path.exists(
    os.path.join(L.HERE, "c14_frozen.json")) else {}


def _cfg(d):
    return {k: tuple(v) for k, v in d.items()}


def part(mode, regime, cell, j, nwin):
    t0 = time.time()
    S = R.Setup(cell, j, regime, oracle=False)
    E, _pr, EF, prF = S.eng["cat"]
    rows = []
    sec = {}

    def run(nm, fn):
        t = time.time()
        a = fn()
        sec[nm] = sec.get(nm, 0.0) + time.time() - t
        return a
    for design, (tag, draw, kmax) in DESIGNS.items():
        st = L.stream(S.nc, tag, draw, parts=(mode == "dev"))
        for wi in range(nwin[design]):
            w = next(st)
            M, k = w["MFn"], len(w["truth"])
            base = run("BASE", lambda: C.mt_loop(E, M, kmax))
            row = {"design": design, "w": wi, "k": k, "truth": [x.tolist() for x in w["truth"]], "taus": w["taus"],
                   "atoms": {"BASE": L.pack(base)}, "twin": {}}
            if mode == "dev":
                dg = dict(dir=w["dir"], coda=w["coda"], arr=w["arr"], cn=w["MFn0"], truth=w["truth"], taus=w["taus"])
                b2 = L.mt_loop14(E, M, kmax, diag=dg)
                assert L.same(base, b2)
                row["dg"] = [a["dg"] for a in b2]
                row["arr"] = w["arr"].tolist()
                for nm, cf in DEV.items():
                    row["atoms"][nm] = L.pack(run(nm, lambda: L.mt_loop14(E, M, kmax, cf)))
                row["atoms"]["T3"] = L.pack(run("T3", lambda: C.mt_loop(E, M, kmax, EF=EF, prF=prF)))
                if k:
                    row["twin"]["BASE"] = L.pack(C.mt_loop(E, w["MFn0"], kmax))
                    ev = w["dir"] + w["coda"]
                    row["atoms"]["NOCODA"] = L.pack(L.mt_loop14(E, M - w["coda"].sum(0), kmax))
                    row["atoms"]["ORAdir"] = L.pack(L.mt_loop14(E, M, kmax, ora={"ev": w["dir"], "truth": w["truth"]}))
                    row["atoms"]["ORAall"] = L.pack(L.mt_loop14(E, M, kmax, ora={"ev": ev, "truth": w["truth"]}))
                else:
                    for nm in ("NOCODA", "ORAdir", "ORAall"):
                        row["atoms"][nm] = row["atoms"]["BASE"]
            else:
                if wi < 2:
                    assert L.same(base, L.mt_loop14(E, M, kmax))
                grid = DEV2 if mode == "dev2" else FROZEN
                for nm, cf in grid.items():
                    if nm == "T3":
                        row["atoms"][nm] = L.pack(run(nm, lambda: C.mt_loop(E, M, kmax, EF=EF, prF=prF)))
                        fn0 = lambda: C.mt_loop(E, w["MFn0"], kmax, EF=EF, prF=prF)  # noqa: E731
                    else:
                        cf = _cfg(cf)
                        row["atoms"][nm] = L.pack(run(nm, lambda: L.mt_loop14(E, M, kmax, cf)))
                        fn0 = lambda: L.mt_loop14(E, w["MFn0"], kmax, cf)  # noqa: E731
                    if mode == "eval" and k and wi < NTWIN and nm != "T3":
                        row["twin"][nm] = L.pack(fn0())
                if mode == "eval" and k and wi < NTWIN:
                    row["twin"]["BASE"] = L.pack(C.mt_loop(E, w["MFn0"], kmax))
            rows.append(row)
    return {"regime": regime, "cell": cell, "net": j, "ns": int(E.ns), "t_setup": S.t_setup,
            "t_total": time.time() - t0, "sec": sec, "nwin": sum(nwin.values()), "rows": rows}


def main(mode, nets, nwin):
    parts = os.path.join(L.HERE, "out", f"c14_{mode}_parts")
    os.makedirs(parts, exist_ok=True)
    for j in nets:
        for regime in REGIMES:
            for cell in CELLS:
                path = os.path.join(parts, f"{regime}_{cell}_{j}.json")
                if os.path.exists(path):
                    continue
                res = part(mode, regime, cell, j, nwin)
                with open(path + ".tmp", "w") as fh:
                    json.dump(res, fh)
                os.replace(path + ".tmp", path)
                print(time.strftime("%H:%M:%S"), regime, cell, j, f"{res['t_total']:.0f} s", flush=True)
    print("ALL DONE", time.strftime("%H:%M:%S"), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], [int(x) for x in sys.argv[2].split(",")], {"STD": int(sys.argv[3]), "DENSE": int(sys.argv[4])})
