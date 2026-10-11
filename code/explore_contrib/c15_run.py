"""C15 runner.  SEED 999 ONLY.  Networks task2a_confirm.network(999, cell, j): development j = 1800-1802, held-out
j = 1810-1821 (never used before).  Streams [999, TAG, ci, j, ri, 41|233|61]: 1851 STD windows (xc.draw_sources, 0-2
events, kmax 5), 1852 DENSE windows (c2_run.draw_dense, 3-5 events, kmax 7).
    python c15_run.py dev|eval j0,j1,... NSTD NDENSE    -> out/c15_<mode>_parts/<regime>_<cell>_<j>.json
"""
import json
import os
import sys
import time

import c2_run as C2
import c10_lib as C
import c10_run as R
import c14_lib as L14
import c15_lib as L
import xc

CELLS = R.CELLS
REGIMES = ("strong", "thin")
DESIGNS = {"STD": (1851, xc.draw_sources, 5), "DENSE": (1852, C2.draw_dense, 7)}


def same(a, b):
    return len(a) == len(b) and all((p["x"] == q["x"]).all() and p["tau"] == q["tau"] and p["s"] == q["s"]
                                    for p, q in zip(a, b))


def part(regime, cell, j, nwin):
    t0 = time.time()
    S = R.Setup(cell, j, regime, oracle=False)
    E = S.eng["cat"][0]
    rows = []
    for design, (tag, draw, kmax) in DESIGNS.items():
        st = L14.stream(S.nc, tag, draw, parts=False)
        for wi in range(nwin[design]):
            w = next(st)
            at = L.mt_loop15(E, w["MFn"], kmax, keep_et=not CYC)
            if wi < 2:
                assert same(at, C.mt_loop(E, w["MFn"], kmax))
            if CYC:
                pk = lambda a_: [[*map(float, a["x"]), float(a["tau"]), float(a["s"])] for a in a_]  # noqa: E731
                atoms = {"BASE": pk(at)}
                for nm, (cy, rs) in CYC.items():
                    atoms[nm] = pk(L.mt_loop_cyc(E, w["MFn"], kmax, cycles=cy, rescore=rs))
                rows.append({"design": design, "w": wi, "k": len(w["truth"]),
                             "truth": [x.tolist() for x in w["truth"]], "taus": w["taus"], "atoms": atoms,
                             "plateau": L.plateau_stats(E, at)})
                continue
            rows.append({"design": design, "w": wi, "k": len(w["truth"]), "truth": [x.tolist() for x in w["truth"]],
                         "taus": w["taus"], "atoms": L.pack_all(E, at), "plateau": L.plateau_stats(E, at)})
    return {"regime": regime, "cell": cell, "net": j, "ns": int(E.ns), "t_total": time.time() - t0, "rows": rows}


CYC = {}


def main(mode, nets, nwin):
    global CYC
    if mode.startswith("cyc"):
        CYC = {"CYC1": (1, True), "CYC2": (2, True), "CYC2g": (2, False), "CYC4": (4, True)}
    parts = os.path.join(L.HERE, "out", f"c15_{mode}_parts")
    os.makedirs(parts, exist_ok=True)
    for j in nets:
        for regime in REGIMES:
            for cell in CELLS:
                path = os.path.join(parts, f"{regime}_{cell}_{j}.json")
                if os.path.exists(path):
                    continue
                res = part(regime, cell, j, nwin)
                with open(path + ".tmp", "w") as fh:
                    json.dump(res, fh)
                os.replace(path + ".tmp", path)
                print(time.strftime("%H:%M:%S"), regime, cell, j, f"{res['t_total']:.0f} s", flush=True)
    print("ALL DONE", time.strftime("%H:%M:%S"), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], [int(x) for x in sys.argv[2].split(",")], {"STD": int(sys.argv[3]), "DENSE": int(sys.argv[4])})
