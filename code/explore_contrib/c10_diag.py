"""C10 diagnostic (POST HOC, designed after the partial look at the 4-station dense windows showed no gain even for
the oracle): first-iteration separability.  On the raw window (no deflation) the node statistic of the MT-FW family is
evaluated under four rules; for every true event, its best value within 5 km (epicentral; any origin) is compared with
the best value of the window among nodes farther than 5 km from EVERY true event (the best false alternative).
  hard      registered statistic, coarse lattice
  T2w-cat   corr prior, +-0.5 s, catalogue level and corrections, coarse lattice
  T2w-ora   corr prior, +-0.5 s, true statics and true model-error level, coarse lattice
  fine-ora  indep objective on the 1 km lattice, true statics, true model-error level + fine cell term
  time-ora  as fine-ora, but the location of the maximum within 10 km of the true event is reported (epicentral error)
Networks 160-169, the main windows (streams 4001 / 4002 regenerated).
    python c10_diag.py DESIGN REGIME j0,... NWIN out.json"""
import json
import sys

import numpy as np

import c2_run as C2
import c9_run as R9
import c10_lib as C
import c10_run as R
import xc


def score_map(E, MFn, mode, pr, W, chunk=800):
    Eh = C._half_energy(MFn, C.CAPV)
    ns = E.ns
    if mode != "fine":
        tpk, pen = C.choose(E, Eh, mode, pr, W, R=MFn)
        vf = C._vectors(MFn, tpk, C.CAPV).reshape(E.N, E.T, 6 * ns)
        y = np.einsum("gtr,grk->gtk", vf, E.Q)
        s = (y * y).sum(-1) / 2.0
        return s if pen is None else s - pen
    out = np.zeros((E.N, E.T))
    tall = np.arange(E.T)
    for c0 in range(0, E.N, chunk):
        sel = np.arange(c0, min(E.N, c0 + chunk))
        tpk, pen = C.choose(E, Eh, "indep", pr, W, sel=sel, tsel=tall)
        vf = C._vectors(MFn, tpk, C.CAPV).reshape(len(sel), E.T, 6 * ns)
        y = np.einsum("gtr,grk->gtk", vf, E.Q[sel])
        out[sel] = (y * y).sum(-1) / 2.0 - pen
    return out


def run(design, regime, nets, nwin, path):
    tag, kind = R.DESIGNS[design]
    draw = C2.draw_dense if kind == "dense" else xc.draw_sources
    rows = []
    for cell in R.CELLS:
        for j in nets:
            S = R.Setup(cell, j, regime)
            E, pr, EF, prF = S.eng["cat"]
            EO, prO, EFO, prFO = S.eng["ora"]
            st = R9.stream(S.nc, tag, draw, None)
            for wi in range(nwin):
                w = next(st)
                tr = w["truth"]
                if not tr:
                    continue
                maps = {"hard": (E, score_map(E, w["MFn"], "hard", None, 4)),
                        "T2w-cat": (E, score_map(E, w["MFn"], "corr", pr, 8)),
                        "T2w-ora": (EO, score_map(EO, w["MFn"], "corr", prO, 8)),
                        "fine-ora": (EFO, score_map(EFO, w["MFn"], "fine", prFO, 4))}
                row = {"cell": cell, "net": j, "w": wi, "k": len(tr), "ev": {}}
                for nm, (En, sc) in maps.items():
                    best = sc.max(axis=1)
                    dd = np.linalg.norm(En.XN[:, None, :2] - np.array(tr)[None, :, :2], axis=2)      # (N, k)
                    far = dd.min(axis=1) > 5.0
                    row["ev"][nm] = {"true": [float(best[dd[:, k] <= 5.0].max()) for k in range(len(tr))],
                                     "false": float(best[far].max())}
                    if nm == "fine-ora":
                        err = []
                        for k in range(len(tr)):
                            near = np.where(dd[:, k] <= 10.0)[0]
                            err.append(float(dd[near[int(np.argmax(best[near]))], k]))
                        row["ev"]["time-ora"] = err
                rows.append(row)
            print("DIAG", design, regime, cell, j, flush=True)
            with open(path, "w") as fh:
                json.dump(rows, fh)


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], [int(x) for x in sys.argv[3].split(",")], int(sys.argv[4]), sys.argv[5])
