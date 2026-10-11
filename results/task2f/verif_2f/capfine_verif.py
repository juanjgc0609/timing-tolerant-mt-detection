"""capfine_verif.py - finer CAPV sweep (4, 5, 8) for MT-FW on the seed-999 base dense windows (exploratory)."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json, sys
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mech_verif as V

def job(a):
    reg, cell, j = a
    out = os.path.join(V.OUT, f"capfine_{reg}_{cell}_net{j}.json")
    if os.path.exists(out):
        return out
    net = V.R.setup_net(999, cell, j, reg); net.update({"j": j, "regime": reg})
    ops = V.L.ModelOps(net["STA"], V.xc.make_grid(), V.L.Homog(), net["corr"])
    st = V.stream(999, net, V.C2.TAG, V.C2.draw_dense, True, None)
    wins = []
    for wi in range(10):
        w = next(st)
        rec = {f"MTFW_cap{c:g}": [list(map(float, x)) for x in V.mtfw(ops, w["MFn"], w["truth"], sel=c)] for c in (4.0, 5.0, 8.0)}
        wins.append({"k": len(w["truth"]), "rec": rec})
    json.dump({"windows": wins}, open(out, "w"))
    return out

if __name__ == "__main__":
    jobs = [(r, c, j) for r in ("thin", "strong") for c in V.C2.CELLS for j in V.NETS]
    with ProcessPoolExecutor(2) as ex:
        for o in ex.map(job, jobs):
            pass
    print("done")
