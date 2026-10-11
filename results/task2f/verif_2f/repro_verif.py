"""repro_verif.py - regenerate a few confirmatory 2F parts from scratch (fresh process, BOMP-PK run BEFORE MT-FW, so a
shared-state or call-order effect would show) and compare with the stored records (verification)."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json, sys
HERE = os.path.dirname(os.path.abspath(__file__)); WORK = os.path.dirname(HERE)
for p in (os.path.join(WORK, "explore_contrib"), WORK):
    sys.path.insert(0, p)
import task2f_run as R, c6_lib as L, c2_ablate as CA, xc, task2d_methods as D2
for arm, seed, reg, cell, j in (("dense", 4242, "thin", "ladder_6", 5), ("dense", 8888, "strong", "unilateral_c40", 77),
                                ("std", 8888, "strong", "ladder_4", 101)):
    spec = R.ARMS[arm]
    net = R.setup_net(seed, cell, j, reg); net.update({"j": j, "regime": reg})
    ops = L.ModelOps(net["STA"], xc.make_grid(), L.Homog(), net["corr"])
    PDc = CA.PKDictCap(ops.BD)
    st = R.window_stream(seed, net, spec["tag"], spec["draw"])
    e = json.load(open(os.path.join(WORK, "task2f_parts", f"{arm}_seed{seed}_{reg}_{cell}_net{j}.json")))
    same = {m: 0 for m in R.METHODS}
    for wi in range(10):
        w = next(st)
        out = {}
        out["BOMP-PK+SC"] = L.bomp_pk_m(ops, w["MFn"], w["truth"], spec["kmax"])[0]
        out["BOMP-PK-CAP+SC"] = [(float(s), int(h)) for s, h in D2.bomp_pk_record(PDc, w["MFn"], w["truth"], kmax=spec["kmax"])]
        out["MT-FW-NOCAP+SC"] = L.mtfw_dense(ops, w["MFn"], w["truth"], spec["kmax"], cap=False)[0]
        out["MT-FW+SC"] = L.mtfw_dense(ops, w["MFn"], w["truth"], spec["kmax"])[0]
        sw = e["windows"][wi]
        assert [list(map(float, t)) for t in w["truth"]] == sw["truth"]
        for m in R.METHODS:
            same[m] += [[float(s), int(h)] for s, h in out[m]] == sw["records"][m]
    print(arm, seed, reg, cell, j, "records identical (of 10):", same, flush=True)
