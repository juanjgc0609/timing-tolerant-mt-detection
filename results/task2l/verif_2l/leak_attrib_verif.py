"""leak_attrib_verif.py - (1) empirical leakage checks of S2 on fresh re-generation: catalogue carries errors, S2 depends
only on (STA, catalogue waveforms, catalogue hypocentres, S0), S2 reproduced exactly; (2) attribution of the contrast
between the terms and the T3 timing level (which the registered design recomputes with each term set): extra arms
MT@S2+T3 with the S0 level and MT@S0+T3 with the S2 level, STANDARD windows.  Writes only under verif_2l/."""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
HERE = os.path.dirname(os.path.abspath(__file__)); WORK = os.path.dirname(HERE)
sys.path.insert(0, WORK)
import numpy as np
import task2l_run as R
from concurrent.futures import ProcessPoolExecutor
OUT = os.path.join(HERE, "attrib_parts_verif")

def job(a):
    seed, regime, cell, j = a
    out = os.path.join(OUT, f"attrib_seed{seed}_{regime}_{cell}_net{j}.json")
    if os.path.exists(out):
        return out
    R.check_kernel(); R.set_seed(seed)
    L, C10, xc, R9, R10 = R.L, R.C10, R.xc, R.R9, R.R10
    net = R.R11.Net(cell, j, regime); nc = net.nc
    st_part = json.load(open(os.path.join(WORK, "task2l_parts", f"main_seed{seed}_{regime}_{cell}_net{j}.json")))
    # catalogue errors actually present in what the estimators see
    dx = np.array([c[1] - t[0] for c, t in zip(nc.cat, nc.cat_truth)]); dt = np.array([c[2] - t[1] for c, t in zip(nc.cat, nc.cat_truth)])
    # S2 reproduced; blind call on copies with the truth objects of the network destroyed
    c2 = L.s2_terms(nc.STA, nc.cat, nc.corr, [L.energy(D) for D, _s, _t in nc.cat], L.P)
    rep = bool(np.array_equal(c2, np.array(st_part["terms"]["S2"])))
    cat_copy = [(D.copy(), np.array(s, float), float(t)) for D, s, t in nc.cat]
    sta, corr = nc.STA.copy(), nc.corr.copy()
    keep_stat, keep_truth = nc.statics, nc.cat_truth
    nc.statics, nc.cat_truth = None, None
    c2b = L.s2_terms(sta, cat_copy, corr, None, dict(L.P))
    nc.statics, nc.cat_truth = keep_stat, keep_truth
    blind = bool(np.array_equal(c2, c2b))
    s0_rep = bool(np.array_equal(nc.corr, np.array(st_part["terms"]["S0"])))
    # engines
    sig0 = C10.sigma_catalogue(nc.STA, nc.cat, nc.corr)[0]; sig2 = C10.sigma_catalogue(nc.STA, nc.cat, c2)[0]
    E2, EF2 = L.clone_eng(net.E0, c2), L.clone_eng(net.EF0, c2)
    E0, EF0 = net.E0, net.EF0
    eng = {"MT@S2+T3": (E2, E2.prior(sig2), EF2, EF2.prior(sig2)),
           "MT@S2+T3|sigS0": (E2, E2.prior(sig0), EF2, EF2.prior(sig0)),
           "MT@S0+T3|sigS2": (E0, E0.prior(sig2), EF0, EF0.prior(sig2))}
    stream = R9.stream(nc, R.TAG_STD, xc.draw_sources, None)
    rows, same_truth, same_s2 = [], True, True
    for wi in range(R.N_WIN):
        w = next(stream)
        sw = st_part["windows"]["std"][wi]
        same_truth &= (np.array(sw["truth"]).tolist() == [t.tolist() for t in w["truth"]]) and sw["taus"] == w["taus"]
        atoms = {}
        for k, (E, pr, EF, prF) in eng.items():
            atoms[k] = R10.pack(C10.mt_loop(E, w["MFn"], R.KMAX["std"], mode="hard", pr=pr, W=4, EF=EF, prF=prF))
        same_s2 &= json.dumps(atoms["MT@S2+T3"]) == json.dumps(sw["atoms"]["MT@S2+T3"])
        rows.append({"w": wi, "atoms": {k: v for k, v in atoms.items() if k != "MT@S2+T3"}})
    obj = {"seed": seed, "regime": regime, "cell": cell, "net": j, "cat_dx_rms": float(np.sqrt((dx ** 2).mean())), "cat_dt_rms": float(np.sqrt((dt ** 2).mean())),
           "cat_dx": dx.tolist(), "cat_dt": dt.tolist(), "s2_reproduced": rep, "s2_blind_identical": blind, "s0_reproduced": s0_rep,
           "same_truth": bool(same_truth), "mt_s2t3_reproduced": bool(same_s2), "sig0": sig0.tolist(), "sig2": sig2.tolist(), "rows": rows}
    json.dump(obj, open(out + ".tmp", "w")); os.replace(out + ".tmp", out)
    return out

if __name__ == "__main__":
    lo, hi = int(sys.argv[1]), int(sys.argv[2])
    jobs = [(seed, reg, cell, j) for j in range(lo, hi) for seed in (7171, 9393) for cell in ("unilateral_c40", "ladder_6", "ladder_4") for reg in ("thin", "strong")]
    t0 = time.time()
    with ProcessPoolExecutor(2) as ex:
        for n, p in enumerate(ex.map(job, jobs)):
            if n % 20 == 0:
                print(n, len(jobs), round(time.time() - t0), os.path.basename(p), flush=True)
    print("done", round(time.time() - t0))
