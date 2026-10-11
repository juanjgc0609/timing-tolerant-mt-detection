"""C11: station terms estimated self-consistently (S0-S4, ORACLE) inside MT-FW (cap 6) and BOMP-PK, with and without the
T3 refinement of C10.  SEED 999 ONLY.  Networks task2a_confirm.network(999, cell, j): DEVELOPMENT j = 200-209
(calibration quality only), EVALUATION j = 210-239; j = 195 validation.
Stream TAGs [999, TAG, ci, j, ri, 41|233|61], all new: 4101 STD evaluation windows (xc.draw_sources, kmax 5), 4102 DENSE
evaluation windows (c2_run.draw_dense, kmax 7), 4103 unlabelled calibration stream (xc.draw_sources; truth never read
by the estimator), 4104 noise-free validation catalogue, 4105 validation windows (network 195).
    python c11_run.py REGIME j0,j1,... NWIN out.json
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
import c10_screloc as SR
import c11_lib as L
import xc

CELLS = R.CELLS
DESIGNS = {"STD": (4101, xc.draw_sources, 5), "DENSE": (4102, C2.draw_dense, 7)}
TAG_UNL = 4103
TERMS_MT = ["S0", "S1", "S2", "S3", "S2@40", "S2@80", "ORA"]
TERMS_BOMP = ["S0", "S2", "ORA"]


def cat_k(nc, K):
    cat, _drawn, truth = xc.catalogue_cfg(nc.STA, K, nc.statics, [L.SEED, xc.T.TAG_CAT, nc.ci, nc.j, nc.ri], nc.nclut,
                                          nc.clut_amp)
    return cat, truth


class Net:
    """engines of every term set of a network (operators shared; only the sample tables differ)."""

    def __init__(self, cell, j, regime):
        self.S = R.Setup(cell, j, regime, oracle=True)
        self.nc = self.S.nc
        self.E0, _pr, self.EF0, _prF = self.S.eng["cat"]
        self.eng = {"S0": self.S.eng["cat"], "ORA": self.S.eng["ora"]}
        self.sig = {"S0": self.S.sig.tolist(), "ORA": self.S.sig_ora.tolist()}
        self.terms = {"S0": self.nc.corr.copy(), "ORA": self.nc.statics.copy()}

    def add(self, name, corr, cat):
        sig, _cnt = C.sigma_catalogue(self.nc.STA, cat, corr)
        E, EF = L.clone_eng(self.E0, corr), L.clone_eng(self.EF0, corr)
        self.eng[name] = (E, E.prior(sig), EF, EF.prior(sig))
        self.sig[name] = sig.tolist()
        self.terms[name] = np.asarray(corr, float).copy()

    def run(self, fam, name, MFn, kmax, t3):
        E, pr, EF, prF = self.eng[name]
        kw = dict(mode="hard", pr=pr, W=4, EF=EF if t3 else None, prF=prF)
        if fam == "MT":
            return C.mt_loop(E, MFn, kmax, **kw)
        return C.bomp_loop(E, None, MFn, kmax, **kw)


def unl_pass(net, name, n, tag=TAG_UNL):
    """MT-FW + T3 with the terms `name` over n windows of the unlabelled stream: atoms and energies per window (and
    the truth, stored for descriptive development diagnostics only; never passed to an estimator)."""
    st = R9.stream(net.nc, tag, xc.draw_sources, None)
    atoms, Ewin, truth = [], [], []
    for _ in range(n):
        w = next(st)
        atoms.append(net.run("MT", name, w["MFn"], 5, True))
        Ewin.append(np.sum(w["MFn"] ** 2, axis=1))
        truth.append((w["truth"], w["taus"]))
    return atoms, Ewin, truth


def calibrate(net, par=L.P, s3=True, ks=(40, 80)):
    """all term sets of a network; returns info (errors, timings)."""
    nc = net.nc
    info = {"t": {}}
    cat80, _truth = cat_k(nc, max(ks) if ks else 20)
    assert all(np.array_equal(a[1], b[1]) and a[2] == b[2] for a, b in zip(nc.cat, cat80[:20]))
    t = time.time()
    E80 = [L.energy(D) for D, _s, _t in cat80]
    info["t"]["energy_per_event"] = (time.time() - t) / len(cat80)
    t = time.time()
    sig0 = C.sigma_catalogue(nc.STA, nc.cat, nc.corr)[0][0]
    net.add("S1", SR.sc_reloc(nc.STA, nc.cat, nc.corr, sig0), nc.cat)
    info["t"]["S1"] = time.time() - t
    t = time.time()
    i2 = {}
    net.add("S2", L.s2_terms(nc.STA, nc.cat, nc.corr, E80[:20], par, info=i2), nc.cat)
    info["t"]["S2"] = time.time() - t
    info["S2"] = {k: i2[k] for k in ("sig", "n_res", "n_ev")}
    for K in ks:
        t = time.time()
        c0 = xc.estimate_station_corrections(nc.STA, cat80[:K])[0]
        net.terms[f"S0@{K}"] = c0
        net.add(f"S2@{K}", L.s2_terms(nc.STA, cat80[:K], c0, E80[:K], par), cat80[:K])
        info["t"][f"S2@{K}"] = time.time() - t
    if s3:
        t = time.time()
        start = "S2"
        for it in range(par["n_iter_s3"]):
            atoms, Ewin, _truth = unl_pass(net, start, par["n_unl"])
            i3 = {}
            c3 = L.s3_terms(nc.STA, nc.cat, nc.corr, net.terms["S2"], Ewin, atoms, E80[:20], par, info=i3)
            net.add("S3", c3, nc.cat)
            start = "S3"
        info["t"]["S3"] = time.time() - t
        info["S3"] = {k: i3[k] for k in ("sig", "n_res", "n_ev", "n_sel", "n_det_used")}
    info["err"] = {k: L.err(v, nc.statics) for k, v in net.terms.items()}
    info["off"] = {k: L.offset(v, nc.statics) for k, v in net.terms.items()}
    info["sigT3"] = dict(net.sig)
    return info


def run(regime, nets, nwin, path):
    res = json.load(open(path)) if os.path.exists(path) else {"rows": [], "nets": []}
    done = {(n["cell"], n["net"]) for n in res["nets"]}
    for j in nets:
        for cell in CELLS:
            if (cell, j) in done:
                continue
            t0 = time.time()
            net = Net(cell, j, regime)
            t_setup = time.time() - t0
            info = calibrate(net)
            t_cal = time.time() - t0 - t_setup
            for design, (tag, draw, kmax) in DESIGNS.items():
                st = R9.stream(net.nc, tag, draw, None)
                for wi in range(nwin):
                    w = next(st)
                    atoms, sec = {}, {}
                    for fam, names in (("MT", TERMS_MT), ("BOMP", TERMS_BOMP)):
                        for nm in names:
                            for t3 in (False, True):
                                if (fam, nm, t3) == ("BOMP", "ORA", False):
                                    continue
                                key = f"{fam}@{nm}" + ("+T3" if t3 else "")
                                t = time.time()
                                atoms[key] = R.pack(net.run(fam, nm, w["MFn"], kmax, t3))
                                sec[key] = time.time() - t
                    res["rows"].append({"design": design, "cell": cell, "net": j, "w": wi, "k": len(w["truth"]),
                                        "truth": [t_.tolist() for t_ in w["truth"]], "taus": w["taus"], "atoms": atoms,
                                        "sec": sec})
            res["nets"].append({"cell": cell, "net": j, "regime": regime, "t_total": time.time() - t0, "t_setup": t_setup,
                                "t_cal": t_cal, "info": info,
                                "terms": {k: np.asarray(v).tolist() for k, v in net.terms.items()}})
            with open(path + ".tmp", "w") as fh:
                json.dump(res, fh)
            os.replace(path + ".tmp", path)
            print(regime, cell, j, "done", time.strftime("%H:%M:%S"), f"{time.time() - t0:.0f} s; cal {t_cal:.0f} s; err",
                  {k: round(v, 3) for k, v in info["err"].items()}, flush=True)


if __name__ == "__main__":
    run(sys.argv[1], [int(x) for x in sys.argv[2].split(",")], int(sys.argv[3]), sys.argv[4])
