"""C11 development: calibration quality (rms error of the terms against the true statics) on the DEVELOPMENT networks
j = 200-209.  No detection recall is computed here.  python c11_dev.py s2 | s3pass REGIME | s3
"""
import itertools
import pickle
import sys
import time

import numpy as np

import c10_lib as C
import c10_screloc as SR
import c11_lib as L
import xc

CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
REGIMES = ("thin", "strong")
DEV = list(range(200, 210))


def cat_k(nc, K):
    cat, _drawn, truth = xc.catalogue_cfg(nc.STA, K, nc.statics, [L.SEED, xc.T.TAG_CAT, nc.ci, nc.j, nc.ri], nc.nclut,
                                          nc.clut_amp)
    return cat, truth


def load_dev():
    out = {}
    for reg in REGIMES:
        for cell in CELLS:
            for j in DEV:
                nc = xc.NetCfg(cell, j, reg)
                cat80, _tr = cat_k(nc, 80)
                assert all(np.array_equal(a[1], b[1]) and a[2] == b[2] for a, b in zip(nc.cat, cat80[:20]))
                out[(reg, cell, j)] = (nc, cat80, [L.energy(D) for D, _s, _t in cat80])
    return out


def table(P_, name, vals):
    """vals[(reg, cell)] = list over networks"""
    s = f"{name:34s}"
    for reg in REGIMES:
        for cell in CELLS:
            s += f" {np.median(vals[(reg, cell)]):.3f}"
        s += f" | {np.median(sum((vals[(reg, c)] for c in CELLS), [])):.3f} ||"
    P_(s)


def dev_s2(out="out/c11_dev_s2.txt"):
    lines = []

    def P_(s=""):
        print(s, flush=True)
        lines.append(s)
    t0 = time.time()
    dev = load_dev()
    P_(f"C11 DEV S2 (networks 200-209; rms error of the terms against the true statics, median removed, s; median over "
       f"10 networks; columns thin: ladder_4 ladder_6 unilateral_c40 | pooled || strong: ... ); load {time.time()-t0:.0f} s")
    res = {}

    def run(name, fn):
        v = {}
        for (reg, cell, j), (nc, cat80, E80) in dev.items():
            v.setdefault((reg, cell), []).append(L.err(fn(nc, cat80, E80), nc.statics))
        res[name] = v
        table(P_, name, v)
    run("S0 (+SC, K=20)", lambda nc, c, E: nc.corr)
    run("S1 (SC-R, K=20)", lambda nc, c, E: SR.sc_reloc(nc.STA, nc.cat, nc.corr, C.sigma_catalogue(nc.STA, nc.cat, nc.corr)[0][0]))
    for K in (40, 80):
        run(f"S0 median, K={K}", lambda nc, c, E, K=K: xc.estimate_station_corrections(nc.STA, c[:K])[0])
    base = dict(L.P)
    grid = [("default", {})]
    for w in (0.3, 0.75, 1.0):
        grid.append((f"w_re={w}", {"w_re": w}))
    for npass in (1, 2, 5):
        grid.append((f"n_pass={npass}", {"n_pass": npass}))
    grid += [("term=huber", {"term": "huber"}), ("sx=2", {"sx": 2.0}), ("sx=0.5", {"sx": 0.5}), ("st=0.2", {"st": 0.2}),
             ("min_arr=5", {"min_arr": 5}), ("min_arr=3", {"min_arr": 3}), ("huber=1.0", {"huber": 1.0}),
             ("huber=2.5", {"huber": 2.5}), ("gate_re=16", {"gate_re": 16.0}), ("gate_re=9", {"gate_re": 9.0}),
             ("gate_re=16 w_re=0.3", {"gate_re": 16.0, "w_re": 0.3}), ("gate_re=9 w_re=0.3", {"gate_re": 9.0, "w_re": 0.3}),
             ("gate_re=16 huber", {"gate_re": 16.0, "term": "huber"})]
    for nm, kw in grid:
        par = dict(base)
        par.update(kw)
        run(f"S2 K=20 {nm}", lambda nc, c, E, par=par: L.s2_terms(nc.STA, c[:20], nc.corr, E[:20], par))
    for K in (40, 80):
        run(f"S2 K={K} default", lambda nc, c, E, K=K: L.s2_terms(
            nc.STA, c[:K], xc.estimate_station_corrections(nc.STA, c[:K])[0], E[:K], base))
    P_(f"total {time.time()-t0:.0f} s")
    open(out, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    if sys.argv[1] == "s2":
        dev_s2(*(sys.argv[2:3]))


def dev_s2b(out="out/c11_dev_s2b.txt"):
    """second grid (after the first showed that the energy gate of the re-read is the lever)."""
    lines = []

    def P_(s=""):
        print(s, flush=True)
        lines.append(s)
    dev = load_dev()
    P_("C11 DEV S2, second grid (same layout as c11_dev_s2.txt)")
    base = dict(L.P)
    grid = []
    for g, m, w in itertools.product((4.0, 9.0, 12.0, 16.0), (3, 4), (0.3, 0.5)):
        grid.append((f"gate_re={g:g} min_arr={m} w_re={w}", {"gate_re": g, "min_arr": m, "w_re": w}))
    for npass in (2, 4):
        grid.append((f"gate_re=12 n_pass={npass}", {"gate_re": 12.0, "n_pass": npass}))
    grid += [("gate_re=12 huber", {"gate_re": 12.0, "term": "huber"}), ("gate_re=12 sx=2", {"gate_re": 12.0, "sx": 2.0}),
             ("gate_re=12 st=0.2", {"gate_re": 12.0, "st": 0.2})]
    for K in (20, 40, 80):
        for nm, kw in (grid if K == 20 else [("gate_re=12", {"gate_re": 12.0}), ("gate_re=9", {"gate_re": 9.0}),
                                             ("gate_re=16", {"gate_re": 16.0})]):
            par = dict(base)
            par.update(kw)
            v = {}
            for (reg, cell, j), (nc, cat80, E80) in dev.items():
                c0 = nc.corr if K == 20 else xc.estimate_station_corrections(nc.STA, cat80[:K])[0]
                v.setdefault((reg, cell), []).append(L.err(L.s2_terms(nc.STA, cat80[:K], c0, E80[:K], par), nc.statics))
            table(P_, f"S2 K={K} {nm}", v)
    open(out, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__" and sys.argv[1] == "s2b":
    dev_s2b()


def dev_s3pass(reg, nwin=40):
    """detection pass (MT-FW + T3 with the S2 terms) over the unlabelled stream of the development networks; stores
    atoms and energies.  No recall is computed."""
    import c11_run as RR
    out = {}
    for cell in CELLS:
        for j in DEV:
            t0 = time.time()
            net = RR.Net(cell, j, reg)
            nc = net.nc
            E20 = [L.energy(D) for D, _s, _t in nc.cat]
            net.add("S2", L.s2_terms(nc.STA, nc.cat, nc.corr, E20), nc.cat)
            atoms, Ewin, truth = RR.unl_pass(net, "S2", nwin)
            out[(reg, cell, j)] = {"S2": net.terms["S2"], "atoms": [[{k: a[k] for k in ("x", "tau", "s")} for a in at]
                                                                    for at in atoms],
                                   "Ewin": [e.astype(np.float32) for e in Ewin], "truth": truth}
            print(reg, cell, j, f"{time.time() - t0:.0f} s", flush=True)
    pickle.dump(out, open(f"out/c11_dev_s3pass_{reg}.pkl", "wb"))


def dev_s3(out="out/c11_dev_s3.txt"):
    lines = []

    def P_(s=""):
        print(s, flush=True)
        lines.append(s)
    dev = load_dev()
    ps = {}
    for reg in REGIMES:
        ps.update(pickle.load(open(f"out/c11_dev_s3pass_{reg}.pkl", "rb")))
    P_("C11 DEV S3 (networks 200-209; unlabelled stream TAG 4103; layout as c11_dev_s2.txt).  purity = share of the "
       "selected atoms within 5 km of a true event (descriptive; uses the truth; not used by any rule)")
    base = dict(L.P)
    v0, v2 = {}, {}
    for key, (nc, cat80, E80) in dev.items():
        v0.setdefault(key[:2], []).append(L.err(nc.corr, nc.statics))
        v2.setdefault(key[:2], []).append(L.err(ps[key]["S2"], nc.statics))
    s2, v2n = {}, {}
    for key, (nc, cat80, E80) in dev.items():
        s2[key] = L.s2_terms(nc.STA, nc.cat, nc.corr, E80[:20])
        v2n.setdefault(key[:2], []).append(L.err(s2[key], nc.statics))
    table(P_, "S0", v0)
    table(P_, "S2 of the detection pass", v2)
    table(P_, "S2 final (start of S3)", v2n)
    grid = [("default", {})]
    for n in (20, 40):
        grid.append((f"n_unl={n}", {"n_unl": n}))
    for q in (0.6, 0.7, 0.9):
        grid.append((f"q_det={q}", {"q_det": q}))
    grid += [("sx_det=1.5", {"sx_det": 1.5}), ("sx_det=6", {"sx_det": 6.0}), ("min_arr_det=4", {"min_arr_det": 4}),
             ("min_arr_det=8", {"min_arr_det": 8}), ("n_pass_s3=1", {"n_pass_s3": 1}), ("n_pass_s3=3", {"n_pass_s3": 3})]
    for nm, kw in grid:
        par = dict(base)
        par.update(kw)
        v, pur, nsel, nuse = {}, [], [], []
        for key, (nc, cat80, E80) in dev.items():
            d = ps[key]
            n = par["n_unl"]
            info = {}
            c = L.s3_terms(nc.STA, nc.cat, nc.corr, s2[key], d["Ewin"][:n], d["atoms"][:n], E80[:20], par, info=info)
            v.setdefault(key[:2], []).append(L.err(c, nc.statics))
            sel = L.select_atoms(d["atoms"][:n], par)
            for wi, a in sel:
                tr = d["truth"][wi][0]
                pur.append(float(any(np.linalg.norm(a["x"][:2] - t[:2]) <= 5.0 for t in tr)))
            nsel.append(info["n_sel"])
            nuse.append(info["n_det_used"])
        table(P_, f"S3 {nm}", v)
        P_(f"{'':34s} selected atoms per network {np.mean(nsel):.1f}, used {np.mean(nuse):.1f}, purity {np.mean(pur):.2f}")
    open(out, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__" and sys.argv[1] == "s3pass":
    dev_s3pass(sys.argv[2])
if __name__ == "__main__" and sys.argv[1] == "s3":
    dev_s3()


def dev_s2c(out="out/c11_dev_s2c.txt"):
    """third grid: narrowing schedules of the re-read window."""
    lines = []

    def P_(s=""):
        print(s, flush=True)
        lines.append(s)
    dev = load_dev()
    P_("C11 DEV S2, third grid: schedules of the re-read half-width (same layout)")
    base = dict(L.P)
    grid = []
    for g in (12.0, 16.0):
        grid += [(f"gate_re={g:g} w_re=0.5,0.3", {"gate_re": g, "w_re": [0.5, 0.3]}),
                 (f"gate_re={g:g} w_re=0.5,0.3,0.2 n_pass=4", {"gate_re": g, "w_re": [0.5, 0.3, 0.2], "n_pass": 4}),
                 (f"gate_re={g:g} w_re=0.75,0.5,0.3 n_pass=4", {"gate_re": g, "w_re": [0.75, 0.5, 0.3], "n_pass": 4}),
                 (f"gate_re={g:g} w_re=0.5", {"gate_re": g})]
    for nm, kw in grid:
        par = dict(base)
        par.update(kw)
        for K in (20, 40):
            v = {}
            for (reg, cell, j), (nc, cat80, E80) in dev.items():
                c0 = nc.corr if K == 20 else xc.estimate_station_corrections(nc.STA, cat80[:K])[0]
                v.setdefault((reg, cell), []).append(L.err(L.s2_terms(nc.STA, cat80[:K], c0, E80[:K], par), nc.statics))
            table(P_, f"S2 K={K} {nm}", v)
    open(out, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__" and sys.argv[1] == "s2c":
    dev_s2c()


def dev_s2d(out="out/c11_dev_s2d.txt"):
    """fourth grid: joint Gauss-Newton stage (added after the noise-free validation showed that the alternation alone
    does not converge: 0.02-0.08 s left with an exact model)."""
    lines = []

    def P_(s=""):
        print(s, flush=True)
        lines.append(s)
    dev = load_dev()
    P_("C11 DEV S2, fourth grid: joint stage (same layout)")
    base = dict(L.P)
    grid = [("alternation only (n_joint=0)", {"n_joint": 0}), ("joint, floor 0.05", {}), ("joint, floor 0.02", {"sig_floor": 0.02}),
            ("joint, floor 0.02, cauchy 1.5", {"sig_floor": 0.02, "cauchy": 1.5}),
            ("joint, floor 0.02, cauchy 4", {"sig_floor": 0.02, "cauchy": 4.0}),
            ("joint, floor 0.02, sx=2", {"sig_floor": 0.02, "sx": 2.0}),
            ("joint, floor 0.02, min_arr=3", {"sig_floor": 0.02, "min_arr": 3}),
            ("joint, floor 0.02, gate_re=9", {"sig_floor": 0.02, "gate_re": 9.0}),
            ("joint, floor 0.02, gate_re=16", {"sig_floor": 0.02, "gate_re": 16.0}),
            ("joint, floor 0.02, w_re=0.5,0.3", {"sig_floor": 0.02, "w_re": [0.5, 0.3]}),
            ("joint, floor 0.02, n_pass=4", {"sig_floor": 0.02, "n_pass": 4})]
    for nm, kw in grid:
        par = dict(base)
        par.update(kw)
        for K in (20, 40, 80) if nm in ("joint, floor 0.02", "alternation only (n_joint=0)") else (20,):
            v, off = {}, []
            for (reg, cell, j), (nc, cat80, E80) in dev.items():
                c0 = nc.corr if K == 20 else xc.estimate_station_corrections(nc.STA, cat80[:K])[0]
                c = L.s2_terms(nc.STA, cat80[:K], c0, E80[:K], par)
                v.setdefault((reg, cell), []).append(L.err(c, nc.statics))
                off.append(abs(L.offset(c, nc.statics)))
            table(P_, f"S2 K={K} {nm}", v)
            P_(f"{'':34s} median |offset| {np.median(off):.3f} s")
    open(out, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__" and sys.argv[1] == "s2d":
    dev_s2d()


def dev_s3b(out="out/c11_dev_s3b.txt"):
    """second S3 grid (first read width of the atoms, self-consistency filter) and descriptive diagnostics of the
    atoms that enter the estimate (purity and location error use the truth; not used by any rule)."""
    lines = []

    def P_(s=""):
        print(s, flush=True)
        lines.append(s)
    dev = load_dev()
    ps = {}
    for reg in REGIMES:
        ps.update(pickle.load(open(f"out/c11_dev_s3pass_{reg}.pkl", "rb")))
    P_("C11 DEV S3, second grid (layout as c11_dev_s2.txt); eff = error against statics + mean model mismatch of the "
       "station (task2a_confirm.mismatch_offsets on the true catalogue sources)")
    base = dict(L.P)
    s2, v2, v2e, eff = {}, {}, {}, {}
    for key, (nc, cat80, E80) in dev.items():
        s2[key] = L.s2_terms(nc.STA, nc.cat, nc.corr, E80[:20])
        eff[key] = nc.statics + xc.T.mismatch_offsets(nc.STA, [s for s, _t in nc.cat_truth])
        v2.setdefault(key[:2], []).append(L.err(s2[key], nc.statics))
        v2e.setdefault(key[:2], []).append(L.err(s2[key], eff[key]))
    table(P_, "S2", v2)
    table(P_, "S2 (eff)", v2e)
    grid = [("default", {}), ("w_det0=1.0", {"w_det0": 1.0}), ("det_consist=1.5", {"det_consist": 1.5}),
            ("det_consist=2.5", {"det_consist": 2.5}), ("w_det0=1.0 det_consist=1.5", {"w_det0": 1.0, "det_consist": 1.5}),
            ("q_det=0.6 det_consist=1.5", {"q_det": 0.6, "det_consist": 1.5}),
            ("q_det=0.6 det_consist=1.5 n_unl=40", {"q_det": 0.6, "det_consist": 1.5, "n_unl": 40}),
            ("det_consist=1.5 sx_det=6 st_det=0.6", {"det_consist": 1.5, "sx_det": 6.0, "st_det": 0.6})]
    for nm, kw in grid:
        par = dict(base)
        par.update(kw)
        v, ve, nuse = {}, {}, []
        for key, (nc, cat80, E80) in dev.items():
            d = ps[key]
            n = par["n_unl"]
            info = {}
            c = L.s3_terms(nc.STA, nc.cat, nc.corr, s2[key], d["Ewin"][:n], d["atoms"][:n], E80[:20], par, info=info)
            v.setdefault(key[:2], []).append(L.err(c, nc.statics))
            ve.setdefault(key[:2], []).append(L.err(c, eff[key]))
            nuse.append(info["n_det_used"])
        table(P_, f"S3 {nm}", v)
        table(P_, f"S3 {nm} (eff)", ve)
        P_(f"{'':34s} atoms used per network {np.mean(nuse):.1f}")
    open(out, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__" and sys.argv[1] == "s3b":
    dev_s3b()
