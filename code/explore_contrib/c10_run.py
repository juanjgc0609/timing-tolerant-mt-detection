"""C10: soft timing tolerance inside BOMP-PK+SC and MT-FW+SC (cap 6).  SEED 999 ONLY.
Networks task2a_confirm.network(999, cell, j), j = 160-179 (never used before: 2A/2E 0-4, C1-C3 10-15, C6 20-29, leak2
30-35, C7 40-59, C8 60-89 and 99, C9 130-149 and 199); j = 195: validation and timing only.
Window streams [999, TAG, ci, j, ri, 41|233|61], all new: 4001 STD (xc.draw_sources, 0-2 events, kmax 5), 4002 DENSE
(c2_run.draw_dense, 3-5 events, kmax 7), 4003 validation windows (network 195), 4004 noise-free validation (network 195).
    python c10_run.py DESIGN REGIME j0,j1,... NWIN out.json     |     python c10_run.py validate
"""
import json
import os
import sys
import time

import numpy as np

import c2_run as C2
import c6_ceil as CZ
import c8_lib as C8
import c9_run as R9
import c10_lib as C
import xc

CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
DESIGNS = {"STD": (4001, "std"), "DENSE": (4002, "dense")}
GRID = xc.make_grid()
GARR = np.asarray(GRID)
_ax = np.arange(-21.0, 21.0 + 1e-9, C.H_FINE) + 0.37
FINE = np.array([[x, y, z] for x in _ax for y in _ax for z in (5.0, 9.0, 13.0)])
TAUS = np.asarray(xc.DEFAULT_TAUS, float)

# name: (family, engine set, mode, W, refine)
VARIANTS = {
    "MT-FW": ("mt", "cat", "hard", 4, False),
    "MT+T1": ("mt", "cat", "indep", 4, False),
    "MT+T1c": ("mt", "cat", "indepc", 4, False),
    "MT+T2": ("mt", "cat", "corr", 4, False),
    "MT+T3": ("mt", "cat", "hard", 4, True),
    "MT+T4": ("mt", "cat", "indep", 8, False),
    "MT+T4h": ("mt", "cat", "hard", 8, False),
    "MT+T2w": ("mt", "cat", "corr", 8, False),
    "MT-ORA": ("mt", "ora", "hard", 4, False),
    "MT-ORA+T2": ("mt", "ora", "corr", 4, False),
    "MT-ORA+T2w": ("mt", "ora", "corr", 8, False),
    "MT-ORA+T3": ("mt", "ora", "hard", 4, True),
    "BOMP-PK": ("bomp", "cat", "hard", 4, False),
    "BOMP+T1": ("bomp", "cat", "indep", 4, False),
    "BOMP+T1c": ("bomp", "cat", "indepc", 4, False),
    "BOMP+T2": ("bomp", "cat", "corr", 4, False),
    "BOMP+T3": ("bomp", "cat", "hard", 4, True),
    "BOMP+T4": ("bomp", "cat", "indep", 8, False),
    "BOMP+T2w": ("bomp", "cat", "corr", 8, False),
    "BOMP-ORA": ("bomp", "ora", "hard", 4, False),
    "BOMP-ORA+T2": ("bomp", "ora", "corr", 4, False),
    "BOMP-ORA+T3": ("bomp", "ora", "hard", 4, True),
}


# main run set (trimmed for the compute budget after the timing of the validation run; the other variants were
# only exercised in the validation)
RUNSET = ["MT-FW", "MT+T1", "MT+T1c", "MT+T2", "MT+T3", "MT+T4", "MT+T4h", "MT+T2w", "MT-ORA", "MT-ORA+T2w",
          "MT-ORA+T3", "BOMP-PK", "BOMP+T1", "BOMP+T2w", "BOMP+T3"]


class Setup:
    def __init__(self, cell, j, regime, oracle=True, corr=None, statics=None, sig=None, cfg_cat=None):
        t0 = time.time()
        nc = xc.NetCfg(cell, j, regime)
        if statics is not None:
            nc.statics = statics
        if corr is not None:
            nc.corr = corr
        self.nc = nc
        self.xops = xc.Ops(nc.STA, nc.corr, GRID)
        self.dops = C8.DOps(self.xops)
        self.sig, self.cnt = C.sigma_catalogue(nc.STA, nc.cat, nc.corr) if sig is None else (np.asarray(sig), None)
        E = C.Eng(GARR, nc.STA, nc.corr, TAUS, C.H_COARSE, G=self.dops.Gs, Q=self.dops.Qst, bomp=True)
        assert (E.IDX == self.dops.IDX).all() and (E.IDX == self.xops.BD.idx).all()
        EF = C.Eng(FINE, nc.STA, nc.corr, TAUS, C.H_FINE, bomp=True)
        self.eng = {"cat": (E, E.prior(self.sig), EF, EF.prior(self.sig))}
        if oracle:
            self.sig_ora = C.sigma_oracle(nc.STA, nc.cat_truth)
            EO = C.Eng(GARR, nc.STA, nc.statics, TAUS, C.H_COARSE, G=E.G, Q=E.Q, bomp=True)
            EFO = C.Eng(FINE, nc.STA, nc.statics, TAUS, C.H_FINE, G=EF.G, Q=EF.Q, bomp=True)
            self.eng["ora"] = (EO, EO.prior(self.sig_ora), EFO, EFO.prior(self.sig_ora))
        self.t_setup = time.time() - t0

    def run(self, name, MFn, kmax, capv=C.CAPV):
        fam, es, mode, W, ref = VARIANTS[name]
        E, pr, EF, prF = self.eng[es]
        kw = dict(mode=mode, pr=pr, W=W, EF=EF if ref else None, prF=prF)
        if fam == "mt":
            return C.mt_loop(E, MFn, kmax, capv=capv, **kw)
        return C.bomp_loop(E, self.xops.PD, MFn, kmax, **kw)


def pack(atoms):
    return [[float(a["x"][0]), float(a["x"][1]), float(a["x"][2]), float(a["tau"]), float(a["s"]), float(a["s0"]),
             float(a.get("sref", np.nan))] for a in atoms]


def run(design, regime, nets, nwin, path, names=None):
    tag, kind = DESIGNS[design]
    draw = C2.draw_dense if kind == "dense" else xc.draw_sources
    kmax = 7 if kind == "dense" else 5
    names = list(RUNSET) if names is None else names
    res = json.load(open(path)) if os.path.exists(path) else {"rows": [], "nets": []}
    done = {(n["cell"], n["net"]) for n in res["nets"]}
    for cell in CELLS:
        for j in nets:
            if (cell, j) in done:
                continue
            t0 = time.time()
            S = Setup(cell, j, regime)
            st = R9.stream(S.nc, tag, draw, None)
            for wi in range(nwin):
                w = next(st)
                tr = w["truth"]
                atoms, sec = {}, {}
                for nm in names:
                    t = time.time()
                    atoms[nm] = pack(S.run(nm, w["MFn"], kmax))
                    sec[nm] = time.time() - t
                t = time.time()
                R_or = CZ.oracle_residual_refined(S.xops.net, w, tr, w["taus"])
                ceil = CZ.record(*CZ.ceil_scores(S.xops.PD, GARR, w["MFn"], R_or, tr))
                sec["CEIL-GD"] = time.time() - t
                res["rows"].append({"cell": cell, "net": j, "w": wi, "k": len(tr), "truth": [t_.tolist() for t_ in tr],
                                    "taus": w["taus"], "atoms": atoms, "ceil": ceil, "sec": sec})
            res["nets"].append({"cell": cell, "net": j, "regime": regime, "design": design, "t_total": time.time() - t0,
                                "t_setup": S.t_setup, "sig": S.sig.tolist(), "sig_ora": S.sig_ora.tolist(),
                                "cnt": S.cnt.tolist(),
                                "cell_sd": np.sqrt(S.eng["cat"][0].cellvar.mean((0, 2))).tolist(),
                                "cell_sd_fine": np.sqrt(S.eng["cat"][2].cellvar.mean((0, 2))).tolist()})
            with open(path + ".tmp", "w") as fh:
                json.dump(res, fh)
            os.replace(path + ".tmp", path)
            print(design, regime, cell, j, "done", time.strftime("%H:%M:%S"), f"{time.time() - t0:.0f} s",
                  "sig", np.round(S.sig, 3), "ora", np.round(S.sig_ora, 3), flush=True)


# ============================================================================================ validation
def validate():
    P = print
    P("C10 VALIDATION (network 195; TAGs 4003, 4004)")
    okb = okf = n = 0
    de = 0.0
    tsec = {}
    for cell in CELLS:
        S = Setup(cell, 195, "strong")
        P(f"{cell}: setup {S.t_setup:.1f} s; sig_cat (P, S) {np.round(S.sig, 3)} from {S.cnt} residuals; sig_ora "
          f"{np.round(S.sig_ora, 3)}; coarse cell sd (P, S) "
          f"{np.round(np.sqrt(S.eng['cat'][0].cellvar.mean((0, 2))), 3)}; fine "
          f"{np.round(np.sqrt(S.eng['cat'][2].cellvar.mean((0, 2))), 3)}")
        st = R9.stream(S.nc, 4003, C2.draw_dense, None)
        E = S.eng["cat"][0]
        for wi in range(3):
            w = next(st)
            tr = w["truth"]
            n += 1
            r0, l0 = xc.bomp_pk(S.xops.PD, w["MFn"], tr, 7)
            a = S.run("BOMP-PK", w["MFn"], 7)
            r1 = C.record(a, C.label(a, tr, w["taus"]))
            okb += int(len(r0) == len(r1) and np.allclose([x[0] for x in r0], [x[0] for x in r1], rtol=1e-9)
                       and [x[1] for x in r0] == [x[1] for x in r1]
                       and np.allclose([x[:3] for x in l0], [y["x"] for y in a])
                       and np.allclose([x[3] for x in l0], [y["tau"] for y in a]))
            r0, l0 = C8.mtfw_general(S.dops, w["MFn"], tr, 7, stat="ls", capv=6.0)
            a = S.run("MT-FW", w["MFn"], 7)
            r1 = C.record(a, C.label(a, tr, w["taus"]))
            okf += int(len(r0) == len(r1) and np.allclose([x[0] for x in r0], [x[0] for x in r1], rtol=1e-10)
                       and [x[1] for x in r0] == [x[1] for x in r1]
                       and np.allclose([x[:3] for x in l0], [y["x"] for y in a]))
            e0, I0, _h = S.xops.PD.energy(w["MFn"], True)
            tpk, _p = C.choose(E, None, "hard", None, 4, R=w["MFn"])
            assert (tpk == I0).all()
            de = max(de, float(np.abs(C.bomp_energy(E, w["MFn"], tpk) - e0).max() / e0.max()))
            for nm in VARIANTS:
                t = time.time()
                S.run(nm, w["MFn"], 7)
                tsec.setdefault((cell, nm), []).append(time.time() - t)
        P(f"  V1 bomp_loop(hard) = xc.bomp_pk: {okb}/{n}; mt_loop(hard) = c8_lib.mtfw_general ls cap 6: {okf}/{n}; "
          f"bomp_energy vs PKDict.energy max rel diff {de:.1e}")
    P("  seconds per dense window (3 windows per cell): ")
    for nm in VARIANTS:
        P(f"    {nm:14s} " + " | ".join(f"{np.mean(tsec[(c, nm)]):.2f}" for c in CELLS))
    # ---------------------------------------------------------------- noise-free, exact statics
    P("V2 exact model and exact statics (straight rays VS0 3.0 g 0, no deviation, no site factors, no coda, 2 Hz "
      "zero-phase pulse, no clutter, statics 0, corrections 0, sig = 0.05 s fixed), single events, kmax 1;")
    P("   noise 1e-4 (noise-free; MT family without cap) and noise 0.02 (nominal waveform noise; cap 6)")
    for noise, capv in ((1e-4, None), (0.02, C.CAPV)):
        cfg = {"g": 0.0, "VS0": 3.0, "dir_deg": 0.0, "site_sd": 0.0, "coda_ratio": 0.0, "fc_range": (2.0, 2.0),
               "phase_max_deg": 0.0, "noise": noise}
        rng = np.random.default_rng([C.SEED, 4004, 0])
        for ci, cell in enumerate(CELLS):
            nst = len(xc.T.network(C.SEED, cell, 195)[0])
            z = np.zeros((nst, 2))
            S = Setup(cell, 195, "strong", oracle=False, corr=z, statics=z, sig=(0.05, 0.05))
            E, pr, EF, prF = S.eng["cat"]
            for kind in ("on-node", "off-node"):
                err = {nm: [] for nm in VARIANTS if VARIANTS[nm][1] == "cat"}
                pen = {nm: [] for nm in err}
                terr, therr = [], []
                for ev in range(10):
                    if kind == "on-node":
                        src = GARR[rng.integers(len(GARR))].copy()
                        while np.abs(src[:2]).max() > 15:
                            src = GARR[rng.integers(len(GARR))].copy()
                        tau = float(TAUS[rng.integers(5, 26)])
                    else:
                        src = np.array([rng.uniform(-15, 15), rng.uniform(-15, 15), rng.uniform(5, 13)])
                        tau = float(rng.uniform(5.5, 7.5))
                    b = [C.SEED, 4004, ci, 195, ev, 0 if kind == "on-node" else 1]
                    D, _pk, _site, _amp = C8.generate_site(S.nc.STA, [src], [tau], 0, np.random.default_rng(b + [41]),
                                                           np.random.default_rng(b + [233]),
                                                           np.random.default_rng(b + [61]), 0.0, cfg=cfg, statics=z)
                    MF, sd = xc.matched(D)
                    MFn = MF / sd
                    for nm in err:
                        a = S.run(nm, MFn, 1, capv=capv)[0]
                        err[nm].append(float(np.linalg.norm(a["x"][:2] - src[:2])))
                        pen[nm].append(a["s0"] - a["s"])
                    # timing-only locator: minimum over the 1 km lattice of the charge of the hard-rule peaks
                    Eh = C._half_energy(MFn, None)
                    tph, _p = C.choose(EF, Eh, "hard", None, 4, R=MFn)
                    dlh = (tph - EF.IDX - EF.frac) * C.DT
                    pf = (dlh * dlh / (2.0 * prF["sc2"][None, None])).sum((-1, -2))
                    kf = int(np.argmin(pf.min(axis=1)))
                    terr.append(float(np.linalg.norm(EF.XN[kf, :2] - src[:2])))
                    # theta of the corr rule at the nearest coarse node: recovered offset
                    kn = int(np.argmin(np.linalg.norm(GARR - src[None], axis=1)))
                    tn = int(np.argmin(np.abs(TAUS - tau)))
                    _t, dl = C.pick(Eh, E.IDX[kn:kn + 1, tn:tn + 1], E.frac[kn:kn + 1, tn:tn + 1],
                                    pr["s2"][kn:kn + 1, None], 8)
                    th = pr["Tm"][kn] @ dl.reshape(-1)
                    therr.append(float(np.linalg.norm(GARR[kn, :2] + th[:2] - src[:2])))
                P(f"  noise {noise:g} {cell:15s} {kind:8s} epicentral error km (median / max), median timing charge of "
                  "the atom:")
                P("      " + "; ".join(f"{nm} {np.median(err[nm]):.2f}/{np.max(err[nm]):.2f} c {np.median(pen[nm]):.2f}"
                                       for nm in err))
                P(f"      timing-only minimum of the charge of the hard-rule peaks on the 1 km lattice: "
                  f"{np.median(terr):.2f} / {np.max(terr):.2f} km; nearest node + theta (corr prior, +-0.5 s): "
                  f"{np.median(therr):.2f} / {np.max(therr):.2f} km")


if __name__ == "__main__":
    if sys.argv[1] == "validate":
        validate()
    else:
        run(sys.argv[1], sys.argv[2], [int(x) for x in sys.argv[3].split(",")], int(sys.argv[4]), sys.argv[5],
            sys.argv[6].split(",") if len(sys.argv) > 6 else None)
