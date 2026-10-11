"""C9 P1-P3: standalone travel-time-free location of single events (thin clutter).  SEED 999 ONLY.
Networks task2a_confirm.network(999, cell, j), j = 130-149; event streams [999, 3901, ci, j, ev, 41|233|61] (new).
Each event is regenerated with identical streams under each generator setting (dir_deg, site_sd), so settings are
paired (same source, mechanism, pulse, clutter, noise; the deviation and site draws are the same normal variates
scaled).  Location on a 1 km grid (|x|, |y| <= 21 km, depths 5, 9, 13 km), epicentral error of the maximum.
    python c9_standalone.py j0,j1,... NEV out.json      |      python c9_standalone.py validate
"""
import json
import sys
import time

import numpy as np

import c9_lib as C
import c8_lib as C8
import task2a_confirm as T
import xc
from task1y_clutter import CLUT_AMP

TAG = 3901
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
SETTINGS = [(0.0, 0.5), (5.0, 0.5), (10.0, 0.5), (20.0, 0.5), (10.0, 0.0), (10.0, 0.25), (0.0, 0.0)]
AX = np.arange(-21.0, 21.0 + 1e-9, 1.0)
FINE = np.array([[x, y, z] for x in AX for y in AX for z in (5.0, 9.0, 13.0)])


def estimators(geo, STA, v, K, statics):
    """epicentres of every travel-time-free estimator (and the arrival-time one) from the vectors v and samples K."""
    out = {}
    out["P1-azim"] = C.azimuth_intersection(STA, v[0]).tolist()
    for name, s in (("P1-grid", C.st_ppol_angle(geo, v)), ("P1-grid-gradmodel", C.st_ppol_angle(geo, v, grad=True)),
                    ("P3-Strans", C.st_pol(geo, v, parts=True)[1]), ("P13-pol", C.st_pol(geo, v)),
                    ("P3-gainfree", C.st_gfree(geo, v)), ("P2-mt", C.st_mt(geo, v)),
                    ("P2-mt-cap6", C.st_mt(geo, C.cap(v))), ("P2b-amp", C.st_amp(geo, v)),
                    ("TIME-raw", C.time_locate(geo, K)), ("TIME-truestatics", C.time_locate(geo, K, statics))):
        out[name] = geo.XN[int(np.argmax(s))][:2].tolist()
    return out


def one_event(STA, statics, nclut, clut_amp, seeds, cfg, geo, src=None, tau=None):
    rng, rng_m, rng_h = (np.random.default_rng(seeds + [41]), np.random.default_rng(seeds + [233]),
                         np.random.default_rng(seeds + [61]))
    s0 = np.array([rng.uniform(-15, 15), rng.uniform(-15, 15), rng.uniform(5, 13)])
    t0 = float(rng.uniform(5.5, 7.5))
    src = s0 if src is None else src
    tau = t0 if tau is None else tau
    D, picks, site, _amp = C8.generate_site(STA, [src], [tau], nclut, rng, rng_m, rng_h, clut_amp, cfg=cfg,
                                            statics=statics)
    MF, sd = xc.matched(D)
    MFn = MF / sd
    ns = len(STA)
    tarr = np.array([[picks[i][b][0] for i in range(ns)] for b in range(2)])
    vt, Kt = C.read_true(MFn, tarr)
    vb, Kb = C.read_blind(MFn)
    res = {"src": src.tolist(), "ampP": np.linalg.norm(vt[0], axis=-1).tolist(),
           "ampS": np.linalg.norm(vt[1], axis=-1).tolist(),
           "blind_ok": (np.abs(Kb - Kt) <= 3).mean(axis=1).tolist(),
           "true": estimators(geo, STA, vt, Kt, statics), "blind": estimators(geo, STA, vb, Kb, statics)}
    return res


def run(nets, nev, path):
    rows = []
    t0 = time.time()
    for cell in CELLS:
        ci, _k, nst, nclut = T.net_setup(C.SEED, cell)
        for j in nets:
            STA, statics = T.network(C.SEED, cell, j)
            geo = C.Geo(FINE, STA, grad=True, need_g=True)
            for ev in range(nev):
                for dd, ss in SETTINGS:
                    r = one_event(STA, statics, nclut, CLUT_AMP["thin"], [C.SEED, TAG, ci, j, ev],
                                  {"dir_deg": dd, "site_sd": ss}, geo)
                    r.update(cell=cell, net=j, ev=ev, dir_deg=dd, site_sd=ss)
                    rows.append(r)
            print(cell, j, f"{time.time() - t0:.0f} s", flush=True)
            json.dump(rows, open(path, "w"))


def validate():
    """noise-free recovery: straight rays (g = 0, VS0 = 3.0), no deviation, no site factor, no coda, 2 Hz zero-phase
    pulse, no clutter, noise 1e-4, no statics; network j = 199, sources ON fine-grid nodes."""
    cfg = {"g": 0.0, "VS0": 3.0, "dir_deg": 0.0, "site_sd": 0.0, "coda_ratio": 0.0, "fc_range": (2.0, 2.0),
           "phase_max_deg": 0.0, "noise": 1e-4}
    rng = np.random.default_rng([C.SEED, 3902])
    for cell in CELLS:
        ci, _k, nst, _nc = T.net_setup(C.SEED, cell)
        STA, _st = T.network(C.SEED, cell, 199)
        geo = C.Geo(FINE, STA, grad=True, need_g=True)
        errs = {}
        for ev in range(12):
            src = np.array([float(rng.integers(-15, 16)), float(rng.integers(-15, 16)), 9.0])
            r = one_event(STA, np.zeros((nst, 2)), 0, 0.0, [C.SEED, 3902, ci, 199, ev], cfg, geo, src=src, tau=6.0)
            for mode in ("true", "blind"):
                for k, xy in r[mode].items():
                    errs.setdefault((mode, k), []).append(float(np.hypot(xy[0] - src[0], xy[1] - src[1])))
        print(cell, "noise-free recovery, max epicentral error over 12 events (km):")
        for (mode, k), e in sorted(errs.items()):
            print(f"   {mode:5s} {k:20s} max {np.nanmax(e):6.2f}  median {np.nanmedian(e):6.2f}")


if __name__ == "__main__":
    if sys.argv[1] == "validate":
        validate()
    else:
        run([int(x) for x in sys.argv[1].split(",")], int(sys.argv[2]), sys.argv[3])
