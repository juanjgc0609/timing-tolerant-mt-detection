"""C9 P4: fusion of travel-time-free evidence with BOMP-PK+SC and MT-FW+SC (cap 6).  SEED 999 ONLY.
Networks task2a_confirm.network(999, cell, j), j = 130-149 (never used before; 2A/2E 0-4, C1-C3 10-15, C6 20-29, leak2
30-35, C7 40-59, C8 60-89 and 99).  Window streams [999, TAG, ci, j, ri, 41|233|61]:
  TAG 3911 STD    xc.draw_sources (0-2 events), kmax 5      TAG 3912 DENSE   c2_run.draw_dense (3-5 events), kmax 7
  Oracle worlds reuse the DENSE / STD streams with a changed generator cfg (paired windows: same sources, mechanisms,
  pulses, clutter, noise): -POL0 (dir_deg 0), -AMP0 (site_sd 0), -BOTH0 (both 0); -POL20, -POL45 (dir_deg 20, 45; post hoc).  The +SC corrections are those of the
  nominal catalogue of the network in every world.
  TAG 3913 validation windows (network 199).
Per window the atoms of each loop are stored with their travel-time-free features; records are built in c9_analyze.py.
    python c9_run.py DESIGN REGIME j0,j1,... NWIN out.json        |      python c9_run.py validate
"""
import json
import os
import sys
import time

import numpy as np

import c2_run as C2
import c6_ceil as CZ
import c8_lib as C8
import c9_lib as C
import xc

CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
DESIGNS = {          # tag, kind, cfg, ceiling, full (BOMP in-loop variants)
    "STD": (3911, "std", None, True, True),
    "DENSE": (3912, "dense", None, True, True),
    "DENSE-POL0": (3912, "dense", {"dir_deg": 0.0}, False, False),
    "DENSE-AMP0": (3912, "dense", {"site_sd": 0.0}, False, False),
    "DENSE-BOTH0": (3912, "dense", {"dir_deg": 0.0, "site_sd": 0.0}, False, False),
    "STD-BOTH0": (3911, "std", {"dir_deg": 0.0, "site_sd": 0.0}, False, False),
    # added after a first look at DENSE thin (post hoc): dependence of the re-ranking on the polarization scatter
    "DENSE-POL20": (3912, "dense", {"dir_deg": 20.0}, False, False),
    "DENSE-POL45": (3912, "dense", {"dir_deg": 45.0}, False, False),
}


def stream(nc, tag, draw, cfg):
    b = [C.SEED, tag, nc.ci, nc.j, nc.ri]
    rng, rng_m, rng_h = (np.random.default_rng(b + [41]), np.random.default_rng(b + [233]),
                         np.random.default_rng(b + [61]))
    while True:
        srcs, taus = draw(rng)
        D, picks, site, _amp = C8.generate_site(nc.STA, srcs, taus, nc.nclut, rng, rng_m, rng_h, nc.clut_amp, cfg=cfg,
                                                statics=nc.statics)
        MF, sd = xc.matched(D)
        yield {"D": D, "MFn": MF / sd, "sd": sd, "truth": [np.asarray(s, float) for s in srcs],
               "taus": [float(t) for t in taus], "picks": picks, "site": site}


def pack(atoms, feats=None):
    out = []
    for k, a in enumerate(atoms):
        d = {"n": a["n"], "n0": a["n0"], "t": a["t"], "s": a["s"]}
        if feats is not None:
            d["f"] = feats[k]
        out.append(d)
    return out


def window(xops, dops, geo, w, kmax, full):
    MFn = w["MFn"]
    _Gw, Qw = dops.scaled(np.asarray(w["site"]).T)
    res, sec = {}, {}
    t = time.time()
    a = C.bomp_c9(xops.PD, geo, MFn, kmax)
    sec["BOMP-PK"] = time.time() - t
    res["BOMP-PK"] = pack(a, [C.atom_features(geo, MFn, x["n"], x["I"], x["v"], Qw) for x in a])
    t = time.time()
    a = C.mtfw_c9(dops, geo, MFn, kmax)
    sec["MT-FW"] = time.time() - t
    t = time.time()
    res["MT-FW"] = pack(a, [C.atom_features(geo, MFn, x["n"], x["I"], x["v"], Qw) for x in a])
    sec["features"] = (time.time() - t) / max(1, len(a))
    for name, kw in (("MT-FW+RLIN-mt", dict(reloc="mt")), ("MT-FW+RLIN-pol", dict(reloc="pol")),
                     ("MT-FW+PEN", dict(pen=C.PEN_LAM))):
        t = time.time()
        res[name] = pack(C.mtfw_c9(dops, geo, MFn, kmax, **kw))
        sec[name] = time.time() - t
    if full:
        t = time.time()
        res["BOMP-PK+RLIN-mt"] = pack(C.bomp_c9(xops.PD, geo, MFn, kmax, reloc="mt"))
        sec["BOMP-PK+RLIN-mt"] = time.time() - t
    return res, sec


def run(design, regime, nets, nwin, path):
    tag, kind, cfg, ceil, full = DESIGNS[design]
    draw = C2.draw_dense if kind == "dense" else xc.draw_sources
    kmax = 7 if kind == "dense" else 5
    res = json.load(open(path)) if os.path.exists(path) else {"rows": [], "nets": []}
    done = {(n["cell"], n["net"]) for n in res["nets"]}
    grid = xc.make_grid()
    garr = np.asarray(grid)
    for cell in CELLS:
        for j in nets:
            if (cell, j) in done:
                continue
            t0 = time.time()
            nc = xc.NetCfg(cell, j, regime)
            xops = xc.Ops(nc.STA, nc.corr, grid)
            dops = C8.DOps(xops)
            geo = C.Geo(garr, nc.STA)
            st = stream(nc, tag, draw, cfg)
            for wi in range(nwin):
                w = next(st)
                tr = w["truth"]
                atoms, sec = window(xops, dops, geo, w, kmax, full)
                row = {"cell": cell, "net": j, "w": wi, "k": len(tr), "truth": [t.tolist() for t in tr],
                       "taus": w["taus"], "atoms": atoms, "sec": sec}
                if ceil:
                    t = time.time()
                    R_or = CZ.oracle_residual_refined(xops.net, w, tr, w["taus"])
                    row["ceil"] = CZ.record(*CZ.ceil_scores(xops.PD, garr, w["MFn"], R_or, tr))
                    sec["CEIL-GD"] = time.time() - t
                res["rows"].append(row)
            res["nets"].append({"cell": cell, "net": j, "regime": regime, "design": design, "t_total": time.time() - t0})
            with open(path + ".tmp", "w") as fh:
                json.dump(res, fh)
            os.replace(path + ".tmp", path)
            print(design, regime, cell, j, "done", time.strftime("%H:%M:%S"), f"{time.time() - t0:.0f} s", flush=True)


def validate():
    """bomp_c9 / mtfw_c9 without relocation reproduce xc.bomp_pk and c8_lib.mtfw_general (ls, cap 6): hits, locations,
    scores; network 199, dense windows of TAG 3913."""
    grid = xc.make_grid()
    garr = np.asarray(grid)
    okb = okf = n = 0
    for cell in CELLS:
        nc = xc.NetCfg(cell, 199, "strong")
        xops = xc.Ops(nc.STA, nc.corr, grid)
        dops = C8.DOps(xops)
        geo = C.Geo(garr, nc.STA)
        st = stream(nc, 3913, C2.draw_dense, None)
        for _ in range(3):
            w = next(st)
            tr = w["truth"]
            n += 1
            r0, l0 = xc.bomp_pk(xops.PD, w["MFn"], tr, 7)
            a = C.bomp_c9(xops.PD, geo, w["MFn"], 7)
            r1 = C.record(a, garr, tr)
            okb += int(len(r0) == len(r1) and np.allclose([x[0] for x in r0], [x[0] for x in r1], rtol=1e-10)
                       and [x[1] for x in r0] == [x[1] for x in r1]
                       and np.allclose([x[:3] for x in l0], [garr[y["n"]] for y in a]))
            r0, l0 = C8.mtfw_general(dops, w["MFn"], tr, 7, stat="ls", capv=6.0)
            a = C.mtfw_c9(dops, geo, w["MFn"], 7)
            r1 = C.record(a, garr, tr)
            okf += int(len(r0) == len(r1) and np.allclose([x[0] for x in r0], [x[0] for x in r1], rtol=1e-10)
                       and [x[1] for x in r0] == [x[1] for x in r1]
                       and np.allclose([x[:3] for x in l0], [garr[y["n"]] for y in a]))
            # vectors stored at the read samples: the MT statistic recomputed from them equals the score
            a0 = a[0]
            s = C.st_mt(geo, C.cap(a0["v"]), np.array([a0["n"]]), Q=dops.Qst)[0] / 2.0
            assert abs(s - a0["s"]) < 1e-8 * max(1.0, a0["s"]), (s, a0["s"])
    print(f"V3 bomp_c9 = xc.bomp_pk (scores, hits, locations): {okb}/{n} windows")
    print(f"V4 mtfw_c9 = c8_lib.mtfw_general ls cap 6 (scores, hits, locations): {okf}/{n} windows; first-atom score "
          f"recomputed from the stored vectors with the c9 basis: equal")


if __name__ == "__main__":
    if sys.argv[1] == "validate":
        validate()
    else:
        run(sys.argv[1], sys.argv[2], [int(x) for x in sys.argv[3].split(",")], int(sys.argv[4]), sys.argv[5])
