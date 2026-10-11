"""C11 validation (network 195; TAGs 4104, 4105).  -> out/c11_validation.txt
V1 identity of S0 with the registered +SC and of the C11 engines / loops with the registered detectors.
V2 exact model, noise-free catalogue with catalogue errors: S2 recovers the true terms up to the gauge."""
import numpy as np

import c2_run as C2
import c8_lib as C8
import c9_run as R9
import c10_lib as C
import c10_screloc as SR
import c11_lib as L
import c11_run as RR
import xc
from station_corrections import peak_delays

lines = []


def P_(s=""):
    print(s, flush=True)
    lines.append(s)


P_("C11 VALIDATION (network 195)")
ok_read = ok_s0 = ok_idx = okm = okb = n = nw = 0
for cell in RR.CELLS:
    for reg in ("thin", "strong"):
        net = RR.Net(cell, 195, reg)
        nc = net.nc
        n += 1
        evs = L.cat_events(nc.cat)
        acc = [[[] for _ in range(2)] for _ in range(len(nc.STA))]
        same = True
        for ev, (D, src, tau) in zip(evs, nc.cat):
            dl, en = peak_delays(nc.STA, D, np.asarray(src, float), float(tau))
            tobs, en2 = L.read_peaks(ev.E, L.tt_hom(nc.STA, ev.x0) + ev.t0, 1.0)
            d2 = tobs - L.tt_hom(nc.STA, ev.x0) - ev.t0
            same &= bool(np.allclose(d2, dl, atol=1e-12) and np.allclose(en, en2))
            for i in range(len(nc.STA)):
                for b in range(2):
                    if en2[i, b] >= L.GATE:
                        acc[i][b].append(d2[i, b])
        c = np.array([[np.median(a) if a else 0.0 for a in row] for row in acc])
        ok_read += int(same)
        ok_s0 += int(np.allclose(c, nc.corr, atol=1e-12))
        E0, _pr, EF0, _prF = net.S.eng["cat"]
        ok_idx += int((L.clone_eng(E0, nc.corr).IDX == E0.IDX).all() and (L.clone_eng(EF0, nc.corr).IDX == EF0.IDX).all()
                      and np.allclose(L.clone_eng(E0, nc.corr).frac, E0.frac))
        # a zero-iteration JHD (no relocation possible: min_arr huge, one pass, one alternation) is the +SC median
        par = dict(L.P)
        par.update({"n_pass": 1, "n_alt": 1, "min_arr": 99, "n_joint": 0})
        cz = L.s2_terms(nc.STA, nc.cat, np.zeros_like(nc.corr), par=par)
        ok_s0 += 0
        P_(f"  {cell:15s} {reg:6s} read = peak_delays: {same}; median of gated reads = registered +SC: "
           f"{np.allclose(c, nc.corr, atol=1e-12)}; JHD without relocation = +SC: {np.allclose(cz, nc.corr, atol=1e-12)}")
        if reg == "strong":
            net.add("S0c", nc.corr, nc.cat)
            st = R9.stream(nc, 4105, C2.draw_dense, None)
            for wi in range(3):
                w = next(st)
                tr = w["truth"]
                nw += 1
                r0, l0 = C8.mtfw_general(net.S.dops, w["MFn"], tr, 7, stat="ls", capv=6.0)
                a = net.run("MT", "S0c", w["MFn"], 7, False)
                r1 = C.record(a, C.label(a, tr, w["taus"]))
                okm += int(np.allclose([x[0] for x in r0], [x[0] for x in r1], rtol=1e-10)
                           and [x[1] for x in r0] == [x[1] for x in r1]
                           and np.allclose([x[:3] for x in l0], [y["x"] for y in a]))
                r0, l0 = xc.bomp_pk(net.S.xops.PD, w["MFn"], tr, 7)
                a = net.run("BOMP", "S0c", w["MFn"], 7, False)
                r1 = C.record(a, C.label(a, tr, w["taus"]))
                okb += int(np.allclose([x[0] for x in r0], [x[0] for x in r1], rtol=1e-9)
                           and [x[1] for x in r0] == [x[1] for x in r1]
                           and np.allclose([x[:3] for x in l0], [y["x"] for y in a]))
                a3 = net.run("MT", "S0c", w["MFn"], 7, True)
                b3 = net.S.run("MT+T3", w["MFn"], 7)
                okm += 0
                assert np.allclose([y["x"] for y in a3], [y["x"] for y in b3]) and np.allclose(
                    [y["s"] for y in a3], [y["s"] for y in b3])
P_(f"V1 read identical to station_corrections.peak_delays {ok_read}/{n}; S0 = registered +SC {ok_s0}/{n}; cloned engine "
   f"tables = c10 engine tables {ok_idx}/{n}; MT loop on a cloned engine = c8_lib.mtfw_general ls cap 6: {okm}/{nw}; BOMP "
   f"loop = xc.bomp_pk: {okb}/{nw}; MT+T3 on a cloned engine = c10 MT+T3 (asserted) {nw}/{nw}")

P_("V2 exact model (straight rays VS 3.0, g 0, no direction or site scatter, no coda, 2 Hz zero-phase pulse, no clutter, "
   "noise 1e-4), statics N(0, 0.3 s) of the network, 20 catalogue events (TAG 4104).  rms error against the true statics "
   "(median removed) and offset (median of term - static), s:")
cfg = {"g": 0.0, "VS0": 3.0, "dir_deg": 0.0, "site_sd": 0.0, "coda_ratio": 0.0, "fc_range": (2.0, 2.0),
       "phase_max_deg": 0.0, "noise": 1e-4}
for cell in RR.CELLS:
    ci = RR.CELLS.index(cell)
    STA, statics = xc.T.network(L.SEED, cell, 195)
    for nm, ce in (("catalogue exact", (0.0, 0.0)), ("catalogue errors 1 km / 0.1 s", (1.0, 0.1))):
        cat, _d, truth = xc.catalogue_cfg(STA, 20, statics, [L.SEED, 4104, ci, 195, 0], 0, 0.0, cfg=cfg, cat_err=ce)
        c0 = xc.estimate_station_corrections(STA, cat)[0]
        s1 = SR.sc_reloc(STA, cat, c0, C.sigma_catalogue(STA, cat, c0)[0][0])
        info = {}
        evs = L.cat_events(cat)
        c2 = L.jhd(STA, evs, c0, info=info)
        dx = np.median([np.linalg.norm(ev.x - np.asarray(t[0])) for ev, t in zip(evs, truth)])
        dx0 = np.median([np.linalg.norm(ev.x0 - np.asarray(t[0])) for ev, t in zip(evs, truth)])
        par = dict(L.P)
        par.update({"sig_floor": 0.01})
        c2b = L.s2_terms(STA, cat, c0, par=par)
        par.update({"sx": 30.0, "st": 10.0})
        c2f = L.s2_terms(STA, cat, c0, par=par)
        par = dict(L.P)
        par.update({"n_joint": 0})
        c2a = L.s2_terms(STA, cat, c0, par=par)
        P_(f"  {cell:15s} {nm:30s} S0 {L.err(c0, statics):.4f} (off {L.offset(c0, statics):+.4f})  S1 "
           f"{L.err(s1, statics):.4f}  S2 {L.err(c2, statics):.4f} (off {L.offset(c2, statics):+.4f}; max abs "
           f"{np.abs(c2 - statics - np.median(c2 - statics)).max():.4f}; residual scale {info['sig']:.3f}; events "
           f"relocated {info['n_ev']}; median 3-D error of the catalogue {dx0:.2f} -> {dx:.2f} km)  S2 alternation only "
           f"{L.err(c2a, statics):.4f}  S2 scale floor 0.01 s {L.err(c2b, statics):.4f}  S2 floor 0.01 s and no prior (sx "
           f"30 km, st 10 s) {L.err(c2f, statics):.4f}")
open("out/c11_validation.txt", "w").write("\n".join(lines) + "\n")
