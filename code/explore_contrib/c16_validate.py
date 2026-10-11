"""C16 validation of the injection code (development network 2400 of each cell).  SEED 999 ONLY.
  V1 GEN-arm catalogue = xc.catalogue_cfg (standard path) bit for bit; +SC terms identical.
  V2 GEN-arm STD / DENSE / noise-only windows = c8_lib.generate_site (standard path) bit for bit; detector atoms of
     every method identical on the standard-path traces; MT@S0 = c8_lib.mtfw_general (ls, cap 6) atom for atom.
  V3 REAL arm: same sources, origin times and signal traces as GEN (pairing); catalogue candidates identical.
  V4 noise slots: eligible, inside the network's block, disjoint between the uses of a network and between the networks
     of a cell; the injected segment is the stored day array.
    python c16_validate.py > out/c16_validation.txt"""
import numpy as np

import c2_run as C2
import c8_lib as C8
import c10_lib as C
import c16_lib as L16
import xc


def main():
    P = print
    _d, meta = L16.noise_store()
    ev = np.array(meta["events_s"])
    ok_all = True
    for cell in L16.CELLS:
        j = L16.NETS_DEV[0]
        g = L16.Net16(cell, j, "GEN")
        ci = g.ci
        cat_s, _dr, tr_s = xc.catalogue_cfg(g.STA, 20, g.statics, [L16.SEED, L16.TAG_CAT, ci, j, L16.RI], g.nclut,
                                            g.clut_amp)
        v1 = all(np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1]) and a[2] == b[2] for a, b in zip(g.cat, cat_s))
        v1 &= np.array_equal(g.corr, xc.estimate_station_corrections(g.STA, cat_s)[0])
        P(f"{cell} net {j}: stations {g.sub.tolist()}; V1 catalogue (20 events, {g.drawn} drawn) and +SC identical to "
          f"xc.catalogue_cfg: {v1}")
        g.setup()
        nwd, nat, nmt = 0, 0, 0
        for design, tag, draw, kmax in (("STD", L16.TAG_STD, xc.draw_sources, 5), ("DENSE", L16.TAG_DENSE, C2.draw_dense, 7)):
            b = [L16.SEED, tag, ci, j, L16.RI]
            rng, rng_m, rng_h = (np.random.default_rng(b + [41]), np.random.default_rng(b + [233]),
                                 np.random.default_rng(b + [61]))
            st = g.stream(design)
            for wi in range(4):
                srcs, taus = draw(rng)
                D0, _p, _s, _a = C8.generate_site(g.STA, srcs, taus, g.nclut, rng, rng_m, rng_h, g.clut_amp,
                                                  statics=g.statics)
                w = next(st)
                same = np.array_equal(D0, w["D"])
                nwd += int(same)
                MF, sd = xc.matched(D0)
                MFn = MF / sd
                okm = True
                for nm in L16.METHODS:
                    a0 = L16.pack(g.run(nm, MFn, kmax))
                    a1 = L16.pack(g.run(nm, w["MFn"], kmax))
                    okm &= np.array_equal(np.array(a0), np.array(a1), equal_nan=True)
                nat += int(okm)
                rec, locs = C8.mtfw_general(g.dops, MFn, [], kmax, stat="ls", capv=6.0)
                a = g.run("MT@S0", MFn, kmax)
                nmt += int(len(a) == len(locs) and np.allclose([x[:3] for x in locs], [y["x"] for y in a])
                           and np.allclose([x[3] for x in locs], [y["tau"] for y in a])
                           and np.allclose([r_[0] for r_ in rec], np.minimum.accumulate([y["s"] for y in a])))
        rng = np.random.default_rng([L16.SEED, L16.TAG_NOISE0, ci, j, L16.RI, 41])
        st = g.noise_windows()
        n0 = 0
        for wi in range(4):
            D0, _p, _s, _a = C8.generate_site(g.STA, [], [], g.nclut, rng, np.random.default_rng(1),
                                              np.random.default_rng(2), g.clut_amp, statics=g.statics)
            n0 += int(np.array_equal(D0, next(st)["D"]))
        v2 = nwd == 8 and nat == 8 and nmt == 8 and n0 == 4
        P(f"  V2 windows bit-identical to c8_lib.generate_site {nwd}/8, noise-only {n0}/4; atoms of all "
          f"{len(L16.METHODS)} methods identical {nat}/8; MT@S0 = c8_lib.mtfw_general ls cap 6 {nmt}/8")
        r = L16.Net16(cell, j, "REAL")
        cand_g, cand_r = [], []
        v3 = True
        for design in ("STD", "DENSE"):
            sg, sr = g.stream(design), r.stream(design)
            for wi in range(6):
                a, b_ = next(sg), next(sr)
                v3 &= (len(a["truth"]) == len(b_["truth"]) and all(np.array_equal(x, y) for x, y in zip(a["truth"], b_["truth"]))
                       and a["taus"] == b_["taus"] and np.array_equal(a["Ds"], b_["Ds"]))
                seg = L16.real_segment(r.sub, b_["slot"])
                v3 &= np.array_equal(b_["N"], seg) and np.allclose(b_["D"], b_["Ds"] + seg, rtol=0, atol=1e-12)
        # catalogue candidate sequence: the accepted truth of each arm is a subsequence of the same candidate list
        seeds = [L16.SEED, L16.TAG_CAT, ci, j, L16.RI]
        rng = np.random.default_rng(list(seeds) + [0])
        cands = []
        for _ in range(max(g.drawn, r.drawn)):
            src = np.array([rng.uniform(-15, 15), rng.uniform(-15, 15), rng.uniform(5, 13)])
            tau = float(rng.uniform(5.5, 7.5))
            C8.generate_site(g.STA, [src], [tau], g.nclut, rng, np.random.default_rng(1), np.random.default_rng(2), 0.0)
            cands.append((src, tau))

        def subseq(tr):
            k = 0
            for s_, t_ in tr:
                while k < len(cands) and not (np.array_equal(cands[k][0], s_) and cands[k][1] == t_):
                    k += 1
                if k == len(cands):
                    return False
                k += 1
            return True
        v3 &= subseq(g.cat_truth) and subseq(r.cat_truth)
        P(f"  V3 REAL arm: sources, origins and signal traces identical to GEN in 12/12 windows, injected noise = stored "
          f"segment, catalogue truth a subsequence of the common candidate list (GEN {g.drawn} drawn, REAL {r.drawn}): {v3}")
        ok_all &= v1 and v2 and v3
    # V4 slots
    el = set(meta["eligible"])
    v4 = True
    allb = []
    for j in L16.NETS_DEV + L16.NETS_EVAL:
        sl = L16.slots(j)
        u = sl["cat"] + sl["STD"] + sl["DENSE"] + sl["NOISE"]
        v4 &= len(set(u)) == len(u) and set(u) <= el
        for k in u:
            a, b_ = k * 40.0, k * 40.0 + 40.0
            v4 &= not np.any((ev > a - 60) & (ev < b_ + 60))
        allb.append(set(u))
    for a in range(len(allb)):
        for b_ in range(a + 1, len(allb)):
            v4 &= not (allb[a] & allb[b_])
    sizes = [len(L16.slots(j)[k]) for j in L16.NETS_DEV + L16.NETS_EVAL for k in ("cat", "STD", "DENSE", "NOISE")]
    P(f"V4 slots eligible, >= 60 s from every catalogued origin, disjoint within and between networks: {v4}; slots per "
      f"use and network min {min(sizes)} max {max(sizes)}; block span (h) " +
      ", ".join(f"{min(L16.slots(j)['cat'] + L16.slots(j)['NOISE']) * 40 / 3600:.2f}-"
                f"{max(L16.slots(j)['cat'] + L16.slots(j)['NOISE']) * 40 / 3600:.2f}" for j in L16.NETS_DEV + L16.NETS_EVAL))
    P("ALL VALIDATION CHECKS PASSED" if (ok_all and v4) else "VALIDATION FAILURE")


if __name__ == "__main__":
    main()
