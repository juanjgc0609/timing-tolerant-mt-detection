"""C8 validation (seed 999; network j = 60 for V1-V5, TAG 3810 validation windows; V6-V7 networks 60-69 catalogues).
    python c8_validate.py > out/c8_validation.txt"""
import time

import numpy as np

import c2_run as C2
import c6_lib as L
import c8_lib as C
import xc

GR = C.GR
TAGV = 3810


def main():
    t00 = time.time()
    print("C8 VALIDATION")
    # V1 generator copy
    w0 = C.World()
    net = C.WNet("ladder_6", 60, "strong", w0)
    b = [999, TAGV, net.ci, 60, net.ri]
    r1 = [np.random.default_rng(b + [s]) for s in (41, 233, 61)]
    r2 = [np.random.default_rng(b + [s]) for s in (41, 233, 61)]
    md = 0.0
    for _ in range(5):
        srcs, taus = C2.draw_dense(r1[0])
        srcs2, taus2 = C2.draw_dense(r2[0])
        D1, p1, _t, _a = GR.generate_real(net.STA, srcs, taus, net.nclut, r1[0], r1[1], r1[2], net.clut_amp,
                                          statics=net.statics)
        D2_, p2, site, _a = C.generate_site(net.STA, srcs2, taus2, net.nclut, r2[0], r2[1], r2[2], net.clut_amp,
                                            statics=net.statics)
        md = max(md, float(np.abs(D1 - D2_).max()))
        assert p1 == p2
    print(f"V1 generate_site == generate_real (5 dense windows, max |dD|) {md:.1e}  {'PASS' if md == 0 else 'FAIL'}")
    # V2 catalogue / corrections
    nc = xc.NetCfg("ladder_6", 60, "strong")
    dmax = max(float(np.abs(a[0] - b_[0]).max()) for a, b_ in zip(nc.cat, net.cat))
    dc = float(np.abs(nc.corr - net.corr).max())
    print(f"V2 WNet(default world) catalogue == xc.NetCfg: max |dD| {dmax:.1e}, max |dcorr| {dc:.1e}  "
          f"{'PASS' if dmax == 0 and dc == 0 else 'FAIL'}")
    # V3 engine == c6_lib.mtfw_dense
    grid = xc.make_grid()
    xops = xc.Ops(net.STA, net.corr, grid)
    dops = C.DOps(xops)
    mops = L.ModelOps(net.STA, grid, L.Homog(), net.corr)
    print(f"   IDX equal to ModelOps: {np.array_equal(dops.IDX, mops.IDX)}; Qst equal: {np.allclose(dops.Qst, mops.Qst)}")
    st = net.stream(TAGV, C2.draw_dense)
    ok3 = ok3b = 0
    wins = [next(st) for _ in range(6)]
    for w in wins:
        a = L.mtfw_dense(mops, w["MFn"], w["truth"], 7)
        b_ = C.mtfw_general(dops, w["MFn"], w["truth"], 7, "ls", 6.0)
        ok3 += int(a[1] == b_[1] and all(h1 == h2 and abs(s1 - s2) <= 1e-9 * abs(s1) for (s1, h1), (s2, h2) in
                                         zip(a[0], b_[0])))
        a = L.mtfw_dense(mops, w["MFn"], w["truth"], 7, cap=False)
        b_ = C.mtfw_general(dops, w["MFn"], w["truth"], 7, "ls", None)
        ok3b += int(a[1] == b_[1] and all(h1 == h2 and abs(s1 - s2) <= 1e-9 * abs(s1) for (s1, h1), (s2, h2) in
                                          zip(a[0], b_[0])))
    print(f"V3 mtfw_general(ls, cap 6) == c6_lib.mtfw_dense: {ok3}/6; uncapped: {ok3b}/6  "
          f"{'PASS' if ok3 == 6 and ok3b == 6 else 'FAIL'}")
    # V4 unit row scaling
    ok4 = 0
    for w in wins[:3]:
        a = C.mtfw_general(dops, w["MFn"], w["truth"], 7, "ls", None)
        b_ = C.mtfw_general(dops, w["MFn"], w["truth"], 7, "ls", None, wsite=np.ones((2, dops.ns)))
        ok4 += int(a[1] == b_[1] and np.allclose([s for s, _ in a[0]], [s for s, _ in b_[0]], rtol=1e-9))
    print(f"V4 unit site scaling == unscaled: {ok4}/3  {'PASS' if ok4 == 3 else 'FAIL'}")
    # V5 limits of the robust statistics
    w = wins[0]
    K, V = C._leaf(w["MFn"], None)
    stt = np.arange(dops.ns)[None, None, None, :]
    tpk = K[stt, dops.IDX]
    vf = V[stt, tpk].reshape(dops.IDX.shape[0], dops.IDX.shape[1], -1)
    y = np.einsum("gtr,grk->gtk", vf, dops.Qst)
    ls = (y * y).sum(-1) / 2
    old = C.C_HUB
    C.C_HUB = 1e9
    hb, _a = C.stat_huber(vf, dops.Qst, False)
    C.C_HUB = old
    gn, _m, _s = C.stat_gain(vf, dops.Gs, dops.Gp, 1e-4)
    gf, _m, _s = C.stat_gain(vf, dops.Gs, dops.Gp, np.inf)
    hn, _a = C.stat_huber(vf, dops.Qst, False)
    gp, _m, _s = C.stat_gain(vf, dops.Gs, dops.Gp, 0.5)
    print(f"V5 Huber(c -> inf) == LS: max rel {np.max(np.abs(hb - ls) / ls):.1e}; GAIN(sp -> 0) == LS: max rel "
          f"{np.max(np.abs(gn - ls) / ls):.1e}")
    print(f"   ordering (should hold per candidate up to IRLS/grid tolerance): GAIN-FREE >= GAIN(0.5): "
          f"{np.mean(gf >= gp - 1e-6):.4f}; GAIN(0.5) >= LS: {np.mean(gp >= ls - 1e-6):.4f}; HUBER-N <= LS (not implied): "
          f"{np.mean(hn <= ls + 1e-6):.4f}")
    # V6 oracle captured fraction, noise-free single events, pers world
    wp = C.World("pers", 0.5)
    fr = {"LS": [], "ORACLE": [], "SITE-CAT": [], "SITE-CAT-ST": []}
    cors = []
    for j in range(60, 64):
        for cell in ("ladder_4", "unilateral_c40"):
            np_ = C.WNet(cell, j, "thin", wp)
            xo = xc.Ops(np_.STA, np_.corr, grid)
            do = C.DOps(xo)
            ok = np.array(np_.u_info["n_obs"]) >= 3
            cors.append(np.corrcoef(np_.u_cat[ok], np.log(np_.site_net.T)[ok])[0, 1])
            Qs = {"LS": do.Qst, "ORACLE": do.scaled(np_.site_net.T)[1], "SITE-CAT": do.scaled(np.exp(np_.u_cat))[1],
                  "SITE-CAT-ST": do.scaled(np.exp(np_.u_cat_st))[1]}
            rng = np.random.default_rng([999, TAGV, 7, j, np_.ci])
            for _ in range(4):
                gi = int(rng.integers(len(do.garr)))
                src, tau = do.garr[gi], 6.0
                D, _p, site, _a = np_.gen([src], [tau], 0, rng, np.random.default_rng(int(rng.integers(1 << 30))),
                                          np.random.default_rng(int(rng.integers(1 << 30))))
                MF, sd = xc.matched(D)
                MFn = MF / sd
                K, V = C._leaf(MFn, None)
                ti = int(np.argmin(np.abs(do.taus - tau)))
                tp = K[np.arange(do.ns)[None, :], do.IDX[gi, ti]]
                v = V[np.arange(do.ns)[None, :], tp].reshape(-1)
                for k_, Q in Qs.items():
                    yy = Q[gi].T @ v
                    fr[k_].append(float(yy @ yy / (v @ v)))
    print("V6 pers world (sigma 0.5), 32 single events with noise (no clutter), captured fraction at the true node "
          "(mean): " + ", ".join(f"{k} {np.mean(v):.3f}" for k, v in fr.items()))
    print(f"V7 catalogue amplitude terms vs true log site (pers, 8 networks): corr per network "
          f"{np.round(cors, 2).tolist()} mean {np.mean(cors):.3f}")
    # V8 bound-ordered exact argmax (GAIN, GAIN-FREE) == full-lattice argmax; one engine run per robust stat
    okc = tot = 0
    for w in wins[:4]:
        K, V = C._leaf(w["MFn"], None)
        tpk = K[stt, dops.IDX]
        vf = V[stt, tpk].reshape(dops.IDX.shape[0], dops.IDX.shape[1], -1)
        exc = np.zeros(len(dops.garr), bool)
        exc[:60] = True
        for sp in (0.5, np.inf):
            gi, tj, s_, _f = C.robust_argmax(vf, exc, "gain", dops.Qst, dops.Gs, dops.Gp, sp)
            full = C.stat_gain(vf, dops.Gs, dops.Gp, sp)[0]
            full[exc] = -np.inf
            tot += 1
            okc += int(gi * vf.shape[1] + tj == int(np.argmax(full)) and abs(s_ - full.max()) <= 1e-8 * full.max())
    print(f"V8 bound-ordered argmax == full-lattice argmax (GAIN, GAIN-FREE, 4 windows): {okc}/{tot}  "
          f"{'PASS' if okc == tot else 'FAIL'}")
    for stat, kw in (("huber_n", {}), ("huber_s", {}), ("gain", {"sp": 0.5}), ("gain", {"sp": np.inf})):
        t = time.time()
        r, _l = C.mtfw_general(dops, wins[0]["MFn"], wins[0]["truth"], 7, stat, **kw)
        print(f"   engine {stat} {kw}: {len(r)} atoms, hits {[h for _s, h in r]}, {time.time() - t:.2f} s")
    print(f"total {time.time() - t00:.0f} s")


if __name__ == "__main__":
    main()
