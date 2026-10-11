"""C8 estimator development for the catalogue amplitude terms (declared: networks j = 60-63, 3 cells, thin, pers
world sigma 0.5; uses catalogue data and the true persistent site factors only, NO detection runs).
    python c8_sitedev.py > out/c8_sitedev.txt"""
import numpy as np
import c8_lib as C
import c6_lib as L
import xc
X = C.X
from task1w_mech import rand_dc


def feats(net, gate=25.0):
    STA, ns = net.STA, len(net.STA)
    rng = np.random.default_rng(5)
    Ms = np.array([[M_[0,0], M_[1,1], M_[2,2], M_[0,1], M_[0,2], M_[1,2]] for M_ in (rand_dc(rng) for _ in range(300))])
    out = []
    for D, src, tau in net.cat:
        MF, sd = xc.matched(D)
        MFn = MF / sd
        E = np.sum(MFn ** 2, axis=1)
        tobs, _tp, en = L.pick_event(E, STA, np.asarray(src, float), float(tau), L.Homog(), center_corr=net.corr)
        k = np.clip(np.rint(tobs / C.DT).astype(int), 0, C.NT - 1)
        v = np.array([[MFn[i, :, k[i, b]] for i in range(ns)] for b in range(2)])
        G = X.radiation_operator(np.asarray(src, float), STA).reshape(2, ns, 3, 6)
        rms = np.sqrt(np.mean(np.einsum("bick,mk->mbic", G, Ms) ** 2, axis=(0, 3)) * 3)
        out.append((G, v, (en >= gate).T, rms))
    return out


def est_geo(F, ns, min_obs=3):
    lr = [[[] for _ in range(ns)] for _ in range(2)]
    for G, v, ok, rms in F:
        if ok.sum() < 2:
            continue
        l = np.log(np.linalg.norm(v, axis=-1) / rms)
        ev = np.median(l[ok])
        for b in range(2):
            for i in range(ns):
                if ok[b, i]:
                    lr[b][i].append(l[b, i] - ev)
    return lr


def fin(lr, ns, min_obs=3):
    u = np.zeros((2, ns)); est = np.zeros((2, ns), bool)
    for b in range(2):
        for i in range(ns):
            if len(lr[b][i]) >= min_obs:
                u[b, i] = np.median(lr[b][i]); est[b, i] = True
    if est.any():
        u[est] -= u[est].mean()
    return u, est


def est_geo_iter(F, ns, rounds=3):
    u = np.zeros((2, ns))
    for _ in range(rounds):
        lr = [[[] for _ in range(ns)] for _ in range(2)]
        for G, v, ok, rms in F:
            if ok.sum() < 2:
                continue
            l = np.log(np.linalg.norm(v, axis=-1) / rms) - u
            ev = np.median(l[ok])
            for b in range(2):
                for i in range(ns):
                    if ok[b, i]:
                        lr[b][i].append(l[b, i] + u[b, i] - ev)
        u, est = fin(lr, ns)
    return u, est


def est_mt(F, ns, nodal=None, rounds=4):
    u = np.zeros((2, ns))
    for _ in range(rounds):
        lr = [[[] for _ in range(ns)] for _ in range(2)]
        for G, v, ok, rms in F:
            if ok.sum() < 2:
                continue
            s = np.exp(u)
            Gw = (G * s[:, :, None, None])[ok].reshape(-1, 6)
            m, *_ = np.linalg.lstsq(Gw, v[ok].reshape(-1), rcond=None)
            g = np.linalg.norm(np.einsum("bick,k->bic", G, m), axis=-1)
            thr = 0 if nodal is None else nodal * np.median(g[ok])
            for b in range(2):
                for i in range(ns):
                    if ok[b, i] and g[b, i] > max(thr, 1e-9):
                        lr[b][i].append(np.log(np.linalg.norm(v[b, i]) / g[b, i]))
        u, est = fin(lr, ns)
    return u, est


res = {}
for j in range(60, 64):
    for cell in xc.T.CELLS:
        net = C.WNet(cell, j, "thin", C.World("pers", 0.5))
        ns = len(net.STA)
        F = feats(net)
        tru = np.log(net.site_net.T)
        cand = {"MT (c8_lib)": (net.u_cat, np.array(net.u_info["n_obs"]) >= 3),
                "MT nodal0.5": est_mt(F, ns, 0.5),
                "GEO": fin(est_geo(F, ns), ns),
                "GEO iter": est_geo_iter(F, ns)}
        for k, (u, est) in cand.items():
            t = tru[est] - tru[est].mean()
            e = u[est]
            res.setdefault(k, []).append((np.corrcoef(t, e)[0, 1], np.sqrt(np.mean((t - e) ** 2)), np.std(t), est.mean()))
print("catalogue amplitude-term estimators, pers world sigma 0.5, networks 60-63 x 3 cells, thin (12 networks)")
for k, v in res.items():
    v = np.array(v)
    print(f"  {k:14s} corr mean {v[:,0].mean():.3f} (min {v[:,0].min():.2f})  rmse {v[:,1].mean():.3f}  (sd of truth {v[:,2].mean():.3f})  fraction estimated {v[:,3].mean():.2f}")
