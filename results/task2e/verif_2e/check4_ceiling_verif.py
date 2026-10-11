"""Check 4: headroom / ceiling on the seed-999 calibration windows (2A networks 0-4, 20 windows each, 3 cells x 2
regimes = 600 windows).  Per window it stores, for oracle ('genie') detectors:
  true candidates (one per true source, hit = its index):
    g_true   max of BOMP-PK's block statistic e/2 (peak reading, +-0.25 s) over all nodes within 5 km epicentral of the
             source and all origins, on the full data (no deflation needed: sources are > 12 km apart);
    c_true   continuous single-atom energy E/2 (TS-SFW refinement, 500 it.) started at the best near-truth block and at
             the exact true (x, y, z, t0) with peak-initialised shifts; max over the starts whose final epicentre stays
             within 5 km of the source (else the start value);
  noise candidates (false):
    g_noise  greedy (5 km exclusion, up to 5) maxima of e/2 on the clutter+noise-only data (D minus the regenerated
             event signal, same sd), over nodes > 5 km from every true source;
    c_noise  the continuous refinement of each of those blocks on the clutter+noise-only data;
    g_resid  greedy maxima of e/2 (> 5 km from every source) on the full data after an ORACLE deflation (joint LS fit
             of the truth-start refined atoms): the false alarms left by event-model mismatch after perfect detection.
    OPENBLAS_CORETYPE=Prescott python3 check4_ceiling_verif.py"""
import json
import os
import sys
import time

import numpy as np

import common_verif as V

E, D2, T = V.E, V.D2, V.T
OUT = os.path.join(V.HERE, "verif_2e", "check4_parts_verif")
TOL = 5.0
EMPTY_A = None


def near_nodes(garr, src):
    return np.where(np.linalg.norm(garr[:, :2] - np.asarray(src)[None, :2], axis=1) <= TOL)[0]


def cont_refine(net, win, R, p_start_xyz, t0, maxiter=500):
    """single-atom refinement (no active atoms) from (x, t0) with shifts from the energy peaks of R."""
    Kt = D2.peak_table(R)
    x = np.clip(np.asarray(p_start_xyz, float), net.lo[:3], net.hi[:3])
    t0 = float(np.clip(t0, net.lo[3], net.hi[3]))
    d = np.linalg.norm(net.STA - x, axis=1)
    T0 = d[None, :] * E.VINV[:, None] + t0 + net.corrT
    idx = np.clip((T0 / E.DT).astype(int), 0, E.NT - 1)
    I = Kt[np.arange(net.ns)[None, :], idx]
    dl = np.clip(I * E.DT - T0, -E.SHIFT_MAX, E.SHIFT_MAX)
    p0 = np.concatenate([x, [t0], dl.ravel()])
    z0 = (np.zeros((0, 2, net.ns, 3)), np.zeros((0, 2, net.ns)))
    E0 = -2 * E.refine_obj(net, win, np.clip(p0, net.lo, net.hi), *z0, grad=False)[0]
    p1, inf = E.refine_atom(net, win, p0, *z0, maxiter=maxiter)
    return p0, E0, p1, inf["E1"]


def process(cell, regime, j, nw=20):
    S = V.setup(cell, regime, j)
    net, PD = S["net"], S["PD"]
    garr = net.garr
    stored = json.load(open(os.path.join(V.HERE, "task2e_parts", f"cal_seed999_{regime}_{cell}_net{j}.json")))
    W = []
    for w_i, w in enumerate(V.windows(S, cell, regime, j, nw, with_noise_only=True)):
        st = stored["windows"][w_i]
        assert [list(map(float, s)) for s in w["srcs"]] == st["truth"]
        assert w["regen_equal"]
        srcs = [np.asarray(s, float) for s in w["srcs"]]
        e_full = PD.energy(w["MFn"], True)[0] / 2.0
        e_noise = PD.energy(w["MFn_noise"], True)[0] / 2.0
        win = E.Win(w["D"], w["MFn"], w["sd"])
        win_n = E.Win(w["D_noise"], w["MFn_noise"], w["sd"])
        g_true, c_true, c_dbg = [], [], []
        P_or = []
        for s_i, src in enumerate(srcs):
            nn = near_nodes(garr, src)
            sub = e_full[nn]
            a = int(np.argmax(sub))
            n_best, t_best = nn[a // net.T], a % net.T
            g_true.append(float(sub.max()))
            best = -np.inf
            dbg = []
            for (xs, t0s) in ((garr[n_best], float(net.taus[t_best])), (src, float(w["taus"][s_i]))):
                p0, E0, p1, E1 = cont_refine(net, win, w["MFn"], xs, t0s)
                ok = np.linalg.norm(p1[:2] - src[:2]) <= TOL
                val = E1 if ok else (E0 if np.linalg.norm(p0[:2] - src[:2]) <= TOL else -np.inf)
                best = max(best, val / 2.0)
                dbg.append([E0 / 2, E1 / 2, bool(ok), float(np.linalg.norm(p1[:3] - src))])
            P_or.append(p1)                      # the truth-start refined atom (last start), used for deflation
            c_true.append(float(best))
            c_dbg.append(dbg)
        # noise candidates: greedy on the clutter+noise-only statistic, > 5 km from every true source
        emax = e_noise.max(axis=1)
        allowed = np.ones(net.N, bool)
        for src in srcs:
            allowed[near_nodes(garr, src)] = False
        g_noise, c_noise, nodes = [], [], []
        for _ in range(5):
            cand = np.where(allowed, emax, -np.inf)
            n = int(np.argmax(cand))
            if not np.isfinite(cand[n]):
                break
            t = int(np.argmax(e_noise[n]))
            g_noise.append(float(cand[n]))
            nodes.append(n)
            p0, E0, p1, E1 = cont_refine(net, win_n, w["MFn_noise"], garr[n], float(net.taus[t]))
            c_noise.append(float(E1 / 2.0))
            allowed[np.linalg.norm(garr[:, :2] - garr[n][None, :2], axis=1) <= TOL] = False
        # oracle deflation: joint LS fit of the truth-start atoms, then greedy on the residual statistic
        R_or = w["MFn"]
        if srcs:
            Po = np.array(P_or)
            mo, _c = E.convex_step(net, win, Po, np.zeros((len(Po), 6)), 0.0)
            R_or = w["MFn"] - E.predict_sampled(net, win, Po, mo)
        e_res = PD.energy(R_or, True)[0] / 2.0
        emr = e_res.max(axis=1)
        allowed = np.ones(net.N, bool)
        for src in srcs:
            allowed[near_nodes(garr, src)] = False
        g_resid = []
        for _ in range(5):
            cand = np.where(allowed, emr, -np.inf)
            n = int(np.argmax(cand))
            if not np.isfinite(cand[n]):
                break
            g_resid.append(float(cand[n]))
            allowed[np.linalg.norm(garr[:, :2] - garr[n][None, :2], axis=1) <= TOL] = False
        W.append({"k": len(srcs), "g_resid": g_resid, "g_true": g_true, "c_true": c_true, "c_dbg": c_dbg, "g_noise": g_noise,
                  "c_noise": c_noise, "noise_nodes": nodes,
                  "rec": {m: st["records"][m] for m in ("TS-SFW+SC f=0.01", "BOMP-PK+SC", "BOMP-PK-SLIDE+SC",
                                                          "MT-FW+SC")}})
    return W


def main():
    os.makedirs(OUT, exist_ok=True)
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    for cell in T.CELLS:
        for regime in T.REGIMES:
            for j in range(T.N_NETS):
                path = os.path.join(OUT, f"c_{regime}_{cell}_net{j}.json")
                if os.path.exists(path):
                    continue
                t0 = time.time()
                W = process(cell, regime, j, nw)
                with open(path + ".tmp", "w") as fh:
                    json.dump({"cell": cell, "regime": regime, "net": j, "windows": W}, fh)
                os.replace(path + ".tmp", path)
                print(f"{regime} {cell} net{j}: {len(W)} windows in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
