"""C5: location accuracy and uncertainty of hits (single-event windows, TAG 3006; thin; 3 km grid).
  NODE: BOMP-PK+SC hit node.  TS: TS-SFW+SC hit (continuous, free shifts).  CR: continuous refinement of the BOMP-PK
  hit with shifts fixed at 0 (x, y, z, t0 only), maximising the projection energy E of 2E (L-BFGS-B, <= 200 it).
  U1 (AA.6 curvature): Cov = inverse Hessian of -E/2 at CR (unit noise in the matched-filter domain).
  U3 (travel-time residuals): peak times (+-0.5 s, energy >= 25) vs CR prediction; Cov = s^2 (J^T J)^-1."""
import json, sys, time
import numpy as np
import xc

E = xc.E
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")


def draw_one(rng):
    s = np.array([rng.uniform(-15, 15), rng.uniform(-15, 15), rng.uniform(5, 13)])
    return [s], [float(rng.uniform(5.5, 7.5))]


def cr_refine(net, win, x0):
    lo, hi = net.lo.copy(), net.hi.copy()
    lo[4:] = 0.0
    hi[4:] = 0.0
    e0 = np.zeros((0, 2, net.ns, 3)), np.zeros((0, 2, net.ns))

    def fun(p):
        v, g, _mu = E.refine_obj(net, win, p, *e0)
        return v, g
    p, r = E._lbfgs(fun, x0, lo, hi, net.scale, 200)
    return p, fun, r


def hessian4(fun, p, h=(0.05, 0.05, 0.05, 0.005)):
    H = np.zeros((4, 4))
    for a in range(4):
        dp = np.zeros_like(p)
        dp[a] = h[a]
        H[a] = (fun(p + dp)[1][:4] - fun(p - dp)[1][:4]) / (2 * h[a])
    return 0.5 * (H + H.T)


def tt_cov(net, MFn, p, gate=25.0, search=0.5):
    G, _dG, Tm, dTdX = net.geom(p[None])
    Tm, dTdX = Tm[0], dTdX[0]
    Et = np.sum(MFn ** 2, axis=1)
    NT = MFn.shape[-1]
    W = int(round(search / xc.DT))
    rows, res = [], []
    for b in range(2):
        for i in range(net.ns):
            c = int(round(Tm[b, i] / xc.DT))
            idx = np.arange(max(1, c - W), min(NT - 2, c + W) + 1)
            k = int(idx[np.argmax(Et[i, idx])])
            if Et[i, k] < gate:
                continue
            e0, e1, e2 = Et[i, k - 1], Et[i, k], Et[i, k + 1]
            den = e0 - 2 * e1 + e2
            dk = 0.5 * (e0 - e2) / den if den < 0 else 0.0
            res.append((k + float(np.clip(dk, -0.5, 0.5))) * xc.DT - Tm[b, i])
            rows.append(np.r_[dTdX[b, i], 1.0])
    n = len(res)
    if n < 6:
        return None, n, None
    J = np.array(rows)
    r = np.array(res)
    beta = np.linalg.lstsq(J, r, rcond=None)[0]
    s2 = float(np.sum((r - J @ beta) ** 2) / (n - 4))
    try:
        C = s2 * np.linalg.inv(J.T @ J)
    except np.linalg.LinAlgError:
        return None, n, s2
    return C, n, s2


def run(regime, nets, nwin, out):
    res = []
    for cell in CELLS:
        for j in nets:
            nc = xc.NetCfg(cell, j, regime)
            ops = xc.Ops(nc.STA, nc.corr, xc.make_grid())
            st = xc.window_stream(nc, 3006, draw=draw_one)
            for wi in range(nwin):
                w = next(st)
                tr = w["truth"]
                dg = {}
                rb = xc.D2.bomp_pk_record(ops.PD, w["MFn"], tr, diag=dg)
                row = {"cell": cell, "net": j, "w": wi, "truth": tr[0].tolist(), "tau": w["taus"][0]}
                hit = [(a, s) for (s, h), a in zip(rb, dg["atoms"]) if h == 0]
                if not hit:
                    row["hit"] = False
                    res.append(row)
                    continue
                a, s = hit[0]
                row["hit"] = True
                row["score"] = float(s)
                row["rank"] = int([h for _s, h in rb].index(0))
                node = ops.BD.garr[a[0]]
                row["node"] = node.tolist() + [float(ops.BD.taus[a[1]])]
                rt, lt, _ = xc.tssfw(ops, w, tr, f=xc.F_FROZEN)
                th = [l for (s_, h), l in zip(rt, lt) if h == 0]
                row["ts"] = th[0] if th else None
                win = E.Win(w["D"], w["MFn"], w["sd"])
                x0 = np.r_[node, ops.BD.taus[a[1]], np.zeros(2 * ops.net.ns)]
                t0 = time.time()
                p, fun, r = cr_refine(ops.net, win, x0)
                row["cr"] = p[:4].tolist()
                row["cr_E"] = -2 * float(r.fun)
                row["cr_nit"] = int(r.nit)
                H = hessian4(fun, p)
                row["H"] = H.tolist()
                row["cr_sec"] = time.time() - t0
                C, n, s2 = tt_cov(ops.net, w["MFn"], p)
                row["C_tt"] = None if C is None else C.tolist()
                row["n_arr"], row["s2_tt"] = n, s2
                res.append(row)
            print(cell, j, "done", time.strftime("%H:%M:%S"), flush=True)
            with open(out, "w") as fh:
                json.dump(res, fh)


if __name__ == "__main__":
    run(sys.argv[1], [int(x) for x in sys.argv[2].split(",")], int(sys.argv[3]), sys.argv[4])
