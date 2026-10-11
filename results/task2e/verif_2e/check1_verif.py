"""Check 1: independent review checks of the continuous atom model of task2e_methods against the grid model (BOMP-PK).
Seed 999 calibration networks only.  OPENBLAS_CORETYPE=Prescott python3 check1_verif.py"""
import time

import numpy as np

import common_verif as V
from task1x_mtglrt import radiation_operator
from task0_common import DT, NT, VP, VS, FC

E, D2, M, T = V.E, V.D2, V.M, V.T
L = []


def P_(s):
    print(s, flush=True)
    L.append(s)


def ricker_own(t):
    x = (np.pi * FC * t) ** 2
    return (1 - 2 * x) * np.exp(-x)


def explicit_A(STA, corr, P, ns):
    """own data-domain dictionary: columns G_k(x) (x) ricker(t - T_kbi) / sqrt(E_R), T = d/V + t0 + corr + delta."""
    tt = np.arange(NT) * DT
    na = len(P)
    A = np.zeros((ns, 3, NT, na, 6))
    for k in range(na):
        x = P[k, :3]
        G = radiation_operator(x, STA).reshape(2, ns, 3, 6)
        d = np.linalg.norm(STA - x, axis=1)
        dl = P[k, 4:].reshape(2, ns)
        for b, v in enumerate((VP, VS)):
            Tb = d / v + P[k, 3] + corr[:, b] + dl[b]
            for i in range(ns):
                A[i, :, :, k, :] += np.einsum("cx,n->cnx", G[b, i], ricker_own(tt - Tb[i])) / np.sqrt(M.E_R)
    return A.reshape(-1, 6 * na)


def J_own(STA, corr, ns, y, P, m, lam):
    A = explicit_A(STA, corr, P, ns)
    r = y - A @ m.ravel()
    Om = sum(np.linalg.norm(A[:, 6 * k:6 * k + 6] @ m[k]) for k in range(len(P)))
    return 0.5 * r @ r + lam * Om


def main():
    t0 = time.time()
    rng = np.random.default_rng(99901)
    P_("# Check 1 (verif_2e): continuous atom model vs grid model, independent gradients, sign conventions")
    # ---------------------------------------------------------------- A/B: grid equivalence on real windows
    P_("")
    P_("## 1A. theta on the grid node, shifts = BOMP-PK peak offsets: TS-SFW atom vs BOMP-PK block (5-atom BOMP-PK "
       "solutions on real seed-999 windows)")
    worst = {"pred": 0.0, "K": 0.0, "b": 0.0, "Jfit": 0.0}
    n_off = n_clip = 0
    n_atoms = 0
    e_rel = []
    s_rel = []
    for cell in ("ladder_4", "ladder_6", "unilateral_c40"):
        for regime in ("thin", "strong"):
            S = V.setup(cell, regime, 1)
            net, PD = S["net"], S["PD"]
            for w in V.windows(S, cell, regime, 1, 4):
                win = E.Win(w["D"], w["MFn"], w["sd"])
                dg = {}
                D2.bomp_pk_record(PD, w["MFn"], w["truth"], diag=dg)
                atoms, U = dg["atoms"], dg["U"]
                Pg = np.array([E.init_atom(net, n, t, I, clip=False) for n, t, I in atoms])
                # sampled prediction: MFn - residual of BOMP-PK vs predict_sampled of TS-SFW
                pred_b = w["MFn"] - PD.residual(w["MFn"], atoms, U)
                pred_t = E.predict_sampled(net, win, Pg, U)
                worst["pred"] = max(worst["pred"], float(np.abs(pred_b - pred_t).max() / np.abs(pred_b).max()))
                Kc, bc = E.gram(net, win, Pg)
                Kd, bd = PD.gram(atoms), PD.corr_atoms(w["MFn"], atoms).ravel()
                worst["K"] = max(worst["K"], float(np.abs(Kc - Kd).max() / np.abs(Kd).max()))
                worst["b"] = max(worst["b"], float(np.abs(bc - bd).max() / np.abs(bd).max()))
                # objective of the joint fit (lambda = 0) at the BOMP-PK solution: -u^T b + 1/2 u^T K u
                Jt = E.joint_obj(net, win, Pg, U, 0.0, grad=False)
                Jb = float(-U.ravel() @ bd + 0.5 * U.ravel() @ Kd @ U.ravel())
                worst["Jfit"] = max(worst["Jfit"], abs(Jt - Jb) / abs(Jb))
                # shifts outside [-0.25, 0.25] (clipped by init_atom(clip=True))
                dl = Pg[:, 4:]
                n_off += dl.size
                n_clip += int(np.sum(np.abs(dl) > E.SHIFT_MAX + 1e-12))
                n_atoms += len(atoms)
                # first atom: refinement energy at the init = oracle value (MT-FW-NOCAP) ; BOMP-PK e/2
                s, n, t, I = E.oracle(net, w["MFn"], np.zeros(net.N, bool))
                p0 = E.init_atom(net, n, t, I, clip=False)
                Eh, _mu = E.refine_obj(net, win, p0, np.zeros((0, 2, net.ns, 3)), np.zeros((0, 2, net.ns)),
                                       grad=False)
                s_rel.append(abs(-Eh - s) / s)
                eb, _Ib, _h = PD.energy(w["MFn"], True)
                e_rel.append(abs(eb[n, t] / 2 - s) / s)
    ok = max(worst.values()) < 1e-10 and max(s_rel) < 1e-10
    P_(f"- {n_atoms} BOMP-PK atoms in 24 windows (network 1, 4 per cell x regime): max rel. difference sampled "
       f"prediction {worst['pred']:.1e}, Gram {worst['K']:.1e}, correlation {worst['b']:.1e}, joint LS objective "
       f"{worst['Jfit']:.1e} -> {'IDENTICAL' if ok else 'DIFFERENT'}")
    P_(f"- refinement energy E/2 at the grid init = oracle value s (MT-FW-NOCAP statistic): max rel. diff "
       f"{max(s_rel):.1e}; BOMP-PK's e/2 = h^T K(I)^+ h / 2 differs from s by the P/S cross term only: median rel. "
       f"{np.median(e_rel):.1e}, max {max(e_rel):.1e}")
    P_(f"- peak offsets outside +-0.25 s in continuous time (BOMP-PK reads idx + [-4, 4] samples with idx = "
       f"floor(T0/DT), i.e. [-0.25 - frac, 0.25 - frac] s; TS-SFW clips to [-0.25, 0.25]): {n_clip}/{n_off} "
       f"({n_clip / n_off:.3f}); a clipped shift moves the read time by < 1 sample (0.0625 s)")
    # ---------------------------------------------------------------- C: independent objective and gradient
    P_("")
    P_("## 1C. objective and gradient vs an independent data-domain implementation (own ricker, own travel times, "
       "task1x_mtglrt.radiation_operator; central differences on the independent J)")
    errs, jerr = [], []
    for trial in range(6):
        cell = ("ladder_4", "ladder_6", "unilateral_c40")[trial % 3]
        regime = ("thin", "strong")[trial % 2]
        S = V.setup(cell, regime, 2)
        net = S["net"]
        w = V.windows(S, cell, regime, 2, 1 + trial % 2)[-1]
        win = E.Win(w["D"], w["MFn"], w["sd"])
        y = win.Dn.ravel() * np.sqrt(M.E_R)
        na = 1 + trial % 3
        P = np.column_stack([rng.uniform(-18, 18, na), rng.uniform(-18, 18, na), rng.uniform(2, 18, na),
                             rng.uniform(5.0, 8.0, na), rng.uniform(-0.25, 0.25, (na, 2 * net.ns))])
        m = rng.normal(size=(na, 6)) * 2.0
        for lam in (0.0, 1.7):
            J, g = E.joint_obj(net, win, P, m, lam)
            Jo = J_own(S["STA"], S["corr"], net.ns, y, P, m, lam)
            jerr.append(abs((J + 0.5 * win.yy) - Jo) / abs(Jo))
            x = P.ravel().copy()
            gf = np.zeros_like(x)
            steps = np.tile(np.array([2e-4] * 3 + [2e-5] * (1 + 2 * net.ns)), na)
            for i in range(len(x)):
                e = np.zeros_like(x)
                e[i] = steps[i]
                gf[i] = (J_own(S["STA"], S["corr"], net.ns, y, (x + e).reshape(na, -1), m, lam) -
                         J_own(S["STA"], S["corr"], net.ns, y, (x - e).reshape(na, -1), m, lam)) / (2 * steps[i])
            rel = float(np.linalg.norm(gf - g.ravel()) / np.linalg.norm(gf))
            errs.append(rel)
            P_(f"- {cell} {regime}, {na} atom(s), lambda {lam}: J(method) + ||y||^2/2 = {J + 0.5 * win.yy:.6f}, "
               f"J(independent) = {Jo:.6f} (rel {jerr[-1]:.1e}); ||g_fd(indep) - g(method)|| / ||g_fd|| = {rel:.1e}")
        # refinement energy: independent h via the explicit dictionary on the residual y - A_old m_old
        if na >= 2:
            A_old = explicit_A(S["STA"], S["corr"], P[:-1], net.ns)
            r = y - A_old @ m[:-1].ravel()
            A_new = explicit_A(S["STA"], S["corr"], P[-1:], net.ns)
            h = A_new.T @ r
            G = radiation_operator(P[-1, :3], S["STA"])
            Eo = float(h @ np.linalg.pinv(G.T @ G) @ h)
            a_old, T_old = E.pred_old(net, win, P[:-1], m[:-1])
            v, _mu = E.refine_obj(net, win, P[-1], a_old, T_old, grad=False)
            jerr.append(abs(-2 * v - Eo) / Eo)
            P_(f"  refinement energy with {na - 1} active atom(s): method {-2 * v:.6f}, independent {Eo:.6f} "
               f"(rel {jerr[-1]:.1e})")
    ok_c = max(errs) < 1e-5 and max(jerr) < 1e-9
    P_(f"- max objective discrepancy {max(jerr):.1e}; max gradient discrepancy {max(errs):.1e} (FD truncation level "
       f"for these steps) -> {'CORRECT' if ok_c else 'CHECK'}")
    # ---------------------------------------------------------------- D: sign conventions (noise-free, own data)
    P_("")
    P_("## 1D. sign conventions: noise-free source ON a grid node, origin on the tau grid, real network corrections, "
       "data built with own code at T = d/V + t0 + corr + delta_true (delta_true = +0.12 s on P of station 0, "
       "-0.08 s on S of station 1, 0 elsewhere)")
    for cell in ("ladder_6", "unilateral_c40"):
        S = V.setup(cell, "thin", 0)
        net, PD, STA, corr = S["net"], S["PD"], S["STA"], S["corr"]
        n_node = int(np.argmin(np.linalg.norm(net.garr - np.array([2.0, -3.0, 9.0]), axis=1)))
        it = 12
        x = net.garr[n_node]
        t0s = float(net.taus[it])
        dtrue = np.zeros((2, net.ns))
        dtrue[0, 0], dtrue[1, 1] = 0.12, -0.08
        m = rng.normal(size=6)
        G = radiation_operator(x, STA).reshape(2, net.ns, 3, 6)
        a = np.einsum("bicx,x->bic", G, m)
        m = m * 12.0 / np.abs(a).max()
        a = a * 12.0 / np.abs(a).max()
        d = np.linalg.norm(STA - x, axis=1)
        tt = np.arange(NT) * DT
        Dn = np.zeros((net.ns, 3, NT))
        Ttrue = np.zeros((2, net.ns))
        for b, v in enumerate((VP, VS)):
            Ttrue[b] = d / v + t0s + corr[:, b] + dtrue[b]
            for i in range(net.ns):
                Dn[i] += np.outer(a[b, i], ricker_own(tt - Ttrue[b, i])) / M.E_R
        MFn, _sd = T.matched(Dn)
        win = E.Win(Dn, MFn, 1.0)
        dg = {}
        D2.bomp_pk_record(PD, MFn, [x], diag=dg, kmax=1)
        n1, t1, I1 = dg["atoms"][0]
        idx = net.BD.idx[n1, t1]
        exp_off = np.round(Ttrue / DT).astype(int) - idx
        P_(f"- {cell}: BOMP-PK first atom node {'= true node' if n1 == n_node else '!= true node'}, origin "
           f"{'= true' if t1 == it else '!= true'}; peak offsets I - idx on (P, st0) {int(I1[0, 0] - idx[0, 0]):+d}, "
           f"(S, st1) {int(I1[1, 1] - idx[1, 1]):+d} samples; expected round(T_true/DT) - idx: {int(exp_off[0, 0]):+d}, "
           f"{int(exp_off[1, 1]):+d}; all offsets equal expected: {bool(np.array_equal(I1 - idx, exp_off))}")
        # TS-SFW refinement from the grid init
        s, n, t, I = E.oracle(net, MFn, np.zeros(net.N, bool))
        p0 = E.init_atom(net, n, t, I)
        p1, inf = E.refine_atom(net, win, p0, np.zeros((0, 2, net.ns, 3)), np.zeros((0, 2, net.ns)), maxiter=500)
        _G, _dG, Tp, _d = net.geom(p1[None])
        dl = p1[4:].reshape(2, net.ns)
        P_(f"  TS-SFW refinement (500 it.): position error {np.linalg.norm(p1[:3] - x):.4f} km, t0 error "
           f"{p1[3] - t0s:+.4f} s, max |T_est - T_true| {np.abs(Tp[0] - Ttrue).max():.4f} s; estimated "
           f"delta(P, st0) - common {dl[0, 0] - np.median(dl):+.4f} s (true +0.12), delta(S, st1) - common "
           f"{dl[1, 1] - np.median(dl):+.4f} s (true -0.08)")
        # gradient sign at delta = 0 for the true atom (lambda = 0, m = true)
        ptrue = np.concatenate([x, [t0s], dtrue.ravel()])
        pz = np.concatenate([x, [t0s], np.zeros(2 * net.ns)])
        _J, gt = E.joint_obj(net, win, ptrue[None], m[None] * np.sqrt(M.E_R) ** 0, 0.0)
        mm, _c = E.convex_step(net, win, ptrue[None], np.zeros((1, 6)), 0.0)
        _J, gt = E.joint_obj(net, win, ptrue[None], mm, 0.0)
        _J, gz = E.joint_obj(net, win, pz[None], mm, 0.0)
        gdz = gz[0, 4:].reshape(2, net.ns)
        P_(f"  joint objective at the truth (m = LS fit): ||grad|| {np.linalg.norm(gt):.1e}; at delta = 0: "
           f"dJ/d delta(P, st0) = {gdz[0, 0]:+.3f} (negative => increasing delta lowers J, toward +0.12), "
           f"dJ/d delta(S, st1) = {gdz[1, 1]:+.3f} (positive => toward -0.08)")
    P_("")
    P_(f"(check 1 runtime {time.time() - t0:.0f} s)")
    open(V.HERE + "/verif_2e/check1_verif.txt", "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
