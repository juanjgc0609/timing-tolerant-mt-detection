"""
task2e_validate.py - validaciones de task2e_methods.py en la semilla 999 (redes de calibracion de 2A, nucleo Prescott),
antes de la calibracion y de cualquier corrida de evaluacion.  Salida: task2e_validation.txt.

V1 identidades: AC continuo en las muestras = task1x_mtglrt.AC; rad_op = radiation_operator; MFc en tiempos de muestra =
   MFn; con los tiempos en muestras enteras, Gram y correlacion continuas = las de PKDict (BOMP-PK).
V2 gradientes analiticos de J (lambda = 0 y lambda > 0) y de la energia del refinamiento contra diferencias finitas
   centrales, en ventanas reales con 1-3 atomos en posiciones y desplazamientos aleatorios.
V3 modelo continuo contra el diccionario explicito en el dominio de los datos (columnas G_k (x) ricker(t - T) / sqrt(E_R),
   y = Dn sqrt(E_R)): lambda = 0 del paso convexo = lstsq explicito; lambda > 0: el objetivo y el KKT del paso convexo
   contra FISTA independiente sobre el diccionario explicito; brecha y KKT al final del paso convexo.
V4 fuente unica fuera de la rejilla sin ruido, con desplazamientos inyectados de 0.1-0.2 s (signo aleatorio) por estacion y
   fase: TS-SFW (lambda = f s_1 para cada f de la calibracion, y lambda = 0) recupera la posicion (< 0.5 km) y los
   desplazamientos (< 0.02 s, modulo el desplazamiento comun, que es indistinguible del tiempo de origen); BOMP-PK no.
V5 sin refinamiento ni deslizamiento y lambda = 0, el primer atomo es el de BOMP-PK (nodo, origen, muestras de pico,
   amplitudes) cuando el oraculo elige el mismo bloque; con el estadistico de BOMP-PK en el oraculo y sin recorte
   de los desplazamientos iniciales, el registro completo es el de BOMP-PK; el primer valor del oraculo es el de
   MT-FW-NOCAP+SC (task2d_methods.fw_record_v).

    OPENBLAS_CORETYPE=Prescott python task2e_validate.py
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import task2d_run as R2  # noqa: E402
import task2e_methods as E  # noqa: E402
from task1x_mtglrt import radiation_operator  # noqa: E402
from task1w_mech import rand_dc  # noqa: E402

T = R2.T
B = R2.B
D2 = R2.D2
M = R2.M
FGRID = (0.01, 0.03, 0.1, 0.3)
OUT = os.path.join(HERE, "task2e_validation.txt")
N_RAND = 30
V4_SUMMARY = {}


def setup(cell, regime="thin", j=0):
    ci, nst, nclut, STA, statics, cat, drawn, corr = B.setup_net(T.CAL_SEED, regime, cell, j)
    grid, Gs, Qs, BD, ST, PD = R2.ops(STA, corr)
    net = E.Net(STA, corr, grid, T.DEFAULT_TAUS, Gs, Qs, BD)
    return net, PD, STA, statics, nclut, corr, grid, Gs, Qs


def windows(STA, statics, nclut, cell, regime="thin", j=0, n=4):
    rngs = T.window_rngs(T.CAL_SEED, cell, j, regime)
    out = []
    for _ in range(n):
        srcs, taus, D, _p, tr, _a = T.generate_window(STA, statics, nclut, T.CLUT_AMP[regime], rngs)
        MF, sd = T.matched(D)
        out.append((srcs, taus, D, MF / sd, sd, tr))
    return out


def fd_grad(fun, x, steps):
    g = np.zeros_like(x)
    for i in range(len(x)):
        e = np.zeros_like(x)
        e[i] = steps[i]
        g[i] = (fun(x + e) - fun(x - e)) / (2 * steps[i])
    return g


def explicit_dict(net, P, sqer):
    """diccionario explicito en el dominio de los datos: (ns*3*NT, 6 na)."""
    G, _dG, Tt, _d = net.geom(P)
    na = len(P)
    tt = np.arange(E.NT) * E.DT
    A = np.zeros((net.ns, 3, E.NT, na, 6))
    for k in range(na):
        for p in range(2):
            for i in range(net.ns):
                r, _rd = E.ricker_d(tt - Tt[k, p, i])
                A[i, :, :, k, :] += np.einsum("cx,n->cnx", G[k, p, i], r) / sqer
    return A.reshape(-1, 6 * na)


def fista_group(A, y, lam, na, iters=20000):
    """lasso de grupos con ||A_k m_k|| por FISTA en coordenadas blanqueadas (independiente del codigo del metodo)."""
    Ws, Tis = [], []
    for k in range(na):
        Ak = A[:, 6 * k:6 * k + 6]
        U, s, Vt = np.linalg.svd(Ak, full_matrices=False)
        Ws.append(U)                                     # A_k = U diag(s) Vt  -> w = diag(s) Vt m
        Tis.append(Vt.T / s)
    W = np.concatenate(Ws, axis=1)
    L = np.linalg.norm(W, 2) ** 2
    w = np.zeros(6 * na)
    z = w.copy()
    tk = 1.0
    for _ in range(iters):
        g = W.T @ (W @ z - y)
        v = z - g / L
        wn = np.zeros_like(v)
        for k in range(na):
            vk = v[6 * k:6 * k + 6]
            nv = np.linalg.norm(vk)
            if nv > lam / L:
                wn[6 * k:6 * k + 6] = (1 - lam / (L * nv)) * vk
        t2 = (1 + np.sqrt(1 + 4 * tk * tk)) / 2
        z = wn + (tk - 1) / t2 * (wn - w)
        w, tk = wn, t2
    m = np.concatenate([Tis[k] @ w[6 * k:6 * k + 6] for k in range(na)])
    return m


def J_explicit(A, y, m, lam, na):
    r = y - A @ m
    return 0.5 * r @ r + lam * sum(np.linalg.norm(A[:, 6 * k:6 * k + 6] @ m[6 * k:6 * k + 6]) for k in range(na))


def main():
    L = []
    P_ = L.append
    ok_all = True
    t_start = time.time()
    P_(f"# Validacion de task2e_methods.py (semilla 999; {R2.C.kernel_info()})")
    rng = np.random.default_rng(20260927)
    # ------------------------------------------------------------------ V1
    net, PD, STA, statics, nclut, corr, grid, Gs, Qs = setup("ladder_6")
    k = np.arange(len(E.X.AC))
    d_ac = float(np.abs(E.acc((k - E.AC_HALF) * E.DT) - E.X.AC).max())
    Xs = np.column_stack([rng.uniform(-20, 20, 5), rng.uniform(-20, 20, 5), rng.uniform(1, 20, 5)])
    Gv = E.rad_op(STA, Xs)
    d_g = max(float(np.abs(Gv[a].reshape(-1, 6) - radiation_operator(Xs[a], STA)).max()) for a in range(5))
    W = windows(STA, statics, nclut, "ladder_6", n=4)
    srcs, taus, D, MFn, sd, tr = W[1]
    win = E.Win(D, MFn, sd)
    Ti = rng.integers(50, E.NT - 50, size=(3, 2, net.ns))
    v, _ = win.mfc(Ti * E.DT)
    d_mf = float(np.abs(v - MFn[np.arange(net.ns)[None, None, :], :, Ti]).max())
    # Gram continua en tiempos enteros = PKDict
    e0, I0, _h = PD.energy(MFn, True)
    picks = np.argsort(e0.ravel())[::-1][[0, 50, 400]]
    atoms = []
    Pint = []
    for jj in picks:
        n, t = divmod(int(jj), PD.T)
        atoms.append((n, t, I0[n, t].copy()))
        Pint.append(E.init_atom(net, n, t, I0[n, t], clip=False))
    Pint = np.array(Pint)
    Kc, bc = E.gram(net, win, Pint)
    Kd = PD.gram(atoms)
    bd = PD.corr_atoms(MFn, atoms).ravel()
    d_K = float(np.abs(Kc - Kd).max() / np.abs(Kd).max())
    d_b = float(np.abs(bc - bd).max() / np.abs(bd).max())
    ok1 = max(d_ac, d_g, d_mf) < 1e-12 and max(d_K, d_b) < 1e-12
    ok_all &= ok1
    P_("")
    P_("## V1 identidades")
    P_(f"- max |AC continuo - AC| en las muestras {d_ac:.1e}; max |rad_op - radiation_operator| {d_g:.1e}; "
       f"max |MFc(n DT) - MFn[n]| {d_mf:.1e}; Gram y correlacion continuas en tiempos enteros contra PKDict (3 atomos de "
       f"pico): dif. rel. max K {d_K:.1e}, b {d_b:.1e} -> {'correcto' if ok1 else 'FALLA'}")
    # ------------------------------------------------------------------ V2
    P_("")
    P_("## V2 gradientes contra diferencias finitas centrales (paso 1e-5 km, 1e-6 s)")
    errs = []
    for trial in range(6):
        cell = ("ladder_4", "ladder_6", "unilateral_c40")[trial % 3]
        netc, PDc, STAc, stc, ncl, _c, *_ = setup(cell)
        Wc = windows(STAc, stc, ncl, cell, n=2)
        _s, _t, Dc, MFc_, sdc, _tr = Wc[trial % 2]
        wc = E.Win(Dc, MFc_, sdc)
        na = 1 + trial % 3
        P = np.column_stack([rng.uniform(-15, 15, na), rng.uniform(-15, 15, na), rng.uniform(3, 15, na),
                             rng.uniform(5.5, 7.5, na), rng.uniform(-0.2, 0.2, (na, 2 * netc.ns))])
        m = rng.normal(size=(na, 6)) * 3.0
        steps = np.array([1e-5] * 3 + [1e-6] * (1 + 2 * netc.ns))
        for lam in (0.0, 2.5):
            J, g = E.joint_obj(netc, wc, P, m, lam)
            gf = fd_grad(lambda z: E.joint_obj(netc, wc, z.reshape(na, -1), m, lam, grad=False), P.ravel(),
                         np.tile(steps, na))
            rel = float(np.linalg.norm(gf - g.ravel()) / np.linalg.norm(gf))
            errs.append(rel)
            P_(f"- {cell} {na} atomo(s), lambda {lam}: J {J:.4f}, ||g_fd - g|| / ||g_fd|| = {rel:.1e}")
        a_old, T_old = E.pred_old(netc, wc, P[:-1], m[:-1])
        p = P[-1]
        v0, g, _mu = E.refine_obj(netc, wc, p, a_old, T_old)
        gf = fd_grad(lambda z: E.refine_obj(netc, wc, z, a_old, T_old, grad=False)[0], p, steps)
        rel = float(np.linalg.norm(gf - g) / np.linalg.norm(gf))
        errs.append(rel)
        P_(f"- {cell} refinamiento con {na - 1} atomo(s) activo(s): -E/2 {v0:.4f}, error relativo {rel:.1e}")
    ok2 = max(errs) < 1e-6
    ok_all &= ok2
    P_(f"- maximo error relativo {max(errs):.1e} -> {'correcto' if ok2 else 'FALLA'}")
    # ------------------------------------------------------------------ V3
    P_("")
    P_("## V3 paso convexo contra el diccionario explicito en el dominio de los datos")
    sqer = np.sqrt(M.E_R)
    y = win.Dn.ravel() * sqer
    worst = []
    for trial in range(3):
        na = 2 + trial
        P = np.column_stack([rng.uniform(-15, 15, na), rng.uniform(-15, 15, na), rng.uniform(3, 15, na),
                             rng.uniform(5.5, 7.5, na), rng.uniform(-0.25, 0.25, (na, 2 * net.ns))])
        A = explicit_dict(net, P, sqer)
        K, b = E.gram(net, win, P)
        dK = float(np.abs(A.T @ A - K).max() / np.abs(K).max())
        db = float(np.abs(A.T @ y - b).max() / np.abs(b).max())
        dyy = abs(y @ y - win.yy) / win.yy
        m0, inf0 = E.convex_step(net, win, P, np.zeros((na, 6)), 0.0)
        mls = np.linalg.lstsq(A, y, rcond=None)[0]
        dls = float(np.abs(m0.ravel() - mls).max() / np.abs(mls).max())
        c0 = np.array([np.sqrt(max(0.0, h @ np.linalg.pinv(K[6 * k:6 * k + 6, 6 * k:6 * k + 6]) @ h))
                       for k, h in enumerate(b.reshape(na, 6))])
        line = (f"- {na} atomos aleatorios: |A^T A - K| rel {dK:.1e}, |A^T y - b| rel {db:.1e}, ||y||^2 rel {dyy:.1e}; "
                f"lambda = 0: |m - lstsq(A, y)| rel {dls:.1e}")
        worst += [dK, db, dls]
        for fr in (0.05, 0.3, 0.8):
            lam = fr * float(c0.max())
            m1, inf1 = E.convex_step(net, win, P, np.zeros((na, 6)), lam)
            mf = fista_group(A, y, lam, na)
            J1 = J_explicit(A, y, m1.ravel(), lam, na)
            Jf = J_explicit(A, y, mf, lam, na)
            line += (f"; lambda = {fr} c_max: brecha rel {inf1['rel_gap']:.1e}, KKT {inf1['kkt']:.1e}, grupos nulos "
                     f"{inf1['n_zero']}, J(metodo) - J(FISTA) = {J1 - Jf:+.1e} (rel {(J1 - Jf) / abs(Jf):+.1e})")
            worst.append(max(0.0, (J1 - Jf) / abs(Jf)))
            worst.append(inf1["rel_gap"] * 1e-3)
            worst.append(inf1["kkt"] * 1e-2)
        m_small, _ = E.convex_step(net, win, P, np.zeros((na, 6)), 1e-9 * float(c0.max()))
        dsm = float(np.abs(m_small - m0).max() / np.abs(m0).max())
        line += f"; lambda = 1e-9 c_max contra lambda = 0: |dm| rel {dsm:.1e}"
        worst.append(dsm * 1e-2)
        P_(line)
    ok3 = max(worst) < 1e-8
    ok_all &= ok3
    P_(f"- peor discrepancia (ponderada) {max(worst):.1e} -> {'correcto' if ok3 else 'FALLA'}")
    # ------------------------------------------------------------------ V4
    P_("")
    P_("## V4 fuente unica fuera de la rejilla, sin ruido, desplazamientos inyectados |delta| en [0.1, 0.2] s (signo "
       "aleatorio)")
    P_("Error de desplazamientos: max |(delta_est - media) - (delta_true - media)| (el desplazamiento comun se confunde "
       "exactamente con t0: T = d/V + t0 + corr + delta); tambien max |T_est - T_true| de los tiempos de llegada. "
       "Criterio de recuperacion: error 3D < 0.5 km y desplazamientos < 0.02 s.")
    labs = [("lambda=0", {"lam": 0.0})] + [(f"f={f}", {"f": f}) for f in FGRID]
    v4 = {}
    basin = {}
    fixed = [("ladder_4", np.array([1.5, -1.5, 7.0]), 6.33), ("ladder_6", np.array([-7.4, 4.6, 11.0]), 6.02),
             ("unilateral_c40", np.array([10.3, 8.8, 6.1]), 5.87), ("ladder_6", np.array([4.4, -10.6, 9.3]), 7.11)]
    rs = np.random.default_rng(4242)
    rand_cases = [(c, np.array([rs.uniform(-15, 15), rs.uniform(-15, 15), rs.uniform(5, 13)]), float(rs.uniform(5.5, 7.5)))
                  for c in ("ladder_4", "ladder_6", "unilateral_c40") for _ in range(N_RAND)]
    nets = {c: setup(c) for c in ("ladder_4", "ladder_6", "unilateral_c40")}
    for ic, (cell, xs, t0s) in enumerate(fixed + rand_cases):
        netc, PDc, STAc = nets[cell][:3]
        Mt = rand_dc(np.random.default_rng([ic, 5]))
        mtrue = np.array([Mt[0, 0], Mt[1, 1], Mt[2, 2], Mt[0, 1], Mt[0, 2], Mt[1, 2]])
        sign = rs.choice([-1.0, 1.0], size=(2, netc.ns))
        dtrue = sign * rs.uniform(0.1, 0.2, size=(2, netc.ns))
        ptrue = np.concatenate([xs, [t0s], dtrue.ravel()])
        G, _dG, Tt, _d = netc.geom(ptrue[None])
        a = np.einsum("pics,s->pic", G[0], mtrue)
        scale = 12.0 / np.abs(a).max()
        mtrue = mtrue * scale
        a = a * scale
        tt = np.arange(E.NT) * E.DT
        Dn = np.zeros((netc.ns, 3, E.NT))
        for p in range(2):
            for i in range(netc.ns):
                Dn[i] += np.outer(a[p, i], E.ricker_d(tt - Tt[0, p, i])[0]) / M.E_R
        MFn, _sd = T.matched(Dn)
        wc = E.Win(Dn, MFn, 1.0)
        truth = [xs]
        rb_d = {}
        D2.bomp_pk_record(PDc, MFn, truth, diag=rb_d)
        n0 = rb_d["atoms"][0][0]
        eb = float(np.linalg.norm(netc.garr[n0] - xs))
        ebh = float(np.linalg.norm(netc.garr[n0][:2] - xs[:2]))

        def errs(p):
            el = float(np.linalg.norm(p[:3] - xs))
            dl = p[4:].reshape(2, -1)
            es = float(np.abs((dl - dl.mean()) - (dtrue - dtrue.mean())).max())
            _G, _dG, Te, _d = netc.geom(p[None])
            return el, es, float(np.abs(Te[0] - Tt[0]).max())
        line = (f"- {cell} fuente {np.round(xs, 2).tolist()} t0 {t0s:.2f}: BOMP-PK primer atomo en el nodo "
                f"{netc.garr[n0].tolist()} (error 3D {eb:.2f} km, epicentral {ebh:.2f} km)")
        for lab, kw in labs:
            dg = {}
            E.tssfw_record(netc, wc, truth, kmax=1, diag=dg, **kw)
            el, es, et = errs(dg["P"][0])
            v4[(ic, lab)] = (el, es, et, eb, ic < len(fixed))
            line += f"; TS-SFW {lab}: {el:.3f} km, {es:.4f} s, tiempos {et:.4f} s"
        # cuenca (prueba del optimizador local, no del metodo): refinamiento + deslizamiento (lambda = 0) desde
        # (A) posicion verdadera + 0.2-0.4 km, t0 verdadero, delta verdadero + U(-0.03, 0.03) s; (B) posicion + 0.5-1 km,
        # delta verdadero; (C) la verdad exacta
        u = rs.normal(size=3)
        u /= np.linalg.norm(u)
        starts = {"A": np.concatenate([xs + u * rs.uniform(0.2, 0.4), [t0s],
                                       (dtrue + rs.uniform(-0.03, 0.03, dtrue.shape)).ravel()]),
                  "B": np.concatenate([xs + u * rs.uniform(0.5, 1.0), [t0s], dtrue.ravel()]),
                  "C": ptrue.copy()}
        for kb, pst in starts.items():
            p1, _inf = E.refine_atom(netc, wc, pst, np.zeros((0, 2, netc.ns, 3)), np.zeros((0, 2, netc.ns)))
            Pb = p1[None]
            mb, _c = E.convex_step(netc, wc, Pb, np.zeros((1, 6)), 0.0)
            for _alt in range(E.MAX_ALT):
                Pb, li = E.slide_step(netc, wc, Pb, mb, 0.0)
                mb, _c = E.convex_step(netc, wc, Pb, mb, 0.0)
                if li["J0"] - li["J1"] <= 1e-10 * max(1.0, abs(li["J0"])):
                    break
            basin.setdefault(kb, []).append(errs(Pb[0]))
        if ic < len(fixed):
            P_(line)
    rec_ = lambda v: v[0] < 0.5 and v[1] < 0.02  # noqa: E731
    P_(f"- {N_RAND} fuentes aleatorias por celda (x, y en U[-15, 15], z en U[5, 13], t0 en U[5.5, 7.5], doble par "
       f"aleatorio), fraccion recuperada (error 3D < 0.5 km y desplazamientos < 0.02 s) y mediana del error 3D:")
    for lab, _kw in labs:
        vv = [v for (ic, l_), v in v4.items() if l_ == lab and not v[4]]
        P_(f"  - TS-SFW {lab}: recuperadas {sum(rec_(v) for v in vv)}/{len(vv)}; error 3D mediana "
           f"{np.median([v[0] for v in vv]):.2f} km, P90 {np.percentile([v[0] for v in vv], 90):.2f} km; tiempos "
           f"de llegada con error < 0.02 s en {sum(v[2] < 0.02 for v in vv)}/{len(vv)}")
    vb = [v[3] for (ic, l_), v in v4.items() if l_ == "lambda=0" and not v[4]]
    P_(f"  - BOMP-PK (primer atomo en la rejilla): error 3D < 0.5 km en {sum(x < 0.5 for x in vb)}/{len(vb)}; "
       f"mediana {np.median(vb):.2f} km")
    for kb, lab in (("C", "desde la verdad exacta"), ("A", "posicion + 0.2-0.4 km, delta verdadero + U(-0.03, 0.03) s"),
                    ("B", "posicion + 0.5-1 km, delta verdadero")):
        vv = basin[kb]
        P_(f"- Cuenca del optimizador local (refinamiento y deslizamiento, lambda = 0), {lab}: recuperadas "
           f"{sum(rec_(v) for v in vv)}/{len(vv)}; error 3D mediana {np.median([v[0] for v in vv]):.3f} km, max "
           f"{max(v[0] for v in vv):.3f} km; desplazamientos max {max(v[1] for v in vv):.4f} s")
    ok4_basin = (all(rec_(v) and v[0] < 1e-3 for v in basin["C"])
                 and sum(rec_(v) for v in basin["A"]) >= 0.8 * len(basin["A"]))
    ok4_fixed = all(rec_(v) for (ic, l_), v in v4.items() if v[4] and l_ == "lambda=0")
    P_(f"- Criterio preespecificado por el encargo (TS-SFW recupera < 0.5 km y < 0.02 s donde BOMP-PK no) en los 4 casos "
       f"fijos: {'si' if ok4_fixed else 'NO'}; en las fuentes aleatorias ver las fracciones. El optimizador local converge "
       f"en su cuenca: {'si' if ok4_basin else 'NO'}.")
    V4_SUMMARY.update({"fixed_ok": ok4_fixed, "basin_ok": ok4_basin,
                       "rand": {lab: [sum(rec_(v) for (ic, l_), v in v4.items() if l_ == lab and not v[4]),
                                      sum(1 for (ic, l_), v in v4.items() if l_ == lab and not v[4])]
                                for lab, _ in labs}})
    ok_all &= ok4_basin
    # ------------------------------------------------------------------ V5
    P_("")
    P_("## V5 sin deslizamiento: primer atomo = BOMP-PK; oraculo = MT-FW-NOCAP+SC")
    n_same, n_tot, dm_max, n_rec_eq, n_nocap = 0, 0, 0.0, 0, 0
    lines = []
    for cell in ("ladder_4", "ladder_6", "unilateral_c40"):
        for regime in ("thin", "strong"):
            netc, PDc, STAc, stc, ncl, corrc, gridc, Gsc, Qsc = setup(cell, regime)
            for (srcs, taus, D, MFn, sd, tr) in windows(STAc, stc, ncl, cell, regime, n=4):
                wc = E.Win(D, MFn, sd)
                n_tot += 1
                dgb = {}
                rb = D2.bomp_pk_record(PDc, MFn, tr, diag=dgb)
                dg = {}
                E.tssfw_record(netc, wc, tr, lam=0.0, kmax=1, refine=False, slide=False, clip_init=False, diag=dg)
                it = dg["iters"][0]
                nb, tb, Ib = dgb["atoms"][0]
                same = (it["node"] == nb and it["tau"] == tb and
                        it["I_shift"] == (Ib - netc.BD.idx[nb, tb]).ravel().tolist())
                if same:
                    n_same += 1
                    # amplitudes del primer atomo de BOMP-PK (antes de entrar el segundo)
                    K1 = PDc.gram([dgb["atoms"][0]])
                    b1 = PDc.corr_atoms(MFn, [dgb["atoms"][0]]).ravel()
                    u1 = np.linalg.lstsq(K1, b1, rcond=None)[0]
                    dm_max = max(dm_max, float(np.abs(dg["m"][0] - u1).max() / np.abs(u1).max()))
                dg2 = {}
                r2 = E.tssfw_record(netc, wc, tr, lam=0.0, refine=False, slide=False, clip_init=False, stat="bomp",
                                    PD=PDc, diag=dg2)
                eq = (len(r2) == len(rb) and all(h1 == h2 and abs(s1 - s2) <= 1e-9 * abs(s2)
                                                  for (s1, h1), (s2, h2) in zip(r2, rb)))
                n_rec_eq += eq
                if not eq:
                    lines.append(f"  registro distinto {cell} {regime}: {r2} vs {rb}")
                rn = D2.fw_record_v(STAc, MFn, gridc, T.DEFAULT_TAUS, Gsc, Qsc, tr, corrc, peak=True, cap=False, kmax=1)
                n_nocap += (abs(rn[0][0] - it["s"]) <= 1e-9 * abs(rn[0][0]))
    ok5 = n_rec_eq == n_tot and n_nocap == n_tot and dm_max < 1e-9
    ok_all &= ok5
    P_(f"- {n_tot} ventanas (4 por celda y regimen, red 0): el oraculo (MT-FW sin tope) elige el mismo bloque y las "
       f"mismas muestras de pico que BOMP-PK en {n_same}/{n_tot} (el resto difiere solo por la metrica: G^T G contra "
       f"K(I) con el solape P/S); en esos, max dif. rel. de las amplitudes del primer atomo {dm_max:.1e}. Con el "
       f"estadistico de BOMP-PK en el oraculo, el registro completo (5 atomos, reajuste por minimos cuadrados con tiempos "
       f"fijos) es el de BOMP-PK en {n_rec_eq}/{n_tot}. Primer valor del oraculo = MT-FW-NOCAP+SC en {n_nocap}/{n_tot} "
       f"-> {'correcto' if ok5 else 'FALLA'}")
    L.extend(lines)
    P_("")
    P_(f"Resultado global de las comprobaciones de implementacion (V1, V2, V3, cuenca de V4, V5): "
       f"{'CORRECTAS' if ok_all else 'HAY FALLAS'}; criterio de recuperacion de V4 en los casos fijos: "
       f"{'se cumple' if V4_SUMMARY.get('fixed_ok') else 'NO se cumple'} ({time.time() - t_start:.0f} s)")
    txt = "\n".join(L)
    print(txt)
    with open(OUT, "w") as fh:
        fh.write(txt + "\n")
    return ok_all


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
