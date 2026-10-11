"""
task2e_methods.py - metodos de la Tarea 2E (PREREGISTRO_1j_1k.txt, seccion "Tarea 2E"): BLASSO de grupos con
desplazamientos temporales por estacion y fase resuelto fuera de la rejilla por Frank-Wolfe deslizante (TS-SFW+SC) y su
control con lambda = 0 (BOMP-PK-SLIDE+SC).  Reutiliza task2d_methods (peak_table, PKDict, bomp_pk_record),
task2c_methods (E_R, EIG_FLOOR, solve_restricted, kkt_restricted, hit_of), certified_fw (certified_oracle, Lattice,
stack_Q, _exclusion_update), task1x_mtglrt (radiation_operator, basis, AC) y task0_common (VP, VS, DT, NT, rejilla).

MODELO CONTINUO (dominio del filtro adaptado, MFn = MF / sd, MF = correlacion de D con el nucleo Ricker KER de +-1 s).
Un atomo tiene theta = (x, y, z, t0) y desplazamientos delta (2, ns) (fase P/S, estacion), acotados a +-0.25 s.  Su
tiempo en la estacion i y la fase b es  T_bi = d_i(x) / V_b + t0 + corr_ib + delta_bi  (continuo, sin redondeo a
muestras; corr = correcciones de estacion de 2A) y su prediccion es  a_bi (x) AC(t - T_bi),  a_bi = G_bi(x) m (G =
task1x_mtglrt.radiation_operator, evaluado en x continuo) y AC(s) = F(s) / F(0), F(s) = sum_{|l|<=16} ricker(s + l DT)
KER_l: la funcion continua cuyas muestras son task1x_mtglrt.AC (identica hasta el redondeo).  Es la imagen exacta por el
filtro adaptado del atomo de datos G m (x) ricker(t - T) / sqrt(E_R) (2C), y la correlacion de los datos con el atomo es
la traza del filtro adaptado evaluada en el tiempo continuo T (interpolacion de banda limitada exacta):
    MFc_i(T) = sum_n Dn_i[n] ricker(n DT - T) (|n DT - T| <= 1 s, Dn = D / sd)   (= MFn_i[n] si T = n DT).
Con y los datos (||y||^2 = E_R ||Dn||^2, como GL-FW de 2C), para atomos k = 1..K con amplitudes m_k:
    b_k  = sum_bi G_kbi^T MFc_i(T_kbi)                                  (A^T y)
    K_kl = sum_i sum_{p,q} G_kpi^T G_lqi AC(T_kpi - T_lqi)              (A^T A, Gram con AC)
    J    = 1/2 ||y - sum_k A_k m_k||^2 + lambda sum_k ||m_k||_G,  ||m_k||_G = sqrt(m_k^T K_kk m_k)  (norma de grupo
           blanqueada por bloque, la de GL-FW de 2C) = 1/2 ||y||^2 - m^T b + 1/2 m^T K m + lambda Omega(m).
Gradientes analiticos respecto de (x, y, z, t0, delta) de todos los atomos (derivada de G por paso complejo, exacta a
precision de maquina; derivadas de MFc y AC por las de ricker), verificados contra diferencias finitas en
task2e_validate.py.

TS-SFW+SC (en cada iteracion k <= 5):
  (a) oraculo: maximo sobre los nodos de la rejilla (excluidos los que estan a <= 5 km epicentrales de la posicion ACTUAL
      de un atomo aceptado) y los origenes DEFAULT_TAUS de la energia ||Q(g)^T v||^2 / 2 de los vectores de tres
      componentes del residuo muestreado en el pico de energia de cada estacion y fase dentro de +-4 muestras (0.25 s)
      de la muestra predicha (con correcciones): el estadistico de MT-FW sin tope (task2d MT-FW-NOCAP+SC), vectorizado
      con task2d_methods.peak_table.  Valor del oraculo s_k.
  (b) el atomo nuevo arranca en el nodo y el origen elegidos con delta_bi = clip(tpk_bi DT - T_bi, -0.25, 0.25) (el pico
      encontrado) y se refina maximizando E(theta, delta) = v^T G (G^T G)^+ G^T v / 2, v = residuo continuo en T_bi
      (MFc menos las predicciones continuas de los atomos activos), con L-BFGS-B acotado (x, y en [-21, 21] km, la
      extension de la rejilla; z en [1, 20] km; t0 en [4.5, 8.5] s; delta en [-0.25, 0.25] s), a lo sumo 50 iteraciones.
  (c) paso deslizante: alternancia de (i) el paso convexo en las amplitudes de todos los atomos con theta fijo (lasso de
      grupos restringido en coordenadas blanqueadas, task2c_methods.solve_restricted, con brecha de dualidad relativa
      <= 1e-6 comprobada y KKT registrado; con lambda = 0, minimos cuadrados) y (ii) L-BFGS-B sobre (theta, delta) de
      todos los atomos con m fijo (J completo, incluida la penalizacion, a lo sumo 50 iteraciones); a lo sumo 5
      alternancias (termina antes si el paso local no baja J en mas de 1e-10 |J|); al final, un paso convexo.
  (d) residuo muestreado = MFn menos la prediccion conjunta (AC continuo evaluado en las muestras, +-48 muestras).
  lambda = f * s_1 (s_1 = valor del oraculo en la primera iteracion de la ventana, en las unidades del estadistico de
  MT-FW, las del puntaje).  Registro: minimo acumulado de s_k; acierto con la posicion final (x, y) de cada atomo tras la
  ultima iteracion (task2c_methods.hit_of, 5 km epicentrales).
BOMP-PK-SLIDE+SC: lo mismo con lambda = 0.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys
import time

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import task1x_mtglrt as X  # noqa: E402
import task2c_methods as M  # noqa: E402
import task2d_methods as D2  # noqa: E402
from certified_fw import Lattice, _exclusion_update, certified_oracle, stack_Q  # noqa: E402
from task0_common import AP, AS, DEFAULT_TAUS, DT, FC, NT, VP, VS  # noqa: E402
from task1o_consistent import _KER, TOL_KM  # noqa: E402
from task1w_mech import CP, CS  # noqa: E402

SHIFT_MAX = 0.25
XY_BOUND = 21.0
Z_BOUND = (1.0, 20.0)
T0_BOUND = (float(DEFAULT_TAUS[0]) - 0.5, float(DEFAULT_TAUS[-1]) + 0.5)
MAXIT_LOCAL = 50
MAX_ALT = 5
GAP_TOL = 1e-6
TSC = 0.2                       # escala interna de las variables de tiempo en L-BFGS-B (s por unidad)
H_CS = 1e-30                    # paso complejo
KER = np.asarray(_KER, float)
LK = (np.arange(len(KER)) - len(KER) // 2) * DT
F0 = float(np.sum(KER ** 2))
LW = np.arange(-17, 19)         # ventana de muestras de MFc (se enmascara a |n DT - T| <= 1 s)
AC_HALF = X.AC_HALF
LAC = np.arange(-AC_HALF, AC_HALF + 2)
VINV = np.array([1.0 / VP, 1.0 / VS])


# ============================================================================================ ondicula
def ricker_d(t):
    """ricker(t) y su derivada."""
    x = (np.pi * FC * t) ** 2
    e = np.exp(-x)
    return (1.0 - 2.0 * x) * e, (2.0 * x - 3.0) * e * 2.0 * (np.pi * FC) ** 2 * t


def acc(s, deriv=False):
    """AC continuo F(s) / F(0) (y su derivada)."""
    s = np.asarray(s, float)
    r, rd = ricker_d(s[..., None] + LK)
    v = r @ KER / F0
    if not deriv:
        return v
    return v, rd @ KER / F0


# ============================================================================================ operador
def rad_op(STA, Xs):
    """task1x_mtglrt.radiation_operator vectorizado sobre fuentes Xs (na, 3), real o complejo: (na, 2, ns, 3, 6)."""
    v = STA[None, :, :] - Xs[:, None, :]
    d = np.sqrt(np.sum(v * v, axis=-1))
    n = v / d[..., None]
    nx, ny, nz = n[..., 0], n[..., 1], n[..., 2]
    c = np.stack([nx * nx, ny * ny, nz * nz, 2 * nx * ny, 2 * nx * nz, 2 * ny * nz], axis=-1)
    gP = n[..., :, None] * c[..., None, :] * (AP / np.sqrt(CP)) / d[..., None, None]
    z = np.zeros_like(nx)
    B = np.stack([np.stack([nx, z, z, ny, nz, z], -1), np.stack([z, ny, z, nx, z, nz], -1),
                  np.stack([z, z, nz, z, nx, ny], -1)], -2)
    Pm = np.eye(3)[None, None] - n[..., :, None] * n[..., None, :]
    gS = np.einsum("...ij,...jk->...ik", Pm, B) * (AS / np.sqrt(CS)) / d[..., None, None]
    return np.stack([gP, gS], axis=1)


def rad_op_cs(STA, Xs):
    """G (na, 2, ns, 3, 6) y dG/dx_j (na, 3, 2, ns, 3, 6) por paso complejo."""
    na = len(Xs)
    Xc = Xs[None, :, :].astype(complex) + 1j * H_CS * np.eye(3)[:, None, :]
    Gc = rad_op(STA, Xc.reshape(3 * na, 3)).reshape(3, na, 2, len(STA), 3, 6)
    return Gc[0].real.copy(), np.moveaxis(Gc.imag / H_CS, 0, 1)


# ============================================================================================ red y ventana
class Net:
    """operadores por red: estaciones, correcciones, rejilla, bases Q, indices predichos (BlockDict/PKDict de 2D)."""

    def __init__(self, STA, corr, grid, taus, Gs, Qs, BD):
        self.STA = np.asarray(STA, float)
        self.ns = len(self.STA)
        self.corrT = np.asarray(corr, float).T.copy()                 # (2, ns)
        self.corr = np.asarray(corr, float)
        self.garr = np.asarray(grid, float)
        self.taus = np.asarray(taus, float)
        self.Gs, self.Qs = Gs, Qs
        self.Qst = stack_Q(Qs)                                        # (N, 6ns, r)
        self.BD = BD
        self.N, self.T = len(self.garr), len(self.taus)
        self.st4 = np.arange(self.ns)[None, None, None, :]
        self.nv = 4 + 2 * self.ns
        lo = [-XY_BOUND, -XY_BOUND, Z_BOUND[0], T0_BOUND[0]] + [-SHIFT_MAX] * (2 * self.ns)
        hi = [XY_BOUND, XY_BOUND, Z_BOUND[1], T0_BOUND[1]] + [SHIFT_MAX] * (2 * self.ns)
        self.lo, self.hi = np.array(lo), np.array(hi)
        self.scale = np.array([1.0, 1.0, 1.0] + [TSC] * (1 + 2 * self.ns))

    def geom(self, P):
        """P (na, 4 + 2ns) -> G, dG, T (na, 2, ns), dT/dx (na, 2, ns, 3)."""
        Xs = P[:, :3]
        na = len(P)
        v = self.STA[None] - Xs[:, None]
        d = np.sqrt(np.sum(v * v, axis=-1))
        T = d[:, None, :] * VINV[None, :, None] + P[:, 3][:, None, None] + self.corrT[None] + \
            P[:, 4:].reshape(na, 2, self.ns)
        dTdX = -(v / d[..., None])[:, None, :, :] * VINV[None, :, None, None]
        G, dG = rad_op_cs(self.STA, Xs)
        return G, dG, T, dTdX


class Win:
    """datos de una ventana: MFn (muestras), Dn = D / sd (para MFc continuo), ||y||^2."""

    def __init__(self, D, MFn, sd):
        self.Dn = np.asarray(D, float) / sd
        self.MFn = MFn
        self.yy = M.E_R * float(np.sum(self.Dn ** 2))
        self.ns = D.shape[0]

    def mfc(self, T):
        """MFc y su derivada en tiempos T (..., 2, ns) (estacion = ultimo eje): (..., 2, ns, 3) cada una."""
        sh = T.shape
        ist = np.broadcast_to(np.arange(self.ns), sh).ravel()
        Tf = T.ravel()
        n0 = np.floor(Tf / DT).astype(int)
        idx = n0[:, None] + LW[None, :]
        s = idx * DT - Tf[:, None]
        ok = (idx >= 0) & (idx < NT) & (np.abs(s) <= 1.0 + 1e-9)
        r, rd = ricker_d(s)
        r = np.where(ok, r, 0.0)
        rd = np.where(ok, rd, 0.0)
        seg = self.Dn[ist[:, None], :, np.clip(idx, 0, NT - 1)]            # (Mt, W, 3)
        val = np.einsum("mlc,ml->mc", seg, r)
        der = -np.einsum("mlc,ml->mc", seg, rd)
        return val.reshape(sh + (3,)), der.reshape(sh + (3,))


# ============================================================================================ oraculo
def oracle(net, R, excluded, stat="mtfw", PD=None):
    """maximo de ||Q^T v||^2 / 2 (v en el pico, sin tope) sobre nodos no excluidos y origenes: (s, n, t, I (2, ns)).
    stat='bomp' usa la energia de BOMP-PK (h^T K(I)^+ h / 2, PKDict.energy), solo para validacion."""
    if stat == "bomp":
        e, I, _h = PD.energy(R, True)
        e = e / 2.0
    else:
        Kt = D2.peak_table(R)
        I = Kt[net.st4, net.BD.idx]                                          # (N, T, 2, ns)
        v = R[np.arange(net.ns)[None, None, None, :], :, I].reshape(net.N, net.T, -1)
        e = np.sum(np.matmul(v, net.Qst) ** 2, axis=-1) / 2.0
    e = np.where(excluded[:, None], -np.inf, e)
    j = int(np.argmax(e))
    n, t = divmod(j, net.T)
    return float(e[n, t]), n, t, I[n, t].copy()


def excl_from(net, P, min_sep):
    mask = np.zeros(net.N, bool)
    for p in P:
        mask = _exclusion_update(mask, net.garr, np.asarray(p[:3]), min_sep)
    return mask


# ============================================================================================ objetivos
def joint_obj(net, win, P, m, lam, grad=True):
    """J - ||y||^2/2 = -m^T b + 1/2 m^T K m + lam Omega y su gradiente respecto de P (na, nv)."""
    na = len(P)
    G, dG, T, dTdX = net.geom(P)
    a = np.einsum("kpicx,kx->kpic", G, m)
    mf, mfd = win.mfc(T)
    Lag = T[:, :, None, None, :] - T[None, None, :, :, :]                    # (k, p, l, q, i)
    A0, A1 = acc(Lag, deriv=True)
    S0 = np.einsum("kplqi,lqic->kpic", A0, a)
    Jfit = -np.sum(a * mf) + 0.5 * np.sum(a * S0)
    kk = np.arange(na)
    A0d = A0[kk, :, kk]                                                      # (k, p, q, i)
    A1d = A1[kk, :, kk]
    Qk = np.einsum("kpqi,kpic,kqic->k", A0d, a, a)
    sq = np.sqrt(np.maximum(Qk, 0.0))
    J = float(Jfit + lam * np.sum(sq))
    if not grad:
        return J
    S1 = np.einsum("kplqi,lqic->kpic", A1, a)
    gA = -(mf - S0)
    gT = -np.sum(a * (mfd - S1), axis=-1)
    if lam > 0:
        pos = sq > 0
        coef = np.where(pos, lam / np.where(pos, sq, 1.0), 0.0)
        gA = gA + coef[:, None, None, None] * np.einsum("kpqi,kqic->kpic", A0d, a)
        gT = gT + coef[:, None, None] * np.einsum("kpqi,kpic,kqic->kpi", A1d, a, a)
    gX = np.einsum("kpic,kjpicx,kx->kj", gA, dG, m) + np.einsum("kpi,kpij->kj", gT, dTdX)
    g = np.concatenate([gX, gT.sum(axis=(1, 2))[:, None], gT.reshape(na, -1)], axis=1)
    return J, g


def pred_old(net, win, P, m):
    """a (na, 2, ns, 3) y T de los atomos activos (fijos durante el refinamiento)."""
    if len(P) == 0:
        return np.zeros((0, 2, net.ns, 3)), np.zeros((0, 2, net.ns))
    G, _dG, T, _d = net.geom(P)
    return np.einsum("kpicx,kx->kpic", G, m), T


def refine_obj(net, win, p, a_old, T_old, grad=True):
    """-E / 2 del atomo p (nv,) sobre el residuo continuo (atomos activos fijos) y su gradiente; devuelve tambien mu."""
    G, dG, T, dTdX = net.geom(p[None])
    G, dG, T, dTdX = G[0], dG[0], T[0], dTdX[0]
    mf, mfd = win.mfc(T)
    if len(a_old):
        Lag = T[:, None, None, :] - T_old[None]                              # (p, l, q, i)
        A0, A1 = acc(Lag, deriv=True)
        v = mf - np.einsum("plqi,lqic->pic", A0, a_old)
        vd = mfd - np.einsum("plqi,lqic->pic", A1, a_old)
    else:
        v, vd = mf, mfd
    h = np.einsum("picx,pic->x", G, v)
    S = np.einsum("picx,picy->xy", G, G)
    mu = np.linalg.lstsq(S, h, rcond=None)[0]
    E = float(h @ mu)
    if not grad:
        return -E / 2.0, mu
    Gmu = np.einsum("picx,x->pic", G, mu)
    gT = -np.sum(Gmu * vd, axis=-1)                                          # d(-E/2)/dT
    gX = -np.einsum("pic,jpicx,x->j", v - Gmu, dG, mu) + np.einsum("pi,pij->j", gT, dTdX)
    g = np.concatenate([gX, [gT.sum()], gT.ravel()])
    return -E / 2.0, g, mu


# ============================================================================================ paso convexo
def gram(net, win, P):
    G, _dG, T, _d = net.geom(P)
    mf, _mfd = win.mfc(T)
    Lag = T[:, :, None, None, :] - T[None, None, :, :, :]
    A0 = acc(Lag)
    na = len(P)
    K = np.einsum("kpicx,lqicy,kplqi->kxly", G, G, A0, optimize=True).reshape(6 * na, 6 * na)
    K = 0.5 * (K + K.T)
    b = np.einsum("kpicx,kpic->kx", G, mf).ravel()
    return K, b


def whiten_block(Kb):
    ev, V = np.linalg.eigh(Kb)
    keep = ev > M.EIG_FLOOR * ev[-1]
    sq = np.where(keep, np.sqrt(np.where(keep, ev, 1.0)), 0.0)
    isq = np.where(keep, 1.0 / np.where(keep, np.sqrt(np.where(keep, ev, 1.0)), 1.0), 0.0)
    return np.einsum("j,kj->jk", sq, V), np.einsum("kj,j->kj", V, isq)


def convex_step(net, win, P, m0, lam):
    """amplitudes optimas con theta fijo. lam = 0: minimos cuadrados (lstsq). lam > 0: lasso de grupos blanqueado
    (solve_restricted) hasta brecha de dualidad relativa <= GAP_TOL.  Devuelve m (na, 6) y diagnosticos."""
    na = len(P)
    K, b = gram(net, win, P)
    if lam == 0:
        u = np.linalg.lstsq(K, b, rcond=None)[0]
        ne = float(np.abs(K @ u - b).max() / max(1e-300, np.abs(b).max()))
        return u.reshape(na, 6), {"gap": 0.0, "rel_gap": 0.0, "kkt": ne, "n_zero": 0, "resolves": 0}
    Tw, Ti = zip(*[whiten_block(K[6 * k:6 * k + 6, 6 * k:6 * k + 6]) for k in range(na)])
    Tw, Ti = np.stack(Tw), np.stack(Ti)
    K4 = K.reshape(na, 6, na, 6)
    Kp = np.einsum("iak,iajb,jbl->ikjl", Ti, K4, Ti, optimize=True).reshape(6 * na, 6 * na)
    Kp = np.ascontiguousarray(0.5 * (Kp + Kp.T))
    bp = np.einsum("iak,ia->ik", Ti, b.reshape(na, 6)).ravel()
    w = np.einsum("kjl,kl->kj", Tw, m0).ravel().copy()
    res = 0
    tol = M.RESTR_TOL
    while True:
        M.solve_restricted(Kp, bp, w, lam, tol=tol)
        u = np.einsum("ikl,il->ik", Ti, w.reshape(na, 6))
        uf = u.ravel()
        h = b - K @ uf
        c = np.sqrt(np.maximum(np.einsum("kj,kj->k", np.einsum("kaj,ka->kj", Ti, h.reshape(na, 6)),
                                         np.einsum("kaj,ka->kj", Ti, h.reshape(na, 6))), 0.0))
        cmax = float(c.max())
        s = min(1.0, lam / cmax) if cmax > 0 else 1.0
        Om = float(np.sum(np.linalg.norm(w.reshape(na, 6), axis=1)))
        rr = win.yy - 2.0 * uf @ b + uf @ K @ uf
        gap = 0.5 * (1.0 - s) ** 2 * rr + lam * Om - s * (uf @ h)
        dec = -(0.5 * uf @ K @ uf - uf @ b + lam * Om)
        rel = gap / dec if dec > 0 else (0.0 if gap <= 1e-12 * max(1.0, abs(rr)) else np.inf)
        if rel <= GAP_TOL or res >= 3:
            break
        res += 1
        tol *= 1e-2
    kkt = M.kkt_restricted(Kp, bp, w, lam)
    nz = int(np.sum(np.linalg.norm(w.reshape(na, 6), axis=1) == 0))
    return u, {"gap": float(gap), "rel_gap": float(rel), "kkt": float(kkt), "n_zero": nz, "resolves": res,
               "cmax_over_lam": cmax / lam}


# ============================================================================================ pasos locales
def _lbfgs(fun, x0, lo, hi, scale, maxiter):
    z0 = np.clip(x0, lo, hi) / scale
    bnds = list(zip(lo / scale, hi / scale))

    def f(z):
        v, g = fun(z * scale)
        return v, g * scale
    r = minimize(f, z0, jac=True, method="L-BFGS-B", bounds=bnds, options={"maxiter": maxiter})
    return r.x * scale, r


def refine_atom(net, win, p0, a_old, T_old, maxiter=MAXIT_LOCAL):
    def fun(p):
        v, g, _mu = refine_obj(net, win, p, a_old, T_old)
        return v, g
    p, r = _lbfgs(fun, p0, net.lo, net.hi, net.scale, maxiter)
    return p, {"nit": int(r.nit), "E0": -2 * fun(np.clip(p0, net.lo, net.hi))[0], "E1": -2 * float(r.fun),
               "status": int(r.status)}


def slide_step(net, win, P, m, lam, maxiter=MAXIT_LOCAL):
    na = len(P)
    lo, hi, sc = np.tile(net.lo, na), np.tile(net.hi, na), np.tile(net.scale, na)

    def fun(z):
        J, g = joint_obj(net, win, z.reshape(na, net.nv), m, lam)
        return J, g.ravel()
    J0 = fun(np.clip(P.ravel(), lo, hi))[0]
    z, r = _lbfgs(fun, P.ravel(), lo, hi, sc, maxiter)
    return z.reshape(na, net.nv), {"nit": int(r.nit), "J0": float(J0), "J1": float(r.fun), "status": int(r.status)}


# ============================================================================================ prediccion muestreada
def predict_sampled(net, win, P, m):
    """prediccion conjunta en las muestras (ns, 3, NT), AC continuo en +-AC_HALF muestras del tiempo de cada atomo."""
    Pr = np.zeros_like(win.MFn)
    if len(P) == 0:
        return Pr
    G, _dG, T, _d = net.geom(P)
    a = np.einsum("kpicx,kx->kpic", G, m)
    n0 = np.floor(T / DT).astype(int)
    idx = n0[..., None] + LAC                                                 # (na, 2, ns, W)
    w = acc(idx * DT - T[..., None])
    ok = (idx >= 0) & (idx < NT)
    w = np.where(ok, w, 0.0)
    W = len(LAC)
    for k in range(len(P)):
        for p in range(2):
            for i in range(net.ns):
                s0 = int(n0[k, p, i]) - AC_HALF
                lo, hi = max(0, s0), min(NT, s0 + W)
                if hi > lo:
                    Pr[i, :, lo:hi] += np.outer(a[k, p, i], w[k, p, i, lo - s0:hi - s0])
    return Pr


# ============================================================================================ TS-SFW
def init_atom(net, n, t, I, clip=True):
    src = net.garr[n]
    tau = float(net.taus[t])
    d = np.linalg.norm(net.STA - src, axis=1)
    T0 = d[None, :] * VINV[:, None] + tau + net.corrT
    dl = I * DT - T0
    if clip:
        dl = np.clip(dl, -SHIFT_MAX, SHIFT_MAX)
    return np.concatenate([src, [tau], dl.ravel()])


def tssfw_record(net, win, truth, f=None, lam=None, kmax=5, min_sep=5.0, refine=True, slide=True, stat="mtfw",
                 PD=None, clip_init=True, diag=None):
    """TS-SFW+SC (lam = f * s_1) o BOMP-PK-SLIDE+SC (lam = 0; f=None, lam=0).  refine/slide/stat/clip_init solo
    para validacion.  Devuelve el registro [(puntaje, acierto)]."""
    t_start = time.time()
    R = win.MFn.copy()
    P = np.zeros((0, net.nv))
    m = np.zeros((0, 6))
    scores, lam_used, it_d = [], None, []
    for k in range(kmax):
        excluded = excl_from(net, P, min_sep)
        if excluded.all():
            break
        s, n, t, I = oracle(net, R, excluded, stat, PD)
        if not np.isfinite(s):
            break
        if k == 0:
            lam_used = float(lam) if lam is not None else float(f) * s
        scores.append(s)
        p0 = init_atom(net, n, t, I, clip_init)
        d_it = {"s": s, "node": int(n), "tau": int(t), "I_shift": (I - net.BD.idx[n, t]).ravel().tolist()}
        a_old, T_old = pred_old(net, win, P, m)
        if refine:
            p1, inf = refine_atom(net, win, p0, a_old, T_old)
            d_it["refine"] = inf
        else:
            p1 = p0
        P = np.vstack([P, p1])
        m = np.vstack([m, np.zeros(6)])
        m, cinf = convex_step(net, win, P, m, lam_used)
        conv = [cinf]
        loc = []
        if slide:
            for _alt in range(MAX_ALT):
                P2, linf = slide_step(net, win, P, m, lam_used)
                loc.append(linf)
                P = P2
                m, cinf = convex_step(net, win, P, m, lam_used)
                conv.append(cinf)
                if linf["J0"] - linf["J1"] <= 1e-10 * max(1.0, abs(linf["J0"])):
                    break
        d_it["alternations"] = len(loc)
        d_it["local_nit"] = [x["nit"] for x in loc]
        d_it["local_status"] = [x["status"] for x in loc]
        d_it["max_rel_gap"] = max(c["rel_gap"] for c in conv)
        d_it["final_rel_gap"] = conv[-1]["rel_gap"]
        d_it["final_kkt"] = conv[-1]["kkt"]
        d_it["J"] = joint_obj(net, win, P, m, lam_used, grad=False)
        R = win.MFn - predict_sampled(net, win, P, m)
        it_d.append(d_it)
    rec, run = [], np.inf
    for k, s in enumerate(scores):
        run = min(run, s)
        rec.append((run, M.hit_of(P[k, :3], truth)))
    if diag is not None:
        diag["P"] = P
        diag["m"] = m
        diag["lam"] = lam_used
        diag["iters"] = it_d
        diag["seconds"] = time.time() - t_start
    return rec


def summarize_diag(dg, net):
    """resumen JSON de los diagnosticos de tssfw_record."""
    P = dg["P"]
    dl = P[:, 4:] if len(P) else np.zeros((0, 2 * net.ns))
    mn = np.sqrt(np.maximum(0.0, [float(np.linalg.norm(x)) ** 2 for x in dg["m"]])) if len(P) else []
    return {"lam": dg["lam"], "pos": P[:, :4].tolist(), "shifts": dl.tolist(), "m_norm": [float(x) for x in mn],
            "frac_shift_at_bound": float(np.mean(np.abs(dl) >= SHIFT_MAX - 1e-6)) if dl.size else 0.0,
            "iters": dg["iters"], "seconds": dg["seconds"]}


# ============================================================================================ MT-FW+SC con posiciones
def mtfw_sc_locs(STA, MFn, grid, taus, Gs, Qs, truth, tcorr, kmax=5, min_sep=5.0):
    """el bucle de certified_fw.certified_fw_record (mode 'enclosure') devolviendo ademas los nodos aceptados."""
    ns = len(STA)
    R = MFn.copy()
    lat = Lattice(grid, taus)
    garr = lat.garr
    taus = np.asarray(taus)
    excluded = np.zeros(len(garr), bool)
    Qst = stack_Q(Qs)
    acc_, rec, run = [], [], np.inf
    for _ in range(kmax):
        s, gi, vf, tpk = certified_oracle(STA, R, lat, taus, Qs, excluded, "enclosure", tcorr, None, Qst)
        if gi is None:
            break
        src = garr[gi]
        run = min(run, s)
        hit = -1
        if truth:
            dd = [float(np.linalg.norm(src[:2] - np.asarray(t_)[:2])) for t_ in truth]
            j = int(np.argmin(dd))
            if dd[j] <= TOL_KM:
                hit = j
        rec.append((run, hit))
        acc_.append(src)
        excluded = _exclusion_update(excluded, garr, src, min_sep)
        G = Gs[gi]
        mm, *_ = np.linalg.lstsq(G, vf, rcond=None)
        pred = (G @ mm).reshape(2, ns, 3)
        for b in range(2):
            for i in range(ns):
                t0 = int(tpk[b, i])
                lo, hi = max(0, t0 - X.AC_HALF), min(NT, t0 + X.AC_HALF + 1)
                R[i, :, lo:hi] -= np.outer(pred[b, i], X.AC[X.AC_HALF - (t0 - lo):X.AC_HALF + (hi - t0)])
    return rec, [a.tolist() for a in acc_]
