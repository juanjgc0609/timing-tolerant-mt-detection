"""
task2c_methods.py - metodos de la Tarea 2C (PREREGISTRO_1j_1k.txt): BOMP+SC, MTSTACK+SC y GL-FW+SC, todos en el
dominio del filtro adaptado (MFn = MF / sd, task1o_consistent.matched), sobre la rejilla default_grid() x DEFAULT_TAUS
con las correcciones de estacion tcorr (ns, 2) sumadas a los tiempos predichos (certified_fw._idx_node, el mismo
indice de muestra que node_vectors_c: int((d / V + tau + corr) / DT), recortado a [0, NT-1]).

MODELO.  Un bloque (atomo) es un par (nodo g, origen tau).  Su columna k (componente k del tensor de momento,
m = (Mxx, Myy, Mzz, Mxy, Mxz, Myz)) predice en el dominio del filtro adaptado, en la estacion i y la fase b, la
columna G_g[(b, i, :), k] de task1x_mtglrt.radiation_operator por la autocorrelacion AC (pico 1) centrada en la
muestra predicha EXACTA t_bi(g, tau) (sin busqueda de pico ni tope CAPV).  Es la imagen, por el filtro adaptado, de
un diccionario en el dominio de los datos: atomo de datos a_k = G[:, k] (x) ricker(t - t_bi) / sqrt(E_r),
E_r = sum KER^2, datos y = D sqrt(E_r) / sd.  Entonces (identidades exactas salvo redondeo; AC coincide con la
autocorrelacion de KER normalizada a 1e-17):
    A^T y        = muestreo de MFn en los tiempos predichos, contraido con G:  h = sum_bi G_bi^T MFn_i(:, t_bi)
    A^T A        = K, K_kl = sum_i sum_{p,q in P,S} G_pi,k . G_qi,l AC(t_pi - t_qi')   (Gram con AC, no identidad)
    A^T (y - Au) = muestreo del residuo R = MFn - sum_bloques sum_bi (G_bi u) (x) AC(t - t_bi)
Asi, correlacionar el residuo con un bloque es muestrear R en los tiempos predichos, y la energia de la proyeccion
del residuo sobre el bloque es e = h^T K_bb^{-1} h (en unidades de datos).  Con K_bb = G^T G (sin solape P/S) y
h = G^T v, e = ||Q^T v||^2: el estadistico de WAVE-MT evaluado en los tiempos exactos.  Los puntajes se dan como
e / 2, la convencion de WAVE-MT / MT-FW.

BOMP+SC (OMP por bloques, Vera Rodriguez et al. 2012): en cada paso, el bloque de mayor e sobre los nodos no
excluidos (exclusion de 5 km de fw_record: nodos a <= 5 km epicentrales de uno aceptado); puntaje = e / 2 del bloque
al entrar (la energia del residuo que explica el bloque por si solo), minimo acumulado como fw_record; luego
reajuste conjunto por minimos cuadrados de TODOS los bloques activos (u = K_SS^{-1} h_S(MFn), ortogonal) y residuo
recalculado desde MFn.  A lo sumo 5 eventos.

MTSTACK+SC (apilamiento corregido por mecanismo, Anikiev et al. 2014): en cada nodo y tau, m = lstsq(G, v) con v las
amplitudes de MFn en los tiempos predichos exactos (pinv de G); valor = ||G m||^2 / 2 (energia de la prediccion);
campo por nodo = maximo sobre tau; catalogo por supresion voraz de maximos a 20 km (froc_exact.greedy_record), a lo
sumo 5 eventos, sin deflacion.

GL-FW+SC (lasso de grupos sobre el diccionario de BOMP, el problema convexo de la tesis en la rejilla):
    min_u  F(u) = 1/2 ||y - A u||^2 + lambda sum_b ||u_b||_K,   ||u_b||_K = sqrt(u_b^T K_bb u_b) = ||A_b u_b||
certificado dual del bloque b: c_b = sqrt(h_b^T K_bb^{-1} h_b) con h = A^T (y - Au) (KKT: c_b <= lambda, = lambda en
el soporte).  Frank-Wolfe totalmente correctivo (conjunto activo): en cada iteracion se resuelve exactamente el
problema restringido a los bloques seleccionados (descenso por bloques en coordenadas blanqueadas w_b = T_b u_b,
T_b^T T_b = K_bb, cuya actualizacion por bloque es el umbral de grupo exacto), se calculan los certificados de todos
los bloques y, si el mayor supera lambda (1 + TOL_FW), se anade ese bloque (paso de Frank-Wolfe) y se repite.
Trayectoria: 30 valores geometricos de lambda_max = max_b c_b(y) a lambda_max / 30, con arranque en caliente.
Brecha de dualidad en cada lambda con el punto dual theta = s r, s = min(1, lambda / max_b c_b(r)):
    gap = 1/2 (1 - s)^2 ||r||^2 + lambda Omega(u) - s u^T h(r),   ||r||^2 = ||y||^2 - 2 u^T A^T y + u^T K u.
Puntaje de un bloque = el mayor lambda de la rejilla en que su norma de grupo es > 0; registro: bloques por lambda
de entrada decreciente (empates: mayor norma de grupo en ese lambda), saltando los que caen a <= 5 km de uno ya
aceptado, a lo sumo 5.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import task1x_mtglrt as X  # noqa: E402
from certified_fw import _idx_node  # noqa: E402
from froc_exact import greedy_record  # noqa: E402
from task0_common import NT  # noqa: E402
from task1o_consistent import _KER, TOL_KM  # noqa: E402

AC, AC_HALF = X.AC, X.AC_HALF
E_R = float(np.sum(_KER ** 2))
N_LAM = 30
LAM_RATIO = 30.0
TOL_FW = 1e-6            # tolerancia relativa del certificado (KKT fuera del soporte)
BCD_TOL = 1e-12          # cambio relativo maximo por barrido en el problema restringido
BCD_MAX = 200000
RESTR_TOL = 1e-9       # violacion KKT relativa del problema restringido (<< TOL_FW)
EIG_FLOOR = 1e-10


def ac_at(lag):
    """AC(lag) en muestras (0 fuera del soporte de AC)."""
    lag = np.asarray(lag)
    ok = np.abs(lag) <= AC_HALF
    return np.where(ok, AC[np.clip(lag + AC_HALF, 0, len(AC) - 1)], 0.0)


def hit_of(src, truth):
    hit = -1
    if truth:
        dd = [float(np.linalg.norm(src[:2] - np.asarray(t_)[:2])) for t_ in truth]
        j = int(np.argmin(dd))
        if dd[j] <= TOL_KM:
            hit = j
    return hit


class BlockDict:
    """Diccionario por bloques de una red (depende de STA, la rejilla, los origenes y tcorr, no de la ventana)."""

    def __init__(self, STA, grid, taus, tcorr=None, Gs=None):
        self.STA = np.asarray(STA, float)
        self.garr = np.asarray(grid, float)
        self.taus = np.asarray(taus, float)
        self.ns = ns = len(STA)
        N, T = len(self.garr), len(self.taus)
        self.N, self.T = N, T
        if Gs is None:
            Gs = [X.radiation_operator(g, self.STA) for g in self.garr]
        self.Gf = np.ascontiguousarray(np.stack(Gs))                       # (N, 6ns, 6), filas (b, i, c)
        self.Gb = self.Gf.reshape(N, 2, ns, 3, 6)
        self.idx = np.stack([_idx_node(self.STA, g, self.taus, tcorr) for g in self.garr])   # (N, T, 2, ns)
        self.st = np.arange(ns)[None, None, None, :]
        # Gram diagonal de cada bloque: G^T G + terminos cruzados P/S de la misma estacion con AC(t_P - t_S)
        GtG = np.einsum("nxk,nxl->nkl", self.Gf, self.Gf)
        C = np.einsum("nick,nicl->nikl", self.Gb[:, 0], self.Gb[:, 1])      # (N, ns, 6, 6): G_P^T G_S
        lag = self.idx[:, :, 1, :] - self.idx[:, :, 0, :]                   # t_S - t_P (N, T, ns)
        wm, wp = ac_at(-lag), ac_at(lag)
        K = (GtG[:, None] + np.einsum("nti,nikl->ntkl", wm, C) + np.einsum("nti,nilk->ntkl", wp, C))
        self.K = 0.5 * (K + np.swapaxes(K, -1, -2))                         # (N, T, 6, 6)
        ev, V = np.linalg.eigh(self.K)
        keep = ev > EIG_FLOOR * ev[..., -1:]
        self.rank = keep.sum(-1)
        inv = np.where(keep, 1.0 / np.where(keep, ev, 1.0), 0.0)
        self.Kinv = np.einsum("ntkj,ntj,ntlj->ntkl", V, inv, V)
        sq = np.where(keep, np.sqrt(np.where(keep, ev, 1.0)), 0.0)
        isq = np.where(keep, 1.0 / np.where(keep, np.sqrt(np.where(keep, ev, 1.0)), 1.0), 0.0)
        self.Tw = np.einsum("ntj,ntkj->ntjk", sq, V)                        # w = Tw u  (Tw^T Tw = K)
        self.Twinv = np.einsum("ntkj,ntj->ntkj", V, isq)                    # u = Twinv w (en el rango)
        self.min_ps_lag = int(np.abs(lag).min())
        self.max_ps_ac = float(np.abs(ac_at(lag)).max())

    # ------------------------------------------------------------------ correlacion y sintesis
    def sample(self, R):
        """amplitudes de R en los tiempos predichos de todos los bloques: (N, T, 6ns), orden (b, i, c)."""
        V = R[self.st, :, self.idx]                                          # (N, T, 2, ns, 3)
        return V.reshape(self.N, self.T, -1)

    def correlate(self, R):
        """h = A^T R para todos los bloques (N, T, 6)."""
        return np.matmul(self.sample(R), self.Gf)

    def energy(self, h):
        """e = h^T K^{-1} h (N, T)."""
        return np.einsum("ntk,ntk->nt", h, np.matmul(self.Kinv, h[..., None])[..., 0])

    def corr_blocks(self, R, blocks):
        """h de una lista de bloques [(n, t)]: (nb, 6)."""
        n = np.array([b[0] for b in blocks], int)
        t = np.array([b[1] for b in blocks], int)
        V = R[np.arange(self.ns)[None, None, :], :, self.idx[n, t]]          # (nb, 2, ns, 3)
        return np.einsum("bx,bxk->bk", V.reshape(len(blocks), -1), self.Gf[n])

    def synth_sub(self, R, n, t, u):
        """R -= contribucion del bloque (n, t) con coeficientes u (la deflacion de fw_record en los tiempos exactos)."""
        pred = (self.Gf[n] @ u).reshape(2, self.ns, 3)
        idx = self.idx[n, t]
        for b in range(2):
            for i in range(self.ns):
                t0 = int(idx[b, i])
                lo, hi = max(0, t0 - AC_HALF), min(NT, t0 + AC_HALF + 1)
                R[i, :, lo:hi] -= np.outer(pred[b, i], AC[AC_HALF - (t0 - lo):AC_HALF + (hi - t0)])

    def residual(self, MFn, blocks, U):
        R = MFn.copy()
        for (n, t), u in zip(blocks, U):
            self.synth_sub(R, n, t, u)
        return R

    def cross_gram(self, A, B):
        """K entre dos listas de bloques: (na, 6, nb, 6)."""
        nA = np.array([a[0] for a in A], int)
        tA = np.array([a[1] for a in A], int)
        nB = np.array([b[0] for b in B], int)
        tB = np.array([b[1] for b in B], int)
        IA, IB = self.idx[nA, tA], self.idx[nB, tB]                           # (na, 2, ns)
        lag = IA[:, None, :, None, :] - IB[None, :, None, :, :]              # (na, nb, p, q, ns)
        W = ac_at(lag)
        return np.einsum("apick,bqicl,abpqi->akbl", self.Gb[nA], self.Gb[nB], W, optimize=True)

    def gram(self, S):
        n = len(S)
        K = self.cross_gram(S, S).reshape(6 * n, 6 * n)
        return 0.5 * (K + K.T)


def excl_mask(mask, garr, src, min_sep):
    """nodos a <= min_sep km epicentrales de src (la prueba de fw_record)."""
    return mask | (np.linalg.norm(garr[:, :2] - src[:2], axis=1) <= min_sep)


# ================================================================================================ BOMP
def bomp_record(BD, MFn, truth, kmax=5, min_sep=5.0, diag=None):
    garr = BD.garr
    excluded = np.zeros(BD.N, bool)
    S, rec, run = [], [], np.inf
    R = MFn.copy()
    U = None
    for _ in range(kmax):
        e = BD.energy(BD.correlate(R))
        e[excluded] = -np.inf
        j = int(np.argmax(e))
        n, t = divmod(j, BD.T)
        s = float(e[n, t]) / 2.0
        if not np.isfinite(s):
            break
        run = min(run, s)
        src = garr[n]
        rec.append((run, hit_of(src, truth)))
        S.append((n, t))
        excluded = excl_mask(excluded, garr, src, min_sep)
        K = BD.gram(S)
        b = BD.corr_blocks(MFn, S).ravel()
        U, *_ = np.linalg.lstsq(K, b, rcond=None)
        U = U.reshape(len(S), 6)
        R = BD.residual(MFn, S, U)
        if diag is not None:
            hS = BD.corr_blocks(R, S)
            diag.setdefault("orth", []).append(float(np.abs(hS).max() / max(1e-300, np.abs(b).max())))
    if diag is not None:
        diag["blocks"], diag["U"] = S, U
    return rec


# ============================================================================================= MTSTACK
class Stacker:
    def __init__(self, BD):
        self.BD = BD
        self.Pinv = np.stack([np.linalg.pinv(G) for G in BD.Gf])             # (N, 6, 6ns)

    def field(self, MFn):
        v = self.BD.sample(MFn)                                               # (N, T, 6ns)
        m = np.matmul(v, np.swapaxes(self.Pinv, 1, 2))                        # (N, T, 6)  m = pinv(G) v
        pred = np.matmul(m, np.swapaxes(self.BD.Gf, 1, 2))                    # (N, T, 6ns)
        return np.sum(pred ** 2, axis=-1) / 2.0                               # (N, T)

    def record(self, MFn, truth, min_sep_km=20.0, max_events=5):
        F = self.field(MFn)
        return greedy_record(list(self.BD.garr), F.max(axis=1), truth, min_sep_km=min_sep_km, max_events=max_events)


# ============================================================================================== GL-FW
try:
    from numba import njit
except Exception:  # noqa: BLE001
    njit = None


def _bcd_py(Kp, bp, w, Kw, lam, ng, max_sweeps, tol):
    for sw in range(max_sweeps):
        md, mw = 0.0, 0.0
        for j in range(ng):
            a = 6 * j
            z = bp[a:a + 6] - Kw[a:a + 6] + w[a:a + 6]
            nz = np.sqrt(z @ z)
            wn = (1.0 - lam / nz) * z if nz > lam else np.zeros(6)
            d = wn - w[a:a + 6]
            dd = np.sqrt(d @ d)
            if dd > 0:
                Kw += Kp[:, a:a + 6] @ d
                w[a:a + 6] = wn
            md = max(md, dd)
            mw = max(mw, np.sqrt(wn @ wn))
        if md <= tol * max(mw, 1e-300):
            return sw + 1
    return max_sweeps


if njit is not None:
    @njit(cache=True)
    def _bcd_nb(Kp, bp, w, Kw, lam, ng, max_sweeps, tol):
        n6 = 6 * ng
        z = np.zeros(6)
        wn = np.zeros(6)
        d = np.zeros(6)
        for sw in range(max_sweeps):
            md = 0.0
            mw = 0.0
            for j in range(ng):
                a = 6 * j
                nz = 0.0
                for k in range(6):
                    z[k] = bp[a + k] - Kw[a + k] + w[a + k]
                    nz += z[k] * z[k]
                nz = np.sqrt(nz)
                f = 1.0 - lam / nz if nz > lam else 0.0
                dd = 0.0
                ww = 0.0
                for k in range(6):
                    wn[k] = f * z[k]
                    d[k] = wn[k] - w[a + k]
                    dd += d[k] * d[k]
                    ww += wn[k] * wn[k]
                dd = np.sqrt(dd)
                if dd > 0.0:
                    for r in range(n6):
                        acc = 0.0
                        for k in range(6):
                            acc += Kp[r, a + k] * d[k]
                        Kw[r] += acc
                    for k in range(6):
                        w[a + k] = wn[k]
                if dd > md:
                    md = dd
                ww = np.sqrt(ww)
                if ww > mw:
                    mw = ww
            if md <= tol * max(mw, 1e-300):
                return sw + 1
        return max_sweeps


def bcd(Kp, bp, w, lam, max_sweeps=BCD_MAX, tol=BCD_TOL):
    Kw = Kp @ w
    f = _bcd_nb if njit is not None else _bcd_py
    return f(Kp, bp, w, Kw, float(lam), len(bp) // 6, int(max_sweeps), float(tol))


def _phi(Kp, bp, w, lam):
    W = w.reshape(-1, 6)
    return 0.5 * w @ Kp @ w - bp @ w + lam * float(np.sum(np.sqrt(np.sum(W * W, axis=1))))


def kkt_restricted(Kp, bp, w, lam):
    """violacion KKT relativa (/lambda) del problema restringido: grupos no nulos |(Kw - b)_g + lam w_g/|w_g||,
    grupos nulos max(0, |b_g - (Kw)_g| - lam)."""
    r = (bp - Kp @ w).reshape(-1, 6)
    W = w.reshape(-1, 6)
    nw = np.sqrt(np.sum(W * W, axis=1))
    nr = np.sqrt(np.sum(r * r, axis=1))
    act = nw > 0
    v = np.zeros(len(nw))
    if act.any():
        v[act] = np.sqrt(np.sum((-r[act] + lam * W[act] / nw[act, None]) ** 2, axis=1))
    v[~act] = np.maximum(0.0, nr[~act] - lam)
    return float(v.max() / lam) if len(v) else 0.0


def solve_restricted(Kp, bp, w, lam, tol=RESTR_TOL, max_outer=50):
    """lasso de grupos restringido (bloques diagonales de Kp = I): barridos de descenso por bloques (umbral de grupo
    exacto) para fijar el soporte y pulido de Newton amortiguado sobre los grupos no nulos (la funcion objetivo es
    suave alli; busqueda lineal de Armijo con holgura de redondeo 1e-14 |F|), hasta que la violacion KKT relativa
    (kkt_restricted) sea <= tol en todos los grupos de S."""
    sweeps = 0
    newton = 0
    for _ in range(max_outer):
        sweeps += bcd(Kp, bp, w, lam, max_sweeps=30, tol=1e-13)
        if kkt_restricted(Kp, bp, w, lam) <= tol:
            return sweeps, newton
        W = w.reshape(-1, 6)
        act = np.where(np.sqrt(np.sum(W * W, axis=1)) > 0)[0]
        if len(act):
            ix = (6 * act[:, None] + np.arange(6)[None, :]).ravel()
            Kaa = Kp[np.ix_(ix, ix)]
            for _it in range(30):
                wa = w[ix]
                Wa = wa.reshape(-1, 6)
                nw = np.sqrt(np.sum(Wa * Wa, axis=1))
                if np.any(nw <= 1e-14 * max(1.0, nw.max())):
                    break
                Nn = Wa / nw[:, None]
                g = Kp[ix] @ w - bp[ix] + lam * Nn.ravel()
                if np.max(np.sqrt(np.sum(g.reshape(-1, 6) ** 2, axis=1))) <= 0.1 * tol * lam:
                    break
                H = Kaa.copy()
                for q in range(len(act)):
                    sl = slice(6 * q, 6 * q + 6)
                    H[sl, sl] += lam * (np.eye(6) - np.outer(Nn[q], Nn[q])) / nw[q]
                try:
                    dw = -np.linalg.solve(H, g)
                except np.linalg.LinAlgError:
                    dw = -np.linalg.lstsq(H, g, rcond=None)[0]
                f0 = _phi(Kp, bp, w, lam)
                t = 1.0
                gd = float(g @ dw)
                ok = False
                for _ls in range(30):
                    w2 = w.copy()
                    w2[ix] = wa + t * dw
                    if _phi(Kp, bp, w2, lam) <= f0 + 1e-4 * t * gd + 1e-14 * abs(f0):
                        ok = True
                        break
                    t *= 0.5
                if not ok:
                    break
                w[ix] = w2[ix]
                newton += 1
        if kkt_restricted(Kp, bp, w, lam) <= tol:
            return sweeps, newton
    return sweeps, newton


class GLPath:
    """trayectoria del lasso de grupos por Frank-Wolfe totalmente correctivo (conjunto activo)."""

    def __init__(self, BD, MFn, yy):
        self.BD, self.MFn, self.yy = BD, MFn, float(yy)          # yy = ||y||^2 en unidades de datos
        self.S = []                                              # bloques activos (n, t)
        self.w = np.zeros(0)                                     # coordenadas blanqueadas (6|S|)
        self.Kp = np.zeros((0, 0))                               # Gram blanqueado de S
        self.bp = np.zeros(0)
        self.KS = np.zeros((0, 0))                               # Gram de S (sin blanquear)
        self.bS = np.zeros(0)

    def _T(self, blk):
        n, t = blk
        return self.BD.Tw[n, t], self.BD.Twinv[n, t]

    def add(self, blk):
        BD = self.BD
        k_new = BD.cross_gram(self.S + [blk], [blk]).reshape(6 * (len(self.S) + 1), 6)
        m = 6 * len(self.S)
        KS = np.zeros((m + 6, m + 6))
        KS[:m, :m] = self.KS
        KS[:, m:] = k_new
        KS[m:, :] = k_new.T
        KS[m:, m:] = 0.5 * (KS[m:, m:] + KS[m:, m:].T)
        self.KS = KS
        self.bS = np.concatenate([self.bS, BD.corr_blocks(self.MFn, [blk])[0]])
        self.S.append(blk)
        self.w = np.concatenate([self.w, np.zeros(6)])
        self._whiten()

    def drop_zeros(self):
        keep = [j for j in range(len(self.S)) if np.any(self.w[6 * j:6 * j + 6] != 0.0)]
        if len(keep) == len(self.S):
            return
        ix = np.concatenate([np.arange(6 * j, 6 * j + 6) for j in keep]) if keep else np.zeros(0, int)
        self.S = [self.S[j] for j in keep]
        self.KS = self.KS[np.ix_(ix, ix)]
        self.bS = self.bS[ix]
        self.w = self.w[ix]
        self._whiten()

    def _whiten(self):
        """Kp = Ti^T KS Ti y bp = Ti^T bS con Ti = diag(Twinv_b) (bloque a bloque)."""
        n = len(self.S)
        if n == 0:
            self.Kp, self.bp = np.zeros((0, 0)), np.zeros(0)
            self.Ti = np.zeros((0, 6, 6))
            return
        self.Ti = np.stack([self._T(blk)[1] for blk in self.S])                    # (n, 6, 6)
        K4 = self.KS.reshape(n, 6, n, 6)
        Kp = np.einsum("iak,iajb,jbl->ikjl", self.Ti, K4, self.Ti, optimize=True).reshape(6 * n, 6 * n)
        self.Kp = np.ascontiguousarray(0.5 * (Kp + Kp.T))
        self.bp = np.einsum("iak,ia->ik", self.Ti, self.bS.reshape(n, 6)).ravel()

    def u(self):
        n = len(self.S)
        return np.einsum("ikl,il->ik", self.Ti, self.w.reshape(n, 6)) if n else np.zeros((0, 6))

    def omega(self):
        return float(sum(np.linalg.norm(self.w[6 * j:6 * j + 6]) for j in range(len(self.S))))

    def certificates(self):
        R = self.BD.residual(self.MFn, self.S, self.u())
        h = self.BD.correlate(R)
        return np.sqrt(np.maximum(self.BD.energy(h), 0.0)), h

    def gap(self, lam, c, h):
        U = self.u().ravel()
        cmax = float(c.max())
        s = min(1.0, lam / cmax) if cmax > 0 else 1.0
        rr = self.yy - 2.0 * U @ self.bS + U @ self.KS @ U
        hS = np.concatenate([h[n, t] for n, t in self.S]) if self.S else np.zeros(0)
        g = 0.5 * (1.0 - s) ** 2 * rr + lam * self.omega() - s * (U @ hS)
        dec = -(0.5 * U @ self.KS @ U - U @ self.bS + lam * self.omega())     # F(0) - F(u)
        return float(g), float(dec), float(rr)

    def solve(self, lam, max_fw=100000, max_stall=5):
        """devuelve certificados, h, pasos FW, barridos y el numero de re-solves del problema restringido por un
        bloque activo con certificado > lambda (1 + TOL_FW) (restringido no convergido; 0 en toda corrida sana)."""
        n_fw, sweeps, stall = 0, 0, 0
        while True:
            if self.S:
                sw, _nn = solve_restricted(self.Kp, self.bp, self.w, lam)
                sweeps += sw
            c, h = self.certificates()
            j = int(np.argmax(c))
            n, t = divmod(j, self.BD.T)
            if c[n, t] <= lam * (1.0 + TOL_FW) or n_fw >= max_fw:
                return c, h, n_fw, sweeps, stall
            if (n, t) in self.S:
                stall += 1
                if stall > max_stall:
                    return c, h, n_fw, sweeps, stall
                continue
            self.add((n, t))
            n_fw += 1


def path_record(BD, entry, truth, kmax=5, min_sep=5.0, with_nv=False):
    order = sorted(entry.items(), key=lambda kv: (-kv[1][0], -kv[1][1], kv[0]))
    acc, rec = [], []
    for (n, t), (lam, nv) in order:
        src = BD.garr[n]
        if any(np.linalg.norm(src[:2] - a[:2]) <= min_sep for a in acc):
            continue
        acc.append(src)
        rec.append((lam, hit_of(src, truth), nv, int(n), int(t)) if with_nv else (lam, hit_of(src, truth)))
        if len(acc) >= kmax:
            break
    return rec, order


def glfw_record(BD, MFn, D, sd, truth, kmax=5, min_sep=5.0, n_lam=N_LAM, ratio=LAM_RATIO, diag=None,
                full_path=False, max_support=None, extra=0):
    """registro de GL-FW+SC y diagnosticos por lambda (brecha, KKT, soporte, iteraciones FW).
    full_path=False: la trayectoria se detiene en el primer lambda de la rejilla en que el registro ya tiene kmax
    eventos (los atomos que entran despues tienen lambda de entrada estrictamente menor y quedan detras de los kmax
    aceptados: el registro es identico al de la trayectoria completa); extra > 0 resuelve ademas esos lambda
    siguientes (solo para comprobar la identidad)."""
    yy = E_R * float(np.sum(D.astype(float) ** 2)) / sd ** 2
    P = GLPath(BD, MFn, yy)
    c0 = np.sqrt(np.maximum(BD.energy(BD.correlate(MFn)), 0.0))
    lam_max = float(c0.max())
    lams = lam_max * ratio ** (-np.arange(n_lam) / (n_lam - 1))
    entry = {}                                    # bloque -> (lambda de entrada, norma de grupo en ese lambda)
    per = []
    stopped = None
    for k, lam in enumerate(lams):
        if max_support is not None and len(P.S) > max_support:
            stopped = ("max_support", k)
            break
        c, h, n_fw, sweeps, stall = P.solve(lam)
        nrm = [float(np.linalg.norm(P.w[6 * j:6 * j + 6])) for j in range(len(P.S))]
        for blk, nv in zip(P.S, nrm):
            if nv > 0.0 and blk not in entry:
                entry[blk] = (float(lam), nv)
        g, dec, rr = P.gap(lam, c, h)
        insup = np.zeros(c.shape, bool)
        onsup = [c[n, t] / lam for (n, t), nv in zip(P.S, nrm) if nv > 0]
        for (n, t), nv in zip(P.S, nrm):
            if nv > 0:
                insup[n, t] = True
        off = c[~insup]
        per.append({"lam": float(lam), "gap": g, "rel_gap": (g / dec if dec > 0 else (0.0 if g == 0 else np.inf)),
                    "gap_over_lam_omega": (g / (lam * P.omega()) if P.omega() > 0 else 0.0),
                    "n_support": int(sum(nv > 0 for nv in nrm)), "n_fw": n_fw, "sweeps": int(sweeps), "stall": int(stall),
                    "kkt_off_max": float(off.max() / lam) if off.size else 0.0,
                    "kkt_on_maxdev": float(max(abs(x - 1.0) for x in onsup)) if onsup else 0.0})
        P.drop_zeros()
        if not full_path and stopped is None and len(path_record(BD, entry, truth, kmax, min_sep)[0]) >= kmax:
            stopped = ("record_complete", k)
        if stopped is not None and k >= stopped[1] + extra:
            break
    rec, order = path_record(BD, entry, truth, kmax, min_sep)
    if diag is not None:
        diag["per_lambda"] = per
        diag["lam_max"] = lam_max
        diag["n_entered"] = len(entry)
        diag["stopped"] = stopped
        diag["first_blocks"] = [[int(n), int(t), lam, nv] for (n, t), (lam, nv) in order[:10]]
        diag["entry"] = entry
        diag["rec_nv"] = path_record(BD, entry, truth, kmax, min_sep, with_nv=True)[0]
    return rec
