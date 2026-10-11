"""
task2d_methods.py - metodos de la Tarea 2D (PREREGISTRO_1j_1k.txt, seccion "Tarea 2D"): igualdad de tolerancia temporal
para los comparadores de 2C y ablaciones de MT-FW.  Reutiliza task2c_methods (BlockDict, GLPath, bcd, solve_restricted),
certified_fw (_idx_node, node_vectors_c) y task1x_mtglrt (OFFS, CAPV, AC).

MUESTREO EN EL PICO (PK).  Para un residuo R (ns, 3, NT) y una muestra predicha idx (int((d / V + tau + corr) / DT),
recortada, la de node_vectors_c), la muestra de pico es el argmax de la energia de tres componentes E_i(k) = ||R_i[:, k]||^2
en la ventana clip(idx + o), |o| <= SH = 4 muestras (0.25 s), con la regla de empate de node_vectors (primer maximo en el
orden de los desplazamientos).  peak_table(R) tabula ese argmax para cada estacion y cada centro (la formula de
certified_fw.LeafTables); en toda la rejilla es identico al tpk de node_vectors_c (verificado en task2d_validate.py).

ATOMOS CON TIEMPOS FIJOS.  Un atomo es (n, t, I): nodo n, origen t y las muestras I (2, ns) de P y S por estacion en que se
lee.  Su columna k predice en el dominio del filtro adaptado G_bi[:, k] AC(s - I_bi).  La Gram entre atomos es la de
BlockDict.cross_gram con I en lugar de los tiempos predichos; la Gram diagonal K(I) = G^T G + terminos cruzados P/S con
AC(I_S - I_P) (la formula de BlockDict con el retardo de pico).  Con I = los tiempos predichos todo se reduce a 2C.

BOMP-PK+SC: en cada paso, para cada bloque (n, t) no excluido, h = G^T v con v el residuo en las muestras de pico (sin
tope), e = h^T K(I_pk)^+ h (seudoinversa con el piso de valores propios de BlockDict); entra el de mayor e con sus
tiempos FIJADOS en esos picos; puntaje e / 2 (minimo acumulado); reajuste ortogonal de todos los atomos activos con sus
tiempos fijos y residuo recalculado desde MFn.  A lo sumo 5, exclusion de 5 km.
MTSTACK-PK+SC: v = MFn en las muestras de pico de MFn (sin tope), m = pinv(G) v, valor ||G m||^2 / 2, maximo sobre tau,
supresion voraz de maximos a 20 km, a lo sumo 5, sin deflacion.
GL-FW-PK+SC: el lasso de grupos de 2C sobre un diccionario que crece: el oraculo de Frank-Wolfe evalua el certificado
c = sqrt(h^T K(I_pk)^+ h) de cada bloque en las muestras de pico del residuo actual; el atomo que entra queda con esos
tiempos fijos; el problema restringido (exacto, task2c_methods.solve_restricted) es sobre los atomos fijos.  Trayectoria,
puntaje (lambda de entrada), registro y parada temprana como 2C.  La brecha de dualidad se calcula respecto del
diccionario realizado (atomos activos con sus tiempos fijos y los atomos de pico del residuo final).
MT-FW-EX+SC: fw_record (dense_fw_record_c) con la busqueda de pico desactivada: cada estacion y fase se lee en la muestra
predicha exacta (con tope CAPV) y la deflacion se centra alli.
MT-FW-NOCAP+SC: fw_record con busqueda de pico y sin tope CAPV.
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
import task2c_methods as M  # noqa: E402
from certified_fw import _idx_node, node_vectors_c  # noqa: E402
from froc_exact import greedy_record  # noqa: E402
from task0_common import NT  # noqa: E402
from task1o_consistent import TOL_KM  # noqa: E402

AC, AC_HALF = X.AC, X.AC_HALF
SH = X.SH
OFFS = X.OFFS
II_TAB = np.clip(np.arange(NT)[:, None] + OFFS[None, :], 0, NT - 1)          # (NT, W)
MAX_FW_PER_LAM = 400


# ============================================================================== MT-FW variants (dense loop)
def node_vectors_v(MFn, STA, src, taus, tcorr, peak=True, cap=True):
    """node_vectors_c con la busqueda de pico y el tope opcionales.  peak=True, cap=True es node_vectors_c."""
    if peak and cap:
        v, tpk, _s, _e = node_vectors_c(MFn, STA, src, taus, tcorr)
        return v, tpk
    ns = len(STA)
    st = np.arange(ns)
    idx = _idx_node(STA, src, taus, tcorr)                                     # (nt, 2, ns)
    if peak:
        ii = np.clip(idx[..., None] + OFFS, 0, NT - 1)
        seg = MFn[st[None, None, :, None], :, ii]
        ener = np.sum(seg ** 2, axis=-1)
        k = ener.argmax(axis=-1)
        v = np.take_along_axis(seg, k[..., None, None], axis=3)[:, :, :, 0, :]
        tpk = np.take_along_axis(ii, k[..., None], axis=3)[..., 0]
    else:
        v = MFn[st[None, None, :], :, idx]                                     # (nt, 2, ns, 3)
        tpk = idx
    if cap:
        nv = np.linalg.norm(v, axis=-1, keepdims=True)
        v = v * np.minimum(1.0, X.CAPV / np.maximum(nv, 1e-12))
    return v, tpk


def fw_record_v(STA, MFn, grid, taus, Gs, Qs, truth, tcorr, peak=True, cap=True, kmax=5, min_sep=5.0):
    """el bucle de certified_fw.dense_fw_record_c (= task1x_mtglrt.fw_record con tcorr) con node_vectors_v."""
    ns = len(STA)
    R = MFn.copy()
    garr = np.asarray(grid)
    acc, rec, run = [], [], np.inf
    for _ in range(kmax):
        best = (-np.inf, None, None, None)
        for gi, src in enumerate(garr):
            if any(np.linalg.norm(src[:2] - a[:2]) <= min_sep for a in acc):
                continue
            v, tpk = node_vectors_v(R, STA, src, taus, tcorr, peak, cap)
            vf = v.reshape(len(taus), -1)
            sc = np.sum((vf @ Qs[gi]) ** 2, axis=1) / 2.0
            j = int(sc.argmax())
            if sc[j] > best[0]:
                best = (float(sc[j]), gi, vf[j], tpk[j])
        s, gi, vf, tpk = best
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
        acc.append(src)
        G = Gs[gi]
        m, *_ = np.linalg.lstsq(G, vf, rcond=None)
        pred = (G @ m).reshape(2, ns, 3)
        for b in range(2):
            for i in range(ns):
                t0 = int(tpk[b, i])
                lo, hi = max(0, t0 - AC_HALF), min(NT, t0 + AC_HALF + 1)
                R[i, :, lo:hi] -= np.outer(pred[b, i], AC[AC_HALF - (t0 - lo):AC_HALF + (hi - t0)])
    return rec


# ============================================================================== peak sampling
def peak_table(R):
    """(ns, NT): para cada estacion y centro t0, la muestra de la ventana clip(t0 + o), |o| <= SH, con la mayor energia
    de tres componentes (primer maximo), como node_vectors (formula de certified_fw.LeafTables)."""
    Et = np.sum(np.ascontiguousarray(R.transpose(0, 2, 1)) ** 2, axis=-1)       # (ns, NT)
    W = Et[:, II_TAB]                                                            # (ns, NT, W)
    a = W.argmax(axis=-1)
    return II_TAB[np.arange(NT)[None, :], a]                                     # (ns, NT)


class PKDict:
    """extension de task2c_methods.BlockDict a atomos con tiempos de lectura arbitrarios I (2, ns)."""

    def __init__(self, BD):
        self.BD = BD
        self.ns, self.N, self.T = BD.ns, BD.N, BD.T
        self.stv = np.arange(BD.ns)
        self.st4 = BD.st                                                         # (1, 1, 1, ns)
        self.GtG = np.einsum("nxk,nxl->nkl", BD.Gf, BD.Gf)
        self.C = np.einsum("nick,nicl->nikl", BD.Gb[:, 0], BD.Gb[:, 1])          # (N, ns, 6, 6): G_P^T G_S
        self.nn = np.arange(BD.N)[:, None]

    # ------------------------------------------------------------------ indices de lectura
    def read_idx(self, R, peak):
        """(N, T, 2, ns): muestras de lectura de todos los bloques (picos del residuo R, o las predichas si peak=False)."""
        if not peak:
            return self.BD.idx
        Kt = peak_table(R)
        return Kt[self.st4, self.BD.idx]

    def sample_idx(self, R, I):
        """amplitudes de R en I (..., 2, ns): (..., 2, ns, 3)."""
        st = self.stv.reshape((1,) * (I.ndim - 1) + (self.ns,))
        return R[st, :, I]

    # ------------------------------------------------------------------ Gram diagonal y energia
    def K_all(self, I):
        """K(I) de todos los bloques (N, T, 6, 6), la formula de BlockDict con el retardo I_S - I_P."""
        lag = I[:, :, 1, :] - I[:, :, 0, :]
        wm, wp = M.ac_at(-lag), M.ac_at(lag)
        K = (self.GtG[:, None] + np.einsum("nti,nikl->ntkl", wm, self.C) + np.einsum("nti,nilk->ntkl", wp, self.C))
        return 0.5 * (K + np.swapaxes(K, -1, -2))

    @staticmethod
    def pinv_sym(K):
        ev, V = np.linalg.eigh(K)
        keep = ev > M.EIG_FLOOR * ev[..., -1:]
        inv = np.where(keep, 1.0 / np.where(keep, ev, 1.0), 0.0)
        return np.einsum("...kj,...j,...lj->...kl", V, inv, V)

    def energy(self, R, peak):
        """e = h^T K(I)^+ h en todos los bloques (N, T) y los indices de lectura I (N, T, 2, ns).  Con peak=False es
        BlockDict.energy(BlockDict.correlate(R)).  Con peak=True, la seudoinversa se recalcula solo en los bloques cuyo
        K(I) difiere del de los tiempos predichos (distinto retardo P/S con AC no nulo); en los demas es BlockDict.Kinv."""
        BD = self.BD
        I = self.read_idx(R, peak)
        h = np.matmul(self.sample_idx(R, I).reshape(self.N, self.T, -1), BD.Gf)
        if not peak:
            return BD.energy(h), I, h
        lag0 = BD.idx[:, :, 1, :] - BD.idx[:, :, 0, :]
        lag = I[:, :, 1, :] - I[:, :, 0, :]
        diff = np.any((M.ac_at(lag) != M.ac_at(lag0)) | (M.ac_at(-lag) != M.ac_at(-lag0)), axis=-1)   # (N, T)
        Kinv = BD.Kinv
        if diff.any():
            Kinv = Kinv.copy()
            n_, t_ = np.nonzero(diff)
            Isub = I[n_, t_]                                                     # (m, 2, ns)
            lg = Isub[:, 1, :] - Isub[:, 0, :]
            wm, wp = M.ac_at(-lg), M.ac_at(lg)
            K = (self.GtG[n_] + np.einsum("mi,mikl->mkl", wm, self.C[n_]) + np.einsum("mi,milk->mkl", wp, self.C[n_]))
            K = 0.5 * (K + np.swapaxes(K, -1, -2))
            Kinv[n_, t_] = self.pinv_sym(K)
        e = np.einsum("ntk,ntk->nt", h, np.matmul(Kinv, h[..., None])[..., 0])
        return e, I, h

    # ------------------------------------------------------------------ atomos (n, t, I)
    def K_atom(self, n, I):
        lg = I[1] - I[0]
        wm, wp = M.ac_at(-lg), M.ac_at(lg)
        K = self.GtG[n] + np.einsum("i,ikl->kl", wm, self.C[n]) + np.einsum("i,ilk->kl", wp, self.C[n])
        return 0.5 * (K + K.T)

    def cross_gram(self, A, B):
        """K entre dos listas de atomos [(n, t, I)]: (na, 6, nb, 6) (BlockDict.cross_gram con I)."""
        nA = np.array([a[0] for a in A], int)
        nB = np.array([b[0] for b in B], int)
        IA = np.stack([a[2] for a in A])
        IB = np.stack([b[2] for b in B])
        lag = IA[:, None, :, None, :] - IB[None, :, None, :, :]
        W = M.ac_at(lag)
        return np.einsum("apick,bqicl,abpqi->akbl", self.BD.Gb[nA], self.BD.Gb[nB], W, optimize=True)

    def gram(self, S):
        n = len(S)
        K = self.cross_gram(S, S).reshape(6 * n, 6 * n)
        return 0.5 * (K + K.T)

    def corr_atoms(self, R, atoms):
        """h = A^T R de una lista de atomos: (na, 6)."""
        n = np.array([a[0] for a in atoms], int)
        I = np.stack([a[2] for a in atoms])                                      # (na, 2, ns)
        V = R[self.stv[None, None, :], :, I]                                     # (na, 2, ns, 3)
        return np.einsum("bx,bxk->bk", V.reshape(len(atoms), -1), self.BD.Gf[n])

    def synth_sub(self, R, n, I, u):
        pred = (self.BD.Gf[n] @ u).reshape(2, self.ns, 3)
        for b in range(2):
            for i in range(self.ns):
                t0 = int(I[b, i])
                lo, hi = max(0, t0 - AC_HALF), min(NT, t0 + AC_HALF + 1)
                R[i, :, lo:hi] -= np.outer(pred[b, i], AC[AC_HALF - (t0 - lo):AC_HALF + (hi - t0)])

    def residual(self, MFn, atoms, U):
        R = MFn.copy()
        for a, u in zip(atoms, U):
            self.synth_sub(R, a[0], a[2], u)
        return R


def atom_key(a):
    return (int(a[0]), int(a[1]), tuple(int(x) for x in np.asarray(a[2]).ravel()))


# ============================================================================== BOMP-PK
def bomp_pk_record(PD, MFn, truth, kmax=5, min_sep=5.0, peak=True, diag=None):
    """BOMP+SC con lectura en el pico (peak=False reproduce task2c_methods.bomp_record)."""
    garr = PD.BD.garr
    excluded = np.zeros(PD.N, bool)
    S, rec, run = [], [], np.inf
    R = MFn.copy()
    U = None
    shifts = []
    for _ in range(kmax):
        e, I, _h = PD.energy(R, peak)
        e[excluded] = -np.inf
        j = int(np.argmax(e))
        n, t = divmod(j, PD.T)
        s = float(e[n, t]) / 2.0
        if not np.isfinite(s):
            break
        run = min(run, s)
        src = garr[n]
        rec.append((run, M.hit_of(src, truth)))
        a = (n, t, I[n, t].copy())
        shifts.append((I[n, t] - PD.BD.idx[n, t]).ravel().tolist())
        S.append(a)
        excluded = M.excl_mask(excluded, garr, src, min_sep)
        K = PD.gram(S)
        b = PD.corr_atoms(MFn, S).ravel()
        U, *_ = np.linalg.lstsq(K, b, rcond=None)
        U = U.reshape(len(S), 6)
        R = PD.residual(MFn, S, U)
        if diag is not None:
            hS = PD.corr_atoms(R, S)
            diag.setdefault("orth", []).append(float(np.abs(hS).max() / max(1e-300, np.abs(b).max())))
    if diag is not None:
        diag["atoms"], diag["U"], diag["shifts"] = S, U, shifts
    return rec


# ============================================================================== MTSTACK-PK
def mtstack_pk_record(PD, ST, MFn, truth, min_sep_km=20.0, max_events=5, peak=True):
    """task2c_methods.Stacker.record con v leido en las muestras de pico de MFn (peak=False: 2C)."""
    I = PD.read_idx(MFn, peak)
    v = PD.sample_idx(MFn, I).reshape(PD.N, PD.T, -1)
    m = np.matmul(v, np.swapaxes(ST.Pinv, 1, 2))
    pred = np.matmul(m, np.swapaxes(PD.BD.Gf, 1, 2))
    F = np.sum(pred ** 2, axis=-1) / 2.0
    return greedy_record(list(PD.BD.garr), F.max(axis=1), truth, min_sep_km=min_sep_km, max_events=max_events)


# ============================================================================== GL-FW-PK
class GLPathPK(M.GLPath):
    """task2c_methods.GLPath con atomos (n, t, I) de tiempos fijos y oraculo con lectura en el pico."""

    def __init__(self, PD, MFn, yy, peak=True):
        super().__init__(PD.BD, MFn, yy)
        self.PD = PD
        self.peak = peak
        self.keys = set()

    def add(self, atom):
        PD = self.PD
        k_new = PD.cross_gram(self.S + [atom], [atom]).reshape(6 * (len(self.S) + 1), 6)
        m = 6 * len(self.S)
        KS = np.zeros((m + 6, m + 6))
        KS[:m, :m] = self.KS
        KS[:, m:] = k_new
        KS[m:, :] = k_new.T
        KS[m:, m:] = 0.5 * (KS[m:, m:] + KS[m:, m:].T)
        self.KS = KS
        self.bS = np.concatenate([self.bS, PD.corr_atoms(self.MFn, [atom])[0]])
        self.S.append(atom)
        self.keys.add(atom_key(atom))
        self.w = np.concatenate([self.w, np.zeros(6)])
        self._whiten()

    def drop_zeros(self):
        super().drop_zeros()
        self.keys = {atom_key(a) for a in self.S}

    def _T(self, atom):
        """(Tw, Twinv) del atomo desde su K(I) (la construccion de BlockDict: w = Tw u, Tw^T Tw = K)."""
        n, t, I = atom
        if not self.peak:
            return self.BD.Tw[n, t], self.BD.Twinv[n, t]
        K = self.PD.K_atom(n, I)
        ev, V = np.linalg.eigh(K)
        keep = ev > M.EIG_FLOOR * ev[-1]
        sq = np.where(keep, np.sqrt(np.where(keep, ev, 1.0)), 0.0)
        isq = np.where(keep, 1.0 / np.where(keep, np.sqrt(np.where(keep, ev, 1.0)), 1.0), 0.0)
        return np.einsum("j,kj->jk", sq, V), np.einsum("kj,j->kj", V, isq)

    def certificates(self):
        R = self.PD.residual(self.MFn, self.S, self.u())
        e, I, _h = self.PD.energy(R, self.peak)
        c = np.sqrt(np.maximum(e, 0.0))
        hS = self.PD.corr_atoms(R, self.S) if self.S else np.zeros((0, 6))
        cS = np.array([np.sqrt(max(0.0, float(h @ np.linalg.pinv(self.PD.K_atom(a[0], a[2])) @ h)))
                       for a, h in zip(self.S, hS)]) if self.S else np.zeros(0)
        return c, I, hS, cS

    def gap(self, lam, cmax, hS):
        U = self.u().ravel()
        s = min(1.0, lam / cmax) if cmax > 0 else 1.0
        rr = self.yy - 2.0 * U @ self.bS + U @ self.KS @ U
        g = 0.5 * (1.0 - s) ** 2 * rr + lam * self.omega() - s * (U @ hS.ravel())
        dec = -(0.5 * U @ self.KS @ U - U @ self.bS + lam * self.omega())
        return float(g), float(dec), float(rr)

    def solve(self, lam, max_fw=MAX_FW_PER_LAM, max_stall=5):
        n_fw, sweeps, stall = 0, 0, 0
        while True:
            if self.S:
                sw, _nn = M.solve_restricted(self.Kp, self.bp, self.w, lam)
                sweeps += sw
            c, I, hS, cS = self.certificates()
            j = int(np.argmax(c))
            n, t = divmod(j, self.BD.T)
            if c[n, t] <= lam * (1.0 + M.TOL_FW) or n_fw >= max_fw:
                return c, hS, cS, n_fw, sweeps, stall
            atom = (n, t, I[n, t].copy())
            if atom_key(atom) in self.keys:
                stall += 1
                if stall > max_stall:
                    return c, hS, cS, n_fw, sweeps, stall
                continue
            self.add(atom)
            n_fw += 1


def path_record_pk(garr, entry, truth, kmax=5, min_sep=5.0, with_nv=False):
    order = sorted(entry.items(), key=lambda kv: (-kv[1][0], -kv[1][1], kv[0]))
    acc, rec = [], []
    for key, (lam, nv) in order:
        n = key[0]
        src = garr[n]
        if any(np.linalg.norm(src[:2] - a[:2]) <= min_sep for a in acc):
            continue
        acc.append(src)
        rec.append((lam, M.hit_of(src, truth), nv, int(n), int(key[1])) if with_nv else (lam, M.hit_of(src, truth)))
        if len(acc) >= kmax:
            break
    return rec, order


def glfw_pk_record(PD, MFn, D, sd, truth, kmax=5, min_sep=5.0, n_lam=M.N_LAM, ratio=M.LAM_RATIO, diag=None,
                   extra=0, peak=True):
    """GL-FW+SC con oraculo de pico y atomos de tiempos fijos (peak=False reproduce task2c_methods.glfw_record)."""
    yy = M.E_R * float(np.sum(D.astype(float) ** 2)) / sd ** 2
    P = GLPathPK(PD, MFn, yy, peak)
    e0, _I, _h = PD.energy(MFn, peak)
    lam_max = float(np.sqrt(np.maximum(e0, 0.0)).max())
    lams = lam_max * ratio ** (-np.arange(n_lam) / (n_lam - 1))
    entry = {}
    per = []
    stopped = None
    garr = PD.BD.garr
    for k, lam in enumerate(lams):
        c, hS, cS, n_fw, sweeps, stall = P.solve(lam)
        nrm = [float(np.linalg.norm(P.w[6 * j:6 * j + 6])) for j in range(len(P.S))]
        for a, nv in zip(P.S, nrm):
            key = atom_key(a)
            if nv > 0.0 and key not in entry:
                entry[key] = (float(lam), nv)
        cmax = max(float(c.max()), float(cS.max()) if len(cS) else 0.0)
        g, dec, rr = P.gap(lam, cmax, hS)
        onsup = [cs / lam for cs, nv in zip(cS, nrm) if nv > 0]
        per.append({"lam": float(lam), "gap": g, "rel_gap": (g / dec if dec > 0 else (0.0 if g == 0 else np.inf)),
                    "gap_over_lam_omega": (g / (lam * P.omega()) if P.omega() > 0 else 0.0),
                    "n_support": int(sum(nv > 0 for nv in nrm)), "n_fw": n_fw, "sweeps": int(sweeps), "stall": int(stall),
                    "kkt_off_max": float(c.max() / lam),
                    "kkt_on_maxdev": float(max(abs(x - 1.0) for x in onsup)) if onsup else 0.0,
                    "fw_capped": int(n_fw >= MAX_FW_PER_LAM)})
        P.drop_zeros()
        if stopped is None and len(path_record_pk(garr, entry, truth, kmax, min_sep)[0]) >= kmax:
            stopped = ("record_complete", k)
        if stopped is not None and k >= stopped[1] + extra:
            break
    rec, order = path_record_pk(garr, entry, truth, kmax, min_sep)
    if diag is not None:
        diag["per_lambda"] = per
        diag["lam_max"] = lam_max
        diag["n_entered"] = len(entry)
        diag["stopped"] = stopped
        diag["entry"] = entry
        diag["rec_nv"] = path_record_pk(garr, entry, truth, kmax, min_sep, with_nv=True)[0]
    return rec
