"""
c6_lib.py - C6 (physics candidate): 1D velocity models inside the pick-free moment-tensor detectors.  SEED 999 ONLY.

Velocity models (all return travel times (N, 2, ns) for sources (N, 3) and stations (ns, 3) at z = 0, and ray
directions at the source / station for the radiation operator):
  Homog(vs, vpvs)           straight rays (the detectors' registered model: vs = 3.0, vpvs = 1.73)
  Gradient(vs0, g, vpvs)    v(z) = v0 + g z, circular rays: gen_realistic.travel_time_gradient / _ray (the generator's
                            own medium when (2.8, 0.04, 1.73))
  TwoLayer(v1, v2, h, vpvs) layer over half-space, direct (refracted) ray by vectorised bisection on the ray parameter
Radiation operator G(src) (2*ns*3 x 6), rows (b, i, c) as task1x_mtglrt.radiation_operator:
  P: n_r (c(n_s) . m) AP / (d sqrt(CP));  S: (e_SV,r e_SV,s^T + e_SH e_SH^T) B(n_s) m AS / (d sqrt(CS)),
the generator's radiation_gradient written as a linear map of m (1/d straight-line spreading as in the generator).
With Homog(3.0, 1.73) both reduce exactly to the registered operators (checked in c6_validate.py).

Detectors with a model:
  make_BD(STA, grid, model, corr)       task2c_methods.BlockDict with the model's sample indices and operators
  bomp_pk_m(PD, MFn, truth, kmax)       task2d_methods.bomp_pk_record (unchanged) -> record, atom locations
  mtfw_dense(ops, MFn, truth, kmax, cap) MT-FW (task1x fw_record logic: capped peak vectors, ||Q^T v||^2/2, sequential
                                        deflation at the peak samples, 5 km exclusion) evaluated densely over the
                                        whole (node, origin) lattice with the model's indices; tie rule of fw_record
                                        (lowest node index, first origin).  cap=False gives MT-FW-NOCAP.
Calibration from the catalogue (K = 20 events, locations +-1 km, origins +-0.1 s):
  pick_catalogue(STA, cat, model)       joint P/S matched-filter energy peaks within +-1 s of the model prediction
                                        (station_corrections.peak_delays logic), gate 25
  jhd_fit(STA, cat, family, ...)        joint fit of the model parameters, station/phase terms, event relocations
                                        (prior 1 km) and origin shifts (prior 0.1 s), soft-L1 loss, 2 pick rounds
  corr_model(STA, cat, model)           station corrections = median delay w.r.t. the model (the registered +SC rule)
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys

import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
for p in (HERE, WORK):
    if p not in sys.path:
        sys.path.insert(0, p)

import gen_realistic as GR  # noqa: E402
import task1x_mtglrt as X  # noqa: E402
import task2c_methods as M  # noqa: E402
import task2d_methods as D2  # noqa: E402
from station_corrections import MIN_PS  # noqa: E402
from task0_common import AP, AS, DEFAULT_TAUS, DT, NT, VP, VS  # noqa: E402
from task1o_consistent import TOL_KM, matched  # noqa: E402
from task1w_mech import CP, CS  # noqa: E402

EZ = np.array([0.0, 0.0, 1.0])


# ============================================================================================ velocity models
class Homog:
    family = "homog"

    def __init__(self, vs=VS, vpvs=VP / VS):
        self.vs, self.vpvs = float(vs), float(vpvs)
        self.params = (self.vs, self.vpvs)

    def v(self, b):
        return self.vs * (self.vpvs if b == 0 else 1.0)

    def tt(self, SRC, STA):
        d = np.linalg.norm(np.asarray(STA)[None] - np.atleast_2d(SRC)[:, None], axis=2)      # (N, ns)
        return np.stack([d / self.v(0), d / self.v(1)], axis=1)

    def dirs(self, src, sta, b):
        n = (np.asarray(sta, float) - src)
        n = n / np.linalg.norm(n)
        dxy = sta[:2] - src[:2]
        X_ = float(np.hypot(*dxy))
        eh = None if X_ < 1e-9 else np.array([dxy[0] / X_, dxy[1] / X_, 0.0])
        return n, n, eh


class Gradient:
    family = "gradient"

    def __init__(self, vs0=2.8, g=0.04, vpvs=1.73):
        self.vs0, self.g, self.vpvs = float(vs0), float(g), float(vpvs)
        self.params = (self.vs0, self.g, self.vpvs)
        self.cfg = {"VS0": self.vs0, "g": self.g, "vpvs": self.vpvs}

    def vg(self, b):
        k = self.vpvs if b == 0 else 1.0
        return k * self.vs0, k * self.g

    def tt(self, SRC, STA):
        SRC = np.atleast_2d(np.asarray(SRC, float))
        STA = np.asarray(STA, float)
        out = np.zeros((len(SRC), 2, len(STA)))
        r = np.linalg.norm(STA[None] - SRC[:, None], axis=2)
        for b in range(2):
            v0, g = self.vg(b)
            q = np.sqrt((v0 + g * SRC[:, 2:3]) * (v0 + g * STA[None, :, 2]))
            x = g * r / (2.0 * q)
            xs = np.where(x > 1e-6, x, 1.0)
            f = np.where(x > 1e-6, np.arcsinh(xs) / xs, 1.0 - x * x / 6.0)
            out[:, b] = r / q * f
        return out

    def dirs(self, src, sta, b):
        v0, g = self.vg(b)
        return GR._ray(src, sta, v0, g)


class TwoLayer:
    """S velocity v1 in 0 <= z < h, v2 >= v1 below; P = vpvs * S.  Direct (refracted) ray only."""
    family = "twolayer"

    def __init__(self, v1=3.0, v2=3.3, h=5.0, vpvs=1.73):
        self.v1, self.v2, self.h, self.vpvs = float(v1), float(max(v2, v1)), float(h), float(vpvs)
        self.params = (self.v1, self.v2, self.h, self.vpvs)

    def _solve(self, X_, zs, b):
        """X_, zs arrays: sin of the angle in layer 2 (u), travel time; straight ray when zs <= h."""
        k = self.vpvs if b == 0 else 1.0
        v1, v2, h = k * self.v1, k * self.v2, self.h
        X_ = np.asarray(X_, float)
        zs = np.asarray(zs, float)
        deep = zs > h
        z2 = np.where(deep, zs - h, 0.0)
        z1 = np.where(deep, h, zs)
        rr = v1 / v2
        lo = np.zeros_like(X_)
        hi = np.full_like(X_, 1.0 - 1e-12)
        for _ in range(60):
            u = 0.5 * (lo + hi)
            s1 = u * rr
            Xu = z2 * u / np.sqrt(1 - u * u) + z1 * s1 / np.sqrt(1 - s1 * s1)
            big = Xu > X_
            hi = np.where(big, u, hi)
            lo = np.where(big, lo, u)
        u = 0.5 * (lo + hi)
        s1 = u * rr
        t_deep = z2 / (v2 * np.sqrt(1 - u * u)) + z1 / (v1 * np.sqrt(1 - s1 * s1))
        d_str = np.sqrt(X_ ** 2 + zs ** 2)
        t = np.where(deep, t_deep, d_str / v1)
        u_str = np.where(d_str > 0, X_ / np.maximum(d_str, 1e-12), 0.0)       # sin of the straight-ray angle
        s_src = np.where(deep, u, u_str)
        s_rec = np.where(deep, s1, u_str)
        return t, s_src, s_rec

    def tt(self, SRC, STA):
        SRC = np.atleast_2d(np.asarray(SRC, float))
        STA = np.asarray(STA, float)
        X_ = np.linalg.norm(STA[None, :, :2] - SRC[:, None, :2], axis=2)
        zs = np.broadcast_to(SRC[:, 2:3], X_.shape)
        return np.stack([self._solve(X_, zs, b)[0] for b in range(2)], axis=1)

    def dirs(self, src, sta, b):
        dxy = sta[:2] - src[:2]
        X_ = float(np.hypot(*dxy))
        if X_ < 1e-9:
            n = np.array([0.0, 0.0, -1.0])
            return n, n.copy(), None
        eh = np.array([dxy[0] / X_, dxy[1] / X_, 0.0])
        _t, ss, sr = self._solve(np.array([X_]), np.array([src[2]]), b)
        ss, sr = float(ss[0]), float(sr[0])
        ns_ = ss * eh - np.sqrt(1 - ss * ss) * EZ
        nr_ = sr * eh - np.sqrt(1 - sr * sr) * EZ
        return ns_, nr_, eh


FAMILIES = {"homog": Homog, "gradient": Gradient, "twolayer": TwoLayer}


# ============================================================================================ radiation operator
def _Bmat(n):
    nx, ny, nz = n
    return np.array([[nx, 0, 0, ny, nz, 0], [0, ny, 0, nx, 0, nz], [0, 0, nz, 0, nx, ny]])


def _cvec(n):
    nx, ny, nz = n
    return np.array([nx * nx, ny * ny, nz * nz, 2 * nx * ny, 2 * nx * nz, 2 * ny * nz])


def radiation_operator_m(src, STA, model):
    src = np.asarray(src, float)
    rows = []
    for b in range(2):
        for i in range(len(STA)):
            sta = np.asarray(STA[i], float)
            d = float(np.linalg.norm(sta - src))
            ns_, nr_, eh = model.dirs(src, sta, b)
            if b == 0:
                g = np.outer(nr_, _cvec(ns_)) * AP / (d * np.sqrt(CP))
            else:
                B = _Bmat(ns_)
                if eh is None:
                    Pm = np.eye(3) - np.outer(ns_, ns_)
                else:
                    esh = np.cross(EZ, eh)
                    esv_s, esv_r = np.cross(esh, ns_), np.cross(esh, nr_)
                    Pm = np.outer(esv_r, esv_s) + np.outer(esh, esh)
                g = Pm @ B * AS / (d * np.sqrt(CS))
            rows.append(g)
    return np.vstack(rows)


# ============================================================================================ operators per network x model
def idx_table(STA, grid, model, corr, taus=DEFAULT_TAUS):
    """(N, T, 2, ns) sample indices int((tt + tau + corr) / DT), clipped (certified_fw._idx_node with the model)."""
    TT = model.tt(np.asarray(grid), STA)                                          # (N, 2, ns)
    tt = TT[:, None] + np.asarray(taus)[None, :, None, None]
    if corr is not None:
        tt = tt + np.asarray(corr).T[None, None]
    return np.clip((tt / DT).astype(int), 0, NT - 1)


class ModelOps:
    def __init__(self, STA, grid, model, corr, taus=DEFAULT_TAUS):
        self.STA = np.asarray(STA, float)
        self.grid = grid
        self.garr = np.asarray(grid, float)
        self.model, self.corr = model, np.asarray(corr, float)
        self.taus = np.asarray(taus, float)
        self.Gs = [radiation_operator_m(g, self.STA, model) for g in self.garr]
        self.Qs = [X.basis(G) for G in self.Gs]
        from certified_fw import stack_Q
        self.Qst = stack_Q(self.Qs)
        self.IDX = idx_table(self.STA, grid, model, self.corr, self.taus)
        tab = {tuple(np.round(g, 9)): k for k, g in enumerate(self.garr)}
        old = M._idx_node

        def _idx(STA_, src, taus_, tcorr_):
            return self.IDX[tab[tuple(np.round(np.asarray(src, float), 9))]]
        M._idx_node = _idx
        try:
            self.BD = M.BlockDict(self.STA, grid, self.taus, self.corr, Gs=self.Gs)
        finally:
            M._idx_node = old
        self.PD = D2.PKDict(self.BD)


def bomp_pk_m(ops, MFn, truth, kmax=5):
    d = {}
    rec = D2.bomp_pk_record(ops.PD, MFn, truth, kmax=kmax, diag=d)
    locs = [ops.garr[a[0]].tolist() + [float(ops.taus[a[1]])] for a in d["atoms"]]
    return [(float(s), int(h)) for s, h in rec], locs


def _leaf(R, cap=True):
    Et = np.sum(np.ascontiguousarray(R.transpose(0, 2, 1)) ** 2, axis=-1)          # (ns, NT)
    W = Et[:, D2.II_TAB]
    a = W.argmax(axis=-1)
    K = D2.II_TAB[np.arange(NT)[None, :], a]                                        # (ns, NT)
    V = np.ascontiguousarray(R.transpose(0, 2, 1))                                  # (ns, NT, 3)
    if cap:
        nv = np.linalg.norm(V, axis=-1, keepdims=True)
        V = V * np.minimum(1.0, X.CAPV / np.maximum(nv, 1e-12))
    return K, V


def mtfw_dense(ops, MFn, truth, kmax=5, min_sep=5.0, cap=True):
    STA, ns = ops.STA, len(ops.STA)
    R = MFn.copy()
    garr, IDX = ops.garr, ops.IDX
    N, T = IDX.shape[:2]
    st = np.arange(ns)[None, None, None, :]
    excluded = np.zeros(N, bool)
    rec, run, locs = [], np.inf, []
    for _ in range(kmax):
        if excluded.all():
            break
        K, V = _leaf(R, cap)
        tpk = K[st, IDX]                                                            # (N, T, 2, ns)
        vf = V[st, tpk].reshape(N, T, 6 * ns)
        y = np.einsum("gtr,grk->gtk", vf, ops.Qst)
        sc = (y * y).sum(-1) / 2.0
        sc[excluded] = -np.inf
        best = sc.max(axis=1)
        gi = int(np.argmax(best))                     # lowest node index among ties
        tj = int(np.argmax(sc[gi]))                   # first origin among ties
        s = float(sc[gi, tj])
        if not np.isfinite(s):
            break
        src = garr[gi]
        run = min(run, s)
        rec.append((float(run), int(M.hit_of(src, truth))))
        locs.append(src.tolist() + [float(ops.taus[tj])])
        excluded = excluded | (np.linalg.norm(garr[:, :2] - src[:2], axis=1) <= min_sep)
        G = ops.Gs[gi]
        m, *_ = np.linalg.lstsq(G, vf[gi, tj], rcond=None)
        pred = (G @ m).reshape(2, ns, 3)
        tp = tpk[gi, tj]
        for b in range(2):
            for i in range(ns):
                t0 = int(tp[b, i])
                lo, hi = max(0, t0 - X.AC_HALF), min(NT, t0 + X.AC_HALF + 1)
                R[i, :, lo:hi] -= np.outer(pred[b, i], X.AC[X.AC_HALF - (t0 - lo):X.AC_HALF + (hi - t0)])
    return rec, locs


# ============================================================================================ catalogue picks, corrections
def cat_energy(cat):
    out = []
    for D, _src, _tau in cat:
        MF, sd = matched(D)
        out.append(np.sum((MF / sd) ** 2, axis=1))
    return out


def pick_event(E, STA, src, tau, model, search=1.0, center_corr=None):
    """station_corrections.peak_delays (joint P/S) around the model prediction (+ center_corr): absolute peak
    times (ns, 2), predicted times without corrections (ns, 2), energies (ns, 2)."""
    tp_all = model.tt(np.asarray(src, float)[None], STA)[0].T + tau                  # (ns, 2)
    W = int(round(search / DT))
    ns = len(STA)
    tobs = np.zeros((ns, 2))
    en = np.zeros((ns, 2))
    for i in range(ns):
        win = []
        for b in range(2):
            c0 = tp_all[i, b] + (0.0 if center_corr is None else center_corr[i, b])
            c = int(round(c0 / DT))
            win.append(np.arange(max(1, c - W), min(NT - 2, c + W) + 1))
        sc = E[i, win[0]][:, None] + E[i, win[1]][None, :]
        sc = np.where(win[1][None, :] >= win[0][:, None] + MIN_PS, sc, -np.inf)
        a, b_ = np.unravel_index(int(np.argmax(sc)), sc.shape)
        ks = (int(win[0][a]), int(win[1][b_]))
        for b in range(2):
            k = ks[b]
            e0, e1, e2 = E[i, k - 1], E[i, k], E[i, k + 1]
            den = e0 - 2 * e1 + e2
            dk = 0.5 * (e0 - e2) / den if den < 0 else 0.0
            tobs[i, b] = (k + float(np.clip(dk, -0.5, 0.5))) * DT
            en[i, b] = e1
    return tobs, tp_all, en


def corr_model(STA, cat, model, Ecat=None, gate=25.0):
    """registered +SC rule re-referenced to the model: median over catalogue events of (peak - model prediction)."""
    Ecat = Ecat if Ecat is not None else cat_energy(cat)
    ns = len(STA)
    obs = [[[] for _ in range(2)] for _ in range(ns)]
    for E, (_D, src, tau) in zip(Ecat, cat):
        tobs, tp, en = pick_event(E, STA, src, tau, model)
        for i in range(ns):
            for b in range(2):
                if en[i, b] >= gate:
                    obs[i][b].append(tobs[i, b] - tp[i, b])
    c = np.zeros((ns, 2))
    for i in range(ns):
        for b in range(2):
            if obs[i][b]:
                c[i, b] = float(np.median(obs[i][b]))
    return c


BOUNDS = {"homog": ([2.0, 1.5], [4.5, 2.1]),
          "gradient": ([2.0, 0.0, 1.5], [4.5, 0.15, 2.1]),
          "twolayer": ([2.0, 0.0, 1.0, 1.5], [4.5, 2.0, 15.0, 2.1])}
INIT = {"homog": [VS, VP / VS], "gradient": [VS, 0.0, VP / VS], "twolayer": [VS, 0.3, 5.0, VP / VS]}


def _mk(family, th):
    if family == "twolayer":
        return TwoLayer(th[0], th[0] + th[1], th[2], th[3])
    return FAMILIES[family](*th)


def jhd_fit(STA, cat, family, Ecat=None, gate=25.0, sig_d=0.05, cat_err=(1.0, 0.1), rounds=2, return_info=False):
    """joint fit (model params, station/phase terms, event dx (3) and dt0) to the catalogue peak times."""
    Ecat = Ecat if Ecat is not None else cat_energy(cat)
    ns, K = len(STA), len(cat)
    npar = len(INIT[family])
    th = np.array(INIT[family], float)
    lo_b, hi_b = BOUNDS[family]
    info = {}
    x = None
    for rnd in range(rounds):
        model = _mk(family, th)
        rows = []                                          # (k, i, b, tobs)
        for k, (E, (_D, src, tau)) in enumerate(zip(Ecat, cat)):
            cc = None
            if x is not None:
                cc = x[npar:npar + 2 * ns].reshape(ns, 2) + x[npar + 2 * ns + 3 * K + k]
            tobs, _tp, en = pick_event(E, STA, np.asarray(src, float), float(tau), model, center_corr=cc)
            for i in range(ns):
                for b in range(2):
                    if en[i, b] >= gate:
                        rows.append((k, i, b, tobs[i, b]))
        rows = np.array(rows)
        kk, ii, bb, tob = rows[:, 0].astype(int), rows[:, 1].astype(int), rows[:, 2].astype(int), rows[:, 3]
        src0 = np.array([np.asarray(s, float) for _D, s, _t in cat])
        tau0 = np.array([float(t) for _D, _s, t in cat])

        def resid(p):
            mdl = _mk(family, p[:npar])
            s = p[npar:npar + 2 * ns].reshape(ns, 2)
            dx = p[npar + 2 * ns:npar + 2 * ns + 3 * K].reshape(K, 3)
            dt = p[npar + 2 * ns + 3 * K:]
            SRC = src0 + dx
            SRC[:, 2] = np.clip(SRC[:, 2], 0.5, 25.0)
            TT = mdl.tt(SRC, STA)                                           # (K, 2, ns)
            pred = TT[kk, bb, ii] + tau0[kk] + dt[kk] + s[ii, bb]
            return np.concatenate([(tob - pred) / sig_d, dx.ravel() / cat_err[0], dt / cat_err[1], s.ravel() / 1.0])

        p0 = np.concatenate([th, np.zeros(2 * ns + 4 * K)]) if x is None else x.copy()
        lo = np.concatenate([lo_b, np.full(2 * ns, -3.0), np.full(3 * K, -6.0), np.full(K, -1.0)])
        hi = np.concatenate([hi_b, np.full(2 * ns, 3.0), np.full(3 * K, 6.0), np.full(K, 1.0)])
        p0 = np.clip(p0, lo + 1e-9, hi - 1e-9)
        sol = least_squares(resid, p0, bounds=(lo, hi), loss="soft_l1", f_scale=2.0, x_scale="jac", max_nfev=200)
        x = sol.x
        th = x[:npar]
        r = resid(x)[:len(tob)] * sig_d
        info[f"round{rnd}"] = {"n_obs": int(len(tob)), "theta": th.tolist(), "rms": float(np.sqrt(np.mean(r ** 2))),
                               "mad": float(np.median(np.abs(r))), "nfev": int(sol.nfev)}
    model = _mk(family, th)
    return (model, info) if return_info else model
