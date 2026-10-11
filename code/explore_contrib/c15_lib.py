"""
c15_lib.py - C15: the origin time of a timing-tolerant detector is set-valued; second-stage estimation of (x, tau)
from the samples the detector itself read.  SEED 999 ONLY.

Fact (exact, from the hard rule of c10_lib.choose): the statistic of a hypothesis (node g, origin tau) depends on the
data only through the samples tp_bi read at the energy peak within +-W samples of the predicted sample of every
station/phase.  Shifting tau by less than the slack of every read window leaves every tp_bi, hence the statistic,
unchanged.  The maximiser in tau is therefore an interval (a plateau), and np.argmax returns its EARLIEST point.  The
registered atom reports that point.  C14 found that 55-58% (STD) and 80-95% (DENSE) of the false atoms at the operating
point under LOC+T are well-placed real events with the origin off by a median 0.8 s, mostly early.

The detection sequence (selection, scores, exclusion, deflation) is NOT changed here: every estimator below only
replaces the reported (x, tau) of an accepted atom.  mt_loop15 with no estimator reproduces c10_lib.mt_loop exactly.

Estimators (truth-free; inputs: the node, the plateau, the read samples tp, their energies):
  MID    tau = midpoint of the plateau {tau_j : s(g, tau_j) = max} (the Chebyshev centre of the argmax set).
  LSO    tau = tau_0 + weighted median of the read delays d_bi = (tp_bi - x_bi(tau_0)) DT, x_bi the continuous predicted
         sample; weights = capped read energy minus the noise mean (|v|^2 - 3, floored at 0).
  LSOs   LSO with sub-sample (parabolic) interpolation of the read peaks.
  LSX    Gauss-Newton step on (x, y, z, tau) from the read delays (sub-sample), robust (Huber, 3 IRLS passes), with a
         Gaussian prior centred on the node: sd h/sqrt(12) on x, y, hz on z, 1 s on tau; delay noise sd SIG_D.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
for p in (HERE, WORK):
    if p not in sys.path:
        sys.path.insert(0, p)

import c10_lib as C  # noqa: E402
from task0_common import DT, NT  # noqa: E402

SIG_D = 0.10                    # s, delay noise level of LSX (fixed a priori, not tuned)
HUB = 1.345


def _plateau(row, s, eps=1e-9):
    return np.where(row >= s - eps * max(1.0, abs(s)))[0]


def _subsample(Et, tp):
    """parabolic offset (samples) of the energy peak at tp, per station-phase; Et (ns, NT) energy."""
    ns = Et.shape[0]
    off = np.zeros(tp.shape)
    for b in range(tp.shape[0]):
        for i in range(ns):
            t = int(tp[b, i])
            if 0 < t < NT - 1:
                a, c, d = Et[i, t - 1], Et[i, t], Et[i, t + 1]
                den = a - 2.0 * c + d
                if den < 0:
                    off[b, i] = float(np.clip(0.5 * (a - d) / den, -0.5, 0.5))
    return off


def _wmedian(x, w):
    o = np.argsort(x)
    cw = np.cumsum(w[o])
    return float(x[o][np.searchsorted(cw, 0.5 * cw[-1])])


def mt_loop15(E, MFn, kmax, capv=C.CAPV, min_sep=5.0, W=C.SH, keep_et=False):
    """c10_lib.mt_loop (hard rule, no refinement) recording, per atom, the node, plateau, read samples and energies."""
    ns = E.ns
    R = MFn.copy()
    excluded = np.zeros(E.N, bool)
    atoms = []
    for _ in range(kmax):
        if excluded.all():
            break
        Eh = C._half_energy(R, capv)
        tpk, pen = C.choose(E, Eh, "hard", None, W, R=R)
        vf = C._vectors(R, tpk, capv).reshape(E.N, E.T, 6 * ns)
        y = np.einsum("gtr,grk->gtk", vf, E.Q)
        s0 = (y * y).sum(-1) / 2.0
        sc = s0.copy()
        sc[excluded] = -np.inf
        gi = int(np.argmax(sc.max(axis=1)))
        tj = int(np.argmax(sc[gi]))
        s = float(sc[gi, tj])
        if not np.isfinite(s):
            break
        src, tau, tp, vsel, G = E.XN[gi], float(E.taus[tj]), tpk[gi, tj].copy(), vf[gi, tj], E.G[gi]
        Et = np.sum(R * R, axis=1)
        v = vsel.reshape(2, ns, 3)
        a = {"s": s, "s0": float(s0[gi, tj]), "x": np.asarray(src, float), "tau": tau, "gi": gi, "tj": tj,
             "plateau": _plateau(sc[gi], s), "tp": tp, "en": (v * v).sum(-1), "sub": _subsample(Et, tp)}
        if keep_et:
            a["Et"] = Et
        atoms.append(a)
        excluded = excluded | (np.linalg.norm(E.XN[:, :2] - src[None, :2], axis=1) <= min_sep)
        m, *_ = np.linalg.lstsq(G, vsel, rcond=None)
        C._deflate(R, (G @ m).reshape(2, ns, 3), tp)
    return atoms


def delays(E, a, sub=False):
    gi, tj = a["gi"], a["tj"]
    xc = E.IDX[gi, tj] + E.frac[gi, tj]                                   # (2, ns) continuous predicted sample
    t = a["tp"] + (a["sub"] if sub else 0.0)
    return (t - xc) * DT


def weights(a):
    return np.maximum(a["en"] - 3.0, 0.0)


def est_mid(E, a):
    return a["x"], float(np.median(E.taus[a["plateau"]]))


def est_lso(E, a, sub=False):
    d, w = delays(E, a, sub).ravel(), weights(a).ravel()
    if w.sum() <= 0:
        return est_mid(E, a)
    return a["x"], a["tau"] + _wmedian(d, w)


def est_lsx(E, a, h=C.H_COARSE, hz=C.HZ):
    d, w = delays(E, a, True).ravel(), weights(a).ravel()
    if w.sum() <= 0:
        return est_mid(E, a)
    J = E.J[a["gi"]]                                                      # (2ns, 4): d t / d (x, y, z, tau)
    P = np.diag(1.0 / np.array([h * h / 12.0, h * h / 12.0, hz * hz, 1.0]))
    wn = w / w.max()
    th = np.zeros(4)
    th[3] = _wmedian(d, w)
    for _ in range(3):
        r = d - J @ th
        u = np.abs(r) / SIG_D
        hw = np.where(u <= HUB, 1.0, HUB / np.maximum(u, 1e-12))
        Wd = wn * hw / SIG_D ** 2
        A = J.T @ (Wd[:, None] * J) + P
        th = np.linalg.solve(A, J.T @ (Wd * d))
    x = a["x"] + th[:3]
    x[2] = float(np.clip(x[2], E.XN[:, 2].min() - hz / 2, E.XN[:, 2].max() + hz / 2))
    return x, a["tau"] + th[3]


# ------------------------------------------------------------------------------------------ EDT relocation
# Equal-differential-time likelihood (Font et al. 2004; Lomax 2005, NonLinLoc EDT): origin-free and robust to
# outlying reads.  score(x) = sum_{i<j} w_i w_j exp(-((r_i - r_j)^2) / (2 SIG_EDT^2)), r = read time - tt(x) - corr.
SIG_EDT = 0.2                   # s, fixed a priori
EDT_XY = np.arange(-6.0, 6.0 + 1e-9, 0.5)
EDT_Z = np.arange(1.0, 17.0 + 1e-9, 1.0)


def engine_corr(E):
    """registered station corrections (2, ns) recovered from the engine tables (interior origin index)."""
    t = E.T // 2
    d = np.linalg.norm(E.STA[None] - E.XN[:1, None], axis=-1)[0]
    tt = np.stack([d * C.VINV[0], d * C.VINV[1]])
    return (E.IDX[0, t] + E.frac[0, t]) * DT - tt - E.taus[t]


def _edt_grid(E, x0):
    gx, gy, gz = np.meshgrid(x0[0] + EDT_XY, x0[1] + EDT_XY, EDT_Z, indexing="ij")
    return np.stack([gx.ravel(), gy.ravel(), gz.ravel()], 1)


def _tt(E, P, corr):
    d = np.linalg.norm(P[:, None, :] - E.STA[None], axis=-1)                    # (n, ns)
    return np.stack([d * C.VINV[0], d * C.VINV[1]], 1) + corr[None]           # (n, 2, ns)


def edt_locate(E, x0, t, w, corr):
    """t (2, ns) read times (s), w (2, ns) weights; returns (x, tau) or None."""
    m = (w > 0).ravel()
    if m.sum() < 3:
        return None
    P = _edt_grid(E, x0)
    r = (t[None] - _tt(E, P, corr)).reshape(len(P), -1)[:, m]
    ww = w.ravel()[m]
    D = r[:, :, None] - r[:, None, :]
    sc = np.einsum("i,j,nij->n", ww, ww, np.exp(-D * D / (2.0 * SIG_EDT ** 2)))
    k = int(np.argmax(sc))
    return P[k], _wmedian(r[k], ww)


def est_edt(E, a):
    corr = engine_corr(E) if not hasattr(E, "_corr15") else E._corr15
    E._corr15 = corr
    t = (a["tp"] + a["sub"]) * DT
    res = edt_locate(E, a["x"], t, weights(a), corr)
    return (a["x"], est_lso(E, a, True)[1]) if res is None else res


def est_edt2(E, a, W=C.SH, rounds=2):
    """EDT, then re-read the peaks of the pre-deflation energy within +-W samples of the times predicted at the EDT
    location and origin, and locate again (rounds times)."""
    x, tau = est_edt(E, a)
    Et = a["Et"]
    ns = Et.shape[0]
    for _ in range(rounds):
        pred = (tau + _tt(E, np.asarray(x, float)[None], E._corr15)[0]) / DT       # (2, ns) continuous sample
        base = np.floor(pred).astype(int)
        tp = np.zeros((2, ns), int)
        for b in range(2):
            for i in range(ns):
                lo, hi = max(0, base[b, i] - W), min(Et.shape[1], base[b, i] + W + 1)
                tp[b, i] = lo + int(np.argmax(Et[i, lo:hi])) if hi > lo else int(np.clip(base[b, i], 0, NT - 1))
        en = np.minimum(Et[np.arange(ns)[None], tp], C.CAPV ** 2)
        b2 = {"tp": tp, "en": en, "sub": _subsample(Et, tp)}
        res = edt_locate(E, np.asarray(x, float), (tp + b2["sub"]) * DT, weights(b2), E._corr15)
        if res is None:
            break
        x, tau = res
    return x, tau


ESTIMATORS = {"MID": est_mid, "LSO": est_lso, "LSOs": lambda E, a: est_lso(E, a, True), "LSX": est_lsx,
              "EDT": est_edt, "EDT2": est_edt2}


def pack_all(E, atoms):
    out = {"BASE": [[*map(float, a["x"]), float(a["tau"]), float(a["s"])] for a in atoms]}
    for nm, f in ESTIMATORS.items():
        out[nm] = []
        for a in atoms:
            x, t = f(E, a)
            out[nm].append([*map(float, x), float(t), float(a["s"])])
    return out


def plateau_stats(E, atoms):
    return [[len(a["plateau"]), float(E.taus[a["plateau"]].min()), float(E.taus[a["plateau"]].max())] for a in atoms]


# ------------------------------------------------------------------------------------------ cyclic re-estimation
def _best(E, R, capv, excluded, W=C.SH):
    Eh = C._half_energy(R, capv)
    tpk, _pen = C.choose(E, Eh, "hard", None, W, R=R)
    vf = C._vectors(R, tpk, capv).reshape(E.N, E.T, 6 * E.ns)
    y = np.einsum("gtr,grk->gtk", vf, E.Q)
    sc = (y * y).sum(-1) / 2.0
    sc[excluded] = -np.inf
    gi = int(np.argmax(sc.max(axis=1)))
    tj = int(np.argmax(sc[gi]))
    return gi, tj, float(sc[gi, tj]), tpk[gi, tj].copy(), vf[gi, tj]


def _contrib(E, gi, vsel, tp, ns):
    """deflation template of an atom as a full array (ns, 3, NT)."""
    G = E.G[gi]
    m, *_ = np.linalg.lstsq(G, vsel, rcond=None)
    Z = np.zeros((ns, 3, NT))
    C._deflate(Z, -(G @ m).reshape(2, ns, 3), tp)
    return Z


def mt_loop_cyc(E, MFn, kmax, cycles=2, capv=C.CAPV, min_sep=5.0, rescore=True):
    """greedy MT-FW (identical to c10_lib.mt_loop for cycles = 0), then block-coordinate re-estimation: every atom in
    turn is put back into the residual and re-selected (node, origin) with the others deflated, outside 5 km of the
    others; its deflation is replaced.  rescore: the score of an atom is its statistic in the last cycle, otherwise
    the greedy score."""
    ns = E.ns
    R = MFn.copy()
    excluded = np.zeros(E.N, bool)
    at = []
    for _ in range(kmax):
        if excluded.all():
            break
        gi, tj, s, tp, v = _best(E, R, capv, excluded)
        if not np.isfinite(s):
            break
        Z = _contrib(E, gi, v, tp, ns)
        at.append({"gi": gi, "tj": tj, "s": s, "s_greedy": s, "Z": Z})
        R -= Z
        excluded = excluded | (np.linalg.norm(E.XN[:, :2] - E.XN[gi][None, :2], axis=1) <= min_sep)
    for _c in range(cycles):
        moved = 0
        for k, a in enumerate(at):
            R += a["Z"]
            exk = np.zeros(E.N, bool)
            for q, b in enumerate(at):
                if q != k:
                    exk |= np.linalg.norm(E.XN[:, :2] - E.XN[b["gi"]][None, :2], axis=1) <= min_sep
            if exk.all():
                R -= a["Z"]
                continue
            gi, tj, s, tp, v = _best(E, R, capv, exk)
            moved += int(gi != a["gi"] or tj != a["tj"])
            a.update(gi=gi, tj=tj, s=s, Z=_contrib(E, gi, v, tp, ns))
            R -= a["Z"]
        if moved == 0:
            break
    return [{"x": np.asarray(E.XN[a["gi"]], float), "tau": float(E.taus[a["tj"]]),
             "s": a["s"] if rescore else a["s_greedy"]} for a in at]
