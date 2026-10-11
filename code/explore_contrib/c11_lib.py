"""
c11_lib.py - C11: self-consistent station terms for the pick-free moment-tensor subspace detectors.  SEED 999 ONLY.

Station terms c (ns, 2) [station, phase], fixed per network, used as the corrections of BOMP-PK / MT-FW (cap 6):
  S0   registered +SC (station_corrections.estimate_station_corrections on the 20-event catalogue).
  S1   SC-R, frozen in c10_screloc.sc_reloc.
  S2   JHD: alternation of (a) relocation of every catalogue event (position and origin time, Huber misfit of its
       peak delays with the current terms, Gaussian prior of the catalogue error around the catalogue hypocentre) and
       (b) station terms = median over events of the residuals; peaks gated by energy and RE-READ around the current
       prediction in a narrower window between passes.  Gauge: median origin-time shift of the catalogue events = 0.
  S3   self-calibration: S2 terms -> MT-FW + T3 over an unlabelled stream -> atoms above a score quantile -> their
       peaks read around the prediction of the refined atom -> joint JHD with the catalogue events -> terms.
  S4   S2 with 40 / 80 catalogue events.
All estimators use only waveforms, catalogue hypocentres and (S3) the detector's own atoms; the true statics enter only
in err() (calibration quality) and in the ORACLE arm.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import copy
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
for p in (HERE, WORK):
    if p not in sys.path:
        sys.path.insert(0, p)

import c10_lib as C  # noqa: E402
import xc  # noqa: E402
from station_corrections import MIN_PS  # noqa: E402
from task0_common import DT, NT, VP, VS  # noqa: E402

SEED = 999
VEL = np.array([VP, VS])
GATE = 25.0

# ---- constants of the estimators (values and the rules that fixed them: C11_PLAN.txt)
P = {
    "search0": 1.0,        # s, first read: the registered +SC read (around the catalogue prediction, no terms)
    "w_re": 0.5,           # s, half-width of the re-read around the current prediction
    "n_pass": 3,           # reads: the first + (n_pass - 1) re-reads
    "n_alt": 8,            # alternations per pass (stop when the terms move < tol)
    "tol": 0.002,          # s
    "huber": 1.345,
    "sig_floor": 0.05,     # s, floor of the residual scale
    "sx": 1.0,             # km, prior sd of the catalogue position per coordinate (the catalogue error)
    "st": 0.1,             # s, prior sd of the catalogue origin time
    "gate_re": 12.0,       # energy gate of the re-read arrivals (first read: the registered gate 25)
    "min_arr": 4,          # gated arrivals needed to relocate an event (fewer: event kept at its prior position)
    "n_joint": 15,         # joint Gauss-Newton iterations after the alternation of each pass (0: alternation only)
    "cauchy": 2.385,       # Cauchy weight constant of the joint stage (95% Gaussian efficiency)
    "term": "median",      # station-term update: median | huber
    # S3
    "n_unl": 30,           # unlabelled windows per network
    "q_det": 0.80,         # atoms with selection score above this quantile of the stream's atom scores
    "sx_det": 3.0,         # km, prior sd of a detected atom's position
    "st_det": 0.3,         # s, prior sd of a detected atom's origin time
    "min_arr_det": 6,      # gated arrivals needed for a detected atom to enter
    "w_det0": 0.5,         # s, half-width of the first read of an atom (around its own prediction)
    "det_consist": 0.0,    # atoms with median |residual| > det_consist x scale after a pass are dropped (0: off)
    "n_pass_s3": 2,        # re-reads of the joint (catalogue + atoms) estimation
    "n_iter_s3": 1,
}


# ============================================================================================ reading peaks
def energy(D):
    MF, sd = xc.matched(D)
    return np.sum((MF / sd) ** 2, axis=1)                                          # (ns, NT)


def tt_hom(STA, x):
    d = np.linalg.norm(STA - np.asarray(x, float)[None], axis=1)
    return d[:, None] / VEL[None, :]                                               # (ns, 2)


def read_peaks(E, centre, half):
    """station_corrections.peak_delays (joint P/S read, P before S by MIN_PS samples, parabolic refinement) around
    arbitrary centres (ns, 2) [s] with half-width half [s].  Returns absolute peak times (ns, 2) and energies."""
    W = int(round(half / DT))
    ns = E.shape[0]
    tobs = np.zeros((ns, 2))
    en = np.zeros((ns, 2))
    for i in range(ns):
        win = []
        for b in range(2):
            c = int(round(centre[i, b] / DT))
            c = min(max(c, 1), NT - 2)
            win.append(np.arange(max(1, c - W), min(NT - 2, c + W) + 1))
        sc = E[i, win[0]][:, None] + E[i, win[1]][None, :]
        sc = np.where(win[1][None, :] >= win[0][:, None] + MIN_PS, sc, -np.inf)
        if not np.isfinite(sc).any():
            ks = (int(win[0][np.argmax(E[i, win[0]])]), int(win[1][np.argmax(E[i, win[1]])]))
        else:
            a, b_ = np.unravel_index(int(np.argmax(sc)), sc.shape)
            ks = (int(win[0][a]), int(win[1][b_]))
        for b in range(2):
            k = ks[b]
            e0, e1, e2 = E[i, k - 1], E[i, k], E[i, k + 1]
            den = e0 - 2 * e1 + e2
            dk = 0.5 * (e0 - e2) / den if den < 0 else 0.0
            tobs[i, b] = (k + float(np.clip(dk, -0.5, 0.5))) * DT
            en[i, b] = e1
    return tobs, en


# ============================================================================================ JHD
class Ev:
    """an event of the joint estimation: energy traces, prior hypocentre and its prior sds, current hypocentre."""

    def __init__(self, E, x0, t0, sx, st, min_arr, kind):
        self.E = E
        self.x0, self.t0 = np.asarray(x0, float).copy(), float(t0)
        self.x, self.t = self.x0.copy(), self.t0
        self.sx, self.st, self.min_arr, self.kind = sx, st, min_arr, kind
        self.tobs = self.ok = None


def relocate(ev, STA, c, sig, par):
    """MAP hypocentre of an event: Huber misfit of the gated arrivals (scale sig) + Gaussian prior around
    (x0, t0).  Iteratively reweighted Gauss-Newton from the current hypocentre."""
    ok = ev.ok
    if ok.sum() < ev.min_arr:
        ev.x, ev.t = ev.x0.copy(), ev.t0
        return
    x, t = ev.x.copy(), ev.t
    tb, cb = ev.tobs[ok], c[ok]
    ib = np.nonzero(ok)                                                            # station index, phase index
    Pinv = np.diag([1 / ev.sx ** 2] * 3 + [1 / ev.st ** 2])
    for _ in range(12):
        dv = STA - x[None]
        d = np.linalg.norm(dv, axis=1)
        n = dv / d[:, None]
        r = tb - d[ib[0]] / VEL[ib[1]] - t - cb
        J = np.concatenate([-n[ib[0]] / VEL[ib[1]][:, None], np.ones((len(r), 1))], axis=1)   # d pred / d (x, t)
        u = np.abs(r) / sig
        w = np.where(u <= par["huber"], 1.0, par["huber"] / np.maximum(u, 1e-12)) / sig ** 2
        dp0 = np.r_[x - ev.x0, t - ev.t0]
        A = J.T @ (J * w[:, None]) + Pinv
        g = J.T @ (w * r) - Pinv @ dp0
        step = np.linalg.solve(A, g)
        x = x + step[:3]
        x[2] = float(np.clip(x[2], 1.0, 20.0))
        t = t + float(step[3])
        if np.abs(step).max() < 1e-4:
            break
    ev.x, ev.t = x, t


def _loc(v, sig, par):
    v = np.asarray(v, float)
    m = float(np.median(v))
    if par["term"] == "median" or len(v) < 3:
        return m
    for _ in range(20):
        u = np.abs(v - m) / sig
        w = np.where(u <= par["huber"], 1.0, par["huber"] / np.maximum(u, 1e-12))
        m2 = float((w * v).sum() / w.sum())
        if abs(m2 - m) < 1e-5:
            m = m2
            break
        m = m2
    return m


def jhd(STA, events, c0, par=P, first_read=True, info=None):
    """events: list of Ev.  c0: starting terms (ns, 2).  Returns the terms; events hold the final hypocentres."""
    ns = len(STA)
    c = np.asarray(c0, float).copy()
    sig = 0.15
    hist = []
    for ps in range(par["n_pass"]):
        for ev in events:
            if ps == 0 and first_read and ev.kind == "cat":
                ev.tobs, en = read_peaks(ev.E, tt_hom(STA, ev.x0) + ev.t0, par["search0"])
            else:
                wr = par["w_re"]
                wr = wr[min(ps - (1 if first_read else 0), len(wr) - 1)] if isinstance(wr, (list, tuple)) else wr
                if ev.kind == "det" and ev.tobs is None:
                    wr = par["w_det0"]
                ev.tobs, en = read_peaks(ev.E, tt_hom(STA, ev.x) + ev.t + c, wr)
            ev.ok = en >= (GATE if (ps == 0 and first_read and ev.kind == "cat") else par["gate_re"])
        for _ in range(par["n_alt"]):
            # residual scale from the current state
            res = np.concatenate([(ev.tobs - tt_hom(STA, ev.x) - ev.t - c)[ev.ok] for ev in events] + [np.zeros(0)])
            if len(res) >= 10:
                sig = max(par["sig_floor"], 1.4826 * float(np.median(np.abs(res - np.median(res)))))
            for ev in events:
                relocate(ev, STA, c, sig, par)
            # gauge: median origin shift of the catalogue events = 0 (the catalogue defines absolute time)
            dts = [ev.t - ev.t0 for ev in events if ev.kind == "cat" and ev.ok.sum() >= ev.min_arr]
            if dts:
                g = float(np.median(dts))
                for ev in events:
                    if ev.ok.sum() >= ev.min_arr:
                        ev.t -= g
            new = c.copy()
            for i in range(ns):
                for b in range(2):
                    v = [ev.tobs[i, b] - np.linalg.norm(STA[i] - ev.x) / VEL[b] - ev.t for ev in events
                         if ev.ok[i, b] and (ev.kind == "cat" or ev.ok.sum() >= ev.min_arr)]
                    if v:
                        new[i, b] = _loc(v, sig, par)
            dc = float(np.abs(new - c).max())
            c = new
            hist.append(dc)
            if dc < par["tol"]:
                break
        if par["n_joint"]:
            c, sig = joint_solve(STA, events, c, sig, par)
        if par["det_consist"]:
            # self-consistency of the atoms: an atom whose gated arrivals do not fit a hypocentre leaves the estimate
            for ev in events:
                if ev.kind == "det" and ev.ok.sum() >= ev.min_arr:
                    r = (ev.tobs - tt_hom(STA, ev.x) - ev.t - c)[ev.ok]
                    if np.median(np.abs(r)) > par["det_consist"] * sig:
                        ev.min_arr = 99
    if info is not None:
        res = np.concatenate([(ev.tobs - tt_hom(STA, ev.x) - ev.t - c)[ev.ok] for ev in events] + [np.zeros(0)])
        info.update({"sig": sig, "hist": hist, "n_res": int(len(res)),
                     "n_ev": int(sum(ev.ok.sum() >= ev.min_arr for ev in events)),
                     "res_mad": float(1.4826 * np.median(np.abs(res))) if len(res) else np.nan})
    return c


def joint_solve(STA, events, c, sig, par):
    """joint iteratively reweighted Gauss-Newton over the terms and the hypocentres of the relocatable events
    (the alternation converges slowly along the directions that couple terms and hypocentres).  Cauchy weights with the
    residual scale re-estimated at every iteration; Gaussian priors of the events; gauge fixed afterwards."""
    ns = len(STA)
    act = [ev for ev in events if ev.ok.sum() >= ev.min_arr]
    fix = [ev for ev in events if ev.ok.sum() < ev.min_arr and ev.kind == "cat" and ev.ok.any()]
    for ev in fix:
        ev.x, ev.t = ev.x0.copy(), ev.t0
    npar = 2 * ns + 4 * len(act)
    c = c.copy()
    for _ in range(par["n_joint"]):
        rr, cols, vals, rows = [], [], [], []
        nrow = 0
        for k, ev in enumerate([None] * 0 + act + fix):
            i_, b_ = np.nonzero(ev.ok)
            dv = STA - ev.x[None]
            d = np.linalg.norm(dv, axis=1)
            n = dv / d[:, None]
            r = ev.tobs[i_, b_] - d[i_] / VEL[b_] - ev.t - c[i_, b_]
            m = len(r)
            ridx = nrow + np.arange(m)
            rows.append(ridx)
            cols.append(2 * i_ + b_)
            vals.append(np.ones(m))
            if k < len(act):
                Jx = -n[i_] / VEL[b_][:, None]
                for a in range(3):
                    rows.append(ridx)
                    cols.append(np.full(m, 2 * ns + 4 * k + a))
                    vals.append(Jx[:, a])
                rows.append(ridx)
                cols.append(np.full(m, 2 * ns + 4 * k + 3))
                vals.append(np.ones(m))
            rr.append(r)
            nrow += m
        r = np.concatenate(rr)
        if len(r) < 10:
            break
        sig = max(par["sig_floor"], 1.4826 * float(np.median(np.abs(r - np.median(r)))))
        J = np.zeros((nrow, npar))
        J[np.concatenate(rows), np.concatenate(cols)] = np.concatenate(vals)
        w = 1.0 / (1.0 + (r / (par["cauchy"] * sig)) ** 2) / sig ** 2
        A = J.T @ (J * w[:, None])
        g = J.T @ (w * r)
        A[np.arange(2 * ns), np.arange(2 * ns)] += 1e-6
        for k, ev in enumerate(act):
            sl = slice(2 * ns + 4 * k, 2 * ns + 4 * k + 4)
            pv = np.array([1 / ev.sx ** 2] * 3 + [1 / ev.st ** 2])
            A[sl, sl] += np.diag(pv)
            g[sl] -= pv * np.r_[ev.x - ev.x0, ev.t - ev.t0]
        step = np.linalg.solve(A, g)
        c = c + step[:2 * ns].reshape(ns, 2)
        for k, ev in enumerate(act):
            sp = step[2 * ns + 4 * k:2 * ns + 4 * k + 4]
            ev.x = ev.x + sp[:3]
            ev.x[2] = float(np.clip(ev.x[2], 1.0, 20.0))
            ev.t = ev.t + float(sp[3])
        if np.abs(step[:2 * ns]).max() < 1e-4:
            break
    dts = [ev.t - ev.t0 for ev in act if ev.kind == "cat"]
    if dts:
        gsh = float(np.median(dts))
        for ev in act:
            ev.t -= gsh
        has = np.zeros((ns, 2), bool)
        for ev in act + fix:
            has |= ev.ok
        c = np.where(has, c + gsh, c)
    return c, sig


def cat_events(cat, Ecat=None, par=P):
    Ecat = Ecat if Ecat is not None else [energy(D) for D, _s, _t in cat]
    return [Ev(E, src, tau, par["sx"], par["st"], par["min_arr"], "cat") for E, (_D, src, tau) in zip(Ecat, cat)]


def s2_terms(STA, cat, c0, Ecat=None, par=P, info=None):
    return jhd(STA, cat_events(cat, Ecat, par), c0, par, info=info)


# ============================================================================================ S3
def select_atoms(stream_atoms, par=P):
    """stream_atoms: list over windows of lists of atoms (dicts x, tau, s).  Returns [(window index, atom)] with
    s >= the q_det quantile of all atom scores of the stream."""
    sc = np.array([a["s"] for at in stream_atoms for a in at])
    if not len(sc):
        return []
    thr = float(np.quantile(sc, par["q_det"]))
    return [(wi, a) for wi, at in enumerate(stream_atoms) for a in at if a["s"] >= thr]


def s3_terms(STA, cat, c_s0, c_s2, Ewin, stream_atoms, Ecat=None, par=P, info=None):
    """joint JHD of the catalogue events and the selected atoms of the unlabelled stream.  Ewin: energies per window.
    The catalogue events first go through S2 (from the S0 terms; this reproduces c_s2, the terms the atoms were detected
    with, and leaves the events relocated); the atoms are then added, read around their own prediction with the S2 terms,
    and the joint estimation continues with n_pass_s3 re-reads of all events."""
    evs = cat_events(cat, Ecat, par)
    c2 = jhd(STA, evs, c_s0, par)
    assert np.allclose(c2, c_s2)
    sel = select_atoms(stream_atoms[:par["n_unl"]], par)
    det = [Ev(Ewin[wi], a["x"], a["tau"], par["sx_det"], par["st_det"], par["min_arr_det"], "det") for wi, a in sel]
    p3 = dict(par)
    p3["n_pass"] = par["n_pass_s3"]
    c = jhd(STA, evs + det, c2, p3, first_read=False, info=info)
    if info is not None:
        info["n_sel"] = len(sel)
        info["n_det_used"] = int(sum(ev.ok.sum() >= ev.min_arr for ev in det))
    return c


# ============================================================================================ quality
def err(c, statics):
    """C10 metric: rms of (c - statics) with its median removed (the common mode is absorbed by the origin time)."""
    d = c - statics
    return float(np.sqrt(((d - np.median(d)) ** 2).mean()))


def offset(c, statics):
    return float(np.median(c - statics))


# ============================================================================================ engines
def clone_eng(E, corr):
    """engine with the same nodes, operators and Jacobians and other corrections (only the sample tables change)."""
    E2 = copy.copy(E)
    d = np.linalg.norm(E.XN[:, None, :] - E.STA[None, :, :], axis=2)
    tt = np.stack([d / VP, d / VS], axis=1)
    x = (tt[:, None] + E.taus[None, :, None, None] + np.asarray(corr, float).T[None, None]) / DT
    E2.IDX = np.clip(x.astype(int), 0, NT - 1)
    E2.frac = x - E2.IDX
    E2.cache = {}
    return E2
