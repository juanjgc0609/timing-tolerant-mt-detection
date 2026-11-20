"""M6 with S-wave splitting as a constraint inside the joint fit (director's request of 20 November 2026; design in
the supervision note of that day, section 3). Exploratory; ``ttmt.v2.m6`` is unchanged.

Model. The S block of atom k at station i is split along two orthogonal axes in the plane normal to the ray of that
atom (fast axis ef at angle phi_i from e1 = z x n, slow axis es):

    S_ki(t) = (Pf g_ki) r(t - T_kSi) + (Ps g_ki) r(t - T_kSi - delta_i),   Pf = ef ef' + n n',  Ps = es es',

with g_ki = (G_k m_k)_{S,i}. The P block is one pulse (splitting cannot explain a second P). phi_i and delta_i, with
0 <= delta_i <= DMAX_S, belong to the station (a property of the path): in ``mode="station"`` they are common to
every atom of the window and estimated in the same joint fit, in BOTH hypotheses of the resolvability test (the two
events, and the one event with free arrivals); delta = 0 is the model of ``ttmt.v2.m6``, so the class is nested.
``mode="event"`` gives every atom (and the free one-event alternative) its own phi and delta per station: the free
splitting of exploration X1, kept as a descriptive arm.

Estimation, per round of the joint fit: tensors in closed form (linear for fixed times and splitting), block times
by bounded scalar search, then (phi_i, delta_i) per station on a grid (PHI_STEP, DELTA_STEP) refined by bounded
scalar searches, the other parameters fixed. The free one-event alternative fits per station one free P 3-vector and
an S whose fast-axis and along-ray components arrive at t_S and slow-axis component at t_S + delta (free amplitudes),
with free times, starting from the unsplit fit of ``m6._free_station``.
"""
import numpy as np
from scipy.optimize import minimize, minimize_scalar

from ..kernel import LEGACY as LEGACY_KERNEL
from .loop import RANK_DEFAULT, default_reader, refine_times
from .m6 import _free_station, _place, _resp_full, _Resp, onset_candidates, select

DMAX_S = 0.2                                    # largest split delay (s), fixed by the design
PHI_GRID = np.arange(12) * np.pi / 12           # fast-axis angle grid (15 degrees)
N_DELTA = 9                                     # delay grid 0 .. DMAX on N_DELTA points


class Frames:
    """per (node, station) ray frame (e1, e2, n), cached."""

    def __init__(self, XN, STA):
        self.XN, self.STA, self.c = np.asarray(XN, float), np.asarray(STA, float), {}

    def __call__(self, gi):
        if gi not in self.c:
            n = self.STA - self.XN[gi][None]
            n /= np.linalg.norm(n, axis=1, keepdims=True)
            e1 = np.cross(np.array([0.0, 0.0, 1.0])[None], n)
            nr = np.linalg.norm(e1, axis=1, keepdims=True)
            e1 = np.where(nr > 1e-9, e1 / np.maximum(nr, 1e-12), np.array([1.0, 0.0, 0.0])[None])
            e2 = np.cross(n, e1)
            self.c[gi] = (e1, e2, n)
        return self.c[gi]


def _proj(fr, i, phi):
    e1, e2, n = fr[0][i], fr[1][i], fr[2][i]
    ef = np.cos(phi) * e1 + np.sin(phi) * e2
    es = -np.sin(phi) * e1 + np.cos(phi) * e2
    return np.outer(ef, ef) + np.outer(n, n), np.outer(es, es), (ef, es, n)


def _place_S(out_i, g, t, Pf, Ps, dl, resp, sign=1.0):
    if dl <= 0:
        _place(out_i, g, t, resp, sign)
        return
    _place(out_i, Pf @ g, t, resp, sign)
    _place(out_i, Ps @ g, t + dl, resp, sign)


def _block(out_i, a, b, i, g_bi, t, resp, frames, sign=1.0):
    if b == 0 or a["dl"][i] <= 0:
        _place(out_i, g_bi, t, resp, sign)
    else:
        Pf, Ps, _ = _proj(frames(a["gi"]), i, a["phi"][i])
        _place_S(out_i, g_bi, t, Pf, Ps, a["dl"][i], resp, sign)


def model(atoms, G, shape, resp, frames, stations=None):
    M = np.zeros(shape)
    ns = shape[0]
    for a in atoms:
        g = (G[a["gi"]] @ a["m"]).reshape(2, ns, 3)
        for b in range(2):
            for i in (range(ns) if stations is None else stations):
                _block(M[i], a, b, i, g[b, i], a["T"][b, i], resp, frames)
    return M


def fit_tensors(Y, atoms, G, resp, frames):
    ns, _, nt = Y.shape
    K = len(atoms)
    A = np.zeros((ns, 3, nt, 6 * K))
    H = resp.H

    def add(i, Mat, t, k):
        k0, r = resp(t)
        lo, hi = max(0, k0 - H), min(nt, k0 + H + 1)
        if hi > lo:
            A[i, :, lo:hi, 6 * k:6 * k + 6] += Mat[:, None, :] * r[H - (k0 - lo):H + (hi - k0)][None, :, None]
    for k, a in enumerate(atoms):
        Gk = G[a["gi"]].reshape(2, ns, 3, 6)
        for i in range(ns):
            add(i, Gk[0, i], a["T"][0, i], k)
            if a["dl"][i] <= 0:
                add(i, Gk[1, i], a["T"][1, i], k)
            else:
                Pf, Ps, _ = _proj(frames(a["gi"]), i, a["phi"][i])
                add(i, Pf @ Gk[1, i], a["T"][1, i], k)
                add(i, Ps @ Gk[1, i], a["T"][1, i] + a["dl"][i], k)
    A = A.reshape(-1, 6 * K)
    rows = np.any(A != 0, axis=1)
    m, *_ = np.linalg.lstsq(A[rows], Y.reshape(-1)[rows], rcond=None)
    for k, a in enumerate(atoms):
        a["m"] = m[6 * k:6 * k + 6]
    return atoms


def refine_block_times(Y, atoms, G, resp, frames):
    ns, _, nt = Y.shape
    M = model(atoms, G, Y.shape, resp, frames)
    H = resp.H
    for a in atoms:
        g = (G[a["gi"]] @ a["m"]).reshape(2, ns, 3)
        for b in range(2):
            for i in range(ns):
                if not np.any(g[b, i]):
                    continue
                base = Y[i] - M[i]
                _block(base, a, b, i, g[b, i], a["T"][b, i], resp, frames, +1.0)
                cands = list(a.get("cands", [[[]] * ns] * 2)[b][i]) + [a["anchor"][b, i]]
                ext = int(np.ceil(a["dl"][i])) if b == 1 else 0
                bd_max = float(np.max(a["bound"]))
                lo = max(0, int(np.floor(min(cands))) - H - int(bd_max) - 1)
                hi = min(nt, int(np.ceil(max(cands))) + H + int(bd_max) + ext + 2)
                seg = base[:, lo:hi]

                def cost(t):
                    tmp = np.zeros((3, nt))
                    _block(tmp, a, b, i, g[b, i], t, resp, frames)
                    d = seg - tmp[:, lo:hi]
                    return float((d * d).sum())
                if len(cands) > 1:
                    a["anchor"][b, i] = min(cands, key=cost)
                t0, bd = a["anchor"][b, i], a["bound"][b, i]
                t_new = minimize_scalar(cost, bounds=(t0 - bd, t0 + bd), method="bounded",
                                        options={"xatol": 1e-3}).x
                _block(M[i], a, b, i, g[b, i], a["T"][b, i], resp, frames, -1.0)
                a["T"][b, i] = float(t_new)
                _block(M[i], a, b, i, g[b, i], a["T"][b, i], resp, frames, +1.0)
    return atoms


def _grid_split(cost, dmax):
    """minimise cost(phi, dl) over the grid, then bounded searches in dl and phi around the best node."""
    dls = np.linspace(0.0, dmax, N_DELTA)
    best = (cost(0.0, 0.0), 0.0, 0.0)
    for ph in PHI_GRID:
        for dl in dls[1:]:
            c = cost(ph, dl)
            if c < best[0]:
                best = (c, ph, dl)
    c0, ph, dl = best
    if dl > 0:
        st = dls[1]
        r = minimize_scalar(lambda d: cost(ph, d), bounds=(max(0.0, dl - st), min(dmax, dl + st)), method="bounded",
                            options={"xatol": 1e-3})
        if r.fun < c0:
            c0, dl = r.fun, float(r.x)
        r = minimize_scalar(lambda p: cost(p, dl), bounds=(ph - np.pi / 12, ph + np.pi / 12), method="bounded",
                            options={"xatol": 1e-3})
        if r.fun < c0:
            c0, ph = r.fun, float(r.x)
    return float(ph % np.pi), float(dl), float(c0)


def estimate_split(Y, atoms, G, resp, frames, mode, dmax, extra=None):
    """(phi_i, delta_i) per station on the residual of the window: one pair shared by all atoms (station) or one per
    atom (event); tensors and times fixed. extra (ns, 3, NT): fixed model terms (subtracted)."""
    ns = Y.shape[0]
    groups = [list(range(len(atoms)))] if mode == "station" else [[k] for k in range(len(atoms))]
    gs = [(G[a["gi"]] @ a["m"]).reshape(2, ns, 3) for a in atoms]
    for i in range(ns):
        for grp in groups:
            rest = [k for k in range(len(atoms)) if k not in grp]
            base = Y[i] - (extra[i] if extra is not None else 0.0)
            tmp = np.zeros_like(base)
            for k in range(len(atoms)):
                _block(tmp, atoms[k], 0, i, gs[k][0, i], atoms[k]["T"][0, i], resp, frames)
            for k in rest:
                _block(tmp, atoms[k], 1, i, gs[k][1, i], atoms[k]["T"][1, i], resp, frames)
            base = base - tmp

            def cost(ph, dl):
                m = np.zeros_like(base)
                for k in grp:
                    a = atoms[k]
                    Pf, Ps, _ = _proj(frames(a["gi"]), i, ph)
                    _place_S(m, gs[k][1, i], a["T"][1, i], Pf, Ps, dl, resp)
                d = base - m
                return float((d * d).sum())
            ph, dl, _c = _grid_split(cost, dmax)
            for k in grp:
                atoms[k]["phi"][i], atoms[k]["dl"][i] = ph, dl
    return atoms


def joint_refit(Y, atoms, G, resp, frames, mode, dmax, rounds=3, extra=None, split=True):
    Z = Y - extra if extra is not None else Y
    for _ in range(rounds):
        fit_tensors(Z, atoms, G, resp, frames)
        refine_block_times(Z, atoms, G, resp, frames)
        if split:
            fit_tensors(Z, atoms, G, resp, frames)
            estimate_split(Z, atoms, G, resp, frames, mode, dmax)
    fit_tensors(Z, atoms, G, resp, frames)
    R = Z - model(atoms, G, Y.shape, resp, frames)
    return atoms, R, float((R * R).sum())


def _free_split_fit(Ri, tP, tS, ax, dl, resp):
    """one free P 3-vector at tP; S: fast-axis and along-ray amplitudes at tS, slow-axis amplitude at tS + dl."""
    nt = Ri.shape[-1]
    Q = np.stack(ax, 1)                                         # columns ef, es, n
    Rr = Q.T @ Ri
    rp, rs = _resp_full(tP, nt, resp), _resp_full(tS, nt, resp)
    rsd = _resp_full(tS + dl, nt, resp) if dl > 0 else rs
    F = np.zeros_like(Rr)
    for c, rS in ((0, rs), (1, rsd), (2, rs)):
        A = np.stack([rp, rS], 1)
        v, *_ = np.linalg.lstsq(A, Rr[c], rcond=None)
        F[c] = A @ v
    Fi = Q @ F
    d = Ri - Fi
    return Fi, float((d * d).sum())


def free_one_event(Y, others, a, G, resp, frames, mode, dmax, rounds=3, half=16):
    """RSS of the one-event alternative with splitting in the same class as the two-event fit: atom a's blocks are one
    free P and one split S per station (free times and amplitudes), the other atoms refitted around it; in "station"
    mode the free S shares (phi_i, delta_i) with the other atoms, estimated here."""
    ns, _, nt = Y.shape
    oth = [dict(o, T=o["T"].copy(), anchor=np.array(o.get("anchor", o["T"]), float).copy(),
                bound=o.get("bound", np.full(o["T"].shape, 2.0)), phi=o["phi"].copy(), dl=o["dl"].copy())
           for o in others]
    fr = frames(a["gi"])
    tc = a["T"].copy()
    ph, dls = a["phi"].copy(), a["dl"].copy()
    F = np.zeros_like(Y)
    for it in range(rounds):
        for i in range(ns):
            if it == 0:                                          # start from the unsplit free fit
                Ri = Y[i] - (model(oth, G, Y.shape, resp, frames, [i])[i] if oth else 0.0)
                _d, _r, tc[0, i], tc[1, i] = _free_station(Ri, tc[0, i], tc[1, i], resp, half, return_times=True)

            def resid(phv, dlv):
                if mode == "station":
                    for o in oth:
                        o["phi"][i], o["dl"][i] = phv, dlv
                Ri = Y[i] - (model(oth, G, Y.shape, resp, frames, [i])[i] if oth else 0.0)
                return Ri

            def cost(phv, dlv):
                Ri = resid(phv, dlv)
                return _free_split_fit(Ri, tc[0, i], tc[1, i], _proj(fr, i, phv)[2], dlv, resp)[1]
            p_, d_, _c = _grid_split(cost, dmax)
            ph[i], dls[i] = p_, d_
            Ri = resid(p_, d_)
            ax = _proj(fr, i, p_)[2]
            th0 = np.array([tc[0, i], tc[1, i]], float)
            res = minimize(lambda th: _free_split_fit(Ri, th[0], th[1], ax, d_, resp)[1], th0, method="Nelder-Mead",
                           options={"xatol": 1e-3, "fatol": 1e-9,
                                    "initial_simplex": th0 + np.array([[0, 0], [0.5, 0], [0, 0.5]])})
            if (res.fun < _free_split_fit(Ri, th0[0], th0[1], ax, d_, resp)[1]
                    and np.all(np.abs(res.x - a["T"][:, i]) <= half)):
                tc[0, i], tc[1, i] = res.x
            F[i] = _free_split_fit(Ri, tc[0, i], tc[1, i], ax, d_, resp)[0]
        if not oth:
            break
        oth, _R, _rss = joint_refit(Y, oth, G, resp, frames, mode, dmax, rounds=1, extra=F, split=(mode == "event"))
    R = Y - F - (model(oth, G, Y.shape, resp, frames) if oth else 0.0)
    return float((R * R).sum())


def adcg_loop(E, MFn, q, kmax, lag, statistic=RANK_DEFAULT, reader=None, thr=None, near_km=10.0, near_s=1.5,
              bound_read=2.0, bound_unread=4.0, rounds=3, cand_half=8, kernel=LEGACY_KERNEL, mode="station",
              dmax_s=DMAX_S):
    """``ttmt.v2.m6.adcg_loop`` with splitting estimated inside the joint fit and inside the one-event alternative
    (mode "station": shared per station; "event": free per atom, descriptive). Atoms carry phi, dl (samples)."""
    reader = reader or default_reader(kernel.dt)
    resp = _Resp(kernel)
    frames = Frames(E.XN, E.STA)
    dmax = dmax_s / kernel.dt
    lagT = np.asarray(lag, float).T
    Y = np.asarray(MFn, float)
    R = Y.copy()
    ns = E.ns
    excluded = np.zeros(E.N, bool)
    atoms = []
    tries = 0
    while len(atoms) < kmax and tries < 3 * kmax:
        tries += 1
        sel = select(E, R, q, lag, statistic, reader, excluded)
        if sel is None:
            break
        gi, tj, s_val, s_on = sel
        t_read = s_on + lagT
        anchor = np.where(np.isfinite(t_read), t_read, E.tpred[gi, tj] + lagT)
        t_ref, _amp = refine_times(R, anchor, kernel)
        bound = np.where(np.isfinite(t_read), bound_read, bound_unread)
        cands = [[(c + lagT[b, i]).tolist() for i, c in enumerate(onset_candidates(q[:, b], E.tpred[gi, tj, b],
                                                                              cand_half))] for b in range(2)]
        if mode == "station" and atoms:
            phi0, dl0 = atoms[0]["phi"].copy(), atoms[0]["dl"].copy()
        else:
            phi0, dl0 = np.zeros(ns), np.zeros(ns)
        new = {"gi": gi, "x": np.asarray(E.XN[gi], float), "tau": float(E.taus[tj]), "s": s_val,
               "anchor": t_ref.copy(), "bound": bound, "T": t_ref.copy(), "m": np.zeros(6), "pair_gain": None,
               "cands": cands, "phi": phi0, "dl": dl0}
        near = [a for a in atoms if np.linalg.norm(a["x"][:2] - new["x"][:2]) <= near_km
                and abs(a["tau"] - new["tau"]) <= near_s]
        trial = [dict(a, T=a["T"].copy(), anchor=a["anchor"].copy(), phi=a["phi"].copy(), dl=a["dl"].copy())
                 for a in atoms] + [new]
        trial, R_t, rss_t = joint_refit(Y, trial, E.G, resp, frames, mode, dmax, rounds)
        if near and thr is not None:
            j = atoms.index(min(near, key=lambda a: np.linalg.norm(a["x"][:2] - new["x"][:2])))
            others = [a for k, a in enumerate(atoms) if k != j]
            rss_one = free_one_event(Y, others, atoms[j], E.G, resp, frames, mode, dmax)
            gain = rss_one - rss_t
            new["pair_gain"] = gain
            if gain <= thr:
                excluded |= np.linalg.norm(E.XN[:, :2] - new["x"][None, :2], axis=1) <= near_km
                continue
        atoms, R = trial, R_t
        atoms[-1]["pair_gain"] = new["pair_gain"]
        excluded[gi] = True
    for a in atoms:
        a.pop("anchor", None)
        a.pop("cands", None)
    return atoms
