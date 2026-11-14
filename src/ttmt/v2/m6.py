"""M6 (first version): joint refit of all accepted atoms, and the "one against two" resolvability test.

Why. A greedy step subtracts one atom from a residual that still holds the other events. When two events overlap,
the arrivals of the second fall under the response of the first and every greedy deflation leaves 6-12 % behind or
eats the neighbour (M6C-SUBSAMPLE, PR #39). The stage-loss diagnosis (PR #38) puts the loss of PyOcto and GaMMA in
the 0.1-1 s regime in merged pairs. M6 explains the window with all accepted atoms at once.

Model. Atom k: node x_k (kept at its lattice node), moment tensor m_k, and a continuous peak time T_kib for every
block (station i, phase b), bounded around its anchor (the refined read time, or the prediction + lag for a block
with no onset). The model of the normalised matched filter is

    Y_i(t) ~ sum_k sum_b (G_k m_k)_{ib} r(t - T_kib),    r = analytic deflation response of the kernel.

For fixed times the model is linear in (m_1 ... m_K): least squares in closed form. Times are refined block by
block, with the others fixed (bounded scalar search); a few alternations. This is the "reoptimise everything" step
of ADCG [Boyd et al.] restricted to times and tensors. Before each block's bounded search, the block may be
re-assigned to any onset maximum of q within the reader's window (the M4 choice, now scored by the full model's
residual, so a block that first locked on the other event of a pair can move to its own arrival) (positions stay on nodes: within the timing tolerance a node
offset is absorbed by the block times, which is the timing-tolerant class of Appendix AB).

Resolvability (replaces the 5 km exclusion of the registered loop). When a new atom falls within `near_km` and
`near_s` of an accepted one, the two-event fit is compared with the best ONE-event fit whose arrivals are free (one P
and one S pulse per station, free times and free 3-vectors fitted jointly, the class of Theorem A.3(a)), with the
other atoms refitted around it (``free_one_event_joint``). The pair is accepted only if the residual sum of squares
drops by more than `thr`. `thr` has no theoretical value: the matched-filter samples are correlated, so a difference
of RSS is not on the scale of the noise-free margin delta^2 of Corollary A.1 (which is a necessary condition on the
signal, not the level of a statistic). The level is set only by a calibration control with neighbouring events
(review T5.2 of 9 November 2026, `teoria-fase5/revision-A1-A3a-M6.md`).

Not yet: EGF / coda dictionary (Ricker + MT only), pruning by S_mix, the a-posteriori gap of AA.4(d).
"""
import numpy as np
from scipy.optimize import minimize_scalar

from ..kernel import LEGACY as LEGACY_KERNEL
from .loop import RANK_DEFAULT, default_reader, refine_times
from .m1 import read_vectors
from .m4 import candidates


class _Resp:
    """cached analytic responses on the support -H..H for a float peak time."""

    def __init__(self, kernel):
        self.k, self.H, self.cache = kernel, kernel.ac_half, {}

    def __call__(self, t):
        k0 = int(np.floor(t + 0.5))
        d = round(float(t - k0), 5)
        if d not in self.cache:
            self.cache[d] = self.k.response(d)
        return k0, self.cache[d]


def _place(out_i, a, t, resp, sign=1.0):
    """add sign * a (3,) x r(. - t) into out_i (3, NT) in place."""
    nt, H = out_i.shape[-1], resp.H
    k0, r = resp(t)
    lo, hi = max(0, k0 - H), min(nt, k0 + H + 1)
    if hi > lo:
        out_i[:, lo:hi] += sign * np.outer(a, r[H - (k0 - lo):H + (hi - k0)])


def model(atoms, G, shape, resp):
    """model traces (ns, 3, NT) of the atoms (dicts gi, T (2, ns), m (6,))."""
    M = np.zeros(shape)
    ns = shape[0]
    for a in atoms:
        g = (G[a["gi"]] @ a["m"]).reshape(2, ns, 3)
        for b in range(2):
            for i in range(ns):
                _place(M[i], g[b, i], a["T"][b, i], resp)
    return M


def fit_tensors(Y, atoms, G, resp):
    """least-squares tensors of all atoms at their current times (closed form). In place; returns the atoms."""
    ns, _, nt = Y.shape
    K = len(atoms)
    A = np.zeros((ns, 3, nt, 6 * K))
    for k, a in enumerate(atoms):
        Gk = G[a["gi"]].reshape(2, ns, 3, 6)
        for b in range(2):
            for i in range(ns):
                k0, r = resp(a["T"][b, i])
                H = resp.H
                lo, hi = max(0, k0 - H), min(nt, k0 + H + 1)
                if hi <= lo:
                    continue
                seg = r[H - (k0 - lo):H + (hi - k0)]
                A[i, :, lo:hi, 6 * k:6 * k + 6] += Gk[b, i][:, None, :] * seg[None, :, None]
    A = A.reshape(-1, 6 * K)
    rows = np.any(A != 0, axis=1)
    m, *_ = np.linalg.lstsq(A[rows], Y.reshape(-1)[rows], rcond=None)
    for k, a in enumerate(atoms):
        a["m"] = m[6 * k:6 * k + 6]
    return atoms


def refine_block_times(Y, atoms, G, resp, sweeps=1):
    """bounded scalar refinement of every block time of every atom, the others fixed. In place."""
    ns = Y.shape[0]
    M = model(atoms, G, Y.shape, resp)
    for _ in range(sweeps):
        for a in atoms:
            g = (G[a["gi"]] @ a["m"]).reshape(2, ns, 3)
            for b in range(2):
                for i in range(ns):
                    if not np.any(g[b, i]):
                        continue
                    base = Y[i] - M[i]
                    _place(base, g[b, i], a["T"][b, i], resp, +1.0)     # residual without this block
                    H = resp.H
                    cands = list(a.get("cands", [[[]] * ns] * 2)[b][i]) + [a["anchor"][b, i]]
                    span = [min(cands), max(cands)]
                    bd_max = float(np.max(a["bound"]))
                    lo = max(0, int(np.floor(span[0])) - H - int(bd_max) - 1)
                    hi = min(Y.shape[-1], int(np.ceil(span[1])) + H + int(bd_max) + 2)
                    seg = base[:, lo:hi]

                    def cost(t):
                        tmp = np.zeros((3, Y.shape[-1]))
                        _place(tmp, g[b, i], t, resp)
                        d = seg - tmp[:, lo:hi]
                        return float((d * d).sum())
                    if len(cands) > 1:                                  # re-assignment among the onset maxima
                        a["anchor"][b, i] = min(cands, key=cost)
                    t0, bd = a["anchor"][b, i], a["bound"][b, i]
                    t_new = minimize_scalar(cost, bounds=(t0 - bd, t0 + bd), method="bounded",
                                            options={"xatol": 1e-3}).x
                    _place(M[i], g[b, i], a["T"][b, i], resp, -1.0)
                    a["T"][b, i] = float(t_new)
                    _place(M[i], g[b, i], a["T"][b, i], resp, +1.0)
    return atoms


def joint_refit(Y, atoms, G, resp, rounds=3):
    """alternate closed-form tensors and block-time refinement; returns (atoms, residual, rss)."""
    for _ in range(rounds):
        fit_tensors(Y, atoms, G, resp)
        refine_block_times(Y, atoms, G, resp)
    fit_tensors(Y, atoms, G, resp)
    R = Y - model(atoms, G, Y.shape, resp)
    return atoms, R, float((R * R).sum())


def free_one_event_rss(Y, others, a, G, resp, free_samples=8):
    """residual sum of squares when atom a's blocks are explained by ONE arrival each with free time
    (within +-free_samples of its current time) and free 3-vector, the other atoms fixed (Theorem A.3(a))."""
    ns, _, nt = Y.shape
    R = Y - model(others, G, Y.shape, resp)
    H = resp.H
    for b in range(2):
        for i in range(ns):
            t0 = a["T"][b, i]
            k0 = int(np.floor(t0 + 0.5))
            lo, hi = max(0, k0 - H - free_samples - 1), min(nt, k0 + H + free_samples + 2)
            seg = R[i, :, lo:hi].copy()

            def fit(t):
                tmp = np.zeros((3, nt))
                _place(tmp, np.ones(3), t, resp)
                r = tmp[0, lo:hi]
                rr = float(r @ r)
                if rr <= 0:
                    return seg, float((seg * seg).sum())
                v = seg @ r / rr
                d = seg - np.outer(v, r)
                return d, float((d * d).sum())
            # the cost oscillates with the response lobes: grid over whole samples, then a local bounded search
            grid = t0 + np.arange(-free_samples, free_samples + 1)
            tg = grid[int(np.argmin([fit(t)[1] for t in grid]))]
            t_best = minimize_scalar(lambda t: fit(t)[1], bounds=(tg - 0.6, tg + 0.6),
                                     method="bounded", options={"xatol": 1e-3}).x
            R[i, :, lo:hi] = fit(t_best)[0]
    return float((R * R).sum())


def _resp_full(t, nt, resp):
    tmp = np.zeros((1, nt))
    _place(tmp, np.ones(1), t, resp)
    return tmp[0]


def _free_station(Ri, tP, tS, resp, half):
    """best joint fit of ONE P and ONE S pulse with free times and free 3-vectors to Ri (3, NT) (Theorem A.3(a): the two
    pulses solved together with their 2 x 2 Gram, not one after the other). Grid of whole samples within +-half of
    (tP, tS), then a local search. Returns (residual (3, NT), rss)."""
    nt = Ri.shape[-1]
    off = np.arange(-half, half + 1)
    RP = np.array([_resp_full(tP + o, nt, resp) for o in off])
    RS = np.array([_resp_full(tS + o, nt, resp) for o in off])
    pp, ss = (RP * RP).sum(1), (RS * RS).sum(1)
    ps = RP @ RS.T
    bP, bS = Ri @ RP.T, Ri @ RS.T                                    # (3, n)
    nP, nS = (bP * bP).sum(0), (bS * bS).sum(0)
    cr = bP.T @ bS                                                   # <bP_k, bS_l>
    det = pp[:, None] * ss[None, :] - ps ** 2
    ok = det > 1e-9 * pp[:, None] * ss[None, :]
    expl = np.where(ok, (ss[None, :] * nP[:, None] - 2 * ps * cr + pp[:, None] * nS[None, :]) / np.where(ok, det, 1),
                    np.maximum(nP[:, None] / np.maximum(pp[:, None], 1e-300), nS[None, :] / np.maximum(ss[None, :], 1e-300)))
    k, l = np.unravel_index(int(np.argmax(expl)), expl.shape)

    def fit(th):
        rp, rs = _resp_full(th[0], nt, resp), _resp_full(th[1], nt, resp)
        A = np.stack([rp, rs], 1)                                    # (NT, 2)
        V, *_ = np.linalg.lstsq(A, Ri.T, rcond=None)                 # (2, 3)
        d = Ri - (A @ V).T
        return d, float((d * d).sum())
    from scipy.optimize import minimize
    th0 = np.array([tP + off[k], tS + off[l]], float)
    res = minimize(lambda th: fit(th)[1], th0, method="Nelder-Mead",
                   options={"xatol": 1e-3, "fatol": 1e-9, "initial_simplex": th0 + np.array([[0, 0], [0.5, 0], [0, 0.5]])})
    th = res.x if res.fun <= fit(th0)[1] else th0
    return fit(th)


def free_one_event_joint(Y, others, a, G, resp, rounds=3, half=16):
    """residual sum of squares of the ONE-event alternative in the same class as the two-event fit (review T5.2 of
    9 November 2026, items b1-b3): atom a's blocks are explained per station by one P and one S pulse with free times
    (+-half samples, 1 s at 16 Hz) and free 3-vectors fitted jointly (A.3(a)), and the OTHER atoms are refitted
    (tensors and block times, ``joint_refit``) around it, alternating ``rounds`` times."""
    ns, _, nt = Y.shape
    oth = [dict(o, T=o["T"].copy(), anchor=np.array(o.get("anchor", o["T"]), float).copy(),
                bound=o.get("bound", np.full(o["T"].shape, 2.0))) for o in others]
    tc = a["T"].copy()
    F = np.zeros_like(Y)
    for _ in range(rounds):
        R = Y - (model(oth, G, Y.shape, resp) if oth else 0.0)
        for i in range(ns):
            d, _r = _free_station(R[i], tc[0, i], tc[1, i], resp, half)
            F[i] = R[i] - d
        if not oth:
            break
        oth, _R, _rss = joint_refit(Y - F, oth, G, resp, rounds=1)
    R = Y - F - (model(oth, G, Y.shape, resp) if oth else 0.0)
    return float((R * R).sum())


def onset_candidates(q_b, tpred_b, half):
    """per station, the sub-sample times (samples) of every onset maximum of q_b (ns, NT) within tpred_b +- half."""
    out = []
    for i in range(q_b.shape[0]):
        pos, _qv = candidates(q_b[i], tpred_b[i], half / 3.0, nsig=3.0)
        out.append(pos)
    return out


def select(E, R, q, lag, statistic, reader, excluded):
    """best (node, origin) on the residual with onset reading (as v2_loop without M4)."""
    ns = E.ns
    lagT = np.asarray(lag, float).T
    s_on = reader.read(q, E.tpred)
    vf = read_vectors(R, s_on + lagT[None, None]).reshape(E.N, E.T, 6 * ns)
    sc, _aux = statistic.score(vf, E.Q)
    sc[excluded] = -np.inf
    if not np.isfinite(sc).any():
        return None
    gi = int(np.argmax(sc.max(axis=1)))
    tj = int(np.argmax(sc[gi]))
    return gi, tj, float(sc[gi, tj]), s_on[gi, tj]


def adcg_loop(E, MFn, q, kmax, lag, statistic=RANK_DEFAULT, reader=None, thr=None, near_km=5.0, near_s=1.5,
              bound_read=2.0, bound_unread=4.0, rounds=3, cand_half=8, kernel=LEGACY_KERNEL, return_residual=False,
              free="joint"):
    """M6 loop: add the best atom on the residual, refit all atoms jointly on MFn, and keep a new atom near an
    accepted one only if the resolvability test passes (thr; None = always keep, i.e. no test and no exclusion).
    A rejected candidate excludes the nodes within near_km of it, so the loop moves on.
    Atoms: dicts x, tau, s (ranking score at its step), gi, T (2, ns) peak times, m (6,), pair_gain (RSS drop of the
    resolvability test, None if not tested). free: "joint" (default since 9 November 2026, review T5.2) computes
    the one-event alternative with the other atoms refitted and P and S fitted jointly; "legacy" is the first version
    (``free_one_event_rss``), kept for the bias measurement."""
    reader = reader or default_reader(kernel.dt)
    resp = _Resp(kernel)
    lagT = np.asarray(lag, float).T
    Y = np.asarray(MFn, float)
    R = Y.copy()
    excluded = np.zeros(E.N, bool)
    atoms, rss = [], float((Y * Y).sum())
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
        new = {"gi": gi, "x": np.asarray(E.XN[gi], float), "tau": float(E.taus[tj]), "s": s_val,
               "anchor": t_ref.copy(), "bound": bound, "T": t_ref.copy(), "m": np.zeros(6), "pair_gain": None,
               "cands": cands}
        near = [a for a in atoms if np.linalg.norm(a["x"][:2] - new["x"][:2]) <= near_km
                and abs(a["tau"] - new["tau"]) <= near_s]
        trial = [dict(a, T=a["T"].copy(), anchor=a["anchor"].copy()) for a in atoms] + [new]
        # the candidate times of every atom are re-read on the residual-free q (q is never deflated)
        trial, R_t, rss_t = joint_refit(Y, trial, E.G, resp, rounds)
        if near and thr is not None:
            j = atoms.index(min(near, key=lambda a: np.linalg.norm(a["x"][:2] - new["x"][:2])))
            others = [a for k, a in enumerate(atoms) if k != j]
            if free == "joint":
                rss_one = free_one_event_joint(Y, others, atoms[j], E.G, resp)
            else:                                      # "legacy": the first version (bounded, sequential, others fixed)
                rss_one = free_one_event_rss(Y, others, atoms[j], E.G, resp)
            gain = rss_one - rss_t
            new["pair_gain"] = gain
            if gain <= thr:
                excluded |= np.linalg.norm(E.XN[:, :2] - new["x"][None, :2], axis=1) <= near_km
                continue
        atoms, R, rss = trial, R_t, rss_t
        atoms[-1]["pair_gain"] = new["pair_gain"]
        excluded[gi] = True
    for a in atoms:
        a.pop("anchor", None)
        a.pop("cands", None)
    return (atoms, R) if return_residual else atoms


def hits_loctv(atoms, sources, tol_km=None, tol_s=None):
    """indices of the true sources hit under LOC+T-V (ADR-0008 item 2, ``ttmt.stats.label_voronoi``), atoms scored in
    acceptance order (the first accepted atom has the highest score). sources: list of (xyz, tau, ...). The one
    function every M6 experiment scores with (CR-29)."""
    from ..stats import TOL_KM, TOL_S, label_voronoi
    if not atoms or not sources:
        return set()
    sc = -np.arange(len(atoms), dtype=float)
    dxy = np.array([np.asarray(a["x"], float)[:2] for a in atoms])
    dt = np.array([float(a["tau"]) for a in atoms])
    exy = np.array([np.asarray(s[0], float)[:2] for s in sources])
    et = np.array([float(s[1]) for s in sources])
    lab = label_voronoi(sc, dxy, dt, exy, et, tol_km or TOL_KM, tol_s or TOL_S)
    return {int(v) for v in lab if v >= 0}


def duplicate_atoms(atoms, half_sep):
    """the duplicator negative control (ADR-0008 item 3, ``ttmt.stats.duplicator``) on an atom list: every atom twice,
    at tau - half_sep and tau + half_sep; a sound pair metric gives it no more pair recall than the atoms."""
    out = []
    for a in atoms:
        out.append(dict(a, tau=float(a["tau"]) - half_sep))
    for a in atoms:
        out.append(dict(a, tau=float(a["tau"]) + half_sep))
    return out
