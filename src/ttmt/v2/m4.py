"""M4: joint reading of one hypothesis (specification v2, M4).

Why. Read block by block, a station can lock on the wrong local maximum (a side lobe, the other event of an
overlapping pair). Appendix AB, Lemma B.2: a block read at the wrong peak does not just lose energy, it SUBTRACTS,
because the vector read there has the wrong sign relative to the moment-tensor subspace. M4 chooses the reading of
every block of a hypothesis theta together, by consistency with that subspace:

    s* = argmax_s  score(v(s))  +  sum_read [ -(s_ib - mu_ib)^2 / (2 sigma_ib^2) + log q~_ib(s_ib) ]  +  n_unread log eps

over s in prod_ib C_ib, C_ib = {local maxima of q_ib within mu_ib +- 3 sigma_ib} U {unread}, subject to C-III in its
per-station form (s_S >= s_P + min_ps when both phases of a station are read). The inter-station Lipschitz part of
C-III is not imposed (inside +-3 sigma windows it rarely binds; declared). v(s) is read at s + onset-to-peak lag
with sub-sample interpolation (M1). Scale: the score is the statistic of the normalised matched filter (unit noise
sd), so it is a log-likelihood ratio and the timing and onset terms are added in the same units.

Solution: exact enumeration when prod |C_ib| <= max_enum, else coordinate ascent from three starts (per-block best
local term, nearest to mu, largest q), at most `rounds` sweeps each.
"""
import itertools

import numpy as np

from .m1 import FLOOR, local_maxima, log_likelihood_ratio, parabolic_offset, read_vectors

MIN_PS = 6


def candidates(q_b, mu, sig, nsig=3.0, floor=FLOOR):
    """sub-sample positions and q values of the local maxima of one block's q within mu +- nsig sig."""
    nt = len(q_b)
    lo, hi = max(1, int(np.floor(mu - nsig * sig))), min(nt - 2, int(np.ceil(mu + nsig * sig)))
    if hi < lo:
        return np.zeros(0), np.zeros(0)
    mx = local_maxima(q_b[None], floor)[0]
    k = np.arange(lo, hi + 1)
    k = k[mx[k]]
    d = parabolic_offset(q_b[k - 1], q_b[k], q_b[k + 1])
    return k + d, q_b[k]


def joint_read(q, R, mu, sig, lag, Q, statistic, eps=0.05, calib=None, nsig=3.0, max_enum=4096, rounds=5,
               min_ps=MIN_PS):
    """joint reading of one hypothesis.

    q (ns, 2, NT) onset evidence; R (ns, 3, NT) residual matched filter; mu, sig (2, ns) predicted onset and its
    tolerance (samples); lag (ns, 2) onset-to-peak lag (samples); Q (6 ns, k) basis of the hypothesis; statistic: a
    detector BlockStatistic (ACAP by default in the v2 loop); calib: an IsotonicCalibration or None (raw q).
    Returns (s_star (2, ns) float samples, NaN = unread; v (6 ns,); objective; diagnostics dict)."""
    ns = R.shape[0]
    blocks = [(b, i) for b in range(2) for i in range(ns)]
    C, local, vecs = [], [], []
    for b, i in blocks:
        pos, qv = candidates(q[i, b], mu[b, i], sig[b, i], nsig)
        qc = calib.transform(qv) if calib is not None else qv
        loc = -(pos - mu[b, i]) ** 2 / (2 * sig[b, i] ** 2) + log_likelihood_ratio(qc)
        t = np.full((2, ns), np.nan)
        vv = []
        for p in pos:
            t[b, i] = p + lag[i, b]
            vv.append(read_vectors(R, t)[b, i])
        C.append(np.r_[pos, np.nan])
        local.append(np.r_[loc, np.log(eps)])
        vecs.append(np.array(vv + [np.zeros(3)]))
    sizes = [len(c) for c in C]

    def assemble(sel):
        v = np.zeros((2, ns, 3))
        for j, (b, i) in enumerate(blocks):
            v[b, i] = vecs[j][sel[j]]
        return v

    def feasible(sel):
        for i in range(ns):
            sp, ss = C[i][sel[i]], C[ns + i][sel[ns + i]]
            if np.isfinite(sp) and np.isfinite(ss) and ss < sp + min_ps:
                return False
        return True

    def objective(sel):
        if not feasible(sel):
            return -np.inf
        v = assemble(sel).reshape(1, 1, -1)
        sc, _ = statistic.score(v, Q[None])
        return float(sc[0, 0]) + float(sum(local[j][sel[j]] for j in range(len(blocks))))

    n_eval = 0
    if np.prod(sizes, dtype=float) <= max_enum:
        best, bsel = -np.inf, None
        for sel in itertools.product(*[range(n) for n in sizes]):
            o = objective(sel)
            n_eval += 1
            if o > best:
                best, bsel = o, sel
        method = "exact"
    else:
        starts = [
            [int(np.argmax(l)) for l in local],
            [int(np.nanargmin(np.abs(c - mu[b, i]))) if len(c) > 1 else 0 for c, (b, i) in zip(C, blocks)],
            [int(np.argmax(np.r_[q[i, b][np.clip(np.round(c[:-1]).astype(int), 0, q.shape[-1] - 1)], -1]))
             for c, (b, i) in zip(C, blocks)],
        ]
        best, bsel = -np.inf, None
        for sel in starts:
            sel = list(sel)
            if not feasible(sel):
                sel = [len(c) - 1 for c in C]                     # all unread is always feasible
            cur = objective(sel)
            for _ in range(rounds):
                changed = False
                for j in range(len(blocks)):
                    for a in range(sizes[j]):
                        if a == sel[j]:
                            continue
                        trial = sel.copy()
                        trial[j] = a
                        o = objective(trial)
                        n_eval += 1
                        if o > cur + 1e-12:
                            sel, cur, changed = trial, o, True
                if not changed:
                    break
            if cur > best:
                best, bsel = cur, tuple(sel)
        method = "coordinate"
    s = np.array([C[j][bsel[j]] for j in range(len(blocks))]).reshape(2, ns)
    return s, assemble(bsel).reshape(-1), best, {"method": method, "evaluations": n_eval, "sizes": sizes}
