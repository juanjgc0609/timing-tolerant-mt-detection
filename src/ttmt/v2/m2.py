"""M2: travel-time model with station terms in a Gaussian process (spec ``especificacion-modelo-v2.md`` M2).

For station i and phase b the arrival time from a source at x is

    t_ib(x) = T_b(x, r_i) + delta_ib(x),    delta_ib ~ GP(m_ib, k),    k Matern 3/2 in (x, y, z),

with T_b a published 1-D medium (``ttmt.traveltime``) and observations delta_ib(x_e) + eps, eps ~ N(0, s_n^2).
Hyperparameters (s_f, ell, s_n) by maximum marginal likelihood on the calibration events only (the mean m_ib in
closed form by generalised least squares); the predictive centre and tolerance at a new position are

    mu_ib(x) = T_b(x, r_i) + E[delta_ib(x) | cal],    sigma_ib(x) = (Var[delta_ib(x) | cal] + s_n^2)^(1/2),

the total tolerance (GP + reading/model noise) that M4 takes (Observation B.8: the smallest tolerance that covers
the timing error, here per position). The posterior variance falls near the calibration events (kriging), which a
constant term cannot do. Evaluation is by blocks of time (``time_blocks``): the terms are never estimated on the
events they predict (E53).

Interface (spec): ``fit_station_gp(cal_events, folds) -> TravelTimeModel``; ``TravelTimeModel.predict(x, station,
phase) -> (mu, sigma)``. Plain NumPy/SciPy; the hyperparameter search runs on at most ``n_hyper`` points (random,
seeded) and the posterior conditions on at most ``n_post`` (the most recent ones are not preferred: random, seeded),
so the cost is bounded (O(n^3) Cholesky).
"""
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize

from ..traveltime import HK, R_VPVS, SourceTimes

SQ3 = np.sqrt(3.0)


def matern32(A, B, ell):
    d = np.sqrt(np.maximum(((A[:, None, :] - B[None, :, :]) ** 2).sum(-1), 0.0))
    u = SQ3 * d / ell
    return (1.0 + u) * np.exp(-u)


class GP1:
    """GP of one (station, phase): constant mean (GLS) + s_f^2 Matern32(ell) + white noise. ``fit(X, y, g)``: with
    integer groups g (e.g. pick-quality classes) the white level is s_n * exp(c_g), c of the most frequent group 0
    and the others fitted with the rest by marginal likelihood (heteroscedastic variant); g None = one level."""

    def __init__(self, n_hyper=600, n_post=3000, seed=0):
        self.n_hyper, self.n_post, self.seed = n_hyper, n_post, seed

    def _noise(self, theta, gi):
        sn = np.exp(theta[2])
        c = np.r_[0.0, theta[3:]]
        return (sn * np.exp(c[gi])) ** 2

    def _nll(self, theta, X, y, gi):
        sf, ell = np.exp(theta[:2])
        K = sf ** 2 * matern32(X, X, ell) + np.diag(self._noise(theta, gi) + 1e-10)
        try:
            c = cho_factor(K, lower=True)
        except np.linalg.LinAlgError:
            return 1e10
        one = np.ones(len(X))
        Ki1 = cho_solve(c, one)
        m = float(Ki1 @ y / (Ki1 @ one))
        r = y - m
        return 0.5 * float(r @ cho_solve(c, r)) + float(np.log(np.diag(c[0])).sum())

    def _gi(self, g, n):
        if g is None or self.levels is None:
            return np.zeros(n, int)
        g = np.asarray(g)
        out = np.zeros(n, int)                      # an unseen group takes the reference level
        for j, v in enumerate(self.levels):
            out[g == v] = j
        return out

    def fit(self, X, y, g=None):
        X, y = np.asarray(X, float), np.asarray(y, float)
        if g is None:
            self.levels = None
        else:
            v, n = np.unique(np.asarray(g), return_counts=True)
            self.levels = list(v[np.argsort(-n, kind="stable")])
        gi = self._gi(g, len(y))
        ng = 1 if self.levels is None else len(self.levels)
        rng = np.random.default_rng(self.seed)
        h = rng.choice(len(X), min(len(X), self.n_hyper), replace=False)
        sd = float(np.std(y)) or 1e-3
        bounds = [(np.log(1e-4), np.log(10.0)), (np.log(0.1), np.log(500.0)), (np.log(1e-3), np.log(10.0))]
        bounds += [(-3.0, 3.0)] * (ng - 1)
        best = None
        for ell0 in (1.0, 5.0, 20.0):                       # a few starts: the likelihood can be multimodal in ell
            th0 = np.r_[np.log([sd, ell0, 0.5 * sd]), np.zeros(ng - 1)]
            res = minimize(self._nll, th0, args=(X[h], y[h], gi[h]), method="L-BFGS-B", bounds=bounds)
            if best is None or res.fun < best.fun:
                best = res
        self.theta = best.x
        self.sf, self.ell, self.sn = np.exp(best.x[:3])
        p = rng.choice(len(X), min(len(X), self.n_post), replace=False)
        self.X = X[p]
        K = self.sf ** 2 * matern32(self.X, self.X, self.ell) + np.diag(self._noise(self.theta, gi[p]) + 1e-10)
        self.c = cho_factor(K, lower=True)
        one = np.ones(len(p))
        Ki1 = cho_solve(self.c, one)
        self.m = float(Ki1 @ y[p] / (Ki1 @ one))
        self.alpha = cho_solve(self.c, y[p] - self.m)
        self.n = len(X)
        return self

    def predict(self, Xs, g=None):
        """posterior mean of delta and total sd (GP posterior + the white level of the reading's group)."""
        Xs = np.atleast_2d(np.asarray(Xs, float))
        ks = self.sf ** 2 * matern32(Xs, self.X, self.ell)
        mu = self.m + ks @ self.alpha
        v = self.sf ** 2 - np.einsum("ij,ji->i", ks, cho_solve(self.c, ks.T))
        return mu, np.sqrt(np.maximum(v, 0.0) + self._noise(self.theta, self._gi(g, len(Xs))))


class Const1:
    """the ablation of the spec: a constant term (mean of the calibration residuals) and its sd as tolerance."""

    def fit(self, X, y, g=None):
        y = np.asarray(y, float)
        self.m, self.s, self.n = float(np.mean(y)), float(np.std(y, ddof=1)) if len(y) > 1 else np.nan, len(y)
        return self

    def predict(self, Xs, g=None):
        n = len(np.atleast_2d(Xs))
        return np.full(n, self.m), np.full(n, self.s)


class TravelTimeModel:
    """1-D times plus one fitted term per (station, phase). Stations: dict code -> (x, y, z) km in the local frame."""

    def __init__(self, stations, terms, model=HK, r=R_VPVS):
        self.stations, self.terms = stations, terms
        self.tt = SourceTimes(model, r)

    def t1d(self, x, station, phase):
        s = np.asarray(self.stations[station], float)[None]
        return float(self.tt(np.asarray(x, float), s)[0 if phase == "P" else 1, 0])

    def predict(self, x, station, phase):
        """(mu, sigma) of the arrival time (s after origin) at station for a source at x (x, y, z km)."""
        d, s = self.terms[(station, phase)].predict(np.asarray(x, float)[None])
        return self.t1d(x, station, phase) + float(d[0]), float(s[0])


def time_blocks(origin_s, k=5):
    """fold index (0..k-1) of every event: k contiguous blocks of equal count in time order."""
    o = np.argsort(np.asarray(origin_s, float), kind="stable")
    f = np.empty(len(o), int)
    for j, blk in enumerate(np.array_split(o, k)):
        f[blk] = j
    return f


def residuals(cal_events, stations, model=HK, r=R_VPVS):
    """cal_events: iterable of (x (3,), {(station, phase): tt_obs}). Returns {(station, phase): (X (n, 3), res (n,),
    event index (n,))}, res = observed travel time - 1-D time."""
    tt = SourceTimes(model, r)
    codes = list(stations)
    S = np.array([stations[c] for c in codes], float)
    out = {}
    for e, (x, obs) in enumerate(cal_events):
        T = tt(np.asarray(x, float), S)
        for (sta, ph), t in obs.items():
            if sta not in stations:
                continue
            out.setdefault((sta, ph), ([], [], []))
            out[(sta, ph)][0].append(x)
            out[(sta, ph)][1].append(t - T[0 if ph == "P" else 1, codes.index(sta)])
            out[(sta, ph)][2].append(e)
    return {k: (np.array(a, float), np.array(b, float), np.array(c, int)) for k, (a, b, c) in out.items()}


def fit_station_gp(cal_events, stations, folds=None, kind="gp", min_n=30, model=HK, r=R_VPVS, **kw):
    """fit one term per (station, phase) with at least min_n calibration readings; ``folds``: None (fit on all) or
    an index array over cal_events with the fold to LEAVE OUT (-1 = keep), used by the block validation."""
    res = residuals(cal_events, stations, model, r)
    terms = {}
    for key, (X, y, ev) in res.items():
        keep = np.ones(len(y), bool) if folds is None else np.asarray(folds)[ev] < 0
        if keep.sum() < min_n:
            continue
        terms[key] = (GP1(**kw) if kind == "gp" else Const1()).fit(X[keep], y[keep])
    return TravelTimeModel(stations, terms, model, r)
