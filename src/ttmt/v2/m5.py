"""M5 decision step: the mixture e-value S_mix with a polynomial-tailed (Student t) background.

Specification ``especificacion-modelo-v2.md`` M5 ("Decidir"), theory ``teoria-fase5.md`` D.1 and D.2.

For one hypothesis with blocks k = 1..K (a block is one phase at one station: a read 3-vector v_k and the radiation
operator A_k, 3 x 6, of the hypothesis), a background density f_k per block, a contamination rate eps_k and a model
error width s_k = s0_k + gamma ||A_k m||,

    l_k(m)  = log[ eps_k f_k(v_k) + (1 - eps_k) N(v_k; A_k m, s_k^2 I_3) ] - log f_k(v_k)
    S_mix   = log( (1/J) sum_j exp( sum_k l_k(m_j) ) ),   m_j ~ pi drawn independently of the data.

Theorem D.2: if under H0 the blocks are independent with densities exactly f_k, E[exp(S_mix)] = 1 for every pi fixed
in advance (also for the J-draw Monte Carlo average), so P(S_mix >= log(1/p)) <= p with no calibration. With a
misspecified null f*_k <= (1 + beta_k) f_k on the region where the model pulse has mass, the level is at most
p prod_k (1 + (1 - eps_k) beta_k). Proposition D.1(a): a gaussian background makes a transient's evidence quadratic
in its amplitude; the background must have polynomial tails, hence the isotropic Student t here.

Exploratory (no registered run uses it). Conventions are arguments; nothing in ``ttmt.detector`` changes.

  TBackground      isotropic 3-D Student t (scale b, degrees of freedom nu), ML fit on read vectors of noise
  RadialBackground isotropic, smooth empirical radial law + Pareto tail (the background used for the registered
                   peak reader, where a single t cannot fit; see its docstring)
  GaussBackground  isotropic 3-D gaussian (the ablation that Proposition D.1(a) predicts fails on transients)
  prior_draws      the J tensors m_j ~ pi: calibration tensors, lognormal amplitude spread, isotropic jitter
  evalue_mix       S_mix for a batch of hypotheses
"""
import math

import numpy as np

LOG2PI = math.log(2.0 * math.pi)


class TBackground:
    """isotropic Student t in R^3: f(v) = c(nu) b^-3 (1 + |v|^2 / (nu b^2))^(-(nu + 3) / 2)."""

    def __init__(self, b, nu):
        self.b, self.nu = float(b), float(nu)

    def logpdf(self, V):
        r2 = np.sum(np.asarray(V, float) ** 2, axis=-1)
        return _t_logpdf_r2(r2, self.b, self.nu)

    @classmethod
    def fit(cls, V, nu_bounds=(0.3, 300.0)):
        """maximum likelihood (b, nu) from read vectors V (n, 3). Profile: for fixed nu the scale solves the EM
        fixed point b^2 = mean(w r^2) / 3, w = (nu + 3) / (nu + r^2 / b^2); nu by golden-section on log nu."""
        r2 = np.sum(np.asarray(V, float) ** 2, axis=-1)
        r2 = r2[np.isfinite(r2)]

        def prof(lnu):
            nu = math.exp(lnu)
            b2 = float(np.mean(r2)) / 3.0 + 1e-300
            for _ in range(200):
                w = (nu + 3.0) / (nu + r2 / b2)
                b2n = float(np.mean(w * r2)) / 3.0 + 1e-300
                if abs(b2n - b2) <= 1e-10 * b2:
                    b2 = b2n
                    break
                b2 = b2n
            return float(np.sum(_t_logpdf_r2(r2, math.sqrt(b2), nu))), math.sqrt(b2)

        lo, hi = math.log(nu_bounds[0]), math.log(nu_bounds[1])
        g = (math.sqrt(5.0) - 1.0) / 2.0
        a, d = lo, hi
        x1, x2 = d - g * (d - a), a + g * (d - a)
        f1, f2 = prof(x1)[0], prof(x2)[0]
        for _ in range(60):
            if f1 >= f2:
                d, x2, f2 = x2, x1, f1
                x1 = d - g * (d - a)
                f1 = prof(x1)[0]
            else:
                a, x1, f1 = x1, x2, f2
                x2 = a + g * (d - a)
                f2 = prof(x2)[0]
        lnu = 0.5 * (a + d)
        return cls(prof(lnu)[1], math.exp(lnu))

    def cdf_norm(self, r):
        """P(|v| <= r): |v|^2 / (3 b^2) is F(3, nu)."""
        from scipy.stats import f as fdist
        return fdist.cdf(np.asarray(r, float) ** 2 / (3.0 * self.b ** 2), 3, self.nu)


class RadialBackground:
    """isotropic density with a smooth empirical radial law and a polynomial (Pareto) tail:

        f(v) = f_r(|v|) / (4 pi |v|^2);  below r0 (the q0 quantile of the fit norms) f_r is a gaussian kernel density
        of log |v| (Silverman bandwidth) holding mass q0; above r0, P(|v| > r) = (1 - q0) (r / r0)^-alpha with the
        Hill estimate of alpha (capped at alpha_max), so f ~ |v|^-(3 + alpha) (Proposition D.1: polynomial tail of index nu = alpha).

    Why not the isotropic t: the registered reader takes the energy peak within +-4 samples, so the norm of a read
    vector in noise has almost no mass near 0 and a concentrated bulk, which no single t (nor a generalised-t radial
    law, whose shape runs to its bound) fits together with a heavy tail (see tests/test_m5.py and the card)."""

    NGRID = 4096

    def __init__(self, y_grid, logfy_grid, r0, alpha, q0):
        self.y, self.lfy = np.asarray(y_grid, float), np.asarray(logfy_grid, float)
        self.r0, self.alpha, self.q0 = float(r0), float(alpha), float(q0)
        self.nu = self.alpha

    @classmethod
    def fit(cls, V, q0=0.99, alpha_max=10.0):
        from scipy.special import ndtr
        r = np.linalg.norm(np.asarray(V, float), axis=-1)
        r = np.sort(r[np.isfinite(r) & (r > 0)])
        r0 = float(np.quantile(r, q0))
        tail = r[r > r0]
        alpha = 1.0 / float(np.mean(np.log(tail / r0))) if len(tail) >= 10 else alpha_max
        alpha = min(alpha, alpha_max)                    # a heavier tail than estimated is conservative
        y = np.log(r[r <= r0])
        sd = min(float(np.std(y)), float(np.subtract(*np.quantile(y, [0.75, 0.25]))) / 1.349)
        h = 0.9 * sd * len(y) ** -0.2
        y0 = math.log(r0)
        grid = np.linspace(y.min() - 6 * h, y0, cls.NGRID)
        ys = y[:: max(1, len(y) // 20000)]                                   # kernel centres (subsampled if large)
        dens = np.zeros_like(grid)
        for c in np.array_split(ys, max(1, len(ys) // 2000)):
            dens += np.exp(-0.5 * ((grid[:, None] - c[None]) / h) ** 2).sum(1)
        dens /= len(ys) * h * math.sqrt(2 * math.pi)
        mass = float(np.mean(ndtr((y0 - ys) / h)))                            # kernel mass below y0
        with np.errstate(divide="ignore"):
            lfy = np.log(dens * q0 / mass)
        return cls(grid, lfy, r0, alpha, q0)

    def _logfr(self, r):
        r = np.maximum(np.asarray(r, float), 1e-300)
        y = np.log(r)
        bulk = np.interp(y, self.y, self.lfy, left=-np.inf) - y
        tail = math.log((1 - self.q0) * self.alpha) + self.alpha * math.log(self.r0) - (self.alpha + 1) * y
        return np.where(r <= self.r0, bulk, tail)

    def logpdf(self, V):
        r = np.linalg.norm(np.asarray(V, float), axis=-1)
        with np.errstate(divide="ignore"):
            return self._logfr(r) - math.log(4 * math.pi) - 2 * np.log(np.maximum(r, 1e-300))

    def cdf_norm(self, r):
        r = np.asarray(r, float)
        fy = np.exp(self.lfy)
        cy = np.concatenate([[0.0], np.cumsum(0.5 * (fy[1:] + fy[:-1]) * np.diff(self.y))])
        with np.errstate(divide="ignore"):
            yb = np.log(np.maximum(r, 1e-300))
        bulk = np.interp(yb, self.y, cy, left=0.0)
        tail = 1 - (1 - self.q0) * (np.maximum(r, self.r0) / self.r0) ** -self.alpha
        return np.where(r <= self.r0, bulk, tail)


class GaussBackground:
    """isotropic gaussian N(0, s^2 I_3); fit: s^2 = mean |v|^2 / 3 (ML)."""

    def __init__(self, s):
        self.s = float(s)

    def logpdf(self, V):
        r2 = np.sum(np.asarray(V, float) ** 2, axis=-1)
        return -1.5 * (LOG2PI + 2.0 * math.log(self.s)) - r2 / (2.0 * self.s ** 2)

    @classmethod
    def fit(cls, V):
        return cls(math.sqrt(float(np.mean(np.sum(np.asarray(V, float) ** 2, axis=-1))) / 3.0))


def _t_logpdf_r2(r2, b, nu):
    c = math.lgamma((nu + 3.0) / 2.0) - math.lgamma(nu / 2.0) - 1.5 * math.log(nu * math.pi) - 3.0 * math.log(b)
    return c - 0.5 * (nu + 3.0) * np.log1p(np.asarray(r2) / (nu * b * b))


def prior_draws(m_cal, J, rng, amp_sd=0.5, jitter=0.25, scale=None):
    """J tensors from the smoothed calibration distribution: m = exp(amp_sd z) (m_c + jitter |m_c| g / sqrt 6),
    m_c drawn uniformly from the rows of m_cal (n, 6), z ~ N(0, 1), g ~ N(0, I_6). scale (n,) optional: each m_c is
    first multiplied by scale[c] (used to move pi to a chosen amplitude level)."""
    m_cal = np.asarray(m_cal, float)
    if scale is not None:
        m_cal = m_cal * np.asarray(scale, float)[:, None]
    c = rng.integers(len(m_cal), size=J)
    mc = m_cal[c]
    nm = np.linalg.norm(mc, axis=1, keepdims=True)
    m = mc + jitter * nm * rng.normal(size=(J, 6)) / math.sqrt(6.0)
    return m * np.exp(amp_sd * rng.normal(size=(J, 1)))


def evalue_mix(v, A, logf, eps, draws, s0, gamma=0.3, chunk=None):
    """S_mix (log e-value) for H hypotheses.

    v (H, K, 3) read vectors; A (H, K, 3, 6) operators; logf (H, K) log f_k(v_k) of the background, or a list of such
    arrays (one S_mix per background, sharing the model term); eps (K,) or (H, K) contamination rate in (0, 1];
    draws (J, 6) tensors m_j ~ pi (independent of v); s0 (K,) or (H, K) model error at zero amplitude; gamma
    proportional model error. Returns (H,) S_mix, or a list of them if logf is a list."""
    v = np.asarray(v, float)
    A = np.asarray(A, float)
    H, K = v.shape[:2]
    many = isinstance(logf, (list, tuple))
    logfs = [np.asarray(x, float) for x in (logf if many else [logf])]
    eps = np.broadcast_to(np.asarray(eps, float), (H, K))
    s0 = np.broadcast_to(np.asarray(s0, float), (H, K))
    with np.errstate(divide="ignore"):
        le = np.log(eps)[..., None]
        l1 = np.log1p(-eps)[..., None]                                                  # -inf where eps = 1
    draws = np.asarray(draws, float)
    J = len(draws)
    step = chunk or max(1, int(4e6 // max(1, J * K * 3)))
    outs = [np.empty(H) for _ in logfs]
    for h0 in range(0, H, step):
        sl = slice(h0, min(H, h0 + step))
        h = sl.stop - sl.start
        P = np.matmul(A[sl].reshape(h, K * 3, 6), draws.T).reshape(h, K, 3, J)       # (h, K, 3, J)
        s = s0[sl][..., None] + gamma * np.sqrt(np.sum(P * P, axis=2))                 # (h, K, J)
        d2 = np.sum((v[sl][..., None] - P) ** 2, axis=2)
        logN = -1.5 * (LOG2PI + 2.0 * np.log(s)) - d2 / (2.0 * s * s) + l1[sl]
        for lf, out in zip(logfs, outs):
            tot = np.logaddexp(le[sl], logN - lf[sl][..., None]).sum(1)                # (h, J)
            mx = tot.max(axis=1, keepdims=True)
            with np.errstate(invalid="ignore"):
                val = mx[:, 0] + np.log(np.mean(np.exp(tot - mx), axis=1))
            out[sl] = np.where(np.isposinf(mx[:, 0]), np.inf, val)               # f_k(v_k) = 0: S = +inf, not NaN
    return outs if many else outs[0]
