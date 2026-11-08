"""M5 decision step (ttmt.v2.m5): the mixture e-value S_mix with a Student t background (theory D.1, D.2)."""
import math

import numpy as np

from ttmt.v2.m5 import GaussBackground, TBackground, evalue_mix, prior_draws


def _t_draw(rng, n, b, nu):
    """isotropic 3-D Student t: b * z / sqrt(chi2_nu / nu)."""
    z = rng.normal(size=(n, 3))
    return b * z / np.sqrt(rng.chisquare(nu, size=(n, 1)) / nu)


def _scene(rng, H, K, b, nu):
    A = rng.normal(size=(H, K, 3, 6)) * 0.3
    v = _t_draw(rng, H * K, b, nu).reshape(H, K, 3)
    return v, A


def test_t_fit_recovers_parameters():
    rng = np.random.default_rng(1)
    for b, nu in ((1.0, 3.0), (2.5, 8.0)):
        f = TBackground.fit(_t_draw(rng, 40000, b, nu))
        assert abs(f.b / b - 1) < 0.03 and abs(math.log(f.nu / nu)) < 0.15, (f.b, f.nu)


def test_t_density_integrates_to_one():
    f = TBackground(1.3, 2.5)
    r = np.linspace(0, 400, 400001)
    dens = 4 * np.pi * r ** 2 * np.exp(f.logpdf(np.stack([r, 0 * r, 0 * r], -1)))
    assert abs(np.trapezoid(dens, r) - 1) < 2e-3
    assert abs(f.cdf_norm(400.0) - 1) < 1e-3


def test_eps_one_is_neutral():
    rng = np.random.default_rng(2)
    v, A = _scene(rng, 5, 6, 1.0, 3.0)
    f = TBackground(1.0, 3.0)
    S = evalue_mix(v, A, f.logpdf(v), 1.0, prior_draws(rng.normal(size=(4, 6)), 50, rng), 1.0)
    assert np.allclose(S, 0.0)


def test_evalue_is_calibrated_under_exact_background():
    """Theorem D.2: E exp(S_mix) = 1 and P(S >= log 1/p) <= p when the blocks follow f exactly. pi is put at the
    noise level (the hardest case for the level; a pi far above the noise makes S very negative)."""
    rng = np.random.default_rng(3)
    H, K, b, nu = 20000, 10, 1.0, 3.0
    v, A = _scene(rng, H, K, b, nu)
    f = TBackground(b, nu)
    draws = prior_draws(rng.normal(size=(8, 6)), 200, rng, amp_sd=0.5)
    S = evalue_mix(v, A, f.logpdf(v), 0.05, draws, s0=b, gamma=0.3)
    e = np.exp(S)
    se = e.std() / math.sqrt(H)
    assert abs(e.mean() - 1) < max(4 * se, 0.05), (e.mean(), se)
    for p in (1e-1, 1e-2):
        assert (S >= math.log(1 / p)).mean() <= p


def test_gaussian_background_inflates_transient_evidence():
    """Proposition D.1(a): one transient block of amplitude a raises S_mix by ~ a^2 / (2 s^2) under a gaussian
    background, but only by ~ (3 + nu) log a under the t background."""
    rng = np.random.default_rng(4)
    K = 8
    A = rng.normal(size=(1, K, 3, 6)) * 0.3
    v0 = _t_draw(rng, K, 1.0, 3.0).reshape(1, K, 3)
    draws = prior_draws(rng.normal(size=(8, 6)), 400, rng)
    ft, fg = TBackground(1.0, 3.0), GaussBackground(1.0)
    gain = {}
    for name, f in (("t", ft), ("g", fg)):
        out = []
        for a in (10.0, 100.0):
            v = v0.copy()
            v[0, 0] = a * np.array([1.0, 0.0, 0.0])
            out.append(evalue_mix(v, A, f.logpdf(v), 0.05, draws, s0=1.0)[0])
        gain[name] = out[1] - out[0]
    assert gain["t"] < (3 + 3.0) * math.log(10) + 1.0
    assert gain["g"] > 100 * max(gain["t"], 1.0)


def test_radial_background_fits_peak_reads():
    """on max-energy-of-9 gaussian reads (the registered reader) with 2 % transients, the radial background fits the
    norm law out of sample (PIT deciles near 0.1, tail ratio near 1), where the isotropic t does not."""
    from ttmt.v2.m5 import RadialBackground
    rng = np.random.default_rng(5)
    X = rng.normal(size=(60000, 9, 3))
    k = np.argmax((X ** 2).sum(-1), axis=1)
    R = X[np.arange(len(X)), k]
    R = R * np.where(rng.random((len(R), 1)) < 0.02, np.exp(rng.normal(1.5, 0.7, (len(R), 1))), 1.0)
    fit, test = R[:30000], R[30000:]
    rb, t = RadialBackground.fit(fit), TBackground.fit(fit)
    rn = np.linalg.norm(test, axis=1)
    dev = {n: np.abs(np.histogram(f.cdf_norm(rn), np.linspace(0, 1, 11))[0] / len(rn) - 0.1).max()
           for n, f in (("r", rb), ("t", t))}
    assert dev["r"] < 0.01 and dev["t"] > 0.05, dev
    tail = float((rb.cdf_norm(rn) > 1 - 1e-2).mean() / 1e-2)
    assert 0.7 < tail < 1.4, tail
    r = np.exp(np.linspace(-8, 8, 400001))
    dens = 4 * np.pi * r ** 3 * np.exp(rb.logpdf(np.stack([r, 0 * r, 0 * r], -1)))      # d r = r d log r
    assert abs(np.trapezoid(dens, np.log(r)) - 1) < 5e-3


def test_evalue_calibrated_with_radial_background_on_peak_reads():
    """Theorem D.2 with the fitted radial background on peak reads: null rates at most p (pi at the noise level)."""
    from ttmt.v2.m5 import RadialBackground
    rng = np.random.default_rng(6)

    def reads(n):
        X = rng.normal(size=(n, 9, 3))
        return X[np.arange(n), np.argmax((X ** 2).sum(-1), axis=1)]

    rb = RadialBackground.fit(reads(40000))
    H, K = 20000, 10
    v = reads(H * K).reshape(H, K, 3)
    A = rng.normal(size=(H, K, 3, 6)) * 0.3
    draws = prior_draws(rng.normal(size=(8, 6)), 200, rng)
    S = evalue_mix(v, A, rb.logpdf(v), 0.05, draws, s0=1.0)
    for p in (1e-1, 1e-2):
        assert (S >= math.log(1 / p)).mean() <= p


def test_zero_background_density_gives_infinite_evidence_not_nan():
    """a block with f_k(v_k) = 0 (outside the fitted support) has l_k = +inf; S_mix must be +inf, never NaN (a NaN
    compares False with every threshold and would be counted as a non-rejection)."""
    rng = np.random.default_rng(7)
    v, A = _scene(rng, 3, 4, 1.0, 3.0)
    lf = TBackground(1.0, 3.0).logpdf(v)
    lf[1, 2] = -np.inf
    S = evalue_mix(v, A, lf, 0.05, prior_draws(rng.normal(size=(4, 6)), 50, rng), 1.0)
    assert np.isposinf(S[1]) and np.all(np.isfinite(S[[0, 2]]))
