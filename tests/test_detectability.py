"""T4.20 ceiling (ttmt.detectability) and static station terms (ttmt.associate.station_terms)."""
import numpy as np
from scipy import stats

from ttmt import detectability as D
from ttmt.associate import station_terms as ST
from ttmt.associate.velocity import Velocity1D
from ttmt.atlas.schema import Picks, Truth


def test_threshold_meets_budget_under_null():
    # Monte Carlo under H0: the per-trial exceedance rate of chi2_p at alpha is alpha
    rng = np.random.default_rng(0)
    a = 1e-3
    for p in (1, 6):
        x = (rng.standard_normal((400_000, p)) ** 2).sum(1)
        rate = (x > D.threshold(p, a)).mean()
        assert abs(rate - a) < 4 * np.sqrt(a / 4e5)


def test_pd_matches_monte_carlo_and_orders_detectors():
    rng = np.random.default_rng(1)
    a, lam = 1e-4, 30.0
    for p in (1, 6):
        mu = np.zeros(p)
        mu[0] = np.sqrt(lam)
        x = ((rng.standard_normal((200_000, p)) + mu) ** 2).sum(1)
        assert abs((x > D.threshold(p, a)).mean() - D.pd(lam, p, a)) < 0.01
    # more unknowns cost detection probability at equal energy: p = 1 >= p = 6 >= p = 12
    v = [float(D.pd(lam, p, a)) for p in (1, 6, 12)]
    assert v[0] >= v[1] >= v[2]
    # p = 1 at lambda = z^2 of the one-sided level is close to one half (two-sided chi2_1)
    z = stats.norm.isf(a)
    assert 0.45 < float(D.pd(z ** 2, 1, 2 * a)) < 0.55


def test_union_bound_alpha_and_lam_lower():
    assert np.isclose(D.trial_alpha(30, 100, 1e4), 30 / 3.6e9)
    assert np.isclose(D.lam_lower([3.0, 4.0, np.nan]), 25.0)
    assert np.isnan(D.lam_lower([np.nan]))
    c = D.ceiling([25.0, 25.0], [4, 8], 1e-6)
    assert c["p=K"][0] > c["p=K"][1]          # same energy over more free amplitudes is harder


def _scene(rng, n_src=8, ns=6, delta=None, sd=0.0):
    vel = Velocity1D.jma2001()
    sta = np.c_[rng.uniform(-30, 30, (ns, 2)), np.zeros(ns)]
    src = np.c_[rng.uniform(-20, 20, (n_src, 2)), rng.uniform(5, 15, n_src)]
    delta = np.zeros((ns, 2)) if delta is None else delta
    ev = np.repeat(np.arange(n_src), 3)                     # three injections per source
    t = np.arange(len(ev)) * 50.0
    arr = np.stack([t[i] + np.asarray(vel.times(src[e], sta)).T for i, e in enumerate(ev)])
    obs = arr + delta[None] + rng.normal(0, sd, arr.shape)
    m = len(ev)
    tr = Truth(np.array([f"x#{i}" for i in range(m)]), t, src[ev], np.zeros(m), arr, np.ones((m, ns)), 0.0,
               float(t[-1] + 60), np.ones(m, bool), -np.ones(m, int), np.full(m, np.nan), obs)
    return tr, sta, vel, ev


def test_station_terms_recover_delays_cross_fitted():
    rng = np.random.default_rng(2)
    ns = 6
    delta = np.c_[rng.normal(0, 0.2, ns), rng.normal(0, 0.4, ns)]
    tr, sta, vel, ev = _scene(rng, ns=ns, delta=delta, sd=0.02)
    R = ST.residuals(tr, sta, vel)
    K = ST.source_keys(tr.xyz)
    assert len(np.unique(K)) == 8
    ks, sres = ST.per_source(R, K)
    d = ST.record_terms(sres, ks, np.array([0, 1]))         # the record holds sources 0 and 1: estimated from 2..7
    assert np.abs(d - delta).max() < 0.05
    # leaving out sources that saw a station drops it below min_sources -> 0
    d0 = ST.record_terms(sres, ks, ks[:-2], min_sources=3)
    assert np.all(d0 == 0.0)
    # applying the terms removes the delay from the picks
    e, s, ph = np.nonzero(np.isfinite(tr.arr_t_obs))
    pk = Picks(s, ph, tr.arr_t_obs[e, s, ph], np.ones(len(s)))
    cp = ST.apply(pk, d)
    assert np.abs(cp.t - tr.arr_t[e, s, ph]).max() < 0.15
    assert np.abs(cp.t - tr.arr_t[e, s, ph]).mean() < np.abs(pk.t - tr.arr_t[e, s, ph]).mean() / 4


def test_station_terms_zero_without_misfit_and_origin_variant():
    rng = np.random.default_rng(3)
    tr, sta, vel, _ = _scene(rng, delta=np.full((6, 2), 0.3))   # a constant delay is an origin shift
    ks, sres = ST.per_source(ST.residuals(tr, sta, vel), ST.source_keys(tr.xyz), remove_origin=True)
    assert np.abs(ST.terms(sres)).max() < 1e-9
    ks, sres = ST.per_source(ST.residuals(tr, sta, vel), ST.source_keys(tr.xyz))
    assert np.allclose(ST.terms(sres), 0.3)


def test_fill_missing_gives_the_common_shift_to_stations_without_a_term():
    d = np.array([[0.2, 0.5], [0.0, 0.3], [0.4, 0.0], [0.3, 0.4]])
    f = ST.fill_missing(d)
    assert f[1, 0] == 0.3 and f[2, 1] == 0.4                 # median of the estimated terms of the phase
    assert np.array_equal(f[d != 0], d[d != 0])              # estimated terms untouched
    assert np.array_equal(ST.fill_missing(np.zeros((3, 2))), np.zeros((3, 2)))   # nothing estimated: nothing to fill
