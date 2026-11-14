"""M2 (ttmt.v2.m2): the GP term recovers a smooth position-dependent delay that a constant term cannot, its
tolerance covers held-out errors, and the time blocks never mix events of one fold."""
import numpy as np

from ttmt.v2.m2 import Const1, GP1, TravelTimeModel, fit_station_gp, time_blocks


def _field(X):
    return 0.3 * np.sin(X[:, 0] / 6.0) + 0.2 * np.cos(X[:, 1] / 8.0) + 0.1


def test_gp_beats_constant_and_covers():
    r = np.random.default_rng(1)
    X = np.c_[r.uniform(-20, 20, (600, 2)), r.uniform(2, 12, 600)]
    y = _field(X) + r.normal(0, 0.05, 600)
    tr, te = slice(0, 400), slice(400, None)
    gp = GP1(seed=0).fit(X[tr], y[tr])
    c = Const1().fit(X[tr], y[tr])
    mg, sg = gp.predict(X[te])
    mc, _ = c.predict(X[te])
    assert np.mean((y[te] - mg) ** 2) < 0.3 * np.mean((y[te] - mc) ** 2)
    cov = np.mean(np.abs(y[te] - mg) <= 2 * sg)
    assert 0.9 <= cov <= 0.99
    assert 0.03 < gp.sn < 0.08                         # the white level is identified


def test_posterior_sd_grows_away_from_data():
    r = np.random.default_rng(2)
    X = np.c_[r.uniform(-5, 5, (200, 2)), np.full(200, 5.0)]
    gp = GP1(seed=0).fit(X, _field(X) + r.normal(0, 0.05, 200))
    _, s_in = gp.predict(np.array([[0.0, 0.0, 5.0]]))
    _, s_out = gp.predict(np.array([[60.0, 60.0, 5.0]]))
    assert s_out[0] > 2 * s_in[0]


def test_time_blocks_contiguous_and_balanced():
    t = np.random.default_rng(3).permutation(103).astype(float)
    f = time_blocks(t, 5)
    assert sorted(np.bincount(f).tolist()) == [20, 20, 21, 21, 21]
    for k in range(4):
        assert t[f == k].max() < t[f == k + 1].min()


def test_fit_station_gp_leaves_fold_out_and_predicts_time():
    r = np.random.default_rng(4)
    sta = {"A": (0.0, 0.0, 0.0)}
    ev = [(np.r_[r.uniform(-10, 10, 2), 6.0], {}) for _ in range(80)]
    m = TravelTimeModel(sta, {})
    ev = [(x, {("A", "P"): m.t1d(x, "A", "P") + 0.2}) for x, _ in ev]
    folds = np.where(np.arange(80) < 20, 0, -1)
    M = fit_station_gp(ev, sta, folds=folds, kind="const")
    assert M.terms[("A", "P")].n == 60
    mu, s = M.predict(ev[0][0], "A", "P")
    assert abs(mu - ev[0][1][("A", "P")]) < 1e-9


def test_heteroscedastic_levels_recovered_and_one_group_equals_plain():
    r = np.random.default_rng(5)
    X = np.c_[r.uniform(-20, 20, (600, 2)), r.uniform(2, 12, 600)]
    g = np.where(r.random(600) < 0.6, 1.0, 0.5)
    y = _field(X) + r.normal(0, np.where(g == 1.0, 0.04, 0.16))
    gp = GP1(seed=0).fit(X, y, g)
    assert gp.levels[0] == 1.0
    ratio = np.exp(gp.theta[3])
    assert 2.5 < ratio < 6.0                            # true ratio 4
    plain = GP1(seed=0).fit(X, y)
    same = GP1(seed=0).fit(X, y, np.ones(600))
    assert np.allclose(plain.predict(X[:5])[0], same.predict(X[:5], np.ones(5))[0])
