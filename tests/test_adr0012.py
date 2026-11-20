"""ADR-0012: operator weighted like the data (Q of diag(1/sd) G) when each station is divided by its noise level."""
import numpy as np
import pytest

from ttmt.radiation import _rand_dc, basis, node_operators, radiation_operator, station_weights

RNG = np.random.default_rng(20261120)


def _setup(ns=7):
    STA = np.c_[RNG.uniform(-20, 20, (ns, 2)), np.zeros(ns)]
    XN = np.array([[0.37, 1.37, 8.0], [-3.63, 4.37, 12.0]])
    sd = np.geomspace(1.0, 10.0, ns)                                      # ratio 10 between stations
    RNG.shuffle(sd)
    return STA, XN, sd


def _f(Q, v):
    y = v @ Q
    return float(y @ y) / float(v @ v)


def test_weighted_operator_keeps_noise_free_event_in_span():
    STA, XN, sd = _setup()
    w = station_weights(1.0 / sd, len(STA))
    _G0, Q0 = node_operators(XN, STA)
    _Gw, Qw = node_operators(XN, STA, 1.0 / sd)
    f_w, f_0 = [], []
    for g, x in enumerate(XN):
        G = radiation_operator(x, STA)
        for _ in range(50):
            m = _rand_dc(RNG)
            m6 = np.array([m[0, 0], m[1, 1], m[2, 2], m[0, 1], m[0, 2], m[1, 2]])
            v = w * (G @ m6)                                              # noise-free data divided by sd per station
            f_w.append(_f(Qw[g], v))
            f_0.append(_f(Q0[g], v))
    assert min(f_w) >= 0.99
    assert np.median(f_0) < 0.95


def test_default_path_unchanged():
    STA, XN, _sd = _setup()
    Gs, Qst = node_operators(XN, STA)
    ref = np.stack([radiation_operator(x, STA) for x in XN])
    assert np.array_equal(Gs, ref)                                        # bit for bit
    for g in range(len(XN)):
        assert np.array_equal(Qst[g][:, :basis(ref[g]).shape[1]], basis(ref[g]))


def test_unit_weights_equal_default_span():
    STA, XN, _sd = _setup()
    _G0, Q0 = node_operators(XN, STA)
    _G1, Q1 = node_operators(XN, STA, np.ones(len(STA)))
    for g in range(len(XN)):
        assert np.allclose(Q0[g] @ Q0[g].T, Q1[g] @ Q1[g].T, atol=1e-12)


@pytest.mark.parametrize("w", [np.ones(3), np.array([1, 0, 1, 1, 1, 1, 1.0]), np.array([1, np.nan, 1, 1, 1, 1, 1.0])])
def test_bad_weights_rejected(w):
    STA, XN, _sd = _setup()
    with pytest.raises(ValueError):
        node_operators(XN, STA, w)
