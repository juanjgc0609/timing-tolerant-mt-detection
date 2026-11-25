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


def _fake_net(monkeypatch):
    """r1.Net with every data access replaced: records the weights it passes to node_operators."""
    import types
    from ttmt import r1
    scale = np.array([1.0, 2.0, 4.0])
    R = types.SimpleNamespace(subnets={"H7": ["a", "b", "c"]}, frame=lambda: None, catalogue_local=lambda c: [],
                              stations=lambda n, c: (np.zeros((3, 3)), None), segment=lambda s: (0.0, 1.0),
                              partition=None, load_waveforms=lambda n, s, t: (None, np.ones((3, 3, 8))),
                              station_scale=lambda d: scale, cal=(0.0, 1.0), m_terms=1.0, depths=None, vp=6.0, vs=3.5)
    seen = []
    monkeypatch.setattr(r1, "R", R)
    monkeypatch.setattr(r1, "legacy_t_end", lambda p, s: 1.0)
    monkeypatch.setattr(r1, "estimate_peak_terms", lambda *a: None)
    monkeypatch.setattr(r1, "coarse_grid", lambda d: np.zeros((1, 3)))
    monkeypatch.setattr(r1, "node_operators", lambda XN, STA, w=None: (seen.append(w), (None, None))[1])
    return r1, scale, seen


def test_net_passes_inverse_scale_only_when_weighted(monkeypatch):
    r1, scale, seen = _fake_net(monkeypatch)
    assert r1.Net("cal").weighted_q is False and seen[-1] is None
    assert r1.Net("cal", weighted_q=True).weighted_q is True
    assert np.array_equal(seen[-1], 1.0 / scale)


def test_hybrid_refuses_a_mixed_weighting():
    from ttmt import r1
    for p1d, wq in (({"weighted_q": True}, False), ({}, True)):
        with pytest.raises(ValueError, match="refusing the mix"):
            r1.run_hybrid("cal", {}, p1d, weighted_q=wq)
