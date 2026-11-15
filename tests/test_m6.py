"""Tests of M6 (ttmt.v2.m6): joint refit of overlapping atoms and the one-against-two resolvability test."""
import numpy as np

from test_v2 import LAG_S, K, VP, VS, _nodes, _q, _scene
from ttmt.radiation import node_operators
from ttmt.station_terms import OnsetDelay
from ttmt.traveltime import homogeneous
from ttmt.v2 import m1, m6
from ttmt.v2.loop import v2_loop


def _engine(STA):
    XN = _nodes()
    taus = np.round(np.arange(4.0, 9.0 + 1e-9, 0.1), 10)
    tt, _ = homogeneous(XN, STA, VP, VS)
    Gs, Qs = node_operators(XN, STA)
    return m1.OnsetEngine(XN, STA, tt, taus, OnsetDelay(np.full((len(STA), 2), -LAG_S)), Gs, Qs, rounding="nearest")


def test_joint_refit_cleans_an_overlapping_pair_better_than_the_greedy_loop():
    STA, srcs, F = _scene(2, 0.5, 3)                       # the scene where the greedy v2 loop leaves the most
    E = _engine(STA)
    Y = F[0][0] + F[1][0]
    q = _q([f[2] - LAG_S for f in F], 8)
    lag = np.full((8, 2), LAG_S / K.dt)
    _, Rg = v2_loop(E, Y, q, 2, lag, return_residual=True)
    atoms, Rj = m6.adcg_loop(E, Y, q, 2, lag, thr=None, return_residual=True)
    assert len(atoms) == 2
    assert np.linalg.norm(Rj) < 0.6 * np.linalg.norm(Rg)


def test_free_one_event_fit_is_at_least_as_good_as_the_tensor_fit():
    STA, srcs, F = _scene(1, 0.0, 11)
    E = _engine(STA)
    Y = F[0][0]
    q = _q([F[0][2] - LAG_S], 8)
    lag = np.full((8, 2), LAG_S / K.dt)
    atoms, R = m6.adcg_loop(E, Y, q, 1, lag, thr=None, return_residual=True)
    resp = m6._Resp(K)
    a = dict(atoms[0])
    assert m6.free_one_event_rss(Y, [], a, E.G, resp) <= float((R * R).sum()) + 1e-6


def test_resolvability_rejects_a_spurious_split_of_one_event():
    STA, srcs, F = _scene(1, 0.0, 11)
    E = _engine(STA)
    r = np.random.default_rng(2)
    Y = F[0][0] / np.abs(F[0][0]).max() * 30 + r.normal(size=F[0][0].shape)
    q = _q([F[0][2] - LAG_S], 8)
    lag = np.full((8, 2), LAG_S / K.dt)
    kept = m6.adcg_loop(E, Y, q, 2, lag, thr=3.29 ** 2)
    near = [a for a in kept[1:] if np.linalg.norm(a["x"][:2] - kept[0]["x"][:2]) <= 5.0
            and abs(a["tau"] - kept[0]["tau"]) <= 1.5]
    assert near == []


def test_joint_free_fit_is_at_least_as_good_as_the_legacy_one():
    """review T5.2 (b1-b3): the one-event alternative with joint P/S, wider times and the neighbour refitted is a
    larger class than the legacy one (bounded, sequential, neighbour fixed), so its RSS is not larger."""
    STA, srcs, F = _scene(2, 0.5, 13)
    E = _engine(STA)
    Y = F[0][0] + F[1][0]
    q = _q([F[0][2] - LAG_S, F[1][2] - LAG_S], 8)
    lag = np.full((8, 2), LAG_S / K.dt)
    atoms = m6.adcg_loop(E, Y, q, 2, lag, thr=None)
    assert len(atoms) == 2
    resp = m6._Resp(K)
    legacy = m6.free_one_event_rss(Y, [atoms[1]], atoms[0], E.G, resp)
    joint = m6.free_one_event_joint(Y, [atoms[1]], atoms[0], E.G, resp)
    assert joint <= legacy * (1 + 1e-6)
