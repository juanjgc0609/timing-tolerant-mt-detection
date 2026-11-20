"""RAY1D card: the ray operator G_ray (ttmt.ray_operator)."""
import numpy as np

from ttmt.radiation import basis, radiation_operator
from ttmt.ray_operator import ray_directions, ray_operator
from ttmt.synth.media import LayeredModel
from ttmt.synth.realdata import pred_model

STA = np.c_[np.random.default_rng(3).uniform(-30, 30, (8, 2)), np.zeros(8)]


def test_homogeneous_limit_is_the_straight_ray():
    m = LayeredModel(thickness=(), vp=(6.0,), vs=(6.0 / 1.73,), rho=(2.7,))
    for src in ([1.0, -2.0, 8.0], [10.0, 5.0, 3.0]):
        assert np.max(np.abs(ray_operator(src, STA, m) - radiation_operator(src, STA))) < 1e-12


def test_directions_are_unit_and_steeper_at_the_surface():
    model = pred_model("ridgecrest_2019")
    dirs = ray_directions([0.0, 0.0, 8.0], STA, model)
    for ph in ("P", "S"):
        p, ns_, nr, h = dirs[ph]
        assert np.allclose(np.linalg.norm(ns_, axis=1), 1) and np.allclose(np.linalg.norm(nr, axis=1), 1)
        assert np.all(nr[:, 2] < 0)                                    # arriving upward
        # velocity grows with depth: the ray arrives steeper than it left (direct rays)
        assert np.all(np.abs(nr[:, 2]) >= np.abs(ns_[:, 2]) - 1e-12)


def test_p_is_along_the_arriving_ray_and_s_is_perpendicular():
    model = pred_model("hinet_2012")
    src = [2.0, 1.0, 12.0]
    G = ray_operator(src, STA, model).reshape(2, len(STA), 3, 6)
    _p, _ns, nrP, _h = ray_directions(src, STA, model)["P"]
    _p, _ns, nrS, _h = ray_directions(src, STA, model)["S"]
    m = np.random.default_rng(0).normal(size=6)
    uP, uS = G[0] @ m, G[1] @ m
    cosP = np.abs(np.sum(uP * nrP, 1)) / np.linalg.norm(uP, axis=1)
    assert np.allclose(cosP, 1.0)
    assert np.allclose(np.sum(uS * nrS, 1), 0.0, atol=1e-12)


def test_ray_events_are_in_the_span_of_g_ray_weighted():
    model = pred_model("ridgecrest_2019")
    src = [0.0, 0.0, 9.0]
    w = np.tile(np.repeat(np.geomspace(1, 10, len(STA)), 3), 2)
    Gw = w[:, None] * ray_operator(src, STA, model, free_surface=True)
    Q = basis(Gw)
    v = Gw @ np.random.default_rng(1).normal(size=6)
    y = v @ Q
    assert float(y @ y) / float(v @ v) > 1 - 1e-10
