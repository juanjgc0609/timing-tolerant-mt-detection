"""Unit and mechanism tests of the v2 components (ttmt.v2: M0, M1, M4, M5 loop). Noise-free or Gaussian synthetics,
fixed seeds; no data files."""
import numpy as np
import pytest

from ttmt.detector import ACAP, Engine, ProjectedEnergy, NOCAP
from ttmt.kernel import LEGACY as K, ricker
from ttmt.radiation import _rand_dc, node_operators, radiation_operator
from ttmt.station_terms import OnsetDelay, PeakDelay
from ttmt.traveltime import homogeneous
from ttmt.v2 import m0, m1
from ttmt.v2.loop import deflate_at, v2_loop
from ttmt.v2.m4 import joint_read

RNG = np.random.default_rng(20261108)
VP, VS = 6.0, 6.0 / 1.73
NT = 640
LAG_S = 0.25                      # onset-to-peak lag of the synthetic arrivals (s), as measured in R1


def _ring(n=8, r=20.0):
    a = np.arange(n) * 2 * np.pi / n + 0.2
    return np.stack([r * np.cos(a), r * np.sin(a), np.zeros(n)], 1)


def _nodes():
    ax = np.arange(-15.0, 15.0 + 1e-9, 1.5) + 0.37
    return np.array([[x, y, 8.0] for x in ax for y in ax])


def _mt(r):
    return _rand_dc(r).ravel()[[0, 4, 8, 1, 2, 5]] * 2000.0


def _source(STA, src, tau, m, stat):
    """matched filter of one source, its block vectors g (2, ns, 3), peak times (2, ns) in s."""
    g = (radiation_operator(src, STA) @ m).reshape(2, len(STA), 3)
    d = np.linalg.norm(STA - src, axis=1)
    tarr = np.stack([d / VP, d / VS]) + tau + stat
    t = np.arange(NT) * K.dt
    D = np.zeros((len(STA), 3, NT))
    for b in range(2):
        for i in range(len(STA)):
            D[i] += g[b, i][:, None] * ricker(t - tarr[b, i])[None]
    MF = np.array([[np.convolve(D[i, c], K.taps[::-1], mode="same") for c in range(3)] for i in range(len(STA))])
    return MF, g, tarr


def _q(onsets_list, ns, sd=0.04, amp=0.9):
    """synthetic onset evidence: a Gaussian bump of width sd (s) at every onset (s)."""
    t = np.arange(NT) * K.dt
    q = np.zeros((ns, 2, NT))
    for on in onsets_list:
        for b in range(2):
            for i in range(ns):
                q[i, b] = np.maximum(q[i, b], amp * np.exp(-0.5 * ((t - on[b, i]) / sd) ** 2))
    return q


# ------------------------------------------------------------------------------------------- M0
def test_tyler_covariance_recovers_anisotropic_noise_and_resists_transients():
    L = np.array([[4.0, 0.5, 0.0], [0.5, 2.0, 0.0], [0.0, 0.0, 1.0]])
    A = np.linalg.cholesky(L)
    X = (A @ RNG.normal(size=(3, 20000)))
    MF = X[None].copy()
    MF[0, :, ::100] *= 30.0                                       # 1 % transients, 30 x
    S = m0.noise_covariance(MF)[0]
    assert np.allclose(S, L, rtol=0.0, atol=0.15 * np.abs(L).max())
    clean = np.ones(MF.shape[-1], bool)
    clean[::100] = False
    Z = m0.whiten(MF[:, :, clean], m0.noise_covariance(MF))       # clean samples, whitened
    C = np.cov(Z[0])
    assert np.allclose(C, np.eye(3), atol=0.1)
    raw = np.cov(MF[0])                                            # the plain covariance is ruined by the transients
    assert np.abs(raw - L).max() > 5


def test_whitened_operators_keep_the_signal_in_the_subspace():
    STA = _ring(5)
    src = np.array([1.0, 2.0, 8.0])
    G = radiation_operator(src, STA)[None]
    L = np.stack([np.diag(RNG.uniform(0.5, 4.0, 3)) for _ in range(5)])
    m = RNG.normal(size=6)
    v = (G[0] @ m).reshape(2, 5, 3)
    W = m0.inv_sqrt(L)
    vw = np.einsum("ijk,bik->bij", W, v).reshape(-1)
    Gw, Qw = m0.whiten_operators(G, L)
    assert np.allclose(Gw[0] @ m, vw)
    assert np.isclose(np.linalg.norm(Qw[0].T @ vw), np.linalg.norm(vw))       # fully in the whitened subspace


def test_orientation_audit_recovers_a_rotation():
    n = RNG.normal(size=(40, 3))
    n /= np.linalg.norm(n, axis=1)[:, None]
    ang = np.deg2rad(30.0)
    pol = RNG.choice([-1.0, 1.0], 40)
    v = n * pol[:, None]
    MF = v.T[None, :, :]                                          # samples = events
    obs = m0.rotate_horizontals(MF, -ang)[0].T                    # misoriented by +30 degrees
    a, r, k = m0.orientation_audit(obs, n)
    assert abs(np.rad2deg(a) - 30.0) < 1e-6 and r > 0.999 and k == 40
    back = m0.rotate_horizontals(obs.T[None], a)[0].T
    assert np.allclose(back, v)


# ------------------------------------------------------------------------------------------- M1
def test_local_maxima_and_isotonic_calibration():
    q = np.array([0.0, 0.5, 0.5, 0.1, 0.3, 0.01, 0.02, 0.0])
    assert np.nonzero(m1.local_maxima(q[None]))[1].tolist() == [1, 4, 6]    # plateau counted once; floor 0.02 kept
    assert np.nonzero(m1.local_maxima(q[None], floor=0.05))[1].tolist() == [1, 4]
    cal = m1.IsotonicCalibration().fit([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], [0, 1, 0, 1, 1, 1])
    y = cal.transform(np.linspace(0, 1, 11))
    assert np.all(np.diff(y) >= 0) and y[-1] == 1.0
    assert np.isclose(cal.transform(0.25), 0.5)


def test_onset_engine_requires_onset_terms():
    STA = _ring(4)
    XN = _nodes()[:5]
    tt, _ = homogeneous(XN, STA, VP, VS)
    Gs, Qs = node_operators(XN, STA)
    with pytest.raises(TypeError):
        m1.OnsetEngine(XN, STA, tt, np.arange(4.0, 5.0, 0.1), PeakDelay.zeros(4), Gs, Qs)
    with pytest.raises(TypeError):
        Engine(XN, STA, tt, np.arange(4.0, 5.0, 0.1), OnsetDelay.zeros(4), Gs, Qs)
    E = m1.OnsetEngine(XN, STA, tt, np.arange(4.0, 5.0, 0.1), OnsetDelay.zeros(4), Gs, Qs)
    assert np.allclose(E.IDX + E.frac, E.tpred)


def test_onset_reader_takes_the_largest_maximum_with_subsample_refinement():
    q = np.zeros((1, 2, 100))
    t = np.arange(100)
    q[0, 0] = 0.8 * np.exp(-0.5 * ((t - 40.3) / 0.8) ** 2) + 0.5 * np.exp(-0.5 * ((t - 44) / 0.8) ** 2)
    s = m1.OnsetReader(4).read(q, np.array([[[42.0], [80.0]]]))
    assert abs(s[0, 0, 0] - 40.3) < 0.1 and np.isnan(s[0, 1, 0])


def test_onset_to_peak_lag_is_recovered():
    STA = _ring(6)
    MFs, ons = [], []
    for _ in range(5):
        src = np.r_[RNG.uniform(-8, 8, 2), 8.0]
        MF, g, tarr = _source(STA, src, 6.0 + RNG.uniform(0, 0.06), _mt(RNG), RNG.uniform(-0.1, 0.1, (2, 6)))
        MFs.append(MF)
        ons.append((tarr - LAG_S).T / K.dt)
    lag = m1.onset_to_peak_lag(MFs, ons)
    assert np.allclose(lag, LAG_S / K.dt, atol=0.15)


# ------------------------------------------------------------------------------------------- M1 mechanism
def test_onset_reading_survives_where_the_energy_peak_jumps():
    """two arrivals 0.3 s apart at one station, the later stronger: the energy argmax of a block predicted on the
    earlier one jumps to the later; the onset reader stays on the earlier onset and reads its vector, sign included."""
    nt, dt = NT, K.dt
    t = np.arange(nt) * dt
    a1, a2 = np.array([1.0, -0.5, 0.3]), np.array([-2.5, 1.5, -1.0])
    t1, t2 = 10.0, 10.3
    D = np.zeros((1, 3, nt))
    for a, tt in ((a1, t1), (a2, t2)):
        D[0] += a[:, None] * ricker(t - tt)[None]
    MF = np.array([[np.convolve(D[0, c], K.taps[::-1], mode="same") for c in range(3)]])
    from ttmt.detector import ENERGY_ARGMAX
    pred = np.array([[[int(round(t1 / dt))], [int(round(t1 / dt))]]])
    k_e = ENERGY_ARGMAX.read(MF, pred)[0, 0, 0]
    v_e = MF[0, :, k_e]
    q = _q([np.array([[t1 - LAG_S], [t1 + 5]]), np.array([[t2 - LAG_S], [t1 + 5]])], 1)
    s = m1.OnsetReader(4).read(q, np.array([[[(t1 - LAG_S) / dt], [(t1 + 5) / dt]]]))
    v_o = m1.read_vectors(MF, s + LAG_S / dt)[0, 0, 0]
    cos_e = v_e @ a1 / np.linalg.norm(v_e) / np.linalg.norm(a1)
    cos_o = v_o @ a1 / np.linalg.norm(v_o) / np.linalg.norm(a1)
    assert abs(k_e * dt - t1) > 0.2 and cos_e < 0                       # the energy peak jumped, wrong sign
    assert abs(s[0, 0, 0] * dt - (t1 - LAG_S)) < 0.02 and cos_o > 0.5   # the onset read keeps the right sign


# ------------------------------------------------------------------------------------------- M4 mechanism
def test_joint_read_rejects_a_stronger_wrong_maximum():
    """one station has a distractor maximum with HIGHER q inside its window, carrying a vector of the wrong sign.
    The isolated onset reader takes it; the joint reading takes the true one because the moment-tensor subspace
    penalises the sign."""
    STA = _ring(6)
    src = np.array([1.37, -2.63, 8.0])
    m = _mt(np.random.default_rng(5))
    MF, g, tarr = _source(STA, src, 6.0, m, np.zeros((2, 6)))
    on = tarr - LAG_S
    # distractor at station 0, P: an opposite-sign copy 0.4 s later (a later event), with a sharper onset
    t = np.arange(NT) * K.dt
    D = -1.2 * g[0, 0][:, None] * ricker(t - (tarr[0, 0] + 0.4))[None]
    MF[0] += np.array([np.convolve(D[c], K.taps[::-1], mode="same") for c in range(3)])
    q = _q([on], 6, amp=0.6)
    on2 = on.copy()
    on2[0, 0] += 0.4
    q[0, 0] = np.maximum(q[0, 0], 0.95 * np.exp(-0.5 * ((t - on2[0, 0]) / 0.04) ** 2))
    G = radiation_operator(src, STA)
    from ttmt.radiation import basis
    Q = basis(G)
    mu = on / K.dt
    lag = np.full((6, 2), LAG_S / K.dt)
    stat = ProjectedEnergy(NOCAP)
    iso = m1.OnsetReader(int(0.5 / K.dt)).read(q, mu[None, None])[0, 0]
    assert abs(iso[0, 0] - on2[0, 0] / K.dt) < 0.5                     # isolated read takes the distractor
    s, v, obj, d = joint_read(q, MF, mu, np.full((2, 6), 0.2 / K.dt), lag, Q, stat)
    assert abs(s[0, 0] - mu[0, 0]) < 0.5                               # joint read takes the true onset
    assert np.all(np.abs(s[np.isfinite(s)] - mu[np.isfinite(s)]) < 0.5)


# ------------------------------------------------------------------------------------------- loop and deflation
def test_deflate_at_removes_an_offgrid_arrival_read_at_its_time():
    STA = _ring(4)
    src = np.array([0.37, 0.37, 8.0])
    MF, g, tarr = _source(STA, src, 6.0 + 0.031, _mt(RNG), RNG.uniform(-0.1, 0.1, (2, 4)))
    # the fitted vector is the matched-filter value at the continuous peak time: g times the filtered pulse peak
    a0max = np.convolve(ricker(np.arange(-3, 3 + 1e-9, K.dt)), K.taps[::-1], "same").max()
    pred = np.zeros((2, 4, 3))
    for b in range(2):
        for i in range(4):
            pred[b, i] = g[b, i] * a0max
    R = deflate_at(MF.copy(), pred, tarr / K.dt)
    assert np.linalg.norm(R) < 1e-6 * np.linalg.norm(MF)


def _scene(n_src, dtau, seed):
    r = np.random.default_rng(seed)
    STA = _ring(8)
    srcs = []
    for k in range(n_src):
        while True:
            x = np.r_[r.uniform(-10, 10, 2), 8.0]
            if all(np.linalg.norm(x[:2] - s[0][:2]) >= 6.0 for s in srcs):
                break
        srcs.append((x, 6.0 + r.uniform(0, K.dt) + k * dtau, _mt(r), r.uniform(-0.15, 0.15, (2, 8))))
    F = [_source(STA, *s) for s in srcs]
    return STA, srcs, F


def test_v2_loop_finds_two_overlapping_sources_and_deflates_cleanly():
    STA, srcs, F = _scene(2, 0.5, 3)
    XN = _nodes()
    taus = np.round(np.arange(4.0, 10.0 + 1e-9, 0.1), 10) - LAG_S * 0       # origin grid
    tt, _ = homogeneous(XN, STA, VP, VS)
    Gs, Qs = node_operators(XN, STA)
    E = m1.OnsetEngine(XN, STA, tt, taus, OnsetDelay(np.full((8, 2), -LAG_S)), Gs, Qs, rounding="nearest")
    Y = F[0][0] + F[1][0]
    q = _q([f[2] - LAG_S for f in F], 8)
    lag = np.full((8, 2), LAG_S / K.dt)
    atoms, R = v2_loop(E, Y, q, 2, lag, reader=m1.OnsetReader(4), return_residual=True)
    hits = set()
    for a in atoms:
        for k, s in enumerate(srcs):
            if np.linalg.norm(a["x"][:2] - s[0][:2]) <= 5.0 and abs(a["tau"] - s[1]) <= 0.5:
                hits.add(k)
    assert hits == {0, 1}
    # the residual after both atoms is below the registered loop's (energy reading, grid deflation) on the same scene
    Eg = Engine(XN, STA, tt, taus, PeakDelay.zeros(8), Gs, Qs, rounding="nearest")
    from ttmt.detector import greedy_loop
    _, Rg = greedy_loop(Eg, Y, 2, ACAP(), return_residual=True)
    assert np.linalg.norm(R) < np.linalg.norm(Rg)


def test_v2_loop_with_joint_reading_runs_and_agrees_on_a_clean_scene():
    STA, srcs, F = _scene(1, 0.0, 11)
    XN = _nodes()
    taus = np.round(np.arange(4.0, 9.0 + 1e-9, 0.1), 10)
    tt, _ = homogeneous(XN, STA, VP, VS)
    Gs, Qs = node_operators(XN, STA)
    E = m1.OnsetEngine(XN, STA, tt, taus, OnsetDelay(np.full((8, 2), -LAG_S)), Gs, Qs, rounding="nearest")
    q = _q([F[0][2] - LAG_S], 8)
    lag = np.full((8, 2), LAG_S / K.dt)
    a_ind = v2_loop(E, F[0][0], q, 1, lag)
    a_joint = v2_loop(E, F[0][0], q, 1, lag, joint={}, n_joint=5, sig=0.15 / K.dt)
    for a in (a_ind[0], a_joint[0]):
        assert np.linalg.norm(a["x"][:2] - srcs[0][0][:2]) <= 3.0 and abs(a["tau"] - srcs[0][1]) <= 0.3
    assert a_joint[0]["obj"] is not None and isinstance(ACAP(), ACAP)


def test_refit_deflation_removes_an_offgrid_atom_even_with_a_wrong_amplitude():
    from ttmt.detector import deflate
    STA = _ring(3)
    MF, g, tarr = _source(STA, np.array([0.37, 0.37, 8.0]), 6.0 + 0.027, _mt(RNG), np.zeros((2, 3)))
    tp = np.floor(tarr / K.dt + 0.5).astype(int)
    pred = np.ones((2, 3, 3))                                        # wrong amplitude: only marks the blocks
    R = deflate(MF.copy(), pred, tp, K, "refit")
    assert np.linalg.norm(R) < 1e-3 * np.linalg.norm(MF)
