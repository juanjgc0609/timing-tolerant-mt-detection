"""The package against the frozen legacy functions, in the same process and on the same inputs: bit-for-bit equality
(no tolerance; both run on the same BLAS and libm). The legacy modules are imported from code/ through ttmt.legacy,
unchanged. Modules that set global state at import (r0_lib, r1_lib) are not imported here."""
import numpy as np
import pytest

from ttmt import legacy

legacy.install()


@pytest.fixture(scope="module")
def L():
    import c8_lib
    import c10_lib
    import froc_exact
    import station_corrections
    import task0_common
    import task1o_consistent
    import task1w_mech
    import task1x_mtglrt
    import task2d_methods
    import xc
    return dict(c8=c8_lib, c10=c10_lib, froc=froc_exact, sc=station_corrections, T0=task0_common,
                t1o=task1o_consistent, t1w=task1w_mech, X=task1x_mtglrt, D2=task2d_methods, xc=xc)


RNG = np.random.default_rng(7)


def _setup(L, ns=6):
    from ttmt.detector import Engine
    from ttmt.radiation import node_operators
    from ttmt.station_terms import PeakDelay
    from ttmt.traveltime import homogeneous
    a = np.linspace(0, 2 * np.pi, ns, endpoint=False) + 0.4
    STA = np.stack([np.linspace(9, 40, ns) * np.cos(a), np.linspace(9, 40, ns) * np.sin(a), np.zeros(ns)], axis=1)
    XN = np.asarray(L["xc"].make_grid(step=6.0, depths=(5.0, 9.0)))
    taus = np.round(np.arange(5.0, 8.01, 0.1), 10)
    corr = RNG.normal(0, 0.2, (ns, 2))
    Gs, Qst = node_operators(XN, STA)
    tt, _ = homogeneous(XN, STA, L["T0"].VP, L["T0"].VS)
    E = Engine(XN, STA, tt, taus, PeakDelay(corr), Gs, Qst, bomp=True)
    Eo = L["c10"].Eng(XN, STA, corr, taus, 6.0, G=Gs, Q=Qst, bomp=True)
    return STA, XN, taus, corr, E, Eo


def _window(L, STA, n_ev=3):
    """a synthetic window of the legacy generator (task0_common.generate) and its normalised matched filter."""
    T0 = L["T0"]
    rng = np.random.default_rng(11)
    srcs = [np.array([rng.uniform(-15, 15), rng.uniform(-15, 15), rng.uniform(5, 13)]) for _ in range(n_ev)]
    taus = [float(rng.uniform(5.5, 7.5)) for _ in range(n_ev)]
    D, _p, _s = T0.generate(STA, srcs, taus, n_clutter=10, rng=rng, noise_scale=0.05)
    MF, sd = L["t1o"].matched(D)
    return D, MF / sd, srcs, taus


def test_constants(L):
    from ttmt.kernel import LEGACY
    from ttmt.radiation import CP, CS
    assert np.array_equal(LEGACY.taps, L["t1o"]._KER)
    assert np.array_equal(LEGACY.ac, L["X"].AC) and LEGACY.ac_half == L["X"].AC_HALF
    assert CP == L["t1w"].CP and CS == L["t1w"].CS
    assert LEGACY.dt == L["T0"].DT


def test_operators_and_matched(L):
    from ttmt.kernel import LEGACY
    from ttmt.radiation import basis, radiation_operator
    STA, *_ = _setup(L)
    src = np.array([1.3, -4.2, 9.0])
    G = radiation_operator(src, STA)
    assert np.array_equal(G, L["X"].radiation_operator(src, STA))
    assert np.array_equal(basis(G), L["X"].basis(G))
    D, _MFn, _s, _t = _window(L, STA)
    a, b = LEGACY.matched(D), L["t1o"].matched(D)
    assert np.array_equal(a[0], b[0]) and a[1] == b[1]


def test_engine_tables_and_peak_table(L):
    from ttmt.detector import peak_table
    STA, XN, taus, corr, E, Eo = _setup(L)
    assert np.array_equal(E.IDX, Eo.IDX) and np.array_equal(E.frac, Eo.frac)
    _D, MFn, _s, _t = _window(L, STA)
    assert np.array_equal(peak_table(MFn), L["D2"].peak_table(MFn))


@pytest.mark.parametrize("capv", [6.0, None])
def test_mt_loop(L, capv):
    from ttmt.detector import mt_loop
    STA, XN, taus, corr, E, Eo = _setup(L)
    _D, MFn, _s, _t = _window(L, STA)
    a = mt_loop(E, MFn, 5, capv=capv)
    b = L["c10"].mt_loop(Eo, MFn, 5, capv=capv)
    assert len(a) == len(b) == 5
    for x, y in zip(a, b):
        assert x["s"] == y["s"] and x["tau"] == y["tau"] and np.array_equal(x["x"], y["x"])


def test_bomp_loop(L):
    from ttmt.detector import bomp_loop
    STA, XN, taus, corr, E, Eo = _setup(L)
    _D, MFn, _s, _t = _window(L, STA)
    a = bomp_loop(E, MFn, 4)
    b = L["c10"].bomp_loop(Eo, None, MFn, 4)
    assert [(x["s"], x["tau"]) for x in a] == [(y["s"], y["tau"]) for y in b]


def test_huber_loop(L):
    import types
    from ttmt.detector import huber_loop
    STA, XN, taus, corr, E, Eo = _setup(L)
    _D, MFn, _s, _t = _window(L, STA)
    dops = types.SimpleNamespace(ns=len(STA), garr=XN, IDX=E.IDX, Gs=E.G, Qst=E.Q, Gp=np.linalg.pinv(E.G),
                                 taus=np.asarray(taus, float))
    rec, locs = L["c8"].mtfw_general(dops, MFn, [np.array([1e6, 1e6, 0.0])], kmax=4, stat="huber_n")
    a = huber_loop(E, MFn, 4)
    assert [x["run"] for x in a] == [r[0] for r in rec]
    assert [list(x["x"]) + [x["tau"]] for x in a] == locs


def test_peak_delays_and_terms(L):
    from ttmt.station_terms import estimate_peak_terms, peak_delays
    STA, *_ = _setup(L)
    T0 = L["T0"]
    rng = np.random.default_rng(3)
    cal = []
    for _ in range(4):
        src = np.array([rng.uniform(-10, 10), rng.uniform(-10, 10), 7.0])
        D, _p, _s = T0.generate(STA, [src], [6.5], n_clutter=3, rng=rng, noise_scale=0.05)
        cal.append((D, src, 6.5))
    for search in (1.0, 2.0):
        a = peak_delays(STA, cal[0][0], cal[0][1], 6.5, T0.VP, T0.VS, search)
        b = L["sc"].peak_delays(STA, cal[0][0], cal[0][1], 6.5, search)
        assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])
    t = estimate_peak_terms(STA, cal, T0.VP, T0.VS)
    c, n = L["sc"].estimate_station_corrections(STA, cal)
    assert np.array_equal(t.values, c) and np.array_equal(t.info["nacc"], n)


def test_hit_rules_and_froc(L):
    from ttmt import scoring as S
    STA, *_ = _setup(L)
    truth = [np.array([0.0, 0.0, 7.0]), np.array([4.0, 3.0, 7.0])]
    taus = [6.0, 6.4]
    atoms = [{"x": np.array([RNG.uniform(-8, 8), RNG.uniform(-8, 8), 7.0]), "tau": float(RNG.uniform(5, 8)),
              "s": float(RNG.uniform(0, 100))} for _ in range(30)]
    for tol in (None, 0.5):
        assert S.label(atoms, truth, taus, tol) == L["c10"].label(atoms, truth, taus, tol)
    recs = [S.record(atoms[i:i + 5], S.label(atoms[i:i + 5], truth, taus, None)) for i in range(0, 30, 5)]
    assert recs == [L["c10"].record(atoms[i:i + 5], L["c10"].label(atoms[i:i + 5], truth, taus, None))
                    for i in range(0, 30, 5)]
    a, b = S.exact_HF(recs), L["froc"].exact_HF(recs)
    assert all(np.array_equal(x, y) for x, y in zip(a, b))
    assert S.recall_caps(recs, 12, (0.1, 0.25, 1.0)) == L["xc"].recall_caps(recs, 12, (0.1, 0.25, 1.0))
    grid = L["xc"].make_grid(step=6.0)
    field = RNG.uniform(0, 1, len(grid))
    assert S.greedy_record(grid, field, truth) == L["froc"].greedy_record(grid, field, truth)


def test_layered_times(L):
    import c18_tt
    from ttmt import traveltime as TT
    assert TT.JMA2001 == c18_tt.JMA2001 and TT.HK == c18_tt.HK
    src = np.array([1.0, 2.0, 7.3])
    STA, *_ = _setup(L)
    assert np.array_equal(TT.times(TT.JMA2001, 1.73, src, STA), c18_tt.times(c18_tt.JMA2001, 1.73, src, STA))
