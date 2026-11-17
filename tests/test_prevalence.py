"""CR-04 prevalence (ttmt.atlas.prevalence): windowed separation equals the card's, and the table on a toy catalogue."""
import numpy as np

from ttmt.atlas import cells as C
from ttmt.atlas import prevalence as P

VP, VS = 6.0, 3.5


def tt(src, sta):
    d = np.sqrt(((np.asarray(sta, float) - np.asarray(src, float)) ** 2).sum(axis=1))
    return np.stack([d / VP, d / VS])


def ring(n=12, r=20.0):
    a = np.arange(n) * 2 * np.pi / n
    return np.c_[r * np.cos(a), r * np.sin(a), np.zeros(n)]


def catalogue(seed=1, m=60):
    rng = np.random.default_rng(seed)
    t = np.sort(rng.uniform(0, 7200, m))
    t[10] = t[9] + 0.4                       # a pair at 0.4 s: in T
    xyz = np.c_[rng.uniform(-5, 5, (m, 2)), rng.uniform(5, 10, m)]
    return t, xyz


def test_windowed_separation_matches_cells():
    t, xyz = catalogue()
    sta = ring()
    arr = P.arrivals(t, xyz, sta, tt)
    st = [np.arange(len(sta))] * len(t)
    full = C.separation(arr, st)
    win = P.separation_windowed(arr, t, st, window=60.0)
    assert np.array_equal(C.bin_of(full, C.SEP_EDGES), C.bin_of(win, C.SEP_EDGES))


def test_table_counts_the_pair_and_shares_sum_to_one():
    t, xyz = catalogue()
    sta = ring()
    codes = [f"S{i:02d}" for i in range(len(sta))]
    avail = np.ones((len(t), len(sta)), bool)
    avail[:, 8:] = False                     # only 8 stations available: k = 12 contributes nothing
    out = P.table(t, xyz, sta, codes, avail, tt, n_boot=200, seed=3)
    assert out["k"][12]["n_counted"] == 0 and out["k"][12]["shares"] is None
    k4 = out["k"][4]
    assert k4["n_counted"] == len(t)
    assert np.isclose(np.sum(k4["shares"]), 1.0)
    assert k4["Y"] >= 2 / len(t)             # both events of the pair
    lo, hi = k4["Y_ci95"]
    assert lo <= k4["Y"] <= hi
    assert 0.0 <= k4["Y_poisson"] <= 1.0
    assert np.isclose(out["dead_time_s"][0], np.diff(np.sort(t)).min())


def test_surface_source_has_finite_times():
    """a source whose depth rounds to 0 km or above the surface (QTM lists depths down to -0.74 km) gets the surface
    direct wave, not infinity, and matches a 0.1 km source within 0.05 s (CR-04 QTM run, 2026-11-15)."""
    import numpy as np
    from ttmt import traveltime as TT
    st = TT.SourceTimes(TT.HK)
    sta = np.array([[1.0, 0, 0], [5.0, 0, 0], [10.0, 0, 0], [20.0, 0, 0], [50.0, 0, 0]])
    ref = st(np.array([0.0, 0.0, 0.1]), sta)
    for z in (-0.74, 0.0, 0.04):
        t = st(np.array([0.0, 0.0, z]), sta)
        assert np.isfinite(t).all()
        assert np.abs(t - ref).max() < 0.05


def test_surface_source_frozen_path_unchanged():
    """the frozen R1 model (traveltime.frozen, E65) keeps its registered behaviour: a source at 0 km has no direct
    wave there, so the surface fix does not change R1's frozen lists."""
    import numpy as np
    from ttmt import traveltime as TT
    X = np.array([1.0, 5.0, 10.0])
    t_new = TT.first_arrival(TT.HK, 0.0, X)
    t_old = TT.first_arrival(TT.frozen(TT.HK), 0.0, X)
    assert np.isfinite(t_new).all()
    assert not np.isfinite(t_old).all()
