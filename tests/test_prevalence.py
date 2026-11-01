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
