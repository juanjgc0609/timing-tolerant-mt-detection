"""Unit tests of the atlas harness (ttmt.atlas): cells, LOC+T-V scoring with the duplicator control, the rule of
ADR-0008 on a toy evaluation, pooled losses and the pick and event diagnostics."""
import numpy as np
import pytest

from ttmt import stats as S
from ttmt.atlas import cells as C
from ttmt.atlas import froc as F
from ttmt.atlas import instrument as I
from ttmt.atlas.schema import Detections, Picks, Truth

VP, VS = 6.0, 3.5


def ring(n=8, r=20.0):
    a = np.arange(n) * 2 * np.pi / n
    return np.c_[r * np.cos(a), r * np.sin(a), np.zeros(n)]


def arrivals(t, xyz, sta):
    d = np.linalg.norm(xyz[:, None, :] - sta[None, :, :], axis=-1)
    return np.stack([t[:, None] + d / VP, t[:, None] + d / VS], axis=-1)


def toy_record(rng, n_pairs=20, sep=0.3, span=600.0, sta=None):
    """n_pairs pairs at origin separation sep (same epicentre, 1 km apart) and n_pairs isolated events."""
    sta = ring() if sta is None else sta
    t1 = np.sort(rng.uniform(20, span - 40, n_pairs))
    t1 = t1 + np.arange(n_pairs) * 0.0
    t_iso = np.sort(rng.uniform(20, span - 40, n_pairs))
    t = np.r_[t1, t1 + sep, t_iso]
    xy = rng.uniform(-8, 8, (len(t), 2))
    xy[n_pairs:2 * n_pairs] = xy[:n_pairs] + [1.0, 0.0]
    xyz = np.c_[xy, np.full(len(t), 8.0)]
    keep = np.ones(len(t), bool)
    # drop isolated events that land within 5 s of another event so "iso" is really isolated
    for i in range(2 * n_pairs, len(t)):
        if (np.abs(np.delete(t, i) - t[i]) < 5.0).any():
            keep[i] = False
    t, xyz = t[keep], xyz[keep]
    m = len(t)
    pair = np.r_[np.arange(n_pairs), np.arange(n_pairs), -np.ones(m - 2 * n_pairs, int)]
    return Truth(np.arange(m), t, xyz, np.full(m, 1.5), arrivals(t, xyz, sta), np.full((m, len(sta)), 10.0),
                 0.0, span, pair=pair)


def perfect(truth, rng, n_false=5):
    t = np.r_[truth.t + rng.normal(0, 0.02, truth.n_events), rng.uniform(truth.t0, truth.t1, n_false)]
    xyz = np.r_[truth.xyz, np.c_[rng.uniform(-30, 30, (n_false, 2)), np.zeros(n_false)]]
    s = np.r_[rng.uniform(5, 10, truth.n_events), rng.uniform(0, 4, n_false)]
    return Detections(t, xyz, s)


def merging(truth, rng, n_false=5):
    """detects every event except the second of a pair closer than 1 s (one detection per merged pair)."""
    keep = np.ones(truth.n_events, bool)
    for e in range(truth.n_events):
        if truth.pair[e] >= 0:
            mate = np.flatnonzero((truth.pair == truth.pair[e]) & (np.arange(truth.n_events) != e))[0]
            if truth.t[e] > truth.t[mate] and truth.t[e] - truth.t[mate] < 1.0:
                keep[e] = False
    sub = Truth(truth.id[keep], truth.t[keep], truth.xyz[keep], truth.mag[keep], truth.arr_t[keep],
                truth.snr[keep], truth.t0, truth.t1)
    return perfect(sub, rng, n_false)


# ------------------------------------------------------------------------------------------------- cells
def test_bins_and_cell_roundtrip():
    assert C.bin_of([0.05, 0.1, 0.49, 0.5, 2.0, 3.0, np.inf, np.nan], C.SEP_EDGES).tolist() == [0, 1, 1, 2, 3, 4, 4, -1]
    for i in [0, 7, 119, C.N_CELLS - 1]:
        assert int(C.cell_index(*C.cell_axes(i))) == i
    assert C.N_CELLS == 240
    assert C.cell_kind(C.cell_index(1, 0, 0, 0)) == "target"
    assert C.cell_kind(C.cell_index(0, 0, 0, 0)) == "info"


def test_station_subset_availability_and_ties():
    xyz = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [5, 5, 0]], float)
    codes = np.array(["D", "C", "B", "A", "E"])
    sub = C.station_subset(xyz, codes, np.ones(5, bool), 2, net_centroid=[0, 0])
    assert sub.tolist() == [2, 3]                      # four stations at distance 1; codes A and B win the tie
    assert C.station_subset(xyz, codes, np.array([1, 0, 0, 0, 1], bool), 3) is None
    av = C.availability([100.0, 595.0], [[(0.0, 600.0)], [(0.0, 660.0)]])
    assert av.tolist() == [[True, True], [False, True]]


def test_separation_median_over_stations():
    sta = ring(4)
    t = np.array([10.0, 10.3])
    xyz = np.array([[0, 0, 8.0], [0, 0, 8.0]])
    sep = C.separation(arrivals(t, xyz, sta))
    assert np.allclose(sep, 0.3)
    assert F.regime_of([0.05, 0.3, 2.0, 5.0, np.nan]).tolist() == ["C", "T", "C", "iso", ""]


# ------------------------------------------------------------------------------------------------- scoring
def test_isolated_rule_is_loct():
    rng = np.random.default_rng(1)
    tr = toy_record(rng, n_pairs=0, sep=0.3)
    det = Detections(tr.t + 0.49, tr.xyz, np.arange(tr.n_events, dtype=float))
    assert (F.label(det, tr) >= 0).all()
    det2 = Detections(tr.t + 0.51, tr.xyz, np.arange(tr.n_events, dtype=float))
    assert (F.label(det2, tr) == F.FALSE).all()


def test_duplicator_gains_nothing_under_loctv():
    rng = np.random.default_rng(2)
    tr = toy_record(rng, n_pairs=30, sep=0.4)
    mer = merging(tr, rng, n_false=0)
    dup = F.duplicate(mer, 0.2)
    T = np.zeros(tr.n_events, bool)
    T[tr.pair >= 0] = True
    lab_m, lab_d = F.label(mer, tr), F.label(dup, tr)
    hit = lambda lab: np.isin(np.flatnonzero(T), lab[lab >= 0]).sum()
    assert hit(lab_d) <= hit(lab_m)
    # under plain LOC+T the duplicator would take both events of most pairs: the hole of CR-02
    dist = np.hypot(dup.xyz[:, None, 0] - tr.xyz[None, :, 0], dup.xyz[:, None, 1] - tr.xyz[None, :, 1])
    ok = (dist <= 5) & (np.abs(dup.t[:, None] - tr.t[None, :]) <= 0.5)
    from ttmt.scoring import one_to_one
    lab_plain = one_to_one(dup.s, dist, ok, order="score")
    assert hit(lab_plain) > hit(lab_m)


def test_weighted_operating_point_matches_stats():
    rng = np.random.default_rng(3)
    for _ in range(50):
        s = np.round(rng.uniform(0, 10, 40), 1)
        lab = np.where(rng.random(40) < 0.4, F.FALSE, rng.integers(0, 30, 40))
        assert F.operating_point(s, lab, np.ones(40), 2.0, 3.0) == S.operating_point(s, lab, 2.0, 3.0)


def _pool(seed, n_records=40, sep=0.4):
    rng = np.random.default_rng(seed)
    sta = ring()
    pool = F.Pool(block_hours=np.full(n_records, 600 / 3600), block_stratum=np.zeros(n_records, int))
    for b in range(n_records):
        tr = toy_record(rng, n_pairs=15, sep=sep, sta=sta)
        subsets = [np.arange(8)] * tr.n_events
        runs = {"perfect": perfect(tr, rng), "merging": merging(tr, rng)}
        F.score_record(pool, tr, 8, subsets, [0, 0], np.full(tr.n_events, b), lambda t, b=b: b, runs,
                       duplicator_of="merging")
    return pool


def test_gate_on_toy_evaluation():
    pool = _pool(4)
    assert (pool.ev_regime == "T").sum() > 0 and (pool.ev_regime == "iso").sum() > 0
    res = F.gate(pool, "perfect", ["merging"], duplicator_of=("merging", "dup:merging"))
    assert res["H1"]["rejected"] and res["H2"]["rejected"] and res["beats"]
    assert not res["control"]["broken"]
    rev = F.gate(pool, "merging", ["perfect"])
    assert not rev["H1"]["rejected"] and not rev["beats"]


def test_describe_and_pooled_loss():
    pool = _pool(5, n_records=30)
    d = F.describe(pool, ["perfect", "merging"], n_boot=50)
    assert d["perfect"]["regime"]["T"][0] == pytest.approx(1.0)
    assert d["merging"]["regime"]["T"][0] == pytest.approx(0.5, abs=0.05)
    lp = F.pooled_loss(pool, "perfect", n_boot=50)
    lm = F.pooled_loss(pool, "merging", n_boot=50)
    assert abs(lp["loss"]) < 1e-9
    assert lm["loss"] == pytest.approx(0.5, abs=0.05) and lm["lo"] > 0.3


# ------------------------------------------------------------------------------------------------- instrumentation
def _one_station_truth(times):
    t = np.asarray(times, float)
    arr = np.stack([t + 2.0, t + 4.0], axis=-1)[:, None, :]
    return Truth(np.arange(len(t)), t, np.zeros((len(t), 3)), np.ones(len(t)), arr, np.ones((len(t), 1)), 0, 100)


def test_pick_states():
    tr = _one_station_truth([10.0, 10.3, 30.0, 50.0])
    # event 0 P picked; event 1 P merged into it (one pick between them); event 2 P picked as S; event 3 lost
    pk = Picks(sta=[0, 0, 0], phase=[0, 1, 0], t=[12.05, 32.0, 70.0], p=[0.9, 0.5, 0.4])
    st, sep = I.arrival_states(tr, pk)
    assert st[:, 0, 0].tolist() == [I.PICKED, I.MERGED, I.WRONG_PHASE, I.LOST]
    assert sep[0, 0, 0] == pytest.approx(0.3)
    frac, tot = I.pick_table(st, sep)
    assert tot.sum() == 8


def test_event_merge_split_and_cause():
    sta = ring(8)
    t = np.array([10.0, 10.4, 40.0, 70.0])
    xyz = np.array([[0, 0, 8], [1, 0, 8], [5, 5, 8], [-5, 5, 8]], float)
    arr = arrivals(t, xyz, sta)
    tr = Truth(np.arange(4), t, xyz, np.ones(4), arr, np.ones((4, 8)), 0, 100)
    sta_i, ph_i, tt, det = [], [], [], []
    for s in range(8):                       # detection 0 gets the P picks of events 0 and 1 (merged)
        for e in (0, 1):
            sta_i.append(s); ph_i.append(0); tt.append(arr[e, s, 0]); det.append(0)
    for s in range(8):                       # event 2 split over detections 1 and 2
        sta_i.append(s); ph_i.append(0); tt.append(arr[2, s, 0]); det.append(1 if s < 4 else 2)
    for s in range(2):                       # event 3: only two picks, never associated
        sta_i.append(s); ph_i.append(0); tt.append(arr[3, s, 0]); det.append(-1)
    pk = Picks(sta_i, ph_i, tt, np.ones(len(tt)), det)
    labels = np.array([0, 2, F.FALSE])
    merged, split, cause = I.event_diagnostics(tr, pk, labels, 3, min_picks=4)
    assert merged.tolist() == [True, False, False]
    assert split.tolist() == [False, False, True, False]
    assert cause.tolist() == ["hit", "associator", "hit", "picker"]


def test_location_failures():
    tr = _one_station_truth([10.0, 30.0])
    det = Detections([10.2, 31.0], np.zeros((2, 3)), [2.0, 1.0])
    out = I.location_failures(det, tr)
    assert out["frac_fail"] == pytest.approx(0.5)


def test_duplicator_shift_per_pair_with_mixed_separations():
    """review of PR #12: one shift per record lets copies land on partners at other separations; the per-detection
    shift of ADR-0008 item 3 must gain no partner at any separation."""
    rng = np.random.default_rng(7)
    sta = ring()
    n = 200
    t1 = 20.0 + 10.0 * np.arange(n)
    sep = rng.uniform(0.1, 1.0, n)
    t = np.r_[t1, t1 + sep]
    xyz = np.tile(np.c_[rng.uniform(-8, 8, (n, 2)), np.full(n, 8.0)], (2, 1))
    tr = Truth(np.arange(2 * n), t, xyz, np.ones(2 * n), arrivals(t, xyz, sta), np.full((2 * n, 8), 10.0),
               0.0, t.max() + 20, pair=np.r_[np.arange(n), np.arange(n)])
    first = Detections(t1, xyz[:n], rng.uniform(1, 2, n))
    _, s, _, _ = C.assign(tr, 8, [np.arange(8)] * tr.n_events, [0, 0])
    assert np.allclose(s[:n], sep)
    dup = F.duplicate(first, F.duplicator_shift(first, tr, s))
    lab = F.label(dup, tr)
    assert not (lab >= n).any()                          # no partner hit
    # a single record-wide shift (the defect) does take partners
    bad = F.duplicate(first, 0.5 * float(np.median(sep)))
    assert (F.label(bad, tr) >= n).sum() > 20


def test_truth_precision_excludes_close_pairs():
    tr = _one_station_truth([10.0, 10.3, 30.0, 30.9, 60.0])
    keep = F.truth_precision_targets(tr, sigma_t=0.05)            # 9 sigma = 0.45 s
    assert keep.tolist() == [False, False, True, True, True]
    assert F.truth_precision_targets(tr, None).all()


def test_anchor_seeds():
    import hashlib
    from ttmt.atlas import seeds as Z
    h = hashlib.sha256(b"7bafb7abd308d0376ff1176d3d7abeb989aa0564" b"2026-10-24T22:37:39Z").digest()
    assert Z.master_seed() == int.from_bytes(h[:8], "big") == 0x8e1a8d5630dc0d52
    t = Z.tuning_seeds()
    assert len(t) == 60 and len(set(t)) == 60
    a, b = Z.rng("inject", 3).random(4), Z.rng("inject", 3).random(4)
    assert np.array_equal(a, b) and not np.array_equal(a, Z.rng("inject", 4).random(4))
    with pytest.raises(KeyError):
        Z.rng("unregistered", 0)
    from pathlib import Path
    assert Z.check_file(Path(__file__).parents[1] / "experiments/exploratory/ATLAS-0/seeds.json")


def test_intervals_use_anchored_streams():
    from ttmt.atlas import seeds as Z
    pool = _pool(4)
    a = F.pooled_loss(pool, "merging", n_boot=50)
    b = F.pooled_loss(pool, "merging", n_boot=50, seed=None)
    c = F.pooled_loss(pool, "merging", n_boot=50, seed=12345)
    assert a == b
    assert Z.rng("bootstrap", 1).integers(0, 10**9) == Z.rng("bootstrap", 1).integers(0, 10**9)
    assert set(a) == set(c)
