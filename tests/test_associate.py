"""Unit tests of the comparator package (ttmt.associate): the shared random search, the append-only ledger, the
tuning objective under LOC+T-V, the picker peak extraction, the 1-D model handed to every program, and (when the
programs are installed) the GaMMA, PyOcto and REAL wrappers on a toy network."""
import json
import os

import numpy as np
import pytest

from ttmt import traveltime as TT
from ttmt.associate import ledger as L
from ttmt.associate import pickers as PK
from ttmt.associate import provenance as PV
from ttmt.associate import space as SP
from ttmt.associate import tune as TU
from ttmt.associate.base import Frame
from ttmt.associate.velocity import Velocity1D
from ttmt.atlas.schema import Detections, Picks, Stations, Truth

V = Velocity1D.hadley_kanamori()


def ring(n=8, r=20.0):
    a = np.arange(n) * 2 * np.pi / n
    return np.c_[r * np.cos(a), r * np.sin(a), np.zeros(n)]


def toy_truth(rng, n_pairs=15, n_iso=15, sep=0.3, span=1200.0, sta=None):
    """pairs at origin separation sep (1 km apart) and isolated events at least 6 s from any other event."""
    sta = ring() if sta is None else sta
    t1 = 20 + np.arange(n_pairs) * (span - 60) / (2 * max(n_pairs, 1))
    start = (t1[-1] + 30) if n_pairs else 20.0
    ti = start + np.arange(n_iso) * (span - 40 - start) / max(n_iso, 1)
    t = np.r_[t1, t1 + sep, ti]
    xy = rng.uniform(-8, 8, (len(t), 2))
    xy[n_pairs:2 * n_pairs] = xy[:n_pairs] + [1.0, 0.0]
    xyz = np.c_[xy, np.full(len(t), 8.0)]
    arr = np.stack([np.stack([t[e] + V.times(xyz[e], sta)[p] for p in (0, 1)], -1) for e in range(len(t))])
    m = len(t)
    pair = np.r_[np.arange(n_pairs), np.arange(n_pairs), -np.ones(n_iso, int)]
    return Truth(np.arange(m), t, xyz, np.full(m, 1.5), arr, np.full((m, len(sta)), 10.0), 0.0, span, pair=pair)


def toy_record(truth, name="toy", n_blocks=4, data=None):
    n = truth.arr_t.shape[1]
    return TU.TuningRecord(name, truth, 8, [np.arange(n)] * truth.n_events, (0.0, 0.0),
                           np.linspace(truth.t0, truth.t1, n_blocks + 1), "tune", data=data)


def toy_picks(truth, rng, sd=0.02, drop_second=False):
    st, ph, tt = [], [], []
    for e in range(truth.n_events):
        for i in range(truth.arr_t.shape[1]):
            for p in (0, 1):
                st.append(i), ph.append(p), tt.append(truth.arr_t[e, i, p] + rng.normal(0, sd))
    o = np.argsort(tt)
    return Picks(np.array(st)[o], np.array(ph)[o], np.array(tt)[o], np.full(len(tt), 0.9))


# ------------------------------------------------------------------------------------------------- search
def test_seeds_and_draws_are_fixed_and_shared():
    assert len(SP.SEEDS) == SP.N_TRIALS == 60 and len(set(SP.SEEDS)) == 60
    a, b = SP.trials("phasenet+gamma"), SP.trials("phasenet+gamma")
    assert a == b
    sp = SP.space_of("phasenet+gamma")
    assert set(sp) == set(SP.PHASENET) | set(SP.GAMMA)
    for _i, _s, cfg in a:
        SP.check(sp, cfg)
    # the same seed list for every method: the picker draws of two pipelines with the same picker coincide in name
    # order only through the seed, and the trial seeds are identical
    assert [s for _i, s, _c in SP.trials("eqnet")] == [s for _i, s, _c in a]
    with pytest.raises(ValueError):
        SP.check(sp, dict(a[0][2], as_min_picks=9))
    with pytest.raises(ValueError):
        SP.check(sp, {k: v for k, v in a[0][2].items() if k != "pk_threshold"})


def test_log_uniform_covers_its_range():
    x = np.array([SP.draw({"x": SP.LogUniform(0.05, 0.5)}, s)["x"] for s in range(4000)])
    assert x.min() >= 0.05 and x.max() <= 0.5
    assert abs(np.mean(x < np.sqrt(0.05 * 0.5)) - 0.5) < 0.04      # half below the geometric midpoint


# ------------------------------------------------------------------------------------------------- ledger
def test_ledger_is_append_only_and_freezes_the_best(tmp_path):
    p = tmp_path / "ledger.csv"
    for i, obj in enumerate([0.4, 0.7, float("nan"), 0.7]):
        L.append(dict(method="m", trial=i, seed=i, params=L.params_json({"a": i}), objective=obj, recall_T=obj,
                      recall_C=obj, status="ok" if not np.isnan(obj) else "error"), p)
    with pytest.raises(ValueError):
        L.append(dict(method="m", trial=1, params="{}"), p)
    h = L.sha256(p)
    assert L.done("m", p) == {0, 1, 2, 3}
    b = L.best("m", p)
    assert b["trial"] == "1"                                         # tie at 0.7: the lower trial
    rec = L.freeze("m", tmp_path / "frozen.json", p)
    assert rec["params"] == {"a": 1} and rec["trials_in_ledger"] == 4 and rec["trials_ok"] == 3
    assert rec["ledger_sha256"] == h
    s = L.sensitivity("m", p)
    assert "a" in s and np.isfinite(s["a"]["spearman"])


def test_ledger_rejects_foreign_columns(tmp_path):
    p = tmp_path / "x.csv"
    p.write_text("a,b\n1,2\n")
    with pytest.raises(ValueError):
        L.read(p)


# ------------------------------------------------------------------------------------------------- objective
def test_objective_and_partition_guard():
    rng = np.random.default_rng(0)
    tr = toy_truth(rng)
    rec = toy_record(tr)
    perfect = Detections(tr.t + rng.normal(0, 0.02, tr.n_events), tr.xyz, np.full(tr.n_events, 10.0))
    obj, rT, rC, nT, nC, hours = TU.objective([rec], [perfect])
    assert rT == 1.0 and rC == 1.0 and obj == 1.0 and nT == 30 and nC == 15 and abs(hours - 1200 / 3600) < 1e-12
    # a merger that reports one detection per pair loses the second event of every pair in T only
    keep = np.r_[np.ones(15, bool), np.zeros(15, bool), np.ones(15, bool)]
    merged = Detections(tr.t[keep], tr.xyz[keep], np.full(keep.sum(), 10.0))
    obj, rT, rC, *_ = TU.objective([rec], [merged])
    assert rT == 0.5 and rC == 1.0 and obj == 0.75
    with pytest.raises(ValueError):
        TU.TuningRecord("r", tr, 8, [np.arange(8)] * tr.n_events, (0, 0), [0, 1200], "read")


def test_tuning_split_is_the_generator_rule():
    from ttmt.synth.realdata import tune_read
    starts = np.arange(2000) * 300.0
    s = [TU.tuning_split("altotiberina_days", a) for a in starts]
    assert s == [tune_read("altotiberina_days", a) for a in starts]
    assert abs(np.mean([x == "tune" for x in s]) - 0.30) < 0.03


def test_run_is_resumable_and_records_failures(tmp_path):
    rng = np.random.default_rng(1)
    tr = toy_truth(rng)
    rec = toy_record(tr)
    p = tmp_path / "ledger.csv"

    def runner(cfg, record):
        if cfg["wf_threshold"] > 0.65:
            raise RuntimeError("boom")
        k = record.truth.t < 600 + 1000 * cfg["wf_threshold"]
        return Detections(record.truth.t[k], record.truth.xyz[k], np.ones(k.sum()))

    rows = TU.run("eqnet", runner, [rec], ledger=p, trials=SP.trials("eqnet")[:10], log=lambda *_: None)
    assert len(rows) == 10
    rows = TU.run("eqnet", runner, [rec], ledger=p, n_trials=12, log=lambda *_: None)
    assert [int(r["trial"]) for r in rows] == list(range(12))
    err = [r for r in rows if r["status"] == "error"]
    assert all(json.loads(r["params"])["wf_threshold"] > 0.65 for r in err)
    assert all(r["partition"] == rows[0]["partition"] and r["code"] for r in rows)
    ok = [r for r in rows if r["status"] == "ok"]
    top = max(float(r["objective"]) for r in ok)
    b = L.best("eqnet", p)
    assert float(b["objective"]) == top
    assert int(b["trial"]) == min(int(r["trial"]) for r in ok if float(r["objective"]) == top)


# ------------------------------------------------------------------------------------------------- pickers
def test_picks_from_probabilities_threshold_and_distance():
    dt, n = 0.01, 2000
    t = np.arange(n) * dt
    pr = np.zeros((2, n, 2), np.float32)
    bump = lambda c, h: h * np.exp(-0.5 * ((t - c) / 0.03) ** 2)
    pr[0, :, 0] = bump(5.0, 0.9) + bump(5.15, 0.6) + bump(12.0, 0.1)
    pr[1, :, 1] = bump(7.0, 0.5)
    pk = PK.picks_from_probabilities(pr, dt, threshold=0.3, peak_distance=0.1, t0=100.0)
    assert sorted(np.round(pk.t, 2).tolist()) == [105.0, 105.15, 107.0]
    assert pk.phase[np.argmin(np.abs(pk.t - 107.0))] == 1 and pk.sta[np.argmin(np.abs(pk.t - 107.0))] == 1
    pk = PK.picks_from_probabilities(pr, dt, threshold=0.3, peak_distance=0.5)
    assert sorted(np.round(pk.t, 2).tolist()) == [5.0, 7.0]       # the 0.15 s neighbour is suppressed
    pk = PK.picks_from_probabilities(pr, dt, threshold=0.05, peak_distance=0.5)
    assert 12.0 in np.round(pk.t, 2)


def test_probability_file_roundtrip(tmp_path):
    pr = np.random.default_rng(0).random((3, 50, 2)).astype(np.float32)
    PK.save_probabilities(tmp_path / "p.npz", pr, 0.01, 5.0, ["A", "B", "C"], "PhaseNet", "original")
    got, dt, t0, meta = PK.load_probabilities(tmp_path / "p.npz", ["C", "A"])
    assert np.array_equal(got, pr[[2, 0]]) and dt == 0.01 and t0 == 5.0 and meta["weights"] == "original"


# ------------------------------------------------------------------------------------------------- velocity
def test_velocity_model_forms():
    one = Velocity1D(((0.0, 6.0),), vpvs=1.75)
    assert np.allclose(one.mean_velocities(10.0), (6.0, 6.0 / 1.75))
    vp, _vs = V.mean_velocities(16.0)
    assert np.isclose(16.0 / vp, 5.5 / 5.5 + 10.5 / 6.3)
    z, p, s = V.profile(20.0)
    assert z[0] == 0 and np.all(np.diff(z) >= 0) and np.allclose(p / s, V.vpvs)
    with pytest.raises(ValueError):
        Velocity1D(((1.0, 6.0),))


def test_real_table_matches_first_arrival(tmp_path):
    path = V.real_ttdb(tmp_path / "ttdb.txt", 0.2, 10.0, dx_deg=0.02, dz=2.0)
    a = np.loadtxt(path, usecols=range(8))
    assert np.isfinite(a).all() and a.shape == (10 * 6, 8)
    row = a[(np.isclose(a[:, 0], 0.1)) & (a[:, 1] == 6.0)][0]
    x = 0.1 * 111.19
    assert np.isclose(row[2], TT.first_arrival(list(V.layers), 6.0, [x])[0], atol=1e-5)
    assert np.isclose(row[3], TT.first_arrival(list(V.layers), 6.0, [x], V.vpvs)[0], atol=1e-5)
    # REAL's linear correction t + dGCarc * prayp + ddep * hslow predicts the neighbouring entries
    nxt = a[(np.isclose(a[:, 0], 0.12)) & (a[:, 1] == 6.0)][0]
    assert abs(row[2] + 0.02 * row[4] - nxt[2]) < 0.01
    deeper = a[(np.isclose(a[:, 0], 0.1)) & (a[:, 1] == 8.0)][0]
    assert abs(row[2] + 2.0 * row[6] - deeper[2]) < 0.05


# ------------------------------------------------------------------------------------------------- provenance
def test_provenance_declares_known_overlaps():
    ov = PV.overlaps()
    assert ov["phasenet/instance"]["altotiberina_days"] == "period"
    assert ov["phasenet/scedc"]["ridgecrest_2019"] == "period"
    assert ov["phasenet/original"]["ridgecrest_2019"] == ""
    assert "phasenet/instance" in PV.table_markdown()


# ------------------------------------------------------------------------------------------------- wrappers
def _frame(tmp_path, **fixed):
    sta = ring()
    return Frame.around(Stations([f"S{i}" for i in range(8)], sta), V, zmax=20.0, workdir=str(tmp_path),
                        fixed=fixed)


def _recall(det, truth, tol=0.5):
    return np.mean([np.min(np.abs(det.t - t)) <= tol if len(det) else False for t in truth.t])


def test_gamma_wrapper_serial_and_deterministic(tmp_path):
    pytest.importorskip("gamma")
    from ttmt.associate import gamma_assoc as G
    rng = np.random.default_rng(3)
    tr = toy_truth(rng, n_pairs=0, n_iso=12, span=400.0)
    pk = toy_picks(tr, rng)
    fr = _frame(tmp_path, gamma={"ncpu": 8})
    c = G.config_for(fr, dict(as_min_picks=5, as_max_sigma11=1.0, as_method="BGMM"))
    assert c["ncpu"] == 1                                              # cannot be overridden
    cfg = dict(as_min_picks=5, as_max_sigma11=1.0, as_method="BGMM")
    d1, p1 = G.associate(fr, pk, cfg)
    d2, p2 = G.associate(fr, pk, cfg)
    assert np.array_equal(d1.t, d2.t) and np.array_equal(p1.det, p2.det)
    assert _recall(d1, tr) >= 0.9
    assert (p1.det >= 0).mean() > 0.8 and p1.det.max() < len(d1)
    assert np.allclose(d1.s, np.bincount(p1.det[p1.det >= 0], minlength=len(d1)))


def test_pyocto_wrapper(tmp_path):
    pytest.importorskip("pyocto")
    pytest.importorskip("pyrocko")                    # PyOcto's 1-D tables need it (comparator env, not core)
    from ttmt.associate import pyocto_assoc as PO
    rng = np.random.default_rng(4)
    tr = toy_truth(rng, n_pairs=0, n_iso=12, span=400.0)
    d, p = PO.associate(_frame(tmp_path), toy_picks(tr, rng),
                        dict(as_min_picks=6, as_tolerance=0.5, as_pick_match_tolerance=1.0))
    assert _recall(d, tr) >= 0.9 and (p.det >= 0).mean() > 0.8
    assert np.median([np.min(np.hypot(*(d.xyz[:, :2] - x[:2]).T)) for x in tr.xyz]) < 2.0


@pytest.mark.skipif(not (os.environ.get("TTMT_REAL_BIN") and os.path.exists(os.environ["TTMT_REAL_BIN"])),
                    reason="REAL binary not built (TTMT_REAL_BIN)")
def test_real_wrapper(tmp_path):
    from ttmt.associate import real_assoc as RE
    rng = np.random.default_rng(5)
    tr = toy_truth(rng, n_pairs=0, n_iso=10, span=400.0)
    d, p = RE.associate(_frame(tmp_path), toy_picks(tr, rng), dict(as_np0=3, as_nps0=6, as_tint=1.5))
    assert _recall(d, tr) >= 0.8 and (p.det >= 0).mean() > 0.7
    assert not [x for x in os.listdir(tmp_path) if x.startswith("real_") and os.path.isdir(tmp_path / x)]


def test_pipeline_trial_end_to_end(tmp_path):
    pytest.importorskip("pyocto")
    pytest.importorskip("pyrocko")
    from ttmt.associate import pipeline as PL
    rng = np.random.default_rng(6)
    tr = toy_truth(rng, n_pairs=4, n_iso=8, span=400.0)
    dt = 0.01
    n = int(tr.t1 / dt)
    t = np.arange(n) * dt
    prob = np.zeros((8, n, 2), np.float32)
    for e in range(tr.n_events):
        for i in range(8):
            for k in (0, 1):
                prob[i, :, k] = np.maximum(prob[i, :, k], 0.9 * np.exp(-0.5 * ((t - tr.arr_t[e, i, k]) / 0.05) ** 2))
    PK.save_probabilities(tmp_path / "orig.npz", prob, dt, 0.0, [f"S{i}" for i in range(8)], "PhaseNet", "original")
    fr = _frame(tmp_path)
    rec = toy_record(tr, data={"frame": fr, "prob": lambda w: tmp_path / "orig.npz"})
    cfg = dict(pk_peak_distance=0.1, pk_threshold=0.3, pk_weights="original", as_min_picks=6, as_tolerance=0.5,
               as_pick_match_tolerance=1.0)
    mr = PL.run_pipeline("phasenet+pyocto", cfg, fr, rec.data["prob"])
    assert mr.min_picks == 6 and len(mr.picks) >= 0.95 * tr.n_events * 16
    obj, rT, rC, *_ = TU.objective([rec], [mr.det])
    assert rC >= 0.9 and 0.0 <= rT <= 1.0


def test_manifest_roundtrip(tmp_path):
    rng = np.random.default_rng(7)
    tr = toy_truth(rng, n_pairs=2, n_iso=3, span=300.0)
    np.savez(tmp_path / "truth.npz", id=tr.id, t=tr.t, xyz=tr.xyz, mag=tr.mag, arr_t=tr.arr_t, snr=tr.snr, t0=tr.t0,
             t1=tr.t1, pair=tr.pair)
    man = {"records": [{"name": "r0", "partition": "tune", "truth": "truth.npz", "k": 8,
                        "stations": {"code": [f"S{i}" for i in range(8)], "xyz": ring().tolist()},
                        "subsets": [list(range(8))] * tr.n_events, "centroid": [0, 0],
                        "block_edges": [0, 150, 300], "prob": {"original": "p.npz"},
                        "frame": {"xlim": [-30, 30], "ylim": [-30, 30], "zlim": [0, 20]}}]}
    (tmp_path / "m.json").write_text(json.dumps(man))
    (rec,) = TU.load_manifest(tmp_path / "m.json", "phasenet+pyocto")
    assert rec.truth.n_events == tr.n_events and np.array_equal(rec.truth.pair, tr.pair)
    assert rec.data["frame"].zlim == (0, 20) and str(rec.data["prob"]("original")).endswith("p.npz")
    assert rec.block(160.0) == 1 and np.allclose(rec.block_hours, 150 / 3600)
    # Hi-net: no noise block hashes below 0.3, so the block at 0 s is "tune" only by the ADR-0009 fallback; the
    # generator's label is accepted as written
    from ttmt.synth.realdata import tune_read
    hinet = [0.0, 600.0, 1200.0, 1500.0]
    assert all(tune_read("hinet_2012", b) == "read" for b in hinet)
    man["records"][0].update(dataset="hinet_2012", noise_blocks=hinet[:1])
    (tmp_path / "m.json").write_text(json.dumps(man))
    assert len(TU.load_manifest(tmp_path / "m.json", "phasenet+pyocto")) == 1
    man["records"][0]["partition"] = "read"
    (tmp_path / "m.json").write_text(json.dumps(man))
    with pytest.raises(ValueError):
        TU.load_manifest(tmp_path / "m.json", "phasenet+pyocto")


def test_seeds_are_the_anchored_atlas_list():
    from ttmt.atlas.seeds import tuning_seeds
    stored = json.load(open(L.REPO_ROOT / "experiments" / "exploratory" / "ATLAS-0" / "seeds.json"))["tuning_seeds"]
    assert list(SP.SEEDS) == tuning_seeds() == stored
    i, s, cfg = SP.trials("phasenet+pyocto")[7]
    rng = np.random.default_rng(stored[7])
    sp = SP.space_of("phasenet+pyocto")
    assert cfg == {k: sp[k].draw(rng) for k in sorted(sp)} and s == stored[7]


def test_amendment_1_spaces_and_objective_guard(tmp_path):
    assert SP.REAL["as_tint"] == SP.LogUniform(0.1, 5.0)
    assert SP.PYOCTO["as_min_interevent_time"] == SP.LogUniform(0.1, 5.0)
    rng = np.random.default_rng(8)
    only_c = toy_record(toy_truth(rng, n_pairs=0, n_iso=6, span=300.0))
    with pytest.raises(ValueError):
        TU.objective([only_c], [Detections(only_c.truth.t, only_c.truth.xyz, np.ones(only_c.truth.n_events))])
    with pytest.raises(ValueError):                 # refused before any trial is spent
        TU.run("eqnet", lambda c, r: None, [only_c], ledger=tmp_path / "l.csv", log=lambda *_: None)
    assert not (tmp_path / "l.csv").exists()
