"""Fast checks of the synthetic generator v2 (ttmt.synth). The long validation is tools/validate_synth.py."""
import ast
import json
from pathlib import Path

import numpy as np
import pytest

from ttmt.partition import PartitionError
from ttmt.synth import analytic, fk, inject, media, picklists, realdata, sets
from ttmt.synth import source as SRC

SYNTH = Path(__file__).resolve().parents[1] / "src" / "ttmt" / "synth"
FORWARD = ("media.py", "fk.py", "analytic.py", "source.py", "fklib.py")
BANNED = {"ttmt.kernel", "ttmt.radiation", "ttmt.traveltime", "ttmt.detector", "kernel", "radiation", "traveltime",
          "detector", "task0_common", "task1w_mech", "gen_realistic", "gen_b"}


def test_forward_modules_share_nothing_with_the_detector():
    for f in FORWARD:
        tree = ast.parse((SYNTH / f).read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                mod = (node.module or "")
                names = {mod, mod.split(".")[-1]} | {a.name for a in node.names}
                assert not names & BANNED, (f, mod, [a.name for a in node.names])
            elif isinstance(node, ast.Import):
                assert not {a.name for a in node.names} & BANNED, f
        assert "ricker" not in (SYNTH / f).read_text().lower().replace("not a ricker", "").replace(
            "no ricker", "").replace("none of them is a ricker", ""), f


def test_engine_matches_the_whole_space_solution():
    m = media.whole_space(6.0, 3.5, 2.7)
    cfg = fk.FKConfig(fs=10.0, nt=300, f_gauss=2.0, t_shift=2.0, L_factor=3.0)
    r, phi = np.array([6.0, 14.0]), np.array([0.4, -2.0])
    x = fk.greens(m, 5.0, r, phi, cfg, free_surface=False)
    t = np.arange(cfg.nt) / cfg.fs
    for i in range(2):
        rec = [r[i] * np.cos(phi[i]), r[i] * np.sin(phi[i]), 0.0]
        a = analytic.whole_space_velocity([0, 0, 5.0], rec, 6.0, 3.5, 2.7, t, 2.0, 2.0)
        for j in range(6):
            assert np.max(np.abs(x[i, j] - a[j])) < 0.01 * np.max(np.abs(a[j])), (i, j)


def test_ray_tracer_agrees_with_an_independent_implementation():
    from ttmt import traveltime as TT
    hk = [(0.0, 5.5), (5.5, 6.3), (16.0, 6.7), (32.0, 7.8)]
    m = media.LayeredModel((5.5, 10.5, 16.0), (5.5, 6.3, 6.7, 7.8), tuple(v / 1.73 for v in (5.5, 6.3, 6.7, 7.8)),
                           (2.5,) * 4)
    X = np.array([0.0, 3.0, 10.0, 30.0, 60.0, 100.0])
    for zs in (3.0, 8.0, 12.0):
        for ph, rr in (("P", 1.0), ("S", 1.73)):
            mine = media.first_arrival(m, zs, X, ph)[0]
            ref = TT.first_arrival(hk, zs, X, rr)
            assert np.max(np.abs(mine - ref)) < 1e-6


def test_model_split_and_velocities():
    m = media.PRESETS["crust_a"]
    m2, k = m.split(2.5)
    assert m2.n == m.n + 1 and abs(m2.tops[k] - 2.5) < 1e-12
    a, b = m.complex_velocities([1.0, 10.0])
    assert np.allclose(a[0].real, m.vp, rtol=1e-3) and np.all(a[1].real > a[0].real) and np.all(a.imag < 0)


def test_sources():
    M = SRC.dc_from_sdr(30, 60, 90)
    assert abs(np.trace(M)) < 1e-12 and abs(np.linalg.norm(M) - np.sqrt(2)) < 1e-12
    assert abs(np.linalg.det(M)) < 1e-12
    rng = np.random.default_rng(0)
    for fam in ("dc", "clvd", "iso", "generic"):
        T = SRC.mechanism(rng, fam)
        assert abs(np.linalg.norm(T) / np.sqrt(2) - 1) < 1e-9 and np.allclose(T, T.T)
    iso, clvd, dc = SRC.decompose(SRC.random_dc(rng))
    assert dc > 0.999
    f = np.array([0.0, 1.0, 10.0])
    s = SRC.stf_spectrum(f, "brune", 5.0)
    assert abs(s[0] - 1) < 1e-12 and abs(abs(s[2]) - 1 / (1 + 4)) < 1e-12
    assert abs(SRC.stf_spectrum(f, "boxcar", duration=0.2)[0] - 1) < 1e-12


def test_measurements():
    rng = np.random.default_rng(1)
    fs, n = 100.0, 4000
    w = rng.normal(size=(3, 3, n)) * np.exp(-((np.arange(n) - 1500) / 50.0) ** 2)
    cc, k = inject.max_cc(w, [np.roll(w, 37, axis=-1), rng.normal(size=w.shape)], fs)
    assert k == 0 and cc > 0.999
    g1 = np.array([[[1.0, 0, 0], [0, 1.0, 0]]])
    assert inject.kmec_over_e(g1, 2 * g1) < 1e-12
    assert inject.kmec_over_e(g1, g1[:, :, [1, 0, 2]]) > 0.1
    sig = np.zeros((2, 3, n))
    sig[:, 2, 1000] = 1.0
    noise = rng.normal(size=(2, 3, n))
    s = inject.card_snr(sig, noise, fs, np.array([10.0, np.nan]))
    assert np.isfinite(s[0]) and np.isnan(s[1])


def test_truth_roundtrip_and_reader(tmp_path, monkeypatch):
    from ttmt.atlas.schema import Truth
    t = Truth(np.array(["a", "b"]), [1.0, 2.0], np.zeros((2, 3)), [1.0, 2.0], np.full((2, 3, 2), np.nan),
              np.ones((2, 3)), 0.0, 10.0, pair=[0, 0], kmec=[0.1, np.nan], meta={"stations": ["x", "y", "z"]})
    d = sets.truth_to_dict(t)
    t2 = sets.truth_from_dict(json.loads(json.dumps(d)))
    assert np.array_equal(t2.pair, t.pair) and np.isnan(t2.kmec[1]) and t2.arr_t.shape == (2, 3, 2)
    v = inject.truth_contaminated(t, [1], 10.0)
    assert np.allclose(v.snr[:, 1], 0.1) and np.allclose(v.snr[:, 0], 1.0)
    man = {"realisations": [{"label": "L1", "partition": "tune"}], "files": {}, "pre_anchor": False}
    p = tmp_path / "m.json"
    p.write_text(json.dumps(man))
    r = sets.SetReader("x", "read", manifest=p)
    assert r.labels == []
    with pytest.raises(PartitionError):
        r._real("L1")
    with pytest.raises(PartitionError):
        sets.SetReader("x", "eval", manifest=p)
    q = tmp_path / "dev.json"                    # a set not seeded from the anchor serves no tuning or reading
    q.write_text(json.dumps({**man, "pre_anchor": True}))
    with pytest.raises(PartitionError):
        sets.SetReader("x", "tune", manifest=q)
    assert sets.SetReader("x", "tune", manifest=q, allow_pre_anchor=True).labels == ["L1"]


def test_contamination_and_subsets():
    xyz = np.array([[0, 0, 0], [1, 0, 0], [5, 5, 0], [-3, 1, 0], [10, 10, 0], [0, 2, 0]], float)
    sub = sets.subset_k(xyz, list("abcdef"), 4)
    assert len(sub) == 4 and 4 not in sub
    var = inject.contamination_variants(sub, np.random.default_rng(0))
    assert [(len(s), f) for _n, s, f in var] == [(1, 3.0), (1, 10.0), (2, 3.0), (2, 10.0)]
    assert all(set(s) <= set(sub) for _n, s, _f in var)


def test_picklists():
    STA = np.array([[0, 0, 0], [20, 0, 0], [0, 20, 0], [20, 20, 0]], float)
    tr, pk = picklists.scenario(STA, 3600.0, np.random.default_rng(2), rate_per_day=500.0, false_frac=1.0,
                                timing_sd=0.2)
    assert tr.arr_t.shape[1:] == (4, 2) and len(pk) > 0
    assert np.all(tr.arr_t[:, :, 1] > tr.arr_t[:, :, 0])
    err = (tr.arr_t_obs - tr.arr_t).ravel()
    assert 0.1 < err.std() < 0.35


def test_partition_rule_adr0009():
    """realdata.partition_blocks is the rule of ADR-0009 item 4, fallback both ways."""
    import hashlib

    def u(ds, a):
        h = hashlib.sha256(f"{ds}|{a:.3f}|tune-read-v1".encode()).digest()[:8]
        return int.from_bytes(h, "little") / 2 ** 64

    starts = [300.0 * i for i in range(12)]
    p = realdata.partition_blocks("ridgecrest_2019", starts)
    raw = {a: ("tune" if u("ridgecrest_2019", a) < 0.3 else "read") for a in starts}
    assert p == raw
    lo = [a for a in starts if raw[a] == "read"]
    p = realdata.partition_blocks("ridgecrest_2019", lo)
    assert [a for a in lo if p[a] == "tune"] == [min(lo, key=lambda a: u("ridgecrest_2019", a))]
    hi = [a for a in starts if raw[a] == "tune"]
    p = realdata.partition_blocks("ridgecrest_2019", hi)
    assert [a for a in hi if p[a] == "read"] == [max(hi, key=lambda a: u("ridgecrest_2019", a))]


def test_tune_read_split_hinet_fallback():
    """Hi-net: no noise block hashes below 0.3, so the noise block of smallest value goes to tuning (ADR-0009)."""
    import types
    rec = types.SimpleNamespace(dataset="hinet_2012", t0=0.0, t1=2000.0)
    sp = realdata.tune_read_split(rec)
    noise = [(bl, r) for bl, r in sp if r != "egf"]
    assert len(sp) == 6 and noise
    assert sum(r == "tune" for _bl, r in noise) >= 1 and sum(r == "read" for _bl, r in noise) >= 1
    u = {bl[0]: realdata.block_hash("hinet_2012", bl[0], "tune-read-v1") for bl, _r in noise}
    assert all(v >= 0.3 for v in u.values())
    assert [bl[0] for bl, r in noise if r == "tune"] == [min(u, key=u.get)]


def test_fklib_cache_resumes(tmp_path, monkeypatch):
    """the per-depth Green's function cache gives the same waveforms (to rounding: two direct runs already differ by
    ~1e-15 relative), and a second run computes nothing."""
    from ttmt.synth import fklib
    c = fklib.config({"fs": 20.0, "fmax": 5.0, "f_gauss": 4.0, "gf_len": 8.0, "medium_seed": None,
                      "depth_range": (4.0, 5.0), "depth_step": 1.0})
    STA = np.array([[0.0, 0.0, 0.0], [6.0, 2.0, 0.0], [-3.0, 5.0, 0.0]])
    rng = np.random.default_rng(3)
    ev = fklib.draw_events(STA, 3, rng, c)
    eff = fklib.network_effects(len(STA), rng, c)
    a = [w.wave.copy() for w in fklib.synthesise(ev, STA, c, eff)]
    b = [w.wave.copy() for w in fklib.synthesise(ev, STA, c, eff, cache=tmp_path)]
    n = len(list(tmp_path.glob("greens_*.npy")))
    assert n == len({float(e.xyz[2]) for e in ev}) and all(np.allclose(x, y, rtol=0, atol=1e-12 * abs(x).max()) for x, y in zip(a, b))
    calls = []
    d = [w.wave.copy() for w in fklib.synthesise(ev, STA, c, eff, cache=tmp_path, progress=lambda i, k: calls.append(i))]
    assert calls == [] and all(np.array_equal(x, y) for x, y in zip(b, d))
    monkeypatch.setattr(fklib, "_engine_digest", lambda: "another engine")   # a changed engine reuses nothing
    fklib.synthesise(ev, STA, c, eff, cache=tmp_path, progress=lambda i, k: calls.append(i))
    assert len(calls) == n
