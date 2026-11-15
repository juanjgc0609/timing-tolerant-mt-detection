"""Seal of the reserved confirmatory sequences (ttmt.registry) and the catalogue checks of ttmt.qc."""
import json

import numpy as np
import pytest

from ttmt import qc, registry
from ttmt.registry import SealedError, check_read, sealed_sha256


def _reservation(tmp_path, unseal=None):
    seqs = [{"id": "conf_x", "networks": ["XX"],
             "window": {"start": "2024-03-06T00:00:00Z", "end": "2024-04-05T00:00:00Z", "lat": 38.8, "lon": -122.8,
                        "radius_km": 60.0}}]
    doc = {"rule": "test", "sequences": seqs, "excluded": [], "unseal": unseal or {}}
    doc["sealed_sha256"] = sealed_sha256(doc)
    p = tmp_path / "reservation.json"
    p.write_text(json.dumps(doc))
    return p


def test_committed_reservation_matches_its_generator():
    import subprocess
    import sys
    r = subprocess.run([sys.executable, str(registry.config.REPO_ROOT / "tools" / "data" / "seal_reservation.py"),
                        "--check"], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    doc, seqs = registry.load_reservation()
    assert doc["sealed_sha256"] == sealed_sha256(doc) and seqs


def test_sealed_read_fails_inside_and_passes_outside(tmp_path):
    p = _reservation(tmp_path)
    with pytest.raises(SealedError):
        check_read("BG", 38.8, -122.7, "2024-03-10T00:00:00", "2024-03-10T01:00:00", reservation=p)
    # any network in the area is sealed, and an unknown station position counts as inside
    with pytest.raises(SealedError):
        check_read("ZZ", None, None, "2024-04-04T23:59:00", "2024-04-05T00:01:00", reservation=p)
    check_read("BG", 38.8, -122.7, "2024-04-05T00:00:00", "2024-04-06T00:00:00", reservation=p)   # after the window
    check_read("BG", 35.0, -118.0, "2024-03-10T00:00:00", "2024-03-11T00:00:00", reservation=p)   # outside the disc
    for purpose in ("fetch", "qc-metadata"):
        check_read("BG", 38.8, -122.7, "2024-03-10T00:00:00", "2024-03-11T00:00:00", purpose=purpose, reservation=p)


def test_tampered_reservation_is_refused(tmp_path):
    p = _reservation(tmp_path)
    doc = json.loads(p.read_text())
    doc["sequences"][0]["window"]["end"] = "2024-03-07T00:00:00Z"
    p.write_text(json.dumps(doc))
    with pytest.raises(SealedError, match="changed after sealing"):
        check_read("BG", 0.0, 0.0, 0, 1, reservation=p)


def test_rule_change_is_refused_and_unseal_is_outside_the_hash(tmp_path):
    p = _reservation(tmp_path)
    doc = json.loads(p.read_text())
    doc["unseal"] = {"conf_x": ["R9"]}
    p.write_text(json.dumps(doc))
    registry.load_reservation(p)                                       # an unseal amendment keeps the seal
    doc["rule"] = "another rule"
    p.write_text(json.dumps(doc))
    with pytest.raises(SealedError, match="changed after sealing"):
        registry.load_reservation(p)


def test_unseal_needs_listed_and_logged_registration(tmp_path, monkeypatch):
    log = tmp_path / "log.txt"
    log.write_text("... R9-confirmatory ...")
    monkeypatch.setenv("TTMT_WORK_ROOT", str(tmp_path / "work"))
    args = ("BG", 38.8, -122.7, "2024-03-10T00:00:00", "2024-03-10T01:00:00")
    monkeypatch.setenv("TTMT_UNSEAL", "conf_x:R9-confirmatory")
    with pytest.raises(SealedError):                                   # not listed in the reservation
        check_read(*args, reservation=_reservation(tmp_path), prereg_log=log)
    p = _reservation(tmp_path, {"conf_x": ["R9-confirmatory"]})
    check_read(*args, reservation=p, prereg_log=log)
    assert (tmp_path / "work" / "unseal_log.jsonl").read_text().count("R9-confirmatory") == 1
    log.write_text("R9-confirmatory-2 only")                           # a longer id does not match
    with pytest.raises(SealedError):
        check_read(*args, reservation=p, prereg_log=log)
    log.write_text("nothing")                                          # listed but not in the registration log
    with pytest.raises(SealedError):
        check_read(*args, reservation=p, prereg_log=log)


def _gr(b, mc, n, rng):
    m = mc - 0.05 + rng.exponential(1 / (b * np.log(10)), n)
    return np.round(np.concatenate([m, mc - 0.05 - rng.exponential(0.3, n // 3)]), 2)


def test_b_value_and_completeness_recover_a_known_law():
    rng = np.random.default_rng(1)
    for b in (0.8, 1.0, 1.2):
        m = _gr(b, 1.0, 5000, rng)
        bh, se, _ = qc.b_value(m, 1.0)
        assert abs(bh - b) < 3 * se + 0.02
        assert 0.8 <= qc.mc_gft(m)[90.0] <= 1.1
        assert 0.9 <= qc.mc_maxc(m) <= 1.3


def test_close_pairs_and_duplicates():
    t = np.array([0.0, 0.4, 5.0, 5.05, 20.0])
    assert qc.close_pairs(t)[0.5] == 2 and qc.close_pairs(t)[2.0] == 2
    x = np.array([0.0, 50.0, 0.0, 0.1, 0.0])
    assert qc.close_pairs(t, x, np.zeros(5), r_km=10.0)[1.0] == 1
    assert qc.duplicates(t, np.full(5, 38.0), np.full(5, -122.0)) == 1


def test_waveform_qc_finds_gaps_rate_changes_and_flat_lines(tmp_path):
    from obspy import Stream, Trace, UTCDateTime
    t0 = UTCDateTime("2020-01-01T00:00:00")
    rng = np.random.default_rng(0)
    a = Trace(rng.normal(size=6000).astype(np.int32), header={"network": "XX", "station": "A", "channel": "HHZ",
                                                               "sampling_rate": 100.0, "starttime": t0})
    b = Trace(rng.normal(size=3000).astype(np.int32), header={"network": "XX", "station": "A", "channel": "HHZ",
                                                               "sampling_rate": 50.0, "starttime": t0 + 120})
    flat = Trace(np.zeros(18000, np.int32), header={"network": "XX", "station": "B", "channel": "HHZ",
                                                     "sampling_rate": 100.0, "starttime": t0})
    p = tmp_path / "x.mseed"
    Stream([a, b, flat]).write(str(p), format="MSEED")
    out = qc.waveform_qc([p], t0, t0 + 180, samples=True)
    qa, qb = out["XX.A..HHZ"], out["XX.B..HHZ"]
    assert qa.n_gaps == 1 and abs(qa.max_gap_s - 60.0) < 0.05 and qa.sampling_rates == [50.0, 100.0]
    assert any("sampling rate" in i for i in qa.issues) and not qa.usable
    assert qb.coverage > 0.99 and qb.flatline_frac > 0.99 and not qb.usable
    head = qc.waveform_qc([p], t0, t0 + 180, samples=False)["XX.B..HHZ"]
    assert head.flatline_frac is None and head.coverage > 0.99


def _sealed_over(tmp_path, monkeypatch, start, end, lat, lon):
    seqs = [{"id": "conf_y", "networks": [],
             "window": {"start": start, "end": end, "lat": lat, "lon": lon, "radius_km": 60.0}}]
    doc = {"rule": "test", "sequences": seqs, "excluded": [], "unseal": {}}
    doc["sealed_sha256"] = sealed_sha256(doc)
    p = tmp_path / "reservation.json"
    p.write_text(json.dumps(doc))
    monkeypatch.setattr(registry, "RESERVATION", p)


def test_dataset_reader_refuses_a_sealed_window(tmp_path, monkeypatch):
    """the reader, not check_read directly: Ridgecrest.load_waveforms on a record inside a sealed window."""
    from obspy import Stream, Trace, UTCDateTime
    from ttmt import datasets
    root = tmp_path / "data"
    d = root / "skillsbench_eqassoc"
    d.mkdir(parents=True)
    t0 = UTCDateTime("2024-03-10T00:00:00")
    Stream([Trace(np.zeros(1000, np.int32), header={"network": "CI", "station": "AAA", "channel": "HH" + c,
                                                   "sampling_rate": 100.0, "starttime": t0}) for c in "ENZ"]) \
        .write(str(d / "wave.mseed"), format="MSEED")
    (d / "stations.csv").write_text("station,channel,latitude,longitude,response\n" +
                                    "".join(f"AAA,HH{c},38.8,-122.75,1.0\n" for c in "ENZ"))
    monkeypatch.setenv("TTMT_DATA_ROOT", str(root))
    ds = datasets.Ridgecrest(name="T", folder="ridgecrest_2019", t0_utc="2024-03-10T00:00:00", cal=(0.0, 5.0),
                             eval=(5.0, 10.0), subnets={}, m_target=1.0)
    sens = np.ones((1, 3))
    _sealed_over(tmp_path, monkeypatch, "2024-03-06T00:00:00Z", "2024-04-05T00:00:00Z", 38.8, -122.8)
    with pytest.raises(SealedError):
        ds.load_waveforms(["AAA"], sens, 10.0)
    _sealed_over(tmp_path, monkeypatch, "2024-05-01T00:00:00Z", "2024-05-02T00:00:00Z", 38.8, -122.8)
    raw, _ = ds.load_waveforms(["AAA"], sens, 10.0)                    # outside the window: reads
    assert raw.shape == (1, 3, 1000)


def test_waveform_qc_with_samples_refuses_a_sealed_window(tmp_path, monkeypatch):
    from obspy import Trace, UTCDateTime
    t0 = UTCDateTime("2024-03-10T00:00:00")
    p = tmp_path / "x.mseed"
    Trace(np.zeros(1000, np.int32), header={"network": "BG", "station": "B", "channel": "DPZ",
                                            "sampling_rate": 100.0, "starttime": t0}).write(str(p), format="MSEED")
    _sealed_over(tmp_path, monkeypatch, "2024-03-06T00:00:00Z", "2024-04-05T00:00:00Z", 38.8, -122.8)
    with pytest.raises(SealedError):
        qc.waveform_qc([p], t0, t0 + 10, samples=True)
    with pytest.raises(SealedError):
        qc.waveform_qc([p], t0, t0 + 10, samples=True, coords={"BG.B": (38.8, -122.7)})
    qc.waveform_qc([p], t0, t0 + 10, samples=True, coords={"BG.B": (30.0, -110.0)})   # outside the disc
    assert qc.waveform_qc([p], t0, t0 + 10, samples=False)["BG.B..DPZ"].coverage > 0.99  # headers only
