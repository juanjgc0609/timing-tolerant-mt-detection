"""Regressions of the package against the frozen outputs (ADR-0002 item 4; ROADMAP phase 1 gate).

(a) 2M rescored from its stored part files reproduces results/task2M/analyze_2m_out.txt (no raw data needed).
(b) R0 and R1 scoring stages from the stored detections reproduce their registered analysis files (need the catalogues
    of the raw data folders).
(c) R1 re-run from the raw waveforms through the package (stored picks as input) reproduces the stored P1 lists, the
    hybrid lists and the registered analysis (slow; raw data).
(d) 2F from regenerated records: pending (its raw records were not stored, E49; they are being regenerated elsewhere).
Comparison: ``--match exact`` on Linux, the declared tolerance of ttmt.testing elsewhere.
"""
import json
import os

import numpy as np
import pytest

from ttmt import config
from ttmt.testing import compare_rows, rows_ok

RES = config.RESULTS_DIR


def _text(p):
    with open(p, newline="") as fh:
        return fh.read()


# ------------------------------------------------------------------------------------------- (a) 2M
@pytest.mark.slow
def test_2m_analysis_from_stored_parts():
    from ttmt.analysis.t2m import analyze, registered_specs
    lines, res = analyze(registered_specs(RES), 2000)
    assert "\n".join(lines) + "\n" == _text(RES / "task2M" / "analyze_2m_out.txt")
    c = res["2718"][("DENSE", "LOC", "MT@S0", "BOMP@S0")]
    assert round(c["point"], 3) == 0.228 and round(c["lo"], 3) == 0.148


# ------------------------------------------------------------------------------------------- (b) R0, R1 scoring
@pytest.mark.rawdata("ridgecrest_2019")
def test_r0_scoring_from_stored_detections(tmp_path):
    from ttmt.analysis.r0 import analyze
    txt, js = analyze("eval", str(RES / "taskR0" / "run_close"), 2000)
    assert txt + "\n" == _text(RES / "taskR0" / "run_close" / "analysis_eval.txt")
    ref = json.load(open(RES / "taskR0" / "run_close" / "analysis_eval.json"))
    assert json.loads(json.dumps(js)) == ref


@pytest.mark.rawdata("ridgecrest_2019")
def test_r0_closure_reproduces_registered_m1_m3():
    """the registered R0 recalls of M1 (all subnets) and M3 (N8) are those of the stored re-run (README of R0)."""
    reg = _text(RES / "taskR0" / "analysis_eval.txt")
    new = _text(RES / "taskR0" / "run_close" / "analysis_eval.txt")

    def rows(t):
        out = {}
        for line in t.splitlines():
            if line.startswith("  LOC") and "/h" in line:
                key = tuple(line.split()[:2])
                out[key] = {f.split(":")[0]: f.split(":")[1].split(",") for f in line.split()[2:] if ":" in f}
        return out
    r, n = rows(reg), rows(new)
    assert r and set(n) <= set(r)
    for k in n:
        assert r[k]["M1"] == n[k]["M1"], k
        assert r[k]["M3"][0] == n[k]["M3"][0], k                 # N8 is the first subnet


@pytest.mark.rawdata("hinet_2012")
def test_r1_scoring_from_stored_detections():
    from ttmt.analysis.r1 import analyze
    d = RES / "taskR1" / "run"
    txt = analyze("eval", json.load(open(d / "hyb_eval.json")), json.load(open(d / "hyb_cal.json")))
    assert txt == _text(d / "analysis_eval.txt")
    assert "R1-a: HYB - BASE, LOC+T 30/h: 23 vs 22 of 54" in txt


# ------------------------------------------------------------------------------------------- (c) R1 end to end
@pytest.mark.slow
@pytest.mark.rawdata("hinet_2012")
def test_r1_end_to_end_from_raw_data(match_mode):
    """every stage is compared and reported (work/r1_e2e/report.json) before the assertions."""
    from ttmt import r1
    from ttmt.analysis.r1 import analyze
    d = RES / "taskR1" / "run"
    out = config.work_root() / "r1_e2e"
    out.mkdir(parents=True, exist_ok=True)
    pk = json.load(open(d / "picks_all.json"))
    jobs = int(os.environ.get("TTMT_JOBS", os.cpu_count() or 1))
    report, ok, hyb = {"mode": match_mode}, {}, {}
    for seg in ("cal", "eval"):
        p1 = r1.run_p1(seg, pk, jobs=jobs, log=lambda s: None)
        json.dump(p1, open(out / f"p1_{seg}.json", "w"))
        ref = json.load(open(d / f"p1_{seg}.json"))
        dif = compare_rows(p1["atoms"]["NOCAP"], ref["atoms"]["NOCAP"], match_mode, label=f"P1 {seg}")
        dif["max_corr_diff"] = float(np.abs(np.array(p1["corr_layered"]) - np.array(ref["corr_layered"])).max())
        report[f"P1run_{seg}"], ok[f"P1run_{seg}"] = dif, rows_ok(dif, match_mode)
        h = r1.run_hybrid(seg, pk, p1)
        json.dump(h, open(out / f"hyb_{seg}.json", "w"))
        href = json.load(open(d / f"hyb_{seg}.json"))
        hj = json.loads(json.dumps(h))
        ok[f"base_{seg}"] = hj["base"] == href["base"] and hj["selection"] == href["selection"]
        for key in ("BASE", "P1", "HYB"):
            dif = compare_rows(h[key], href[key], match_mode, label=f"{key} {seg}")
            report[f"{key}_{seg}"], ok[f"{key}_{seg}"] = dif, rows_ok(dif, match_mode)
        hyb[seg] = h
    txt = analyze("eval", hyb["eval"], hyb["cal"])
    (out / "analysis_eval.txt").write_text(txt)
    ok["analysis_text"] = txt == _text(d / "analysis_eval.txt")
    report["ok"] = ok
    json.dump(report, open(out / "report.json", "w"), indent=1)
    print(json.dumps(report, indent=1))
    assert all(ok.values()), report


# ------------------------------------------------------------------------------------------- (d) 2F
@pytest.mark.skip(reason="2F records not stored (E49); regenerated on Linux by the duplicate-scoring work (PR #4). "
                         "The test reads results/task2f/records/ once that PR is merged.")
def test_2f_from_regenerated_records():
    path = RES / "task2f" / "records"
    assert path.is_dir()
