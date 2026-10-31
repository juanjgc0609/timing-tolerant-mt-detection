"""Reproduction gate of ATLAS-0 (card, "Gate"): write experiments/exploratory/ATLAS-0/gate.json from the evidence.

    python gate.py write [--puente puente.json] [--ethz ethz.json] [--quakemigrate qm.json]

gate.json is {"targets": {name: {"passed": bool, "evidence": str}}}; the atlas reads data only when at least two
targets passed. Targets (card): puente_pyocto and puente_gamma (Puente et al. 2025, crustal scenario, Table S10,
parameters of Tables S3-S5, from repro_puente.py run), ethz_phasenet (SeisBench pick-benchmark, in-domain ETHZ) and
quakemigrate_example (Winder et al. 2026, stated event counts). A target passes when every published condition it was
run on is within 5 per cent (relative) of the published value. A Puente run on fewer than MIN_STATIONS stations is the
cloud stand-in, not the published network, and never passes. ethz and quakemigrate evidence files are
{"conditions": [{"name", "value", "published"}], "source": str}. A target without evidence is written as not passed.
"""
import argparse
import json
import sys
from pathlib import Path

GATE = Path(__file__).resolve().parents[1] / "ATLAS-0" / "gate.json"
TOL = 0.05
MIN_STATIONS = 100            # the SCSN list of the scenario has several hundred; the cloud stand-in has 23
NAMES = ("puente_pyocto", "puente_gamma", "ethz_phasenet", "quakemigrate_example")


def puente(path, tol=TOL):
    """{target: (passed, evidence)} from a repro_puente.py run file."""
    from repro_puente import S10
    res = json.load(open(path))
    out = {}
    for m in ("pyocto", "gamma"):
        rows = [r for r in res["runs"] if r["method"] == m and r.get("status") == "ok"
                and S10[m].get((r["n_events"], r["noise"])) is not None]
        rel = [(r["recall"] - S10[m][(r["n_events"], r["noise"])][1]) / S10[m][(r["n_events"], r["noise"])][1]
               for r in rows]
        ok = sum(abs(x) <= tol for x in rel)
        stand_in = res["stations"] < MIN_STATIONS
        passed = bool(rows) and ok == len(rows) and not stand_in
        ev = (f"{path}: {ok}/{len(rows)} conditions within {tol:.0%} of Table S10 recall on {res['stations']} stations"
              f" ({res['station_file']})" + ("; stand-in network, not the published one" if stand_in else ""))
        out[f"puente_{m}"] = (passed, ev)
    return out


def generic(path, tol=TOL):
    res = json.load(open(path))
    c = res["conditions"]
    ok = sum(abs(x["value"] - x["published"]) <= tol * abs(x["published"]) for x in c)
    return bool(c) and ok == len(c), f"{path}: {ok}/{len(c)} conditions within {tol:.0%} ({res.get('source', '')})"


def build(puente_path=None, ethz=None, quakemigrate=None):
    t = {n: {"passed": False, "evidence": "not run"} for n in NAMES}
    if puente_path:
        for n, (p, ev) in puente(puente_path).items():
            t[n] = {"passed": p, "evidence": ev}
    for n, path in (("ethz_phasenet", ethz), ("quakemigrate_example", quakemigrate)):
        if path:
            p, ev = generic(path)
            t[n] = {"passed": p, "evidence": ev}
    return {"targets": t, "passed": sum(v["passed"] for v in t.values()), "required": 2}


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    w = sub.add_parser("write")
    w.add_argument("--puente")
    w.add_argument("--ethz")
    w.add_argument("--quakemigrate")
    w.add_argument("--out", default=str(GATE))
    a = ap.parse_args(argv)
    g = build(a.puente, a.ethz, a.quakemigrate)
    Path(a.out).write_text(json.dumps(g, indent=1) + "\n")
    print(json.dumps(g, indent=1))


if __name__ == "__main__":
    sys.exit(main())
