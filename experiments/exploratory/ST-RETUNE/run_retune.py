"""The retuning of ATLAS-0 amendment 6 (station terms), in the order of the card.

    python run_retune.py check <tune>            # no trial: environment, data, terms file, ledger, tag
    python run_retune.py smoke <tune>            # the first 3 trials of each method (they count in the 60)
    python run_retune.py full <tune> [--methods phasenet+pyocto]
    python run_retune.py read                    # freeze, reading rule, pivot and the paired comparison

<tune> holds ss_ridgecrest_2019_egf_v1/tune_phasenet.json and ss_hinet_2012_egf_v1/tune_phasenet.json (written by
COMPARATORS/tune_cli.py records). "smoke", "full" and "read" refuse to start unless the annotated tag
reg/ATLAS-0/amendment-6 exists and is an ancestor of HEAD, so no trial with terms can precede the anchoring. "check"
runs no association and writes nothing.
"""
import argparse
import json
import pathlib
import subprocess
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "experiments" / "exploratory" / "COMPARATORS"))

import tune_cli  # noqa: E402
from ttmt.associate import ledger as L  # noqa: E402
from ttmt.associate import pipeline as PL  # noqa: E402
from ttmt.associate import space as SP  # noqa: E402
from ttmt.associate import tune as TU  # noqa: E402

TAG = "reg/ATLAS-0/amendment-6"
TERMS = HERE / "terms_ridgecrest_2019.json"
TERMS_SHA = "1ba1719e07a6193ca490009a7d488efba8469a881d7007bdd91a50360ab898bc"
PARTITION = "b60d99b3d0bbbbfa"
SETS = ("ss_ridgecrest_2019_egf_v1", "ss_hinet_2012_egf_v1")
FROZEN = REPO / "experiments" / "exploratory" / "COMPARATORS" / "frozen"
# method: (ledger of its frozen search without terms, that frozen file)
METHODS = {"phasenet+pyocto": (L.LEDGER, FROZEN / "phasenet+pyocto.json"),
           "phasenet+gamma": (L.LEDGER_A4, FROZEN / "phasenet+gamma_a4.json")}
GAIN = 0.01                                                 # amendment 6, item 5 (pivot), and the reading of H
RESULT = REPO / "results" / "st_retune" / "retune_a6.json"


def git(*a):
    return subprocess.run(("git", "-C", str(REPO)) + a, capture_output=True, text=True)


def anchored():
    """the amendment tag is annotated and an ancestor of HEAD."""
    if git("cat-file", "-t", TAG).stdout.strip() != "tag":
        return False
    return git("merge-base", "--is-ancestor", TAG + "^{commit}", "HEAD").returncode == 0


def same_trial(row, seed, params):
    """the ledger row (seed, params) is this trial; floats may differ in the last digit between machines (seen: the
    Windows ledger against Linux, 1 part in 1e16)."""
    if row is None or row[0] != seed or set(row[1]) != set(params):
        return False
    return all(np.isclose(row[1][k], v, rtol=1e-12, atol=0.0) if isinstance(v, float) else row[1][k] == v
               for k, v in params.items())


def manifests(tune):
    return [str(pathlib.Path(tune) / s / "tune_phasenet.json") for s in SETS]


def check(tune, methods):
    """every condition of the run that can be verified without a trial; returns the list of failures."""
    bad = []
    if L.sha256(TERMS) != TERMS_SHA:
        bad.append("terms file differs from the hash of amendment 6, item 1")
    for m in methods:
        if not PL.STATUS.get(m, "").startswith("ready"):
            bad.append(f"{m}: {PL.STATUS.get(m, 'unknown method')}")
            continue
        try:
            PL.runner(m)
        except Exception as e:                                     # noqa: BLE001
            bad.append(f"{m}: runner not importable ({type(e).__name__}: {e})")
        recs = []
        for path in manifests(tune):
            for r in TU.load_manifest(path, m):
                if "ridgecrest" in path:
                    r.data["terms"] = tune_cli.load_terms(TERMS, r.data["frame"])
                recs.append(r)
        TU.check_records(recs)
        with_terms = sum("terms" in r.data for r in recs)
        part = TU.partition_hash(recs)
        print(f"{m}: {len(recs)} records ({with_terms} with terms), partition {part}, "
              f"{len(L.done(m, L.LEDGER_A6)) if L.LEDGER_A6.exists() else 0} trials already in the a6 ledger")
        if part != PARTITION:
            bad.append(f"{m}: partition {part}, registered {PARTITION}")
        if not with_terms or with_terms == len(recs):
            bad.append(f"{m}: terms must be on the Ridgecrest records only")
        led = {int(r["trial"]): (int(r["seed"]), json.loads(r["params"])) for r in L.read(METHODS[m][0])
               if r["method"] == m}
        if any(not same_trial(led.get(i), s, json.loads(L.params_json(c))) for i, s, c in SP.trials(m, SP.N_TRIALS)):
            bad.append(f"{m}: the 60 trials differ from those of its ledger without terms (pairing lost)")
        L.require_ledger(m, L.LEDGER_A6, "run", True)
    print("tag", TAG, "anchored" if anchored() else "NOT anchored: smoke, full and read will refuse")
    return bad


def search(tune, methods, n):
    for m in methods:
        tune_cli.main(["run", m, *manifests(tune), "--terms", f"ridgecrest={TERMS}", "--ledger", str(L.LEDGER_A6),
                       "--trials", str(n), "--env", "cloud-cpu"])


def paired(m):
    """objective with terms minus without, trial by trial (the same seeds give the same parameters)."""
    a = {int(r["trial"]): float(r["objective"]) for r in L.read(L.LEDGER_A6) if r["method"] == m
         and r["status"] == "ok"}
    b = {int(r["trial"]): float(r["objective"]) for r in L.read(METHODS[m][0]) if r["method"] == m
         and r["status"] == "ok"}
    d = np.array([a[i] - b[i] for i in sorted(set(a) & set(b))])
    return dict(pairs=len(d), up=int(np.sum(d > 0)), down=int(np.sum(d < 0)), equal=int(np.sum(d == 0)),
                median=float(np.median(d)), mean=float(np.mean(d)))


def read(methods):
    out = {}
    for m in methods:
        if len(L.done(m, L.LEDGER_A6)) < SP.N_TRIALS:
            raise SystemExit(f"{m}: the 60 trials are not all in the ledger")
        new = L.freeze(m, FROZEN / (m + "_a6.json"), L.LEDGER_A6)
        old = json.load(open(METHODS[m][1]))
        gain = new["objective"] - old["objective"]
        out[m] = dict(frozen=old["objective"], retuned=new["objective"], gain=gain,
                      H="holds" if gain >= GAIN else "no change" if gain > -GAIN else "lower",
                      used_against_model="with terms" if new["objective"] > old["objective"] else "without terms",
                      paired=paired(m))
    out["pivot"] = ("closed as negative" if all(v["gain"] < GAIN for v in out.values())
                    else "at least one method gains")
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(RESULT, "w"), indent=1, sort_keys=True)
    print(json.dumps(out, indent=1, sort_keys=True))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cmd", choices=("check", "smoke", "full", "read"))
    ap.add_argument("tune", nargs="?")
    ap.add_argument("--methods", nargs="+", default=list(METHODS), choices=list(METHODS))
    a = ap.parse_args(argv)
    if a.cmd != "read" and not a.tune:
        ap.error("the tune folder is required")
    if a.cmd == "check":
        bad = check(a.tune, a.methods)
        print("\n".join(["FAIL: " + b for b in bad]) if bad else "check passed: nothing was run")
        return 1 if bad else 0
    if not anchored():
        raise SystemExit(f"{TAG} is not anchored on this checkout: no trial with station terms may run")
    if a.cmd == "read":
        read(a.methods)
        return 0
    bad = check(a.tune, a.methods)
    if bad:
        raise SystemExit("\n".join(bad))
    search(a.tune, a.methods, 3 if a.cmd == "smoke" else SP.N_TRIALS)
    return 0


if __name__ == "__main__":
    sys.exit(main())
