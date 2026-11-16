"""Run the ATLAS-0 random search of one comparator on the tuning partition, append to the ledger, freeze the best.

    python tune_cli.py records ss_ridgecrest_2019_egf_v1 ss_hinet_2012_egf_v1 --out C:/Users/juanj/ttmt/tune
    python tune_cli.py run phasenet+pyocto <out>/*/tune_phasenet.json [--env ttmt-gpu] [--trials 60]
    python tune_cli.py freeze phasenet+pyocto frozen/phasenet+pyocto.json
    python tune_cli.py flags            # methods whose best development recall in T or C is outside [0.2, 0.85]
    python tune_cli.py sensitivity phasenet+pyocto
GaMMA runs with ncpu = 1 (fixed in the wrapper). Set TTMT_REAL_BIN for REAL. The manifest format is
ttmt.associate.tune.load_manifest; "records" writes it from the generator's sets (ttmt.associate.records), reading
only their tuning partition (fingerprint check by default, since the waveforms are rebuilt on each machine).
"""
import argparse
import json
import sys

from ttmt.runtime import pin_blas

pin_blas(strict=False)

from ttmt.associate import ledger as L  # noqa: E402
from ttmt.associate import pipeline as PL  # noqa: E402
from ttmt.associate import tune as TU  # noqa: E402


def load_terms(path, frame):
    """(ns, 2) terms of a file written by ST-RETUNE/catalogue_terms.py --terms-file; the station codes must be the
    frame's, in the same order."""
    import numpy as np
    d = json.load(open(path))
    if [str(c) for c in d["codes"]] != [str(c) for c in frame.codes]:
        raise SystemExit(f"{path}: station codes differ from the frame's")
    return np.asarray(d["terms_P_S_s"], float)


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("records")
    b.add_argument("set_ids", nargs="+")
    b.add_argument("--out", required=True)
    b.add_argument("--verify", default="fingerprint", choices=("sha256", "fingerprint"))
    b.add_argument("--device", default=None)
    r = sub.add_parser("run")
    r.add_argument("method")
    r.add_argument("manifests", nargs="+")
    r.add_argument("--env", default="")
    r.add_argument("--trials", type=int, default=60)
    r.add_argument("--ledger", default=str(L.LEDGER))
    r.add_argument("--terms", action="append", default=[], metavar="KEY=FILE",
                   help="station terms (ATLAS-0 amendment 6) for the manifests whose path holds KEY; needs "
                        "--ledger tuning_ledger_a6.csv")
    f = sub.add_parser("freeze")
    f.add_argument("method")
    f.add_argument("out")
    f.add_argument("--ledger", default=str(L.LEDGER))
    s = sub.add_parser("sensitivity")
    s.add_argument("method")
    s.add_argument("--ledger", default=str(L.LEDGER))
    g = sub.add_parser("flags")
    g.add_argument("--ledger", default=str(L.LEDGER))
    a = ap.parse_args(argv)
    if a.cmd == "records":
        from ttmt.associate import records as R
        from ttmt.synth.sets import SetReader
        for sid in a.set_ids:
            print(R.write_manifest(SetReader(sid, "tune", verify=a.verify), a.out, R.seisbench_annotate(a.device)))
    elif a.cmd == "run":
        if not PL.STATUS.get(a.method, "").startswith("ready"):
            raise SystemExit(f"{a.method}: {PL.STATUS.get(a.method, 'unknown method')}")
        recs, env = [], a.env
        terms = dict(t.split("=", 1) for t in a.terms)
        for key, path in sorted(terms.items()):
            env += f" terms:{key}:{L.sha256(path)[:12]}"
        for m in a.manifests:
            keys = [k for k in terms if k in m]
            if len(keys) > 1:
                raise SystemExit(f"{m}: more than one --terms key matches")
            for r in TU.load_manifest(m, a.method):
                if keys:
                    r.data["terms"] = load_terms(terms[keys[0]], r.data["frame"])
                recs.append(r)
        if terms and not any("terms" in r.data for r in recs):
            raise SystemExit("--terms given but no manifest matches a key")
        TU.run(a.method, PL.runner(a.method), recs, ledger=a.ledger, n_trials=a.trials, env=env.strip(),
               terms=bool(terms))
    elif a.cmd == "freeze":
        print(json.dumps(L.freeze(a.method, a.out, a.ledger), indent=1))
    elif a.cmd == "sensitivity":
        print(json.dumps(L.sensitivity(a.method, a.ledger), indent=1))
    else:
        print(json.dumps(TU.realism_flags(a.ledger), indent=1))


if __name__ == "__main__":
    sys.exit(main())
