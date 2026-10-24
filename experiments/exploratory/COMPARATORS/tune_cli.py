"""Run the ATLAS-0 random search of one comparator on the tuning partition, append to the ledger, freeze the best.

    python tune_cli.py run phasenet+pyocto manifest.json [--env ttmt-gpu] [--trials 60]
    python tune_cli.py freeze phasenet+pyocto frozen/phasenet+pyocto.json
    python tune_cli.py flags            # methods whose best development recall in T or C is outside [0.2, 0.85]
    python tune_cli.py sensitivity phasenet+pyocto
GaMMA runs with ncpu = 1 (fixed in the wrapper). Set TTMT_REAL_BIN for REAL. The manifest format is
ttmt.associate.tune.load_manifest.
"""
import argparse
import json
import sys

from ttmt.runtime import pin_blas

pin_blas(strict=False)

from ttmt.associate import ledger as L  # noqa: E402
from ttmt.associate import pipeline as PL  # noqa: E402
from ttmt.associate import tune as TU  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("method")
    r.add_argument("manifest")
    r.add_argument("--env", default="")
    r.add_argument("--trials", type=int, default=60)
    r.add_argument("--ledger", default=str(L.LEDGER))
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
    if a.cmd == "run":
        if not PL.STATUS.get(a.method, "").startswith("ready"):
            raise SystemExit(f"{a.method}: {PL.STATUS.get(a.method, 'unknown method')}")
        recs = TU.load_manifest(a.manifest, a.method)
        TU.run(a.method, PL.runner(a.method), recs, ledger=a.ledger, n_trials=a.trials, env=a.env)
    elif a.cmd == "freeze":
        print(json.dumps(L.freeze(a.method, a.out, a.ledger), indent=1))
    elif a.cmd == "sensitivity":
        print(json.dumps(L.sensitivity(a.method, a.ledger), indent=1))
    else:
        print(json.dumps(TU.realism_flags(a.ledger), indent=1))


if __name__ == "__main__":
    sys.exit(main())
