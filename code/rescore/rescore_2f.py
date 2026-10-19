"""rescore_2f.py - Task 2F analysis (analyze_2f.py, unchanged) under a chosen hit rule.
    python rescore/rescore_2f.py RULE OUT.txt
RULE  legacy        registered labels as stored in the parts (must reproduce results/task2f/analyze_2f_out.txt)
      legacy-loct   registered nearest-event rule with origin time (|tau error| <= 0.5 s), duplicates neutral
      strict        froc_strict.relabel, LOC: one match per event, later detections of a matched event are false
      strict-loct   the same with origin time within 0.5 s
BOMP-PK-CAP+SC (descriptive) stores no locations: under strict it gets froc_strict.drop_duplicates of its stored
labels (no reassignment); under the LOC+T rules it is dropped from the analysis."""
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)

import analyze_2f as A  # noqa: E402
import froc_strict as FS  # noqa: E402

TOL_T = 0.5
NOLOC = "BOMP-PK-CAP+SC"
rule = sys.argv.pop(1)
if rule not in ("legacy", "legacy-loct", "strict", "strict-loct"):
    raise SystemExit(__doc__)
tol_t = TOL_T if rule.endswith("loct") else None
one = rule.startswith("strict")
if tol_t is not None:
    A.METHODS = tuple(m for m in A.METHODS if m != NOLOC)
    A.PAIRS = tuple(p for p in A.PAIRS if NOLOC not in p)
_load = A.load_arm


def load_arm(seed, arm, parts=A.PARTS, nets=None, prefix=None):
    if rule == "legacy":
        return _load(seed, arm, parts, nets, prefix)
    import json
    nets = A.ARMS[arm]["nets"] if nets is None else nets
    prefix = arm if prefix is None else prefix
    D = _load(seed, arm, parts, nets, prefix)
    for reg in A.REGIMES:
        for cell in A.CELLS:
            W = D[(reg, cell)]
            wi = 0
            for j in nets:
                with open(os.path.join(parts, f"{prefix}_seed{seed}_{reg}_{cell}_net{j}.json")) as fh:
                    e = json.load(fh)
                for w in e["windows"]:
                    for m in A.METHODS:
                        rec = W["rec"][m][wi]
                        if m in w["loc"]:
                            assert len(w["loc"][m]) == len(rec)
                            assert FS.relabel(w["loc"][m], w["truth"], w["taus"], None, one_to_one=False) == \
                                [h for _s, h in rec], "registered labels not reproduced from the locations"
                            lab = FS.relabel(w["loc"][m], w["truth"], w["taus"], tol_t, one_to_one=one)
                        else:
                            lab = FS.drop_duplicates([h for _s, h in rec])
                        W["rec"][m][wi] = [[s, h] for (s, _h), h in zip(rec, lab)]
                    wi += 1
            assert wi == len(W["k"])
    return D


A.load_arm = load_arm
sys.argv = [sys.argv[0], "--out", sys.argv[1]]
A.main()
