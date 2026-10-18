"""rescore_2i.py - Task 2I analysis (analyze_2i.py, unchanged) under a chosen hit rule.
    python rescore/rescore_2i.py RULE OUT.txt
RULE  legacy | legacy-loct | strict | strict-loct, as in rescore_2f.py.  With only the prim arm regenerated, legacy
must reproduce results/task2i/analyze_2i_primary_out.txt."""
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)

import analyze_2i as A  # noqa: E402
import froc_strict as FS  # noqa: E402

rule = sys.argv.pop(1)
if rule not in ("legacy", "legacy-loct", "strict", "strict-loct"):
    raise SystemExit(__doc__)
tol_t = A.TOL_T if rule.endswith("loct") else None
one = rule.startswith("strict")
_load = A.load


def load(prefix, seed, nets, parts=A.PARTS):
    D, meta = _load(prefix, seed, nets, parts)
    if D is None or rule == "legacy":
        return D, meta
    for Dk in D.values():
        for W in Dk.values():
            for m in W["rec"]:
                for wi in range(len(W["k"])):
                    rec, loc = W["rec"][m][wi], W["loc"][m][wi]
                    assert len(loc) == len(rec)
                    assert FS.relabel(loc, W["truth"][wi], W["taus"][wi], None, one_to_one=False) == \
                        [h for _s, h in rec], "registered labels not reproduced from the locations"
                    lab = FS.relabel(loc, W["truth"][wi], W["taus"][wi], tol_t, one_to_one=one)
                    W["rec"][m][wi] = [[s, h] for (s, _h), h in zip(rec, lab)]
    return D, meta


A.load = load
sys.argv = [sys.argv[0], "--out", sys.argv[1]]
A.main()
