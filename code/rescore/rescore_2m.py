"""rescore_2m.py - Task 2M analysis (t2m_analyze.py, unchanged) under the registered or the one-to-one hit rule.
    python rescore/rescore_2m.py RULE OUT.txt B LABEL=PARTSDIR [...]
RULE: legacy (registered labels; must reproduce results/task2M/analyze_2m_out.txt) | strict (froc_strict.relabel:
one match per true event, later detections of a matched event are false; LOC and LOC+T as registered)."""
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for p in (os.path.join(HERE, "explore_contrib"), HERE):
    sys.path.insert(0, p)

import c10_lib as C  # noqa: E402
import froc_strict as FS  # noqa: E402
import t2m_analyze  # noqa: E402

rule = sys.argv.pop(1)
if rule == "strict":
    def label(atoms, truth, taus, tol_t=None, tol_km=5.0):
        locs = [list(a["x"][:3]) + [a["tau"]] for a in atoms]
        return FS.relabel(locs, truth, taus, tol_t, tol_km, one_to_one=True)
    C.label = label
elif rule != "legacy":
    raise SystemExit("RULE = legacy | strict")
t2m_analyze.main(sys.argv[1], int(sys.argv[2]), sys.argv[3:])
