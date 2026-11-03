"""rescore_2l.py - Task 2L (analyze_2l.py, unchanged) under the one-to-one hit rule of froc_strict (ADR-0007).
    python rescore/rescore_2l.py RULE CODE_COPY OUT
RULE       legacy  the registered rule (must reproduce results/task2l/analyze_2l_out.txt from the regenerated parts)
           strict  froc_strict.relabel(one_to_one=True): atoms matched in record order to the nearest eligible event
                   still unmatched; an atom with none is false.  Applied to both rules of the analysis, LOC+T (5 km and
                   0.5 s) and LOC, and to the descriptive geometry block, since all three go through c10_lib.label.
CODE_COPY  a copy of code/ whose task2l_parts/ holds the parts regenerated with the frozen code (task2l_code_hashes.txt
           checked, OPENBLAS_CORETYPE=Prescott, run_queue_2l.sh).  The 2L parts store every atom as
           [x, y, z, origin, s, s0, sref], so, unlike 2A to 2D, the full one-to-one reassignment can be applied.
OUT        output text file.
Only the labelling function is replaced, in memory; the frozen files keep their hashes."""
import os
import sys

import numpy as np


def main(rule, code, out):
    if rule not in ("legacy", "strict"):
        raise SystemExit(__doc__)
    code = os.path.abspath(code)
    sys.path[:0] = [os.path.join(code, "explore_contrib"), code]
    import c10_lib as C10
    import froc_strict as FS
    if rule == "strict":
        def label(atoms, truth, taus, tol_t=None, tol_km=5.0):
            locs = [list(np.asarray(a["x"], float)[:3]) + [float(a["tau"])] for a in atoms]
            return FS.relabel(locs, [np.asarray(t, float) for t in truth], taus, tol_t, tol_km, one_to_one=True)
        C10.label = label
    import analyze_2l as A
    assert A.C10 is C10
    sys.argv = ["analyze_2l.py", "--out", os.path.abspath(out)]
    A.main()


if __name__ == "__main__":
    main(*sys.argv[1:4])
