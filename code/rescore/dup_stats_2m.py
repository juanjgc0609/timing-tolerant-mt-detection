"""dup_stats_2m.py - Task 2M REAL arm: share of atoms that are repeat hits under the registered LOC rule, atoms that
froc_strict.relabel reassigns to another unmatched event, and origin-time errors of LOC hits.
    python rescore/dup_stats_2m.py [RESULTS_TASK2M_DIR]"""
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for p in (os.path.join(HERE, "explore_contrib"), HERE):
    sys.path.insert(0, p)

import numpy as np  # noqa: E402

import c10_analyze as A  # noqa: E402
import c10_lib as C  # noqa: E402
import c16_analyze as Z  # noqa: E402
import froc_strict as FS  # noqa: E402

base = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "results", "task2M")
for lab in ("2718", "3141"):
    _nets, rows = Z.load(os.path.join(base, f"parts_{lab}"), "REAL")
    for des in ("DENSE", "STD"):
        sub = [r for r in rows if r["design"] == des]
        for m in ("MT@S0", "BOMP@S0", "HUBER-N@S0"):
            n = dup = re_ = 0
            err = []
            for r in sub:
                at = A.atoms_of(r, m)
                h = C.label(at, [np.array(t) for t in r["truth"]], r["taus"], None)
                d = FS.drop_duplicates(h)
                s = FS.relabel([list(a["x"]) + [a["tau"]] for a in at], r["truth"], r["taus"], None)
                n += len(h)
                dup += sum(a >= 0 > b for a, b in zip(h, d))
                re_ += sum(a < 0 <= b for a, b in zip(d, s))
                err += [abs(a["tau"] - r["taus"][j]) for a, j in zip(at, h) if j >= 0]
            e = np.array(err)
            print(f"{lab} {des:5s} {m:11s} atoms {n:5d}  repeat hits {dup:4d} ({dup / n:.3f})  reassigned {re_:3d}  "
                  f"LOC-hit |origin error| median {np.median(e):.2f} s, <=0.5 s {np.mean(e <= .5):.2f}, "
                  f"<=1 s {np.mean(e <= 1):.2f}")
