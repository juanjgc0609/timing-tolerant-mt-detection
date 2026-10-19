"""
froc_strict.py - one-to-one hit rule for the exact FROC (proposed correction to the rule of froc_exact.py).

Registered rule (froc_exact.greedy_record, task2c_methods.hit_of, c10_lib.label): each detection is labelled with the
nearest true event within 5 km (LOC), or the nearest one within 5 km whose origin is within tol_t (LOC+T).  The FROC
counts an event once (newhit in analyze_2f.Pool and c10_analyze.Boot), so a second detection of an event already hit
is NEUTRAL: it is neither a hit nor a false event, and a method that emits several atoms for one earthquake pays
nothing for them.

Corrected rule (this module): detections are matched in score order (the order of the record, whose scores are
non-increasing) to the nearest eligible true event that is still UNMATCHED; a detection with no unmatched eligible
event is FALSE.  This is the convention of JAFROC (marks after the first on a lesion are non-lesion localisations) and
of COCO-style detection scoring.  Because the matching only looks at earlier detections, every threshold still gives
a prefix of the record, so the exact FROC (all scores as thresholds) and the bootstrap code stay unchanged: only the
labels change.  With one_to_one=False the function reproduces the registered labels (used as a check).
"""
import numpy as np

TOL_KM = 5.0


def relabel(locs, truth, taus, tol_t=None, tol_km=TOL_KM, one_to_one=True):
    """locs: per detection [x, y, (z), tau, ...] in record order (tau = locs[i][3]); truth: list of [x, y, z];
    taus: true origin times.  Returns a label per detection: index of the matched event, or -1 (false)."""
    tr = np.asarray(truth, float).reshape(-1, 3) if len(truth) else np.zeros((0, 3))
    ta = np.asarray(taus, float)
    used = np.zeros(len(tr), bool)
    out = []
    for a in locs:
        h = -1
        if len(tr):
            d = np.linalg.norm(tr[:, :2] - np.asarray(a[:2], float), axis=1)
            if tol_t is None and not one_to_one:
                j = int(np.argmin(d))           # registered LOC: nearest event, then the distance test
                h = j if d[j] <= tol_km else -1
            else:
                ok = d <= tol_km
                if tol_t is not None:
                    ok &= np.abs(ta - float(a[3])) <= tol_t + 1e-9
                if one_to_one:
                    ok &= ~used
                if ok.any():
                    h = int(np.where(ok, d, np.inf).argmin())
        if h >= 0 and one_to_one:
            used[h] = True
        out.append(h)
    return out


def drop_duplicates(hits):
    """the corrected rule when only the registered labels are stored (no locations): a repeated event index becomes
    false.  Weaker than relabel (no reassignment to another unmatched event within tolerance)."""
    seen, out = set(), []
    for h in hits:
        if h >= 0 and h in seen:
            out.append(-1)
        else:
            if h >= 0:
                seen.add(h)
            out.append(h)
    return out
