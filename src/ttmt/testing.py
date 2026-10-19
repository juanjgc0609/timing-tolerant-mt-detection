"""Comparison rules of the regression suite.

Discrete outputs are exact on every platform: counts, locations, origin times, window indices, labels, detection
lists, recall and false-detection counts, and the analysis texts.

Continuous detector scores (and the station terms they are built from):
``exact``  bit-for-bit equality. Linux, the platform of the frozen runs: CI and the cloud check.
``tol``    elsewhere, relative ``SCORE_RTOL`` = 1e-12. The libm and OpenBLAS builds of other platforms round
           differently. On Windows the differences are already visible upstream of the detector, in the station scale
           after filtfilt and resample (4.6e-12), and reach 4.8e-13 relative in the R1 P1 scores, while the legacy
           code on the same machine gives values identical to the package's (environments/DEVIATIONS.md).
The tolerance is a named constant and is never loosened silently.
"""
import numpy as np

SCORE_RTOL = 1e-12


def compare_rows(got, ref, mode, score_cols=(4,), label=""):
    """rows of detections [x, y, z, t, s, ...]: counts, the largest absolute difference of the non-score columns and
    the largest absolute and relative differences of the score columns."""
    out = {"label": label, "n_got": len(got), "n_ref": len(ref), "max_score_diff": 0.0, "max_score_rel": 0.0,
           "max_other_diff": 0.0}
    if len(got) != len(ref):
        out["count_mismatch"] = True
        return out
    if not len(got):
        return out
    g, r = np.asarray(got, float), np.asarray(ref, float)
    d = np.abs(g - r)
    sc = list(score_cols)
    other = [c for c in range(g.shape[1]) if c not in sc]
    out["max_score_diff"] = float(d[:, sc].max())
    out["max_score_rel"] = float((d[:, sc] / np.maximum(np.abs(r[:, sc]), 1e-300)).max())
    out["max_other_diff"] = float(d[:, other].max()) if other else 0.0
    return out


def rows_ok(diff, mode):
    """discrete columns exact in both modes; scores exact (mode exact) or within SCORE_RTOL relative (mode tol)."""
    if diff.get("count_mismatch") or diff["max_other_diff"] != 0.0:
        return False
    if mode == "exact":
        return diff["max_score_diff"] == 0.0
    return diff["max_score_rel"] <= SCORE_RTOL


def close(got, ref, mode):
    """continuous arrays: equal (exact) or within SCORE_RTOL relative (tol)."""
    got, ref = np.asarray(got, float), np.asarray(ref, float)
    if mode == "exact":
        return np.array_equal(got, ref)
    return bool(np.all(np.abs(got - ref) <= SCORE_RTOL * np.abs(ref)))
