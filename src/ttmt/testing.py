"""Comparison rules of the regression suite.

``exact``: bit-for-bit equality (Linux, the platform of the frozen runs: CI and the cloud check).
``tol``:   the declared tolerance of other platforms (Windows libm and OpenBLAS builds differ): absolute ``SCORE_ATOL``
           on detector scores and every other float field, identical counts, identical labels, identical locations
           and window indices, and identical analysis text (recalls, false counts and intervals are printed to three
           decimals and must not move). Never loosened silently: the tolerance is a named constant.
"""
import numpy as np

SCORE_ATOL = 1e-9


def compare_rows(got, ref, mode, score_cols=(4,), label=""):
    """rows of detections [x, y, z, t, s, ...]: same count; non-score columns exactly equal (mode exact) or within
    SCORE_ATOL; score columns equal (exact) or within SCORE_ATOL (tol). Returns a dict of differences."""
    out = {"label": label, "n_got": len(got), "n_ref": len(ref), "max_score_diff": 0.0, "max_other_diff": 0.0}
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
    out["max_other_diff"] = float(d[:, other].max()) if other else 0.0
    return out


def rows_ok(diff, mode):
    if diff.get("count_mismatch"):
        return False
    if mode == "exact":
        return diff["max_score_diff"] == 0.0 and diff["max_other_diff"] == 0.0
    return diff["max_score_diff"] <= SCORE_ATOL and diff["max_other_diff"] <= SCORE_ATOL
