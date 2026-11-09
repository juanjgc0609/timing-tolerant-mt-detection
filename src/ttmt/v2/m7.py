"""M7 false-alarm control, first step (specification M7.1; theory D.3): the time-slide null.

Exploratory: no registered run uses it. ``time_slide_pvalue`` is the specification's interface
``time_slide_pvalue(window, detector, B) -> float``: each station's trace block is rolled by an independent uniform
offset, the detector's statistic is recomputed, and p = (1 + #{T_b >= T}) / (B + 1). It is exact (super-uniform)
when the stations are independent and each one is circularly shift-invariant in law within the window (D.3). Coherent
noise across stations (transients seen by several stations at once) breaks the premise; experiment M7-TIMESLIDE
measures that on real noise before anything is built on it.
"""
import numpy as np


def slide(window, rng):
    """roll each station (axis 0) of window (ns, ..., nt) along the last axis by an independent uniform offset."""
    off = rng.integers(window.shape[-1], size=window.shape[0])
    return np.stack([np.roll(window[i], int(off[i]), axis=-1) for i in range(window.shape[0])]), off


def slide_stats(window, detector, B, rng):
    """the observed statistic and B slid statistics (the detector maps a window to a float)."""
    t = float(detector(window))
    return t, np.array([float(detector(slide(window, rng)[0])) for _ in range(B)])


def pvalue(t, tb):
    """p = (1 + #{T_b >= T}) / (B + 1); a non-finite observed T counts as the smallest p (conservative for the level
    test: a NaN never hides a rejection)."""
    tb = np.asarray(tb, float)
    if not np.isfinite(t):
        return 1.0 / (len(tb) + 1.0) if np.isnan(t) or t > 0 else 1.0
    return (1.0 + np.sum(tb >= t)) / (len(tb) + 1.0)


def time_slide_pvalue(window, detector, B, rng):
    t, tb = slide_stats(window, detector, B, rng)
    return pvalue(t, tb)
