"""Conventions that change numbers, as explicit, versioned parameters (protocol 5.2; ADR-0002 item 5).

Sample index (E52). A predicted arrival time t (s) is read at a sample of the 1/DT grid:
  ``"floor"``   (legacy, every registered run): ``int(t / dt)``, truncation, so the read sample is on average half a
                sample early. For t >= 0 truncation equals the floor.
  ``"nearest"`` (corrected): the nearest sample, ties (exact half samples) rounded up: ``floor(t / dt + 0.5)``.
The fraction returned with the index is ``t / dt - index`` (in [0, 1) for floor, in [-0.5, 0.5) for nearest), which
is what the timing charge of ``c10_lib`` uses.
"""
import numpy as np

ROUNDING = ("floor", "nearest")
LEGACY_ROUNDING = "floor"


def sample_index(t_over_dt, rounding=LEGACY_ROUNDING):
    """index (int array) of the samples of x = t / dt under the rounding convention; no clipping."""
    x = np.asarray(t_over_dt, float)
    if rounding == "floor":
        return x.astype(int)                       # legacy expression (truncation), kept bit for bit
    if rounding == "nearest":
        return np.floor(x + 0.5).astype(int)
    raise ValueError(f"rounding must be one of {ROUNDING}, got {rounding!r}")


def sample_index_clipped(t_over_dt, n, rounding=LEGACY_ROUNDING):
    """(index clipped to [0, n - 1], fraction x - index): the engine tables of c10_lib.Eng and certified_fw."""
    x = np.asarray(t_over_dt, float)
    idx = np.clip(sample_index(x, rounding), 0, n - 1)
    return idx, x - idx
