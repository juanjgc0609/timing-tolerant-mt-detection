"""Development / evaluation partitions enforced by the reader (protocol 2.2; error E56).

A :class:`Partition` names disjoint time segments of one continuous record. A :class:`SegmentReader` serves analysis
windows of one segment and raises :class:`PartitionError` when any sample it would read, context margins included,
lies outside that segment. The legacy R1 development read its segment plus a 40 s margin and so read 40 s of the
evaluation segment (E56); :func:`legacy_t_end` keeps that behaviour available for the regressions only.

Per-station time shifts (for time-slide null calibration): each station's trace may be shifted by an independent
number of samples, applied as a circular shift WITHIN the segment, so a shifted read never leaves the segment.
"""
from dataclasses import dataclass

import numpy as np


class PartitionError(ValueError):
    """a read crosses the boundary of its segment."""


@dataclass(frozen=True)
class Partition:
    """segments: {name: (start s, end s)}, half-open, in record time; t0_sample: record time of sample 0 is 0."""
    segments: dict

    def __post_init__(self):
        iv = sorted(self.segments.values())
        for (a0, a1), (b0, b1) in zip(iv, iv[1:]):
            if b0 < a1:
                raise ValueError(f"segments overlap: {(a0, a1)} and {(b0, b1)}")
        for a, b in iv:
            if not b > a:
                raise ValueError(f"empty segment {(a, b)}")

    def bounds(self, segment):
        return self.segments[segment]

    def check(self, segment, t0, t1):
        """raise unless [t0, t1) lies inside the segment."""
        lo, hi = self.segments[segment]
        if t0 < lo - 1e-9 or t1 > hi + 1e-9:
            raise PartitionError(f"read [{t0:.3f}, {t1:.3f}) s crosses segment {segment!r} [{lo}, {hi}) s")


def legacy_t_end(partition, segment, margin=40.0):
    """r1_lib.t_end / r0_run.main: read up to the end of the segment plus the window margin (crosses into the next
    segment: E56). Regressions only."""
    return partition.bounds(segment)[1] + margin


class SegmentReader:
    """windows of one segment of a continuous array X (ns, 3, n) sampled at dt with sample 0 at record time t_origin.

    strict (default): every window must lie inside the segment, margins included. shifts: per-station circular
    shifts in samples, applied within the segment (None or zeros: identity)."""

    def __init__(self, X, dt, partition, segment, t_origin=0.0, shifts=None, strict=True):
        self.X, self.dt, self.partition, self.segment = X, float(dt), partition, segment
        self.t_origin, self.strict = float(t_origin), strict
        lo, hi = partition.bounds(segment)
        self.k_lo = int(round((lo - self.t_origin) / self.dt))
        self.k_hi = min(int(round((hi - self.t_origin) / self.dt)), X.shape[-1])
        self.shifts = np.zeros(X.shape[0], int) if shifts is None else np.asarray(shifts, int)
        if self.shifts.shape != (X.shape[0],):
            raise ValueError("one shift per station")

    def segment_view(self):
        """the segment's samples with the per-station circular shifts applied (copy)."""
        if getattr(self, "_view", None) is None:
            seg = self.X[:, :, self.k_lo:self.k_hi]
            self._view = seg.copy() if not self.shifts.any() else \
                np.stack([np.roll(seg[i], int(self.shifts[i]), axis=-1) for i in range(seg.shape[0])])
        return self._view

    def window(self, t_start, nt):
        """(ns, 3, nt) window starting at record time t_start (nearest sample, as r0_lib.window)."""
        k = int(round((t_start - self.t_origin) / self.dt))
        if self.strict:
            self.partition.check(self.segment, self.t_origin + k * self.dt, self.t_origin + (k + nt) * self.dt)
        a, b = k - self.k_lo, k - self.k_lo + nt
        if self.strict or self.shifts.any():
            if a < 0 or b > self.k_hi - self.k_lo:
                raise PartitionError(f"window [{a}, {b}) outside segment {self.segment!r}")
            return self.segment_view()[:, :, a:b]
        W = np.zeros((self.X.shape[0], 3, nt))
        seg = self.X[:, :, k:k + nt]
        W[:, :, :seg.shape[-1]] = seg
        return W
