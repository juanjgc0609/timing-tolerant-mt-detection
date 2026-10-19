"""Matched filter and deflation response, tied in one object (E51).

Legacy (``task1o_consistent._KER``, ``task1x_mtglrt.AC``): the matched-filter kernel is a Ricker wavelet of 2 Hz sampled
at DT = 1/16 s over [-1, 1] s; the deflation response ``AC`` is the matched-filter output of a 2 Hz Ricker pulse over
[-3, 3] s, normalised to 1 at its peak. In the legacy code the two lived in different modules and ``AC`` was computed
once at import, so ``c18h_kernel.set_kernel`` could change the filter without changing the deflation (E51). A
:class:`Kernel` computes both from the same parameters and cannot be split.
"""
from dataclasses import dataclass, field

import numpy as np

FC_LEGACY = 2.0
DT_LEGACY = 1.0 / (8 * FC_LEGACY)            # task0_common.DT = 1/16 s


def ricker(t, fc=FC_LEGACY):
    x = (np.pi * fc * t) ** 2
    return (1 - 2 * x) * np.exp(-x)


@dataclass(frozen=True)
class Kernel:
    """fc: centre frequency of the matched filter (Hz); pulse_fc: of the pulse whose filter output is the deflation
    response (default: fc, the autocorrelation of the kernel as in every registered run); dt: sample interval (s)."""
    fc: float = FC_LEGACY
    dt: float = DT_LEGACY
    pulse_fc: float = None
    taps: np.ndarray = field(init=False, repr=False, compare=False)
    ac: np.ndarray = field(init=False, repr=False, compare=False)
    ac_half: int = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        pf = self.fc if self.pulse_fc is None else self.pulse_fc
        object.__setattr__(self, "pulse_fc", float(pf))
        taps = ricker(np.arange(-1.0, 1.0 + 1e-9, self.dt), self.fc)
        a = np.convolve(ricker(np.arange(-3.0, 3.0 + 1e-9, self.dt), pf), taps[::-1], mode="same")
        ac = a / a.max()
        object.__setattr__(self, "taps", taps)
        object.__setattr__(self, "ac", ac)
        object.__setattr__(self, "ac_half", len(ac) // 2)

    @property
    def energy(self):
        """sum of squared taps (task2c_methods.E_R)."""
        return float(np.sum(self.taps ** 2))

    def matched(self, D):
        """task1o_consistent.matched: per station and component, convolution with the reversed kernel ('same'), and
        the robust sd of the whole window (median |MF| / 0.6745, floor 1e-12)."""
        MF = np.array([[np.convolve(D[i, c], self.taps[::-1], mode="same") for c in range(3)]
                       for i in range(D.shape[0])])
        sd = float(np.median(np.abs(MF)) / 0.6745)
        return MF, max(sd, 1e-12)

    def response(self, delta):
        """deflation response of a pulse delayed by delta samples (fractional), on the same support and with the
        same normalisation as ``ac`` (response(0.0) is ac bit for bit): the filter output of the pulse sampled at
        (k - delta) dt, k = -ac_half .. ac_half. Exact (analytic pulse), no interpolation."""
        a0 = np.convolve(ricker(np.arange(-3.0, 3.0 + 1e-9, self.dt), self.pulse_fc), self.taps[::-1], mode="same")
        a = np.convolve(ricker(np.arange(-3.0, 3.0 + 1e-9, self.dt) - delta * self.dt, self.pulse_fc),
                        self.taps[::-1], mode="same")
        return a / a0.max()

    def ac_at(self, lag):
        """AC(lag) in samples, 0 outside its support (task2c_methods.ac_at)."""
        lag = np.asarray(lag)
        ok = np.abs(lag) <= self.ac_half
        return np.where(ok, self.ac[np.clip(lag + self.ac_half, 0, len(self.ac) - 1)], 0.0)


LEGACY = Kernel()


def matched(D, kernel=LEGACY):
    return kernel.matched(D)
