"""
c18h_kernel.py - exploration C18h (Ridgecrest CAL only; second candidate for the registered arrival-reading gate).
A higher-frequency matched filter: the Ricker kernel of the matched filter (2 Hz in every task) is replaced by one of
centre frequency FC, everywhere the project computes a matched filter (peak readings, station corrections, scores).
Then (a) the leave-one-event-out timing ladder of C18a and (b) the oracle scores of C18c (cap 6, LOG6, no cap) are
recomputed.
    python c18h_kernel.py FC SUBNET OUTPREFIX
"""
import sys

import numpy as np

sys.path.insert(0, "/home/claude/work/r0")
import r0_lib as R  # noqa: E402,F401  (sets the R0 medium first)

import task1o_consistent as T1O  # noqa: E402
from task0_common import DT  # noqa: E402


def set_kernel(fc):
    t = np.arange(-1.0, 1.0 + 1e-9, DT)
    x = (np.pi * fc * t) ** 2
    T1O._KER = (1 - 2 * x) * np.exp(-x)


if __name__ == "__main__":
    fc, sub, pre = float(sys.argv[1]), sys.argv[2], sys.argv[3]
    set_kernel(fc)
    import c18a_timing as TA  # noqa: E402
    import c18c_oracle as OC  # noqa: E402
    TA.main(sub, 2.0, pre + "_timing.json")
    OC.VARIANTS = ("CAP6", "LOG6", "NOCAP")
    OC.main(sub, pre + "_oracle.json")
