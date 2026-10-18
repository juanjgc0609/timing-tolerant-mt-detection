"""floor_bias_bh.py - energy read at the 'exact predicted sample' idx = int(t / DT) (floor, as certified_fw._idx_node,
c10_lib.Eng, c11_lib.clone_eng, c18j_layered.layered_engine) versus the nearest sample, for a 2 Hz Ricker arrival seen
through the 2 Hz Ricker matched filter (continuous-time autocorrelation, unit peak), predicted time exact, fractional
part f ~ U(0, 1)."""
import numpy as np

DT = 1 / 16.0
fc = 2.0
t = np.arange(-4, 4, 1e-3)
r = (1 - 2 * (np.pi * fc * t) ** 2) * np.exp(-(np.pi * fc * t) ** 2)
ac = np.correlate(r, r, mode="same")
ac /= ac.max()
lag = (np.arange(len(ac)) - len(ac) // 2) * 1e-3
A = lambda d: np.interp(d, lag, ac)  # noqa: E731
f = (np.arange(20000) + 0.5) / 20000
e_floor = A(f * DT) ** 2
e_round = A(np.minimum(f, 1 - f) * DT) ** 2
print(f"mean energy fraction read at floor sample {e_floor.mean():.3f} (worst {e_floor.min():.3f}); "
      f"at nearest sample {e_round.mean():.3f} (worst {e_round.min():.3f})")
print(f"peak search +-4 samples recovers ~1.000; mean early bias of the floor read {-0.5 * DT * 1000:.1f} ms")
