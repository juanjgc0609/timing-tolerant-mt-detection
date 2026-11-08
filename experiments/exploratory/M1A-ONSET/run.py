"""M1A-ONSET: gate G-M1a of the v2 specification (plan T4.13), on the R1 calibration events (Hi-net, CAL segment).

Question. With the reading moved from the matched-filter energy peak to the onset (the largest local maximum of the
PhaseNet probability, no threshold, no peak distance), do >= 90 % of the P readings lie within 0.25 s of the
reference time, with leave-one-event-out station terms? The registered energy-peak reading gives 55 % for P and 76 %
for S (``results/taskR1/posthoc/timing_cal_out.txt``).

Same events, windows and reference as ``code/r1/posthoc_timing_cal.py``: the 17 CAL truth events with M >= 1.0,
the 7 stations of H7, JMA2001 times from the catalogue location, search +-min(2 s, half the predicted S-P) around
the prediction. Reference time = origin + JMA2001 time + term, term = median over the OTHER events of the reading
delay (leave one event out). Two denominators: (a) the blocks that pass the registered energy gate (25), the same
blocks as the registered figure; (b) every block with an onset maximum in its window.
It also prints the onset-to-peak lag per phase (median of energy-peak reading - onset reading), the input of the v2
loop.

PhaseNet: the original model of the repository (``code/real_picker.py`` via ``r0_picks.phasenet_hour``, the same call
as the registered R1 picks), on the 100 Hz traces of the CAL segment. Its probabilities are cached next to this
script (``q_cal.npz``) so that a re-run does not repeat the inference.

Run on the machine with the raw data and the picks environment (TensorFlow PhaseNet):
    python experiments/exploratory/M1A-ONSET/run.py > experiments/exploratory/M1A-ONSET/run_out.txt
"""
import sys
from pathlib import Path

import numpy as np

from ttmt import legacy
from ttmt.r1 import Net, model_times
from ttmt.station_terms import peak_delays
from ttmt.v2.m1 import FLOOR, local_maxima, parabolic_offset

HERE = Path(__file__).resolve().parent
FS = 100.0
GATE = 25.0


def phasenet_q(net):
    cache = HERE / "q_cal.npz"
    if cache.exists():
        return np.load(cache)["P"]
    legacy.install()
    import r0_picks
    P = r0_picks.phasenet_hour(net.raw).astype(np.float32)          # (ns, n, 3): noise, P, S
    np.savez_compressed(cache, P=P)
    return P


def onset_read(qb, t_pred, half):
    """largest local maximum of qb (100 Hz curve) within t_pred +- half (s), sub-sample by parabola; NaN if none."""
    lo, hi = max(1, int(np.floor((t_pred - half) * FS))), min(len(qb) - 2, int(np.ceil((t_pred + half) * FS)))
    if hi <= lo:
        return np.nan, 0.0
    k = np.arange(lo, hi + 1)
    mx = local_maxima(qb[None, lo - 1:hi + 2])[0][1:-1]
    k = k[mx]
    if not len(k):
        return np.nan, 0.0
    j = int(k[np.argmax(qb[k])])
    d = float(parabolic_offset(qb[j - 1], qb[j], qb[j + 1]))
    return (j + d) / FS, float(qb[j])


def summary(name, delays, ok):
    """delays (events, stations, 2) reading - prediction (s), ok mask; leave-one-event-out residuals."""
    out = []
    for b in range(2):
        res = []
        for e in range(delays.shape[0]):
            for i in range(delays.shape[1]):
                if not ok[e, i, b]:
                    continue
                oth = [delays[k, i, b] for k in range(delays.shape[0]) if k != e and ok[k, i, b]]
                res.append(delays[e, i, b] - (np.median(oth) if oth else 0.0))
        r = np.array(res)
        if not len(r):
            out.append(f"   {name:34s} {'PS'[b]}: no readings")
            continue
        mad = 1.4826 * np.median(np.abs(r - np.median(r)))
        out.append(f"   {name:34s} {'PS'[b]}: n {len(r):3d}, robust sd {mad:.3f} s, within 0.25 s "
                   f"{np.mean(np.abs(r) <= 0.25):.2f}, within 0.5 s {np.mean(np.abs(r) <= 0.5):.2f}")
    return out


def main():
    net = Net("cal")
    P = phasenet_q(net)
    ne, ns = len(net.calev), len(net.STA)
    d_peak = np.full((ne, ns, 2), np.nan)
    en = np.zeros((ne, ns, 2))
    d_on = np.full((ne, ns, 2), np.nan)
    qv = np.zeros((ne, ns, 2))
    from ttmt.datasets import HINET as R
    used = 0
    for e, ev in enumerate(net.calev):
        if ev["t"] - 6.5 < 0:
            continue
        used += 1
        D = R.window(net.d16, ev["t"] - 6.5)
        src = np.array([ev["x"], ev["y"], ev["z"]], float)
        tl = model_times(src, net.STA)                                  # (2, ns) JMA2001
        d = np.linalg.norm(net.STA - src, axis=1)
        th = np.stack([d / R.vp, d / R.vs], axis=1)
        for i in range(ns):
            half = min(2.0, 0.5 * float(tl[1, i] - tl[0, i]))
            a, b_ = peak_delays(net.STA, D, src, 6.5, R.vp, R.vs, half)
            d_peak[e, i] = a[i] + th[i] - tl[:, i]
            en[e, i] = b_[i]
            for b in range(2):
                tp = ev["t"] + float(tl[b, i])
                t_on, qq = onset_read(P[i, :, 1 + b], tp, half)
                d_on[e, i, b] = t_on - tp
                qv[e, i, b] = qq
    gated = (en >= GATE) & np.isfinite(d_peak)
    lines = [f"M1A-ONSET  R1 CAL, {used} calibration events (truth M >= 1.0), {ns} stations, PhaseNet original model,",
             f"onset = largest local maximum of q >= {FLOOR} in +-min(2 s, S-P / 2); leave-one-event-out terms.",
             f"blocks above the energy gate: P {int(gated[..., 0].sum())}, S {int(gated[..., 1].sum())}", ""]
    lines += summary("(a) energy peak (registered rule)", d_peak, gated)
    lines += summary("(a) onset, same blocks", d_on, gated & np.isfinite(d_on))
    lines += summary("(b) onset, every block with a maximum", d_on, np.isfinite(d_on))
    lag = [float(np.nanmedian((d_peak - d_on)[..., b][gated[..., b] & np.isfinite(d_on[..., b])])) for b in range(2)]
    lines += ["", f"onset-to-peak lag (median energy-peak reading - onset reading, gated blocks): P {lag[0]:+.3f} s, "
                  f"S {lag[1]:+.3f} s",
              "G-M1a passes if line '(a) onset, same blocks' P shows within 0.25 s >= 0.90."]
    print("\n".join(lines))


if __name__ == "__main__":
    sys.exit(main())
