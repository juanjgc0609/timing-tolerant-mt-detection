"""Typed station terms (E54; protocol 3.2).

Two kinds of per-(station, phase) delays exist and must not be mixed:
  :class:`PeakDelay`   delay of the matched-filter ENERGY PEAK relative to the model prediction. Read by detectors
                       that read energy peaks (``ttmt.detector.Engine``).
  :class:`OnsetDelay`  delay of the picker ONSET relative to the model prediction. Read by pick-based methods
                       (associators, the pick support of the R1 hybrid).
An onset precedes the energy peak by a phase-dependent lag (0.24 to 0.27 s in R1), so subtracting a peak term from
picks pulls them early (E54). Every reader checks the type and raises ``TypeError`` on the wrong one.

Estimators ported unchanged: ``station_corrections.peak_delays`` and ``estimate_station_corrections`` (peak terms,
homogeneous medium), ``hyb_dev2.pick_terms`` (onset terms), ``r1_run.mf_terms`` (peak terms anchored on onset terms).
"""
from dataclasses import dataclass, field

import numpy as np

from .kernel import LEGACY as LEGACY_KERNEL

MIN_PS = 6        # samples: the P peak must precede the S peak by at least this much (station_corrections.MIN_PS)


@dataclass(frozen=True)
class _Terms:
    values: np.ndarray                     # (ns, 2): P, S delays in s
    info: dict = field(default_factory=dict, compare=False)

    def __post_init__(self):
        v = np.asarray(self.values, float)
        if v.ndim != 2 or v.shape[1] != 2:
            raise ValueError(f"station terms must have shape (ns, 2), got {v.shape}")
        object.__setattr__(self, "values", v)

    @classmethod
    def zeros(cls, ns):
        return cls(np.zeros((ns, 2)))

    def tolist(self):
        return self.values.tolist()


class PeakDelay(_Terms):
    """delays of matched-filter energy peaks (s), shape (ns, 2)."""
    kind = "peak"


class OnsetDelay(_Terms):
    """delays of picker onsets (s), shape (ns, 2)."""
    kind = "onset"


def require(terms, cls):
    """reader-side type check (E54)."""
    if not isinstance(terms, cls):
        raise TypeError(f"expected {cls.__name__} station terms, got {type(terms).__name__} (E54: energy-peak and onset "
                        "delays are not interchangeable)")
    return terms.values


# ============================================================================================ energy-peak readings
def peak_delays(STA, D, src, tau, vp, vs, search=1.0, joint=True, kernel=LEGACY_KERNEL):
    """station_corrections.peak_delays with the medium explicit (legacy: task0_common.VP, VS at import time).
    Per (station, phase): delay of the matched-filter energy peak relative to the homogeneous prediction d / V + tau,
    and the peak energy. joint: P and S peaks chosen together maximising E(kP) + E(kS) with kP + MIN_PS <= kS.
    Parabolic sub-sample refinement, clipped to +-0.5 sample."""
    dt = kernel.dt
    MF, sd = kernel.matched(D)
    E = np.sum((MF / sd) ** 2, axis=1)
    nt = E.shape[-1]
    d = np.linalg.norm(STA - src, axis=1)
    W = int(round(search / dt))
    out = np.zeros((len(STA), 2))
    en = np.zeros((len(STA), 2))
    for i in range(len(STA)):
        tp = [d[i] / vp + tau, d[i] / vs + tau]
        win = []
        for b in range(2):
            c = int(round(tp[b] / dt))
            win.append(np.arange(max(1, c - W), min(nt - 2, c + W) + 1))
        if joint:
            sc = E[i, win[0]][:, None] + E[i, win[1]][None, :]
            sc = np.where(win[1][None, :] >= win[0][:, None] + MIN_PS, sc, -np.inf)
            a, b_ = np.unravel_index(int(np.argmax(sc)), sc.shape)
            ks = (int(win[0][a]), int(win[1][b_]))
        else:
            ks = tuple(int(w[np.argmax(E[i, w])]) for w in win)
        for b in range(2):
            k = ks[b]
            e0, e1, e2 = E[i, k - 1], E[i, k], E[i, k + 1]
            den = e0 - 2 * e1 + e2
            dk = 0.5 * (e0 - e2) / den if den < 0 else 0.0
            out[i, b] = (k + float(np.clip(dk, -0.5, 0.5))) * dt - tp[b]
            en[i, b] = e1
    return out, en


def estimate_peak_terms(STA, cal, vp, vs, search=1.0, gate=25.0, joint=True, kernel=LEGACY_KERNEL):
    """station_corrections.estimate_station_corrections: median over the calibration windows cal = [(D, src, tau)] of
    the gated peak delays. Returns PeakDelay (info: nacc, the accepted readings per station and phase)."""
    ds, ens = [], []
    for D, src, tau in cal:
        dl, en = peak_delays(STA, D, np.asarray(src, float), float(tau), vp, vs, search, joint, kernel)
        ds.append(dl)
        ens.append(en)
    ds, ens = np.array(ds), np.array(ens)
    ok = ens >= gate
    corr = np.zeros(ds.shape[1:])
    nacc = ok.sum(0)
    for i in range(ds.shape[1]):
        for b in range(2):
            if nacc[i, b]:
                corr[i, b] = float(np.median(ds[ok[:, i, b], i, b]))
    return PeakDelay(corr, {"nacc": nacc})


# ============================================================================================ onset readings
def pick_arrays(picks):
    """hyb_dev2.pick_arrays: per station {phase: array of pick times}. picks: per station [(t, phase, prob)]."""
    return [{b: np.array([t for t, ph, _h in pl if ph == b]) for b in "PS"} for pl in picks]


def onset_terms(cat, STA, pa, times, lo, hi, mmin):
    """hyb_dev2.pick_terms: per station and phase, median over the catalogue events with origin in [lo, hi) and
    magnitude >= mmin of (nearest pick of the phase within +-2 s of the predicted arrival - predicted arrival); 0 with
    fewer than 2 readings. times(src, STA) -> (2, ns) predicted travel times. Returns OnsetDelay."""
    res = [[[], []] for _ in range(len(STA))]
    for e in cat:
        if not (lo <= e["t"] < hi and e["M"] >= mmin):
            continue
        tt = times([e["x"], e["y"], e["z"]], STA)
        for i in range(len(STA)):
            for b, ph in enumerate("PS"):
                p = pa[i][ph]
                d = p - (e["t"] + tt[b, i])
                d = d[np.abs(d) <= 2.0]
                if len(d):
                    res[i][b].append(float(d[np.argmin(np.abs(d))]))
    return OnsetDelay(np.array([[np.median(x) if len(x) >= 2 else 0.0 for x in r] for r in res]))


def anchored_peak_terms(calwin, STA, onset, model_times, vp, vs, cap=2.0, frac=0.5, gate=25.0, nmin=3,
                        kernel=LEGACY_KERNEL):
    """r1_run.mf_terms (amendment of R1): energy-peak terms anchored on onset terms. The energy-peak readings of the
    calibration windows calwin = [(D, src, tau)] are searched within +-min(cap, frac x predicted S-P) of the
    homogeneous prediction, converted to delays against the layered times (model_times(src, STA) -> (2, ns)); per
    phase the onset-to-peak lag is the median over station-phases with >= nmin gated readings of (median reading -
    onset term); term = onset term + lag. Returns PeakDelay (info: lag (2,), readings per phase)."""
    pick_terms = require(onset, OnsetDelay)
    ds, ens = [], []
    for D, src, tau in calwin:
        src = np.asarray(src, float)
        tl = model_times(src, STA)
        d = np.linalg.norm(STA - src, axis=1)
        th = np.stack([d / vp, d / vs], axis=1)
        dl, en = np.zeros((len(d), 2)), np.zeros((len(d), 2))
        for i in range(len(d)):
            a, b = peak_delays(STA, D, src, float(tau), vp, vs, min(cap, frac * float(tl[1, i] - tl[0, i])),
                               kernel=kernel)
            dl[i], en[i] = a[i], b[i]
        ds.append(dl + th - tl.T)
        ens.append(en)
    ds, ens = np.array(ds), np.array(ens)
    lag, info = np.zeros(2), {}
    for b in range(2):
        v = [float(np.median(ds[ens[:, i, b] >= gate, i, b])) - pick_terms[i, b] for i in range(ds.shape[1])
             if (ens[:, i, b] >= gate).sum() >= nmin]
        lag[b] = float(np.median(v)) if v else 0.0
        info["PS"[b]] = v
    return PeakDelay(pick_terms + lag[None, :], {"lag": lag, "lag_readings": info})
