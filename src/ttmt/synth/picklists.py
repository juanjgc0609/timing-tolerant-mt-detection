"""Pick-list scenarios for the associators (ATLAS-0, associator scenarios A and B; protocol of Puente et al. 2025).

No waveforms: events at a Poisson rate inside (and around) a network, true P and S arrivals from the true layered
medium (``media.first_arrival``), and picks made from them with
- timing error = station static (common to every event, per station and phase) + an event-specific term correlated
  in source position (a Gaussian process over the hypocentre, squared-exponential kernel of length ``corr_km``, one
  independent process per station and phase) + independent jitter; the three standard deviations split ``timing_sd``
  as ``split`` (variances);
- false picks: a fraction ``false_frac`` of the number of true picks, uniform in time over the stations, phase at
  random, probability drawn like the true picks';
- missing picks: each true arrival is kept with probability ``p_detect``.
Returns ``ttmt.atlas.schema`` Truth (arr_t true, arr_t_obs = picked times where picked) and Picks.
"""
import numpy as np

from ..atlas.schema import Picks, Truth
from .media import PRESETS, first_arrival

RATES_PER_DAY = (None, 100.0, 500.0, 2000.0)          # None: isolated (one event per 60 s, no overlap)
FALSE_FRACS = (0.0, 0.3, 1.0, 3.0)
TIMING_SDS = (0.0, 0.1, 0.2, 0.4)


def _gp_factor(xyz, corr_km):
    """Cholesky factor of the unit-variance squared-exponential covariance at the points xyz (m, 3)."""
    d2 = np.sum((xyz[:, None] - xyz[None]) ** 2, -1)
    K = np.exp(-0.5 * d2 / corr_km ** 2) + 1e-6 * np.eye(len(xyz))
    return np.linalg.cholesky(K)


def scenario(STA, duration_s, rng, rate_per_day=None, false_frac=0.0, timing_sd=0.0, split=(0.4, 0.4, 0.2),
             corr_km=5.0, p_detect=0.9, medium="crust_a", depth=(2.0, 14.0), out_factor=1.5, mag=(0.5, 3.0)):
    STA = np.asarray(STA, float)
    ns = len(STA)
    med = PRESETS[medium]
    if rate_per_day is None:
        t = np.arange(30.0, duration_s - 30.0, 60.0)
    else:
        n = rng.poisson(rate_per_day * duration_s / 86400.0)
        t = np.sort(rng.uniform(0.0, duration_s - 30.0, n))
    m = len(t)
    cen = STA[:, :2].mean(0)
    ap = float(np.max(np.linalg.norm(STA[:, :2] - cen, axis=1)))
    rr = ap * out_factor * np.sqrt(rng.random(m))
    az = rng.uniform(0, 2 * np.pi, m)
    xyz = np.c_[cen[0] + rr * np.cos(az), cen[1] + rr * np.sin(az), rng.uniform(*depth, m)]
    arr = np.zeros((m, ns, 2))
    for j in range(m):
        X = np.hypot(STA[:, 0] - xyz[j, 0], STA[:, 1] - xyz[j, 1])
        arr[j, :, 0] = t[j] + first_arrival(med, xyz[j, 2], X, "P")[0]
        arr[j, :, 1] = t[j] + first_arrival(med, xyz[j, 2], X, "S")[0]
    sd = np.sqrt(np.asarray(split, float) * timing_sd ** 2)
    err = np.zeros_like(arr)
    if timing_sd > 0:
        err += rng.normal(0.0, sd[0], (1, ns, 2))
        if m:
            L = _gp_factor(xyz, corr_km)
            err += sd[1] * np.einsum("jk,kib->jib", L, rng.normal(size=(m, ns, 2)))
        err += rng.normal(0.0, sd[2], err.shape)
    obs = arr + err
    keep = rng.random(arr.shape) < p_detect
    sta_i, ph_i = np.meshgrid(np.arange(ns), np.arange(2), indexing="ij")
    ps = np.broadcast_to(sta_i, arr.shape)[keep]
    pp = np.broadcast_to(ph_i, arr.shape)[keep]
    pt = obs[keep]
    pr = rng.uniform(0.3, 1.0, len(pt))
    nf = rng.poisson(false_frac * len(pt)) if false_frac > 0 else 0
    if nf:
        ps = np.r_[ps, rng.integers(0, ns, nf)]
        pp = np.r_[pp, rng.integers(0, 2, nf)]
        pt = np.r_[pt, rng.uniform(0.0, duration_s, nf)]
        pr = np.r_[pr, rng.uniform(0.3, 1.0, nf)]
    o = np.argsort(pt)
    picks = Picks(ps[o], pp[o], pt[o], pr[o])
    truth = Truth(np.array([f"k{j}" for j in range(m)]), t, xyz, rng.uniform(*mag, m), arr, np.full((m, ns), np.nan),
                  0.0, float(duration_s), arr_t_obs=obs,
                  meta={"rate_per_day": rate_per_day, "false_frac": false_frac, "timing_sd": timing_sd,
                        "split": list(split), "corr_km": corr_km, "p_detect": p_detect, "medium": medium,
                        "n_false": int(nf), "picked": keep.tolist()})
    return truth, picks
