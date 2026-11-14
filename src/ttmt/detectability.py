"""Network detectability ceiling (T4.20): the detection probability of an ideal network subspace detector for one
event, at a false-alarm budget per hour, as the regime-C analogue of the resolution limit of Corollary A.7.

Model. Gaussian, stationary noise, independent between stations, and a signal s that lies in a known p-dimensional
subspace of the whitened data (Green's functions known: the detector is clairvoyant on the propagation, not on the
source). The generalised likelihood ratio statistic is ||P y||^2, chi^2_p under H0 and noncentral chi^2_p(lambda)
under H1, lambda = s' C^-1 s. Three subspaces bracket the detectors of the project:

    p = 1    coherent oracle: the waveform is known up to a scalar (template matching with the true template)
    p = 6    moment-tensor GLRT: the source is any moment tensor (the model's detector, M7 stage)
    p = K    per-station amplitude: one free scalar per station (stacking without a source model)

Signal energy. lambda is not observed in the truth files. The card SNR of ttmt.synth.inject.card_snr is the peak of
the band-passed vertical P window over the standard deviation of the band-passed noise. Bandpass-then-sample is one
linear filter, so its output SNR cannot exceed that of the whitened matched filter: SNR_card,s^2 <= lambda_s for each
station, and lambda >= sum_s SNR_card,s^2 over independent stations (S waves and horizontals ignored). The ceiling
computed with lambda = kappa^2 * sum SNR_card^2, kappa = 1, is therefore conservative in the signal energy; it is
optimistic in that the propagation is known and the hit rule's location test (LOC+T) is not imposed.

Threshold. The budget B false detections per hour is spread over N = 3600 * f_trial * n_loc trials per hour (origin
times on an f_trial grid, n_loc resolvable source positions). By the union bound the per-trial false-alarm
probability alpha = B / N keeps the expected false rate at or below B, whatever the correlation between trials.
"""
import numpy as np
from scipy import stats


def trial_alpha(budget_per_hour=30.0, f_trial=100.0, n_loc=1e4):
    """per-trial false-alarm probability for a budget of false detections per hour (union bound)."""
    return float(budget_per_hour) / (3600.0 * float(f_trial) * float(n_loc))


def threshold(p, alpha):
    """chi^2_p threshold for the per-trial false-alarm probability alpha."""
    return float(stats.chi2.isf(alpha, p))


def pd(lam, p, alpha):
    """detection probability of the p-dimensional subspace detector at noncentrality lam (array) and level alpha."""
    lam = np.asarray(lam, float)
    c = threshold(p, alpha)
    out = stats.ncx2.sf(c, p, np.maximum(lam, 1e-12))
    return np.where(np.isfinite(lam), out, np.nan)


def lam_lower(snr_card, kappa=1.0):
    """lower bound of lambda from per-station card SNR (nan stations ignored)."""
    v = np.asarray(snr_card, float)
    return float(kappa ** 2 * np.nansum(v ** 2)) if np.isfinite(v).any() else float("nan")


def ceiling(lams, ks, alpha, ps=(1, 6, "K")):
    """per-target detection probability for each detector in ps; "K" uses each target's number of stations."""
    lams, ks = np.asarray(lams, float), np.asarray(ks, int)
    out = {}
    for p in ps:
        if p == "K":
            v = np.full(len(lams), np.nan)
            for k in np.unique(ks):
                m = ks == k
                v[m] = pd(lams[m], max(int(k), 1), alpha)
            out["p=K"] = v
        else:
            out[f"p={p}"] = pd(lams, int(p), alpha)
    return out


def snr_tot_for_pd(p, alpha, target=0.5):
    """sqrt(lambda) at which the p-dimensional detector reaches detection probability target (bisection)."""
    lo, hi = 0.0, 100.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if pd(mid ** 2, p, alpha) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
