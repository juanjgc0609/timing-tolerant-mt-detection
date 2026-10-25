"""Pickers of the pick-based pipelines: probability traces once per (weights, record), picks per trial.

A picker's network output does not depend on the tuned parameters (peak distance, threshold), so the probability
traces of a record are computed once per weights (:func:`annotate_seisbench`, on the GPU machine) and stored
(:func:`save_probabilities`); every trial then only extracts picks (:func:`picks_from_probabilities`), as R0 and R1
did (``code/r0/r0_picks.py``: scipy ``find_peaks`` with a height and a minimum distance). PhaseNO writes the same
file layout from its own inference (it has no SeisBench port).

Probability file (.npz): prob (n_sta, n_samples, 2) float32, P and S; dt (s); t0 (record time of sample 0, s);
codes (n_sta,) station codes in the order of the frame; model, weights (str).
"""
import numpy as np

from ..atlas.schema import Picks


def picks_from_probabilities(prob, dt, threshold, peak_distance, t0=0.0):
    """Picks from probability traces (n_sta, n, 2): local maxima above threshold at least peak_distance s apart,
    per station and phase (scipy.signal.find_peaks)."""
    from scipy.signal import find_peaks
    mpd = max(1, int(round(float(peak_distance) / float(dt))))
    st, ph, tt, pp = [], [], [], []
    for i in range(prob.shape[0]):
        for k in (0, 1):
            idx, prop = find_peaks(prob[i, :, k], height=float(threshold), distance=mpd)
            st += [i] * len(idx)
            ph += [k] * len(idx)
            tt += list(t0 + idx * dt)
            pp += list(prop["peak_heights"])
    o = np.lexsort((np.asarray(ph), np.asarray(st), np.asarray(tt))) if tt else np.zeros(0, int)
    return Picks(np.asarray(st, int)[o], np.asarray(ph, int)[o], np.asarray(tt, float)[o], np.asarray(pp, float)[o])


def save_probabilities(path, prob, dt, t0, codes, model, weights):
    np.savez_compressed(path, prob=np.asarray(prob, np.float32), dt=float(dt), t0=float(t0),
                        codes=np.asarray(codes, dtype=str), model=str(model), weights=str(weights))


def load_probabilities(path, codes=None):
    """(prob, dt, t0, meta); with codes, the stations are reordered to codes (missing ones raise)."""
    z = np.load(path, allow_pickle=False)
    prob, have = z["prob"], z["codes"].tolist()
    if codes is not None:
        prob = prob[[have.index(c) for c in codes]]
    return prob, float(z["dt"]), float(z["t0"]), dict(model=str(z["model"]), weights=str(z["weights"]), codes=have)


def annotate_seisbench(stream, codes, model="PhaseNet", weights="original", device=None, t0=None):
    """probability traces (n_sta, n, 2) of a SeisBench model on an obspy Stream, one station per code, at the
    model's sampling rate; returns (prob, dt, t0). t0: record time (s) of the first sample, relative to the
    earliest start of the stream when None."""
    import seisbench.models as sbm
    import torch
    net = getattr(sbm, model).from_pretrained(weights)
    if device or torch.cuda.is_available():
        net.to(device or "cuda")
    ann = net.annotate(stream)
    start = min(tr.stats.starttime for tr in ann)
    dt = float(ann[0].stats.delta)
    n = max(int(round((tr.stats.endtime - start) / dt)) + 1 for tr in ann)
    prob = np.zeros((len(codes), n, 2), np.float32)
    for i, c in enumerate(codes):
        for k, ph in enumerate("PS"):
            sel = [tr for tr in ann if tr.stats.station == c and tr.stats.channel.endswith("_" + ph)]
            for tr in sel:
                a = int(round((tr.stats.starttime - start) / dt))
                prob[i, a:a + tr.stats.npts, k] = np.maximum(prob[i, a:a + tr.stats.npts, k], tr.data[:n - a])
    return prob, dt, (0.0 if t0 is None else float(t0))
