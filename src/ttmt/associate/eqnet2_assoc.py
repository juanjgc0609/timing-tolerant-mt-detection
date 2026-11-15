"""EQNet2 (Wang and Zhu, 2026, arXiv 2610.17404; AI4EPS/EQNet branch develop 40a48ed1, weights EQNet2-v0
``eqnet_base_noprompt_v0.pth``) as an end-to-end waveform comparator, in two stages, as the pickers are.

Stage 1 (:func:`annotate`, GPU machine, once per record: the network output does not depend on any tuned parameter).
EQNet2 reads 3-component waveforms at 100 Hz in windows of 4096 samples (40.96 s) for a set of stations with their
longitude, latitude and depth, and returns per station an event-centre score at 1/16 of the sampling rate (0.16 s)
with, at each such sample, the travel time to that station, the epicentral distance and azimuth of the event from
the station and its depth below the station (``predict.py`` and ``eqnet.utils.postprocess.extract_events`` of that
commit). Windows overlap by half and each output sample is taken from the window where it is most central; stations
are padded to a multiple of 16 with silent copies of the first one (``predict.py --is_pad``), and the padding is
dropped. The scores are saved (:func:`save_scores`) as the picker probabilities are.

Stage 2 (:func:`detect`, CPU, per trial). As published: per station, peaks of the event score above the threshold,
origin time = peak time - travel time, epicentre = station + (distance, azimuth) on the local plane of the frame,
depth = station depth + predicted depth; then EQNet2's "time" clustering (``metrics_continuous.py``): a histogram of
the origin times weighted by the scores, local maxima over three bins whose summed score is at least MIN_STATION (3),
each origin assigned to the nearest such peak within one bin. A detection is the score-weighted mean origin and
location of its members, scored by their summed score.

Tuned (space.EQNET2): wf_threshold [0.05, 0.7] (the card's EQNet range); wf_min_interevent log-uniform [0.1, 5] s,
both the bin width of the clustering (published 2 s) and the minimum peak distance per station (published kernel 16
samples, 2.56 s). As published the two would keep EQNet2 from reporting two events inside regime T by construction, so
they are tuned on the range of ATLAS-0 amendment 1 (REAL tint, PyOcto min_interevent_time), as one parameter because
both are the method's minimum interevent time.

Score file (.npz): score, tt, dist, az, dz (n_sta, m) float32; dt (s, 0.16); t0 (record time of sample 0); codes;
model, weights (str).
"""
import numpy as np

from ..atlas.schema import Detections

FS = 100.0
NT = 4096
SCALE = 16
MIN_STATION = 3.0
WEIGHTS = "eqnet_base_noprompt_v0.pth"


def load_model(checkpoint, eqnet_path=None, device=None):
    """EQNetPredictor swin2 + upernet2 without prompt, strict load (keys stripped of torch.compile's prefix)."""
    import sys

    import torch
    if eqnet_path and eqnet_path not in sys.path:
        sys.path.insert(0, eqnet_path)
    import eqnet.models.eqnet as E
    m = E.build_model(backbone="swin2", head="upernet2", use_station_location=True, prediction=True, add_prompt=False)
    ck = torch.load(checkpoint, map_location="cpu", weights_only=False)
    state = {k.replace("_orig_mod.", ""): v for k, v in ck["model"].items()}
    m.load_state_dict(state, strict=True)
    dev = device or ("cuda" if torch.cuda.is_available() else "cpu")
    return m.to(dev).eval(), dev


def windows(n, nt=NT):
    """starts of half-overlapping windows covering n samples (the last one ends at n), at least one."""
    if n <= nt:
        return [0]
    hop = nt // 2
    s = list(range(0, n - nt + 1, hop))
    if s[-1] != n - nt:
        s.append(n - nt)
    return s


def annotate(model, device, wave, lld, t0=0.0):
    """wave (n_sta, 3, n) at 100 Hz, lld (n_sta, 3) longitude, latitude (deg), depth (km, positive down).
    Returns dict(score, tt, dist, az, dz) of shape (n_sta, m) with m = ceil(n / 16), dt, t0."""
    import torch
    ns, _, n = wave.shape
    n_pad = -(-max(n, NT) // SCALE) * SCALE           # whole output samples, so every window starts on one
    npad = -ns % 16
    w = np.zeros((ns + npad, 3, n_pad), np.float32)
    w[:ns, :, :n] = wave
    loc = np.concatenate([lld, np.repeat(lld[:1], npad, 0)]).astype(np.float32)
    m = -(-n_pad // SCALE)
    out = {k: np.zeros((ns, m), np.float32) for k in ("score", "tt", "dist", "az", "dz")}
    best = np.full(m, np.inf)                                    # distance of each output sample to its window centre
    me = NT // SCALE
    for s in windows(n_pad):
        x = torch.from_numpy(w[None, :, :, s:s + NT].transpose(0, 2, 3, 1).copy()).to(device)   # (1, 3, nt, ns)
        batch = {"data": x, "amplitude": x.clone(), "station_location": torch.from_numpy(loc[None]).to(device)}
        with torch.inference_mode():
            o = model(batch)
        ev = torch.sigmoid(o["event"][0]).float().cpu().numpy()[:, :ns].T        # (ns, 256)
        hyp = o["hypocenter"][0].float().cpu().numpy()[:, :, :ns]                # (5, 256, ns)
        j0 = s // SCALE
        idx = np.arange(me)
        dc = np.abs(idx + 0.5 - me / 2)
        take = dc < best[j0:j0 + me]
        cols = j0 + idx[take]
        best[cols] = dc[take]
        out["score"][:, cols] = ev[:, take]
        out["tt"][:, cols] = hyp[0][take].T
        out["dist"][:, cols] = hyp[1][take].T
        out["az"][:, cols] = hyp[2][take].T
        out["dz"][:, cols] = hyp[3][take].T
    keep = -(-n // SCALE)
    return {**{k: v[:, :keep] for k, v in out.items()}, "dt": SCALE / FS, "t0": float(t0)}


def save_scores(path, sc, codes, weights=WEIGHTS):
    np.savez_compressed(path, **{k: sc[k] for k in ("score", "tt", "dist", "az", "dz")}, dt=float(sc["dt"]),
                        t0=float(sc["t0"]), codes=np.asarray(codes, dtype=str), model="EQNet2", weights=str(weights))


def load_scores(path, codes=None):
    z = np.load(path, allow_pickle=False)
    have = z["codes"].tolist()
    o = [have.index(c) for c in codes] if codes is not None else list(range(len(have)))
    return {**{k: z[k][o] for k in ("score", "tt", "dist", "az", "dz")}, "dt": float(z["dt"]), "t0": float(z["t0"])}


def candidates(sc, sta_xyz, threshold, min_interevent):
    """per-station event candidates: (origin t, x, y, z, score), peaks above threshold at least min_interevent apart."""
    from scipy.signal import find_peaks
    dt, rows = sc["dt"], []
    dist = max(1, int(round(float(min_interevent) / dt)))
    for i in range(sc["score"].shape[0]):
        idx, prop = find_peaks(sc["score"][i], height=float(threshold), distance=dist)
        for j, h in zip(idx, prop["peak_heights"]):
            tt = max(0.0, float(sc["tt"][i, j]))
            d, a = float(sc["dist"][i, j]), float(sc["az"][i, j])
            rows.append((sc["t0"] + j * dt - tt, sta_xyz[i, 0] + d * np.sin(a), sta_xyz[i, 1] + d * np.cos(a),
                         sta_xyz[i, 2] + float(sc["dz"][i, j]), float(h)))
    return np.array(rows, float).reshape(-1, 5)


def cluster(c, bin_s):
    """EQNet2's "time" clustering of candidates c (t, x, y, z, w): Detections."""
    if len(c) == 0:
        return Detections(np.zeros(0), np.zeros((0, 3)), np.zeros(0))
    t, w = c[:, 0], c[:, 4]
    edges = np.arange(t.min() - bin_s, t.max() + 2 * bin_s, bin_s)
    hist, _ = np.histogram(t, bins=edges, weights=w)
    pad = np.r_[-np.inf, hist, -np.inf]
    peak = np.where((hist >= pad[:-2]) & (hist >= pad[2:]) & (hist >= MIN_STATION))[0]
    if len(peak) == 0:
        return Detections(np.zeros(0), np.zeros((0, 3)), np.zeros(0))
    b = np.clip(np.searchsorted(edges, t, side="right") - 1, 0, len(hist) - 1)
    k = np.abs(b[:, None] - peak[None, :])
    near = np.argmin(k, axis=1)
    ok = k[np.arange(len(t)), near] <= 1
    T, X, S = [], [], []
    for p in range(len(peak)):
        m = ok & (near == p)
        if not m.any():
            continue
        ww = w[m] / w[m].sum()
        T.append(float(ww @ t[m]))
        X.append(ww @ c[m, 1:4])
        S.append(float(w[m].sum()))
    return Detections(np.array(T), np.array(X).reshape(-1, 3), np.array(S))


def detect(sc, sta_xyz, cfg):
    """one trial: Detections from stored scores."""
    c = candidates(sc, np.asarray(sta_xyz, float), cfg["wf_threshold"], cfg["wf_min_interevent"])
    return cluster(c, float(cfg["wf_min_interevent"]))


def runner():
    """tune.run runner: the record's data must be {"frame": Frame, "eqnet2": score-file path}."""
    def run(cfg, record):
        fr = record.data["frame"]
        return detect(load_scores(record.data["eqnet2"], fr.codes), fr.stations.xyz, cfg)
    return run
