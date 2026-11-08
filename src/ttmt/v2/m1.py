"""M1: per-station evidence without a threshold, and reading at the onset (specification v2, M1).

Why. The registered rule reads every block at the matched-filter energy peak within +-0.25 s of the prediction. The
energy peak lies 0.24-0.27 s after the onset and can jump to a side lobe or to another arrival (Appendix AB, B.1;
R1: 36 % of P and 20 % of S readings more than 0.5 s off). M1 reads the TIME at a local maximum of an onset
probability q (a picker's continuous output, never thresholded, never ``find_peaks`` with a distance), and the
AMPLITUDE vector on the signed 3-component matched filter at onset + a calibrated onset-to-peak lag, with sub-sample
interpolation.

Shapes. q (ns, 2, NT) on the detector grid (P, S), values in [0, 1]. Times are in samples of the detector grid
(floats allowed). The picker itself is not part of the package: ``seisbench_onset_prob`` adapts a SeisBench model
when that library is installed (picks environment), and any other picker can fill q directly.
"""
import numpy as np

from ..detector import Engine
from ..station_terms import OnsetDelay

FLOOR = 0.02            # technical floor of candidate maxima (computation only, not a decision threshold)


# ============================================================================================ onset evidence
def seisbench_onset_prob(stream, model, dt, t0, nt, phases=("P", "S")):
    """q (ns, 2, nt) from a SeisBench picker (``model.annotate``), resampled to the detector grid t0 + k dt by linear
    interpolation of each probability trace. stream: obspy Stream with one station per group of traces; the station
    order is the order of first appearance. Lazy import: the core environment does not need SeisBench."""
    ann = model.annotate(stream)
    stas = []
    for tr in stream:
        if tr.stats.station not in stas:
            stas.append(tr.stats.station)
    tg = t0 + np.arange(nt) * dt
    q = np.zeros((len(stas), 2, nt))
    for tr in ann:
        ph = tr.stats.channel.split("_")[-1]
        if ph not in phases or tr.stats.station not in stas:
            continue
        ts = tr.stats.starttime.timestamp + np.arange(tr.stats.npts) * tr.stats.delta
        q[stas.index(tr.stats.station), phases.index(ph)] = np.interp(tg, ts, tr.data, left=0.0, right=0.0)
    return q, stas


def local_maxima(q, floor=FLOOR):
    """boolean mask of the strict-left / weak-right local maxima of q along time above the technical floor (a plateau
    counts once, at its first sample). No distance, no threshold beyond the floor."""
    q = np.asarray(q, float)
    m = np.zeros(q.shape, bool)
    m[..., 1:-1] = (q[..., 1:-1] > q[..., :-2]) & (q[..., 1:-1] >= q[..., 2:])
    return m & (q >= floor)


def parabolic_offset(y0, y1, y2):
    """sub-sample offset in [-0.5, 0.5] of the vertex of the parabola through (-1, y0), (0, y1), (1, y2)."""
    den = y0 - 2.0 * y1 + y2
    with np.errstate(divide="ignore", invalid="ignore"):
        d = np.where(den < 0, 0.5 * (y0 - y2) / den, 0.0)
    return np.clip(d, -0.5, 0.5)


def peak_basin(q, k):
    """[lo, hi) of the basin of the local maximum at sample k of a 1-D curve: walk down both sides while q decreases.
    Used to suppress one accepted onset without touching the neighbouring maxima."""
    n = len(q)
    lo = k
    while lo > 0 and q[lo - 1] < q[lo]:
        lo -= 1
    hi = k + 1
    while hi < n and q[hi] <= q[hi - 1]:
        hi += 1
    return lo, hi


# ============================================================================================ calibration
class IsotonicCalibration:
    """isotonic regression (pool-adjacent-violators) of hit indicators on raw probabilities, so that q can be used as
    a likelihood and not only as a ranking (M1 step 2). fit(q_raw, hit) on candidate maxima of a development set
    with hit = 1 when the maximum lies within tolerance of a reference arrival; transform by step interpolation."""

    def fit(self, q, hit):
        q, y = np.asarray(q, float), np.asarray(hit, float)
        o = np.argsort(q, kind="stable")
        xs, ys = q[o], y[o]
        vals, wts, xr = [], [], []
        for x, v in zip(xs, ys):
            vals.append(v)
            wts.append(1.0)
            xr.append(x)
            while len(vals) > 1 and vals[-2] > vals[-1]:
                w = wts[-2] + wts[-1]
                v2 = (vals[-2] * wts[-2] + vals[-1] * wts[-1]) / w
                vals[-2:], wts[-2:], xr[-2:] = [v2], [w], [xr[-2]]
        self.x_, self.y_ = np.array(xr), np.array(vals)
        return self

    def transform(self, q):
        i = np.searchsorted(self.x_, np.asarray(q, float), side="right") - 1
        return self.y_[np.clip(i, 0, len(self.y_) - 1)]


def log_likelihood_ratio(qcal, base=FLOOR, eps=1e-6):
    """log q~: calibrated probability as a likelihood ratio against the technical floor (0 at the floor)."""
    return np.log(np.clip(qcal, eps, 1 - eps)) - np.log(base)


# ============================================================================================ onset engine
class OnsetEngine(Engine):
    """an Engine whose tables point at predicted ONSETS: travel time + origin + OnsetDelay terms (E54: onset terms
    are not energy-peak terms). ``tpred`` (N, T, 2, ns) holds the float sample times of every block."""
    terms_type = OnsetDelay


# ============================================================================================ reading
class OnsetReader:
    """M1 reader: in every block, among the local maxima of q within +-half_width samples of the predicted onset,
    the one with the largest q (sigma None; ties: the nearest to the prediction) or, with sigma (2,) per phase in
    samples, the largest log q - (s - mu)^2 / (2 sigma^2), the per-block local term of M4 ("prior" reading; on the
    R1 EVAL events P 0.97 within 0.25 s against 0.92 for the largest q, M1A-ONSET step 2). Refined to sub-sample by
    a parabola; a block with no maximum in its window is unread (NaN). Returns float onset samples, shape of the
    prediction."""

    def __init__(self, half_width=4, floor=FLOOR, sigma=None):
        self.half_width, self.floor = int(half_width), float(floor)
        self.sigma = None if sigma is None else np.asarray(sigma, float).reshape(2)

    def read(self, q, tpred):
        q = np.asarray(q, float)
        ns, _, nt = q.shape
        W = self.half_width
        mx = local_maxima(q, self.floor)
        qq = np.where(mx, q, -1.0)
        c = np.floor(np.asarray(tpred) + 0.5).astype(int)                 # (..., 2, ns)
        offs = np.arange(-W, W + 1)
        ii = np.clip(c[..., None] + offs, 0, nt - 1)                     # (..., 2, ns, 2W+1)
        b = np.arange(2).reshape((1,) * (c.ndim - 2) + (2, 1, 1))
        st = np.arange(ns).reshape((1,) * (c.ndim - 2) + (1, ns, 1))
        val = qq[st, b, ii]
        dist = ii - np.asarray(tpred)[..., None]
        if self.sigma is None:                                           # largest q, ties towards the prediction
            key = val - 1e-9 * np.abs(dist)
        else:
            sg = self.sigma.reshape((1,) * (c.ndim - 2) + (2, 1, 1))
            key = np.where(val >= self.floor, np.log(np.maximum(val, 1e-12)) - dist ** 2 / (2 * sg ** 2), -np.inf)
        a = key.argmax(-1)
        k = np.take_along_axis(ii, a[..., None], -1)[..., 0]
        ok = np.take_along_axis(val, a[..., None], -1)[..., 0] >= self.floor
        km, kp = np.clip(k - 1, 0, nt - 1), np.clip(k + 1, 0, nt - 1)
        b2, st2 = b[..., 0], st[..., 0]
        d = parabolic_offset(q[st2, b2, km], q[st2, b2, k], q[st2, b2, kp])
        return np.where(ok, k + d, np.nan)


def read_vectors(R, t):
    """(..., 2, ns, 3) vectors of R (ns, 3, NT) at float samples t (..., 2, ns) by linear interpolation; NaN times
    (unread blocks) give zero vectors."""
    ns, _, nt = R.shape
    ok = np.isfinite(t)
    tt = np.clip(np.where(ok, t, 0.0), 0, nt - 1)
    k0 = np.floor(tt).astype(int)
    k1 = np.minimum(k0 + 1, nt - 1)
    w = (tt - k0)[..., None]
    V = np.ascontiguousarray(R.transpose(0, 2, 1))                         # (ns, NT, 3)
    st = np.arange(ns).reshape((1,) * (t.ndim - 1) + (ns,))
    v = (1 - w) * V[st, k0] + w * V[st, k1]
    return np.where(ok[..., None], v, 0.0)


def onset_to_peak_lag(MF_list, onset_list, search=8, nmin=3):
    """(ns, 2) onset-to-peak lag in samples: per station and phase, the median over calibration windows of (energy
    peak of the matched filter within [onset, onset + search] samples - onset). MF_list: matched-filter outputs
    (ns, 3, NT); onset_list: reference onsets (ns, 2) in samples (NaN = not available). A block with fewer than nmin
    readings takes the median lag of its phase over the other stations (as anchored_peak_terms does per phase)."""
    ns = MF_list[0].shape[0]
    lags = [[[] for _ in range(2)] for _ in range(ns)]
    for MF, on in zip(MF_list, onset_list):
        E = (MF ** 2).sum(1)
        nt = E.shape[-1]
        for i in range(ns):
            for b in range(2):
                o = on[i][b]
                if not np.isfinite(o):
                    continue
                k0 = int(np.floor(o))
                w = np.arange(max(k0, 1), min(k0 + search + 1, nt - 1))
                if not len(w):
                    continue
                k = int(w[np.argmax(E[i, w])])
                d = parabolic_offset(E[i, k - 1], E[i, k], E[i, k + 1])
                lags[i][b].append(k + float(d) - o)
    out = np.full((ns, 2), np.nan)
    for b in range(2):
        for i in range(ns):
            if len(lags[i][b]) >= nmin:
                out[i, b] = np.median(lags[i][b])
        fill = np.nanmedian(out[:, b]) if np.isfinite(out[:, b]).any() else 0.0
        out[np.isnan(out[:, b]), b] = fill
    return out
