"""CR-04: prevalence of the ATLAS-0 strata in a reviewed real catalogue (docs/protocol/CR04_PREVALENCE.md).

Reads only catalogue origins and hypocentres, station coordinates and data availability: never a pick, a detection
or a waveform. The separation is the one of the card (`cells.separation`), computed only over the events whose origin
lies within `window` s of the event, which gives the same bins as the full computation whenever window exceeds the
largest travel time plus 3 s (events farther apart cannot have arrivals closer than 3 s)."""
import numpy as np

from . import cells as C

BLOCK_S = 3600.0          # bootstrap block: 1 hour of origin time
RATE_HALF_S = 300.0       # local rate of the Poisson reference: events in a centred 10-minute window


def arrivals(origins, src_xyz, sta_xyz, tt):
    """(m, n, 2) absolute P and S arrival times. tt(src (3,), STA (n, 3)) -> (2, n) travel times
    (e.g. ttmt.traveltime.SourceTimes)."""
    origins, src_xyz, sta_xyz = np.asarray(origins, float), np.asarray(src_xyz, float), np.asarray(sta_xyz, float)
    out = np.empty((len(origins), len(sta_xyz), 2))
    for e in range(len(origins)):
        out[e] = origins[e] + np.asarray(tt(src_xyz[e], sta_xyz), float).T
    return out


def separation_windowed(arr_t, origins, sta_idx, window=60.0):
    """cells.separation restricted to the other events with |origin difference| <= window (inf when there is none)."""
    arr_t, origins = np.asarray(arr_t, float), np.asarray(origins, float)
    order = np.argsort(origins, kind="stable")
    so = origins[order]
    out = np.full(len(origins), np.nan)
    for e in range(len(origins)):
        st = np.asarray(sta_idx[e], int)
        if st.size == 0:
            continue
        lo, hi = np.searchsorted(so, origins[e] - window), np.searchsorted(so, origins[e] + window, side="right")
        cand = order[lo:hi]
        cand = cand[cand != e]
        d = []
        for s in st:
            tp = arr_t[e, s, 0]
            if np.isnan(tp):
                continue
            other = arr_t[cand, s, :].ravel()
            other = other[~np.isnan(other)]
            d.append(np.abs(other - tp).min() if other.size else np.inf)
        if d:
            out[e] = float(np.median(d))
    return out


def subsets(sta_xyz, codes, available, k, net_centroid=None):
    """per event the cells.station_subset of size k (None when fewer than k stations are available)."""
    return [C.station_subset(sta_xyz, codes, row, k, net_centroid) for row in np.asarray(available, bool)]


def _shares(sep, dist, ok):
    cnt = np.zeros((len(C.SEP_LABELS), len(C.DIST_LABELS)))
    if ok.any():
        np.add.at(cnt, (C.bin_of(sep[ok], C.SEP_EDGES), C.bin_of(dist[ok], C.DIST_EDGES)), 1.0)
    return cnt


def target_share(sep, ok):
    """Y: share of the counted events whose separation is in the target regime T = [0.1, 1) s."""
    if not ok.any():
        return np.nan
    b = C.bin_of(sep[ok], C.SEP_EDGES)
    return float(np.isin(b, C.TARGET_SEP_BINS).mean())


def block_bootstrap(origins, sep, ok, n_boot=2000, seed=0, block=BLOCK_S):
    """95 % percentile interval of Y over blocks of `block` s of origin time resampled with replacement."""
    origins = np.asarray(origins, float)
    blk = np.floor((origins - origins.min()) / block).astype(int)
    ids = np.unique(blk)
    rng = np.random.default_rng(seed)
    inT = np.isin(C.bin_of(np.where(ok, sep, np.inf), C.SEP_EDGES), C.TARGET_SEP_BINS) & ok
    num = np.array([inT[blk == b].sum() for b in ids], float)
    den = np.array([ok[blk == b].sum() for b in ids], float)
    ys = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(ids), len(ids))
        d = den[pick].sum()
        if d > 0:
            ys.append(num[pick].sum() / d)
    if not ys:
        return (np.nan, np.nan)
    return (float(np.quantile(ys, 0.025)), float(np.quantile(ys, 0.975)))


def poisson_reference(origins, src_xyz, sta_xyz, subs, tt, window=60.0, seed=0, half=RATE_HALF_S):
    """Y under a Poisson reference: for every counted event, partners drawn at its local catalogue rate (events in
    origin +- half, itself excluded) over +- window s, each with the hypocentre of a random catalogue event of that
    10-minute window; the event is in T when the nearest partner arrival (cells.separation rule) is in [0.1, 1) s."""
    origins, src_xyz, sta_xyz = np.asarray(origins, float), np.asarray(src_xyz, float), np.asarray(sta_xyz, float)
    rng = np.random.default_rng(seed)
    order = np.argsort(origins, kind="stable")
    so = origins[order]
    hit, n = 0, 0
    for e, st in enumerate(subs):
        if st is None:
            continue
        n += 1
        lo, hi = np.searchsorted(so, origins[e] - half), np.searchsorted(so, origins[e] + half, side="right")
        near = order[lo:hi]
        near = near[near != e]
        lam = len(near) / (2.0 * half)
        k = rng.poisson(lam * 2.0 * window)
        if k == 0 or near.size == 0:
            continue
        te = origins[e] + np.asarray(tt(src_xyz[e], sta_xyz[st]), float)[0]
        d = np.full(len(st), np.inf)
        for _ in range(k):
            j = near[rng.integers(0, near.size)]
            tj = origins[e] + rng.uniform(-window, window) + np.asarray(tt(src_xyz[j], sta_xyz[st]), float)
            d = np.minimum(d, np.abs(tj - te[None, :]).min(axis=0))
        if 0.1 <= np.median(d) < 1.0:
            hit += 1
    return hit / n if n else np.nan


def dead_time(origins):
    """(minimum, 1st percentile) of the inter-origin times of the catalogue, s."""
    dt = np.diff(np.sort(np.asarray(origins, float)))
    return (float(dt.min()), float(np.quantile(dt, 0.01))) if dt.size else (np.nan, np.nan)


def table(origins, src_xyz, sta_xyz, codes, available, tt, window=60.0, n_boot=2000, seed=0, poisson=True):
    """per station count k: stratum counts (separation x distance), shares, Y with its block-bootstrap interval and
    the Poisson reference, plus the catalogue dead time. Distance is to the catalogue centroid (card, axis 4)."""
    origins, src_xyz = np.asarray(origins, float), np.asarray(src_xyz, float)
    arr = arrivals(origins, src_xyz, sta_xyz, tt)
    dist = np.hypot(*(src_xyz[:, :2] - src_xyz[:, :2].mean(axis=0)).T)
    out = {"n_events": int(len(origins)), "dead_time_s": dead_time(origins), "window_s": window, "seed": seed,
           "sep_labels": C.SEP_LABELS, "dist_labels": C.DIST_LABELS, "k": {}}
    for k in C.STATION_COUNTS:
        subs = subsets(sta_xyz, codes, available, k)
        ok = np.array([s is not None for s in subs], bool)
        st = [s if s is not None else np.zeros(0, int) for s in subs]
        sep = separation_windowed(arr, origins, st, window)
        ok &= ~np.isnan(sep)
        cnt = _shares(sep, dist, ok)
        out["k"][k] = {
            "n_counted": int(ok.sum()),
            "counts": cnt.tolist(),
            "shares": (cnt / cnt.sum()).tolist() if cnt.sum() else None,
            "Y": target_share(sep, ok),
            "Y_ci95": block_bootstrap(origins, sep, ok, n_boot, seed),
            "Y_poisson": poisson_reference(origins, src_xyz, sta_xyz, subs, tt, window, seed) if poisson else None,
        }
    return out
