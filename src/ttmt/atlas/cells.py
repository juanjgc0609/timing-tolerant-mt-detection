"""Regime axes of the ATLAS-0 card and the cell of every event.

Axes (card, "Regime grid"): separation of arrivals, network signal-to-noise ratio, station subset size, distance to
the centroid of the development catalogue. Every function here is deterministic and reads only the truth and the
data-availability mask: never a pick or a detection."""
import numpy as np

SEP_EDGES = (0.1, 0.5, 1.0, 3.0)                 # s; bins < 0.1, 0.1-0.5, 0.5-1, 1-3, > 3
SNR_EDGES = (4.0, 8.0, 20.0)                     # linear amplitude ratio; bins < 4, 4-8, 8-20, >= 20
STATION_COUNTS = (4, 6, 8, 12)
DIST_EDGES = (10.0, 30.0)                        # km; bins < 10, 10-30, > 30
SEP_LABELS = ("<0.1", "0.1-0.5", "0.5-1", "1-3", ">3")
SNR_LABELS = ("<4", "4-8", "8-20", ">=20")
DIST_LABELS = ("<10", "10-30", ">30")
TARGET_SEP_BINS = (1, 2)                         # 0.1-0.5 and 0.5-1 s
INFO_LIMITED_SEP_BINS = (0,)                     # < 0.1 s: reported, never counted
ISOLATED_SEP_BIN = 4                             # > 3 s, the reference of the loss
AVAIL_BEFORE, AVAIL_AFTER = 10.0, 60.0           # s around the origin that the data must cover
N_CELLS = len(SEP_LABELS) * len(SNR_LABELS) * len(STATION_COUNTS) * len(DIST_LABELS)


def bin_of(x, edges, right_closed=False):
    """bin index of x for increasing edges; values equal to an edge go to the upper bin (lower-closed bins), or to
    the lower one when right_closed. nan -> -1."""
    x = np.asarray(x, float)
    b = np.searchsorted(np.asarray(edges, float), x, side="left" if right_closed else "right")
    return np.where(np.isnan(x), -1, b)


def separation(arr_t, sta_idx=None):
    """median over stations of |t_P(e) - nearest arrival (P or S) of any other event| at the same station.
    arr_t (m, n, 2); stations restricted to sta_idx (per-event list of station indices) when given. inf when there is
    no other arrival; nan when the event has no P at any station of its subset."""
    arr_t = np.asarray(arr_t, float)
    m, n, _ = arr_t.shape
    out = np.full(m, np.nan)
    for e in range(m):
        st = np.arange(n) if sta_idx is None else np.asarray(sta_idx[e], int)
        d = []
        for s in st:
            tp = arr_t[e, s, 0]
            if np.isnan(tp):
                continue
            other = np.delete(arr_t[:, s, :], e, axis=0).ravel()
            other = other[~np.isnan(other)]
            d.append(np.abs(other - tp).min() if len(other) else np.inf)
        if d:
            out[e] = float(np.median(d))
    return out


def network_snr(snr, sta_idx=None):
    """median over the subset stations of the per-station ratio (card, axis 2)."""
    snr = np.asarray(snr, float)
    if sta_idx is None:
        return np.nanmedian(snr, axis=1) if snr.size else np.zeros(0)
    out = np.full(len(snr), np.nan)
    for e, st in enumerate(sta_idx):
        v = snr[e, np.asarray(st, int)]
        if np.isfinite(v).any():
            out[e] = float(np.nanmedian(v))
    return out


def centroid(xyz):
    """epicentral centroid (x, y) of a development catalogue."""
    return np.asarray(xyz, float)[:, :2].mean(axis=0)


def station_subset(sta_xyz, codes, available, k, net_centroid=None):
    """the k stations nearest the network centroid among the available ones, ties broken by station code; None when
    fewer than k are available. available (n,) bool: data cover [origin - 10 s, origin + 60 s] without a gap."""
    sta_xyz = np.asarray(sta_xyz, float)
    c = sta_xyz[:, :2].mean(axis=0) if net_centroid is None else np.asarray(net_centroid, float)
    idx = np.flatnonzero(np.asarray(available, bool))
    if len(idx) < k:
        return None
    d = np.hypot(*(sta_xyz[idx, :2] - c).T)
    order = sorted(range(len(idx)), key=lambda i: (round(float(d[i]), 9), str(codes[idx[i]])))
    return np.sort(idx[order[:k]])


def availability(origins, coverage):
    """(m, n) bool: station n has data over [t - 10 s, t + 60 s] for every origin t. coverage: per station a list of
    gap-free intervals (t_start, t_end)."""
    origins = np.asarray(origins, float)
    out = np.zeros((len(origins), len(coverage)), bool)
    for s, ivs in enumerate(coverage):
        for a, b in ivs:
            out[:, s] |= (origins - AVAIL_BEFORE >= a) & (origins + AVAIL_AFTER <= b)
    return out


def cell_index(sep_bin, snr_bin, k_bin, dist_bin):
    """flat cell id; -1 when any axis is undefined."""
    sep_bin, snr_bin, k_bin, dist_bin = map(np.asarray, (sep_bin, snr_bin, k_bin, dist_bin))
    bad = (sep_bin < 0) | (snr_bin < 0) | (k_bin < 0) | (dist_bin < 0)
    i = ((sep_bin * len(SNR_LABELS) + snr_bin) * len(STATION_COUNTS) + k_bin) * len(DIST_LABELS) + dist_bin
    return np.where(bad, -1, i)


def cell_axes(i):
    """(sep_bin, snr_bin, k_bin, dist_bin) of a flat cell id."""
    i = int(i)
    i, d = divmod(i, len(DIST_LABELS))
    i, k = divmod(i, len(STATION_COUNTS))
    s, q = divmod(i, len(SNR_LABELS))
    return s, q, k, d


def cell_name(i):
    s, q, k, d = cell_axes(i)
    return f"sep {SEP_LABELS[s]} s | SNR {SNR_LABELS[q]} | {STATION_COUNTS[k]} sta | dist {DIST_LABELS[d]} km"


def isolated_cell(i):
    """the matched isolated cell of a cell (same SNR, stations and distance, separation > 3 s)."""
    _, q, k, d = cell_axes(i)
    return int(cell_index(ISOLATED_SEP_BIN, q, k, d))


def cell_kind(i):
    """"target", "info" (information-limited, reported only) or "other"."""
    s = cell_axes(i)[0]
    if s in TARGET_SEP_BINS:
        return "target"
    if s in INFO_LIMITED_SEP_BINS:
        return "info"
    return "other"


def assign(truth, k, subsets, net_centroid):
    """cells of the events of truth for subset size k. subsets: per event the station indices of its subset (None
    when the network has fewer than k available stations). Returns (cells (m,), sep (m,), snr (m,), dist (m,))."""
    m = truth.n_events
    kb = STATION_COUNTS.index(k)
    ok = np.array([s is not None for s in subsets], bool)
    st = [s if s is not None else np.zeros(0, int) for s in subsets]
    sep = separation(truth.arr_t, st)
    snr = network_snr(truth.snr, st)
    dist = np.hypot(*(truth.xyz[:, :2] - np.asarray(net_centroid, float)).T)
    cells = cell_index(bin_of(sep, SEP_EDGES), bin_of(snr, SNR_EDGES), np.full(m, kb), bin_of(dist, DIST_EDGES))
    cells = np.where(ok & truth.target, cells, -1)
    return cells, sep, snr, dist
