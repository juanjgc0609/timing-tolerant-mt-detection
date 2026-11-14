"""Static station terms (T4.24 / T4.25): per-station, per-phase delays of the observed arrivals relative to the 1-D
velocity model of the frame, applied to picks before association. PyOcto 0.2.0 and GaMMA take no station terms in
their API, so correcting the picks (t - delta[station, phase]) is the equivalent, and the same correction is given to
every comparator and to our model.

Estimation without leakage: the terms come from the injected sources' own arrival times (what the network operator
knows from its catalogue), and no target ever informs its own correction: the terms of a record are estimated from
the sources that do not appear in that record (cross-fitting by record), and one term set is applied to all its
picks. A source is identified by its catalogue position (EGF sources are fixed catalogue events). Residuals are
relative to the catalogue origin, as an operator's terms are; a constant over stations trades off with the origin
time, so remove_origin=True (per-source median over stations and phases removed) is kept as a sensitivity variant.
"""
import numpy as np

from ..atlas.schema import Picks


def source_keys(xyz, decimals=3):
    """integer label of each event's source: events at the same catalogue position share a label."""
    r = np.round(np.asarray(xyz, float), decimals)
    _, lab = np.unique(r, axis=0, return_inverse=True)
    return lab.ravel()


def residuals(truth, stations_xyz, velocity, events=None, observed=True):
    """(n_events, ns, 2) arrival minus origin minus 1-D travel time, nan where absent; events: indices to use."""
    arr = truth.arr_t_obs if (observed and truth.arr_t_obs is not None) else truth.arr_t
    ev = np.arange(truth.n_events) if events is None else np.asarray(events, int)
    out = np.full((truth.n_events, arr.shape[1], 2), np.nan)
    sta = np.asarray(stations_xyz, float)
    for e in ev:
        if not np.isfinite(arr[e]).any():
            continue
        tt = np.asarray(velocity.times(np.asarray(truth.xyz[e], float), sta))  # (2, ns)
        out[e] = arr[e] - truth.t[e] - tt.T
    return out


def per_source(res, keys, remove_origin=False):
    """(n_sources, ns, 2) median residual of each source (its repeated injections agree up to timing noise); with
    remove_origin the source's median over stations and phases is removed."""
    ks = np.unique(keys)
    out = np.full((len(ks), res.shape[1], 2), np.nan)
    for i, k in enumerate(ks):
        m = keys == k
        with np.errstate(all="ignore"):
            r = np.nanmedian(res[m], axis=0) if np.isfinite(res[m]).any() else np.full(res.shape[1:], np.nan)
        out[i] = r - np.nanmedian(r) if (remove_origin and np.isfinite(r).any()) else r
    return ks, out


def terms(src_res, exclude=None, min_sources=3):
    """(ns, 2) station terms: median over sources of the per-source residuals, leaving out source index exclude;
    0 where fewer than min_sources sources see the station and phase."""
    use = np.ones(len(src_res), bool)
    if exclude is not None:
        use[exclude] = False
    r = src_res[use]
    n = np.isfinite(r).sum(axis=0)
    with np.errstate(all="ignore"):
        d = np.nanmedian(r, axis=0) if len(r) else np.zeros(src_res.shape[1:])
    return np.where((n >= min_sources) & np.isfinite(d), d, 0.0)


def apply(picks, delta):
    """picks with t - delta[station, phase]."""
    return Picks(picks.sta, picks.phase, picks.t - delta[picks.sta, picks.phase], picks.p)


def record_terms(src_res, ks, keys_in_record, min_sources=3):
    """(ns, 2) terms for one record, estimated from every source that does not appear in the record (cross-fitting by
    record: one term set per record, as an operator would apply it, and no target informs its own correction)."""
    excl = np.isin(ks, np.unique(keys_in_record))
    return terms(src_res, exclude=np.where(excl)[0], min_sources=min_sources)


def fill_missing(delta):
    """terms with the stations and phases that have no estimate (0) set to the median of the estimated ones of the
    same phase: a network-wide shift trades off with the origin time, but a shift applied to some stations only is a
    relative error between stations, so a missing term must take the common shift and not zero."""
    d = np.array(delta, float)
    for b in range(d.shape[1]):
        have = d[:, b] != 0
        if have.any():
            d[~have, b] = np.median(d[have, b])
    return d
