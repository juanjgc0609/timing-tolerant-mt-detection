"""Internal diagnostics of the atlas (card, "Instrumentation"): where each pipeline loses an event.

Picks against true arrivals, per station and phase: picked, merged, wrong phase or lost, against the separation to the
nearest other arrival at the station. Events: merged and split detections traced through the pick matching, and
misses attributed to the picker or to the associator. Location: LOC hits that fail LOC+T."""
import numpy as np

from .. import stats as S
from .cells import SEP_EDGES, bin_of

PICK_TOL = 0.1                  # s, primary pick tolerance (secondary 0.2)
MERGE_WINDOW = 1.0              # s, a neighbouring arrival of the same phase within it can absorb the pick
MIN_EVENT_PICKS = 3             # picks a true event must contribute to count in a merge or a split
PICKED, MERGED, WRONG_PHASE, LOST = 0, 1, 2, 3
PICK_STATES = ("picked", "merged", "wrong_phase", "lost")


def match_picks(truth, picks, tol=PICK_TOL, observed=True):
    """one-to-one matching of picks to true arrivals at the same station and phase, nearest first.
    Returns arr_pick (m, n, 2): index of the pick of each arrival or -1, and pick_ev (k,): the event of each pick or
    -1. observed: match against arrival times with the injected timing error (what a picker can see)."""
    at = truth.arr_t_obs if observed else truth.arr_t
    m, n, _ = at.shape
    arr_pick = -np.ones((m, n, 2), int)
    pick_ev = -np.ones(len(picks), int)
    for s in range(n):
        for ph in (0, 1):
            pk = np.flatnonzero((picks.sta == s) & (picks.phase == ph))
            ev = np.flatnonzero(~np.isnan(at[:, s, ph]))
            if not len(pk) or not len(ev):
                continue
            d = np.abs(picks.t[pk][:, None] - at[ev, s, ph][None, :])
            pairs = np.argwhere(d <= tol + 1e-12)
            for i, j in sorted(pairs.tolist(), key=lambda q: (d[q[0], q[1]], q[0], q[1])):
                if pick_ev[pk[i]] < 0 and arr_pick[ev[j], s, ph] < 0:
                    pick_ev[pk[i]] = ev[j]
                    arr_pick[ev[j], s, ph] = pk[i]
    return arr_pick, pick_ev


def arrival_states(truth, picks, tol=PICK_TOL, observed=True):
    """state (m, n, 2) of every true arrival (PICKED, MERGED, WRONG_PHASE, LOST; -1 where the arrival is absent) and
    the separation (m, n, 2) to the nearest other arrival at the same station (any phase, any other event)."""
    at = truth.arr_t_obs if observed else truth.arr_t
    m, n, _ = at.shape
    arr_pick, pick_ev = match_picks(truth, picks, tol, observed)
    state = np.full((m, n, 2), -1, int)
    sep = np.full((m, n, 2), np.nan)
    for e in range(m):
        for s in range(n):
            for ph in (0, 1):
                t = at[e, s, ph]
                if np.isnan(t):
                    continue
                other = np.delete(at[:, s, :], e, axis=0).ravel()
                other = other[~np.isnan(other)]
                sep[e, s, ph] = np.abs(other - t).min() if len(other) else np.inf
                if arr_pick[e, s, ph] >= 0:
                    state[e, s, ph] = PICKED
                    continue
                state[e, s, ph] = LOST
                same = (picks.sta == s) & (picks.phase == ph)
                # merged: a pick of this phase, matched to another arrival within MERGE_WINDOW, lies within the span
                # of the two arrivals widened by the tolerance
                for k in np.flatnonzero(same & (pick_ev >= 0) & (pick_ev != e)):
                    tb = at[pick_ev[k], s, ph]
                    if abs(tb - t) <= MERGE_WINDOW and min(t, tb) - tol <= picks.t[k] <= max(t, tb) + tol:
                        state[e, s, ph] = MERGED
                        break
                if state[e, s, ph] == LOST:
                    near = (picks.sta == s) & (picks.phase == 1 - ph) & (np.abs(picks.t - t) <= tol + 1e-12)
                    if near.any():
                        state[e, s, ph] = WRONG_PHASE
    return state, sep


def pick_table(state, sep, edges=SEP_EDGES):
    """fractions of the four pick states per separation bin: (n_bins, 4) and the counts (n_bins,)."""
    ok = state >= 0
    b = bin_of(sep[ok], edges)
    st = state[ok]
    nb = len(edges) + 1
    cnt = np.zeros((nb, 4))
    np.add.at(cnt, (b, st), 1)
    tot = cnt.sum(axis=1)
    return np.divide(cnt, tot[:, None], out=np.full_like(cnt, np.nan), where=tot[:, None] > 0), tot


def event_diagnostics(truth, picks, det_labels, n_det, min_picks, tol=PICK_TOL):
    """event-level diagnostics of a pick-based pipeline.
    det_labels (d,): LOC+T-V labels of the detections (ttmt.atlas.froc.label); picks.det: detection of each pick.
    Returns per detection merged (d,) bool, per event split (m,) bool, and per missed target its cause (m,) as
    "hit", "picker" (fewer of its arrivals picked than min_picks) or "associator"."""
    _, pick_ev = match_picks(truth, picks, tol)
    m = truth.n_events
    contrib = np.zeros((n_det, m), int)
    k = (picks.det >= 0) & (pick_ev >= 0)
    np.add.at(contrib, (picks.det[k], pick_ev[k]), 1)
    merged = (contrib >= MIN_EVENT_PICKS).sum(axis=1) >= 2
    split = (contrib >= MIN_EVENT_PICKS).sum(axis=0) >= 2
    npick = np.bincount(pick_ev[pick_ev >= 0], minlength=m)
    hit = np.zeros(m, bool)
    hit[det_labels[det_labels >= 0]] = True
    cause = np.where(hit, "hit", np.where(npick < min_picks, "picker", "associator")).astype(object)
    cause[~truth.target] = ""
    return merged, split, cause


def waveform_merges(det, det_labels, truth, dt=2.0, dkm=5.0):
    """for methods without picks: a missed target is "merged" when a detection hits another event within dt s and
    dkm km of it."""
    hit = np.zeros(truth.n_events, bool)
    hit[det_labels[det_labels >= 0]] = True
    out = np.zeros(truth.n_events, bool)
    for e in np.flatnonzero(truth.target & ~hit):
        d = np.hypot(*(truth.xyz[hit, :2] - truth.xyz[e, :2]).T)
        out[e] = bool(((np.abs(truth.t[hit] - truth.t[e]) <= dt) & (d <= dkm)).any())
    return out


def location_failures(det, truth):
    """fraction of LOC hits (epicentre within 5 km, any origin, one to one) that fail LOC+T, with the origin and
    epicentre errors of the LOC hits. Uses the nearest-epicentre one-to-one rule of ADR-0007 with no time test."""
    from ..scoring import one_to_one
    if len(det) == 0 or truth.n_events == 0:
        return dict(frac_fail=float("nan"), dt=np.zeros(0), dx=np.zeros(0))
    dx = np.hypot(det.xyz[:, None, 0] - truth.xyz[None, :, 0], det.xyz[:, None, 1] - truth.xyz[None, :, 1])
    lab = one_to_one(det.s, dx, dx <= S.TOL_KM, order="score")
    ok = lab >= 0
    err_t = np.abs(det.t[ok] - truth.t[lab[ok]])
    err_x = dx[np.flatnonzero(ok), lab[ok]]
    return dict(frac_fail=float(np.mean(err_t > S.TOL_S + 1e-9)) if ok.any() else float("nan"), dt=err_t, dx=err_x)
