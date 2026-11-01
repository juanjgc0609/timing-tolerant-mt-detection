"""Scoring of the atlas under the single rule of ADR-0008.

Labels: LOC+T-V (ttmt.stats.label_voronoi) against every catalogued event of the record, targets and non-targets;
a detection matched to a non-target (a catalogued event of the noise record, or an event outside the station subset
of the run) is neutral. Operating point: exact FROC at a false-detection budget, chosen per method on the evaluated
data (ADR-0008, item 4); in the block bootstrap it is re-derived in every replicate from the resampled blocks.
Decisions: the fixed sequence H1 (superiority in T) then H2 (non-inferiority in C, margin 0.05), each an
intersection-union test over the comparators of block sign-flip tests (ttmt.stats). Cells are described only.
"""
import warnings
from dataclasses import dataclass, field

import numpy as np

from .. import stats as S
from . import cells as C
from . import seeds as Z

FALSE, NEUTRAL = S.FALSE, S.NEUTRAL
BUDGETS = (30.0, 10.0, 100.0)            # false detections per hour; the first is primary
ALPHA = 0.025                            # one-sided, ADR-0008 items 5 and 7
NI_MARGIN = 0.05                         # H2
NI_MARGIN_SECONDARY = 0.03
MIN_BLOCKS = 30                          # below it an interval is indicative (SAP S7.5)
TRUTH_PRECISION_FACTOR = 9.0             # a real-data pair enters T only if its separation >= 9 sigma_t (SAP 3.4)


def label(det, truth, target=None):
    """LOC+T-V labels (d,) of a Detections object against truth: pooled-event index of a target hit, NEUTRAL when the
    matched event is not a target, FALSE otherwise."""
    target = truth.target if target is None else np.asarray(target, bool)
    lab = S.label_voronoi(det.s, det.xyz[:, :2], det.t, truth.xyz[:, :2], truth.t)
    return np.where((lab >= 0) & ~target[np.maximum(lab, 0)], NEUTRAL, lab)


def duplicate(det, half_sep):
    """the duplicator control of ADR-0008 item 3, as a Detections object; half_sep scalar or one value per
    detection (duplicator_shift)."""
    from .schema import Detections
    s, xy, t = S.duplicator(det.s, det.xyz[:, :2], det.t, half_sep)
    xyz = np.c_[xy, np.r_[det.xyz[:, 2], det.xyz[:, 2]]]
    return Detections(t, xyz, s, dict(det.meta, control="duplicator"))


def duplicator_shift(det, truth, sep, cap=1.0):
    """half-shift (d,) of the duplicator of ADR-0008 item 3: for every detection, half the separation (cells.separation,
    the arrival separation that defines T) of the catalogued event nearest to it in the normalised cost of
    ttmt.stats.cost_matrix. A copy then lands on the bisector of that event's pair, where LOC+T-V never makes it
    eligible; one shift per record would land copies on partners at other separations (review of PR #12).
    Separations above cap (isolated events, or none) use cap: such copies cannot reach another event."""
    if len(det) == 0 or truth.n_events == 0:
        return np.zeros(len(det))
    c = S.cost_matrix(det.xyz[:, :2], det.t, truth.xyz[:, :2], truth.t)
    near = np.argmin(c, axis=1)
    sp = np.asarray(sep, float)[near]
    return 0.5 * np.where(np.isfinite(sp), np.minimum(sp, cap), cap)


def isolated_duplicates(labels, det, truth, regime):
    """number of false detections within LOC+T of a hit isolated target (separation > 3 s): the per-method duplicate
    count of isolated events that ADR-0008 item 3 asks to report."""
    hit = np.zeros(truth.n_events, bool)
    hit[labels[labels >= 0]] = True
    iso = hit & (regime == "iso")
    if not iso.any() or len(det) == 0:
        return 0
    c = S.cost_matrix(det.xyz[:, :2], det.t, truth.xyz[iso, :2], truth.t[iso])
    return int(((labels == FALSE) & (c <= 1.0).any(axis=1)).sum())


def regime_of(sep):
    """"T" for separation in [0.1, 1) s, "C" otherwise (ADR-0008 item 1); "iso" marks the C events above 3 s, the
    reference of the atlas losses. nan separation (no P in the subset) -> ""."""
    sep = np.asarray(sep, float)
    out = np.where((sep >= 0.1) & (sep < 1.0), "T", "C").astype(object)
    out[sep > 3.0] = "iso"
    out[np.isnan(sep)] = ""
    return out


@dataclass
class Pool:
    """every record of one evaluation, pooled. Events (E,): cell, regime ("T", "C", "iso" or "" for none), block,
    stratum (sequence or network, for the stratified bootstrap), weight (prevalence weight, ADR-0008 item 4).
    block_hours (nb,); block_cluster (nb,) optional, the unit of the sensitivity bootstrap (_resampler). Per method:
    scores, labels into the pooled events (or FALSE / NEUTRAL), detection blocks."""
    block_hours: np.ndarray
    block_stratum: np.ndarray
    ev_cell: np.ndarray = field(default_factory=lambda: np.zeros(0, int))
    ev_regime: np.ndarray = field(default_factory=lambda: np.zeros(0, object))
    ev_block: np.ndarray = field(default_factory=lambda: np.zeros(0, int))
    ev_weight: np.ndarray = field(default_factory=lambda: np.zeros(0))
    runs: dict = field(default_factory=dict)
    block_cluster: np.ndarray = None

    def __post_init__(self):
        self.block_hours = np.asarray(self.block_hours, float)
        self.block_stratum = np.asarray(self.block_stratum, int)

    @property
    def n_blocks(self):
        return len(self.block_hours)

    def add_record(self, cells, regime, ev_block, method_runs, weights=None):
        """append one record: its target cells, regimes and blocks, and per method (scores, labels into this
        record's events, detection blocks)."""
        off = len(self.ev_cell)
        n = len(cells)
        self.ev_cell = np.r_[self.ev_cell, np.asarray(cells, int)]
        self.ev_regime = np.r_[self.ev_regime, np.asarray(regime, object)]
        self.ev_block = np.r_[self.ev_block, np.asarray(ev_block, int)]
        self.ev_weight = np.r_[self.ev_weight, np.ones(n) if weights is None else np.asarray(weights, float)]
        for name, (s, lab, blk) in method_runs.items():
            lab = np.asarray(lab, int)
            lab = np.where(lab >= 0, lab + off, lab)
            r = self.runs.setdefault(name, [np.zeros(0), np.zeros(0, int), np.zeros(0, int)])
            r[0], r[1], r[2] = np.r_[r[0], np.asarray(s, float)], np.r_[r[1], lab], np.r_[r[2], np.asarray(blk, int)]


def truth_precision_targets(truth, sigma_t, factor=TRUTH_PRECISION_FACTOR):
    """target flags after the truth-precision rule of the phase-6 plan (SAP 3.4): with real truth whose origin-time
    standard deviation is sigma_t, an event closer in origin time than factor * sigma_t to another catalogued event
    is excluded (it becomes a non-target, so detections of it are neutral). sigma_t None: unchanged (semi-synthetic
    truth is exact)."""
    if sigma_t is None or truth.n_events < 2:
        return truth.target.copy()
    o = np.argsort(truth.t, kind="stable")
    gap = np.diff(truth.t[o])
    near = np.empty(truth.n_events)
    near[o] = np.minimum(np.r_[np.inf, gap], np.r_[gap, np.inf])
    return truth.target & (near >= factor * sigma_t)


def score_record(pool, truth, k, subsets, net_centroid, ev_block, det_block, runs, duplicator_of=None,
                 weights=None, sigma_t=None):
    """label every method's detections of one record and add them to the pool. k: subset size; subsets: per event
    its station indices (cells.station_subset); ev_block (m,) and det_block: callable t -> block id; runs: name ->
    Detections. duplicator_of: name of the comparator whose duplicator is added as "dup:<name>", with half_sep half
    the separation of the event nearest to each detection (duplicator_shift; ADR-0008 item 3). sigma_t: origin-time
    standard deviation of a real truth catalogue (measured on T3.V); pairs closer than 9 sigma_t are excluded."""
    if sigma_t is not None:
        from dataclasses import replace
        truth = replace(truth, target=truth_precision_targets(truth, sigma_t))
    cells, sep, _, _ = C.assign(truth, k, subsets, net_centroid)
    reg = regime_of(sep)
    reg[~truth.target] = ""
    if duplicator_of is not None:
        det = runs[duplicator_of]
        runs = dict(runs, **{"dup:" + duplicator_of: duplicate(det, duplicator_shift(det, truth, sep))})
    entry = {}
    for name, det in runs.items():
        entry[name] = (det.s, label(det, truth), np.array([det_block(x) for x in det.t], int))
    pool.add_record(cells, reg, ev_block, entry, weights)
    return cells, reg


def operating_point(scores, labels, mult, hours, budget):
    """ttmt.stats.operating_point with a multiplicity per detection (its block's count in a bootstrap replicate):
    the lowest threshold whose detections (score >= threshold) hold at most floor(budget * hours) weighted false
    detections, tied scores together. +inf when no detection can be emitted."""
    cap = np.floor(budget * hours + 1e-9)
    if len(scores) == 0:
        return np.inf
    o = np.argsort(-scores, kind="stable")
    s = scores[o]
    f = np.cumsum((labels[o] == FALSE) * mult[o])
    last = np.r_[s[1:] != s[:-1], True]
    ok = last & (f <= cap + 1e-9)
    return float(s[np.flatnonzero(ok)[-1]]) if ok.any() else np.inf


def hits(pool, method, budget, mult=None):
    """bool (E,): pooled events hit at the method's operating point; mult (nb,) block multiplicities (default 1)."""
    mult = np.ones(pool.n_blocks) if mult is None else np.asarray(mult, float)
    s, lab, blk = pool.runs[method]
    m = mult[blk] if len(blk) else np.zeros(0)
    keep = m > 0
    hours = float((mult * pool.block_hours).sum())
    thr = operating_point(s[keep], lab[keep], m[keep], hours, budget)
    got = np.zeros(len(pool.ev_cell), bool)
    sel = (s >= thr) & (lab >= 0)
    got[lab[sel]] = True
    return got


def _multiplicity(pool, idx):
    return np.bincount(np.asarray(idx, int), minlength=pool.n_blocks).astype(float)


def _strata(pool):
    return [np.flatnonzero(pool.block_stratum == k) for k in np.unique(pool.block_stratum)]


def _resampler(pool, cluster=None):
    """draw() -> block multiplicities (nb,) of one stratified bootstrap replicate. cluster None: the registered unit,
    blocks drawn with replacement within each stratum. cluster (nb,) ints: a sensitivity unit, clusters (e.g. the
    noise block shared by several realisations) drawn with replacement within each stratum, each drawn cluster
    bringing all its blocks. Strata that share a cluster (e.g. two sets injected into the same noise) are drawn
    together, so a shared cluster is never split."""
    strata = _strata(pool)
    if cluster is None:
        return lambda rng: _multiplicity(pool, np.concatenate([s[rng.integers(0, s.size, s.size)] for s in strata]))
    cluster = np.asarray(cluster, int)
    if cluster.shape != (pool.n_blocks,):
        raise ValueError(f"cluster must give one id per block ({pool.n_blocks}), got shape {cluster.shape}")
    lab = np.arange(len(strata))                                # merge strata linked by a shared cluster
    owner = {}
    for i, s in enumerate(strata):
        for c in np.unique(cluster[s]):
            j = owner.setdefault(int(c), i)
            a, b = lab[i], lab[j]
            lab[lab == max(a, b)] = min(a, b)
    groups = []
    for g in np.unique(lab):
        s = np.concatenate([strata[i] for i in np.flatnonzero(lab == g)])
        ids = np.unique(cluster[s])
        groups.append((ids, [s[cluster[s] == c] for c in ids]))

    def draw(rng):
        mu = np.zeros(pool.n_blocks)
        for ids, members in groups:
            for j in rng.integers(0, ids.size, ids.size):
                mu[members[j]] += 1.0
        return mu
    return draw


def n_clusters(pool, cluster):
    """number of distinct clusters (e.g. noise blocks) of the pool."""
    return int(np.unique(np.asarray(cluster, int)).size)


def recall(pool, got, mask, mult=None):
    """weighted recall of the events of mask (bool (E,)), with block multiplicities."""
    mult = np.ones(pool.n_blocks) if mult is None else mult
    w = pool.ev_weight * mult[pool.ev_block] * mask
    return float((w * got).sum() / w.sum()) if w.sum() > 0 else float("nan")


def cell_recalls(pool, got, mult=None):
    """weighted recall of every cell (N_CELLS,), nan for empty cells."""
    mult = np.ones(pool.n_blocks) if mult is None else mult
    ok = pool.ev_cell >= 0
    w = pool.ev_weight * mult[pool.ev_block] * ok
    c = np.where(ok, pool.ev_cell, 0)
    n = np.bincount(c, weights=w, minlength=C.N_CELLS)
    h = np.bincount(c, weights=w * got, minlength=C.N_CELLS)
    return np.divide(h, n, out=np.full(C.N_CELLS, np.nan), where=n > 0)


def _gen(seed, stream, index):
    """an integer seed gives its own generator; None gives the anchored stream (ttmt.atlas.seeds)."""
    return Z.rng(stream, index) if seed is None else np.random.default_rng(seed)


def describe(pool, methods, budget=BUDGETS[0], n_boot=2000, seed=None, cluster=None):
    """the descriptive tables of the atlas: per method, the recall of every cell and of the regimes T, C and iso
    with stratified block-bootstrap 95 % intervals (operating points re-derived per replicate). No test.
    seed None: the anchored stream ("bootstrap", 0), or ("bootstrap", 2) with a sensitivity cluster (_resampler)."""
    rng = _gen(seed, "bootstrap", 0 if cluster is None else 2)
    draw = _resampler(pool, cluster)
    regmask = {r: pool.ev_regime == r for r in ("T", "C", "iso")}
    out = {}
    for m in methods:
        got = hits(pool, m, budget)
        pt_c = cell_recalls(pool, got)
        pt_r = {r: recall(pool, got, mk) for r, mk in regmask.items()}
        rc, rr = np.empty((n_boot, C.N_CELLS)), {r: np.empty(n_boot) for r in regmask}
        for b in range(n_boot):
            mu = draw(rng)
            g = hits(pool, m, budget, mu)
            rc[b] = cell_recalls(pool, g, mu)
            for r, mk in regmask.items():
                rr[r][b] = recall(pool, g, mk, mu)
        n_cell = np.bincount(pool.ev_cell[pool.ev_cell >= 0], minlength=C.N_CELLS)
        with warnings.catch_warnings():                 # empty cells have no replicate: their bounds stay nan
            warnings.simplefilter("ignore", RuntimeWarning)
            lo, hi = np.nanquantile(rc, 0.025, axis=0), np.nanquantile(rc, 0.975, axis=0)
        out[m] = dict(cell=pt_c, cell_lo=lo, cell_hi=hi,
                      n_cell=n_cell, regime={r: (pt_r[r], float(np.nanquantile(rr[r], 0.025)),
                                                 float(np.nanquantile(rr[r], 0.975)), int(regmask[r].sum()))
                                             for r in regmask},
                      indicative=pool.n_blocks < MIN_BLOCKS)
    return out


def _pairwise(pool, a, b, mask, budget, margin, signs):
    ga, gb = hits(pool, a, budget), hits(pool, b, budget)
    w = pool.ev_weight * mask
    d, n = S.block_sums(ga, gb, pool.ev_block, pool.n_blocks, w)
    if margin:      # non-inferiority: clustered score test (ADR-0008 item 7), not the margin sign-flip
        bb, cc, nn = S.block_counts(ga, gb, pool.ev_block, pool.n_blocks, w)
        p = S.ni_score_pvalue(bb, cc, nn, margin)
    elif signs is not None:
        p = S.signflip_pvalue(d, n, signs=signs)
    else:
        p = S.signflip_pvalue(d, n, rng=Z.rng("signflip", 0))
    diff = S.paired_difference(d, n) if n.sum() > 0 else float("nan")
    return dict(diff=diff, p=float(p), n=float(n.sum()), blocks=int((n > 0).sum()))


def gate(pool, model, comparators, budget=BUDGETS[0], duplicator_of=None, signs=None, alpha_h2=ALPHA):
    """the rule of ADR-0008 item 7 on a pooled evaluation.
    H1: superiority in T against every comparator (IUT of sign-flip p-values); H2, only if H1 is rejected:
    non-inferiority in C at margin 0.05. Secondary family (only if H1 and H2): superiority in [0.1, 0.5) s, in
    [0.5, 1) s, non-inferiority in C at 0.03, Holm at 0.025. duplicator_of: (comparator, its duplicator run name) for
    the negative control of item 3; if the duplicator beats its comparator in T at 0.025 the metric is broken and
    no claim is made. alpha_h2: level of H2, 0.025 unless the registered design check (sample_size.py calibrate)
    lowered it; the Holm threshold of the secondary non-inferiority is scaled by the same factor."""
    T = pool.ev_regime == "T"
    Cm = (pool.ev_regime == "C") | (pool.ev_regime == "iso")
    sep_bin = np.array([C.cell_axes(c)[0] if c >= 0 else -1 for c in pool.ev_cell])
    res = dict(budget=budget, comparators=list(comparators))
    if duplicator_of is not None:
        comp, dup = duplicator_of
        ctl = _pairwise(pool, dup, comp, T, budget, 0.0, signs)
        res["control"] = dict(ctl, broken=bool(ctl["p"] <= ALPHA))
    h1 = {c: _pairwise(pool, model, c, T, budget, 0.0, signs) for c in comparators}
    h2 = {c: _pairwise(pool, model, c, Cm, budget, NI_MARGIN, signs) for c in comparators}
    p1, p2 = S.iut([r["p"] for r in h1.values()]), S.iut([r["p"] for r in h2.values()])
    rej1 = bool(p1 <= ALPHA)
    rej2 = bool(rej1 and p2 <= alpha_h2)
    res.update(H1=dict(p=p1, rows=h1, rejected=rej1), H2=dict(p=p2, rows=h2, rejected=rej2))
    if rej1 and rej2:
        fam = {"T[0.1,0.5)": (T & (sep_bin == 1), 0.0), "T[0.5,1)": (T & (sep_bin == 2), 0.0),
               "C NI 0.03": (Cm, NI_MARGIN_SECONDARY)}
        ps = {k: S.iut([_pairwise(pool, model, c, mk, budget, mg, signs)["p"] for c in comparators])
              for k, (mk, mg) in fam.items()}
        ps["C NI 0.03"] = min(1.0, ps["C NI 0.03"] * ALPHA / alpha_h2)
        rej = S.holm(list(ps.values()), ALPHA)
        res["secondary"] = {k: dict(p=p, rejected=r) for (k, p), r in zip(ps.items(), rej)}
    broken = res.get("control", {}).get("broken", False)
    res["beats"] = bool(rej1 and rej2 and not broken)
    res["indicative"] = pool.n_blocks < MIN_BLOCKS
    return res


def pooled_loss(pool, method, budget=BUDGETS[0], n_boot=2000, seed=None, mask=None, cluster=None):
    """the loss of the atlas hypothesis (card, H-ATLAS): recall of the isolated events (separation > 3 s),
    standardised to the (SNR, stations, distance) mix of the events of mask (default regime T), minus the recall of
    the events of mask; point and stratified block-bootstrap 95 % interval. seed None: the anchored stream
    ("bootstrap", 1), or ("bootstrap", 3) with a sensitivity cluster (_resampler)."""
    mask = (pool.ev_regime == "T") if mask is None else mask
    iso = pool.ev_regime == "iso"
    key = np.array([C.cell_axes(c)[1:] if c >= 0 else (-1, -1, -1) for c in pool.ev_cell])
    kid = key[:, 0] * 100 + key[:, 1] * 10 + key[:, 2]

    def stat(mu):
        g = hits(pool, method, budget, mu)
        w = pool.ev_weight * mu[pool.ev_block]
        wt = w * mask
        r_t = (wt * g).sum() / wt.sum() if wt.sum() > 0 else np.nan
        acc = tot = 0.0
        for k in np.unique(kid[mask]):
            share = wt[kid == k].sum()
            wi = w * iso * (kid == k)
            if wi.sum() > 0 and share > 0:
                acc += share * (wi * g).sum() / wi.sum()
                tot += share
        r_iso = acc / tot if tot > 0 else np.nan
        return r_iso - r_t

    rng = _gen(seed, "bootstrap", 1 if cluster is None else 3)
    draw = _resampler(pool, cluster)
    pt = stat(np.ones(pool.n_blocks))
    reps = np.array([stat(draw(rng)) for _ in range(n_boot)])
    return dict(loss=float(pt), lo=float(np.nanquantile(reps, 0.025)), hi=float(np.nanquantile(reps, 0.975)),
                indicative=pool.n_blocks < MIN_BLOCKS)
