"""Decision rules of the single "beat" rule (ADR-0008) and of the phase 6 analysis plan (docs/protocol/PHASE6_SAP.md).

Every decision rule of the plan is a function here, so that the power simulation
(experiments/exploratory/P6-SAP/sample_size.py) runs the same code the registered analysis will run. numpy only.

Input convention. Detections of each method are labelled once, on the full list, by the LOC+T-V rule of ADR-0008
(label_voronoi: ADR-0007 matching in decreasing score order, restricted to the event each detection is closest to).
That matching is prefix-stable: the label of a detection
depends only on detections with a higher score, so the labels of the full list are the labels at every threshold.
A label is a target index (hit), -1 (false), or -2 (neutral: matched to an excluded event or inside an excluded
zone, SAP Section 5).

Sections of the plan implemented here:
  S6   LOC+T-V hit rule and the duplicator negative control           label_voronoi, duplicator
  S7.1 operating point at a false-detection budget (exact FROC)        operating_point, hits_at
  S7.2 paired block statistics                                         block_sums
  S7.3 block sign-flip test (exact randomisation form of the clustered
       McNemar test of Durkalski et al. 2003)                          signflip_pvalue, clustered_mcnemar_z
  S7.4 intersection-union over comparators, fixed sequence, Holm        iut, fixed_sequence, holm
  S7.5 stratified block bootstrap interval with atom diagnostics       block_bootstrap
  S7.6 false-alarm budget check                                        budget_check
"""
from math import comb, erf, sqrt
from statistics import NormalDist

import numpy as np

FALSE, NEUTRAL = -1, -2


# ---------------------------------------------------------------- S6 hit rule
TOL_KM, TOL_S = 5.0, 0.5
DOMINANCE = 0.5      # ADR-0008 item 2: c(d, e) <= DOMINANCE * c(d, e') for every other catalogued event e'


def cost_matrix(det_xy, det_t, ev_xy, ev_t, tol_km=TOL_KM, tol_s=TOL_S):
    """Normalised cost c(d, e) = max(epicentral distance / tol_km, |origin difference| / tol_s); LOC+T is c <= 1."""
    det_xy, ev_xy = np.atleast_2d(np.asarray(det_xy, float)), np.atleast_2d(np.asarray(ev_xy, float))
    dx = np.hypot(det_xy[:, None, 0] - ev_xy[None, :, 0], det_xy[:, None, 1] - ev_xy[None, :, 1])
    dt = np.abs(np.asarray(det_t, float)[:, None] - np.asarray(ev_t, float)[None, :])
    return np.maximum(dx / tol_km, dt / tol_s)


def label_voronoi(scores, det_xy, det_t, ev_xy, ev_t, tol_km=TOL_KM, tol_s=TOL_S, dominance=DOMINANCE):
    """LOC+T-V labels (ADR-0008, item 2). A detection is eligible for event e only if c(d, e) <= 1 (LOC+T) and
    c(d, e) <= dominance * c(d, e') for every other catalogued event e' (with dominance 0.5 and two events at the same
    epicentre, the detection must lie within a third of the separation of e). Eligible detections are matched one to
    one in decreasing score order (ttmt.scoring.one_to_one, order "score"); the rest are false (-1). Two detections of
    one event then cannot be counted as the two events of a close pair, which LOC+T alone allows when the tolerance
    exceeds the separation (CR-02). Plain nearest-event (Voronoi) eligibility, dominance 1, is not enough: a copy
    of a detection moved by half the separation lands on the bisector and falls on the partner's side half of the time
    (tests/test_stats.py, duplicator control). For an isolated event the rule is LOC+T. Prefix-stable, so labels of the
    full list hold at every threshold."""
    from ttmt.scoring import one_to_one
    n = len(scores)
    if n == 0 or len(ev_t) == 0:
        return np.full(n, FALSE, int)
    c = cost_matrix(det_xy, det_t, ev_xy, ev_t, tol_km, tol_s)
    if c.shape[1] > 1:
        part = np.partition(c, 1, axis=1)
        best, second = part[:, :1], part[:, 1:2]
        ok = (c <= 1.0 + 1e-12) & (c == best) & (best <= dominance * second + 1e-12)
    else:
        ok = c <= 1.0 + 1e-12
    return one_to_one(np.asarray(scores, float), c, ok, order="score")


def duplicator(scores, det_xy, det_t, half_sep):
    """Negative control of ADR-0008 (CR-02): every detection emitted twice, at origin -half_sep and +half_sep, with the
    same epicentre and score. A sound pair metric gives it no more pair recall than the original detector."""
    det_t = np.asarray(det_t, float)
    return (np.r_[scores, scores], np.r_[det_xy, det_xy], np.r_[det_t - half_sep, det_t + half_sep])


# ---------------------------------------------------------------- S7.1 operating point
def operating_point(scores, labels, hours, budget):
    """Score threshold of the exact-FROC operating point: the lowest threshold whose detections (score >= threshold)
    contain at most floor(budget * hours) false detections. Tied scores enter or leave together. Returns +inf when
    even the highest-scored group exceeds the budget (no detection emitted)."""
    scores, labels = np.asarray(scores, float), np.asarray(labels)
    cap = int(np.floor(budget * hours + 1e-9))
    if scores.size == 0:
        return np.inf
    order = np.argsort(-scores, kind="stable")
    s, f = scores[order], np.cumsum(labels[order] == FALSE)
    last = np.r_[s[1:] != s[:-1], True]           # last index of each tie group
    ok = last & (f <= cap)
    if not ok.any():
        return np.inf
    return float(s[np.flatnonzero(ok)[-1]])


def hits_at(scores, labels, threshold, n_targets):
    """Boolean vector (n_targets,) of targets hit by detections with score >= threshold."""
    scores, labels = np.asarray(scores, float), np.asarray(labels)
    hit = np.zeros(n_targets, bool)
    sel = (scores >= threshold) & (labels >= 0)
    hit[labels[sel].astype(int)] = True
    return hit


def false_count(scores, labels, threshold):
    scores, labels = np.asarray(scores, float), np.asarray(labels)
    return int(((scores >= threshold) & (labels == FALSE)).sum())


# ---------------------------------------------------------------- S7.2 paired block statistics
def block_sums(hit_a, hit_b, block, n_blocks=None, weights=None):
    """Per block j: d_j = sum w_i (hit_a - hit_b) over its targets and n_j = sum w_i. block: integer ids 0..B-1.
    weights (default 1): the prevalence weights of ADR-0008 item 4, w = real share of the target's cell / share of the
    cell among the evaluated targets."""
    block = np.asarray(block, int)
    nb = (int(block.max()) + 1 if block.size else 0) if n_blocks is None else n_blocks
    w = np.ones(block.size) if weights is None else np.asarray(weights, float)
    d = np.bincount(block, weights=w * (np.asarray(hit_a, float) - np.asarray(hit_b, float)), minlength=nb)
    n = np.bincount(block, weights=w, minlength=nb)
    return d, n


def paired_difference(d, n):
    """Recall difference A - B pooled over blocks: sum d / sum n."""
    return float(np.sum(d) / np.sum(n))


# ---------------------------------------------------------------- S7.3 tests
def exact_mcnemar(b, c):
    """One-sided exact McNemar p-value P(X >= b | X ~ Bin(b + c, 1/2)); b = targets hit only by A, c = only by B."""
    m = b + c
    return sum(comb(m, j) for j in range(b, m + 1)) / 2 ** m if m else 1.0


def clustered_mcnemar_z(d, n=None, margin=0.0):
    """Durkalski et al. (2003) clustered McNemar statistic, one-sided: sum d' / sqrt(sum d'^2), d' = d + margin * n.
    With one target per block it is the McNemar z statistic."""
    d = np.asarray(d, float) + (margin * np.asarray(n, float) if margin else 0.0)
    ss = np.sum(d ** 2)
    return float(np.sum(d) / np.sqrt(ss)) if ss > 0 else 0.0


def normal_sf(z):
    return 0.5 * (1.0 - erf(z / sqrt(2.0)))


def signflip_pvalue(d, n=None, margin=0.0, n_flip=100_000, rng=None, exact_max=16, signs=None, studentize=None):
    """One-sided block sign-flip p-value for H0: recall(A) - recall(B) <= -margin, against A better.

    Margin 0 (superiority): T = sum_j d_j. Under the exchangeability null (the two methods' labels are exchangeable
    within each target) every d_j is symmetric about 0 and independent across blocks, so flipping the sign of each
    block's d_j gives the exact null distribution of T: the randomisation form of the clustered McNemar test, valid for
    any dependence inside a block.

    Margin > 0 (non-inferiority): d'_j = d_j + margin * n_j has mean zero at the boundary but is not symmetric, so the
    raw sum is not exact (its size reached 0.032 at 30 blocks, review of PR #11). The statistic is then studentised by
    the standard deviation of the flipped sample, t(s) = mean(s d') / (sd(s d') / sqrt(B)), which is asymptotically
    exact under asymmetry (Janssen 1999). Asymptotics are slow when discordances are rare (high method correlation,
    few targets per block), so H2 is run at a level calibrated by simulation for each design (ADR-0008 item 7,
    sample_size.calibrate_alpha_h2), never at a nominal 0.025 unchecked. studentize=None picks raw for margin 0 and
    studentised otherwise.

    Exact enumeration when the number of non-zero blocks is <= exact_max, otherwise Monte Carlo with
    p = (1 + #{T* >= T}) / (1 + n_flip). A precomputed (n_flip, B) matrix of i.i.d. +-1 can be passed as signs (the
    power simulation does); then the Monte Carlo form is always used."""
    d = np.asarray(d, float) + (margin * np.asarray(n, float) if margin else 0.0)
    d = d[d != 0]
    if d.size == 0:
        return 1.0
    stud = (margin != 0) if studentize is None else studentize
    if signs is None and d.size <= exact_max:
        s = ((np.arange(2 ** d.size)[:, None] >> np.arange(d.size)) & 1) * 2.0 - 1.0
        include_observed = False
    else:
        if signs is None:
            rng = np.random.default_rng() if rng is None else rng
            s = rng.integers(0, 2, size=(n_flip, d.size), dtype=np.int8) * 2.0 - 1.0
        else:
            s = signs[:, : d.size]
        include_observed = True
    if stud and d.size > 1:
        stat = lambda x: _t_stat(x)
        t, tstar = stat(d[None, :])[0], stat(s * d)
    else:
        t, tstar = d.sum(), s @ d
    tol = 1e-9 * max(1.0, abs(t)) if np.isfinite(t) else 0.0
    k = np.sum(tstar >= t - tol)
    return float((1 + k) / (1 + s.shape[0])) if include_observed else float(k / s.shape[0])


def _t_stat(x):
    """Row-wise one-sample t statistic of the flipped block sums (0 where the spread is 0)."""
    b = x.shape[1]
    m = x.mean(axis=1)
    sd = x.std(axis=1, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(sd > 0, m / (sd / np.sqrt(b)), np.sign(m) * np.inf)


# ---------------------------------------------------------------- S7.4 multiplicity
def iut(pvalues):
    """Intersection-union test: 'A beats every comparator' is rejected-null only if every pairwise null is rejected;
    its p-value is the largest pairwise p-value and needs no multiplicity adjustment (Berger 1982)."""
    return float(max(pvalues))


def fixed_sequence(pvalues, alpha):
    """Hypotheses tested in the registered order at full alpha; testing stops at the first non-rejection.
    Returns a list of booleans (rejected)."""
    out, go = [], True
    for p in pvalues:
        go = go and p <= alpha
        out.append(go)
    return out


def holm(pvalues, alpha):
    """Holm step-down; returns booleans in the input order."""
    p = np.asarray(pvalues, float)
    order = np.argsort(p, kind="stable")
    rej = np.zeros(p.size, bool)
    for k, i in enumerate(order):
        if p[i] > alpha / (p.size - k):
            break
        rej[i] = True
    return rej.tolist()


# ---------------------------------------------------------------- S7.5 interval
def block_bootstrap(stat_fn, strata, n_boot=10_000, level=0.95, rng=None):
    """Stratified block bootstrap: within each stratum (a sequence) resample its blocks with replacement; stat_fn
    receives the concatenated array of resampled block ids (with repetitions) and recomputes the statistic, including
    any operating point. strata: list of arrays of block ids. Returns a dict with the percentile interval, the number
    of blocks and the atom diagnostics (share of replicates equal to each bound): an interval whose bound carries an
    atom of more than 5 % of the replicates, or one built on fewer than 30 blocks, is reported as indicative (SAP
    S7.5)."""
    rng = np.random.default_rng() if rng is None else rng
    strata = [np.asarray(s, int) for s in strata]
    reps = np.empty(n_boot)
    for r in range(n_boot):
        idx = np.concatenate([s[rng.integers(0, s.size, s.size)] for s in strata])
        reps[r] = stat_fn(idx)
    a = (1 - level) / 2
    lo, hi = np.quantile(reps, [a, 1 - a], method="inverted_cdf")
    n_blocks = int(sum(s.size for s in strata))
    atom_lo, atom_hi = float(np.mean(reps == lo)), float(np.mean(reps == hi))
    return {"lo": float(lo), "hi": float(hi), "n_blocks": n_blocks, "atom_lo": atom_lo, "atom_hi": atom_hi,
            "indicative": n_blocks < 30 or max(atom_lo, atom_hi) > 0.05, "replicates": reps}


def recall_diff_stat(d, n):
    """stat_fn for block_bootstrap at fixed operating points: pooled recall difference of the resampled blocks."""
    d, n = np.asarray(d, float), np.asarray(n, float)
    return lambda idx: float(d[idx].sum() / n[idx].sum()) if n[idx].sum() > 0 else 0.0


# ---------------------------------------------------------------- S7.6 false-alarm budget
def budget_check(false_per_block, hours_per_block, budget, tol=0.20, level=0.975):
    """False-alarm budget check (SAP S7.6). Realised rate r = sum f / sum h over blocks. The check passes when a
    quasi-Poisson score test rejects 'rate >= (1 + tol) * budget' at one-sided level 1 - level: the variance is
    phi * mu0 / H at the null mean mu0, with phi the Pearson dispersion of the block counts (floored at 1). The
    score form and the 0.975 level were chosen because the Wald bound r + z se under-covers for right-skewed counts
    (simulated in the tests: size 0.03 to 0.06 with 20 to 100 blocks and negative-binomial counts). 'ub' is the
    score-inverted upper bound (pass iff ub <= mu0); 'calibrated' is descriptive: the symmetric score interval lies
    inside [(1 - tol), (1 + tol)] * budget."""
    f, h = np.asarray(false_per_block, float), np.asarray(hours_per_block, float)
    H, nb = h.sum(), f.size
    r = f.sum() / H
    phi = max(1.0, float(np.sum((f - r * h) ** 2 / (r * h)) / (nb - 1))) if r > 0 and nb > 1 else 1.0
    z = NormalDist().inv_cdf(level)
    a = z * z * phi / H
    ub = r + a / 2 + sqrt(a * a / 4 + a * r)
    lb = r + a / 2 - sqrt(a * a / 4 + a * r)
    mu0 = (1 + tol) * budget
    return {"rate": float(r), "ub": float(ub), "lb": float(lb), "phi": phi, "hours": float(H), "n_blocks": nb,
            "pass": bool(ub <= mu0), "calibrated": bool(lb >= (1 - tol) * budget and ub <= mu0)}
