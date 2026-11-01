"""The ATLAS-0 analysis on the reading partition: pools the runs of the plan (ttmt.atlas.plan) and reports, per method,
the recall of the regimes and strata with stratified block-bootstrap intervals (froc.describe), the pooled loss of
H-ATLAS (froc.pooled_loss) and the duplicator control (ADR-0008 item 3). The atlas makes no beat claim.

Pooling (card, "Regimes and strata" and "The single rule"): one bootstrap block is one realisation (one noise block of
one set with its injected events), holding all its station-subset runs; its hours are the sum of those runs' hours,
since every run emits its own false detections. The bootstrap is stratified by set (sequence or network). The
operating point of a method is set on the whole pooled evaluation at the budget.

Sensitivity to the bootstrap unit (review of PR #21): realisations of one set on the same 300 s noise block share
every noise sample, so the false detections driven by the noise are not independent across them; so do the EGF
and known-tensor sets of one network, injected into the same noise blocks (the noise block is keyed by dataset).
Next to the registered intervals, every scope reports the number of distinct noise blocks and a sensitivity row whose
bootstrap draws noise blocks (all realisations of a noise block together, strata sharing a noise block drawn
together), with H-ATLAS read against both. The registered unit stays the realisation; the sensitivity row is
descriptive and is stated in the run log.

Scopes reported: every set pooled; each set alone; the stratum of injected events whose largest correlation to any
template or deflation atom is below 0.5 (CR-06; truth.meta["max_cc"] per event); each contamination variant of the
4-station subsets (CR-24), apart from the main pool. A run is read only from the "read" partition, and only once the
reproduction gate of the card has passed for at least two comparators (gate.json)."""
import json
from collections import OrderedDict
from dataclasses import replace
from pathlib import Path

import numpy as np

from . import froc as F

ATLAS_DIR = Path(__file__).resolve().parents[3] / "experiments" / "exploratory" / "ATLAS-0"
GATE_FILE = ATLAS_DIR / "gate.json"
CC_LIMIT = 0.5


class GateError(RuntimeError):
    pass


def require_gate(path=GATE_FILE):
    """the reproduction gate (card, "Gate"): at least two comparators reproduce their published numbers. gate.json:
    {"targets": {name: {"passed": bool, "evidence": str}}}."""
    p = Path(path)
    if not p.exists():
        raise GateError(f"{p} is missing: the reproduction gate has not been reported, so no atlas data are read")
    g = json.loads(p.read_text())
    passed = sorted(k for k, v in g.get("targets", {}).items() if v.get("passed"))
    if len(passed) < 2:
        raise GateError(f"the reproduction gate passed for {passed or 'no target'}; the card needs at least two")
    return passed


def build_pool(rows, truth_of, dets_of, methods, centroid_of, target_filter=None, duplicator_of=None):
    """Pool of the runs `rows` (plan rows, all of the "read" partition). truth_of(row) -> Truth of the run (with the
    variant's SNR for a contamination run); dets_of(method, row) -> Detections; centroid_of(set_id) -> (x, y) of the
    development catalogue; target_filter(truth) -> bool (m,) ANDed with the target flags (events outside it become
    neutral); duplicator_of: comparator whose duplicator joins the pool as "dup:<name>"."""
    bad = [r["run"] for r in rows if r["partition"] != "read"]
    if bad:
        raise ValueError(f"the atlas reads the 'read' partition only; {len(bad)} runs are not, e.g. {bad[0]}")
    blocks, sets, noise, hours, stratum, cluster = OrderedDict(), OrderedDict(), OrderedDict(), [], [], []
    truths = [truth_of(r) for r in rows]
    for r, tr in zip(rows, truths):
        key = (r["set_id"], r["label"])
        if key not in blocks:
            blocks[key] = len(blocks)
            sets.setdefault(r["set_id"], len(sets))
            hours.append(0.0)
            stratum.append(sets[r["set_id"]])
            nkey = (r.get("dataset") or r["set_id"], int(r["block_index"]))   # sets on one noise share it
            cluster.append(noise.setdefault(nkey, len(noise)))
        hours[blocks[key]] += tr.hours
    pool = F.Pool(np.array(hours), np.array(stratum, int), block_cluster=np.array(cluster, int))
    for r, tr in zip(rows, truths):
        if target_filter is not None:
            tr = replace(tr, target=tr.target & np.asarray(target_filter(tr), bool))
        b = blocks[(r["set_id"], r["label"])]
        idx = np.array([int(i) for i in str(r["station_idx"]).split()], int)
        runs = {m: dets_of(m, r) for m in methods}
        F.score_record(pool, tr, int(r["k"]), [idx] * tr.n_events, centroid_of(r["set_id"]),
                       np.full(tr.n_events, b), lambda t, b=b: b, runs, duplicator_of=duplicator_of)
    return pool


def _regimes(d):
    return {m: {g: dict(zip(("recall", "lo", "hi", "n"), v)) for g, v in x["regime"].items()} for m, x in d.items()}


LOSS_FLOOR = 0.2                                            # card, H-ATLAS: every core pipeline loses >= 20 %


def h_atlas_reading(losses):
    """H-ATLAS on the losses (pooled_loss) of the core pipelines: refuted if the upper bound of any is below 0.20,
    confirmed if every lower bound is at least 0.20, otherwise undecided."""
    refuted_by = sorted(m for m, x in losses.items() if np.isfinite(x["hi"]) and x["hi"] < LOSS_FLOOR)
    if refuted_by:
        verdict = "refuted"
    elif losses and all(np.isfinite(x["lo"]) and x["lo"] >= LOSS_FLOOR for x in losses.values()):
        verdict = "confirmed"
    else:
        verdict = "undecided"
    return dict(verdict=verdict, refuted_by=refuted_by)


def scope_report(pool, methods, comparators, n_boot=2000):
    """the regime table, H-ATLAS loss and duplicator control of one pool, with the noise-block sensitivity row."""
    d = F.describe(pool, methods, n_boot=n_boot)
    core = [m for m in methods if not m.startswith("dup:")]
    loss = {m: F.pooled_loss(pool, m, n_boot=n_boot) for m in core}
    pipes = list(comparators) or core                       # H-ATLAS is read on the core pipelines (comparators)
    out = dict(blocks=int(pool.n_blocks), hours=float(pool.block_hours.sum()),
               targets={g: int((pool.ev_regime == g).sum()) for g in ("T", "C", "iso")},
               indicative=bool(pool.n_blocks < F.MIN_BLOCKS), regimes=_regimes(d),
               h_atlas=loss, h_atlas_reading=h_atlas_reading({m: loss[m] for m in pipes}))
    if pool.block_cluster is not None:
        nc = F.n_clusters(pool, pool.block_cluster)
        ds = F.describe(pool, methods, n_boot=n_boot, cluster=pool.block_cluster)
        ls = {m: F.pooled_loss(pool, m, n_boot=n_boot, cluster=pool.block_cluster) for m in core}
        out["noise_blocks"] = nc
        out["sensitivity_noise_block"] = dict(unit="noise block", blocks=nc, indicative=bool(nc < F.MIN_BLOCKS),
                                              regimes=_regimes(ds), h_atlas=ls,
                                              h_atlas_reading=h_atlas_reading({m: ls[m] for m in pipes}))
    ctl = {}
    for c in comparators:
        if "dup:" + c in pool.runs:
            r = F._pairwise(pool, "dup:" + c, c, pool.ev_regime == "T", F.BUDGETS[0], 0.0, None)
            ctl[c] = dict(r, broken=bool(r["p"] <= F.ALPHA))
    out["control"] = ctl
    return out, d


def max_cc_below(limit=CC_LIMIT):
    def f(truth):
        cc = np.asarray(truth.meta.get("max_cc", np.full(truth.n_events, np.nan)), float)
        return cc < limit                                   # nan (not injected / unknown) is outside the stratum
    return f


def report(rows, truth_of, dets_of, methods, centroid_of, comparators=(), n_boot=2000, gate=GATE_FILE):
    """every scope of the atlas on the reading rows of the plan. Returns (report dict, cell tables per scope). The
    duplicator control runs on the first comparator listed (the card's primary comparator)."""
    passed = require_gate(gate) if gate is not None else []
    main = [r for r in rows if not r["variant"]]
    dup = comparators[0] if comparators else None
    names = list(methods) + (["dup:" + dup] if dup else [])
    scopes = OrderedDict([("all", (main, None))])
    for s in sorted({r["set_id"] for r in main}):
        scopes[s] = ([r for r in main if r["set_id"] == s], None)
    scopes["cc<0.5"] = (main, max_cc_below())
    for v in sorted({r["variant"] for r in rows if r["variant"]}):
        scopes["contamination:" + v] = ([r for r in rows if r["variant"] == v], None)
    rep, cells = dict(gate=passed, methods=list(methods), comparators=list(comparators), scopes={}), {}
    for name, (rs, filt) in scopes.items():
        if not rs:
            continue
        pool = build_pool(rs, truth_of, dets_of, methods, centroid_of, filt, dup)
        rep["scopes"][name], cells[name] = scope_report(pool, names, comparators, n_boot)
    return rep, cells


def to_json(obj):
    def conv(x):
        if isinstance(x, dict):
            return {str(k): conv(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [conv(v) for v in x]
        if isinstance(x, np.ndarray):
            return conv(x.tolist())
        if isinstance(x, (np.floating, float)):
            return None if not np.isfinite(x) else float(x)
        if isinstance(x, (np.integer, np.bool_)):
            return x.item()
        return x
    return json.dumps(conv(obj), indent=1)
