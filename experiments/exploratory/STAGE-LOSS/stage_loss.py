"""Stage-loss diagnostic of the frozen pick-based comparators on the tuning partition (b60d99b3d0bbbbfa).

For a frozen configuration (experiments/exploratory/COMPARATORS/frozen/*.json) it runs counterfactual pick sources
through the same associator and the same objective as the registered tuning (ttmt.associate.tune.objective, 30 false
detections per hour, LOC+T-V):

  PN        PhaseNet picks of the frozen picker parameters (reproduces the frozen objective)
  PN-clean  the same picks with every pick that matches no true arrival removed (false picks off)
  OR        oracle picks: every finite injected arrival (arr_t_obs) of every catalogued event, p = 1; its recall on
            isolated events (regime "iso", reported apart as recall_iso) is the check that the oracle is near 1
  OR-snrX   oracle picks restricted to arrivals with per-station SNR >= X
  OR-mit    oracle picks with the associator's min_interevent_time set to 0.1 s (PyOcto only)
  OR-tight  oracle picks, min_interevent_time, tolerance and pick_match_tolerance all 0.1 s (PyOcto only)
  PN-mit    PhaseNet picks with min_interevent_time 0.1 s (PyOcto only)

and attributes every target missed by PN to one stage, in this order: picker (fewer than min_picks true arrivals
picked), association (no detection within LOC+T, c <= 1), resolution (a detection within LOC+T but not eligible under
the dominance rule or taken by the partner, one to one), budget (an eligible detection below the operating point).
It also reports, for T targets, the information limit of Corollary A.7 (scenario S2, same mechanism, error 0.1:
dtau_min = 0.345 / sqrt(SNR_tot), SNR_tot = sqrt(sum of squared per-station SNR over the event's subset)), the
conservative case because K_mec is unknown (nan) in the EGF sets.

Inputs are read only; nothing here touches the ledgers. Usage:
    python stage_loss.py <frozen.json> <tune_phasenet.json>... [--out results.json] [--only PN,OR]
"""
import argparse
import json
import sys
import time

import numpy as np

from ttmt.associate import pickers, tune as TU
from ttmt.associate.pipeline import associator
from ttmt.atlas import cells as C
from ttmt.atlas import froc as F
from ttmt.atlas.schema import Picks
from ttmt import stats as S

MATCH_S = 0.5           # a pick matches a true arrival of its station and phase within this many seconds
A7_S2 = 0.345           # s * sqrt(SNR): 0.098 s at SNR 12.4, 0.061 at 30, 0.033 at 100 (teoria-fase5 A.7 table)


def phasenet_picks(rec, cfg):
    prob, dt, t0, _ = pickers.load_probabilities(rec.data["prob"](cfg.get("pk_weights")), rec.data["frame"].codes)
    return pickers.picks_from_probabilities(prob, dt, cfg["pk_threshold"], cfg["pk_peak_distance"], t0)


def oracle_picks(truth, snr_min=None):
    a = truth.arr_t_obs
    e, s, ph = np.nonzero(np.isfinite(a))
    if snr_min is not None:
        keep = truth.snr[e, s] >= snr_min
        e, s, ph = e[keep], s[keep], ph[keep]
    t = a[e, s, ph]
    o = np.lexsort((ph, s, t))
    return Picks(s[o], ph[o], t[o], np.ones(len(o)))


def match(picks, truth, tol=MATCH_S):
    """(event of each pick or -1, number of true arrivals picked per event): each pick to its nearest arrival of the
    same station and phase within tol, each arrival at most once (closest pick wins)."""
    a = truth.arr_t_obs
    ev = -np.ones(len(picks), int)
    best = {}
    for i, (s, ph, t) in enumerate(zip(picks.sta, picks.phase, picks.t)):
        d = np.abs(a[:, s, ph] - t)
        if not np.isfinite(d).any():
            continue
        e = int(np.nanargmin(d))
        if d[e] <= tol and (best.get((e, s, ph)) is None or d[e] < best[(e, s, ph)][1]):
            if (e, s, ph) in best:
                ev[best[(e, s, ph)][0]] = -1
            best[(e, s, ph)] = (i, d[e])
            ev[i] = e
    n = np.bincount(ev[ev >= 0], minlength=truth.n_events)
    return ev, n


def subset_picks(p, keep):
    return Picks(p.sta[keep], p.phase[keep], p.t[keep], p.p[keep])


def evaluate(records, dets):
    """objective plus per-target arrays over the pooled events (same order as the objective's pool)."""
    hours = np.concatenate([r.block_hours for r in records])
    pool = F.Pool(hours, np.concatenate([np.full(len(r.block_hours), r.stratum) for r in records]))
    off, meta = 0, []
    for r, d in zip(records, dets):
        ev_block = np.array([r.block(t) for t in r.truth.t], int) + off
        _cells, reg = F.score_record(pool, r.truth, r.k, r.subsets, r.centroid, ev_block,
                                     lambda t, r=r, o=off: r.block(t) + o, {"m": d})
        _c, sep, snr_med, _ = C.assign(r.truth, r.k, r.subsets, r.centroid)
        for e in range(r.truth.n_events):
            st = r.subsets[e] if r.subsets[e] is not None else np.zeros(0, int)
            v = r.truth.snr[e, st]
            meta.append((r, e, sep[e], snr_med[e], float(np.sqrt(np.nansum(v ** 2)))))
        off += len(r.block_hours)
    got = F.hits(pool, "m", TU.BUDGET)
    s, lab, _ = pool.runs["m"]
    thr = F.operating_point(s, lab, np.ones(len(s)), float(hours.sum()), TU.BUDGET)
    T = pool.ev_regime == "T"
    Cm = (pool.ev_regime == "C") | (pool.ev_regime == "iso")
    iso = pool.ev_regime == "iso"
    return dict(recall_T=F.recall(pool, got, T), recall_C=F.recall(pool, got, Cm),
                recall_iso=F.recall(pool, got, iso) if iso.any() else float("nan"), got=got, T=T, Cm=Cm, iso=iso,
                thr=thr, meta=meta, pool=pool)


def attribute(ev, records, dets, picks_of, min_picks):
    """stage of every pooled target (hit, picker, association, resolution, budget)."""
    stage = np.array(["hit"] * len(ev["meta"]), object)
    for j, (r, e, *_rest) in enumerate(ev["meta"]):
        if ev["got"][j]:
            continue
        d = dets[records.index(r)]
        _ev, n = picks_of[r.name]
        if n[e] < min_picks:
            stage[j] = "picker"
            continue
        if len(d) == 0:
            stage[j] = "association"
            continue
        c = S.cost_matrix(d.xyz[:, :2], d.t, r.truth.xyz[:, :2], r.truth.t)
        near = c[:, e] <= 1.0 + 1e-12
        if not near.any():
            stage[j] = "association"
            continue
        lab = F.label(d, r.truth, np.ones(r.truth.n_events, bool))
        own = lab == e
        stage[j] = "budget" if (own & (d.s < ev["thr"])).any() else "resolution"
    return stage


def summarise(ev, stage=None):
    out = {"recall_T": round(ev["recall_T"], 4), "recall_C": round(ev["recall_C"], 4),
           "objective": round((ev["recall_T"] + ev["recall_C"]) / 2, 4), "n_T": int(ev["T"].sum()),
           "n_C": int(ev["Cm"].sum()), "recall_iso": round(ev["recall_iso"], 4), "n_iso": int(ev["iso"].sum()),
           "threshold": float(ev["thr"])}
    if stage is not None:
        for name, m in (("T", ev["T"]), ("C", ev["Cm"]), ("iso", ev["iso"])):
            out["stages_" + name] = {k: int((stage[m] == k).sum())
                                     for k in ("hit", "picker", "association", "resolution", "budget")}
    return out


def a7_ceiling(ev):
    sep = np.array([m[2] for m in ev["meta"]])
    snr_tot = np.array([m[4] for m in ev["meta"]])
    T = ev["T"]
    lim = A7_S2 / np.sqrt(np.maximum(snr_tot, 1e-9))
    return {"T_resolvable_fraction": round(float(np.mean(sep[T] >= lim[T])), 4),
            "dtau_min_median_s": round(float(np.median(lim[T])), 4),
            "snr_tot_median_T": round(float(np.median(snr_tot[T])), 2),
            "snr_median_station_T": round(float(np.nanmedian([m[3] for m, t in zip(ev["meta"], T) if t])), 3),
            "sep_T_quartiles_s": [round(float(x), 3) for x in np.percentile(sep[T], [25, 50, 75])]}


def pick_recall_by_snr(records, picks_of):
    edges = [0, 0.5, 1, 2, 4, 8, np.inf]
    tot, got = np.zeros(len(edges) - 1), np.zeros(len(edges) - 1)
    seen = set()
    for r in records:
        if r.name.split("__")[0] in seen:
            continue
        seen.add(r.name.split("__")[0])
        ev_of_pick, _ = picks_of[r.name]
        tr = r.truth
        hit = np.zeros(tr.arr_t_obs.shape, bool)
        p = picks_of[r.name + "#picks"]
        for i in np.flatnonzero(ev_of_pick >= 0):
            hit[ev_of_pick[i], p.sta[i], p.phase[i]] = True
        snr = np.repeat(tr.snr[:, :, None], 2, axis=2)
        fin = np.isfinite(tr.arr_t_obs) & np.isfinite(snr)
        b = np.digitize(snr[fin], edges) - 1
        np.add.at(tot, b, 1)
        np.add.at(got, b, hit[fin].astype(float))
    return {"snr_bins": [f"{edges[i]}-{edges[i + 1]}" for i in range(len(edges) - 1)],
            "arrivals": tot.astype(int).tolist(), "picked_fraction": np.round(got / np.maximum(tot, 1), 3).tolist()}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("frozen")
    ap.add_argument("manifests", nargs="+")
    ap.add_argument("--out")
    ap.add_argument("--only", default="")
    a = ap.parse_args(argv)
    fz = json.load(open(a.frozen))
    method, cfg = fz["method"], dict(fz["params"])
    asc = method.split("+", 1)[1]
    recs = [r for m in a.manifests for r in TU.load_manifest(m, method)]
    part = TU.partition_hash(recs)
    if part != fz["partition"]:
        sys.exit(f"partition {part} is not the frozen one {fz['partition']}")
    acfg = {k: v for k, v in cfg.items() if k.startswith("as_")}
    min_picks = int(cfg["as_min_picks"])
    labels = {}
    for r in recs:
        labels.setdefault(r.name.split("__")[0], []).append(r)

    def run(picks_fn, acfg_):
        cache, dets, pk = {}, [], {}
        for lab, rs in labels.items():
            r0 = rs[0]
            p = picks_fn(r0)
            cache[lab] = associator(asc).associate(r0.data["frame"], p, acfg_)[0], p, match(p, r0.truth)
        for r in recs:
            d, p, m = cache[r.name.split("__")[0]]
            dets.append(d)
            pk[r.name], pk[r.name + "#picks"] = m, p
        return dets, pk

    conds = {
        "PN": (lambda r: phasenet_picks(r, cfg), acfg),
        "PN-clean": (lambda r: (lambda p: subset_picks(p, match(p, r.truth)[0] >= 0))(phasenet_picks(r, cfg)), acfg),
        "OR": (lambda r: oracle_picks(r.truth), acfg),
        "OR-snr1": (lambda r: oracle_picks(r.truth, 1.0), acfg),
        "OR-snr2": (lambda r: oracle_picks(r.truth, 2.0), acfg),
    }
    if asc == "pyocto":
        conds["OR-mit"] = (lambda r: oracle_picks(r.truth), dict(acfg, as_min_interevent_time=0.1))
        conds["OR-tight"] = (lambda r: oracle_picks(r.truth), dict(acfg, as_min_interevent_time=0.1, as_tolerance=0.1,
                                                                    as_pick_match_tolerance=0.1))
        conds["PN-mit"] = (lambda r: phasenet_picks(r, cfg), dict(acfg, as_min_interevent_time=0.1))
    only = [x for x in a.only.split(",") if x]
    res = {"method": method, "frozen": a.frozen, "partition": part, "frozen_objective": fz["objective"],
           "frozen_recall_T": fz["recall_T"], "frozen_recall_C": fz["recall_C"], "match_s": MATCH_S, "conditions": {}}
    for name, (fn, ac) in conds.items():
        if only and name not in only:
            continue
        t = time.time()
        dets, pk = run(fn, ac)
        ev = evaluate(recs, dets)
        st = attribute(ev, recs, dets, pk, min_picks)
        res["conditions"][name] = dict(summarise(ev, st), wall_s=round(time.time() - t, 1))
        if name == "PN":
            res["pick_recall_by_snr"] = pick_recall_by_snr(recs, pk)
            res["a7"] = a7_ceiling(ev)
        print(name, json.dumps(res["conditions"][name]), flush=True)
    print(json.dumps({k: v for k, v in res.items() if k != "conditions"}, indent=1))
    if a.out:
        json.dump(res, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
