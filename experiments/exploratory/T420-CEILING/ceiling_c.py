"""T4.20 detectability ceiling in regime C, and station terms for Hi-net (T4.24 / T4.25), on the tuning partition.

Part 1 (ceiling): for every C target (regime C or iso, as the objective counts them) of the frozen comparator, the
detection probability of three ideal network subspace detectors (ttmt.detectability: p = 1 coherent oracle, p = 6
moment-tensor GLRT, p = K per-station amplitude) at the objective's budget of 30 false detections per hour, with the
conservative signal energy lambda = sum of squared card SNR over the target's station subset. Output: the ceiling
recall (mean detection probability over C), and the recoverable fraction of the comparator's loss (mean detection
probability over the C targets the comparator misses), overall, per network and in the picker-loss stratum.

Part 2 (station terms): the frozen associator with PhaseNet picks (PN) and oracle picks (OR), with and without static
station terms (ttmt.associate.station_terms, cross-fitted by realisation: the terms of a realisation come only from
sources that do not appear in it). Recall per network, isolated recall included (the diagnosis of PR #38: PyOcto with
oracle picks reaches 0.96 isolated on Ridgecrest and 0.50 on Hi-net).

Inputs are read only (the same files as STAGE-LOSS); nothing touches a ledger or the read partition. Usage:
    python ceiling_c.py <frozen.json> <tune_phasenet.json>... [--out results.json] [--no-terms]
"""
import argparse
import importlib.util
import json
import os
import sys
import time

import numpy as np

from ttmt import detectability as D
from ttmt.associate import station_terms as ST
from ttmt.associate import tune as TU
from ttmt.associate.pipeline import associator

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("stage_loss_mod", os.path.join(HERE, "..", "STAGE-LOSS", "stage_loss.py"))
SL = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(SL)

BUDGET = 30.0
F_TRIAL = 100.0
N_LOC = (1e2, 1e4, 1e6)        # 1e4 is the primary value (declared in the card before the run)
PS = (1, 6, "K")
KAPPA = (1.0, 2.0)             # 1 is the bound; 2 is a sensitivity value only


def network_label(path):
    """network of a manifest, from its path (ss_hinet_2012_egf_v1/... -> hinet)."""
    p = os.path.abspath(path).lower()
    return "hinet" if "hinet" in p else ("ridgecrest" if "ridgecrest" in p else os.path.basename(os.path.dirname(p)))


def network_of(r):
    return r.meta["net"]


def group_of(r):
    return (r.meta["net"], r.name.split("__")[0])


def run_condition(recs, labels, asc, acfg, picks_fn):
    cache, dets, pk = {}, [], {}
    for lab, rs in labels.items():
        r0 = rs[0]
        p = picks_fn(r0, lab)
        cache[lab] = associator(asc).associate(r0.data["frame"], p, acfg)[0], p, SL.match(p, r0.truth)
    for r in recs:
        d, p, m = cache[group_of(r)]
        dets.append(d)
        pk[r.name], pk[r.name + "#picks"] = m, p
    return dets, pk


def by_network(ev):
    out = {}
    nets = np.array([network_of(m[0]) for m in ev["meta"]])
    for net in sorted(set(nets)):
        row = {}
        for nm in ("T", "Cm", "iso"):
            m = ev[nm] & (nets == net)
            row["recall_" + ("C" if nm == "Cm" else nm)] = round(float(ev["got"][m].mean()), 4) if m.any() else None
            row["n_" + ("C" if nm == "Cm" else nm)] = int(m.sum())
        out[net] = row
    return out


def ceiling_part(ev, stage, pk):
    meta = ev["meta"]
    Cm = ev["Cm"]
    lam = np.array([D.lam_lower(m[0].truth.snr[m[1], m[0].subsets[m[1]]]) if m[0].subsets[m[1]] is not None
                    else np.nan for m in meta])
    ks = np.array([int(np.isfinite(m[0].truth.snr[m[1], m[0].subsets[m[1]]]).sum()) if m[0].subsets[m[1]] is not None
                   else 0 for m in meta])
    nets = np.array([network_of(m[0]) for m in meta])
    got = ev["got"]
    groups = {"all": Cm, "loss": Cm & ~got, "loss_picker": Cm & ~got & (stage == "picker")}
    for net in sorted(set(nets)):
        groups[f"{net}"] = Cm & (nets == net)
        groups[f"{net}_loss"] = Cm & ~got & (nets == net)
    res = {"budget_per_hour": BUDGET, "f_trial_hz": F_TRIAL,
           "snr_tot_C_quartiles": [round(float(q), 2) for q in np.nanpercentile(np.sqrt(lam[Cm]), [25, 50, 75])],
           "k_C_counts": {int(k): int((ks[Cm] == k).sum()) for k in np.unique(ks[Cm])},
           "comparator_recall_C": round(float(got[Cm].mean()), 4), "n": {g: int(m.sum()) for g, m in groups.items()},
           "ceilings": {}}
    hits = Cm & got
    for kappa in KAPPA:
        for nl in N_LOC:
            a = D.trial_alpha(BUDGET, F_TRIAL, nl)
            c = D.ceiling(kappa ** 2 * lam, ks, a, PS)
            row = {"alpha": a, "snr_tot_at_pd_half": {f"p={p}": round(D.snr_tot_for_pd(p, a), 3) for p in (1, 6)}}
            for name, v in c.items():
                row[name] = {g: (round(float(np.nanmean(v[m])), 4) if m.any() else None) for g, m in groups.items()}
                # consistency check of the conservative energy: comparator hits the ideal detector would miss
                row[name]["hits_with_pd_below_half"] = round(float(np.mean(v[hits] < 0.5)), 4) if hits.any() else None
            res["ceilings"][f"kappa={kappa:g},n_loc={nl:.0e}"] = row
    return res, lam, ks


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("frozen")
    ap.add_argument("manifests", nargs="+")
    ap.add_argument("--out")
    ap.add_argument("--no-terms", action="store_true")
    a = ap.parse_args(argv)
    fz = json.load(open(a.frozen))
    method, cfg = fz["method"], dict(fz["params"])
    asc = method.split("+", 1)[1]
    recs = []
    for m in a.manifests:
        for r in TU.load_manifest(m, method):
            r.meta["net"] = network_label(m)
            recs.append(r)
    part = TU.partition_hash(recs)
    if part != fz["partition"]:
        sys.exit(f"partition {part} is not the frozen one {fz['partition']}")
    acfg = {k: v for k, v in cfg.items() if k.startswith("as_")}
    min_picks = int(cfg["as_min_picks"])
    labels = {}
    for r in recs:
        labels.setdefault(group_of(r), []).append(r)
    res = {"method": method, "frozen": a.frozen, "partition": part, "frozen_recall_T": fz["recall_T"],
           "frozen_recall_C": fz["recall_C"], "conditions": {}}

    # ---- PN: reproduces the frozen objective, then the ceiling
    t = time.time()
    dets, pk = run_condition(recs, labels, asc, acfg, lambda r, lab: SL.phasenet_picks(r, cfg))
    ev = SL.evaluate(recs, dets)
    stage = SL.attribute(ev, recs, dets, pk, min_picks)
    res["conditions"]["PN"] = dict(SL.summarise(ev, stage), by_network=by_network(ev), wall_s=round(time.time() - t, 1))
    print("PN", json.dumps(res["conditions"]["PN"]), flush=True)
    res["ceiling"], _lam, _ks = ceiling_part(ev, stage, pk)
    print("ceiling", json.dumps(res["ceiling"], indent=1), flush=True)

    # ---- station terms, cross-fitted by realisation
    if not a.no_terms:
        src = {}
        for net in sorted({network_of(r) for r in recs}):
            rs = [labels[lab][0] for lab in labels if lab[0] == net]
            res_all, keys_all = [], []
            for r in rs:
                tr = r.truth
                inj = np.array([not str(i).startswith("cat:") for i in tr.id])
                rr = ST.residuals(tr, r.data["frame"].stations.xyz, r.data["frame"].velocity, np.where(inj)[0])[inj]
                res_all.append(rr)
                keys_all.append(np.round(np.asarray(tr.xyz, float)[inj], 3))
            R = np.concatenate(res_all)
            K = ST.source_keys(np.concatenate(keys_all))
            uk = np.unique(np.concatenate(keys_all), axis=0)
            for variant, ro in (("", False), ("-origin", True)):
                ks, sres = ST.per_source(R, K, remove_origin=ro)
                src[(net, variant)] = (uk, ks, sres)
            ks, sres = ST.per_source(R, K)
            full = ST.terms(sres)
            res.setdefault("terms", {})[net] = {
                "n_sources": int(len(ks)),
                "median_residual_P_S_s": [round(float(np.nanmedian(sres[..., 0])), 3), round(float(np.nanmedian(sres[..., 1])), 3)],
                "iqr_term_P_S_s": [round(float(np.subtract(*np.nanpercentile(full[:, b], [75, 25]))), 3) for b in (0, 1)],
                "rms_residual_after_terms_in_sample_P_S_s": [round(float(np.sqrt(np.nanmean((sres[..., b] - full[None, :, b]) ** 2))), 3) for b in (0, 1)],
                "rms_residual_before_P_S_s": [round(float(np.sqrt(np.nanmean(sres[..., b] ** 2))), 3) for b in (0, 1)]}

        def corrected(base_fn, variant):
            def f(r, lab):
                uk, ks, sres = src[(lab[0], variant)]
                tr = r.truth
                inj = np.array([not str(i).startswith("cat:") for i in tr.id])
                here = np.round(np.asarray(tr.xyz, float)[inj], 3)
                idx = [int(np.where((uk == x).all(axis=1))[0][0]) for x in here]
                delta = ST.record_terms(sres, ks, np.array(idx, int))
                return ST.apply(base_fn(r), delta)
            return f

        conds = {"OR": lambda r, lab: SL.oracle_picks(r.truth),
                 "PN+ST": corrected(lambda r: SL.phasenet_picks(r, cfg), ""),
                 "OR+ST": corrected(lambda r: SL.oracle_picks(r.truth), ""),
                 "OR+ST-origin": corrected(lambda r: SL.oracle_picks(r.truth), "-origin")}
        for name, fn in conds.items():
            t = time.time()
            d2, pk2 = run_condition(recs, labels, asc, acfg, fn)
            ev2 = SL.evaluate(recs, d2)
            st2 = SL.attribute(ev2, recs, d2, pk2, min_picks)
            res["conditions"][name] = dict(SL.summarise(ev2, st2), by_network=by_network(ev2),
                                           wall_s=round(time.time() - t, 1))
            print(name, json.dumps(res["conditions"][name]), flush=True)
            if name == "PN+ST":                         # rule (c) of the card: the loss measured against PN+ST
                res["ceiling_vs_PN+ST"] = ceiling_part(ev2, st2, pk2)[0]
    if a.out:
        json.dump(res, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
