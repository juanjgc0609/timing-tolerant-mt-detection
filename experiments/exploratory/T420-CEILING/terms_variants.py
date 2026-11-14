"""Station terms for the comparators (T4.24 / T4.25): how the estimate of the terms and the association tolerances
change the result. Added after the first run of ceiling_c.py (14 Nov 2026), which showed that the registered
cross-fitting leaves Hi-net without any term and that the terms lower the recall with PhaseNet picks. Exploratory,
tuning partition only, nothing written to a ledger.

Estimates of the terms (one set per realisation, applied to the picks as t - delta):
  none      no terms (PN reproduces the frozen objective)
  cf        registered: sources that appear in the realisation are left out, at least 3 sources per term, 0 otherwise
  cf-fill   cf, missing terms take the median of the estimated terms of the phase (station_terms.fill_missing)
  cf1-fill  cf-fill with at least 1 source per term
  all       every source of the network, filled (in-sample: optimistic, a target informs its own correction)
  loro      the other realisations, filled (the scheme of the supervision check X7; the same sources recur, optimistic)

Usage:
    python terms_variants.py <frozen.json> <variants,comma> <tune_phasenet.json>... [--tol T --mit M] [--picks PN,OR]
--tol sets tolerance and pick_match_tolerance, --mit min_interevent_time (PyOcto only; a sensitivity, not a tuning).
"""
import argparse
import importlib.util
import json
import os
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("ceiling_c_mod", os.path.join(HERE, "ceiling_c.py"))
CC = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(CC)
SL, ST, TU = CC.SL, CC.ST, CC.TU


def residual_tables(labels):
    """per network: residuals of the injected events, source label and realisation of each."""
    tab = {}
    for net in sorted({lab[0] for lab in labels}):
        R, X, own = [], [], []
        for lab in (x for x in labels if x[0] == net):
            r = labels[lab][0]
            tr = r.truth
            inj = np.array([not str(i).startswith("cat:") for i in tr.id])
            R.append(ST.residuals(tr, r.data["frame"].stations.xyz, r.data["frame"].velocity, np.where(inj)[0])[inj])
            X.append(np.round(np.asarray(tr.xyz, float)[inj], 3))
            own += [lab[1]] * int(inj.sum())
        R = np.concatenate(R)
        K = ST.source_keys(np.concatenate(X))
        ks, sres = ST.per_source(R, K)
        tab[net] = (R, K, ks, sres, np.array(own))
    return tab


def delta(tab, lab, kind):
    R, K, ks, sres, own = tab[lab[0]]
    here = np.unique(K[own == lab[1]])
    if kind == "cf":
        return ST.record_terms(sres, ks, here)
    if kind == "cf-fill":
        return ST.fill_missing(ST.record_terms(sres, ks, here))
    if kind == "cf1-fill":
        return ST.fill_missing(ST.record_terms(sres, ks, here, min_sources=1))
    if kind == "all":
        return ST.fill_missing(ST.terms(sres))
    if kind == "loro":
        m = own != lab[1]
        return ST.fill_missing(ST.terms(ST.per_source(R[m], K[m])[1], min_sources=1))
    raise ValueError(kind)


def out_of_sample(tab):
    """median |residual| of the injected arrivals before and after the registered (cf) terms, and how many terms
    the cross-fitting leaves, per network."""
    out = {}
    for net, (R, K, ks, sres, own) in tab.items():
        b, a, n = [], [], []
        for lab in np.unique(own):
            d = ST.record_terms(sres, ks, np.unique(K[own == lab]))
            x = R[own == lab]
            b.append(x[np.isfinite(x)])
            a.append((x - d[None])[np.isfinite(x)])
            n.append(int((d != 0).sum()))
        b, a = np.concatenate(b), np.concatenate(a)
        out[net] = {"sources": int(len(ks)), "terms_possible": int(R.shape[1] * 2), "terms_estimated_per_realisation": n,
                    "median_abs_residual_before_s": round(float(np.median(np.abs(b))), 3),
                    "median_abs_residual_after_cf_s": round(float(np.median(np.abs(a))), 3)}
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("frozen")
    ap.add_argument("variants")
    ap.add_argument("manifests", nargs="+")
    ap.add_argument("--tol", type=float)
    ap.add_argument("--mit", type=float)
    ap.add_argument("--picks", default="PN,OR")
    a = ap.parse_args(argv)
    warnings.simplefilter("ignore", RuntimeWarning)          # all-NaN medians of unseen stations
    fz = json.load(open(a.frozen))
    method, cfg = fz["method"], dict(fz["params"])
    asc = method.split("+", 1)[1]
    acfg = {k: v for k, v in cfg.items() if k.startswith("as_")}
    if a.tol is not None:
        acfg.update(as_tolerance=a.tol, as_pick_match_tolerance=a.tol)
    if a.mit is not None:
        acfg["as_min_interevent_time"] = a.mit
    recs = []
    for m in a.manifests:
        for r in TU.load_manifest(m, method):
            r.meta["net"] = CC.network_label(m)
            recs.append(r)
    if TU.partition_hash(recs) != fz["partition"]:
        raise SystemExit("not the frozen partition")
    labels = {}
    for r in recs:
        labels.setdefault(CC.group_of(r), []).append(r)
    tab = residual_tables(labels)
    print(method, "tol", a.tol, "mit", a.mit, json.dumps(out_of_sample(tab)), flush=True)
    base = {"PN": lambda r: SL.phasenet_picks(r, cfg), "OR": lambda r: SL.oracle_picks(r.truth)}
    for pk in a.picks.split(","):
        for kind in a.variants.split(","):
            fn = (lambda r, lab: base[pk](r)) if kind == "none" else \
                (lambda r, lab: ST.apply(base[pk](r), delta(tab, lab, kind)))
            dets, _ = CC.run_condition(recs, labels, asc, acfg, fn)
            ev = SL.evaluate(recs, dets)
            bn = CC.by_network(ev)
            row = f"{pk}+{kind:9s} T {ev['recall_T']:.3f} C {ev['recall_C']:.3f} iso {ev['recall_iso']:.3f}"
            for net, v in bn.items():
                row += f" | {net} T {v['recall_T']:.3f} C {v['recall_C']:.3f} iso {v['recall_iso']:.3f}"
            print(row, flush=True)


if __name__ == "__main__":
    main()
