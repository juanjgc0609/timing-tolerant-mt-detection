"""Station terms from an independent calibration set (catalogue picks) and their check on the tuning partition.

    python catalogue_terms.py <scedc_dir> <network> <frozen.json> <tune_phasenet.json>... [--out terms.json] [--assoc]

<scedc_dir> holds events.csv and picks.csv of tools/fetch_scedc_phases.py; <network> is the label of the network the
picks belong to (ridgecrest). The terms are estimated with ttmt.associate.station_terms.from_catalogue against the
1-D model of the frame, filled, and compared with the in-sample terms of the injected sources (T420-CEILING). With
--assoc the frozen PyOcto is run with PhaseNet and oracle picks, with and without the terms (the other networks keep
their picks unchanged). Exploratory; reads the tuning partition only; nothing is written to a ledger.
"""
import argparse
import csv
import importlib.util
import json
import os
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("terms_variants_mod",
                                               os.path.join(HERE, "..", "T420-CEILING", "terms_variants.py"))
TV = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(TV)
CC, SL, ST, TU = TV.CC, TV.SL, TV.ST, TV.TU


def stats(x):
    x = x[np.isfinite(x)]
    return {"median_abs_s": round(float(np.median(np.abs(x))), 3), "rms_s": round(float(np.sqrt(np.mean(x ** 2))), 3)}


def free_origin(x):
    """residuals with each event's median removed (what an associator that fits the origin time sees)."""
    return np.concatenate([(e - np.nanmedian(e))[np.isfinite(e)] for e in x if np.isfinite(e).sum() >= 2])


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("scedc_dir")
    ap.add_argument("network")
    ap.add_argument("frozen")
    ap.add_argument("manifests", nargs="+")
    ap.add_argument("--out")
    ap.add_argument("--assoc", action="store_true")
    a = ap.parse_args(argv)
    warnings.simplefilter("ignore", RuntimeWarning)
    fz = json.load(open(a.frozen))
    method, cfg = fz["method"], dict(fz["params"])
    acfg = {k: v for k, v in cfg.items() if k.startswith("as_")}
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
    frame = next(r for r in recs if r.meta["net"] == a.network).data["frame"]
    codes = [str(c) for c in frame.stations.code]
    events = list(csv.DictReader(open(os.path.join(a.scedc_dir, "events.csv"))))
    depth = {e["evid"]: float(e["depth_km"]) for e in events}
    with open(os.path.join(a.scedc_dir, "picks.csv")) as f:
        cat, n = ST.from_catalogue(codes, frame.velocity, csv.DictReader(f), depth)
    catf = ST.fill_missing(cat)
    R, _K, _ks, sres, _own = TV.residual_tables(labels)[a.network]
    ins = ST.terms(sres)
    both = (cat != 0) & (ins != 0)
    out = {"network": a.network, "events": len(events), "origin_range": [min(e["origin"] for e in events),
                                                                         max(e["origin"] for e in events)],
           "codes": codes, "terms_P_S_s": np.round(catf, 4).tolist(), "n_picks_P_S": n.tolist(),
           "in_sample_injected_P_S_s": np.round(ins, 4).tolist(),
           "terms_in_both": int(both.sum()),
           "correlation_with_injected": round(float(np.corrcoef(cat[both], ins[both])[0, 1]), 3),
           "median_difference_s": round(float(np.median(cat[both] - ins[both])), 3),
           "injected_residual": {"none": stats(R), "catalogue": stats(R - catf[None]),
                                 "in_sample": stats(R - ST.fill_missing(ins)[None])},
           "injected_residual_free_origin": {"none": stats(free_origin(R)), "catalogue": stats(free_origin(R - catf[None])),
                                             "in_sample": stats(free_origin(R - ST.fill_missing(ins)[None]))}}
    print(json.dumps(out, indent=1), flush=True)
    if a.assoc:
        base = {"PN": lambda r: SL.phasenet_picks(r, cfg), "OR": lambda r: SL.oracle_picks(r.truth)}
        out["association"] = {}
        for pk in ("PN", "OR"):
            for kind in ("none", "catalogue"):
                fn = (lambda r, lab: base[pk](r)) if kind == "none" else \
                    (lambda r, lab: ST.apply(base[pk](r), catf) if lab[0] == a.network else base[pk](r))
                dets, _ = CC.run_condition(recs, labels, method.split("+", 1)[1], acfg, fn)
                ev = SL.evaluate(recs, dets)
                row = {"T": round(ev["recall_T"], 4), "C": round(ev["recall_C"], 4), a.network: CC.by_network(ev)[a.network]}
                out["association"][f"{pk}+{kind}"] = row
                print(f"{pk}+{kind}", json.dumps(row), flush=True)
    if a.out:
        json.dump(out, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
