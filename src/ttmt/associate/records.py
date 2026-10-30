"""Tuning records from the generator v2 semi-synthetic sets (thread "Generador sintético v2 y datos").

``write_manifest`` reads one set through a reader bound to the tuning partition (``ttmt.synth.sets.SetReader(set_id,
"tune")``, which refuses every other partition and every set not seeded from the ATLAS-0 anchor), runs the pickers on
each realisation, and writes the manifest that :func:`ttmt.associate.tune.load_manifest` reads: per realisation and
per station count k of the card (``cells.STATION_COUNTS``, those the network has), one record with its truth, the k
stations nearest the network centroid (``cells.station_subset``; the block covers every station, so availability is
the 70 s window around each origin inside the block), the realisation as the one bootstrap block (card, "Analysis":
the block is an injected network realisation) and the partition label the generator wrote, never recomputed.

The distance axis uses the centroid of the network's development catalogue (card, axis 4), read from the set
manifest's "catalogue_centroid"; a manifest without it is refused, never replaced by another centroid. Record names
are the atlas run names (label__k<k>). Both follow ttmt.atlas.plan (run_name, catalogue_centroid), so tuning and the
atlas share the same axes.
"""
import json
import os

import numpy as np

from ..atlas import cells as C
from . import pickers

VELOCITY = {"ridgecrest_2019": "hadley_kanamori", "hinet_2012": "jma2001"}
TRUTH_FIELDS = ("id", "t", "xyz", "mag", "arr_t", "snr", "target", "pair", "kmec", "arr_t_obs")


def run_name(label, k, variant=""):
    """the atlas run name (ttmt.atlas.plan.run_name)."""
    return f"{label}__k{k}" + (f"__{variant}" if variant else "")


def catalogue_centroid(manifest):
    """(x, y) local km of the development catalogue (ttmt.atlas.plan.catalogue_centroid); KeyError if absent."""
    c = manifest.get("catalogue_centroid")
    if c is None:
        raise KeyError(f"manifest of {manifest.get('set_id')!r} has no catalogue_centroid (card item 4: distance to "
                       "the centroid of the network's development catalogue)")
    return tuple(map(float, c[:2]))


def save_truth(path, truth):
    """truth.npz with the fields load_manifest reads (None fields left out)."""
    f = {k: getattr(truth, k) for k in TRUTH_FIELDS if getattr(truth, k) is not None}
    np.savez_compressed(path, t0=float(truth.t0), t1=float(truth.t1), **f)


def stream_of(data, codes, fs, t0, net="XX", band="HH"):
    """obspy Stream of a (ns, 3, n) array, components E, N, Z, first sample at record time t0 (s after base.EPOCH)."""
    from obspy import Stream, Trace, UTCDateTime

    from .base import to_datetime
    start = UTCDateTime(to_datetime(t0))
    st = Stream()
    for i, c in enumerate(codes):
        for j, comp in enumerate("ENZ"):
            st.append(Trace(np.asarray(data[i, j], np.float32), header=dict(
                network=net, station=str(c), channel=band + comp, sampling_rate=fs, starttime=start)))
    return st


def subsets_for(truth, sta_xyz, codes, k):
    """per event the k stations nearest the network centroid among those whose data cover [t - 10 s, t + 60 s]."""
    cov = [[(truth.t0, truth.t1)]] * len(codes)
    avail = C.availability(truth.t, cov)
    return [C.station_subset(sta_xyz, codes, avail[e], k) for e in range(len(truth.t))]


def write_manifest(reader, out, annotate, pickers_=(("PhaseNet", "original"), ("PhaseNet", "instance")), log=print):
    """write out/<set id>/ (truth and probability files) and out/<set id>/tune_<model>.json, one manifest per picker
    model with its weights as "prob" entries. reader: SetReader of partition "tune" (or anything with .man, .set_id,
    .partition, iteration over labels, .truth(label) and .data(label)); annotate(stream, codes, model, weights, t0)
    -> (prob, dt, t0), e.g. pickers.annotate_seisbench. Returns {model: manifest path}."""
    if reader.partition != "tune":
        raise ValueError(f"tuning records come from the 'tune' partition, not {reader.partition!r}")
    man = reader.man
    root = os.path.join(out, reader.set_id)
    os.makedirs(root, exist_ok=True)
    codes, sta_xyz, fs = list(man["stations"]), np.asarray(man["xyz"], float), float(man["fs"])
    ks = [k for k in C.STATION_COUNTS if k <= len(codes)]
    centroid = catalogue_centroid(man)
    labels = list(reader)
    truths = {lab: reader.truth(lab) for lab in labels}
    models = sorted({m for m, _w in pickers_})
    recs = {m: [] for m in models}
    real = {r["label"]: r for r in man["realisations"]}
    for i, lab in enumerate(labels):
        truth, r = truths[lab], real[lab]
        save_truth(os.path.join(root, f"truth_{lab}.npz"), truth)
        st = stream_of(reader.data(lab), codes, fs, truth.t0)
        prob = {}
        for m, w in pickers_:
            p, dt, t0 = annotate(st, codes, m, w, truth.t0)
            name = f"prob_{lab}_{m}_{w}.npz"
            pickers.save_probabilities(os.path.join(root, name), p, dt, t0, codes, m, w)
            prob.setdefault(m, {})[w] = name
        for k in ks:
            sub = subsets_for(truth, sta_xyz, codes, k)
            for m in models:
                recs[m].append({
                    "name": run_name(lab, k), "partition": r["partition"], "truth": f"truth_{lab}.npz", "k": k,
                    "stations": {"code": codes, "xyz": sta_xyz.tolist()},
                    "subsets": [None if s is None else s.tolist() for s in sub],
                    "centroid": [float(centroid[0]), float(centroid[1])],
                    "block_edges": [float(truth.t0), float(truth.t1)], "stratum": i,
                    "prob": prob[m], "velocity": VELOCITY.get(man["dataset"], "hadley_kanamori"),
                    "dataset": man["dataset"], "set_id": reader.set_id, "label": lab,
                    "block_index": r.get("block_index"), "r": r.get("r")})
        log(f"{lab}: {int(truth.target.sum()) if truth.target is not None else 0} targets, k {ks}")
    paths = {}
    for m in models:
        p = os.path.join(root, f"tune_{m.lower()}.json")
        with open(p, "w") as fh:
            json.dump({"set_id": reader.set_id, "inject_config_hash": man.get("inject_config_hash"),
                       "records": recs[m]}, fh, indent=1)
        paths[m] = p
    return paths


def seisbench_annotate(device=None):
    """annotate callable for write_manifest backed by pickers.annotate_seisbench."""
    def f(stream, codes, model, weights, t0):
        return pickers.annotate_seisbench(stream, codes, model, weights, device=device, t0=t0)
    return f
