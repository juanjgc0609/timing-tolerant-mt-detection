"""EQNet2 stage 1 on the tuning partition (ATLAS-0 amendment 5): network scores of every tuning record, saved next
to the PhaseNet manifest, and a manifest tune_eqnet2.json that points each record to its score file.

    python eqnet2_stage1.py C:/Users/juanj/ttmt/tune/<set>/tune_phasenet.json ... \
        --checkpoint C:/Users/juanj/ttmt/eqnet2/eqnet_base_noprompt_v0.pth --eqnet-path C:/Users/juanj/ttmt/eqnet2/EQNet \
        --scale ridgecrest_2019=1e6 [--scale hinet_2012=<um/s per count>]

Inputs as the published continuous run of EQNet2 (README of AI4EPS/EQNet at 40a48ed1, mseed command; loader
``eqnet.data.seismic_network`` lines 1800-1856):
- components E, N, Z (the sets' order), 100 Hz;
- amplitude in micrometres per second: the published loader removes the response to m/s and multiplies by 1e6, and
  the hypocentre branch reads log10 of that amplitude (``uper_head.py`` lines 649-677), so the units matter for travel
  time, distance, azimuth and depth (not for the event score, whose branch is moving-normalised). ``--scale`` gives
  the factor from the set's units to um/s per dataset; a dataset without one is refused (Hi-net sets are in counts);
- station location (longitude, latitude, depth in km positive down) from the frame by inverting
  ``ttmt.datasets.Dataset.to_xy`` about the dataset's centre (the centre the sets were built with);
- stations ordered as with ``--sort``: three groups by ``optimal_transport_clustering`` of (lon, lat), then
  lexsort by (group, lon, lat), when there are more than 5 stations. The scores are saved with the station codes in
  that order and ``eqnet2_assoc.load_scores`` restores the frame's order.
Reads only the 'tune' partition (SetReader), as tune_cli.py records does. Uses the GPU when available.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

from ttmt.associate import eqnet2_assoc as EQ


def lonlat_depth(xyz, centre):
    """(n, 3) longitude, latitude (deg), depth (km, positive down) of local frame coordinates (x east, y north, z
    down, km) about centre (lat0, lon0): the inverse of ttmt.datasets.Dataset.to_xy."""
    lat0, lon0 = centre
    xyz = np.asarray(xyz, float)
    lat = lat0 + xyz[:, 1] / 111.19
    lon = lon0 + xyz[:, 0] / (111.19 * np.cos(np.radians(lat0)))
    return np.stack([lon, lat, xyz[:, 2]], 1)


def station_order(lld, cluster=None):
    """published ``--sort`` order (seismic_network.py lines 1838-1845); identity for 5 stations or fewer."""
    n = len(lld)
    if n <= 5:
        return np.arange(n)
    if cluster is None:
        from eqnet.data.seismic_network import optimal_transport_clustering as cluster
    groups = cluster(lld[:, :2], k=3)
    return np.lexsort([lld[:, 1], lld[:, 0], groups])


def dataset_centre(name):
    from ttmt import datasets as DS
    return {"ridgecrest_2019": DS.RIDGECREST, "hinet_2012": DS.HINET}[name].frame()


def run_manifest(path, model, device, scales, centre_of=dataset_centre, reader_of=None, cluster=None, log=print):
    man = json.load(open(path))
    recs = man["records"]
    if not recs:
        raise SystemExit(f"{path}: no records")
    ds = {r["dataset"] for r in recs}
    if len(ds) != 1:
        raise SystemExit(f"{path}: more than one dataset")
    ds = ds.pop()
    if ds not in scales:
        raise SystemExit(f"{path}: no --scale for {ds} (factor from the set's units to um/s)")
    if reader_of is None:
        from ttmt.synth.sets import SetReader

        def reader_of(set_id):
            return SetReader(set_id, "tune", verify="fingerprint")
    reader = reader_of(man["set_id"])
    if reader.partition != "tune":
        raise SystemExit("EQNet2 stage 1 reads the tuning partition only")
    if float(reader.man["fs"]) != EQ.FS:
        raise SystemExit(f"{path}: sampling rate {reader.man['fs']} Hz, EQNet2 needs {EQ.FS}")
    codes = list(reader.man["stations"])
    lld = lonlat_depth(reader.man["xyz"], centre_of(ds))
    order = station_order(lld, cluster)
    root = os.path.dirname(os.path.abspath(path))
    files = {}
    for lab in sorted({r["label"] for r in recs}):
        if lab not in reader.labels:
            raise SystemExit(f"{path}: record {lab} is not in the tuning partition")
        t = time.time()
        wave = np.asarray(reader.data(lab), np.float32) * float(scales[ds])
        t0 = float(reader.truth(lab).t0)
        sc = EQ.annotate(model, device, wave[order], lld[order], t0)
        name = f"eqnet2_{lab}.npz"
        EQ.save_scores(os.path.join(root, name), sc, [codes[i] for i in order])
        files[lab] = name
        log(f"{lab}: {wave.shape[-1] / EQ.FS:.0f} s, {time.time() - t:.1f} s, max score {float(sc['score'].max()):.3f}")
    out = dict(man, records=[dict(r, eqnet2=files[r["label"]]) for r in recs],
               eqnet2={"model": "EQNet2", "weights": EQ.WEIGHTS, "scale_um_per_s": float(scales[ds]),
                       "station_order": [codes[i] for i in order]})
    p = os.path.join(root, "tune_eqnet2.json")
    with open(p, "w") as fh:
        json.dump(out, fh, indent=1)
    return p


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("manifests", nargs="+", help="tune_phasenet.json of each set (written by tune_cli.py records)")
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--eqnet-path", required=True)
    ap.add_argument("--device", default=None)
    ap.add_argument("--scale", action="append", default=[], metavar="DATASET=FACTOR")
    a = ap.parse_args(argv)
    scales = {k: float(v) for k, v in (s.split("=", 1) for s in a.scale)}
    if a.eqnet_path not in sys.path:
        sys.path.insert(0, a.eqnet_path)
    model, device = EQ.load_model(a.checkpoint, a.eqnet_path, a.device)
    for m in a.manifests:
        print(run_manifest(m, model, device, scales))


if __name__ == "__main__":
    sys.exit(main())
