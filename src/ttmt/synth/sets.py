"""Semi-synthetic sets on disk: build, manifest with checksums, and a reader bound to one partition.

Layout of a set (under ``$TTMT_WORK_ROOT/synth/<set id>/``, never in git):
    noise_<block start>.npz      the real noise block, float32 (ns, 3, n), components (E, N, Z up)
    sig_<label>.npz              the injected signal of one realisation, float32, same shape
    truth_<label>.json           its ttmt.atlas.schema.Truth (nan as null)
and the manifest ``data/manifests/synth/<set id>.json`` (in git): set id, dataset, source kind, generator config and
hashes, the anchor string the seeds derive from, and per realisation its label, seed, noise block, partition, the
4-station contamination variants (CR-24), and per file its SHA-256 and a float fingerprint (sum of squares per
station, so a rebuild on another platform can be checked to a tolerance when the bytes differ).

Partitions (ATLAS-0, Tuning): 'tune' and 'read' split the noise blocks 30/70 by a hash of the block start
(``realdata.partition_blocks``, at least one block each); every realisation on a block inherits its partition, so no noise sample is shared. A
:class:`SetReader` serves one partition and raises :class:`ttmt.partition.PartitionError` on any other.

Seeds (ADR-0004 item 4, ADR-0009 item 5): with ``anchor='ATLAS-0'`` (the default) realisation r on noise block i
(the block starting at rec.t0 + 300 i) of set s draws its injections, timing errors and amplitudes from
``ttmt.atlas.seeds.realisation_rng(s, i, r)`` and its contamination stations from ``realisation_rng(s, i, r,
'subset')``; s, i, r and the master seed go into Truth.meta. Any other anchor string hashes into a development seed
with ``derive_seed``, is flagged ``pre_anchor`` in the manifest and must not be used for tuning or reading.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

from .. import config as CFG
from ..atlas.schema import Truth
from ..partition import PartitionError
from . import inject as INJ
from . import realdata as RD

MANIFEST_DIR = CFG.REPO_ROOT / "data" / "manifests" / "synth"
GENERATOR_VERSION = "ttmt.synth v2.0"


def set_root(set_id):
    return Path(CFG.work_root()) / "synth" / set_id


def derive_seed(anchor, set_id, block, r):
    h = hashlib.sha256(f"{anchor}|{set_id}|{block:.3f}|{r}".encode()).digest()
    return int.from_bytes(h[:8], "little") & ((1 << 63) - 1)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fingerprint(a):
    a = np.asarray(a, np.float64)
    return np.sum(a * a, axis=tuple(range(1, a.ndim))).tolist()


# ------------------------------------------------------------------------------------------------ truth json
def _clean(x):
    if isinstance(x, dict):
        return {k: _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return _clean(x.tolist())
    if isinstance(x, (np.floating, float)):
        return None if not np.isfinite(x) else float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def truth_to_dict(t):
    return _clean({"id": t.id.tolist(), "t": t.t, "xyz": t.xyz, "mag": t.mag, "arr_t": t.arr_t, "snr": t.snr,
                   "t0": t.t0, "t1": t.t1, "target": t.target, "pair": t.pair, "kmec": t.kmec,
                   "arr_t_obs": t.arr_t_obs, "meta": t.meta})


def _f(x):
    return np.array(x, dtype=float) if x is not None else None


def _nan(x):
    """nested lists with None -> float array with nan."""
    return np.array(json.loads(json.dumps(x).replace("null", "NaN")), float)


def truth_from_dict(d):
    m = len(d["t"])
    ns = len(d["meta"]["stations"])
    shp = (m, ns, 2)
    return Truth(np.array(d["id"]), _nan(d["t"]), _nan(d["xyz"]).reshape(m, 3), _nan(d["mag"]),
                 _nan(d["arr_t"]).reshape(shp), _nan(d["snr"]).reshape(m, ns), d["t0"], d["t1"],
                 np.array(d["target"], bool), np.array(d["pair"], int), _nan(d["kmec"]),
                 _nan(d["arr_t_obs"]).reshape(shp), d["meta"])


# ------------------------------------------------------------------------------------------------ build
def subset_k(xyz, codes, k):
    """the k stations nearest the network centroid, ties broken by code (``ttmt.atlas.cells.station_subset`` with every
    station available: the injected records have no gaps)."""
    from ..atlas.cells import station_subset
    return station_subset(xyz, codes, np.ones(len(codes), bool), k).tolist()


def build(set_id, rec, sources, per_block, anchor="ATLAS-0", inj_cfg=None, atoms=None, atom_ids=None,
          source_meta=None, blocks=None, log=print):
    """build a set: per_block realisations on every noise-role block of rec (or the given blocks); writes files and the
    manifest; returns the manifest dict."""
    root = set_root(set_id)
    root.mkdir(parents=True, exist_ok=True)
    c = INJ.config(inj_cfg)
    split = {float(bl[0]): r for bl, r in RD.tune_read_split(rec)}
    if blocks is None:
        blocks = [bl for bl, r in RD.tune_read_split(rec) if r != "egf"]
    if any(split.get(float(a)) not in ("tune", "read") for a, _b in blocks):
        raise ValueError("injections go only into noise-role blocks")
    sub4 = subset_k(rec.xyz, rec.codes, 4)
    parts = {float(a): split[float(a)] for a, _b in blocks}
    reals, files = [], {}
    for a, b in blocks:
        part = parts[float(a)]
        nf = root / f"noise_{a:07.1f}.npz"
        noise = rec.cut(a, b).astype(np.float32)
        np.savez_compressed(nf, noise=noise)
        files[nf.name] = {"sha256": sha256(nf), "fingerprint": fingerprint(noise)}
        for r in range(per_block):
            label = f"{set_id}_{a:07.1f}_{r:03d}"
            i = int(round((a - rec.t0) / RD.BLOCK_S))
            if anchor == "ATLAS-0":
                from ..atlas import seeds as SEEDS
                seed = SEEDS.realisation_rng(set_id, i, r)
                rng_sub = SEEDS.realisation_rng(set_id, i, r, "subset")
                seed_rec = {"set_id": set_id, "block_index": i, "r": r,
                            "master_seed_hex": f"{SEEDS.master_seed():016x}", "streams": ["inject", "subset"]}
            else:
                seed = derive_seed(anchor, set_id, a, r)
                rng_sub = np.random.default_rng(seed ^ 0x5EED)
                seed_rec = {"seed": seed}
            sig, _noise, truth = INJ.realise(rec, (a, b), sources, seed, c, atoms, atom_ids, label)
            variants = INJ.contamination_variants(sub4, rng_sub)
            truth.meta.update(seed_rec)
            truth.meta.update({"partition": part, "subset4": sub4, "contamination_variants": variants,
                               "generator": GENERATOR_VERSION, "anchor": anchor, "set_id": set_id,
                               "source_meta": source_meta or {}})
            sf = root / f"sig_{label}.npz"
            np.savez_compressed(sf, sig=sig.astype(np.float32))
            tf = root / f"truth_{label}.json"
            tf.write_text(json.dumps(truth_to_dict(truth)))
            files[sf.name] = {"sha256": sha256(sf), "fingerprint": fingerprint(sig.astype(np.float32))}
            files[tf.name] = {"sha256": sha256(tf)}
            reals.append({"label": label, **seed_rec, "block": [a, b], "partition": part,
                          "noise": nf.name, "sig": sf.name, "truth": tf.name, "variants": variants,
                          "n_targets": int(truth.target.sum()), "n_pairs": int(len(set(truth.pair[truth.target]) - {-1}))})
            log(f"{label} {part} targets {reals[-1]['n_targets']} pairs {reals[-1]['n_pairs']}")
    man = {"set_id": set_id, "generator": GENERATOR_VERSION, "dataset": rec.dataset, "fs": rec.fs,
           "components": "E, N, Z up", "stations": rec.codes, "xyz": np.asarray(rec.xyz).tolist(),
           "anchor": anchor, "pre_anchor": anchor != "ATLAS-0", "inject_config": c,
           "inject_config_hash": INJ.config_hash(c), "source_meta": source_meta or {},
           "partitions": {"rule": "realdata.partition_blocks: sha256(dataset|block start|tune-read-v1) < 0.3 -> tune, "
                                  "at least one block per partition",
                          "tune": [x["label"] for x in reals if x["partition"] == "tune"],
                          "read": [x["label"] for x in reals if x["partition"] == "read"]},
           "realisations": reals, "files": files}
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    (MANIFEST_DIR / f"{set_id}.json").write_text(json.dumps(_clean(man), indent=1))
    return man


# ------------------------------------------------------------------------------------------------ read
class SetReader:
    """realisations of one partition of a set; a set not built from the ATLAS-0 anchor (manifest ``pre_anchor``) is
    refused unless ``allow_pre_anchor`` (development only). verify: 'sha256' (exact bytes), 'fingerprint' (rel. tolerance tol, for
    a rebuild on another platform) or None."""

    def __init__(self, set_id, partition, verify="sha256", tol=1e-6, manifest=None, allow_pre_anchor=False):
        self.man = json.loads(Path(manifest or MANIFEST_DIR / f"{set_id}.json").read_text())
        if self.man.get("pre_anchor", True) and not allow_pre_anchor:
            raise PartitionError(f"set {set_id!r} is not seeded from the ATLAS-0 anchor (pre_anchor); it serves no "
                                 "tuning or reading (pass allow_pre_anchor=True for development only)")
        self.set_id, self.partition, self.verify, self.tol = set_id, partition, verify, tol
        if partition not in ("tune", "read"):
            raise PartitionError(f"unknown partition {partition!r}")
        self.root = set_root(set_id)
        self.labels = [r["label"] for r in self.man["realisations"] if r["partition"] == partition]

    def _check(self, name, arr=None):
        ent = self.man["files"][name]
        if self.verify == "sha256":
            if sha256(self.root / name) != ent["sha256"]:
                raise OSError(f"checksum mismatch: {name}")
        elif self.verify == "fingerprint" and arr is not None:
            fp = np.array(fingerprint(arr))
            ref = np.array(ent["fingerprint"])
            if not np.allclose(fp, ref, rtol=self.tol, atol=0):
                raise OSError(f"fingerprint mismatch: {name}")

    def _real(self, label):
        r = next((x for x in self.man["realisations"] if x["label"] == label), None)
        if r is None:
            raise KeyError(label)
        if r["partition"] != self.partition:
            raise PartitionError(f"{label} belongs to partition {r['partition']!r}, reader serves {self.partition!r}")
        return r

    def truth(self, label, variant=None):
        r = self._real(label)
        self._check(r["truth"])
        t = truth_from_dict(json.loads((self.root / r["truth"]).read_text()))
        if variant is not None:
            name, st, f = self._variant(r, variant)
            t = INJ.truth_contaminated(t, st, f)
        return t

    def _variant(self, r, variant):
        v = next((v for v in r["variants"] if v[0] == variant), None)
        if v is None:
            raise KeyError(variant)
        return v

    def arrays(self, label, variant=None):
        """(signal, noise) float64 (ns, 3, n); with a variant the noise of its stations is scaled (CR-24)."""
        r = self._real(label)
        noise = np.load(self.root / r["noise"])["noise"]
        sig = np.load(self.root / r["sig"])["sig"]
        self._check(r["noise"], noise)
        self._check(r["sig"], sig)
        noise, sig = noise.astype(np.float64), sig.astype(np.float64)
        if variant is not None:
            _n, st, f = self._variant(r, variant)
            noise = INJ.contaminate(noise, st, f)
        return sig, noise

    def data(self, label, variant=None):
        s, n = self.arrays(label, variant)
        return s + n

    def __iter__(self):
        return iter(self.labels)
