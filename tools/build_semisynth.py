"""Build a semi-synthetic set for ATLAS-0 (ttmt.synth.sets) from a registered real record.

    TTMT_DATA_ROOT=... TTMT_WORK_ROOT=... python tools/build_semisynth.py --dataset ridgecrest_2019 --source egf \\
        --per-block 10 [--fk-lib work/synth/fk_ridgecrest_2019_11.pkl]

source egf: recorded events of the EGF-role blocks (tensors unknown); source fk: full-wave synthetic events with known
tensors (tools/build_fk_library.py). Writes the files under $TTMT_WORK_ROOT/synth/<set id>/ and the manifest
data/manifests/synth/<set id>.json. Seeds: ttmt.atlas.seeds.realisation_rng(set id, block, r) (ADR-0009 item 5).
"""
import argparse
import hashlib
import pickle

from ttmt.runtime import pin_blas

pin_blas()

from ttmt.synth import fklib  # noqa: E402
from ttmt.synth import inject as INJ  # noqa: E402
from ttmt.synth import realdata as RD  # noqa: E402
from ttmt.synth import sets  # noqa: E402

EGF_KW = {"ridgecrest_2019": {}, "hinet_2012": {"iso_pre": 4.0, "iso_post": 15.0, "iso_dm": 1.0, "pre": 4.0, "post": 20.0}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=sorted(EGF_KW))
    ap.add_argument("--source", required=True, choices=("egf", "fk"))
    ap.add_argument("--per-block", type=int, default=10)
    ap.add_argument("--anchor", default="ATLAS-0", help="ATLAS-0 (anchored streams) or a development string")
    ap.add_argument("--fk-lib")
    ap.add_argument("--set-id")
    a = ap.parse_args()
    rec = RD.load_record(a.dataset)
    atoms, atom_ids = RD.template_candidates(rec)
    if a.source == "egf":
        lib = RD.egf_library(rec, **EGF_KW[a.dataset])
        src = [INJ.Source.from_egf(g, rec.fs) for g in lib]
        meta = {"kind": "egf", "n": len(src), "ids": [g.id for g in lib], "egf_blocks": sorted({g.block for g in lib}),
                "egf_kwargs": EGF_KW[a.dataset], "atoms": "realdata.template_candidates", "n_atoms": len(atoms)}
    else:
        d = pickle.load(open(a.fk_lib, "rb"))
        lib_sha = hashlib.sha256(open(a.fk_lib, "rb").read()).hexdigest()
        assert d["codes"] == rec.codes, "fk library built for another station list"
        src = [INJ.Source.from_fk(e, rec.fs) for e in d["events"]]
        meta = {"kind": "fk", "n": len(src), "fk_config": d["cfg"], "fk_seed": d["seed"],
                "library": d["library"], "library_stream": d["stream"], "library_file_sha256": lib_sha,
                # engine that computed the library (libraries built before it was stored: the engine at set build,
                # unchanged since; engine_digest_source says which)
                "engine_digest": d.get("engine_digest", fklib._engine_digest()),
                "engine_digest_source": "library" if "engine_digest" in d else "code at set build",
                "medium": d["events"][0].meta, "receiver_effects": {k: v.tolist() for k, v in d["eff"].items()},
                "atoms": "realdata.template_candidates", "n_atoms": len(atoms)}
    set_id = a.set_id or f"ss_{a.dataset}_{a.source}_v1"
    man = sets.build(set_id, rec, src, a.per_block, a.anchor, None, atoms, atom_ids, meta)
    print(set_id, len(man["realisations"]), "realisations;", len(man["partitions"]["tune"]), "tune,",
          len(man["partitions"]["read"]), "read")


if __name__ == "__main__":
    main()
