"""Synthetic seeds of ATLAS-0, derived after the anchor (ADR-0004 item 4).

The master seed is the first 8 bytes, read big-endian, of SHA-256(<tag commit hash> || <Release created_at>): the
40-character hexadecimal commit the registration tag points to, followed with no separator by the `created_at`
string of its GitHub Release exactly as the API returns it. Every random stream of the atlas is a child of that
master by a fixed name and integer indices (numpy SeedSequence spawn keys), so a stream never depends on the order in
which other streams are drawn.

Named streams (card, "Tuning" and "Truth"):
- "tuning", trial: the configuration drawn in trial 0..59 of the random search, the same list for every method;
- "inject", (set key, block, r): the injected events and timing errors of realisation r on noise block `block` of a
  semi-synthetic set (`realisation_rng`); "inject", 1000000 + j: event library j, drawn whole
  (`library_rng`). A realisation key has three indices and a library key one, so they never coincide;
- "subset", (set key, block, r): the contamination station draws of that realisation (`realisation_rng`);
- "bootstrap", analysis: resampling of blocks for the intervals (0: froc.describe, 1: froc.pooled_loss);
- "signflip", analysis: the Monte Carlo signs of ttmt.stats.signflip_pvalue when the exact test is too large
  (0: froc.gate, when no signs are passed).
The layout of these streams was fixed in PR #16, before any atlas data were drawn.
Realisation keys, library keys and the partition rule were added by amendment 1 (ADR-0009). The partition of noise
blocks into tuning and reading (30/70) is `ttmt.synth.realdata.tune_read` and does not use these streams."""
import hashlib
import json
from pathlib import Path

import numpy as np

TAG = "reg/ATLAS-0/registration-1"
TAG_COMMIT = "7bafb7abd308d0376ff1176d3d7abeb989aa0564"
RELEASE_CREATED_AT = "2026-10-24T22:37:39Z"
STREAMS = ("tuning", "inject", "subset", "bootstrap", "signflip")
N_TRIALS = 60
SEEDS_FILE = Path(__file__).resolve().parents[3] / "experiments" / "exploratory" / "ATLAS-0" / "seeds.json"


def master_seed(commit=TAG_COMMIT, created_at=RELEASE_CREATED_AT):
    """first 8 bytes, big-endian, of SHA-256(commit || created_at) (ADR-0004 item 4)."""
    if len(commit) != 40 or any(c not in "0123456789abcdef" for c in commit):
        raise ValueError("commit must be the full 40-character lowercase hexadecimal hash")
    return int.from_bytes(hashlib.sha256((commit + created_at).encode("ascii")).digest()[:8], "big")


def _name_key(name):
    if name not in STREAMS:
        raise KeyError(f"unknown stream {name!r}; registered streams: {STREAMS}")
    return int.from_bytes(hashlib.sha256(name.encode("ascii")).digest()[:4], "big")


def seed_sequence(name, *index, master=None):
    """SeedSequence of stream `name` at the integer indices `index` (e.g. ("inject", record))."""
    master = master_seed() if master is None else master
    return np.random.SeedSequence(master, spawn_key=(_name_key(name), *(int(i) for i in index)))


def rng(name, *index, master=None):
    """numpy Generator (PCG64) of stream `name` at `index`."""
    return np.random.default_rng(seed_sequence(name, *index, master=master))


def set_key(set_id):
    """integer key of a semi-synthetic set: the first 4 bytes, big-endian, of sha256(set_id)."""
    return int.from_bytes(hashlib.sha256(str(set_id).encode("utf-8")).digest()[:4], "big")


def realisation_rng(set_id, block, r, stream="inject", master=None):
    """generator of realisation r on noise block number `block` (the block starting at record t0 + 300 block s) of
    set `set_id` (ADR-0009 item 5): stream "inject" for its events, "subset" for its contamination stations."""
    if stream not in ("inject", "subset"):
        raise KeyError("a realisation draws only from the streams 'inject' and 'subset'")
    if int(block) < 0 or int(r) < 0:
        raise ValueError("block and r are non-negative")
    return rng(stream, set_key(set_id), int(block), int(r), master=master)


LIBRARY_OFFSET = 1_000_000


def library_rng(j, master=None):
    """generator of event library j, drawn whole (medium, event sequence, receiver effects; ADR-0009 item 5):
    stream ("inject", 1000000 + j)."""
    if int(j) < 0:
        raise ValueError("library index is non-negative")
    return rng("inject", LIBRARY_OFFSET + int(j), master=master)


def tuning_seeds(n=N_TRIALS, master=None):
    """the seed list of the random search: one 32-bit seed per trial, shared by every method."""
    return [int(seed_sequence("tuning", t, master=master).generate_state(1)[0]) for t in range(n)]


def record(master=None):
    """the derivation as stored in seeds.json."""
    m = master_seed() if master is None else master
    return dict(adr="ADR-0004 item 4", tag=TAG, tag_commit=TAG_COMMIT, release_created_at=RELEASE_CREATED_AT,
                rule="int.from_bytes(sha256((tag_commit + release_created_at).encode('ascii')).digest()[:8], 'big')",
                master_seed=m, master_seed_hex=f"{m:016x}", streams=list(STREAMS),
                stream_rule="numpy.random.SeedSequence(master_seed, spawn_key=(first 4 bytes big-endian of "
                            "sha256(name), *indices))",
                tuning_seeds=tuning_seeds(master=m))


def check_file(path=SEEDS_FILE):
    """True when the stored seeds.json matches the derivation."""
    return json.loads(Path(path).read_text()) == record()


if __name__ == "__main__":
    import sys
    if "--check" in sys.argv:
        ok = check_file()
        print("seeds.json matches the derivation" if ok else "seeds.json does NOT match the derivation")
        sys.exit(0 if ok else 1)
    SEEDS_FILE.write_text(json.dumps(record(), indent=1) + "\n")
    print(f"wrote {SEEDS_FILE} (master {record()['master_seed_hex']})")
