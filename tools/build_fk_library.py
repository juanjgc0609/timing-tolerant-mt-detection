"""Build a full-wave synthetic event library with known tensors at a network's stations (ttmt.synth.fklib).

    TTMT_DATA_ROOT=... python tools/build_fk_library.py DATASET N LIBRARY ZR_KM OUT.pkl [WORKERS]

LIBRARY: library index j; every draw (medium perturbation, events, receiver effects) comes from the anchored stream
ttmt.atlas.seeds.library_rng(j) (stream "inject", 1000000 + j; ADR-0009 item 5): Ridgecrest j = 0, Hi-net j = 1.

Medium: crust_a perturbed for Ridgecrest, crust_b perturbed for Hi-net (none is a detector model). ZR_KM: receiver
depth (Hi-net sensors are in boreholes; the depth is not in the registered SAC headers, 0.1 km is assumed).
"""
from ttmt.runtime import pin_blas

pin_blas()

import pickle  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402


from ttmt.synth import fklib as F  # noqa: E402
from ttmt.synth import realdata as R  # noqa: E402


def main():
    from ttmt.atlas import seeds as SEEDS
    name, n, lib, zr, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
    rng = SEEDS.library_rng(lib)
    seed = int(rng.integers(2 ** 31))                       # medium perturbation seed
    workers = int(sys.argv[6]) if len(sys.argv) > 6 else 1
    rec = R.load_record(name)
    c = F.config({"zr": zr, "medium": "crust_a" if name == "ridgecrest_2019" else "crust_b", "medium_seed": seed})
    ev = F.draw_events(rec.xyz, n, rng, c, prefix=f"fk{lib}_")
    eff = F.network_effects(len(rec.xyz), rng, c)
    t0 = time.time()
    F.synthesise(ev, rec.xyz, c, eff, workers=workers,
                 progress=lambda i, k: print(f"depth {i}/{k} {time.time() - t0:.0f} s", flush=True))
    print(name, len(ev), "events", f"{time.time() - t0:.0f} s", flush=True)
    pickle.dump({"events": ev, "eff": eff, "cfg": c, "seed": seed, "library": lib, "stream": ["inject", SEEDS.LIBRARY_OFFSET + lib],
                 "master_seed_hex": f"{SEEDS.master_seed():016x}", "codes": rec.codes, "xyz": rec.xyz},
                open(out, "wb"))


if __name__ == "__main__":
    main()
