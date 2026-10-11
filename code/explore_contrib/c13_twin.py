"""C13 supplement: clutter-and-noise twin of every STD / DENSE window of c13_run.py (post hoc; planned after the first
look at one network).  The window W of the stream is regenerated (same streams, same draws) together with its event
part W_sig (the same generator call with no clutter and no noise, from clones of the random generators that keep
the seed sequence and its spawn counter; every window is also regenerated in full from clones and compared); the twin is
W0 = W - W_sig: the clutter and the noise of W, sample by sample, with no event.  MT-FW (cap 6) is run on W0 with the
normalisation of W.  The first two windows of every part are also re-run on W and checked against the stored target
atoms.  Used to measure how much the presence of events (sidelobes, coda, deflation residue) raises the number of false
atoms away from the events.
    python c13_twin.py [max parts]     (all parts present in out/c13_parts) -> out/c13_twin_parts/REGIME_CELL_J.json
"""
import copy
import glob
import json
import os
import time

import numpy as np

import c8_lib as C8
import c10_lib as C
import c10_run as R
import c13_lib as L
import c13_run as RR
import xc

PARTS = os.path.join(L.HERE, "out", "c13_twin_parts")


def clone(g, entropy):
    """copy of a generator INCLUDING its seed sequence and spawn counter (copy.deepcopy loses them, and the generator
    of gen_realistic spawns a child of rng_h in every window)."""
    ss = np.random.SeedSequence(entropy, n_children_spawned=g.bit_generator.seed_seq.n_children_spawned)
    bg = np.random.PCG64(ss)
    bg.state = copy.deepcopy(g.bit_generator.state)
    return np.random.Generator(bg)


def twin_stream(nc, tag, draw):
    b = [L.SEED, tag, nc.ci, nc.j, nc.ri]
    rng, rng_m, rng_h = (np.random.default_rng(b + [41]), np.random.default_rng(b + [233]),
                         np.random.default_rng(b + [61]))
    while True:
        srcs, taus = draw(rng)
        cp, cq = ([clone(g, b + [e]) for g, e in zip((rng, rng_m, rng_h), (41, 233, 61))] for _ in range(2))
        D, _p, _s, _a = C8.generate_site(nc.STA, srcs, taus, nc.nclut, rng, rng_m, rng_h, nc.clut_amp, cfg=None,
                                         statics=nc.statics)
        Dq, _p, _s, _a = C8.generate_site(nc.STA, srcs, taus, nc.nclut, cq[0], cq[1], cq[2], nc.clut_amp, cfg=None,
                                          statics=nc.statics)
        assert np.array_equal(D, Dq)                 # the clones reproduce the window exactly
        Ds, _p, _s, _a = C8.generate_site(nc.STA, srcs, taus, 0, cp[0], cp[1], cp[2], nc.clut_amp, cfg={"noise": 0.0},
                                          statics=nc.statics)
        MF, sd = xc.matched(D)
        MFs, _sd = xc.matched(Ds)
        yield {"MFn": MF / sd, "MFn0": (MF - MFs) / sd, "truth": [np.asarray(s, float) for s in srcs],
               "taus": [float(t) for t in taus], "sig_rms": float(np.sqrt((Ds ** 2).mean())),
               "tw_rms": float(np.sqrt(((D - Ds) ** 2).mean()))}


def main(limit=None):
    os.makedirs(PARTS, exist_ok=True)
    for f in sorted(glob.glob(os.path.join(RR.PARTS, "*.json")))[:limit]:
        out = os.path.join(PARTS, os.path.basename(f))
        if os.path.exists(out):
            continue
        t0 = time.time()
        p = json.load(open(f))
        S = R.Setup(p["cell"], p["net"], p["regime"], oracle=False)
        E = S.eng["cat"][0]
        rows = []
        for design in ("STD", "DENSE"):
            tag, draw, kmax, _cl = RR.DESIGNS[design]
            st = twin_stream(S.nc, tag, draw)
            stored = [r for r in p["rows"] if r["design"] == design]
            for wi, r in enumerate(stored):
                w = next(st)
                assert len(w["truth"]) == r["k"] and np.allclose(w["taus"], r["taus"])
                if wi < 2:
                    a = L.pack(C.mt_loop(E, w["MFn"], kmax))
                    assert np.allclose(np.array(a), np.array(r["atoms"]["T"])), (f, design, wi)
                rows.append({"design": design, "w": wi, "k": r["k"], "truth": r["truth"], "taus": r["taus"],
                             "atoms0": L.pack(C.mt_loop(E, w["MFn0"], kmax)) if r["k"] else None,
                             "sig_rms": w["sig_rms"], "tw_rms": w["tw_rms"]})
        with open(out + ".tmp", "w") as fh:
            json.dump({"regime": p["regime"], "cell": p["cell"], "net": p["net"], "rows": rows}, fh)
        os.replace(out + ".tmp", out)
        print(time.strftime("%H:%M:%S"), os.path.basename(f), f"{time.time() - t0:.0f} s", flush=True)
    print("ALL DONE", time.strftime("%H:%M:%S"), flush=True)


if __name__ == "__main__":
    import sys
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
