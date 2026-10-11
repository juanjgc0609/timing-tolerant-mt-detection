"""Shared setup for the 2E verification (seed 999 calibration networks only). Imports the frozen 2E code read-only."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import copy
import sys

import numpy as np

HERE = "/home/claude/work"
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import task2d_run as R2  # noqa: E402
import task2e_methods as E  # noqa: E402
import gen_realistic as GR  # noqa: E402
from gen_realistic import generate_real  # noqa: E402
_ORIG_CHILD = GR._child
from task1o_consistent import draw_sources  # noqa: E402

T = R2.T
B = R2.B
D2 = R2.D2
M = R2.M
SEED = T.CAL_SEED


def setup(cell, regime, j):
    ci, nst, nclut, STA, statics, cat, drawn, corr = B.setup_net(SEED, regime, cell, j)
    grid, Gs, Qs, BD, ST, PD = R2.ops(STA, corr)
    net = E.Net(STA, corr, grid, T.DEFAULT_TAUS, Gs, Qs, BD)
    return dict(net=net, PD=PD, STA=STA, statics=statics, nclut=nclut, corr=corr, grid=grid, Gs=Gs, Qs=Qs)


def windows(S, cell, regime, j, n, with_noise_only=False):
    """windows 0..n-1 of network j (same streams as 2A/2E). With with_noise_only, also the clutter+noise-only part
    D - D_signal, where D_signal (events + coda, no clutter, no noise) is regenerated from deep copies of the streams."""
    rngs = T.window_rngs(SEED, cell, j, regime)
    out = []
    for _ in range(n):
        cp = copy.deepcopy(rngs) if with_noise_only else None
        kids = []

        def rec_child(r):
            c = _ORIG_CHILD(r)
            kids.append(copy.deepcopy(c))
            return c
        GR._child = rec_child
        try:
            srcs, taus, D, _p, tr, _a = T.generate_window(S["STA"], S["statics"], S["nclut"], T.CLUT_AMP[regime],
                                                          rngs)
        finally:
            GR._child = _ORIG_CHILD
        MF, sd = T.matched(D)
        w = {"srcs": srcs, "taus": list(taus), "D": D, "MFn": MF / sd, "sd": sd, "truth": tr}
        if with_noise_only:
            rng, rng_m, rng_h = cp[:3]
            s2, t2 = draw_sources(rng)
            assert all(np.array_equal(a, b) for a, b in zip(s2, srcs)) and list(t2) == list(taus)
            GR._child = lambda r: copy.deepcopy(kids[0])
            try:
                Df, _pp, _tt, _aa = generate_real(S["STA"], s2, t2, S["nclut"], copy.deepcopy(rng),
                                                  copy.deepcopy(rng_m), copy.deepcopy(rng_h), T.CLUT_AMP[regime],
                                                  statics=S["statics"])
                Ds, _pp, _tt, _aa = generate_real(S["STA"], s2, t2, 0, rng, rng_m, rng_h, T.CLUT_AMP[regime],
                                                  cfg={"noise": 0.0}, statics=S["statics"])
            finally:
                GR._child = _ORIG_CHILD
            w["regen_equal"] = bool(np.array_equal(Df, D))
            Dn = D - Ds
            MFs, _ = T.matched(Dn)
            w["D_noise"] = Dn
            w["MFn_noise"] = MFs / sd
            w["D_signal"] = Ds
        out.append(w)
    return out
