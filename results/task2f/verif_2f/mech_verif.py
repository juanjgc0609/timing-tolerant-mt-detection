"""mech_verif.py - exploratory mechanism checks for Task 2F on seed-999 dense windows (verification; not registered).

Windows: the dense pilot windows of seed 999 (networks 10-15, TAG 3002, c2_run.draw_dense, kmax 7), regenerated with
task2f_run.setup_net and a copy of task2f_run.window_stream that also accepts a generator cfg and a clutter switch.
Conditions (generator):  base (as registered), noclut (n_clutter = 0), noclut_nosite (n_clutter = 0, site_sd = 0).
Methods per window:
  MT-FW at CAPV in {2, 3, 6, 12, 24} and uncapped (own re-implementation of c6_lib.mtfw_dense with the cap level as a
  parameter; checked to reproduce c6_lib.mtfw_dense at CAPV 6 and cap=False);
  swap variants: SELcap/SCOREnocap and SELnocap/SCOREcap (selection with one statistic, the recorded score from the
  other statistic at the selected (node, origin), on the same residual; deflation follows the selection statistic);
  BOMP-PK+SC (c6_lib.bomp_pk_m).
Output: verif_2f/mech_parts/{cond}_{regime}_{cell}_net{j}.json
    OPENBLAS_CORETYPE=Prescott python verif_2f/mech_verif.py run [--workers 1]
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
for p in (os.path.join(WORK, "explore_contrib"), WORK):
    if p not in sys.path:
        sys.path.insert(0, p)

import task2f_run as R  # noqa: E402
import c6_lib as L  # noqa: E402
import c2_run as C2  # noqa: E402
import xc  # noqa: E402
from gen_realistic import generate_real  # noqa: E402
from task1o_consistent import matched  # noqa: E402
from task1y_clutter import CLUT_AMP  # noqa: E402

X, M, D2 = L.X, L.M, L.D2
NT = L.NT
OUT = os.path.join(HERE, "mech_parts")
CAPS = (2.0, 3.0, 6.0, 12.0, 24.0)
CONDS = {"base": (True, None), "noclut": (False, None), "noclut_nosite": (False, {"site_sd": 0.0})}
NETS = tuple(range(10, 16))
KMAX = 7


def leaf(R_, capv):
    Et = np.sum(np.ascontiguousarray(R_.transpose(0, 2, 1)) ** 2, axis=-1)
    W = Et[:, D2.II_TAB]
    a = W.argmax(axis=-1)
    K = D2.II_TAB[np.arange(NT)[None, :], a]
    V = np.ascontiguousarray(R_.transpose(0, 2, 1))
    if capv is not None:
        nv = np.linalg.norm(V, axis=-1, keepdims=True)
        V = V * np.minimum(1.0, capv / np.maximum(nv, 1e-12))
    return K, V


def mtfw(ops, MFn, truth, kmax=KMAX, sel=6.0, score="same", min_sep=5.0):
    """sel: cap level of the selection statistic (None = uncapped); score: 'same' or a cap level / None for the
    recorded score.  Deflation uses the selection statistic's (capped) peak vector, as c6_lib.mtfw_dense."""
    ns = len(ops.STA)
    R_ = MFn.copy()
    garr, IDX = ops.garr, ops.IDX
    N, T = IDX.shape[:2]
    st = np.arange(ns)[None, None, None, :]
    excluded = np.zeros(N, bool)
    rec, run = [], np.inf
    for _ in range(kmax):
        if excluded.all():
            break
        K, V = leaf(R_, sel)
        tpk = K[st, IDX]
        vf = V[st, tpk].reshape(N, T, 6 * ns)
        y = np.einsum("gtr,grk->gtk", vf, ops.Qst)
        sc = (y * y).sum(-1) / 2.0
        sc[excluded] = -np.inf
        best = sc.max(axis=1)
        gi = int(np.argmax(best))
        tj = int(np.argmax(sc[gi]))
        s = float(sc[gi, tj])
        if not np.isfinite(s):
            break
        if score != "same":
            _K2, V2 = leaf(R_, score)
            v2 = V2[np.arange(ns)[None, :], tpk[gi, tj]].reshape(6 * ns)       # same peak samples (K identical)
            y2 = v2 @ ops.Qst[gi]
            s = float(y2 @ y2) / 2.0
        src = garr[gi]
        run = min(run, s)
        rec.append((float(run), int(M.hit_of(src, truth))))
        excluded = excluded | (np.linalg.norm(garr[:, :2] - src[:2], axis=1) <= min_sep)
        G = ops.Gs[gi]
        m, *_ = np.linalg.lstsq(G, vf[gi, tj], rcond=None)
        pred = (G @ m).reshape(2, ns, 3)
        tp = tpk[gi, tj]
        for b in range(2):
            for i in range(ns):
                t0 = int(tp[b, i])
                lo, hi = max(0, t0 - X.AC_HALF), min(NT, t0 + X.AC_HALF + 1)
                R_[i, :, lo:hi] -= np.outer(pred[b, i], X.AC[X.AC_HALF - (t0 - lo):X.AC_HALF + (hi - t0)])
    return rec


def stream(seed, net, tag, draw, clutter, cfg):
    b = [seed, tag, net["ci"], net["j"], net["ri"]]
    rng, rng_m, rng_h = (np.random.default_rng(b + [41]), np.random.default_rng(b + [233]),
                         np.random.default_rng(b + [61]))
    while True:
        srcs, taus = draw(rng)
        D, _p, _t, _a = generate_real(net["STA"], srcs, taus, net["nclut"] if clutter else 0, rng, rng_m, rng_h,
                                      CLUT_AMP[net["regime"]], cfg=cfg, statics=net["statics"])
        MF, sd = matched(D)
        yield {"MFn": MF / sd, "truth": [np.asarray(s, float) for s in srcs]}


def job(args):
    cond, regime, cell, j = args
    out = os.path.join(OUT, f"{cond}_{regime}_{cell}_net{j}.json")
    if os.path.exists(out):
        return out
    clutter, cfg = CONDS[cond]
    net = R.setup_net(999, cell, j, regime)
    net.update({"j": j, "regime": regime})
    ops = L.ModelOps(net["STA"], xc.make_grid(), L.Homog(), net["corr"])
    st = stream(999, net, C2.TAG, C2.draw_dense, clutter, cfg)
    wins = []
    for wi in range(10):
        w = next(st)
        MFn, tr = w["MFn"], w["truth"]
        rec = {}
        for c in CAPS:
            rec[f"MTFW_cap{c:g}"] = mtfw(ops, MFn, tr, sel=c)
        rec["MTFW_nocap"] = mtfw(ops, MFn, tr, sel=None)
        rec["SELcap6_SCOREnocap"] = mtfw(ops, MFn, tr, sel=6.0, score=None)
        rec["SELnocap_SCOREcap6"] = mtfw(ops, MFn, tr, sel=None, score=6.0)
        rec["BOMP-PK"] = L.bomp_pk_m(ops, MFn, tr, KMAX)[0]
        if cond == "base" and wi < 3:
            ref6 = L.mtfw_dense(ops, MFn, tr, KMAX)[0]
            refn = L.mtfw_dense(ops, MFn, tr, KMAX, cap=False)[0]
            rec["_check"] = [bool(np.allclose(np.array(ref6, float), np.array(rec["MTFW_cap6"], float), rtol=0,
                                              atol=1e-9)),
                             bool(np.allclose(np.array(refn, float), np.array(rec["MTFW_nocap"], float), rtol=0,
                                              atol=1e-9))]
        wins.append({"k": len(tr), "rec": {k: [list(map(float, x)) if not isinstance(x, bool) else x for x in v]
                                           for k, v in rec.items()}})
    with open(out, "w") as fh:
        json.dump({"cond": cond, "regime": regime, "cell": cell, "net": j, "windows": wins}, fh)
    return out


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    workers = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 1
    conds = sys.argv[sys.argv.index("--conds") + 1].split(",") if "--conds" in sys.argv else list(CONDS)
    jobs = [(c, r, cell, j) for c in conds for r in ("thin", "strong") for cell in C2.CELLS for j in NETS]
    t0 = time.time()
    with ProcessPoolExecutor(workers) as ex:
        for i, o in enumerate(ex.map(job, jobs)):
            print(f"{i + 1}/{len(jobs)} {os.path.basename(o)} {time.time() - t0:.0f} s", flush=True)
