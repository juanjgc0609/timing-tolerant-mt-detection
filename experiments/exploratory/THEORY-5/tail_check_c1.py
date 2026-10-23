"""
tail_check_c1.py - code/theory_tail_check.py with the bound of Proposition AB.C.1 printed next to the bound of
AA.9.5(b).  Same kernel (task1o_consistent._KER, 2 Hz Ricker at 16 Hz), configuration, seed and simulation as the
frozen script, so the empirical quantiles must equal results/theory/tail_check_out.txt.  code/ is append-only, so the
corrected check lives here rather than as an edit of the frozen script.
    AA.9.5(b):  sqrt(6) + kappa (sqrt(2 n_b log(2W+1)) + sqrt(2 log(1/p)))
    AB.C.1:     sqrt(6) +        sqrt(2 n_b log(2W+1)) + sqrt(2 log(1/p))
    python tail_check_c1.py OUT.txt
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "code"))
from task1o_consistent import _KER  # noqa: E402

NB, W, C, REPS = 14, 4, 6.0, 200000


def kappa(k, n=1 << 16):
    K = np.abs(np.fft.rfft(k, n))
    return float(K.max() / np.linalg.norm(k))


def radial(v, kind):
    t = np.linalg.norm(v, axis=-1, keepdims=True)
    if kind == "identity":
        return v
    if kind == "cap":
        return v * np.minimum(1.0, C / np.maximum(t, 1e-12))
    r = C * np.sqrt(np.log1p(t * t / (C * C)))
    return v * r / np.maximum(t, 1e-12)


def main(outp):
    rng = np.random.default_rng(20261017)
    k = _KER / np.linalg.norm(_KER)
    kap = kappa(k)
    L = 2 * W + 1
    Q = np.linalg.qr(rng.normal(size=(3 * NB, 6)))[0]
    out = {m: [] for m in ("identity", "cap", "log")}
    batch = 5000
    for _ in range(REPS // batch):
        z = rng.normal(size=(batch, NB, 3, L + len(k) - 1))
        x = np.stack([np.apply_along_axis(lambda s: np.convolve(s, k, mode="valid"), -1, z[:, :, c]) for c in range(3)],
                     axis=-1)                                                     # (batch, NB, L, 3)
        e = (x * x).sum(-1)
        sel = np.take_along_axis(x, e.argmax(-1)[..., None, None].repeat(3, -1), axis=2)[:, :, 0]   # (batch, NB, 3)
        for m in out:
            v = radial(sel, m).reshape(batch, -1)
            out[m].append(np.linalg.norm(v @ Q, axis=1))
    lines = [f"kappa (max |K(f)| / ||K||) = {kap:.3f}; n_b {NB}, W {W}, cap {C}, {REPS} replicates; bounds: AA.9.5(b) "
             f"sqrt(6) + kappa (sqrt(2 n_b log(2W+1)) + sqrt(2 log(1/p))), AB.C.1 the same with kappa = 1"]
    lines.append("   p        AA.9.5b  AB.C.1   empirical quantile of sigma: identity   cap   log   ratio AB.C.1/identity")
    for p in (1e-1, 1e-2, 1e-3, 1e-4):
        a = np.sqrt(2 * NB * np.log(L)) + np.sqrt(2 * np.log(1 / p))
        b_old, b_new = np.sqrt(6) + kap * a, np.sqrt(6) + a
        q = {m: float(np.quantile(np.concatenate(v), 1 - p)) for m, v in out.items()}
        lines.append(f"   {p:.0e}   {b_old:6.2f}   {b_new:6.2f}   {q['identity']:6.2f}  {q['cap']:6.2f}  {q['log']:6.2f}"
                     f"   {b_new / q['identity']:5.2f}")
    s = {m: np.concatenate(v) for m, v in out.items()}
    lines.append("mean sigma: " + "  ".join(f"{m} {s[m].mean():.3f}" for m in s) + f"; sqrt(6) = {np.sqrt(6):.3f}")
    open(outp, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1])
