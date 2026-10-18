"""
tail_check.py - numerical check of Proposition AA.9.5 (false-alarm tail of the peak-searched, radially transformed
moment-tensor score under Gaussian noise).  White Gaussian noise -> the 2 Hz Ricker matched filter of the project
(task1o_consistent._KER, 16 Hz) -> unit variance; n_b = 14 blocks (7 stations x P, S) of 3 components, each read at
the sample of largest energy within +-W = 4 samples of its own predicted sample (blocks 40 samples apart, so
independent); a random 6-dimensional orthonormal Q.  Compares the empirical upper quantiles of sigma = ||Q^T phi(v)||
with the bound  sqrt(6) + kappa * (sqrt(2 n_b log(2W+1)) + sqrt(2 log(1/p))).
    python tail_check.py OUT.txt
"""
import sys

import numpy as np

sys.path.insert(0, "/home/claude/work")
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
    lines = [f"kappa (max |K(f)| / ||K||) = {kap:.3f}; n_b {NB}, W {W}, cap {C}, {REPS} replicates; per-hypothesis "
             f"bound sigma_p = sqrt(6) + kappa (sqrt(2 n_b log(2W+1)) + sqrt(2 log(1/p)))"]
    lines.append("   p        bound   empirical quantile of sigma: identity   cap   log")
    for p in (1e-1, 1e-2, 1e-3, 1e-4):
        b = np.sqrt(6) + kap * (np.sqrt(2 * NB * np.log(L)) + np.sqrt(2 * np.log(1 / p)))
        q = {m: float(np.quantile(np.concatenate(v), 1 - p)) for m, v in out.items()}
        lines.append(f"   {p:.0e}   {b:6.2f}   {q['identity']:6.2f}  {q['cap']:6.2f}  {q['log']:6.2f}")
    s = {m: np.concatenate(v) for m, v in out.items()}
    lines.append("mean sigma: " + "  ".join(f"{m} {s[m].mean():.3f}" for m in s) + f"; sqrt(6) = {np.sqrt(6):.3f}")
    open(outp, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1])
