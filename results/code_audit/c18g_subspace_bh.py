"""c18g_subspace_bh.py - C18g down-weights the S peak vectors (w_S = 0.5 or 0) but keeps the orthonormal basis Q of the
FULL P+S operator G.  For w_S < 1 the matching statistic would project onto the span of the equally weighted operator
[G_P; w_S G_S] (P-only: span(G_P)).  This script measures, on the R0 N8 geometry, for nodes near the CAL target events
and random double couples, the signal energy kept and the noise mean of the two statistics (no cap, no noise draw):
  signal_full = |Q^T [v_P; w v_S]|^2, signal_ok = |Q_w^T [v_P; w v_S]|^2, v = G m
  noise mean  = trace of the projector restricted to the weighted rows (unit white noise on every component).
Ratio of deflections (signal / noise mean) of the C18g statistic to the matched one."""
import sys

import numpy as np

sys.path[:0] = ["/home/claude/work/r0", "/home/claude/work/explore_contrib", "/home/claude/work"]
import r0_lib as R  # noqa: E402
import task1x_mtglrt as X  # noqa: E402
from task1w_mech import rand_dc  # noqa: E402


def orth(G):
    U, s, _ = np.linalg.svd(G, full_matrices=False)
    return U[:, s > 1e-8 * s.max()]


def main():
    centre = R.frame()
    cat = R.catalogue_local(centre)
    out = []
    for sub in ("N8", "U8"):
        STA, _ = R.stations(R.SUBNETS[sub], centre)
        ns = len(STA)
        tg = [e for e in cat if e["t"] < 900 and e["M"] >= 2.0]
        rng = np.random.default_rng(1)
        for w in (0.5, 0.0):
            rs, rn, rd = [], [], []
            for e in tg:
                src = np.array([e["x"], e["y"], max(e["z"], 1.0)])
                G = X.radiation_operator(src, STA)                       # (6 ns, 6), rows (b, i, c)
                wr = np.r_[np.ones(3 * ns), np.full(3 * ns, w)]
                Q = orth(G)
                Gw = G * wr[:, None]
                Qw = orth(Gw[wr > 0]) if w == 0 else orth(Gw)
                for _ in range(50):
                    Mt = rand_dc(rng)
                    m = np.array([Mt[0, 0], Mt[1, 1], Mt[2, 2], Mt[0, 1], Mt[0, 2], Mt[1, 2]])
                    v = (G @ m) * wr
                    sf = float(np.sum((Q.T @ v) ** 2))
                    so = float(np.sum((Qw.T @ (v[wr > 0] if w == 0 else v)) ** 2))
                    nf = float(np.sum((Q * wr[:, None]) ** 2))         # E|Q^T (w n)|^2
                    no = float(Qw.shape[1]) if w == 0 else float(np.sum((Qw * wr[:, None]) ** 2))
                    rs.append(sf / so)
                    rn.append(nf / no)
                    rd.append((sf / nf) / (so / no))
            out.append(f"{sub} w_S {w}: signal kept by C18g statistic / matched {np.median(rs):.3f} "
                       f"[{np.percentile(rs, 10):.3f}, {np.percentile(rs, 90):.3f}]; noise mean ratio {np.median(rn):.3f}; "
                       f"deflection ratio {np.median(rd):.3f} [{np.percentile(rd, 10):.3f}, {np.percentile(rd, 90):.3f}]")
    print("\n".join(out))


if __name__ == "__main__":
    main()
