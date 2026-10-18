"""tt_check_bh.py - independent check of explore_contrib/c18_tt.first_arrival.
(1) one-layer model must equal straight-ray time; (2) Hadley-Kanamori model against a shortest-path (Dijkstra) first
arrival on a 2-D (x, z) graph with a wide stencil (an independent, upper-bound approximation of the eikonal first
arrival that includes head waves as rays travelling just below an interface)."""
import sys

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra

sys.path.insert(0, "/home/claude/work/explore_contrib")
import c18_tt as TT  # noqa: E402


def vel(model, z):
    tops = np.array([m[0] for m in model])
    v = np.array([m[1] for m in model])
    k = np.searchsorted(tops, z, side="right") - 1
    return v[np.clip(k, 0, len(v) - 1)]


def graph_times(model, zs, xmax=60.0, zmax=20.0, h=0.25, rad=6, nsub=24):
    xs = np.arange(0.0, xmax + 1e-9, h)
    zz = np.arange(0.0, zmax + 1e-9, h)
    nx, nz = len(xs), len(zz)
    idx = lambda i, k: i * nz + k  # noqa: E731
    offs = [(a, b) for a in range(-rad, rad + 1) for b in range(-rad, rad + 1)
            if (a, b) != (0, 0) and np.gcd(abs(a), abs(b)) == 1]
    I, J, W = [], [], []
    ii, kk = np.meshgrid(np.arange(nx), np.arange(nz), indexing="ij")
    ii, kk = ii.ravel(), kk.ravel()
    u = (np.arange(nsub) + 0.5) / nsub
    for a, b in offs:
        i2, k2 = ii + a, kk + b
        ok = (i2 >= 0) & (i2 < nx) & (k2 >= 0) & (k2 < nz)
        i1, k1, i2, k2 = ii[ok], kk[ok], i2[ok], k2[ok]
        L = h * np.hypot(a, b)
        zpath = zz[k1][:, None] + (zz[k2] - zz[k1])[:, None] * u[None]
        if b == 0:  # horizontal edge on a grid line: if the line is an interface, the faster (lower) layer
            s = 1.0 / vel(model, zpath)
        else:
            s = 1.0 / vel(model, zpath)
        I.append(idx(i1, k1))
        J.append(idx(i2, k2))
        W.append(L * s.mean(1))
    A = coo_matrix((np.concatenate(W), (np.concatenate(I), np.concatenate(J))), shape=(nx * nz, nx * nz)).tocsr()
    ks = int(round(zs / h))
    assert abs(zz[ks] - zs) < 1e-9
    d = dijkstra(A, indices=idx(0, ks))
    return xs, d.reshape(nx, nz)[:, 0]


def main():
    out = []
    one = [(0.0, 6.0)]
    X = np.array([0.0, 1.0, 5.0, 20.0, 60.0])
    for zs in (3.0, 7.0):
        t = TT.first_arrival(one, zs, X)
        ref = np.hypot(X, zs) / 6.0
        out.append(f"one-layer zs {zs}: max |c18_tt - straight| = {np.abs(t - ref).max():.2e} s")
    for zs in (3.0, 7.0, 11.0):
        xs, tg = graph_times(TT.HK, zs)
        tc = TT.first_arrival(TT.HK, zs, xs)
        sel = xs <= 60.0
        d = tc[sel] - tg[sel]
        out.append(f"HK zs {zs}: c18_tt - graph over X 0-60 km: min {d.min():+.3f} max {d.max():+.3f} s "
                   f"(graph is an upper bound; c18_tt should be <= graph and close)")
        for x0 in (0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0):
            j = int(round(x0 / 0.25))
            out.append(f"   X {x0:5.1f}  c18_tt {tc[j]:7.3f}  graph {tg[j]:7.3f}  diff {tc[j] - tg[j]:+.3f}")
    print("\n".join(out))


if __name__ == "__main__":
    main()
