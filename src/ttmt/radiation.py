"""Far-field radiation operator of a point moment tensor and its signal subspace (task1x_mtglrt, c9_lib, certified_fw).

For a hypothesis (node x), the P and S displacement vectors at every station are linear in the moment tensor
m = (Mxx, Myy, Mzz, Mxy, Mxz, Myz): v = G(x) m, with G of shape (2 ns 3, 6), rows ordered (phase, station, component).
The GLRT with m as a linear nuisance is the energy of the data projected on span(G): ||Q^T v||^2 with Q an orthonormal
basis of span(G) (Scharf and Friedlander 1994).
"""
import numpy as np

AP, AS = 1.0, 1.5                         # task0_common: relative P and S amplitudes


def _rand_dc(r):
    s = r.normal(size=3)
    s /= np.linalg.norm(s)
    d = r.normal(size=3)
    d -= d.dot(s) * s
    d /= np.linalg.norm(d)
    return np.outer(s, d) + np.outer(d, s)


def _radiation_constants():
    """task1w_mech.CP, CS: mean squared P and S radiation of random double couples over random directions (seed 0)."""
    r = np.random.default_rng(0)
    Ms = [_rand_dc(r) for _ in range(300)]
    n = r.normal(size=(3000, 3))
    n /= np.linalg.norm(n, axis=1)[:, None]
    cp = float(np.mean([np.mean(np.einsum("ki,ij,kj->k", n, M, n) ** 2) for M in Ms]))
    cs = float(np.mean([np.mean(np.sum((np.einsum("ij,kj->ki", M, n)
                                        - np.einsum("ki,ij,kj->k", n, M, n)[:, None] * n) ** 2, 1)) for M in Ms]))
    return cp, cs


CP, CS = _radiation_constants()


def radiation_operator(src, STA):
    """task1x_mtglrt.radiation_operator: G (2 ns 3, 6) with 1/d spreading and the relative P and S amplitudes."""
    rows = []
    for b in ("P", "S"):
        for i in range(len(STA)):
            v = STA[i] - src
            d = np.linalg.norm(v)
            n = v / d
            nx, ny, nz = n
            if b == "P":
                c = np.array([nx * nx, ny * ny, nz * nz, 2 * nx * ny, 2 * nx * nz, 2 * ny * nz])
                g = np.outer(n, c) * AP / (d * np.sqrt(CP))
            else:
                B = np.array([[nx, 0, 0, ny, nz, 0], [0, ny, 0, nx, 0, nz], [0, 0, nz, 0, nx, ny]])
                g = (np.eye(3) - np.outer(n, n)) @ B * AS / (d * np.sqrt(CS))
            rows.append(g)
    return np.vstack(rows)


def rays(XN, STA):
    """straight rays (c9_lib.rays_hom): unit vectors node -> station (N, ns, 3) and distances (N, ns)."""
    dv = np.asarray(STA, float)[None] - np.asarray(XN, float)[:, None]
    d = np.linalg.norm(dv, axis=2)
    return dv / d[..., None], d


def operators(XN, STA):
    """c9_lib.operators_hom: radiation_operator vectorised over nodes, G (N, 6 ns, 6). Equal to radiation_operator up
    to rounding (different order of operations); the registered engines use :func:`radiation_operator` per node."""
    n, d = rays(XN, STA)
    N, ns = d.shape
    nx, ny, nz = n[..., 0], n[..., 1], n[..., 2]
    c = np.stack([nx * nx, ny * ny, nz * nz, 2 * nx * ny, 2 * nx * nz, 2 * ny * nz], axis=-1)
    GP = n[..., :, None] * c[..., None, :] * (AP / np.sqrt(CP)) / d[..., None, None]
    z = np.zeros_like(nx)
    B = np.stack([np.stack([nx, z, z, ny, nz, z], -1), np.stack([z, ny, z, nx, z, nz], -1),
                  np.stack([z, z, nz, z, nx, ny], -1)], axis=-2)
    Pm = np.eye(3) - n[..., :, None] * n[..., None, :]
    GS = np.matmul(Pm, B) * (AS / np.sqrt(CS)) / d[..., None, None]
    return np.concatenate([GP.reshape(N, ns * 3, 6), GS.reshape(N, ns * 3, 6)], axis=1)


def basis(G, rtol=1e-8):
    """task1x_mtglrt.basis: orthonormal basis of span(G) (left singular vectors above rtol x the largest)."""
    U, s, _ = np.linalg.svd(G, full_matrices=False)
    return U[:, s > rtol * s.max()]


def bases(G, rtol=1e-8):
    """c9_lib.bases: basis of every node of a stack G (N, R, 6), columns below the rank zeroed (fixed width 6)."""
    U, s, _ = np.linalg.svd(G, full_matrices=False)
    keep = s > rtol * s.max(axis=1, keepdims=True)
    return U * keep[:, None, :]


def stack_Q(Qs):
    """certified_fw.stack_Q: (N, 6 ns, r) stack of per-node bases, zero-padded to the largest rank."""
    r = max(Q.shape[1] for Q in Qs)
    out = np.zeros((len(Qs), Qs[0].shape[0], r))
    for g, Q in enumerate(Qs):
        out[g, :, :Q.shape[1]] = Q
    return out


def node_operators(XN, STA):
    """per-node operators of the registered engines (r0_run.Net): Gs (N, 6 ns, 6) from radiation_operator and the
    stacked bases Qst (N, 6 ns, 6)."""
    STA = np.asarray(STA, float)
    Gs = np.stack([radiation_operator(g, STA) for g in np.asarray(XN, float)])
    return Gs, stack_Q([basis(g) for g in Gs])


def project(Q, v):
    """orthogonal projection of v (..., R) on span(Q) (Q (R, k) orthonormal): Q Q^T v."""
    return (v @ Q) @ Q.T


def subspace_energy(Q, v):
    """the GLRT statistic ||Q^T v||^2 / 2 (record units)."""
    y = v @ Q
    return (y * y).sum(-1) / 2.0
