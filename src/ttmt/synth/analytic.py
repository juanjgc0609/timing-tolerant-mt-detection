"""Closed-form whole-space displacement of a point moment tensor (Aki and Richards 2002, eq. 4.29): the reference of
the engine test. Near, intermediate and far field; moment-rate function a Gaussian of unit area."""
import numpy as np

from .fk import MT_ORDER

_PQ = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))


def _gauss(t, sig):
    return np.exp(-0.5 * (t / sig) ** 2) / (np.sqrt(2 * np.pi) * sig)


def whole_space_velocity(src, rec, vp, vs, rho, t, f_gauss, t_shift):
    """(6, 3, nt) particle velocity at rec from unit basis tensors at src (MT_ORDER, symmetric pairs), moment rate
    Gaussian exp(-(f / f_gauss)^2) centred at t_shift. Units as ttmt.synth.fk (km, s, g/cm^3)."""
    sig = 1.0 / (np.sqrt(2.0) * np.pi * f_gauss)
    d = np.asarray(rec, float) - np.asarray(src, float)
    r = float(np.linalg.norm(d))
    g = d / r
    dl = np.eye(3)
    mdot = lambda tt: _gauss(tt - t_shift, sig)                                    # noqa: E731  M'(t)
    mddot = lambda tt: -(tt - t_shift) / sig ** 2 * _gauss(tt - t_shift, sig)       # noqa: E731  M''(t)
    # d/dt int_{r/a}^{r/b} tau M(t - tau) dtau = int tau M'(t - tau) dtau (Gaussian quadrature, fine)
    taus = np.linspace(r / vp, r / vs, 4001)
    wq = np.full(len(taus), taus[1] - taus[0])
    wq[0] = wq[-1] = 0.5 * wq[0]
    near = np.array([np.sum(wq * taus * mdot(tt - taus)) for tt in t])
    out = np.zeros((6, 3, len(t)))
    for j, (p0, q0) in enumerate(_PQ):
        for (p, q) in ({(p0, q0), (q0, p0)}):
            for n in range(3):
                gn, gp, gq = g[n], g[p], g[q]
                cn = 15 * gn * gp * gq - 3 * gn * dl[p, q] - 3 * gp * dl[n, q] - 3 * gq * dl[n, p]
                ip = 6 * gn * gp * gq - gn * dl[p, q] - gp * dl[n, q] - gq * dl[n, p]
                is_ = 6 * gn * gp * gq - gn * dl[p, q] - gp * dl[n, q] - 2 * gq * dl[n, p]
                fp = gn * gp * gq
                fs = gn * gp * gq - dl[n, p] * gq
                out[j, n] += (cn / r ** 4 * near
                              + ip / (vp ** 2 * r ** 2) * mdot(t - r / vp)
                              - is_ / (vs ** 2 * r ** 2) * mdot(t - r / vs)
                              + fp / (vp ** 3 * r) * mddot(t - r / vp)
                              - fs / (vs ** 3 * r) * mddot(t - r / vs)) / (4 * np.pi * rho)
    return out


__all__ = ["whole_space_velocity", "MT_ORDER"]
