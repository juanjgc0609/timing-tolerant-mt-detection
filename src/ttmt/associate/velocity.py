"""One 1-D velocity model, handed to every comparator that takes one (card: "the same 1-D model for every method
that takes one").

The model is the layered model of :mod:`ttmt.traveltime` (layers (top depth km, VP km/s), the last a half-space, VS =
VP / vpvs). Each program receives it in its own form: GaMMA as the layer profile of its eikonal solver, PyOcto as a
travel-time table built by ``VelocityModel1D.create_model`` from the same profile, REAL as a ``ttdb`` table computed
here with :func:`ttmt.traveltime.first_arrival`. Methods that only take a homogeneous medium get the slowness-averaged
velocities of :meth:`Velocity1D.mean_velocities`.
"""
from dataclasses import dataclass

import numpy as np

from .. import traveltime as TT

KM_PER_DEG = 111.19            # REAL's conversion (station.dat in degrees)


@dataclass(frozen=True)
class Velocity1D:
    layers: tuple
    vpvs: float = TT.R_VPVS
    name: str = ""

    def __post_init__(self):
        object.__setattr__(self, "layers", tuple((float(z), float(v)) for z, v in self.layers))
        if self.layers[0][0] != 0.0 or any(b[0] <= a[0] for a, b in zip(self.layers, self.layers[1:])):
            raise ValueError("layers must start at 0 km and increase in depth")

    @classmethod
    def hadley_kanamori(cls):
        return cls(tuple(TT.HK), name="Hadley-Kanamori 1977")

    @classmethod
    def jma2001(cls):
        return cls(tuple(TT.JMA2001), name="JMA2001 (Maeda 28-point approximation)")

    def profile(self, zmax):
        """depth, VP, VS samples of the step profile down to zmax (two samples per interface)."""
        z, vp = [], []
        for k, (top, v) in enumerate(self.layers):
            if top >= zmax:
                break
            bot = min(self.layers[k + 1][0] if k + 1 < len(self.layers) else zmax, zmax)
            z += [top, max(top, bot - 1e-3)]
            vp += [v, v]
        vp = np.array(vp)
        return np.array(z), vp, vp / self.vpvs

    def mean_velocities(self, zmax):
        """VP and VS that give the vertical travel time through [0, zmax] (harmonic mean by thickness)."""
        tops = [z for z, _ in self.layers] + [np.inf]
        t = 0.0
        for k, (top, v) in enumerate(self.layers):
            h = max(0.0, min(tops[k + 1], zmax) - top)
            t += h / v
        vp = zmax / t
        return float(vp), float(vp / self.vpvs)

    def times(self, src, sta):
        """(2, ns) P and S first arrivals from src (x, y, z km) to the stations sta (ns, 3) at the surface."""
        return TT.times(list(self.layers), self.vpvs, src, sta)

    def gamma_eikonal(self, xlim, ylim, zlim, h=0.5):
        """the "eikonal" entry of a GaMMA configuration."""
        z, vp, vs = self.profile(float(zlim[1]) + 1.0)
        return {"vel": {"z": z.tolist(), "p": vp.tolist(), "s": vs.tolist()}, "h": float(h),
                "xlim": list(map(float, xlim)), "ylim": list(map(float, ylim)), "zlim": list(map(float, zlim))}

    def pyocto_model(self, path, xdist, zdist, delta=0.25):
        """write PyOcto's binary 1-D table to path (once; reused when it exists) and return the path."""
        import os

        import pandas as pd
        import pyocto
        if not os.path.exists(path):
            z, vp, vs = self.profile(float(zdist) + 1.0)
            pyocto.VelocityModel1D.create_model(pd.DataFrame({"depth": z, "vp": vp, "vs": vs}), float(delta),
                                                float(xdist), float(zdist), str(path))
        return path

    def real_ttdb(self, path, xmax_deg, zmax, dx_deg=0.01, dz=1.0):
        """write REAL's ttdb.txt: for every depth step dz in [0, zmax] km and epicentral distance step dx_deg in
        (0, xmax_deg], the P and S first arrivals, their derivatives in distance (s/deg) and in source depth (s/km), by
        central differences of :func:`ttmt.traveltime.first_arrival`. REAL uses the two derivatives as the linear
        correction t(GCarc, dep) = t + dGCarc * prayp + ddep * hslow (REAL.c), so they must be the true derivatives."""
        xs = np.arange(1, int(round(xmax_deg / dx_deg)) + 1) * dx_deg
        X = xs * KM_PER_DEG
        model = list(self.layers)
        e = 1e-3
        with open(path, "w") as fh:
            for zs in np.arange(0.0, zmax + 1e-9, dz):
                zq = max(zs, 0.01)          # first_arrival has no direct ray from a source exactly at the surface
                cols = []
                for r in (1.0, self.vpvs):
                    t = TT.first_arrival(model, zq, X, r)
                    p = (TT.first_arrival(model, zq, X + e, r) - TT.first_arrival(model, zq, X - e, r)) / (2 * e)
                    z1, z0 = zq + e, zq - e
                    hz = (TT.first_arrival(model, z1, X, r) - TT.first_arrival(model, z0, X, r)) / (z1 - z0)
                    cols.append((t, p * KM_PER_DEG, hz))
                for i, x in enumerate(xs):
                    (tp, pp, hp), (ts, ps, hs) = [(c[0][i], c[1][i], c[2][i]) for c in cols]
                    fh.write(f"{x:.4f} {zs:g} {tp:.6f} {ts:.6f} {pp:.6f} {ps:.6f} {hp:.6f} {hs:.6f} P S\n")
        return path
