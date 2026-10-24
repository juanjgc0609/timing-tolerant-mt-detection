"""Flat layered elastic media with a free surface and constant-Q attenuation (synthetic generator v2).

A :class:`LayeredModel` lists layers from the surface down; the last is a half-space. Depth z is positive downward
(km), velocities km/s, density g/cm^3. Quality factors are frequency independent (Kjartansson 1979 constant Q), with
the velocities given at the reference frequency ``f_ref``; the complex velocity at frequency f is
    c(f) = c_ref [1 + ln(f / f_ref) / (pi Q)] (1 - i / (2 Q))
(first order in 1/Q; time convention exp(-i omega t), so a wave exp(i omega z / c) decays with distance).

Independence (v2 rule): this module and the rest of ``ttmt.synth`` import nothing from the detector (``ttmt.kernel``,
``ttmt.radiation``, ``ttmt.traveltime``, ``ttmt.detector``) nor from ``code/``. The presets below are not the
detector's models (homogeneous VP 6 / VP/VS 1.73, Hadley-Kanamori, JMA2001).
"""
import hashlib
import json
from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class LayeredModel:
    """thickness: km for every layer but the last (the half-space); vp, vs, rho, qp, qs: one value per layer."""
    thickness: tuple
    vp: tuple
    vs: tuple
    rho: tuple
    qp: tuple = None
    qs: tuple = None
    f_ref: float = 1.0
    name: str = ""
    tops: np.ndarray = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        n = len(self.vp)
        if not (len(self.thickness) == n - 1 and len(self.vs) == n and len(self.rho) == n):
            raise ValueError("thickness needs one value per layer except the half-space")
        if any(h <= 0 for h in self.thickness):
            raise ValueError("thicknesses must be positive")
        if any(b <= 0 or a <= b * np.sqrt(4.0 / 3.0) for a, b in zip(self.vp, self.vs)):
            raise ValueError("need vp > sqrt(4/3) vs > 0 (positive bulk modulus)")
        for k in ("qp", "qs"):
            v = getattr(self, k)
            object.__setattr__(self, k, tuple([1e9] * n) if v is None else tuple(float(x) for x in v))
        for k in ("thickness", "vp", "vs", "rho"):
            object.__setattr__(self, k, tuple(float(x) for x in getattr(self, k)))
        object.__setattr__(self, "tops", np.r_[0.0, np.cumsum(self.thickness)])

    @property
    def n(self):
        return len(self.vp)

    def layer_of(self, z):
        """index of the layer holding depth z (a depth on an interface belongs to the layer below)."""
        return int(np.searchsorted(self.tops, z, side="right") - 1)

    def to_dict(self):
        return {"name": self.name, "thickness": list(self.thickness), "vp": list(self.vp), "vs": list(self.vs),
                "rho": list(self.rho), "qp": list(self.qp), "qs": list(self.qs), "f_ref": self.f_ref}

    @classmethod
    def from_dict(cls, d):
        return cls(**{k: (tuple(v) if isinstance(v, list) else v) for k, v in d.items()})

    def digest(self):
        """sha256 of the canonical JSON of the model (keys of GF caches and of dataset manifests)."""
        return hashlib.sha256(json.dumps(self.to_dict(), sort_keys=True).encode()).hexdigest()

    def split(self, z):
        """the same medium with an extra interface at depth z (no-op if z is already an interface) and the index of
        the layer that starts at z. Used to place the source on an interface of the global matrix."""
        if z <= 0:
            raise ValueError("source must be below the free surface")
        tops = list(self.tops)
        k = self.layer_of(z)
        if abs(z - tops[k]) < 1e-9:
            return self, k
        th = list(self.thickness)
        if k < self.n - 1:
            th[k:k + 1] = [z - tops[k], tops[k + 1] - z]
        else:
            th.append(z - tops[k])
        dup = lambda v: tuple(v[:k + 1]) + tuple(v[k:])          # noqa: E731
        m = LayeredModel(tuple(th), dup(self.vp), dup(self.vs), dup(self.rho), dup(self.qp), dup(self.qs),
                         self.f_ref, self.name)
        return m, k + 1

    def complex_velocities(self, f):
        """(alpha, beta) complex velocities, arrays (len(f), n), at real frequencies f (Hz, > 0)."""
        f = np.maximum(np.abs(np.atleast_1d(np.asarray(f, float))), 1e-6)[:, None]
        lf = np.log(f / self.f_ref)
        out = []
        for v, q in ((self.vp, self.qp), (self.vs, self.qs)):
            v, q = np.asarray(v)[None], np.asarray(q)[None]
            out.append(v * (1.0 + lf / (np.pi * q)) * (1.0 - 0.5j / q))
        return out[0], out[1]

    def perturbed(self, rng, v_sd=0.03, h_sd=0.10, vpvs_sd=0.02, name=None):
        """a random neighbour of the model: layer velocities x (1 + N(0, v_sd)), thicknesses x (1 + N(0, h_sd)),
        VP/VS x (1 + N(0, vpvs_sd)), per layer; the ordering vp > sqrt(4/3) vs is kept."""
        n = self.n
        f = 1.0 + rng.normal(0.0, v_sd, n)
        r = (np.asarray(self.vp) / np.asarray(self.vs)) * (1.0 + rng.normal(0.0, vpvs_sd, n))
        r = np.maximum(r, 1.45)
        vp = np.asarray(self.vp) * f
        th = np.asarray(self.thickness) * np.clip(1.0 + rng.normal(0.0, h_sd, n - 1), 0.5, 1.5)
        return LayeredModel(tuple(th), tuple(vp), tuple(vp / r), self.rho, self.qp, self.qs, self.f_ref,
                            name or (self.name + "~"))


# Presets. None is a detector model. Values are generic upper-crust profiles chosen to differ from the detector's
# homogeneous VP 6.0 / VS 3.47 and from Hadley-Kanamori; Q from typical crustal values (Qs ~ 50-600, Qp ~ 2 Qs).
PRESETS = {
    # moderately soft top layer, two crustal interfaces
    "crust_a": LayeredModel(thickness=(1.0, 3.0, 8.0, 12.0), vp=(3.6, 5.2, 5.9, 6.4, 7.0),
                            vs=(2.0, 3.0, 3.4, 3.7, 4.0), rho=(2.30, 2.55, 2.70, 2.80, 2.95),
                            qp=(80.0, 300.0, 600.0, 900.0, 1200.0), qs=(40.0, 150.0, 300.0, 450.0, 600.0),
                            name="crust_a"),
    # carbonate-like stiff upper crust with a velocity inversion (Apennine-type), higher VP/VS at depth
    "crust_b": LayeredModel(thickness=(1.5, 4.5, 3.0, 9.0), vp=(4.2, 6.1, 5.6, 6.3, 6.9),
                            vs=(2.3, 3.3, 3.0, 3.5, 3.9), rho=(2.45, 2.70, 2.60, 2.78, 2.90),
                            qp=(100.0, 500.0, 400.0, 800.0, 1000.0), qs=(50.0, 250.0, 200.0, 400.0, 500.0),
                            name="crust_b"),
    # strong gradient approximated by many layers (no low-velocity zone)
    "crust_c": LayeredModel(thickness=(0.8, 1.2, 2.0, 3.0, 4.0, 6.0), vp=(3.9, 4.6, 5.3, 5.8, 6.1, 6.5, 7.2),
                            vs=(2.1, 2.6, 3.0, 3.3, 3.5, 3.75, 4.1), rho=(2.35, 2.45, 2.60, 2.70, 2.75, 2.85, 3.0),
                            qp=(90.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0),
                            qs=(45.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0), name="crust_c"),
}


def whole_space(vp=6.0, vs=3.5, rho=2.7):
    """homogeneous medium for the analytic tests (with free_surface=False in the engine it is a whole space)."""
    return LayeredModel(thickness=(), vp=(vp,), vs=(vs,), rho=(rho,), name="whole_space")


# ---------------------------------------------------------------------------------------------- ray theory (truth)
def _xt(p, h, v):
    """horizontal distance and time of a ray with slowness p through thicknesses h with velocities v (real)."""
    q = np.sqrt(np.maximum(1.0 / v ** 2 - p ** 2, 1e-30))
    return np.sum(h * p / q), np.sum(h * q) + p * np.sum(h * p / q)


def _path_layers(model, zs, v):
    """thicknesses of the layers crossed between depth zs and the surface, and their velocities."""
    tops = np.r_[model.tops, np.inf]
    h, vv = [], []
    for i in range(model.n):
        a, b = tops[i], min(tops[i + 1], zs)
        if b > a:
            h.append(b - a)
            vv.append(v[i])
    return np.array(h), np.array(vv)


def first_arrival(model, zs, X, phase="P"):
    """first-arrival time (s), ray parameter (s/km) and kind ('direct' or 'head<i>') from a source at depth zs to
    surface receivers at epicentral distances X (array), in the model at its reference frequency.

    Own implementation (the v2 rule forbids ttmt.traveltime): the direct ray is found by Newton iteration on the
    take-off sine with a bracketing safeguard; head waves are the refractions along every interface below the
    source whose velocity exceeds every velocity above it."""
    v = np.asarray(model.vp if phase == "P" else model.vs, float)
    X = np.atleast_1d(np.asarray(X, float))
    h, va = _path_layers(model, zs, v)
    T = np.empty(len(X))
    P = np.empty(len(X))
    kind = []
    vmax = va.max()
    for j, x in enumerate(X):
        # direct: x(p) increasing on [0, 1/vmax); bisection in u = p vmax in [0, 1)
        if x < 1e-9:
            p = 0.0
        else:
            lo, hi = 0.0, 1.0 - 1e-15
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                if _xt(mid / vmax, h, va)[0] < x:
                    lo = mid
                else:
                    hi = mid
                if hi - lo < 1e-15:
                    break
            p = 0.5 * (lo + hi) / vmax
        best, bp, bk = _xt(p, h, va)[1] if x > 1e-9 else float(np.sum(h / va)), p, "direct"
        # head waves on interfaces at or below the source
        for i in range(1, model.n):
            zi = model.tops[i]
            if zi < zs - 1e-12:
                continue
            vi = v[i]
            above = v[:i]
            if vi <= above.max():
                continue
            pi_ = 1.0 / vi
            hd = np.array([max(0.0, min(model.tops[k + 1], zi) - max(model.tops[k], zs)) for k in range(i)])
            hu = np.array([model.tops[k + 1] - model.tops[k] for k in range(i)])
            hh = hd + hu
            q = np.sqrt(np.maximum(1.0 / above ** 2 - pi_ ** 2, 0.0))
            xc = np.sum(hh * pi_ / np.maximum(q, 1e-30))
            if x < xc:
                continue
            t = x * pi_ + np.sum(hh * q)
            if t < best:
                best, bp, bk = t, pi_, f"head{i}"
        T[j], P[j] = best, bp
        kind.append(bk)
    return T, P, kind
