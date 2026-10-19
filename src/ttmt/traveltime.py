"""Travel times: homogeneous straight rays and first arrivals in flat layered media (c18_tt, c18j_layered).

Layered model: a list of layers (top depth km, VP km/s), the last a half-space; VS = VP / r. The first arrival is the
minimum of the direct up-going ray (ray parameter by bisection) and the head waves on the interfaces below the source
whose velocity exceeds all those above. A linear gradient is represented by thin layers (0.25 km). Pure functions,
identical to ``code/explore_contrib/c18_tt.py``.
"""
import numpy as np

HK = [(0.0, 5.5), (5.5, 6.3), (16.0, 6.7), (32.0, 7.8)]      # Hadley and Kanamori (1977), southern California
# JMA2001 P velocities, 28-point approximation of Y. Maeda (Nagoya University opentools documentation,
# create_1DvelocityStructure_JMA2001_Vp), first 13 points (0 to 65 km), linear between points
JMA2001_POINTS = [(0.0, 4.800), (1.5, 5.220), (2.5, 5.450), (4.0, 5.720), (5.0, 5.820), (15.0, 6.220), (22.0, 6.544),
                  (34.0, 7.188), (41.0, 7.494), (44.0, 7.598), (48.0, 7.688), (52.5, 7.750), (65.0, 7.845)]


def points_model(points, zmax=60.0, dz=0.25):
    """thin-layer representation of a piecewise-linear velocity profile (velocity at the layer midpoint)."""
    zp, vp = np.array([p[0] for p in points]), np.array([p[1] for p in points])
    z = np.arange(0.0, zmax, dz)
    return [(float(zz), float(np.interp(zz + dz / 2, zp, vp))) for zz in z]


JMA2001 = points_model(JMA2001_POINTS)


def gradient_model(v0, g, zmax=40.0, dz=0.25, vmax=7.8):
    z = np.arange(0.0, zmax, dz)
    return [(float(zz), float(min(v0 + g * (zz + dz / 2), vmax))) for zz in z]


def _split(model, zs):
    """layers between the surface and the source: thicknesses and velocities; and the layers below the source."""
    tops = np.array([m[0] for m in model] + [np.inf])
    v = np.array([m[1] for m in model])
    above_h, above_v = [], []
    for i in range(len(model)):
        a, b = tops[i], tops[i + 1]
        if a >= zs:
            break
        above_h.append(min(b, zs) - a)
        above_v.append(v[i])
    return np.array(above_h), np.array(above_v), tops, v


def first_arrival(model, zs, X, r=1.0):
    """first-arrival time (s) of the phase with velocities model / r, source depth zs, epicentral distances X (array)."""
    X = np.atleast_1d(np.asarray(X, float))
    h, va, tops, vall = _split(model, zs)
    va = va / r
    vall = vall / r
    out = np.full(X.shape, np.inf)
    # direct wave
    if len(h):
        vmax = va.max()
        for k, x in enumerate(X):
            if x < 1e-6:
                out[k] = float(np.sum(h / va))
                continue
            lo, hi = 0.0, (1.0 - 1e-12) / vmax
            for _ in range(80):
                p = 0.5 * (lo + hi)
                q = np.sqrt(np.maximum(1.0 - (p * va) ** 2, 1e-30))
                xp = np.sum(h * p * va / q)
                if xp < x:
                    lo = p
                else:
                    hi = p
            p = 0.5 * (lo + hi)
            q = np.sqrt(np.maximum(1.0 - (p * va) ** 2, 1e-30))
            out[k] = float(np.sum(h / (va * q)))
    # head waves on interfaces below the source
    for j in range(len(vall)):
        ztop = tops[j]
        if ztop <= zs or j == 0:
            continue
        vn = vall[j]
        vabove = vall[:j]
        if vn <= vabove.max():
            continue
        p = 1.0 / vn
        th = np.diff(np.r_[tops[:j], ztop])                  # thicknesses of layers 0..j-1
        # path: from the source down to ztop, then up from ztop to the surface
        down = np.clip(np.minimum(tops[1:j + 1], ztop) - np.maximum(tops[:j], zs), 0.0, None)
        eta = np.sqrt(np.maximum(1.0 / vabove ** 2 - p * p, 0.0))
        tint = np.sum((th + down) * eta)
        xc = np.sum((th + down) * p * vabove / np.sqrt(np.maximum(1.0 - (p * vabove) ** 2, 1e-30)))
        t = X * p + tint
        out = np.where(X >= xc, np.minimum(out, t), out)
    return out


def table(model, r, depths, xmax=120.0, dx=0.25):
    """travel-time tables: dict depth -> (x grid, tP, tS)."""
    xs = np.arange(0.0, xmax + dx, dx)
    return {float(z): (xs, first_arrival(model, z, xs, 1.0), first_arrival(model, z, xs, r)) for z in depths}


def times(model, r, src, STA):
    """(2, ns) P and S first-arrival times from src (x, y, z km) to stations at the surface."""
    STA = np.asarray(STA, float)
    X = np.hypot(STA[:, 0] - src[0], STA[:, 1] - src[1])
    return np.stack([first_arrival(model, float(src[2]), X, 1.0), first_arrival(model, float(src[2]), X, r)])


R_VPVS = 1.73                               # c18j_layered.R_VPVS


def tt_table(model, depths, r=R_VPVS, xmax=120.0, dx=0.25):
    """c18j_layered.tt_table with the model explicit (legacy: module constant c18j_layered.MODEL)."""
    return table(model, r, depths, xmax=xmax, dx=dx)


def tt_nodes(XN, STA, tab):
    """(N, 2, ns) layered first-arrival times, linear interpolation in epicentral distance at the node depth
    (c18j_layered.tt_nodes)."""
    XN, STA = np.asarray(XN, float), np.asarray(STA, float)
    out = np.zeros((len(XN), 2, len(STA)))
    for z, (xs, tp, ts) in tab.items():
        m = np.abs(XN[:, 2] - z) < 1e-6
        if not m.any():
            continue
        X = np.hypot(XN[m, None, 0] - STA[None, :, 0], XN[m, None, 1] - STA[None, :, 1])
        out[m, 0] = np.interp(X, xs, tp)
        out[m, 1] = np.interp(X, xs, ts)
    assert np.all(out > 0), "node depth not in the table"
    return out


class SourceTimes:
    """hyb_dev2.times: P and S times from a source to the stations, from a 0.25 km table at the source depth rounded to
    0.1 km (one table per rounded depth, cached per instance; legacy: a module-level cache)."""

    def __init__(self, model, r=R_VPVS):
        self.model, self.r, self.cache = model, r, {}

    def __call__(self, src, STA):
        key = round(float(src[2]), 1)
        if key not in self.cache:
            self.cache[key] = table(self.model, self.r, [key], xmax=120.0, dx=0.25)[key]
        xs, tp, ts = self.cache[key]
        X = np.hypot(STA[:, 0] - src[0], STA[:, 1] - src[1])
        return np.stack([np.interp(X, xs, tp), np.interp(X, xs, ts)])


def homogeneous(XN, STA, vp, vs):
    """(N, 2, ns) straight-ray times and the distances (N, ns) (c10_lib.Eng: d from c9_lib.rays_hom)."""
    from .radiation import rays
    _n, d = rays(XN, STA)
    return np.stack([d / vp, d / vs], axis=1), d
