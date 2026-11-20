"""Radiation operator along the ray of a 1-D layered model (RAY1D card; not used by the model unless an ADR adopts it).

Same layout and constants as radiation.radiation_operator: G (2 ns 3, 6), rows (phase, station, component), columns
m = (Mxx, Myy, Mzz, Mxy, Mxz, Myz), z down, stations at z = 0, spreading 1/d with d the straight source-station
distance. What changes is the direction of the ray:

- the ray parameter p of every station and phase comes from synth.media.first_arrival (the same call that predicts
  the arrival times of the sets); the take-off direction at the source has sine p*v_source along the source-station
  azimuth, upward for a direct ray and downward for a head wave; it sets the P, SV and SH radiation;
- the incidence at the surface has sine p*v_surface; P is polarised along the arriving ray, SV perpendicular to it in
  the vertical plane of the ray, SH horizontal and transverse (the same vector at source and station). The S of the
  source is projected on SV and SH at the source and carried to the station on the same two directions.

Optional (descriptive): free-surface displacement of an upgoing P and SV (Aki and Richards 2002, section 5.2), SH
doubled. In a homogeneous model with free_surface False the operator equals radiation_operator.
"""
import numpy as np

from ttmt.radiation import AP, AS, CP, CS


def _c(n):
    nx, ny, nz = n
    return np.array([nx * nx, ny * ny, nz * nz, 2 * nx * ny, 2 * nx * nz, 2 * ny * nz])


def _B(n):
    nx, ny, nz = n
    return np.array([[nx, 0, 0, ny, nz, 0], [0, ny, 0, nx, 0, nz], [0, 0, nz, 0, nx, ny]])


def layer_velocity(model, z, phase):
    """velocity of the layer that holds depth z (km, z down); the top of a layer belongs to that layer."""
    v = np.asarray(model.vp if phase == "P" else model.vs, float)
    k = int(np.searchsorted(model.tops, z, side="right") - 1)
    return float(v[max(k, 0)])


def _eta(p, v):
    return np.real(np.sqrt(complex(1.0 / v ** 2 - p * p)))


def _fs_p(p, a, b, h):
    """surface displacement of a unit upgoing P: radial along h, vertical up = -z."""
    ea, eb = _eta(p, a), _eta(p, b)
    D = (1 / b ** 2 - 2 * p * p) ** 2 + 4 * p * p * ea * eb
    ux = 4 * a * p * ea * eb / (b * b * D)
    uz = 2 * a * ea * (1 / b ** 2 - 2 * p * p) / (b * b * D)
    return ux * h + uz * np.array([0.0, 0.0, -1.0])


def _fs_sv(p, a, b, h, svr):
    """surface displacement of a unit upgoing SV with the polarisation svr, signs taken from svr."""
    ea, eb = _eta(p, a), _eta(p, b)
    D = (1 / b ** 2 - 2 * p * p) ** 2 + 4 * p * p * ea * eb
    ux = 2 * b * eb * (1 / b ** 2 - 2 * p * p) / (b * b * D)
    uz = 4 * b * p * ea * eb / (b * b * D)
    sx = np.sign(svr @ h) or 1.0
    sz = np.sign(-svr[2]) or 1.0
    return ux * sx * h + uz * sz * np.array([0.0, 0.0, -1.0])


def ray_directions(src, STA, model):
    """per phase: (p, take-off unit vectors n_s (ns, 3), arrival unit vectors n_r (ns, 3), horizontal unit h (ns, 3))."""
    from ttmt.synth.media import first_arrival
    src = np.asarray(src, float); STA = np.asarray(STA, float)
    zs = max(float(src[2]), 0.5)                       # the same floor as the sets' arrival predictions
    dv = STA - src
    X = np.hypot(dv[:, 0], dv[:, 1])
    h = np.zeros((len(STA), 3))
    far = X > 1e-9
    h[far, :2] = dv[far, :2] / X[far, None]
    h[~far, 0] = 1.0
    out = {}
    for ph in ("P", "S"):
        _T, p, kind = first_arrival(model, zs, X, ph)
        p = np.asarray(p, float)
        si = np.minimum(p * layer_velocity(model, zs, ph), 1.0); ci = np.sqrt(1.0 - si * si)
        down = np.array([k != "direct" for k in kind])
        ns_ = si[:, None] * h; ns_[:, 2] = np.where(down, ci, -ci)
        sj = np.minimum(p * layer_velocity(model, 0.0, ph), 1.0); cj = np.sqrt(1.0 - sj * sj)
        nr = sj[:, None] * h; nr[:, 2] = -cj
        out[ph] = (p, ns_, nr, h)
    return out


def clipped(src, model, STA):
    """number of (phase, station) rays with p*v > 1 at the source or at the surface (sine clipped to 1)."""
    from ttmt.synth.media import first_arrival
    src = np.asarray(src, float); STA = np.asarray(STA, float)
    zs = max(float(src[2]), 0.5)
    X = np.hypot(STA[:, 0] - src[0], STA[:, 1] - src[1])
    k = 0
    for ph in ("P", "S"):
        p = np.asarray(first_arrival(model, zs, X, ph)[1], float)
        k += int(np.sum((p * layer_velocity(model, zs, ph) > 1.0) | (p * layer_velocity(model, 0.0, ph) > 1.0)))
    return k


def ray_operator(src, STA, model, free_surface=False):
    """G_ray (2 ns 3, 6) of a source at src (km, z down) for stations STA (ns, 3), in the 1-D model."""
    src = np.asarray(src, float); STA = np.asarray(STA, float)
    d = np.linalg.norm(STA - src, axis=1)
    dirs = ray_directions(src, STA, model)
    a0, b0 = layer_velocity(model, 0.0, "P"), layer_velocity(model, 0.0, "S")
    blocks = []
    for ph in ("P", "S"):
        p, NS, NR, H = dirs[ph]
        for i in range(len(STA)):
            ns_, nr, h = NS[i], NR[i], H[i]
            sh = np.array([-h[1], h[0], 0.0])
            if ph == "P":
                pol = _fs_p(p[i], a0, b0, h) if free_surface else nr
                g = np.outer(pol, _c(ns_)) * AP / (d[i] * np.sqrt(CP))
            else:
                Bm = _B(ns_)
                svs, svr = np.cross(sh, ns_), np.cross(sh, nr)
                if free_surface:
                    g = np.outer(_fs_sv(p[i], a0, b0, h, svr), svs @ Bm) + 2.0 * np.outer(sh, sh @ Bm)
                else:
                    g = np.outer(svr, svs @ Bm) + np.outer(sh, sh @ Bm)
                g = g * AS / (d[i] * np.sqrt(CS))
            blocks.append(g)
    return np.vstack(blocks)
