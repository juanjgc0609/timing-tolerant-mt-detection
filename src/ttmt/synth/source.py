"""Source mechanisms and source-time functions of the synthetic generator v2.

Moment tensors are 3x3 symmetric arrays in the (x east, y north, z down) frame; ``as_vector`` gives the order of
``ttmt.synth.fk.MT_ORDER`` (Mxx, Myy, Mzz, Mxy, Mxz, Myz). Every family is normalised to the scalar moment
M0 = ||M||_F / sqrt(2) (Silver and Jordan 1982), so a double couple of moment M0 has Frobenius norm sqrt(2) M0.

Source-time functions are spectra of unit-area moment-rate functions, evaluated at real frequencies f (Hz) in the
numpy FFT sign convention X(f) = sum x(t) exp(-2 pi i f t). None of them is a Ricker wavelet (v2 rule).
"""
import numpy as np


# ----------------------------------------------------------------------------------------------- mechanisms
def scalar_moment(mw):
    """M0 in N m (Hanks and Kanamori 1979)."""
    return 10.0 ** (1.5 * np.asarray(mw, float) + 9.1)


def dc_from_sdr(strike, dip, rake):
    """double couple of unit scalar moment from strike, dip, rake (degrees; Aki and Richards 2002 box 4.4), converted
    from their (north, east, down) frame to (east, north, down)."""
    s, d, r = np.radians([strike, dip, rake])
    mnn = -(np.sin(d) * np.cos(r) * np.sin(2 * s) + np.sin(2 * d) * np.sin(r) * np.sin(s) ** 2)
    mee = np.sin(d) * np.cos(r) * np.sin(2 * s) - np.sin(2 * d) * np.sin(r) * np.cos(s) ** 2
    mdd = np.sin(2 * d) * np.sin(r)
    mne = np.sin(d) * np.cos(r) * np.cos(2 * s) + 0.5 * np.sin(2 * d) * np.sin(r) * np.sin(2 * s)
    mnd = -(np.cos(d) * np.cos(r) * np.cos(s) + np.cos(2 * d) * np.sin(r) * np.sin(s))
    med = -(np.cos(d) * np.cos(r) * np.sin(s) - np.cos(2 * d) * np.sin(r) * np.cos(s))
    return np.array([[mee, mne, med], [mne, mnn, mnd], [med, mnd, mdd]])


def random_dc(rng):
    """double couple with a uniformly random orientation (random rotation of a reference couple), unit M0."""
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    a, b, c, d = q
    R = np.array([[a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
                  [2 * (b * c + a * d), a * a - b * b + c * c - d * d, 2 * (c * d - a * b)],
                  [2 * (b * d - a * c), 2 * (c * d + a * b), a * a - b * b - c * c + d * d]])
    M0 = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    return R @ M0 @ R.T


def clvd_like(M_dc, rng, frac):
    """the DC plus a CLVD along its null axis with weight frac in [0, 1] of the Frobenius norm (same M0)."""
    w, V = np.linalg.eigh(M_dc)
    b = V[:, int(np.argmin(np.abs(w)))]
    clvd = (3.0 * np.outer(b, b) - np.eye(3)) / np.sqrt(6.0) * np.sqrt(2.0)
    sgn = 1.0 if rng.random() < 0.5 else -1.0
    return np.sqrt(1.0 - frac ** 2) * M_dc + frac * sgn * clvd


def isotropic(sign=1.0):
    return sign * np.eye(3) * np.sqrt(2.0 / 3.0)


def lune_mt(rng, gamma=None, delta=None):
    """moment tensor uniform in orientation with lune longitude gamma in [-30, 30] and latitude delta in [-90, 90]
    degrees (Tape and Tape 2012); uniform on the lune when both are None. Unit M0."""
    if gamma is None:
        # uniform on the lune: gamma, beta with density proportional to sin^4(beta)... sampled by rejection
        while True:
            g = rng.uniform(-np.pi / 6, np.pi / 6)
            u = rng.uniform(-1, 1)
            beta = np.arccos(u)
            lam = _lune_eig(g, beta)
            if rng.random() < 1.0:          # area element of the lune in (gamma, cos beta) is uniform
                break
    else:
        g, beta = np.radians(gamma), np.radians(90.0 - delta)
        lam = _lune_eig(g, beta)
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    a, b, c, d = q
    R = np.array([[a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
                  [2 * (b * c + a * d), a * a - b * b + c * c - d * d, 2 * (c * d - a * b)],
                  [2 * (b * d - a * c), 2 * (c * d + a * b), a * a - b * b - c * c + d * d]])
    M = R @ np.diag(lam) @ R.T
    return M / (np.linalg.norm(M) / np.sqrt(2.0))


def _lune_eig(gamma, beta):
    """eigenvalues of a unit-norm moment tensor at lune longitude gamma and colatitude beta (Tape and Tape 2012 eq. 7)."""
    U = np.array([[np.sqrt(3), 0, -np.sqrt(3)], [-1, 2, -1], [np.sqrt(2), np.sqrt(2), np.sqrt(2)]]) / np.sqrt(6)
    v = np.array([np.cos(gamma) * np.sin(beta), np.sin(gamma) * np.sin(beta), np.cos(beta)])
    return U.T @ v


def as_vector(M):
    M = np.asarray(M, float)
    return np.array([M[0, 0], M[1, 1], M[2, 2], M[0, 1], M[0, 2], M[1, 2]])


def decompose(M):
    """(iso fraction, clvd fraction, dc fraction) of the Frobenius energy of a moment tensor (Vavrycuk 2001 style,
    by energy, which sums to one)."""
    M = np.asarray(M, float)
    tr = np.trace(M) / 3.0
    D = M - tr * np.eye(3)
    e_iso = 3 * tr * tr
    w = np.linalg.eigvalsh(D)
    w = w[np.argsort(np.abs(w))]
    eps = -w[0] / max(abs(w[2]), 1e-30)
    e_dev = float(np.sum(w * w))
    f_clvd = min(1.0, 2 * abs(eps)) if e_dev > 0 else 0.0
    tot = e_iso + e_dev
    return e_iso / tot, (1 - e_iso / tot) * f_clvd, (1 - e_iso / tot) * (1 - f_clvd)


def mechanism(rng, family, frac=None):
    """one unit-M0 tensor of a family: 'dc', 'clvd' (DC + CLVD with frac in [0.2, 0.6]), 'iso' (DC + ISO with frac in
    [0.2, 0.6]), 'generic' (uniform on the lune)."""
    if family == "dc":
        return random_dc(rng)
    if family == "clvd":
        return clvd_like(random_dc(rng), rng, rng.uniform(0.2, 0.6) if frac is None else frac)
    if family == "iso":
        f = rng.uniform(0.2, 0.6) if frac is None else frac
        return np.sqrt(1 - f * f) * random_dc(rng) + f * isotropic(1.0 if rng.random() < 0.5 else -1.0)
    if family == "generic":
        return lune_mt(rng)
    raise ValueError(family)


# ----------------------------------------------------------------------------------------------- source-time functions
def brune_fc(mw, stress_mpa=3.0, beta=3.5):
    """corner frequency (Hz) of Brune (1970): 0.49 beta (stress / M0)^(1/3), beta in km/s, stress in MPa."""
    return 0.4906 * beta * 1e3 * (stress_mpa * 1e6 / scalar_moment(mw)) ** (1.0 / 3.0)


def stf_spectrum(f, kind="brune", fc=5.0, duration=None):
    """unit-area moment-rate spectrum at real frequencies f (numpy sign).
    brune: omega-squared causal pulse t exp(-2 pi fc t) normalised, 1 / (1 + i f / fc)^2.
    boxcar: duration T (extended source, Haskell 1964 one-dimensional rupture), sinc with delay T / 2.
    triangle: duration T, sinc^2 with delay T / 2.
    brune_box: the Brune pulse convolved with a boxcar of duration T (corner and rupture length together)."""
    f = np.asarray(f, float)
    br = 1.0 / (1.0 + 1j * f / fc) ** 2
    if kind == "brune":
        return br
    T = float(duration)
    box = np.sinc(f * T) * np.exp(-1j * np.pi * f * T)
    if kind == "boxcar":
        return box
    if kind == "triangle":
        return np.sinc(f * T / 2) ** 2 * np.exp(-1j * np.pi * f * T)
    if kind == "brune_box":
        return br * box
    raise ValueError(kind)


def rupture_duration(mw, stress_mpa=3.0, beta=3.5):
    """duration of an extended source proportional to its size: T = 0.37 / fc (circular crack, Madariaga 1976)."""
    return 0.37 / brune_fc(mw, stress_mpa, beta) * 2.0


# ----------------------------------------------------------------------------------------------- composite sources
def subevents(rng, M, mw, kind="single", dt_range=(0.1, 0.5), dx_km=0.3):
    """list of subevents (dt s, dx (3,) km, tensor, moment fraction) of one event.
    single: one; double: two subevents dt apart in dt_range, the second 0..dx_km away, with a second mechanism
    rotated from the first (C19 'DC2'), moment split 60/40; extended: one subevent with a rupture-duration flag."""
    if kind == "single" or kind == "extended":
        return [(0.0, np.zeros(3), M, 1.0)]
    if kind == "double":
        dt = rng.uniform(*dt_range)
        d = rng.normal(size=3)
        d *= rng.uniform(0, dx_km) / np.linalg.norm(d)
        M2 = random_dc(rng) * 0.5 + M * 0.5
        M2 /= np.linalg.norm(M2) / np.sqrt(2.0)
        return [(0.0, np.zeros(3), M, 0.6), (dt, d, M2, 0.4)]
    raise ValueError(kind)
