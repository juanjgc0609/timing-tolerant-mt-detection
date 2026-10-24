"""Full-wave Green's functions of a point moment tensor in a flat layered medium (discrete wavenumber method).

Synthetic generator v2. Computes the complete elastic response (near, intermediate and far field, direct and head
waves, reflections, multiples and conversions, surface waves, free-surface effect at the receiver) of a point moment
tensor at depth zs, at receivers on the free surface, by wavenumber integration (Bouchon 1981 discrete wavenumber
sum with complex frequency) of plane-wave responses obtained with a global matrix (Chin, Hedstrom and Thigpen 1984;
Schmidt and Tango 1986): one linear system per (frequency, wavenumber) for the up- and down-going P, SV and SH
amplitudes of every layer. No ray approximation and no operator shared with the detector.

Conventions. Time dependence exp(-i omega t); horizontal transform u(x) = (2 pi)^-2 int u^(k) exp(i k.x) d2k;
coordinates x = east, y = north, z = down (the detector's frame, components E, N, -Z); moment tensor order
(Mxx, Myy, Mzz, Mxy, Mxz, Myz). The source is a stress glut sigma = C:eps - M delta(x - xs), whose plane-wave form is
a jump across z = zs of the motion-stress vector b = (u, t), t = sigma . e_z:
    [u_x] = Mxz / mu, [u_y] = Myz / mu, [u_z] = Mzz / (lambda + 2 mu),
    [t_a] = i k_b M_ab - i k_a lambda Mzz / (lambda + 2 mu)  (a, b horizontal),   [t_z] = 0.
Each horizontal wavevector k (k cos th, k sin th) is solved in its local frame (L along k, T horizontal normal), where
P-SV (u_L, u_z, t_L, t_z) and SH (u_T, t_T) decouple. The azimuthal integral is exact: the Cartesian response is a
trigonometric polynomial of degree <= 3 in th, sampled at 8 angles and integrated with
int exp(i m th) exp(i k r cos(th - phi)) dth = 2 pi i^m J_m(k r) exp(i m phi).

Units: km, s, g/cm^3, so stresses are in GPa and a unit moment is 1e18 N m (1 GPa km^3); displacement in km.
"""
import numpy as np
from scipy.special import jv

MT_ORDER = ("Mxx", "Myy", "Mzz", "Mxy", "Mxz", "Myz")
_TH = np.arange(8) * (2 * np.pi / 8)


def _sqrt_up(z):
    """square root with non-negative imaginary part (decaying / outgoing vertical wavenumber)."""
    s = np.sqrt(z)
    return np.where(s.imag < 0, -s, s)


def _psv_cols(k, w, a, b, rho):
    """motion-stress columns (..., 4, 4) of [P down, S down, P up, S up] for b = (u_L, u_z, t_L, t_z), and the vertical
    wavenumbers (gamma, nu). k, w, a, b broadcast; a, b complex velocities."""
    mu = rho * b * b
    lam = rho * a * a - 2 * mu
    g = _sqrt_up((w / a) ** 2 - k * k)
    n = _sqrt_up((w / b) ** 2 - k * k)
    cols = []
    for pl, pz, s in ((k, g, g), (n, -k, n), (k, -g, -g), (n, k, -n)):
        pl, pz, s = np.broadcast_arrays(pl * np.ones_like(g), pz * np.ones_like(g), s * np.ones_like(g))
        tl = mu * (1j * s * pl + 1j * k * pz)
        tz = lam * (1j * k * pl + 1j * s * pz) + 2 * mu * 1j * s * pz
        cols.append(np.stack([pl, pz, tl, tz], axis=-1))
    return np.stack(cols, axis=-1), g, n


def _sh_cols(k, w, b, rho):
    """columns (..., 2, 2) of [SH down, SH up] for b = (u_T, t_T), and nu."""
    mu = rho * b * b
    n = _sqrt_up((w / b) ** 2 - k * k)
    one = np.ones_like(n)
    down = np.stack([one, mu * 1j * n], axis=-1)
    up = np.stack([one, -mu * 1j * n], axis=-1)
    return np.stack([down, up], axis=-1), n


def _assemble(k, w, alpha, beta, rho, h, src_if, free_surface, nwave, rec_layer=0):
    """global matrix and the map from unknowns to the surface field, for P-SV (nwave = 4) or SH (nwave = 2).

    Layers 0..N-1 (N-1 the lower half-space); h thicknesses (N-1,). Unknown amplitudes: down waves referenced at
    the layer top, up waves at the layer bottom (every exponential is then <= 1 in modulus). Layer 0 is bounded by the
    free surface, or (free_surface=False) is an upper half-space with up-going waves only. Rows: free-surface
    tractions (if any), then the continuity of b at each interface (b_below - b_above = jump).
    Returns A (B, n, n), the interface row offset of the source, and S (B, nwave/2, n) the surface displacement."""
    N = len(rho)
    B = k.shape[0]
    half = nwave // 2
    has_down = [free_surface or i > 0 for i in range(N)]
    has_up = [i < N - 1 for i in range(N)]
    off, cnt = [], 0
    for i in range(N):
        off.append(cnt)
        cnt += half * (has_down[i] + has_up[i])
    nrow = (half if free_surface else 0) + nwave * (N - 1)
    if nrow != cnt:
        raise RuntimeError("global matrix not square")
    A = np.zeros((B, cnt, cnt), complex)
    cols, vz = [], []
    for i in range(N):
        if nwave == 4:
            c, g, n = _psv_cols(k, w, alpha[:, i], beta[:, i], rho[i])
            vz.append(np.stack([g, n], -1))
        else:
            c, n = _sh_cols(k, w, beta[:, i], rho[i])
            vz.append(n[:, None])
        cols.append(c)

    def block(i, z_rel_top):
        """(B, nwave, n_i) motion-stress at depth z_rel_top below the top of layer i, and the column slice."""
        parts = []
        if has_down[i]:
            e = np.exp(1j * vz[i] * z_rel_top) if np.isfinite(z_rel_top) else None
            parts.append(cols[i][..., :half] * e[:, None, :])
        if has_up[i]:
            hb = h[i] - z_rel_top
            parts.append(cols[i][..., half:] * np.exp(1j * vz[i] * hb)[:, None, :])
        return np.concatenate(parts, axis=-1)

    r = 0
    surf = None
    if free_surface:
        b0 = block(0, 0.0)
        A[:, :half, off[0]:off[0] + b0.shape[-1]] = b0[:, half:, :]
        surf = (0, b0)
        r = half
    else:
        # upper half-space: only up waves, referenced at its bottom (the first interface); receiver at z = 0
        surf = (0, block(0, 0.0))
    src_row = None
    for i in range(N - 1):
        above = block(i, h[i])
        below = block(i + 1, 0.0)
        A[:, r:r + nwave, off[i]:off[i] + above.shape[-1]] = -above
        A[:, r:r + nwave, off[i + 1]:off[i + 1] + below.shape[-1]] = below
        if i + 1 == src_if:
            src_row = r
        r += nwave
    S = np.zeros((B, half, cnt), complex)
    if rec_layer == 0:
        S[:, :, off[0]:off[0] + surf[1].shape[-1]] = surf[1][:, :half, :]
    else:                                   # receiver at the top of layer rec_layer (a borehole sensor)
        br = block(rec_layer, 0.0)
        S[:, :, off[rec_layer]:off[rec_layer] + br.shape[-1]] = br[:, :half, :]
    return A, src_row, S


class FKConfig:
    """numerical parameters. fs: sampling rate of the output (Hz); nt: samples (window length nt / fs from the
    origin time); fmax: highest frequency computed (<= fs / 2); f_gauss: the moment-rate function is a Gaussian
    exp(-(f / f_gauss)^2) in amplitude (an anti-alias low-pass that is exact under complex frequency); t_shift: delay
    of that Gaussian (s); eps_t: damping exp(-eps t) with eps = eps_t / T; L: period of the wavenumber sum (km);
    k_margin: extra wavenumbers beyond omega / vs_min (1/km), plus kz / zs for evanescent waves."""

    def __init__(self, fs=25.0, nt=1200, fmax=None, f_gauss=8.0, t_shift=1.0, eps_t=3.0, L=None, L_factor=3.0, k_over=1.15,
                 kz=12.0, taper=(0.8, 1.0)):
        self.fs, self.nt = float(fs), int(nt)
        self.T = self.nt / self.fs
        self.fmax = float(fmax or 0.5 * fs)
        self.f_gauss, self.t_shift, self.eps_t = float(f_gauss), float(t_shift), float(eps_t)
        self.L, self.L_factor = L, float(L_factor)
        self.k_over, self.kz, self.taper = float(k_over), float(kz), taper

    def as_dict(self):
        return {k: v for k, v in vars(self).items()}


def greens(model, zs, r, phi, cfg=None, free_surface=True, quantity="velocity", chunk=6000, zr=0.0):
    """time series (nrec, 6, 3, nt) at surface receivers at epicentral distance r (km) and azimuth phi (radians,
    counter-clockwise from x = east, i.e. phi = atan2(dy, dx)), for unit moment-tensor basis elements in MT_ORDER
    (Mxy means Mxy = Myx = 1) with moment-rate function the Gaussian of cfg (unit area, centred at cfg.t_shift).
    Components (x, y, z-down). quantity: 'displacement' (km per 1e18 N m) or 'velocity' (km/s per 1e18 N m).
    zr: receiver depth (km, 0 = free surface; > 0 a borehole sensor, which records the down-going free-surface
    reflection, the 'ghost'); requires zr < zs."""
    cfg = cfg or FKConfig()
    r = np.atleast_1d(np.asarray(r, float))
    phi = np.atleast_1d(np.asarray(phi, float))
    rec_layer = 0
    m = model
    if zr > 0:
        if zr >= zs:
            raise ValueError("receiver must be above the source")
        m, rec_layer = model.split(zr)
    m, src_if = m.split(zs)
    N = m.n
    rho = np.asarray(m.rho, float)
    h = np.asarray(m.thickness, float)
    df = 1.0 / cfg.T
    f = np.arange(0.0, cfg.fmax + 1e-9, df)
    eps = cfg.eps_t / cfg.T
    w = 2 * np.pi * f + 1j * eps
    alpha, beta = m.complex_velocities(np.maximum(f, df * 0.5))
    vmax = float(np.max(m.vp))
    vsmin = float(np.min(m.vs))
    L = cfg.L or (cfg.L_factor * (r.max() + vmax * cfg.T) + 20.0)
    dk = 2 * np.pi / L
    # source-layer elastic constants (complex, per frequency)
    a_s, b_s = alpha[:, src_if], beta[:, src_if]
    mu_s = rho[src_if] * b_s * b_s
    lam_s = rho[src_if] * a_s * a_s - 2 * mu_s
    # moment function spectrum M(w) = G(w) / (-i w), G the Gaussian moment rate
    G = np.exp(1j * w * cfg.t_shift - (w / (2 * np.pi * cfg.f_gauss)) ** 2)
    if cfg.taper:
        t0, t1 = cfg.taper[0] * cfg.fmax, cfg.taper[1] * cfg.fmax
        tp = np.where(f <= t0, 1.0, np.where(f >= t1, 0.0, np.cos(0.5 * np.pi * (f - t0) / max(t1 - t0, 1e-12)) ** 2))
        G = G * tp
    Mw = G / (-1j * w)
    U = np.zeros((len(f), len(r), 6, 3), complex)
    kmax_f = np.abs(w.real) / vsmin * cfg.k_over + cfg.kz / max(zs, 0.5) + 5 * dk
    nks = np.ceil(kmax_f / dk).astype(int)
    fi = np.repeat(np.arange(len(f)), nks)
    kk = np.concatenate([np.arange(1, n + 1) * dk for n in nks])
    c8, s8 = np.cos(_TH), np.sin(_TH)
    # unit MT basis tensors
    Mb = np.zeros((6, 3, 3))
    for j, (p, q) in enumerate(((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))):
        Mb[j, p, q] = Mb[j, q, p] = 1.0
    for a0 in range(0, len(kk), chunk):
        sl = slice(a0, a0 + chunk)
        k = kk[sl]
        ii = fi[sl]
        ww = w[ii]
        al, be = alpha[ii], beta[ii]
        Ap, rp, Sp = _assemble(k, ww, al, be, rho, h, src_if, free_surface, 4, rec_layer)
        As, rs, Ss = _assemble(k, ww, al, be, rho, h, src_if, free_surface, 2, rec_layer)
        Bn = len(k)
        rhs_p = np.zeros((Bn, Ap.shape[1], 3), complex)
        for j, comp in enumerate((0, 1, 2)):               # unit jumps in u_L, u_z, t_L
            rhs_p[:, rp + comp, j] = 1.0
        rhs_s = np.zeros((Bn, As.shape[1], 2), complex)
        rhs_s[:, rs + 0, 0] = 1.0
        rhs_s[:, rs + 1, 1] = 1.0
        Rp = Sp @ np.linalg.solve(Ap, rhs_p)                 # (B, 2 [u_L, u_z], 3 jumps)
        Rs = Ss @ np.linalg.solve(As, rhs_s)                 # (B, 1 [u_T], 2 jumps)
        mu, lam = mu_s[ii][:, None], lam_s[ii][:, None]      # (B, 1) against angles
        c, s = c8[None], s8[None]
        kc = k[:, None]
        resp = np.zeros((Bn, 6, 3, 8), complex)
        for j in range(6):
            M = Mb[j]
            juL = (c * M[0, 2] + s * M[1, 2]) / mu
            juT = (-s * M[0, 2] + c * M[1, 2]) / mu
            juz = M[2, 2] / (lam + 2 * mu) * np.ones_like(c)
            jtL = 1j * kc * (c * c * M[0, 0] + 2 * c * s * M[0, 1] + s * s * M[1, 1]) \
                - 1j * kc * lam / (lam + 2 * mu) * M[2, 2]
            jtT = 1j * kc * (c * s * (M[1, 1] - M[0, 0]) + (c * c - s * s) * M[0, 1]) * np.ones_like(mu)
            uL = Rp[:, 0, 0, None] * juL + Rp[:, 0, 1, None] * juz + Rp[:, 0, 2, None] * jtL
            uz = Rp[:, 1, 0, None] * juL + Rp[:, 1, 1, None] * juz + Rp[:, 1, 2, None] * jtL
            uT = Rs[:, 0, 0, None] * juT + Rs[:, 0, 1, None] * jtT
            resp[:, j, 0] = c * uL - s * uT
            resp[:, j, 1] = s * uL + c * uT
            resp[:, j, 2] = uz
        # harmonics a_m, m = -3..3, of resp(th) = sum_m a_m exp(i m th)
        # np.fft.fft gives X_q = sum_n x_n exp(-i q th_n) = 8 a_q, so a_m = X[m mod 8] / 8
        A_m = np.fft.fft(resp, axis=-1) / 8.0
        ms = np.arange(-3, 4)
        a = np.stack([A_m[..., mm % 8] for mm in ms], axis=-1)      # (B, 6, 3, 7)
        # Bessel terms per (k, receiver): (1/2pi) sum_k dk k a_m i^m J_m(k r) exp(i m phi)
        J = np.stack([jv(abs(mm), k[:, None] * r[None]) * (1.0 if mm >= 0 or mm % 2 == 0 else -1.0) for mm in ms],
                     axis=-1)                                 # J_{m} with J_{-m} = (-1)^m J_m: (B, nrec, 7)
        ph = (1j ** ms)[None] * np.exp(1j * ms[None] * phi[:, None])   # (nrec, 7)
        Kr = J * ph[None]                                      # (B, nrec, 7)
        contrib = np.einsum("bjcm,brm->brjc", a, Kr) * (k * dk / (2 * np.pi))[:, None, None, None]
        np.add.at(U, ii, contrib)
    U *= Mw[:, None, None, None]
    if quantity == "velocity":
        U *= (-1j * w)[:, None, None, None]
    elif quantity != "displacement":
        raise ValueError(quantity)
    # inverse transform: x(t) = e^{eps t} df N irfft(conj U)
    nt = cfg.nt
    nf_full = nt // 2 + 1
    Uf = np.zeros((nf_full,) + U.shape[1:], complex)
    nuse = min(len(f), nf_full)
    Uf[:nuse] = U[:nuse]
    x = np.fft.irfft(np.conj(Uf), n=nt, axis=0) * nt * df
    t = np.arange(nt) / cfg.fs
    x *= np.exp(eps * t)[:, None, None, None]
    return np.moveaxis(x, 0, -1)                               # (nrec, 6, 3, nt)
