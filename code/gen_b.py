"""
gen_b.py - generador B: segundo generador sintetico, independiente de gen_realistic (generador A), con la misma
interfaz de llamada que gen_realistic.generate_real.

    generate_real(STA, sources, taus, n_clutter, rng, rng_m, rng_h, clut_amp, cfg=None, meta=None, statics=None)
        -> D (ns, 3, NT), picks, fuentes, amp
    generate_site(...)   igual, devolviendo los factores de sitio (ns, 2) en el lugar de las fuentes (interfaz de
                         c8_lib.generate_site, la que usa c9_run.stream)
    install()            sustituye xc.generate_real y c8_lib.generate_site por los de B (bandera de los corredores)

Diferencias de modelado frente a A (constantes y justificacion: DEFAULT_CFG y explore_contrib/REPORT_GENB.txt):
  medio        capas planas con discontinuidades (3 capas sobre un semiespacio, capa superficial lenta) en lugar de un
               gradiente suave; tiempos, parametro de rayo y angulos de salida y llegada con un trazador de rayos de
               capas propio (rayo directo, ondas conicas en cada interfaz mas rapida y reflexion en la interfaz
               profunda; primera llegada). Vp/Vs distinto por capa.
  heterogeneidad lateral   el medio verdadero es v(x, y, z) = v_capa(z) (1 + eps(x, y)), eps suave (suma de modos de
               Fourier al azar, fija por red); el retardo se integra a primer orden (Fermat) sobre el rayo 1D, de modo
               que el error de tiempo depende del trayecto y no es un termino de estacion.
  estaticas    espesor de sedimento por estacion (mezcla roca / cuenca, sesgada y positiva): retardos P y S
               correlacionados y no gaussianos; en A son N(0, sd) independientes por fase.
  fuente       doble par mas, en una fraccion de los eventos, componentes CLVD e isotropa; pulso de Brune (omega
               cuadrado) en velocidad, con frecuencia de esquina dependiente de la magnitud y de la fase, convolucionado
               con el operador de atenuacion causal t* (Q dependiente de la frecuencia); magnitudes con ley de
               Gutenberg-Richter truncada; razon de amplitudes S/P fisica (vp/vs)^3.
  propagacion  atenuacion, onda conica debil antes del rayo directo donde corresponde, reflexion en la interfaz
               profunda con coeficiente complejo (poscritica), conversiones Sp y Ps en la base del sedimento de las
               estaciones de cuenca, desdoblamiento de la onda S en un subconjunto de estaciones, dispersion de la
               polarizacion fija por estacion (mas fuerte en un subconjunto), coda con envolvente de dispersion simple.
  ruido        coloreado (pico microsismico mas fondo rojo), nivel distinto por estacion y componente, modulacion lenta
               no estacionaria; espigas y huecos en un canal.
  clutter      rafagas locales de banda estrecha, duracion y frecuencia variables, polarizacion eliptica, amplitud de
               cola pesada, tasa heterogenea entre estaciones (en A: pulsos de Ricker de 2 Hz de amplitud fija).

cfg["mode"] = "matched": configuracion igualada (medio homogeneo VP, VS de los detectores, solo doble par, sin extras)
que reproduce el esqueleto estocastico de A con el motor de B (trazador, radiacion, sintesis en frecuencia); sirve de
validacion contra gen_realistic con cfg {"g": 0, "VS0": VS}.

    python gen_b_validate.py       (validacion; no dibuja nada)
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import hashlib

import numpy as np
from scipy.signal import butter, freqs
from scipy.special import dawsn

from task0_common import AP, AS, DT, FC, NT, SNR, T, VP, VS, ricker
from task1w_mech import CP, CS, rand_dc

EZ = np.array([0.0, 0.0, 1.0])
N2 = 2 * NT                                   # rejilla de sintesis (80 s) para que la FFT no envuelva
FREQ = np.fft.rfftfreq(N2, DT)
NF = len(FREQ)

DEFAULT_CFG = {
    "mode": "full",
    # ---- medio 1D: techo de capa (km), VP, VS (km/s), densidad (g/cm3), Q_S a 1 Hz
    "ztop": (0.0, 1.2, 7.5, 22.0),
    "vp": (3.50, 5.20, 5.80, 6.70),
    "vs": (1.90, 3.00, 3.36, 3.85),
    "rho": (2.38, 2.63, 2.70, 2.80),
    "qs": (60.0, 250.0, 500.0, 800.0),
    "qp_over_qs": 2.25,
    "q_eta": 0.3,                             # Q(f) = Q0 f^eta
    "f_ref": 1.0,                             # Hz, frecuencia de referencia de la dispersion causal
    "vel_scale": 1.0,                         # factor comun de las velocidades (calibracion geometrica, ver informe)
    # ---- heterogeneidad lateral
    "eps_sd": 0.03, "eps_modes": 24, "eps_lambda": 30.0, "eps_lambda_logsd": 0.35,
    # ---- sedimento por estacion (estaticas, sitio, conversiones)
    "basin_frac": 0.4, "dh_rock_mean": 0.15, "dh_basin_shape": 2.0, "dh_basin_scale": 0.4, "dh_max": 1.5,
    "sed_vp": 2.6, "sed_vs": 1.2, "sed_rho": 2.0, "sed_qs": 25.0,
    "static_extra_sd": (0.03, 0.05),
    "site_dh0": 0.3, "site_fix_logsd": 0.25, "site_arr_logsd": 0.2,
    "conv_range": (0.15, 0.35),
    # ---- polarizacion
    "pol_fix_deg": 6.0, "pol_rough_deg": 25.0, "pol_rough_frac": 0.3, "pol_arr_deg": 4.0,
    "split_frac": 0.3, "split_dt": (0.05, 0.20),
    # ---- fuente
    "mw_range": (3.5, 4.3), "b_value": 1.0, "mw_ref": 3.9,
    "stress_drop": 3.0e6, "stress_logsd": 0.58, "fcp_over_fcs": 1.5,
    "nondc_frac": 0.2, "nondc_range": (0.1, 0.3),
    "aa_corner": 6.0, "aa_order": 4,
    "tstar_ref": 0.04,
    # ---- ondas conicas y reflexion
    "head_gamma": 1.0, "head_lmin": 1.0,
    # ---- coda
    "coda_ratio": 0.25, "coda_basin_factor": 2.5, "coda_logsd": 0.3, "coda_p_factor": 0.5, "coda_t0": 1.0,
    "coda_tc": 8.0, "coda_len": 20.0,
    # ---- ruido
    "noise": 0.02, "micro_f": 0.2, "micro_logw": 0.35, "micro_amp": 30.0, "noise_slope": 1.5,
    "noise_sta_logsd": 0.35, "noise_basin_factor": 1.6, "noise_h_over_z": 1.3, "noise_mod_sd": 0.25,
    "noise_mod_tc": 8.0,
    "spike_rate": 1.0, "spike_amp": (10.0, 40.0), "gap_rate": 0.5, "gap_len": (0.5, 3.0),
    # ---- clutter
    "burst_ref_dist": 30.0, "burst_tau": 0.5, "burst_tau_logsd": 0.5, "burst_f": (1.0, 5.0), "burst_pareto": 3.0,
    "burst_floor": 0.7, "burst_dirichlet": 0.7,
    # ---- escala global senal / ruido (el unico ajuste permitido tras ver detectores; ver informe)
    "snr_scale": 1.0,
    # ---- modo matched (constantes de A)
    "static_sd": 0.15, "dir_deg": 10.0, "site_sd": 0.5, "fc_range": (1.4, 2.8), "phase_max_deg": 45.0,
    "coda_tau_A": 3.0, "coda_ratio_A": 0.3, "coda_band_A": (1.0, 4.0),
}


def make_cfg(cfg=None):
    c = dict(DEFAULT_CFG)
    c.update(cfg or {})
    return c


# ============================================================================================ medio y trazador
class Model:
    """capas planas homogeneas; la ultima es un semiespacio."""

    def __init__(self, ztop, vp, vs, rho=None, qs=None, qp_over_qs=2.25):
        self.ztop = np.asarray(ztop, float)
        self.nl = len(self.ztop)
        self.zbot = np.r_[self.ztop[1:], np.inf]
        self.v = {"P": np.asarray(vp, float), "S": np.asarray(vs, float)}
        self.rho = np.ones(self.nl) if rho is None else np.asarray(rho, float)
        qs = np.full(self.nl, 1e9) if qs is None else np.asarray(qs, float)
        self.q = {"S": qs, "P": qs * qp_over_qs}

    def layer(self, z):
        return int(np.clip(np.searchsorted(self.ztop, z, side="right") - 1, 0, self.nl - 1))

    def thick(self, za, zb):
        lo, hi = min(za, zb), max(za, zb)
        return np.clip(np.minimum(self.zbot, hi) - np.maximum(self.ztop, lo), 0.0, None)


def model_from_cfg(c):
    if c["mode"] == "matched":
        return Model([0.0], [VP], [VS])
    s = c["vel_scale"]
    return Model(c["ztop"], np.asarray(c["vp"]) * s, np.asarray(c["vs"]) * s, c["rho"], c["qs"], c["qp_over_qs"])


def _xt(h, v, p):
    """distancia horizontal y tiempo de un rayo de parametro p (n,) que cruza espesores h (nl,)."""
    act = h > 0
    c = np.sqrt(np.clip(1.0 - (p[:, None] * v[None, act]) ** 2, 1e-300, None))
    return (h[act] * p[:, None] * v[None, act] / c).sum(1), (h[act] / (v[None, act] * c)).sum(1)


def _solve_p(h, v, X):
    """parametro de rayo p (n,) con X(p) = X para espesores h; biseccion en el angulo de la capa mas rapida."""
    X = np.asarray(X, float)
    vmax = v[h > 0].max()
    lo = np.zeros_like(X)
    hi = np.full_like(X, np.pi / 2)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        x, _t = _xt(h, v, np.sin(mid) / vmax)
        big = x > X
        hi = np.where(big, mid, hi)
        lo = np.where(big, lo, mid)
    return np.sin(0.5 * (lo + hi)) / vmax


def _legs(m, v, p, order, h):
    """tramos (capa, dx, dt, longitud) de un rayo de parametro p (n,) por las capas `order` con espesores h."""
    out = []
    for k in order:
        if h[k] <= 0:
            continue
        c = np.sqrt(np.clip(1.0 - (p * v[k]) ** 2, 1e-300, None))
        out.append((k, h[k] * p * v[k] / c, h[k] / (v[k] * c), h[k] / c))
    return out


def trace(m, zs, zr, X, phase, kind="direct", iface=None):
    """rayo entre profundidades zs (fuente) y zr (receptor) a distancias horizontales X (n,).
    kind: direct | head (onda conica sobre el techo de la capa iface) | refl (reflexion en el techo de iface).
    Devuelve dict: T, p, valid (n,), down (salida hacia abajo), segs [(capa, dx, dt, len)] en orden de viaje,
    ks, kr (capas de fuente y receptor), run (longitud deslizada, solo head)."""
    v = m.v[phase]
    X = np.atleast_1d(np.asarray(X, float))
    ks, kr = m.layer(zs), m.layer(zr)
    n = len(X)
    if kind == "direct":
        h = m.thick(zs, zr)
        if h.sum() <= 0:
            p = np.full(n, 1.0 / v[ks])
            return {"T": X / v[ks], "p": p, "valid": np.ones(n, bool), "down": False, "ks": ks, "kr": kr,
                    "segs": [(ks, X, X / v[ks], X)], "run": np.zeros(n), "horizontal": True}
        p = _solve_p(h, v, X)
        up = zs > zr
        order = range(ks, kr - 1, -1) if up else range(ks, kr + 1)
        segs = _legs(m, v, p, order, h)
        return {"T": sum(s[2] for s in segs), "p": p, "valid": np.ones(n, bool), "down": not up, "ks": ks, "kr": kr,
                "segs": segs, "run": np.zeros(n)}
    zi = m.ztop[iface]
    if zi < max(zs, zr) - 1e-12:
        return None
    hd, hu = m.thick(zs, zi), m.thick(zr, zi)
    od, ou = range(ks, iface), range(iface - 1, kr - 1, -1)
    if kind == "head":
        if v[iface] <= v[(hd + hu) > 0].max() if ((hd + hu) > 0).any() else False:
            return None
        p = np.full(n, 1.0 / v[iface])
        segs_d, segs_u = _legs(m, v, p, od, hd), _legs(m, v, p, ou, hu)
        xc = sum(s[1] for s in segs_d + segs_u) if (segs_d or segs_u) else np.zeros(n)
        run = X - xc
        segs = segs_d + [(iface, run, run / v[iface], run)] + segs_u
        return {"T": sum(s[2] for s in segs), "p": p, "valid": run >= 0, "down": True, "ks": ks, "kr": kr,
                "segs": segs, "run": run, "xc": xc}
    h = hd + hu
    if h.sum() <= 0:
        return None
    p = _solve_p(h, v, X)
    segs = _legs(m, v, p, od, hd) + _legs(m, v, p, ou, hu)
    return {"T": sum(s[2] for s in segs), "p": p, "valid": np.ones(n, bool), "down": True, "ks": ks, "kr": kr,
            "segs": segs, "run": np.zeros(n)}


def first_arrival(m, zs, zr, X, phase):
    """tiempo de la primera llegada (directa u onda conica) y su parametro de rayo, medio 1D sin perturbar."""
    r = trace(m, zs, zr, X, phase)
    Tm, pm = r["T"].copy(), r["p"].copy()
    for i in range(1, m.nl):
        h = trace(m, zs, zr, X, phase, "head", i)
        if h is None:
            continue
        better = h["valid"] & (h["T"] < Tm)
        Tm, pm = np.where(better, h["T"], Tm), np.where(better, h["p"], pm)
    return Tm, pm


def ray_dirs(m, ray, phase, eh):
    """direcciones unitarias de salida (fuente) y de llegada (receptor), z positivo hacia abajo; eh (n, 3)."""
    v = m.v[phase]
    ss = np.clip(ray["p"] * v[ray["ks"]], 0.0, 1.0)
    sr = np.clip(ray["p"] * v[ray["kr"]], 0.0, 1.0)
    cs, cr = np.sqrt(1.0 - ss ** 2), np.sqrt(1.0 - sr ** 2)
    if ray.get("horizontal"):
        return eh.copy(), eh.copy()
    sgn_s = 1.0 if ray["down"] else -1.0
    sgn_r = -1.0 if (ray["down"] or ray["ks"] >= ray["kr"]) else 1.0
    if not ray["down"]:
        sgn_r = -1.0
    return ss[:, None] * eh + sgn_s * cs[:, None] * EZ, sr[:, None] * eh + sgn_r * cr[:, None] * EZ


# ============================================================================================ propiedades de la red
class NetProps:
    """propiedades fijas de una red (geologia y estaciones); todo sale de la semilla."""

    def __init__(self, STA, seed, cfg=None):
        c = make_cfg(cfg)
        r = np.random.default_rng(seed)
        ns = len(STA)
        self.ns = ns
        lam = c["eps_lambda"] * np.exp(r.normal(0.0, c["eps_lambda_logsd"], c["eps_modes"]))
        az = r.uniform(0, 2 * np.pi, c["eps_modes"])
        self.k = (2 * np.pi / lam)[:, None] * np.stack([np.cos(az), np.sin(az)], axis=1)
        self.ph = r.uniform(0, 2 * np.pi, c["eps_modes"])
        self.ea = c["eps_sd"] * np.sqrt(2.0 / c["eps_modes"])
        self.basin = r.random(ns) < c["basin_frac"]
        dh = np.where(self.basin, r.gamma(c["dh_basin_shape"], c["dh_basin_scale"], ns),
                      r.exponential(c["dh_rock_mean"], ns))
        self.dh = np.minimum(dh, c["dh_max"])
        vp1, vs1 = c["vp"][0] * c["vel_scale"], c["vs"][0] * c["vel_scale"]
        self.lead = self.dh * (1.0 / c["sed_vs"] - 1.0 / c["sed_vp"])          # adelanto Sp / retardo Ps (s)
        st = np.stack([self.dh * (1.0 / c["sed_vp"] - 1.0 / vp1), self.dh * (1.0 / c["sed_vs"] - 1.0 / vs1)], axis=1)
        self.statics = st + r.normal(0.0, 1.0, (ns, 2)) * np.asarray(c["static_extra_sd"])[None]
        amax = np.array([np.sqrt(c["rho"][0] * vp1 / (c["sed_rho"] * c["sed_vp"])),
                         np.sqrt(c["rho"][0] * vs1 / (c["sed_rho"] * c["sed_vs"]))])
        sat = 1.0 - np.exp(-self.dh / c["site_dh0"])
        self.site = (1.0 + (amax[None] - 1.0) * sat[:, None]) * np.exp(r.normal(0.0, c["site_fix_logsd"], (ns, 2)))
        self.tstar_sed = np.stack([self.dh / c["sed_vp"] / (c["sed_qs"] * c["qp_over_qs"]),
                                   self.dh / c["sed_vs"] / c["sed_qs"]], axis=1)
        self.conv = np.where(self.basin, r.uniform(*c["conv_range"], ns), 0.0)
        self.rough = r.random(ns) < c["pol_rough_frac"]
        sd = np.tan(np.radians(np.where(self.rough, c["pol_rough_deg"], c["pol_fix_deg"])))
        self.dn = r.normal(0.0, 1.0, (ns, 2, 3)) * sd[:, None, None]
        self.split = r.random(ns) < c["split_frac"]
        self.split_az = r.uniform(0, np.pi, ns)
        self.split_dt = np.where(self.split, r.uniform(*c["split_dt"], ns), 0.0)
        self.noise = np.exp(r.normal(0.0, c["noise_sta_logsd"], ns)) * np.where(self.basin, c["noise_basin_factor"], 1.0)
        self.burst_w = r.dirichlet(np.full(ns, c["burst_dirichlet"]))
        self.coda = c["coda_ratio"] * np.where(self.basin, c["coda_basin_factor"], 1.0) * np.exp(
            r.normal(0.0, c["coda_logsd"], ns))

    def eps(self, xy):
        """perturbacion relativa de velocidad en puntos xy (..., 2)."""
        return self.ea * np.cos(np.einsum("...j,mj->...m", xy, self.k) + self.ph).sum(-1)


_PROPS = {}


def network_props(STA, cfg=None):
    """propiedades de la red: cfg['net_seed'] (lista de enteros) o, si falta, un resumen de las coordenadas."""
    c = make_cfg(cfg)
    STA = np.asarray(STA, float)
    seed = c.get("net_seed")
    if seed is None:
        seed = [int.from_bytes(hashlib.sha256(np.round(STA, 6).tobytes()).digest()[:8], "little"), 7001]
    key = (tuple(int(x) for x in seed), STA.tobytes(), c["vel_scale"], c["eps_sd"])
    if key not in _PROPS:
        _PROPS[key] = NetProps(STA, list(key[0]), c)
    return _PROPS[key]


_GL = np.array([0.5 - np.sqrt(0.15), 0.5, 0.5 + np.sqrt(0.15)])
_GW = np.array([5.0, 8.0, 5.0]) / 18.0


def lateral_dt(props, src, STA, ray, X):
    """retardo de primer orden del rayo por la heterogeneidad lateral: - integral de eps dt sobre el trayecto."""
    n = len(X)
    dxy = STA[:, :2] - src[None, :2]
    pos = np.zeros(n)
    out = np.zeros(n)
    Xs = np.maximum(X, 1e-9)
    for _k, dx, dt, _l in ray["segs"]:
        dx = np.broadcast_to(dx, (n,))
        f = (pos[:, None] + dx[:, None] * _GL[None]) / Xs[:, None]
        pts = src[None, None, :2] + f[..., None] * dxy[:, None, :]
        out -= dt * (props.eps(pts) * _GW[None]).sum(1)
        pos = pos + dx
    return out


def ray_tstar(m, ray, phase):
    return sum(dt / m.q[phase][k] for k, _dx, dt, _l in ray["segs"])


def ray_length(ray):
    return sum(l for _k, _dx, _dt, l in ray["segs"])


# ============================================================================================ radiacion
def rad_vectors(M, ns_, nr_, eh, phase):
    """vector de desplazamiento (n, 3) del tensor M en el receptor, sin esparcimiento ni normalizacion (CP, CS): P segun
    la direccion de llegada con amplitud n_s.M.n_s; S = parte transversal de M n_s transportada conservando SV y SH."""
    if phase == "P":
        return nr_ * np.einsum("ni,ij,nj->n", ns_, M, ns_)[:, None]
    u = ns_ @ M.T
    ut = u - np.einsum("ni,ni->n", u, ns_)[:, None] * ns_
    esh = np.cross(EZ[None], eh)
    sv_s, sv_r = np.cross(esh, ns_), np.cross(esh, nr_)
    return (np.einsum("ni,ni->n", ut, sv_s)[:, None] * sv_r + np.einsum("ni,ni->n", ut, esh)[:, None] * esh)


def nondc_tensor(Mdc, f, psi, sgn_c, sgn_i):
    """sqrt(1 - f^2) DC + f (cos(psi) CLVD + sin(psi) ISO), con la norma de Frobenius del doble par (sqrt 2). El CLVD
    tiene su eje de simetria en el eje B del doble par (ortogonal a el en el producto de Frobenius)."""
    w, V = np.linalg.eigh(Mdc)
    b = V[:, int(np.argmin(np.abs(w)))]
    clvd = (3.0 * np.outer(b, b) - np.eye(3)) / np.sqrt(6.0) * np.sqrt(2.0)
    iso = np.eye(3) * np.sqrt(2.0 / 3.0)
    return np.sqrt(1.0 - f * f) * Mdc + f * (np.cos(psi) * sgn_c * clvd + np.sin(psi) * sgn_i * iso)


def refl_coef(m, phase, p, iface):
    """coeficiente de reflexion en el techo de la capa iface (complejo en regimen poscritico). S: formula SH;
    P: formula acustica (aproximacion documentada)."""
    v, rho = m.v[phase], m.rho
    v1, v2, r1, r2 = v[iface - 1], v[iface], rho[iface - 1], rho[iface]
    c1 = np.sqrt(np.clip(1.0 - (p * v1) ** 2, 0.0, None)).astype(complex)
    s2 = (p * v2) ** 2
    c2 = np.where(s2 <= 1.0, np.sqrt(np.clip(1.0 - s2, 0.0, None)), -1j * np.sqrt(np.clip(s2 - 1.0, 0.0, None)))
    if phase == "S":
        a, b = r1 * v1 * c1, r2 * v2 * c2
        return (a - b) / (a + b)
    a, b = r2 * v2 * c1, r1 * v1 * c2
    return (a - b) / (a + b)


# ============================================================================================ pulsos (frecuencia)
_AA = {}


def _aa(c):
    key = (c["aa_corner"], c["aa_order"])
    if key not in _AA:
        b, a = butter(c["aa_order"], 2 * np.pi * c["aa_corner"], analog=True)
        _AA[key] = freqs(b, a, 2 * np.pi * FREQ)[1]
    return _AA[key]


def brune_fc(mw, stress, beta_km):
    """frecuencia de esquina de Brune de la onda S (Hz): 0.49 beta (stress / M0)^(1/3), SI."""
    m0 = 10.0 ** (1.5 * mw + 9.1)
    return 0.49 * beta_km * 1e3 * (stress / m0) ** (1.0 / 3.0)


def disp_spectrum(fc, tstar, c):
    """espectro del pulso de desplazamiento de campo lejano (area unidad): Brune causal 1 / (1 + i f / fc)^2, operador
    de atenuacion causal exp(-pi f t*(f)) exp(2 i f t*(f) ln(f / f_ref)) con t*(f) = t* f^-eta, y filtro antialias."""
    f = FREQ
    fe = np.maximum(f, 0.05)
    ts = tstar * fe ** (-c["q_eta"])
    att = np.exp(-np.pi * f * ts + 2j * f * ts * np.log(fe / c["f_ref"]))
    return att * _aa(c) / (1.0 + 1j * f / fc) ** 2


_NORM = {}


def pulse_norm(c):
    """pico del pulso de velocidad de referencia (onda S de un evento mw_ref con caida de esfuerzo mediana y t*_ref)."""
    key = (c["mw_ref"], c["stress_drop"], c["tstar_ref"], c["aa_corner"], c["q_eta"], c["vs"][2], c["vel_scale"])
    if key not in _NORM:
        fc = brune_fc(c["mw_ref"], c["stress_drop"], c["vs"][2] * c["vel_scale"])
        w = np.fft.irfft(2j * np.pi * FREQ * disp_spectrum(fc, c["tstar_ref"], c), N2)
        _NORM[key] = float(np.abs(w).max())
    return _NORM[key]


def ricker_spectrum(fc, phi):
    """espectro (rejilla FREQ, muestreo DT) del Ricker de frecuencia central fc con pico 1 en t = 0, rotado en fase phi
    como gen_realistic.source_wavelet: cos(phi) r - sin(phi) H[r]."""
    f = FREQ
    R = (2.0 / (np.sqrt(np.pi) * fc)) * (f / fc) ** 2 * np.exp(-(f / fc) ** 2) / DT
    return R * np.exp(1j * phi)


def _shift(t):
    return np.exp(-2j * np.pi * FREQ * t)


# ============================================================================================ ruido
_NS = {}


def noise_shape(c):
    """raiz de la densidad espectral del ruido de fondo, normalizada para que su salida por el filtro adaptado de los
    detectores (Ricker de 2 Hz) tenga la varianza de un ruido blanco de desviacion 1."""
    key = (c["micro_f"], c["micro_logw"], c["micro_amp"], c["noise_slope"])
    if key not in _NS:
        f = np.maximum(FREQ, 0.05)
        S = c["micro_amp"] * np.exp(-np.log(f / c["micro_f"]) ** 2 / (2 * c["micro_logw"] ** 2)) + (f / FC) ** (
            -c["noise_slope"])
        S[0] = 0.0
        ker = np.zeros(N2)
        k = ricker(np.arange(-1.0, 1.0 + 1e-9, DT))
        ker[:len(k)] = k
        K2 = np.abs(np.fft.rfft(ker)) ** 2
        wgt = np.r_[1.0, np.full(NF - 2, 2.0), 1.0]
        S *= (wgt * K2).sum() / (wgt * K2 * S).sum()
        _NS[key] = np.sqrt(S)
    return _NS[key]


def _slow_modulation(r, shape, tc):
    nk = int(np.ceil(T[-1] / tc)) + 3
    g = r.normal(size=shape + (nk,))
    x = T / tc
    i = x.astype(int)
    w = x - i
    w = w * w * (3 - 2 * w)
    return g[..., i] * (1 - w) + g[..., i + 1] * w


# ============================================================================================ generador
def _child(rng):
    try:
        return rng.spawn(1)[0]
    except (AttributeError, TypeError, ValueError):
        return np.random.default_rng(int(rng.integers(2 ** 63)))


def _horiz(src, STA):
    dxy = STA[:, :2] - src[None, :2]
    X = np.hypot(dxy[:, 0], dxy[:, 1])
    eh = np.zeros((len(STA), 3))
    ok = X > 1e-9
    eh[ok, :2] = dxy[ok] / X[ok, None]
    eh[~ok, 0] = 1.0
    return X, eh


def _perturb(vec, dn):
    a0 = np.linalg.norm(vec)
    if a0 <= 0:
        return vec
    w = vec / a0 + dn
    return a0 * w / np.linalg.norm(w)


def event_arrivals(m, props, src, STA, c):
    """llegadas de un evento por fase: lista de dict por rama con T (incluye heterogeneidad lateral), direcciones,
    t*, factor geometrico y tipo; y el tiempo de la primera llegada por estacion (ns, 2)."""
    X, eh = _horiz(src, STA)
    out = {}
    first = np.zeros((len(STA), 2))
    for bi, ph in enumerate("PS"):
        v = m.v[ph]
        br = []
        r = trace(m, src[2], 0.0, X, ph)
        r["kind"] = "direct"
        br.append(r)
        for i in range(1, m.nl):
            for kind in ("head", "refl"):
                if kind == "refl" and i != m.nl - 1:
                    continue
                q = trace(m, src[2], 0.0, X, ph, kind, i)
                if q is not None:
                    q["kind"], q["iface"] = kind, i
                    br.append(q)
        for q in br:
            q["Tt"] = q["T"] + (lateral_dt(props, src, STA, q, X) if props is not None else 0.0)
            q["ns"], q["nr"] = ray_dirs(m, q, ph, eh)
            q["tstar"] = ray_tstar(m, q, ph)
            L = ray_length(q)
            if q["kind"] == "head":
                run = np.maximum(q["run"], c["head_lmin"])
                q["geom"] = c["head_gamma"] * v[q["iface"]] / (2 * np.pi) / (np.sqrt(np.maximum(X, 1.0)) * run ** 1.5)
            elif q["kind"] == "refl":
                q["geom"] = refl_coef(m, ph, q["p"], q["iface"]) / L
            else:
                q["geom"] = 1.0 / L
        cand = np.stack([np.where(q["valid"], q["Tt"], np.inf) for q in br if q["kind"] != "refl"])
        first[:, bi] = cand.min(0)
        out[ph] = br
    return out, first, X, eh


def true_times(STA, src, cfg=None, props=None, statics=None):
    """tiempos de la primera llegada P y S (ns, 2) del generador B desde src, con estaticas."""
    c = make_cfg(cfg)
    m = model_from_cfg(c)
    STA = np.asarray(STA, float)
    props = network_props(STA, c) if props is None else props
    _a, first, _X, _eh = event_arrivals(m, props, np.asarray(src, float), STA, c)
    return first + (props.statics if statics is None else np.asarray(statics, float))


def _generate_full(STA, sources, taus, n_clutter, rng, rng_m, rng_h, clut_amp, c, meta, statics):
    ns = len(STA)
    m = model_from_cfg(c)
    props = network_props(STA, c)
    st = props.statics if statics is None else np.asarray(statics, float).reshape(ns, 2)
    rx = _child(rng_h)
    Dh = np.zeros((ns, 3, NF), complex)
    Dt = np.zeros((ns, 3, N2))
    picks = [[] for _ in range(ns)]
    amp = [[] for _ in range(ns)]
    norm = pulse_norm(c)
    kappa = SNR * c["snr_scale"]
    coef_s = AS / np.sqrt(CS)
    tarr = np.tan(np.radians(c["pol_arr_deg"]))
    ncoda = int(c["coda_len"] / DT)
    ucoda = np.arange(ncoda) * DT
    env = (1.0 + ucoda / c["coda_t0"]) ** -1.5 * np.exp(-ucoda / c["coda_tc"])
    for src, tau in zip(sources, taus):
        src = np.asarray(src, float)
        Mt = rand_dc(rng_m)
        lo, hi = c["mw_range"]
        u = rx.random()
        b = c["b_value"] * np.log(10.0)
        mw = lo - np.log(1.0 - u * (1.0 - np.exp(-b * (hi - lo)))) / b
        stress = c["stress_drop"] * float(np.exp(rx.normal(0.0, c["stress_logsd"])))
        nondc = 0.0
        if rx.random() < c["nondc_frac"]:
            nondc = float(rx.uniform(*c["nondc_range"]))
            Mt = nondc_tensor(Mt, nondc, float(rx.uniform(-np.pi / 2, np.pi / 2)), 1.0 if rx.random() < 0.5 else -1.0,
                              1.0)
        ks = m.layer(src[2])
        fcs = brune_fc(mw, stress, m.v["S"][ks])
        mrel = 10.0 ** (1.5 * (mw - c["mw_ref"]))
        if meta is not None:
            meta.append({"mw": float(mw), "stress": stress, "fc": float(fcs), "nondc": nondc})
        arr, first, X, eh = event_arrivals(m, props, src, STA, c)
        dna = rx.normal(0.0, tarr, (ns, 2, 3))
        sca = np.exp(rx.normal(0.0, c["site_arr_logsd"], (ns, 2)))
        for bi, ph in enumerate("PS"):
            fc = fcs * (c["fcp_over_fcs"] if ph == "P" else 1.0)
            cph = coef_s * (1.0 if ph == "S" else (m.v["S"][ks] / m.v["P"][ks]) ** 3)
            imp = np.sqrt(m.rho[ks] * m.v[ph][ks] / (m.rho[0] * m.v[ph][0]))
            peak = np.zeros(ns)
            energy = np.zeros(ns)
            pulse_dir = [None] * ns
            for q in arr[ph]:
                vec = rad_vectors(Mt, q["ns"], q["nr"], eh, ph)
                for i in range(ns):
                    if not q["valid"][i]:
                        continue
                    a = kappa * cph * mrel * imp * props.site[i, bi] * sca[i, bi] * q["geom"][i] / norm
                    B = disp_spectrum(fc, float(q["tstar"][i]) + props.tstar_sed[i, bi], c)
                    W = B if q["kind"] == "head" else 2j * np.pi * FREQ * B
                    ta = float(q["Tt"][i]) + tau + st[i, bi]
                    v3 = _perturb(vec[i], props.dn[i, bi] + dna[i, bi])
                    S_ = a * W * _shift(ta)
                    if ph == "S" and props.split[i]:
                        fa = np.array([np.cos(props.split_az[i]), np.sin(props.split_az[i]), 0.0])
                        fa = fa - (fa @ q["nr"][i]) * q["nr"][i]
                        nf = np.linalg.norm(fa)
                        if nf > 1e-6:
                            fa /= nf
                            vf = (v3 @ fa) * fa
                            Dh[i] += np.outer(vf, S_) + np.outer(v3 - vf, S_ * _shift(props.split_dt[i]))
                        else:
                            Dh[i] += np.outer(v3, S_)
                    else:
                        Dh[i] += np.outer(v3, S_)
                    if q["kind"] == "direct":
                        w = np.fft.irfft(W, N2)
                        pk = float(np.abs(a)) * float(np.linalg.norm(v3)) * float(np.abs(w).max())
                        peak[i] = pk
                        energy[i] = (float(np.abs(a)) * float(np.linalg.norm(v3))) ** 2 * float((w * w).sum()) * DT
                        pulse_dir[i] = (W, ta)
                        if props.conv[i] > 0:
                            # conversion en la base del sedimento: Sp adelantada a la S (polarizacion P), Ps retrasada
                            # respecto de la P (polarizacion SV)
                            esh = np.cross(EZ, eh[i])
                            if ph == "S":
                                pol_p = arr["P"][0]["nr"][i]
                                sv = np.cross(esh, q["nr"][i])
                                Dh[i] += np.outer(pol_p * props.conv[i] * (v3 @ sv), S_ * _shift(-props.lead[i]))
                            else:
                                sv = np.cross(esh, arr["S"][0]["nr"][i])
                                Dh[i] += np.outer(sv * props.conv[i] * (v3 @ q["nr"][i]), S_ * _shift(props.lead[i]))
            for i in range(ns):
                picks[i].append((float(first[i, bi]) + tau + float(st[i, bi]), ph))
                amp[i].append(float(peak[i]))
                # coda: ruido con el espectro del pulso directo y envolvente de dispersion simple desde la llegada
                ratio = props.coda[i] * (c["coda_p_factor"] if ph == "P" else 1.0)
                if ratio <= 0 or energy[i] <= 0:
                    continue
                W, ta = pulse_dir[i]
                k0 = int(np.ceil(ta / DT))
                if k0 >= N2 - 2:
                    continue
                L = min(ncoda, N2 - k0)
                x = np.zeros((3, N2))
                x[:, k0:k0 + L] = rx.normal(size=(3, L)) * env[:L]
                x = np.fft.irfft(np.fft.rfft(x, axis=1) * W[None], N2, axis=1)
                e = float((x * x).sum()) * DT
                if e > 0:
                    Dt[i] += x * np.sqrt(ratio * energy[i] / e)
    D = (np.fft.irfft(Dh, N2, axis=2) + Dt)[:, :, :NT]
    # ---- clutter: rafagas locales
    nb = int(rng.poisson(n_clutter)) if n_clutter > 0 else 0
    a0 = SNR * clut_amp / c["burst_ref_dist"]
    for _ in range(nb):
        i = int(rng.choice(ns, p=props.burst_w))
        tc = float(rng.uniform(1.0, T[-1] - 2.0))
        tb = c["burst_tau"] * float(np.exp(rng.normal(0.0, c["burst_tau_logsd"])))
        fb = float(rng.uniform(*c["burst_f"]))
        q1 = rng.normal(size=3)
        q1 /= np.linalg.norm(q1)
        q2 = rng.normal(size=3)
        q2 -= (q2 @ q1) * q1
        q2 /= np.linalg.norm(q2)
        ell, phi = float(rng.uniform()), float(rng.uniform(0, 2 * np.pi))
        a = a0 * c["burst_floor"] * float(rng.random()) ** (-1.0 / c["burst_pareto"])
        g = np.exp(-0.5 * ((T - tc) / tb) ** 2)
        arg = 2 * np.pi * fb * (T - tc) + phi
        D[i] += a * g * (np.outer(q1, np.cos(arg)) + ell * np.outer(q2, np.sin(arg)))
        picks[i].append((tc, "P" if rng.random() < 0.5 else "S"))
        amp[i].append(a)
    # ---- ruido coloreado, no estacionario, heterogeneo
    sig = SNR * c["noise"]
    if sig > 0:
        w = np.fft.irfft(np.fft.rfft(rng.normal(size=(ns, 3, N2)), axis=2) * noise_shape(c), N2, axis=2)[:, :, :NT]
        comp = np.array([c["noise_h_over_z"], c["noise_h_over_z"], 1.0])
        comp = comp / np.sqrt((comp ** 2).mean())
        mod = np.exp(c["noise_mod_sd"] * _slow_modulation(rng, (ns,), c["noise_mod_tc"]))
        noise = sig * w * props.noise[:, None, None] * comp[None, :, None] * mod[:, None, :]
        D += noise
        raw_sd = noise.std(axis=2)
        for _ in range(int(rx.poisson(c["spike_rate"]))):
            i, k, j = int(rx.integers(ns)), int(rx.integers(3)), int(rx.integers(NT))
            D[i, k, j] += (1.0 if rx.random() < 0.5 else -1.0) * float(rx.uniform(*c["spike_amp"])) * raw_sd[i, k]
    for _ in range(int(rx.poisson(c["gap_rate"]))):
        i, k = int(rx.integers(ns)), int(rx.integers(3))
        j0 = int(rx.integers(NT))
        D[i, k, j0:j0 + int(float(rx.uniform(*c["gap_len"])) / DT)] = 0.0
    return D, picks, [np.asarray(s, float) for s in sources], amp, props.site


def _generate_matched(STA, sources, taus, n_clutter, rng, rng_m, rng_h, clut_amp, c, meta, statics):
    """esqueleto estocastico de gen_realistic.generate_real (mismo consumo de rng, rng_m, rng_h) con el motor de B en
    el limite homogeneo: trazador de capas con una sola capa, radiacion de B, sintesis en frecuencia."""
    import gen_realistic as GR
    ns = len(STA)
    m = model_from_cfg(c)
    Dh = np.zeros((ns, 3, NF), complex)
    D = np.zeros((ns, 3, NT))
    picks = [[] for _ in range(ns)]
    amp = [[] for _ in range(ns)]
    st = rng_h.normal(0.0, c["static_sd"], size=(ns, 2))
    if statics is not None:
        st = np.asarray(statics, float).reshape(ns, 2)
    dn = rng_h.normal(0.0, float(np.tan(np.radians(c["dir_deg"]))), size=(ns, 2, 3))
    site = np.exp(rng_h.normal(0.0, c["site_sd"], size=(ns, 2)))
    rx = _child(rng_h)
    ca = {"coda_tau": c["coda_tau_A"], "coda_ratio": c["coda_ratio_A"]}
    sos = GR._coda_filter(c["coda_band_A"])
    for src, tau in zip(sources, taus):
        src = np.asarray(src, float)
        Mt = rand_dc(rng_m)
        fc = float(rx.uniform(*c["fc_range"]))
        phi = float(np.radians(rx.uniform(-c["phase_max_deg"], c["phase_max_deg"])))
        ew = GR.wavelet_energy(fc)
        if meta is not None:
            meta.append({"fc": fc, "phase_deg": float(np.degrees(phi))})
        X, eh = _horiz(src, STA)
        rays = {}
        for ph in "PS":
            r = trace(m, src[2], 0.0, X, ph)
            r["ns"], r["nr"] = ray_dirs(m, r, ph, eh)
            r["vec"] = rad_vectors(Mt, r["ns"], r["nr"], eh, ph) / np.sqrt(CP if ph == "P" else CS)
            r["L"] = ray_length(r)
            rays[ph] = r
        W = ricker_spectrum(fc, phi)
        for i in range(ns):
            for bi, (b, A) in enumerate((("P", AP), ("S", AS))):
                r = rays[b]
                ta = float(r["T"][i]) + tau + float(st[i, bi])
                picks[i].append((ta, b))
                vec = r["vec"][i]
                a0 = np.linalg.norm(vec)
                if a0 > 0:
                    w = vec / a0 + dn[i, bi]
                    vec = a0 * site[i, bi] * w / np.linalg.norm(w)
                av = SNR * A * float(np.linalg.norm(vec)) / float(r["L"][i])
                amp[i].append(av)
                Dh[i] += SNR * A * np.outer(vec, W * _shift(ta)) / float(r["L"][i])
                GR._add_coda(D[i], ta, av, ew, ca, sos, rx)
    D += np.fft.irfft(Dh, N2, axis=2)[:, :, :NT]
    ref = np.asarray(sources[0], float) if len(sources) else np.zeros(3)
    scale = float(np.median(np.linalg.norm(STA - ref, axis=1)))
    for _ in range(n_clutter):
        i = int(rng.integers(ns))
        tc = float(rng.uniform(2, T[-1] - 4))
        q = rng.normal(size=3)
        q /= np.linalg.norm(q)
        D[i] += SNR * clut_amp * np.outer(q, ricker(T - tc)) / scale
        picks[i].append((tc, "P" if rng.random() < 0.5 else "S"))
        amp[i].append(SNR * clut_amp / scale)
    D += rng.normal(scale=SNR * c["noise"], size=D.shape)
    return D, picks, [np.asarray(s, float) for s in sources], amp, site


def _generate(STA, sources, taus, n_clutter, rng, rng_m, rng_h, clut_amp, cfg, meta, statics):
    c = make_cfg(cfg)
    STA = np.asarray(STA, float)
    f = _generate_matched if c["mode"] == "matched" else _generate_full
    return f(STA, sources, taus, n_clutter, rng, rng_m, rng_h, clut_amp, c, meta, statics)


def generate_real(STA, sources, taus, n_clutter, rng, rng_m, rng_h, clut_amp, cfg=None, meta=None, statics=None):
    """ventana (ns, 3, NT) del generador B; misma interfaz y formato de salida que gen_realistic.generate_real.
    statics (ns, 2): retardos de estacion que sustituyen a los de la red (None: los de NetProps)."""
    D, picks, srcs, amp, _site = _generate(STA, sources, taus, n_clutter, rng, rng_m, rng_h, clut_amp, cfg, meta,
                                           statics)
    return D, picks, srcs, amp


def generate_site(STA, sources, taus, n_clutter, rng, rng_m, rng_h, clut_amp, cfg=None, meta=None, statics=None,
                  **_kw):
    """interfaz de c8_lib.generate_site: devuelve D, picks, factores de sitio (ns, 2), amp."""
    D, picks, _srcs, amp, site = _generate(STA, sources, taus, n_clutter, rng, rng_m, rng_h, clut_amp, cfg, meta,
                                           statics)
    return D, picks, site, amp


def install(cfg=None):
    """los corredores existentes (xc.NetCfg / xc.catalogue_cfg / xc.window_stream / c9_run.stream) pasan a usar B.
    cfg: configuracion por defecto que se aplica cuando el corredor pasa cfg=None."""
    import sys
    here = os.path.join(os.path.dirname(os.path.abspath(__file__)), "explore_contrib")
    if here not in sys.path:
        sys.path.insert(0, here)
    import c8_lib
    import xc
    base = dict(cfg or {})

    def gr(*a, cfg=None, **k):
        return generate_real(*a, cfg={**base, **(cfg or {})}, **k)

    def gs(*a, cfg=None, **k):
        return generate_site(*a, cfg={**base, **(cfg or {})}, **k)

    xc.generate_real = gr
    c8_lib.generate_site = gs
    return gr, gs
