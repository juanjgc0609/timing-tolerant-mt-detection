"""
gen_realistic.py - generador de formas de onda mas realista que generate_y (Tarea 1y), con la misma interfaz.

generate_y supone medio homogeneo (tiempos d/VP y d/VS, rayos rectos), el mismo pulso de Ricker de 2 Hz que
el filtro adaptado, ruido blanco y ninguna energia despues del arribo directo. Aqui:
  (a) medio 1D con gradiente lineal v(z) = v0 + g z para P y S (VS0 = 2.8 km/s, g = 0.04 1/s, VP = 1.73 VS,
      z profundidad positiva hacia abajo como en task0_common). Tiempo analitico
          t = (1/g) arccosh(1 + g^2 r^2 / (2 v(z_s) v(z_r))),  r = distancia fuente-estacion,
      evaluado como (2/g) asinh(g r / (2 sqrt(v_s v_r))), estable cuando g -> 0. Los rayos son arcos de
      circulo con centro a la profundidad -v0/g (donde v se anula); la P se polariza segun la tangente del
      rayo en la estacion (angulo de incidencia del rayo curvo). La S es la parte transversal de M n_s con
      la direccion de salida n_s y se transporta a la estacion conservando sus componentes SV y SH: rotacion
      alrededor de la normal al plano vertical del rayo, que en un medio 1D es el transporte paralelo exacto
      del marco de polarizacion (el rayo no sale de su plano). Divergencia geometrica 1/r como en generate_y.
      P y S comparten el cociente g/v0, asi que sus rayos son el mismo circulo.
  (b) retardos por estacion y fase N(0, static_sd), desviacion de polarizacion de ~10 grados y factor de sitio
      lognormal sigma 0.5, con el mismo consumo de rng_h que generate_y.
  (c) funcion fuente por evento: Ricker con frecuencia central U[1.4, 2.8] Hz (el filtro adaptado de los
      detectores sigue en 2 Hz) rotado en fase U[-45, 45] grados, w = cos(phi) r - sin(phi) H[r], con la
      transformada de Hilbert del Ricker en forma cerrada: H[r](t) = (2 s + (2 - 4 s^2) F(s)) / sqrt(pi),
      s = pi fc t, F la funcion de Dawson. La rotacion conserva la energia del pulso.
  (d) coda dispersada despues de cada arribo P y S: ruido gaussiano de tres componentes con envolvente
      exp(-(t - t_a)/coda_tau) desde el arribo, filtro Butterworth causal de orden 4 en 1-4 Hz, escalado a
      coda_ratio veces la energia del arribo directo (energia de la coda sin truncar por el fin de la ventana).
  (e) clutter identico al de generate_y (amplitud clut_amp, estacion, tiempo y polarizacion al azar, Ricker de
      2 Hz, sin coda) y, si cfg['spikes'] > 0, espigas de una muestra en una componente al azar, de amplitud
      spike_amp veces la desviacion del ruido. Las espigas no son arribos y no entran en picks ni en amp.
Semillas: el mecanismo consume rng_m y el clutter y el ruido consumen rng exactamente como generate_y, de modo
que con las mismas semillas las ventanas comparten fuentes, mecanismos, retardos, clutter y ruido; la
aleatoriedad nueva (funcion fuente, coda, espigas) sale de un hijo de rng_h (Generator.spawn), que no altera
la secuencia de rng_h. Con cfg = {"g": 0, "VS0": 3.0} los tiempos y amplitudes de los picks coinciden con los
de generate_y (autoprueba).

picks tiene el formato de generate_y: por estacion, los 2*nsrc arribos verdaderos (P y S de cada fuente, en
orden) seguidos de los picks de clutter; amp es la amplitud por pick en unidades de SNR*A/r, de modo que
degrade_y (task1y_clutter) se aplica sin cambios.

    python gen_realistic.py        (autoprueba; no dibuja nada)
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
from scipy.signal import butter, sosfilt
from scipy.special import dawsn

from task0_common import AP, AS, DT, NT, SNR, T, ricker
from task1w_mech import CP, CS, rand_dc

DEFAULT_CFG = {
    "VS0": 2.8, "g": 0.04, "vpvs": 1.73,            # medio: v(z) = v0 + g z, VP = vpvs * VS
    "static_sd": 0.15, "dir_deg": 10.0, "site_sd": 0.5,
    "fc_range": (1.4, 2.8), "phase_max_deg": 45.0,   # funcion fuente por evento
    "coda_tau": 3.0, "coda_ratio": 0.3, "coda_band": (1.0, 4.0),
    "spikes": 0, "spike_amp": 10.0,                  # espigas por ventana y su amplitud en desviaciones del ruido
    "noise": 0.02,                                   # desviacion del ruido en unidades de SNR (la de generate_y)
}
EZ = np.array([0.0, 0.0, 1.0])
_SOS = {}


def make_cfg(cfg=None):
    c = dict(DEFAULT_CFG)
    c.update(cfg or {})
    return c


def phase_velocity(phase, cfg=None):
    """(v0, g) de la fase: la P es vpvs veces la S a toda profundidad."""
    c = make_cfg(cfg)
    k = c["vpvs"] if phase == "P" else 1.0
    return k * c["VS0"], k * c["g"]


def travel_time_gradient(src, STA, phase, cfg=None):
    """tiempo de viaje de la fase en v(z) = v0 + g z desde src (3,) a STA (ns, 3) o a una estacion (3,).
    (1/g) arccosh(1 + g^2 r^2 / (2 v_s v_r)) = (2/g) asinh(x), x = g r / (2 sqrt(v_s v_r)); se evalua como
    r / sqrt(v_s v_r) * asinh(x)/x, que tiende a r/v0 sin cancelacion cuando g -> 0."""
    v0, g = phase_velocity(phase, cfg)
    src = np.asarray(src, float)
    R = np.atleast_2d(np.asarray(STA, float))
    r = np.linalg.norm(R - src, axis=1)
    q = np.sqrt((v0 + g * src[2]) * (v0 + g * R[:, 2]))
    x = g * r / (2.0 * q)
    xs = np.where(x > 1e-6, x, 1.0)
    f = np.where(x > 1e-6, np.arcsinh(xs) / xs, 1.0 - x * x / 6.0)
    t = r / q * f
    return t if np.ndim(STA) == 2 else float(t[0])


def _ray(src, sta, v0, g):
    """tangentes unitarias del rayo circular en la fuente y en la estacion, en el sentido de propagacion, y el
    versor horizontal fuente -> estacion (None si el rayo es vertical). En el plano vertical del rayo, con h la
    distancia horizontal desde la fuente, el centro esta en (h_c, -v0/g) y
        g h_c = (g X^2 + (z_r - z_s)(v_s + v_r)) / (2 X),   g R_c = hypot(g h_c, v_s),
        n_s = (v_s e_h + g h_c e_z) / (g R_c),   n_r = (v_r e_h - (g X - g h_c) e_z) / (g R_c).
    Todo va multiplicado por g, de modo que g = 0 da exactamente el rayo recto."""
    src, sta = np.asarray(src, float), np.asarray(sta, float)
    dxy = sta[:2] - src[:2]
    X = float(np.hypot(dxy[0], dxy[1]))
    dz = float(sta[2] - src[2])
    if X < 1e-9:
        n = np.array([0.0, 0.0, 1.0 if dz >= 0 else -1.0])
        return n, n.copy(), None
    eh = np.array([dxy[0] / X, dxy[1] / X, 0.0])
    vs, vr = v0 + g * src[2], v0 + g * sta[2]
    gxc = (g * X * X + dz * (vs + vr)) / (2.0 * X)
    gR = float(np.hypot(gxc, vs))
    return (vs * eh + gxc * EZ) / gR, (vr * eh - (g * X - gxc) * EZ) / gR, eh


def ray_directions(src, sta, phase="P", cfg=None):
    """(n_s, n_r): direccion de salida en la fuente y de llegada a la estacion (vectores unitarios)."""
    v0, g = phase_velocity(phase, cfg)
    ns_, nr_, _eh = _ray(src, sta, v0, g)
    return ns_, nr_


def arrival_direction(src, sta, phase="P", cfg=None):
    """direccion de propagacion del rayo al llegar a la estacion; es la polarizacion P."""
    return ray_directions(src, sta, phase, cfg)[1]


def radiation_gradient(M, src, sta, phase, cfg=None):
    """vector de desplazamiento del doble par M en la estacion, sin 1/r ni amplitud de fase, con la misma
    normalizacion (CP, CS) que generate_y. P: n_r (n_s.M.n_s). S: la parte transversal de M n_s descompuesta
    en SV y SH en la fuente y recompuesta con los versores SV y SH de la llegada."""
    v0, g = phase_velocity(phase, cfg)
    ns_, nr_, eh = _ray(src, sta, v0, g)
    if phase == "P":
        return nr_ * float(ns_ @ M @ ns_) / np.sqrt(CP)
    u = M @ ns_
    ut = u - float(u @ ns_) * ns_
    if eh is None:
        return ut / np.sqrt(CS)
    esh = np.cross(EZ, eh)
    return (float(ut @ np.cross(esh, ns_)) * np.cross(esh, nr_) + float(ut @ esh) * esh) / np.sqrt(CS)


def source_wavelet(t, fc, phi=0.0):
    """Ricker de frecuencia central fc (pico 1 en t = 0) rotado en fase phi (radianes)."""
    s = np.pi * fc * np.asarray(t, float)
    s2 = s * s
    r = (1.0 - 2.0 * s2) * np.exp(-s2)
    if phi == 0.0:
        return r
    h = (2.0 * s + (2.0 - 4.0 * s2) * dawsn(s)) / np.sqrt(np.pi)
    return np.cos(phi) * r - np.sin(phi) * h


def wavelet_energy(fc):
    """integral de w^2: 3/4 sqrt(pi/2) / (pi fc), independiente de la rotacion de fase."""
    return 0.75 * np.sqrt(np.pi / 2.0) / (np.pi * fc)


def _coda_filter(band):
    key = tuple(float(b) for b in band)
    if key not in _SOS:
        _SOS[key] = butter(4, key, btype="bandpass", fs=1.0 / DT, output="sos")
    return _SOS[key]


def _child(rng):
    try:
        return rng.spawn(1)[0]
    except (AttributeError, TypeError, ValueError):
        return np.random.default_rng(int(rng.integers(2 ** 63)))


def _add_coda(Di, ta, a, ew, c, sos, rx):
    """coda de un arribo de amplitud vectorial a (norma) y energia de pulso ew, sumada en Di (3, NT)."""
    L = int(np.ceil((6.0 * c["coda_tau"] + 1.0) / DT))
    x = sosfilt(sos, rx.normal(size=(3, L)) * np.exp(-np.arange(L) * DT / c["coda_tau"]), axis=1)
    e = float(np.sum(x ** 2) * DT)
    if c["coda_ratio"] <= 0 or e <= 0 or a <= 0:
        return
    x *= np.sqrt(c["coda_ratio"] * a * a * ew / e)
    k0 = int(np.ceil(ta / DT))
    lo, hi = max(k0, 0), min(k0 + L, NT)
    if hi > lo:
        Di[:, lo:hi] += x[:, lo - k0:hi - k0]


def generate_real(STA, sources, taus, n_clutter, rng, rng_m, rng_h, clut_amp, cfg=None, meta=None, statics=None):
    """ventana de forma de onda (ns, 3, NT) con los efectos (a)-(e); devuelve D, picks, fuentes, amp.
    meta, si es una lista, recibe por evento {'fc', 'phase_deg'}. statics (ns, 2), en s, columnas P y S:
    retardos de estacion fijos (por red) que sustituyen al sorteo por ventana; el sorteo de rng_h se hace
    igual, de modo que la secuencia de rng_h no cambia. statics=None reproduce el sorteo por ventana."""
    c = make_cfg(cfg)
    ns = len(STA)
    D = np.zeros((ns, 3, NT))
    picks = [[] for _ in range(ns)]
    amp = [[] for _ in range(ns)]
    st = rng_h.normal(0.0, c["static_sd"], size=(ns, 2))
    if statics is not None:
        st = np.asarray(statics, float).reshape(ns, 2)
    dn = rng_h.normal(0.0, float(np.tan(np.radians(c["dir_deg"]))), size=(ns, 2, 3))
    site = np.exp(rng_h.normal(0.0, c["site_sd"], size=(ns, 2)))
    rx = _child(rng_h)
    sos = _coda_filter(c["coda_band"])
    for src, tau in zip(sources, taus):
        src = np.asarray(src, float)
        d = np.linalg.norm(STA - src, axis=1)
        M = rand_dc(rng_m)
        fc = float(rx.uniform(*c["fc_range"]))
        phi = float(np.radians(rx.uniform(-c["phase_max_deg"], c["phase_max_deg"])))
        ew = wavelet_energy(fc)
        if meta is not None:
            meta.append({"fc": fc, "phase_deg": float(np.degrees(phi))})
        tt = {b: travel_time_gradient(src, STA, b, c) for b in ("P", "S")}
        for i in range(ns):
            for bi, (b, A) in enumerate((("P", AP), ("S", AS))):
                ta = float(tt[b][i]) + tau + float(st[i, bi])
                picks[i].append((ta, b))
                vec = radiation_gradient(M, src, STA[i], b, c)
                a0 = np.linalg.norm(vec)
                if a0 > 0:
                    w = vec / a0 + dn[i, bi]
                    vec = a0 * site[i, bi] * w / np.linalg.norm(w)
                av = SNR * A * float(np.linalg.norm(vec)) / d[i]
                amp[i].append(av)
                D[i] += SNR * A * np.outer(vec, source_wavelet(T - ta, fc, phi)) / d[i]
                _add_coda(D[i], ta, av, ew, c, sos, rx)
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
    for _ in range(int(c["spikes"])):
        i, k, j = int(rx.integers(ns)), int(rx.integers(3)), int(rx.integers(NT))
        D[i, k, j] += (1.0 if rx.random() < 0.5 else -1.0) * c["spike_amp"] * SNR * c["noise"]
    return D, picks, [np.asarray(s, float) for s in sources], amp


def shoot_ray(src, sta, phase="P", cfg=None):
    """trazado numerico de rayos (control independiente de las formulas): ecuaciones cinematicas
    dx/dt = v^2 p, dp/dt = -grad(v)/v integradas en 3D con DOP853 y tiro sobre el angulo de salida en el
    plano vertical fuente-estacion. Devuelve (tiempo, direccion de salida, direccion de llegada)."""
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq
    v0, g = phase_velocity(phase, cfg)
    src, sta = np.asarray(src, float), np.asarray(sta, float)
    dxy = sta[:2] - src[:2]
    X = float(np.hypot(dxy[0], dxy[1]))
    eh = np.array([dxy[0] / X, dxy[1] / X, 0.0])

    def rhs(_t, y):
        v = v0 + g * y[2]
        return np.concatenate([v * v * y[3:], [0.0, 0.0, -g / v]])

    def reach(_t, y):
        return float((y[:2] - src[:2]) @ eh[:2]) - X
    reach.terminal, reach.direction = True, 1

    def top(_t, y):
        return v0 + g * y[2] - 0.05 * v0
    top.terminal = True

    def run(th):
        n0 = np.cos(th) * eh + np.sin(th) * EZ
        s = solve_ivp(rhs, (0.0, 300.0), np.concatenate([src, n0 / (v0 + g * src[2])]), method="DOP853",
                      events=(reach, top), rtol=1e-12, atol=1e-12)
        if len(s.t_events[0]) == 0:
            return -1e3, None, None, n0
        y = s.y_events[0][0]
        return float(y[2] - sta[2]), float(s.t_events[0][0]), y, n0

    th = brentq(lambda a: run(a)[0], -1.3, 1.3, xtol=1e-14)
    _m, t, y, n0 = run(th)
    return t, n0, (v0 + g * y[2]) * y[3:]


def _selftest():
    from task0_common import VP, VS
    pairs = [(np.array([0.0, 0.0, 5.0]), np.array([3.0, 4.0, 0.0])),
             (np.array([2.0, -3.0, 13.0]), np.array([40.0, 10.0, 0.0])),
             (np.array([-5.0, 5.0, 9.0]), np.array([-30.0, -35.0, 0.0]))]
    ok = True

    def check(name, val, tol):
        nonlocal ok
        ok &= bool(val < tol)
        print(f"   {name:<62} {val:.3e}  (< {tol:g})  {'ok' if val < tol else 'FALLA'}")

    # 1. limite g -> 0: tiempo r/v y rayo recto
    et = ed = 0.0
    for g in (1e-9, 1e-7, 0.0):
        cg = make_cfg({"g": g})
        for src, sta in pairs:
            r = float(np.linalg.norm(sta - src))
            for ph in "PS":
                v0 = phase_velocity(ph, cg)[0]
                et = max(et, abs(travel_time_gradient(src, sta[None], ph, cg)[0] - r / v0) / (r / v0))
                for n in ray_directions(src, sta, ph, cg):
                    ed = max(ed, 1.0 - float(n @ (sta - src)) / r)
    check("g->0: error relativo del tiempo frente a r/v0", et, 1e-6)
    check("g->0: 1 - cos(direccion, rayo recto) en salida y llegada", ed, 1e-9)
    # 2. forma estable contra la formula arccosh; vectores unitarios con el gradiente por defecto
    c = make_cfg()
    ea = eu = 0.0
    for src, sta in pairs:
        r = float(np.linalg.norm(sta - src))
        for ph in "PS":
            v0, g = phase_velocity(ph, c)
            ref = np.arccosh(1 + g * g * r * r / (2 * (v0 + g * src[2]) * (v0 + g * sta[2]))) / g
            ea = max(ea, abs(travel_time_gradient(src, sta, ph, c) - ref) / ref)
            eu = max(eu, max(abs(float(np.linalg.norm(n)) - 1.0) for n in ray_directions(src, sta, ph, c)))
    check("forma asinh frente a (1/g) arccosh(...), error relativo", ea, 1e-12)
    check("| |n| - 1 | de las direcciones de salida y llegada", eu, 1e-12)
    # 3. trazado numerico de rayos: tiempo y direcciones
    print("   trazado numerico (DOP853 + tiro), gradiente por defecto:")
    wt = wa = 0.0
    for src, sta in pairs:
        for ph in "PS":
            tn, ns_n, nr_n = shoot_ray(src, sta, ph, c)
            ta_ = travel_time_gradient(src, sta, ph, c)
            ns_a, nr_a = ray_directions(src, sta, ph, c)
            ang = max(np.degrees(np.arccos(np.clip(ns_a @ ns_n / np.linalg.norm(ns_n), -1, 1))),
                      np.degrees(np.arccos(np.clip(nr_a @ nr_n / np.linalg.norm(nr_n), -1, 1))))
            wt, wa = max(wt, abs(tn - ta_) / tn), max(wa, ang)
            inc = np.degrees(np.arccos(-nr_a[2]))
            print(f"      {ph} X={np.hypot(*(sta - src)[:2]):5.1f} km z={src[2]:4.1f}: t {ta_:8.5f} s (numerico "
                  f"{tn:8.5f}, recto {np.linalg.norm(sta - src) / phase_velocity(ph, c)[0]:8.5f}), incidencia "
                  f"{inc:5.2f} grados, diferencia angular {ang:.1e} grados")
    check("error relativo del tiempo frente al trazado numerico", wt, 1e-2)
    check("diferencia angular de las direcciones (grados)", wa, 1e-3)
    # 4. transporte de S: conserva la norma y la transversalidad; con g = 0 es la formula de generate_y
    rng = np.random.default_rng(3)
    en = et_ = e0 = 0.0
    for _ in range(50):
        M = rand_dc(rng)
        src = np.array([rng.uniform(-15, 15), rng.uniform(-15, 15), rng.uniform(5, 13)])
        sta = np.array([rng.uniform(-45, 45), rng.uniform(-45, 45), 0.0])
        ns_, nr_ = ray_directions(src, sta, "S", c)
        u = M @ ns_
        vS = radiation_gradient(M, src, sta, "S", c)
        en = max(en, abs(np.linalg.norm(vS) - np.linalg.norm(u - (u @ ns_) * ns_) / np.sqrt(CS)))
        et_ = max(et_, abs(float(vS @ nr_)))
        n = (sta - src) / np.linalg.norm(sta - src)
        u = M @ n
        c0 = make_cfg({"g": 0.0})
        e0 = max(e0, float(np.linalg.norm(radiation_gradient(M, src, sta, "S", c0) - (u - (u @ n) * n) / np.sqrt(CS))),
                 float(np.linalg.norm(radiation_gradient(M, src, sta, "P", c0) - n * float(n @ M @ n) / np.sqrt(CP))))
    check("S transportada: cambio de norma", en, 1e-12)
    check("S transportada: componente sobre la direccion de llegada", et_, 1e-12)
    check("g = 0: radiacion P y S igual a la de generate_y", e0, 1e-12)
    # 5. funcion fuente: Hilbert en forma cerrada frente a la numerica, energia
    from scipy.signal import hilbert
    tf_ = np.arange(-40, 40, 0.001)
    fc = 1.7
    eh_ = float(np.max(np.abs(np.imag(hilbert(source_wavelet(tf_, fc)))[np.abs(tf_) < 3]
                              - (source_wavelet(tf_, fc, np.pi / 2) * -1.0)[np.abs(tf_) < 3])))
    check("Hilbert del Ricker: forma cerrada frente a scipy.signal.hilbert", eh_, 1e-6)
    ee = abs(np.sum(source_wavelet(tf_, fc, 0.6) ** 2) * 0.001 - wavelet_energy(fc)) / wavelet_energy(fc)
    check("energia del pulso rotado frente a la analitica (relativa)", ee, 1e-6)
    check("Ricker de 2 Hz sin rotar igual a ricker() de task0_common", float(np.max(np.abs(source_wavelet(T - 7.3, 2.0) - ricker(T - 7.3)))), 1e-12)
    # 6. consumo de semillas: con g = 0 y VS0 = VS los picks y amplitudes son los de generate_y
    from task0_common import make_network
    from task1o_consistent import draw_sources
    from task1y_clutter import generate_y
    assert abs(1.73 * VS - VP) < 1e-12
    ep = ea_ = 0.0
    nw = 0
    for seed in (5, 6, 7, 8):
        S = make_network("random", 6, seed=seed)
        gens = [(np.random.default_rng(seed + 1), np.random.default_rng(seed + 2), np.random.default_rng(seed + 3))
                for _ in range(2)]
        srcs, taus = draw_sources(gens[0][0])
        _s2, _t2 = draw_sources(gens[1][0])
        _D1, p1, _tr1, a1 = generate_y(S, srcs, taus, 35, *gens[0], 1.75)
        D2, p2, _tr2, a2 = generate_real(S, srcs, taus, 35, *gens[1], 1.75, cfg={"g": 0.0, "VS0": VS})
        nw += len(srcs) > 0
        for x1, x2, b1, b2 in zip(p1, p2, a1, a2):
            assert len(x1) == len(x2) and [b for _t, b in x1] == [b for _t, b in x2]
            ep = max(ep, max(abs(u1[0] - u2[0]) for u1, u2 in zip(x1, x2)))
            ea_ = max(ea_, max(abs(u1 - u2) for u1, u2 in zip(b1, b2)))
        assert D2.shape == _D1.shape and np.all(np.isfinite(D2))
    check(f"g = 0: tiempos de pick iguales a generate_y ({nw} ventanas con fuente)", ep, 1e-9)
    check("g = 0: amplitudes por pick iguales a generate_y", ea_, 1e-9)
    print("autoprueba", "correcta" if ok else "CON FALLAS")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if _selftest() else 1)
