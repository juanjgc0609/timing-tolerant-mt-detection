"""check_identifiability.py - verificacion numerica de la Parte A (identificabilidad de eventos solapados).

  A.1  residuo exacto de un arribo (dos pulsos frente a uno) contra los terminos principales del
       Teorema A.4 (primer orden con diversidad de polarizacion, segundo orden con polarizaciones
       paralelas) y contra la cota de separacion grande (Proposicion A.5');
  A.2  margen de identificabilidad de red delta_1 (distancia exacta de la senal de dos eventos al modelo
       de un evento, optimizado sobre hipocentro, tiempo de origen y tensor) frente a la cota inferior
       por relajacion por estacion (Teorema A.3) y la formula linealizada (Corolario A.6), en tres
       escenarios: colocados con mecanismos distintos, colocados con el mismo mecanismo, y separados
       2 km con el mismo mecanismo;
  A.3  pendientes log-log y separacion minima resoluble para un error total de 0.1 (Teorema A.1).
    python3 check_identifiability.py > check_identifiability_out.txt
"""
import numpy as np
from scipy.optimize import minimize, minimize_scalar

from lib5 import OMEGA2, PSI, acorr, network, radiation, rand_dc, ttimes

STA = network(8)
NS = len(STA)


# ---------------------------------------------------------------- A.1: un arribo
def R_exact(g1, g2, D):
    h2 = g1 @ g1 + g2 @ g2 + 2 * (g1 @ g2) * acorr(D)
    f = lambda e: -np.sum((g1 * acorr(e) + g2 * acorr(e - D)) ** 2)
    grid = np.linspace(-0.6 + min(0, D), 0.6 + max(0, D), 4001)
    vals = np.array([f(e) for e in grid])
    e0 = grid[vals.argmin()]
    r = minimize_scalar(f, bounds=(e0 - 0.002, e0 + 0.002), method="bounded", options={"xatol": 1e-12})
    return h2 + min(r.fun, vals.min())


def part_a1(out, rng):
    out.append("A.1  un arribo: residuo exacto R frente a los terminos principales")
    g1, g2 = rng.normal(size=3), rng.normal(size=3)
    cross = np.sum(np.cross(g1, g2) ** 2) / np.sum((g1 + g2) ** 2)
    out.append(f"   polarizaciones distintas: |g1| {np.linalg.norm(g1):.3f} |g2| {np.linalg.norm(g2):.3f}; "
               f"termino de primer orden omega^2 D^2 |g1 x g2|^2/|g1+g2|^2")
    for D in (0.1, 0.03, 0.01, 0.003, 0.001):
        out.append(f"     D = {D:6.3f} s: R = {R_exact(g1, g2, D):.4e}; principal {OMEGA2 * D * D * cross:.4e}; "
                   f"cociente {R_exact(g1, g2, D) / (OMEGA2 * D * D * cross):.5f}")
    e = rng.normal(size=3)
    e /= np.linalg.norm(e)
    al, be = 1.3, 0.6
    c4 = (PSI - OMEGA2 ** 2) / 4 * (al * be / (al + be)) ** 2
    out.append(f"   polarizaciones paralelas (|g1| = {al}, |g2| = {be}): termino (psi - omega^4)/4 D^4 "
               f"(|g1||g2|/(|g1|+|g2|))^2")
    for D in (0.1, 0.03, 0.01, 0.003):
        r = R_exact(al * e, be * e, D)
        out.append(f"     D = {D:6.3f} s: R = {r:.4e}; principal {c4 * D ** 4:.4e}; cociente {r / (c4 * D ** 4):.5f}")
    out.append("   separacion grande (Prop. A.5'): R >= m^2 (1 - alpha^2) - 4 alpha M m, alpha = sup_{|u|>=D/2}|a(u)|")
    uu = np.linspace(0, 5, 50001)
    for D in (0.5, 0.75, 1.0, 1.5, 2.0):
        alpha = np.abs(acorr(uu[uu >= D / 2])).max()
        worst = np.inf
        for _ in range(200):
            a1, a2 = rng.normal(size=3), rng.normal(size=3) * rng.uniform(0.1, 3)
            M, m = max(np.linalg.norm(a1), np.linalg.norm(a2)), min(np.linalg.norm(a1), np.linalg.norm(a2))
            worst = min(worst, R_exact(a1, a2, D) - (m * m * (1 - alpha ** 2) - 4 * alpha * M * m))
        out.append(f"     D = {D:4.2f} s: alpha = {alpha:.2e}; min sobre 200 pares de (R - cota) = {worst:.3e}")


# ---------------------------------------------------------------- A.2: red
def event(x, tau, m):
    A = radiation(x, STA).reshape(2, NS, 3, 6)
    g = np.einsum("bicj,j->bic", A, m)                      # (2, NS, 3)
    return ttimes(x, tau, STA), g


def one_event_residual(theta, evs, E):
    x, tau = theta[:3], theta[3]
    t = ttimes(x, tau, STA)                                 # (2, NS)
    A = radiation(x, STA).reshape(2, NS, 3, 6)
    # Gram del modelo: sum_i sum_{b,b'} A_ib^T A_ib' a(t_ib - t_ib')
    C = acorr(t[:, None, :] - t[None, :, :])                # (b, b', i)
    Gm = np.einsum("bicj,Bick,bBi->jk", A, A, C)
    bvec = np.zeros(6)
    for (tk, gk) in evs:
        Ck = acorr(t[:, None, :] - tk[None, :, :])          # (b, b', i)
        bvec += np.einsum("bicj,Bic,bBi->j", A, gk, Ck)
    try:
        val = bvec @ np.linalg.solve(Gm, bvec)
    except np.linalg.LinAlgError:
        val = 0.0
    return E - val


def total_energy(evs):
    E = 0.0
    for (t1, g1) in evs:
        for (t2, g2) in evs:
            C = acorr(t1[:, None, :] - t2[None, :, :])      # (b, b', i)
            E += np.einsum("bic,Bic,bBi->", g1, g2, C)
    return E


def delta1(evs, starts):
    E = total_energy(evs)
    best = np.inf
    for s in starts:
        r = minimize(one_event_residual, s, args=(evs, E), method="Nelder-Mead",
                     options={"xatol": 1e-6, "fatol": 1e-12, "maxiter": 4000, "maxfev": 8000})
        best = min(best, r.fun)
    return np.sqrt(max(best, 0.0)), E


def station_relaxation(evs):
    """sum_i min_{tP, tS, gP, gS} ||d_i - gP s_tP - gS s_tS||^2 (Teorema A.3), sin suponer P/S separadas."""
    tot = 0.0
    for i in range(NS):
        pul = [(tk[b, i], gk[b, i]) for (tk, gk) in evs for b in range(2)]
        di2 = sum(g @ h * acorr(t - u) for (t, g) in pul for (u, h) in pul)

        def neg(tt):
            tp, ts = tt
            c = acorr(ts - tp)
            Gi = np.array([[1, c], [c, 1]])
            bp = sum(g * acorr(tp - t) for (t, g) in pul)
            bs = sum(g * acorr(ts - t) for (t, g) in pul)
            B = np.stack([bp, bs])                          # (2, 3)
            if abs(c) > 0.999999:
                return -np.sum(bp ** 2)
            return -np.trace(B.T @ np.linalg.solve(Gi, B))
        # busqueda global vectorizada sobre tiempos candidatos (rejillas locales de +-0.3 s alrededor de
        # cada pulso de los datos y de los tiempos medios ponderados), luego refinamiento local
        w = [g @ g for (_t, g) in pul]
        tp0 = sum(w[j] * pul[j][0] for j in (0, 2)) / (w[0] + w[2])
        ts0 = sum(w[j] * pul[j][0] for j in (1, 3)) / (w[1] + w[3])
        cand = np.unique(np.concatenate([c + np.linspace(-0.3, 0.3, 61)
                                         for c in [p[0] for p in pul] + [tp0, ts0]]))
        Bc = sum(np.outer(acorr(cand - t), g) for (t, g) in pul)          # (nc, 3)
        nb2 = np.sum(Bc ** 2, 1)
        cc = acorr(cand[:, None] - cand[None, :])
        dots = Bc @ Bc.T
        with np.errstate(divide="ignore", invalid="ignore"):
            val = (nb2[:, None] + nb2[None, :] - 2 * cc * dots) / (1 - cc ** 2)
        val[np.abs(cc) > 0.999999] = -np.inf
        np.fill_diagonal(val, nb2)
        j1, j2 = np.unravel_index(np.nanargmax(val), val.shape)
        x = np.array([cand[j1], cand[j2]])
        bestv = neg(x)
        for _ in range(30):
            for k in (0, 1):
                f1 = lambda u: neg((u, x[1])) if k == 0 else neg((x[0], u))
                r1 = minimize_scalar(f1, bounds=(x[k] - 0.02, x[k] + 0.02), method="bounded",
                                     options={"xatol": 1e-12})
                if r1.fun < bestv:
                    bestv = r1.fun
                    x[k] = r1.x
        r = minimize(neg, x, method="Nelder-Mead",
                     options={"xatol": 1e-12, "fatol": 1e-16, "initial_simplex": [x, x + [1e-3, 0], x + [0, 1e-3]]})
        tot += di2 + min(bestv, r.fun)
    return np.sqrt(max(tot, 0.0))


def linearized(evs):
    (t1, g1), (t2, g2) = evs
    D = t2 - t1
    lin1 = OMEGA2 * np.sum(D ** 2 * np.sum(np.cross(g1, g2) ** 2, -1) / np.sum((g1 + g2) ** 2, -1))
    n1, n2 = np.linalg.norm(g1, axis=-1), np.linalg.norm(g2, axis=-1)
    lin2 = (PSI - OMEGA2 ** 2) / 4 * np.sum(D ** 4 * (n1 * n2 / (n1 + n2)) ** 2)
    return np.sqrt(lin1), np.sqrt(lin2)


def part_a2(out, rng):
    x0 = np.array([2.0, -3.0, 9.0])
    m1 = rand_dc(rng)
    m2 = rand_dc(rng)
    # normalizar: cada evento con energia de red 1
    def norm_m(x, m):
        _t, g = event(x, 0.0, m)
        return m / np.sqrt(np.sum(g ** 2))
    scen = {
        "S1 colocados, mecanismos distintos": (np.zeros(3), m1, m2),
        "S2 colocados, mismo mecanismo": (np.zeros(3), m1, m1),
        "S3 mismo mecanismo y origen, separacion horizontal dx (columna = dx en km)": (np.array([1.0, 0.0, 0.0]), m1, m1),
    }
    dts = (0.003, 0.01, 0.03, 0.1, 0.3, 1.0)
    out.append("A.2  margen de red delta_1 / sqrt(E1) (E1 = E2 = 1), red 'spread' de 8 estaciones, fuente "
               "en (2, -3, 9) km")
    res = {}
    for name, (dx, ma, mb) in scen.items():
        ma, mb = norm_m(x0, ma), norm_m(x0 + dx, mb)
        out.append(f"   {name}")
        out.append("     dtau(s)   delta_1 exacto   cota relajacion   lineal 1er orden   lineal 2o orden")
        rows = []
        for dt in (dts if name[:2] != "S3" else (0.03, 0.1, 0.3, 1.0, 2.0)):
            if name[:2] == "S3":
                dx = np.array([dt, 0.0, 0.0])
                mb = norm_m(x0 + dx, m1)
                evs = [event(x0, 5.0, ma), event(x0 + dx, 5.0, mb)]
                mid = np.concatenate([x0 + dx / 2, [5.0]])
                starts = [np.concatenate([x0, [5.0]]), np.concatenate([x0 + dx, [5.0]]), mid,
                          mid + np.array([0.3, 0.3, 0.3, 0.03]), mid - np.array([0.3, -0.3, 0.3, 0.03])]
                d1, _E = delta1(evs, starts)
                rl = station_relaxation(evs)
                l1, l2 = linearized(evs)
                rows.append((dt, d1, rl, l1, l2))
                out.append(f"     {dt:7.3f}   {d1:12.4e}   {rl:13.4e}   {l1:14.4e}   {l2:14.4e}")
                continue
            evs = [event(x0, 5.0, ma), event(x0 + dx, 5.0 + dt, mb)]
            mid = np.concatenate([x0 + dx / 2, [5.0 + dt / 2]])
            starts = [np.concatenate([x0, [5.0]]), np.concatenate([x0 + dx, [5.0 + dt]]), mid,
                      mid + np.array([0.5, 0.5, 0.5, 0.05]), mid - np.array([0.5, -0.5, 0.5, 0.05])]
            d1, _E = delta1(evs, starts)
            rl = station_relaxation(evs)
            l1, l2 = linearized(evs)
            rows.append((dt, d1, rl, l1, l2))
            out.append(f"     {dt:7.3f}   {d1:12.4e}   {rl:13.4e}   {l1:14.4e}   {l2:14.4e}")
        res[name] = rows
    return res


def part_a3(out, res):
    out.append("A.3  pendientes log-log de delta_1 en dtau pequeno y separacion minima resoluble")
    z = 3.29                      # 2 Phi(-z/2) = 0.10
    for name, rows in res.items():
        r = [x for x in rows if x[0] > 0]
        dt = np.array([x[0] for x in r])
        d = np.array([x[1] for x in r])
        k = dt <= (0.3 if name[:2] == "S3" else 0.03)
        un = "km" if name[:2] == "S3" else "s"
        if k.sum() >= 2 and np.all(d[k] > 0):
            sl = np.polyfit(np.log(dt[k]), np.log(d[k]), 1)[0]
        else:
            sl = np.nan
        msg = f"   {name}: pendiente {sl:.2f}"
        for snr in (12.4, 30.0, 100.0):
            need = z / snr
            # interpolacion log-log sobre la curva exacta
            ok = np.nonzero(d >= need)[0]
            if len(ok) == 0:
                msg += f"; SNR {snr:g}: > {dt.max():g} s"
                continue
            j = ok[0]
            if j == 0:
                msg += f"; SNR {snr:g}: <= {dt[0]:g} s"
                continue
            f = (np.log(need) - np.log(d[j - 1])) / (np.log(d[j]) - np.log(d[j - 1]))
            msg += f"; SNR {snr:g}: {np.exp(np.log(dt[j - 1]) + f * (np.log(dt[j]) - np.log(dt[j - 1]))):.3f} {un}"
        out.append(msg)


def main():
    rng = np.random.default_rng(5)
    out = [f"Ricker {2.0} Hz: omega^2 = ||s'||^2 = {OMEGA2:.2f} s^-2, psi = ||s''||^2 = {PSI:.1f} s^-4, "
           f"psi - omega^4 = {PSI - OMEGA2 ** 2:.1f} s^-4"]
    part_a1(out, rng)
    res = part_a2(out, rng)
    part_a3(out, res)
    print("\n".join(out))


if __name__ == "__main__":
    main()
