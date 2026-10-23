"""check_adaptive.py - verificacion numerica de la Parte C (estadistico robusto adaptativo).

Configuracion de code/theory_tail_check.py (Proposicion AA.9.5): ruido blanco -> filtro adaptado
Ricker 2 Hz a 16 Hz (varianza unitaria), n_b = 14 bloques de 3 componentes, lectura en la muestra de
mayor energia dentro de +-W = 4 muestras, Q aleatoria ortonormal de 6 columnas.

Comprueba:
  (C.1) la cota corregida (constante de Lipschitz 1 en lugar de kappa) frente a la de AA.9.5(b);
  (C.2) el umbral condicional (Teorema C.2, Hsu-Kakade-Zhang) bajo ruido gaussiano: validez y holgura;
  (C.2) validez del mismo umbral bajo ruido de colas pesadas, no estacionario y con transitorios
       isotropos, frente a umbrales fijos calibrados en ruido gaussiano; y su fallo con transitorios
       polarizados (supuesto violado);
  (C.4) equivariancia de escala, influencia acotada relativa y no saturacion de ACAP.
    python3 check_adaptive.py > check_adaptive_out.txt
"""
import numpy as np

from lib5 import DT, FC

NB, W, C, C1 = 14, 4, 6.0, 3.0
L = 2 * W + 1
JMED = int(np.ceil(NB / 2))           # estadistico de orden usado por ACAP (7.o de 14)
_t = np.arange(-1.0, 1.0 + 1e-9, DT)
_x = (np.pi * FC * _t) ** 2
KER = (1 - 2 * _x) * np.exp(-_x)
KER /= np.linalg.norm(KER)
KAPPA = float(np.abs(np.fft.rfft(KER, 1 << 16)).max())


def mf_reads(z):
    """z (n, NB, 3, L + len(KER) - 1) blanco -> lecturas de pico (n, NB, 3)."""
    n = z.shape[0]
    zz = z.reshape(-1, z.shape[-1])
    x = np.stack([np.convolve(s, KER, mode="valid") for s in zz]).reshape(n, NB, 3, L)
    x = np.moveaxis(x, 2, 3)                                    # (n, NB, L, 3)
    e = (x * x).sum(-1)
    return np.take_along_axis(x, e.argmax(-1)[..., None, None].repeat(3, -1), axis=2)[:, :, 0]


def phi(v, kind):
    t = np.linalg.norm(v, axis=-1, keepdims=True)
    tt = np.maximum(t, 1e-300)
    if kind == "identity":
        return v
    if kind == "cap":
        return v * np.minimum(1.0, C / tt)
    if kind == "log":
        return v * C * np.sqrt(np.log1p(t * t / C ** 2)) / tt
    if kind == "acap":
        q = np.sort(t[..., 0], axis=-1)[..., JMED - 1][..., None, None]
        return v * np.minimum(1.0, C1 * q / tt)
    raise ValueError(kind)


def sigma(v, Q, kind):
    return np.linalg.norm(phi(v, kind).reshape(v.shape[0], -1) @ Q, axis=1)


def cond_threshold(v, Q, kind, t):
    """Teorema C.2: tau^2 = (tr M + 2 sqrt(tr M^2 t) + 2 ||M|| t) / 3, M = sum_k rho_k^2 Q_k^T Q_k."""
    rho2 = (phi(v, kind) ** 2).sum(-1)                          # (n, NB)
    Qb = Q.reshape(NB, 3, 6)
    Gk = np.einsum("kci,kcj->kij", Qb, Qb)                       # (NB, 6, 6)
    M = np.einsum("nk,kij->nij", rho2, Gk)
    ev = np.linalg.eigvalsh(M)
    tr, tr2, op = ev.sum(1), (ev ** 2).sum(1), ev[:, -1]
    return np.sqrt((tr + 2 * np.sqrt(tr2 * t) + 2 * op * t) / 3.0)


def null_reads(rng, n, kind="gauss"):
    shape = (n, NB, 3, L + len(KER) - 1)
    if kind == "gauss":
        return mf_reads(rng.normal(size=shape))
    # colas pesadas: t de Student isotropa (nu = 2) por muestra, escala lognormal por bloque
    g = rng.normal(size=shape)
    chi = rng.chisquare(2, size=(n, NB, 1, shape[-1]))
    z = g * np.sqrt(2.0 / chi)
    scale = np.exp(rng.normal(0, 1.0, size=(n, NB, 1, 1)))
    # transitorios: probabilidad 0.15 por bloque, amplitud log-uniforme 1e2..1e4
    hit = rng.random((n, NB)) < 0.15
    amp = 10 ** rng.uniform(2, 4, size=(n, NB))
    pos = rng.integers(0, shape[-1], size=(n, NB))
    if kind == "heavy_iso":
        d = rng.normal(size=(n, NB, 3))
        d /= np.linalg.norm(d, axis=-1, keepdims=True)
    elif kind == "heavy_pol":
        d = np.zeros((n, NB, 3))
        d[..., 2] = rng.choice([-1.0, 1.0], size=(n, NB))           # transitorio vertical
    else:
        raise ValueError(kind)
    ii, kk = np.nonzero(hit)
    z[ii, kk, :, pos[ii, kk]] += (amp[ii, kk, None] * d[ii, kk])
    # escala por bloque: la estadistica condicional es invariante a ella; los umbrales fijos se evaluan
    # sobre lecturas divididas por la escala verdadera (lo mas favorable a ellos)
    v = mf_reads(z)
    return v * scale[:, :, 0, 0, None], v


def main():
    rng = np.random.default_rng(20261018)
    Q = np.linalg.qr(rng.normal(size=(3 * NB, 6)))[0]
    kinds = ("identity", "cap", "log", "acap")
    out = []
    pr = out.append
    pr(f"kappa = {KAPPA:.3f}; n_b = {NB}; W = {W}; L = {L}; cap c = {C}; ACAP c1 = {C1}, orden {JMED}/{NB}")

    # ---------- C.1 y C.2: ruido gaussiano
    N = 200000
    V = np.concatenate([null_reads(rng, 5000) for _ in range(N // 5000)])
    S = {k: sigma(V, Q, k) for k in kinds}
    nrm = np.linalg.norm(V, axis=-1)
    pr(f"ruido gaussiano, {N} hipotesis: mediana de la norma de la lectura de pico "
       f"{np.median(nrm):.3f}; cuantil {JMED}/{NB} medio {np.sort(nrm, 1)[:, JMED - 1].mean():.3f}"
       f" -> radio efectivo de ACAP {C1 * np.sort(nrm, 1)[:, JMED - 1].mean():.2f}")
    pr("C.1  cota por hipotesis de sigma (gaussiana):")
    pr("   p       AA.9.5(b) (kappa)   corregida (Lip 1)   cuantil empirico: identidad  tope  log  acap")
    for p in (1e-2, 1e-3, 1e-4):
        old = np.sqrt(6) + KAPPA * (np.sqrt(2 * NB * np.log(L)) + np.sqrt(2 * np.log(1 / p)))
        new = np.sqrt(6) + np.sqrt(2 * NB * np.log(L)) + np.sqrt(2 * np.log(1 / p))
        q = "  ".join(f"{np.quantile(S[k], 1 - p):5.2f}" for k in kinds)
        pr(f"   {p:.0e}   {old:6.2f}              {new:6.2f}              {q}")
    pr("C.2  umbral condicional del Teorema C.2 bajo ruido gaussiano")
    pr("   p       tasa realizada (identidad tope log acap)          mediana del umbral (id tope log acap)")
    for p in (1e-2, 1e-3, 1e-4):
        t = np.log(1 / p)
        rates, meds = [], []
        for k in kinds:
            tau = cond_threshold(V, Q, k, t)
            rates.append(f"{np.mean(S[k] > tau):.1e}")
            meds.append(f"{np.median(tau):5.2f}")
        pr(f"   {p:.0e}   {'  '.join(rates)}     {'  '.join(meds)}")
    # umbrales fijos de referencia: cuantil gaussiano 1e-3 de cada mapa
    fixed = {k: np.quantile(S[k], 1 - 1e-3) for k in kinds}

    # cuantil condicional exacto por Monte Carlo de direcciones (holgura de la cota HKZ)
    sub = V[:400]
    rho = np.linalg.norm(phi(sub, "identity"), axis=-1)
    Qb = Q.reshape(NB, 3, 6)
    ratio = []
    for i in range(len(sub)):
        d = rng.normal(size=(20000, NB, 3))
        d /= np.linalg.norm(d, axis=-1, keepdims=True)
        s = np.linalg.norm(np.einsum("bk,bkc,kci->bi", np.broadcast_to(rho[i], (20000, NB)), d, Qb), axis=1)
        ratio.append(cond_threshold(sub[i:i + 1], Q, "identity", np.log(1e3))[0] / np.quantile(s, 1 - 1e-3))
    pr(f"   holgura de la cota HKZ frente al cuantil condicional exacto (p = 1e-3, identidad, 400 casos): "
       f"cociente mediano {np.median(ratio):.3f}, rango {np.min(ratio):.3f}..{np.max(ratio):.3f}")

    # ---------- C.3: ruido de colas pesadas
    for kind in ("heavy_iso", "heavy_pol"):
        Nh = 50000
        pairs = [null_reads(rng, 5000, kind) for _ in range(Nh // 5000)]
        Vh = np.concatenate([p[0] for p in pairs])
        Vn = np.concatenate([p[1] for p in pairs])
        pr(f"C.3  ruido {kind}, {Nh} hipotesis (t de Student nu = 2, escala lognormal por bloque, "
           f"transitorios 15 % de bloques, 1e2..1e4):")
        pr("   p       umbral condicional: tasa (id tope log acap)    umbral fijo gaussiano 1e-3: tasa")
        Sh = {k: sigma(Vh, Q, k) for k in kinds}
        Sn = {k: sigma(Vn, Q, k) for k in kinds}
        for p in (1e-2, 1e-3):
            rates = [f"{np.mean(Sh[k] > cond_threshold(Vh, Q, k, np.log(1 / p))):.1e}" for k in kinds]
            fr = "" if p != 1e-3 else "  ".join(f"{np.mean(Sn[k] > fixed[k]):.1e}" for k in kinds)
            pr(f"   {p:.0e}   {'  '.join(rates)}          {fr}")

    # ---------- C.4: propiedades deterministas de ACAP
    v = V[:2000].copy()
    a = 37.0
    eq = np.max(np.abs(sigma(a * v, Q, "acap") - a * sigma(v, Q, "acap")) / (a * sigma(v, Q, "acap")))
    pr(f"C.4  equivariancia de escala de ACAP: error relativo maximo {eq:.1e}")
    worst, bound = 0.0, np.inf
    viol = 0
    for m in (1, 3, 6):
        for rep in range(300):
            x = V[rep].copy()
            ks = rng.choice(NB, m, replace=False)
            xr = x.copy()
            xr[ks] = rng.normal(size=(m, 3)) * 10 ** rng.uniform(0, 5, size=(m, 1))
            s0 = sigma(x[None], Q, "acap")[0]
            s1 = sigma(xr[None], Q, "acap")[0]
            clean = np.sort(np.linalg.norm(np.delete(x, ks, 0), axis=1))
            qm, qp = clean[max(JMED - m, 1) - 1], clean[min(JMED, NB - m) - 1]
            b = C1 * (2 * np.sqrt(m) * qp + np.sqrt(NB - m) * (qp - qm))
            viol += abs(s1 - s0) > b + 1e-9
            worst = max(worst, abs(s1 - s0) / b)
        pr(f"     m = {m} bloques reemplazados (amplitud hasta 1e5): violaciones de la Prop. C.4(b) = {viol}, "
           f"max |cambio| / cota = {worst:.3f}")
    ev = rng.normal(size=6)
    w = Q @ ev
    w = w.reshape(NB, 3)
    pr("     no saturacion: evento a*w sobre un fondo saturado (todas las normas del fondo = 20)")
    bg = rng.normal(size=(NB, 3))
    bg *= 20 / np.linalg.norm(bg, axis=1, keepdims=True)
    for amp in (0, 1, 3, 10, 30, 100):
        x = (bg + amp * w / np.linalg.norm(w) * np.sqrt(NB) * 20)[None]
        pr(f"       a = {amp:4d}:  identidad {sigma(x, Q, 'identity')[0]:9.1f}  tope {sigma(x, Q, 'cap')[0]:6.2f}"
           f"  log {sigma(x, Q, 'log')[0]:7.2f}  acap {sigma(x, Q, 'acap')[0]:9.1f}   (cota del tope "
           f"c sqrt(n_b) = {C * np.sqrt(NB):.2f})")
    print("\n".join(out))


if __name__ == "__main__":
    main()
