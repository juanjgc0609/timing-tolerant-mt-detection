"""check_extensions.py - verificacion numerica de la Parte D (componentes del modelo nuevo).

  D.1  verosimilitud de contaminacion por bloque: influencia de un transitorio frente a su amplitud
       (crece como log(1/f_fondo)) y evidencia de un evento frente a su amplitud sobre un fondo denso, con
       error de modelo proporcional (crece como nu log a);
  D.2  el estadistico de mezcla sobre el tensor es un e-valor: P(S_mix >= log(1/p)) <= p para cualquier
       fondo correctamente especificado; comparacion con el GLR (maximo sobre m) y con un fondo mal
       especificado;
  D.4  angulo principal entre los subespacios de dos eventos solapados en la red frente a la separacion:
       sigma_min(Gamma)^2 = 1 - cos(theta_1) escala como dtau^2, de modo que la condicion suficiente de AA.6(a)
       pide SNR ~ dtau^-3 frente al dtau^-1 necesario del Corolario A.7.
    python3 check_extensions.py > check_extensions_out.txt
"""
import numpy as np
from scipy.special import gammaln, logsumexp

from lib5 import acorr, network, radiation, rand_dc, ttimes

NB = 14
NU = 3.0


def logt3(v, b, nu=NU):
    """log densidad de la t isotropa en R^3 con escala b y nu grados de libertad."""
    r2 = np.sum(v * v, -1) / (b * b)
    return (gammaln((nu + 3) / 2) - gammaln(nu / 2) - 1.5 * np.log(nu * np.pi) - 3 * np.log(b)
            - (nu + 3) / 2 * np.log1p(r2 / nu))


def lognorm3(v, mu, s):
    return -1.5 * np.log(2 * np.pi * s * s) - np.sum((v - mu) ** 2, -1) / (2 * s * s)


def ell(v, mu, s, b, eps, nu=NU):
    """log[(eps f + (1-eps) N) / f] por bloque."""
    lf = logt3(v, b, nu)
    return np.logaddexp(np.log(eps), np.log1p(-eps) + lognorm3(v, mu, s) - lf)


def sample_t3(rng, n, b, nu):
    g = rng.normal(size=n + (3,))
    w = rng.chisquare(nu, size=n + (1,))
    return b[..., None] * g * np.sqrt(nu / w) if np.ndim(b) else b * g * np.sqrt(nu / w)


def main():
    rng = np.random.default_rng(7)
    out = []
    eps = 0.05
    Q = np.linalg.qr(rng.normal(size=(3 * NB, 6)))[0]
    Qb = Q.reshape(NB, 3, 6)

    # ---------------- D.1 influencia y no saturacion
    out.append(f"D.1  verosimilitud de contaminacion, eps = {eps}, fondo t isotropa nu = {NU}")
    out.append("   influencia maxima de un bloque (sup_m ell - inf ell) frente a la amplitud M del transitorio:")
    s, b = 1.0, 1.0
    for M in (1e1, 1e2, 1e3, 1e4):
        v = np.array([M, 0, 0])
        inf_l = ell(v, v + 1e9, s, b, eps)            # prediccion lejana: ell -> log eps
        sup_l = ell(v, v, s, b, eps)                  # prediccion exacta
        bound = np.log(eps + (1 - eps) * np.exp(lognorm3(v, v, s) - logt3(v, b)))
        out.append(f"     M = {M:7.0e}: sup ell = {sup_l:7.2f}, cota log(eps+(1-eps) N_max/f) = {bound:7.2f},"
                   f" inf ell = {inf_l:6.2f}; prediccion (nu+3) log M = {(NU + 3) * np.log(M):6.2f}")
    out.append("   evidencia de un evento exacto a*w por bloque sobre fondo denso b = 20, error de modelo "
               "s = 1 + gamma |a w| (gamma = 0.2):")
    w = np.array([1.0, 0, 0])
    prev = None
    for a in (30, 100, 300, 1000, 3000):
        v = a * w
        sk = 1 + 0.2 * a
        e = ell(v, v, sk, 20.0, eps)
        sl = "" if prev is None else f"; pendiente en log a = {(e - prev[1]) / np.log(a / prev[0]):.2f}"
        out.append(f"     a = {a:5d}: ell = {e:7.3f}{sl}")
        prev = (a, e)
    out.append(f"     (prediccion de la Proposicion D.1(b): pendiente nu = {NU})")

    # ---------------- D.2 e-valor
    out.append("D.2  S_mix = log E_pi[prod_k exp(ell_k(m))], pi = N(0, 10^2 I), 4000 muestras de m; "
               "nulo: bloques t isotropos")
    ms = rng.normal(scale=10.0, size=(4000, 6))
    mu = np.einsum("kci,ji->jkc", Qb, ms)                    # (4000, NB, 3)
    scales = np.exp(rng.normal(0, 1, NB))
    for label, nu_true, sc in (("fondo bien especificado", NU, 1.0),
                               ("fondo real con colas mas pesadas (nu 1.5)", 1.5, 1.0),
                               ("fondo real 2x mas ancho", NU, 2.0)):
        Smix, Sglr = [], []
        for rep in range(4000):
            v = sample_t3(rng, (NB,), scales * sc, nu_true)
            L = ell(v[None], mu, 1.0, scales[None], eps).sum(1)     # (4000,)
            Smix.append(logsumexp(L) - np.log(len(L)))
            Sglr.append(L.max())
        Smix, Sglr = np.array(Smix), np.array(Sglr)
        rows = []
        for p in (0.1, 0.01, 1e-3):
            rows.append(f"p={p:g}: mezcla {np.mean(Smix >= np.log(1 / p)):.4f}, GLR {np.mean(Sglr >= np.log(1 / p)):.4f}")
        out.append(f"   {label}: media de exp(S_mix) = {np.mean(np.exp(np.minimum(Smix, 50))):.3f}; " + "; ".join(rows))

    # ---------------- D.4 angulos principales
    STA = network(8)
    ns = len(STA)
    x0 = np.array([2.0, -3.0, 9.0])

    def atoms_gram(ta, Aa, tb, Ab):
        C = acorr(ta[:, None, :] - tb[None, :, :])          # (b, b', i)
        return np.einsum("bicj,Bick,bBi->jk", Aa, Ab, C)

    A = radiation(x0, STA).reshape(2, ns, 3, 6)
    out.append("D.4  colocados, origen separado dtau: 1 - cos(theta_1) entre los subespacios de las dos hipotesis")
    prev = None
    for dt in (0.003, 0.01, 0.03, 0.1, 0.3):
        t1, t2 = ttimes(x0, 5.0, STA), ttimes(x0, 5.0 + dt, STA)
        G11, G22, G12 = atoms_gram(t1, A, t1, A), atoms_gram(t2, A, t2, A), atoms_gram(t1, A, t2, A)
        w1, V1 = np.linalg.eigh(G11)
        w2, V2 = np.linalg.eigh(G22)
        R1 = V1 @ np.diag(w1 ** -0.5) @ V1.T
        R2 = V2 @ np.diag(w2 ** -0.5) @ V2.T
        c = np.linalg.svd(R1 @ G12 @ R2, compute_uv=False)
        gap = 1 - c.max()
        sl = "" if prev is None else f"; pendiente {np.log(gap / prev[1]) / np.log(dt / prev[0]):.2f}"
        out.append(f"     dtau = {dt:5.3f} s: 1 - cos theta_1 = {gap:.3e}, sigma_min(Gamma) = {np.sqrt(gap):.3e}{sl}")
        prev = (dt, gap)
    print("\n".join(out))


if __name__ == "__main__":
    main()
