"""check_greedy.py - verificacion numerica de la Parte B (busqueda codiciosa con diccionario perturbado).

Detector desplegado sin tope (MT-FW NOCAP, igual a code/task1x_mtglrt.fw_record con tope infinito):
lectura de pico de energia dentro de +-W = 4 muestras del tiempo predicho (redondeado hacia abajo),
puntaje ||Pi_j v||, maximo global, ajuste por minimos cuadrados, resta de la autocorrelacion muestreada en
los tiempos leidos, exclusion de 5 km, dos iteraciones. Dos fuentes de doble par EN nodos de la grilla
(3 km, 15 x 15 x 3, origenes cada 0.1 s) con estaticas e_ib ~ U[-eps, eps] comunes a ambas fuentes (el
error de tiempo del diccionario), arribos fuera de la rejilla de muestreo, ruido gaussiano filtrado.

Comprueba:
  B.1  Lema B.2 (lectura en el nodo correcto): cota inferior LB_k <= puntaje del nodo verdadero (100 %);
  B.2  Lema B.3 (fuga hacia nodos equivocados) en registros sin ruido: puntaje <= sum_k Lambda_jk (100 %);
  B.3  Teorema B.4 (certificado de seleccion correcta en el primer paso): si certifica, la seleccion es
       correcta (100 %); tasa de certificacion frente a la tasa real de acierto;
  B.4  Proposicion B.6 (residuo de la deflacion desplegada por desfase submuestral) y exito en K = 2.
    python3 check_greedy.py > check_greedy_out.txt
"""
import itertools

import numpy as np

from lib5 import DT, FC, acorr, basis, network, radiation, rand_dc, ttimes

STA = network(8)
NS = len(STA)
NB = 2 * NS
W = 4
L = 2 * W + 1
OFF = np.arange(-W, W + 1)
NT = 640
TN = np.arange(NT) * DT
AMIN = float(-acorr(np.linspace(0, 2, 200001)).min())       # a^- = 0.618
XS = np.arange(-21, 21.01, 3.0)
ZS = np.array([5.0, 9.0, 13.0])
NODES = np.array([[x, y, z] for x in XS for y in XS for z in ZS])
TAUS = np.round(np.arange(3.0, 12.01, 0.1), 3)
NJ = len(NODES)
_t = np.arange(-1.0, 1.0 + 1e-9, DT)
_x = (np.pi * FC * _t) ** 2
KER = (1 - 2 * _x) * np.exp(-_x)
KER /= np.linalg.norm(KER)

Qs = np.stack([basis(radiation(x, STA)) for x in NODES])      # (NJ, 48, 6)
TT = np.stack([np.stack([ttimes(x, 0.0, STA) for x in NODES])])[0]   # (NJ, 2, NS) sin origen
IDX = np.floor((TT[:, None] + TAUS[None, :, None, None]) / DT + 1e-9).astype(int)   # (NJ, ntau, 2, NS)
WIN = np.clip(IDX[..., None] + OFF, 0, NT - 1)                 # (NJ, ntau, 2, NS, L)


def noise(rng):
    z = rng.normal(size=(NS, 3, NT + len(KER) - 1))
    return np.stack([[np.convolve(z[i, c], KER, mode="valid") for c in range(3)] for i in range(NS)])


def field(srcs):
    """campo sin ruido del filtro adaptado: sum_k sum_b g a(t - t_arr) -> (NS, 3, NT)."""
    Y = np.zeros((NS, 3, NT))
    for (tarr, g) in srcs:
        for b in range(2):
            for i in range(NS):
                Y[i] += np.outer(g[b, i], acorr(TN - tarr[b, i]))
    return Y


def score_all(Y, excl):
    """puntaje ||Pi_j v||, lectura y tiempos leidos para todas las hipotesis admisibles."""
    best = (-np.inf, None)
    st = np.arange(NS)
    for j in range(NJ):
        if excl[j]:
            continue
        seg = Y[st[None, None, :, None], :, WIN[j]]           # (ntau, 2, NS, L, 3)
        e = (seg ** 2).sum(-1)
        k = e.argmax(-1)
        v = np.take_along_axis(seg, k[..., None, None], axis=3)[..., 0, :]     # (ntau, 2, NS, 3)
        vf = v.reshape(len(TAUS), -1)
        sc = np.linalg.norm(vf @ Qs[j], axis=1)
        it = int(sc.argmax())
        if sc[it] > best[0]:
            tpk = np.take_along_axis(WIN[j, it], k[it][..., None], axis=-1)[..., 0]
            best = (float(sc[it]), (j, it, vf[it], tpk))
    return best


def deflate(Y, j, vf, tpk):
    A = radiation(NODES[j], STA)
    m, *_ = np.linalg.lstsq(A, vf, rcond=None)
    pred = (A @ m).reshape(2, NS, 3)
    Y = Y.copy()
    for b in range(2):
        for i in range(NS):
            Y[i] -= np.outer(pred[b, i], acorr(TN - tpk[b, i] * DT))
    return Y


def is_hit(j, it, src_nodes, src_taus):
    return [k for k in range(len(src_nodes)) if np.linalg.norm(NODES[j, :2] - NODES[src_nodes[k], :2]) <= 5.0
            and abs(TAUS[it] - src_taus[k]) <= 1.0]


def lower_bound(j, it, g, tarr, Yother):
    """Lema B.2 en el nodo verdadero (j, it): g (2, NS, 3) y tarr (2, NS) de la fuente, Yother el resto."""
    gn = np.linalg.norm(g)
    lb = 0.0
    for b in range(2):
        for i in range(NS):
            w = WIN[j, it, b, i]
            gi = np.linalg.norm(g[b, i])
            # interferencia: ruido + otras fuentes + la otra fase de la misma fuente
            other = Yother[i][:, w] + np.outer(g[1 - b, i], acorr(TN[w] - tarr[1 - b, i]))
            Z = np.linalg.norm(other, axis=0).max()
            off = np.abs(TN[w] - tarr[b, i])
            if off.min() <= DT / 2 + 1e-12:
                lam = acorr(off.min()) - 2 * Z / max(gi, 1e-300)
            else:
                lam = -np.inf
            if lam > AMIN:
                lb += lam * gi * gi - gi * Z
            else:
                y = np.outer(g[b, i], acorr(TN[w] - tarr[b, i])) + other
                lb -= gi * np.linalg.norm(y, axis=0).max()
    return lb / gn


def alpha_tables(srcs):
    """alpha[k] (NJ, ntau, 2[b], NS, 2[b']) = max_ventana |a(t - t^k_{i b'})|."""
    out = []
    tw = WIN * DT                                               # (NJ, ntau, 2, NS, L)
    for (tarr, _g) in srcs:
        d = tw[..., None] - tarr.T[None, None, None, :, None, :]  # (NJ, ntau, 2, NS, L, 2)
        out.append(np.abs(acorr(d)).max(axis=4))
    return out


def leak_cheap(al, g):
    """cota de Cauchy-Schwarz de Lambda_jk ||g_k||: sqrt(sum_ib (sum_b' alpha ||g_ib'||)^2)."""
    gn = np.linalg.norm(g, axis=-1)                             # (2, NS) por b'
    s = np.einsum("...bic,ci->...bi", al, gn)                   # (NJ, ntau, 2, NS)
    return np.sqrt((s ** 2).sum((-1, -2)))


def leak_exact(j, al, g, cap=16):
    """max sobre signos de ||Pi_j sum_r sigma_r W_r|| (relajacion exacta del Lema B.3), mas un resto de C-S."""
    Wr, mag = [], []
    for b in range(2):
        for i in range(NS):
            for bp in range(2):
                if al[b, i, bp] > 1e-9:
                    v = np.zeros((2, NS, 3))
                    v[b, i] = al[b, i, bp] * g[bp, i]
                    Wr.append(v.reshape(-1))
                    mag.append(np.linalg.norm(v))
    if not Wr:
        return 0.0
    order = np.argsort(mag)[::-1]
    Wr = np.array(Wr)[order]
    mag = np.array(mag)[order]
    head, rest = Wr[:cap], mag[cap:]
    P = head @ Qs[j]                                            # (r, 6)
    S = _signs(len(head))
    best = float(np.linalg.norm(S @ P, axis=1).max())
    return best + np.sqrt((rest ** 2).sum())


_SIGNS = {}


def _signs(r):
    if r not in _SIGNS:
        _SIGNS[r] = np.array([(1.0,) + sg for sg in itertools.product((1.0, -1.0), repeat=r - 1)])
    return _SIGNS[r]


def main():
    rng = np.random.default_rng(31)
    out = [f"red 'spread' 8 estaciones, n_b = {NB}, W = {W} muestras, a^- = {AMIN:.4f}, a(DT/2) = "
           f"{acorr(DT / 2):.4f}; grilla {NJ} nodos x {len(TAUS)} origenes"]
    # cota de ruido del Lema C.1 corregido, union sobre todas las hipotesis, p = 0.05
    Nh = NJ * len(TAUS)
    NW = np.sqrt(6) + np.sqrt(2 * NB * np.log(L)) + np.sqrt(2 * np.log(Nh / 0.05))
    out.append(f"cota de ruido N_W (p = 0.05, {Nh} hipotesis) = {NW:.2f}")
    rows = []
    viol_lb = viol_leak = viol_cert = 0
    nlb = nleak = 0
    left = []
    for snr in (30.0, 60.0, 120.0):
        for dtau in (0.0, 1.0, 3.0):
            for eps in (0.0, 0.15):
                ncert = ncorr1 = nsucc = 0
                nrep = 12
                for rep in range(nrep):
                    # fuentes en nodos, separadas >= 9 km en planta
                    while True:
                        a, b = rng.choice(NJ, 2, replace=False)
                        if np.linalg.norm(NODES[a, :2] - NODES[b, :2]) >= 9.0:
                            break
                    t1 = 5.0 + 0.1 * rng.integers(0, 10)
                    taus = [round(t1, 3), round(t1 + dtau, 3)]
                    stat = rng.uniform(-eps, eps, size=(2, NS))
                    srcs = []
                    for k, nd in enumerate((a, b)):
                        m = rand_dc(rng)
                        g = (radiation(NODES[nd], STA) @ m).reshape(2, NS, 3)
                        g *= snr / np.linalg.norm(g)
                        srcs.append((ttimes(NODES[nd], taus[k], STA) + stat, g))
                    z = noise(rng)
                    Yk = [field([s]) for s in srcs]
                    Y = Yk[0] + Yk[1] + z
                    itrue = [int(np.argmin(np.abs(TAUS - t))) for t in taus]
                    # B.1: cota inferior en el nodo verdadero
                    lbs = []
                    for k, nd in enumerate((a, b)):
                        lb = lower_bound(nd, itrue[k], srcs[k][1], srcs[k][0], Yk[1 - k] + z)
                        st = np.arange(NS)
                        seg = Y[st[None, :, None], :, WIN[nd, itrue[k]]]
                        e = (seg ** 2).sum(-1)
                        v = np.take_along_axis(seg, e.argmax(-1)[..., None, None], axis=2)[..., 0, :]
                        sc = np.linalg.norm(v.reshape(-1) @ Qs[nd])
                        nlb += 1
                        viol_lb += lb > sc + 1e-9
                        lbs.append(lb)
                    # B.3: certificado del primer paso
                    al = alpha_tables(srcs)
                    hitmask = np.zeros((NJ, len(TAUS)), bool)
                    for k, nd in enumerate((a, b)):
                        near = np.linalg.norm(NODES[:, :2] - NODES[nd, :2], axis=1) <= 5.0
                        hitmask |= near[:, None] & (np.abs(TAUS - taus[k]) <= 1.0)[None, :]
                    cheap = leak_cheap(al[0], srcs[0][1]) + leak_cheap(al[1], srcs[1][1])
                    cheap[hitmask] = -np.inf
                    LB = max(lbs)
                    cert = False
                    if LB > NW:
                        flat = np.argsort(cheap, axis=None)[::-1]
                        ub = -np.inf
                        for r_, f in enumerate(flat[:40]):
                            j, it = np.unravel_index(f, cheap.shape)
                            ex = sum(leak_exact(j, al[k][j, it], srcs[k][1]) for k in range(2))
                            ub = max(ub, ex)
                            # B.2: fuga sin ruido en ese nodo equivocado
                            st = np.arange(NS)
                            Ys = Yk[0] + Yk[1]
                            seg = Ys[st[None, :, None], :, WIN[j, it]]
                            ee = (seg ** 2).sum(-1)
                            v = np.take_along_axis(seg, ee.argmax(-1)[..., None, None], axis=2)[..., 0, :]
                            nleak += 1
                            viol_leak += np.linalg.norm(v.reshape(-1) @ Qs[j]) > ex + 1e-9
                            nxt = cheap.flat[flat[r_ + 1]]
                            if max(ub, nxt) + NW < LB:
                                cert = True
                                break
                            if ub + NW >= LB:
                                break
                    # algoritmo desplegado
                    s1, (j1, it1, vf1, tp1) = score_all(Y, np.zeros(NJ, bool))
                    h1 = is_hit(j1, it1, (a, b), taus)
                    ok1 = len(h1) > 0
                    ncert += cert
                    ncorr1 += ok1
                    viol_cert += cert and not ok1
                    Y2 = deflate(Y, j1, vf1, tp1)
                    if ok1 and dtau >= 3.0:
                        kk = h1[0]
                        # residuo de la fuente hallada en sus arribos (sin ruido ni la otra fuente)
                        R = deflate(Yk[kk], j1, vf1, tp1)
                        g = srcs[kk][1]
                        fr = []
                        for bb in range(2):
                            for i in range(NS):
                                tt = int(round(srcs[kk][0][bb, i] / DT))
                                lo, hi = max(tt - 8, 0), min(tt + 9, NT)
                                fr.append(np.linalg.norm(R[i][:, lo:hi], axis=0).max()
                                          / np.linalg.norm(g[bb, i]))
                        left.append((np.median(fr), np.max(fr)))
                    excl = np.linalg.norm(NODES[:, :2] - NODES[j1, :2], axis=1) <= 5.0
                    s2, (j2, it2, _v, _t) = score_all(Y2, excl)
                    h2 = is_hit(j2, it2, (a, b), taus)
                    nsucc += any(k1 != k2 for k1 in h1 for k2 in h2)
                rows.append((snr, dtau, eps, ncert / nrep, ncorr1 / nrep, nsucc / nrep))
    out.append(f"B.1  cota inferior del Lema B.2 por encima del puntaje real: {viol_lb} de {nlb}")
    out.append(f"B.2  puntaje sin ruido de nodos equivocados por encima de la fuga del Lema B.3: "
               f"{viol_leak} de {nleak}")
    out.append(f"B.3  certificados del Teorema B.4 con seleccion equivocada: {viol_cert}")
    out.append("     SNR   dtau(s)  eps(s)   certificado   acierto paso 1   exito K=2 (2 aciertos distintos)")
    for r in rows:
        out.append(f"     {r[0]:5.0f}  {r[1]:6.1f}   {r[2]:5.2f}     {r[3]:6.2f}         {r[4]:6.2f}           {r[5]:6.2f}")
    if left:
        left = np.array(left)
        out.append(f"B.4  amplitud residual de la fuente hallada tras la deflacion desplegada (fraccion de "
                   f"||g_ib||, max en +-0.5 s del arribo): mediana sobre bloques {np.median(left[:, 0]):.3f}, "
                   f"maximo {left[:, 1].max():.3f} ({len(left)} casos)")
    # B.4 teorica: residuo de un arribo leido a desfase d con amplitud a(d) g, sin ruido
    dd = np.linspace(0, DT / 2, 6)
    tt = np.linspace(-1, 1, 20001)
    out.append("     residuo de un arribo aislado sup_t |a(t) - a(d) a(t - d)| (prediccion de la Prop. B.6):")
    out.append("       " + "  ".join(f"d = {d * 1000:4.1f} ms: {np.abs(acorr(tt) - acorr(d) * acorr(tt - d)).max():.3f}"
                             for d in dd))
    print("\n".join(out))


if __name__ == "__main__":
    main()
