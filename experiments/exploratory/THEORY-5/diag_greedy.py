"""diag_greedy.py - diagnostico de la holgura del certificado del Teorema B.4 y del residuo de la
deflacion (Proposicion B.6), con las funciones de check_greedy.py. SNR 120, dtau 3 s, eps 0.15 s, 8 replicas.
    python3 diag_greedy.py > diag_greedy_out.txt
"""
import numpy as np

import check_greedy as cg
from lib5 import radiation, rand_dc, ttimes


def read(Y, j, it):
    st = np.arange(cg.NS)
    seg = Y[st[None, :, None], :, cg.WIN[j, it]]
    e = (seg ** 2).sum(-1)
    return np.take_along_axis(seg, e.argmax(-1)[..., None, None], axis=2)[..., 0, :].reshape(-1)


def main():
    rng = np.random.default_rng(99)
    rows, fr = [], []
    for rep in range(8):
        while True:
            a, b = rng.choice(cg.NJ, 2, replace=False)
            if np.linalg.norm(cg.NODES[a, :2] - cg.NODES[b, :2]) >= 9.0:
                break
        taus = [5.0, 8.0]
        stat = rng.uniform(-0.15, 0.15, size=(2, cg.NS))
        srcs = []
        for k, nd in enumerate((a, b)):
            g = (radiation(cg.NODES[nd], cg.STA) @ rand_dc(rng)).reshape(2, cg.NS, 3)
            g *= 120.0 / np.linalg.norm(g)
            srcs.append((ttimes(cg.NODES[nd], taus[k], cg.STA) + stat, g))
        z = cg.noise(rng)
        Yk = [cg.field([s]) for s in srcs]
        Y = Yk[0] + Yk[1] + z
        it = [int(np.argmin(np.abs(cg.TAUS - t))) for t in taus]
        k = 0 if np.linalg.norm(srcs[0][1]) >= np.linalg.norm(srcs[1][1]) else 1
        nd = (a, b)[k]
        T_true = np.linalg.norm(read(Y, nd, it[k]) @ cg.Qs[nd])
        LB = cg.lower_bound(nd, it[k], srcs[k][1], srcs[k][0], Yk[1 - k] + z)
        hit = np.zeros((cg.NJ, len(cg.TAUS)), bool)
        for kk, n2 in enumerate((a, b)):
            near = np.linalg.norm(cg.NODES[:, :2] - cg.NODES[n2, :2], axis=1) <= 5.0
            hit |= near[:, None] & (np.abs(cg.TAUS - taus[kk]) <= 1.0)[None, :]
        # puntaje real maximo entre nodos equivocados y su cota en ese nodo
        best, arg = -1, None
        for j in range(cg.NJ):
            st = np.arange(cg.NS)
            seg = Y[st[None, None, :, None], :, cg.WIN[j]]
            e = (seg ** 2).sum(-1)
            v = np.take_along_axis(seg, e.argmax(-1)[..., None, None], axis=3)[..., 0, :].reshape(len(cg.TAUS), -1)
            sc = np.linalg.norm(v @ cg.Qs[j], axis=1)
            sc[hit[j]] = -1
            if sc.max() > best:
                best, arg = sc.max(), (j, int(sc.argmax()))
        al = cg.alpha_tables(srcs)
        j, t = arg
        UBj = sum(cg.leak_exact(j, al[q][j, t], srcs[q][1]) for q in range(2)) + cg.NW if hasattr(cg, "NW") else None
        leak = sum(cg.leak_exact(j, al[q][j, t], srcs[q][1]) for q in range(2))
        dkm = np.linalg.norm(cg.NODES[j, :2] - cg.NODES[nd, :2])
        rows.append((T_true, LB, best, leak, dkm, abs(cg.TAUS[t] - taus[k])))
        # residuo de la deflacion en el nodo verdadero, sin ruido
        st = np.arange(cg.NS)
        seg = Yk[k][st[None, :, None], :, cg.WIN[nd, it[k]]]
        e = (seg ** 2).sum(-1)
        kk_ = e.argmax(-1)
        v = np.take_along_axis(seg, kk_[..., None, None], axis=2)[..., 0, :].reshape(-1)
        tpk = np.take_along_axis(cg.WIN[nd, it[k]], kk_[..., None], axis=-1)[..., 0]
        R = cg.deflate(Yk[k], nd, v, tpk)
        g = srcs[k][1]
        gmax = np.linalg.norm(g, axis=-1).max()
        for bb in range(2):
            for i in range(cg.NS):
                gi = np.linalg.norm(g[bb, i])
                if gi < 0.1 * gmax:
                    continue
                tt = int(round(srcs[k][0][bb, i] / cg.DT))
                fr.append(np.linalg.norm(R[i][:, max(tt - 8, 0):tt + 9], axis=0).max() / gi)
    rows = np.array(rows)
    print("SNR 120, dtau 3 s, eps 0.15 s; fuente mas fuerte en su nodo verdadero")
    print("  rep  T_nodo_verdadero  LB(Lema B.2)  max T equivocado  fuga Lambda en ese nodo  dist(km)  |dtau|(s)")
    for r_, r in enumerate(rows):
        print(f"  {r_:3d}  {r[0]:14.1f}  {r[1]:12.1f}  {r[2]:15.1f}  {r[3]:20.1f}  {r[4]:8.1f}  {r[5]:6.1f}")
    print(f"  medianas: LB/T_verdadero = {np.median(rows[:, 1] / rows[:, 0]):.2f}; "
          f"max T equivocado/T_verdadero = {np.median(rows[:, 2] / rows[:, 0]):.2f}; "
          f"(fuga + N_W)/T_verdadero = {np.median((rows[:, 3] + 16.13) / rows[:, 0]):.2f}")
    fr = np.array(fr)
    print(f"residuo de la deflacion desplegada en bloques con ||g_ib|| >= 0.1 max ({len(fr)} bloques): "
          f"mediana {np.median(fr):.3f}, p90 {np.quantile(fr, 0.9):.3f}, max {fr.max():.3f} de ||g_ib||")


if __name__ == "__main__":
    main()
