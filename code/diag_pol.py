"""Diagnostico del canal de polarizacion sintetico: cuanto separa |cos| arribos verdaderos de clutter."""
import numpy as np
import task0_common as T
from task0_common import make_network, generate, polarization, DT, NT
S1 = np.array([6.0, -8.0, 9.0])

def pd_at(D, t, W):
    x = D; k = int(np.clip(round(t / DT), 0, NT - 1)); h = max(1, int(W / DT) // 2)
    seg = x[:, max(0, k - h):k + h]
    C = seg @ seg.T
    w, v = np.linalg.eigh(C); return v[:, -1], (w[-1] - w[-2]) / max(w[-1], 1e-12)

def auc(a, b):
    a = np.asarray(a); b = np.asarray(b)
    return float((a[:, None] > b[None, :]).mean() + 0.5 * (a[:, None] == b[None, :]).mean())

rng = np.random.default_rng(1)
for W in (0.6, 1.2, 2.0):
    ct, cc, azok = [], [], []
    for k in range(40):
        S = make_network("random", 8, seed=100 + k)
        D, picks, tr = generate(S, [S1], [6.0], n_clutter=35, rng=rng)
        ns = len(S)
        for i in range(ns):
            ntrue = 2  # los dos primeros picks de cada estacion son P y S verdaderos
            for j, (t, b) in enumerate(picks[i]):
                v, _ = pd_at(D[i], t, W)
                if j < ntrue:
                    p = polarization(S1, S[i], b); ct.append(abs(v @ p))
                    if b == "P":
                        hp = np.arctan2(p[1], p[0]); hv = np.arctan2(v[1], v[0])
                        d = abs((np.degrees(hp - hv) + 90) % 180 - 90); azok.append(d <= 20)
                else:
                    # clutter visto como si lo explicara la fuente verdadera
                    p = polarization(S1, S[i], b); cc.append(abs(v @ p))
    print(f"ventana {W:.1f} s: |cos| verdaderos media {np.mean(ct):.3f} mediana {np.median(ct):.3f}; "
          f"clutter media {np.mean(cc):.3f}; AUC {auc(ct, cc):.3f}; P azimut<=20 deg {np.mean(azok):.2f}")
