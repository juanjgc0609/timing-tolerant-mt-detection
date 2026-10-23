"""lib5.py - modelo minimo autocontenido para las verificaciones de la teoria de la fase 5.

Replica las convenciones del repositorio (code/task0_common.py, code/task1x_mtglrt.py):
Ricker de 2 Hz, muestreo de 16 Hz, VS = 3 km/s, VP = 1.73 VS, AP = 1, AS = 1.5, red 'spread' de
make_network, operador de radiacion de momento tensorial (m = Mxx, Myy, Mzz, Mxy, Mxz, Myz) con 1/d.
Las constantes CP y CS del repositorio solo reescalan P y S; aqui se fijan en 1 (no cambia ninguna
conclusion cualitativa y las cantidades se reportan normalizadas por la energia del evento).
"""
import numpy as np

FC = 2.0
DT = 1.0 / 16.0
VS = 3.0
VP = 1.73 * VS
AP, AS = 1.0, 1.5


def acorr(t):
    """Autocorrelacion normalizada (a(0) = 1) del Ricker de pico FC, forma cerrada:
    a(t) = (1 - 2u + u^2/3) exp(-u/2), u = (pi FC t)^2."""
    u = (np.pi * FC * np.asarray(t, float)) ** 2
    return (1 - 2 * u + u * u / 3) * np.exp(-u / 2)


def dacorr(t):
    """a'(t)."""
    t = np.asarray(t, float)
    k = (np.pi * FC) ** 2
    u = k * t * t
    du = 2 * k * t
    return ((-2 + 2 * u / 3) - 0.5 * (1 - 2 * u + u * u / 3)) * np.exp(-u / 2) * du


OMEGA2 = 5 * np.pi ** 2 * FC ** 2          # ||s'||^2 = -a''(0)
PSI = 35 * np.pi ** 4 * FC ** 4            # ||s''||^2 = a''''(0)


def network(n=8):
    a = np.linspace(0, 2 * np.pi, n, endpoint=False) + 0.3
    r = np.linspace(9, 58, n)
    return np.stack([r * np.cos(a), r * np.sin(a), np.zeros(n)], axis=1)


def radiation(src, STA):
    """A (2 ns x 3, 6): bloques P (estaciones 0..ns-1) y luego S, como radiation_operator."""
    rows = []
    for b in ("P", "S"):
        for i in range(len(STA)):
            v = STA[i] - src
            d = np.linalg.norm(v)
            n = v / d
            nx, ny, nz = n
            if b == "P":
                c = np.array([nx * nx, ny * ny, nz * nz, 2 * nx * ny, 2 * nx * nz, 2 * ny * nz])
                g = np.outer(n, c) * AP / d
            else:
                B = np.array([[nx, 0, 0, ny, nz, 0], [0, ny, 0, nx, 0, nz], [0, 0, nz, 0, nx, ny]])
                g = (np.eye(3) - np.outer(n, n)) @ B * AS / d
            rows.append(g)
    return np.vstack(rows)


def ttimes(src, tau, STA):
    """tiempos (2, ns): P y S."""
    d = np.linalg.norm(STA - src, axis=1)
    return np.stack([tau + d / VP, tau + d / VS])


def rand_dc(rng):
    s = rng.normal(size=3)
    s /= np.linalg.norm(s)
    d = rng.normal(size=3)
    d -= d.dot(s) * s
    d /= np.linalg.norm(d)
    M = np.outer(s, d) + np.outer(d, s)
    return np.array([M[0, 0], M[1, 1], M[2, 2], M[0, 1], M[0, 2], M[1, 2]])


def basis(G):
    U, s, _ = np.linalg.svd(G, full_matrices=False)
    return U[:, s > 1e-8 * s.max()]
