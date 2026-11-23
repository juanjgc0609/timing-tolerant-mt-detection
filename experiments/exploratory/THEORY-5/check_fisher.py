"""check_fisher.py - numerical result AB.A.8 (Fisher information of a colocated pair, 20 parameters).

Model of lib5 (straight rays, homogeneous, Ricker 2 Hz unit energy, white noise sd 1), pair of point sources at the
same hypocentre with origins tau and tau + dtau. Parameters per event: x, y, z (km), tau (s), m (6). Jacobian by
central differences on a 400 Hz grid. Reports the smallest eigenvalue of J (normalised by the event energy E) and
the Cramer-Rao standard deviation of dtau at network SNR sqrt(E) = 12.4, against dtau. Seed 20261120.
"""
import numpy as np
from lib5 import network, radiation, rand_dc, ttimes, FC

rng = np.random.default_rng(20261120)
STA = network(8)
SRC = np.array([2.0, -3.0, 9.0])
fs = 400.0
t = np.arange(0, 22, 1 / fs)


def pulse(x):
    u = (np.pi * FC * x) ** 2
    return (1 - 2 * u) * np.exp(-u)


norm = np.sqrt(np.sum(pulse(t - 5) ** 2) / fs)


def signal(theta):
    """theta (20,): [x1 y1 z1 tau1 m1(6) x2 y2 z2 tau2 m2(6)] -> stacked waveforms (2 x 8 x 3 x nt) flattened."""
    out = np.zeros((2, 8, 3, len(t)))
    for k in range(2):
        th = theta[10 * k:10 * k + 10]
        src, tau, m = th[:3], th[3], th[4:]
        G = radiation(src, STA).reshape(2, 8, 3, 6)
        T = ttimes(src, tau, STA)
        for b in range(2):
            for i in range(8):
                out[b, i] += np.outer(G[b, i] @ m, pulse(t - T[b, i]) / norm)
    return out.ravel() / np.sqrt(fs)                      # so that sum of squares = continuous energy


def fisher(theta, h):
    cols = []
    for j in range(20):
        e = np.zeros(20); e[j] = h[j]
        cols.append((signal(theta + e) - signal(theta - e)) / (2 * h[j]))
    D = np.stack(cols, 1)
    return D.T @ D


lines = []
for lab, same in (("different mechanisms", False), ("same mechanism (m2 = 0.7 m1)", True)):
    m1 = rand_dc(rng); m2 = 0.7 * m1 if same else rand_dc(rng)
    for dtau in (0.005, 0.01, 0.02, 0.05, 0.1):
        th = np.r_[SRC, 5.0, m1, SRC, 5.0 + dtau, m2]
        E = np.sum(signal(np.r_[SRC, 5.0, m1, SRC, 50.0, 0 * m2]) ** 2)   # energy of event 1 alone
        h = np.r_[[1e-3] * 3, 1e-4, [1e-3] * 6, [1e-3] * 3, 1e-4, [1e-3] * 6]
        J = fisher(th, h) / E
        ev = np.linalg.eigvalsh(J)
        c = np.zeros(20); c[13], c[3] = 1, -1                # dtau = tau2 - tau1
        crb = np.sqrt(c @ np.linalg.solve(J, c)) / 12.4       # J was per unit E; SNR sqrt(E) = 12.4
        lines.append(f"[A.8] {lab:30s} dtau={dtau:.3f} s: lambda_min/E = {ev[0]:.3e}, lambda_2/E = {ev[1]:.3e}, "
                     f"CRB sd(dtau) at SNR 12.4 = {crb:.4f} s")
        print(lines[-1])
# log-log slopes of lambda_min
open("check_fisher_out.txt", "w").write("\n".join(lines) + "\n")
