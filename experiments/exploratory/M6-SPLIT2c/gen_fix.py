"""Exact (Fourier) delay of the slow S component for the split data generator. The M6-SPLIT* generators delay the
sampled pulse by linear interpolation (np.interp) at a sub-sample delay (0.02-0.2 s at 0.0625 s sampling), which low-
pass filters the slow component: the data then differ from any split model (noise-free check: 1.1e-3 relative RSS at
the true parameters). install(M) replaces M.S.field by the same generator with a band-limited exact shift."""
import numpy as np


def shift(w, dt_samples):
    n = w.size
    f = np.fft.rfftfreq(2 * n)
    W = np.fft.rfft(np.r_[w, np.zeros(n)])
    return np.fft.irfft(W * np.exp(-2j * np.pi * f * dt_samples), 2 * n)[:n]


def install(M):
    S, K, SPLIT, F0 = M.S, M.K, M.SPLIT, M._FIELD0
    from ttmt.radiation import radiation_operator

    def field(STA, src, tau, m, stat):
        if not SPLIT["on"]:
            return F0(STA, src, tau, m, stat)
        g = (radiation_operator(src, STA) @ m).reshape(2, len(STA), 3)
        d = np.linalg.norm(STA - src, axis=1)
        tarr = np.stack([d / S.VP, d / S.VS]) + tau + stat
        t = np.arange(S.NT) * S.DT
        D = np.zeros((len(STA), 3, S.NT))
        for i in range(len(STA)):
            D[i] += g[0, i][:, None] * S.ricker(t - tarr[0, i], K.pulse_fc)[None]
            n = (STA[i] - src) / d[i]
            e1 = np.cross([0, 0, 1.0], n)
            e1 /= np.linalg.norm(e1)
            e2 = np.cross(n, e1)
            ph, dt = SPLIT["phi"][i], SPLIT["dt"][i]
            ef = np.cos(ph) * e1 + np.sin(ph) * e2
            es = -np.sin(ph) * e1 + np.cos(ph) * e2
            w = S.ricker(t - tarr[1, i], K.pulse_fc)
            D[i] += (g[1, i] @ ef) * ef[:, None] * w[None] + (g[1, i] @ n) * n[:, None] * w[None]
            D[i] += (g[1, i] @ es) * es[:, None] * shift(w, dt / S.DT)[None]
        MF = np.array([[np.convolve(D[i, c], K.taps[::-1], mode="same") for c in range(3)] for i in range(len(STA))])
        return MF, g, tarr
    S.field = field
