"""C12 mechanism check (POST HOC, vector level, no waveforms): loss s6 - sdc of an exact double-couple signal read at
the true node, with unit Gaussian noise, with and without the generator's amplitude errors (site factor exp(N(0, sd))
per station and phase, polarisation scatter of dir_deg degrees).  Network 1290 (validation), ladder_6 and ladder_4.
    python c12_mech.py > out/c12_mech.txt"""
import numpy as np

import c12_lib as L
import gen_realistic as GR

rng = np.random.default_rng(12)
for cell in ("ladder_4", "ladder_6", "unilateral_c40"):
    net = L.Net(cell, 1290, "thin")
    E, tab = net.E, net.tab
    ns = E.ns
    for target in (30.0, 60.0):
        for site_sd, dir_deg in ((0.0, 0.0), (0.0, 10.0), (0.25, 0.0), (0.5, 0.0), (0.5, 10.0)):
            n = 3000
            g = rng.integers(0, E.N, n)
            md = np.array([L.vec6(GR.rand_dc(rng)) for _ in range(n)])
            v = np.einsum("brk,bk->br", E.G[g], md).reshape(n, 2, ns, 3)
            a0 = np.linalg.norm(v, axis=-1, keepdims=True)
            w = v / np.maximum(a0, 1e-12) + rng.normal(0.0, np.tan(np.radians(dir_deg)), size=v.shape) if dir_deg else v / np.maximum(a0, 1e-12)
            v = a0 * w / np.linalg.norm(w, axis=-1, keepdims=True)
            v = v * np.exp(rng.normal(0.0, site_sd, size=(n, 2, ns, 1))) if site_sd else v
            v = v.reshape(n, 6 * ns)
            y0 = np.einsum("br,brk->bk", v, E.Q[g])
            v = v * np.sqrt(2 * target / (y0 * y0).sum(-1))[:, None]              # noise-free s6 = target
            y = np.einsum("br,brk->bk", v + rng.normal(size=v.shape), E.Q[g])
            s6 = (y * y).sum(-1) / 2
            s5 = L.s5_of(y, tab.nvec[g])
            sd = np.concatenate([L.dc_max(tab.Rm[g[i:i + 1000]], tab.Ri[g[i:i + 1000]], y[i:i + 1000]) for i in range(0, n, 1000)])
            print(f"{cell:15s} noise-free s6 {target:4.0f}  site sd {site_sd:4.2f} dir {dir_deg:4.1f} deg:  mean s6 {s6.mean():6.2f}  "
                  f"loss s6 - s5 {np.mean(s6 - s5):5.2f}  loss s6 - sdc mean {np.mean(s6 - sd):5.2f} median {np.median(s6 - sd):5.2f}  "
                  f"sdc / s6 median {np.median(sd / s6):.3f}", flush=True)
