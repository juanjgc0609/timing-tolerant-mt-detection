"""C10 post hoc: station corrections re-estimated with the catalogue events relocated by their own arrival times
(SC-R).  Motivated by the oracle result (true statics alone give +0.05 dense / +0.10 standard).  Catalogue data only.
  corr0 = registered +SC (median of the gated peak delays at the catalogue positions).
  repeat NIT times: for every catalogue event with >= 4 gated arrivals, r = delay - corr (gated, |r - median| <= 0.5 s),
     MAP position shift d under d ~ N(0, 1 km^2 per coordinate) (the catalogue error) and residual sd sig (catalogue
     rule), origin shift free: d = (Jc^T Jc / sig^2 + I)^-1 Jc^T rc / sig^2, Jc, rc centred over the used arrivals,
     J = -n / V; corr = median over events of (delay - J d - origin shift) over the gated arrivals, then re-referenced
     so that its mean equals the mean of corr0 (the common mode is not identifiable and is absorbed by the origin).
    python c10_screloc.py check     -> out/c10_screloc_check.txt (correction error against the true statics)"""
import sys

import numpy as np

import c10_lib as C
import xc
from station_corrections import peak_delays
from task0_common import VP, VS

NIT = 3
VINV = np.array([1.0 / VP, 1.0 / VS])


def sc_reloc(STA, cat, corr0, sig, gate=25.0, nit=NIT):
    obs = []
    for D, src, tau in cat:
        src = np.asarray(src, float)
        dl, en = peak_delays(STA, D, src, float(tau))
        dv = STA - src[None]
        n = dv / np.linalg.norm(dv, axis=1, keepdims=True)
        J = -n[:, None, :] * VINV[None, :, None]                      # (ns, 2, 3)
        obs.append((dl, en >= gate, J))
    corr = corr0.copy()
    for _ in range(nit):
        acc = [[[] for _ in range(2)] for _ in range(len(STA))]
        for dl, ok, J in obs:
            r = dl - corr
            adj = np.zeros_like(dl)
            mu = 0.0
            if ok.sum() >= 2:
                mu = float(np.median(r[ok]))
            use = ok & (np.abs(r - mu) <= 0.5)
            if use.sum() >= 4:
                Jc = J[use] - J[use].mean(axis=0, keepdims=True)
                rc = r[use] - r[use].mean()
                d = np.linalg.solve(Jc.T @ Jc / sig ** 2 + np.eye(3), Jc.T @ rc / sig ** 2)
                adj = J @ d
                mu = float((r - adj)[use].mean())
            for i in range(len(STA)):
                for b in range(2):
                    if ok[i, b]:
                        acc[i][b].append(dl[i, b] - adj[i, b] - mu)
        new = np.zeros_like(corr0)
        est = np.zeros(corr0.shape, bool)
        for i in range(len(STA)):
            for b in range(2):
                if acc[i][b]:
                    new[i, b] = float(np.median(acc[i][b]))
                    est[i, b] = True
        if est.any():
            new[est] += corr0[est].mean() - new[est].mean()
        corr = new
    return corr


def err(c, statics):
    d = c - statics
    return float(np.sqrt(((d - np.median(d)) ** 2).mean()))


if __name__ == "__main__":
    lines = []
    for reg in ("thin", "strong"):
        e0, e1 = [], []
        for cell in ("ladder_4", "ladder_6", "unilateral_c40"):
            a, b = [], []
            for j in range(160, 180):
                nc = xc.NetCfg(cell, j, reg)
                sig = C.sigma_catalogue(nc.STA, nc.cat, nc.corr)[0][0]
                a.append(err(nc.corr, nc.statics))
                b.append(err(sc_reloc(nc.STA, nc.cat, nc.corr, sig), nc.statics))
            lines.append(f"{reg:6s} {cell:15s} rms correction error (s; median removed), median over 20 networks: "
                         f"+SC {np.median(a):.3f}  SC-R {np.median(b):.3f};  mean {np.mean(a):.3f} -> {np.mean(b):.3f}; "
                         f"SC-R better in {np.mean(np.array(b) < np.array(a)):.2f} of the networks")
            print(lines[-1], flush=True)
    open("out/c10_screloc_check.txt", "w").write("\n".join(lines) + "\n")
