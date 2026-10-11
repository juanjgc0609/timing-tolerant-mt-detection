"""C12 validation (network 1290, stream TAG 4305): identities of the statistics, the optimiser against a dense search,
the loop against c10_lib.mt_loop, the tube formula against Monte Carlo.   python c12_validate.py > out/c12_validation.txt"""
import time

import numpy as np
from scipy import stats

import c2_run as C2
import c10_lib as C
import c12_lib as L
import xc

rng = np.random.default_rng(5)
print("V0 curvature code on the Frobenius metric: L3, L1, mc-se(L3) =", L.lk_curvatures(L.ISO, 2000),
      " expected L3 = 4 pi^2 =", 4 * np.pi ** 2)

t0 = time.time()
net = L.Net("ladder_6", 1290, "thin")
E, tab = net.E, net.tab
print(f"setup {time.time() - t0:.1f} s; N {E.N} T {E.T} ns {E.ns}")
# V1 identities
g = rng.integers(0, E.N, 40)
v = rng.normal(size=(40, 6 * E.ns))
y = np.einsum("br,brk->bk", v, E.Q[g])
h = np.einsum("br,brk->bk", v, E.G[g])
A = np.einsum("brk,brl->bkl", E.G[g], E.G[g])
s6 = (y * y).sum(-1) / 2
s6b = np.einsum("bk,bk->b", h, np.linalg.solve(A, h[..., None])[..., 0]) / 2
print("V1a |Q^T v|^2 = h^T A^-1 h: max abs diff", np.abs(s6 - s6b).max())
Bd = np.linalg.svd(L.TRACE[None])[2][1:].T                                     # (6, 5) basis of trace-free
AB = np.einsum("ka,bkl,lc->bac", Bd, A, Bd)
hB = h @ Bd
s5b = np.einsum("ba,ba->b", hB, np.linalg.solve(AB, hB[..., None])[..., 0]) / 2
s5 = L.s5_of(y, tab.nvec[g])
print("V1b DEV5 closed form = constrained LS: max abs diff", np.abs(s5 - s5b).max())
sdc, Rr, fcl, fgr = L.dc_max(tab.Rm[g], tab.Ri[g], y, full=True)
m = L.vec6(L.dc_of_rot(Rr))
sdc_h = np.einsum("bk,bk->b", m, h) ** 2 / np.einsum("bk,bkl,bl->b", m, A, m) / 2
print("V1c sdc in Q coordinates = (m^T h)^2 / (2 m^T A m): max abs diff", np.abs(sdc - sdc_h).max(),
      "; trace of m:", np.abs(m[:, :3].sum(1)).max(), "; det:", np.abs(np.linalg.det(L.mat3(m))).max())
br = np.array([L.dc_brute(tab.Rm[k], yy, 300000, 3) for k, yy in zip(g, y)])
print("V2 optimiser vs dense search (300000 rotations), noise vectors: sdc - brute: min %.4f median %.4f max %.4f; "
      "sdc >= brute in %d/40; ordering sdc <= s5 <= s6: %s" % ((sdc - br).min(), np.median(sdc - br), (sdc - br).max(),
                                                                int((sdc >= br - 1e-9).sum()),
                                                                bool(np.all(sdc <= s5 + 1e-9) and np.all(s5 <= s6 + 1e-9))))
print("   closed-form start alone / best grid point alone / refined, mean ratio to refined: %.3f / %.3f / 1" %
      ((fcl / sdc).mean(), (fgr / sdc).mean()))
# exact double-couple signal
import gen_realistic as GR
ok = []
for k in g[:20]:
    md = L.vec6(GR.rand_dc(rng))
    vv = E.G[k] @ md * 50
    yy = vv @ E.Q[k]
    ok.append(float(L.dc_max(tab.Rm[[k]], tab.Ri[[k]], yy[None])[0] / ((yy * yy).sum() / 2)))
print("V3 exact DC signal: sdc / s6 min %.6f max %.6f" % (min(ok), max(ok)))
# V4 loop identity
for design, draw, kmax in (("STD", xc.draw_sources, 5), ("DENSE", C2.draw_dense, 7)):
    st = L.stream(net.nc, 4305, draw)
    same, n, tt = 0, 0, {"c10": 0.0, "full6": 0.0, "dev5": 0.0, "dc": 0.0}
    stat = {}
    for _ in range(6):
        w = next(st)
        t = time.time(); a0 = C.mt_loop(E, w["MFn"], kmax); tt["c10"] += time.time() - t
        res = {}
        for rk in ("full6", "dev5", "dc"):
            t = time.time(); res[rk] = L.loop(E, tab, w["MFn"], kmax, rk, stat=stat); tt[rk] += time.time() - t
        n += 1
        same += int(len(a0) == len(res["full6"]) and all(np.allclose(p["x"], q["x"]) and p["tau"] == q["tau"] and
                    abs(p["s"] - q["s6"]) < 1e-9 for p, q in zip(a0, res["full6"])))
    print(f"V4 {design}: loop(full6) == c10_lib.mt_loop in {same}/{n} windows; s per window:",
          {k: round(v / n, 3) for k, v in tt.items()}, "; DC evaluations per greedy iteration",
          round(stat["ndc"] / stat["nit"], 1))
# V5 tube formula against Monte Carlo, z ~ N(0, I6): isotropic metric and 6 nodes of the network
print("V5 single candidate, z ~ N(0, I_6): P(2 s > x); MC 100000 draws (se at 1e-3: 1e-4)")
nmc = 100000
z = rng.normal(size=(nmc, 6))
for name, Rm in [("Frobenius", L.ISO)] + [(f"node {k}", tab.Rm[k]) for k in (10, 200, 337, 500, 650)]:
    L3, L1, se = L.lk_curvatures(Rm, 4000)
    Ri = np.linalg.pinv(Rm)
    s = np.concatenate([L.dc_max(np.broadcast_to(Rm, (5000, 6, 6)), np.broadcast_to(Ri, (5000, 6, 6)), z[i:i + 5000],
                                 ngn=5, nstart=2) for i in range(0, nmc, 5000)])
    line = f"   {name:10s} L3 {L3:7.2f} (x vol S3 = {L3 / (2 * np.pi ** 2):.2f}) L1 {L1:7.2f} |"
    for x in (12.0, 16.0, 20.0, 25.0):
        mc = float(np.mean(2 * s > x))
        th = float(L.tube_tail(x, L3, L1))
        line += f" x={x:g}: MC {mc:.2e} tube {th:.2e} (chi2_4 {stats.chi2.sf(x, 4):.2e} chi2_5 {stats.chi2.sf(x, 5):.2e} chi2_6 {stats.chi2.sf(x, 6):.2e}) |"
    q = np.quantile(2 * s, [0.99, 0.999])
    line += f" q99 MC {q[0]:.2f} tube {L.tube_quantile(0.01, L3, L1):.2f}; q999 MC {q[1]:.2f} tube {L.tube_quantile(0.001, L3, L1):.2f}"
    print(line, flush=True)
