"""check_assumptions.py - numerical checks of the assumptions of Appendix AB (21 Nov 2026 revision: AB.A.9, AB.A.10, AB.C.2, AB.D.2, AB.D.3, AB.D.6, AB.D.4').

N1  inter-station correlated noise: Corollary C.2 conditional threshold, and time-slide p-values (Prop. D.3)
N2  data-dependent joint reading (M4) vs the rotation-equivariant energy peak: validity of Corollary C.2
N3  operator misspecification (P/S ratio, free surface, site gains, misorientation): signal fraction f in span(G),
    and the 12-D (separate P and S tensors) relaxation with its detection cost
N4  S splitting with free fast axis: order of the one-event residual of one S block (Delta^2 -> Delta^4) and the P
    share of the near-pair margin
Synthetic only; nothing here is real data. Seed 20261120. numpy only.
"""
import numpy as np
from lib5 import acorr, network, radiation, rand_dc, basis, VP, VS, AP, AS, OMEGA2

rng = np.random.default_rng(20261120)
out = []
P = lambda *a: (out.append(" ".join(str(x) for x in a)), print(out[-1]))
STA = network(8)
SRC = np.array([2.0, -3.0, 9.0])
ns = len(STA)
nb = 2 * ns


def c2_threshold(Q, rho, p):
    """Corollary C.2 threshold on sigma^2 for block weights rho (nb,), Q (3 nb, k)."""
    Qb = Q.reshape(nb, 3, -1)
    M = np.einsum("k,kcj,kcl->jl", rho ** 2, Qb, Qb)
    t = np.log(1 / p)
    return (np.trace(M) + 2 * np.sqrt(np.trace(M @ M) * t) + 2 * np.linalg.norm(M, 2) * t) / 3


# ---------------- N1a: correlated stations, Corollary C.2 -------------------------------------------------------
G = radiation(SRC, STA)
Q = basis(G)
R = 20000
for r in (0.0, 0.1, 0.3, 0.5):
    # one common 3-vector per window and phase group (zero lag), shared by all stations: marginals stay N(0, I3)
    e = rng.normal(size=(R, 2, ns, 3))
    c = rng.normal(size=(R, 2, 1, 3))
    v = (np.sqrt(1 - r) * e + np.sqrt(r) * c).reshape(R, nb, 3)
    s2 = np.sum((v.reshape(R, -1) @ Q) ** 2, 1)
    rates = []
    for p in (1e-2, 1e-3):
        thr = np.array([c2_threshold(Q, np.linalg.norm(v[j], axis=1), p) for j in range(R)])
        rates.append(np.mean(s2 > thr))
    # spatial whitening with the known equicorrelation: Sigma = (1-r) I + r 11^T per group, W = Sigma^{-1/2}
    if r > 0:
        S = (1 - r) * np.eye(ns) + r * np.ones((ns, ns))
        ev, V = np.linalg.eigh(S)
        W = V @ np.diag(ev ** -0.5) @ V.T
        Wf = np.kron(np.eye(2), np.kron(W, np.eye(3)))
        vw = (v.reshape(R, -1) @ Wf.T).reshape(R, nb, 3)
        Qw = basis(Wf @ G)
        s2w = np.sum((vw.reshape(R, -1) @ Qw) ** 2, 1)
        rw = [np.mean(s2w > np.array([c2_threshold(Qw, np.linalg.norm(vw[j], axis=1), p) for j in range(R)]))
              for p in (1e-2, 1e-3)]
    else:
        rw = rates
    P(f"[N1a] inter-station corr {r}: C.2 exceedance at p=1e-2, 1e-3: {rates[0]:.4f}, {rates[1]:.5f}"
      f" | after spatial whitening: {rw[0]:.4f}, {rw[1]:.5f}")

# ---------------- N1b: time-slides with correlated stations ------------------------------------------------------
T, B, trials = 1024, 99, 300
k = acorr(np.arange(-12, 13) / 16.0)
k /= np.linalg.norm(k)
def traces(r):
    e = rng.normal(size=(ns, T + 24))
    c = rng.normal(size=(1, T + 24))
    z = np.sqrt(1 - r) * e + np.sqrt(r) * c
    return np.array([np.convolve(x, k, "valid") for x in z])[:, :T]
def stat(y):                                     # coherent zero-moveout beam, maximum over origins
    return np.max(np.sum(y, 0) ** 2) / ns
for r in (0.0, 0.05, 0.2):
    pv = []
    for _ in range(trials):
        y = traces(r)
        m0 = stat(y)
        sh = rng.integers(0, T, size=(B, ns))
        ms = np.array([stat(np.array([np.roll(y[i], sh[b, i]) for i in range(ns)])) for b in range(B)])
        pv.append((1 + np.sum(ms >= m0)) / (B + 1))
    pv = np.array(pv)
    P(f"[N1b] inter-station corr {r}: time-slide p-values, P(p<=0.05) = {np.mean(pv<=0.05):.3f}, "
      f"P(p<=0.01) = {np.mean(pv<=0.01):.3f} (valid: 0.05, 0.01; {trials} null windows, B={B})")

# ---------------- N2: joint reading (M4-like) breaks C.2 --------------------------------------------------------
Kc, R2 = 4, 3000
def run_sel(joint):
    viol = {1e-2: 0, 1e-3: 0}
    for _ in range(R2):
        C = rng.normal(size=(nb, Kc, 3))                      # candidate readings: isotropic, independent
        nrm = np.linalg.norm(C, axis=2)
        sel = np.argmax(nrm, 1)                               # energy peak: rotation-equivariant
        if joint:                                             # coordinate ascent on ||Q^T v||^2
            for _ in range(3):
                for kk in range(nb):
                    best, bj = -1, sel[kk]
                    for j in range(Kc):
                        sel[kk] = j
                        val = np.sum((C[np.arange(nb), sel].ravel() @ Q) ** 2)
                        if val > best:
                            best, bj = val, j
                    sel[kk] = bj
        v = C[np.arange(nb), sel]
        s2 = np.sum((v.ravel() @ Q) ** 2)
        for p in viol:
            viol[p] += s2 > c2_threshold(Q, np.linalg.norm(v, axis=1), p)
    return {p: viol[p] / R2 for p in viol}
a = run_sel(False)
P(f"[N2] C.2 exceedance, energy-peak reading: p=1e-2 {a[1e-2]:.4f}, p=1e-3 {a[1e-3]:.4f}")
for Kc in (4, 16):
    b = run_sel(True)
    P(f"[N2] joint reading ({Kc} candidates/block): p=1e-2 {b[1e-2]:.4f}, p=1e-3 {b[1e-3]:.4f}")

# ---------------- N3: operator misspecification ------------------------------------------------------------------
EZ = np.array([0.0, 0.0, 1.0])                                  # lib5: z grows with depth (stations z=0, source z=9)


def fs_coeffs(p, a=VP, b=VS):
    """surface displacement (x radial, z down) per unit incident upgoing P and SV of slowness p in a Poisson
    half-space, from the traction-free conditions; complex beyond the SV critical angle."""
    lm = (a / b) ** 2 - 2                                        # lambda / mu
    ea = np.sqrt(complex(1 / a ** 2 - p ** 2)); eb = np.sqrt(complex(1 / b ** 2 - p ** 2))
    if ea.imag < 0:
        ea = -ea                                                 # evanescent reflected P decays with depth
    trac = lambda U, q: np.array([q * U[0] + p * U[1], lm * (p * U[0] + q * U[1]) + 2 * q * U[1]])
    UPd, USd = np.array([p, ea]) * a, np.array([eb, -p]) * b      # reflected P and SV (downgoing)
    A_ = np.stack([trac(UPd, ea), trac(USd, eb)], 1)
    res = {}
    for kind, Ui, qi in (("P", np.array([p, -ea]) * a, -ea), ("SV", np.array([eb, p]) * b, -eb)):
        cf = np.linalg.solve(A_, -trac(Ui, qi))
        res[kind] = Ui + cf[0] * UPd + cf[1] * USd
    return res


def station_rows(Gi, n_src, n_arr, fs):
    """P and S rows (2, 3, 6) of one station: radiation with the take-off direction n_src, polarisation at the
    arrival direction n_arr, optional free surface. Returns (in-phase rows, out-of-phase rows)."""
    hh = np.r_[n_src[:2] / np.linalg.norm(n_src[:2]), 0.0]
    sh = np.cross(EZ, hh)
    sv_src = np.cross(sh, n_src); sv_arr = np.cross(sh, n_arr)
    aP, aSH, aSV = n_src @ Gi[0], sh @ Gi[1], sv_src @ Gi[1]
    re, im = np.zeros((2, 3, 6)), np.zeros((2, 3, 6))
    if not fs:
        re[0] = np.outer(n_arr, aP); re[1] = np.outer(sh, aSH) + np.outer(sv_arr, aSV)
        return re, im
    sin_i = np.linalg.norm(n_arr[:2])
    uP = fs_coeffs(sin_i / VP)["P"]; uS = fs_coeffs(sin_i / VS)["SV"]
    # incident polarisations in (x, z): P (sin, -cos) = n_arr; SV (cos, sin) = sign-matched to sv_arr
    svx = np.array([np.sqrt(1 - sin_i ** 2), sin_i]); sgn = np.sign((hh * svx[0] + EZ * svx[1]) @ sv_arr)
    vP = hh * uP[0] + EZ * uP[1]; vS = sgn * (hh * uS[0] + EZ * uS[1])
    re[0], im[0] = np.outer(vP.real, aP), np.outer(vP.imag, aP)
    re[1] = np.outer(sh, 2 * aSH) + np.outer(vS.real, aSV); im[1] = np.outer(vS.imag, aSV)
    return re, im


def rad(src, sta, ps=AS / AP, fs=False, gain=None, mis=None, steep=None):
    """'truth' (or model) operator variants: S/P ratio ps, free surface, site gains, horizontal misorientation,
    steepened arrival (sin of incidence multiplied by `steep`, as a slower near-surface layer would do).
    Returns (G_in_phase, G_out_of_phase), each (6 ns, 6), rows P stations then S stations."""
    G0 = radiation(src, sta).reshape(2, ns, 3, 6)
    G0[1] *= ps / (AS / AP)
    re, im = np.zeros_like(G0), np.zeros_like(G0)
    for i in range(ns):
        n = (sta[i] - src) / np.linalg.norm(sta[i] - src)
        na = n
        if steep is not None:
            si = np.linalg.norm(n[:2]) * steep
            na = np.r_[n[:2] / np.linalg.norm(n[:2]) * si, np.sign(n[2]) * np.sqrt(1 - si ** 2)]
        r_, i_ = station_rows(G0[:, i], n, na, fs)
        re[:, i], im[:, i] = r_, i_
    if gain is not None:
        re *= gain[None, :, None, None]; im *= gain[None, :, None, None]
    if mis is not None:
        for i in range(ns):
            c, s_ = np.cos(mis[i]), np.sin(mis[i])
            Rm = np.array([[c, -s_, 0], [s_, c, 0], [0, 0, 1]])
            re[:, i] = np.einsum("ab,pbj->paj", Rm, re[:, i]); im[:, i] = np.einsum("ab,pbj->paj", Rm, im[:, i])
    return re.reshape(-1, 6), im.reshape(-1, 6)


def frac(Qm, Gt, m):
    g, h = Gt[0] @ m, Gt[1] @ m
    return np.sum((g @ Qm) ** 2) / (np.sum(g * g) + np.sum(h * h))


dcs = [rand_dc(rng) for _ in range(300)]
G0 = radiation(SRC, STA)
def bases(Gm):
    G12 = np.zeros((6 * ns, 12)); G12[:3 * ns, :6] = Gm[:3 * ns]; G12[3 * ns:, 6:] = Gm[3 * ns:]
    return basis(Gm), basis(G12)
Q6, Q12 = bases(G0)
Gr1d = rad(SRC, STA, steep=0.45)[0]                              # 'RAY1D-like' model: steep arrival, no free surface
Q6r, Q12r = bases(Gr1d)
cases = [
    ("exact model", dict(), (Q6, Q12)),
    ("S/P ratio 3 (model 1.5)", dict(ps=3.0), (Q6, Q12)),
    ("S/P ratio 5.2 (model 1.5)", dict(ps=5.2), (Q6, Q12)),
    ("free surface, straight-ray incidence", dict(fs=True), (Q6, Q12)),
    ("steep arrival (x0.45), straight model", dict(steep=0.45), (Q6, Q12)),
    ("steep arrival + free surface, straight model", dict(steep=0.45, fs=True), (Q6, Q12)),
    ("steep arrival + free surface, RAY1D model", dict(steep=0.45, fs=True), (Q6r, Q12r)),
    ("site gains lognormal sd 0.5", dict(gain=rng.lognormal(0, 0.5, ns)), (Q6, Q12)),
    ("misorientation 2 of 8 stations, 30 deg", dict(mis=np.r_[np.radians([30, -30]), np.zeros(ns - 2)]), (Q6, Q12)),
    ("ALL: steep+FS+S/P 5.2+gains+misor., straight", dict(steep=0.45, fs=True, ps=5.2, gain=rng.lognormal(0, 0.5, ns),
      mis=np.r_[np.radians([30, -30]), np.zeros(ns - 2)]), (Q6, Q12)),
    ("ALL, RAY1D model", dict(steep=0.45, fs=True, ps=5.2, gain=rng.lognormal(0, 0.5, ns),
      mis=np.r_[np.radians([30, -30]), np.zeros(ns - 2)]), (Q6r, Q12r)),
]
for lab, kw, (qa, qb) in cases:
    Gt = rad(SRC, STA, **kw)
    f6 = np.mean([frac(qa, Gt, m) for m in dcs]); f12 = np.mean([frac(qb, Gt, m) for m in dcs])
    P(f"[N3] {lab:46s} f(6-D) = {f6:.3f}  f(12-D, separate P and S tensors) = {f12:.3f}")
inc = [np.degrees(np.arcsin(np.linalg.norm((s - SRC)[:2]) / np.linalg.norm(s - SRC))) for s in STA]
P("[N3] straight-ray incidence angles (deg):", np.round(inc, 0).tolist(), "| SV critical angle",
  f"{np.degrees(np.arcsin(VS / VP)):.1f}")
# detection cost of the 12-D relaxation: per-hypothesis level 1e-5, power at network energy E
q6 = np.quantile(rng.chisquare(6, 4_000_000), 1 - 1e-5); q12 = np.quantile(rng.chisquare(12, 4_000_000), 1 - 1e-5)
for E in (40.0, 80.0, 154.0):
    for lab, f6, f12 in (("exact", 1.0, 1.0), ("S/P 5.2", None, None)):
        if f6 is None:
            Gt = rad(SRC, STA, ps=5.2)
            f6 = np.mean([frac(Q6, Gt, m) for m in dcs]); f12 = np.mean([frac(Q12, Gt, m) for m in dcs])
        p6 = np.mean(rng.noncentral_chisquare(6, f6 * E, 200000) > q6)
        p12 = np.mean(rng.noncentral_chisquare(12, f12 * E, 200000) > q12)
        P(f"[N3] power at level 1e-5, E={E:.0f} ({lab}): 6-D {p6:.3f}, 12-D {p12:.3f}")

# ---------------- N4: S splitting with a free fast axis ----------------------------------------------------------
eps_grid = np.linspace(-0.3, 0.3, 6001)
def r_vec(g1, g2, D):
    tot = g1 @ g1 + g2 @ g2 + 2 * (g1 @ g2) * acorr(D)
    f = np.sum((np.outer(acorr(eps_grid), g1) + np.outer(acorr(eps_grid - D), g2)) ** 2, 1)
    return tot - f.max(), eps_grid[f.argmax()]
P("[N4] one S block, data g1 s(t) + g2 s(t-D) in the plane normal to the ray")
for D in (0.005, 0.01, 0.02, 0.04):
    rt, rs, dmax = [], [], []
    for _ in range(40):
        g1, g2 = rng.normal(size=2), rng.normal(size=2)
        rt.append(r_vec(g1, g2, D)[0])
        best = np.inf
        for ph in np.linspace(0, np.pi / 2, 91):             # fast axis; components f, s each get their own time
            f_ = np.array([np.cos(ph), np.sin(ph)]); s_ = np.array([-f_[1], f_[0]])
            a1, e1 = r_vec(np.array([g1 @ f_]), np.array([g2 @ f_]), D)
            a2, e2 = r_vec(np.array([g1 @ s_]), np.array([g2 @ s_]), D)
            if a1 + a2 < best:
                best, dl = a1 + a2, abs(e1 - e2)
        rs.append(best); dmax.append(dl)
    rt, rs = np.median(rt), np.median(rs)
    P(f"[N4] D={D:.3f} s: median R_TT/D^2 = {rt/D**2:9.2f}  median R_split/D^4 = {rs/D**4:12.1f}  "
      f"R_split/R_TT = {rs/rt:.4f}  max fast-slow delay used {max(dmax):.3f} s")
# P polarisation is along the ray for every mechanism, so kappa_P = 0: the first-order (K_mec) margin of a colocated
# pair comes from S alone. Network margin delta/sqrt(E), colocated pair, different mechanisms, with and without a
# free fast axis per station (|delay| <= 0.2 s), and the SNR needed for delta >= z_e = 3.29 (error 0.1).
GA = radiation(SRC, STA).reshape(2, ns, 3, 6)
kP = []
for _ in range(50):
    m1, m2 = rand_dc(rng), rand_dc(rng)
    kP.append(max(np.sum(np.cross(GA[0, i] @ m1, GA[0, i] @ m2) ** 2) for i in range(ns)))
P(f"[N4] max over stations of |g1 x g2|^2 for P (50 random pairs): {max(kP):.2e} (P carries no first-order term)")
def split_res(g1, g2, D, dmax=0.2):
    best = np.inf
    for ph in np.linspace(0, np.pi / 2, 46):
        f_ = np.array([np.cos(ph), np.sin(ph)]); s_ = np.array([-f_[1], f_[0]])
        a1, e1 = r_vec(np.array([g1 @ f_]), np.array([g2 @ f_]), D)
        a2, e2 = r_vec(np.array([g1 @ s_]), np.array([g2 @ s_]), D)
        if abs(e1 - e2) <= dmax and a1 + a2 < best:
            best = a1 + a2
    return best
pairs = [(rand_dc(rng), rand_dc(rng)) for _ in range(12)]
for D in (0.02, 0.05, 0.1, 0.2):
    r_tt, r_sp = [], []
    for m1, m2 in pairs:
        E = sum(np.sum((GA[b_, i] @ m1) ** 2) for b_ in range(2) for i in range(ns))
        tt = sp = 0.0
        for i in range(ns):
            n = (STA[i] - SRC) / np.linalg.norm(STA[i] - SRC)
            gP1, gP2 = GA[0, i] @ m1, GA[0, i] @ m2
            rp = r_vec(gP1, gP2, D)[0]
            e1 = np.cross(n, [0, 0, 1.0]); e1 /= np.linalg.norm(e1); e2 = np.cross(n, e1)
            B2 = np.stack([e1, e2])
            gS1, gS2 = B2 @ (GA[1, i] @ m1), B2 @ (GA[1, i] @ m2)
            tt += rp + r_vec(gS1, gS2, D)[0]
            sp += rp + split_res(gS1, gS2, D)
        r_tt.append(np.sqrt(max(tt, 0) / E)); r_sp.append(np.sqrt(max(sp, 0) / E))
    P(f"[N4] colocated pair, D={D:.2f} s: median delta_TT/sqrtE = {np.median(r_tt):.4f}, with free splitting "
      f"{np.median(r_sp):.4f} (ratio {np.median(np.array(r_sp)/np.array(r_tt)):.3f}); SNR for delta>=3.29: "
      f"{3.29/np.median(r_tt):.1f} vs {3.29/max(np.median(r_sp),1e-9):.1f}")
open("check_assumptions_out.txt", "w").write("\n".join(out) + "\n")
