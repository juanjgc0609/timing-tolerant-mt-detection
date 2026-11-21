"""check_rederive_bcd.py - independent re-derivation checks of Parts AB.B, AB.C and AB.D (21 Nov 2026 revision).

Written from the statements, not from the original scripts (only lib5's closed-form Ricker autocorrelation is reused
and itself re-checked against a numerical autocorrelation). Seed 20261120. numpy only.
"""
import numpy as np
from lib5 import acorr, FC, DT, OMEGA2, PSI

rng = np.random.default_rng(20261120)
out = []
P = lambda *a: out.append(" ".join(str(x) for x in a))

# ---- Ricker pulse and its autocorrelation, numerically -------------------------------------------------------
dt = 1e-4
t = np.arange(-3, 3, dt)
u = (np.pi * FC * t) ** 2
s = (1 - 2 * u) * np.exp(-u)
s /= np.sqrt((s * s).sum() * dt)
lags = np.array([0.0, 0.031, 0.1, 0.216, 0.5])
num = [np.sum(s * np.interp(t - L, t, s, left=0, right=0)) * dt for L in lags]
P("[R] a(t) closed form vs numeric at", lags.tolist(), "max abs diff", f"{np.max(np.abs(np.array(num) - acorr(lags))):.2e}")
ds = np.gradient(s, dt)
dds = np.gradient(ds, dt)
P("[R] omega2 closed", f"{OMEGA2:.2f}", "numeric", f"{(ds*ds).sum()*dt:.2f}", "| psi-omega^4 closed",
  f"{PSI-OMEGA2**2:.0f}", "numeric", f"{(dds*dds).sum()*dt-((ds*ds).sum()*dt)**2:.0f}")

# ---- B.1 constants ---------------------------------------------------------------------------------------------
tt = np.linspace(0, 1, 200001)
a = acorr(tt)
i0 = np.argmax(a < 0)
P("[B.1] first zero of a:", f"{tt[i0]:.4f} s (claim 0.118)", "| min a:", f"{a.min():.4f} at {tt[a.argmin()]:.4f} s (claim -0.618 at 0.216)")
P("[B.1] a(Delta/2) =", f"{acorr(DT/2):.4f} (claim 0.906)", "-> Z/||g|| bound (a(D/2)-a^-)/2 =",
  f"{(acorr(DT/2) + a.min())/2:.4f} (claim 0.144)")
# largest positive side lobe (must be < 0.618 for the lemma's 'main lobe' conclusion)
P("[B.1] largest positive side lobe beyond first zero:", f"{a[i0:].max():.4f}")

# Monte Carlo of Lemma B.1: one block, offset delta in [0, D/2], interference Z (noise + 2nd arrival)
viol = 0; ncase = 0; jump = 0; njump = 0
for _ in range(20000):
    g = rng.normal(size=3); g *= rng.uniform(1, 10) / np.linalg.norm(g)
    t0 = rng.uniform(-DT/2, DT/2)
    ts = np.arange(-8, 9) * DT                                  # window of 17 samples centred on nearest sample
    g2 = rng.normal(size=3) * rng.uniform(0, 1.0) * np.linalg.norm(g)
    D2 = rng.uniform(-0.4, 0.4)
    zeta = rng.normal(size=(len(ts), 3)) * rng.uniform(0, 0.3)
    interf = g2[None] * acorr(ts - t0 - D2)[:, None] + zeta
    y = g[None] * acorr(ts - t0)[:, None] + interf
    Z = np.linalg.norm(interf, axis=1).max()
    k = np.argmax(np.linalg.norm(y, axis=1))
    lam = acorr(np.min(np.abs(ts - t0))) - 2 * Z / np.linalg.norm(g)
    if lam > -a.min():
        ncase += 1
        if acorr(ts[k] - t0) < lam - 1e-12:
            viol += 1
    else:
        njump += 1
        jump += acorr(ts[k] - t0) < 0
P("[B.1] MC: cases meeting the condition", ncase, "violations", viol,
  "| cases outside the condition", njump, "of which read on a negative lobe", int(jump))

# ---- B.6: c(delta) and the deflation residue -------------------------------------------------------------------
xs = np.linspace(-1, 1, 40001)
for dms in (6.25, 12.5, 18.75, 25, 31.25):
    d = dms / 1000
    c = np.max(np.abs(acorr(xs) - acorr(d) * acorr(xs - d)))
    P(f"[B.6] c({dms} ms) = {c:.3f}", "(claims 0.077 0.154 0.232 0.310 0.387)")
P("[B.6] bound (1-a(D/2))/2 =", f"{(1-acorr(DT/2))/2:.4f}")
# residue: subtract a(delta) g a(t - t0 - delta) (read at the sample) vs refined sub-sample subtraction
d = DT / 2
res_grid = np.max(np.abs(acorr(xs) - acorr(d) * acorr(xs - d)))
dh = 0.002                                                       # refined time with 2 ms error, amplitude refit
amp = np.sum(acorr(xs) * acorr(xs - dh)) / np.sum(acorr(xs - dh) ** 2)
res_ref = np.max(np.abs(acorr(xs) - amp * acorr(xs - dh)))
P(f"[B.6] worst-offset residue: read-sample subtraction {res_grid:.3f}, sub-sample (2 ms error, refit) {res_ref:.3f}")

# ---- C.2 subgaussian constant: sinh(y)/y <= exp(y^2/6) ---------------------------------------------------------
y = np.linspace(1e-6, 50, 100000)
P("[C.2] max of log(sinh y / y) - y^2/6 over (0,50]:", f"{np.max(np.log(np.sinh(y[y<700])/y[y<700]) - y[y<700]**2/6):.2e} (must be <= 0)")
# conditional mean: E[sigma^2 | norms] = tr M / 3 , MC
nb = 14
Q, _ = np.linalg.qr(rng.normal(size=(3 * nb, 6)))
rho = rng.lognormal(0, 1, nb)
M = sum(rho[k] ** 2 * Q[3*k:3*k+3].T @ Q[3*k:3*k+3] for k in range(nb))
th = rng.normal(size=(200000, nb, 3)); th /= np.linalg.norm(th, axis=2, keepdims=True)
sig2 = np.sum((np.einsum("nkc,kcj->nj", th * rho[None, :, None], Q.reshape(nb, 3, 6))) ** 2, 1)
P("[C.2] E[sigma^2|norms]: MC", f"{sig2.mean():.4f}", "tr M/3", f"{np.trace(M)/3:.4f}")
for p in (1e-2, 1e-3):
    tl = np.log(1 / p)
    thr = (np.trace(M) + 2 * np.sqrt(np.trace(M @ M) * tl) + 2 * np.linalg.norm(M, 2) * tl) / 3
    P(f"[C.2] p={p}: exceedance {np.mean(sig2 > thr):.5f} (must be <= p); threshold/q_emp(1-p) = "
      f"{np.sqrt(thr/np.quantile(sig2, 1-p)):.3f}")

# ---- C.4(a) Lipschitz constant of ACAP, worst ratio on random pairs --------------------------------------------
def acap(v, c1=3.0):
    n = np.linalg.norm(v, axis=1)
    q = np.sort(n)[int(np.ceil(len(n) / 2)) - 1]
    return v * np.minimum(1, c1 * q / np.maximum(n, 1e-300))[:, None]
worst = 0
for _ in range(20000):
    v = rng.standard_t(2, size=(nb, 3)) * rng.lognormal(0, 1, (nb, 1))
    w = v + rng.normal(size=(nb, 3)) * 10 ** rng.uniform(-3, 1)
    sa = np.linalg.norm(Q.T @ acap(v).ravel()); sb = np.linalg.norm(Q.T @ acap(w).ravel())
    worst = max(worst, abs(sa - sb) / np.linalg.norm(v - w))
P("[C.4] worst |dsigma_A|/||dv|| on 20000 pairs:", f"{worst:.3f}", "bound 1+c1 sqrt(nb) =", f"{1+3*np.sqrt(nb):.3f}")

# ---- D.2: E[exp S_mix] = 1 (well specified, independent blocks) and what dependence does -------------------------
def tlogpdf(v, nu, b):
    from math import lgamma, log, pi
    c = lgamma((nu + 3) / 2) - lgamma(nu / 2) - 1.5 * log(nu * pi) - 3 * log(b)
    return c - (nu + 3) / 2 * np.log1p(np.sum(v * v, -1) / (nu * b * b))
nu, b, eps, sk, nb2, J = 3.0, 1.0, 0.05, 1.0, 14, 400
A = rng.normal(size=(nb2, 3, 6)) * 0.4
mj = rng.normal(size=(J, 6)) * 0.5                                # prior draws, fixed before the data
mu = np.einsum("kcj,nj->nkc", A, mj)                              # (J, nb, 3)
def emix(V):                                                      # V (R, nb, 3) -> e-values (R,)
    lf = tlogpdf(V, nu, b)                                        # (R, nb)
    lg = -1.5 * np.log(2 * np.pi * sk**2) - np.sum((V[:, None] - mu[None]) ** 2, -1) / (2 * sk**2)  # (R,J,nb)
    ell = np.log(eps + (1 - eps) * np.exp(lg - lf[:, None]))
    return np.exp(ell.sum(-1)).mean(1)
def tdraw(R, common=0.0):
    g = rng.normal(size=(R, nb2, 3))
    if common:                                                    # shared component across blocks (correlated stations)
        g = np.sqrt(1 - common) * g + np.sqrt(common) * rng.normal(size=(R, 1, 3))
    w = rng.chisquare(nu, size=(R, nb2 if not common else 1, 1)) / nu
    return b * g / np.sqrt(w)
for lab, cm in (("independent", 0.0), ("common scale per window + direction corr 0.3", 0.3), ("common scale + direction corr 0.6", 0.6)):
    E = np.concatenate([emix(tdraw(2000, cm)) for _ in range(20)])
    P(f"[D.2] {lab}: mean e = {E.mean():.3f} (se {E.std()/np.sqrt(len(E)):.3f}); P(e>=10) = {np.mean(E>=10):.4f} (<=0.1); P(e>=100) = {np.mean(E>=100):.4f} (<=0.01)")

# ---- D.4: 1-cos theta1 for colocated hypotheses vs omega^2 dtau^2 / 2 ---------------------------------------------
for dtau in (0.003, 0.01, 0.03):
    P(f"[D.4] 1-a({dtau}) = {1-acorr(dtau):.3e}, omega^2 dtau^2/2 = {OMEGA2*dtau**2/2:.3e}")

open("check_rederive_bcd_out.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
