"""AUDIT-NEG: known-truth controls of the code behind the recorded negatives (card.md, declared before running).

    python experiments/exploratory/AUDIT-NEG/controls.py [m5 m2 bl gen] > experiments/exploratory/AUDIT-NEG/controls_out.txt

Synthetic data only, fixed seeds (base 20261124); no project data, no laptop. Each control prints its declared
criterion and PASS / FAIL as written in the card.
"""
import math
import sys
import time

import numpy as np

SEED = 20261124


def verdict(ok):
    return "PASS" if ok else "FAIL"


# ------------------------------------------------------------------------------------------------------------ C-M5
def m5():
    from ttmt.v2.m5 import GaussBackground, RadialBackground, evalue_mix, prior_draws
    print("== C-M5: mixture e-value with known truth (ttmt.v2.m5.evalue_mix)")
    rng = np.random.default_rng(SEED)
    H, K, J = 200_000, 16, 500
    groups = np.repeat(np.arange(K // 2), 2)                    # station = (P, S) pair
    A = rng.normal(size=(H, K, 3, 6)) * 0.3
    draws = prior_draws(rng.normal(size=(8, 6)), J, rng, amp_sd=0.5)
    kw = dict(eps=0.05, draws=draws, s0=1.0, gamma=0.3)
    v = rng.normal(size=(H, K, 3))
    fx, fl = GaussBackground(1.0), GaussBackground(0.8)

    def rates(S):
        return {p: float(np.mean(S >= math.log(1 / p))) for p in (1e-2, 1e-3)}

    Sp, Sa = evalue_mix(v, A, fx.logpdf(v), **kw), evalue_mix(v, A, fx.logpdf(v), groups=groups, **kw)
    rp, ra = rates(Sp), rates(Sa)
    me = float(np.mean(np.exp(np.minimum(Sa, 700))))
    ok = all(rp[p] <= p and ra[p] <= p for p in rp) and 0.9 <= me <= 1.1
    print(f"M5-a exact background, H {H}: product {rp}, S_avg {ra}, mean exp(S_avg) {me:.3f} "
          f"(criterion: all <= p, mean in [0.9, 1.1]) -> {verdict(ok)}")
    rb = rates(evalue_mix(v, A, fl.logpdf(v), **kw))
    print(f"M5-b too-light background (sd 0.8): product {rb} (criterion: rate at 1e-3 > 1e-3) -> "
          f"{verdict(rb[1e-3] > 1e-3)}")
    Hs = 20_000
    m = draws[rng.integers(J, size=Hs)] * 3.0
    vs = np.einsum("hkcj,hj->hkc", A[:Hs], m) + rng.normal(size=(Hs, K, 3))
    Ss = evalue_mix(vs, A[:Hs], fx.logpdf(vs), groups=groups, **kw)
    pw = float(np.mean(Ss >= math.log(1000)))
    print(f"M5-c power, signal 3 x pi, H {Hs}: P(S_avg >= log 1000) {pw:.3f}, median S_avg {np.median(Ss):.1f} "
          f"(criterion >= 0.5) -> {verdict(pw >= 0.5)}")

    def reads(n, r):
        out = np.empty((n, 3))
        for a in range(0, n, 500_000):
            b = min(n, a + 500_000)
            X = r.normal(size=(b - a, 9, 3))
            out[a:b] = X[np.arange(b - a), np.argmax((X ** 2).sum(-1), axis=1)]
        return out

    vpk = reads(H * K, rng).reshape(H, K, 3)
    for nfit in (30_000, 2_000_000):
        bg = RadialBackground.fit(reads(nfit, np.random.default_rng(SEED + nfit)), left="flat")
        lf = bg.logpdf(vpk)
        r_p = rates(evalue_mix(vpk, A, lf, **kw))
        r_a = rates(evalue_mix(vpk, A, lf, groups=groups, **kw))
        print(f"M5-d peak reads, radial background fit on {nfit} reads (alpha {bg.alpha:.2f}): product {r_p}, "
              f"S_avg {r_a}")
    print()


# ------------------------------------------------------------------------------------------------------------ C-M2
def m2():
    from ttmt.v2.m2 import GP1, Const1, matern32, time_blocks
    print("== C-M2: GP station term with a known field (ttmt.v2.m2)")
    rng = np.random.default_rng(SEED + 2)
    ne, kf = 3000, 5
    X = np.c_[rng.uniform(0, 60, ne), rng.uniform(-1.5, 1.5, ne), rng.uniform(2, 12, ne)]
    t = rng.uniform(0, 7 * 86400, ne)
    fold = time_blocks(t, kf)

    def field(sf, ell, r):
        if sf == 0:
            return np.zeros(ne)
        C = sf ** 2 * matern32(X, X, ell) + 1e-8 * np.eye(ne)
        return np.linalg.cholesky(C) @ r.normal(size=ne)

    def run(sf, noise, npairs=6):
        rows = []
        for j in range(npairs):
            r = np.random.default_rng(SEED + 100 + j)
            f = 0.3 * r.normal() + field(sf, 8.0, r)
            e = 0.15 * (r.normal(size=ne) if noise == "gauss" else r.standard_t(3, ne) / math.sqrt(3.0))
            y = f + e
            for k in range(kf):
                tr, te = fold != k, fold == k
                c = Const1().fit(X[tr], y[tr])
                gp = GP1(seed=k).fit(X[tr], y[tr])
                mc, _ = c.predict(X[te])
                mg, sg = gp.predict(X[te])
                rows.append((y[te] - mc, y[te] - mg, sg, e[te], f[te] - np.mean(f[tr])))
        ec, eg, sg, ee, fc = (np.concatenate([x[i] for x in rows]) for i in range(5))
        msec, mseg = np.mean(ec ** 2), np.mean(eg ** 2)
        reach = 1 - np.mean(ee ** 2) / msec
        return 1 - mseg / msec, float(np.mean(np.abs(eg) <= 2 * sg)), reach

    t0 = time.time()
    d, cov, reach = run(0.15, "gauss")
    ok = d >= 0.8 * reach and 0.90 <= cov <= 0.98
    print(f"M2-a field sf 0.15 ell 8 km + gaussian noise 0.15: MSE drop {d:.3f}, reachable {reach:.3f} "
          f"(ratio {d / reach:.2f}), coverage {cov:.3f} (criterion drop >= 0.8 x reachable, cov 0.90-0.98) -> {verdict(ok)}")
    d, cov, _ = run(0.0, "gauss")
    ok = abs(d) <= 0.03 and 0.90 <= cov <= 0.98
    print(f"M2-b noise only: MSE drop {d:.3f}, coverage {cov:.3f} (criterion |drop| <= 0.03, cov 0.90-0.98) -> {verdict(ok)}")
    d, cov, reach = run(0.15, "t3")
    print(f"M2-c field + t3 noise of the same variance: MSE drop {d:.3f}, reachable {reach:.3f} "
          f"(ratio {d / reach:.2f}), coverage {cov:.3f} (descriptive; real data ratio about 0.6)")
    print(f"   {time.time() - t0:.0f} s\n")


# ------------------------------------------------------------------------------------------------------------ C-BL
def bl():
    from ttmt.radiation import basis, radiation_operator
    print("== C-BL: whitening, per-station sd with an unweighted subspace (pattern of PR #68)")
    rng = np.random.default_rng(SEED + 3)
    fa, fb, fc = [], [], []
    for _ in range(200):
        ns = 7
        STA = np.c_[rng.uniform(-30, 30, (ns, 2)), np.zeros(ns)]
        src = np.r_[rng.uniform(-10, 10, 2), rng.uniform(3, 12)]
        G = radiation_operator(src, STA)                            # (2 ns 3, 6), rows (phase, station, comp)
        v0 = G @ rng.normal(size=6)
        sd = np.exp(rng.uniform(0, math.log(7.6), ns))
        sd[np.argmin(sd)], sd[np.argmax(sd)] = 1.0, 7.6
        w = np.tile(np.repeat(1 / sd, 3), 2)
        f = lambda Q, v: float(np.sum((Q.T @ v) ** 2) / np.sum(v ** 2))   # noqa: E731
        Q0 = basis(G)
        fa.append(f(Q0, v0 / 2.0))                                  # equal sd
        fb.append(f(basis(w[:, None] * G), w * v0))                 # fix: Q of diag(1/sd) G
        fc.append(f(Q0, w * v0))                                    # defect: data weighted, Q not
    print(f"BL-a equal sd, unweighted Q: min f {min(fa):.4f} (criterion >= 0.99) -> {verdict(min(fa) >= 0.99)}")
    print(f"BL-b unequal sd, Q of diag(1/sd) G: min f {min(fb):.4f} (criterion >= 0.99) -> {verdict(min(fb) >= 0.99)}")
    q = np.quantile(fc, [0.25, 0.5, 0.75])
    print(f"BL-c unequal sd (max/min 7.6), unweighted Q: median f {q[1]:.3f} [{q[0]:.3f}, {q[2]:.3f}] "
          f"(descriptive; N-1 gcheck without noise 0.705 Hi-net, 0.771 Ridgecrest)\n")


# ------------------------------------------------------------------------------------------------------------ C-GEN
def gen():
    from ttmt.synth.inject import _shift
    print("== C-GEN: fractional delay of the generator (ttmt.synth.inject._shift)")
    fs, f0, n = 100.0, 4.0, 2048
    tt = np.arange(n) / fs

    def ricker(t):
        a = (np.pi * f0 * (t - 5.0)) ** 2
        return (1 - 2 * a) * np.exp(-a)
    dt = 0.37 / fs
    x, ref = ricker(tt), ricker(tt - dt)
    y = _shift(x, fs, dt)
    lin = np.interp(tt - dt, tt, x, left=0, right=0)
    err = float(np.abs(y - ref).max() / np.abs(ref).max())
    en = float(np.sum(y ** 2) / np.sum(ref ** 2) - 1)
    print(f"GEN-a _shift (frequency domain): max error {err:.2e} of the peak, energy {en:+.2e} "
          f"(criterion both <= 1e-3) -> {verdict(err <= 1e-3 and abs(en) <= 1e-3)}")
    el = float(np.sum(lin ** 2) / np.sum(ref ** 2) - 1)
    print(f"GEN-b np.interp (linear, the M6-SPLIT pattern): max error "
          f"{np.abs(lin - ref).max() / np.abs(ref).max():.2e}, energy {el:+.2e} (criterion |energy| > 1e-3) -> "
          f"{verdict(abs(el) > 1e-3)}\n")


if __name__ == "__main__":
    which = sys.argv[1:] or ["gen", "bl", "m2", "m5"]
    t0 = time.time()
    print(f"AUDIT-NEG controls, seed {SEED}, numpy {np.__version__}\n")
    for w in which:
        globals()[w]()
    print(f"total {time.time() - t0:.0f} s")
    # (added after the run, review of 78c9bab) controls_out.txt has 9 PASS lines = 13 criteria (M5-a holds 5)
