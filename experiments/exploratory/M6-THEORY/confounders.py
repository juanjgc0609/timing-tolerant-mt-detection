"""M6-THEORY, Part E.4: confounders and source (critique CR-09, CR-10, CR-25).

  (a) S-wave splitting. The null class of a split S (two orthogonal pulses, free axes, delay <= dt_max) absorbs the
      first-order term of an S pair; check of Proposition E.4(a) on random S pairs (2-D, plane normal to the ray).
      Then the operative question for the current test, whose null class has NO splitting: does a single event
      with split S make ttmt.v2.m6 accept a false pair? M6-JOINT scenes, one source, S split by dt per station
      (random fast axis), adcg_loop with the test.
  (b) Free source time function (rank-1 bound per block): sigma_2^2 = (T - sqrt(T^2 - 4D))/2 against the SVD, and its
      first order against Theorem A.5(i).
  (c) Opposite polarity (G = 0): R = Delta^2 |g1|^2 (omega^2 - sup a'^2) + O(Delta^4) against Lemma A.4.

    python experiments/exploratory/M6-THEORY/confounders.py [NREP] > experiments/exploratory/M6-THEORY/confounders_out.txt
"""
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
NREP = int(sys.argv[1]) if len(sys.argv) > 1 else 20
sys.argv = sys.argv[:1]                       # M6-JOINT parses argv at import


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


J = _load("m6joint", HERE.parent / "M6-JOINT" / "run.py")
S = J.S
from ttmt.kernel import LEGACY as K, ricker  # noqa: E402
from ttmt.radiation import node_operators, radiation_operator  # noqa: E402
from ttmt.station_terms import OnsetDelay  # noqa: E402
from ttmt.traveltime import homogeneous  # noqa: E402
from ttmt.v2.m1 import OnsetEngine  # noqa: E402
from ttmt.v2.m6 import adcg_loop  # noqa: E402

SEED = 20261108
THR = 10.8
FC = 2.0
OM2 = 5 * np.pi ** 2 * FC ** 2
PSI = 35 * np.pi ** 4 * FC ** 4


def a_ricker(t):
    u = (np.pi * FC * t) ** 2
    return (1 - 2 * u + u ** 2 / 3) * np.exp(-u / 2)


def da_ricker(t, h=1e-6):
    return (a_ricker(t + h) - a_ricker(t - h)) / (2 * h)


EPS = np.linspace(-0.6, 0.6, 24001)


def R_lemma(g1, g2, D):
    """Lemma A.4 (whitened, continuous pulse): one free pulse against g1 s + g2 s(. - D)."""
    g1, g2 = np.atleast_1d(g1), np.atleast_1d(g2)
    E = g1 @ g1 + g2 @ g2 + 2 * (g1 @ g2) * a_ricker(D)
    f = np.linalg.norm(np.outer(a_ricker(EPS), g1) + np.outer(a_ricker(EPS - D), g2), axis=1) ** 2
    return E - f.max()


def R_split(g1, g2, D, nphi=181):
    """2-D S pair against the split class (orthogonal axes at angle phi, one free pulse per axis; delay free)."""
    best = np.inf
    for phi in np.linspace(0, np.pi / 2, nphi):
        ef, es = np.array([np.cos(phi), np.sin(phi)]), np.array([-np.sin(phi), np.cos(phi)])
        best = min(best, R_lemma(g1 @ ef, g2 @ ef, D) + R_lemma(g1 @ es, g2 @ es, D))
    return best


def part_a_theory(r):
    print("(a) S pairs (2-D) against the split class; first-order (A.5(i)) term vs split residual")
    print(f"  {'angle':>7s} {'Delta':>6s} | {'R_TT (one pulse)':>16s} {'A.5(i)':>9s} | {'R_split':>10s} "
          f"{'A.5(ii) per axis':>16s} {'slope':>6s}")
    for ang in (30, 60, 90, 120, 150):
        th = np.radians(ang)
        g1, g2 = np.array([1.0, 0.0]), 0.8 * np.array([np.cos(th), np.sin(th)])
        prev = None
        for D in (0.04, 0.02, 0.01):
            rt = R_lemma(g1, g2, D)
            G = g1 + g2
            first = OM2 * D ** 2 * ((g1[0] * g2[1] - g1[1] * g2[0]) ** 2) / (G @ G)
            rs = R_split(g1, g2, D)
            # A.5(ii)-type prediction on the best axes (scalar problems, valid for G_c != 0)
            best = np.inf
            for phi in np.linspace(0, np.pi / 2, 181):
                ef, es = np.array([np.cos(phi), np.sin(phi)]), np.array([-np.sin(phi), np.cos(phi)])
                v = sum((PSI - OM2 ** 2) / 4 * D ** 4 * ((g1 @ e) * (g2 @ e)) ** 2 / ((g1 @ e + g2 @ e) ** 2 + 1e-300)
                        for e in (ef, es))
                best = min(best, v)
            slope = np.log(prev / rs) / np.log(2) if prev else float("nan")
            prev = rs
            print(f"  {ang:5d}deg {D:6.3f} | {rt:16.3e} {first:9.3e} | {rs:10.3e} {best:16.3e} {slope:6.2f}")


def split_field(STA, src, tau, m, stat, dt_split, phis):
    """M6C-SUBSAMPLE field with the S pulse split: fast axis at angle phi in the plane normal to the ray."""
    g = (radiation_operator(src, STA) @ m).reshape(2, len(STA), 3)
    d = np.linalg.norm(STA - src, axis=1)
    tarr = np.stack([d / S.VP, d / S.VS]) + tau + stat
    t = np.arange(S.NT) * S.DT
    D = np.zeros((len(STA), 3, S.NT))
    for i in range(len(STA)):
        D[i] += g[0, i][:, None] * ricker(t - tarr[0, i], K.pulse_fc)[None]
        n = (STA[i] - src) / d[i]
        e1 = np.cross([0, 0, 1.0], n)
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(n, e1)
        ef = np.cos(phis[i]) * e1 + np.sin(phis[i]) * e2
        es = -np.sin(phis[i]) * e1 + np.cos(phis[i]) * e2
        D[i] += (g[1, i] @ ef) * ef[:, None] * ricker(t - tarr[1, i], K.pulse_fc)[None]
        D[i] += (g[1, i] @ es) * es[:, None] * ricker(t - tarr[1, i] - dt_split, K.pulse_fc)[None]
    MF = np.array([[np.convolve(D[i, c], K.taps[::-1], mode="same") for c in range(3)] for i in range(len(STA))])
    return MF, g, tarr


def part_a_sim():
    print("\n(a') single source with split S through the current test (null class without splitting); M6-JOINT setup")
    STA, XN = S.network(), S.nodes()
    taus = np.round(np.arange(4.0, 12.0 + 1e-9, 0.1), 10)
    tt, _ = homogeneous(XN, STA, S.VP, S.VS)
    Gs, Qs = node_operators(XN, STA)
    E = OnsetEngine(XN, STA, tt, taus, OnsetDelay(np.full((8, 2), -S.LAG)), Gs, Qs, rounding="nearest")
    lag = np.full((8, 2), S.LAG / S.DT)
    zsd = float(np.median(np.abs(J.noise_mf(np.random.default_rng(0), 8))) / 0.6745)
    print(f"  {'dt_split':>8s} | {'tested':>6s} {'gain median':>11s} {'gain max':>9s} | false pairs at thr {THR}")
    for dts in (0.0, 0.05, 0.1, 0.2):
        gains, fp = [], 0
        for rep in range(NREP):
            r = np.random.default_rng([SEED, 7, int(1000 * dts), rep])
            s = J.source(r)
            phis = r.uniform(0, np.pi, 8)
            MF, g, tarr = split_field(STA, *s, dts, phis)
            pk = np.linalg.norm(MF, axis=1).max(axis=-1)
            s[2] = s[2] * (J.SNR * zsd / np.median(np.sort(pk)[-6:]))
            MF, g, tarr = split_field(STA, *s, dts, phis)
            Y = MF + J.noise_mf(r, 8)
            Y = Y / float(np.median(np.abs(Y)) / 0.6745)
            q = S.synth_q([tarr], 8, r)
            at = adcg_loop(E, Y, q, 2, lag, thr=-np.inf)
            gg = [x["pair_gain"] for x in at if x["pair_gain"] is not None]
            gains += gg
            fp += sum(1 for x in gg if x > THR)
        gs = np.array(gains) if gains else np.array([np.nan])
        print(f"  {dts:8.2f} | {len(gains):6d} {np.median(gs):11.1f} {np.max(gs):9.1f} | {fp} of {NREP}")


def part_b(r):
    print("\n(b) free STF: rank-1 bound per block, formula vs SVD, and first order vs A.5(i)")
    worst = 0.0
    t = np.arange(-3, 3, 1 / 400)
    s = ricker(t, FC)
    s /= np.sqrt(np.sum(s ** 2))
    for _ in range(200):
        g1, g2 = r.normal(size=3), r.normal(size=3)
        D = r.uniform(0.005, 0.5)
        sd = np.interp(t - D, t, s, left=0, right=0)
        X = np.outer(g1, s) + np.outer(g2, sd)
        sv = np.linalg.svd(X, compute_uv=False)
        a = float(s @ sd)
        T = g1 @ g1 + g2 @ g2 + 2 * a * (g1 @ g2)
        Dt = np.sum(np.cross(g1, g2) ** 2) * (1 - a ** 2)
        s2 = (T - np.sqrt(max(T * T - 4 * Dt, 0))) / 2
        worst = max(worst, abs(s2 - sv[1] ** 2) / max(sv[1] ** 2, 1e-12))
    print(f"  200 random blocks: max relative error formula vs SVD = {worst:.1e}")
    g1, g2 = np.array([1.0, 0.2, 0.0]), np.array([0.3, 0.9, 0.1])
    G = g1 + g2
    for D in (0.04, 0.02, 0.01):
        a = a_ricker(D)
        T = g1 @ g1 + g2 @ g2 + 2 * a * (g1 @ g2)
        Dt = np.sum(np.cross(g1, g2) ** 2) * (1 - a ** 2)
        s2 = (T - np.sqrt(T * T - 4 * Dt)) / 2
        first = OM2 * D ** 2 * np.sum((g1[0] * g2[1] - g1[1] * g2[0]) ** 2) / (G @ G)
        print(f"  Delta {D:5.3f}: rank-1 residual {s2:.4e}, A.5(i) {first:.4e}, Lemma A.4 {R_lemma(g1, g2, D):.4e}, "
              f"ratio rank-1/A.5(i) {s2 / first:.4f}")
    print("  same mechanism (g2 = 0.7 g1): D = 0, rank-1 residual = 0 at every Delta (exactly)")


def part_c():
    print("\n(c) opposite polarity g2 = -g1: Lemma A.4 vs Delta^2 |g1|^2 (omega^2 - sup a'^2)")
    sup2 = float(np.max(da_ricker(EPS) ** 2))
    print(f"  omega^2 = {OM2:.1f} s^-2, sup a'^2 = {sup2:.1f} s^-2, ratio (omega^2 - sup a'^2)/omega^2 = "
          f"{(OM2 - sup2) / OM2:.4f}")
    g1 = np.array([1.0, 0.0, 0.0])
    for D in (0.04, 0.02, 0.01):
        R = R_lemma(g1, -g1, D)
        print(f"  Delta {D:5.3f}: Lemma A.4 {R:.4e}, prediction {D ** 2 * (OM2 - sup2):.4e}, ratio "
              f"{R / (D ** 2 * (OM2 - sup2)):.4f}")


def main():
    t0 = time.time()
    r = np.random.default_rng([SEED, 8])
    print(f"M6-THEORY confounders  seed {SEED}  numpy {np.__version__}  NREP {NREP}\n")
    part_a_theory(r)
    part_b(r)
    part_c()
    part_a_sim()
    print(f"\n{time.time() - t0:.0f} s")


if __name__ == "__main__":
    sys.exit(main())
