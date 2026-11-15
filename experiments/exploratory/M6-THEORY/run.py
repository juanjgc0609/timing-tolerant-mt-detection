"""M6-THEORY: numerical checks of the theory of the M6 resolvability test ("one against two", ttmt.v2.m6).

Checks Part E of the Phase 5 theory (teoria-fase5/teoria-m6-prueba.md) on the M6-JOINT simulation:

  E.1  units: the test's RSS is the matched-filter (MF) norm, not the whitened norm of Appendix AB. Spectral
       constants of the LEGACY kernel: A2 = sum_k a(k)^2, first- and second-order factors kappa1, kappa2 that turn
       whitened margins into MF margins, and the noise weight lambda of one amplitude direction.
  E.2  null: with one event, at fixed times the gain is eps^T (P2 - P1) eps, eps ~ N(0, C) (C = MF noise covariance).
       Its mean tr(C (P2 - P1)) is negative (the free one-event class has 6 ns amplitude dofs, the pair 12) and its
       tail has the Laurent-Massart bound of Proposition E.2. Compared with the gains of the code on the M6-JOINT
       control scenes (same seeds as M6-JOINT).
  E.3  power: test-only operating characteristic. Two sources at the same lattice node, origin separation dtau,
       atoms started at the truth. Noise-free gain of the code against the MF-norm margin delta_TT,MF^2 (Lemma A.4
       with the autocorrelation of the response), then acceptance at SNR 12 against the Gaussian prediction.
  E.5  the near 0.5 s cell of M6-JOINT: hit scoring with the greedy matching of M6-JOINT against an optimal
       (bipartite) matching on the same atoms.

    python experiments/exploratory/M6-THEORY/run.py [NCTRL] [NREP] > experiments/exploratory/M6-THEORY/run_out.txt
"""
import importlib.util
import itertools
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import norm

HERE = Path(__file__).resolve().parent


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


J = _load("m6joint", HERE.parent / "M6-JOINT" / "run.py")          # scenes, source(), scene(), hits() of M6-JOINT
S = J.S                                                             # helpers of M6C-SUBSAMPLE

from ttmt.kernel import LEGACY as K  # noqa: E402
from ttmt.radiation import _rand_dc, node_operators  # noqa: E402
from ttmt.station_terms import OnsetDelay  # noqa: E402
from ttmt.traveltime import homogeneous  # noqa: E402
from ttmt.v2.m1 import OnsetEngine  # noqa: E402
from ttmt.v2.m6 import _place, _Resp, adcg_loop, free_one_event_rss, joint_refit  # noqa: E402

NCTRL = int(sys.argv[1]) if len(sys.argv) > 1 else 40
NREP = int(sys.argv[2]) if len(sys.argv) > 2 else 40
THR = 10.8                     # thr of M6-JOINT (z_e^2 + max(0, q95) with q95 < 0)
SEED = 20261108
NPROC = int(os.environ.get("NPROC", "4"))   # E.3 cells in parallel, one process per core (BLAS 1 thread)
_CELL = None
RESP = _Resp(K)
NS, NT = 8, S.NT


# ---------------------------------------------------------------- E.1 spectral constants
def spectral_constants():
    taps = K.taps / np.linalg.norm(K.taps)
    a = np.convolve(taps, taps[::-1])                       # normalised autocorrelation of the MF noise, a(0) = 1
    A2 = float(a @ a)
    b = np.convolve(a, a[::-1])                             # autocorrelation of the response r = a
    # spectral moments by the DTFT (exact for the sampled sequences)
    nu = np.linspace(-0.5, 0.5, 20001)
    P = np.abs(np.exp(-2j * np.pi * np.outer(nu, np.arange(len(taps)))) @ taps) ** 2   # |S(nu)|^2, mean 1
    w = 2 * np.pi * nu
    m = lambda f: float(np.trapezoid(f, nu))
    om2 = m(w ** 2 * P) / m(P)
    psi = m(w ** 4 * P) / m(P)
    om2r = m(w ** 2 * P ** 2) / m(P ** 2)
    psir = m(w ** 4 * P ** 2) / m(P ** 2)
    lam_r = m(P ** 3) / m(P ** 2)                           # r^T C r / r^T r
    lam_d = m(w ** 2 * P ** 3) / m(w ** 2 * P ** 2)         # same for the derivative direction r'
    k1 = A2 * om2r / om2
    k2 = A2 * (psir - om2r ** 2) / (psi - om2 ** 2)
    P2n = P ** 2 / m(P ** 2)
    bfun = lambda t: np.cos(2 * np.pi * np.multiply.outer(np.asarray(t, float), nu)) @ P2n * (nu[1] - nu[0])
    return dict(a=a, A2=A2, bfun=bfun, b0=float(b[len(b) // 2]), om2=om2, psi=psi, om2r=om2r,
                psir=psir, lam_r=lam_r, lam_d=lam_d, k1=k1, k2=k2, Pmax=float(P.max()))


# ---------------------------------------------------------------- designs, covariance
def C_apply(x, a):
    """MF noise covariance applied to a field (ns, 3, NT) (Toeplitz of a, per trace)."""
    return np.array([[np.convolve(x[i, c], a, mode="same") for c in range(3)] for i in range(x.shape[0])])


def tt_design(T):
    """free one-event (TT) class at fixed block times T (2, ns): 3 columns per block."""
    cols = []
    for b in range(2):
        for i in range(NS):
            for c in range(3):
                x = np.zeros((NS, 3, NT))
                e = np.zeros(3)
                e[c] = 1.0
                _place(x[i], e, T[b, i], RESP)
                cols.append(x)
    return cols


def pair_design(G, atoms):
    """two-event class at fixed times: 6 tensor columns per atom."""
    cols = []
    for at in atoms:
        Gk = G[at["gi"]].reshape(2, NS, 3, 6)
        for j in range(6):
            x = np.zeros((NS, 3, NT))
            for b in range(2):
                for i in range(NS):
                    _place(x[i], Gk[b, i, :, j], at["T"][b, i], RESP)
            cols.append(x)
    return cols


def orth(cols):
    M = np.stack([c.ravel() for c in cols], 1)
    Q, R = np.linalg.qr(M)
    keep = np.abs(np.diag(R)) > 1e-9 * np.abs(np.diag(R)).max()
    return Q[:, keep]


def null_eigs(Q1, Q2, a):
    """nonzero eigenvalues of C^1/2 (P2 - P1) C^1/2 (= those of A_w C_w on W = span(Q1, Q2))."""
    W = orth([q.reshape(NS, 3, NT) for q in np.concatenate([Q1, Q2], 1).T])
    CW = np.stack([C_apply(w.reshape(NS, 3, NT), a).ravel() for w in W.T], 1)
    Cw = W.T @ CW
    A = (W.T @ Q2) @ (Q2.T @ W) - (W.T @ Q1) @ (Q1.T @ W)
    L = np.linalg.cholesky(Cw + 1e-12 * np.eye(len(Cw)))
    return np.linalg.eigvalsh(L.T @ A @ L)


def lm_bound(lam, t):
    """Proposition E.2: P(gain >= m0 + 2 sqrt(t sum lp^2) + 2 t max lp + 2 sqrt(t sum ln^2)) <= 2 e^-t."""
    lp, ln = lam[lam > 0], -lam[lam < 0]
    return lam.sum() + 2 * np.sqrt(t * (lp ** 2).sum()) + 2 * t * (lp.max() if len(lp) else 0) \
        + 2 * np.sqrt(t * (ln ** 2).sum())


# ---------------------------------------------------------------- E.3 MF-norm margin (Lemma A.4 with b)
def margin_mf(g1, g2, t1, t2, sc):
    """sum over blocks of the MF-norm residual of the best single pulse (Lemma A.4 with the response's
    autocorrelation b; blocks treated as orthogonal, as in Corollary A.6)."""
    b0, bf = sc["b0"], sc["bfun"]          # exact b(t)/b(0) at fractional lags (spectral formula)
    tot = 0.0
    eps = np.linspace(-8, 8, 1601)
    for bb in range(2):
        for i in range(NS):
            u, v, D = g1[bb, i], g2[bb, i], t2[bb, i] - t1[bb, i]
            E = b0 * (u @ u + v @ v + 2 * (u @ v) * float(bf(D)))
            f = np.linalg.norm(np.outer(bf(eps), u) + np.outer(bf(eps - D), v), axis=1) ** 2 * b0
            j = int(np.argmax(f))                                   # refine the maximum on a fine local grid
            e2 = np.linspace(eps[max(j - 1, 0)], eps[min(j + 1, len(eps) - 1)], 201)
            f2 = np.linalg.norm(np.outer(bf(e2), u) + np.outer(bf(e2 - D), v), axis=1) ** 2 * b0
            tot += E - max(f.max(), f2.max())
    return tot


def _run_cell(cd):
    return _CELL(cd)


def main():
    t0 = time.time()
    print(f"M6-THEORY  seed {SEED}  numpy {np.__version__}  NCTRL {NCTRL}  NREP {NREP}\n")
    sc = spectral_constants()
    print("E.1 spectral constants of the LEGACY kernel (per sample; Ricker 2 Hz at 16 Hz)")
    print(f"  A2 = sum a(k)^2 = {sc['A2']:.3f} (closed form 0.4114/(fc dt) = {0.4114 / (2.0 * K.dt):.3f})")
    print(f"  omega^2 = {sc['om2']:.4f} /sample^2 = {sc['om2'] / K.dt ** 2:.1f} s^-2 (5 pi^2 fc^2 = "
          f"{5 * np.pi ** 2 * 4:.1f}); omega_r^2/omega^2 = {sc['om2r'] / sc['om2']:.4f} (0.9)")
    print(f"  (psi_r - omega_r^4)/(psi - omega^4) = {(sc['psir'] - sc['om2r'] ** 2) / (sc['psi'] - sc['om2'] ** 2):.4f}"
          f" (0.45)")
    print(f"  kappa1 = A2 omega_r^2/omega^2 = {sc['k1']:.3f}; kappa2 = {sc['k2']:.3f}")
    print(f"  lambda_r = r'Cr/r'r = {sc['lam_r']:.3f}; lambda for r' = {sc['lam_d']:.3f}; max |S|^2 = {sc['Pmax']:.3f}\n")

    STA, XN = S.network(), S.nodes()
    taus = np.round(np.arange(4.0, 12.0 + 1e-9, 0.1), 10)
    tt, _ = homogeneous(XN, STA, S.VP, S.VS)
    Gs, Qs = node_operators(XN, STA)
    E = OnsetEngine(XN, STA, tt, taus, OnsetDelay(np.full((8, 2), -S.LAG)), Gs, Qs, rounding="nearest")
    lag = np.full((8, 2), S.LAG / S.DT)
    zsd = float(np.median(np.abs(J.noise_mf(np.random.default_rng(0), 8))) / 0.6745)
    a = sc["a"]

    # ------------------------------------------------------------ E.2 null
    print("E.2 null: one source + noise (M6-JOINT control seeds); gain of the code vs fixed-time linear theory")
    rows, gains_all, mc_all = [], [], []
    rng = np.random.default_rng([SEED, 9])
    for c in range(NCTRL):
        r = np.random.default_rng([SEED, 1, c])
        srcs = [J.source(r)]
        Y, q, _F = J.scene(r, STA, srcs, zsd)
        at = adcg_loop(E, Y, q, 2, lag, thr=-np.inf)
        if len(at) < 2 or at[1]["pair_gain"] is None:
            continue
        g = at[1]["pair_gain"]
        Q1 = orth(tt_design(at[0]["T"]))
        Q2 = orth(pair_design(E.G, at))
        lam = null_eigs(Q1, Q2, a)
        # fixed-time Monte Carlo with the code's noise (raw white noise through the MF, divided by zsd)
        mc = []
        for _ in range(50):
            e = (J.noise_mf(rng, NS) / zsd).ravel()
            mc.append(float(((Q2.T @ e) ** 2).sum() - ((Q1.T @ e) ** 2).sum()))
        rows.append((g, lam.sum(), np.mean(mc), lm_bound(lam, np.log(2 / 0.05)), Q1.shape[1], Q2.shape[1],
                     (lam > 0).sum()))
        gains_all.append(g)
        mc_all += mc
    R = np.array([x[:4] for x in rows])
    print(f"  scenes with a tested spurious atom: {len(rows)} of {NCTRL}; dims TT {rows[0][4]}, pair {rows[0][5]}, "
          f"positive eigenvalues {rows[0][6]}")
    print(f"  code gain     median {np.median(R[:, 0]):8.1f}  q95 {np.quantile(R[:, 0], .95):8.1f}  max {R[:, 0].max():8.1f}")
    print(f"  theory mean   median {np.median(R[:, 1]):8.1f}  (tr C(P2-P1)); MC fixed-time mean median "
          f"{np.median(R[:, 2]):8.1f}")
    print(f"  fixed-time MC pooled q95 {np.quantile(mc_all, .95):8.1f}  max {np.max(mc_all):8.1f}")
    print(f"  E.2 bound (level 0.05) median {np.median(R[:, 3]):8.1f}  max {R[:, 3].max():8.1f}")
    d = R[:, 0] - R[:, 1]
    print(f"  search inflation (code gain - theory mean): median {np.median(d):6.1f}, q95 {np.quantile(d, .95):6.1f}")
    print(f"  corr(code gain, theory mean) = {np.corrcoef(R[:, 0], R[:, 1])[0, 1]:.2f}\n")
    m0 = float(np.median(R[:, 0]))          # offset under H0 of the code's gain (search included)
    s_null = float(np.std(R[:, 0]))

    # ------------------------------------------------------------ E.3 power, test only
    print("E.3 power of the test alone: two sources at one lattice node, atoms started at the truth, SNR 12")
    print(f"  {'case':5s} {'dtau':>5s} | {'margin MF':>9s} {'code g0':>8s} {'delta_TT^2 (wh.)':>16s} | "
          f"{'accept':>6s} {'pred':>5s} {'gain med':>8s}")
    grids = {"diff": (0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10), "same": (0.06, 0.08, 0.09, 0.10, 0.12, 0.15, 0.25)}
    def cell(cd):
        case, dtau = cd
        acc, gs, g0s, mfs, whs, preds = 0, [], [], [], [], []
        for rep in range(NREP):
            r = np.random.default_rng([SEED, 5, int(1000 * dtau), len(case), rep])
            gi = int(r.integers(len(XN)))
            while np.linalg.norm(XN[gi, :2]) > 9:
                gi = int(r.integers(len(XN)))
            m1 = _rand_dc(r).ravel()[[0, 4, 8, 1, 2, 5]]
            m2 = m1.copy() if case == "same" else _rand_dc(r).ravel()[[0, 4, 8, 1, 2, 5]]
            s1 = [XN[gi].copy(), 6.0, m1, np.zeros((2, 8))]
            s2 = [XN[gi].copy(), 6.0 + dtau, m2, np.zeros((2, 8))]
            F = []
            for s in (s1, s2):
                MF, g, tarr = S.field(STA, *s)
                pk = np.linalg.norm(MF, axis=1).max(axis=-1)
                s[2] = s[2] * (J.SNR * zsd / np.median(np.sort(pk)[-6:]))
                F.append(S.field(STA, *s))
            sig = sum(f[0] for f in F) / zsd
            gv = [f[1] * np.linalg.norm(K.taps) ** 2 / zsd for f in F]      # amplitudes in Y units
            tt_ = [f[2] / K.dt for f in F]                                     # peak times (samples)
            mfs.append(margin_mf(gv[0], gv[1], tt_[0], tt_[1], sc))
            whs.append(mfs[-1] / sc["k1" if case == "diff" else "k2"])

            def test(Y):
                atoms = [{"gi": gi, "T": tt_[k].copy(), "anchor": tt_[k].copy(), "bound": np.full((2, NS), 2.0),
                          "m": np.zeros(6)} for k in range(2)]
                atoms, _R, rss2 = joint_refit(Y, atoms, E.G, RESP, 3)
                one = [{"gi": gi, "T": tt_[0].copy(), "anchor": tt_[0].copy(), "bound": np.full((2, NS), 2.0),
                        "m": np.zeros(6)}]
                joint_refit(Y, one, E.G, RESP, 3)          # the accepted atom as the loop holds it (fit alone)
                return free_one_event_rss(Y, [], one[0], E.G, RESP) - rss2
            g0 = test(sig)
            gn = test(sig + J.noise_mf(r, NS) / zsd)
            g0s.append(g0)
            gs.append(gn)
            acc += gn > THR
            sd = np.sqrt(4 * sc["lam_d"] * max(g0, 0) + s_null ** 2)
            preds.append(norm.sf((THR - g0 - m0) / sd))
        return (f"  {case:5s} {dtau:5.2f} | {np.median(mfs):9.1f} {np.median(g0s):8.1f} {np.median(whs):16.1f} | "
              f"{acc / NREP:6.2f} {np.mean(preds):5.2f} {np.median(gs):8.1f}")
    global _CELL
    _CELL = cell
    cells = [(c, d) for c in ("diff", "same") for d in grids[c]]
    with mp.get_context("fork").Pool(NPROC) as pool:
        for line in pool.map(_run_cell, cells):
            print(line, flush=True)
    print(f"  pred = mean over scenes of P(N(g0 + m0, 4 lambda' g0 + s_null^2) > thr), m0 = median code gain under H0 "
          f"{m0:.1f}, s_null {s_null:.1f}\n")

    # ------------------------------------------------------------ E.5 the near 0.5 s cell
    print("E.5 near 0.5 s cell of M6-JOINT (same seeds): greedy matching vs optimal matching of atoms to sources")

    def hits_opt(atoms, srcs):
        best = 0
        for perm in itertools.permutations(range(len(atoms)), min(len(atoms), len(srcs))):
            best = max(best, sum(1 for k, j in enumerate(perm) if np.linalg.norm(atoms[j]["x"][:2] - srcs[k][0][:2])
                                 <= 5.0 and abs(atoms[j]["tau"] - srcs[k][1]) <= 0.5))
        return best
    for dtau in (0.25, 0.5):
        n_g = n_o = 0
        for rep in range(20):
            r = np.random.default_rng([SEED, 3, int(100 * dtau), len("near"), rep])
            s1 = J.source(r)
            s2 = J.source(r, s1[0], J.STRATA["near"], 6.0 + dtau)
            srcs = [s1, s2]
            Y, q, _F = J.scene(r, STA, srcs, zsd)
            at = adcg_loop(E, Y, q, 2, lag, thr=None)
            n_g += len(J.hits(at, srcs)) == 2
            n_o += hits_opt(at, srcs) == 2
        print(f"  near {dtau:4.2f}: m6-none pair recall greedy matching {n_g / 20:.2f}, optimal {n_o / 20:.2f}")
    print(f"\n{time.time() - t0:.0f} s")


if __name__ == "__main__":
    sys.exit(main())
