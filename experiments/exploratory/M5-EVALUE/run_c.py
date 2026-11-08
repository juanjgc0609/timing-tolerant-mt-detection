"""M5-EVALUE, gate G-M5c (calibration) of the v2 specification (M5 "Decidir"; theory D.1, D.2; ADR-0010 pivot).

Question. On real Task 2M noise, with the Student t background estimated per station and noise block, is the null rate
of the mixture e-value S_mix >= log(1/p) at most p, for p in {1e-2, 1e-3}?

Design (fixed before the run; see card.md):
  networks   the Task 2M REAL-arm networks of results/task2M/cfg_SEED.json (station subsets, catalogue, station
             corrections: code/t2m_lib.py unchanged). 2718: day 101231 (transients), 3141: days 100210 + 100211.
  slots      each network uses its noise block (about 50 min of consecutive 40 s slots). The slots of the catalogue
             (net.cat_slots) are removed; the rest, in time order, alternate FIT (even position) / TEST (odd).
  reading    the registered hard reading: per station and phase, the matched-filter vector (MFn = MF / sd of the
             window, as the pipeline) at the energy peak within +-4 samples of the predicted sample.
  background per station, fit on the DISTINCT peak read vectors (duplicates removed) at the predicted samples of
             400 random lattice hypotheses per FIT slot (time slides on FIT data only), both phases pooled:
               R  (primary) ttmt.v2.m5.RadialBackground: isotropic, kernel density of log|v| below the 0.99 quantile,
                  Pareto tail above (Hill index, capped at 10: a
                  heavier tail than estimated is conservative, a lighter one is not). Declared deviation from the specification's isotropic t, made before
                  any real-data S_mix: the peak reader leaves no mass near 0, which a single t cannot fit (card.md).
               T  isotropic Student t by ML (the specification's background; declared ablation).
               G  isotropic gaussian (ablation; D.1(a) predicts it exceeds p on the transient day).
  eps        per station 1 - r_s, r_s the polarisation reliability (ttmt.v2.m0.orientation_audit) of the P vectors of
             the network's catalogue events, clipped to [0.05, 0.5].
  model      s_k = s0_station + 0.3 |A_k m|, s0 = median FIT read norm / sqrt(3) (same for R, T, G).
  pi         CAL: the catalogue tensors (least squares at the catalogue location), smoothed (prior_draws: lognormal
             amplitude sd 0.5, jitter 0.25). LOW (stress): the same directions rescaled so that the median predicted
             block norm over 200 lattice nodes is 2 x the median station s0, the hardest place for the level.
             J = 2000 draws per network, seeded, independent of the TEST data.
  hypotheses 40 random (node, origin) of the registered coarse lattice per TEST slot, seeded, independent of data.
Gate: with the R background, for pi CAL and pi LOW, the pooled TEST rate of S_mix >= log(1/p) is <= p for p = 1e-2
and 1e-3, on seed 2718; seed 3141 is the replicate and must pass too. Reported: 95 % upper bounds by a bootstrap of
time clusters (network block, slot), the mean of exp(S_mix), per cell rates, the gaussian ablation (predicted by D.1(a)
to exceed p on the transient day), and the out-of-sample fit of the t (PIT deciles of the TEST read norms).

    python run_c.py measure SEED NOISE_DIR[,NOISE_DIR] [j,j,...|all] [cell,...]     (step 0: the fit, seconds)
    python run_c.py run SEED NOISE_DIR[,NOISE_DIR] PARTSDIR [j,j,...|all] [cell,...]
    python run_c.py analyze PARTSDIR[,PARTSDIR] [B]
"""
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
for p in (os.path.join(ROOT, "code"), os.path.join(ROOT, "code", "explore_contrib"), os.path.join(ROOT, "src")):
    if p not in sys.path:
        sys.path.insert(0, p)
import numpy as np  # noqa: E402

import t2m_lib as L2  # noqa: E402
import c10_lib as C  # noqa: E402
import c10_run as R  # noqa: E402
import xc  # noqa: E402
from task1o_consistent import matched  # noqa: E402
from ttmt.v2.m0 import orientation_audit  # noqa: E402
from ttmt.v2.m5 import GaussBackground, RadialBackground, TBackground, evalue_mix, prior_draws  # noqa: E402

T2M = os.path.join(ROOT, "results", "task2M")
TAG_FIT, TAG_PI, TAG_HYP = 7301, 7302, 7303
N_FIT_READS, H_TEST, J_DRAWS = 400, 40, 2000
GAMMA, EPS_LO, EPS_HI, LOW_LEVEL = 0.3, 0.05, 0.5, 2.0
PS = (1e-2, 1e-3)
BGS = ("R", "T", "G")


def configure(seed, ndirs):
    with open(os.path.join(T2M, f"cfg_{seed}.json")) as fh:
        c = json.load(fh)
    L2.configure(c["seed"], ndirs, c["j0"], c["nblocks"], c["nwin"], ["MT@S0"])
    return c


def mfn(net, slot):
    MF, sd = matched(L2.real_segment(net.di, net.sub, slot))
    return MF / sd


def split_slots(net):
    _di, blk = L2.block(net.j)
    used = set(int(s) for s in net.cat_slots)
    rest = [int(s) for s in blk if int(s) not in used]
    return rest[0::2], rest[1::2]


def read(MFn, IDX):
    """hard reading: vectors (..., 2, ns, 3) at the energy peak within +-4 samples of IDX (..., 2, ns)."""
    Kt = xc.peak_table_w(MFn, C.SH)
    ns = MFn.shape[0]
    st = np.arange(ns).reshape((1,) * (IDX.ndim - 1) + (ns,))
    return C._vectors(MFn, Kt[st, IDX], None)


def catalogue_terms(net):
    """per station r_s from the P vectors of the catalogue events, and the catalogue tensors m_c (n, 6)."""
    vP, nP, mc = [[] for _ in range(net.ns)], [[] for _ in range(net.ns)], []
    for D, src, tau in net.cat:
        MF, sd = matched(D)
        E1 = C.Eng(np.asarray(src, float)[None], net.STA, net.corr, np.array([tau]), C.H_COARSE)
        v = read(MF / sd, E1.IDX[0, 0])                                        # (2, ns, 3)
        mc.append(np.linalg.lstsq(E1.G[0], v.reshape(-1), rcond=None)[0])
        n = (net.STA - np.asarray(src, float)[None]) / np.linalg.norm(net.STA - np.asarray(src, float)[None],
                                                                        axis=1, keepdims=True)
        for i in range(net.ns):
            vP[i].append(v[0, i])
            nP[i].append(n[i])
    r = np.array([orientation_audit(np.array(vP[i]), np.array(nP[i]))[1] for i in range(net.ns)])
    return r, np.array(mc)


def setup(cell, j, seed):
    net = L2.Net2M(cell, j, "REAL")
    E = C.Eng(R.GARR, net.STA, net.corr, R.TAUS, C.H_COARSE)
    fit, test = split_slots(net)
    rng = np.random.default_rng([seed, TAG_FIT, net.ci, j])
    V = [[] for _ in range(net.ns)]
    for s in fit:
        g, t = rng.integers(E.N, size=N_FIT_READS), rng.integers(E.T, size=N_FIT_READS)
        M = mfn(net, s)
        tp = xc.peak_table_w(M, C.SH)[np.arange(net.ns)[None, None], E.IDX[g, t]]   # (n, 2, ns) peak samples
        for i in range(net.ns):
            V[i].append(M[i][:, np.unique(tp[:, :, i])].T)                    # distinct peaks only (no duplicates)
    V = [np.concatenate(x) for x in V]
    bg = {"R": [RadialBackground.fit(V[i]) for i in range(net.ns)],
          "T": [TBackground.fit(V[i]) for i in range(net.ns)],
          "G": [GaussBackground.fit(V[i]) for i in range(net.ns)]}
    s0 = np.array([np.median(np.linalg.norm(x, axis=-1)) for x in V]) / math.sqrt(3.0)
    r, mc = catalogue_terms(net)
    eps = np.clip(1.0 - r, EPS_LO, EPS_HI)
    rp = np.random.default_rng([seed, TAG_PI, net.ci, j])
    gsel = rp.integers(E.N, size=200)
    Gs = E.G[gsel].reshape(200, 2 * net.ns, 3, 6)
    med = np.array([np.median(np.linalg.norm(np.einsum("nkcl,l->nkc", Gs, m), axis=-1)) for m in mc])
    pis = {"CAL": prior_draws(mc, J_DRAWS, rp),
           "LOW": prior_draws(mc, J_DRAWS, rp, scale=LOW_LEVEL * np.median(s0) / np.maximum(med, 1e-300))}
    return net, E, fit, test, bg, s0, eps, r, pis, med


def measure(seed, ndirs, nets, cells):
    """step 0 (before any S_mix): the t fit per station, its out-of-sample PIT on TEST reads, isotropy."""
    configure(seed, ndirs)
    nets = L2.CFG["nets"] if nets == "all" else [int(x) for x in nets.split(",")]
    for j in nets:
        for cell in cells:
            net, E, fit, test, bg, s0, eps, r, pis, med = setup(cell, j, seed)
            rng = np.random.default_rng([seed, TAG_HYP, net.ci, j, 0])
            Vt = []
            for s in test[:10]:
                g, t = rng.integers(E.N, size=N_FIT_READS), rng.integers(E.T, size=N_FIT_READS)
                Vt.append(read(mfn(net, s), E.IDX[g, t]))
            Vt = np.concatenate(Vt).transpose(2, 0, 1, 3).reshape(net.ns, -1, 3)
            print(f"== {cell} {j}: {len(fit)} FIT / {len(test)} TEST slots, {len(net.cat)} catalogue events")
            for i in range(net.ns):
                rn = np.linalg.norm(Vt[i], axis=1)
                cv = Vt[i].var(axis=0)
                print(f"   st {net.stations[i]:5s} s0 {s0[i]:.3f} eps {eps[i]:.2f} |v| q50 {np.median(rn):.2f} "
                      f"q99 {np.quantile(rn, .99):.2f} max {rn.max():.1f}  R alpha {bg['R'][i].alpha:.2f}  "
                      f"T nu {bg['T'][i].nu:.2f}  comp var {np.round(cv / cv.mean(), 2).tolist()}")
                for k in ("R", "T"):
                    u = bg[k][i].cdf_norm(rn)
                    dec = np.histogram(u, bins=np.linspace(0, 1, 11))[0] / len(u)
                    tail = [float((u > 1 - q).mean() / q) for q in (1e-2, 1e-3)]
                    print(f"      {k}: PIT deciles x10 {np.round(dec * 10, 2).tolist()}  max dev {np.abs(dec - .1).max():.3f}"
                          f"  tail ratio @1e-2 {tail[0]:.2f} @1e-3 {tail[1]:.2f}")
            print(f"   pi CAL median block norm {np.median(med):.2f}  (LOW rescaled to {LOW_LEVEL} x median s0 "
                  f"{np.median(s0):.3f})", flush=True)


def run_part(cell, j, seed, pdir):
    fn = os.path.join(pdir, f"{cell}_{j}.npz")
    if os.path.exists(fn):
        return
    net, E, fit, test, bg, s0, eps, r, pis, med = setup(cell, j, seed)
    ns = net.ns
    st = np.tile(np.arange(ns), 2)                                               # block order (b, i)
    out = {f"{k}_{pk}": [] for k in BGS for pk in pis}
    for s in test:
        rng = np.random.default_rng([seed, TAG_HYP, net.ci, j, s])
        g, t = rng.integers(E.N, size=H_TEST), rng.integers(E.T, size=H_TEST)
        v = read(mfn(net, s), E.IDX[g, t]).reshape(H_TEST, 2 * ns, 3)
        A = E.G[g].reshape(H_TEST, 2 * ns, 3, 6)
        lf = [np.stack([bg[k][i].logpdf(v[:, q]) for q, i in enumerate(st)], axis=1) for k in BGS]
        for pk, dr in pis.items():
            for k, S in zip(BGS, evalue_mix(v, A, lf, eps[st], dr, s0[st], GAMMA)):
                out[f"{k}_{pk}"].append(S)
    np.savez(fn, slots=np.array(test), cell=cell, j=j, seed=seed, eps=eps, r=r, med=med, s0=s0,
             alpha=np.array([f.alpha for f in bg["R"]]), nu=np.array([f.nu for f in bg["T"]]),
             **{k: np.array(v) for k, v in out.items()})
    print(f"{cell} {j}: {len(test)} TEST slots, R alpha median {np.median([f.alpha for f in bg['R']]):.2f}, "
          f"max S R_LOW {np.max(out['R_LOW']):.2f}", flush=True)


def run(seed, ndirs, pdir, nets, cells):
    configure(seed, ndirs)
    nets = L2.CFG["nets"] if nets == "all" else [int(x) for x in nets.split(",")]
    os.makedirs(pdir, exist_ok=True)
    for j in nets:
        for cell in cells:
            run_part(cell, j, seed, pdir)


def analyze(pdirs, B):
    parts = []
    for d in pdirs:
        parts += [dict(np.load(os.path.join(d, f), allow_pickle=False)) for f in sorted(os.listdir(d)) if f.endswith(".npz")]
    for seed in sorted({int(p["seed"]) for p in parts}):
        P = [p for p in parts if int(p["seed"]) == seed]
        n_eval = sum(p["R_CAL"].size for p in P)
        print(f"== seed {seed}: {len(P)} parts, {n_eval} TEST hypotheses, numpy {np.__version__}")
        # time clusters: (network block j, slot); the three cells of a block share the noise
        keys = sorted({(int(p["j"]), int(s)) for p in P for s in p["slots"]})
        kid = {k: i for i, k in enumerate(keys)}
        verdict = {}
        for name in [f"{k}_{pk}" for k in BGS for pk in ("CAL", "LOW")]:
            S = np.concatenate([p[name].reshape(-1) for p in P])
            n_nan = int(np.isnan(S).sum())
            S = np.where(np.isnan(S), np.inf, S)       # NaN came from inf - inf with f_k(v_k) = 0: it is S = +inf
            cl = np.concatenate([np.repeat([kid[(int(p["j"]), int(s))] for s in p["slots"]], p[name].shape[1])
                                 for p in P])
            e = np.exp(np.minimum(S, 700))
            line = f"   {name}: non-finite S {int((~np.isfinite(S)).sum())} (NaN read as +inf {n_nan})  mean exp(S) {e.mean():.3g}  median S {np.median(S):.2f}"
            rng = np.random.default_rng(0)
            for p_ in PS:
                x = (S >= math.log(1 / p_)).astype(float)
                ks, ns_ = np.bincount(cl, x, len(keys)), np.bincount(cl, None, len(keys))
                bt = [ks[i].sum() / ns_[i].sum() for i in (rng.integers(len(keys), size=len(keys)) for _ in range(B))]
                rate, ub = x.mean(), float(np.quantile(bt, 0.975))
                line += f"  | p {p_:g}: rate {rate:.2e} (95 % ub {ub:.2e}) {'<=' if rate <= p_ else '>'} p"
                verdict[(name, p_)] = rate <= p_
            print(line)
            for cell in L2.CELLS:
                if not any(str(p["cell"]) == cell for p in P):
                    continue
                Sc = np.concatenate([p[name].reshape(-1) for p in P if str(p["cell"]) == cell])
                Sc = np.where(np.isnan(Sc), np.inf, Sc)
                print(f"      {cell:15s} " + "  ".join(f"p {p_:g}: {(Sc >= math.log(1 / p_)).mean():.2e}" for p_ in PS))
        nu = np.concatenate([p["alpha"] for p in P])
        eps = np.concatenate([p["eps"] for p in P])
        print(f"   R tail index quantiles 5/50/95: {np.round(np.quantile(nu, [.05, .5, .95]), 2).tolist()};"
              f" eps 5/50/95: {np.round(np.quantile(eps, [.05, .5, .95]), 2).tolist()}")
        ok = all(verdict[(n, p_)] for n in ("R_CAL", "R_LOW") for p_ in PS)
        print(f"   G-M5c seed {seed}: {'PASS' if ok else 'FAIL'} (R background, pi CAL and LOW, p 1e-2 and 1e-3)\n")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "measure":
        measure(int(sys.argv[2]), sys.argv[3].split(","), sys.argv[4] if len(sys.argv) > 4 else "all",
                sys.argv[5].split(",") if len(sys.argv) > 5 else list(L2.CELLS))
    elif cmd == "run":
        run(int(sys.argv[2]), sys.argv[3].split(","), sys.argv[4], sys.argv[5] if len(sys.argv) > 5 else "all",
            sys.argv[6].split(",") if len(sys.argv) > 6 else list(L2.CELLS))
    elif cmd == "analyze":
        analyze(sys.argv[2].split(","), int(sys.argv[3]) if len(sys.argv) > 3 else 2000)
    else:
        raise SystemExit(__doc__)
