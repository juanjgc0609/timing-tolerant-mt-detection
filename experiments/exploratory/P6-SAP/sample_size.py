"""sample_size.py - required sample size for the phase 6 tests (docs/protocol/PHASE6_SAP.md, Section 9; rule of ADR-0008).

It simulates the registered decision rules themselves (ttmt.stats): the block sign-flip test of each pairwise
contrast, the intersection-union over comparators and the fixed sequence H1 -> H2. Targets are clustered in time
blocks and the two methods' hits are correlated, so the answer is a number of blocks and of targets, not a number of
independent pairs.

Model of one design (all parameters are inputs; defaults are the R1 / P6-power values):
  - B blocks; targets per block ~ negative binomial with mean m and dispersion k (k = 0 means Poisson);
  - each target is in the target regime R* with probability f, otherwise in the complement;
  - method m hits target i of block j when z_mi < Phi^-1(p_m), z_mi = sqrt(icc) u_j + sqrt(1 - icc) e_mi,
    u_j ~ N(0, 1) shared by all methods (a noisy hour hurts every method), corr(e_mi, e_m'i) chosen so that the total
    latent correlation of two methods on the same target is rho;
  - the proposed method's latent threshold gets a block shift v_j ~ N(0, tau^2) (heterogeneity of the effect);
  - the best comparator has recall p0; the proposed method has p0 + delta_r in R* and p0 + delta_c in the complement;
    K comparators, all at p0 by default (the worst case for the intersection-union test).

Subcommands
  power   power and size of one design;
  solve   smallest number of blocks (and expected targets) reaching the target power;
  calibrate  level of H2 whose simulated size is <= 0.025 for a design (ADR-0008 item 7);
  hours   hours of record the false-alarm budget check needs;
  connor  closed-form number of independent targets (Connor 1987), a lower bound used to bracket `solve`;
  estimate  p, rho and icc from a CSV of development hits (columns: block, regime, then one 0/1 column per method);
  table   the default design table written to sample_size_out.txt.
Run `python3 sample_size.py <subcommand> -h` for the options. numpy only.
"""
import argparse
import csv
import json
from math import sqrt
from pathlib import Path
from statistics import NormalDist

import numpy as np

from ttmt import stats as st

ND = NormalDist()
OUT = Path(__file__).with_name("sample_size_out.txt")


# ---------------------------------------------------------------- bivariate normal and the latent model
def phi2(a, b, rho):
    """P(Z1 < a, Z2 < b) for standard bivariate normal with correlation rho (numerical integration)."""
    if rho >= 1 - 1e-12:
        return ND.cdf(min(a, b))
    x = np.linspace(-9.0, a, 6001)
    pdf = np.exp(-x ** 2 / 2) / sqrt(2 * np.pi)
    cond = np.array([ND.cdf((b - rho * xi) / sqrt(1 - rho ** 2)) for xi in x])
    trapezoid = getattr(np, "trapezoid", None) or np.trapz      # numpy < 2 (the core lock pins 1.26.4)
    return float(trapezoid(pdf * cond, x))


def discordance(p_a, p_b, rho):
    """(p10, p01): P(A hits, B misses), P(A misses, B hits) under the Gaussian copula."""
    both = phi2(ND.inv_cdf(p_a), ND.inv_cdf(p_b), rho)
    return p_a - both, p_b - both


def connor_n(p0, delta, rho, alpha=0.025, power=0.8, margin=0.0):
    """Independent targets for a one-sided paired test (Connor 1987; with a margin, its non-inferiority form)."""
    p10, p01 = discordance(p0 + delta, p0, rho)
    psi, dd = p10 + p01, p10 - p01
    za, zb = ND.inv_cdf(1 - alpha), ND.inv_cdf(power)
    return int(np.ceil((za * sqrt(psi) + zb * sqrt(max(psi - dd ** 2, 1e-12))) ** 2 / (dd + margin) ** 2))


def simulate_hits(rng, n_blocks, m, k, f, p0, delta_r, delta_c, rho, icc, tau, n_comp):
    """One simulated design. Returns block ids, regime mask (True = R*), hits of A (n,) and of comparators (K, n)."""
    if k and k > 0:
        lam = rng.gamma(k, m / k, n_blocks)
        sizes = rng.poisson(lam)
    else:
        sizes = rng.poisson(m, n_blocks)
    block = np.repeat(np.arange(n_blocks), sizes)
    n = block.size
    reg = rng.random(n) < f
    rho_e = (rho - icc) / (1 - icc)
    u = rng.standard_normal(n_blocks)[block]
    common = rng.standard_normal(n)                     # equicorrelated e: sqrt(rho_e) common + sqrt(1 - rho_e) own
    def latent():
        e = sqrt(rho_e) * common + sqrt(1 - rho_e) * rng.standard_normal(n)
        return sqrt(icc) * u + sqrt(1 - icc) * e
    v = rng.normal(0.0, tau, n_blocks)[block] if tau > 0 else 0.0
    # the block shift v ~ N(0, tau^2) is added to the threshold, so the threshold is scaled by sqrt(1 + tau^2) to keep
    # the marginal recall at p0 + delta (review of PR #11)
    thr_a = np.where(reg, ND.inv_cdf(p0 + delta_r), ND.inv_cdf(p0 + delta_c)) * sqrt(1 + tau ** 2) + v
    hit_a = latent() < thr_a
    thr0 = ND.inv_cdf(p0)
    hit_c = np.array([latent() < thr0 for _ in range(n_comp)])
    return block, reg, hit_a, hit_c


def run_design(n_blocks, m=3.0, k=1.0, f=0.3, p0=0.40, delta_r=0.10, delta_c=0.0, rho=0.7, icc=0.10, tau=0.0,
               n_comp=2, alpha=0.025, margin=0.05, sims=1000, n_flip=1999, seed=20261019, alpha_h2=None, keep_p=False):
    """Power of H1 (superiority over every comparator in R*), H2 (non-inferiority by margin in the complement) and of
    the fixed sequence H1 -> H2; mean targets in R* and in total."""
    rng = np.random.default_rng(seed)
    # one sign matrix for the whole design: signflip_pvalue drops zero blocks and uses the first columns, and the
    # columns are i.i.d., so any subset is a valid sign matrix for the non-zero blocks
    signs = rng.integers(0, 2, size=(n_flip, n_blocks), dtype=np.int8).astype(np.float64) * 2 - 1
    r1 = r2 = both = 0
    n_r = n_t = 0
    p2_all = []
    for _ in range(sims):
        block, reg, ha, hc = simulate_hits(rng, n_blocks, m, k, f, p0, delta_r, delta_c, rho, icc, tau, n_comp)
        n_r += reg.sum(); n_t += reg.size
        p1, p2 = [], []
        for c in hc:
            d, n = st.block_sums(ha[reg], c[reg], block[reg], n_blocks)
            p1.append(st.signflip_pvalue(d, n, signs=signs))
            d, n = st.block_sums(ha[~reg], c[~reg], block[~reg], n_blocks)
            p2.append(st.signflip_pvalue(d, n, margin=margin, signs=signs))
        a2 = alpha if alpha_h2 is None else alpha_h2
        h1, h2 = st.iut(p1) <= alpha, st.iut(p2) <= a2
        p2_all.append(st.iut(p2))
        r1 += h1
        r2 += h2
        both += h1 and h2
    out = {"blocks": n_blocks, "targets_R": n_r / sims, "targets": n_t / sims, "power_H1": r1 / sims,
           "power_H2": r2 / sims, "power_chain": both / sims}
    if keep_p:
        out["p_H2"] = np.array(p2_all)
    return out


def solve(target="power_chain", power=0.8, lo=25, hi=6000, **kw):
    """Smallest number of blocks with power >= target: doubling from lo, then bisection (common random numbers).
    Returns (None, last result) when hi blocks do not reach the power."""
    n, res = lo, run_design(lo, **kw)
    while res[target] < power:
        if n >= hi:
            return None, res
        lo, n = n, min(2 * n, hi)
        res = run_design(n, **kw)
    best, top = res, n
    while top - lo > max(2, lo // 20):
        mid = (lo + top) // 2
        r = run_design(mid, **kw)
        if r[target] >= power:
            top, best = mid, r
        else:
            lo = mid
    return top, best


ALPHA_GRID = (0.025, 0.02, 0.015, 0.01, 0.005)


def h2_sizes(n_blocks, rhos=(0.5, 0.7, 0.9), iccs=(0.0, 0.2), sims=4000, margin=0.05, **kw):
    """Simulated size of H2 at every level of ALPHA_GRID, worst case over a parameter range, at the least favourable
    configuration of the intersection-union test: one comparator exactly at the margin (n_comp 1; review of PR #11).
    Returns {level: worst size}."""
    worst = {a: 0.0 for a in ALPHA_GRID}
    for rho in rhos:
        for icc in iccs:
            if icc > rho:
                continue
            p = run_design(n_blocks, delta_r=0.10, delta_c=-margin, margin=margin, rho=rho, icc=icc, n_comp=1,
                           sims=sims, n_flip=999, keep_p=True, **kw)["p_H2"]
            for a in ALPHA_GRID:
                worst[a] = max(worst[a], float(np.mean(p <= a)))
    return worst


def calibrate_alpha_h2(n_blocks, target=0.025, sims=4000, **kw):
    """Level at which H2 is run for a design (ADR-0008 item 7): the largest level of ALPHA_GRID whose worst simulated
    size plus 1.645 Monte Carlo standard errors is at most target. Returns (level or None, its size, all sizes)."""
    sizes = h2_sizes(n_blocks, sims=sims, **kw)
    for a in ALPHA_GRID:
        s = sizes[a]
        if s + 1.645 * sqrt(max(s, 1e-3) * (1 - s) / sims) <= target:
            return a, s, sizes
    return None, sizes[ALPHA_GRID[-1]], sizes


def hours_needed(budget=30.0, phi=5.0, tol=0.20, level=0.975, pass_prob=0.80):
    """Hours of record for the budget check (SAP S7.6) to pass with probability pass_prob when the true rate equals
    the budget, with quasi-Poisson overdispersion phi of the false counts:
    H >= phi * (z_level * sqrt(1 + tol) + z_pass)^2 / (tol^2 * budget). The estimated dispersion floored at 1 makes
    the realised pass probability a little lower (0.75 for a nominal 0.80 at phi 1, tests/test_stats.py); plan 1.25 times
    these hours."""
    z = ND.inv_cdf(level) * sqrt(1 + tol) + ND.inv_cdf(pass_prob)
    return phi * z * z / (tol * tol * budget)


# ---------------------------------------------------------------- estimation from development hits
def _bisect(fn, target, lo=0.0, hi=0.999):
    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if fn(mid) < target else (lo, mid)
    return (lo + hi) / 2


def estimate(path):
    """Recall per method, tetrachoric rho of each pair, and latent icc per method (one-way ANOVA ICC of the binary
    hits inverted through the latent model) from development hits; blocks as in the registration (SAP S7.5)."""
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    methods = [c for c in rows[0] if c not in ("block", "regime")]
    blocks = np.array([int(r["block"]) for r in rows])
    H = {m: np.array([int(r[m]) for r in rows], float) for m in methods}
    out = {"n_targets": len(rows), "n_blocks": int(np.unique(blocks).size), "recall": {}, "rho": {}, "icc": {}}
    for m in methods:
        p = H[m].mean(); out["recall"][m] = float(p)
        icc_bin = _anova_icc(H[m], blocks)
        a = ND.inv_cdf(min(max(p, 1e-6), 1 - 1e-6))
        out["icc"][m] = float(_bisect(lambda r: (phi2(a, a, r) - p * p) / (p * (1 - p)), max(icc_bin, 0.0)))
    for i, m1 in enumerate(methods):
        for m2 in methods[i + 1:]:
            p1, p2, p11 = H[m1].mean(), H[m2].mean(), (H[m1] * H[m2]).mean()
            a, b = ND.inv_cdf(p1), ND.inv_cdf(p2)
            out["rho"][f"{m1}|{m2}"] = float(_bisect(lambda r: phi2(a, b, r), p11, -0.999, 0.999))
    return out


def _anova_icc(y, g):
    ids, inv = np.unique(g, return_inverse=True)
    n_i = np.bincount(inv).astype(float)
    mean_i = np.bincount(inv, weights=y) / n_i
    N, G, grand = y.size, ids.size, y.mean()
    msb = np.sum(n_i * (mean_i - grand) ** 2) / (G - 1)
    msw = np.sum((y - mean_i[inv]) ** 2) / (N - G)
    n0 = (N - np.sum(n_i ** 2) / N) / (G - 1)
    return float((msb - msw) / (msb + (n0 - 1) * msw))


# ---------------------------------------------------------------- per-cell rule against the pooled rule (K2)
def per_cell_vs_pooled(n_cells=10, per_cell=100, p0=0.40, delta=0.05, rho=0.7, alpha=0.025, sims=2000, seed=20261023):
    """Independent targets, one comparator. Per-cell rule (ATLAS-0 draft): a cell is won when its exact one-sided
    McNemar test rejects at alpha (the paired 95 % interval above zero); the method wins when it wins a majority of
    cells. Pooled rule (ADR-0008): one exact McNemar test on all target-regime targets. Returns (mean share of cells
    won, P(majority of cells won), power of the pooled test)."""
    rng = np.random.default_rng(seed)
    ta, tb = ND.inv_cdf(p0 + delta), ND.inv_cdf(p0)
    share = maj = pooled = 0.0
    for _ in range(sims):
        z1 = rng.standard_normal((n_cells, per_cell))
        z2 = rho * z1 + sqrt(1 - rho ** 2) * rng.standard_normal((n_cells, per_cell))
        ha, hb = z1 < ta, z2 < tb
        b, c = (ha & ~hb).sum(1), (~ha & hb).sum(1)
        won = np.array([st.exact_mcnemar(int(x), int(y)) <= alpha for x, y in zip(b, c)])
        share += won.mean(); maj += won.sum() > n_cells / 2
        pooled += st.exact_mcnemar(int(b.sum()), int(c.sum())) <= alpha
    return share / sims, maj / sims, pooled / sims


# ---------------------------------------------------------------- default table
def table(sims=400):
    lines = ["Sample size of the single beat rule (ADR-0008) and of the phase 6 tests (docs/protocol/PHASE6_SAP.md, S9).",
             "Block sign-flip tests, intersection-union over 2 comparators at equal recall 0.40, one-sided alpha 0.025,",
             "fixed sequence H1 (superiority in the target regime) -> H2 (non-inferiority, margin 0.05, true difference 0,",
             f"in the complement). Targets per block: negative binomial, mean 3, dispersion 1. {sims} simulations per",
             "design, seed 20261019.", ""]
    lines.append("1. Independent targets, closed form (Connor 1987), one comparator, alpha 0.025, power 0.80:")
    lines.append("   superiority   rho 0.5  rho 0.7  rho 0.9")
    for d in (0.05, 0.10, 0.15):
        lines.append(f"   +{d:.2f}       " + "".join(f"{connor_n(0.40, d, r):>9d}" for r in (0.5, 0.7, 0.9)))
    for mg in (0.03, 0.05):
        lines.append(f"   NI {mg:.2f}, true 0" + "".join(f"{connor_n(0.40, 0.0, r, margin=mg):>9d}" for r in (0.5, 0.7, 0.9)))
    lines.append("   (P6-power used alpha 0.05 one-sided; at the registered 0.025 the numbers are about 25 % larger.)")
    lines.append("")
    lines.append("2. K2: per-cell rule against one pooled test, true gain +0.05 in every target cell, rho 0.7, 2000 sims:")
    lines.append("   cells x targets   share of cells won   P(majority of cells won)   power of the pooled test")
    for nc, pc in ((10, 100), (10, 400), (24, 100), (24, 400)):
        a, b, c = per_cell_vs_pooled(nc, pc)
        lines.append(f"   {nc:>3d} x {pc:<4d}        {a:>8.2f}             {b:>10.3f}                 {c:>8.3f}")
    lines.append("")
    lines.append("3. Size at the least favourable configuration of the intersection-union test: one comparator at the")
    lines.append("   boundary (n_comp 1), 3 targets per block on average, f 0.3, 4000 simulations. H1 at gain 0 (exact test);")
    lines.append("   H2 at true difference -0.05 with the studentised sign-flip, worst case over rho 0.5/0.7/0.9 and icc 0/0.2,")
    lines.append("   and the calibrated level of H2 (largest level whose size + 1.645 MC se <= 0.025):")
    for nb in (30, 50, 100, 300):
        h1 = max(run_design(nb, delta_r=0.0, rho=r, icc=i, n_comp=1, sims=4000, n_flip=999)["power_H1"]
                 for r, i in ((0.5, 0.0), (0.9, 0.2)))
        a, sa, sizes = calibrate_alpha_h2(nb)
        lines.append(f"   {nb:>4d} blocks: size H1 {h1:.3f}   size H2 at 0.025 {sizes[0.025]:.3f}, at 0.02 {sizes[0.02]:.3f}, "
                     f"at 0.015 {sizes[0.015]:.3f}   calibrated level {a}")
    lines.append("")
    lines.append("4. Blocks for 80 % power, share of targets in the target regime f = 0.3: blocks (targets in regime / total)")
    lines.append("   gain  rho  icc   H1 alone                    chain H1 -> H2 (H2 at 0.015, calibrated at 50 blocks)")
    for d in (0.05, 0.10):
        for rho, icc in ((0.5, 0.1), (0.7, 0.1), (0.9, 0.1), (0.7, 0.3)):
            cells = []
            for tgt in ("power_H1", "power_chain"):
                b, r = solve(tgt, delta_r=d, rho=rho, icc=icc, f=0.3, sims=sims, hi=6400, alpha_h2=0.015)
                cells.append(f"{b:>5d} ({r['targets_R']:.0f} / {r['targets']:.0f})" if b else "> 6400 blocks")
            lines.append(f"   +{d:.2f} {rho}  {icc}   " + "".join(f"{c:<28s}" for c in cells))
    lines.append("")
    lines.append("5. Hours of record for the false-alarm budget check (30/h, tolerance 20 %, level 0.975, pass probability 0.80):")
    lines.append("   " + "   ".join(f"phi {p}: {hours_needed(phi=p):.1f} h" for p in (1, 5, 10, 20)))
    OUT.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


# ---------------------------------------------------------------- CLI
def _design_args(p):
    p.add_argument("--m", type=float, default=3.0, help="mean targets per block")
    p.add_argument("--k", type=float, default=1.0, help="negative-binomial dispersion of targets per block (0 = Poisson)")
    p.add_argument("--f", type=float, default=0.3, help="share of targets in the target regime R*")
    p.add_argument("--p0", type=float, default=0.40, help="recall of the best comparator")
    p.add_argument("--delta-r", type=float, default=0.10, help="recall gain of the proposed method in R*")
    p.add_argument("--delta-c", type=float, default=0.0, help="true recall difference in the complement")
    p.add_argument("--rho", type=float, default=0.7, help="latent correlation of two methods on the same target")
    p.add_argument("--icc", type=float, default=0.10, help="latent intra-block correlation (shared by methods)")
    p.add_argument("--tau", type=float, default=0.0, help="sd of the block shift of the proposed method (latent)")
    p.add_argument("--n-comp", type=int, default=2, help="number of comparators in the intersection-union test")
    p.add_argument("--alpha", type=float, default=0.025)
    p.add_argument("--alpha-h2", type=float, default=None, help="level of H2 (default: alpha; use `calibrate`)")
    p.add_argument("--margin", type=float, default=0.05, help="non-inferiority margin of H2 (ADR-0008)")
    p.add_argument("--sims", type=int, default=1000)
    p.add_argument("--seed", type=int, default=20261019)


def _kw(a):
    return dict(m=a.m, k=a.k, f=a.f, p0=a.p0, delta_r=a.delta_r, delta_c=a.delta_c, rho=a.rho, icc=a.icc, tau=a.tau,
                n_comp=a.n_comp, alpha=a.alpha, margin=a.margin, sims=a.sims, seed=a.seed, alpha_h2=a.alpha_h2)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("power"); p.add_argument("--blocks", type=int, required=True); _design_args(p)
    p = sub.add_parser("solve"); _design_args(p)
    p.add_argument("--target", default="power_chain", choices=["power_H1", "power_H2", "power_chain"])
    p.add_argument("--power", type=float, default=0.8); p.add_argument("--max-blocks", type=int, default=6000)
    p = sub.add_parser("calibrate"); p.add_argument("--blocks", type=int, required=True)
    p.add_argument("--m", type=float, default=3.0); p.add_argument("--k", type=float, default=1.0)
    p.add_argument("--f", type=float, default=0.3); p.add_argument("--margin", type=float, default=0.05)
    p.add_argument("--sims", type=int, default=4000)
    p = sub.add_parser("hours"); p.add_argument("--budget", type=float, default=30.0)
    p.add_argument("--phi", type=float, default=5.0); p.add_argument("--tol", type=float, default=0.2)
    p.add_argument("--pass-prob", type=float, default=0.8)
    p = sub.add_parser("connor"); p.add_argument("--p0", type=float, default=0.40)
    p.add_argument("--delta", type=float, default=0.10); p.add_argument("--rho", type=float, default=0.7)
    p.add_argument("--alpha", type=float, default=0.025)
    p.add_argument("--alpha-h2", type=float, default=None, help="level of H2 (default: alpha; use `calibrate`)"); p.add_argument("--power", type=float, default=0.8)
    p.add_argument("--margin", type=float, default=0.0)
    p = sub.add_parser("estimate"); p.add_argument("csv")
    p = sub.add_parser("table"); p.add_argument("--sims", type=int, default=400)
    a = ap.parse_args(argv)
    if a.cmd == "power":
        print(json.dumps(run_design(a.blocks, **_kw(a)), indent=1))
    elif a.cmd == "solve":
        b, r = solve(a.target, a.power, hi=a.max_blocks, **_kw(a))
        print(json.dumps({"blocks": b, **r}, indent=1))
    elif a.cmd == "calibrate":
        lvl, s, sizes = calibrate_alpha_h2(a.blocks, sims=a.sims, m=a.m, k=a.k, f=a.f, margin=a.margin)
        print(json.dumps({"alpha_h2": lvl, "size": s, "sizes": {str(k): v for k, v in sizes.items()}}, indent=1))
    elif a.cmd == "hours":
        print(f"{hours_needed(a.budget, a.phi, a.tol, pass_prob=a.pass_prob):.1f} hours")
    elif a.cmd == "connor":
        print(connor_n(a.p0, a.delta, a.rho, a.alpha, a.power, a.margin))
    elif a.cmd == "estimate":
        print(json.dumps(estimate(a.csv), indent=1))
    else:
        table(a.sims)


if __name__ == "__main__":
    main()
