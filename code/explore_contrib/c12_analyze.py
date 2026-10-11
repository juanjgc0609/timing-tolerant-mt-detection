"""C12 analysis.
    python c12_analyze.py eval DESIGN [B]   -> out/c12_analysis_DESIGN.txt   (recall, paired differences, weak events, cost)
    python c12_analyze.py null [B]          -> out/c12_null.txt              (null calibration and tube formula)
Hierarchical bootstrap as c10_analyze.Boot (strata regime x cell, networks then windows, operating point re-chosen per
replicate), extended with an event subset: hits on events outside the subset are ignored (neither hit nor false), the
false count and the operating point use every window."""
import glob
import json
import os
import sys

import numpy as np
from scipy import stats

import c10_lib as C
import c12_lib as L
import xc

REGIMES = ("thin", "strong")
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
CAPS = (0.1, 0.25)
PARTS = os.path.join(L.HERE, "out", "c12_parts")
# method: (atom set, score index in the packed atom [x, y, z, tau, s6, s5, sdc])
METHODS = {"FULL6": ("F6", 4), "DEV5": ("D5", 5), "DC": ("DC", 6), "MIX-DC": ("F6", 6), "MIX-DEV": ("F6", 5),
           "DCrank/s6": ("DC", 4)}


def load(kind):
    rows = []
    for f in sorted(glob.glob(os.path.join(PARTS, f"{kind}_*.json"))):
        rows += json.load(open(f))["rows"]
    return rows


def atoms_of(row, aset, k):
    return [{"x": np.array(a[:3]), "tau": a[3], "s": a[k]} for a in row["atoms"][aset]]


def build(rows, tol_t):
    recs = {}
    for nm, (aset, k) in METHODS.items():
        out = []
        for r in rows:
            at = atoms_of(r, aset, k)
            out.append(C.record(at, C.label(at, [np.array(t) for t in r["truth"]], r["taus"], tol_t), "s"))
        recs[nm] = out
    return recs


class Boot:
    def __init__(self, rows, recs, sel, B, sub=None, seed=12345):
        self.methods = list(recs)
        idx = [i for i, r in enumerate(rows) if sel(r)]
        self.idx = idx
        self.ks = np.array([rows[i]["k"] if sub is None else len(sub[i]) for i in idx])
        g = {}
        for p, i in enumerate(idx):
            r = rows[i]
            g.setdefault((r["regime"], r["cell"], r["net"]), []).append(p)
        self.groups = list(g.values())
        self.strata = [k[:2] for k in g]
        self.S = {m: [np.array([x[0] for x in recs[m][i]], float) for i in idx] for m in self.methods}
        self.H = {}
        for m in self.methods:
            hh = []
            for i in idx:
                h = np.array([x[1] for x in recs[m][i]], int)
                if sub is not None and len(h):
                    h = np.where((h >= 0) & ~np.isin(h, list(sub[i])), -2, h)
                hh.append(h)
            self.H[m] = hh
        allw = np.arange(len(idx))
        self.pt = {m: self.recall(m, allw) for m in self.methods}
        if sub is None:
            for m in self.methods[:3]:
                ref = xc.recall_caps([recs[m][i] for i in idx], int(self.ks.sum()), CAPS)[0]
                assert np.allclose(ref, self.pt[m]), (m, ref, self.pt[m])
        rng = np.random.default_rng(seed)
        sidx = {}
        for gi, s in enumerate(self.strata):
            sidx.setdefault(s, []).append(gi)
        self.rep = {m: np.zeros((B, len(CAPS))) for m in self.methods}
        for b in range(B):
            ws = []
            for s, gl in sidx.items():
                for gi in np.asarray(gl)[rng.integers(0, len(gl), len(gl))]:
                    gg = np.asarray(self.groups[gi])
                    ws.extend(gg[rng.integers(0, len(gg), len(gg))].tolist())
            ws = np.asarray(ws)
            for m in self.methods:
                self.rep[m][b] = self.recall(m, ws)

    def recall(self, m, ws):
        nte = int(self.ks[ws].sum())
        S, H = self.S[m], self.H[m]
        sc = np.concatenate([S[w] for w in ws])
        if not len(sc) or nte == 0:
            return [0.0] * len(CAPS)
        hit = np.concatenate([H[w] for w in ws])
        pos = np.repeat(np.arange(len(ws)), [len(S[w]) for w in ws])
        key = np.where(hit >= 0, pos * 64 + hit, -1)
        o = np.argsort(-sc, kind="stable")
        sc, hit, key = sc[o], hit[o], key[o]
        _u, first = np.unique(key, return_index=True)
        newhit = np.zeros(len(sc), bool)
        newhit[first] = True
        newhit &= key >= 0
        Hc, Fc = np.cumsum(newhit), np.cumsum(hit == -1)
        last = np.r_[np.where(np.diff(sc) != 0)[0], len(sc) - 1]
        Hc, Fc = Hc[last], Fc[last]
        out = []
        for c in CAPS:
            ok = Fc / len(ws) <= c
            out.append(float(Hc[ok].max() / nte) if ok.any() else 0.0)
        return out

    def ci(self, m):
        return np.percentile(self.rep[m], [2.5, 97.5], axis=0).T

    def diff(self, a, b):
        d = self.rep[a] - self.rep[b]
        return np.subtract(self.pt[a], self.pt[b]), np.percentile(d, [2.5, 97.5], axis=0).T, d.std(axis=0)


def fmt(p, ci):
    return f"{p:+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}]"


def terciles(rows):
    """event subsets by amplitude tercile within cell (rss arrival amplitude in noise sd)."""
    cut = {}
    for cell in CELLS:
        a = np.array([x for r in rows if r["cell"] == cell for x in r["amp"]])
        cut[cell] = np.quantile(a, [1 / 3, 2 / 3]) if len(a) else np.array([0.0, 0.0])
    subs = [dict(), dict(), dict()]
    for i, r in enumerate(rows):
        t = np.searchsorted(cut[r["cell"]], r["amp"], side="right") if r["k"] else np.array([], int)
        for q in range(3):
            subs[q][i] = set(np.where(t == q)[0].tolist())
    return subs, cut


def block(P, rows, recs, title, sel, B, sub=None):
    bt = Boot(rows, recs, sel, B, sub)
    P("")
    P(f"== {title}: {len(bt.idx)} windows, {int(bt.ks.sum())} events")
    P("   recall @0.1 / @0.25 [95%]                                  paired difference to FULL6 @0.1 / @0.25 [95%] (sd)")
    for m in bt.methods:
        c = bt.ci(m)
        s = f"   {m:10s} {bt.pt[m][0]:.3f} [{c[0][0]:.3f},{c[0][1]:.3f}]  {bt.pt[m][1]:.3f} [{c[1][0]:.3f},{c[1][1]:.3f}]"
        if m != "FULL6":
            p, ci, sd = bt.diff(m, "FULL6")
            s += f"   {fmt(p[0], ci[0])} / {fmt(p[1], ci[1])} (sd {sd[0]:.4f} / {sd[1]:.4f})"
        P(s)
    return bt


def main_eval(design, B=2000):
    rows = [r for r in load("eval") if r["design"] == design]
    lines = []

    def P(s=""):
        print(s, flush=True)
        lines.append(s)
    nets = sorted({r["net"] for r in rows})
    P(f"C12 ANALYSIS {design}: {len(rows)} windows, {sum(r['k'] for r in rows)} events, networks {nets[0]}-{nets[-1]} "
      f"({len(nets)}), bootstrap B = {B}")
    subs, cut = terciles(rows)
    P("   amplitude tercile cuts per cell (rss arrival amplitude over stations and phases, noise sd): " +
      "; ".join(f"{c} {cut[c][0]:.1f} / {cut[c][1]:.1f}" for c in CELLS))
    for rule, tol in (("LOC+T0.5", 0.5), ("LOC", None)):
        recs = build(rows, tol)
        block(P, rows, recs, f"hit rule {rule}, pooled", lambda r: True, B)
        for reg in REGIMES:
            block(P, rows, recs, f"hit rule {rule}, regime {reg}", lambda r, g=reg: r["regime"] == g, B)
        for cell in CELLS:
            block(P, rows, recs, f"hit rule {rule}, cell {cell}", lambda r, g=cell: r["cell"] == g, B)
        for q, nm in enumerate(("lowest", "middle", "highest")):
            block(P, rows, recs, f"hit rule {rule}, pooled, events of the {nm} amplitude tercile", lambda r: True, B,
                  sub=subs[q])
        # realised thresholds
        P("")
        P(f"   realised thresholds (statistic units) at <= 0.1 / <= 0.25 false per window, rule {rule}: " + "; ".join(
            "%s %.2f / %.2f" % ((m,) + tuple(xc.recall_caps(recs[m], sum(r["k"] for r in rows), CAPS)[1])) for m in recs))
    # retention of the statistic on hits and on false atoms (FULL6 atoms, LOC+T rule)
    P("")
    P("== FULL6 atoms: ratio of the restricted statistics to s6, by label (rule LOC+T) [median, quartiles]; and by "
      "rank of the atom")
    lab = {"hit": [], "false": []}
    for r in rows:
        at = atoms_of(r, "F6", 4)
        hh = C.label(at, [np.array(t) for t in r["truth"]], r["taus"], 0.5)
        for a, h in zip(r["atoms"]["F6"], hh):
            lab["hit" if h >= 0 else "false"].append(a[4:7])
    for k, v in lab.items():
        v = np.array(v)
        for nm, col in (("s5/s6", 1), ("sdc/s6", 2)):
            q = np.percentile(v[:, col] / v[:, 0], [25, 50, 75])
            P(f"   {k:5s} n {len(v):5d}  {nm}: {q[1]:.3f} [{q[0]:.3f},{q[2]:.3f}]   mean {np.mean(v[:, col] / v[:, 0]):.3f}"
              f"   mean loss s6 - stat: {np.mean(v[:, 0] - v[:, col]):.2f}")
    for lo, hi in ((0, 15), (15, 25), (25, 50), (50, 1e9)):
        s = f"   s6 in [{lo:g},{hi:g}):"
        for k, v in lab.items():
            v = np.array(v)
            mk = (v[:, 0] >= lo) & (v[:, 0] < hi)
            s += f" {k} n {int(mk.sum())} median sdc/s6 {np.median(v[mk, 2] / v[mk, 0]) if mk.any() else np.nan:.3f}"
            s += f" mean loss {np.mean(v[mk, 0] - v[mk, 2]) if mk.any() else np.nan:.2f};"
        P(s)
    # AUC-like separation of hits and false atoms by each score among FULL6 atoms
    h_, f_ = np.array(lab["hit"]), np.array(lab["false"])
    P("   separation of hit and false FULL6 atoms, P(score hit > score false): " + "; ".join(
        f"{nm} {stats.mannwhitneyu(h_[:, c], f_[:, c]).statistic / (len(h_) * len(f_)):.4f}"
        for nm, c in (("s6", 0), ("s5", 1), ("sdc", 2))))
    # agreement of the atom sets
    same = [np.mean([np.allclose(a[:4], b[:4]) for a, b in zip(r["atoms"]["F6"], r["atoms"][k])]) for r in rows
            for k in ("DC",) if r["atoms"]["F6"]]
    same1 = [np.allclose(r["atoms"]["F6"][0][:4], r["atoms"]["DC"][0][:4]) for r in rows if r["atoms"]["F6"]]
    P(f"   DC-ranked loop picks the same (node, origin) as FULL6: first atom {np.mean(same1):.3f} of the windows; "
      f"all atoms position by position {np.mean(same):.3f}")
    P("")
    P("== cost, s per window (one process): " + "; ".join(f"{k} {np.mean([r['sec'][k] for r in rows]):.3f}"
                                                          for k in ("F6", "D5", "DC")) +
      f"; DC evaluations per greedy iteration of the DC loop {sum(r['ndc'] for r in rows) / max(1, sum(r['nit'] for r in rows)):.1f}"
      f" of {675 * 31} candidates")
    with open(f"out/c12_analysis_{design}.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ============================================================================================ null
def thr_at(acc, ws, caps=CAPS):
    """thresholds with <= cap false per window: all atoms are false; acc = running-min scores per window (as in the
    records); scores strictly above the returned value are accepted."""
    sc = np.sort(np.concatenate([acc[w] for w in ws]))[::-1]
    out = []
    for cap in caps:
        n = int(np.floor(cap * len(ws)))
        out.append(float(sc[n]) if n < len(sc) else float("nan"))
    return out


def main_null(B=1000):
    lines = []

    def P(s=""):
        print(s, flush=True)
        lines.append(s)
    P("C12 NULL CALIBRATION")
    # ---- tube constants per cell
    P("")
    P("== tube constants of the DC manifold in the G-metric (30 random nodes x 2 networks per cell, 1500 Haar rotations "
      "per node): L3 = volume, in units of vol(S^3) = 2 pi^2; L1 = (1 / 4 pi) int scalar curvature")
    Lc = {}
    rngn = np.random.default_rng(3)
    for cell in CELLS:
        v = []
        for j in (1226, 1227):
            net = L.Net(cell, j, "thin")
            for g in rngn.integers(0, net.E.N, 30):
                v.append(L.lk_curvatures(net.tab.Rm[g], 1500)[:2])
        v = np.array(v)
        Lc[cell] = v
        P(f"   {cell:15s} L3 / vol(S^3): median {np.median(v[:, 0]) / (2 * np.pi ** 2):.3f} range "
          f"[{v[:, 0].min() / (2 * np.pi ** 2):.3f}, {v[:, 0].max() / (2 * np.pi ** 2):.3f}];  L1: median "
          f"{np.median(v[:, 1]):+.3f} range [{v[:, 1].min():+.3f}, {v[:, 1].max():+.3f}]")
    P(f"   Frobenius (isotropic) metric: L3 = 4 pi^2 = 2 vol(S^3), L1 = 0 exactly (the orbit is flat on average: Milnor's "
      "formula for the left-invariant metric (1, 1, 4) on SO(3) gives scalar curvature 0)")
    Lall = np.concatenate(list(Lc.values()))

    def tube_mix(x):
        return float(np.mean([L.tube_tail(x, a, b) for a, b in Lall]))

    def tube_q(alpha):
        from scipy.optimize import brentq
        return brentq(lambda x: np.log(tube_mix(x)) - np.log(alpha), 3.0, 400.0)
    P("")
    P("== theory, single candidate, z ~ N(0, I_6): upper quantile x of 2 s; ratio to FULL6; amplitude gain sqrt(x6 / x) - 1;"
      " effective dof = k with P(chi2_k > x_DC) = alpha")
    for al in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7):
        x6, x5, xd = stats.chi2.isf(al, 6), stats.chi2.isf(al, 5), tube_q(al)
        P(f"   alpha {al:.0e}: FULL6 {x6:6.2f}  DEV5 {x5:6.2f} (ratio {x5 / x6:.3f}, gain {np.sqrt(x6 / x5) - 1:+.3f})  "
          f"DC {xd:6.2f} (ratio {xd / x6:.3f}, gain {np.sqrt(x6 / xd) - 1:+.3f}); effective dof of DC {L.eff_dof(xd, al):.2f}")
    # ---- single-candidate simulation
    for kind, title in (("nulln", "event-free windows WITHOUT clutter (noise only), networks 1226-1229"),
                        ("nullc", "event-free windows WITH the clutter of the regime, networks 1220-1225")):
        rows = load(kind)
        if not rows:
            continue
        P("")
        P(f"== {title}: {len(rows)} windows")
        for read in ("fixed", "peak"):
            for reg in sorted({r["regime"] for r in rows}):
                c = 2.0 * np.array([x for r in rows if r["regime"] == reg for x in r["cand"][read]])
                P(f"   single candidate, read = {read} ({'predicted sample, no cap' if read == 'fixed' else 'registered peak rule, cap 6'}), "
                  f"regime {reg}: n = {len(c)}; mean of 2 s: FULL6 {c[:, 0].mean():.2f} DEV5 {c[:, 1].mean():.2f} DC {c[:, 2].mean():.2f}")
                for al in (0.1, 0.01, 0.001):
                    q = np.quantile(c, 1 - al, axis=0)
                    P(f"      upper quantile {al:g}: sim FULL6 {q[0]:6.2f} DEV5 {q[1]:6.2f} DC {q[2]:6.2f} | theory "
                      f"{stats.chi2.isf(al, 6):6.2f} {stats.chi2.isf(al, 5):6.2f} {tube_q(al):6.2f} | ratio DC/FULL6 sim "
                      f"{q[2] / q[0]:.3f} theory {tube_q(al) / stats.chi2.isf(al, 6):.3f}; DEV5/FULL6 sim {q[1] / q[0]:.3f} "
                      f"theory {stats.chi2.isf(al, 5) / stats.chi2.isf(al, 6):.3f}")
                if read == "fixed":
                    P("      tail probabilities sim / theory: " + "; ".join(
                        f"x={x:g}: FULL6 {np.mean(c[:, 0] > x):.2e}/{stats.chi2.sf(x, 6):.2e} DEV5 {np.mean(c[:, 1] > x):.2e}/"
                        f"{stats.chi2.sf(x, 5):.2e} DC {np.mean(c[:, 2] > x):.2e}/{tube_mix(x):.2e}" for x in (10.0, 14.0, 18.0)))
        # ---- window level
        P("   window level (greedy loop with deflation, 5 atoms; statistic units s): threshold with <= 0.1 / <= 0.25 false "
          "per window; ratio to FULL6 [95% bootstrap over networks then windows]; amplitude gain sqrt(thr6 / thr) - 1;"
          " theory = DC (or DEV5) threshold predicted from the FULL6 threshold by equal single-candidate tail")
        for reg in sorted({r["regime"] for r in rows}) + (["both"] if kind == "nullc" else []):
            rr = [r for r in rows if reg in ("both", r["regime"])]
            sc = {m: [np.array([a[k] for a in r["atoms"][aset]]) for r in rr] for m, (aset, k) in METHODS.items()}
            acc = {m: [np.minimum.accumulate(x) if len(x) else x for x in sc[m]] for m in METHODS}
            g = {}
            for i, r in enumerate(rr):
                g.setdefault((r["regime"], r["cell"], r["net"]), []).append(i)
            gl = list(g.values())
            rng = np.random.default_rng(99)
            reps = []
            for _ in range(B):
                ws = [w for gi in rng.integers(0, len(gl), len(gl)) for w in np.asarray(gl[gi])[rng.integers(0, len(gl[gi]), len(gl[gi]))]]
                reps.append([thr_at(acc[m], ws) for m in METHODS])
            reps = np.array(reps)                                         # (B, M, 2)
            pt = np.array([thr_at(acc[m], range(len(rr))) for m in METHODS])
            P(f"      regime {reg}: {len(rr)} windows, {len(gl)} network x cell x regime setups")
            for mi, m in enumerate(METHODS):
                s = f"      {m:10s}"
                for ci_, cap in enumerate(CAPS):
                    lo, hi = np.percentile(reps[:, mi, ci_], [2.5, 97.5])
                    s += f"  @{cap:g}: thr {pt[mi, ci_]:6.2f} [{lo:.2f},{hi:.2f}]"
                    if m != "FULL6":
                        rt = reps[:, mi, ci_] / reps[:, 0, ci_]
                        lo, hi = np.percentile(rt, [2.5, 97.5])
                        s += (f" ratio {pt[mi, ci_] / pt[0, ci_]:.3f} [{lo:.3f},{hi:.3f}] gain "
                              f"{np.sqrt(pt[0, ci_] / pt[mi, ci_]) - 1:+.3f} [{np.sqrt(1 / hi) - 1:+.3f},{np.sqrt(1 / lo) - 1:+.3f}]")
                P(s)
            c6 = 2.0 * np.array([x[0] for r in rr for x in r["cand"]["peak"]])
            csc = c6.var() / (2.0 * c6.mean())
            kk = c6.mean() / csc
            P(f"      peak-read single-candidate law of 2 s6 fitted by moments: {csc:.3f} x chi2 with {kk:.2f} dof (check: "
              f"upper 0.001 quantile fit {csc * stats.chi2.isf(1e-3, kk):.2f}, simulated {np.quantile(c6, 0.999):.2f})")
            for ci_, cap in enumerate(CAPS):
                x6 = 2 * pt[0, ci_]
                p6 = float(stats.chi2.sf(x6 / csc, kk))
                xd, x5, xf = tube_q(p6), stats.chi2.isf(p6, 5), stats.chi2.isf(p6, 6)
                P(f"      theory @{cap:g}: FULL6 thr {x6 / 2:.2f} = per-candidate tail {p6:.1e} under the fitted law (implied number "
                  f"of independent candidates {cap / p6:.0f}); Gaussian theory at that tail: DC / FULL6 threshold ratio "
                  f"{xd / xf:.3f} (gain {np.sqrt(xf / xd) - 1:+.3f}, effective dof of DC {L.eff_dof(xd, p6):.2f}), DEV5 / FULL6 "
                  f"{x5 / xf:.3f}; simulated ratios DC {pt[2, ci_] / pt[0, ci_]:.3f}, MIX-DC {pt[3, ci_] / pt[0, ci_]:.3f}, "
                  f"DEV5 {pt[1, ci_] / pt[0, ci_]:.3f}")
            mx = np.array([[s_[0] if len(s_) else np.nan for s_ in sc[m]] for m in ("FULL6", "DEV5", "DC")])
            P(f"      first atom (window maximum of the statistic): mean FULL6 {np.nanmean(mx[0]):.2f} DEV5 {np.nanmean(mx[1]):.2f} "
              f"DC {np.nanmean(mx[2]):.2f}; median {np.nanmedian(mx[0]):.2f} {np.nanmedian(mx[1]):.2f} {np.nanmedian(mx[2]):.2f}; "
              f"mean paired difference FULL6 - DC {np.nanmean(mx[0] - mx[2]):.2f} (sd {np.nanstd(mx[0] - mx[2]):.2f})")
    with open("out/c12_null.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    if sys.argv[1] == "eval":
        main_eval(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 2000)
    else:
        main_null(int(sys.argv[2]) if len(sys.argv) > 2 else 1000)
