"""C8 analysis.  Exact FROC recall at <= 0.1 / <= 0.25 realised false events per window (xc.recall_caps), hierarchical
bootstrap (strata = regime x cell; networks then windows; operating point re-chosen per replicate; B = 1000), paired
differences, fraction of the CEIL-GD headroom captured, catalogue amplitude-term quality, runtime.
    python c8_analyze.py DENSE-IID [STD-IID ...] > out/c8_analysis.txt"""
import json
import os
import sys

import numpy as np

import xc

B = int(os.environ.get("NBOOT", 1000))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
METHODS = ["BOMP-PK", "MT-FW", "MT-FW-NOCAP", "MT-FW-C3", "MT-FW-C8", "MT-FW-C12", "ORACLE", "ORACLE-CAP", "SITE-CAT",
           "SITE-CAT-CAP", "SITE-CAT-ST", "SITE-CAT-GAIN", "HUBER-N", "HUBER-N90", "HUBER-N99", "HUBER-C1.5", "HUBER-C5",
           "HUBER-C8", "HUBER-S", "GAIN", "GAIN-FREE", "NORM-TR", "CEIL-GD"]
CAND = ["ORACLE", "ORACLE-CAP", "SITE-CAT", "SITE-CAT-CAP", "SITE-CAT-ST", "SITE-CAT-GAIN", "HUBER-N", "HUBER-N90",
        "HUBER-N99", "HUBER-C1.5", "HUBER-C5", "HUBER-C8", "HUBER-S", "GAIN", "GAIN-FREE", "NORM-TR"]
PAIRS = ([("MT-FW", "BOMP-PK"), ("MT-FW-NOCAP", "BOMP-PK"), ("MT-FW-C3", "BOMP-PK"), ("MT-FW-C8", "BOMP-PK"),
          ("MT-FW-C12", "BOMP-PK"), ("CEIL-GD", "BOMP-PK")]
         + [(c, "MT-FW") for c in CAND] + [(c, "MT-FW-NOCAP") for c in CAND] + [(c, "BOMP-PK") for c in CAND]
         + [("SITE-CAT", "ORACLE"), ("GAIN", "HUBER-S"), ("HUBER-N90", "HUBER-N"), ("HUBER-N99", "HUBER-N"),
            ("HUBER-C1.5", "HUBER-N"), ("HUBER-C5", "HUBER-N"), ("HUBER-C8", "HUBER-N"),
            ("HUBER-N", "ORACLE-CAP"), ("HUBER-N", "SITE-CAT-CAP")])
FRAC = ["MT-FW", "MT-FW-NOCAP"] + CAND


def load(design):
    if design in ("DENSE-IID+H", "DENSE-IID+W"):
        rows, nets = load("DENSE-IID")
        hr, _n = load("DENSE-IID-HSWEEP" if design.endswith("H") else "DENSE-IID-HWIDE")
        extra = ("HUBER-N90", "HUBER-N99") if design.endswith("H") else ("HUBER-C1.5", "HUBER-C5", "HUBER-C8")
        key = {(r["regime"], r["cell"], r["net"], r["w"]): r for r in hr}
        same = tot = 0
        out = []
        for r in rows:
            h = key.get((r["regime"], r["cell"], r["net"], r["w"]))
            if h is None:
                continue
            assert h["truth"] == r["truth"]
            for m in ("MT-FW", "HUBER-N"):
                tot += 1
                same += int(h["loc"][m] == r["loc"][m] and np.allclose([x[0] for x in h["rec"][m]],
                                                                         [x[0] for x in r["rec"][m]], rtol=1e-9))
            for m in extra:
                r["rec"][m], r["loc"][m], r["sec"][m] = h["rec"][m], h["loc"][m], h["sec"][m]
            out.append(r)
        print(f"  HSWEEP merge: {len(out)} of {len(rows)} windows; identity of regenerated MT-FW / HUBER-N records "
              f"{same}/{tot}")
        return out, nets
    rows, nets = [], []
    for reg in ("thin", "strong"):
        p = os.path.join(OUT, f"c8_{design}_{reg}.json")
        if not os.path.exists(p):
            continue
        d = json.load(open(p))
        for r in d["rows"]:
            r["regime"] = reg
            rows.append(r)
        nets.extend(d["nets"])
    return rows, nets


def boot(rows, meths, seed=17):
    keys = sorted({(r["regime"], r["cell"], r["net"]) for r in rows})
    idx = {k: [] for k in keys}
    for i, r in enumerate(rows):
        idx[(r["regime"], r["cell"], r["net"])].append(i)
    strata = sorted({k[:2] for k in keys})
    recs = {m: [[tuple(x) for x in r["rec"][m]] for r in rows] for m in meths}
    ks = [r["k"] for r in rows]

    def rc(ws):
        n = sum(ks[w] for w in ws)
        return {m: xc.recall_caps([recs[m][w] for w in ws], n)[0] for m in meths}
    pt = rc(list(range(len(rows))))
    rng = np.random.default_rng(seed)
    reps = {m: [] for m in meths}
    for _ in range(B):
        ws = []
        for s in strata:
            nets = [k for k in keys if k[:2] == s]
            for g in rng.integers(0, len(nets), len(nets)):
                gg = idx[nets[g]]
                ws.extend(np.asarray(gg)[rng.integers(0, len(gg), len(gg))].tolist())
        r_ = rc(ws)
        for m in meths:
            reps[m].append(r_[m])
    return pt, {m: np.array(v) for m, v in reps.items()}


def ci(a):
    return np.percentile(a, [2.5, 97.5], axis=0)


def f3(p, lo, hi, sign=True):
    return f"{p:+.3f} [{lo:+.3f}, {hi:+.3f}]" if sign else f"{p:.3f} [{lo:.3f}, {hi:.3f}]"


def block(rows, title, meths, pairs=PAIRS):
    nw, ne = len(rows), sum(r["k"] for r in rows)
    nn = len({(r["regime"], r["cell"], r["net"]) for r in rows})
    print(f"\n{title}: windows {nw}, events {ne}, network-regime units {nn}")
    pt, reps = boot(rows, meths)
    print("  recall @0.1 / @0.25 with 95% CI")
    for m in meths:
        c = ci(reps[m])
        print(f"    {m:14s} {f3(pt[m][0], c[0][0], c[1][0], False)}   {f3(pt[m][1], c[0][1], c[1][1], False)}")
    print("  paired differences @0.1 / @0.25 (95% CI, P(diff <= 0), boot SD @0.25)")
    for a, b in pairs:
        if a not in meths or b not in meths:
            continue
        d = reps[a] - reps[b]
        c = ci(d)
        p0 = np.mean(d <= 0, axis=0)
        pd_ = np.subtract(pt[a], pt[b])
        c90 = np.percentile(d[:, 1], [5, 95])
        print(f"    {a:13s} - {b:11s} @0.1 {f3(pd_[0], c[0][0], c[1][0])} P {p0[0]:.3f}   @0.25 "
              f"{f3(pd_[1], c[0][1], c[1][1])} P {p0[1]:.3f}  SD {d.std(axis=0)[1]:.4f}  90% [{c90[0]:+.3f}, "
              f"{c90[1]:+.3f}]")
    if "CEIL-GD" in meths:
        print("  fraction of the CEIL-GD headroom captured at 0.25, (m - BOMP-PK) / (CEIL-GD - BOMP-PK): point, 95% "
              "interval over replicates")
        H = reps["CEIL-GD"][:, 1] - reps["BOMP-PK"][:, 1]
        Hp = pt["CEIL-GD"][1] - pt["BOMP-PK"][1]
        for m in FRAC:
            if m not in meths:
                continue
            G = reps[m][:, 1] - reps["BOMP-PK"][:, 1]
            with np.errstate(divide="ignore", invalid="ignore"):
                fr = np.where(np.abs(H) > 1e-9, G / H, np.nan)
            v = fr[np.isfinite(fr)]
            print(f"    {m:14s} {(pt[m][1] - pt['BOMP-PK'][1]) / Hp:+.2f} [{np.percentile(v, 2.5):+.2f}, "
                  f"{np.percentile(v, 97.5):+.2f}]")
    return pt, reps


def site_quality(nets):
    """catalogue amplitude terms vs the true persistent site (pers worlds only), per network-regime."""
    cs, cst, frac = [], [], []
    for n in nets:
        if n.get("site_net") is None:
            continue
        tru = np.log(np.asarray(n["site_net"])).T                      # (2, ns)
        u = np.asarray(n["u_cat"])
        ok = np.asarray(n["n_obs"]) >= 3
        frac.append(ok.mean())
        if ok.sum() >= 3:
            cs.append(np.corrcoef(u[ok], tru[ok])[0, 1])
        ust = np.asarray(n["u_cat_st"])[0]
        tst = tru.mean(0)
        if np.ptp(ust) > 0:
            cst.append(np.corrcoef(ust, tst)[0, 1])
    if cs:
        print(f"  catalogue amplitude terms vs true log site (per network-regime, blocks with >= 3 picks): corr median "
              f"{np.median(cs):.2f} [IQR {np.percentile(cs, 25):.2f}, {np.percentile(cs, 75):.2f}], n {len(cs)}; "
              f"fraction of blocks estimated {np.mean(frac):.2f}; station-level (-ST vs mean of P,S log site) median "
              f"{np.median(cst):.2f}")


def runtime(rows, meths):
    print("  runtime s/window (2 processes on 2 cores): " + ", ".join(
        f"{m} {np.mean([r['sec'][m] for r in rows]):.3f}" for m in meths if m in rows[0]["sec"]))


def main():
    for design in sys.argv[1:]:
        rows, nets = load(design)
        if not rows:
            continue
        meths = [m for m in METHODS if m in rows[0]["rec"]]
        print("=" * 118)
        print(f"DESIGN {design}")
        site_quality(nets)
        block(rows, f"{design} pooled (thin + strong)", meths)
        for reg in ("thin", "strong"):
            sub = [r for r in rows if r["regime"] == reg]
            if sub:
                block(sub, f"{design} {reg}", meths,
                      pairs=[p for p in PAIRS if p[1] in ("MT-FW", "BOMP-PK") and p[0] in
                             ("MT-FW", "ORACLE", "ORACLE-CAP", "SITE-CAT", "SITE-CAT-CAP", "HUBER-N", "HUBER-N90",
                              "HUBER-N99", "HUBER-S", "GAIN", "GAIN-FREE", "NORM-TR", "MT-FW-NOCAP")])
        runtime(rows, meths)


if __name__ == "__main__":
    main()
