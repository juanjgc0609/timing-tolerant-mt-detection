"""C7 analysis.  Exact FROC recall at <= 0.1 / <= 0.25 realised false events per window (xc.recall_caps), hierarchical
bootstrap (strata = regime x cell; networks then windows; operating point re-chosen per replicate), paired differences,
fraction of the CEIL-GD headroom captured (ratio per replicate), empty-window score distributions, recall by pulse
parameters, runtime.
    python c7_analyze.py STD [DENSE ...] > out/c7_analysis_STD.txt"""
import json
import os
import sys

import numpy as np

import xc

B = int(os.environ.get("NBOOT", 1000))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
METHODS = ["BOMP-PK", "BOMP-PK-NX", "BOMP-PK-CAP", "BOMP-PK-PH", "BOMP-PK-PHF", "BOMP-PK-PH-FULL", "BOMP-PK-PHF-FULL",
           "BOMP-PK-PHF3", "BOMP-PK-PH-CAP", "BOMP-PK-PHF-CAP", "MT-FW", "MT-FW-PH", "MT-FW-PHF", "CEIL-GD"]
PAIRS = [("BOMP-PK-PH", "BOMP-PK"), ("BOMP-PK-PHF", "BOMP-PK"),
         ("BOMP-PK-PH-FULL", "BOMP-PK"), ("BOMP-PK-PHF-FULL", "BOMP-PK"), ("BOMP-PK-PHF3", "BOMP-PK"),
         ("BOMP-PK-NX", "BOMP-PK"), ("BOMP-PK-CAP", "BOMP-PK"),
         ("BOMP-PK-PH-CAP", "BOMP-PK"), ("BOMP-PK-PHF-CAP", "BOMP-PK"),
         ("BOMP-PK-PH-CAP", "BOMP-PK-CAP"), ("BOMP-PK-PHF-CAP", "BOMP-PK-CAP"),
         ("BOMP-PK-PH", "MT-FW"), ("BOMP-PK-PHF", "MT-FW"), ("BOMP-PK-PH-CAP", "MT-FW"), ("BOMP-PK-PHF-CAP", "MT-FW"),
         ("MT-FW", "BOMP-PK"), ("MT-FW-PH", "MT-FW"), ("MT-FW-PHF", "MT-FW"), ("MT-FW-PH", "BOMP-PK"),
         ("MT-FW-PHF", "BOMP-PK"), ("CEIL-GD", "BOMP-PK")]
FRAC = ["BOMP-PK-PH", "BOMP-PK-PHF", "BOMP-PK-PH-CAP", "BOMP-PK-PHF-CAP", "MT-FW", "MT-FW-PH", "MT-FW-PHF"]


def load(design):
    rows = []
    for reg in ("thin", "strong"):
        p = os.path.join(OUT, f"c7_{design}_{reg}.json")
        if not os.path.exists(p):
            continue
        for r in json.load(open(p))["rows"]:
            r["regime"] = reg
            rows.append(r)
    return rows


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


def block(rows, title, meths):
    nw, ne = len(rows), sum(r["k"] for r in rows)
    nn = len({(r["regime"], r["cell"], r["net"]) for r in rows})
    print(f"\n{title}: windows {nw}, events {ne}, networks {nn}")
    pt, reps = boot(rows, meths)
    print("  recall @0.1 / @0.25 with 95% CI")
    for m in meths:
        c = ci(reps[m])
        print(f"    {m:18s} {f3(pt[m][0], c[0][0], c[1][0], False)}   {f3(pt[m][1], c[0][1], c[1][1], False)}")
    print("  paired differences @0.1 / @0.25 (95% CI, P(diff <= 0))")
    out = {}
    for a, b in PAIRS:
        if a not in meths or b not in meths:
            continue
        d = reps[a] - reps[b]
        c = ci(d)
        p0 = np.mean(d <= 0, axis=0)
        pd_ = np.subtract(pt[a], pt[b])
        out[(a, b)] = (pd_, c, d.std(axis=0))
        print(f"    {a:16s} - {b:12s} @0.1 {f3(pd_[0], c[0][0], c[1][0])} P {p0[0]:.3f}   @0.25 "
              f"{f3(pd_[1], c[0][1], c[1][1])} P {p0[1]:.3f}   (boot SD {d.std(axis=0)[0]:.4f} / {d.std(axis=0)[1]:.4f})")
    if "CEIL-GD" in meths:
        print("  fraction of the CEIL-GD headroom captured, (m - BOMP-PK) / (CEIL-GD - BOMP-PK): point, median and 95% "
              "interval over replicates")
        H = reps["CEIL-GD"] - reps["BOMP-PK"]
        Hp = np.subtract(pt["CEIL-GD"], pt["BOMP-PK"])
        for m in FRAC:
            if m not in meths:
                continue
            G = reps[m] - reps["BOMP-PK"]
            with np.errstate(divide="ignore", invalid="ignore"):
                fr = np.where(np.abs(H) > 1e-9, G / H, np.nan)
                fp = np.subtract(pt[m], pt["BOMP-PK"]) / Hp
            s = []
            for c_ in range(2):
                v = fr[:, c_][np.isfinite(fr[:, c_])]
                s.append(f"{fp[c_]:+.2f} (median {np.median(v):+.2f}, [{np.percentile(v, 2.5):+.2f}, {np.percentile(v, 97.5):+.2f}])")
            print(f"    {m:16s} @0.1 {s[0]}   @0.25 {s[1]}")
    return pt, reps, out


def empty_scores(rows, meths):
    print("  empty windows (k = 0): top-atom score quantiles (median / P90 / max), n windows; and the per-method "
          "threshold at the pooled cap-0.25 operating point")
    E = [r for r in rows if r["k"] == 0]
    if not E:
        print("    (no empty windows in this design)")
        return
    ne = sum(r["k"] for r in rows)
    for m in meths:
        if m == "CEIL-GD":
            continue
        v = np.array([r["rec"][m][0][0] for r in E if r["rec"][m]])
        _rc, thr = xc.recall_caps([[tuple(x) for x in r["rec"][m]] for r in rows], ne)
        nx = np.array([r["rec"]["BOMP-PK-NX"][0][0] for r in E])
        print(f"    {m:18s} {np.median(v):7.1f} {np.percentile(v, 90):7.1f} {v.max():7.1f}  n {len(v):3d}   "
              f"mean ratio to NX {np.mean(v / nx):.3f}   thr@0.25 {thr[1]:7.1f}  thr@0.1 {thr[0]:7.1f}")


def false_split(rows, meths):
    """at each method's pooled cap-0.25 threshold: false atoms per window in empty vs event windows, events hit, and
    the first-atom hit rate / coverage (event with an atom within 5 km at any score)."""
    ne = sum(r["k"] for r in rows)
    E0 = [r for r in rows if r["k"] == 0]
    E1 = [r for r in rows if r["k"] > 0]
    print("  at the pooled cap-0.25 threshold: false atoms/window in empty | event windows; first-atom hit rate (event "
          "windows); coverage; median epicentral error of hits (km)")
    for m in meths:
        if m == "CEIL-GD":
            continue
        th = xc.recall_caps([[tuple(x) for x in r["rec"][m]] for r in rows], ne)[1][1]
        f0 = np.mean([sum(1 for s, h in r["rec"][m] if h < 0 and s >= th) for r in E0]) if E0 else np.nan
        f1 = np.mean([sum(1 for s, h in r["rec"][m] if h < 0 and s >= th) for r in E1])
        fh = np.mean([r["rec"][m][0][1] >= 0 for r in E1])
        cov, err = [], []
        for r in rows:
            L_ = [np.asarray(l[:2]) for l in r["loc"][m]]
            for t in r["truth"]:
                d = min(np.linalg.norm(l - np.asarray(t[:2])) for l in L_)
                cov.append(d <= 5)
                if d <= 5:
                    err.append(d)
        print(f"    {m:18s} {f0:6.3f} | {f1:6.3f}   first-hit {fh:.3f}   coverage {np.mean(cov):.3f}   err {np.median(err):.2f}")


def by_pulse(rows, meths):
    """per-event hit at each method's pooled cap-0.25 operating threshold, by |phase| and fc bins."""
    ne = sum(r["k"] for r in rows)
    thr = {m: xc.recall_caps([[tuple(x) for x in r["rec"][m]] for r in rows], ne)[1][1] for m in meths}
    ev = []
    for r in rows:
        for k in range(r["k"]):
            hit = {m: any(h == k and s >= thr[m] for s, h in r["rec"][m]) for m in meths}
            ev.append((abs(r["meta"][k]["phase_deg"]), r["meta"][k]["fc"], hit))
    print("  per-event recall at each method's pooled cap-0.25 threshold, by pulse parameters (n events)")
    bins = [("|phase| 0-15", lambda p, f: p < 15), ("|phase| 15-30", lambda p, f: 15 <= p < 30),
            ("|phase| 30-45", lambda p, f: 30 <= p < 45), ("|phase| 45-90", lambda p, f: p >= 45),
            ("fc < 1.75", lambda p, f: f < 1.75), ("fc 1.75-2.25", lambda p, f: 1.75 <= f < 2.25),
            ("fc >= 2.25", lambda p, f: f >= 2.25)]
    for name, sel in bins:
        sub = [e for e in ev if sel(e[0], e[1])]
        if not sub:
            continue
        print(f"    {name:14s} n {len(sub):4d}  " + "  ".join(f"{m} {np.mean([e[2][m] for e in sub]):.3f}" for m in meths))


def runtime(rows, meths):
    print("  runtime s/window (2 processes on 2 cores): " + ", ".join(
        f"{m} {np.mean([r['sec'][m] for r in rows]):.3f}" for m in meths))


def main():
    for design in sys.argv[1:]:
        rows = load(design)
        if not rows:
            continue
        meths = [m for m in METHODS if m in rows[0]["rec"]]
        print("=" * 118)
        print(f"DESIGN {design}")
        res = {}
        for reg in ("thin", "strong"):
            sub = [r for r in rows if r["regime"] == reg]
            if sub:
                res[reg] = block(sub, f"{design} {reg}", meths)
        if len(res) == 2:
            res["both"] = block(rows, f"{design} both regimes pooled (6 strata)", meths)
        print(f"\n{design} diagnostics (both regimes)")
        empty_scores(rows, meths)
        false_split(rows, meths)
        if all("meta" in r for r in rows):
            by_pulse(rows, ["BOMP-PK", "BOMP-PK-PH", "BOMP-PK-PHF", "BOMP-PK-PH-CAP", "MT-FW", "MT-FW-PH"])
        runtime(rows, [m for m in meths])
        pw = {reg: {f"{a}|{b}": [float(x) for x in v[2]] for (a, b), v in res[reg][2].items()} for reg in res}
        json.dump(pw, open(os.path.join(OUT, f"c7_bootsd_{design}.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
