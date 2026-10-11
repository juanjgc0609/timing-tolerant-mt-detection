"""C6 analysis: recall at 0.1 / 0.25 false per window (exact FROC, xc.recall_caps), hierarchical bootstrap (cells as
strata, networks then windows; 1000 replicates), differences vs HOM, fraction of the CEIL-GD headroom captured
(ratio computed per replicate), location errors of hits at the cap-0.25 operating point, atom coverage.
    python c6_analyze.py out/c6_thin.json out/c6_strong.json > out/c6_analysis.txt"""
import json
import os
import sys

import numpy as np

import xc

MODELS = ("HOM", "ORA", "ORA-TS", "CAL-H", "CAL-G", "CAL-2L")
DETS = ("BOMP-PK", "MT-FW")
CEILS = ("CEIL-GD", "CEIL-GD-O", "CEIL-GD-TT")
B = int(os.environ.get("NBOOT", 1000))


def methods():
    return [f"{d}:{m}" for d in DETS for m in MODELS] + list(CEILS)


def boot(rows, meths, B=B, seed=11):
    keys = sorted({(r["cell"], r["net"]) for r in rows})
    idx = {k: [] for k in keys}
    for i, r in enumerate(rows):
        idx[(r["cell"], r["net"])].append(i)
    cells = sorted({k[0] for k in keys})
    recs = {m: [[tuple(x) for x in r["rec"][m]] for r in rows] for m in meths}
    ks = [r["k"] for r in rows]

    def rc(ws):
        n = sum(ks[w] for w in ws)
        return {m: xc.recall_caps([recs[m][w] for w in ws], n)[0] for m in meths}
    allw = list(range(len(rows)))
    pt = rc(allw)
    rng = np.random.default_rng(seed)
    reps = {m: [] for m in meths}
    for _ in range(B):
        ws = []
        for c in cells:
            nets = [k for k in keys if k[0] == c]
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
    return (f"{p:+.3f} [{lo:+.3f}, {hi:+.3f}]" if sign else f"{p:.3f} [{lo:.3f}, {hi:.3f}]")


def loc_errors(rows, m):
    recs = [[tuple(x) for x in r["rec"][m]] for r in rows]
    n = sum(r["k"] for r in rows)
    _rc, thr = xc.recall_caps(recs, n)
    th = thr[1]
    ep, dp = [], []
    for r in rows:
        for (s, h), loc in zip(r["rec"][m], r["loc"][m]):
            if h >= 0 and s >= th:
                t = np.asarray(r["truth"][h])
                ep.append(float(np.linalg.norm(np.asarray(loc[:2]) - t[:2])))
                dp.append(float(abs(loc[2] - t[2])))
    return (np.median(ep) if ep else np.nan, np.median(dp) if dp else np.nan, len(ep))


def coverage(rows, m):
    """fraction of events with an atom (any score) within 5 km epicentral; median nearest-atom distance of the rest."""
    cov, far = [], []
    for r in rows:
        L = [np.asarray(l[:2]) for l in r["loc"][m]]
        for t in r["truth"]:
            dd = min([np.linalg.norm(l - np.asarray(t[:2])) for l in L] or [np.inf])
            cov.append(dd <= 5.0)
            if dd > 5.0:
                far.append(dd)
    return float(np.mean(cov)) if cov else np.nan, float(np.median(far)) if far else np.nan


def block(rows, title, out, per_cell=True):
    meths = methods()
    pt, reps = boot(rows, meths)
    nev = sum(r["k"] for r in rows)
    keys = sorted({(r["cell"], r["net"]) for r in rows})
    out.append(f"\n{title}: windows {len(rows)}, events {nev}, networks {len(keys)}")
    out.append("  recall @0.1 / @0.25 with 95% CI")
    for m in meths:
        c = ci(reps[m])
        out.append(f"    {m:18s} {f3(pt[m][0], c[0][0], c[1][0], False)}   {f3(pt[m][1], c[0][1], c[1][1], False)}")
    out.append("  differences (paired, same bootstrap replicates)")
    pairs = [(f"{d}:{m}", f"{d}:HOM") for d in DETS for m in MODELS if m != "HOM"]
    pairs += [("MT-FW:HOM", "BOMP-PK:HOM"), ("MT-FW:CAL-G", "BOMP-PK:HOM"), ("CEIL-GD", "BOMP-PK:HOM"),
              ("CEIL-GD-O", "BOMP-PK:HOM"), ("CEIL-GD", "MT-FW:HOM")]
    for a, b in pairs:
        d = reps[a] - reps[b]
        c = ci(d)
        p = np.subtract(pt[a], pt[b])
        pl = (d <= 0).mean(axis=0)
        out.append(f"    {a:14s} - {b:12s} @0.1 {f3(p[0], c[0][0], c[1][0])} P(<=0) {pl[0]:.3f}   "
                   f"@0.25 {f3(p[1], c[0][1], c[1][1])} P(<=0) {pl[1]:.3f}")
    out.append("  fraction of the CEIL-GD headroom (CEIL-GD - HOM of the same detector) captured by the model")
    for d in DETS:
        for m in MODELS[1:]:
            a, b = f"{d}:{m}", f"{d}:HOM"
            num = reps[a] - reps[b]
            den = reps["CEIL-GD"] - reps[b]
            fr = np.where(np.abs(den) > 1e-9, num / np.where(np.abs(den) > 1e-9, den, 1.0), np.nan)
            c = np.nanpercentile(fr, [2.5, 97.5], axis=0)
            dd = np.subtract(pt["CEIL-GD"], pt[b])
            p = np.where(dd > 1e-9, np.subtract(pt[a], pt[b]) / np.where(dd > 1e-9, dd, 1.0), np.nan)
            out.append(f"    {a:14s} @0.1 {f3(p[0], c[0][0], c[1][0])}   @0.25 {f3(p[1], c[0][1], c[1][1])}")
    out.append("  hits at the cap-0.25 operating point: median epicentral / depth error (km), n hits;"
               " atom coverage (event has an atom within 5 km at any score), median nearest-atom distance of the rest")
    for m in [f"{d}:{mm}" for d in DETS for mm in MODELS]:
        e, dz, n = loc_errors(rows, m)
        cv, far = coverage(rows, m)
        out.append(f"    {m:18s} epi {e:5.2f}  depth {dz:5.2f}  n {n:4d}   coverage {cv:.3f}  far-median {far:5.2f}")
    if per_cell:
        out.append("  per cell (point estimates @0.1/@0.25):")
        for cell in sorted({r["cell"] for r in rows}):
            sub = [r for r in rows if r["cell"] == cell]
            n = sum(r["k"] for r in sub)
            s = []
            for m in meths:
                v = xc.recall_caps([[tuple(x) for x in r["rec"][m]] for r in sub], n)[0]
                s.append(f"{m} {v[0]:.2f}/{v[1]:.2f}")
            out.append(f"    {cell} ({len(sub)} windows, {n} events): " + "  ".join(s))
    return pt, reps


def load(path):
    """c6_run json merged with the c6_ceil json (refined-deflation CEIL-GD / CEIL-GD-O, verif_2e definition); the
    true-arrival-time ceilings of c6_run are kept as CEIL-GD-TT / CEIL-GD-O-TT."""
    d = json.load(open(path))
    cp = path.replace("c6_", "c6_ceil_")
    cz = {(r["cell"], r["net"], r["w"]): r for r in json.load(open(cp))["rows"]}
    rows = []
    for r in d["rows"]:
        key = (r["cell"], r["net"], r["w"])
        if key not in cz:
            continue
        assert cz[key]["k"] == r["k"]
        r["rec"]["CEIL-GD-TT"], r["rec"]["CEIL-GD-O-TT"] = r["rec"]["CEIL-GD"], r["rec"]["CEIL-GD-O"]
        r["rec"]["CEIL-GD"], r["rec"]["CEIL-GD-O"] = cz[key]["rec"]["CEIL-GD"], cz[key]["rec"]["CEIL-GD-O"]
        rows.append(r)
    d["rows"] = rows
    return d


def main():
    out = []
    allrows = []
    for path in sys.argv[1:]:
        d = load(path)
        rows = d["rows"]
        allrows.extend(rows)
        reg = d["nets"][0]["regime"]
        out.append("=" * 118)
        out.append(f"REGIME {reg}  ({path})")
        # calibration summary
        for fam in ("CAL-H", "CAL-G", "CAL-2L"):
            P = np.array([n["params"][fam] for n in d["nets"]])
            med = np.median(P, axis=0)
            q1, q3 = np.percentile(P, [25, 75], axis=0)
            out.append(f"  {fam} params median [IQR]: " + ", ".join(f"{a:.3f} [{b:.3f}, {c:.3f}]"
                                                                    for a, b, c in zip(med, q1, q3)))
        ce = {k: np.median([n["corr_err"][k] for n in d["nets"]]) for k in d["nets"][0]["corr_err"]}
        out.append("  median |corrections - true statics| (s): " + ", ".join(f"{k} {v:.3f}" for k, v in ce.items()))
        rms = {k: np.median([n["calib"][k]["round1"]["rms"] for n in d["nets"]]) for k in ("CAL-H", "CAL-G", "CAL-2L")}
        out.append("  median JHD pick rms (s): " + ", ".join(f"{k} {v:.4f}" for k, v in rms.items()))
        sec = {}
        for r in rows:
            for k, v in r["sec"].items():
                sec.setdefault(k, []).append(v)
        out.append("  s/window: " + ", ".join(f"{k} {np.mean(v):.3f}" for k, v in sec.items()))
        out.append("  setup s/network (6 models incl. JHD): %.1f; JHD s: H %.2f G %.2f 2L %.2f" % (
            np.mean([n["t_setup"] for n in d["nets"]]),
            *[np.mean([n["calib"][k]["sec"] for n in d["nets"]]) for k in ("CAL-H", "CAL-G", "CAL-2L")]))
        block(rows, f"{reg} pooled (3 cells)", out)
    if len(sys.argv) > 2:
        out.append("=" * 118)
        for r in allrows:
            pass
        rows2 = []
        for path in sys.argv[1:]:
            d = load(path)
            reg = d["nets"][0]["regime"]
            for r in d["rows"]:
                r = dict(r)
                r["cell"] = reg + "/" + r["cell"]
                rows2.append(r)
        block(rows2, "BOTH REGIMES pooled (6 strata)", out, per_cell=False)
    print("\n".join(out))


if __name__ == "__main__":
    main()
