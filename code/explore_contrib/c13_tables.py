"""C13: narrow tables for REPORT_C13.txt from out/c13_summary.json  ->  out/c13_report_tables.txt"""
import json
import os

import c13_lib as L

OUT = os.path.join(L.HERE, "out")
DEC = tuple(json.load(open(os.path.join(OUT, "c13_summary.json")))["dec"])
DECT = tuple(d for d in DEC if d != "D1o")          # D1o is shown in TB, TC and TH only (same construction class as D3)
RULES = ("LOC+T", "LOC")
S = json.load(open(os.path.join(OUT, "c13_summary.json")))
cal = S["cal"]
LINES = list(cal)
STRATA = [n for n in LINES if not n.endswith("ALL")]
fh = open(os.path.join(OUT, "c13_report_tables.txt"), "w")


def g(name, bs, d, m, q, r, k):
    return cal[name]["cal"]["|".join(map(str, (bs, d, m, float(q), r)))][k]


def sh(n):
    return n.replace("unilateral_c40", "unil_c40")


def w(s=""):
    fh.write(s + "\n")


w("TA. Sample and benchmark (oracle) operating points of MT-FW: recall at <= 0.1 / <= 0.25 realised false events per")
w("    window, and the false discovery proportion (FDP) of the accepted atoms at those points")
w(f"  {'line':27s} {'win':>4s} {'ev':>5s}   LOC+T: rec@.1 FDP  rec@.25 FDP     LOC: rec@.1 FDP  rec@.25 FDP")
for n in LINES:
    b = cal[n]["bench"]
    w(f"  {n:27s} {cal[n]['nw']:4d} {cal[n]['ne']:5d}          {b[0][0][0]:.3f} {b[0][0][1]:.2f}   {b[0][1][0]:.3f} {b[0][1][1]:.2f}"
      f"          {b[1][0][0]:.3f} {b[1][0][1]:.2f}   {b[1][1][0]:.3f} {b[1][1][1]:.2f}")

w()
w("TB. Windows without events: targets against decoys (first 5 atoms, record scores). KS = Kolmogorov-Smirnov distance")
w("    [p by swapping target and decoy lists within windows]; r80, r95 = #targets / #decoys above the pooled 0.80 and")
w("    0.95 quantile [95% bootstrap]")
ex = S["exchange"]
for lab, title in (("CLUT", "clutter-only windows (CLUT + empty STD)"), ("NOIS", "noise-only windows (no clutter)")):
    w(f"  {title}")
    w(f"  {'regime cell':23s} {'n':>4s} dec   KS top [p]     KS all [p]     r80               r95")
    for reg in ("thin", "strong"):
        for cell in ("ladder_4", "ladder_6", "unilateral_c40"):
            for d in DEC:
                k = "|".join((lab, reg, cell, d))
                if k not in ex:
                    continue
                e = ex[k]
                w(f"  {reg + ' ' + cell:23s} {e['n']:4d} {d:4s}  {e['ks_top']:.3f} [{e['p_top']:.3f}]  {e['ks_all']:.3f} "
                  f"[{e['p_all']:.3f}]  {e['r80'][0]:.2f} [{e['r80'][1]:.2f},{e['r80'][2]:.2f}]  {e['r95'][0]:.2f} "
                  f"[{e['r95'][1]:.2f},{e['r95'][2]:.2f}]")

w()
w("TC. Windows with events: false targets against decoys. Atoms per event window above t80 | t95 (0.80 and 0.95 quantile")
w("    of the target scores of the clutter-only windows of the stratum). false = false targets under LOC+T, of which")
w("    tmiss (within 5 km, origin error > 0.5 s), near (5-10 km), far (> 10 km); fLOC = false under LOC; hit = atoms")
w("    labelled as hits; decoy columns = decoy atoms above the same thresholds")
ev = S["event_windows"]
for ti, tn in ((0, "t80"), (1, "t95")):
    w(f"  at {tn}")
    w(f"  {'line':27s} false tmiss  near   far  fLOC   hit | " + " ".join(f"{d:>5s}" for d in DEC))
    for n in STRATA:
        if n not in ev:
            continue
        e = ev[n]
        w(f"  {n:27s} " + " ".join(f"{e[k][ti][0]:5.2f}" for k in ("false", "tmiss", "near", "far", "falseLOC", "hit"))
          + " | " + " ".join(f"{e[d][ti][0]:5.2f}" for d in DEC))
w("  ratio false targets (LOC+T) / decoys at t80 [95%]; a valid decoy needs <= 1")
w(f"  {'line':22s}" + "".join(f"{d:>17s}" for d in DEC))
for n in STRATA:
    if n in ev:
        w(f"  {sh(n):22s}" + "".join(f" {ev[n]['false/' + d][0][0]:4.2f} [{ev[n]['false/' + d][0][1]:4.2f},"
                                 f"{min(ev[n]['false/' + d][0][2], 9.99):4.2f}]" for d in DEC))
w("  ratio false targets under LOC / decoys at t80 and t95 [95%], and KS distance of the false-target scores (LOC+T)")
w(f"  {'line':27s}  fLOC/D1p t80      fLOC/D1p t95      fLOC/D3 t80       KS " + " ".join(f"{d:>5s}" for d in DEC))
for n in STRATA:
    if n in ev:
        a, b = ev[n]["falseLOC/D1p"], ev[n]["falseLOC/D3"]
        w(f"  {n:27s}  {a[0][0]:4.2f} [{a[0][1]:4.2f},{a[0][2]:4.2f}]  {a[1][0]:4.2f} [{a[1][1]:4.2f},{min(a[1][2], 99):5.2f}]"
          f"  {b[0][0]:4.2f} [{b[0][1]:4.2f},{min(b[0][2], 99):5.2f}]   " + " ".join(f"{ev[n]['ks_false'][d]:.3f}" for d in DEC))

w()
w("TD. Event windows W against their clutter-and-noise twins W0 (same clutter and noise, no event): atoms per event")
w("    window above t80 | t95 by distance to the nearest event epicentre, and ratios [95%]")
tw = S.get("twin", {})
for n in tw:
    if not n.endswith("ALL"):
        continue
    t = tw[n]
    w(f"  {n}")
    for k in ("W_false", "W_in5", "W_near", "W_far", "W_far_pre", "W_far_post", "W0_all", "W0_in5", "W0_near", "W0_far"):
        w(f"      {k:11s} {t[k][0][0]:.3f} [{t[k][0][1]:.3f},{t[k][0][2]:.3f}] | {t[k][1][0]:.3f} [{t[k][1][1]:.3f},{t[k][1][2]:.3f}]")
    for k in ("W_false/W0_all", "W_far/W0_far", "W_near/W0_near", "W_in5/W0_in5", "W_far_post/W0_far",
              "k1h_far_post/k1h_W0_far"):
        w(f"      ratio {k:24s} {t[k][0][0]:.2f} [{t[k][0][1]:.2f},{t[k][0][2]:.2f}] | {t[k][1][0]:.2f} [{t[k][1][1]:.2f},"
          f"{min(t[k][1][2], 99):.2f}]")
w("  per stratum, ratio W_far/W0_far and W_false/W0_all at t80 [95%]")
for n in tw:
    if n.endswith("ALL"):
        continue
    a, b = tw[n]["W_far/W0_far"][0], tw[n]["W_false/W0_all"][0]
    w(f"      {n:27s} {a[0]:.2f} [{a[1]:.2f},{a[2]:.2f}]   {b[0]:.2f} [{b[1]:.2f},{b[2]:.2f}]")

w()
w("TE. Target-decoy threshold (tdc): realised FDR = mean FDP over batches [95% hierarchical bootstrap, percentile].")
w("    * = interval entirely above the nominal q")
for r in RULES:
    for bs in (120, 20):
        for q in (0.05, 0.1, 0.2):
            w(f"  rule {r}, batch {bs} windows, q = {q}")
            w(f"  {'line':22s}" + "".join(f"{d:>19s}" for d in DECT))
            for n in LINES:
                s = f"  {sh(n):22s}"
                for d in DECT:
                    v = g(n, bs, d, "tdc", q, r, "fdr")
                    s += f"  {v[0]:.2f} [{v[1]:.2f},{v[2]:.2f}]{'*' if v[1] > q else ' '}"
                w(s)
w()
w("TE2. The same FDR for batches of 20 windows with the interval mean +- 1.96 standard errors over the independent")
w("    batches (one per network; 12 per stratum, 36 or 72 on the pooled lines). * = interval entirely above q")
for r in RULES:
    for q in (0.05, 0.1, 0.2):
        w(f"  rule {r}, batch 20 windows, q = {q}")
        w(f"  {'line':22s}" + "".join(f"{d:>19s}" for d in DECT))
        for n in LINES:
            s = f"  {sh(n):22s}"
            for d in DECT:
                v, se = g(n, 20, d, "tdc", q, r, "fdr")[0], g(n, 20, d, "tdc", q, r, "se20")[0]
                lo, hi = max(0.0, v - 1.96 * se), min(1.0, v + 1.96 * se)
                s += f"  {v:.2f} [{lo:.2f},{hi:.2f}]{'*' if lo > q else ' '}"
            w(s)

w()
w("TF. tdc: recall | realised false events per window | recall of the oracle threshold at the same realised rate |")
w("    loss = oracle - tdc [95%] | share of batches that accept anything")
for d in ("D1p", "D2r", "D2s"):
    for r in RULES:
        for bs in (120, 20):
            w(f"  decoy {d}, rule {r}, batch {bs}")
            w(f"  {'line':22s}" + "".join(f"{'q = ' + str(q):>44s}" for q in (0.1, 0.2)))
            for n in LINES:
                s = f"  {sh(n):22s}"
                for q in (0.1, 0.2):
                    a, b, c, e, y = (g(n, bs, d, "tdc", q, r, k) for k in ("rec", "few", "orec", "loss", "any"))
                    s += f"  {a[0]:.3f}|{b[0]:.3f}|{c[0]:.3f}|{e[0]:+.3f} [{e[1]:+.3f},{e[2]:+.3f}]|{y[0]:.2f}"
                w(s)

w()
w("TF0. tdc with D1p at q = 0.05 and q = 0.5: recall | realised false events per window | loss against the oracle")
w("    threshold at the same realised rate [95%] | share of batches that accept anything")
for r in RULES:
    for bs in (120, 20):
        w(f"  rule {r}, batch {bs}")
        w(f"  {'line':22s}" + "".join(f"{'q = ' + str(q):>44s}" for q in (0.05, 0.5)))
        for n in LINES:
            s = f"  {sh(n):22s}"
            for q in (0.05, 0.5):
                a, b, c, e, y = (g(n, bs, "D1p", "tdc", q, r, k) for k in ("rec", "few", "orec", "loss", "any"))
                s += f"  {a[0]:.3f}|{b[0]:.3f}|{c[0]:.3f}|{e[0]:+.3f} [{e[1]:+.3f},{e[2]:+.3f}]|{y[0]:.2f}"
            w(s)

w()
w("TG. fpw (threshold = decoy score exceeded by at most c x windows decoys): realised false events per window [95%]")
w("    (* = interval entirely above c) | recall | recall of the benchmark oracle at <= c | loss = oracle - fpw [95%]")
for d in DEC:
    for r in RULES:
        for bs in (120, 20):
            w(f"  decoy {d}, rule {r}, batch {bs}")
            w(f"  {'line':22s}" + "".join(f"{'c = ' + str(c):>52s}" for c in (0.1, 0.25)))
            for n in LINES:
                s = f"  {sh(n):22s}"
                for ci, c in ((0, 0.1), (1, 0.25)):
                    a, b, e = (g(n, bs, d, "fpw", c, r, k) for k in ("few", "rec", "lossn"))
                    ob = cal[n]["bench"][RULES.index(r)][ci][0]
                    s += f"  {a[0]:.2f} [{a[1]:.2f},{a[2]:.2f}]{'*' if a[1] > c else ' '}|{b[0]:.3f}|{ob:.3f}|{e[0]:+.3f} [{e[1]:+.2f},{e[2]:+.2f}]"
                w(s)

w()
w("TH. Variants on the pooled lines: FDR [95%] of tdc, bh (conformal p-value + Benjamini-Hochberg) and comp (rank-paired")
w("    competition), batch 120 and 20")
for r in RULES:
    for n in [x for x in LINES if x.endswith("both ALL")]:
        for bs in (120, 20):
            w(f"  {n}, rule {r}, batch {bs}")
            w(f"    {'q':>5s} {'meth':5s}" + "".join(f"{d:>19s}" for d in DECT))
            for q in (0.05, 0.1, 0.2, 0.3, 0.5):
                for m in ("tdc", "bh", "comp"):
                    s = f"    {q:5.2f} {m:5s}"
                    for d in DECT:
                        v = g(n, bs, d, m, q, r, "fdr")
                        s += f"  {v[0]:.2f} [{v[1]:.2f},{v[2]:.2f}]{'*' if v[1] > q else ' '}"
                    w(s)
w()
w("TI. Cost: seconds per window (target run | five decoy runs)")
for n, v in S["sec"].items():
    w(f"  {n:27s} {v[0]:4d} windows  {v[2]:.3f} | {v[3]:.3f}")
fh.close()
