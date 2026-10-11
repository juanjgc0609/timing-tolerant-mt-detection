"""Task 2M analysis (C16 endpoints and bootstrap; decision rules of the 2M registration).
    python t2m_analyze.py OUT.txt B LABEL=PARTSDIR [LABEL=PARTSDIR ...]
One LABEL=PARTSDIR per seed (e.g. 2718=t2m_out/parts_2718 3141=t2m_out/parts_3141).  B = bootstrap replicates (2000).
Endpoints exactly as C16 (c16_analyze.py, c10_analyze.Boot): exact FROC, recall at <= 0.1 and <= 0.25 realised false
events per window, operating point re-chosen in every replicate; hit rules LOC (epicentre <= 5 km) and LOC+T (and
|origin error| <= 0.5 s); hierarchical bootstrap, strata = cell, networks then windows; paired differences.
Registered readings (per seed, REAL arm, cells pooled):
  2M-a  DENSE, LOC, @0.25: MT@S0 - BOMP@S0, 2.5th percentile (one-sided 97.5 lower bound) > 0 and point >= +0.03.
  2M-b  DENSE, LOC, @0.25: HUBER-N@S0 - BOMP@S0, the same.
  2M-c  secondary: the same two contrasts in DENSE LOC+T and in STD LOC (lower bound > 0).
  2M-d  non-inferiority, DENSE, LOC, @0.25: HUBER-N@S0 - MT@S0, lower bound > -0.03.
  Strict reading: Holm over the four primary tests (2M-a, 2M-b x seeds), one-sided alpha 0.025, bootstrap p-value
  p = (1 + #{d <= 0}) / (B + 1).
Descriptive: @0.1, per cell, recall tables, MT@S2+T3 - MT@S0+T3, noise statistics, station terms, transfer of the
operating point to noise-only windows, REAL - GEN if GEN parts exist, runtime.
"""
import os
import sys

import numpy as np

WORK = os.path.dirname(os.path.abspath(__file__))
EXPL = os.path.join(WORK, "explore_contrib")
for p in (EXPL, WORK):
    if p not in sys.path:
        sys.path.insert(0, p)

import c10_analyze as A  # noqa: E402
import c16_analyze as Z  # noqa: E402

CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
PRIMARY = {"2M-a": ("MT@S0", "BOMP@S0"), "2M-b": ("HUBER-N@S0", "BOMP@S0")}
NI = ("HUBER-N@S0", "MT@S0")
NI_MARGIN = -0.03
POINT_MIN = 0.03
ALPHA = 0.025
DESCR = [("MT@S2+T3", "MT@S0+T3"), ("MT@S0+T3", "MT@S0"), ("MT@S2+T3", "MT@S0")]


def holm(pvals, alpha):
    keys = sorted(pvals, key=lambda k_: pvals[k_])
    m = len(keys)
    adj, run, rej, stop = {}, 0.0, {}, False
    for i, k_ in enumerate(keys):
        run = max(run, min(1.0, (m - i) * pvals[k_]))
        adj[k_] = run
        if not stop and pvals[k_] <= alpha / (m - i):
            rej[k_] = True
        else:
            stop = True
            rej[k_] = False
    return rej, adj


def contrast(bt, a, b, k=1):
    d = bt.rep[a][:, k] - bt.rep[b][:, k]
    p = bt.pt[a][k] - bt.pt[b][k]
    return {"point": float(p), "lo": float(np.percentile(d, 2.5)), "hi": float(np.percentile(d, 97.5)),
            "sd": float(d.std()), "p_le0": float((1 + np.sum(d <= 0)) / (len(d) + 1)),
            "p_le_ni": float((1 + np.sum(d <= NI_MARGIN)) / (len(d) + 1))}


def f3(c):
    return f"{c['point']:+.3f} [{c['lo']:+.3f}, {c['hi']:+.3f}] (sd {c['sd']:.4f}, p {c['p_le0']:.4f})"


def main(outp, B, specs):
    lines = []

    def P(s=""):
        print(s, flush=True)
        lines.append(s)
    seeds = [s.split("=", 1) for s in specs]
    res = {}
    P(f"TASK 2M ANALYSIS: seeds {[s for s, _ in seeds]}; bootstrap B = {B}; REAL arm primary")
    for lab, pdir in seeds:
        data = {arm: Z.load(pdir, arm) for arm in ("REAL", "GEN")}
        nets, rows = data["REAL"]
        P("")
        P("=" * 118)
        P(f"SEED {lab} ({pdir}): REAL {len(nets)} parts (" + ", ".join(f"{c} {sum(n['cell'] == c for n in nets)}"
                                                                     for c in CELLS) + "); days " +
          ", ".join(f"{d} {sum(n.get('day') == d for n in nets)}" for d in sorted({n.get('day') for n in nets})) +
          f"; windows STD {sum(r['design'] == 'STD' for r in rows)} ({sum(r['k'] for r in rows if r['design'] == 'STD')} "
          f"events), DENSE {sum(r['design'] == 'DENSE' for r in rows)} ({sum(r['k'] for r in rows if r['design'] == 'DENSE')}"
          f" events), noise-only {sum(r['design'] == 'NOISE' for r in rows)}; GEN parts {len(data['GEN'][0])}")
        names = list(rows[0]["atoms"])
        res[lab] = {}
        for des, rule, tol in (("DENSE", "LOC", None), ("DENSE", "LOC+T", 0.5), ("STD", "LOC", None),
                               ("STD", "LOC+T", 0.5)):
            sub = [r for r in rows if r["design"] == des]
            if not sub:
                continue
            recs = Z.records(sub, names, tol)
            for grp in ("all",) + CELLS:
                sel = (lambda r: True) if grp == "all" else (lambda r, g=grp: r["cell"] == g)
                bt = A.Boot(sub, recs, sel, B if grp == "all" else max(500, B // 2))
                P("")
                P(f"== {lab} REAL {des} rule {rule} cells {grp}: {len(bt.idx)} windows, {int(bt.ks.sum())} events")
                if grp == "all":
                    for m in bt.methods:
                        c = bt.ci(m)
                        P(f"   {m:12s} {bt.pt[m][0]:.3f} [{c[0][0]:.3f},{c[0][1]:.3f}]  {bt.pt[m][1]:.3f} "
                          f"[{c[1][0]:.3f},{c[1][1]:.3f}]")
                for a, b in list(PRIMARY.values()) + [NI] + DESCR:
                    if a in names and b in names:
                        c1, c0 = contrast(bt, a, b, 1), contrast(bt, a, b, 0)
                        P(f"   {a:12s} - {b:12s} @0.1 {c0['point']:+.3f} [{c0['lo']:+.3f},{c0['hi']:+.3f}]  @0.25 {f3(c1)}")
                        if grp == "all":
                            res[lab][(des, rule, a, b)] = c1
        # noise and terms (descriptive)
        P("")
        lev = np.array([x for r in rows for x in r["lev"]])
        ku = np.array([r["kurt"] for r in rows])
        mx = np.array([r["maxabs"] for r in rows])
        snr = {d: np.array([x for r in rows if r["design"] == d for x in r["snr"]]) for d in ("STD", "DENSE")}
        P(f"   noise level / generator (5/25/50/75/95%) " + " ".join(f"{x:.2f}" for x in np.percentile(lev, [5, 25, 50, 75, 95]))
          + f"; excess kurtosis median {np.median(ku):.2f} p90 {np.percentile(ku, 90):.2f}; max|MF|/sd median "
          f"{np.median(mx):.1f} p90 {np.percentile(mx, 90):.1f}")
        for d in ("STD", "DENSE"):
            if len(snr[d]):
                P(f"   SNR per arrival {d}: quantiles 10/25/50/75/90 " +
                  " ".join(f"{x:.2f}" for x in np.percentile(snr[d], [10, 25, 50, 75, 90])))
        for cell in CELLS + ("all",):
            ns_ = [n for n in nets if cell == "all" or n["cell"] == cell]
            e0 = np.array([n["err"]["S0"] for n in ns_])
            e2 = np.array([n["err"]["S2"] for n in ns_])
            P(f"   terms {cell:15s} S0 {np.median(e0):.3f}  S2 {np.median(e2):.3f}  better {np.mean(e2 < e0):.2f}  "
              f"catalogue drawn mean {np.mean([n['cat_drawn'] for n in ns_]):.1f} max {max(n['cat_drawn'] for n in ns_)}")
        no = [r for r in rows if r["design"] == "NOISE"]
        if no:
            P("   TRANSFER (descriptive): threshold at <= 0.1 / <= 0.25 on DENSE LOC event windows; realised false events "
              "per noise-only window [95%]")
            ev = [r for r in rows if r["design"] == "DENSE"]
            tr = Z.transfer(ev, no, names, None, 7, max(200, B // 4), 4242)
            for nm in names:
                P(f"     {nm:12s} " + "  ".join(f"@{cap}: t {tr[nm][i][0][0]:6.2f} noise {tr[nm][i][0][1]:.3f} "
                                                 f"[{tr[nm][i][1][0][0]:.3f},{tr[nm][i][1][1][0]:.3f}]"
                                                 for i, cap in enumerate(A.CAPS)))
        if data["GEN"][0]:
            rows_p = Z.combine(rows, data["GEN"][1], "R", "G")
            sub = [r for r in rows_p if r["design"] == "DENSE"]
            if sub:
                nm2 = list(sub[0]["atoms"])
                bt = A.Boot(sub, Z.records(sub, nm2, None), lambda r: True, B)
                P(f"   REAL - GEN (DENSE, LOC, joint bootstrap, {len(bt.idx)} windows):")
                for a, b in PRIMARY.values():
                    d = (bt.rep["R:" + a] - bt.rep["R:" + b]) - (bt.rep["G:" + a] - bt.rep["G:" + b])
                    p = (bt.pt["R:" + a][1] - bt.pt["R:" + b][1]) - (bt.pt["G:" + a][1] - bt.pt["G:" + b][1])
                    gp = bt.pt["G:" + a][1] - bt.pt["G:" + b][1]
                    P(f"     {a} - {b}: GEN {gp:+.3f}; REAL minus GEN @0.25 {p:+.3f} "
                      f"[{np.percentile(d[:, 1], 2.5):+.3f},{np.percentile(d[:, 1], 97.5):+.3f}]")
        P(f"   runtime s per window: " + "; ".join(f"{m} {np.mean([r['sec'][m] for r in rows]):.2f}" for m in names) +
          f"; part mean {np.mean([n['t_total'] for n in nets]):.0f} s")
    # ------------------------------------------------------------------ decision
    P("")
    P("=" * 118)
    P("DECISION (registered rules; REAL arm, cells pooled, @0.25)")
    pv = {}
    verdict = {}
    for rule_id, (a, b) in PRIMARY.items():
        ok_all = True
        for lab, _ in seeds:
            c = res[lab].get(("DENSE", "LOC", a, b))
            ok = c is not None and c["lo"] > 0 and c["point"] >= POINT_MIN
            ok_all &= ok
            pv[(rule_id, lab)] = c["p_le0"] if c else 1.0
            P(f"   {rule_id} seed {lab}: {a} - {b} DENSE LOC {f3(c)} -> {'MET' if ok else 'NOT MET'}")
        verdict[rule_id] = ok_all and len(seeds) >= 2
        P(f"   {rule_id}: {'MET in all seeds' if verdict[rule_id] else 'NOT MET'}")
    rej, adj = holm(pv, ALPHA)
    P("   Holm (4 tests, one-sided 0.025): " + "; ".join(f"{k[0]} {k[1]} p {pv[k]:.4f} adj {adj[k]:.4f} "
                                                       f"{'reject' if rej[k] else 'keep'}" for k in sorted(pv)))
    for lab, _ in seeds:
        for (des, rule) in (("DENSE", "LOC+T"), ("STD", "LOC")):
            for a, b in PRIMARY.values():
                c = res[lab].get((des, rule, a, b))
                if c:
                    P(f"   2M-c seed {lab} {des} {rule}: {a} - {b} {f3(c)} -> {'lower bound > 0' if c['lo'] > 0 else 'not'}")
        c = res[lab].get(("DENSE", "LOC") + NI)
        if c:
            P(f"   2M-d seed {lab}: {NI[0]} - {NI[1]} DENSE LOC {f3(c)} -> "
              f"{'non-inferior (lower bound > -0.03)' if c['lo'] > NI_MARGIN else 'non-inferiority not shown'}")
    P(f"   G1-R under real noise confirmed (2M-a and 2M-b in every seed): {verdict.get('2M-a') and verdict.get('2M-b')}")
    with open(outp, "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3:])
