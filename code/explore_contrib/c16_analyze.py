"""C16 analysis.  SEED 999 ONLY.
    python c16_analyze.py PARTSDIR B OUT.txt
Endpoints as C11 / C14: exact FROC, recall at <= 0.1 and <= 0.25 realised false events per window (operating point
re-chosen in every replicate), hit rules LOC+T (epicentre <= 5 km and |origin error| <= 0.5 s) and LOC (<= 5 km),
hierarchical bootstrap (strata = cell; networks, then windows; c10_analyze.Boot), 95% percentile intervals, paired
differences.  Arms REAL and GEN share networks and signals: the difference between arms of a contrast is bootstrapped
jointly (same resampled networks and windows).
Primary readings (fixed in C16_FREEZE.txt):
  G1-R  DENSE, LOC, @0.25: MT@S0 - BOMP@S0 and HUBER-N@S0 - BOMP@S0 (REAL arm); LOC+T and @0.1 secondary.
  G1-C  STD, LOC+T, @0.25: MT@S2+T3 - MT@S0+T3 (REAL arm); LOC, @0.1, DENSE secondary; oracle MT@ORA+T3.
  TRANSFER  per method and design, the score threshold of the operating point (<= 0.1 / <= 0.25 false events per
        window) on the event windows of an arm (DENSE: rule LOC; STD: rule LOC+T), and the realised false events per
        window at that threshold on the arm's noise-only windows (STD: first 5 atoms; DENSE: 7).
"""
import glob
import json
import os
import sys

import numpy as np

import c10_analyze as A
import c10_lib as C

CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
CAPS = A.CAPS
G1R = [("MT@S0", "BOMP@S0"), ("HUBER-N@S0", "BOMP@S0"), ("HUBER-N@S0", "MT@S0")]
G1C = [("MT@S2+T3", "MT@S0+T3"), ("MT@ORA+T3", "MT@S0+T3"), ("MT@S2+T3", "MT@S0"), ("MT@S0+T3", "MT@S0"),
       ("MT@S2", "MT@S0"), ("MT@ORA+T3", "MT@S2+T3")]
RULES = (("LOC+T", 0.5), ("LOC", None))


def load(pdir, arm):
    parts = sorted(glob.glob(os.path.join(pdir, f"{arm}_*.json")))
    nets = [json.load(open(p)) for p in parts]
    rows = []
    for n in nets:
        for r in n["rows"]:
            r = dict(r)
            r.update({"regime": arm, "cell": n["cell"], "net": n["net"], "ceil": []})
            rows.append(r)
    return nets, rows


def combine(rows_a, rows_b, la, lb):
    key = {(r["cell"], r["net"], r["design"], r["w"]): r for r in rows_b}
    out = []
    for r in rows_a:
        q = key.get((r["cell"], r["net"], r["design"], r["w"]))
        if q is None:
            continue
        assert r["k"] == q["k"] and np.allclose(np.array(r["truth"]).reshape(-1), np.array(q["truth"]).reshape(-1))
        at = {f"{la}:{m}": v for m, v in r["atoms"].items()}
        at.update({f"{lb}:{m}": v for m, v in q["atoms"].items()})
        out.append({"regime": "pair", "cell": r["cell"], "net": r["net"], "design": r["design"], "w": r["w"],
                    "k": r["k"], "truth": r["truth"], "taus": r["taus"], "atoms": at})
    return out


def records(rows, names, tol_t, kcut=None):
    recs = {}
    for nm in names:
        out = []
        for r in rows:
            at = A.atoms_of(r, nm)[:kcut]
            out.append(C.record(at, C.label(at, [np.array(t) for t in r["truth"]], r["taus"], tol_t), "s"))
        recs[nm] = out
    return recs


def fmt(p, ci):
    return f"{p:+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}]"


def thr_at(recs_w, nwin, cap):
    """score threshold of the exact-FROC operating point at <= cap false events per window (lowest score level whose
    realised false count per window is <= cap; +inf if none)."""
    sc = np.concatenate([np.array([x[0] for x in r]) for r in recs_w] + [np.zeros(0)])
    hit = np.concatenate([np.array([x[1] for x in r], int) for r in recs_w] + [np.zeros(0, int)])
    if not len(sc):
        return np.inf
    o = np.argsort(-sc, kind="stable")
    sc, hit = sc[o], hit[o]
    Fc = np.cumsum(hit < 0)
    last = np.r_[np.where(np.diff(sc) != 0)[0], len(sc) - 1]
    ok = Fc[last] / nwin <= cap
    return float(sc[last][ok].min()) if ok.any() else np.inf


def fp_noise(recs_w, thr):
    return float(np.mean([sum(1 for s, _h in r if s >= thr) for r in recs_w])) if len(recs_w) else np.nan


def transfer(ev_rows, no_rows, names, tol, kcut, B, seed, cross_rows=None):
    """point and bootstrap (networks within cell, then windows; same network draw for event and noise windows)."""
    rec_e = records(ev_rows, names, tol)
    rec_n = records(no_rows, names, None, kcut)
    rec_c = records(cross_rows, names, None, kcut) if cross_rows is not None else None
    ge, gn = {}, {}
    for i, r in enumerate(ev_rows):
        ge.setdefault((r["cell"], r["net"]), []).append(i)
    for i, r in enumerate(no_rows):
        gn.setdefault((r["cell"], r["net"]), []).append(i)
    gc = {}
    if cross_rows is not None:
        for i, r in enumerate(cross_rows):
            gc.setdefault((r["cell"], r["net"]), []).append(i)
    keys = sorted(set(ge) & set(gn))
    rng = np.random.default_rng(seed)
    res = {}
    draws = []
    for _b in range(B):
        we, wn, wc = [], [], []
        for cell in CELLS:
            kc = [k for k in keys if k[0] == cell]
            for k in [kc[i] for i in rng.integers(0, len(kc), len(kc))]:
                e, n_ = np.asarray(ge[k]), np.asarray(gn[k])
                we += e[rng.integers(0, len(e), len(e))].tolist()
                wn += n_[rng.integers(0, len(n_), len(n_))].tolist()
                if gc:
                    c_ = np.asarray(gc[k])
                    wc += c_[rng.integers(0, len(c_), len(c_))].tolist()
        draws.append((we, wn, wc))
    for nm in names:
        out = []
        for cap in CAPS:
            t0 = thr_at(rec_e[nm], len(ev_rows), cap)
            p = (t0, fp_noise(rec_n[nm], t0), fp_noise(rec_c[nm], t0) if rec_c else np.nan)
            bs = []
            for we, wn, wc in draws:
                t = thr_at([rec_e[nm][i] for i in we], len(we), cap)
                bs.append((fp_noise([rec_n[nm][i] for i in wn], t),
                           fp_noise([rec_c[nm][i] for i in wc], t) if rec_c else np.nan))
            bs = np.array(bs)
            out.append((p, np.nanpercentile(bs, [2.5, 97.5], axis=0)))
        res[nm] = out
    return res


def main(pdir, B=1000, outp="out/c16_analysis.txt"):
    lines = []

    def P(s=""):
        print(s, flush=True)
        lines.append(s)
    data = {arm: load(pdir, arm) for arm in ("REAL", "GEN")}
    arms = [a for a in ("REAL", "GEN") if data[a][0]]
    P(f"C16 ANALYSIS ({pdir}): arms {arms}; bootstrap B = {B}")
    for arm in arms:
        nets, rows = data[arm]
        P(f"  {arm}: {len(nets)} networks (" + ", ".join(f"{c} {sum(n['cell'] == c for n in nets)}" for c in CELLS) +
          f"); windows STD {sum(r['design'] == 'STD' for r in rows)} ({sum(r['k'] for r in rows if r['design'] == 'STD')}"
          f" events), DENSE {sum(r['design'] == 'DENSE' for r in rows)} ({sum(r['k'] for r in rows if r['design'] == 'DENSE')}"
          f" events), noise-only {sum(r['design'] == 'NOISE' for r in rows)}")
    # ------------------------------------------------------------------ SNR and noise
    P("")
    P("== SNR per true arrival (peak 3-component |MF(signal)| within +-0.25 s of the arrival / window noise level sd of "
      "the detector), STD and DENSE event windows; quantiles 10/25/50/75/90; share > 3 and > 6")
    snr = {}
    for arm in arms:
        _n, rows = data[arm]
        for des in ("STD", "DENSE"):
            v = np.array([x for r in rows if r["design"] == des for x in r["snr"]])
            snr[(arm, des)] = v
            q = np.percentile(v, [10, 25, 50, 75, 90])
            P(f"   {arm:4s} {des:5s} n {len(v):6d}  " + " ".join(f"{x:6.2f}" for x in q) +
              f"   >3 {np.mean(v > 3):.3f}  >6 {np.mean(v > 6):.3f}")
    if len(arms) == 2:
        for des in ("STD", "DENSE"):
            a, b = snr[("REAL", des)], snr[("GEN", des)]
            if len(a) == len(b):
                r_ = np.log(a / b)
                P(f"   paired ratio REAL / GEN of the same arrivals, {des}: median {np.exp(np.median(r_)):.3f}, quartiles "
                  f"{np.exp(np.percentile(r_, 25)):.3f} - {np.exp(np.percentile(r_, 75)):.3f}, 5-95% "
                  f"{np.exp(np.percentile(r_, 5)):.3f} - {np.exp(np.percentile(r_, 95)):.3f}")
    P("   noise part of the windows: per-station robust in-band level / generator level (quantiles 5/25/50/75/95), "
      "excess kurtosis of the normalised MF noise (median, 90%), max |MF| / robust sd (median, 90%), window sd (median)")
    for arm in arms:
        _n, rows = data[arm]
        lev = np.array([x for r in rows for x in r["lev"]])
        ku = np.array([r["kurt"] for r in rows])
        mx = np.array([r["maxabs"] for r in rows])
        sd = np.array([r["sd"] for r in rows])
        P(f"   {arm:4s} level " + " ".join(f"{x:.2f}" for x in np.percentile(lev, [5, 25, 50, 75, 95])) +
          f"; kurtosis {np.median(ku):.2f} / {np.percentile(ku, 90):.2f}; max {np.median(mx):.1f} / "
          f"{np.percentile(mx, 90):.1f}; sd {np.median(sd):.3f}")
    # ------------------------------------------------------------------ terms
    P("")
    P("== station terms: rms error against the true statics, median removed (s), median over networks [quartiles]; "
      "better = share of networks where S2 < S0; catalogue candidates drawn for 20 events (mean)")
    for arm in arms:
        nets, _r = data[arm]
        for cell in CELLS + ("all",):
            ns_ = [n for n in nets if cell == "all" or n["cell"] == cell]
            if not ns_:
                continue
            e0 = np.array([n["err"]["S0"] for n in ns_])
            e2 = np.array([n["err"]["S2"] for n in ns_])
            P(f"   {arm:4s} {cell:15s} S0 {np.median(e0):.3f} [{np.percentile(e0, 25):.3f},{np.percentile(e0, 75):.3f}]  "
              f"S2 {np.median(e2):.3f} [{np.percentile(e2, 25):.3f},{np.percentile(e2, 75):.3f}]  better {np.mean(e2 < e0):.2f}"
              f"  drawn {np.mean([n['cat_drawn'] for n in ns_]):.1f}  T3 level S0 "
              f"{np.median([n['sig']['S0'][0] for n in ns_]):.3f} S2 {np.median([n['sig']['S2'][0] for n in ns_]):.3f}")
    # ------------------------------------------------------------------ recall and contrasts per arm
    for arm in arms:
        _n, rows = data[arm]
        names = list(rows[0]["atoms"])
        for des, pairs in (("DENSE", G1R + G1C), ("STD", G1C + G1R)):
            sub = [r for r in rows if r["design"] == des]
            for rule, tol in RULES:
                recs = records(sub, names, tol)
                for grp in ("all",) + CELLS:
                    sel = (lambda r: True) if grp == "all" else (lambda r, g=grp: r["cell"] == g)
                    if not any(sel(r) for r in sub):
                        continue
                    bt = A.Boot(sub, recs, sel, B if grp == "all" else max(200, B // 2))
                    P("")
                    P(f"== {arm} {des} rule {rule} cells {grp}: {len(bt.idx)} windows, {int(bt.ks.sum())} events")
                    if grp == "all":
                        for m in bt.methods:
                            c = bt.ci(m)
                            P(f"   {m:12s} {bt.pt[m][0]:.3f} [{c[0][0]:.3f},{c[0][1]:.3f}]  {bt.pt[m][1]:.3f} "
                              f"[{c[1][0]:.3f},{c[1][1]:.3f}]")
                    for a, b in pairs:
                        p, c, sd = bt.diff(a, b)
                        P(f"   {a:12s} - {b:12s} @0.1 {fmt(p[0], c[0])}  @0.25 {fmt(p[1], c[1])}  (sd {sd[1]:.4f})")
                    if grp == "all" and des == "STD" or (grp == "all" and des == "DENSE"):
                        den = bt.rep["MT@ORA+T3"][:, 1] - bt.rep["MT@S0+T3"][:, 1]
                        num = bt.rep["MT@S2+T3"][:, 1] - bt.rep["MT@S0+T3"][:, 1]
                        d0 = bt.pt["MT@ORA+T3"][1] - bt.pt["MT@S0+T3"][1]
                        f0 = (bt.pt["MT@S2+T3"][1] - bt.pt["MT@S0+T3"][1]) / d0 if abs(d0) > 1e-9 else np.nan
                        fr = num / np.where(np.abs(den) < 1e-9, np.nan, den)
                        P(f"   fraction of the oracle gain (S2+T3 - S0+T3) / (ORA+T3 - S0+T3) @0.25: {f0:+.2f} "
                          f"[{np.nanpercentile(fr, 2.5):+.2f},{np.nanpercentile(fr, 97.5):+.2f}] (oracle gain <= 0 in "
                          f"{np.mean(den <= 0):.3f} of replicates)")
    # ------------------------------------------------------------------ REAL vs GEN, jointly bootstrapped
    if len(arms) == 2:
        rows_p = combine(data["REAL"][1], data["GEN"][1], "R", "G")
        for des, pairs in (("DENSE", G1R), ("STD", G1C[:2]), ("DENSE", G1C[:2]), ("STD", G1R)):
            sub = [r for r in rows_p if r["design"] == des]
            names = list(sub[0]["atoms"])
            for rule, tol in RULES:
                recs = records(sub, names, tol)
                bt = A.Boot(sub, recs, lambda r: True, B)
                P("")
                P(f"== REAL - GEN (same networks and signals) {des} rule {rule}: {len(bt.idx)} windows")
                for m in [x[2:] for x in names if x.startswith("R:")]:
                    p, c, _sd = bt.diff("R:" + m, "G:" + m)
                    P(f"   recall {m:12s} REAL - GEN @0.1 {fmt(p[0], c[0])}  @0.25 {fmt(p[1], c[1])}")
                for a, b in pairs:
                    d = (bt.rep["R:" + a] - bt.rep["R:" + b]) - (bt.rep["G:" + a] - bt.rep["G:" + b])
                    p = (np.subtract(bt.pt["R:" + a], bt.pt["R:" + b]) - np.subtract(bt.pt["G:" + a], bt.pt["G:" + b]))
                    c = np.percentile(d, [2.5, 97.5], axis=0).T
                    P(f"   contrast {a} - {b}: REAL minus GEN @0.1 {fmt(p[0], c[0])}  @0.25 {fmt(p[1], c[1])}")
    # ------------------------------------------------------------------ transfer of the operating point
    P("")
    P("== TRANSFER of the operating point: threshold t fixed on the arm's event windows (DENSE: LOC; STD: LOC+T; and "
      "STD: LOC) at <= 0.1 / <= 0.25 false events per window; realised false events per window at t on the arm's "
      "noise-only windows [95%]; X = REAL noise-only windows at the GEN threshold (synthetic operating point on real noise)")
    for des, kcut, rule, tol in (("DENSE", 7, "LOC", None), ("STD", 5, "LOC+T", 0.5), ("STD", 5, "LOC", None)):
        for arm in arms:
            _n, rows = data[arm]
            ev = [r for r in rows if r["design"] == des]
            no = [r for r in rows if r["design"] == "NOISE"]
            cross = [r for r in data["REAL"][1] if r["design"] == "NOISE"] if (arm == "GEN" and "REAL" in arms) else None
            names = list(ev[0]["atoms"])
            res = transfer(ev, no, names, tol, kcut, max(200, B // 2), 4242, cross)
            P(f"   {arm} {des} ({rule}; noise-only windows {len(no)}, first {kcut} atoms)")
            for nm in names:
                s = f"     {nm:12s}"
                for ci_, cap in enumerate(CAPS):
                    (t0, f0, x0), q = res[nm][ci_]
                    s += f"  @{cap}: t {t0:6.2f} noise {f0:.3f} [{q[0][0]:.3f},{q[1][0]:.3f}]"
                    if cross is not None:
                        s += f" X {x0:.3f} [{q[0][1]:.3f},{q[1][1]:.3f}]"
                P(s)
    # zero-event STD windows of the event stream at the same thresholds (internal reference)
    P("   reference: false events per window of the zero-event STD windows of the event stream at the same STD "
      "threshold (LOC+T), point")
    for arm in arms:
        _n, rows = data[arm]
        ev = [r for r in rows if r["design"] == "STD"]
        z = [r for r in ev if r["k"] == 0]
        names = list(ev[0]["atoms"])
        re_ = records(ev, names, 0.5)
        rz = records(z, names, None, 5)
        P(f"     {arm} ({len(z)} zero-event windows): " + "; ".join(
            f"{nm} " + "/".join(f"{fp_noise(rz[nm], thr_at(re_[nm], len(ev), cap)):.3f}" for cap in CAPS) for nm in names))
    # ------------------------------------------------------------------ runtime
    P("")
    for arm in arms:
        nets, rows = data[arm]
        names = list(rows[0]["atoms"])
        P(f"== runtime {arm}, s per window (mean over designs): " + "; ".join(
            f"{m} {np.mean([r['sec'][m] for r in rows]):.2f}" for m in names) +
          f"; catalogue {np.mean([n['t_cat'] for n in nets]):.1f} s, setup {np.mean([n['t_setup'] for n in nets]):.1f} s, "
          f"network total {np.mean([n['t_total'] for n in nets]):.0f} s")
    with open(outp, "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 1000, sys.argv[3] if len(sys.argv) > 3 else
         "out/c16_analysis.txt")
