"""C14 analysis.   python c14_analyze.py dev|dev2|eval [B]   -> out/c14_analysis_<mode>.txt
Exact FROC (recall at realised false events per window <= 0.1 and <= 0.25, operating point re-chosen in every
bootstrap replicate), hit rules LOC+T (primary) and LOC, hierarchical bootstrap (strata regime x cell; networks, then
windows), paired differences against BASE; false atoms against the event-free twin; close pairs; diagnosis (dev)."""
import glob
import json
import os
import sys

import numpy as np

import c10_lib as C
import c14_lib as L
import xc

OUT = os.path.join(L.HERE, "out")
RULES = {"LOC+T": 0.5, "LOC": None}
CAPS = (0.1, 0.25)
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
KE = 5                                             # max events per window
NINF = -np.inf


def load(mode):
    rows = []
    for f in sorted(glob.glob(os.path.join(OUT, f"c14_{mode}_parts", "*.json"))):
        p = json.load(open(f))
        for r in p["rows"]:
            r.update(regime=p["regime"], cell=p["cell"], net=p["net"])
            rows.append(r)
    return rows


def close_flags(r, dkm=10.0, dts=2.0):
    """events with another event of the window within dkm (epicentre) and dts (origin).  Primary definition of the brief
    (10 km, 2 s); TIGHT (7 km, 1 s) is a secondary descriptor added after the freeze, before any held-out result, because
    90% of the DENSE events satisfy the primary one."""
    k = r["k"]
    out = np.zeros(KE, bool)
    for a in range(k):
        for b in range(k):
            if a != b and np.linalg.norm(np.array(r["truth"][a][:2]) - np.array(r["truth"][b][:2])) <= dkm \
                    and abs(r["taus"][a] - r["taus"][b]) <= dts:
                out[a] = True
    return out


class Tab:
    """per design: matrices over windows.  F[m][rule] (nw, kmax) false record scores (-inf padded); H[m][rule] (nw, KE)
    best record score of the atoms that hit event e (-inf: not hit; nan: no such event)."""

    def __init__(self, rows, methods):
        self.rows, self.methods = rows, methods
        nw = len(rows)
        self.k = np.array([r["k"] for r in rows])
        self.close = np.array([close_flags(r) for r in rows])
        self.tight = np.array([close_flags(r, 7.0, 1.0) for r in rows])
        self.exist = np.arange(KE)[None, :] < self.k[:, None]
        kmax = max(len(r["atoms"]["BASE"]) for r in rows)
        self.F = {m: {u: np.full((nw, kmax), NINF) for u in RULES} for m in methods}
        self.H = {m: {u: np.full((nw, KE), NINF) for u in RULES} for m in methods}
        self.dmin = {m: np.full((nw, kmax), np.inf) for m in methods}
        self.S = {m: np.full((nw, kmax), NINF) for m in methods}
        self.S0 = {m: np.full((nw, kmax), NINF) for m in methods}
        self.tw = {m: np.zeros(nw, bool) for m in methods}
        self.dmin0 = {m: np.full((nw, kmax), np.inf) for m in methods}
        for i, r in enumerate(rows):
            truth = [np.array(t) for t in r["truth"]]
            for m in methods:
                at = L.unpack(r["atoms"][m])
                if not at:
                    continue
                s = np.minimum.accumulate([a["s"] for a in at])
                self.S[m][i, :len(s)] = s
                if truth:
                    self.dmin[m][i, :len(at)] = [min(np.linalg.norm(a["x"][:2] - t[:2]) for t in truth) for a in at]
                for u, tol in RULES.items():
                    lab = C.label(at, truth, r["taus"], tol)
                    for p, h in enumerate(lab):
                        if h < 0:
                            self.F[m][u][i, p] = s[p]
                        else:
                            self.H[m][u][i, h] = max(self.H[m][u][i, h], s[p])
                if m in r.get("twin", {}):
                    self.tw[m][i] = True
                    a0 = L.unpack(r["twin"][m])
                    if a0:
                        s0 = np.minimum.accumulate([a["s"] for a in a0])
                        self.S0[m][i, :len(s0)] = s0
                        self.dmin0[m][i, :len(a0)] = [min(np.linalg.norm(a["x"][:2] - t[:2]) for t in truth)
                                                      for a in a0]
        self.strat = np.array([f"{r['regime']}|{r['cell']}" for r in rows])
        self.net = np.array([f"{r['regime']}|{r['cell']}|{r['net']}" for r in rows])

    def groups(self, sel):
        """strata -> list of arrays of window indices (one per network)."""
        g = {}
        for i in np.where(sel)[0]:
            g.setdefault(self.strat[i], {}).setdefault(self.net[i], []).append(i)
        return {s: [np.array(v) for v in d.values()] for s, d in g.items()}


def froc(F, H, exist, nw, caps=CAPS):
    """exact FROC point: largest number of events hit by one threshold with at most cap * nw false atoms.
    Returns recalls and thresholds (a hit needs a score strictly above the (m + 1)-th largest false score)."""
    f = np.sort(F[np.isfinite(F)])[::-1]
    h = H[exist]
    nte = max(1, int(exist.sum()))
    rec, thr = [], []
    for c in caps:
        m = int(np.floor(c * nw + 1e-9))
        t = f[m] if len(f) > m else NINF
        rec.append(float((h > t).sum()) / nte)
        thr.append(t)
    return rec, thr


def metrics(T, w, methods):
    """all scalar endpoints of a multiset of windows w (index array)."""
    out = {}
    nw = len(w)
    ex = T.exist[w]
    cl = T.close[w] & ex
    for m in methods:
        for u in RULES:
            rec, thr = froc(T.F[m][u][w], T.H[m][u][w], ex, nw)
            for c, r_, t in zip(CAPS, rec, thr):
                out[(m, u, c)] = r_
                h = T.H[m][u][w]
                out[(m, u, c, "close")] = float((h[cl] > t).sum()) / max(1, int(cl.sum()))
                out[(m, u, c, "far")] = float((h[ex & ~cl] > t).sum()) / max(1, int((ex & ~cl).sum()))
                if c == 0.25:
                    out[(m, u, "thr25")] = t
            h = T.H[m][u][w]
            out[(m, u, "any")] = float(np.isfinite(h[ex]).sum()) / max(1, int(ex.sum()))
            if u == "LOC+T":
                Fl, Fo = T.F[m]["LOC+T"][w], T.F[m]["LOC"][w]
                a1 = np.isfinite(Fl) & ~np.isfinite(Fo)
                out[(m, "a1", 0.25)] = float((a1 & (Fl > out[(m, u, "thr25")])).sum()) / nw
                out[(m, "a1", "all")] = float(a1.sum()) / nw
                out[(m, "fal", 0.25)] = float((Fl > out[(m, u, "thr25")]).sum()) / nw
            tg = T.tight[w] & ex
            out[(m, u, 0.25, "tight")] = float((h[tg] > out[(m, u, "thr25")]).sum()) / max(1, int(tg.sum()))
            out[(m, u, "any", "tight")] = float(np.isfinite(h[tg]).sum()) / max(1, int(tg.sum()))
            out[(m, u, "any", "close")] = float(np.isfinite(h[cl]).sum()) / max(1, int(cl.sum()))
    return out


def boot(T, sel, methods, B, seed):
    rng = np.random.default_rng(seed)
    G = T.groups(sel)
    allw = np.concatenate([g for gs in G.values() for g in gs])
    pt = metrics(T, allw, methods)
    reps = {k: np.zeros(B) for k in pt}
    for b in range(B):
        ws = []
        for gs in G.values():
            for gi in rng.integers(0, len(gs), len(gs)):
                g = gs[gi]
                ws.append(g[rng.integers(0, len(g), len(g))])
        mb = metrics(T, np.concatenate(ws), methods)
        for k in pt:
            reps[k][b] = mb[k]
    return pt, reps


def ci(x):
    return np.percentile(x, [2.5, 97.5])


def fmt(p, c, sign=True):
    return (f"{p:+.3f} [{c[0]:+.3f},{c[1]:+.3f}]" if sign else f"{p:.3f} [{c[0]:.3f},{c[1]:.3f}]")


# ============================================================================================ twin ratio
def twin_table(T, methods, B, seed, fh, name):
    """false atoms (LOC+T) per event window above the q80 | q95 quantile of the method's OWN twin record scores of the
    stratum (regime x cell), against the twin atoms above the same threshold; all, and farther than 10 km from every
    event."""
    rng = np.random.default_rng(seed)
    sel = (T.k > 0) & T.tw["BASE"]
    for m in methods:
        if T.tw[m][T.k > 0].any():
            assert np.array_equal(T.tw[m] & (T.k > 0), sel), m
    G = T.groups(sel)
    fh.write(f"\n  {name}: {int(sel.sum())} event windows with twins.  W = false atoms per window (LOC+T), W0 = twin atoms per "
             "window, at q80 | q95 of the method's own twin scores\n")
    fh.write(f"  {'method':10s} {'W q80':>18s} {'W0 q80':>7s} {'ratio q80':>18s} | {'W q95':>18s} {'W0 q95':>7s} "
             f"{'ratio q95':>20s} | {'far ratio q95':>20s}\n")
    res = {}
    for m in methods:
        if not T.tw[m][sel].any():
            continue
        cnt = np.zeros((len(T.k), 2, 4))                       # W, W0, Wfar, W0far
        for s_ in G:
            idx = np.concatenate(G[s_])
            s0 = T.S0[m][idx]
            tq = np.quantile(s0[np.isfinite(s0)], [0.80, 0.95])
            for q, t in enumerate(tq):
                f = T.F[m]["LOC+T"][idx] >= t
                cnt[idx, q, 0] = f.sum(1)
                cnt[idx, q, 2] = (f & (T.dmin[m][idx] > 10.0)).sum(1)
                cnt[idx, q, 1] = (s0 >= t).sum(1)
                cnt[idx, q, 3] = ((s0 >= t) & (T.dmin0[m][idx] > 10.0)).sum(1)
        # twin far: distance of twin atoms to events is not stored per method here; use the far fraction of W0 = all
        gs = [g for v in G.values() for g in v]
        strat_of = [s_ for s_, v in G.items() for _ in v]
        allw = np.concatenate(gs)
        pt = cnt[allw].mean(0)
        reps = np.zeros((B, 2, 4))
        for b in range(B):
            ws = []
            for s_ in G:
                gl = G[s_]
                for gi in rng.integers(0, len(gl), len(gl)):
                    g = gl[gi]
                    ws.append(g[rng.integers(0, len(g), len(g))])
            reps[b] = cnt[np.concatenate(ws)].mean(0)
        del strat_of
        line = f"  {m:10s}"
        res[m] = {}
        for q in range(2):
            rr = reps[:, q, 0] / np.maximum(reps[:, q, 1], 1e-9)
            r0 = pt[q, 0] / max(pt[q, 1], 1e-9)
            line += f" {fmt(pt[q, 0], ci(reps[:, q, 0]), False):>18s} {pt[q, 1]:7.3f} {r0:6.2f} [{ci(rr)[0]:.2f},{ci(rr)[1]:.2f}] |"
            res[m][q] = (pt[q, 0], pt[q, 1], r0, reps[:, q, 0], reps[:, q, 1])
        rr = reps[:, 1, 2] / np.maximum(reps[:, 1, 3], 1e-9)
        line += f" {pt[1, 2] / max(pt[1, 3], 1e-9):6.2f} [{ci(rr)[0]:.2f},{ci(rr)[1]:.2f}]"
        fh.write(line + "\n")
    if "BASE" in res:
        fh.write("  paired change of W (false atoms per event window) against BASE, q80 | q95; and of the ratio W/W0\n")
        for m in res:
            if m == "BASE":
                continue
            line = f"  {m:10s}"
            for q in range(2):
                d = res[m][q][3] - res["BASE"][q][3]
                dr = res[m][q][3] / np.maximum(res[m][q][4], 1e-9) - res["BASE"][q][3] / np.maximum(res["BASE"][q][4], 1e-9)
                line += (f" dW {fmt(res[m][q][0] - res['BASE'][q][0], ci(d))}  dratio "
                         f"{res[m][q][2] - res['BASE'][q][2]:+.2f} [{ci(dr)[0]:+.2f},{ci(dr)[1]:+.2f}] |")
            fh.write(line + "\n")


# ============================================================================================ diagnosis (dev)
CATS = ("a1 right place, wrong origin time", "a2 copy of a detected event elsewhere (same arrivals)",
        "b  cross-phase ghost (P read on S or S on P)", "c  deflation residue of claimed arrivals", "d  coda",
        "e  only detection of an otherwise missed event (mislocated)", "f  clutter + noise (not event-induced)")


def classify(r, p, at, lab_lt, lab_loc):
    dg = r["dg"][p]
    if lab_loc[p] >= 0:
        return 0
    ls, lc = np.array(dg["live_same"]), np.array(dg["live_cross"])
    b = [ls.sum(), lc.sum(), float(np.sum(dg["resid"])) + dg["defl"], float(np.sum(dg["coda"])), dg["cn"]]
    j = int(np.argmax(b))
    if j == 0:
        return 1 if int(ls.argmax()) in set(h for h in lab_lt if h >= 0) else 5
    return {1: 2, 2: 3, 3: 4, 4: 6}[j]


def diagnosis(rows, fh):
    fh.write("\n" + "=" * 118 + "\nDIAGNOSIS (development networks, BASE detector; truth used only here)\n")
    for des in ("STD", "DENSE"):
        rs = [r for r in rows if r["design"] == des]
        T = Tab(rs, ["BASE"])
        w = np.arange(len(rs))
        _rec, thr = froc(T.F["BASE"]["LOC+T"][w], T.H["BASE"]["LOC+T"][w], T.exist[w], len(w))
        s0 = T.S0["BASE"][T.k > 0]
        t95 = float(np.quantile(s0[np.isfinite(s0)], 0.95))
        nev = int((T.k > 0).sum())
        for tname, t in (("operating threshold at <= 0.25 false per window", thr[1]),
                         ("operating threshold at <= 0.1 false per window", thr[0]), ("twin q95", t95),
                         ("all atoms (no threshold)", NINF)):
            cnt = np.zeros(len(CATS))
            shares = np.zeros((len(CATS), 6))
            lit = np.zeros((len(CATS), 5))
            n_tot = 0
            a1d = []
            for i, r in enumerate(rs):
                if not r["k"]:
                    continue
                at = L.unpack(r["atoms"]["BASE"])
                truth = [np.array(x) for x in r["truth"]]
                lt, lo = C.label(at, truth, r["taus"], 0.5), C.label(at, truth, r["taus"], None)
                for p in range(len(at)):
                    if lt[p] >= 0 or not (T.S["BASE"][i, p] >= t if np.isfinite(t) else True):
                        continue
                    c = classify(r, p, at, lt, lo)
                    if c == 0:
                        e0 = lo[p]
                        a1d.append([abs(at[p]["tau"] - r["taus"][e0]), at[p]["tau"] - r["taus"][e0], e0 in lt,
                                    float(np.linalg.norm(at[p]["x"][:2] - truth[e0][:2])), p])
                    dg = r["dg"][p]
                    s = dg["s0"]
                    cnt[c] += 1
                    n_tot += 1
                    shares[c] += np.array([np.sum(dg["live_same"]), np.sum(dg["live_cross"]),
                                           np.sum(dg["resid"]) + dg["defl"], np.sum(dg["coda"]), dg["cn"], s]) / s
                    tot_e = np.array(dg["live_same"]) + np.array(dg["live_cross"]) + np.array(dg["resid"]) \
                        + np.array(dg["coda"])
                    edom = int(np.argmax(tot_e))
                    lit[c] += [dg["f_arr"] / s, dg["f_arr"] / s >= 0.5, dg["f_acc"] / s, edom in dg["acc"],
                               dg["f_cross"] / s]
            fh.write(f"\n  {des}: {nev} event windows; false atoms (LOC+T) above {tname} (t = {t:.1f}): {n_tot} "
                     f"({n_tot / max(1, nev):.2f} per event window)\n")
            fh.write("    class                                                          n  share  | mean share of the "
                     "statistic: live-same live-cross residue  coda  clut+noise | f_arr  f_arr>=.5  f_acc  dom.acc  f_cross\n")
            for c, nm in enumerate(CATS):
                n = max(1, cnt[c])
                fh.write(f"    {nm:58s} {int(cnt[c]):5d}  {cnt[c] / max(1, n_tot):5.2f}  |   "
                         + "  ".join(f"{x:8.2f}" for x in shares[c, :5] / n) + "   | "
                         + "  ".join(f"{x:6.2f}" for x in lit[c] / n) + "\n")
            n = max(1, n_tot)
            if a1d:
                a1d_ = np.array(a1d)
                fh.write(f"    a1 detail: n {len(a1d_)}; |origin error| median {np.median(a1d_[:, 0]):.2f} s, quartiles "
                         f"{np.quantile(a1d_[:, 0], 0.25):.2f}-{np.quantile(a1d_[:, 0], 0.75):.2f}; sign + (late) "
                         f"{np.mean(a1d_[:, 1] > 0):.2f}; same event also hit with the right time by another atom "
                         f"{np.mean(a1d_[:, 2]):.2f}; epicentral error median {np.median(a1d_[:, 3]):.2f} km; rank 0 "
                         f"in the window {np.mean(a1d_[:, 4] == 0):.2f}\n")
            fh.write(f"    {'ALL':58s} {n_tot:5d}  {1.0:5.2f}  |   "
                     + "  ".join(f"{x:8.2f}" for x in shares.sum(0)[:5] / n) + "   | "
                     + "  ".join(f"{x:6.2f}" for x in lit.sum(0) / n) + "\n")


# ============================================================================================ main
def lines_of(T):
    out = [("ALL", np.ones(len(T.k), bool))]
    for c in CELLS:
        out.append((c, np.array([r["cell"] == c for r in T.rows])))
    for g in ("strong", "thin"):
        out.append((g, np.array([r["regime"] == g for r in T.rows])))
    return out


def main(mode, B):
    rows = load(mode)
    methods = list(rows[0]["atoms"].keys())
    fh = open(os.path.join(OUT, f"c14_analysis_{mode}.txt"), "w")
    nets = sorted({r["net"] for r in rows})
    fh.write(f"C14 analysis, mode {mode}: {len(rows)} windows, networks {nets}, methods {methods}, B = {B}\n")
    # check of the FROC against xc.recall_caps
    for des in ("STD", "DENSE"):
        rs = [r for r in rows if r["design"] == des]
        T = Tab(rs, methods)
        w = np.arange(len(rs))
        for m in methods[:3]:
            for u, tol in RULES.items():
                recs = []
                for r in rs:
                    at = L.unpack(r["atoms"][m])
                    recs.append(C.record(at, C.label(at, [np.array(t) for t in r["truth"]], r["taus"], tol), "s"))
                a, _t = xc.recall_caps(recs, int(T.k.sum()))
                b, _t = froc(T.F[m][u][w], T.H[m][u][w], T.exist[w], len(w))
                assert np.allclose(a, b), (des, m, u, a, b)
        fh.write(f"\n{'=' * 118}\n{des}: {len(rs)} windows, {int(T.k.sum())} events, {int((T.close & T.exist).sum())} "
                 "events in a close pair (another event within 10 km and 2 s)\n")
        for lname, sel in lines_of(T):
            if mode != "eval" and lname != "ALL":
                continue
            pt, reps = boot(T, sel, methods, B, 20261005)
            for u in RULES:
                fh.write(f"\n  {des} {lname} rule {u}: recall at <= 0.1 | <= 0.25 false events per window; paired "
                         "difference against BASE\n")
                for m in methods:
                    line = f"    {m:10s}"
                    for c in CAPS:
                        line += f"  {fmt(pt[(m, u, c)], ci(reps[(m, u, c)]), False)}"
                    if m != "BASE":
                        for c in CAPS:
                            d = reps[(m, u, c)] - reps[("BASE", u, c)]
                            line += f"  d {fmt(pt[(m, u, c)] - pt[('BASE', u, c)], ci(d))} sd {d.std():.3f}"
                    fh.write(line + "\n")
            if lname == "ALL":
                fh.write(f"\n  {des} ALL: 'right place, wrong origin time' atoms (within 5 km, |origin error| > 0.5 s) per "
                         "window: above the LOC+T operating threshold at <= 0.25 (and all false atoms there) | at any score\n")
                for m in methods:
                    k1, k2, k3 = (m, "a1", 0.25), (m, "fal", 0.25), (m, "a1", "all")
                    line = (f"    {m:10s}  {fmt(pt[k1], ci(reps[k1]), False)} of {pt[k2]:.3f}   |  "
                            f"{fmt(pt[k3], ci(reps[k3]), False)}")
                    if m != "BASE":
                        line += f"   d(any) {fmt(pt[k3] - pt[('BASE',) + k3[1:]], ci(reps[k3] - reps[('BASE',) + k3[1:]]))}"
                    fh.write(line + "\n")
            if des == "DENSE" and lname == "ALL":
                for u in RULES:
                    fh.write(f"\n  DENSE ALL rule {u}: recall of events in a CLOSE PAIR at the operating point (<= 0.25 "
                             "| <= 0.1), of the other events (<= 0.25), at any score (close | all), TIGHT pairs (<= 0.25 | "
                             "any score)\n")
                    for m in methods:
                        line = f"    {m:10s}"
                        for k in ((m, u, 0.25, "close"), (m, u, 0.1, "close"), (m, u, 0.25, "far"),
                                  (m, u, "any", "close"), (m, u, "any"), (m, u, 0.25, "tight"), (m, u, "any", "tight")):
                            line += f"  {fmt(pt[k], ci(reps[k]), False)}"
                            if m != "BASE":
                                kb = ("BASE",) + k[1:]
                                line += f" d {fmt(pt[k] - pt[kb], ci(reps[k] - reps[kb]))}"
                        fh.write(line + "\n")
        if mode == "eval":
            fh.write("\n" + "-" * 118 + f"\nFALSE ATOMS AGAINST THE EVENT-FREE TWIN ({des})\n")
            twin_table(T, methods, B, 79, fh, f"{des} ALL")
        elif mode == "dev":
            twin_table(T, ["BASE"], B, 79, fh, f"{des} ALL")
    if mode == "dev":
        diagnosis(rows, fh)
    fh.close()
    print(open(os.path.join(OUT, f"c14_analysis_{mode}.txt")).read())


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 500)
