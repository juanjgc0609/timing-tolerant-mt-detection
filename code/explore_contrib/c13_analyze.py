"""C13 analysis: exchangeability of decoys and targets, realised false discovery proportion / false events per window /
recall of the decoy-calibrated thresholds, loss against the oracle operating point at the same realised false-event
rate.  Hierarchical bootstrap (networks, then windows; thresholds, batches and the oracle point recomputed per
replicate).
    python c13_analyze.py [B] [NPERM]   ->  out/c13_analysis.txt, out/c13_summary.json
"""
import glob
import json
import os
import sys

import numpy as np

import c10_lib as C
import c13_lib as L

REGIMES = ("thin", "strong")
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
RULES = ("LOC+T", "LOC")
TOL = {"LOC+T": 0.5, "LOC": None}
METHODS = ("tdc", "bh", "comp", "fpw")
CS = np.array((0.05, 0.1, 0.25, 0.5, 1.0))   # fpw: nominal false events per window (same slots as QS)
BS = (1, 6)                      # networks per batch: 20 and 120 windows
QS = np.array(L.QS)
D3C = os.path.join(L.HERE, "out", "c13_d3c_parts")
HAS_D3C = len(glob.glob(os.path.join(D3C, "*.json"))) == len(glob.glob(os.path.join(L.HERE, "out", "c13_parts", "*.json")))
DEC = L.DECOYS + (("D3c",) if HAS_D3C else ())
OUT = os.path.join(L.HERE, "out")


# ============================================================================================ loading
class Unit:
    """windows of one (regime, cell, network, design)."""

    def __init__(self, rows):
        W = len(rows)
        K = len(rows[0]["atoms"]["T"])
        self.W, self.K = W, K
        self.k = np.array([r["k"] for r in rows])
        self.sT = np.zeros((W, K))
        self.sD = {d: np.zeros((W, K)) for d in DEC}
        self.false = {u: np.zeros((W, K), bool) for u in RULES}
        self.first = {u: np.zeros((W, K), bool) for u in RULES}
        self.dmin = np.full((W, K), np.inf)          # epicentral distance of the atom to the nearest event
        self.dtn = np.full((W, K), np.inf)           # |origin error| against that event
        self.dhit = {d: np.zeros((W, K), bool) for d in DEC}     # decoy atom within 5 km of an event (LOC)
        self.sec = np.array([[r["sec_t"], r["sec_d"]] for r in rows])
        for i, r in enumerate(rows):
            at = L.unpack(r["atoms"]["T"])
            assert len(at) == K
            truth = [np.array(t) for t in r["truth"]]
            self.sT[i] = L.runmin([a["s"] for a in at])
            for u in RULES:
                lab = np.array(C.label(at, truth, r["taus"], TOL[u]))
                self.false[u][i] = lab < 0
                seen = set()
                for p, h in enumerate(lab):
                    if h >= 0 and h not in seen:
                        seen.add(h)
                        self.first[u][i, p] = True
            if truth:
                for p, a in enumerate(at):
                    dd = np.array([np.linalg.norm(a["x"][:2] - t[:2]) for t in truth])
                    j = int(dd.argmin())
                    self.dmin[i, p], self.dtn[i, p] = dd[j], abs(a["tau"] - r["taus"][j])
            for d in DEC:
                da = L.unpack(r["atoms"][d])
                self.sD[d][i] = L.runmin([a["s"] for a in da])
                self.dhit[d][i] = np.array(C.label(da, truth, r["taus"], None)) >= 0


def load():
    data, nets = {}, {}
    for f in sorted(glob.glob(os.path.join(OUT, "c13_parts", "*.json"))):
        p = json.load(open(f))
        if HAS_D3C:
            x = json.load(open(os.path.join(D3C, os.path.basename(f))))
            for r, rx in zip(p["rows"], x["rows"]):
                assert (r["design"], r["w"]) == (rx["design"], rx["w"])
                r["atoms"]["D3c"] = rx["D3c"]
                r["sec_d3c"] = rx["sec"]
        by = {}
        for r in p["rows"]:
            by.setdefault(r["design"], []).append(r)
        for des, rows in by.items():
            data.setdefault((p["regime"], p["cell"], des), {})[p["net"]] = Unit(rows)
        nets.setdefault((p["regime"], p["cell"]), []).append(p["net"])
    return data, nets


# ============================================================================================ calibration metrics
def froc(sc, fal, fst):
    """sorted false / hit counts of a pooled set (ties resolved as in the exact FROC)."""
    o = np.argsort(-sc, kind="stable")
    sc, Fc, Hc = sc[o], np.cumsum(fal[o]), np.cumsum(fst[o])
    last = np.r_[np.where(np.diff(sc) != 0)[0], len(sc) - 1]
    return Fc[last], Hc[last]


def oracle_hits(Fc, Hc, F):
    """largest number of events hit by one threshold with at most F false atoms."""
    i = np.searchsorted(Fc, F, side="right") - 1
    return np.where(i >= 0, Hc[np.maximum(i, 0)], 0)


def batch_eval(units, widx):
    """units: list of Unit; widx: list of window index arrays.  Returns acc arrays [decoy, method, q, rule, metric]
    with metric = (fdp, n_false, n_acc, n_hits, any)."""
    sT = np.concatenate([u.sT[w].ravel() for u, w in zip(units, widx)])
    fal = {r: np.concatenate([u.false[r][w].ravel() for u, w in zip(units, widx)]) for r in RULES}
    fst = {r: np.concatenate([u.first[r][w].ravel() for u, w in zip(units, widx)]) for r in RULES}
    out = np.zeros((len(DEC), len(METHODS), len(QS), len(RULES), 5))
    nwin = sum(len(w) for w in widx)
    for di, d in enumerate(DEC):
        sD = np.concatenate([u.sD[d][w].ravel() for u, w in zip(units, widx)])
        win = sT > sD
        sDd = np.sort(sD)[::-1]
        for mi, m in enumerate(METHODS):
            if m == "fpw":
                # decoy estimate of the false events per window: at most floor(c * windows) decoys above the threshold
                na_ = np.floor(CS * nwin + 1e-9).astype(int)
                thr = [sDd[k] if k < len(sDd) else -np.inf for k in na_]
            else:
                thr = {"tdc": L.thr_tdc, "bh": L.thr_bh, "comp": L.thr_comp}[m](sT, sD, QS)
            for qi, t in enumerate(thr):
                acc = sT > t if m == "fpw" else sT >= t
                if m == "comp":
                    acc &= win
                na = int(acc.sum())
                for ri, r in enumerate(RULES):
                    nf = int((acc & fal[r]).sum())
                    out[di, mi, qi, ri] = (nf / max(1, na), nf, na, int((acc & fst[r]).sum()), na > 0)
    return out


def stratum_eval(units, widx, rng, nrep6):
    """one stratum (regime, cell, design): batches of 1 network and of 6 networks (random partitions, nrep6 times).
    Returns res[bs] = dict(sum over batches of the metrics, number of batches) and the pooled arrays for the oracle."""
    res = {}
    outs = [batch_eval([u], [w]) for u, w in zip(units, widx)]
    acc = sum(outs)
    res[1] = (acc, len(units))
    res["b1"] = np.array([o_[..., 0] for o_ in outs])            # FDP of every 20-window batch
    n = len(units)
    acc6, nb = 0, 0
    for _ in range(nrep6):
        perm = rng.permutation(n)
        for g in range(0, n - n % 6 if n >= 6 else n, 6):
            ii = perm[g:g + 6]
            acc6 = acc6 + batch_eval([units[i] for i in ii], [widx[i] for i in ii])
            nb += 1
    res[6] = (acc6 / nrep6, nb / nrep6)
    pooled = {"sT": np.concatenate([u.sT[w].ravel() for u, w in zip(units, widx)]),
              "nw": int(sum(len(w) for w in widx)), "ne": int(sum(u.k[w].sum() for u, w in zip(units, widx)))}
    for r in RULES:
        pooled["fal" + r] = np.concatenate([u.false[r][w].ravel() for u, w in zip(units, widx)])
        pooled["fst" + r] = np.concatenate([u.first[r][w].ravel() for u, w in zip(units, widx)])
    return res, pooled


def summarise(parts):
    """parts: list of (res, pooled) of the strata pooled into one line.  Returns arrays [bs, decoy, method, q, rule] of
    FDR (mean FDP over batches), pooled FDP, false events per window, recall, oracle recall at the same number of false
    atoms, share of batches accepting something; and the benchmark recalls / FDPs at <= 0.1 and <= 0.25."""
    nw = sum(p["nw"] for _r, p in parts)
    ne = max(1, sum(p["ne"] for _r, p in parts))
    sT = np.concatenate([p["sT"] for _r, p in parts])
    shape = (len(BS), len(DEC), len(METHODS), len(QS), len(RULES))
    o = {k: np.zeros(shape) for k in ("fdr", "pfdp", "few", "rec", "orec", "any")}
    bench = np.zeros((len(RULES), 2, 3))
    for ri, r in enumerate(RULES):
        Fc, Hc = froc(sT, np.concatenate([p["fal" + r] for _r, p in parts]),
                      np.concatenate([p["fst" + r] for _r, p in parts]))
        for ci, cap in enumerate((0.1, 0.25)):
            i = np.searchsorted(Fc, cap * nw + 1e-9, side="right") - 1
            if i >= 0:
                tp = Hc[i]
                bench[ri, ci] = (tp / ne, Fc[i] / max(1, Fc[i] + tp), Fc[i] / nw)
        for bi, bs in enumerate(BS):
            acc = sum(res[bs][0] for res, _p in parts)
            nb = sum(res[bs][1] for res, _p in parts)
            a = acc[..., ri, :]
            o["fdr"][bi, ..., ri] = a[..., 0] / nb
            o["pfdp"][bi, ..., ri] = a[..., 1] / np.maximum(1, a[..., 2])
            o["few"][bi, ..., ri] = a[..., 1] / nw
            o["rec"][bi, ..., ri] = a[..., 3] / ne
            o["orec"][bi, ..., ri] = oracle_hits(Fc, Hc, np.floor(a[..., 1] + 1e-9)) / ne
            o["any"][bi, ..., ri] = a[..., 4] / nb
    o["loss"] = o["orec"] - o["rec"]
    # standard error of the mean FDP over the independent 20-window batches (one per network); batch 20 only
    b1 = np.concatenate([res["b1"] for res, _p in parts])
    o["se20"] = np.zeros(shape)
    o["se20"][0] = b1.std(axis=0, ddof=1) / np.sqrt(len(b1)) if len(b1) > 1 else 0.0
    # fpw only: loss against the oracle point at the NOMINAL false-event rate (the benchmark operating point)
    o["lossn"] = np.zeros(shape)
    fi = METHODS.index("fpw")
    for ri, r in enumerate(RULES):
        Fc, Hc = froc(sT, np.concatenate([p["fal" + r] for _r, p in parts]),
                      np.concatenate([p["fst" + r] for _r, p in parts]))
        on = oracle_hits(Fc, Hc, np.floor(CS * nw + 1e-9)) / ne
        o["lossn"][:, :, fi, :, ri] = on[None, None, :] - o["rec"][:, :, fi, :, ri]
    return o, bench, nw, ne


def lines_of(data, nets):
    """reporting lines: name -> list of strata (regime, cell, design)."""
    ln = {}
    for des in ("STD", "DENSE"):
        for reg in REGIMES:
            for cell in CELLS:
                ln[f"{des} {reg} {cell}"] = [(reg, cell, des)]
            ln[f"{des} {reg} ALL"] = [(reg, c, des) for c in CELLS]
        for cell in CELLS:
            ln[f"{des} both {cell}"] = [(reg, cell, des) for reg in REGIMES]
        ln[f"{des} both ALL"] = [(reg, c, des) for reg in REGIMES for c in CELLS]
    return {k: [s for s in v if s in data] for k, v in ln.items()}


def run_calibration(data, nets, B, fh, seed=20261005):
    rng = np.random.default_rng(seed)
    ln = lines_of(data, nets)
    strata = sorted({s for v in ln.values() for s in v})

    def one(resample):
        ev = {}
        for s in strata:
            js = sorted(data[s])
            if resample:
                js = [js[i] for i in rng.integers(0, len(js), len(js))]
            units = [data[s][j] for j in js]
            widx = [rng.integers(0, u.W, u.W) if resample else np.arange(u.W) for u in units]
            ev[s] = stratum_eval(units, widx, rng, 2 if resample else 20)
        return {name: summarise([ev[s] for s in ss]) for name, ss in ln.items() if ss}

    pt = one(False)
    reps = [one(True) for _ in range(B)]
    summ = {}
    for name in pt:
        o, bench, nw, ne = pt[name]
        rb = np.array([r[name][1] for r in reps])
        fh.write(f"\n=== {name}: {nw} windows, {ne} events\n")
        for ri, r in enumerate(RULES):
            fh.write(f"  benchmark oracle {r}: recall@0.1 {bench[ri, 0, 0]:.3f} [{np.percentile(rb[:, ri, 0, 0], 2.5):.3f},"
                     f" {np.percentile(rb[:, ri, 0, 0], 97.5):.3f}] (FDP there {bench[ri, 0, 1]:.2f})   recall@0.25 "
                     f"{bench[ri, 1, 0]:.3f} [{np.percentile(rb[:, ri, 1, 0], 2.5):.3f}, "
                     f"{np.percentile(rb[:, ri, 1, 0], 97.5):.3f}] (FDP there {bench[ri, 1, 1]:.2f})\n")
        summ[name] = {"nw": nw, "ne": ne, "bench": bench.tolist(), "cal": {}}
        R_ = {k: np.array([r[name][0][k] for r in reps]) for k in o}
        for bi, bs in enumerate(BS):
            for mi, m in enumerate(METHODS):
                for ri, r in enumerate(RULES):
                    fh.write(f"  -- batch {bs * 20} windows, {m}, {r}: per decoy and q: FDR [95%] | false/win | recall "
                             f"[95%] | oracle recall at the same false/win | loss [95%] | share of batches accepting\n")
                    for di, d in enumerate(DEC):
                        for qi, q in enumerate(QS):
                            ix = (bi, di, mi, qi, ri)
                            c = {k: np.percentile(R_[k][(slice(None),) + ix], [2.5, 97.5]) for k in o}
                            q = CS[qi] if m == "fpw" else q
                            fh.write(f"     {d:4s} {'c' if m == 'fpw' else 'q'}={q:.2f}  FDR {o['fdr'][ix]:.3f} [{c['fdr'][0]:.3f}, {c['fdr'][1]:.3f}]"
                                     f"  pooledFDP {o['pfdp'][ix]:.3f}  f/w {o['few'][ix]:.3f} [{c['few'][0]:.3f}, "
                                     f"{c['few'][1]:.3f}]  rec {o['rec'][ix]:.3f} [{c['rec'][0]:.3f}, {c['rec'][1]:.3f}]"
                                     f"  orec {o['orec'][ix]:.3f}  loss {o['loss'][ix]:+.3f} [{c['loss'][0]:+.3f}, "
                                     f"{c['loss'][1]:+.3f}]  any {o['any'][ix]:.2f}\n")
                            summ[name]["cal"]["|".join(map(str, (bs * 20, d, m, float(q), r)))] = {
                                k: [float(o[k][ix]), float(c[k][0]), float(c[k][1])] for k in o}
    return summ


# ============================================================================================ exchangeability
def ks(a, b):
    a, b = np.sort(a), np.sort(b)
    x = np.concatenate([a, b])
    return float(np.abs(np.searchsorted(a, x, side="right") / len(a) - np.searchsorted(b, x, side="right") / len(b)).max())


def swap_test(T, D, col, nperm, rng):
    """T, D (W, K) record scores of the same windows.  KS of T[:, col] against D[:, col]; p-value by swapping the two
    lists of a window (valid under exchangeability of the pair, whatever the dependence between windows of a network is
    NOT assumed away: swaps are per window)."""
    obs = ks(T[:, col].ravel(), D[:, col].ravel())
    cnt = 0
    for _ in range(nperm):
        sw = rng.random(len(T)) < 0.5
        a = np.where(sw[:, None], D, T)[:, col].ravel()
        b = np.where(sw[:, None], T, D)[:, col].ravel()
        cnt += ks(a, b) >= obs - 1e-12
    return obs, (1 + cnt) / (1 + nperm)


def null_windows(data, reg, cell, kinds):
    """units restricted to the windows without events: list of (net, sT, {d: sD})."""
    out = []
    for des in kinds:
        for j, u in sorted(data.get((reg, cell, des), {}).items()):
            w = np.where(u.k == 0)[0]
            if len(w):
                out.append((j, u.sT[w][:, :5], {d: u.sD[d][w][:, :5] for d in DEC}))
    return out


def boot_ratio(groups, fn, B, rng):
    """hierarchical bootstrap (networks, then windows) of a statistic fn(list of (T, D) per network)."""
    nets = sorted({g[0] for g in groups})
    by = {j: [g for g in groups if g[0] == j] for j in nets}
    vals = []
    for _ in range(B):
        sel = []
        for j in [nets[i] for i in rng.integers(0, len(nets), len(nets))]:
            for g in by[j]:
                w = rng.integers(0, len(g[1]), len(g[1]))
                sel.append((g[1][w], g[2][w]))
        vals.append(fn(sel))
    return np.nanpercentile(vals, [2.5, 97.5])


def tail_ratio(sel, t):
    nT = sum(int((a >= t).sum()) for a, _b in sel)
    nD = sum(int((b >= t).sum()) for _a, b in sel)
    return nT / nD if nD else np.nan


def run_exchange(data, nets, B, nperm, fh, seed=77):
    rng = np.random.default_rng(seed)
    summ = {}
    fh.write("\n" + "=" * 110 + "\nEXCHANGEABILITY, WINDOWS WITHOUT EVENTS (first 5 atoms; record scores)\n")
    for label, kinds in (("CLUT + empty STD (clutter only)", ("CLUT", "STD")), ("NOISE (no clutter)", ("NOISE",))):
        fh.write(f"\n-- {label}: KS distance target vs decoy [swap-permutation p]: top atom | all atoms; "
                 f"ratio #T>=t / #D>=t at t = pooled 0.80 and 0.95 quantile [95% bootstrap]\n")
        for reg in REGIMES:
            for cell in CELLS:
                nw_ = null_windows(data, reg, cell, kinds)
                if not nw_:
                    continue
                T = np.concatenate([g[1] for g in nw_])
                fh.write(f"  {reg:6s} {cell:15s} {len(T):4d} windows; target top-atom quantiles 50/75/90/95: "
                         + " ".join(f"{v:.1f}" for v in np.quantile(T[:, 0], [0.5, 0.75, 0.9, 0.95])) + "\n")
                for d in DEC:
                    D = np.concatenate([g[2][d] for g in nw_])
                    k1, p1 = swap_test(T, D, [0], nperm, rng)
                    ka, pa = swap_test(T, D, slice(None), nperm, rng)
                    allv = np.concatenate([T.ravel(), D.ravel()])
                    rr = []
                    for qq in (0.80, 0.95):
                        t = float(np.quantile(allv, qq))
                        grp = [(g[0], g[1], g[2][d]) for g in nw_]
                        rr.append((tail_ratio([(g[1], g[2]) for g in grp], t),
                                   boot_ratio(grp, lambda s, t=t: tail_ratio(s, t), B, rng)))
                    qd = np.quantile(D[:, 0], [0.5, 0.75, 0.9, 0.95])
                    fh.write(f"      {d:4s} top q {' '.join(f'{v:.1f}' for v in qd)}  KS top {k1:.3f} [p {p1:.3f}]  KS all "
                             f"{ka:.3f} [p {pa:.3f}]  ratio@.80 {rr[0][0]:.2f} [{rr[0][1][0]:.2f}, {rr[0][1][1]:.2f}]"
                             f"  ratio@.95 {rr[1][0]:.2f} [{rr[1][1][0]:.2f}, {rr[1][1][1]:.2f}]\n")
                    summ["|".join((label[:4], reg, cell, d))] = {
                        "n": int(len(T)), "ks_top": k1, "p_top": p1, "ks_all": ka, "p_all": pa,
                        "r80": [rr[0][0]] + rr[0][1].tolist(), "r95": [rr[1][0]] + rr[1][1].tolist()}
    return summ


def run_event_windows(data, nets, B, fh, seed=78):
    """windows WITH events: false targets against decoys, and against the targets of the clutter-only windows of the same
    networks, at absolute thresholds t80 / t95 = the 0.80 / 0.95 quantile of the target atom scores of the CLUT windows of
    the stratum (1.0 and 0.25 atoms per window there)."""
    rng = np.random.default_rng(seed)
    summ = {}
    fh.write("\n" + "=" * 110 + "\nWINDOWS WITH EVENTS: FALSE TARGETS AGAINST DECOYS (rates per window at t80 | t95)\n"
             "cats of false targets (LOC+T): tmiss = within 5 km of an event, origin error > 0.5 s; near = 5-10 km of "
             "the nearest event; far = > 10 km; post = false atom after the first LOC hit of its window\n")
    for des in ("STD", "DENSE"):
        for reg in REGIMES:
            for cell in CELLS:
                if (reg, cell, des) not in data or (reg, cell, "CLUT") not in data:
                    continue
                cl = data[(reg, cell, "CLUT")]
                tq = np.quantile(np.concatenate([u.sT.ravel() for u in cl.values()]), [0.80, 0.95])
                U = data[(reg, cell, des)]

                def stats(js, widx, clw):
                    o = {}
                    nw = sum(len(w) for w in widx)
                    ncl = sum(len(w) for w in clw)
                    for ti, t in enumerate(tq):
                        c = {k: 0 for k in ("false", "falseLOC", "tmiss", "near", "far", "post", "hit", "clutT")}
                        c.update({d: 0 for d in DEC})
                        c.update({"clut" + d: 0 for d in DEC})
                        c.update({"hit" + d: 0 for d in DEC})
                        for j, w, wc in zip(js, widx, clw):
                            u = U[j]
                            ev = w[u.k[w] > 0]
                            a = u.sT[ev] >= t
                            f = u.false["LOC+T"][ev] & a
                            c["false"] += f.sum()
                            c["falseLOC"] += (u.false["LOC"][ev] & a).sum()
                            c["hit"] += (~u.false["LOC+T"][ev] & a).sum()
                            dm = u.dmin[ev]
                            c["tmiss"] += (f & (dm <= 5.0)).sum()
                            c["near"] += (f & (dm > 5.0) & (dm <= 10.0)).sum()
                            c["far"] += (f & (dm > 10.0)).sum()
                            hl = ~u.false["LOC"][ev]
                            after = (np.cumsum(hl, axis=1) - hl) > 0
                            c["post"] += (f & after).sum()
                            for d in DEC:
                                c[d] += (u.sD[d][ev] >= t).sum()
                                c["hit" + d] += ((u.sD[d][ev] >= t) & u.dhit[d][ev]).sum()
                                c["clut" + d] += (cl[j].sD[d][wc] >= t).sum()
                            c["clutT"] += (cl[j].sT[wc] >= t).sum()
                        nev = sum(int((U[j].k[w] > 0).sum()) for j, w in zip(js, widx))
                        for k, v in c.items():
                            o[(k, ti)] = v / max(1, ncl if k.startswith("clut") else nev)
                    return o

                js = sorted(set(U) & set(cl))
                pt = stats(js, [np.arange(U[j].W) for j in js], [np.arange(cl[j].W) for j in js])
                reps = []
                for _ in range(B):
                    jb = [js[i] for i in rng.integers(0, len(js), len(js))]
                    reps.append(stats(jb, [rng.integers(0, U[j].W, U[j].W) for j in jb],
                                      [rng.integers(0, cl[j].W, cl[j].W) for j in jb]))
                name = f"{des} {reg} {cell}"
                fh.write(f"\n  {name}: t80 {tq[0]:.1f}, t95 {tq[1]:.1f}\n")
                summ[name] = {}
                fs = np.concatenate([U[j].sT[U[j].k > 0][U[j].false["LOC+T"][U[j].k > 0]] for j in js])
                kk = {d: ks(fs, np.concatenate([U[j].sD[d][U[j].k > 0].ravel() for j in js])) for d in DEC}
                fh.write("      KS distance, scores of false targets (LOC+T) against decoy scores: "
                         + "  ".join(f"{d} {v:.3f}" for d, v in kk.items()) + "\n")
                summ[name]["ks_false"] = kk
                keys = ["clutT", "false", "falseLOC", "tmiss", "near", "far", "post", "hit"] + list(DEC) + \
                       ["hit" + d for d in DEC] + ["clut" + d for d in DEC]
                for k in keys:
                    v = [pt[(k, ti)] for ti in (0, 1)]
                    ci = [np.percentile([r[(k, ti)] for r in reps], [2.5, 97.5]) for ti in (0, 1)]
                    fh.write(f"      {k:9s} {v[0]:.3f} [{ci[0][0]:.3f}, {ci[0][1]:.3f}] | {v[1]:.3f} [{ci[1][0]:.3f}, "
                             f"{ci[1][1]:.3f}]\n")
                    summ[name][k] = [[v[ti], float(ci[ti][0]), float(ci[ti][1])] for ti in (0, 1)]
                for num, den in (("false", "D1p"), ("false", "D1o"), ("false", "D2s"), ("false", "D2r"), ("false", "D3"),
                                 ("far", "clutT"), ("far", "D3"), ("falseLOC", "D3"), ("falseLOC", "D1p")) + \
                        ((("false", "D3c"), ("falseLOC", "D3c")) if HAS_D3C else ()):
                    v = [pt[(num, ti)] / max(pt[(den, ti)], 1e-9) for ti in (0, 1)]
                    ci = [np.percentile([r[(num, ti)] / max(r[(den, ti)], 1e-9) for r in reps], [2.5, 97.5])
                          for ti in (0, 1)]
                    fh.write(f"      ratio {num}/{den:6s} {v[0]:.2f} [{ci[0][0]:.2f}, {ci[0][1]:.2f}] | {v[1]:.2f} "
                             f"[{ci[1][0]:.2f}, {ci[1][1]:.2f}]\n")
                    summ[name][f"{num}/{den}"] = [[v[ti], float(ci[ti][0]), float(ci[ti][1])] for ti in (0, 1)]
    return summ


def compact(cal, fh):
    """compact tables of the calibration summary."""
    names = list(cal)

    def g(name, bs, d, m, q, r, k):
        return cal[name]["cal"]["|".join(map(str, (bs, d, m, float(q), r)))][k]

    fh.write("\n" + "=" * 110 + "\nCOMPACT TABLES\n")
    for r in RULES:
        for bs in (20, 120):
            for m in ("tdc", "comp"):
                for q in (0.05, 0.1, 0.2):
                    fh.write(f"\nFDR (mean FDP over batches) [95%], {m}, batch {bs} windows, nominal q = {q}, rule {r}; "
                             f"* = interval entirely above q\n")
                    fh.write(f"  {'line':28s}" + "".join(f"{d:>24s}" for d in DEC) + "\n")
                    for n in names:
                        fh.write(f"  {n:28s}")
                        for d in DEC:
                            v = g(n, bs, d, m, q, r, "fdr")
                            fh.write(f"   {v[0]:.3f} [{v[1]:.3f},{v[2]:.3f}]{'*' if v[1] > q else ' '}")
                        fh.write("\n")
    for r in RULES:
        for bs in (20, 120):
            for q in (0.1, 0.2, 0.5):
                fh.write(f"\nRecall of tdc | realised false events per window | loss against the oracle threshold at the "
                         f"same realised rate [95%]; batch {bs}, q = {q}, rule {r}\n")
                fh.write(f"  {'line':28s}" + "".join(f"{d:>38s}" for d in DEC) + "\n")
                for n in names:
                    fh.write(f"  {n:28s}")
                    for d in DEC:
                        a, b, c = (g(n, bs, d, "tdc", q, r, k) for k in ("rec", "few", "loss"))
                        fh.write(f"   {a[0]:.3f} | {b[0]:.3f} | {c[0]:+.3f} [{c[1]:+.3f},{c[2]:+.3f}]")
                    fh.write("\n")
    for r in RULES:
        for bs in (20, 120):
            for ci_, c_ in ((1, 0.1), (2, 0.25)):
                fh.write(f"\nfpw (decoy estimate of the false events per window), nominal c = {c_}: realised false events "
                         f"per window [95%] | recall | loss against the benchmark oracle point at <= c [95%]; batch {bs}, "
                         f"rule {r}; * = interval of the realised rate entirely above c\n")
                fh.write(f"  {'line':28s}" + "".join(f"{d:>50s}" for d in DEC) + f"   oracle recall at <= c\n")
                for n in names:
                    fh.write(f"  {n:28s}")
                    for d in DEC:
                        a, b, c = (g(n, bs, d, "fpw", c_, r, k) for k in ("few", "rec", "lossn"))
                        fh.write(f"   {a[0]:.3f} [{a[1]:.3f},{a[2]:.3f}]{'*' if a[1] > c_ else ' '} | {b[0]:.3f} | "
                                 f"{c[0]:+.3f} [{c[1]:+.3f},{c[2]:+.3f}]")
                    fh.write(f"   {cal[n]['bench'][RULES.index(r)][ci_ - 1][0]:.3f}\n")


def run_twin(data, B, fh, seed=79):
    """event windows W against their clutter-and-noise twins W0 (c13_twin.py): atoms per window above t80 | t95 (CLUT
    quantiles of the stratum) by distance to the nearest event epicentre (in5 <= 5 km, near 5-10 km, far > 10 km).
    W: false atoms under LOC+T (in5 = origin-time misses), far atoms before / after the first LOC hit of the window."""
    rng = np.random.default_rng(seed)
    tw = {}
    for f in sorted(glob.glob(os.path.join(OUT, "c13_twin_parts", "*.json"))):
        p = json.load(open(f))
        for r in p["rows"]:
            tw.setdefault((p["regime"], p["cell"], r["design"], p["net"]), []).append(r)
    if not tw:
        return {}
    summ = {}
    fh.write("\n" + "=" * 110 + "\nEVENT WINDOWS AGAINST THEIR CLUTTER-AND-NOISE TWINS (atoms per event window at t80 | t95)\n")
    keys = ("W_false", "W_in5", "W_near", "W_far", "W_far_pre", "W_far_post", "W0_all", "W0_in5", "W0_near", "W0_far",
            "k1h_far_post", "k1h_W0_far")
    ratios = (("W_false", "W0_all"), ("W_far", "W0_far"), ("W_near", "W0_near"), ("W_in5", "W0_in5"),
              ("W_far_post", "W0_far"), ("k1h_far_post", "k1h_W0_far"))
    acc_lines = {}
    for des in ("STD", "DENSE"):
        for reg in REGIMES:
            for cell in CELLS:
                if (reg, cell, des) not in data or (reg, cell, "CLUT") not in data:
                    continue
                U, cl = data[(reg, cell, des)], data[(reg, cell, "CLUT")]
                js = sorted(j for j in U if (reg, cell, des, j) in tw)
                if not js:
                    continue
                tq = np.quantile(np.concatenate([u.sT.ravel() for u in cl.values()]), [0.80, 0.95])
                per = {}                                     # net -> array (windows with events, 2 thresholds, keys)
                for j in js:
                    u, rows = U[j], tw[(reg, cell, des, j)]
                    arr = []
                    for i, r in enumerate(rows):
                        if not r["k"]:
                            continue
                        a0 = L.unpack(r["atoms0"])
                        s0 = L.runmin([a["s"] for a in a0])
                        d0 = np.array([min(np.linalg.norm(a["x"][:2] - np.array(t)[:2]) for t in r["truth"]) for a in a0])
                        hl = ~u.false["LOC"][i]
                        after = (np.cumsum(hl) - hl) > 0
                        f, dm = u.false["LOC+T"][i], u.dmin[i]
                        row = []
                        for t in tq:
                            a, b = u.sT[i] >= t, s0 >= t
                            far = f & a & (dm > 10.0)
                            k1h = float(r["k"] == 1 and hl.any())
                            row.append([(f & a).sum(), (f & a & (dm <= 5.0)).sum(), (f & a & (dm > 5.0) & (dm <= 10.0)).sum(),
                                        far.sum(), (far & ~after).sum(), (far & after).sum(), b.sum(), (b & (d0 <= 5.0)).sum(),
                                        (b & (d0 > 5.0) & (d0 <= 10.0)).sum(), (b & (d0 > 10.0)).sum(),
                                        k1h * (far & after).sum(), k1h * (b & (d0 > 10.0)).sum()])
                        arr.append(row)
                    per[j] = np.array(arr, float).reshape(-1, 2, len(keys))
                name = f"{des} {reg} {cell}"
                acc_lines[name] = per
                acc_lines.setdefault(f"{des} {reg} ALL", {}).update({(cell, j): v for j, v in per.items()})
                acc_lines.setdefault(f"{des} both ALL", {}).update({(reg, cell, j): v for j, v in per.items()})
    for name, per in acc_lines.items():
        ks_ = sorted(per, key=str)

        def stat(sel):
            tot = sum(a.sum(0) for a in sel)
            n = sum(len(a) for a in sel)
            return tot / max(1, n)
        pt = stat([per[k] for k in ks_])
        reps = []
        for _ in range(B):
            sel = []
            for i in rng.integers(0, len(ks_), len(ks_)):
                a = per[ks_[i]]
                if len(a):
                    sel.append(a[rng.integers(0, len(a), len(a))])
            reps.append(stat(sel))
        reps = np.array(reps)
        fh.write(f"\n  {name}: {sum(len(per[k]) for k in ks_)} event windows\n")
        summ[name] = {}
        for ki, k in enumerate(keys):
            ci = np.percentile(reps[:, :, ki], [2.5, 97.5], axis=0)
            fh.write(f"      {k:13s} {pt[0, ki]:.3f} [{ci[0, 0]:.3f}, {ci[1, 0]:.3f}] | {pt[1, ki]:.3f} [{ci[0, 1]:.3f}, "
                     f"{ci[1, 1]:.3f}]\n")
            summ[name][k] = [[float(pt[t, ki]), float(ci[0, t]), float(ci[1, t])] for t in (0, 1)]
        for a, b in ratios:
            ia, ib = keys.index(a), keys.index(b)
            v = pt[:, ia] / np.maximum(pt[:, ib], 1e-9)
            ci = np.percentile(reps[:, :, ia] / np.maximum(reps[:, :, ib], 1e-9), [2.5, 97.5], axis=0)
            fh.write(f"      ratio {a}/{b}: {v[0]:.2f} [{ci[0, 0]:.2f}, {ci[1, 0]:.2f}] | {v[1]:.2f} [{ci[0, 1]:.2f}, "
                     f"{ci[1, 1]:.2f}]\n")
            summ[name][f"{a}/{b}"] = [[float(v[t]), float(ci[0, t]), float(ci[1, t])] for t in (0, 1)]
    return summ


def main(B=300, nperm=500):
    data, nets = load()
    with open(os.path.join(OUT, "c13_analysis.txt"), "w") as fh:
        fh.write(f"C13 analysis; B = {B}, permutations = {nperm}\n")
        for k in sorted(nets):
            fh.write(f"  {k}: networks {sorted(nets[k])}\n")
        sec = {}
        for (reg, cell, des), U in sorted(data.items()):
            s = np.concatenate([u.sec for u in U.values()])
            sec[f"{des} {reg} {cell}"] = [len(s), int(sum(u.k.sum() for u in U.values()))] + s.mean(0).tolist()
            fh.write(f"  {des:5s} {reg:6s} {cell:15s} {len(s):4d} windows {sec[f'{des} {reg} {cell}'][1]:4d} events; "
                     f"s per window: target {s[:, 0].mean():.3f}, five decoys {s[:, 1].mean():.3f}\n")
        ex = run_exchange(data, nets, B, nperm, fh)
        evw = run_event_windows(data, nets, B, fh)
        twn = run_twin(data, B, fh)
        fh.write("\n" + "=" * 110 + "\nCALIBRATION\n")
        cal = run_calibration(data, nets, B, fh)
        compact(cal, fh)
    json.dump({"dec": list(DEC), "sec": sec, "exchange": ex, "event_windows": evw, "twin": twn, "cal": cal},
              open(os.path.join(OUT, "c13_summary.json"), "w"))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 300, int(sys.argv[2]) if len(sys.argv) > 2 else 500)
