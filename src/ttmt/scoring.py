"""Hit rules, exact FROC and the bootstrap of the registered analyses.

Hit rules (``match``):
  ``"legacy"``      registered behaviour. Synthetic tasks (c10_lib.label): the nearest true event within 5 km
                    (LOC; no origin-time test) or, with tol_t, within 5 km and tol_t s (LOC+T). Real-data tasks
                    (r0_analyze.label): within 5 km and tol s, events of the target magnitude preferred, nearest. A
                    second detection of an event already hit gets the same label and is neutral in the FROC (it adds
                    neither a hit nor a false detection).
  ``"one_to_one"``  ADR-0007: detections are taken in order; each is matched to the nearest UNMATCHED event within
                    tolerance (5 km and 0.5 s for LOC+T); a detection with no unmatched event left is false (label -1),
                    so duplicates count as false detections.
Matching order (``order``): ``"record"`` for the records of a greedy loop, whose FROC thresholds the running minimum
of the score, so every threshold admits a prefix in loop order; matching in that order keeps the prefix property (a
later atom with a higher raw score must not take an event from an earlier one). ``"score"`` for detection lists
thresholded on their raw scores (R0 and R1), where every threshold admits the detections above it.
Ported from c10_lib (label, record), froc_exact (greedy_record, exact_HF), xc.recall_caps, c10_analyze.Boot,
c16_analyze (records, thr_at, fp_noise, transfer, combine, load), r0_analyze (dedupe, label, Set, froc).
"""
import glob
import json
import os

import numpy as np

MATCH = ("legacy", "one_to_one")
TOL_KM = 5.0
CAPS = (0.1, 0.25)
LOC_T = 0.5                                # s, LOC+T origin tolerance


# ============================================================================================ synthetic hit rules
def label(atoms, truth, taus, tol_t=None, tol_km=TOL_KM, match="legacy", order="record"):
    """labels of atoms (dicts with x, tau, s; c10_lib.label for match legacy), in loop order: index of the true event
    hit, or -1. truth: true sources (3,); taus: their origin times; tol_t None = LOC (epicentre only); order: the
    matching order of one_to_one (default "record": the loop order of a greedy record)."""
    if match == "legacy":
        out = []
        for a in atoms:
            h = -1
            if truth:
                dd = np.array([np.linalg.norm(a["x"][:2] - np.asarray(t)[:2]) for t in truth])
                if tol_t is None:
                    j = int(np.argmin(dd))
                    h = j if dd[j] <= tol_km else -1
                else:
                    ok = (dd <= tol_km) & (np.abs(a["tau"] - np.asarray(taus)) <= tol_t + 1e-9)
                    if ok.any():
                        h = int(np.where(ok, dd, np.inf).argmin())
            out.append(h)
        return out
    if match == "one_to_one":
        xy = np.array([a["x"][:2] for a in atoms], float).reshape(-1, 2)
        t = np.array([a["tau"] for a in atoms], float)
        s = np.array([a["s"] for a in atoms], float)
        ex = np.array([np.asarray(q, float)[:2] for q in truth]).reshape(-1, 2)
        et = np.asarray(taus, float)
        dist = np.linalg.norm(xy[:, None, :] - ex[None, :, :], axis=-1)
        ok = dist <= tol_km
        if tol_t is not None:
            ok &= np.abs(t[:, None] - et[None, :]) <= tol_t + 1e-9
        return one_to_one(s, dist, ok, order=order).tolist()
    raise ValueError(f"match must be one of {MATCH}")


ORDERS = ("record", "score")


def match_order(scores, order):
    """processing order of the detections: "record" = decreasing running minimum of the scores, stable (the loop
    order); "score" = decreasing raw score, stable."""
    s = np.asarray(scores, float)
    if order == "record":
        return np.argsort(-np.minimum.accumulate(s), kind="stable") if len(s) else np.zeros(0, int)
    if order == "score":
        return np.argsort(-s, kind="stable")
    raise ValueError(f"order must be one of {ORDERS}, got {order!r}")


def one_to_one(scores, dist, ok, prefer=None, order="score"):
    """ADR-0007 matching. scores (n,), dist (n, m) distances, ok (n, m) eligibility, prefer (m,) optional bool: among
    eligible unmatched events the preferred ones are taken first (the target-magnitude preference of r0_analyze).
    Detections in match_order(scores, order); returns labels (n,) with -1 for false."""
    n = len(scores)
    lab = np.full(n, -1, int)
    if n == 0 or dist.shape[1] == 0:
        return lab
    used = np.zeros(dist.shape[1], bool)
    for i in match_order(scores, order):
        cand = ok[i] & ~used
        if prefer is not None and (cand & prefer).any():
            cand = cand & prefer
        if cand.any():
            j = int(np.where(cand, dist[i], np.inf).argmin())
            lab[i] = j
            used[j] = True
    return lab


def record(atoms, labels, key="s"):
    """c10_lib.record: (running minimum of the score, label) per atom, in loop order."""
    rec, run = [], np.inf
    for a, h in zip(atoms, labels):
        run = min(run, a[key])
        rec.append((float(run), int(h)))
    return rec


# ============================================================================================ exact FROC
def greedy_record(grid, field, truth, min_sep_km=5.0, max_events=5, tol_km=TOL_KM):
    """froc_exact.greedy_record: the greedy catalogue of a score field with its hits (prefix property)."""
    f = np.where(np.isfinite(field), field, -np.inf)
    order = np.argsort(-f)
    acc, rec = [], []
    for k in order:
        if not np.isfinite(f[k]):
            break
        p = np.asarray(grid[k], float)
        if all(np.linalg.norm(p[:2] - e[:2]) > min_sep_km for e in acc):
            acc.append(p)
            hit = -1
            if truth:
                dd = [float(np.linalg.norm(p[:2] - np.asarray(s)[:2])) for s in truth]
                j = int(np.argmin(dd))
                if dd[j] <= tol_km:
                    hit = j
            rec.append((float(f[k]), hit))
        if len(acc) >= max_events:
            break
    return rec


def exact_HF(records):
    """froc_exact.exact_HF: descending thresholds and (window x threshold) counts of distinct hits H and false F."""
    thr = np.unique(np.concatenate([np.array([s for s, _h in r], float) for r in records if r]
                                   or [np.array([0.0])]))[::-1]
    H = np.zeros((len(records), len(thr)), int)
    F = np.zeros((len(records), len(thr)), int)
    for w, r in enumerate(records):
        hits, nf, pos = set(), 0, 0
        for j, t in enumerate(thr):
            while pos < len(r) and r[pos][0] >= t:
                if r[pos][1] >= 0:
                    hits.add(r[pos][1])
                else:
                    nf += 1
                pos += 1
            H[w, j], F[w, j] = len(hits), nf
    return thr, H, F


def froc_curve(recs):
    """pooled exact FROC of per-window records: score levels (descending), cumulative distinct hits and false
    detections at each level (ties grouped). Duplicate hits of the same (window, event) count once."""
    sc, hit, key = [], [], []
    for wi, r in enumerate(recs):
        for s_, h in r:
            sc.append(s_)
            hit.append(h)
            key.append(wi * 64 + h if h >= 0 else -1)
    if not sc:
        return np.zeros(0), np.zeros(0, int), np.zeros(0, int)
    sc, hit, key = np.array(sc), np.array(hit), np.array(key)
    o = np.argsort(-sc, kind="stable")
    sc, hit, key = sc[o], hit[o], key[o]
    _u, first = np.unique(key, return_index=True)
    newhit = np.zeros(len(sc), bool)
    newhit[first] = True
    newhit &= key >= 0
    H, F = np.cumsum(newhit), np.cumsum(hit < 0)
    last = np.r_[np.where(np.diff(sc) != 0)[0], len(sc) - 1]
    return sc[last], H[last], F[last]


def recall_caps(recs, nte, caps=CAPS):
    """xc.recall_caps: recall at <= cap realised false events per window, and the score threshold."""
    if not any(len(r) for r in recs) or nte == 0:
        return [0.0 for _ in caps], [np.inf for _ in caps]
    th, H, F = froc_curve(recs)
    out, thr = [], []
    for c in caps:
        ok = F / len(recs) <= c
        if ok.any():
            j = int(np.argmax(np.where(ok, H, -1)))
            out.append(float(H[j] / nte))
            thr.append(float(th[j]))
        else:
            out.append(0.0)
            thr.append(np.inf)
    return out, thr


class Boot:
    """c10_analyze.Boot: point recall at the caps and hierarchical bootstrap (strata (regime, cell); networks, then
    windows), operating point re-chosen in every replicate. rows: windows with keys k, regime, cell, net; recs:
    {method: per-window records}; sel: row filter; B replicates; seed 12345."""

    def __init__(self, rows, recs, sel, B, seed=12345, caps=CAPS):
        self.caps = caps
        self.methods = list(recs)
        idx = [i for i, r in enumerate(rows) if sel(r)]
        self.idx = idx
        self.ks = np.array([rows[i]["k"] for i in idx])
        g = {}
        for p, i in enumerate(idx):
            r = rows[i]
            g.setdefault((r["regime"], r["cell"], r["net"]), []).append(p)
        self.groups = list(g.values())
        self.strata = [k[:2] for k in g]
        self.S = {m: [np.array([x[0] for x in recs[m][i]], float) for i in idx] for m in self.methods}
        self.H = {m: [np.array([x[1] for x in recs[m][i]], int) for i in idx] for m in self.methods}
        allw = np.arange(len(idx))
        self.pt = {m: self.recall(m, allw) for m in self.methods}
        for m in self.methods[:3]:
            ref = recall_caps([recs[m][i] for i in idx], int(self.ks.sum()), caps)[0]
            assert np.allclose(ref, self.pt[m]), (m, ref, self.pt[m])
        rng = np.random.default_rng(seed)
        sidx = {}
        for gi, s in enumerate(self.strata):
            sidx.setdefault(s, []).append(gi)
        self.rep = {m: np.zeros((B, len(caps))) for m in self.methods}
        self.ws = []
        for b in range(B):
            ws = []
            for s, gl in sidx.items():
                for gi in np.asarray(gl)[rng.integers(0, len(gl), len(gl))]:
                    gg = np.asarray(self.groups[gi])
                    ws.extend(gg[rng.integers(0, len(gg), len(gg))].tolist())
            ws = np.asarray(ws)
            self.ws.append(ws)
            for m in self.methods:
                self.rep[m][b] = self.recall(m, ws)

    def recall(self, m, ws):
        nte = int(self.ks[ws].sum())
        S, H = self.S[m], self.H[m]
        sc = np.concatenate([S[w] for w in ws])
        if not len(sc) or nte == 0:
            return [0.0] * len(self.caps)
        hit = np.concatenate([H[w] for w in ws])
        pos = np.repeat(np.arange(len(ws)), [len(S[w]) for w in ws])
        key = np.where(hit >= 0, pos * 64 + hit, -1)
        o = np.argsort(-sc, kind="stable")
        sc, hit, key = sc[o], hit[o], key[o]
        _u, first = np.unique(key, return_index=True)
        newhit = np.zeros(len(sc), bool)
        newhit[first] = True
        newhit &= key >= 0
        Hc, Fc = np.cumsum(newhit), np.cumsum(hit < 0)
        last = np.r_[np.where(np.diff(sc) != 0)[0], len(sc) - 1]
        Hc, Fc = Hc[last], Fc[last]
        out = []
        for c in self.caps:
            ok = Fc / len(ws) <= c
            out.append(float(Hc[ok].max() / nte) if ok.any() else 0.0)
        return out

    def ci(self, m):
        return np.percentile(self.rep[m], [2.5, 97.5], axis=0).T


# ============================================================================================ part files (C16, 2M)
def load_parts(pdir, arm):
    """c16_analyze.load: the part files ARM_*.json of a folder; rows carry regime = arm, cell, net."""
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
    """c16_analyze.combine: pair the windows of two arms (same cell, network, design, window)."""
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


def atoms_of(row, name):
    """c10_analyze.atoms_of: the stored atoms [x, y, z, tau, s, s0, sref] of one method in one window."""
    return [{"x": np.array(a[:3]), "tau": a[3], "s": a[4], "s0": a[5], "sref": a[6] if a[6] == a[6] else a[4]}
            for a in row["atoms"][name]]


def records(rows, names, tol_t, kcut=None, match="legacy"):
    """c16_analyze.records: per method, the per-window records under the hit rule."""
    recs = {}
    for nm in names:
        out = []
        for r in rows:
            at = atoms_of(r, nm)[:kcut]
            out.append(record(at, label(at, [np.array(t) for t in r["truth"]], r["taus"], tol_t, match=match), "s"))
        recs[nm] = out
    return recs


def thr_at(recs_w, nwin, cap):
    """c16_analyze.thr_at: score threshold of the exact-FROC operating point at <= cap false events per window."""
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


def transfer(ev_rows, no_rows, names, tol, kcut, B, seed, cells, caps=CAPS, match="legacy"):
    """c16_analyze.transfer (without the cross arm): threshold of the operating point on event windows and realised
    false events per noise-only window; bootstrap of networks within cell, then windows."""
    rec_e = records(ev_rows, names, tol, match=match)
    rec_n = records(no_rows, names, None, kcut, match=match)
    ge, gn = {}, {}
    for i, r in enumerate(ev_rows):
        ge.setdefault((r["cell"], r["net"]), []).append(i)
    for i, r in enumerate(no_rows):
        gn.setdefault((r["cell"], r["net"]), []).append(i)
    keys = sorted(set(ge) & set(gn))
    rng = np.random.default_rng(seed)
    draws = []
    for _b in range(B):
        we, wn = [], []
        for cell in cells:
            kc = [k for k in keys if k[0] == cell]
            for k in [kc[i] for i in rng.integers(0, len(kc), len(kc))]:
                e, n_ = np.asarray(ge[k]), np.asarray(gn[k])
                we += e[rng.integers(0, len(e), len(e))].tolist()
                wn += n_[rng.integers(0, len(n_), len(n_))].tolist()
        draws.append((we, wn, []))
    res = {}
    for nm in names:
        out = []
        for cap in caps:
            t0 = thr_at(rec_e[nm], len(ev_rows), cap)
            p = (t0, fp_noise(rec_n[nm], t0), np.nan)
            bs = []
            for we, wn, _wc in draws:
                t = thr_at([rec_e[nm][i] for i in we], len(we), cap)
                bs.append((fp_noise([rec_n[nm][i] for i in wn], t), np.nan))
            bs = np.array(bs)
            out.append((p, np.nanpercentile(bs, [2.5, 97.5], axis=0)))
        res[nm] = out
    return res


# ============================================================================================ catalogue hit rules
def dedupe(atoms, dkm=5.0, dts=1.0):
    """r0_analyze.dedupe: pooled atoms [x, y, z, t, s, ...] in decreasing score; one within dkm and dts s of a kept
    atom is removed."""
    a = sorted(atoms, key=lambda r: -r[4])
    keep = []
    for r in a:
        if all(np.hypot(r[0] - k[0], r[1] - k[1]) > dkm or abs(r[3] - k[3]) > dts for k in keep):
            keep.append(r)
    return keep


def label_catalogue(det, cat, tol, m_target, tol_km=TOL_KM, match="legacy", order="score"):
    """detections [x, y, z, t, s] against catalogue events {x, y, t, M}: index of the event hit or -1. Eligible:
    epicentre within tol_km and |origin error| <= tol s; events of M >= m_target preferred; nearest epicentre
    (r0_analyze.label for match legacy; ADR-0007 one-to-one matching otherwise, by default in raw-score order since
    these lists are thresholded on their raw scores)."""
    ex = np.array([[e["x"], e["y"], e["t"]] for e in cat])
    big = np.array([e["M"] >= m_target for e in cat])
    if match == "legacy":
        out = []
        for d in det:
            dd = np.hypot(ex[:, 0] - d[0], ex[:, 1] - d[1])
            ok = (dd <= tol_km) & (np.abs(ex[:, 2] - d[3]) <= tol + 1e-9)
            if (ok & big).any():
                ok = ok & big
            out.append(int(np.where(ok, dd, np.inf).argmin()) if ok.any() else -1)
        return np.array(out, int)
    if match == "one_to_one":
        if not len(det):
            return np.zeros(0, int)
        D = np.asarray([d[:5] for d in det], float)
        dist = np.hypot(D[:, None, 0] - ex[None, :, 0], D[:, None, 1] - ex[None, :, 1])
        ok = (dist <= tol_km) & (np.abs(ex[None, :, 2] - D[:, None, 3]) <= tol + 1e-9)
        return one_to_one(D[:, 4], dist, ok, prefer=big, order=order)
    raise ValueError(f"match must be one of {MATCH}")


class DetSet:
    """r0_analyze.Set: one detection list of one method in a segment [lo, hi): scores, 60 s block of each detection,
    labels per rule (rules: {name: tol s})."""

    def __init__(self, det, cat, lo, hi, m_target, rules, block=60.0, match="legacy"):
        det = [d for d in det if lo <= d[3] < hi]
        self.s = np.array([d[4] for d in det], float)
        self.blk = np.array([int((d[3] - lo) // block) for d in det], int)
        self.lab = {u: label_catalogue(det, cat, tol, m_target, match=match) if det else np.zeros(0, int)
                    for u, tol in rules.items()}
        self.det = det


def froc_blocks(sets, targets_mask, M, wb, rate, rule, m_target, block=60.0):
    """r0_analyze.froc: recall at <= rate false detections per hour for a method given as a list of DetSets (the
    envelope over configurations); wb: weight (multiplicity) of every block in this replicate. A detection hitting
    only events below m_target is neutral."""
    hours = wb.sum() * block / 3600.0
    cap = int(np.floor(rate * hours + 1e-9))
    tw = wb[targets_mask[1]]
    ntar = tw.sum()
    best = 0.0
    for S in sets:
        if not len(S.s):
            continue
        w = wb[S.blk]
        lab = S.lab[rule]
        neutral = (lab >= 0) & (M[np.maximum(lab, 0)] < m_target)
        fal = lab < 0
        fs, fw = S.s[fal], w[fal]
        o = np.argsort(-fs)
        cum = np.cumsum(fw[o])
        k = np.searchsorted(cum, cap, side="right")
        thr = fs[o][k] if k < len(fs) else -np.inf
        hit = (lab >= 0) & ~neutral & (S.s > thr)
        got = np.zeros(len(M), bool)
        got[lab[hit]] = True
        rec = float((got[targets_mask[0]] * tw).sum() / max(ntar, 1e-12))
        best = max(best, rec)
    return best
