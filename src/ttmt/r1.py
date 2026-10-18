"""Task R1 through the package: network setup, anchored energy-peak terms, layered JMA2001 engine, P1 detector
(MT-FW without cap), hybrid detection lists and analysis. Port of code/r1/r1_lib.py, r1_run.py, r1_hybrid.py,
hyb_dev.py (pct, waveform_scores, evaluate) and hyb_dev2.py (score_candidates).

The pick stage (PhaseNet, GaMMA, PyOcto) is not part of the core package: its output ``picks_all.json`` is an input.

    python -m ttmt.r1 p1 SEGMENT PICKS.json OUT.json [--jobs 16] [--strict-partition]
    python -m ttmt.r1 hybrid SEGMENT PICKS.json P1.json OUT.json
    python -m ttmt.r1 all PICKS.json OUTDIR [--jobs 16]          (p1 and hybrid of cal and eval, then the analysis)

Partition (E56): the legacy run read each segment up to its end plus 40 s (``legacy_t_end``), so the CAL reads reached
the first 40 s of EVAL. That is the default here, because it changes the filtered waveforms and so the registered
numbers. ``--strict-partition`` ends every read at the segment boundary; its numbers are a new, separately labelled
result.
"""
import json
import os
import sys
import time

import numpy as np

from . import scoring as S
from . import traveltime as TT
from .datasets import DT, HINET, NT
from .detector import NOCAP, Engine, mt_loop, score_map
from .kernel import LEGACY as KERNEL
from .partition import legacy_t_end
from .radiation import node_operators
from .station_terms import OnsetDelay, PeakDelay, anchored_peak_terms, estimate_peak_terms, onset_terms, pick_arrays

R = HINET
MODEL = TT.JMA2001
TOL = 0.5          # s, pick support of the hybrid
W_WAVE = 0.5       # weight of the waveform-score percentile
DKM, DTS = 5.0, 1.5
RULE_TOL = {"LOC+T": 0.5, "LOC": 3.0}


def coarse_grid(depths, step=3.0, half=21.0):
    """xc.make_grid: epicentral lattice offset by 0.37 km (truth never on a node), at the given depths."""
    ax = np.arange(-half, half + 1e-9, step) + 0.37
    return np.array([[x, y, z] for x in ax for y in ax for z in depths], float)


class Net:
    """r1_lib.make_net / r0_run.Net, the parts the R1 detector uses: stations, normalised 16 Hz traces, calibration
    windows of the CAL truth events with M >= 1.0 (origin at 6.5 s of window time), homogeneous peak terms (used only
    by the descriptive R0-style comparators), coarse lattice with its radiation operators and bases."""

    def __init__(self, segment, strict_partition=False):
        self.names = R.subnets["H7"]
        self.centre = R.frame()
        self.cat = R.catalogue_local(self.centre)
        self.STA, sens = R.stations(self.names, self.centre)
        lo, hi = R.segment(segment)
        t_end = hi if strict_partition else legacy_t_end(R.partition, segment)
        self.t_end = t_end
        self.raw, d16 = R.load_waveforms(self.names, sens, t_end)
        self.scale = R.station_scale(d16)
        self.d16 = d16 / self.scale[:, None, None]
        self.calev = [e for e in self.cat if R.cal[0] <= e["t"] < R.cal[1] and e["M"] >= R.m_terms]
        self.calwin = []
        for e in self.calev:
            ts = e["t"] - 6.5
            if ts < 0:
                continue
            self.calwin.append((R.window(self.d16, ts), np.array([e["x"], e["y"], e["z"]]), 6.5))
        self.corr_sc = estimate_peak_terms(self.STA, self.calwin, R.vp, R.vs)
        self.XN = coarse_grid(R.depths)
        self.Gs, self.Qst = node_operators(self.XN, self.STA)


def model_times(src, STA):
    return TT.times(MODEL, TT.R_VPVS, src, STA)


def engine(net, onset):
    """r1_run.engine: layered JMA2001 first-arrival tables of the coarse nodes and the energy-peak terms anchored on
    the onset terms."""
    tab = TT.tt_table(MODEL, sorted({float(z) for z in net.XN[:, 2]}))
    corr = anchored_peak_terms(net.calwin, net.STA, onset, model_times, R.vp, R.vs)
    tt = TT.tt_nodes(net.XN, net.STA, tab)
    E = Engine(net.XN, net.STA, tt, R.taus, corr, net.Gs, net.Qst, dt=DT, nt=NT)
    return E, corr


def p1_windows(net, E, ws, idx):
    """P1 atoms owned by the windows idx: [x, y, z, origin (record s), score, window index]."""
    out = []
    for wi in idx:
        ts = ws[wi]
        MF, sd = KERNEL.matched(R.window(net.d16, ts))
        for at in mt_loop(E, MF / sd, R.kmax, capv=None):
            if R.own[0] - 1e-9 <= at["tau"] < R.own[1] - 1e-9:
                x = at["x"]
                out.append([float(x[0]), float(x[1]), float(x[2]), float(ts + at["tau"]), float(at["s"]), wi])
    return out


_W = {}


def _worker_init(segment, onset_vals, strict):
    from .runtime import pin_blas
    pin_blas(strict=False)
    net = Net(segment, strict)
    E, _ = engine(net, OnsetDelay(np.array(onset_vals)))
    _W.update(net=net, E=E, ws=R.starts(segment))


def _worker_chunk(rng):
    t0 = time.time()
    return rng, p1_windows(_W["net"], _W["E"], _W["ws"], range(*rng)), time.time() - t0


def run_p1(segment, picks, jobs=1, chunk=40, strict=False, log=print):
    """r1_run.run + merge: the P1 detection list of a segment (dict as p1_<segment>.json)."""
    onset = OnsetDelay(np.array(picks["pick_terms"]))
    net = Net(segment, strict)
    E, corr = engine(net, onset)
    ws = R.starts(segment)
    meta = {"lag": corr.info["lag"].tolist(), "lag_readings": corr.info["lag_readings"], "segment": segment,
            "names": net.names, "centre": net.centre, "STA": net.STA.tolist(), "scale": net.scale.tolist(),
            "corr_layered": corr.tolist(), "n_calwin": len(net.calwin), "n_windows": len(ws)}
    todo = [(a, min(a + chunk, len(ws))) for a in range(0, len(ws), chunk)]
    res = {}
    if jobs <= 1:
        for a, b in todo:
            t0 = time.time()
            res[(a, b)] = (p1_windows(net, E, ws, range(a, b)), time.time() - t0)
            log(f"{time.strftime('%H:%M:%S')} {segment} chunk {a} {b} of {len(ws)}")
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(jobs, initializer=_worker_init, initargs=(segment, onset.values.tolist(), strict)) as ex:
            for rng, atoms, sec in ex.map(_worker_chunk, todo):
                res[rng] = (atoms, sec)
                log(f"{time.strftime('%H:%M:%S')} {segment} chunk {rng[0]} {rng[1]} of {len(ws)}")
    out = dict(meta)
    out["atoms"] = {"NOCAP": [a for k in sorted(res) for a in res[k][0]]}
    out["sec"] = float(sum(v[1] for v in res.values()))
    return out


# ============================================================================================ hybrid
def pct(x):
    """hyb_dev.pct: within-list percentile rank (1 / n ... 1), ties broken by position."""
    x = np.asarray(x, float)
    if len(x) == 0:
        return x
    order = np.argsort(np.argsort(x, kind="stable"), kind="stable")
    return (order + 1) / len(x)


def waveform_scores(net, E, events, radial=NOCAP):
    """hyb_dev.waveform_scores: per event, the maximum first-iteration statistic over coarse nodes within 3 km
    (epicentre) and origins within 0.5 s, in the 40 s window starting 6.5 s before its origin."""
    out = []
    for ev in events:
        ts = ev[3] - 6.5
        W = R.window(net.d16, max(ts, 0.0))
        MF, sd = KERNEL.matched(W)
        s = score_map(E, MF / sd, radial)
        d = np.hypot(E.XN[:, 0] - ev[0], E.XN[:, 1] - ev[1])
        near = d <= 3.0
        tt = np.abs(max(ts, 0.0) + E.taus - ev[3]) <= 0.5
        out.append(float(s[np.ix_(near, tt)].max()) if near.any() and tt.any() else 0.0)
    return np.array(out)


def score_candidates(cands, STA, pa, onset, tol, times):
    """hyb_dev2.score_candidates: every candidate [x, y, z, t, provisional score, waveform pct], in decreasing
    provisional score, claims the unused picks within tol s of its predicted arrivals (layered times + ONSET terms);
    returns [x, y, z, t, number of supporting picks, waveform pct]."""
    terms = onset.values if isinstance(onset, OnsetDelay) else None
    if terms is None:
        raise TypeError(f"pick support needs OnsetDelay terms, got {type(onset).__name__} (E54)")
    used = [{ph: np.zeros(len(pa[i][ph]), bool) for ph in "PS"} for i in range(len(STA))]
    out = []
    for c in sorted(cands, key=lambda r: -r[4]):
        tt = times(c[:3], STA)
        n = 0
        for i in range(len(STA)):
            for b, ph in enumerate("PS"):
                p = pa[i][ph]
                if not len(p):
                    continue
                d = np.abs(p - (c[3] + tt[b, i] + terms[i, b]))
                d[used[i][ph]] = np.inf
                k = int(np.argmin(d))
                if d[k] <= tol:
                    used[i][ph][k] = True
                    n += 1
        out.append([c[0], c[1], c[2], c[3], n, c[5]])
    return out


def merge_hybrid(base, p1, wb, wp, STA, pa, onset, times):
    """r1_hybrid.hybrid: PICKSCORE-UNION of the associator events (provisional score + 100) and the P1 atoms, final
    score = supporting picks + 0.5 x waveform percentile, then a candidate within 5 km and 1.5 s of a higher one is
    dropped."""
    cands = [e[:4] + [e[4] + 100.0, wb[i]] for i, e in enumerate(base)] + [a[:4] + [wp[i], wp[i]] for i, a in enumerate(p1)]
    sc = score_candidates(cands, STA, pa, onset, TOL, times)
    dets = sorted([[c[0], c[1], c[2], c[3], c[4] + W_WAVE * c[5]] for c in sc], key=lambda r: -r[4])
    keep = []
    for r in dets:
        if all(np.hypot(r[0] - k[0], r[1] - k[1]) > DKM or abs(r[3] - k[3]) > DTS for k in keep):
            keep.append(r)
    return keep


def in_seg(evs, lo, hi):
    return [list(e[:5]) for e in evs if lo <= e[3] < hi]


def evaluate(dets, cat, lo, hi):
    """hyb_dev.evaluate: recall at 10, 30, 60 false per hour, LOC+T and LOC, full weights."""
    M = np.array([e["M"] for e in cat])
    tg = np.array([i for i, e in enumerate(cat) if lo <= e["t"] < hi and e["M"] >= R.m_target])
    tb = np.array([int((cat[i]["t"] - lo) // 60.0) for i in tg])
    nb = int(np.ceil((hi - lo) / 60.0))
    St = [S.DetSet(dets, cat, lo, hi, R.m_target, RULE_TOL)]
    return {(u, r): S.froc_blocks(St, (tg, tb), M, np.ones(nb), r, u, R.m_target) for u in ("LOC+T", "LOC")
            for r in (10, 30, 60)}


def select_base(pk, cat):
    """r1_hybrid.select_base: the 1-D associator configuration with the highest CAL recall at LOC+T 30/h; ties: fewer
    CAL detections, then PyOcto, then the first of the grid."""
    from .associators import GAMMA_GRID, PYOCTO_GRID
    rows = []
    for fam, (key, grid) in enumerate((("pyocto_1d", PYOCTO_GRID), ("gamma_1d", GAMMA_GRID))):
        for gi, cfg in enumerate(grid):
            c = json.dumps(cfg)
            ev = in_seg(pk[key][c], *R.cal)
            rec = evaluate(ev, cat, *R.cal)[("LOC+T", 30)]
            rows.append((-rec, len(ev), fam, gi, key, c))
    rows.sort()
    return rows[0][4], rows[0][5], [(k, c, -r, n) for r, n, _f, _g, k, c in rows]


def run_hybrid(segment, pk, p1d, strict=False):
    """r1_hybrid.main: the detection lists of a segment (dict as hyb_<segment>.json)."""
    lo, hi = R.segment(segment)
    net = Net(segment, strict)
    onset = OnsetDelay(np.array(pk["pick_terms"]))
    E = engine(net, onset)[0]
    cat = net.cat
    bkey, bcfg, table = select_base(pk, cat)
    base = in_seg(pk[bkey][bcfg], lo, hi)
    p1 = S.dedupe(in_seg(p1d["atoms"]["NOCAP"], lo, hi))
    wb = pct(waveform_scores(net, E, base)).tolist()
    wp = pct([a[4] for a in p1]).tolist()
    pa = pick_arrays([[tuple(p) for p in pl] for pl in pk["picks"]])
    times = TT.SourceTimes(MODEL)
    hyb = merge_hybrid(base, p1, wb, wp, net.STA, pa, onset, times)
    return {"segment": segment, "base": [bkey, bcfg], "selection": table, "BASE": base, "HYB": hyb, "P1": p1,
            "BESTPICK": {f"{k}|{c}": in_seg(v, lo, hi) for k in ("pyocto_1d", "gamma_1d") for c, v in pk[k].items()},
            "R0STYLE": {f"{k}|{c}": in_seg(v, lo, hi) for k in ("gamma", "pyocto", "gamma_sc", "pyocto_sc")
                        for c, v in pk[k].items()}}


def recompute_onset_terms(pk, segment="eval", strict=False):
    """r1_picks: the onset terms from the stored picks (hyb_dev2.pick_terms over the CAL truth events, M >= 1.0)."""
    net = Net(segment, strict)
    pa = pick_arrays([[tuple(p) for p in pl] for pl in pk["picks"]])
    return onset_terms(net.cat, net.STA, pa, TT.SourceTimes(MODEL), R.cal[0], R.cal[1], R.m_terms), net


def _dump(obj, path):
    with open(path + ".tmp", "w") as fh:
        json.dump(obj, fh)
    os.replace(path + ".tmp", path)


def main(argv=None):
    from .runtime import pin_blas
    pin_blas(strict=False)
    a = list(sys.argv[1:] if argv is None else argv)
    jobs, strict = 1, False
    if "--jobs" in a:
        i = a.index("--jobs")
        jobs = int(a[i + 1])
        del a[i:i + 2]
    if "--strict-partition" in a:
        a.remove("--strict-partition")
        strict = True
    cmd = a[0]
    if cmd == "p1":
        _dump(run_p1(a[1], json.load(open(a[2])), jobs, strict=strict), a[3])
    elif cmd == "hybrid":
        _dump(run_hybrid(a[1], json.load(open(a[2])), json.load(open(a[3])), strict), a[4])
    elif cmd == "all":
        from .analysis.r1 import analyze
        pk = json.load(open(a[1]))
        out = a[2]
        os.makedirs(out, exist_ok=True)
        for seg in ("cal", "eval"):
            t0 = time.time()
            p = os.path.join(out, f"p1_{seg}.json")
            if not os.path.exists(p):
                _dump(run_p1(seg, pk, jobs, strict=strict), p)
            print(f"P1 {seg} {time.time() - t0:.0f} s", flush=True)
            h = os.path.join(out, f"hyb_{seg}.json")
            if not os.path.exists(h):
                _dump(run_hybrid(seg, pk, json.load(open(p)), strict), h)
            print(f"HYB {seg} {time.time() - t0:.0f} s", flush=True)
        txt = analyze("eval", json.load(open(os.path.join(out, "hyb_eval.json"))),
                      json.load(open(os.path.join(out, "hyb_cal.json"))))
        with open(os.path.join(out, "analysis_eval.txt"), "w", newline="\n") as fh:
            fh.write(txt)
        print(txt, end="")
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
