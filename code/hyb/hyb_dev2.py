"""
hyb_dev2.py - second hybrid on Ridgecrest CAL (development only): pick-free candidates, pick evidence.

Candidates = associated events of a pick-based associator (PyOcto, fixed configuration) and pick-free atoms of the C18j
detector (layered medium, station terms over +-2 s, 2 Hz, no cap), merged (a candidate within 5 km and 1.5 s of a
higher-ranked one is dropped).  Every candidate is scored by the PhaseNet picks it explains: for each station and phase
the predicted arrival (layered first-arrival time from the candidate's location and origin, plus a pick-based station
term) is compared with the PhaseNet picks of that phase; a pick within TOL s counts once (each pick may support only
one candidate: candidates are processed in decreasing order of their provisional score and claim their picks).
Score = number of supporting picks + 0.5 x percentile of the waveform score of the candidate (tie-breaker).
Pick-based station terms: median over the CAL catalogue events (M >= mmin, 1.5 on Ridgecrest) of (nearest PhaseNet pick of the phase within
+-2 s of the layered prediction - layered prediction), per station and phase.
    python hyb_dev2.py SUBNET P1_ATOMS.json OUT.txt [TOL [BASE_FILE:KEY BASE_CONFIG]]
"""
import json
import sys

import numpy as np

sys.path[:0] = ["/home/claude/work/r0", "/home/claude/work/explore_contrib", "/home/claude/work"]
import r0_lib as R  # noqa: E402
import r0_analyze as A  # noqa: E402
import r0_run as RR  # noqa: E402

import c18_tt as TT  # noqa: E402
import c18j_layered as LJ  # noqa: E402
import hyb_dev as H  # noqa: E402

CFG = '{"n_picks": 3, "tol": 0.5}'


def times(src, STA, cache={}):
    key = round(float(src[2]), 1)
    if key not in cache:
        cache[key] = TT.table(LJ.MODEL, LJ.R_VPVS, [key], xmax=120.0, dx=0.25)[key]
    xs, tp, ts = cache[key]
    X = np.hypot(STA[:, 0] - src[0], STA[:, 1] - src[1])
    return np.stack([np.interp(X, xs, tp), np.interp(X, xs, ts)])          # (2, ns)


def pick_arrays(picks):
    out = []
    for pl in picks:
        out.append({b: np.array([t for t, ph, _h in pl if ph == b]) for b in "PS"})
    return out


def pick_terms(net, pa, lo, hi, mmin=1.5):
    res = [[[], []] for _ in range(len(net.STA))]
    for e in net.cat:
        if not (lo <= e["t"] < hi and e["M"] >= mmin):
            continue
        tt = times([e["x"], e["y"], e["z"]], net.STA)
        for i in range(len(net.STA)):
            for b, ph in enumerate("PS"):
                p = pa[i][ph]
                d = p - (e["t"] + tt[b, i])
                d = d[np.abs(d) <= 2.0]
                if len(d):
                    res[i][b].append(float(d[np.argmin(np.abs(d))]))
    return np.array([[np.median(x) if len(x) >= 2 else 0.0 for x in r] for r in res])   # (ns, 2)


def score_candidates(cands, STA, pa, terms, tol):
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


def main(sub, p1path, outp, tol=0.5, base_src=None, base_cfg=None):
    lo, hi = R.T_CAL
    net = RR.Net(sub, 940.0)
    tab = LJ.tt_table(sorted({float(z) for z in net.eng["SC"][0].XN[:, 2]}))
    E, _corr = LJ.layered_engine(net, tab)
    cat = net.cat
    pk = json.load(open(f"/home/claude/work/hyb/picks_cal_{sub}.json"))
    pa = pick_arrays(pk["picks"])
    terms = pick_terms(net, pa, lo, hi)
    if base_src:
        bs, bk = base_src.split(":")
        base = [e[:5] for e in json.load(open(bs))[bk][base_cfg] if lo <= e[3] < hi]
    else:
        base = [e[:5] for e in pk["pyocto"][CFG] if lo <= e[3] < hi]
    p1 = A.dedupe([a[:5] for a in json.load(open(p1path))["atoms"]["NOCAP"] if lo <= a[3] < hi])
    wb = H.pct(H.waveform_scores(net, E, base, "NOCAP"))
    wp = H.pct([a[4] for a in p1])
    lines = [f"{sub}: {len(base)} associated, {len(p1)} pick-free; pick terms P {np.round(terms[:, 0], 2).tolist()} "
             f"S {np.round(terms[:, 1], 2).tolist()}; TOL {tol}"]
    res = {"BASE (PyOcto)": H.evaluate(base, cat, lo, hi), "P1": H.evaluate(p1, cat, lo, hi)}
    for tag, cands in (("PICKSCORE-P1", [a + [wp[i]] for i, a in enumerate(p1)]),
                       ("PICKSCORE-BASE", [e[:4] + [e[4], wb[i]] for i, e in enumerate(base)]),
                       ("PICKSCORE-UNION", [e[:4] + [e[4] + 100.0, wb[i]] for i, e in enumerate(base)] +
                        [a[:4] + [wp[i], wp[i]] for i, a in enumerate(p1)])):
        sc = score_candidates(cands, net.STA, pa, terms, tol)
        dets = [[c[0], c[1], c[2], c[3], c[4] + 0.5 * c[5]] for c in sc]
        dets.sort(key=lambda r: -r[4])
        keep = []
        for r_ in dets:
            if all(np.hypot(r_[0] - k[0], r_[1] - k[1]) > 5.0 or abs(r_[3] - k[3]) > 1.5 for k in keep):
                keep.append(r_)
        res[tag] = H.evaluate(keep, cat, lo, hi)
    for k, v in res.items():
        lines.append(f"   {k:18s} " + "  ".join(f"{u} {r}/h {v[(u, r)]:.3f}" for u in ("LOC+T", "LOC") for r in (10, 30, 60)))
    open(outp, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]) if len(sys.argv) > 4 else 0.5,
         sys.argv[5] if len(sys.argv) > 5 else None, sys.argv[6] if len(sys.argv) > 6 else None)
