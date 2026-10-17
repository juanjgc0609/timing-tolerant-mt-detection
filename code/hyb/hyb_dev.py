"""
hyb_dev.py - development of a hybrid detector on the Ridgecrest CAL segment (origins before 19:15; EVAL not read).

The hybrid combines a pick-based associator (PhaseNet picks + PyOcto or GaMMA, frozen r0_picks.py outputs) with the
pick-free moment-tensor evidence of the C18j engine (layered medium, station terms read over +-2 s, 2 Hz):
  1. Rescoring: every associated event gets a waveform score, the maximum of the first-iteration subspace statistic
     over coarse nodes within 3 km (epicentre) of its location and origins within 0.5 s of its origin time, in the
     40 s window that starts 6.5 s before its origin.  Its fused score is the sum of two percentile ranks within the
     segment: of its number of associated picks and of its waveform score.
  2. Fill-in: pick-free atoms (C18j detector) farther than 5 km or 1.5 s from every associated event are added with
     the percentile rank of their own pick-free score (no pick term), so they rank below associated events with the
     same waveform evidence.
Variants: RESCORE (1 only), FILL (2 only, associator score = percentile of picks), HYB (1 and 2); statistic of the
waveform score NOCAP or LOG6.  Output: recall on CAL at 10/30/60 false detections per hour, LOC+T and LOC, for the
base associator configuration, BESTPICK-like envelope, P1 alone and the hybrids.
    python hyb_dev.py SUBNET ASSOC CONFIG_JSON P1_ATOMS.json OUT.txt
"""
import json
import sys

import numpy as np

sys.path[:0] = ["/home/claude/work/r0", "/home/claude/work/explore_contrib", "/home/claude/work"]
import r0_lib as R  # noqa: E402
import r0_analyze as A  # noqa: E402
import r0_run as RR  # noqa: E402

import c18j_layered as LJ  # noqa: E402


def pct(x):
    x = np.asarray(x, float)
    if len(x) == 0:
        return x
    order = np.argsort(np.argsort(x, kind="stable"), kind="stable")
    return (order + 1) / len(x)


def waveform_scores(net, E, events, stat):
    import c10_lib as C
    import c18c_oracle as OC
    import xc
    out = []
    for ev in events:
        ts = ev[3] - 6.5
        W = R.window(net.d16, max(ts, 0.0))
        MF, sd = xc.matched(W)
        kind, c = {"NOCAP": ("none", 0.0), "LOG6": ("log", 6.0)}[stat]
        s = OC.score_map(E, MF / sd, kind, c)
        d = np.hypot(E.XN[:, 0] - ev[0], E.XN[:, 1] - ev[1])
        near = d <= 3.0
        tt = np.abs(max(ts, 0.0) + E.taus - ev[3]) <= 0.5
        out.append(float(s[np.ix_(near, tt)].max()) if near.any() and tt.any() else 0.0)
    return np.array(out)


def evaluate(dets, cat, lo, hi):
    M = np.array([e["M"] for e in cat])
    tg = np.array([i for i, e in enumerate(cat) if lo <= e["t"] < hi and e["M"] >= R.M_TARGET])
    tb = np.array([int((cat[i]["t"] - lo) // A.BLOCK) for i in tg])
    nb = int(np.ceil((hi - lo) / A.BLOCK))
    S = [A.Set(dets, cat, lo, hi)]
    return {(u, r): A.froc(S, (tg, tb), M, nb, np.ones(nb), r, u) for u in ("LOC+T", "LOC") for r in (10, 30, 60)}


def main(sub, assoc, cfg, p1path, outp):
    lo, hi = R.T_CAL
    net = RR.Net(sub, 940.0)
    tab = LJ.tt_table(sorted({float(z) for z in net.eng["SC"][0].XN[:, 2]}))
    E, _corr = LJ.layered_engine(net, tab)
    cat = net.cat
    pk = json.load(open(f"/home/claude/work/hyb/picks_cal_{sub}.json"))
    base = [e for e in pk[assoc][cfg] if lo <= e[3] < hi]
    p1 = json.load(open(p1path))["atoms"]["NOCAP"]
    p1 = A.dedupe([a for a in p1 if lo <= a[3] < hi])
    lines = [f"{sub} {assoc} {cfg}: {len(base)} associated events, {len(p1)} pick-free atoms (after de-duplication)"]
    res = {"BASE": evaluate([[e[0], e[1], e[2], e[3], e[4]] for e in base], cat, lo, hi),
           "P1": evaluate(p1, cat, lo, hi)}
    env = {}
    for key in ("gamma", "pyocto", "gamma_sc", "pyocto_sc"):
        for c, v in pk[key].items():
            r_ = evaluate([e[:5] for e in v if lo <= e[3] < hi], cat, lo, hi)
            for k, x in r_.items():
                env[k] = max(env.get(k, 0.0), x)
    res["ENVELOPE(all assoc, cfgs)"] = env
    far = [a for a in p1 if all(np.hypot(a[0] - e[0], a[1] - e[1]) > 5.0 or abs(a[3] - e[3]) > 1.5 for e in base)]
    pp = pct([e[4] for e in base])
    for stat in ("NOCAP", "LOG6"):
        ws = waveform_scores(net, E, base, stat)
        pw = pct(ws)
        resc = [[e[0], e[1], e[2], e[3], float(pp[i] + pw[i])] for i, e in enumerate(base)]
        res[f"RESCORE-{stat}"] = evaluate(resc, cat, lo, hi)
        pf = pct([a[4] for a in far]) if far else []
        fill = [[a[0], a[1], a[2], a[3], float(pf[i])] for i, a in enumerate(far)]
        res[f"FILL-{stat}"] = evaluate([[e[0], e[1], e[2], e[3], float(1.0 + pp[i])] for i, e in enumerate(base)] + fill,
                                       cat, lo, hi)
        res[f"HYB-{stat}"] = evaluate(resc + fill, cat, lo, hi)
    # rank union: each source scored by its own within-segment percentile; cross-source duplicates (5 km, 1.5 s) keep
    # the higher-scored member
    pp1 = pct([a[4] for a in p1])
    for tag, bs in (("UNION-RANK", [float(x) for x in pp]),):
        cand = [[e[0], e[1], e[2], e[3], bs[i]] for i, e in enumerate(base)] + \
               [[a[0], a[1], a[2], a[3], float(pp1[i])] for i, a in enumerate(p1)]
        cand.sort(key=lambda r: -r[4])
        keep = []
        for r_ in cand:
            if all(np.hypot(r_[0] - k[0], r_[1] - k[1]) > 5.0 or abs(r_[3] - k[3]) > 1.5 for k in keep):
                keep.append(r_)
        res[tag] = evaluate(keep, cat, lo, hi)
    lines.append(f"   fill-in candidates: {len(far)}")
    for k, v in res.items():
        lines.append(f"   {k:26s} " + "  ".join(f"{u} {r}/h {v[(u, r)]:.3f}" for u in ("LOC+T", "LOC") for r in (10, 30, 60)))
    open(outp, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:6])
