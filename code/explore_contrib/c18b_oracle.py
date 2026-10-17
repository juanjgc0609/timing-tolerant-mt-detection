"""
c18b_oracle.py - exploration C18b (after R0; Ridgecrest CAL segment only, origins in [0, 900) s; EVAL is not read).
Oracle diagnostic of the R0 failure on real data: separate the score (ranking) from the location (timing) failure.
    python c18b_oracle.py SUBNET OUT.json
For every 40 s window of CAL (start every 3 s, as R0) the first-iteration score map of MT-FW (cap 6, +SC, hard
+-0.25 s read; the registered engine of r0_run.Net) is computed over all coarse nodes and origins, for four input
normalisations / statistics:
  REG     registered: MF divided by ONE robust sd per window (all stations), cap 6, subspace energy
  STA     MF of each station divided by its own robust sd in the window, cap 6, subspace energy
  STA-COH as STA, score = subspace energy x (subspace energy / capped energy at the read samples) (coherence weight)
  NOCAP   registered normalisation without cap (BOMP-like unbounded energies), subspace energy
For every target event (catalogue, M >= 2.0, origin in CAL) and the window whose owned range [5, 8) s contains its
origin: s_true = best score over nodes within 3 km (epicentre) of the catalogue hypocentre and origins within 0.5 s;
s_map = window maximum; whether the window argmax lies within 5 km of the event.  Null: window maxima of the windows
whose owned range contains no catalogued origin (any magnitude) within +-2 s.
Outputs per variant: fraction of events whose s_true exceeds the 90th / 99th percentile of the null maxima, AUC of s_true
against null maxima, fraction of event windows whose argmax is within 5 km.
"""
import json
import sys
import time

import numpy as np

sys.path.insert(0, "/home/claude/work/r0")
import r0_lib as R  # noqa: E402
import r0_run as RR  # noqa: E402

import c10_lib as C  # noqa: E402
import xc  # noqa: E402

VARIANTS = ("REG", "STA", "STA-COH", "NOCAP")


def score_map(E, MFn, capv, coh=False):
    Eh = C._half_energy(MFn, capv)
    tpk, _ = C.choose(E, Eh, "hard", None, C.SH, R=MFn)
    vf = C._vectors(MFn, tpk, capv).reshape(E.N, E.T, -1)
    y = np.einsum("gtr,grk->gtk", vf, E.Q)
    s0 = (y * y).sum(-1) / 2.0
    if coh:
        tot = (vf * vf).sum(-1) / 2.0
        s0 = s0 * s0 / np.maximum(tot, 1e-9)
    return s0


def main(sub, outp):
    t0 = time.time()
    net = RR.Net(sub, 940.0)
    E = net.eng["SC"][0]
    cat = net.cat
    tg = [e for e in cat if R.T_CAL[0] <= e["t"] < R.T_CAL[1] and e["M"] >= R.M_TARGET]
    allt = np.array([e["t"] for e in cat])
    starts = RR.starts("cal")
    out = {v: {"null": [], "events": []} for v in VARIANTS}
    ownt = (E.taus >= R.OWN[0] - 1e-9) & (E.taus < R.OWN[1] - 1e-9)
    for wi, ts in enumerate(starts):
        W = R.window(net.d16, ts)
        MF, sd = xc.matched(W)
        sds = np.median(np.abs(MF), axis=(1, 2)) / 0.6745
        inp = {"REG": (MF / sd, C.CAPV, False), "STA": (MF / sds[:, None, None], C.CAPV, False),
               "STA-COH": (MF / sds[:, None, None], C.CAPV, True), "NOCAP": (MF / sd, None, False)}
        lo, hi = ts + R.OWN[0], ts + R.OWN[1]
        own = [e for e in tg if lo <= e["t"] < hi]
        quiet = not np.any((allt > lo - 2.0) & (allt < hi + 2.0))
        for v, (M_, capv, coh) in inp.items():
            s = score_map(E, M_, capv, coh)
            s_own = np.where(ownt[None, :], s, -np.inf)
            smax = float(s_own.max())
            if quiet:
                out[v]["null"].append(smax)
            gi, tj = np.unravel_index(int(np.argmax(s_own)), s.shape)
            for e in own:
                d = np.hypot(E.XN[:, 0] - e["x"], E.XN[:, 1] - e["y"])
                near = d <= 3.0
                tt = np.abs(ts + E.taus - e["t"]) <= 0.5
                st = float(s[np.ix_(near, tt)].max()) if near.any() and tt.any() else float("nan")
                out[v]["events"].append({"t": e["t"], "M": e["M"], "s_true": st, "s_map": smax,
                                         "argmax_km": float(d[gi]), "argmax_dt": float(ts + E.taus[tj] - e["t"])})
        if wi % 50 == 0:
            print(time.strftime("%H:%M:%S"), sub, wi, len(starts), flush=True)
    summ = {}
    for v in VARIANTS:
        nul = np.array(out[v]["null"])
        ev = out[v]["events"]
        st = np.array([e["s_true"] for e in ev])
        q90, q99 = np.percentile(nul, [90, 99])
        auc = float(np.mean([(np.sum(nul < x) + 0.5 * np.sum(nul == x)) / len(nul) for x in st]))
        summ[v] = {"n_events": len(ev), "n_null": int(len(nul)), "frac_above_q90": float(np.mean(st > q90)),
                   "frac_above_q99": float(np.mean(st > q99)), "auc": auc,
                   "frac_argmax_5km": float(np.mean([e["argmax_km"] <= 5.0 for e in ev])),
                   "frac_strue_is_max": float(np.mean([e["s_true"] >= 0.999 * e["s_map"] for e in ev])),
                   "median_strue_over_q90": float(np.median(st / q90))}
        print(sub, v, summ[v], flush=True)
    json.dump({"subnet": sub, "summary": summ, "detail": out, "t": time.time() - t0}, open(outp, "w"))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
