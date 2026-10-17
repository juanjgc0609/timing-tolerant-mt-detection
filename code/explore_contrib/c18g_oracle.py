"""
c18c_oracle.py - exploration C18g (Ridgecrest CAL only; first step toward the registered arrival-reading gate).  Same oracle
diagnostic as c18b_oracle.py and c18c_oracle.py, for readings that down-weight the S blocks, whose peaks C18a found
least predictable (weight w_S on the S peak vectors; w_S = 0 is a P-only score):
  CAP6 / LOG6            registered cap 6 / logarithmic compression, both phases (as C18c)
  CAP6-S05 / LOG6-S05    S blocks weighted 0.5
  CAP6-P / LOG6-P        P blocks only
    python c18g_oracle.py SUBNET OUT.json
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

VARIANTS = ("CAP6", "LOG6", "CAP6-S05", "LOG6-S05", "CAP6-P", "LOG6-P")


def compress(v, kind, c):
    if kind == "cap":
        nv = np.linalg.norm(v, axis=-1, keepdims=True)
        return v * np.minimum(1.0, c / np.maximum(nv, 1e-12))
    if kind == "log":
        e = (v * v).sum(-1, keepdims=True)
        return v * np.sqrt(c * c * np.log1p(e / (c * c)) / np.maximum(e, 1e-24))
    return v


def score_map(E, MFn, kind, c, ws=1.0):
    tpk, _ = C.choose(E, None, "hard", None, C.SH, R=MFn)
    vf = C._vectors(MFn, tpk, None)
    vf = compress(vf, kind, c)
    vf[..., 1, :, :] *= ws
    vf = vf.reshape(E.N, E.T, -1)
    y = np.einsum("gtr,grk->gtk", vf, E.Q)
    return (y * y).sum(-1) / 2.0


def main(sub, outp):
    t0 = time.time()
    net = RR.Net(sub, 940.0)
    E = net.eng["SC"][0]
    tg = [e for e in net.cat if R.T_CAL[0] <= e["t"] < R.T_CAL[1] and e["M"] >= R.M_TARGET]
    allt = np.array([e["t"] for e in net.cat])
    out = {v: {"null": [], "events": []} for v in VARIANTS}
    ownt = (E.taus >= R.OWN[0] - 1e-9) & (E.taus < R.OWN[1] - 1e-9)
    starts = RR.starts("cal")
    for wi, ts in enumerate(starts):
        W = R.window(net.d16, ts)
        MF, sd = xc.matched(W)
        inp = {"CAP6": (MF / sd, "cap", 6.0, 1.0), "LOG6": (MF / sd, "log", 6.0, 1.0),
               "CAP6-S05": (MF / sd, "cap", 6.0, 0.5), "LOG6-S05": (MF / sd, "log", 6.0, 0.5),
               "CAP6-P": (MF / sd, "cap", 6.0, 0.0), "LOG6-P": (MF / sd, "log", 6.0, 0.0)}
        lo, hi = ts + R.OWN[0], ts + R.OWN[1]
        own = [e for e in tg if lo <= e["t"] < hi]
        quiet = not np.any((allt > lo - 2.0) & (allt < hi + 2.0))
        for v, (M_, kind, c, ws) in inp.items():
            s = score_map(E, M_, kind, c, ws)
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
        if wi % 100 == 0:
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
                   "frac_strue_is_max": float(np.mean([e["s_true"] >= 0.999 * e["s_map"] for e in ev]))}
        print(sub, v, {k: (round(x, 3) if isinstance(x, float) else x) for k, x in summ[v].items()}, flush=True)
    json.dump({"subnet": sub, "summary": summ, "detail": out, "t": time.time() - t0}, open(outp, "w"))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
