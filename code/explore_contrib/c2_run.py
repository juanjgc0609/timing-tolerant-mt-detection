"""C2: dense overlaps. 3-5 events per window, epicentral separations >= 5 km (clustered within 9 km of a centre),
origin times in [5.2, 7.8] s with pairwise gaps >= 0.5 s. kmax = 7 for all methods (so false atoms remain possible)."""
import json, sys, time
import numpy as np
import xc

CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
TAG = 3002
KMAX = 7


def draw_dense(rng):
    k = int(rng.integers(3, 6))
    c = rng.uniform(-8, 8, 2)
    srcs = []
    while len(srcs) < k:
        r = 9.0 * np.sqrt(rng.uniform())
        a = rng.uniform(0, 2 * np.pi)
        xy = c + r * np.array([np.cos(a), np.sin(a)])
        if np.any(np.abs(xy) > 15):
            continue
        if all(np.linalg.norm(xy - q[:2]) >= 5.0 for q in srcs):
            srcs.append(np.array([xy[0], xy[1], rng.uniform(5, 13)]))
    slack = 2.6 - 0.5 * (k - 1)
    u = np.sort(rng.uniform(0, slack, k))
    taus = (5.2 + u + 0.5 * np.arange(k))
    perm = rng.permutation(k)
    return srcs, [float(taus[p]) for p in perm]


def run(regime, nets, nwin, out):
    res = []
    for cell in CELLS:
        for j in nets:
            nc = xc.NetCfg(cell, j, regime)
            st = xc.window_stream(nc, TAG, draw=draw_dense)
            ops = xc.Ops(nc.STA, nc.corr, xc.make_grid())
            for wi in range(nwin):
                w = next(st)
                tr = w["truth"]
                osc = xc.oracle_scores(ops, w)
                row = {"cell": cell, "net": j, "w": wi, "k": len(tr), "truth": [t.tolist() for t in tr],
                       "taus": w["taus"], "oracle": osc, "rec": {}, "loc": {}, "sec": {}}
                t = time.time(); r, l = xc.bomp_pk(ops.PD, w["MFn"], tr, kmax=KMAX); row["sec"]["BOMP-PK"] = time.time() - t
                row["rec"]["BOMP-PK"], row["loc"]["BOMP-PK"] = r, l
                A, Bc = xc.ceilings(r, osc, tr, l)
                row["rec"]["CEIL-A"], row["rec"]["CEIL-B"] = A, Bc
                t = time.time(); r, l, _ = xc.tssfw(ops, w, tr, f=xc.F_FROZEN, kmax=KMAX); row["sec"]["TS-SFW"] = time.time() - t
                row["rec"]["TS-SFW"], row["loc"]["TS-SFW"] = r, l
                t = time.time(); r, l, _ = xc.tssfw(ops, w, tr, lam=0.0, kmax=KMAX); row["sec"]["SLIDE"] = time.time() - t
                row["rec"]["SLIDE"], row["loc"]["SLIDE"] = r, l
                t = time.time()
                r = [(float(s), int(h)) for s, h in xc.certified_fw_record(ops.STA, w["MFn"], ops.grid, xc.DEFAULT_TAUS, ops.Gs, ops.Qs, tr, kmax=KMAX, tcorr=ops.corr)]
                row["sec"]["MT-FW"] = time.time() - t
                row["rec"]["MT-FW"] = r
                t = time.time(); r, l = xc.bomp_pk_mask(ops.PD, w["MFn"], tr, kmax=KMAX, mw=6, amp_gate=3.0)
                row["sec"]["MASK6g3"] = time.time() - t
                row["rec"]["MASK6g3"], row["loc"]["MASK6g3"] = r, l
                res.append(row)
            print(cell, j, "done", time.strftime("%H:%M:%S"), flush=True)
            with open(out, "w") as fh:
                json.dump(res, fh)


if __name__ == "__main__":
    run(sys.argv[1], [int(x) for x in sys.argv[2].split(",")], int(sys.argv[3]), sys.argv[4])
