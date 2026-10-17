"""
c18d_detect.py - exploration C18d (Ridgecrest CAL only, origins in [0, 900) s).  Full MT-FW loop (registered R0 M1
path: +SC terms, hard +-0.25 s read, kmax 6, coarse grid, 40 s windows every 3 s owning origins in [5, 8) s) with the
per-arrival statistic changed: CAP6 (= R0 M1), CAP24, LOG6 (energy c^2 log(1 + E / c^2), c = 6), NOCAP; and BOMP-PK
(= R0 M3) for reference.  Writes the owned atoms; c18d_analyze.py computes recall at 30 false detections per hour on CAL.
    python c18d_detect.py SUBNET VARIANTS OUT.json
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

_vectors_orig = C._vectors


def _vectors_var(Rr, tpk, capv):
    """capv > 0: cap (registered); capv None: no cap; capv < 0: logarithmic compression with c = -capv."""
    if capv is not None and capv < 0:
        v = _vectors_orig(Rr, tpk, None)
        c = -capv
        e = (v * v).sum(-1, keepdims=True)
        return v * np.sqrt(c * c * np.log1p(e / (c * c)) / np.maximum(e, 1e-24))
    return _vectors_orig(Rr, tpk, capv)


C._vectors = _vectors_var
CAPS = {"CAP6": 6.0, "CAP24": 24.0, "LOG6": -6.0, "NOCAP": None}


def main(sub, variants, outp):
    net = RR.Net(sub, 940.0)
    E = net.eng["SC"][0]
    res = {"subnet": sub, "atoms": {v: [] for v in variants}, "sec": {v: 0.0 for v in variants}}
    ws = RR.starts("cal")
    for wi, ts in enumerate(ws):
        MF, sd = xc.matched(R.window(net.d16, ts))
        MFn = MF / sd
        for v in variants:
            t0 = time.time()
            if v == "BOMP":
                at = C.bomp_loop(E, None, MFn, R.KMAX)
            else:
                at = C.mt_loop(E, MFn, R.KMAX, capv=CAPS[v])
            for a in at:
                if R.OWN[0] - 1e-9 <= a["tau"] < R.OWN[1] - 1e-9:
                    x = a["x"]
                    res["atoms"][v].append([float(x[0]), float(x[1]), float(x[2]), float(ts + a["tau"]), float(a["s"]), wi])
            res["sec"][v] += time.time() - t0
        if wi % 50 == 0:
            print(time.strftime("%H:%M:%S"), sub, wi, len(ws), {k: round(x, 1) for k, x in res["sec"].items()}, flush=True)
    json.dump(res, open(outp, "w"))
    print("DONE", sub, flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2].split(","), sys.argv[3])
