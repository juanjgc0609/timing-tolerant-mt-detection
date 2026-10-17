"""
r0_run_chunks.py - resumable driver for the frozen R0 detectors (closure of R0, post hoc).  It imports r0_run.py
unchanged and runs exactly its per-window loop (same Net, same windows, same methods, same owned-origin rule) on a range
of window indices, writing one chunk file per range atomically, so that a container restart loses at most one chunk.
    python r0_run_chunks.py SEGMENT SUBNET METHODS OUTDIR [CHUNK]          (CHUNK windows per file, default 50)
    python r0_run_chunks.py merge SEGMENT SUBNET METHODS OUTDIR OUT.json    joins the chunks in window order
"""
import glob
import json
import os
import sys
import time

import r0_lib as R
import r0_run as RR
import xc  # noqa: E402


def run(segment, subnet, methods, outdir, chunk=50):
    os.makedirs(outdir, exist_ok=True)
    ws = RR.starts(segment)
    todo = [(a, min(a + chunk, len(ws))) for a in range(0, len(ws), chunk)]
    todo = [(a, b) for a, b in todo if not os.path.exists(os.path.join(outdir, f"chunk_{a:04d}_{b:04d}.json"))]
    if not todo:
        return
    t_end = (R.T_CAL[1] if segment == "cal" else R.T_EVAL[1]) + 40.0
    net = RR.Net(subnet, t_end)
    for a, b in todo:
        res = {"atoms": {m: [] for m in methods}, "sec": {m: 0.0 for m in methods}, "range": [a, b]}
        for wi in range(a, b):
            ts = ws[wi]
            W = R.window(net.d16, ts)
            MF, sd = xc.matched(W)
            MFn = MF / sd
            for m in methods:
                t0 = time.time()
                for x, tau, s in net.run(m, MFn):
                    if R.OWN[0] - 1e-9 <= tau < R.OWN[1] - 1e-9:
                        res["atoms"][m].append([float(x[0]), float(x[1]), float(x[2]), float(ts + tau), float(s), wi])
                res["sec"][m] += time.time() - t0
        fn = os.path.join(outdir, f"chunk_{a:04d}_{b:04d}.json")
        json.dump(res, open(fn + ".tmp", "w"))
        os.replace(fn + ".tmp", fn)
        print(time.strftime("%H:%M:%S"), subnet, segment, methods, "chunk", a, b, "of", len(ws), flush=True)
    meta = {"segment": segment, "subnet": subnet, "names": net.names, "centre": net.centre, "STA": net.STA.tolist(),
            "scale": net.scale.tolist(), "corr_sc": net.corr_sc.tolist(), "corr_s2": net.corr_s2.tolist(),
            "n_calwin": len(net.calwin), "t_setup": net.t_setup, "n_windows": len(ws)}
    json.dump(meta, open(os.path.join(outdir, "meta.json"), "w"))


def merge(segment, subnet, methods, outdir, outp):
    ws = RR.starts(segment)
    files = sorted(glob.glob(os.path.join(outdir, "chunk_*.json")))
    rng = [json.load(open(f))["range"] for f in files]
    covered = sorted(i for a, b in rng for i in range(a, b))
    assert covered == list(range(len(ws))), "chunks do not cover every window exactly once"
    out = json.load(open(os.path.join(outdir, "meta.json")))
    out.update({"atoms": {m: [] for m in methods}, "sec": {m: 0.0 for m in methods}})
    for f in files:
        r = json.load(open(f))
        for m in methods:
            out["atoms"][m].extend(r["atoms"][m])
            out["sec"][m] += r["sec"][m]
    json.dump(out, open(outp + ".tmp", "w"))
    os.replace(outp + ".tmp", outp)
    print("MERGED", outp, {m: len(v) for m, v in out["atoms"].items()}, flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "merge":
        merge(sys.argv[2], sys.argv[3], sys.argv[4].split(","), sys.argv[5], sys.argv[6])
    else:
        run(sys.argv[1], sys.argv[2], sys.argv[3].split(","), sys.argv[4], int(sys.argv[5]) if len(sys.argv) > 5 else 50)
