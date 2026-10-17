"""
r1_run.py - Task R1, P1 (registered): MT-FW without cap on the layered JMA2001 engine (c18j_layered: first-arrival
tables, matched-filter station terms over the CAL truth events with magnitude >= 1.0, +-2 s search), 2 Hz matched
filter, R0 windows (40 s every 3 s, owned origins [5, 8) s), kmax 6.  Resumable: one file per chunk of windows.
    python r1_run.py run SEGMENT PICKS.json OUTDIR [CHUNK]
    python r1_run.py merge SEGMENT OUTDIR OUT.json
"""
import glob
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r1_lib as R  # noqa: E402

R.setup()
import c10_lib as C  # noqa: E402
import c18j_layered as LJ  # noqa: E402
import r0_run as RR  # noqa: E402
import xc  # noqa: E402

import numpy as np  # noqa: E402


def mf_terms(net, pick_terms, cap=2.0, frac=0.5, gate=25.0, nmin=3):
    """matched-filter station terms anchored on the pick terms (amendment of R1, alignment): the energy-peak readings of
    the calibration events are searched within +-min(2 s, half the predicted S-P time) of the homogeneous prediction
    (a P reading must not reach the S arrival), converted to delays against the JMA2001 times; per phase, the lag of
    the energy peak behind the onset is the median over station-phases with >= nmin gated readings of (median reading
    - pick term); term = pick term + lag of its phase."""
    import c18_tt as TT
    from station_corrections import peak_delays
    from task0_common import VP, VS
    ds, ens = [], []
    for D, src, tau in net.calwin:
        src = np.asarray(src, float)
        tl = TT.times(LJ.MODEL, LJ.R_VPVS, src, net.STA)
        d = np.linalg.norm(net.STA - src, axis=1)
        th = np.stack([d / VP, d / VS], axis=1)
        dl, en = np.zeros((len(d), 2)), np.zeros((len(d), 2))
        for i in range(len(d)):
            a, b = peak_delays(net.STA, D, src, float(tau), min(cap, frac * float(tl[1, i] - tl[0, i])))
            dl[i], en[i] = a[i], b[i]
        ds.append(dl + th - tl.T)
        ens.append(en)
    ds, ens = np.array(ds), np.array(ens)
    pick_terms = np.asarray(pick_terms, float)
    lag, info = np.zeros(2), {}
    for b in range(2):
        v = [float(np.median(ds[ens[:, i, b] >= gate, i, b])) - pick_terms[i, b] for i in range(ds.shape[1])
             if (ens[:, i, b] >= gate).sum() >= nmin]
        lag[b] = float(np.median(v)) if v else 0.0
        info["PS"[b]] = v
    return pick_terms + lag[None, :], lag, info


def engine(net, pick_terms):
    """layered JMA2001 engine (c18j_layered tables and timing) with the anchored matched-filter terms."""
    E = net.eng["SC"][0]
    tab = LJ.tt_table(sorted({float(z) for z in E.XN[:, 2]}))
    corr, lag, info = mf_terms(net, pick_terms)
    tt = LJ.tt_nodes(E.XN, E.STA, tab)
    x = (tt[:, None] + E.taus[None, :, None, None] + corr.T[None, None]) / R.DT
    E.IDX = np.clip(x.astype(int), 0, R.NT - 1)
    E.frac = x - E.IDX
    E.cache = {}
    return E, corr, lag, info


def run(segment, pkp, outdir, chunk=40):
    os.makedirs(outdir, exist_ok=True)
    ws = RR.starts(segment)
    todo = [(a, min(a + chunk, len(ws))) for a in range(0, len(ws), chunk)]
    todo = [(a, b) for a, b in todo if not os.path.exists(os.path.join(outdir, f"chunk_{a:04d}_{b:04d}.json"))]
    if not todo:
        return
    net = R.make_net(segment)
    E, corr, lag, info = engine(net, json.load(open(pkp))["pick_terms"])
    meta = {"lag": lag.tolist(), "lag_readings": info, "segment": segment, "names": net.names, "centre": net.centre, "STA": net.STA.tolist(),
            "scale": net.scale.tolist(), "corr_layered": corr.tolist(), "n_calwin": len(net.calwin),
            "n_windows": len(ws)}
    json.dump(meta, open(os.path.join(outdir, "meta.json"), "w"))
    print("corrections P", [round(float(c), 2) for c in corr[:, 0]], "S", [round(float(c), 2) for c in corr[:, 1]],
          "calwin", len(net.calwin), flush=True)
    for a, b in todo:
        res = {"atoms": [], "sec": 0.0, "range": [a, b]}
        t0 = time.time()
        for wi in range(a, b):
            ts = ws[wi]
            MF, sd = xc.matched(R.window(net.d16, ts))
            for at in C.mt_loop(E, MF / sd, R.KMAX, capv=None):
                if R.OWN[0] - 1e-9 <= at["tau"] < R.OWN[1] - 1e-9:
                    x = at["x"]
                    res["atoms"].append([float(x[0]), float(x[1]), float(x[2]), float(ts + at["tau"]), float(at["s"]), wi])
        res["sec"] = time.time() - t0
        fn = os.path.join(outdir, f"chunk_{a:04d}_{b:04d}.json")
        json.dump(res, open(fn + ".tmp", "w"))
        os.replace(fn + ".tmp", fn)
        print(time.strftime("%H:%M:%S"), segment, "chunk", a, b, "of", len(ws), round(res["sec"]), "s", flush=True)


def merge(segment, outdir, outp):
    ws = RR.starts(segment)
    files = sorted(glob.glob(os.path.join(outdir, "chunk_*.json")))
    rng = [json.load(open(f))["range"] for f in files]
    assert sorted(i for a, b in rng for i in range(a, b)) == list(range(len(ws))), "chunks do not cover the windows"
    out = json.load(open(os.path.join(outdir, "meta.json")))
    out.update({"atoms": {"NOCAP": []}, "sec": 0.0})
    for f in files:
        r = json.load(open(f))
        out["atoms"]["NOCAP"].extend(r["atoms"])
        out["sec"] += r["sec"]
    json.dump(out, open(outp + ".tmp", "w"))
    os.replace(outp + ".tmp", outp)
    print("MERGED", outp, len(out["atoms"]["NOCAP"]), flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]) if len(sys.argv) > 5 else 40)
    else:
        merge(sys.argv[2], sys.argv[3], sys.argv[4])
