"""FE-FIX (card.md): fixing the model's reading before pilot P1.

Stage `sp` (step 1, P/S overlap, descriptive, no new waveform run): for each noise-free analytic-G event of
N1-STACK gcheck2.json, the S-P time of every station from the set's prediction model at the event's true location
(the same first_arrival call gcheck used), the minimum S-P of the event, and the median noise-free f after removing
the blocks with S-P < 0.6 s. When no block of an event is below 0.6 s nothing is removed and its f is unchanged by
construction; when some are, the stage stops, because gcheck2.json holds no per-block reads and that case would need
a new run (laptop).
"""
import argparse, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
SP_CUT = 0.6          # s, card step 1
WIN_END = 0.75        # s, descriptive: end of the read window (+8 samples at 16 Hz) plus the largest common shift (+4)


def sp_stage(a):
    from ttmt.synth import realdata as RD
    from ttmt.synth.realdata import first_arrival
    rows = json.load(open(a.gcheck))["rows"]
    out = []
    for sid in sorted({r["set_id"] for r in rows}):
        man = json.loads((ROOT / "data" / "manifests" / "synth" / f"{sid}.json").read_text())
        STA = np.asarray(man["xyz"], float)
        model = RD.pred_model(man["dataset"])
        R = [r for r in rows if r["set_id"] == sid]
        for r in R:
            x = np.asarray(r["xyz"], float)
            X = np.hypot(STA[:, 0] - x[0], STA[:, 1] - x[1])
            z = max(float(x[2]), 0.5)
            sp = first_arrival(model, z, X, "S")[0] - first_arrival(model, z, X, "P")[0]
            n_cut = int((sp < SP_CUT).sum())
            if n_cut:
                sys.exit(f"{sid}: {n_cut} blocks with S-P < {SP_CUT} s; removal needs per-block reads (new run): stop")
            out.append({"set_id": sid, "f_clean": r["f_clean"], "sp_min": float(sp.min()), "n_sta": len(sp),
                        "n_below_cut": n_cut, "n_below_window": int((sp < WIN_END).sum()),
                        "f_clean_without": r["f_clean"]})
        fc = np.array([o["f_clean"] for o in out if o["set_id"] == sid])
        fw = np.array([o["f_clean_without"] for o in out if o["set_id"] == sid])
        sm = np.array([o["sp_min"] for o in out if o["set_id"] == sid])
        nb = sum(o["n_below_cut"] for o in out if o["set_id"] == sid)
        nw = sum(o["n_below_window"] for o in out if o["set_id"] == sid)
        rise = float(np.median(fw) - np.median(fc))
        zone = ("P/S overlap is the main cause" if np.median(fw) >= 0.90 else
                ("P/S overlap is not the cause" if rise < 0.05 else "P/S overlap is part of the cause"))
        rho = float(np.corrcoef(np.argsort(np.argsort(sm)), np.argsort(np.argsort(fc)))[0, 1])
        print(f"{sid}: {len(fc)} events; S-P minimum per event median {np.median(sm):.2f} s [min {sm.min():.2f}, "
              f"max {sm.max():.2f}]; blocks with S-P < {SP_CUT} s: {nb}; noise-free f median {np.median(fc):.3f} -> "
              f"{np.median(fw):.3f} without them (rise {rise:+.3f}) -> {zone} (card rule: >= 0.90 main / < +0.05 not)")
        print(f"  descriptive: blocks with S-P < {WIN_END} s (read window end + largest shift): {nw}; "
              f"Spearman(S-P minimum, noise-free f) {rho:+.2f}")
    json.dump({"rows": out}, open(a.out, "w"), indent=1)


# ============================================================================================ step 2: chain
N1 = None


def _n1():
    global N1
    if N1 is None:
        import importlib.util
        spec = importlib.util.spec_from_file_location("n1", ROOT / "experiments/exploratory/N1-STACK/run.py")
        N1 = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(N1)
    return N1


FAKE_SEED, FAKE_NS = 0, 7       # the cloud check of gcheck (N1-STACK, a7cea2c): 7 stations uniform in +-20 km, z 0
REAL = {"geometry": "real", "model": "real", "sd": "real", "resample": True, "band": (0.5, 6.0), "ricker": True,
        "zflip": True, "window": (-1, 8), "peak": "max", "shift": 4, "qweight": False}
CLOUD = dict(REAL, geometry="fake", model="ridgecrest_2019", sd="one")
# one element at a time; the alternative of every element is fixed here, before running
TO_REAL = [("geometry and stations of the set", {"geometry": "real"}),
           ("velocity model of the set", {"model": "real"}),
           ("per-station noise normalisation (robust sd of the set's noise block)", {"sd": "real"})]
FROM_REAL = [("no per-station noise normalisation (sd = 1)", {"sd": "one"}),
             ("no resampling (read at 100 Hz, same window and shifts in seconds)", {"resample": False}),
             ("no 0.5-6 Hz band-pass (demean and detrend only)", {"band": None}),
             ("no 2 Hz Ricker matched filter", {"ricker": False}),
             ("no Z flip (sign check: must not raise f)", {"zflip": False}),
             ("wider read window [-0.25, +1.0] s", {"window": (-4, 16)}),
             ("peak rule: principal polarisation of the window instead of the energy peak", {"peak": "pca"}),
             ("no common shift (read at the predicted onset only)", {"shift": 0})]


def _front(x, cfg, fs=100.0):
    """the R1 front end (N1-STACK pre16 + matched) with each element switchable; returns (R, fs of R)."""
    from scipy.signal import butter, resample_poly, sosfiltfilt
    from ttmt.kernel import ricker
    x = np.array(x, float, copy=True)
    x -= x.mean(-1, keepdims=True)
    tt = np.arange(x.shape[-1], dtype=float); tc = tt - tt.mean()
    x -= (x * tc).sum(-1, keepdims=True) / float((tc * tc).sum()) * tc
    if cfg["band"] is not None:
        x = sosfiltfilt(butter(4, list(cfg["band"]), btype="band", fs=fs, output="sos"), x, axis=-1)
    fr = fs
    if cfg["resample"]:
        x = resample_poly(x, 4, 25, axis=-1); fr = 16.0
    if cfg["zflip"]:
        x[:, 2] *= -1.0
    if cfg["ricker"]:
        taps = ricker(np.arange(-1.0, 1.0 + 1e-9, 1.0 / fr), 2.0)[::-1]
        x = np.array([[np.convolve(x[i, c], taps, mode="same") for c in range(3)] for i in range(x.shape[0])])
    return x, fr


def _read(R, fr, arr_s, t0, cfg):
    """(2, ns, 3) read: per block, within the window (16 Hz samples, scaled to fr) after the predicted onset plus a
    common shift; energy peak (N1 rule) or principal polarisation (sign of the energy-peak sample, norm sqrt of the
    largest eigenvalue of the window's 3x3 second moment); the common shift with the largest ||v||^2."""
    ns, _, n = R.shape
    k16 = fr / 16.0
    lo, hi = int(round(cfg["window"][0] * k16)), int(round(cfg["window"][1] * k16))
    offs = np.arange(lo, hi + 1)
    k = np.floor((arr_s.T - t0) * fr + 0.5)                                # (2, ns)
    best, out = -1.0, None
    for d in range(-cfg["shift"], cfg["shift"] + 1):
        v = np.zeros((2, ns, 3))
        for b in range(2):
            for i in range(ns):
                if not np.isfinite(k[b, i]):
                    continue
                ii = (k[b, i] + int(round(d * k16)) + offs).astype(int)
                ii = ii[(ii >= 0) & (ii < n)]
                if not len(ii):
                    continue
                W = R[i][:, ii]                                            # (3, w)
                j = int((W * W).sum(0).argmax())
                if cfg["peak"] == "max":
                    v[b, i] = W[:, j]
                else:
                    lam, U = np.linalg.eigh(W @ W.T)
                    u = U[:, -1] * np.sign(U[:, -1] @ W[:, j] or 1.0)
                    v[b, i] = u * np.sqrt(lam[-1])
        e = float((v * v).sum())
        if e > best:
            best, out = e, v
    return out


def chain_f(cfg, sid, rows_set, rng, man, sd_real):
    """noise-free f = ||Q^T v||^2 / ||v||^2 of the gcheck events of one set under one chain configuration. The moment
    tensors and origins are redrawn with gcheck's seed and order (checked against gcheck2.json); locations from
    gcheck2.json."""
    from ttmt.kernel import ricker
    from ttmt.synth import inject as INJ
    from ttmt.synth import realdata as RD
    from ttmt.synth.realdata import first_arrival
    from ttmt.radiation import radiation_operator
    if cfg["geometry"] == "real":
        STA = np.asarray(man["xyz"], float)
    else:
        g0 = np.random.default_rng(FAKE_SEED)
        STA = np.c_[g0.uniform(-20, 20, (FAKE_NS, 2)), np.zeros(FAKE_NS)]
    ns = len(STA)
    model = RD.pred_model(man["dataset"] if cfg["model"] == "real" else cfg["model"])
    sd = sd_real if cfg["sd"] == "real" else np.ones(ns)
    if sd is None:
        return None
    fs, n = float(man["fs"]), 30000
    fvals = []
    for r in rows_set:
        x = np.asarray(r["xyz"], float)
        m = rng.normal(size=6)
        t0 = float(rng.uniform(20.0, 300.0 - 40.0))
        if abs(t0 - r["t0"]) > 1e-9:
            sys.exit(f"{sid}: redrawn origin {t0} differs from gcheck2.json {r['t0']}: stop")
        X = np.hypot(STA[:, 0] - x[0], STA[:, 1] - x[1]); z = max(float(x[2]), 0.5)
        tt = np.stack([first_arrival(model, z, X, "P")[0], first_arrival(model, z, X, "S")[0]], 1)
        g = (radiation_operator(x, STA) @ m).reshape(2, ns, 3)
        pre, L = 2.0, 40.0
        tl = np.arange(int(L * fs)) / fs - pre
        w = np.zeros((ns, 3, len(tl)))
        for b in range(2):
            for i in range(ns):
                w[i] += g[b, i][:, None] * ricker(tl - tt[i, b] - 0.1, 4.0)[None]
        w[:, 2] *= -1.0
        w *= INJ._taper(w.shape[-1], fs, 0.5)
        i0 = int(np.floor((t0 - pre) * fs)); frac = (t0 - pre) * fs - i0
        w = INJ._shift(w, fs, frac / fs)
        unit = np.zeros((ns, 3, n)); unit[:, :, i0:i0 + w.shape[-1]] = w
        R, fr = _front(unit, cfg, fs)
        R = R / sd[:, None, None]
        v = _read(R, fr, t0 + tt, 0.0, cfg).reshape(-1)
        if cfg["qweight"]:                                                 # Q of diag(1/sd) G (supervisor, 17 Nov)
            from ttmt.radiation import basis
            wt = np.tile(np.repeat(1.0 / sd, 3), 2)
            Q = basis(wt[:, None] * radiation_operator(x, STA))
        else:
            Q = _n1().operator_basis(x, STA)
        y = v @ Q
        fvals.append(float(y @ y) / max(float(v @ v), 1e-300))
    return fvals


def chain_stage(a):
    rows = json.load(open(a.gcheck))["rows"]
    sids = sorted({r["set_id"] for r in rows})
    sdr = {}
    if a.laptop:
        from ttmt.synth import sets
        for sid in sids:
            rd = sets.SetReader(sid, "tune", verify="fingerprint")
            _s, noise = rd.arrays(rd.labels[0])
            sdr[sid] = _n1().station_sd(_n1().matched(_n1().pre16(noise)))
    configs = [("cloud check (start)", CLOUD)]
    acc = dict(CLOUD)
    for name, ch in TO_REAL:
        configs.append((f"cloud + {name}", dict(CLOUD, **ch)))
        acc.update(ch)
        configs.append((f"  cumulative up to: {name}", dict(acc)))
    configs.append(("real chain (gcheck) = (a): sd normalisation, Q of unweighted G", REAL))
    configs.append(("(b) real chain with Q of diag(1/sd) G", dict(REAL, qweight=True)))
    for name, ch in FROM_REAL:
        configs.append((f"real - {name}", dict(REAL, **ch)))
    out = []
    for sid in sids:
        man = json.loads((ROOT / "data" / "manifests" / "synth" / f"{sid}.json").read_text())
        R = [r for r in rows if r["set_id"] == sid]
        print(f"{sid}: {len(R)} analytic events, noise-free f median [min, max] per configuration")
        for name, cfg in configs:
            rng = np.random.default_rng(20261125)
            for s2 in sids[:sids.index(sid)]:          # gcheck draws all sets from one generator, in sorted order
                for _ in [r for r in rows if r["set_id"] == s2]:
                    rng.normal(size=6); rng.uniform(20.0, 260.0)
            f = chain_f(cfg, sid, R, rng, man, sdr.get(sid))
            if f is None:
                print(f"  {name:<95s} needs the set's noise (laptop)")
                continue
            f = np.array(f)
            print(f"  {name:<95s} {np.median(f):.3f} [{f.min():.3f}, {f.max():.3f}]")
            out.append({"set_id": sid, "config": name, "cfg": {k: v for k, v in cfg.items()}, "f": f.tolist()})
        if sid in sdr:
            s = sdr[sid]
            print(f"  descriptive: per-station robust sd, max/min {s.max() / s.min():.2f}")
    json.dump({"rows": out, "sd": {k: v.tolist() for k, v in sdr.items()}}, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="stage", required=True)
    s = sub.add_parser("sp"); s.add_argument("gcheck"); s.add_argument("--out", default="sp.json")
    c = sub.add_parser("chain"); c.add_argument("gcheck"); c.add_argument("--out", default="chain.json")
    c.add_argument("--laptop", action="store_true", help="read each set's noise block for the per-station sd")
    a = p.parse_args()
    {"sp": sp_stage, "chain": chain_stage}[a.stage](a)
