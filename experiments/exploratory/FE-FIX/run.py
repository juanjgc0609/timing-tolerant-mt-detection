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
    if len(sd) != ns:          # fake geometry (7 stations) with a set's sd: the set's first 7 stations' sd, in order
        sd = np.asarray(sd)[:ns]
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


# ============================================================================================ step 3: band
BANDS = [(0.5, 6.0), (1.0, 10.0), (2.0, 20.0), (2.0, 30.0)]
SNR_MIN, SEED_SPLIT, SEED_BOOT3, B3 = 5.0, 20261123, 20261124, 10000


def _bandfront(x, band):
    """candidate front end: demean, detrend, Butterworth 4 zero-phase band-pass at 100 Hz, Z down; no resampling and
    no matched filter (a 16 Hz rate cannot hold bands above 8 Hz; the Ricker of 2 Hz belongs to the 0.5-6 Hz band)."""
    return _front(x, dict(REAL, band=band, resample=False, ricker=False), 100.0)[0]


def _split(labels, sid):
    """half for choosing, half for confirming, by realisation, per set (stratified), fixed seed."""
    rng = np.random.default_rng([SEED_SPLIT, sum(map(ord, sid))])
    p = list(rng.permutation(sorted(labels)))
    h = (len(p) + 1) // 2
    return set(p[:h]), set(p[h:])


QNOISE, GAIN_PASS = 0.999, 0.20
CAP = 6.0


def _peaks(R, lo, hi):
    """(ns, n, 3): for every station and centre sample c, the 3-vector at the energy peak in c + [lo, hi]."""
    ns, _, n = R.shape
    E = (R * R).sum(1)
    best = np.full((ns, n), -1.0); k = np.zeros((ns, n), int)
    for o in range(lo, hi + 1):
        j = np.clip(np.arange(n) + o, 0, n - 1)
        upd = E[:, j] > best
        best = np.where(upd, E[:, j], best); k = np.where(upd, j[None, :], k)
    return R.transpose(0, 2, 1)[np.arange(ns)[:, None], k]


def _smt(Vp, fr, origins, off_s, Qw, cap):
    """S_mt = 1/2 ||Qw^T phi(v)||^2 for each origin sample (array), arrivals at origin + off_s (2, ns) seconds; the read
    of N1-STACK: common shift -0.25..+0.25 s in steps of 1/16 s, the one with the largest ||v||^2; phi = CAP6 or none."""
    ns, n, _ = Vp.shape
    off = off_s * fr                                                          # (2, ns) samples
    bestE = np.full(len(origins), -1.0); bestv = np.zeros((len(origins), 2, ns, 3))
    for d in range(-4, 5):
        c = np.floor(origins[:, None, None] + off[None] + d * fr / 16.0 + 0.5)
        ok = np.isfinite(c) & (c >= 0) & (c < n)
        cc = np.where(ok, c, 0).astype(int)
        v = np.where(ok[..., None], Vp[np.arange(ns)[None, None, :], cc], 0.0)
        e = (v * v).sum((1, 2, 3))
        upd = e > bestE
        bestE = np.where(upd, e, bestE); bestv[upd] = v[upd]
    if cap:
        nv = np.linalg.norm(bestv, axis=-1, keepdims=True)
        bestv = bestv * np.minimum(1.0, CAP / np.maximum(nv, 1e-12))
    y = bestv.reshape(len(origins), -1) @ Qw
    return 0.5 * (y * y).sum(1)


def band_stage(a):
    """card.md, step 3 as revised by the supervisor on 17 November (before running). Population: every target of the
    tune partition with snr_med >= 5. For each candidate (four bands at 100 Hz without resampling or Ricker, and the
    operative reference 0.5-6 Hz at 16 Hz with the Ricker), every station divided by the robust sd of its noise
    block in that front end, catalogue mask of N1-STACK, Q of diag(1/sd) G at the true location (the step-2 fix):
    margin m = S_mt of the target (true location, observed arrivals) / the 0.999 quantile of the same S_mt over every
    origin of the unslid masked noise block at the same location (predicted arrivals of the set's model). CAP6 is the
    main statistic, uncapped descriptive. f_emp with g/sd as a reading check of the chosen band."""
    from ttmt.synth import realdata as RD
    from ttmt.synth import sets
    from ttmt.synth.realdata import first_arrival
    from ttmt.radiation import basis, radiation_operator
    n1 = _n1()
    cands = [f"{b[0]:g}-{b[1]:g}" for b in BANDS]
    REF = "R1 0.5-6 16Hz Ricker"
    out = []
    for sid in a.sets:
        rd = sets.SetReader(sid, "tune", verify="fingerprint")
        STA = np.asarray(rd.man["xyz"], float)
        model = RD.pred_model(rd.man["dataset"])
        tr0 = rd.truth(rd.labels[0])
        _s, noise = rd.arrays(rd.labels[0])
        mfn = n1.matched(n1.pre16(noise))
        m16 = n1.catalogue_mask(tr0.t0, mfn.shape[-1], tr0, STA, model)
        idx = np.minimum((np.arange(noise.shape[-1]) * 16 / 100).astype(int), m16.shape[1] - 1)
        k16, k100 = (~m16)[:, None, :].astype(float), (~m16[:, idx])[:, None, :].astype(float)
        fe = {}
        for b, name in zip(BANDS, cands):
            rn = _bandfront(noise, b)
            fe[name] = (lambda x, b=b: _bandfront(x, b), 100.0, np.median(np.abs(rn), axis=(1, 2)) / 0.6745, k100, rn)
        fe[REF] = (lambda x: n1.matched(n1.pre16(x)), 16.0, n1.station_sd(mfn), k16, mfn)
        Vn = {}
        for name, (fn, fr, sd, kk, rn) in fe.items():
            k = fr / 16.0
            Vn[name] = _peaks(rn * kk / sd[:, None, None], int(round(-1 * k)), int(round(8 * k)))
        choose, confirm = _split(rd.labels, sid)
        for lab in rd.labels:
            sig, _n = rd.arrays(lab)
            truth = rd.truth(lab)
            n1.check_station_order(rd, truth)
            snr_med = np.nanmedian(truth.snr, axis=1)
            tg = [e for e in range(len(truth.t)) if truth.target[e] and np.isfinite(snr_med[e]) and snr_med[e] >= SNR_MIN]
            if not tg:
                continue
            Vs, Rs = {}, {}
            for name, (fn, fr, sd, kk, rn) in fe.items():
                k = fr / 16.0
                Rs[name] = (fn(sig) + rn * kk) / sd[:, None, None]
                Vs[name] = _peaks(Rs[name], int(round(-1 * k)), int(round(8 * k)))
            for e in tg:
                x = truth.xyz[e]
                G = radiation_operator(x, STA)
                X = np.hypot(STA[:, 0] - x[0], STA[:, 1] - x[1]); z = max(float(x[2]), 0.5)
                tt = np.stack([first_arrival(model, z, X, "P")[0], first_arrival(model, z, X, "S")[0]], 0)  # (2, ns)
                g = np.array(truth.meta["targets"][e]["g"], float)
                g = np.nan_to_num(g.transpose(1, 0, 2)); g[~np.isfinite(truth.arr_t_obs[e].T)] = 0.0; g[..., 2] *= -1.0
                row = {"set_id": sid, "label": lab, "event": int(e), "half": "choose" if lab in choose else "confirm",
                       "snr_med": float(snr_med[e]), "pair": int(truth.pair[e])}
                for name, (fn, fr, sd, kk, rn) in fe.items():
                    w = np.tile(np.repeat(1.0 / sd, 3), 2)
                    Qw = basis(w[:, None] * G)
                    off_t = (truth.arr_t_obs[e].T - truth.t0)                  # (2, ns) s from the record start
                    o_n = np.arange(0, Vn[name].shape[1], fr / 16.0)
                    o_t = np.array([(truth.t[e] - truth.t0) * fr])                 # true origin, record samples
                    for cap, tag in ((True, "m"), (False, "m_nocap")):
                        sn = max(float(np.quantile(_smt(Vn[name], fr, o_n, tt, Qw, cap), QNOISE)), 1e-300)
                        # main: predicted arrivals at the true location and origin (no oracle times; supervisor)
                        row[f"{tag}|{name}"] = float(_smt(Vs[name], fr, o_t, tt, Qw, cap)[0]) / sn
                        # descriptive: observed arrivals (cost of the timing error in each band)
                        row[f"{tag}_obs|{name}"] = float(_smt(Vs[name], fr, np.array([0.0]), off_t, Qw, cap)[0]) / sn
                    # reading check: f_emp with g/sd, the read at the best common shift of the target
                    gw = g.reshape(-1) * w
                    Rt = Rs[name]
                    vv = (n1.read_at_best(Rt, truth.arr_t_obs[e], truth.t0) if fr == 16.0 else
                          _read(Rt, 100.0, truth.arr_t_obs[e], truth.t0, REAL))
                    vv = vv.reshape(-1)
                    row[f"f_emp|{name}"] = float((vv @ gw) ** 2 / (max(float(vv @ vv), 1e-300) * max(float(gw @ gw), 1e-300)))
                out.append(row)
    def med(sid, half, key):
        return float(np.median([x[key] for x in out if x["set_id"] == sid and x["half"] == half]))
    best = {}
    for sid in a.sets:
        meds = {c: med(sid, "choose", f"m|{c}") for c in cands}
        best[sid] = max(meds, key=meds.get)
        n = sum(1 for x in out if x["set_id"] == sid and x["half"] == "choose")
        print(f"{sid}: choosing half n {n}; median margin m (CAP6) per band " + ", ".join(f"{c} Hz {v:.3f}" for c, v in meds.items())
              + f"; reference {REF} {med(sid, 'choose', f'm|{REF}'):.3f} -> best {best[sid]} Hz")
        print("  descriptive: uncapped " + ", ".join(f"{c} {med(sid, 'choose', f'm_nocap|{c}'):.3f}" for c in cands + [REF]))
        print("  descriptive: with observed arrivals, CAP6 " + ", ".join(
            f"{c} {med(sid, 'choose', f'm_obs|{c}'):.3f}" for c in cands + [REF]))
    key = "m"
    pooled = {c: float(np.mean([med(s, "choose", f"m|{c}") for s in a.sets])) for c in cands}
    if min(pooled.values()) >= 0.95 * max(pooled.values()):
        key = "m_nocap"
        pooled = {c: float(np.mean([med(s, "choose", f"m_nocap|{c}") for s in a.sets])) for c in cands}
        print("no difference with the cap (every candidate within 5 % of the best CAP6 median margin): the choice and "
              "the confirmation use the uncapped margin")
    chosen = max(pooled, key=pooled.get)
    best = {s: max(cands, key=lambda c: med(s, "choose", f"{key}|{c}")) for s in a.sets}
    note = ""
    for s in a.sets:
        if best[s] != chosen and med(s, "choose", f"{key}|{best[s]}") > 1.10 * med(s, "choose", f"{key}|{chosen}"):
            note = "; a set prefers another band by more than 10 %: reported, the owner decides"
    print(f"chosen band (largest mean of the per-set median margins, {key}): {chosen} Hz" + note)
    rng = np.random.default_rng(SEED_BOOT3)
    for sid in a.sets:
        r = [x for x in out if x["set_id"] == sid and x["half"] == "confirm"]
        labs = sorted({x["label"] for x in r})
        mc = {l: np.array([x[f"{key}|{chosen}"] for x in r if x["label"] == l]) for l in labs}
        mr = {l: np.array([x[f"{key}|{REF}"] for x in r if x["label"] == l]) for l in labs}
        gain = float(np.median(np.concatenate(list(mc.values()))) / np.median(np.concatenate(list(mr.values()))) - 1)
        bs = []
        for _ in range(B3):
            ll = rng.choice(labs, len(labs))
            bs.append(np.median(np.concatenate([mc[l] for l in ll])) / np.median(np.concatenate([mr[l] for l in ll])) - 1)
        lo, hi = np.quantile(bs, [0.025, 0.975])
        zone = ("passes: adopt the band" if gain >= GAIN_PASS else
                ("useful smaller gain: the owner decides" if gain > 0 else "no gain: keep the current band (with weights)"))
        fe_c = float(np.median([x[f"f_emp|{chosen}"] for x in r]))
        print(f"{sid}: confirming half n {len(r)} targets in {len(labs)} realisations; median margin {chosen} Hz "
              f"{np.median(np.concatenate(list(mc.values()))):.3f} vs reference {np.median(np.concatenate(list(mr.values()))):.3f}; "
              f"gain {gain:+.1%} [95 % {lo:+.1%}, {hi:+.1%}] -> {zone} (card rule: >= +20 % / 0-20 % / <= 0)")
        print(f"  reading check: f_emp (g/sd) {chosen} Hz median {fe_c:.3f}"
              + (" -> BELOW 0.4: flagged, not adopted without a look" if fe_c < 0.4 else " (>= 0.4)")
              + f"; reference f_emp {np.median([x[f'f_emp|{REF}'] for x in r]):.3f}")
    json.dump({"rows": out, "chosen": chosen, "measure": key}, open(a.out, "w"), indent=1)


def gdiagw_stage(a):
    """supervisor, 17 Nov, item 4 (descriptive): N1-STACK gdiag recomputed with g and Q weighted by 1/sd (the step-2
    fix); same hits, front end, mask and reading as gdiag (028aa3f / 2508ddf)."""
    from ttmt.synth import realdata as RD
    from ttmt.synth import sets
    from ttmt.radiation import basis, radiation_operator
    n1 = _n1()
    pops = json.load(open(a.pops))
    ctl = [dict(x) for x in pops["control"]]
    out, cache = [], {}
    for sid in sorted({x["set_id"] for x in ctl}):
        rd = sets.SetReader(sid, "tune", verify="fingerprint")
        STA = np.asarray(rd.man["xyz"], float)
        model = RD.pred_model(rd.man["dataset"])
        tr0 = rd.truth(rd.labels[0])
        _s, noise = rd.arrays(rd.labels[0])
        mfn = n1.matched(n1.pre16(noise))
        sd = n1.station_sd(mfn)
        keep = (~n1.catalogue_mask(tr0.t0, mfn.shape[-1], tr0, STA, model))[:, None, :].astype(float)
        cache[sid] = (rd, STA, sd, keep, mfn)
    for x in ctl:
        rd, STA, sd, keep, mfn = cache[x["set_id"]]
        sig, _n = rd.arrays(x["label"])
        truth = rd.truth(x["label"])
        n1.check_station_order(rd, truth)
        R = (n1.matched(n1.pre16(sig)) + mfn * keep) / sd[:, None, None]
        e = x["event"]
        v = n1.read_at_best(R, truth.arr_t_obs[e], truth.t0).reshape(-1)
        g = np.array(truth.meta["targets"][e]["g"], float)
        g = np.nan_to_num(g.transpose(1, 0, 2)); g[~np.isfinite(truth.arr_t_obs[e].T)] = 0.0; g[..., 2] *= -1.0
        w = np.tile(np.repeat(1.0 / sd, 3), 2)
        vv = float(v @ v)
        y = v @ basis(w[:, None] * radiation_operator(truth.xyz[e], STA))
        gw = g.reshape(-1) * w
        out.append(dict(x, f_emp_w=float((v @ gw) ** 2 / (vv * max(float(gw @ gw), 1e-300))), f_an_w=float(y @ y) / vv))
    for sid in sorted({x["set_id"] for x in out}):
        r = [x for x in out if x["set_id"] == sid]
        fe = np.array([x["f_emp_w"] for x in r]); fa = np.array([x["f_an_w"] for x in r])
        print(f"{sid}: n {len(r)}; weighted f_emp median {np.median(fe):.3f} [q25 {np.quantile(fe, .25):.3f}, q75 "
              f"{np.quantile(fe, .75):.3f}]; weighted f_an median {np.median(fa):.3f} (descriptive; unweighted in gdiag2)")
    json.dump({"rows": out}, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="stage", required=True)
    s = sub.add_parser("sp"); s.add_argument("gcheck"); s.add_argument("--out", default="sp.json")
    c = sub.add_parser("chain"); c.add_argument("gcheck"); c.add_argument("--out", default="chain.json")
    c.add_argument("--laptop", action="store_true", help="read each set's noise block for the per-station sd")
    b3 = sub.add_parser("band"); b3.add_argument("--sets", nargs="+",
                                                  default=["ss_hinet_2012_egf_v1", "ss_ridgecrest_2019_egf_v1"])
    b3.add_argument("--out", default="band.json")
    gw = sub.add_parser("gdiagw"); gw.add_argument("pops"); gw.add_argument("--out", default="gdiagw.json")
    a = p.parse_args()
    {"sp": sp_stage, "chain": chain_stage, "band": band_stage, "gdiagw": gdiagw_stage}[a.stage](a)
