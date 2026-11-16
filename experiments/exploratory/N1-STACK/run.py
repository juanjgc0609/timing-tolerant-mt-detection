"""N1-STACK: pre-test of the signed stack for weak events (proposal N-1; card.md, declared before running).

Question: of the T targets that the picker does not see even alone (N2-REPICK: fewer than 6 picked arrivals after
the exact subtraction of the partner), how many rise above the noise when the evidence of the whole network is
summed with the signs a common moment tensor predicts, at the registered budget of 30 false detections per hour?
Ideal case: true location, injected arrival times (arr_t_obs), origin tolerance +-0.25 s.

Statistics (one front end, the one of R1: demean, detrend, Butterworth 4 zero-phase 0.5-6 Hz, 16 Hz, Z down, the
2 Hz Ricker matched filter of ttmt.kernel, every station divided by the robust sd of its filtered NOISE block):
  S_mt  0.5 ||Q^T v||^2, Q the orthonormal basis of the moment-tensor radiation operator (ttmt.radiation) at the
        hypothesis location and all stations of the record: the projected energy of the v2 detector without cap
        (MT-FW NOCAP, ttmt.detector.ProjectedEnergy), read without the picker.
  S_en  0.5 ||v||^2 on the same reads: energy summed over stations and phases, no radiation pattern (QuakeMigrate's
        idea, no sign).
  v     the 3-component vector of every block (station, phase) at the energy peak within [-1, +8] samples of its
        predicted onset (-0.06 to +0.5 s, as inject.radiated_vectors); absent arrival -> zero vector. Origin
        tolerance: maximum over common shifts of -4..+4 samples (+-0.25 s), each statistic its own maximum.

Threshold (card.md, design revision of 16 November 2026, before running): the tuning partition has ONE 300 s noise
block per set. Two nulls on that block, the same reading rule, n_loc = 100 positions uniform in the box of the set's
target sources, travel times of the set's prediction model (realdata.pred_model); P and S arrival windows of the
block's catalogued events masked per station BEFORE anything else ([tP - 1 s, tS + 8 s], fraction reported):
  slides     100 time slides: station i rolled by 20 s x m_i, the m_i distinct among stations (slots of 20 s in
             300 s), seed 20261123; local maxima of max-over-positions separated >= 1 s; threshold = the value
             exceeded 30 times per hour of slid noise (check: recount within 20-45 per hour).
  unslid     the block as it is, origins whose reads stay inside it: threshold = its 3rd highest local maximum.
  main       the larger of the two (the zone is read with it). Both are reported, with their ratio.
Sensitivity: n_loc = 1000 against n_loc = 100 on the first 5 slides.

Stages:
  baseline   tuning manifests only: reproduces the frozen objective of PhaseNet+PyOcto (stops otherwise), lists P-C (C
             targets with stage "picker", STAGE-LOSS rule) and the positive control (up to 60 targets the frozen
             comparator detects, label order).  python run.py baseline <tune_phasenet.json>... --out pops.json
  stats      set files and the real records (laptop). P-T from N2-REPICK/repick.json (n_repick < 6), measured on the
             record with the partner subtracted (N2-REPICK rebuild, which checks the events add up to the stored
             signal); P-C and control on the unmodified record.
             python run.py stats pops.json <N2-REPICK/repick.json> [--limit 2] --out stats.json
             --limit L: the small version, first L labels in label order across both sets, plus the set noise of
             those sets; checks only, no R and no reading printed.

Tuning partition b60d99b3d0bbbbfa only; nothing reads a read-partition file, a reserved window or a ledger.
"""
import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
FROZEN = REPO / "experiments/exploratory/COMPARATORS/frozen/phasenet+pyocto.json"

FS16 = 16
READ = (-1, 8)            # samples around the predicted onset where the energy peak is read
SHIFT = 4                 # +- samples of origin tolerance (0.25 s)
SEP = 16                  # samples between local maxima (1 s)
BUDGET = 30.0             # false detections per hour
N_LOC, N_LOC_SENS, N_SENS_SLIDES = 100, 1000, 5
N_SLIDES, SLOT_S = 100, 20.0
MASK = (-1.0, 8.0)        # s: catalogued-event mask [tP + MASK[0], tS + MASK[1]] per station
UNSLID_RANK = 3           # the unslid threshold is the 3rd highest local maximum
CONTROL_MIN = 0.90
N_CONTROL = 60
SEED_GRID, SEED_SLIDES, SEED_BOOT = 20261124, 20261123, 20261122
MIN_PICKS = 6


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ============================================================================================ front end
def pre16(x):
    """(ns, 3, n) at 100 Hz -> (ns, 3, n * 4 / 25) at 16 Hz, components (E, N, Z down): datasets._preprocess. Every
    step is linear, so pre16(a + b) = pre16(a) + pre16(b)."""
    from scipy.signal import butter, resample_poly, sosfiltfilt
    x = np.array(x, float, copy=True)
    x -= x.mean(-1, keepdims=True)
    tt = np.arange(x.shape[-1], dtype=float)
    tc = tt - tt.mean()
    slope = (x * tc).sum(-1, keepdims=True) / float((tc * tc).sum())       # least-squares line, as np.polyfit deg 1
    x -= slope * tc
    sos = butter(4, [0.5, 6.0], btype="band", fs=100.0, output="sos")
    d = resample_poly(sosfiltfilt(sos, x, axis=-1), 4, 25, axis=-1)
    d[:, 2] *= -1.0
    return d


def matched(d16):
    from ttmt.kernel import LEGACY
    return LEGACY.matched(d16)[0]


def station_sd(mf_noise):
    """robust sd of every station's filtered noise (all components): median |x| / 0.6745."""
    return np.median(np.abs(mf_noise), axis=(1, 2)) / 0.6745


def peak_vectors(R):
    """V (ns, n, 3): for every station and centre sample c, the 3-vector of R at the energy peak in c + [READ]."""
    ns, _, n = R.shape
    offs = np.arange(READ[0], READ[1] + 1)
    ii = np.clip(np.arange(n)[:, None] + offs[None, :], 0, n - 1)
    E = (R * R).sum(1)                                                     # (ns, n)
    k = ii[np.arange(n)[None, :], E[:, ii].argmax(-1)]                     # (ns, n)
    Rt = R.transpose(0, 2, 1)
    return Rt[np.arange(ns)[:, None], k]


# ============================================================================================ statistics
def operator_basis(xyz, STA):
    from ttmt.radiation import basis, radiation_operator
    return basis(radiation_operator(np.asarray(xyz, float), np.asarray(STA, float)))


def target_stats(R, arr_s, t0, Q):
    """(S_mt, S_en) of one target: R (ns, 3, n) normalised filter output from record time t0; arr_s (ns, 2) onsets
    (record s, nan absent); Q (6 ns, r). Rows of v ordered (phase, station, component) as radiation_operator."""
    ns, _, n = R.shape
    V = peak_vectors(R)
    k = np.floor((arr_s.T - t0) * FS16 + 0.5)                              # (2, ns)
    best = [-np.inf, -np.inf]
    for d in range(-SHIFT, SHIFT + 1):
        c = k + d
        ok = np.isfinite(c) & (c >= 0) & (c < n)
        cc = np.where(ok, c, 0).astype(int)
        v = np.where(ok[..., None], V[np.arange(ns)[None, :], cc], 0.0).reshape(-1)
        y = v @ Q
        best[0] = max(best[0], 0.5 * float(y @ y))
        best[1] = max(best[1], 0.5 * float(v @ v))
    return best


def grid(set_xyz, n, seed):
    """n positions uniform in the box of the set's target sources (depth at least 0.5 km)."""
    lo, hi = set_xyz.min(0), set_xyz.max(0)
    p = np.random.default_rng(seed).uniform(lo, hi, (n, 3))
    p[:, 2] = np.maximum(p[:, 2], 0.5)
    return p


def lags(model, P, STA):
    """(N, 2, ns) P and S travel times in samples from every position to every station (surface)."""
    from ttmt.synth.realdata import first_arrival
    STA = np.asarray(STA, float)
    out = np.empty((len(P), 2, len(STA)))
    for g, x in enumerate(P):
        X = np.hypot(STA[:, 0] - x[0], STA[:, 1] - x[1])
        for b, ph in enumerate(("P", "S")):
            out[g, b] = first_arrival(model, float(x[2]), X, ph)[0] * FS16
    return out


def scan(V, L, Qs, circular):
    """max over positions of S_mt and S_en for every origin sample: V (ns, n, 3) peak vectors, L (N, 2, ns) lags in
    samples, Qs (N, 6 ns, r). circular: indices taken modulo n (time slides); else origins whose reads stay inside."""
    ns, n, _ = V.shape
    Lr = np.floor(L + 0.5).astype(int)
    if circular:
        t = np.arange(n)
    else:
        t = np.arange(max(0, -int(Lr.min()) - READ[0]), n - int(Lr.max()) - READ[1])
    mt = np.full(len(t), -np.inf)
    en = np.full(len(t), -np.inf)
    st = np.arange(ns)[None, None, :]
    for g in range(len(L)):
        idx = t[:, None, None] + Lr[g][None]                               # (nt, 2, ns)
        if circular:
            idx = idx % n
        v = V[st, idx].reshape(len(t), -1)                                 # (nt, 6 ns)
        y = v @ Qs[g]
        mt = np.maximum(mt, 0.5 * (y * y).sum(1))
        en = np.maximum(en, 0.5 * (v * v).sum(1))
    return t, mt, en


def local_maxima(x, sep=SEP):
    """values of the greedy non-maximum suppression of x: largest first, none within sep samples of a kept one."""
    order = np.argsort(-x, kind="stable")
    taken = np.zeros(len(x), bool)
    out = []
    for i in order:
        if taken[i]:
            continue
        out.append(x[i])
        taken[max(0, i - sep + 1):i + sep] = True
    return np.array(out)


def rate_threshold(maxima, hours, budget=BUDGET):
    """the value exceeded budget x hours times (the k-th largest maximum, k = floor(budget x hours))."""
    m = np.sort(np.asarray(maxima))[::-1]
    k = int(np.floor(budget * hours))
    return float(m[min(k, len(m) - 1)])


def catalogue_mask(noise_t0, n, truth, STA, model):
    """(ns, n) bool, True where a catalogued event of the block can be on that station: [tP - 1 s, tS + 8 s]."""
    from ttmt.synth.realdata import first_arrival
    STA = np.asarray(STA, float)
    m = np.zeros((len(STA), n), bool)
    for e in np.flatnonzero(~truth.target):
        x = truth.xyz[e]
        X = np.hypot(STA[:, 0] - x[0], STA[:, 1] - x[1])
        z = max(float(x[2]), 0.5)
        tp = truth.t[e] + first_arrival(model, z, X, "P")[0]
        ts = truth.t[e] + first_arrival(model, z, X, "S")[0]
        for i in range(len(STA)):
            a = int(np.floor((tp[i] + MASK[0] - noise_t0) * FS16))
            b = int(np.ceil((ts[i] + MASK[1] - noise_t0) * FS16))
            m[i, max(a, 0):max(min(b, n), 0)] = True
    return m


def slide_offsets(ns, n, rng):
    slots = int(round(300.0 / SLOT_S))
    m = rng.permutation(slots)[:ns]
    return (m * int(SLOT_S * FS16)) % n


def noise_null(Rn, mask, P, Lg, Qs, rng, n_slides):
    """thresholds of one set from its masked noise. Returns dict with the slide maxima per slide, hours, unslid."""
    Rm = np.where(mask[:, None, :], 0.0, Rn)
    V = peak_vectors(Rm)
    ns, n, _ = V.shape
    out = {"slides_mt": [], "slides_en": [], "offsets": []}
    for s in range(n_slides):
        off = slide_offsets(ns, n, rng)
        Vs = np.stack([np.roll(V[i], off[i], axis=0) for i in range(ns)])
        _t, mt, en = scan(Vs, Lg, Qs, circular=True)
        out["slides_mt"].append(local_maxima(mt))
        out["slides_en"].append(local_maxima(en))
        out["offsets"].append(off.tolist())
    t, mt, en = scan(V, Lg, Qs, circular=False)
    out["unslid_mt"], out["unslid_en"] = local_maxima(mt), local_maxima(en)
    out["hours_slide"] = n / FS16 / 3600.0
    out["hours_unslid"] = len(t) / FS16 / 3600.0
    return out


def thresholds(nl, rng_boot, B=2000):
    """main, slide and unslid thresholds per statistic, recount check and a bootstrap interval over slides."""
    res = {}
    S = len(nl["slides_mt"])
    H = S * nl["hours_slide"]
    for k in ("mt", "en"):
        sl = nl["slides_" + k]
        allm = np.concatenate(sl)
        th_s = rate_threshold(allm, H)
        um = nl["unslid_" + k]
        th_u = float(np.sort(um)[::-1][UNSLID_RANK - 1]) if len(um) >= UNSLID_RANK else float("nan")
        bs = []
        for _ in range(B):
            pick = rng_boot.integers(0, S, S)
            bs.append(rate_threshold(np.concatenate([sl[j] for j in pick]), H))
        res[k] = {"slides": th_s, "unslid": th_u, "main": float(np.nanmax([th_s, th_u])),
                  "ratio_unslid_over_slides": th_u / th_s, "recount_per_h": float((allm > th_s).sum() / H),
                  "slides_interval": [float(np.quantile(bs, .025)), float(np.quantile(bs, .975))],
                  "hours_slides": H, "hours_unslid": nl["hours_unslid"],
                  "n_unslid_maxima": int(len(um))}
    return res


# ============================================================================================ baseline stage
def baseline(a):
    SL = _load("sl", REPO / "experiments/exploratory/STAGE-LOSS/stage_loss.py")
    TU = SL.TU
    from ttmt.associate.pipeline import associator
    fz = json.load(open(FROZEN))
    cfg = dict(fz["params"])
    recs = [r for m in a.manifests for r in TU.load_manifest(m, fz["method"])]
    part = TU.partition_hash(recs)
    if part != fz["partition"]:
        sys.exit(f"partition {part} is not the frozen one {fz['partition']}")
    acfg = {k: v for k, v in cfg.items() if k.startswith("as_")}
    cache, dets, pk = {}, [], {}
    for r in recs:
        lab = r.name.split("__")[0]
        if lab not in cache:
            p = SL.phasenet_picks(r, cfg)
            cache[lab] = associator("pyocto").associate(r.data["frame"], p, acfg)[0], p, SL.match(p, r.truth)
        d, p, m = cache[lab]
        dets.append(d)
        pk[r.name], pk[r.name + "#picks"] = m, p
    ev = SL.evaluate(recs, dets)
    st = SL.attribute(ev, recs, dets, pk, int(cfg["as_min_picks"]))
    print(f"partition {part}  T {ev['recall_T']:.4f} (frozen {fz['recall_T']})  C {ev['recall_C']:.4f} "
          f"(frozen {fz['recall_C']})")
    if abs(ev["recall_T"] - fz["recall_T"]) > 5e-4 or abs(ev["recall_C"] - fz["recall_C"]) > 5e-4:
        sys.exit("the frozen objective is not reproduced: stop")
    pc, hit = [], {}
    for j, (r, e, sep, _snr, _tot) in enumerate(ev["meta"]):
        lab = r.name.split("__")[0]
        row = {"record": r.name, "label": lab, "set_id": r.name.split("_0")[0], "event": int(e)}
        if ev["Cm"][j] and st[j] == "picker":
            pc.append(row)
        if ev["got"][j]:
            hit.setdefault((lab, int(e)), dict(row, regime="T" if ev["T"][j] else "C"))
    control = [hit[k] for k in sorted(hit)][:N_CONTROL]
    print(f"P-C (C targets lost at the picker): {len(pc)} rows, {len({(x['label'], x['event']) for x in pc})} unique; "
          f"targets detected: {len(hit)} unique, control takes {len(control)}")
    json.dump({"partition": part, "manifests": [str(m) for m in a.manifests], "pc": pc, "control": control},
              open(a.out, "w"), indent=1)


# ============================================================================================ stats stage
def boot_ratio(rows, key, seed=SEED_BOOT, B=10000):
    """R and its 95 % interval, bootstrap over realisations (labels)."""
    labs = sorted({x["label"] for x in rows})
    num = np.array([sum(x[key] for x in rows if x["label"] == l) for l in labs], float)
    den = np.array([sum(1 for x in rows if x["label"] == l) for l in labs], float)
    i = np.random.default_rng(seed).integers(0, len(labs), (B, len(labs)))
    bs = num[i].sum(1) / np.maximum(den[i].sum(1), 1)
    return float(num.sum() / den.sum()), [float(np.quantile(bs, .025)), float(np.quantile(bs, .975))]


def stats(a):
    from ttmt.synth import realdata as RD
    from ttmt.synth import sets
    N2 = _load("n2", REPO / "experiments/exploratory/N2-REPICK/run.py")
    pops = json.load(open(a.pops))
    rp = json.load(open(a.repick))["rows"]
    pt = [dict(x) for x in rp if x["n_repick"] < MIN_PICKS]
    print(f"P-T: {len(pt)} rows (n_repick < {MIN_PICKS}) of {len(rp)}; P-C: {len(pops['pc'])}; "
          f"control: {len(pops['control'])}", flush=True)
    labels = sorted({x["label"] for x in pt + pops["pc"] + pops["control"]})
    if a.limit:
        firsts = [next(l for l in labels if l.startswith(s)) for s in sorted({l.split("_0")[0] for l in labels})]
        labels = sorted(set(firsts[:a.limit]) | set(labels[:max(0, a.limit - len(firsts))]))
    set_ids = sorted({l.split("_0")[0] for l in labels})
    t00 = time.time()

    # ---------------------------------------------------------------- per set: front end of the noise, null
    S = {}
    for sid in set_ids:
        rd = sets.SetReader(sid, "tune", verify="fingerprint")
        man = rd.man
        STA = np.asarray(man["xyz"], float)
        model = RD.pred_model(man["dataset"])
        tl = [l for l in rd.labels]
        tr0 = rd.truth(tl[0])
        _s, noise = rd.arrays(tl[0])
        mfn = matched(pre16(noise))
        sd = station_sd(mfn)
        Rn = mfn / sd[:, None, None]
        n = Rn.shape[-1]
        mask = catalogue_mask(tr0.t0, n, tr0, STA, model)
        trs = [rd.truth(l) for l in tl]
        txyz = np.concatenate([t.xyz[t.target] for t in trs])
        P = grid(txyz, N_LOC, SEED_GRID)
        Lg = lags(model, P, STA)
        Qs = np.stack([operator_basis(x, STA) for x in P])
        nl = noise_null(Rn, mask, P, Lg, Qs, np.random.default_rng(SEED_SLIDES), N_SLIDES)
        th = thresholds(nl, np.random.default_rng(SEED_BOOT))
        # sensitivity n_loc = 1000 against 100 on the first N_SENS_SLIDES slides (same offsets)
        P2 = grid(txyz, N_LOC_SENS, SEED_GRID)
        L2 = lags(model, P2, STA)
        Q2 = np.stack([operator_basis(x, STA) for x in P2])
        sens = {}
        Vm = peak_vectors(np.where(mask[:, None, :], 0.0, Rn))
        for nm, (LL, QQ) in (("100", (Lg, Qs)), ("1000", (L2, Q2))):
            mx = {"mt": [], "en": []}
            for off in nl["offsets"][:N_SENS_SLIDES]:
                Vs = np.stack([np.roll(Vm[i], off[i], axis=0) for i in range(len(STA))])
                _t, mt, en = scan(Vs, LL, QQ, circular=True)
                mx["mt"].append(local_maxima(mt))
                mx["en"].append(local_maxima(en))
            H5 = N_SENS_SLIDES * nl["hours_slide"]
            sens[nm] = {k: rate_threshold(np.concatenate(mx[k]), H5) for k in mx}
        S[sid] = {"rd": rd, "STA": STA, "sd": sd, "noise_t0": tr0.t0, "th": th, "sens": sens,
                  "mask_fraction": float(mask.mean()), "n_catalogued": int((~tr0.target).sum())}
        for k in ("mt", "en"):
            t = th[k]
            print(f"{sid} S_{k}: threshold slides {t['slides']:.2f} [{t['slides_interval'][0]:.2f}, "
                  f"{t['slides_interval'][1]:.2f}] ({t['hours_slides']:.2f} h, recount {t['recount_per_h']:.1f}/h), "
                  f"unslid {t['unslid']:.2f} ({t['hours_unslid']:.3f} h, {t['n_unslid_maxima']} maxima), "
                  f"ratio {t['ratio_unslid_over_slides']:.2f}, main {t['main']:.2f}; "
                  f"n_loc 1000/100 on {N_SENS_SLIDES} slides {sens['1000'][k] / sens['100'][k]:.3f}", flush=True)
        print(f"{sid}: catalogued events in the block {S[sid]['n_catalogued']}, masked fraction "
              f"{S[sid]['mask_fraction']:.3f}; {time.time() - t00:.0f} s", flush=True)
        if not 20.0 <= th["mt"]["recount_per_h"] <= 45.0 or not 20.0 <= th["en"]["recount_per_h"] <= 45.0:
            sys.exit("the slide null does not recount to 20-45 per hour: stop")

    # ---------------------------------------------------------------- per label: targets
    rows = {"P-T": pt, "P-C": [dict(x) for x in pops["pc"]], "control": [dict(x) for x in pops["control"]]}
    for li, lab in enumerate(labels):
        sid = lab.split("_0")[0]
        s = S[sid]
        rd, STA = s["rd"], s["STA"]
        sig, noise = rd.arrays(lab)
        truth = rd.truth(lab)
        norm = lambda x: matched(pre16(x)) / s["sd"][:, None, None]   # noqa: E731
        Rsig, Rnoi = norm(sig), norm(noise)
        R0 = Rsig + Rnoi                                                   # linear front end: unmodified record
        need_pt = [x for x in pt if x["label"] == lab]
        comps = None
        if need_pt:
            _data, comps, _t, _codes, _fs = N2.rebuild(sid)(lab)          # checks the events add up to the signal
        for pop, rr in rows.items():
            for x in rr:
                if x["label"] != lab:
                    continue
                e = x["event"]
                R = R0
                if pop == "P-T":
                    R = R0 - norm(comps[x["partner"]].astype(float))
                Q = operator_basis(truth.xyz[e], STA)
                x["S_mt"], x["S_en"] = target_stats(R, truth.arr_t_obs[e], truth.t0, Q)
                x["n_blocks"] = int(np.isfinite(truth.arr_t_obs[e]).sum())
        print(f"{lab}: done ({li + 1}/{len(labels)}, {time.time() - t00:.0f} s)", flush=True)
    for pop in rows:
        rows[pop] = [x for x in rows[pop] if "S_mt" in x]
        for x in rows[pop]:
            th = S[x["set_id"]]["th"]
            x["above_mt"] = bool(x["S_mt"] > th["mt"]["main"])
            x["above_en"] = bool(x["S_en"] > th["en"]["main"])
            x["above_mt_slides"] = bool(x["S_mt"] > th["mt"]["slides"])
            x["above_mt_unslid"] = bool(x["S_mt"] > th["mt"]["unslid"])
            x["diff"] = int(x["above_mt"]) - int(x["above_en"])

    ctl = rows["control"]
    fc = float(np.mean([x["above_mt"] for x in ctl])) if ctl else float("nan")
    print(f"\npositive control: {sum(x['above_mt'] for x in ctl)} of {len(ctl)} above the main S_mt threshold "
          f"({fc:.3f}; needed {CONTROL_MIN})")
    out = {"limit": a.limit, "labels": labels, "rows": rows,
           "noise": {sid: {"th": S[sid]["th"], "sens": S[sid]["sens"], "sd": S[sid]["sd"].tolist(),
                           "mask_fraction": S[sid]["mask_fraction"], "n_catalogued": S[sid]["n_catalogued"]}
                     for sid in S},
           "control_fraction": fc}
    json.dump(out, open(a.out, "w"), indent=1)
    if a.limit:
        print("small version: checks only, no R and no reading")
        return
    if not fc >= CONTROL_MIN:
        print("STOP: the positive control fails; the statistic or the threshold is wrong and nothing is read")
        return

    def report(name, rr):
        if not rr:
            return
        R, ci = boot_ratio(rr, "above_mt")
        Re, cie = boot_ratio(rr, "above_en")
        Rs, _ = boot_ratio(rr, "above_mt_slides")
        Ru, _ = boot_ratio(rr, "above_mt_unslid")
        d, cid = boot_ratio(rr, "diff")
        uq = {(x["label"], x["event"]): x for x in rr}
        print(f"{name}: n {len(rr)} ({len(uq)} unique); R(S_mt) {R:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]; "
              f"R(S_en) {Re:.3f} [{cie[0]:.3f}, {cie[1]:.3f}]; paired S_mt - S_en {d:+.3f} [{cid[0]:+.3f}, "
              f"{cid[1]:+.3f}]; R with the slide threshold {Rs:.3f}, with the unslid {Ru:.3f}; unique "
              f"{sum(x['above_mt'] for x in uq.values())} of {len(uq)}")
        return R, ci, d, cid

    res = {}
    for pop in ("P-T", "P-C", "control"):
        res[pop] = report(pop, rows[pop])
        for sid in set_ids:
            res[(pop, sid)] = report(f"  {pop} {sid}", [x for x in rows[pop] if x["set_id"] == sid])
    R, ci, d, cid = res["P-T"]
    zone = ("PASSES the target (>= 0.30)" if R >= 0.30 else
            ("MIDDLE zone: target not met, useful effect (0.10-0.30)" if R >= 0.10 else "LOW zone (< 0.10)"))
    print(f"\nREADING (card.md, rule 36, point estimate): R = {R:.3f} on P-T -> {zone}"
          + ("; the interval crosses a boundary" if (ci[0] < 0.30 <= ci[1]) or (ci[0] < 0.10 <= ci[1]) else ""))
    print("sign adds: " + ("yes (interval above 0)" if cid[0] > 0 else "NOT shown (interval includes 0 or below)"))
    rr = {sid: res[("P-T", sid)][0] for sid in set_ids if res[("P-T", sid)]}
    hi = [k for k in rr if "hinet" in k]
    rc = [k for k in rr if "ridgecrest" in k]
    if hi and rc:
        print(f"Hi-net R {rr[hi[0]]:.3f} vs Ridgecrest R {rr[rc[0]]:.3f}: "
              + ("Hi-net below half: look at Green's functions / time model of Hi-net before going on"
                 if rr[hi[0]] < 0.5 * rr[rc[0]] else "Hi-net not below half"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="stage", required=True)
    b = sub.add_parser("baseline")
    b.add_argument("manifests", nargs="+")
    b.add_argument("--out", default="pops.json")
    s = sub.add_parser("stats")
    s.add_argument("pops")
    s.add_argument("repick")
    s.add_argument("--limit", type=int, default=0)
    s.add_argument("--out", default="stats.json")
    a = ap.parse_args()
    {"baseline": baseline, "stats": stats}[a.stage](a)
