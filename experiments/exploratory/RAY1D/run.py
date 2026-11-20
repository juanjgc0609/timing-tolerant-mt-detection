"""RAY1D (card.md): radiation along the ray of the sets' 1-D models instead of the straight ray.

Stage `check` (cloud, before everything): noise-free events made with G_ray and read with G_ray through the operative
chain of N1-STACK (pre16, matched, energy peak at the best common shift), in the real geometry and model of each set,
at the event locations of N1-STACK gcheck2.json. Rule: median f >= 0.99 in both sets, otherwise nothing else runs.

Stage `measure` (laptop, tune partition only): every target with snr_med >= 5, the FE-FIX step-3 reference front end
(0.5-6 Hz, 16 Hz, Ricker of 2 Hz), weighted Q (ADR-0012); margin m with CAP6 (main) and uncapped, paired G_ray against
the straight G; f_emp, f_an and f_g (card); bootstrap by realisation. The zone is read by the supervisor.
"""
import argparse, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
SEED_GCHECK = 20261125        # gcheck's generator (FE-FIX chain stage)
SEED_FAKE_SD = 20261126       # cloud check: fake per-station sd for the weighted variant
SEED_BOOT, B = 20261125, 10000
SNR_MIN, QNOISE, CAP = 5.0, 0.999, 6.0
F_MIN = 0.99
GAIN_PASS = 0.20
F_UP, F_FLAT = 0.20, 0.05
N1 = None


def _n1():
    global N1
    if N1 is None:
        import importlib.util
        spec = importlib.util.spec_from_file_location("n1", ROOT / "experiments/exploratory/N1-STACK/run.py")
        N1 = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(N1)
    return N1


def _weights(sd):
    return np.tile(np.repeat(1.0 / np.asarray(sd, float), 3), 2)


def _ops(x, STA, model, sd):
    """weighted operators and bases: straight, ray, ray with free surface (descriptive)."""
    from ttmt.radiation import basis, radiation_operator
    from ttmt.ray_operator import ray_operator
    w = _weights(sd)
    G = {"straight": radiation_operator(x, STA), "ray": ray_operator(x, STA, model),
         "ray_fs": ray_operator(x, STA, model, free_surface=True)}
    return {k: (w[:, None] * g, basis(w[:, None] * g)) for k, g in G.items()}


# ============================================================================================ cloud check
def _event(Gx, m, tt, t0, ns, fs, n):
    """noise-free 100 Hz record of one event (FE-FIX chain_f): Ricker of 4 Hz at every arrival + 0.1 s, Z down."""
    from ttmt.kernel import ricker
    from ttmt.synth import inject as INJ
    g = (Gx @ m).reshape(2, ns, 3)
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
    return unit


def check_stage(a):
    from ttmt.radiation import basis, radiation_operator
    from ttmt.ray_operator import ray_operator
    from ttmt.synth import realdata as RD
    from ttmt.synth.realdata import first_arrival
    n1 = _n1()
    rows = json.load(open(a.gcheck))["rows"]
    sids = sorted({r["set_id"] for r in rows})
    rng = np.random.default_rng(SEED_GCHECK)
    out, ok = [], True
    for sid in sids:
        man = json.loads((ROOT / "data" / "manifests" / "synth" / f"{sid}.json").read_text())
        STA = np.asarray(man["xyz"], float); ns = len(STA)
        model = RD.pred_model(man["dataset"])
        fs, n = float(man["fs"]), 30000
        sdf = np.random.default_rng([SEED_FAKE_SD, sum(map(ord, sid))]).permutation(np.geomspace(1.0, 10.0, ns))
        res = {k: [] for k in ("ray|sd1", "ray|fake_sd_weighted", "ray|fake_sd_unweighted", "ray_fs|sd1",
                               "straight_read_of_ray|sd1")}
        for r in [r for r in rows if r["set_id"] == sid]:
            x = np.asarray(r["xyz"], float)
            m = rng.normal(size=6)
            t0 = float(rng.uniform(20.0, 300.0 - 40.0))
            if abs(t0 - r["t0"]) > 1e-9:
                sys.exit(f"{sid}: redrawn origin {t0} differs from gcheck2.json {r['t0']}: stop")
            X = np.hypot(STA[:, 0] - x[0], STA[:, 1] - x[1]); z = max(float(x[2]), 0.5)
            tt = np.stack([first_arrival(model, z, X, "P")[0], first_arrival(model, z, X, "S")[0]], 1)
            Gr, Gf = ray_operator(x, STA, model), ray_operator(x, STA, model, free_surface=True)
            for G, tag in ((Gr, "ray"), (Gf, "ray_fs")):
                R0 = n1.matched(n1.pre16(_event(G, m, tt, t0, ns, fs, n)))
                cases = [("sd1", np.ones(ns), True)]
                if tag == "ray":
                    cases += [("fake_sd_weighted", sdf, True), ("fake_sd_unweighted", sdf, False)]
                for name, sd, wq in cases:
                    R = R0 / sd[:, None, None]
                    v = n1.read_at_best(R, t0 + tt, 0.0).reshape(-1)
                    w = _weights(sd) if wq else np.ones(6 * ns)
                    y = v @ basis(w[:, None] * G)
                    res[f"{tag}|{name}"].append(float(y @ y) / max(float(v @ v), 1e-300))
                    if tag == "ray" and name == "sd1":
                        y = v @ basis(radiation_operator(x, STA))
                        res["straight_read_of_ray|sd1"].append(float(y @ y) / max(float(v @ v), 1e-300))
        print(f"{sid}: {len(res['ray|sd1'])} noise-free events, f median [min, max]")
        for k, f in res.items():
            f = np.array(f)
            rule = k in ("ray|sd1", "ray|fake_sd_weighted")
            if rule and np.median(f) < F_MIN:
                ok = False
            print(f"  {k:<28s} {np.median(f):.4f} [{f.min():.4f}, {f.max():.4f}]" + ("  (rule: median >= 0.99)" if rule else "  (descriptive)"))
        out.append({"set_id": sid, "f": res})
    print("cloud check " + ("PASSES: the laptop measures may be run" if ok else "FAILS: nothing else runs"))
    json.dump({"rows": out, "passes": ok}, open(a.out, "w"), indent=1)


# ============================================================================================ laptop measures
def _smt(Vp, origins, off_s, Qw, cap):
    """FE-FIX _smt at 16 Hz: S_mt = 1/2 ||Qw^T phi(v)||^2 per origin sample, arrivals at origin + off_s (2, ns) s;
    common shift -4..+4 samples, the one with the largest ||v||^2; phi = CAP6 or none."""
    fr = 16.0
    ns, n, _ = Vp.shape
    off = off_s * fr
    bestE = np.full(len(origins), -1.0); bestv = np.zeros((len(origins), 2, ns, 3))
    for d in range(-4, 5):
        c = np.floor(origins[:, None, None] + off[None] + d + 0.5)
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


def measure_stage(a):
    from ttmt.synth import realdata as RD
    from ttmt.synth import sets
    from ttmt.synth.realdata import first_arrival
    n1 = _n1()
    out = []
    for sid in a.sets:
        rd = sets.SetReader(sid, "tune", verify="fingerprint")
        STA = np.asarray(rd.man["xyz"], float)
        model = RD.pred_model(rd.man["dataset"])
        tr0 = rd.truth(rd.labels[0])
        _s, noise = rd.arrays(rd.labels[0])
        mfn = n1.matched(n1.pre16(noise))
        sd = n1.station_sd(mfn)
        keep = (~n1.catalogue_mask(tr0.t0, mfn.shape[-1], tr0, STA, model))[:, None, :].astype(float)
        Vn = n1.peak_vectors(mfn * keep / sd[:, None, None])
        o_n = np.arange(Vn.shape[1], dtype=float)
        w = _weights(sd)
        for lab in rd.labels:
            sig, _n = rd.arrays(lab)
            truth = rd.truth(lab)
            n1.check_station_order(rd, truth)
            snr_med = np.nanmedian(truth.snr, axis=1)
            tg = [e for e in range(len(truth.t)) if truth.target[e] and np.isfinite(snr_med[e]) and snr_med[e] >= SNR_MIN]
            if not tg:
                continue
            R = (n1.matched(n1.pre16(sig)) + mfn * keep) / sd[:, None, None]
            Vs = n1.peak_vectors(R)
            for e in tg:
                x = truth.xyz[e]
                X = np.hypot(STA[:, 0] - x[0], STA[:, 1] - x[1]); z = max(float(x[2]), 0.5)
                tt = np.stack([first_arrival(model, z, X, "P")[0], first_arrival(model, z, X, "S")[0]], 0)
                g = np.array(truth.meta["targets"][e]["g"], float)
                g = np.nan_to_num(g.transpose(1, 0, 2)); g[~np.isfinite(truth.arr_t_obs[e].T)] = 0.0; g[..., 2] *= -1.0
                gw = g.reshape(-1) * w
                v = n1.read_at_best(R, truth.arr_t_obs[e], truth.t0).reshape(-1)
                vv, gg = max(float(v @ v), 1e-300), max(float(gw @ gw), 1e-300)
                row = {"set_id": sid, "label": lab, "event": int(e), "snr_med": float(snr_med[e]),
                       "pair": int(truth.pair[e]), "f_emp": float((v @ gw) ** 2 / (vv * gg))}
                o_t = np.array([(truth.t[e] - truth.t0) * 16.0])
                for k, (_Gw, Qw) in _ops(x, STA, model, sd).items():
                    for cap, tag in ((True, "m"), (False, "m_nocap")):
                        sn = max(float(np.quantile(_smt(Vn, o_n, tt, Qw, cap), QNOISE)), 1e-300)
                        row[f"{tag}|{k}"] = float(_smt(Vs, o_t, tt, Qw, cap)[0]) / sn
                    y, yg = v @ Qw, gw @ Qw
                    row[f"f_an|{k}"] = float(y @ y) / vv
                    row[f"f_g|{k}"] = float(yg @ yg) / gg
                out.append(row)
    rng = np.random.default_rng(SEED_BOOT)
    summary = {}
    for sid in a.sets:
        r = [x for x in out if x["set_id"] == sid]
        labs = sorted({x["label"] for x in r})
        by = {l: [x for x in r if x["label"] == l] for l in labs}

        def ratio(key, ll):
            num = np.concatenate([[x[f"{key}|ray"] for x in by[l]] for l in ll])
            den = np.concatenate([[x[f"{key}|straight"] for x in by[l]] for l in ll])
            return float(np.median(num) / np.median(den) - 1)

        def dmed(key, ll):
            return float(np.median(np.concatenate([[x[f"{key}|ray"] for x in by[l]] for l in ll]))
                         - np.median(np.concatenate([[x[f"{key}|straight"] for x in by[l]] for l in ll])))

        bs = {k: [] for k in ("m", "m_nocap", "f_g", "f_an")}
        for _ in range(B):
            ll = rng.choice(labs, len(labs))
            for k in ("m", "m_nocap"):
                bs[k].append(ratio(k, ll))
            for k in ("f_g", "f_an"):
                bs[k].append(dmed(k, ll))
        s = {"n": len(r), "realisations": len(labs)}
        for k in ("m", "m_nocap"):
            s[k] = {"gain": ratio(k, labs), "ci95": np.quantile(bs[k], [0.025, 0.975]).tolist(),
                    "median_ray": float(np.median([x[f"{k}|ray"] for x in r])),
                    "median_straight": float(np.median([x[f"{k}|straight"] for x in r])),
                    "median_ray_fs": float(np.median([x[f"{k}|ray_fs"] for x in r]))}
        for k in ("f_g", "f_an"):
            s[k] = {"delta": dmed(k, labs), "ci95": np.quantile(bs[k], [0.025, 0.975]).tolist(),
                    "median_ray": float(np.median([x[f"{k}|ray"] for x in r])),
                    "median_straight": float(np.median([x[f"{k}|straight"] for x in r])),
                    "median_ray_fs": float(np.median([x[f"{k}|ray_fs"] for x in r]))}
        s["f_emp"] = float(np.median([x["f_emp"] for x in r]))
        summary[sid] = s
        print(f"{sid}: n {s['n']} targets in {s['realisations']} realisations")
        for k in ("m", "m_nocap"):
            q = s[k]
            print(f"  {k} median: ray {q['median_ray']:.3f} vs straight {q['median_straight']:.3f} (ray+free surface "
                  f"{q['median_ray_fs']:.3f}, descriptive); ratio of medians - 1 = {q['gain']:+.1%} "
                  f"[95 % {q['ci95'][0]:+.1%}, {q['ci95'][1]:+.1%}]" + ("  (main)" if k == "m" else "  (descriptive)"))
        for k in ("f_g", "f_an"):
            q = s[k]
            print(f"  {k} median: ray {q['median_ray']:.3f} vs straight {q['median_straight']:.3f} (ray+free surface "
                  f"{q['median_ray_fs']:.3f}); difference {q['delta']:+.3f} [95 % {q['ci95'][0]:+.3f}, {q['ci95'][1]:+.3f}]")
        print(f"  f_emp median {s['f_emp']:.3f} (does not depend on G; anchor)")
    print("zones (card; read by the supervisor): m gain >= +20 % in both sets / between, or only one set / <= 0 in both;"
          " polarisation: f rises >= 0.20 supported, < 0.05 not")
    json.dump({"rows": out, "summary": summary}, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="stage", required=True)
    c = sub.add_parser("check"); c.add_argument("gcheck"); c.add_argument("--out", default="check.json")
    m = sub.add_parser("measure"); m.add_argument("--sets", nargs="+",
                                                  default=["ss_hinet_2012_egf_v1", "ss_ridgecrest_2019_egf_v1"])
    m.add_argument("--out", default="measure.json")
    a = p.parse_args()
    {"check": check_stage, "measure": measure_stage}[a.stage](a)
