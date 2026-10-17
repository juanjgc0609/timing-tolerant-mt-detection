"""
c18a_timing.py - exploration C18a (after R0; Ridgecrest CAL segment only, origins in [0, 900) s; the EVAL segment is
not read).  Timing ladder for the R0 failure: how well do different travel-time models, each with a constant term per
station and phase, predict the matched-filter energy peaks of catalogued events that were left out?
    python c18a_timing.py SUBNET MINMAG OUT.json
Readings: station_corrections.peak_delays at the catalogue hypocentre (joint P/S search +-1 s around the homogeneous
prediction, energy gate 25), as in C17a; observed peak time = homogeneous prediction + delay.
Models (leave one event out; every fitted quantity uses only the other events):
  HOM6     registered homogeneous medium VP 6.0, VP/VS 1.73
  HOMfit   homogeneous, VP and VP/VS fitted (grid search, robust misfit after removing the station-phase constant)
  HK       Hadley-Kanamori (1977) layered model of southern California, VP/VS 1.73 (published; no fitting)
  HKr      HK with VP/VS fitted
  GRADfit  VP = v0 + g z (thin layers), v0, g and VP/VS fitted
Prediction of a held-out reading = model time + median residual of the same station and phase over the other events
(nothing when fewer than 3).  Metrics: robust sd, rms, fraction |e| <= 0.25 s (the read half-width of MT-FW).
"""
import json
import sys

import numpy as np

sys.path.insert(0, "/home/claude/work/r0")
import r0_lib as R  # noqa: E402  (sets VP 6.0, VS = VP / 1.73 first)

import c18_tt as TT  # noqa: E402
from station_corrections import peak_delays  # noqa: E402
from task0_common import VP, VS  # noqa: E402


def readings(sub, minmag):
    centre = R.frame()
    cat = R.catalogue_local(centre)
    STA, sens = R.stations(R.SUBNETS[sub], centre)
    _raw, d16 = R.load_waveforms(R.SUBNETS[sub], sens, 940.0)
    d16 = d16 / R.station_scale(d16)[:, None, None]
    rows = []
    evs = [e for e in cat if R.T_CAL[0] <= e["t"] < R.T_CAL[1] and e["M"] >= minmag and e["t"] >= 6.5]
    for k, e in enumerate(evs):
        src = np.array([e["x"], e["y"], e["z"]])
        D = R.window(d16, e["t"] - 6.5)
        dl, en = peak_delays(STA, D, src, 6.5)
        d = np.linalg.norm(STA - src, axis=1)
        for i in range(len(STA)):
            for b, V in ((0, VP), (1, VS)):
                if en[i, b] >= 25:
                    rows.append({"ev": k, "M": e["M"], "src": src.tolist(), "i": i, "b": b,
                                 "tobs": float(d[i] / V + dl[i, b]), "en": float(en[i, b])})
    return STA, evs, rows


def model_times(name, par, STA, src_list):
    """dict ev -> (2, ns) times."""
    out = {}
    for k, src in src_list.items():
        if name.startswith("HOM"):
            vp, r = par
            d = np.linalg.norm(STA - np.asarray(src), axis=1)
            out[k] = np.stack([d / vp, d / (vp / r)])
        elif name.startswith("HK"):
            out[k] = TT.times(TT.HK, par[0], src, STA)
        else:
            v0, g, r = par
            out[k] = TT.times(TT.gradient_model(v0, g), r, src, STA)
    return out


def residuals(rows, tm):
    return np.array([r["tobs"] - tm[r["ev"]][r["b"], r["i"]] for r in rows])


def misfit(rows, res):
    """robust misfit after removing the median per station-phase."""
    key = np.array([r["i"] * 2 + r["b"] for r in rows])
    e = res.copy()
    for k in np.unique(key):
        e[key == k] -= np.median(res[key == k])
    return float(np.median(np.abs(e)))


GRIDS = {
    "HOMfit": [(vp, r) for vp in np.arange(5.4, 6.61, 0.1) for r in np.arange(1.65, 1.86, 0.02)],
    "HKr": [(r,) for r in np.arange(1.65, 1.86, 0.02)],
    "GRADfit": [(v0, g, r) for v0 in np.arange(4.6, 6.01, 0.2) for g in (0.0, 0.03, 0.06, 0.09, 0.12)
                for r in np.arange(1.67, 1.84, 0.03)],
}
FIXED = {"HOM6": (6.0, 1.73), "HK": (1.73,)}


def main(sub, minmag, outp):
    STA, evs, rows = readings(sub, minmag)
    srcs = {k: e and [e["x"], e["y"], e["z"]] for k, e in enumerate(evs)}
    ev_ids = sorted({r["ev"] for r in rows})
    cache = {}

    def tm(name, par):
        key = (name, tuple(np.round(par, 4)))
        if key not in cache:
            cache[key] = model_times(name, par, STA, srcs)
        return cache[key]

    res = {}
    for name in ("HOM6", "HOMfit", "HK", "HKr", "GRADfit"):
        errs, pars, evlab = [], [], []
        for e0 in ev_ids:
            tr = [r for r in rows if r["ev"] != e0]
            te = [r for r in rows if r["ev"] == e0]
            if name in FIXED:
                par = FIXED[name]
            else:
                best = None
                for p in GRIDS[name]:
                    m = misfit(tr, residuals(tr, tm(name, p)))
                    if best is None or m < best[0]:
                        best = (m, p)
                par = best[1]
            pars.append([float(x) for x in par])
            rt = residuals(tr, tm(name, par))
            for r in te:
                m = np.array([(q["i"], q["b"]) == (r["i"], r["b"]) for q in tr])
                if m.sum() < 3:
                    continue
                pred = tm(name, par)[r["ev"]][r["b"], r["i"]] + np.median(rt[m])
                errs.append(r["tobs"] - pred)
                evlab.append(r["ev"])
        v = np.array(errs)
        # event term: the held-out event's residuals minus their own median (a detector estimates origin time itself)
        ev_e = {}
        for r_, e_ in zip(evlab, errs):
            ev_e.setdefault(r_, []).append(e_)
        ve = np.concatenate([np.array(x) - np.median(x) for x in ev_e.values() if len(x) >= 3]) if ev_e else np.zeros(0)
        res.setdefault("_evterm", {})[name] = {"n": int(len(ve)), "robust_sd": float(1.4826 * np.median(np.abs(ve - np.median(ve)))),
                                              "frac_025": float(np.mean(np.abs(ve) <= 0.25))}
        print(f"      with event term: n {len(ve)} robust sd {res['_evterm'][name]['robust_sd']:.3f} |e|<=0.25 "
              f"{res['_evterm'][name]['frac_025']:.2f}", flush=True)
        res[name] = {"n": int(len(v)), "robust_sd": float(1.4826 * np.median(np.abs(v - np.median(v)))),
                     "rms": float(np.sqrt(np.mean(v ** 2))), "frac_025": float(np.mean(np.abs(v) <= 0.25)),
                     "frac_050": float(np.mean(np.abs(v) <= 0.5)), "pars_median": np.median(np.array(pars), 0).tolist()}
        print(f"{sub} M>={minmag} {name:8s} n {len(v):4d}  robust sd {res[name]['robust_sd']:.3f}  rms {res[name]['rms']:.3f}"
              f"  |e|<=0.25 {res[name]['frac_025']:.2f}  |e|<=0.5 {res[name]['frac_050']:.2f}  params {res[name]['pars_median']}",
              flush=True)
    json.dump({"subnet": sub, "minmag": minmag, "n_events": len(ev_ids), "n_readings": len(rows), "models": res, "rows": rows},
              open(outp, "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]), sys.argv[3])
