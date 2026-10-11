"""
task1r_baselines.py - GaMMA y PyOcto sobre las ventanas de la Tarea 1r, con el mismo dominio de
busqueda que los metodos propios (x, y en [-21, 22] km, profundidad 5-13 km, tiempo de origen
5-8 s), minimo de picks desde 3, y la configuracion elegida con ORACULO sobre la evaluacion (la mejor
configuracion para cada tope de falsos realizados), lo que favorece a los baselines.

    python task1r_baselines.py --seed 2024 --workers 2
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import argparse
import itertools
import json
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "gamma_src"))
warnings.filterwarnings("ignore")

from task0_common import EPOCH, VP, VS, generate, make_network, picks_to_dataframe, stations_to_dataframe
from task1o_consistent import draw_sources, score_events
from task1r_fair import CAPS, COND, SCEN, degrade2

XR, YR, ZR, TR = (-21.0, 22.0), (-21.0, 22.0), (5.0, 13.0), (5.0, 8.0)
GAMMA_GRID = [dict(min_picks=m, sigma=s) for m, s in itertools.product([3, 4, 5, 6, 8], [1.0, 2.0])]
METHODS = ["GAMMA", "PYOCTO"]
PYOCTO_GRID = [dict(n_picks=n, tol=t) for n, t in itertools.product([3, 4, 6, 8, 10], [0.5, 1.0, 2.0])]


def in_domain(x, y, z, t):
    return XR[0] <= x <= XR[1] and YR[0] <= y <= YR[1] and ZR[0] - 1 <= z <= ZR[1] + 1 and TR[0] <= t <= TR[1]


def run_gamma(STA, picks, min_picks, sigma):
    from gamma.utils import association
    pk = picks_to_dataframe(picks, "gamma")
    if len(pk) < min_picks:
        return []
    stdf = stations_to_dataframe(STA, "gamma")
    cfg = {"dims": ["x(km)", "y(km)", "z(km)"], "x(km)": list(XR), "y(km)": list(YR), "z(km)": list(ZR),
           "vel": {"p": VP, "s": VS}, "method": "BGMM", "oversample_factor": 5,
           "use_amplitude": False, "use_dbscan": True, "dbscan_eps": 10.0,
           "dbscan_min_samples": 3, "min_picks_per_eq": int(min_picks),
           "max_sigma11": float(sigma), "max_sigma22": 1.0, "max_sigma12": 1.0,
           "bfgs_bounds": ((XR[0] - 1, XR[1] + 1), (YR[0] - 1, YR[1] + 1), ZR, (None, None)),
           "covariance_prior": [5.0, 5.0], "initial_points": [1, 1, 1]}
    try:
        np.random.seed(42)  # GaMMA inicializa con el RNG global de numpy
        cat, _a = association(pk, stdf, cfg, method="BGMM")
    except Exception:                                                   # noqa: BLE001
        return []
    out = []
    for e in (cat or []):
        t = (pd.Timestamp(e["time"]).tz_localize(None) - pd.Timestamp(EPOCH).tz_localize(None)).total_seconds() \
            if pd.Timestamp(e["time"]).tzinfo else (pd.Timestamp(e["time"]) - pd.Timestamp(EPOCH.replace(tzinfo=None))).total_seconds()
        x, y, z = float(e["x(km)"]), float(e["y(km)"]), float(e["z(km)"])
        if in_domain(x, y, z, t):
            out.append(np.array([x, y, z]))
    return out


def run_pyocto(STA, picks, n_picks, tol):
    import pyocto
    pk = picks_to_dataframe(picks, "pyocto")
    stdf = stations_to_dataframe(STA, "pyocto")
    vm = pyocto.VelocityModel0D(p_velocity=VP, s_velocity=VS, tolerance=float(tol))
    k = max(1, int(n_picks) // 3)
    assoc = pyocto.OctoAssociator(xlim=XR, ylim=YR, zlim=ZR, velocity_model=vm, time_before=60.0,
                                  n_picks=int(n_picks), n_p_picks=k, n_s_picks=k,
                                  n_p_and_s_picks=k, n_threads=1)
    try:
        events, _a = assoc.associate(pk, stdf)
    except Exception:                                                   # noqa: BLE001
        return []
    out = []
    if events is None or len(events) == 0:
        return out
    for _, r in events.iterrows():
        t = float(r["time"]) - EPOCH.timestamp()
        if in_domain(float(r["x"]), float(r["y"]), float(r["z"]), t):
            out.append(np.array([r["x"], r["y"], r["z"]], float))
    return out


def run_cell(args):
    scen, cond, reps, ncal, seed = args
    kind, nst, nclut = SCEN[scen]
    rng = np.random.default_rng(seed + 41)
    rng_d = np.random.default_rng(seed + 43)
    t0 = time.time()

    def dataset(n, base):
        out = []
        for k in range(n):
            S = make_network(kind, nst, seed=seed + base + k)
            srcs, taus = draw_sources(rng)
            D, p, tr = generate(S, srcs, taus, n_clutter=nclut, rng=rng)
            p2, _l = degrade2(S, srcs, p, COND[cond], rng_d)
            out.append((S, [[(t, b) for (t, b) in pl if t > -100] for pl in p2], tr))
        return out

    _cal = dataset(ncal, 1000)
    ev = dataset(reps, 5000)
    nte = sum(len(e[2]) for e in ev)
    res = {"scenario": scen, "cond": cond, "seed": seed, "k_eval": [len(e[2]) for e in ev], "methods": {}}
    for name, runner, grid in (("GAMMA", run_gamma, GAMMA_GRID), ("PYOCTO", run_pyocto, PYOCTO_GRID)):
        if name not in METHODS:
            continue
        pts = []
        per = {}
        for ci, cfg in enumerate(grid):
            H, F = [], []
            for S, p, tr in ev:
                a, b = score_events(runner(S, p, **cfg), tr)
                H.append(a)
                F.append(b)
            pts.append((sum(F) / len(ev), sum(H) / max(1, nte)))
            per[ci] = (H, F)
        ra = {}
        for c in CAPS:
            ok = [(rc, ci) for ci, (fr, rc) in enumerate(pts) if fr <= c]
            ra[str(c)] = max(ok)[0] if ok else 0.0
        res["methods"][name] = {"points": pts, "recall_at": ra,
                                "configs": grid, "per": {str(k): v for k, v in per.items()}}
    res["seconds"] = time.time() - t0
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--reps", type=int, default=60)
    ap.add_argument("--ncal", type=int, default=60)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--methods", default="GAMMA,PYOCTO")
    args = ap.parse_args()
    global METHODS
    METHODS = args.methods.split(",")
    cells = [(s, "consistente", args.reps, args.ncal, args.seed) for s in SCEN]
    out = []
    with ProcessPoolExecutor(args.workers) as ex:
        for r in ex.map(run_cell, cells):
            out.append(r)
            print(f"   hecho {r['scenario']} en {r['seconds']:.0f} s :: "
                  + " | ".join(f"{m} " + "/".join(f"{r['methods'][m]['recall_at'][str(c)]:.2f}" for c in CAPS)
                               for m in r["methods"]), flush=True)
    with open(f"task1r_baselines_seed{args.seed}" + ("" if args.reps == 60 else f"_r{args.reps}")
              + ("" if len(METHODS) == 2 else "_" + "_".join(METHODS).lower()) + ".json", "w") as fh:
        json.dump(out, fh)


if __name__ == "__main__":
    main()
