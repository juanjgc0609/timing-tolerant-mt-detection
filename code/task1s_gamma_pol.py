"""
task1s_gamma_pol.py - Tarea 1s: la fisica encima del estado del arte. GaMMA (cinematico) contra el
MISMO GaMMA con un filtro de polarizacion por estacion sobre sus eventos: cada evento emitido se
puntua con la suma de log f/g sobre sus picks asociados (razon de la ecuacion 9d, estimada en 40
ventanas propias) y se emite solo si supera un umbral. Control negativo: la misma suma con cosenos
permutados entre picks. Mismas ventanas, dominio y ventana de origen que la Tarea 1r; configuracion y
umbral con oraculo sobre la evaluacion para los tres, lo que trata igual a GaMMA y a GaMMA+POL.

    python task1s_gamma_pol.py --seed 2024 --workers 2
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

from task0_common import (DT, EPOCH, NT, VP, VS, _precompute, default_grid, generate, make_network,
                          picks_to_dataframe, polarization, stations_to_dataframe)
from task1o_consistent import NB, draw_sources, pol_hist, score_events
from task1r_baselines import XR, YR, ZR, in_domain
from task1r_fair import CAPS, COND, SCEN, degrade2

def tsec(ts):
    ts = pd.Timestamp(ts)
    ep = pd.Timestamp(EPOCH)
    if ts.tzinfo is None and ep.tzinfo is not None:
        ep = ep.tz_localize(None)
    elif ts.tzinfo is not None and ep.tzinfo is None:
        ts = ts.tz_localize(None)
    return (ts - ep).total_seconds()


GRID = [dict(min_picks=m, sigma=s) for m, s in itertools.product([3, 4, 5, 6, 8], [1.0, 2.0])]


def gamma_with_assign(STA, picks, min_picks, sigma):
    from gamma.utils import association
    pk = picks_to_dataframe(picks, "gamma")
    if len(pk) < min_picks:
        return [], pk
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
        cat, asg = association(pk, stdf, cfg, method="BGMM")
    except Exception:                                                   # noqa: BLE001
        return [], pk
    ev = []
    for e in (cat or []):
        t = tsec(e["time"])
        x, y, z = float(e["x(km)"]), float(e["y(km)"]), float(e["z(km)"])
        if in_domain(x, y, z, t):
            members = [int(pi) for (pi, ei, _p) in asg if ei == e["event_index"]]
            ev.append((np.array([x, y, z]), members))
    return ev, pk


def pol_score(ev_xyz, members, pk, STA, PD, ratio):
    s = 0.0
    for m in members:
        row = pk.iloc[m]
        i = int(str(row["id"])[1:])
        t = tsec(row["timestamp"])
        k = int(np.clip(int(t / DT), 0, NT - 1))
        b = "P" if str(row["type"]).lower() == "p" else "S"
        c = abs(float(PD[i, k] @ polarization(ev_xyz, STA[i], b)))
        s += float(np.log(ratio[min(NB - 1, int(c * NB))]))
    return s


def run_cell(args):
    scen, cond, reps, ncal, seed = args
    kind, nst, nclut = SCEN[scen]
    rng = np.random.default_rng(seed + 41)
    rng_d = np.random.default_rng(seed + 43)
    rng_p = np.random.default_rng(seed + 53)
    t0 = time.time()

    def dataset(n, base):
        out = []
        for k in range(n):
            S = make_network(kind, nst, seed=seed + base + k)
            srcs, taus = draw_sources(rng)
            D, p, tr = generate(S, srcs, taus, n_clutter=nclut, rng=rng)
            p2, lab = degrade2(S, srcs, p, COND[cond], rng_d)
            out.append((S, D, p2, tr, lab))
        return out

    _cal = dataset(ncal, 1000)
    ev = dataset(reps, 5000)
    pol = dataset(40, 3000)
    ratio = np.asarray(pol_hist(pol, default_grid())[0])
    nte = sum(len(e[3]) for e in ev)
    # pool de cosenos observados en picks, para el control negativo
    pool = []
    records = {ci: [] for ci in range(len(GRID))}
    for S, D, p, tr, _l in ev:
        pc = [[(t, b) for (t, b) in pl if t > -100] for pl in p]
        pdp = [pl if pl else [(-1000.0, "P")] for pl in p]
        PD, _H, _N = _precompute(D, pdp)
        for ci, cfg in enumerate(GRID):
            evs, pk = gamma_with_assign(S, pc, **cfg)
            items = []
            for xyz, mem in evs:
                sp = pol_score(xyz, mem, pk, S, PD, ratio)
                items.append((xyz, sp, mem))
                for m in mem:
                    row = pk.iloc[m]
                    i = int(str(row["id"])[1:])
                    t = tsec(row["timestamp"])
                    k = int(np.clip(int(t / DT), 0, NT - 1))
                    b = "P" if str(row["type"]).lower() == "p" else "S"
                    pool.append(abs(float(PD[i, k] @ polarization(xyz, S[i], b))))
            records[ci].append((items, tr))
    pool = np.asarray(pool) if pool else np.array([0.5])
    # puntaje del control negativo: mismos eventos, cosenos tomados al azar del pool
    for ci in records:
        for items, _tr in records[ci]:
            for j, (xyz, sp, mem) in enumerate(items):
                ssh = float(sum(np.log(ratio[min(NB - 1, int(c * NB))]) for c in rng_p.choice(pool, len(mem))))
                items[j] = (xyz, sp, ssh)
    res = {"scenario": scen, "cond": cond, "seed": seed, "k_eval": [len(e[3]) for e in ev], "methods": {}}
    # control positivo del canal: el puntaje de polarizacion debe separar eventos verdaderos de falsos
    from task1o_consistent import TOL_KM
    lab, sco, ssh_ = [], [], []
    for ci in records:
        for items, tr in records[ci]:
            for xyz, sp, ssh in items:
                ok = bool(tr) and min(np.linalg.norm(xyz[:2] - np.asarray(s_)[:2]) for s_ in tr) <= TOL_KM
                lab.append(ok)
                sco.append(sp)
                ssh_.append(ssh)

    def auc(y, x):
        y, x = np.asarray(y, bool), np.asarray(x, float)
        if y.all() or (~y).all():
            return float("nan")
        r = pd.Series(x).rank().to_numpy()
        return float((r[y].sum() - y.sum() * (y.sum() + 1) / 2) / (y.sum() * (~y).sum()))
    res["auc_pol"] = auc(lab, sco)
    res["auc_sh"] = auc(lab, ssh_)
    res["n_true_ev"] = int(np.sum(lab))
    res["n_false_ev"] = int(len(lab) - np.sum(lab))
    for name, key in (("GAMMA", None), ("GAMMA+POL", 1), ("GAMMA+SH", 2)):
        pts, per = [], []
        for ci in records:
            if key is None:
                thr_list = [-np.inf]
            else:
                allv = [it[key] for items, _ in records[ci] for it in items]
                thr_list = [-np.inf] + (list(np.quantile(allv, np.linspace(0, 1, 21))) if allv else [])
            for th in thr_list:
                Hw, Fw = [], []
                for items, tr in records[ci]:
                    kept = [it[0] for it in items if key is None or it[key] >= th]
                    a, b = score_events(kept, tr)
                    Hw.append(a)
                    Fw.append(b)
                pts.append((sum(Fw) / len(ev), sum(Hw) / max(1, nte)))
                per.append((Hw, Fw))
        ra, arg = {}, {}
        for c in CAPS:
            ok = [(rc, j) for j, (fr, rc) in enumerate(pts) if fr <= c]
            ra[str(c)] = max(ok)[0] if ok else 0.0
            arg[str(c)] = max(ok)[1] if ok else -1
        res["methods"][name] = {"recall_at": ra, "points": pts,
                                "chosen": {c: (per[j] if j >= 0 else None) for c, j in arg.items()}}
    res["seconds"] = time.time() - t0
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--reps", type=int, default=60)
    ap.add_argument("--ncal", type=int, default=60)
    ap.add_argument("--workers", type=int, default=2)
    args = ap.parse_args()
    cells = [(s, "consistente", args.reps, args.ncal, args.seed) for s in SCEN]
    out = []
    with ProcessPoolExecutor(args.workers) as ex:
        for r in ex.map(run_cell, cells):
            out.append(r)
            print(f"   hecho {r['scenario']} en {r['seconds']:.0f} s :: " + " | ".join(
                f"{m} " + "/".join(f"{r['methods'][m]['recall_at'][str(c)]:.2f}" for c in CAPS)
                for m in r["methods"]), flush=True)
    with open(f"task1s_seed{args.seed}_r{args.reps}.json", "w") as fh:
        json.dump(out, fh)


if __name__ == "__main__":
    main()
