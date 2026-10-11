"""
task1j2_noise.py - Tarea 1j-b: la polarizacion frente al mejor puntaje cinematico continuo, con
picks imperfectos.

1j mostro que, con picks perfectos y modelo de velocidad exacto, un puntaje continuo de residuos
de tiempo sin polarizacion (TOL-K) iguala o supera al modelo restringido PHYS. Con picks exactos el
residuo de la fuente verdadera es cero en todos los arribos, lo que hace al residuo un
discriminante casi ideal. Esta tarea mide si la polarizacion anade algo cuando el tiempo deja de
ser exacto: error gaussiano en el tiempo de cada pick y error de velocidad en la inferencia.

Metodos (protocolo P3, mismo que 1j):
    TOL-K    (1/2N) sum exp(-r^2/2s^2), sin polarizacion
    PHYS     el prototipo de la Tarea 0 (fraccion de aciertos x coseno medio al cuadrado, con C-III)
    PHYS-S   la forma sumada (24c)
    PHYS-K   (1/2N) sum exp(-r^2/2s^2) cos^2, con C-III: la polarizacion anadida al mejor puntaje
             cinematico continuo
La escala s de TOL-K y PHYS-K se elige en calibracion entre {0.15, 0.3, 0.6} s.

    python task1j2_noise.py --seed 2024 --workers 2
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from task0_common import (DEFAULT_TAUS, DT, NT, VMIN, VP, VS,
                          _precompute, calibrate_threshold, catalog_from_score_field,
                          default_grid, evaluate_catalog, generate, make_network, polarization)
from task1jk_controls import S1, cp_ci, mcnemar

TOL_KM = 5.0
TARGET = 0.98
SIGMAS = (0.15, 0.3, 0.6)
SCEN = {"unilateral_c40": ("onesided", 8, 40), "ladder_4": ("random", 4, 35),
        "ladder_6": ("random", 6, 35)}
NOISE = {"jit0.1": (0.1, 1.0), "jit0.2": (0.2, 1.0), "jit0.1_v5": (0.1, 1.05)}


def fields(STA, D, picks, grid, taus, vscale):
    ns = len(STA)
    PD, HIT, NEAR = _precompute(D, picks)
    DIJ = np.linalg.norm(STA[:, None, :] - STA[None, :, :], axis=2)
    st = np.arange(ns)
    names = ["PHYS", "PHYS-S"] + [f"TOL-K@{s}" for s in SIGMAS] + [f"PHYS-K@{s}" for s in SIGMAS]
    out = {m: np.zeros(len(grid)) for m in names}
    vp, vs = VP * vscale, VS * vscale
    for gi, src in enumerate(grid):
        d = np.linalg.norm(STA - src, axis=1)
        pred = [[polarization(src, STA[i], bb) for i in range(ns)] for bb in ("P", "S")]
        best = dict.fromkeys(names, 0.0)
        for tau in taus:
            tt = np.stack([d / vp, d / vs]) + tau
            idx = np.clip((tt / DT).astype(int), 0, NT - 1)
            cl = NEAR[st[None, :], idx]
            res = np.abs(cl - tt)
            cos2 = np.array([[float(np.dot(pred[b][i], PD[i, idx[b, i]])) ** 2 for i in range(ns)]
                             for b in range(2)])
            for s in SIGMAS:
                k = np.exp(-res ** 2 / (2 * s ** 2))
                best[f"TOL-K@{s}"] = max(best[f"TOL-K@{s}"], float(k.sum()) / (2 * ns))
            hit = HIT[st[None, :], idx]
            if hit.sum() == 0:
                continue
            bad = False
            for b in range(2):
                ii = np.where(hit[b])[0]
                if len(ii) > 1 and np.any(np.abs(cl[b][ii][:, None] - cl[b][ii][None, :])
                                          > DIJ[np.ix_(ii, ii)] / VMIN + 1e-9):
                    bad = True
                    break
            if bad:
                continue
            s_hit = hit.sum() / (2 * ns)
            c = np.sqrt(cos2[hit])
            best["PHYS"] = max(best["PHYS"], s_hit * float(np.mean(c)) ** 2)
            best["PHYS-S"] = max(best["PHYS-S"], float(np.sum(cos2[hit])) / (2 * ns))
            for s in SIGMAS:
                k = np.exp(-res ** 2 / (2 * s ** 2))
                best[f"PHYS-K@{s}"] = max(best[f"PHYS-K@{s}"], float((k * cos2).sum()) / (2 * ns))
        for m in names:
            out[m][gi] = best[m]
    return out


def jitter(picks, sj, rng):
    return [[(t + float(rng.normal(0, sj)), b) for (t, b) in pl] for pl in picks]


def run_cell(args):
    scen, noise, reps, ncal, seed = args
    kind, nst, nclut = SCEN[scen]
    sj, vscale = NOISE[noise]
    grid = default_grid()
    rng = np.random.default_rng(seed)
    rng_j = np.random.default_rng(seed + 99)
    t0 = time.time()

    def dataset(n, base):
        out = []
        for k in range(n):
            S = make_network(kind, nst, seed=seed + base + k)
            D, p, tr = generate(S, [S1], [6.0], n_clutter=nclut, rng=rng)
            out.append((S, D, jitter(p, sj, rng_j), tr))
        return out

    cal, ev = dataset(ncal, 1000), dataset(reps, 5000)
    fcal = [fields(S, D, p, grid, DEFAULT_TAUS, vscale) for S, D, p, _t in cal]
    fev = [fields(S, D, p, grid, DEFAULT_TAUS, vscale) for S, D, p, _t in ev]
    truth_cal = [c[3] for c in cal]
    res = {"scenario": scen, "noise": noise, "seed": seed, "methods": {}}
    groups = {"PHYS": ["PHYS"], "PHYS-S": ["PHYS-S"],
              "TOL-K": [f"TOL-K@{s}" for s in SIGMAS], "PHYS-K": [f"PHYS-K@{s}" for s in SIGMAS]}
    for g, variants in groups.items():
        best = None
        for v in variants:
            thr, rc, ok = calibrate_threshold([f[v] for f in fcal], truth_cal, grid,
                                              target_recall=TARGET, tol_km=TOL_KM)
            fal = 0
            for f, tr in zip(fcal, truth_cal):
                fal += int(evaluate_catalog(catalog_from_score_field(grid, f[v], thr), tr, TOL_KM)[1] > 0)
            key = (ok, -fal if ok else rc)
            if best is None or key > best[0]:
                best = (key, v, thr, rc, ok)
        _k, v, thr, rc, ok = best
        rec, fal = [], []
        for f, (_S, _D, _p, tr) in zip(fev, ev):
            r, nf = evaluate_catalog(catalog_from_score_field(grid, f[v], thr), tr, TOL_KM)
            rec.append(int(r))
            fal.append(int(nf > 0))
        res["methods"][g] = {"variant": v, "thr": thr, "recall_cal": rc, "reachable": bool(ok),
                             "recall": rec, "false": fal}
    res["seconds"] = time.time() - t0
    return res


def report(results):
    L = []
    for r in results:
        L.append(f"\n=== {r['scenario']} | {r['noise']} | semilla {r['seed']} ===")
        ms = r["methods"]
        for m, d in ms.items():
            n = len(d["recall"])
            kf = sum(d["false"])
            lf, hf = cp_ci(kf, n)
            L.append(f"   {m:<7} recall {sum(d['recall']):>2}/{n}  falsos {kf:>2}/{n} [{lf:.2f},{hf:.2f}]"
                     f"  {d['variant']}  {'ok' if d['reachable'] else 'NO alcanza'}")
        for a, b in (("TOL-K", "PHYS-K"), ("TOL-K", "PHYS"), ("TOL-K", "PHYS-S")):
            n10, n01, p = mcnemar(ms[a]["false"], ms[b]["false"])
            L.append(f"      falsos {a} vs {b}: {n10}/{n01}  McNemar p={p:.4f}")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--ncal", type=int, default=30)
    ap.add_argument("--workers", type=int, default=2)
    args = ap.parse_args()
    cells = [(s, n, args.reps, args.ncal, args.seed) for s in SCEN for n in NOISE]
    results = []
    with ProcessPoolExecutor(args.workers) as ex:
        for r in ex.map(run_cell, cells):
            results.append(r)
            print(f"   hecho {r['scenario']}/{r['noise']} en {r['seconds']:.0f} s", flush=True)
    out = f"task1j2_seed{args.seed}.json"
    with open(out, "w") as fh:
        json.dump(results, fh)
    print(report(results))
    print(f"\nescrito {out}")


if __name__ == "__main__":
    main()
