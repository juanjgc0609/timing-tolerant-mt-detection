"""
task1jk_controls.py - Tareas 1j y 1k: los dos controles sinteticos que deciden la forma de G1.

1j (empates). Bajo el protocolo P3 de la Tarea 0 (red por realizacion, umbral relativo al maximo
de la realizacion, calibracion a recall 0.98), la ablacion TOL puntua una fraccion de aciertos
cuyo maximo vale exactamente 1 en varios nodos a la vez; todos se emiten. El modelo restringido
PHYS multiplica por un termino continuo y su maximo es unico. 1j anade ablaciones SIN
polarizacion que desempatan con informacion cinematica continua:

    TOL-TB   fraccion de aciertos - 0.01 * RMS(residuo de los aciertos)/TOL   (solo desempata:
             0.01 < 1/(2N) para N <= 12, asi que no altera ningun orden entre fracciones distintas)
    TOL-K    (1/2N) * sum exp(-r^2 / 2 sigma^2), sigma = TOL/2: puntaje continuo de residuos
    CIII-TB  TOL-TB con el filtro cinematico C-III

y compara contra PHYS (Tarea 0) y PHYS-S (el estimador sumado (20c) de la propuesta).

1k (ventanas vacias). Cada ventana contiene K fuentes, K en {0,1,2} con probabilidades
0.4/0.4/0.2, en posiciones y tiempos de origen aleatorios. El umbral es ABSOLUTO y se calibra en
un conjunto aparte. Se reportan recall por evento, eventos falsos por ventana y la fraccion de
ventanas vacias que emiten algun evento. Es el endpoint que existe en datos continuos reales.

Baselines publicados en ambas: GaMMA y PyOcto, afinados sobre el conjunto de CALIBRACION con una
rejilla de hiperparametros (no solo min_picks) y congelados antes de la evaluacion.

    python task1jk_controls.py --part j --seed 2024 --workers 2
    python task1jk_controls.py --part k --seed 2024 --workers 2
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import argparse
import itertools
import json
import os
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.stats import beta, binomtest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "gamma_src"))
warnings.filterwarnings("ignore")

from task0_common import (DEFAULT_TAUS, DT, NT, TOL, VMIN, VP, VS, _catalog_impl, _precompute,
                          calibrate_threshold, catalog_from_score_field, default_grid,
                          evaluate_catalog, generate, make_network, picks_to_dataframe,
                          polarization, stations_to_dataframe)

TOL_KM = 5.0
TARGET_J = 0.98
TARGET_K = 0.90
S1 = np.array([6.0, -8.0, 9.0])
OWN = ("TOL", "TOL-TB", "TOL-K", "CIII-TB", "PHYS", "PHYS-S", "PHYS-K")
EXT = ("GAMMA", "PYOCTO")

SCEN_J = {
    "unilateral_c15": ("onesided", 8, 15),
    "unilateral_c40": ("onesided", 8, 40),
    "dispersa_c40": ("spread", 8, 40),
    "ladder_3": ("random", 3, 35),
    "ladder_4": ("random", 4, 35),
    "ladder_6": ("random", 6, 35),
}
SCEN_K = {
    "unilateral_c40": ("onesided", 8, 40),
    "dispersa_c40": ("spread", 8, 40),
    "ladder_4": ("random", 4, 35),
    "ladder_6": ("random", 6, 35),
}

GAMMA_GRID = [dict(min_picks=m, sigma=s, eps=10.0)
              for m, s in itertools.product([5, 6, 8, 10, 12, 14], [1.0, 2.0])]
GAMMA_GRID_K = GAMMA_GRID
PYOCTO_GRID = [dict(n_picks=n, tol=t)
               for n, t in itertools.product([4, 6, 8, 10, 12, 14], [0.5, 1.0, 2.0])]


# --------------------------------------------------------------------------- campos propios
def fields_own(STA, D, picks, grid, taus):
    ns = len(STA)
    # una estacion sin ningun pick (posible en ventanas vacias) recibe un pick ficticio muy lejos
    # de la ventana: no produce aciertos y su residuo es enorme, que es lo que significa no tener pick
    picks = [pl if len(pl) else [(-1000.0, "P")] for pl in picks]
    PD, HIT, NEAR = _precompute(D, picks)
    DIJ = np.linalg.norm(STA[:, None, :] - STA[None, :, :], axis=2)
    st = np.arange(ns)
    sig = TOL / 2.0
    out = {m: np.zeros(len(grid)) for m in OWN}
    for gi, src in enumerate(grid):
        d = np.linalg.norm(STA - src, axis=1)
        pred = [[polarization(src, STA[i], bb) for i in range(ns)] for bb in ("P", "S")]
        best = dict.fromkeys(OWN, 0.0)
        for tau in taus:
            tt = np.stack([d / VP, d / VS]) + tau
            idx = np.clip((tt / DT).astype(int), 0, NT - 1)
            cl = NEAR[st[None, :], idx]
            res = np.abs(cl - tt)
            best["TOL-K"] = max(best["TOL-K"], float(np.exp(-res ** 2 / (2 * sig ** 2)).sum()) / (2 * ns))
            hit = HIT[st[None, :], idx]
            if hit.sum() == 0:
                continue
            s = hit.sum() / (2 * ns)
            rms = float(np.sqrt(np.mean(res[hit] ** 2)))
            s_tb = s - 0.01 * rms / TOL
            best["TOL"] = max(best["TOL"], s)
            best["TOL-TB"] = max(best["TOL-TB"], s_tb)
            bad = False
            for b in range(2):
                ii = np.where(hit[b])[0]
                if len(ii) > 1 and np.any(np.abs(cl[b][ii][:, None] - cl[b][ii][None, :])
                                          > DIJ[np.ix_(ii, ii)] / VMIN + 1e-9):
                    bad = True
                    break
            if bad:
                continue
            best["CIII-TB"] = max(best["CIII-TB"], s_tb)
            cs = [abs(float(np.dot(pred[b][i], PD[i, idx[b, i]])))
                  for b in range(2) for i in range(ns) if hit[b, i]]
            if cs:
                best["PHYS"] = max(best["PHYS"], s * float(np.mean(cs)) ** 2)
                best["PHYS-S"] = max(best["PHYS-S"], float(np.sum(np.square(cs))) / (2 * ns))
            cos2 = np.array([[float(np.dot(pred[b][i], PD[i, idx[b, i]])) ** 2 for i in range(ns)]
                             for b in range(2)])
            best["PHYS-K"] = max(best["PHYS-K"],
                                 float((np.exp(-res ** 2 / (2 * sig ** 2)) * cos2).sum()) / (2 * ns))
        for m in OWN:
            out[m][gi] = best[m]
    return out


# --------------------------------------------------------------------------- baselines
def _domain(STA):
    return ([float(STA[:, 0].min() - 30), float(STA[:, 0].max() + 30)],
            [float(STA[:, 1].min() - 30), float(STA[:, 1].max() + 30)], [0.0, 30.0])


def run_gamma(STA, picks, min_picks, sigma, eps):
    from gamma.utils import association
    import pandas as pd
    pk = picks_to_dataframe(picks, "gamma")
    if len(pk) < min_picks:
        return []
    stdf = stations_to_dataframe(STA, "gamma")
    xr, yr, zr = _domain(STA)
    cfg = {"dims": ["x(km)", "y(km)", "z(km)"], "x(km)": xr, "y(km)": yr, "z(km)": zr,
           "vel": {"p": VP, "s": VS}, "method": "BGMM", "oversample_factor": 5,
           "use_amplitude": False, "use_dbscan": True, "dbscan_eps": float(eps),
           "dbscan_min_samples": 3, "min_picks_per_eq": int(min_picks),
           "max_sigma11": float(sigma), "max_sigma22": 1.0, "max_sigma12": 1.0,
           "bfgs_bounds": ((xr[0] - 1, xr[1] + 1), (yr[0] - 1, yr[1] + 1), (zr[0], zr[1]),
                           (None, None)),
           "covariance_prior": [5.0, 5.0], "initial_points": [1, 1, 1]}
    try:
        cat, _a = association(pk, stdf, cfg, method="BGMM")
    except Exception:                                                   # noqa: BLE001
        return []
    if not cat:
        return []
    cdf = pd.DataFrame(cat)
    return [np.array([r["x(km)"], r["y(km)"], r["z(km)"]], float) for _, r in cdf.iterrows()]


def run_pyocto(STA, picks, n_picks, tol):
    import pyocto
    pk = picks_to_dataframe(picks, "pyocto")
    stdf = stations_to_dataframe(STA, "pyocto")
    xr, yr, zr = _domain(STA)
    diag = float(np.hypot(xr[1] - xr[0], yr[1] - yr[0]))
    vm = pyocto.VelocityModel0D(p_velocity=VP, s_velocity=VS, tolerance=float(tol))
    k = max(1, int(n_picks) // 3)
    assoc = pyocto.OctoAssociator(xlim=tuple(xr), ylim=tuple(yr), zlim=tuple(zr), velocity_model=vm,
                                  time_before=float(2.0 * (diag + zr[1]) / VS + 60.0),
                                  n_picks=int(n_picks), n_p_picks=k, n_s_picks=k,
                                  n_p_and_s_picks=k, n_threads=1)
    try:
        events, _a = assoc.associate(pk, stdf)
    except Exception:                                                   # noqa: BLE001
        return []
    if events is None or len(events) == 0:
        return []
    return [np.array([r["x"], r["y"], r["z"]], float) for _, r in events.iterrows()]


def ext_runner(name, part="j"):
    if name == "GAMMA":
        return run_gamma, (GAMMA_GRID if part == "j" else GAMMA_GRID_K)
    return run_pyocto, PYOCTO_GRID


# --------------------------------------------------------------------------- evaluacion
def score_events(events, truth):
    """(eventos verdaderos recuperados, eventos emitidos a mas de TOL_KM de toda verdad)."""
    hit = set()
    nf = 0
    for e in events:
        if not truth:
            nf += 1
            continue
        dd = [float(np.linalg.norm(np.asarray(e)[:2] - np.asarray(s)[:2])) for s in truth]
        j = int(np.argmin(dd))
        if dd[j] <= TOL_KM:
            hit.add(j)
        else:
            nf += 1
    return len(hit), nf


def mcnemar(a, b):
    a, b = np.asarray(a, int), np.asarray(b, int)
    n10 = int(np.sum((a == 1) & (b == 0)))
    n01 = int(np.sum((a == 0) & (b == 1)))
    p = 1.0 if n10 + n01 == 0 else float(binomtest(n10, n10 + n01, 0.5).pvalue)
    return n10, n01, p


def sign_test(x, y):
    d = np.asarray(x, float) - np.asarray(y, float)
    pos, neg = int(np.sum(d > 0)), int(np.sum(d < 0))
    p = 1.0 if pos + neg == 0 else float(binomtest(pos, pos + neg, 0.5).pvalue)
    return pos, neg, p


def cp_ci(k, n, alpha=0.05):
    lo = beta.ppf(alpha / 2, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(1 - alpha / 2, k + 1, n - k) if k < n else 1.0
    return float(lo), float(hi)


# --------------------------------------------------------------------------- parte j
def run_cell_j(args):
    scen, reps, ncal, seed, with_ext = args
    kind, nst, nclut = SCEN_J[scen]
    grid = default_grid()
    rng = np.random.default_rng(seed)
    t0 = time.time()

    def dataset(n, base):
        out = []
        for k in range(n):
            S = make_network(kind, nst, seed=seed + base + k)
            D, p, tr = generate(S, [S1], [6.0], n_clutter=nclut, rng=rng)
            out.append((S, D, p, tr))
        return out

    cal, ev = dataset(ncal, 1000), dataset(reps, 5000)
    fcal = [fields_own(S, D, p, grid, DEFAULT_TAUS) for S, D, p, _t in cal]
    fev = [fields_own(S, D, p, grid, DEFAULT_TAUS) for S, D, p, _t in ev]
    res = {"part": "j", "scenario": scen, "seed": seed, "reps": reps, "ncal": ncal, "methods": {}}
    for m in OWN:
        thr, rc, ok = calibrate_threshold([f[m] for f in fcal], [c[3] for c in cal], grid,
                                          target_recall=TARGET_J, tol_km=TOL_KM)
        rec, fal, nev = [], [], []
        for f, (_S, _D, _p, tr) in zip(fev, ev):
            cat = catalog_from_score_field(grid, f[m], thr)
            r, nf = evaluate_catalog(cat, tr, TOL_KM)
            rec.append(int(r))
            fal.append(int(nf > 0))
            nev.append(len(cat))
        res["methods"][m] = {"thr": thr, "recall_cal": rc, "reachable": bool(ok), "recall": rec,
                             "false": fal, "n_emitted": nev}
    if with_ext:
        for name in EXT:
            runner, gridp = ext_runner(name)
            best = None
            for cfg in gridp:
                rec = fal = 0
                for S, _D, p, tr in cal:
                    r, nf = evaluate_catalog(runner(S, p, **cfg), tr, TOL_KM)
                    rec += int(r)
                    fal += int(nf > 0)
                key = (rec / ncal >= TARGET_J, -fal if rec / ncal >= TARGET_J else rec, )
                if best is None or key > best[0]:
                    best = (key, cfg, rec / ncal)
            cfg = best[1]
            rec, fal, nev = [], [], []
            for S, _D, p, tr in ev:
                cat = runner(S, p, **cfg)
                r, nf = evaluate_catalog(cat, tr, TOL_KM)
                rec.append(int(r))
                fal.append(int(nf > 0))
                nev.append(len(cat))
            res["methods"][name] = {"config": cfg, "recall_cal": best[2],
                                    "reachable": bool(best[0][0]), "recall": rec, "false": fal,
                                    "n_emitted": nev}
    res["seconds"] = time.time() - t0
    return res


# --------------------------------------------------------------------------- parte k
def draw_sources(rng):
    k = int(rng.choice([0, 1, 2], p=[0.4, 0.4, 0.2]))
    srcs, taus = [], []
    while len(srcs) < k:
        s = np.array([rng.uniform(-15, 15), rng.uniform(-15, 15), rng.uniform(5, 13)])
        if all(np.linalg.norm(s[:2] - q[:2]) > 12.0 for q in srcs):
            srcs.append(s)
            taus.append(float(rng.uniform(5.5, 7.5)))
    return srcs, taus


JITTER = float(os.environ.get("TASK1K_JITTER", "0"))


def run_cell_k(args):
    scen, reps, ncal, seed, with_ext = args
    kind, nst, nclut = SCEN_K[scen]
    grid = default_grid()
    rng = np.random.default_rng(seed + 31)
    rng_j = np.random.default_rng(seed + 99)
    t0 = time.time()

    def dataset(n, base):
        out = []
        for k in range(n):
            S = make_network(kind, nst, seed=seed + base + k)
            srcs, taus = draw_sources(rng)
            D, p, tr = generate(S, srcs, taus, n_clutter=nclut, rng=rng)
            if JITTER > 0:
                p = [[(t + float(rng_j.normal(0, JITTER)), b) for (t, b) in pl] for pl in p]
            out.append((S, D, p, tr))
        return out

    cal, ev = dataset(ncal, 1000), dataset(reps, 5000)
    fcal = [fields_own(S, D, p, grid, DEFAULT_TAUS) for S, D, p, _t in cal]
    fev = [fields_own(S, D, p, grid, DEFAULT_TAUS) for S, D, p, _t in ev]
    n_true_cal = sum(len(c[3]) for c in cal)
    res = {"part": "k", "scenario": scen, "seed": seed, "reps": reps, "ncal": ncal, "jitter": JITTER,
           "k_eval": [len(e[3]) for e in ev], "methods": {}}

    def summarize_ops(cats, data):
        hits, nfs = [], []
        for cat, (_S, _D, _p, tr) in zip(cats, data):
            h, nf = score_events(cat, tr)
            hits.append(h)
            nfs.append(nf)
        return hits, nfs

    for m in OWN:
        allmax = np.unique(np.concatenate([f[m][f[m] > 0] for f in fcal]))
        cand = np.unique(np.quantile(allmax, np.linspace(0.0, 1.0, 400)))
        best = None
        for thr in cand[::-1]:
            cats = [_catalog_impl(grid, f[m], float(thr)) for f in fcal]
            hits, nfs = summarize_ops(cats, cal)
            rec = sum(hits) / max(1, n_true_cal)
            if rec >= TARGET_K:
                key = (sum(nfs), -thr)
                if best is None or key < best[0]:
                    best = (key, float(thr), rec)
        reach = best is not None
        if not reach:
            thr = float(cand[0])
            rc = 0.0
        else:
            thr, rc = best[1], best[2]
        cats = [_catalog_impl(grid, f[m], thr) for f in fev]
        hits, nfs = summarize_ops(cats, ev)
        res["methods"][m] = {"thr": thr, "recall_cal": rc, "reachable": reach, "hits": hits,
                             "false": nfs}
    if with_ext:
        for name in EXT:
            runner, gridp = ext_runner(name, "k")
            best = None
            for cfg in gridp:
                cats = [runner(S, p, **cfg) for S, _D, p, _t in cal]
                hits, nfs = summarize_ops(cats, cal)
                rec = sum(hits) / max(1, n_true_cal)
                key = (rec >= TARGET_K, -sum(nfs) if rec >= TARGET_K else rec)
                if best is None or key > best[0]:
                    best = (key, cfg, rec)
            cats = [runner(S, p, **best[1]) for S, _D, p, _t in ev]
            hits, nfs = summarize_ops(cats, ev)
            res["methods"][name] = {"config": best[1], "recall_cal": best[2],
                                    "reachable": bool(best[0][0]), "hits": hits, "false": nfs}
    res["seconds"] = time.time() - t0
    return res


# --------------------------------------------------------------------------- reporte
def report(results):
    L = []
    for r in results:
        L.append(f"\n=== parte {r['part']} | {r['scenario']} | semilla {r['seed']} ===")
        ms = r["methods"]
        if r["part"] == "j":
            n = r["reps"]
            for m, d in ms.items():
                kr, kf = sum(d["recall"]), sum(d["false"])
                lf, hf = cp_ci(kf, n)
                extra = f"  config={d['config']}" if "config" in d else f"  umbral={d['thr']:.3f}"
                L.append(f"   {m:<8} recall {kr:>2}/{n}  falsos {kf:>2}/{n} [{lf:.2f},{hf:.2f}]"
                         f"  emitidos medios {np.mean(d['n_emitted']):.2f}"
                         f"  {'ok' if d['reachable'] else 'NO alcanza'}{extra}")
            for a in [x for x in ms if x != "PHYS"]:
                n10, n01, p = mcnemar(ms[a]["false"], ms["PHYS"]["false"])
                L.append(f"      falsos {a} vs PHYS: {n10}/{n01}  McNemar p={p:.4f}")
        else:
            kev = r["k_eval"]
            ntrue = sum(kev)
            empty = [i for i, k in enumerate(kev) if k == 0]
            L.append(f"   ventanas {len(kev)} (vacias {len(empty)}), eventos verdaderos {ntrue}")
            for m, d in ms.items():
                h, f = sum(d["hits"]), sum(d["false"])
                fa_empty = sum(1 for i in empty if d["false"][i] > 0)
                extra = f"  config={d['config']}" if "config" in d else f"  umbral={d['thr']:.4f}"
                L.append(f"   {m:<8} recall {h:>2}/{ntrue}={h / max(1, ntrue):.2f}"
                         f"  falsos/ventana {f / len(kev):.2f}"
                         f"  vacias con emision {fa_empty}/{len(empty)}"
                         f"  {'ok' if d['reachable'] else 'NO alcanza'}{extra}")
            for a in [x for x in ms if x != "PHYS"]:
                pos, neg, p = sign_test(ms[a]["false"], ms["PHYS"]["false"])
                L.append(f"      falsos por ventana {a} vs PHYS: {pos}/{neg}  signo p={p:.4f}")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", choices=["j", "k"], required=True)
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--reps", type=int, default=None)
    ap.add_argument("--ncal", type=int, default=None)
    ap.add_argument("--scenarios", default=None)
    ap.add_argument("--no-ext", action="store_true")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    scen_all = SCEN_J if args.part == "j" else SCEN_K
    reps = args.reps or (30 if args.part == "j" else 60)
    ncal = args.ncal or (30 if args.part == "j" else 60)
    scens = args.scenarios.split(",") if args.scenarios else list(scen_all)
    fn = run_cell_j if args.part == "j" else run_cell_k
    cells = [(s, reps, ncal, args.seed, not args.no_ext) for s in scens]
    results = []
    with ProcessPoolExecutor(args.workers) as ex:
        for r in ex.map(fn, cells):
            results.append(r)
            print(f"   hecho {r['scenario']} en {r['seconds']:.0f} s", flush=True)
    out = args.out or f"task1{args.part}_seed{args.seed}.json"
    with open(out, "w") as fh:
        json.dump(results, fh)
    print(report(results))
    print(f"\nescrito {out}")


if __name__ == "__main__":
    main()
