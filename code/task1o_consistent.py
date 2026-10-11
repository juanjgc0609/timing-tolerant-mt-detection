"""
task1o_consistent.py - Tarea 1o: la polarizacion y la forma de onda frente al mejor puntaje
cinematico, en un generador fisicamente consistente.

Motivo. El generador de las Tareas 0 y 1j-1k detecta TODOS los arribos verdaderos con tiempo exacto,
aunque su amplitud sobre el ruido sea 1.1-1.7 en las estaciones lejanas, mientras la polarizacion se
mide sobre esa misma forma de onda ruidosa. Eso regala al canal cinematico informacion que ningun
picker tiene. Aqui la deteccion y el error de tiempo salen de la misma relacion senal/ruido que
degrada la polarizacion:

    deteccion de un arribo verdadero con probabilidad Phi((a - 2.0) / 0.5), a = amplitud / ruido
    error de tiempo gaussiano con sigma = 3 x cota de Cramer-Rao del pulso (minimo 0.02 s)

Condicion de control positivo ("jit0.5"): todos los arribos detectados y error de tiempo de 0.5 s.
Ahi el tiempo casi no informa, asi que si la polarizacion esta bien implementada DEBE ganar; si no
gana, la implementacion del canal es sospechosa y no se concluye nada sobre la fisica.

Metodos, todos con el mismo cuidado (escala elegida en calibracion, umbral absoluto):
    TOL-TB   referencia de la Tarea 1j
    KIN-LR   sum_i log(eps + (1-eps) N(r_i; sigma) / u): verosimilitud cinematica robusta
    PHYS-LR  sum_i log(eps + (1-eps) N(r_i; sigma)/u * f(c_i)/g(c_i)): la misma mas la razon de
             verosimilitud de polarizacion, con f y g estimadas en calibracion (ecuacion 9d)
    WAVE-E   sin picks: energia de las tres componentes del filtro adaptado en el tiempo predicho
    WAVE-P   sin picks: la misma proyectada sobre la polarizacion predicha
Ventanas con 0, 1 o 2 fuentes; endpoint: eventos falsos por ventana a recall por evento objetivo 0.90.

    python task1o_consistent.py --seed 2024 --workers 2
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.stats import binomtest, norm

from task0_common import (AP, AS, DEFAULT_TAUS, DT, NT, SNR, T, TOL, VMIN, VP, VS, _catalog_impl,
                          _precompute, default_grid, generate, make_network, polarization, ricker)

TOL_KM = 5.0
TARGET = 0.90
EPS = 0.1
SIGMAS = (0.1, 0.2, 0.4)
SCEN = {"unilateral_c40": ("onesided", 8, 40), "ladder_4": ("random", 4, 35),
        "ladder_6": ("random", 6, 35)}
COND = {"consistente": "consistent", "jit0.5": "control"}
OWN = ["TOL-TB"] + [f"KIN-LR@{s}" for s in SIGMAS] + [f"PHYS-LR@{s}" for s in SIGMAS] + ["WAVE-E", "WAVE-P"]
GROUPS = {"TOL-TB": ["TOL-TB"], "KIN-LR": [f"KIN-LR@{s}" for s in SIGMAS],
          "PHYS-LR": [f"PHYS-LR@{s}" for s in SIGMAS], "WAVE-E": ["WAVE-E"], "WAVE-P": ["WAVE-P"]}
NOISE_SD = SNR * 0.02
_t = np.arange(-2, 2, DT)
_GD2 = float(np.sum(np.gradient(ricker(_t), DT) ** 2))
_KER = ricker(np.arange(-1.0, 1.0 + 1e-9, DT))
NB = 10


def draw_sources(rng):
    k = int(rng.choice([0, 1, 2], p=[0.4, 0.4, 0.2]))
    srcs, taus = [], []
    while len(srcs) < k:
        s = np.array([rng.uniform(-15, 15), rng.uniform(-15, 15), rng.uniform(5, 13)])
        if all(np.linalg.norm(s[:2] - q[:2]) > 12.0 for q in srcs):
            srcs.append(s)
            taus.append(float(rng.uniform(5.5, 7.5)))
    return srcs, taus


def degrade(STA, srcs, picks, cond, rng):
    """Deteccion y error de tiempo consistentes con la amplitud; devuelve picks y etiquetas."""
    ns = len(STA)
    out, lab = [], []
    k = len(srcs)
    for i in range(ns):
        pl, ll = [], []
        for j, (t, b) in enumerate(picks[i]):
            if j < 2 * k:
                s = srcs[j // 2]
                d = float(np.linalg.norm(STA[i] - s))
                A = AP if b == "P" else AS
                a = SNR * A / d / NOISE_SD
                if cond == "control":
                    pl.append((t + float(rng.normal(0, 0.5)), b))
                    ll.append(j // 2 + 1)
                    continue
                if rng.random() > norm.cdf((a - 2.0) / 0.5):
                    continue
                crb = NOISE_SD / (SNR * A / d * np.sqrt(_GD2))
                pl.append((t + float(rng.normal(0, max(0.02, 3.0 * crb))), b))
                ll.append(j // 2 + 1)
            else:
                pl.append((t, b))
                ll.append(0)
        if not pl:
            pl, ll = [(-1000.0, "P")], [0]
        out.append(pl)
        lab.append(ll)
    return out, lab


def matched(D):
    MF = np.array([[np.convolve(D[i, c], _KER[::-1], mode="same") for c in range(3)]
                   for i in range(D.shape[0])])
    sd = float(np.median(np.abs(MF)) / 0.6745)
    return MF, max(sd, 1e-12)


def pol_hist(cal, grid):
    """f: |cos| de arribos verdaderos detectados; g: |cos| de picks de clutter contra una direccion
    predicha por un nodo al azar de la grilla. Histogramas de NB clases con suavizado."""
    ft, gt = [], []
    rng = np.random.default_rng(0)
    for S, D, p, tr, lab in cal:
        PD, _H, _N = _precompute(D, p)
        for i in range(len(S)):
            for (t, b), l in zip(p[i], lab[i]):
                if t < 0:
                    continue
                k = int(np.clip(int(t / DT), 0, NT - 1))
                if l >= 1:
                    ft.append(abs(float(PD[i, k] @ polarization(tr[l - 1], S[i], b))))
                else:
                    g = grid[int(rng.integers(len(grid)))]
                    gt.append(abs(float(PD[i, k] @ polarization(g, S[i], b))))
    edges = np.linspace(0, 1, NB + 1)
    f = (np.histogram(ft, edges)[0] + 1.0)
    g = (np.histogram(gt, edges)[0] + 1.0)
    f /= f.sum()
    g /= g.sum()
    return (f / g).tolist(), len(ft), len(gt)


def fields(STA, D, picks, grid, taus, ratio, u):
    ns = len(STA)
    PD, HIT, NEAR = _precompute(D, picks)
    MF, sd = matched(D)
    DIJ = np.linalg.norm(STA[:, None, :] - STA[None, :, :], axis=2)
    st = np.arange(ns)
    ratio = np.asarray(ratio)
    out = {m: np.full(len(grid), -np.inf) for m in OWN}
    sh = int(round(0.25 / DT))
    offs = np.arange(-sh, sh + 1)
    MFn = MF / sd
    for gi, src in enumerate(grid):
        d = np.linalg.norm(STA - src, axis=1)
        pred = np.array([[polarization(src, STA[i], bb) for i in range(ns)] for bb in ("P", "S")])
        best = {m: -np.inf for m in OWN}
        for tau in taus:
            tt = np.stack([d / VP, d / VS]) + tau
            idx = np.clip((tt / DT).astype(int), 0, NT - 1)
            cl = NEAR[st[None, :], idx]
            res = np.abs(cl - tt)
            # polarizacion en el pick mas cercano (no en el tiempo predicho)
            kpk = np.clip((cl / DT).astype(int), 0, NT - 1)
            cosv = np.abs(np.einsum("bic,bic->bi", pred, PD[st[None, :], kpk]))
            rr = ratio[np.clip((cosv * NB).astype(int), 0, NB - 1)]
            for s in SIGMAS:
                nk = np.exp(-res ** 2 / (2 * s ** 2)) / (np.sqrt(2 * np.pi) * s) / u
                best[f"KIN-LR@{s}"] = max(best[f"KIN-LR@{s}"], float(np.log(EPS + (1 - EPS) * nk).sum()))
                best[f"PHYS-LR@{s}"] = max(best[f"PHYS-LR@{s}"],
                                           float(np.log(EPS + (1 - EPS) * nk * rr).sum()))
            # sin picks: filtro adaptado en +-0.25 s del tiempo predicho, vectorizado
            ii = np.clip(idx[:, :, None] + offs[None, None, :], 0, NT - 1)          # (2, ns, W)
            seg = MFn[st[None, :, None], :, ii]                                       # (2, ns, W, 3)
            proj = np.einsum("bic,biwc->biw", pred, seg) ** 2
            ener = np.sum(seg ** 2, axis=-1)
            best["WAVE-P"] = max(best["WAVE-P"], float(np.minimum(proj.max(axis=2), 36.0).sum()) / 2.0)
            best["WAVE-E"] = max(best["WAVE-E"], float(np.minimum(ener.max(axis=2) - 2.0, 36.0).sum()) / 2.0)
            hit = HIT[st[None, :], idx]
            if hit.sum() == 0:
                continue
            ok = True
            for b in range(2):
                ii = np.where(hit[b])[0]
                if len(ii) > 1 and np.any(np.abs(cl[b][ii][:, None] - cl[b][ii][None, :])
                                          > DIJ[np.ix_(ii, ii)] / VMIN + 1e-9):
                    ok = False
                    break
            s_hit = hit.sum() / (2 * ns)
            rms = float(np.sqrt(np.mean(res[hit] ** 2)))
            best["TOL-TB"] = max(best["TOL-TB"], s_hit - 0.01 * rms / TOL)
        for m in OWN:
            out[m][gi] = best[m]
    return out


def score_events(events, truth):
    hit, nf = set(), 0
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


def run_cell(args):
    scen, cond, reps, ncal, seed = args
    kind, nst, nclut = SCEN[scen]
    grid = default_grid()
    rng = np.random.default_rng(seed + 41)
    rng_d = np.random.default_rng(seed + 43)
    t0 = time.time()

    def dataset(n, base):
        out = []
        for k in range(n):
            S = make_network(kind, nst, seed=seed + base + k)
            srcs, taus = draw_sources(rng)
            D, p, tr = generate(S, srcs, taus, n_clutter=nclut, rng=rng)
            p2, lab = degrade(S, srcs, p, COND[cond], rng_d)
            out.append((S, D, p2, tr, lab))
        return out

    cal, ev = dataset(ncal, 1000), dataset(reps, 5000)
    ratio, nf_, ng_ = pol_hist(cal, grid)
    u = max(1e-3, sum(sum(1 for l in lab_i if l == 0) for (_S, _D, _p, _t, lab) in cal
                      for lab_i in lab) / (ncal * nst * float(T[-1])))
    fcal = [fields(S, D, p, grid, DEFAULT_TAUS, ratio, u) for S, D, p, _t, _l in cal]
    fev = [fields(S, D, p, grid, DEFAULT_TAUS, ratio, u) for S, D, p, _t, _l in ev]
    ntc = sum(len(c[3]) for c in cal)
    res = {"scenario": scen, "cond": cond, "seed": seed, "k_eval": [len(e[3]) for e in ev],
           "pol_ratio": ratio, "n_f": nf_, "n_g": ng_, "u": u,
           "det_frac_cal": float(sum(1 for c in cal for lab_i in c[4] for l in lab_i if l >= 1)
                                 / max(1, sum(2 * len(c[3]) * len(c[0]) for c in cal))),
           "methods": {}}
    BUDGETS = (0.05, 0.1, 0.25, 0.5)

    def cat(f, thr):
        return _catalog_impl(grid, np.where(np.isfinite(f), f, -1e9), thr)

    def curve(fl, data, thr):
        h = nf = 0
        for f, d in zip(fl, data):
            a_, b_ = score_events(cat(f, thr), d[3])
            h += a_
            nf += b_
        return h / max(1, sum(len(d[3]) for d in data)), nf / len(data)

    def thr_for_budget(fl, data, cand, budget):
        # el umbral mas bajo (recorriendo desde arriba) cuyos falsos por ventana en calibracion no
        # exceden el presupuesto
        lo = None
        for thr in cand[::-1]:
            _rc, fr = curve(fl, data, float(thr))
            if fr <= budget:
                lo = float(thr)
            else:
                break
        return lo if lo is not None else float(cand[-1]) + 1e9

    for gname, variants in GROUPS.items():
        best = None
        for v in variants:
            fl = [f[v] for f in fcal]
            vals = np.concatenate([f[np.isfinite(f)] for f in fl])
            cand = np.unique(np.quantile(vals, np.linspace(0.8, 1.0, 160)))
            thrs = [thr_for_budget(fl, cal, cand, bgt) for bgt in BUDGETS]
            rcs = [curve(fl, cal, t_)[0] for t_ in thrs]
            key = float(np.mean(rcs))
            if best is None or key > best[0]:
                best = (key, v, thrs, rcs)
        _k, v, thrs, rcs = best
        entry = {"variant": v, "budgets": list(BUDGETS), "thr": thrs, "recall_cal": rcs}
        for bgt, t_ in zip(BUDGETS, thrs):
            hits, fal, det = [], [], []
            for f, e in zip(fev, ev):
                evs = cat(f[v], t_)
                a_, b_ = score_events(evs, e[3])
                hits.append(a_)
                fal.append(b_)
                for s_ in e[3]:
                    det.append(int(any(np.linalg.norm(np.asarray(x)[:2] - s_[:2]) <= TOL_KM for x in evs)))
            entry[f"b{bgt}"] = {"hits": hits, "false": fal, "det": det}
        res["methods"][gname] = entry
    res["seconds"] = time.time() - t0
    return res


def sign_test(x, y):
    d = np.asarray(x, float) - np.asarray(y, float)
    pos, neg = int(np.sum(d > 0)), int(np.sum(d < 0))
    return pos, neg, (1.0 if pos + neg == 0 else float(binomtest(pos, pos + neg, 0.5).pvalue))


def mcnemar(a, b):
    a, b = np.asarray(a), np.asarray(b)
    n10, n01 = int(np.sum((a == 1) & (b == 0))), int(np.sum((a == 0) & (b == 1)))
    return n10, n01, (1.0 if n10 + n01 == 0 else float(binomtest(n10, n10 + n01, 0.5).pvalue))


def report(results):
    L = []
    for r in results:
        kev = r["k_eval"]
        nt, n = sum(kev), len(kev)
        L.append(f"\n=== {r['scenario']} | {r['cond']} | semilla {r['seed']} | eventos {nt}, ventanas {n}, "
                 f"deteccion de arribos {r['det_frac_cal']:.2f}")
        for m, d in r["methods"].items():
            parts = []
            for bgt in d["budgets"]:
                x = d[f"b{bgt}"]
                parts.append(f"{bgt}: rec {sum(x['hits']) / max(1, nt):.2f} fal {sum(x['false']) / n:.2f}")
            L.append(f"   {m:<8} {d['variant']:<12} " + " | ".join(parts))
        for a, b in (("KIN-LR", "PHYS-LR"), ("KIN-LR", "WAVE-P"), ("PHYS-LR", "WAVE-P"), ("WAVE-E", "WAVE-P")):
            n10, n01, p = mcnemar(r["methods"][a]["b0.1"]["det"], r["methods"][b]["b0.1"]["det"])
            L.append(f"      presupuesto 0.1, eventos detectados solo por {a} / solo por {b}: {n10}/{n01}  McNemar p={p:.4f}")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--reps", type=int, default=60)
    ap.add_argument("--ncal", type=int, default=60)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--scenarios", default=",".join(SCEN))
    ap.add_argument("--conds", default=",".join(COND))
    args = ap.parse_args()
    cells = [(s, c, args.reps, args.ncal, args.seed) for c in args.conds.split(",") for s in args.scenarios.split(",")]
    results = []
    with ProcessPoolExecutor(args.workers) as ex:
        for r in ex.map(run_cell, cells):
            results.append(r)
            print(f"   hecho {r['scenario']}/{r['cond']} en {r['seconds']:.0f} s", flush=True)
    out = f"task1o_seed{args.seed}.json"
    with open(out, "w") as fh:
        json.dump(results, fh)
    print(report(results))
    print(f"\nescrito {out}")


if __name__ == "__main__":
    main()
