"""
task1y_clutter.py - Tarea 1y: la comparacion de 1x con clutter consistente.

La verificacion de 1x hallo que en todos los generadores desde 1o los picks de clutter se entregaban
siempre, con tiempo exacto, aunque sus pulsos estan por debajo del ruido (amplitud sobre el ruido
mediana 0.52): un picker real los detectaria con probabilidad 0.003. GaMMA recibia alrededor de 90 por
ciento de picks de clutter indetectable. Aqui el clutter se trata con la misma regla que los arribos
verdaderos (deteccion con probabilidad Phi((a-2)/0.5), error de tiempo 3 veces la cota de Cramer-Rao),
en dos regimenes consistentes:
    thin    clutter debil (amplitud 0.35, la de 1o a 1x): casi nunca se detecta
    strong  clutter fuerte (amplitud 1.75, a ~ 2.6): transitorios visibles en la forma de onda y
            detectados, el caso de falsos picks por ruido cultural o espigas
Condicion pert de 1x (retardos N(0, 0.15 s), ~10 grados, factor de sitio lognormal sigma 0.5).
Metodos: WAVE-PT y WAVE-MT con supresion de vecinos a 5 y a 20 km, MT-FW, MT-FW-NODE (Frank-Wolfe con
el subespacio de otro nodo, control del metodo completo) y GaMMA con semilla y oraculo.

    python task1y_clutter.py --seed 2024 --reps 100 --clutter thin
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.stats import norm

import task1x_mtglrt as X
from task0_common import AP, AS, DEFAULT_TAUS, NT, SNR, T, VP, VS, default_grid, make_network, ricker
from task1o_consistent import _GD2, NOISE_SD, draw_sources, score_events
from task1r_baselines import GAMMA_GRID, run_gamma
from task1r_fair import CAPS, SCEN, recall_at
from task1w_mech import CP, CS, rand_dc
from froc_exact import exact_HF, greedy_record

CLUT_AMP = {"thin": 0.35, "strong": 1.75}


def generate_y(STA, sources, taus, n_clutter, rng, rng_m, rng_h, clut_amp):
    """generate_x (pert) con amplitud de clutter configurable; devuelve tambien la amplitud por pick."""
    ns = len(STA)
    D = np.zeros((ns, 3, NT))
    picks = [[] for _ in range(ns)]
    amp = [[] for _ in range(ns)]                  # amplitud en las unidades de SNR*A/d
    st = rng_h.normal(0.0, X.STATIC_SD, size=(ns, 2))
    dn = rng_h.normal(0.0, X.DIR_SD, size=(ns, 2, 3))
    site = np.exp(rng_h.normal(0.0, X.SITE_SD, size=(ns, 2)))
    for src, tau in zip(sources, taus):
        d = np.linalg.norm(STA - src, axis=1)
        M = rand_dc(rng_m)
        for i in range(ns):
            n = (STA[i] - src) / d[i]
            for bi, (b, A, V) in enumerate((("P", AP, VP), ("S", AS, VS))):
                ta = d[i] / V + tau + float(st[i, bi])
                picks[i].append((ta, b))
                if b == "P":
                    vec = n * float(n @ M @ n) / np.sqrt(CP)
                else:
                    u = M @ n
                    vec = (u - float(u @ n) * n) / np.sqrt(CS)
                a0 = np.linalg.norm(vec)
                if a0 > 0:
                    w = vec / a0 + dn[i, bi]
                    vec = a0 * site[i, bi] * w / np.linalg.norm(w)
                amp[i].append(SNR * A * float(np.linalg.norm(vec)) / d[i])
                D[i] += SNR * A * np.outer(vec, ricker(T - ta)) / d[i]
    ref = np.asarray(sources[0], float) if len(sources) else np.zeros(3)
    scale = float(np.median(np.linalg.norm(STA - ref, axis=1)))
    for _ in range(n_clutter):
        i = int(rng.integers(ns))
        tc = float(rng.uniform(2, T[-1] - 4))
        q = rng.normal(size=3)
        q /= np.linalg.norm(q)
        D[i] += SNR * clut_amp * np.outer(q, ricker(T - tc)) / scale
        picks[i].append((tc, "P" if rng.random() < 0.5 else "S"))
        amp[i].append(SNR * clut_amp / scale)
    D += rng.normal(scale=SNR * 0.02, size=D.shape)
    return D, picks, [np.asarray(s, float) for s in sources], amp


def degrade_y(picks, amp, nsrc, rng):
    """la misma regla de deteccion y de error de tiempo para arribos verdaderos y para clutter."""
    out, lab = [], []
    for pl0, am in zip(picks, amp):
        pl, ll = [], []
        for j, ((t, b), A) in enumerate(zip(pl0, am)):
            a = A / NOISE_SD
            if rng.random() > norm.cdf((a - 2.0) / 0.5):
                continue
            crb = NOISE_SD / (max(A, 1e-9) * np.sqrt(_GD2))
            pl.append((t + float(rng.normal(0, max(0.02, 3.0 * crb))), b))
            ll.append(j // 2 + 1 if j < 2 * nsrc else 0)
        if not pl:
            pl, ll = [(-1000.0, "P")], [0]
        out.append(pl)
        lab.append(ll)
    return out, lab


def run_cell(args):
    scen, reps, seed, clut, evbase, moff = args
    kind, nst, nclut = SCEN[scen]
    grid = default_grid()
    garr = np.asarray(grid)
    rng = np.random.default_rng(seed + 41 + evbase)
    rng_d = np.random.default_rng(seed + 43 + evbase)
    rng_m = np.random.default_rng(seed + moff)
    rng_h = np.random.default_rng(seed + 61 + evbase)
    rng_n = np.random.default_rng(seed + 71 + evbase)
    t0 = time.time()
    ev = []
    for k in range(reps):
        S = make_network(kind, nst, seed=seed + evbase + k)
        srcs, taus = draw_sources(rng)
        D, p, tr, amp = generate_y(S, srcs, taus, nclut, rng, rng_m, rng_h, CLUT_AMP[clut])
        p2, lab = degrade_y(p, amp, len(srcs), rng_d)
        ev.append((S, D, p2, tr, lab))

    def nodemap():
        m = rng_n.permutation(len(grid))
        for gi in range(len(grid)):
            tries = 0
            while np.linalg.norm(garr[m[gi], :2] - garr[gi, :2]) < 9.0 and tries < 50:
                m[gi] = rng_n.integers(len(grid))
                tries += 1
        return m

    R = {m: [] for m in ("WAVE-PT@5", "WAVE-PT@20", "WAVE-MT@5", "WAVE-MT@20", "MT-FW", "MT-FW-NODE")}
    for S, D, _p, tr, _l in ev:
        Gs = [X.radiation_operator(g, S) for g in garr]
        Qs = [X.basis(G) for G in Gs]
        nm = nodemap()
        f, MFn = X.fields_x(S, D, grid, DEFAULT_TAUS, nm, Qs)
        for r_ in (5, 20):
            R[f"WAVE-PT@{r_}"].append(greedy_record(grid, f["WAVE-PT"], tr, min_sep_km=float(r_)))
            R[f"WAVE-MT@{r_}"].append(greedy_record(grid, f["WAVE-MT"], tr, min_sep_km=float(r_)))
        R["MT-FW"].append(X.fw_record(S, MFn, grid, DEFAULT_TAUS, Gs, Qs, tr))
        R["MT-FW-NODE"].append(X.fw_record(S, MFn, grid, DEFAULT_TAUS, [Gs[j] for j in nm], [Qs[j] for j in nm], tr))
    nte = sum(len(e[3]) for e in ev)
    res = {"scenario": scen, "seed": seed, "clutter": clut, "k_eval": [len(e[3]) for e in ev], "methods": {}}
    for m, recs in R.items():
        _t, H, F = exact_HF(recs)
        res["methods"][m] = {"H": H.tolist(), "F": F.tolist(),
                             "recall_at": {str(c): recall_at(H, F, nte, c) for c in CAPS}}
    pts, per = [], {}
    for ci, cfg in enumerate(GAMMA_GRID):
        Hs, Fs = [], []
        for S, _D, p, tr, _l in ev:
            a, b = score_events(run_gamma(S, [[(t, bb) for (t, bb) in pl if t > -100] for pl in p], **cfg), tr)
            Hs.append(a)
            Fs.append(b)
        pts.append((sum(Fs) / len(ev), sum(Hs) / max(1, nte)))
        per[str(ci)] = (Hs, Fs)
    ra = {}
    for c in CAPS:
        ok = [(rc, ci) for ci, (fr, rc) in enumerate(pts) if fr <= c]
        ra[str(c)] = max(ok)[0] if ok else 0.0
    res["methods"]["GAMMA"] = {"points": pts, "per": per, "recall_at": ra}
    ntrue_arr = sum(2 * len(e[3]) * nst for e in ev)
    res["detected_fraction"] = float(sum(l >= 1 for e in ev for lab_i in e[4] for l in lab_i) / max(1, ntrue_arr))
    res["clutter_picks_per_window"] = float(np.mean([sum(l == 0 for lab_i in e[4] for l in lab_i) for e in ev]))
    res["true_picks_per_window"] = float(np.mean([sum(l >= 1 for lab_i in e[4] for l in lab_i) for e in ev]))
    res["seconds"] = time.time() - t0
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--clutter", default="thin", choices=["thin", "strong"])
    ap.add_argument("--evbase", type=int, default=25000)
    ap.add_argument("--moff", type=int, default=233)
    ap.add_argument("--workers", type=int, default=1)
    args = ap.parse_args()
    cells = [(s, args.reps, args.seed, args.clutter, args.evbase, args.moff) for s in SCEN]
    out = []
    with ProcessPoolExecutor(args.workers) as ex:
        for r in ex.map(run_cell, cells):
            out.append(r)
            print(f"   hecho {r['scenario']}/{r['clutter']} en {r['seconds']:.0f} s (clutter/true picks "
                  f"{r['clutter_picks_per_window']:.1f}/{r['true_picks_per_window']:.1f}) :: " + " | ".join(
                      f"{m} " + "/".join(f"{r['methods'][m]['recall_at'][str(c)]:.2f}" for c in CAPS)
                      for m in r["methods"]), flush=True)
    with open(f"task1y_seed{args.seed}_{args.clutter}_r{args.reps}.json", "w") as fh:
        json.dump(out, fh)


if __name__ == "__main__":
    main()
