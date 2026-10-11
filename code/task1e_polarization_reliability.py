"""
task1e_polarization_reliability.py — Tarea 1e: cuanto de la ventaja de G1 sobrevive con la
confiabilidad de polarizacion medida en la Tarea 1, y si un gate por estacion la recupera.

La Tarea 1 encontro que, en la red CM, la polarizacion P de una estacion o conserva la direccion del
evento (ZAR) o esta al nivel del azar (CHI, RUS). Los pre-estudios y la Tarea 0 suponen polarizacion
fiel en TODAS las estaciones. Este experimento repite el protocolo final de la Tarea 0 (Apendice X:
red sorteada por realizacion, umbral relativo al maximo de cada realizacion, calibracion a recall
0.98, McNemar exacto pareado) cambiando una sola cosa: la polarizacion de los arribos VERDADEROS en
cada estacion, segun su perfil.

Perfiles de estacion
  IDEAL  polarizacion exacta (el supuesto de los pre-estudios).
  BUENA  azimut conservado con error gaussiano de SIGMA_BUENA grados (tipo ZAR).
  MALA   azimut sorteado uniforme, incidencia conservada (tipo CHI: azimut al nivel del azar).
  MALA3D direccion sorteada uniforme en la esfera (tipo CHI si tambien la incidencia es ruido).

Metodos, todos sobre el mismo campo y las mismas realizaciones
  TOL     sin restricciones (la ablacion ABL de la propuesta).
  CIII    solo la restriccion cinematica C-III: el piso si la polarizacion no sirve en ningun lado.
  PHYS    C-III + C-IV en todas las estaciones (lo que asumen los pre-estudios).
  PHYS-G  C-III + C-IV solo en las estaciones que pasan el gate de la Tarea 1; los arribos de las
          demas no entran al factor de polarizacion, y sin ninguna estacion en el gate el factor
          vale 1. Es la forma ingenua de "apagar C-IV donde no sirve".
  PHYS-N  C-III + C-IV donde cada arribo de una estacion fuera del gate aporta un valor NEUTRO
          c_neutro en vez de su coseno medido. c_neutro es el coseno medio de los arribos
          verdaderos en las estaciones del gate, medido sobre las realizaciones de CALIBRACION:
          lo mismo que el control de la Tarea 1 mide con eventos catalogados. Un arribo sin
          informacion de polarizacion no premia ni castiga a la hipotesis que lo explica.

La perturbacion se sortea con un generador aparte, asi que para una misma semilla todas las
condiciones comparten redes, clutter y ruido: las comparaciones entre condiciones tambien son
pareadas.

    python task1e_polarization_reliability.py --smoke
    python task1e_polarization_reliability.py --reps 30 --ncal 30 --workers 2
"""
import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.stats import beta, binomtest

from task0_common import (AP, AS, DEFAULT_TAUS, DT, NT, SNR, T, VMIN, VP, VS,
                          _precompute, calibrate_threshold, catalog_from_score_field,
                          default_grid, evaluate_catalog, make_network, polarization, ricker,
                          score_grid, validate_positive_control)

NI_MARGIN = 0.02
TARGET = 1.0 - NI_MARGIN
TOL_KM = 5.0
SIGMA_BUENA = 10.0
S1 = np.array([6.0, -8.0, 9.0])

SCENARIOS = {
    "unilateral_c40": ("unilateral, clutter 40", "onesided", 8, 40),
    "ladder_3": ("3 estaciones, clutter 35", "random", 3, 35),
    "ladder_4": ("4 estaciones, clutter 35", "random", 4, 35),
}
# fraccion de estaciones buenas y perfil de las malas
CONDITIONS = {
    "ideal": (None, None),
    "todas_buenas": (1.0, "MALA"),
    "mitad_buenas": (0.5, "MALA"),
    "cuarto_buenas": (0.25, "MALA"),
    "ninguna_buena": (0.0, "MALA"),
    "cuarto_buenas_3d": (0.25, "MALA3D"),
    "ninguna_buena_3d": (0.0, "MALA3D"),
}
METHODS = ("TOL", "CIII", "PHYS", "PHYS-G", "PHYS-N")


def rot_z(v, deg):
    a = np.radians(deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([c * v[0] - s * v[1], s * v[0] + c * v[1], v[2]])


def perturb(pol, kind, rng_pol):
    if kind == "IDEAL":
        return pol
    if kind == "BUENA":
        return rot_z(pol, rng_pol.normal(0.0, SIGMA_BUENA))
    if kind == "MALA":
        return rot_z(pol, rng_pol.uniform(0.0, 360.0))
    if kind == "MALA3D":
        q = rng_pol.normal(size=3)
        return q / np.linalg.norm(q)
    raise ValueError(kind)


def generate_rel(STA, sources, taus, n_clutter, rng, kinds, rng_pol, noise_scale=0.02):
    """
    Identico a task0_common.generate (sin senuelo) salvo que la polarizacion de cada arribo
    verdadero en la estacion i pasa por perturb(kinds[i]). Con todos los perfiles IDEAL el
    resultado es identico bit a bit al de generate con la misma semilla (verificado en --selftest).
    """
    ns = len(STA)
    D = np.zeros((ns, 3, NT))
    picks = [[] for _ in range(ns)]
    for src, tau in zip(sources, taus):
        d = np.linalg.norm(STA - src, axis=1)
        for i in range(ns):
            for b, A, V in (("P", AP, VP), ("S", AS, VS)):
                ta = d[i] / V + tau
                picks[i].append((ta, b))
                pol = perturb(polarization(src, STA[i], b), kinds[i], rng_pol)
                D[i] += SNR * A * np.outer(pol, ricker(T - ta)) / d[i]
    ref = np.asarray(sources[0], float)
    scale = float(np.median(np.linalg.norm(STA - ref, axis=1)))
    for _ in range(n_clutter):
        i = int(rng.integers(ns))
        tc = float(rng.uniform(2, T[-1] - 4))
        q = rng.normal(size=3)
        q /= np.linalg.norm(q)
        D[i] += SNR * 0.35 * np.outer(q, ricker(T - tc)) / scale
        picks[i].append((tc, "P" if rng.random() < 0.5 else "S"))
    if noise_scale > 0:
        D += rng.normal(scale=SNR * noise_scale, size=D.shape)
    return D, picks, [np.asarray(s, float) for s in sources]


def neutral_cosine(cal):
    """
    Coseno medio |p_pred . p_obs| de los arribos VERDADEROS (P y S) en estaciones del gate, sobre
    las realizaciones de calibracion. Se evalua en la fuente verdadera y en el tiempo del arribo,
    que es lo que un operador obtiene con eventos catalogados en esas estaciones.
    """
    vals = []
    for S, D, p, tr, gate in cal:
        PD, _H, _N = _precompute(D, p)
        src = tr[0]
        for i in range(len(S)):
            if not gate[i]:
                continue
            for bb in ("P", "S"):
                ta = [t for (t, b) in p[i] if b == bb][0]
                k = int(np.clip(int(ta / DT), 0, NT - 1))
                vals.append(abs(float(np.dot(polarization(src, S[i], bb), PD[i, k]))))
    return float(np.mean(vals)) if vals else 1.0


def fields_all(STA, D, picks, grid, taus, gate, c_neutral=1.0):
    """
    Los cuatro campos en una pasada. TOL y PHYS replican task0_common.score_grid linea a linea
    (verificado en --selftest); CIII es PHYS sin el factor de polarizacion; PHYS-G calcula ese
    factor solo con las estaciones gate=True, y sin ninguna vale 1 (sin informacion no hay
    penalizacion), con lo que se reduce a CIII.
    """
    ns = len(STA)
    PD, HIT, NEAR = _precompute(D, picks)
    DIJ = np.linalg.norm(STA[:, None, :] - STA[None, :, :], axis=2)
    st = np.arange(ns)
    out = {m: np.zeros(len(grid)) for m in METHODS}
    for gi, src in enumerate(grid):
        d = np.linalg.norm(STA - src, axis=1)
        pred = [[polarization(src, STA[i], bb) for i in range(ns)] for bb in ("P", "S")]
        best = dict.fromkeys(METHODS, 0.0)
        for tau in taus:
            tt = np.stack([d / VP, d / VS]) + tau
            idx = np.clip((tt / DT).astype(int), 0, NT - 1)
            hit = HIT[st[None, :], idx]
            if hit.sum() == 0:
                continue
            s = hit.sum() / (2 * ns)
            best["TOL"] = max(best["TOL"], s)
            cl = NEAR[st[None, :], idx]
            bad = False
            for b in range(2):
                ii = np.where(hit[b])[0]
                if len(ii) > 1 and np.any(
                        np.abs(cl[b][ii][:, None] - cl[b][ii][None, :])
                        > DIJ[np.ix_(ii, ii)] / VMIN + 1e-9):
                    bad = True
                    break
            if bad:
                continue
            best["CIII"] = max(best["CIII"], s)
            cs, cs_g, cs_n = [], [], []
            for b in range(2):
                for i in range(ns):
                    if hit[b, i]:
                        c = abs(float(np.dot(pred[b][i], PD[i, idx[b, i]])))
                        cs.append(c)
                        if gate[i]:
                            cs_g.append(c)
                        cs_n.append(c if gate[i] else c_neutral)
            best["PHYS"] = max(best["PHYS"], s * (float(np.mean(cs)) ** 2 if cs else 0.0))
            best["PHYS-G"] = max(best["PHYS-G"], s * (float(np.mean(cs_g)) ** 2 if cs_g else 1.0))
            best["PHYS-N"] = max(best["PHYS-N"], s * (float(np.mean(cs_n)) ** 2 if cs_n else 0.0))
        for m in METHODS:
            out[m][gi] = best[m]
    return out


def station_kinds(ns, cond, rng_assign):
    frac, bad_kind = CONDITIONS[cond]
    if frac is None:
        return ["IDEAL"] * ns, [True] * ns
    n_good = int(round(frac * ns))
    if frac > 0 and n_good == 0:
        n_good = 1
    good = np.zeros(ns, bool)
    good[rng_assign.choice(ns, n_good, replace=False)] = True
    return [("BUENA" if g else bad_kind) for g in good], list(good)


def cp_ci(k, n, alpha=0.05):
    lo = beta.ppf(alpha / 2, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(1 - alpha / 2, k + 1, n - k) if k < n else 1.0
    return float(lo), float(hi)


def mcnemar(a, b):
    a, b = np.asarray(a, int), np.asarray(b, int)
    n10 = int(np.sum((a == 1) & (b == 0)))
    n01 = int(np.sum((a == 0) & (b == 1)))
    p = 1.0 if n10 + n01 == 0 else float(binomtest(n10, n10 + n01, 0.5).pvalue)
    return n10, n01, p


def run_cell(args_tuple):
    """Una celda escenario x condicion: calibracion y evaluacion, todos los metodos."""
    scen, cond, reps, ncal, seed = args_tuple
    label, kind, nst, nclut = SCENARIOS[scen]
    grid = default_grid()
    t0 = time.time()
    rng = np.random.default_rng(seed)
    rng_pol = np.random.default_rng(seed + 777_777)
    rng_assign = np.random.default_rng(seed + 555_555)

    def dataset(n, base):
        data = []
        for k in range(n):
            S = make_network(kind, nst, seed=seed + base + k)
            kinds, gate = station_kinds(nst, cond, rng_assign)
            D, p, tr = generate_rel(S, [S1], [6.0], nclut, rng, kinds, rng_pol)
            validate_positive_control(S, p, tr)
            data.append((S, D, p, tr, gate))
        return data

    cal = dataset(ncal, 1000)
    ev = dataset(reps, 5000)
    c_neu = neutral_cosine(cal)
    fcal = [fields_all(S, D, p, grid, DEFAULT_TAUS, g, c_neu) for S, D, p, _t, g in cal]
    fev = [fields_all(S, D, p, grid, DEFAULT_TAUS, g, c_neu) for S, D, p, _t, g in ev]
    res = {"scenario": scen, "label": label, "condition": cond, "reps": reps, "ncal": ncal,
           "seed": seed, "c_neutral": c_neu, "methods": {}}
    for m in METHODS:
        thr, rc, ok = calibrate_threshold([f[m] for f in fcal], [c[3] for c in cal], grid,
                                          target_recall=TARGET, tol_km=TOL_KM)
        rec, fal = [], []
        for f, (_S, _D, _p, tr, _g) in zip(fev, ev):
            cat = catalog_from_score_field(grid, f[m], thr)
            r, nf = evaluate_catalog(cat, tr, TOL_KM)
            rec.append(int(r))
            fal.append(int(nf > 0))
        res["methods"][m] = {"thr": thr, "recall_cal": rc, "reachable": bool(ok),
                             "recall": rec, "false": fal}
    res["seconds"] = time.time() - t0
    return res


def summarize(results):
    lines = []
    by_scen = {}
    for r in results:
        by_scen.setdefault(r["scenario"], []).append(r)
    for scen, rows in by_scen.items():
        lines.append(f"\n=== {SCENARIOS[scen][0]} ===")
        lines.append(f"   {'condicion':<18}{'metodo':<8}{'recall':>22}{'falsos':>22}   calib")
        order = list(CONDITIONS)
        for r in sorted(rows, key=lambda x: order.index(x["condition"])):
            for m in METHODS:
                d = r["methods"][m]
                n = len(d["recall"])
                kr, kf = int(sum(d["recall"])), int(sum(d["false"]))
                lr, hr = cp_ci(kr, n)
                lf, hf = cp_ci(kf, n)
                cal = "ok" if d["reachable"] else "NO alcanza 0.98"
                lines.append(f"   {r['condition']:<18}{m:<8}{kr:>3}/{n}={kr / n:.2f} [{lr:.2f},{hr:.2f}]"
                             f"{kf:>5}/{n}={kf / n:.2f} [{lf:.2f},{hf:.2f}]   {cal}")
            ph, tol, g, nn = (r["methods"][k] for k in ("PHYS", "TOL", "PHYS-G", "PHYS-N"))
            lines.append(f"      c_neutro={r.get('c_neutral', float('nan')):.3f}")
            for name, a, b in (("TOL vs PHYS", tol, ph), ("TOL vs PHYS-N", tol, nn),
                               ("PHYS vs PHYS-N", ph, nn), ("PHYS-G vs PHYS-N", g, nn)):
                for what in ("false", "recall"):
                    n10, n01, p = mcnemar(a[what], b[what])
                    lab = "falsos" if what == "false" else "recall"
                    lines.append(f"      {lab} {name:<17} discordantes {n10}/{n01}  McNemar p={p:.4f}")
    return "\n".join(lines)


def selftest():
    """Con perfiles IDEAL, generate_rel y fields_all deben reproducir la Tarea 0 exactamente."""
    from task0_common import generate
    grid = default_grid()
    for k in range(3):
        S = make_network("onesided", 8, seed=42 + k)
        D1, p1, _ = generate(S, [S1], [6.0], n_clutter=40, rng=np.random.default_rng(k))
        D2, p2, _ = generate_rel(S, [S1], [6.0], 40, np.random.default_rng(k), ["IDEAL"] * 8,
                                 np.random.default_rng(99))
        assert np.array_equal(D1, D2) and p1 == p2, "generate_rel IDEAL difiere de generate"
        f = fields_all(S, D2, p2, grid, DEFAULT_TAUS, [True] * 8)
        ref_tol = score_grid(S, D1, p1, set(), grid, DEFAULT_TAUS)
        ref_phys = score_grid(S, D1, p1, {"CIII", "CIV"}, grid, DEFAULT_TAUS)
        assert np.allclose(f["TOL"], ref_tol, atol=0, rtol=0), "TOL difiere de score_grid"
        assert np.allclose(f["PHYS"], ref_phys, atol=0, rtol=0), "PHYS difiere de score_grid"
        assert np.array_equal(f["PHYS-G"], f["PHYS"]), "PHYS-G con todas en gate debe ser PHYS"
        assert np.array_equal(f["PHYS-N"], f["PHYS"]), "PHYS-N con todas en gate debe ser PHYS"
        f0 = fields_all(S, D2, p2, grid, DEFAULT_TAUS, [False] * 8, c_neutral=0.8)
        assert np.array_equal(f0["PHYS-G"], f0["CIII"]), "PHYS-G sin gate debe ser CIII"
        # sin estaciones en el gate, PHYS-N es CIII escalado por una constante: mismo ordenamiento
        # y, con umbral relativo al maximo, exactamente el mismo catalogo
        assert np.allclose(f0["PHYS-N"], 0.64 * f0["CIII"]), "PHYS-N sin gate debe ser CIII escalado"
    print("selftest OK: generate_rel(IDEAL) == generate, TOL/PHYS == score_grid, PHYS-G y PHYS-N "
          "en sus dos limites")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--ncal", type=int, default=30)
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--scenarios", default=",".join(SCENARIOS))
    ap.add_argument("--conditions", default=",".join(CONDITIONS))
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", default="task1e_results.json")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    reps, ncal = (3, 4) if args.smoke else (args.reps, args.ncal)
    cells = [(s, c, reps, ncal, args.seed) for s in args.scenarios.split(",")
             for c in args.conditions.split(",")]
    print(f"[protocolo] Tarea 0 final: red por realizacion, umbral relativo, calibracion a "
          f"recall {TARGET:.2f}, {reps} evaluacion + {ncal} calibracion por celda, "
          f"sigma_buena={SIGMA_BUENA} grados, {len(cells)} celdas")
    results = []
    if args.workers > 1:
        with ProcessPoolExecutor(args.workers) as ex:
            for r in ex.map(run_cell, cells):
                results.append(r)
                print(f"   hecho {r['scenario']}/{r['condition']} en {r['seconds']:.0f} s", flush=True)
    else:
        for c in cells:
            r = run_cell(c)
            results.append(r)
            print(f"   hecho {r['scenario']}/{r['condition']} en {r['seconds']:.0f} s", flush=True)
    with open(args.out, "w") as fh:
        json.dump(results, fh)
    print(summarize(results))
    print(f"\nescrito {args.out}")


if __name__ == "__main__":
    main()
