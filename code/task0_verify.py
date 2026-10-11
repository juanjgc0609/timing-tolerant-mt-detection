"""
task0_verify.py — suite de verificación de implementación. EJECUTAR ANTES DE task0_run.py.

Cinco niveles, con criterios de aceptación fijados AQUI y no después de ver los números.
El propósito es separar dos preguntas que no deben mezclarse:

    (a) ¿la implementación reproduce la propiedad matemática que pretende imponer?
    (b) ¿esa propiedad mejora el problema?

Solo si (a) pasa tiene sentido interpretar (b). Si algún nivel FALLA, el resultado de
task0_run.py es sobre el código, no sobre la formulación, y no debe interpretarse.

    python task0_verify.py
    python task0_verify.py --quick     # menos réplicas, para iterar
"""
import argparse
import sys
import numpy as np

import task0_common as tc
from task0_common import (make_network, generate, field_own, catalog_from_score_field,
                          evaluate_catalog, calibrate_threshold, default_grid,
                          DEFAULT_TAUS, VP, VS, VMIN, polarization)

SRC = np.array([6.0, -8.0, 9.0])      # deliberadamente fuera de los nodos de la malla
TAU = [6.0]
PASS, FAIL = "PASS", "FAIL"
results = []


def check(level, name, ok, detail=""):
    results.append((level, name, ok))
    print(f"   [{PASS if ok else FAIL}] {name}" + (f"  — {detail}" if detail else ""))
    return ok


# ====================================================================== NIVEL 1
def level1_oracle(reps=5):
    """
    Sin clutter, sin ruido, picks exactos, modelo idéntico al generador.
    Criterio: recall = 1.00 y falsos = 0.00 para TOL y PHYS, usando SOLO el máximo del campo.
    Si esto falla, la implementación está rota y no hay nada más que discutir.
    """
    print("\nNIVEL 1 — Oracle test (sin clutter, sin ruido)")
    STA = make_network("random", 8, seed=11)
    grid = default_grid()
    rng = np.random.default_rng(11)
    ok_all = True
    for name, use in (("TOL", set()), ("PHYS", {"CIII", "CIV"})):
        rec = fal = 0
        for _ in range(reps):
            D, pk, tr = generate(STA, [SRC], TAU, n_clutter=0, rng=rng)
            D_clean = D  # el generador añade ruido 0.02*SNR; se acepta como piso instrumental
            f = field_own(STA, D_clean, pk, use, grid, DEFAULT_TAUS)
            ev = catalog_from_score_field(grid, f, 0.999 * f.max(), max_events=1)
            r, nf = evaluate_catalog(ev, tr)
            rec += int(r)
            fal += int(nf > 0)
        ok = (rec == reps) and (fal == 0)
        ok_all &= check(1, f"oracle {name}", ok, f"recall={rec}/{reps} falsos={fal}/{reps}")
    return ok_all


# ====================================================================== NIVEL 2
def level2_gradual(reps=6):
    """
    Degradación gradual. Criterio: el recall debe ser MONOTONO NO CRECIENTE al aumentar el
    clutter, sin caídas abruptas (>0.4 entre niveles consecutivos). Una caída abrupta con
    picks perfectos indica una restricción mal formulada, no dificultad del problema.
    """
    print("\nNIVEL 2 — Degradación gradual con el clutter")
    STA = make_network("random", 8, seed=12)
    grid = default_grid()
    levels = (0, 1, 5, 10, 20, 40)
    ok_all = True
    for name, use in (("TOL", set()), ("PHYS", {"CIII", "CIV"})):
        rng = np.random.default_rng(12)
        rr = []
        for nc in levels:
            rec = 0
            for _ in range(reps):
                D, pk, tr = generate(STA, [SRC], TAU, n_clutter=nc, rng=rng)
                f = field_own(STA, D, pk, use, grid, DEFAULT_TAUS)
                ev = catalog_from_score_field(grid, f, 0.999 * f.max(), max_events=1)
                rec += int(evaluate_catalog(ev, tr)[0])
            rr.append(rec / reps)
        drops = [rr[i] - rr[i + 1] for i in range(len(rr) - 1)]
        ok = max(drops) <= 0.40
        ok_all &= check(2, f"sin caída abrupta {name}", ok,
                        "recall=" + " ".join(f"{v:.2f}" for v in rr))
    return ok_all


# ====================================================================== NIVEL 3
def level3_ablation(reps=8):
    """
    Ablación por restricción. Criterio: cada condición debe ser EJECUTABLE y producir recall
    > 0. No se exige que una gane a otra; se exige que ninguna esté inactiva o rota.
    """
    print("\nNIVEL 3 — Ablación por restricción")
    STA = make_network("onesided", 8, seed=13)
    grid = default_grid()
    conds = [("TOL", set()), ("CIII", {"CIII"}), ("CIV", {"CIV"}),
             ("CIII+CIV", {"CIII", "CIV"})]
    ok_all = True
    for name, use in conds:
        rng = np.random.default_rng(13)
        rec = fal = 0
        for _ in range(reps):
            D, pk, tr = generate(STA, [SRC], TAU, n_clutter=20, rng=rng)
            f = field_own(STA, D, pk, use, grid, DEFAULT_TAUS)
            ev = catalog_from_score_field(grid, f, 0.999 * f.max(), max_events=1)
            r, nf = evaluate_catalog(ev, tr)
            rec += int(r)
            fal += int(nf > 0)
        ok = rec > 0
        ok_all &= check(3, f"condición activa {name}", ok,
                        f"recall={rec}/{reps} falsos={fal}/{reps}")
    return ok_all


# ====================================================================== NIVEL 4
def level4_invariance():
    """
    Invariancias sobre DATO DETERMINISTA (sin clutter ni ruido) y evaluadas EN EL PUNTO
    CORRESPONDIENTE. Con clutter el test es inválido: las direcciones aleatorias se sortean con
    la misma semilla pero no se rotan, así que representan direcciones físicas distintas.
    Evaluadas en el punto correspondiente, no como máximo sobre una malla:
    rotar o trasladar la malla remuestrea los nodos y produce diferencias que son de la malla,
    no de la física. La tolerancia se fija por la cuantización del índice temporal, int(t/DT),
    que es del orden de 1e-5 relativo, no 1e-9.
      4a traslación  4b rotación  4c escala cinemática  4d P no se mezcla con S
    """
    print("\nNIVEL 4 — Invariancias físicas (evaluadas en el punto correspondiente)")
    RTOL = 1e-5
    STA = make_network("random", 8, seed=14)
    use = {"CIII", "CIV"}
    rng = np.random.default_rng(14)
    D, pk, _tr = generate(STA, [SRC], TAU, n_clutter=0, rng=rng, noise_scale=0.0)
    base = float(field_own(STA, D, pk, use, [SRC], DEFAULT_TAUS)[0])
    ok_all = True

    a = np.array([13.0, -7.0, 0.0])
    rng2 = np.random.default_rng(14)
    D2, pk2, _ = generate(STA + a, [SRC + a], TAU, n_clutter=0, rng=rng2, noise_scale=0.0)
    v2 = float(field_own(STA + a, D2, pk2, use, [SRC + a], DEFAULT_TAUS)[0])
    ok_all &= check(4, "invariancia por traslación",
                    abs(v2 - base) <= RTOL * max(1.0, abs(base)),
                    f"{base:.10f} vs {v2:.10f}  (rtol {RTOL:g})")

    th = 0.7
    R = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    rng3 = np.random.default_rng(14)
    D3, pk3, _ = generate(STA @ R.T, [R @ SRC], TAU, n_clutter=0, rng=rng3, noise_scale=0.0)
    v3 = float(field_own(STA @ R.T, D3, pk3, use, [R @ SRC], DEFAULT_TAUS)[0])
    ok_all &= check(4, "invariancia por rotación",
                    abs(v3 - base) <= RTOL * max(1.0, abs(base)),
                    f"{base:.10f} vs {v3:.10f}  (rtol {RTOL:g})")

    s = 2.0
    tprev = np.linalg.norm(STA[0] - SRC) / VP
    tnew = np.linalg.norm(s * STA[0] - s * SRC) / (s * VP)
    ok_all &= check(4, "invariancia de escala cinemática",
                    abs(tnew - tprev) <= 1e-12, f"{tprev:.12f} vs {tnew:.12f}")

    pP = polarization(SRC, STA[0], "P")
    pS = polarization(SRC, STA[0], "S")
    ok_all &= check(4, "P y S ortogonales (no se mezclan fases)",
                    abs(float(np.dot(pP, pS))) <= 1e-12,
                    f"|P·S| = {abs(float(np.dot(pP, pS))):.2e}")
    return ok_all


# ====================================================================== NIVEL 5
def level5_independent(n_cases=2000, seed=15):
    """
    Implementación independiente de C-III. La versión vectorizada del score se compara contra
    una implementación explícita con dobles bucles sobre miles de casos aleatorios.
    Criterio: acuerdo exacto.
    """
    print("\nNIVEL 5 — Implementación independiente de C-III")
    rng = np.random.default_rng(seed)
    STA = make_network("random", 8, seed=seed)
    DIJ = np.linalg.norm(STA[:, None, :] - STA[None, :, :], axis=2)
    disagree = 0
    for _ in range(n_cases):
        k = int(rng.integers(2, 9))
        idx = rng.choice(8, k, replace=False)
        t = rng.uniform(0, 30, k)
        # A: vectorizada
        a_ok = not np.any(np.abs(t[:, None] - t[None, :]) >
                          DIJ[np.ix_(idx, idx)] / VMIN + 1e-9)
        # B: explícita, lenta
        b_ok = True
        for i in range(k):
            for j in range(i + 1, k):
                if abs(t[i] - t[j]) > DIJ[idx[i], idx[j]] / VMIN + 1e-9:
                    b_ok = False
        if a_ok != b_ok:
            disagree += 1
    return check(5, "C-III vectorizada == C-III explícita", disagree == 0,
                 f"{n_cases - disagree}/{n_cases} de acuerdo")


# ====================================================================== extra
def level0_calibration(reps=10):
    """
    Chequeo de la calibración, que fue la fuente del último artefacto: el umbral calibrado
    debe dar simultáneamente recall alto Y falsos bajos. Optimizar solo el recall produce
    umbrales que dejan pasar medio mapa.
    """
    print("\nNIVEL 0 — Sanidad de la calibración")
    STA = make_network("random", 12, seed=16)
    grid = default_grid()
    rng = np.random.default_rng(16)
    data = [generate(STA, [SRC], TAU, n_clutter=35, rng=rng) for _ in range(reps)]
    truth = [t for (_D, _p, t) in data]
    ok_all = True
    for name, use in (("TOL", set()), ("PHYS", {"CIII", "CIV"})):
        fs = [field_own(STA, D, p, use, grid, DEFAULT_TAUS) for (D, p, _t) in data]
        thr, rc, reach = calibrate_threshold(fs, truth, grid, target_recall=0.98)
        r, f = tc._recall_false_at(fs, truth, grid, thr, 5.0, 5, 5.0)
        ok = reach and r >= 0.98 and f <= 0.20
        ok_all &= check(0, f"calibración sana {name}", ok,
                        f"umbral={thr:.3f} recall={r:.2f} falsos={f:.2f}")
    return ok_all


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    q = args.quick
    print("=" * 72)
    print("SUITE DE VERIFICACION DE IMPLEMENTACION")
    print("Criterios fijados por adelantado. Si algo FALLA, los resultados de task0_run.py")
    print("son sobre el codigo y no sobre la formulacion fisica.")
    print("=" * 72)

    level0_calibration(reps=6 if q else 10)
    level1_oracle(reps=3 if q else 5)
    level2_gradual(reps=3 if q else 6)
    level3_ablation(reps=4 if q else 8)
    level4_invariance()
    level5_independent(n_cases=400 if q else 2000)

    n_fail = sum(1 for _l, _n, ok in results if not ok)
    print("\n" + "=" * 72)
    if n_fail == 0:
        print("TODOS LOS NIVELES PASAN. Los resultados de task0_run.py pueden interpretarse")
        print("como evidencia sobre la formulacion fisica, no sobre la implementacion.")
    else:
        print(f"{n_fail} COMPROBACION(ES) FALLAN. NO ejecutes la Tarea 0 completa todavia:")
        for l, n, ok in results:
            if not ok:
                print(f"   nivel {l}: {n}")
    print("=" * 72)
    sys.exit(0 if n_fail == 0 else 1)


if __name__ == "__main__":
    main()
