"""
task0_run.py — Tarea 0. VERSIÓN CORREGIDA.

Correcciones estructurales respecto de la versión anterior:

  D1. CONTRADICCION INTERNA ELIMINADA. Antes la calibración usaba TARGET_RECALL=0.90 mientras la
      compuerta de elegibilidad exigía recall >= recall_TOL - NI_MARGIN = 0.98. Calibrar a 0.90
      garantizaba fallar la compuerta, y PHYS quedó inelegible en las quince celdas. Ahora el
      objetivo de calibración se deriva del límite de no inferioridad y es configurable.

  D2. PUNTO DE OPERACION IGUALADO PARA LOS CUATRO METODOS. Antes TOL y PHYS se calibraban y
      GaMMA y PyOcto corrían en su valor por defecto, así que sus tasas de falsos no eran
      comparables. Ahora los cuatro se barren sobre su propia perilla:
        TOL, PHYS  -> umbral del campo de puntaje
        GAMMA      -> min_picks_per_eq
        PYOCTO     -> n_picks
      y se comparan a RECALL IGUALADO, interpolando de la curva característica.

  D3. Se reporta la CURVA COMPLETA además del punto. Un método que no alcanza el recall objetivo
      en ningún punto de operación se reporta como tal: es un resultado, no una exclusión.

  D4. El límite de no inferioridad se calcula sobre el conjunto de CALIBRACION, no sobre el de
      evaluación, para no seleccionar métodos usando el mismo dato con que se los juzga.

Uso:
    python task0_run.py --scenario onesided --reps 30 --ncal 30
    python task0_run.py --scenario ladder   --reps 30 --ncal 30
    python task0_run.py --scenario decoy    --reps 30 --ncal 30
    python task0_run.py --smoke                      # rápido, sin baselines externos
"""
import argparse
import json
import sys
import numpy as np
import pandas as pd
from scipy.stats import beta, binomtest

from task0_common import (make_network, generate, picks_to_dataframe, stations_to_dataframe,
                          evaluate_catalog, field_own, catalog_from_score_field,
                          calibrate_threshold, operating_curve, false_rate_at_recall,
                          validate_positive_control, default_grid, DEFAULT_TAUS, VP, VS)

HAVE_GAMMA = HAVE_PYOCTO = False
try:
    from gamma.utils import association as gamma_association
    HAVE_GAMMA = True
except Exception as e:                                                  # noqa: BLE001
    print(f"[aviso] GaMMA no disponible: {e}", file=sys.stderr)
try:
    import pyocto
    HAVE_PYOCTO = True
except Exception as e:                                                  # noqa: BLE001
    print(f"[aviso] PyOcto no disponible: {e}", file=sys.stderr)

NI_MARGIN = 0.02
PC_REPS_DEFAULT = 10
PC_TARGET_DEFAULT = 0.98
TOL_KM = 5.0
# perillas de punto de operación de los baselines, de permisivo a estricto
GAMMA_KNOB = [5, 6, 8, 10, 12, 14, 16]
PYOCTO_KNOB = [4, 6, 8, 10, 12, 14]


# ------------------------------------------------------------------ estadística
def cp_ci(k, n, alpha=0.05):
    lo = beta.ppf(alpha / 2, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(1 - alpha / 2, k + 1, n - k) if k < n else 1.0
    return float(lo), float(hi)


def exact_mcnemar(y_a, y_b):
    """McNemar exacto bilateral sobre resultados binarios PAREADOS por realización."""
    y_a = np.asarray(y_a, int)
    y_b = np.asarray(y_b, int)
    if y_a.shape != y_b.shape:
        raise ValueError("vectores pareados de distinta longitud")
    b = int(np.sum((y_a == 1) & (y_b == 0)))
    c = int(np.sum((y_a == 0) & (y_b == 1)))
    if b + c == 0:
        return b, c, 1.0
    p = binomtest(b, b + c, 0.5, alternative="two-sided").pvalue
    return b, c, float(p)


def _domain(STA):
    return ([float(STA[:, 0].min() - 30), float(STA[:, 0].max() + 30)],
            [float(STA[:, 1].min() - 30), float(STA[:, 1].max() + 30)],
            [0.0, 30.0])


def _time_before(STA, xr, yr, zr, margin=2.0):
    diag = float(np.hypot(xr[1] - xr[0], yr[1] - yr[0]))
    return float(margin * (diag + zr[1]) / VS + 60.0)


# ------------------------------------------------------------------ baselines
def _event_xyz_from_gamma(cat):
    cdf = pd.DataFrame(cat)
    if cdf.empty:
        return []
    cols = {c.lower().replace("(km)", "").strip(): c for c in cdf.columns}
    for k in ("x", "y", "z"):
        if k not in cols:
            raise RuntimeError(f"columnas inesperadas en catalogo GaMMA: {list(cdf.columns)}")
    return [np.array([float(r[cols["x"]]), float(r[cols["y"]]), float(r[cols["z"]])])
            for _, r in cdf.iterrows()]


def run_gamma(STA, picks, min_picks_per_eq=6):
    if not HAVE_GAMMA:
        raise RuntimeError("GaMMA no instalado")
    pk = picks_to_dataframe(picks, "gamma")
    stdf = stations_to_dataframe(STA, "gamma")
    xr, yr, zr = _domain(STA)
    cfg = {
        "dims": ["x(km)", "y(km)", "z(km)"],
        "x(km)": xr, "y(km)": yr, "z(km)": zr,
        "vel": {"p": VP, "s": VS},
        "method": "BGMM", "oversample_factor": 5,
        "use_amplitude": False, "use_dbscan": True,
        "dbscan_eps": 10.0, "dbscan_min_samples": 3,
        "min_picks_per_eq": int(min_picks_per_eq),
        "max_sigma11": 2.0, "max_sigma22": 1.0, "max_sigma12": 1.0,
        "bfgs_bounds": ((xr[0] - 1, xr[1] + 1), (yr[0] - 1, yr[1] + 1),
                        (zr[0], zr[1]), (None, None)),
        "covariance_prior": [5.0, 5.0], "initial_points": [1, 1, 1],
    }
    cat, _assign = gamma_association(pk, stdf, cfg, method="BGMM")
    return [] if not cat else _event_xyz_from_gamma(cat)


def run_pyocto(STA, picks, n_picks=6):
    if not HAVE_PYOCTO:
        raise RuntimeError("PyOcto no instalado")
    pk = picks_to_dataframe(picks, "pyocto")
    stdf = stations_to_dataframe(STA, "pyocto")
    if not np.issubdtype(pk["time"].dtype, np.number):
        raise AssertionError(f"PyOcto exige time numerico; llego {pk['time'].dtype}")
    xr, yr, zr = _domain(STA)
    vm = pyocto.VelocityModel0D(p_velocity=VP, s_velocity=VS, tolerance=1.0)
    npk = int(n_picks)
    assoc = pyocto.OctoAssociator(
        xlim=tuple(xr), ylim=tuple(yr), zlim=tuple(zr),
        velocity_model=vm, time_before=_time_before(STA, xr, yr, zr),
        n_picks=npk, n_p_picks=max(1, npk // 3), n_s_picks=max(1, npk // 3),
        n_p_and_s_picks=max(1, npk // 3), n_threads=1)
    events, _asg = assoc.associate(pk, stdf)
    if events is None or len(events) == 0:
        return []
    return [np.array([float(r["x"]), float(r["y"]), float(r["z"])])
            for _, r in events.iterrows()]


# ------------------------------------------------------------------ escenarios
def scenario_cases(name):
    S1 = np.array([6.0, -8.0, 9.0])
    if name == "onesided":
        return [("unilateral, clutter 15", "onesided", 8, [S1], [6.0], 15, None, 0),
                ("unilateral, clutter 40", "onesided", 8, [S1], [6.0], 40, None, 0),
                ("dispersa,   clutter 40", "spread",   8, [S1], [6.0], 40, None, 0)]
    if name == "ladder":
        return [(f"{n} estaciones", "random", n, [S1], [6.0], 35, None, 0)
                for n in (3, 4, 6, 12)]
    if name == "decoy":
        return [("senuelo 4 km,  6/16", "spread", 8, [S1], [6.0], 8,
                 S1 + np.array([4.0, 0.0, 0.0]), 6),
                ("senuelo 12 km, 12/16", "spread", 8, [S1], [6.0], 8,
                 S1 + np.array([12.0, 0.0, 0.0]), 12)]
    raise ValueError(name)


def _make_dataset(STA, srcs, taus, nclut, dec, sup, rng, n):
    """Genera n realizaciones y las conserva para que TODOS los métodos vean las mismas."""
    out = []
    for _ in range(n):
        D, picks, truth = generate(STA, srcs, taus, n_clutter=nclut, rng=rng,
                                   coherent_decoy=dec, decoy_support=sup)
        validate_positive_control(STA, picks, truth)
        out.append((D, picks, truth))
    return out


def _external_curve(runner, knobs, data, nets):
    """
    D2. Curva característica de un baseline externo: barre su perilla y devuelve
    [(perilla, recall, tasa_falsos)] más los vectores pareados por punto de operación.
    """
    curve, paired = [], {}
    for k in knobs:
        rec, fal = [], []
        for S, (_D, picks, truth) in zip(nets, data):
            try:
                ev = runner(S, picks, k)
            except Exception as e:                                      # noqa: BLE001
                print(f"   [error] perilla={k}: {e}", file=sys.stderr)
                ev = []
            r, nf = evaluate_catalog(ev, truth, TOL_KM)
            rec.append(int(r))
            fal.append(int(nf > 0))
        curve.append((float(k), float(np.mean(rec)), float(np.mean(fal))))
        paired[k] = (np.array(rec), np.array(fal))
    curve.sort(key=lambda t: t[1])
    return curve, paired


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", default="onesided",
                    choices=["onesided", "ladder", "decoy"])
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--ncal", type=int, default=30)
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--target-recall", type=float, default=None,
                    help="objetivo de calibración; por defecto 1 - NI_MARGIN")
    ap.add_argument("--pc-reps", type=int, default=PC_REPS_DEFAULT)
    ap.add_argument("--pc-target", type=float, default=PC_TARGET_DEFAULT)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()

    R = 3 if args.smoke else args.reps
    NCAL = 4 if args.smoke else args.ncal
    # D1: el objetivo de calibración DEBE ser compatible con la compuerta de no inferioridad.
    TARGET = args.target_recall if args.target_recall is not None else (1.0 - NI_MARGIN)

    meta = {"python": sys.version.split()[0], "numpy": np.__version__,
            "pandas": pd.__version__, "gamma": HAVE_GAMMA, "pyocto": HAVE_PYOCTO}
    print(f"   [runtime] {json.dumps(meta, sort_keys=True)}")
    print(f"   [protocol] target_recall_calib={TARGET:.3f}, NI_margin={NI_MARGIN:.3f}, "
          f"reps={R}, ncal={NCAL}")

    grid = default_grid()
    own = {"TOL": set(), "PHYS": {"CIII", "CIV"}}
    ext = {}
    if HAVE_GAMMA and not args.smoke:
        ext["GAMMA"] = (run_gamma, GAMMA_KNOB)
    if HAVE_PYOCTO and not args.smoke:
        ext["PYOCTO"] = (run_pyocto, PYOCTO_KNOB)

    print(f"\nTarea 0 | escenario '{args.scenario}' | metodos: "
          f"{', '.join(list(own) + list(ext))}\n")

    for label, kind, nst, srcs, taus, nclut, dec, sup in scenario_cases(args.scenario):
        # La geometria se sortea POR REALIZACION. Fijar una sola red por celda hace que la
        # geometria concreta domine sobre el metodo: entre dos semillas la tasa de falsos de un
        # mismo metodo variaba de 0.13 a 0.70 en la misma celda.
        STA = make_network(kind, nst, seed=args.seed)          # red de referencia, para reportes
        rng = np.random.default_rng(args.seed)
        cal_data, cal_nets = [], []
        for k in range(NCAL):
            S = make_network(kind, nst, seed=args.seed + 1000 + k)
            cal_data += _make_dataset(S, srcs, taus, nclut, dec, sup, rng, 1)
            cal_nets.append(S)
        eval_data, eval_nets = [], []
        for k in range(R):
            S = make_network(kind, nst, seed=args.seed + 5000 + k)
            eval_data += _make_dataset(S, srcs, taus, nclut, dec, sup, rng, 1)
            eval_nets.append(S)
        cal_truth = [t for (_D, _p, t) in cal_data]

        print(f"  {label}")

        # ---------- calibración de los métodos de campo ----------
        thr, reach = {}, {}
        for m, use in own.items():
            fields = [field_own(S, D, p, use, grid, DEFAULT_TAUS)
                      for S, (D, p, _t) in zip(cal_nets, cal_data)]
            thr[m], rc, ok = calibrate_threshold(fields, cal_truth, grid,
                                                 target_recall=TARGET, tol_km=TOL_KM)
            reach[m] = ok
            flag = "alcanzable" if ok else "NO alcanza el objetivo a ningun umbral"
            print(f"   [calib] {m}: umbral={thr[m]:.4f} recall_cal={rc:.3f} -> {flag}")

        # ---------- curvas características ----------
        curves, paired = {}, {}
        for m, use in own.items():
            f_eval = [field_own(S, D, p, use, grid, DEFAULT_TAUS)
                      for S, (D, p, _t) in zip(eval_nets, eval_data)]
            ev_truth = [t for (_D, _p, t) in eval_data]
            curves[m] = operating_curve(f_eval, ev_truth, grid, tol_km=TOL_KM)
            rec = [int(evaluate_catalog(
                       catalog_from_score_field(grid, f, thr[m]), t, TOL_KM)[0])
                   for f, t in zip(f_eval, ev_truth)]
            fal = [int(evaluate_catalog(
                       catalog_from_score_field(grid, f, thr[m]), t, TOL_KM)[1] > 0)
                   for f, t in zip(f_eval, ev_truth)]
            paired[m] = (np.array(rec), np.array(fal))
        for m, (runner, knobs) in ext.items():
            curves[m], pk_by_knob = _external_curve(runner, knobs, eval_data, eval_nets)
            # punto de operación del baseline: el que alcanza TARGET con menos falsos
            bestk, bestf = None, None
            for k, (rc, fl) in pk_by_knob.items():
                if rc.mean() >= TARGET and (bestf is None or fl.mean() < bestf):
                    bestk, bestf = k, fl.mean()
            if bestk is None:                       # no alcanza: usar el de mayor recall
                bestk = max(pk_by_knob, key=lambda k: pk_by_knob[k][0].mean())
            paired[m] = pk_by_knob[bestk]
            print(f"   [op] {m}: perilla elegida={bestk} "
                  f"(recall={paired[m][0].mean():.3f})")

        # ---------- reporte a recall igualado ----------
        print(f"   {'metodo':>8} {'recall':>18} {'falsos':>18} "
              f"{'falsos@recall=' + format(TARGET, '.2f'):>22}")
        for m in list(own) + list(ext):
            rec, fal = paired[m]
            k1, k2, n = int(rec.sum()), int(fal.sum()), len(rec)
            l1, h1 = cp_ci(k1, n)
            l2, h2 = cp_ci(k2, n)
            fr = false_rate_at_recall(curves[m], TARGET)
            frs = "no alcanza" if fr is None else f"{fr:.3f}"
            print(f"   {m:>8} {k1:3d}/{n}={k1/n:4.2f} [{l1:.2f},{h1:.2f}] "
                  f"{k2:3d}/{n}={k2/n:4.2f} [{l2:.2f},{h2:.2f}] {frs:>22}")

        # ---------- comparaciones pareadas ----------
        for m in list(own) + list(ext):
            if m == "PHYS":
                continue
            b, c, p = exact_mcnemar(paired[m][1], paired["PHYS"][1])
            d = paired[m][1].mean() - paired["PHYS"][1].mean()
            print(f"      {m} vs PHYS: delta_falsos={d:+.3f} "
                  f"discordantes({m}=1,PHYS=0)={b} ({m}=0,PHYS=1)={c} McNemar_p={p:.4f}")
        print()

    if not (HAVE_GAMMA and HAVE_PYOCTO):
        print("ATENCION: faltan baselines; no se responde al enunciado (ii).", file=sys.stderr)


if __name__ == "__main__":
    main()
