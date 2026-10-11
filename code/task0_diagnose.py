"""
task0_diagnose.py — dos diagnósticos, ejecutar ANTES de volver a correr task0_run.py

A) PyOcto devuelve cero eventos en todas las celdas. Se prueban variantes del formato de entrada
   hasta encontrar cuál devuelve eventos con datos limpios y 12 estaciones.
B) GaMMA corre en su punto de operación por defecto mientras TOL y PHYS están calibrados, así que
   sus tasas de falsos no son comparables. Se barre `min_picks_per_eq` para localizar el punto que
   iguala el recall, que es el que hay que usar en la comparación.

    python task0_diagnose.py
"""
import sys
import numpy as np
import pandas as pd
from task0_common import (make_network, generate, picks_to_dataframe, stations_to_dataframe,
                          evaluate_catalog, VP, VS, EPOCH)

TRUTH = [np.array([6.0, -8.0, 9.0])]
TAUS = [6.0]


def caso_facil(n_sta=12, n_clutter=5, seed=7):
    """Caso deliberadamente fácil: muchas estaciones, poco clutter. Si un asociador falla aquí,
    está roto; no es una limitación del método."""
    STA = make_network("random", n_sta, seed=seed)
    rng = np.random.default_rng(seed)
    D, picks, truth = generate(STA, TRUTH, TAUS, n_clutter=n_clutter, rng=rng)
    return STA, D, picks, truth


# ===================================================================== A) PyOcto
def diagnose_pyocto():
    try:
        import pyocto
    except Exception as e:                                        # noqa: BLE001
        print(f"[A] PyOcto no importable: {e}")
        return
    print("\n" + "=" * 70)
    print("[A] DIAGNOSTICO PYOCTO — caso facil: 12 estaciones, 5 arribos de clutter")
    print("=" * 70)
    STA, D, picks, truth = caso_facil()
    pk = picks_to_dataframe(picks, "pyocto")
    stdf = stations_to_dataframe(STA, "pyocto")

    print(f"\npicks: {len(pk)} filas")
    print(pk.head(4).to_string(index=False))
    print(f"  dtype de 'time': {pk['time'].dtype}")
    print(f"\nestaciones: {len(stdf)} filas")
    print(stdf.head(3).to_string(index=False))

    vm = pyocto.VelocityModel0D(p_velocity=VP, s_velocity=VS, tolerance=1.0)
    xlim = (float(STA[:, 0].min() - 30), float(STA[:, 0].max() + 30))
    ylim = (float(STA[:, 1].min() - 30), float(STA[:, 1].max() + 30))

    variantes = []

    # V1: tal como está en task0_run.py
    variantes.append(("V1 datetime + n_picks=6", pk.copy(), dict(
        n_picks=6, n_p_picks=2, n_s_picks=2, n_p_and_s_picks=2)))

    # V2: tiempo como float de segundos UNIX (PyOcto acepta ambos segun version)
    pk2 = pk.copy()
    pk2["time"] = pd.to_datetime(pk2["time"]).astype("int64") / 1e9
    variantes.append(("V2 time float UNIX", pk2, dict(
        n_picks=6, n_p_picks=2, n_s_picks=2, n_p_and_s_picks=2)))

    # V3: pandas datetime64 explicito
    pk3 = pk.copy()
    pk3["time"] = pd.to_datetime(pk3["time"])
    variantes.append(("V3 datetime64 explicito", pk3, dict(
        n_picks=6, n_p_picks=2, n_s_picks=2, n_p_and_s_picks=2)))

    # V4: umbrales de picks relajados al minimo
    variantes.append(("V4 n_picks=4 relajado", pk3.copy(), dict(
        n_picks=4, n_p_picks=1, n_s_picks=1, n_p_and_s_picks=1)))

    # V5: fase en minusculas
    pk5 = pk3.copy(); pk5["phase"] = pk5["phase"].str.lower()
    variantes.append(("V5 fase en minusculas", pk5, dict(
        n_picks=4, n_p_picks=1, n_s_picks=1, n_p_and_s_picks=1)))

    print("\nvariantes probadas:")
    for nombre, pkv, kw in variantes:
        try:
            assoc = pyocto.OctoAssociator(
                xlim=xlim, ylim=ylim, zlim=(0.0, 30.0),
                time_before=300.0, velocity_model=vm, **kw)
            ev, asg = assoc.associate(pkv, stdf)
            n = 0 if ev is None else len(ev)
            print(f"   {nombre:28s} -> {n} eventos", end="")
            if n:
                print(f"   columnas: {list(ev.columns)}")
                print(ev.head(2).to_string(index=False))
            else:
                print()
        except Exception as e:                                    # noqa: BLE001
            print(f"   {nombre:28s} -> EXCEPCION: {type(e).__name__}: {e}")

    print("\nSi ninguna variante devuelve eventos, prueba el constructor geografico, que es el"
          "\ndocumentado y el que la mayoria de ejemplos usa:")
    print("   assoc = pyocto.OctoAssociator.from_area(lat=(..), lon=(..), zlim=(0,30),"
          " time_before=300, velocity_model=vm)")
    print("   assoc.transform_stations(stations)   # anade x,y,z desde lat/lon")
    print("Eso exige convertir tus coordenadas km a lat/lon con un origen arbitrario.")


# ===================================================================== B) GaMMA
def diagnose_gamma():
    try:
        from gamma.utils import association
    except Exception as e:                                        # noqa: BLE001
        print(f"[B] GaMMA no importable: {e}")
        return
    print("\n" + "=" * 70)
    print("[B] CALIBRACION DE PUNTO DE OPERACION DE GAMMA")
    print("=" * 70)
    print("Se barre min_picks_per_eq. Hay que elegir el valor que iguale el RECALL de los")
    print("metodos propios (~1.00) con el menor numero de eventos emitidos. Sin esto la")
    print("comparacion de tasas de falsos no significa nada.\n")

    STA = make_network("onesided", 8, seed=2024)
    rng = np.random.default_rng(2024)
    xr = [float(STA[:, 0].min() - 30), float(STA[:, 0].max() + 30)]
    yr = [float(STA[:, 1].min() - 30), float(STA[:, 1].max() + 30)]

    print(f"{'min_picks':>10} {'eventos/real':>13} {'recall':>8} {'falsos':>8}")
    print("-" * 44)
    for mpe in (5, 6, 8, 10, 12, 14):
        nev = rec = fal = 0
        R = 10
        for _ in range(R):
            D, picks, truth = generate(STA, TRUTH, TAUS, n_clutter=40, rng=rng)
            pkg = picks_to_dataframe(picks, "gamma")
            stg = stations_to_dataframe(STA, "gamma")
            cfg = {
                "dims": ["x(km)", "y(km)", "z(km)"],
                "x(km)": xr, "y(km)": yr, "z(km)": [0.0, 30.0],
                "vel": {"p": VP, "s": VS}, "method": "BGMM",
                "oversample_factor": 5, "use_amplitude": False,
                "use_dbscan": True, "dbscan_eps": 10.0, "dbscan_min_samples": 3,
                "min_picks_per_eq": mpe,
                "max_sigma11": 2.0, "max_sigma22": 1.0, "max_sigma12": 1.0,
                "bfgs_bounds": ((xr[0] - 1, xr[1] + 1), (yr[0] - 1, yr[1] + 1),
                                (0, 30), (None, None)),
                "covariance_prior": [5.0, 5.0], "initial_points": [1, 1, 1],
            }
            try:
                cat, _asg = association(pkg, stg, cfg, method="BGMM")
            except Exception as e:                                # noqa: BLE001
                print(f"   min_picks={mpe}: EXCEPCION {e}")
                break
            if cat is None or len(cat) == 0:
                continue
            cdf = pd.DataFrame(cat)
            cols = {c.lower().replace("(km)", "").strip(): c for c in cdf.columns}
            try:
                ev = [np.array([float(r[cols["x"]]), float(r[cols["y"]]), float(r[cols["z"]])])
                      for _, r in cdf.iterrows()]
            except KeyError:
                print(f"   columnas del catalogo: {list(cdf.columns)}")
                return
            nev += len(ev)
            r, nf = evaluate_catalog(ev, truth)
            rec += int(r); fal += int(nf > 0)
        print(f"{mpe:10d} {nev/R:13.1f} {rec/R:8.2f} {fal/R:8.2f}")

    print("\nElige el min_picks_per_eq mas alto que conserve recall >= 0.97 y usalo en")
    print("task0_run.py. Declara el valor elegido y el criterio: es sintonizacion del")
    print("baseline a FAVOR, que es legitima; sintonizarlo en contra no lo seria.")


if __name__ == "__main__":
    diagnose_pyocto()
    diagnose_gamma()
