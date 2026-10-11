"""
task2b_ablation.py - Tarea 2B (preregistrada en PREREGISTRO_1j_1k.txt): ablacion del patron de radiacion dentro del
mismo solver de 2A.

MT-FW-FREE+SC es la llamada de MT-FW+SC de 2A (certified_fw_record, kmax 5, min_sep 5 km, CAPV 6, tcorr = las
correcciones de estacion de 2A) con el operador de radiacion reemplazado por la identidad: Gs = [I] y Qs = [I] en
cada nodo, I de tamano 2*ns*3 (vf = v.reshape(len(taus), -1) en orden (fase, estacion, componente) y
pred = (G @ m).reshape(2, ns, 3) en task1x_mtglrt.fw_record). Su estadistico es ||vt||^2 / 2, la mitad de la energia
total de los vectores de pico (con el tope CAPV por estacion y fase), invariante al signo; m = lstsq(I, vf) = vf, asi
que la deflacion resta los vectores observados.

Las ventanas, redes, catalogos y correcciones se regeneran con las funciones de task2a_confirm (importadas, no
copiadas), con los mismos flujos aleatorios que task2a_confirm.run_job: las ventanas solo consumen rng, rng_m y rng_h
(window_rngs[:3]); el picker, degrade_y (rng_d) y el mapa de nodos (rng_n) no tocan esos flujos.

Control por red: STA, retardos, catalogo (localizacion y origen) y correcciones iguales a los guardados en
task2a_parts; en cada ventana la verdad y los origenes iguales; en las ventanas 0 y 1, MT-FW+SC recalculado igual
EXACTAMENTE a los registros guardados, y MT-FW-FREE+SC certificado igual al denso (dense_fw_record_c). Cualquier
diferencia detiene la corrida (excepcion y archivo task2b_CONTROL_FAILED.txt).

    python task2b_ablation.py selfcheck                     # comprobaciones de la identidad en una ventana
    python task2b_ablation.py run --workers 2               # semillas 2024 y 7331, ambos regimenes (reanudable)
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import task2a_confirm as T  # noqa: E402
from certified_fw import dense_fw_record_c, node_vectors_c  # noqa: E402

X = T.X
SEEDS = (2024, 7331)
PARTS_2A = T.PARTS
PARTS = os.path.join(HERE, "task2b_parts")
FAIL_FILE = os.path.join(HERE, "task2b_CONTROL_FAILED.txt")
N_CONTROL = 2
METHODS = ("MT-FW+SC", "MT-FW-FREE+SC")


class ControlError(RuntimeError):
    pass


def as_rec(r):
    """formato de los registros guardados en task2a_parts."""
    return [[float(s), int(h)] for s, h in r]


def identity_ops(n_nodes, nst):
    I = np.eye(2 * nst * 3)
    return [I] * n_nodes, [I] * n_nodes


def setup_net(seed, regime, cell, j):
    """red, retardos, catalogo y correcciones de 2A (mismas llamadas y semillas que task2a_confirm.run_job)."""
    ci, _kind, nst, nclut = T.net_setup(seed, cell)
    ri = T.REGIMES.index(regime)
    STA, statics = T.network(seed, cell, j)
    cat, drawn, _truth = T.catalogue_real(STA, T.K_CAT, statics, [seed, T.TAG_CAT, ci, j, ri], nclut, T.CLUT_AMP[regime])
    corr, _nacc = T.estimate_station_corrections(STA, cat)
    return ci, nst, nclut, STA, statics, cat, drawn, corr


def fail(msg):
    with open(FAIL_FILE, "a") as fh:
        fh.write(time.strftime("%Y-%m-%d %H:%M:%S ") + msg + "\n")
    raise ControlError(msg)


def blas_arch():
    try:
        import threadpoolctl
        return [i.get("architecture") for i in threadpoolctl.threadpool_info() if i.get("user_api") == "blas"]
    except Exception:  # noqa: BLE001
        return None


def run_net(args):
    seed, regime, cell, j, out_path = args
    t_start = time.time()
    tag = f"seed{seed}_{regime}_{cell}_net{j}"
    with open(os.path.join(PARTS_2A, tag + ".json")) as fh:
        P2A = json.load(fh)
    ci, nst, nclut, STA, statics, cat, drawn, corr = setup_net(seed, regime, cell, j)
    ctl = {"STA": STA.tolist() == P2A["STA"], "statics": statics.tolist() == P2A["statics"],
           "cat_drawn": drawn == P2A["cat_drawn"], "cat_src": [c[1].tolist() for c in cat] == P2A["cat_src"],
           "cat_tau": [c[2] for c in cat] == P2A["cat_tau"], "corr": corr.tolist() == P2A["corr"]}
    if not all(ctl.values()):
        fail(f"{tag}: red/catalogo/correcciones distintos de 2A: {ctl}")
    grid = T.default_grid()
    garr = np.asarray(grid)
    taus_grid = T.DEFAULT_TAUS
    Gs = [X.radiation_operator(g, STA) for g in garr]
    Qs = [X.basis(G) for G in Gs]
    GsF, QsF = identity_ops(len(garr), nst)
    assert Gs[0].shape == (2 * nst * 3, 6) and GsF[0].shape == (Gs[0].shape[0],) * 2
    rngs = T.window_rngs(seed, cell, j, regime)
    n_win = P2A["n_win"]
    W, n_sc_ok, n_free_dense_ok, t_free = [], 0, 0, []
    for w in range(n_win):
        srcs, taus, D, _picks, tr, _amp = T.generate_window(STA, statics, nclut, T.CLUT_AMP[regime], rngs)
        st = P2A["windows"][w]
        truth = [list(map(float, s)) for s in srcs]
        taus_l = [float(t) for t in taus]
        if len(srcs) != st["k"] or truth != st["truth"] or taus_l != st["taus"]:
            fail(f"{tag} ventana {w}: fuentes u origenes distintos de 2A")
        MF, sd = T.matched(D)
        MFn = MF / sd
        t0 = time.time()
        rf = as_rec(T.certified_fw_record(STA, MFn, grid, taus_grid, GsF, QsF, tr, tcorr=corr))
        t_free.append(time.time() - t0)
        if w < N_CONTROL:
            rs = as_rec(T.certified_fw_record(STA, MFn, grid, taus_grid, Gs, Qs, tr, tcorr=corr))
            if rs != st["records"]["MT-FW+SC"]:
                fail(f"{tag} ventana {w}: MT-FW+SC recalculado {rs} != guardado {st['records']['MT-FW+SC']}")
            n_sc_ok += 1
            rd = as_rec(dense_fw_record_c(STA, MFn, grid, taus_grid, GsF, QsF, tr, tcorr=corr))
            if rd != rf:
                fail(f"{tag} ventana {w}: MT-FW-FREE+SC certificado {rf} != denso {rd}")
            n_free_dense_ok += 1
        W.append({"k": len(srcs), "truth": truth, "taus": taus_l,
                  "records": {"MT-FW+SC": st["records"]["MT-FW+SC"], "MT-FW-FREE+SC": rf}})
    res = {"seed": seed, "regime": regime, "cell": cell, "net": j, "n_win": n_win, "nst": nst,
           "control": {**ctl, "windows_truth_taus": n_win, "mtfw_sc_exact": n_sc_ok, "free_cert_eq_dense": n_free_dense_ok},
           "corr": corr.tolist(), "blas_arch": blas_arch(),
           "openblas_coretype": os.environ.get("OPENBLAS_CORETYPE"), "windows": W, "free_seconds_mean": float(np.mean(t_free)),
           "seconds": time.time() - t_start}
    if out_path:
        tmp = out_path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(res, fh)
        os.replace(tmp, out_path)
    return res


# ------------------------------------------------------------------------------------ aggregation
def aggregate(seed, regime):
    out = []
    for cell in T.CELLS:
        ps = []
        for j in range(T.N_NETS):
            with open(os.path.join(PARTS, f"seed{seed}_{regime}_{cell}_net{j}.json")) as fh:
                ps.append(json.load(fh))
        W = [w for p in ps for w in p["windows"]]
        k = [w["k"] for w in W]
        methods = {m: {"records": [w["records"][m] for w in W],
                       "recall_at": T.summarize_records([w["records"][m] for w in W], k)} for m in METHODS}
        out.append({"scenario": cell, "seed": seed, "regime": regime, "clut_amp": T.CLUT_AMP[regime],
                    "n_windows": len(W), "k_eval": k, "net_of_window": [p["net"] for p in ps for _ in p["windows"]],
                    "truth": [w["truth"] for w in W], "taus": [w["taus"] for w in W],
                    "control": [p["control"] for p in ps], "methods": methods,
                    "seconds": float(sum(p["seconds"] for p in ps))})
    return out


def run_all(workers, only=None, do_aggregate=True):
    os.makedirs(PARTS, exist_ok=True)
    if os.path.exists(FAIL_FILE):
        raise SystemExit(f"{FAIL_FILE} existe: el control fallo antes; revisar antes de reanudar")
    jobs = []
    for seed in SEEDS:
        for regime in T.REGIMES:
            for cell in T.CELLS:
                for j in range(T.N_NETS):
                    path = os.path.join(PARTS, f"seed{seed}_{regime}_{cell}_net{j}.json")
                    if only is not None and f"{seed}_{regime}_{cell}_{j}" not in only:
                        continue
                    if not os.path.exists(path):
                        jobs.append((seed, regime, cell, j, path))
    jobs.sort(key=lambda a: (-T.SCEN[a[2]][1], a[0], a[1], a[3]))
    print(f"{len(jobs)} redes pendientes", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(workers) as ex:
        for r in ex.map(run_net, jobs):
            k = [w["k"] for w in r["windows"]]
            ra = {m: T.summarize_records([w["records"][m] for w in r["windows"]], k) for m in METHODS}
            print(f"   hecho {r['seed']} {r['regime']} {r['cell']} red {r['net']} en {r['seconds']:.0f} s "
                  f"(total {time.time() - t0:.0f} s); control {r['control']} :: "
                  + " | ".join(f"{m} " + "/".join(f"{ra[m][str(c)]:.2f}" for c in T.CAPS) for m in ra), flush=True)
    if not do_aggregate:
        return
    for seed in SEEDS:
        for regime in T.REGIMES:
            agg = aggregate(seed, regime)
            fn = os.path.join(HERE, f"task2b_seed{seed}_{regime}.json")
            with open(fn, "w") as fh:
                json.dump(agg, fh)
            for c in agg:
                print(f"   {seed} {regime} {c['scenario']}: " + " | ".join(
                    f"{m} " + "/".join(f"{c['methods'][m]['recall_at'][str(x)]:.2f}" for x in T.CAPS) for m in METHODS),
                    flush=True)
            print(f"escrito {fn}", flush=True)


# ------------------------------------------------------------------------------------------ selfcheck
def selfcheck(seed=2024, regime="thin", cell="ladder_6", j=0, out="task2b_selfcheck.txt"):
    """en una ventana con al menos un evento: (1) formas de vf, G y pred; (2) el puntaje FREE de cada hipotesis es la
    mitad de la energia total de los vectores de pico (con tope), y el primer puntaje del registro es su maximo sobre
    todas las hipotesis; m = lstsq(I, vf) = vf; (3) cambiar el signo de la forma de onda de una estacion deja los
    registros FREE identicos y cambia los de MT-FW+SC."""
    L = []
    P = L.append
    ci, nst, nclut, STA, statics, cat, drawn, corr = setup_net(seed, regime, cell, j)
    grid = T.default_grid()
    garr = np.asarray(grid)
    taus_grid = np.asarray(T.DEFAULT_TAUS)
    Gs = [X.radiation_operator(g, STA) for g in garr]
    Qs = [X.basis(G) for G in Gs]
    GsF, QsF = identity_ops(len(garr), nst)
    rngs = T.window_rngs(seed, cell, j, regime)
    for w in range(T.N_WIN):
        srcs, taus, D, _p, tr, _a = T.generate_window(STA, statics, nclut, T.CLUT_AMP[regime], rngs)
        if len(srcs) >= 1:
            break
    P(f"ventana: semilla {seed}, {regime}, {cell}, red {j}, ventana {w}, {len(srcs)} fuente(s), ns = {nst}")
    MF, sd = T.matched(D)
    MFn = MF / sd
    ok = True
    # (1) formas
    v, tpk, _s, _e = node_vectors_c(MFn, STA, garr[0], taus_grid, corr)
    vf = v.reshape(len(taus_grid), -1)
    shp = (v.shape == (len(taus_grid), 2, nst, 3) and vf.shape[1] == 2 * nst * 3 == Gs[0].shape[0] == GsF[0].shape[0]
           and np.array_equal((GsF[0] @ vf[0]).reshape(2, nst, 3), v[0]))
    P(f"(1) v {v.shape}, vf {vf.shape}, G {Gs[0].shape}, I {GsF[0].shape}; (I @ vf).reshape(2, ns, 3) == v: {shp}")
    ok &= shp
    # (2) estadistico = energia / 2 en todas las hipotesis; registro[0] = maximo
    emax, dmax, sign_inv = -np.inf, 0.0, True
    for gi, src in enumerate(garr):
        v, _t, _s, _e = node_vectors_c(MFn, STA, src, taus_grid, corr)
        vf = v.reshape(len(taus_grid), -1)
        sc = np.sum((vf @ QsF[gi]) ** 2, axis=1) / 2.0
        en = 0.5 * np.einsum("tbic,tbic->t", v, v)
        dmax = max(dmax, float(np.max(np.abs(sc - en) / np.maximum(en, 1e-300))))
        flip = v.copy()
        flip[:, :, 0, :] *= -1.0
        sign_inv &= bool(np.array_equal(np.sum((flip.reshape(len(taus_grid), -1) @ QsF[gi]) ** 2, axis=1) / 2.0, sc))
        emax = max(emax, float(en.max()))
    rF = T.certified_fw_record(STA, MFn, grid, taus_grid, GsF, QsF, tr, tcorr=corr)
    rel0 = abs(rF[0][0] - emax) / emax
    P(f"(2) max_hip |FREE - ||vt||^2/2| / (||vt||^2/2) = {dmax:.2e} (675 nodos x 31 origenes); primer puntaje del "
      f"registro FREE {rF[0][0]:.6f}, maximo de ||vt||^2/2 {emax:.6f} (dif. rel. {rel0:.1e}); puntaje invariante al "
      f"cambiar el signo del vector de la estacion 0: {sign_inv}")
    ok &= dmax < 1e-12 and rel0 < 1e-12 and sign_inv
    # m = lstsq(I, vf) en el primer paso (el vector del ganador)
    best = (-np.inf, None, None)
    for gi, src in enumerate(garr):
        v, _t, _s, _e = node_vectors_c(MFn, STA, src, taus_grid, corr)
        vf = v.reshape(len(taus_grid), -1)
        sc = np.sum((vf @ QsF[gi]) ** 2, axis=1) / 2.0
        jj = int(sc.argmax())
        if sc[jj] > best[0]:
            best = (float(sc[jj]), gi, vf[jj])
    m, *_ = np.linalg.lstsq(GsF[best[1]], best[2], rcond=None)
    dm = float(np.max(np.abs(GsF[best[1]] @ m - best[2])))
    P(f"    deflacion: max |I lstsq(I, vf) - vf| = {dm:.1e} (se restan los vectores observados)")
    ok &= dm < 1e-12
    # (3) cambio de signo de una estacion
    rS = T.certified_fw_record(STA, MFn, grid, taus_grid, Gs, Qs, tr, tcorr=corr)
    P(f"    registros sin cambio de signo: MT-FW+SC {as_rec(rS)}")
    P(f"                                   FREE     {as_rec(rF)}")
    n_changed = 0
    for i in range(nst):
        D2 = D.copy()
        D2[i] *= -1.0
        MF2, sd2 = T.matched(D2)
        M2 = MF2 / sd2
        rF2 = T.certified_fw_record(STA, M2, grid, taus_grid, GsF, QsF, tr, tcorr=corr)
        rS2 = T.certified_fw_record(STA, M2, grid, taus_grid, Gs, Qs, tr, tcorr=corr)
        same_free = as_rec(rF2) == as_rec(rF)
        changed_mt = as_rec(rS2) != as_rec(rS)
        n_changed += changed_mt
        ok &= same_free
        P(f"(3) estacion {i} con signo cambiado (sd {sd2 == sd}): FREE identico {same_free}; MT-FW+SC cambia {changed_mt} "
          f"-> {as_rec(rS2)[:2]}")
    ok &= n_changed >= 1
    P(f"selfcheck {'correcto' if ok else 'CON FALLAS'} (MT-FW+SC cambia en {n_changed}/{nst} estaciones)")
    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(HERE, out), "w") as fh:
        fh.write(txt + "\n")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["selfcheck", "run", "aggregate"])
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--nets", default=None, help="solo estas redes: seed_regime_cell_net separadas por comas")
    ap.add_argument("--no-aggregate", action="store_true")
    a = ap.parse_args()
    if a.stage == "selfcheck":
        raise SystemExit(0 if selfcheck() else 1)
    if a.stage == "aggregate":
        for seed in SEEDS:
            for regime in T.REGIMES:
                with open(os.path.join(HERE, f"task2b_seed{seed}_{regime}.json"), "w") as fh:
                    json.dump(aggregate(seed, regime), fh)
        return
    run_all(a.workers, None if a.nets is None else set(a.nets.split(",")), not a.no_aggregate)


if __name__ == "__main__":
    main()
