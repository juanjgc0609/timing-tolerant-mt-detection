"""
task2c_run.py - Tarea 2C (preregistrada en PREREGISTRO_1j_1k.txt): comparadores cercanos (BOMP+SC, MTSTACK+SC), el
metodo convexo de la tesis en la rejilla (GL-FW+SC) y replicacion por red.  Metodos en task2c_methods.py.

Parte (i): las 1200 ventanas de 2A (semillas 2024 y 7331, regimenes thin y strong, 3 celdas, redes 0-4, 20 ventanas),
regeneradas con las funciones de task2a_confirm (importadas; mismos flujos aleatorios, ver task2b_ablation.py), con
control exacto por red: STA, retardos, catalogo, correcciones, verdad y origenes iguales a task2a_parts, y MT-FW+SC
recalculado (certified_fw_record, tcorr = correcciones) EXACTAMENTE igual al guardado en las ventanas 0 y 1.
Metodos BOMP+SC, MTSTACK+SC, GL-FW+SC.  Nucleo de OpenBLAS, fijado explicitamente en cada proceso con
OPENBLAS_CORETYPE (el nucleo por defecto de esta maquina cambio de Prescott en 2A/2B a SkylakeX): las 18 redes de 7331
que 2A produjo con SkylakeX (task2b_report.txt) se corren con SkylakeX, las demas redes de 2A y todas las de la parte
(ii) y de la semilla 999 con Prescott.  Cada red comprueba al empezar que el nucleo cargado (threadpoolctl) es el
requerido y lo guarda en su parte; run_all solo lanza las redes cuyo nucleo requerido es el del proceso.
Parte (ii): 25 redes nuevas por semilla, celda y regimen (indices 5-29), con la construccion de semillas de 2A
(task2a_confirm.network / catalogue_real / window_rngs con j = 5..29), 12 ventanas por red (las ventanas 0-11 de
esa red), metodos MT-FW+SC (certified_fw_record; en las ventanas 0 y 1 igual al denso dense_fw_record_c), QM-3C+SC
(configuracion congelada qm_recal_2a.json con las correcciones sumadas al LUT, task2a_confirm.run_qm_variants),
BOMP+SC y MTSTACK+SC.  Nucleo por defecto.

    python task2c_run.py validate                       # validaciones de los metodos (semilla 999) -> task2c_validation.txt
    python task2c_run.py smoke --nwin 3                 # tiempos en la semilla 999 y estimacion -> task2c_smoke.txt
    OPENBLAS_CORETYPE=SkylakeX python task2c_run.py run --part 1 --workers 2
    OPENBLAS_CORETYPE=Prescott python task2c_run.py run --part 1 --workers 2
    OPENBLAS_CORETYPE=Prescott python task2c_run.py run --part 2 --workers 2 [--nets-ii 25]
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

import task2b_ablation as B  # noqa: E402
import task2c_methods as M  # noqa: E402

T = B.T
X = T.X
QM = T.QM
SEEDS = (2024, 7331)
PARTS = os.path.join(HERE, "task2c_parts")
PARTS_2A = T.PARTS
FAIL_FILE = os.path.join(HERE, "task2c_CONTROL_FAILED.txt")
N_WIN_II = 12
NETS_II_FULL = 25
N_CONTROL = 2
SKX = {(7331, "strong", "ladder_6", j) for j in (2, 3, 4)} | {(7331, "thin", "ladder_6", j) for j in range(5)} | \
      {(7331, r, "ladder_4", j) for r in ("strong", "thin") for j in range(5)}
METHODS_I = ("BOMP+SC", "MTSTACK+SC", "GL-FW+SC")
METHODS_II = ("MT-FW+SC", "QM-3C+SC", "BOMP+SC", "MTSTACK+SC")


class ControlError(RuntimeError):
    pass


def fail(msg):
    with open(FAIL_FILE, "a") as fh:
        fh.write(time.strftime("%Y-%m-%d %H:%M:%S ") + msg + "\n")
    raise ControlError(msg)


def as_rec(r):
    return [[float(s), int(h)] for s, h in r]


def blas_libs():
    """[(arquitectura, biblioteca)] de todas las OpenBLAS cargadas (la de numpy y, si scipy.linalg esta cargado, la
    de scipy, que no se usa en MT-FW, SC ni en los metodos de 2C)."""
    import threadpoolctl
    return [(i.get("architecture"), os.path.basename(i.get("filepath", ""))) for i in threadpoolctl.threadpool_info()
            if i.get("user_api") == "blas"]


def numpy_blas_arch():
    return [a for a, f in blas_libs() if "scipy" not in f]


def kernel_info():
    return {"OPENBLAS_CORETYPE": os.environ.get("OPENBLAS_CORETYPE"), "blas": blas_libs()}


def required_kernel(part, seed, regime, cell, j):
    return "SkylakeX" if (part == 1 and (seed, regime, cell, j) in SKX) else "Prescott"


def check_kernel(part, seed, regime, cell, j, tag):
    req = required_kernel(part, seed, regime, cell, j)
    if numpy_blas_arch() != [req] or os.environ.get("OPENBLAS_CORETYPE") != req:
        fail(f"{tag}: nucleo OpenBLAS {kernel_info()} distinto del requerido {req}")


def write_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    os.replace(tmp, path)


def gl_summary(dg):
    pl = dg["per_lambda"]
    return {"lam_max": dg["lam_max"], "stopped": dg["stopped"], "n_entered": dg["n_entered"],
            "n_lam_solved": len(pl), "per_lambda": [[p["lam"], p["gap"], p["rel_gap"], p["gap_over_lam_omega"],
                                                     p["n_support"], p["n_fw"], p["kkt_off_max"], p["kkt_on_maxdev"],
                                                     p["stall"]]
                                                    for p in pl]}


def ops(STA, corr):
    grid = T.default_grid()
    garr = np.asarray(grid)
    Gs = [X.radiation_operator(g, STA) for g in garr]
    Qs = [X.basis(G) for G in Gs]
    BD = M.BlockDict(STA, grid, T.DEFAULT_TAUS, corr, Gs=Gs)
    return grid, Gs, Qs, BD, M.Stacker(BD)


# =========================================================================================== part (i)
def run_part1_net(args):
    seed, regime, cell, j, out_path, n_win_max = args
    t_start = time.time()
    tag = f"seed{seed}_{regime}_{cell}_net{j}"
    check_kernel(1, seed, regime, cell, j, tag)
    with open(os.path.join(PARTS_2A, tag + ".json")) as fh:
        P2A = json.load(fh)
    ci, nst, nclut, STA, statics, cat, drawn, corr = B.setup_net(seed, regime, cell, j)
    ctl = {"STA": STA.tolist() == P2A["STA"], "statics": statics.tolist() == P2A["statics"],
           "cat_drawn": drawn == P2A["cat_drawn"], "cat_src": [c[1].tolist() for c in cat] == P2A["cat_src"],
           "cat_tau": [c[2] for c in cat] == P2A["cat_tau"], "corr": corr.tolist() == P2A["corr"]}
    if not all(ctl.values()):
        fail(f"{tag}: red/catalogo/correcciones distintos de 2A: {ctl}")
    t0 = time.time()
    grid, Gs, Qs, BD, ST = ops(STA, corr)
    t_ops = time.time() - t0
    rngs = T.window_rngs(seed, cell, j, regime)
    n_win = min(P2A["n_win"], n_win_max or P2A["n_win"])
    W, n_ok = [], 0
    tim = {"BOMP+SC": [], "MTSTACK+SC": [], "GL-FW+SC": [], "control MT-FW+SC": []}
    for w in range(n_win):
        srcs, taus, D, _p, tr, _a = T.generate_window(STA, statics, nclut, T.CLUT_AMP[regime], rngs)
        st = P2A["windows"][w]
        truth = [list(map(float, s)) for s in srcs]
        taus_l = [float(t) for t in taus]
        if len(srcs) != st["k"] or truth != st["truth"] or taus_l != st["taus"]:
            fail(f"{tag} ventana {w}: fuentes u origenes distintos de 2A")
        MF, sd = T.matched(D)
        MFn = MF / sd
        if w < N_CONTROL:
            t0 = time.time()
            rs = as_rec(T.certified_fw_record(STA, MFn, grid, T.DEFAULT_TAUS, Gs, Qs, tr, tcorr=corr))
            tim["control MT-FW+SC"].append(time.time() - t0)
            if rs != st["records"]["MT-FW+SC"]:
                fail(f"{tag} ventana {w}: MT-FW+SC recalculado {rs} != guardado {st['records']['MT-FW+SC']}")
            n_ok += 1
        rec = {}
        t0 = time.time()
        rec["BOMP+SC"] = as_rec(M.bomp_record(BD, MFn, tr))
        tim["BOMP+SC"].append(time.time() - t0)
        t0 = time.time()
        rec["MTSTACK+SC"] = as_rec(ST.record(MFn, tr))
        tim["MTSTACK+SC"].append(time.time() - t0)
        t0 = time.time()
        dg = {}
        rec["GL-FW+SC"] = as_rec(M.glfw_record(BD, MFn, D, sd, tr, diag=dg))
        tim["GL-FW+SC"].append(time.time() - t0)
        W.append({"k": len(srcs), "truth": truth, "taus": taus_l, "records": rec, "gl": gl_summary(dg),
                  "gl_rec_nv": [[float(a), int(b), float(c), int(d), int(e)] for a, b, c, d, e in dg["rec_nv"]]})
    res = {"part": 1, "seed": seed, "regime": regime, "cell": cell, "net": j, "n_win": n_win, "nst": nst,
           "control": {**ctl, "windows_truth_taus": n_win, "mtfw_sc_exact": n_ok}, "corr": corr.tolist(),
           "kernel": kernel_info(), "windows": W, "timing": tim, "t_ops": t_ops, "seconds": time.time() - t_start}
    if out_path:
        write_json(out_path, res)
    return res


# ========================================================================================== part (ii)
def run_part2_net(args):
    seed, regime, cell, j, n_win, out_path, check2a = args
    t_start = time.time()
    tag = f"seed{seed}_{regime}_{cell}_net{j}"
    check_kernel(2, seed, regime, cell, j, tag)
    ci, nst, nclut, STA, statics, cat, drawn, corr = B.setup_net(seed, regime, cell, j)
    t0 = time.time()
    grid, Gs, Qs, BD, ST = ops(STA, corr)
    qmc = T.load_qm_cfgs()
    lut = QM.build_lut(STA, T.VP, T.VS, {**QM.QM_CFG, **qmc["QM-3C"]})
    t_ops = time.time() - t0
    tmpdir = os.path.join(T.SCRATCH, f"t2c_s{seed}_{regime}_{cell}_{j}")
    rngs = T.window_rngs(seed, cell, j, regime)
    P2A = None
    if check2a:
        with open(os.path.join(PARTS_2A, tag + ".json")) as fh:
            P2A = json.load(fh)
    W = []
    tim = {m: [] for m in METHODS_II}
    tim["identity check"] = []
    n_id = 0
    n_2a = 0
    for w in range(n_win):
        srcs, taus, D, _p, tr, _a = T.generate_window(STA, statics, nclut, T.CLUT_AMP[regime], rngs)
        MF, sd = T.matched(D)
        MFn = MF / sd
        rec = {}
        t0 = time.time()
        rec["MT-FW+SC"] = as_rec(T.certified_fw_record(STA, MFn, grid, T.DEFAULT_TAUS, Gs, Qs, tr, tcorr=corr))
        tim["MT-FW+SC"].append(time.time() - t0)
        t0 = time.time()
        evs = T.run_qm_variants(STA, D, tmpdir, {"QM-3C+SC": (qmc["QM-3C"], corr)}, lut=lut)
        rec["QM-3C+SC"] = as_rec(QM.qm_records(evs["QM-3C+SC"], tr))
        tim["QM-3C+SC"].append(time.time() - t0)
        t0 = time.time()
        rec["BOMP+SC"] = as_rec(M.bomp_record(BD, MFn, tr))
        tim["BOMP+SC"].append(time.time() - t0)
        t0 = time.time()
        rec["MTSTACK+SC"] = as_rec(ST.record(MFn, tr))
        tim["MTSTACK+SC"].append(time.time() - t0)
        if w < N_CONTROL:
            t0 = time.time()
            rd = as_rec(T.dense_fw_record_c(STA, MFn, grid, T.DEFAULT_TAUS, Gs, Qs, tr, tcorr=corr))
            if rd != rec["MT-FW+SC"]:
                fail(f"{tag} ventana {w}: MT-FW+SC certificado {rec['MT-FW+SC']} != denso {rd}")
            n_id += 1
            tim["identity check"].append(time.time() - t0)
            if P2A is not None:
                st = P2A["windows"][w]
                same = (rec["MT-FW+SC"] == st["records"]["MT-FW+SC"] and rec["QM-3C+SC"] == st["records"]["QM-3C+SC"])
                if not same:
                    fail(f"{tag} ventana {w}: la tuberia de la parte (ii) no reproduce 2A: "
                         f"{rec['MT-FW+SC']} / {st['records']['MT-FW+SC']} ; {rec['QM-3C+SC']} / {st['records']['QM-3C+SC']}")
                n_2a += 1
        W.append({"k": len(srcs), "truth": [list(map(float, s)) for s in srcs], "taus": [float(t) for t in taus],
                  "records": rec, "qm_events": [[float(x) for x in e] for e in evs["QM-3C+SC"]]})
    res = {"part": 2, "seed": seed, "regime": regime, "cell": cell, "net": j, "n_win": n_win, "nst": nst,
           "seeds": {"network": [seed, T.TAG_NET, ci, j], "statics": [seed, T.TAG_STAT, ci, j],
                     "catalogue": [seed, T.TAG_CAT, ci, j, T.REGIMES.index(regime)],
                     "windows": [seed, T.TAG_WIN, ci, j, T.REGIMES.index(regime)]},
           "STA": STA.tolist(), "statics": statics.tolist(), "corr": corr.tolist(), "cat_drawn": drawn,
           "resid_rms": float(np.sqrt(np.mean((corr - statics) ** 2))),
           "control": {"mtfw_cert_eq_dense": n_id, "reproduces_2a": n_2a}, "kernel": kernel_info(),
           "windows": W, "timing": tim, "t_ops": t_ops, "seconds": time.time() - t_start}
    if out_path:
        write_json(out_path, res)
    return res


# ============================================================================================ runner
def summarize(r, methods):
    k = [w["k"] for w in r["windows"]]
    ra = {m: T.summarize_records([w["records"][m] for w in r["windows"]], k) for m in methods}
    return " | ".join(f"{m} " + "/".join(f"{ra[m][str(c)]:.2f}" for c in T.CAPS) for m in methods)


def run_all(part, workers, nets_ii=NETS_II_FULL, only=None):
    os.makedirs(PARTS, exist_ok=True)
    if os.path.exists(FAIL_FILE):
        raise SystemExit(f"{FAIL_FILE} existe: el control fallo antes; revisar antes de reanudar")
    jobs = []
    for seed in SEEDS:
        for regime in T.REGIMES:
            for cell in T.CELLS:
                nets = range(T.N_NETS) if part == 1 else range(T.N_NETS, T.N_NETS + nets_ii)
                for j in nets:
                    if required_kernel(part, seed, regime, cell, j) != os.environ.get("OPENBLAS_CORETYPE"):
                        continue
                    if only is not None and f"{seed}_{regime}_{cell}_{j}" not in only:
                        continue
                    path = os.path.join(PARTS, f"p{part}_seed{seed}_{regime}_{cell}_net{j}.json")
                    if os.path.exists(path):
                        continue
                    if part == 1:
                        jobs.append((seed, regime, cell, j, path, None))
                    else:
                        jobs.append((seed, regime, cell, j, N_WIN_II, path, False))
    jobs.sort(key=lambda a: (-T.SCEN[a[2]][1], a[0], a[1], a[3]))
    print(f"parte {part}: {len(jobs)} redes pendientes; {kernel_info()}", flush=True)
    t0 = time.time()
    fn = run_part1_net if part == 1 else run_part2_net
    methods = METHODS_I if part == 1 else METHODS_II
    with ProcessPoolExecutor(workers) as ex:
        for r in ex.map(fn, jobs):
            print(f"   hecho p{part} {r['seed']} {r['regime']} {r['cell']} red {r['net']} en {r['seconds']:.0f} s "
                  f"(total {time.time() - t0:.0f} s); control {r['control']} :: {summarize(r, methods)}", flush=True)


# ========================================================================================= smoke test
def smoke(nwin, workers):
    """semilla 999: parte (i) en la red 0 (con el control de 2A: la semilla 999 tiene partes de 2A) y parte (ii) en
    la red 5 (y en la red 0 con el control de reproduccion de 2A), nwin ventanas por celda y regimen."""
    jobs1 = [(T.CAL_SEED, reg, cell, 0, None, nwin) for reg in T.REGIMES for cell in T.CELLS]
    jobs2 = [(T.CAL_SEED, reg, cell, 5, nwin, None, False) for reg in T.REGIMES for cell in T.CELLS]
    jobs2c = [(T.CAL_SEED, reg, cell, 0, 2, None, True) for reg in T.REGIMES for cell in T.CELLS]
    t0 = time.time()
    with ProcessPoolExecutor(workers) as ex:
        r1 = list(ex.map(run_part1_net_999, jobs1))
        r2 = list(ex.map(run_part2_net, jobs2))
        r2c = list(ex.map(run_part2_net, jobs2c))
    wall = time.time() - t0
    L = [f"prueba de humo 2C: semilla 999, {nwin} ventanas por celda y regimen, {workers} procesos, {wall:.0f} s de pared; "
         f"{kernel_info()}"]
    est1 = est2 = 0.0
    for r in r1:
        per = {m: float(np.mean(v)) for m, v in r["timing"].items() if v}
        pw = sum(per[m] for m in METHODS_I)
        ctl = float(np.sum(r["timing"]["control MT-FW+SC"]))
        est1 += len(SEEDS) * T.N_NETS * (T.N_WIN * pw + ctl + r["t_ops"] + 2.0)
        gl = [w["gl"] for w in r["windows"]]
        L.append(f"   (i)  {r['regime']:6s} {r['cell']:15s} por ventana {pw:6.2f} s [" + ", ".join(f"{m} {v:.2f}" for m, v in per.items())
                 + f"]; operadores {r['t_ops']:.1f} s; control {r['control']}; GL: lambdas resueltos "
                 + ",".join(str(g["n_lam_solved"]) for g in gl) + "; soporte final "
                 + ",".join(str(g["per_lambda"][-1][4]) for g in gl) + f"; max brecha rel {max(p[2] for g in gl for p in g['per_lambda']):.1e}"
                 + f"; max KKT fuera {max(p[6] for g in gl for p in g['per_lambda']):.7f}")
    for r in r2:
        per = {m: float(np.mean(v)) for m, v in r["timing"].items() if v}
        pw = sum(per[m] for m in METHODS_II)
        idc = float(np.sum(r["timing"]["identity check"]))
        est2 += len(SEEDS) * NETS_II_FULL * (N_WIN_II * pw + idc + r["t_ops"] + 2.0)
        L.append(f"   (ii) {r['regime']:6s} {r['cell']:15s} por ventana {pw:6.2f} s [" + ", ".join(f"{m} {v:.2f}" for m, v in per.items())
                 + f"]; operadores+LUT {r['t_ops']:.1f} s; control {r['control']}")
    for r in r2c:
        L.append(f"   (ii) tuberia sobre la red 0 de 2A (semilla 999) {r['regime']} {r['cell']}: control {r['control']}")
    L.append(f"   estimacion parte (i): {est1 / 3600:.2f} h de CPU; parte (ii) con 25 redes: {est2 / 3600:.2f} h de CPU; "
             f"total {(est1 + est2) / 3600:.2f} h de CPU, {(est1 + est2) / 7200:.2f} h en 2 nucleos (umbral de la regla: 8 h)")
    est2b = est2 * 15 / 25
    L.append(f"   con 15 redes nuevas: total {(est1 + est2b) / 7200:.2f} h en 2 nucleos")
    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(HERE, "task2c_smoke.txt"), "w") as fh:
        fh.write(txt + "\n")


def run_part1_net_999(args):
    """la parte (i) en la semilla 999 (redes 0-4 de 2A): el nucleo requerido es el por defecto."""
    return run_part1_net(args)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["validate", "smoke", "run"])
    ap.add_argument("--part", type=int, choices=[1, 2])
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--nwin", type=int, default=3)
    ap.add_argument("--nets-ii", type=int, default=NETS_II_FULL)
    ap.add_argument("--nets", default=None)
    a = ap.parse_args()
    os.makedirs(T.SCRATCH, exist_ok=True)
    if a.stage == "validate":
        import task2c_validate as V
        raise SystemExit(0 if V.main() else 1)
    if a.stage == "smoke":
        smoke(a.nwin, a.workers)
        return
    run_all(a.part, a.workers, a.nets_ii, None if a.nets is None else set(a.nets.split(",")))


if __name__ == "__main__":
    main()
