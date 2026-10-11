"""
task2d_run.py - Tarea 2D (preregistrada en PREREGISTRO_1j_1k.txt, seccion "Tarea 2D").  Metodos en task2d_methods.py.

Parte (i): las 1200 ventanas de 2A (semillas 2024 y 7331, regimenes thin y strong, 3 celdas, redes 0-4, 20 ventanas),
regeneradas con task2b_ablation.setup_net y task2a_confirm.window_rngs / generate_window (como task2c_run.run_part1_net),
con el nucleo de OpenBLAS de 2C (task2c_run.required_kernel: SkylakeX en las 18 redes de 7331 listadas en
task2b_report.txt, Prescott en las demas; comprobado con threadpoolctl y guardado en cada parte).  Control exacto por
red: STA, retardos, catalogo y correcciones iguales a task2a_parts; verdad y origenes iguales en cada ventana; MT-FW+SC
recalculado (certified_fw_record, tcorr = correcciones) EXACTAMENTE igual al guardado en las ventanas 0 y 1.  Metodos:
BOMP-PK+SC, MTSTACK-PK+SC, GL-FW-PK+SC, MT-FW-EX+SC, MT-FW-NOCAP+SC.
Parte (ii): 200 ventanas vacias por celda y regimen de la semilla 999: las 5 redes de calibracion de 2A (mismas
estaciones, retardos y correcciones: task2b_ablation.setup_net(999, ...)), 40 ventanas vacias por red generadas con
generate_real sin fuentes (clutter y ruido del regimen) con flujos nuevos [999, 2701, ci, j, ri, 41 | 233 | 61];
MT-FW+SC (certified_fw_record, registro completo; el maximo por ventana es su primer puntaje) y QM-3C+SC (configuracion
congelada qm_recal_2a.json con las correcciones en el LUT, task2a_confirm.run_qm_variants).  Nucleo Prescott.
Control: STA, retardos y correcciones iguales a las partes de 2A de la semilla 999, y la tuberia reproduce MT-FW+SC y
QM-3C+SC guardados en las ventanas 0 y 1 de 2A de cada red (se registra; ver D5 de task2a_report.txt).
Los umbrales split-conformal se calculan y aplican en analyze_2d.py.

    python task2d_run.py validate                          # task2d_validate.py -> task2d_validation.txt
    OPENBLAS_CORETYPE=Prescott python task2d_run.py smoke --nwin 2
    OPENBLAS_CORETYPE=SkylakeX python task2d_run.py run --part 1 --workers 2
    OPENBLAS_CORETYPE=Prescott python task2d_run.py run --part 1 --workers 2
    OPENBLAS_CORETYPE=Prescott python task2d_run.py run --part 2 --workers 2 [--n-empty 40]
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
import task2c_run as C  # noqa: E402
import task2d_methods as D2  # noqa: E402

T = B.T
X = T.X
QM = T.QM
SEEDS = (2024, 7331)
PARTS = os.path.join(HERE, "task2d_parts")
PARTS_2A = T.PARTS
FAIL_FILE = os.path.join(HERE, "task2d_CONTROL_FAILED.txt")
N_CONTROL = 2
N_EMPTY_PER_NET = 40               # 200 ventanas vacias por celda y regimen (regla de presupuesto: 20 -> 100)
TAG_EMPTY = 2701
METHODS_I = ("BOMP-PK+SC", "MTSTACK-PK+SC", "GL-FW-PK+SC", "MT-FW-EX+SC", "MT-FW-NOCAP+SC")


class ControlError(RuntimeError):
    pass


def fail(msg):
    with open(FAIL_FILE, "a") as fh:
        fh.write(time.strftime("%Y-%m-%d %H:%M:%S ") + msg + "\n")
    raise ControlError(msg)


def as_rec(r):
    return [[float(s), int(h)] for s, h in r]


def check_kernel(part, seed, regime, cell, j, tag):
    req = C.required_kernel(1 if part == 1 else 2, seed, regime, cell, j)
    if C.numpy_blas_arch() != [req] or os.environ.get("OPENBLAS_CORETYPE") != req:
        fail(f"{tag}: nucleo OpenBLAS {C.kernel_info()} distinto del requerido {req}")
    return req


def write_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    os.replace(tmp, path)


def gl_summary(dg):
    pl = dg["per_lambda"]
    return {"lam_max": dg["lam_max"], "stopped": dg["stopped"], "n_entered": dg["n_entered"], "n_lam_solved": len(pl),
            "per_lambda": [[p["lam"], p["gap"], p["rel_gap"], p["gap_over_lam_omega"], p["n_support"], p["n_fw"],
                            p["kkt_off_max"], p["kkt_on_maxdev"], p["stall"], p["fw_capped"]] for p in pl]}


def ops(STA, corr):
    grid, Gs, Qs, BD, ST = C.ops(STA, corr)
    return grid, Gs, Qs, BD, ST, D2.PKDict(BD)


def methods_window(STA, grid, Gs, Qs, BD, ST, PD, D, MFn, sd, tr, corr, tim):
    rec = {}
    t0 = time.time()
    dgb = {}
    rec["BOMP-PK+SC"] = as_rec(D2.bomp_pk_record(PD, MFn, tr, diag=dgb))
    tim["BOMP-PK+SC"].append(time.time() - t0)
    t0 = time.time()
    rec["MTSTACK-PK+SC"] = as_rec(D2.mtstack_pk_record(PD, ST, MFn, tr))
    tim["MTSTACK-PK+SC"].append(time.time() - t0)
    t0 = time.time()
    dg = {}
    rec["GL-FW-PK+SC"] = as_rec(D2.glfw_pk_record(PD, MFn, D, sd, tr, diag=dg))
    tim["GL-FW-PK+SC"].append(time.time() - t0)
    t0 = time.time()
    rec["MT-FW-EX+SC"] = as_rec(D2.fw_record_v(STA, MFn, grid, T.DEFAULT_TAUS, Gs, Qs, tr, corr, peak=False, cap=True))
    tim["MT-FW-EX+SC"].append(time.time() - t0)
    t0 = time.time()
    rec["MT-FW-NOCAP+SC"] = as_rec(D2.fw_record_v(STA, MFn, grid, T.DEFAULT_TAUS, Gs, Qs, tr, corr, peak=True, cap=False))
    tim["MT-FW-NOCAP+SC"].append(time.time() - t0)
    extra = {"gl": gl_summary(dg),
             "gl_rec_nv": [[float(a), int(b), float(c), int(d), int(e)] for a, b, c, d, e in dg["rec_nv"]],
             "bomp_shifts": dgb["shifts"], "bomp_orth": dgb.get("orth", [])}
    return rec, extra


# =========================================================================================== part (i)
def run_part1_net(args):
    seed, regime, cell, j, out_path, n_win_max = args
    t_start = time.time()
    tag = f"seed{seed}_{regime}_{cell}_net{j}"
    kern = check_kernel(1, seed, regime, cell, j, tag)
    with open(os.path.join(PARTS_2A, tag + ".json")) as fh:
        P2A = json.load(fh)
    ci, nst, nclut, STA, statics, cat, drawn, corr = B.setup_net(seed, regime, cell, j)
    ctl = {"STA": STA.tolist() == P2A["STA"], "statics": statics.tolist() == P2A["statics"],
           "cat_drawn": drawn == P2A["cat_drawn"], "cat_src": [c[1].tolist() for c in cat] == P2A["cat_src"],
           "cat_tau": [c[2] for c in cat] == P2A["cat_tau"], "corr": corr.tolist() == P2A["corr"]}
    if not all(ctl.values()):
        fail(f"{tag}: red/catalogo/correcciones distintos de 2A: {ctl}")
    t0 = time.time()
    grid, Gs, Qs, BD, ST, PD = ops(STA, corr)
    t_ops = time.time() - t0
    rngs = T.window_rngs(seed, cell, j, regime)
    n_win = min(P2A["n_win"], n_win_max or P2A["n_win"])
    W, n_ok = [], 0
    tim = {m: [] for m in METHODS_I}
    tim["control MT-FW+SC"] = []
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
        rec, extra = methods_window(STA, grid, Gs, Qs, BD, ST, PD, D, MFn, sd, tr, corr, tim)
        W.append({"k": len(srcs), "truth": truth, "taus": taus_l, "records": rec, **extra})
    res = {"part": 1, "seed": seed, "regime": regime, "cell": cell, "net": j, "n_win": n_win, "nst": nst,
           "control": {**ctl, "windows_truth_taus": n_win, "mtfw_sc_exact": n_ok}, "corr": corr.tolist(),
           "kernel_required": kern, "kernel": C.kernel_info(), "windows": W, "timing": tim, "t_ops": t_ops,
           "seconds": time.time() - t_start}
    if out_path:
        write_json(out_path, res)
    return res


# ========================================================================================== part (ii)
def empty_rngs(cell, j, regime):
    ci = T.CELLS.index(cell)
    ri = T.REGIMES.index(regime)
    b = [T.CAL_SEED, TAG_EMPTY, ci, j, ri]
    return (np.random.default_rng(b + [41]), np.random.default_rng(b + [233]), np.random.default_rng(b + [61]))


def run_part2_net(args):
    regime, cell, j, n_empty, out_path = args
    seed = T.CAL_SEED
    t_start = time.time()
    tag = f"seed{seed}_{regime}_{cell}_net{j}"
    kern = check_kernel(2, seed, regime, cell, j, tag)
    with open(os.path.join(PARTS_2A, tag + ".json")) as fh:
        P2A = json.load(fh)
    ci, nst, nclut, STA, statics, cat, drawn, corr = B.setup_net(seed, regime, cell, j)
    ctl = {"STA": STA.tolist() == P2A["STA"], "statics": statics.tolist() == P2A["statics"],
           "cat_drawn": drawn == P2A["cat_drawn"], "corr_exact": corr.tolist() == P2A["corr"],
           "corr_maxabs": float(np.max(np.abs(corr - np.asarray(P2A["corr"]))))}
    if not (ctl["STA"] and ctl["statics"] and ctl["cat_drawn"] and ctl["corr_maxabs"] < 1e-12):
        fail(f"{tag}: red de calibracion distinta de 2A: {ctl}")
    t0 = time.time()
    grid = T.default_grid()
    garr = np.asarray(grid)
    Gs = [X.radiation_operator(g, STA) for g in garr]
    Qs = [X.basis(G) for G in Gs]
    qmc = T.load_qm_cfgs()
    lut = QM.build_lut(STA, T.VP, T.VS, {**QM.QM_CFG, **qmc["QM-3C"]})
    t_ops = time.time() - t0
    tmpdir = os.path.join(T.SCRATCH, f"t2d_s{seed}_{regime}_{cell}_{j}")
    tim = {"MT-FW+SC": [], "QM-3C+SC": [], "reproduction check": []}
    # control: la tuberia reproduce MT-FW+SC y QM-3C+SC de 2A en las ventanas 0-1 de la red
    rngs = T.window_rngs(seed, cell, j, regime)
    rep = []
    for w in range(N_CONTROL):
        t0 = time.time()
        srcs, taus, D, _p, tr, _a = T.generate_window(STA, statics, nclut, T.CLUT_AMP[regime], rngs)
        st = P2A["windows"][w]
        same_truth = [list(map(float, s)) for s in srcs] == st["truth"]
        MF, sd = T.matched(D)
        MFn = MF / sd
        r_fw = as_rec(T.certified_fw_record(STA, MFn, grid, T.DEFAULT_TAUS, Gs, Qs, tr, tcorr=corr))
        evs = T.run_qm_variants(STA, D, tmpdir, {"QM-3C+SC": (qmc["QM-3C"], corr)}, lut=lut)
        r_qm = as_rec(QM.qm_records(evs["QM-3C+SC"], tr))
        rep.append({"truth": same_truth, "MT-FW+SC": r_fw == st["records"]["MT-FW+SC"],
                    "QM-3C+SC": r_qm == st["records"]["QM-3C+SC"]})
        tim["reproduction check"].append(time.time() - t0)
    ctl["reproduces_2a"] = rep
    rngs_e = empty_rngs(cell, j, regime)
    W = []
    for w in range(n_empty):
        D, _p, _tr, _a = T.generate_real(STA, [], [], nclut, rngs_e[0], rngs_e[1], rngs_e[2], T.CLUT_AMP[regime],
                                         statics=statics)
        MF, sd = T.matched(D)
        MFn = MF / sd
        t0 = time.time()
        r_fw = as_rec(T.certified_fw_record(STA, MFn, grid, T.DEFAULT_TAUS, Gs, Qs, [], tcorr=corr))
        tim["MT-FW+SC"].append(time.time() - t0)
        t0 = time.time()
        evs = T.run_qm_variants(STA, D, tmpdir, {"QM-3C+SC": (qmc["QM-3C"], corr)}, lut=lut)
        r_qm = as_rec(QM.qm_records(evs["QM-3C+SC"], []))
        tim["QM-3C+SC"].append(time.time() - t0)
        W.append({"records": {"MT-FW+SC": r_fw, "QM-3C+SC": r_qm},
                  "max": {"MT-FW+SC": (r_fw[0][0] if r_fw else None), "QM-3C+SC": (r_qm[0][0] if r_qm else None)},
                  "qm_events": [[float(x) for x in e] for e in evs["QM-3C+SC"]]})
    res = {"part": 2, "seed": seed, "regime": regime, "cell": cell, "net": j, "n_empty": n_empty, "nst": nst,
           "seeds_empty": [seed, TAG_EMPTY, T.CELLS.index(cell), j, T.REGIMES.index(regime)],
           "STA": STA.tolist(), "statics": statics.tolist(), "corr": corr.tolist(), "control": ctl,
           "kernel_required": kern, "kernel": C.kernel_info(), "windows": W, "timing": tim, "t_ops": t_ops,
           "seconds": time.time() - t_start}
    if out_path:
        write_json(out_path, res)
    return res


# ============================================================================================ runner
def summarize(r):
    k = [w["k"] for w in r["windows"]]
    ra = {m: T.summarize_records([w["records"][m] for w in r["windows"]], k) for m in METHODS_I}
    return " | ".join(f"{m} " + "/".join(f"{ra[m][str(c)]:.2f}" for c in T.CAPS) for m in METHODS_I)


def run_all(part, workers, n_empty=N_EMPTY_PER_NET):
    os.makedirs(PARTS, exist_ok=True)
    if os.path.exists(FAIL_FILE):
        raise SystemExit(f"{FAIL_FILE} existe: el control fallo antes; revisar antes de reanudar")
    jobs = []
    if part == 1:
        for seed in SEEDS:
            for regime in T.REGIMES:
                for cell in T.CELLS:
                    for j in range(T.N_NETS):
                        if C.required_kernel(1, seed, regime, cell, j) != os.environ.get("OPENBLAS_CORETYPE"):
                            continue
                        path = os.path.join(PARTS, f"p1_seed{seed}_{regime}_{cell}_net{j}.json")
                        if not os.path.exists(path):
                            jobs.append((seed, regime, cell, j, path, None))
        jobs.sort(key=lambda a: (-T.SCEN[a[2]][1], a[0], a[1], a[3]))
    else:
        if os.environ.get("OPENBLAS_CORETYPE") != "Prescott":
            raise SystemExit("la parte (ii) (semilla 999) se corre con OPENBLAS_CORETYPE=Prescott")
        for regime in T.REGIMES:
            for cell in T.CELLS:
                for j in range(T.N_NETS):
                    path = os.path.join(PARTS, f"p2_seed999_{regime}_{cell}_net{j}.json")
                    if not os.path.exists(path):
                        jobs.append((regime, cell, j, n_empty, path))
        jobs.sort(key=lambda a: (-T.SCEN[a[1]][1], a[0], a[2]))
    print(f"parte {part}: {len(jobs)} redes pendientes; {C.kernel_info()}", flush=True)
    t0 = time.time()
    fn = run_part1_net if part == 1 else run_part2_net
    with ProcessPoolExecutor(workers) as ex:
        for r in ex.map(fn, jobs):
            if part == 1:
                print(f"   hecho p1 {r['seed']} {r['regime']} {r['cell']} red {r['net']} en {r['seconds']:.0f} s "
                      f"(total {time.time() - t0:.0f} s); control {r['control']} :: {summarize(r)}", flush=True)
            else:
                mx = [w["max"]["MT-FW+SC"] for w in r["windows"]]
                nq = sum(w["max"]["QM-3C+SC"] is not None for w in r["windows"])
                print(f"   hecho p2 999 {r['regime']} {r['cell']} red {r['net']} en {r['seconds']:.0f} s "
                      f"(total {time.time() - t0:.0f} s); control {r['control']}; {len(mx)} ventanas vacias; "
                      f"QM-3C+SC con eventos en {nq}", flush=True)


# ========================================================================================= smoke test
def smoke(nwin, workers):
    """semilla 999: parte (i) en la red 0 (con el control contra las partes de 2A de la semilla 999) y parte (ii) con
    nwin ventanas vacias en la red 0, por celda y regimen; estimacion del tiempo total."""
    jobs1 = [(T.CAL_SEED, reg, cell, 0, None, nwin) for reg in T.REGIMES for cell in T.CELLS]
    jobs2 = [(reg, cell, 0, nwin, None) for reg in T.REGIMES for cell in T.CELLS]
    t0 = time.time()
    with ProcessPoolExecutor(workers) as ex:
        r1 = list(ex.map(run_part1_net, jobs1))
        r2 = list(ex.map(run_part2_net, jobs2))
    wall = time.time() - t0
    L = [f"prueba de humo 2D: semilla 999, {nwin} ventanas por celda y regimen (red 0), {workers} procesos, {wall:.0f} s "
         f"de pared; {C.kernel_info()}"]
    est1 = est2 = 0.0
    for r in r1:
        per = {m: float(np.mean(v)) for m, v in r["timing"].items() if v}
        pw = sum(per[m] for m in METHODS_I)
        ctl = float(np.sum(r["timing"]["control MT-FW+SC"]))
        est1 += len(SEEDS) * T.N_NETS * (T.N_WIN * pw + ctl + r["t_ops"] + 3.0)
        gl = [w["gl"] for w in r["windows"]]
        L.append(f"   (i)  {r['regime']:6s} {r['cell']:15s} por ventana {pw:6.2f} s [" + ", ".join(
            f"{m} {v:.2f}" for m, v in per.items()) + f"]; operadores {r['t_ops']:.1f} s; control {r['control']}; "
            f"GL-PK: lambdas " + ",".join(str(g["n_lam_solved"]) for g in gl) + "; pasos FW "
            + ",".join(str(sum(p[5] for p in g["per_lambda"])) for g in gl)
            + f"; max brecha rel {max(p[2] for g in gl for p in g['per_lambda']):.1e}"
            f"; max c/lambda {max(p[6] for g in gl for p in g['per_lambda']):.9f} :: {summarize(r)}")
    for r in r2:
        per = {m: float(np.mean(v)) for m, v in r["timing"].items() if v and m != "reproduction check"}
        pw = sum(per.values())
        rc = float(np.sum(r["timing"]["reproduction check"]))
        est2 += T.N_NETS * (N_EMPTY_PER_NET * pw + rc + r["t_ops"] + 3.0)
        L.append(f"   (ii) {r['regime']:6s} {r['cell']:15s} por ventana vacia {pw:6.2f} s [" + ", ".join(
            f"{m} {v:.2f}" for m, v in per.items()) + f"]; control {r['control']}; maximos MT-FW+SC "
            + ",".join(f"{w['max']['MT-FW+SC']:.1f}" for w in r["windows"]) + "; QM-3C+SC "
            + ",".join("-" if w["max"]["QM-3C+SC"] is None else f"{w['max']['QM-3C+SC']:.2f}" for w in r["windows"]))
    L.append(f"   estimacion parte (i): {est1 / 3600:.2f} h de CPU; parte (ii) con 200 ventanas vacias por celda y regimen: "
             f"{est2 / 3600:.2f} h de CPU; total {(est1 + est2) / 3600:.2f} h de CPU, {(est1 + est2) / 7200:.2f} h en 2 "
             f"nucleos (umbral de la regla: 6 h)")
    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(HERE, "task2d_smoke.txt"), "w") as fh:
        fh.write(txt + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["validate", "smoke", "run"])
    ap.add_argument("--part", type=int, choices=[1, 2])
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--nwin", type=int, default=2)
    ap.add_argument("--n-empty", type=int, default=N_EMPTY_PER_NET)
    a = ap.parse_args()
    os.makedirs(T.SCRATCH, exist_ok=True)
    if a.stage == "validate":
        import task2d_validate as V
        raise SystemExit(0 if V.main() else 1)
    if a.stage == "smoke":
        smoke(a.nwin, a.workers)
        return
    run_all(a.part, a.workers, a.n_empty)


if __name__ == "__main__":
    main()
