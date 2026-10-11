"""
task2e_run.py - Tarea 2E (PREREGISTRO_1j_1k.txt, seccion "Tarea 2E").  Metodos en task2e_methods.py.

Datos: por semilla (2024, 7331), regimen y celda, las 30 redes de 2C: redes 0-4 de 2A (20 ventanas) y redes 5-29 de 2C
parte (ii) (12 ventanas), regeneradas bit a bit (task2b_ablation.setup_net, task2a_confirm.window_rngs /
generate_window) con el nucleo de OpenBLAS de 2C (task2c_run.required_kernel: SkylakeX en las 18 redes de 7331 de 2A
listadas en task2b_report.txt, Prescott en las demas), comprobado con threadpoolctl y guardado en cada parte.
Control por red: STA, retardos, (catalogo) y correcciones iguales a las partes guardadas (task2a_parts para las redes
0-4, task2c_parts/p2_* para las 5-29); verdad, origenes y k iguales en cada ventana; MT-FW+SC recalculado
(certified_fw_record) EXACTAMENTE igual al guardado en las ventanas 0 y 1 (si no, la corrida se detiene).
Metodos por ventana: TS-SFW+SC (f congelado en task2e_frozen.json), BOMP-PK-SLIDE+SC (lambda = 0), BOMP-PK+SC
(task2d_methods.bomp_pk_record sin cambios; en las redes 0-4 se compara con task2d_parts) y MT-FW+SC con posiciones
(task2e_methods.mtfw_sc_locs, el bucle de certified_fw_record; su registro se compara con el guardado en cada ventana;
el analisis usa los registros guardados y las posiciones solo donde coinciden).
Calibracion (semilla 999, redes 0-4 de 2A, 20 ventanas, nucleo Prescott): TS-SFW+SC con f en {0.01, 0.03, 0.1, 0.3},
BOMP-PK-SLIDE+SC, BOMP-PK+SC y MT-FW+SC con posiciones (para la estimacion de tiempo); f = el de mayor recall a tope
0.25 sobre las 600 ventanas agrupadas (una curva FROC exacta); empate a favor del menor.

    OPENBLAS_CORETYPE=Prescott python task2e_run.py validate
    OPENBLAS_CORETYPE=Prescott python task2e_run.py calibrate --workers 2
    python task2e_run.py cal-summary                       # -> task2e_calibration.txt, task2e_frozen.json (f)
    python task2e_run.py estimate                          # -> task2e_smoke.txt (tiempos y regla de presupuesto)
    OPENBLAS_CORETYPE=SkylakeX python task2e_run.py run --workers 2
    OPENBLAS_CORETYPE=Prescott python task2e_run.py run --workers 2 [--nets-ii 25]
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse
import hashlib
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import task2c_run as C  # noqa: E402
import task2d_run as R2  # noqa: E402
import task2e_methods as E  # noqa: E402
from analyze_2c import HF, point  # noqa: E402

T = R2.T
B = R2.B
D2 = R2.D2
SEEDS = (2024, 7331)
FGRID = (0.01, 0.03, 0.1, 0.3)
PARTS = os.path.join(HERE, "task2e_parts")
PARTS_2A = T.PARTS
PARTS_2C = os.path.join(HERE, "task2c_parts")
PARTS_2D = os.path.join(HERE, "task2d_parts")
FAIL_FILE = os.path.join(HERE, "task2e_CONTROL_FAILED.txt")
FROZEN = os.path.join(HERE, "task2e_frozen.json")
REGISTRY = os.path.join(HERE, "PREREGISTRO_1j_1k.txt")
HASHES = os.path.join(HERE, "PREREGISTRO_hashes.txt")
N_CONTROL = 2
NETS_II_FULL = 25
N_WIN_II = 12
METHODS = ("TS-SFW+SC", "BOMP-PK-SLIDE+SC", "BOMP-PK+SC", "MT-FW+SC")


class ControlError(RuntimeError):
    pass


def fail(msg):
    with open(FAIL_FILE, "a") as fh:
        fh.write(time.strftime("%Y-%m-%d %H:%M:%S ") + msg + "\n")
    raise ControlError(msg)


def as_rec(r):
    return [[float(s), int(h)] for s, h in r]


def write_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    os.replace(tmp, path)


def load(path):
    with open(path) as fh:
        return json.load(fh)


def part_kind(j):
    return 1 if j < T.N_NETS else 2


def check_kernel(seed, regime, cell, j, tag):
    req = C.required_kernel(part_kind(j), seed, regime, cell, j) if seed != T.CAL_SEED else "Prescott"
    if C.numpy_blas_arch() != [req] or os.environ.get("OPENBLAS_CORETYPE") != req:
        fail(f"{tag}: nucleo OpenBLAS {C.kernel_info()} distinto del requerido {req}")
    return req


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


# ============================================================================================ por ventana
def compact_iters(dg):
    out = []
    for it in dg["iters"]:
        r = it.get("refine", {})
        out.append([it["s"], it["node"], it["tau"], it["alternations"], it["local_nit"], it["local_status"],
                    it["max_rel_gap"], it["final_rel_gap"], it["final_kkt"], it["J"], r.get("nit"), r.get("E0"),
                    r.get("E1")])
    return out


def run_ts(net, win, tr, f=None, lam=None):
    dg = {}
    t0 = time.time()
    rec = as_rec(E.tssfw_record(net, win, tr, f=f, lam=lam, diag=dg))
    sec = time.time() - t0
    P = dg["P"]
    dl = P[:, 4:]
    return rec, {"lam": dg["lam"], "pos": P[:, :4].tolist(), "shifts": np.round(dl, 6).tolist(),
                 "m_norm": [float(np.linalg.norm(x)) for x in dg["m"]], "iters": compact_iters(dg),
                 "frac_shift_at_bound": float(np.mean(np.abs(dl) >= E.SHIFT_MAX - 1e-6)) if dl.size else 0.0}, sec


def methods_window(net, PD, win, tr, STA, grid, Gs, Qs, corr, fvals, tim):
    rec, loc, dia = {}, {}, {}
    for f in fvals:
        name = "TS-SFW+SC" if len(fvals) == 1 else f"TS-SFW+SC f={f}"
        rec[name], dia[name], sec = run_ts(net, win, tr, f=f)
        loc[name] = dia[name]["pos"]
        tim.setdefault(name, []).append(sec)
    rec["BOMP-PK-SLIDE+SC"], dia["BOMP-PK-SLIDE+SC"], sec = run_ts(net, win, tr, lam=0.0)
    loc["BOMP-PK-SLIDE+SC"] = dia["BOMP-PK-SLIDE+SC"]["pos"]
    tim.setdefault("BOMP-PK-SLIDE+SC", []).append(sec)
    t0 = time.time()
    dgb = {}
    rec["BOMP-PK+SC"] = as_rec(D2.bomp_pk_record(PD, win.MFn, tr, diag=dgb))
    tim.setdefault("BOMP-PK+SC", []).append(time.time() - t0)
    loc["BOMP-PK+SC"] = [net.garr[a[0]].tolist() + [float(net.taus[a[1]])] for a in dgb["atoms"]]
    t0 = time.time()
    r_m, l_m = E.mtfw_sc_locs(STA, win.MFn, grid, T.DEFAULT_TAUS, Gs, Qs, tr, corr)
    rec["MT-FW+SC(loc)"] = as_rec(r_m)
    loc["MT-FW+SC(loc)"] = l_m
    tim.setdefault("MT-FW+SC(loc)", []).append(time.time() - t0)
    return rec, loc, dia


# ============================================================================================ por red
def run_net(args):
    seed, regime, cell, j, fvals, out_path, n_win_max = args
    t_start = time.time()
    tag = f"seed{seed}_{regime}_{cell}_net{j}"
    kern = check_kernel(seed, regime, cell, j, tag)
    ci, nst, nclut, STA, statics, cat, drawn, corr = B.setup_net(seed, regime, cell, j)
    if j < T.N_NETS:
        ref = load(os.path.join(PARTS_2A, tag + ".json"))
        ctl = {"STA": STA.tolist() == ref["STA"], "statics": statics.tolist() == ref["statics"],
               "cat_drawn": drawn == ref["cat_drawn"], "cat_src": [c[1].tolist() for c in cat] == ref["cat_src"],
               "cat_tau": [c[2] for c in cat] == ref["cat_tau"], "corr": corr.tolist() == ref["corr"]}
        n_win = ref["n_win"]
        ref_d = load(os.path.join(PARTS_2D, f"p1_{tag}.json")) if seed != T.CAL_SEED else None
    else:
        ref = load(os.path.join(PARTS_2C, f"p2_{tag}.json"))
        ctl = {"STA": STA.tolist() == ref["STA"], "statics": statics.tolist() == ref["statics"],
               "cat_drawn": drawn == ref["cat_drawn"], "corr": corr.tolist() == ref["corr"]}
        n_win = ref["n_win"]
        ref_d = None
    if not all(ctl.values()):
        fail(f"{tag}: red/catalogo/correcciones distintos de lo guardado: {ctl}")
    n_win = min(n_win, n_win_max or n_win)
    t0 = time.time()
    grid, Gs, Qs, BD, ST, PD = R2.ops(STA, corr)
    net = E.Net(STA, corr, grid, T.DEFAULT_TAUS, Gs, Qs, BD)
    t_ops = time.time() - t0
    rngs = T.window_rngs(seed, cell, j, regime)
    W = []
    tim = {"control MT-FW+SC": []}
    n_ok, n_loc_eq, n_bpk_eq = 0, 0, 0
    for w in range(n_win):
        srcs, taus, D, _p, tr, _a = T.generate_window(STA, statics, nclut, T.CLUT_AMP[regime], rngs)
        st = ref["windows"][w]
        truth = [list(map(float, s)) for s in srcs]
        taus_l = [float(t) for t in taus]
        if len(srcs) != st["k"] or truth != st["truth"] or taus_l != st["taus"]:
            fail(f"{tag} ventana {w}: fuentes u origenes distintos de lo guardado")
        MF, sd = T.matched(D)
        MFn = MF / sd
        stored = st["records"]["MT-FW+SC"]
        if w < N_CONTROL:
            t0 = time.time()
            rs = as_rec(T.certified_fw_record(STA, MFn, grid, T.DEFAULT_TAUS, Gs, Qs, tr, tcorr=corr))
            tim["control MT-FW+SC"].append(time.time() - t0)
            if rs != stored:
                fail(f"{tag} ventana {w}: MT-FW+SC recalculado {rs} != guardado {stored}")
            n_ok += 1
        win = E.Win(D, MFn, sd)
        rec, loc, dia = methods_window(net, PD, win, tr, STA, grid, Gs, Qs, corr, fvals, tim)
        loc_eq = rec["MT-FW+SC(loc)"] == stored
        n_loc_eq += loc_eq
        rec["MT-FW+SC"] = stored
        bpk_eq = None
        if ref_d is not None:
            bpk_eq = rec["BOMP-PK+SC"] == ref_d["windows"][w]["records"]["BOMP-PK+SC"]
            n_bpk_eq += bpk_eq
        W.append({"k": len(srcs), "truth": truth, "taus": taus_l, "records": rec, "loc": loc, "diag": dia,
                  "mtfw_loc_eq_stored": loc_eq, "bomp_pk_eq_2d": bpk_eq})
    res = {"seed": seed, "regime": regime, "cell": cell, "net": j, "n_win": n_win, "nst": nst, "fvals": list(fvals),
           "control": {**ctl, "windows_truth_taus": n_win, "mtfw_sc_exact": n_ok, "mtfw_loc_eq_stored": n_loc_eq,
                       "bomp_pk_eq_2d": (n_bpk_eq if ref_d is not None else None)},
           "corr": corr.tolist(), "kernel_required": kern, "kernel": C.kernel_info(), "windows": W, "timing": tim,
           "t_ops": t_ops, "seconds": time.time() - t_start}
    if out_path:
        write_json(out_path, res)
    return res


# ============================================================================================ corridas
def summarize(r):
    k = [w["k"] for w in r["windows"]]
    ms = [m for m in r["windows"][0]["records"] if m != "MT-FW+SC(loc)"]
    ra = {m: T.summarize_records([w["records"][m] for w in r["windows"]], k) for m in ms}
    return " | ".join(f"{m} " + "/".join(f"{ra[m][str(c)]:.2f}" for c in T.CAPS) for m in ms)


def run_jobs(jobs, workers, label):
    print(f"{label}: {len(jobs)} redes pendientes; {C.kernel_info()}", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(workers) as ex:
        for r in ex.map(run_net, jobs):
            print(f"   hecho {r['seed']} {r['regime']} {r['cell']} red {r['net']} en {r['seconds']:.0f} s (total "
                  f"{time.time() - t0:.0f} s); control {r['control']} :: {summarize(r)}", flush=True)


def calibrate(workers):
    os.makedirs(PARTS, exist_ok=True)
    if os.path.exists(FAIL_FILE):
        raise SystemExit(f"{FAIL_FILE} existe")
    if os.environ.get("OPENBLAS_CORETYPE") != "Prescott":
        raise SystemExit("la calibracion (semilla 999) se corre con OPENBLAS_CORETYPE=Prescott")
    jobs = []
    for regime in T.REGIMES:
        for cell in T.CELLS:
            for j in range(T.N_NETS):
                path = os.path.join(PARTS, f"cal_seed999_{regime}_{cell}_net{j}.json")
                if not os.path.exists(path):
                    jobs.append((T.CAL_SEED, regime, cell, j, FGRID, path, None))
    jobs.sort(key=lambda a: (-T.SCEN[a[2]][1], a[1], a[3]))
    run_jobs(jobs, workers, "calibracion 999")


def cal_parts():
    out = []
    for regime in T.REGIMES:
        for cell in T.CELLS:
            for j in range(T.N_NETS):
                out.append(load(os.path.join(PARTS, f"cal_seed999_{regime}_{cell}_net{j}.json")))
    return out


def cal_summary():
    parts = cal_parts()
    L = ["# Calibracion de f de TS-SFW+SC (Tarea 2E), semilla 999, redes 0-4 de 2A, 20 ventanas por red: 600 ventanas",
         "(6 celdas x regimen x 100). Recall a <= 0.1 / <= 0.25 falsos por ventana (curva FROC exacta, punto de operacion "
         "sobre la propia calibracion). Regla: f de mayor recall a 0.25 sobre las 600 ventanas agrupadas (una curva); "
         "empate a favor del menor.", ""]
    names = [f"TS-SFW+SC f={f}" for f in FGRID] + ["BOMP-PK-SLIDE+SC", "BOMP-PK+SC", "MT-FW+SC"]
    ctl_ok = all(p["control"]["mtfw_sc_exact"] == 2 for p in parts)
    L.append(f"Controles: {len(parts)} redes; red, catalogo y correcciones iguales a 2A y MT-FW+SC exacto en las ventanas "
             f"0-1: {'correcto' if ctl_ok else 'FALLA'}; MT-FW+SC con posiciones igual al guardado en "
             f"{sum(p['control']['mtfw_loc_eq_stored'] for p in parts)}/{sum(p['n_win'] for p in parts)} ventanas; "
             f"nucleos {sorted(set(tuple(tuple(x) for x in p['kernel']['blas']) for p in parts))}")
    L.append("")
    L.append("| Celda (regimen) | " + " | ".join(names) + " |")
    L.append("|:--|" + ":-:|" * len(names))
    pooled = {m: [] for m in names}
    kall = []
    per_cell = {m: [] for m in names}
    for regime in T.REGIMES:
        for cell in T.CELLS:
            ws = [w for p in parts if p["regime"] == regime and p["cell"] == cell for w in p["windows"]]
            k = np.array([w["k"] for w in ws], float)
            kall.extend(k.tolist())
            row = []
            for m in names:
                recs = [w["records"][m] for w in ws]
                pooled[m].extend(recs)
                pt = point(HF(recs), k)
                per_cell[m].append(pt["0.25"])
                row.append(f"{pt['0.1']:.3f}/{pt['0.25']:.3f}")
            L.append(f"| {cell} ({regime}) | " + " | ".join(row) + " |")
    kall = np.array(kall)
    pts = {m: point(HF(pooled[m]), kall) for m in names}
    L.append(f"| agrupado (600) | " + " | ".join(f"{pts[m]['0.1']:.3f}/{pts[m]['0.25']:.3f}" for m in names) + " |")
    L.append(f"| media de las 6 celdas a 0.25 | " + " | ".join(f"{np.mean(per_cell[m]):.3f}" for m in names) + " |")
    best = max(FGRID, key=lambda f: (round(pts[f"TS-SFW+SC f={f}"]["0.25"], 12), -f))
    L.append("")
    L.append(f"f elegido: {best} (recall agrupado a 0.25 = {pts[f'TS-SFW+SC f={best}']['0.25']:.4f}; "
             + ", ".join(f"f={f}: {pts[f'TS-SFW+SC f={f}']['0.25']:.4f}" for f in FGRID) + ")")
    # lambda / c_1 implicito
    for f in FGRID:
        r = [w["diag"][f"TS-SFW+SC f={f}"]["lam"] / np.sqrt(2 * w["diag"][f"TS-SFW+SC f={f}"]["lam"] / f)
             for p in parts for w in p["windows"]]
        L.append(f"- f = {f}: lambda / sqrt(2 s_1) (fraccion del certificado del primer atomo) mediana {np.median(r):.3f}, "
                 f"P10 {np.percentile(r, 10):.3f}, P90 {np.percentile(r, 90):.3f}; fraccion de atomos con amplitud "
                 f"nula {np.mean([x == 0.0 for p in parts for w in p['windows'] for x in w['diag'][f'TS-SFW+SC f={f}']['m_norm']]):.3f}")
    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(HERE, "task2e_calibration.txt"), "w") as fh:
        fh.write(txt + "\n")
    frozen = {"f": best, "pooled_recall_025": {str(f): pts[f"TS-SFW+SC f={f}"]["0.25"] for f in FGRID},
              "written_utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())}
    write_json(FROZEN, frozen)


def estimate(nets_ii=NETS_II_FULL):
    parts = cal_parts()
    fr = load(FROZEN)
    fname = f"TS-SFW+SC f={fr['f']}"
    L = [f"# Estimacion de tiempo de la evaluacion de 2E desde la calibracion completa en la semilla 999 "
         f"(f = {fr['f']}; {len(parts)} redes, {sum(p['n_win'] for p in parts)} ventanas)"]
    tot = 0.0
    tot15 = 0.0
    for regime in T.REGIMES:
        for cell in T.CELLS:
            ps = [p for p in parts if p["regime"] == regime and p["cell"] == cell]
            per = {m: float(np.mean([x for p in ps for x in p["timing"][m]]))
                   for m in (fname, "BOMP-PK-SLIDE+SC", "BOMP-PK+SC", "MT-FW+SC(loc)")}
            pw = sum(per.values())
            ctl = float(np.mean([sum(p["timing"]["control MT-FW+SC"]) for p in ps]))
            # 'seconds' incluye los 4 valores de f; el costo por red sin metodos = seconds - todos los tiempos medidos
            setup = float(np.mean([p["seconds"] - sum(sum(v) for v in p["timing"].values()) for p in ps]))
            nets_ev = T.N_NETS * T.N_WIN * pw + nets_ii * N_WIN_II * pw
            per_net = ctl + setup
            c = len(SEEDS) * (nets_ev + (T.N_NETS + nets_ii) * per_net)
            c15 = len(SEEDS) * (T.N_NETS * T.N_WIN * pw + 10 * N_WIN_II * pw + (T.N_NETS + 10) * per_net)
            tot += c
            tot15 += c15
            L.append(f"   {regime:6s} {cell:15s} por ventana {pw:6.2f} s [" + ", ".join(f"{m} {v:.2f}" for m, v in per.items())
                     + f"]; por red: control {ctl:.1f} s, preparacion {setup:.1f} s; celda {c / 3600:.2f} h de CPU")
    L.append(f"   estimacion total (30 redes por celda): {tot / 3600:.2f} h de CPU, {tot / 7200:.2f} h en 2 nucleos (umbral de "
             f"la regla de presupuesto: 10 h); con 15 redes: {tot15 / 7200:.2f} h en 2 nucleos")
    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(HERE, "task2e_smoke.txt"), "w") as fh:
        fh.write(txt + "\n")


def run_eval(workers, nets_ii):
    os.makedirs(PARTS, exist_ok=True)
    if os.path.exists(FAIL_FILE):
        raise SystemExit(f"{FAIL_FILE} existe: el control fallo antes; revisar antes de reanudar")
    fr = load(FROZEN)
    reg_hash = sha256(REGISTRY)
    last = [ln.split()[0] for ln in open(HASHES) if ln.strip().endswith("PREREGISTRO_1j_1k.txt")][-1]
    if reg_hash != last or "Nota de calibracion de 2E" not in open(REGISTRY).read():
        raise SystemExit("el registro no tiene la nota de calibracion de 2E con su hash: no se corre la evaluacion")
    if f"f = {fr['f']}" not in open(REGISTRY).read().split("Nota de calibracion de 2E")[1][:600]:
        raise SystemExit("el f congelado no coincide con la nota del registro")
    jobs = []
    for seed in SEEDS:
        for regime in T.REGIMES:
            for cell in T.CELLS:
                for j in list(range(T.N_NETS)) + list(range(T.N_NETS, T.N_NETS + nets_ii)):
                    if C.required_kernel(part_kind(j), seed, regime, cell, j) != os.environ.get("OPENBLAS_CORETYPE"):
                        continue
                    path = os.path.join(PARTS, f"ev_seed{seed}_{regime}_{cell}_net{j}.json")
                    if not os.path.exists(path):
                        jobs.append((seed, regime, cell, j, (fr["f"],), path, None))
    jobs.sort(key=lambda a: (-T.SCEN[a[2]][1], a[0], a[1], a[3]))
    run_jobs(jobs, workers, f"evaluacion (f = {fr['f']})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["validate", "calibrate", "cal-summary", "estimate", "run"])
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--nets-ii", type=int, default=NETS_II_FULL)
    a = ap.parse_args()
    os.makedirs(T.SCRATCH, exist_ok=True)
    if a.stage == "validate":
        import task2e_validate as V
        raise SystemExit(0 if V.main() else 1)
    if a.stage == "calibrate":
        calibrate(a.workers)
    elif a.stage == "cal-summary":
        cal_summary()
    elif a.stage == "estimate":
        estimate(a.nets_ii)
    else:
        run_eval(a.workers, a.nets_ii)


if __name__ == "__main__":
    main()
