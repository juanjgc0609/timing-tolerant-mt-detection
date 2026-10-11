"""
task2l_run.py - Tarea 2L (PREREGISTRO_1j_1k.txt, seccion "Tarea 2L"): terminos de estacion por determinacion conjunta
de hipocentros y terminos sin picks (S2, JHD sobre el catalogo registrado de 20 eventos) dentro del detector refinado
(T3).  Corredor de las semillas 7171 y 9393 y de los controles en la semilla 999.

Diseno implementado (las definiciones son las del piloto C11, sin cambios: explore_contrib/c11_run.py, c11_lib.py,
c10_lib.py, c10_run.py, c9_run.py, c2_run.py, xc.py; la semilla se pasa fijando en el proceso de cada trabajo las
variables SEED de xc, c9_lib, c10_lib, c11_lib y c8_lib, de las que dependen xc.NetCfg (red, retardos, catalogo, +SC) y
c9_run.stream (ventanas)):
  redes         task2a_confirm.network(seed, cell, j); retardos N(0, 0.3 s) por estacion y fase, fijos por red
  catalogo      xc.NetCfg: 20 eventos, flujos [seed, 2003, ci, j, ri, 0..3] (= task2a_confirm.catalogue_real), errores
                1 km / 0.1 s; gen_realistic por defecto
  ventanas      c9_run.stream (c8_lib.generate_site = gen_realistic.generate_real, factores de sitio sorteados por
                ventana): flujos [seed, TAG, ci, j, ri, 41|233|61]
                  TAG 6101  ventanas estandar (xc.draw_sources = task1o_consistent.draw_sources, 0-2 eventos), kmax 5
                  TAG 6102  ventanas densas (c2_run.draw_dense, 3-5 eventos), kmax 7
                6101 y 6102 no aparecen en ningun archivo anterior (grep de todos los TAG y entropias de rng: usados
                2001-2005, 2097, 2098, 24xx, 26xx, 27xx, 2801, 3000-3008, 3701-3706, 3801-3812, 3901-3913, 4001-4004,
                4101-4105, 4201, 4202, 5201, 5202).
  terminos      S0 = +SC registrado (xc.NetCfg.corr = station_corrections.estimate_station_corrections sobre el
                catalogo); S2 = c11_lib.s2_terms(STA, catalogo, S0, energias, c11_lib.P) (JHD, constantes c11_lib.P del
                plan congelado C11_PLAN.txt); ORA = retardos verdaderos (brazo oraculo)
  detectores    c11_run.Net.run (bucles c10_lib.mt_loop tope 6 y c10_lib.bomp_loop; rejilla 3 km x 3 profundidades,
                origenes 5.0-8.0 s paso 0.1, regla dura +-0.25 s, exclusion 5 km); T3 = refinamiento de C10 (rejilla
                de 1 km, radio 10 km, origenes a 1 s); nivel de tiempo de T3 = c10_lib.sigma_catalogue recalculado con
                cada conjunto de terminos; ORA: c10_lib.sigma_oracle
  brazos (por red j, celda, regimen y semilla; una parte por combinacion)
                estandar   j = 0..159: MT@S0, MT@S0+T3, MT@S2+T3           (5 ventanas)
                           j = 0..59 : ademas BOMP@S0, BOMP@S0+T3, BOMP@S2+T3
                           j = 0..29 : ademas MT@ORA+T3 y BOMP@ORA+T3 (brazo oraculo, descriptivo)
                densas     j = 0..99 : MT@S0, MT@S0+T3, MT@S2+T3           (5 ventanas)
                           j = 0..29 : ademas MT@ORA+T3
                (regla de presupuesto, si se activa: estandar j = 0..119, densas j = 0..79; se fija con --budget en la
                orden de corrida y se anota en el registro antes de correr)
  por atomo se guarda [x, y, z, origen, s, s0, sref] (c10_run.pack); los registros (regla de acierto LOC y LOC+T) se
  construyen en analyze_2l.py con c10_lib.label / c10_lib.record.
  nucleo        OPENBLAS_CORETYPE=Prescott en todo proceso; cada trabajo comprueba con threadpoolctl que la OpenBLAS de
                numpy cargo Prescott y lo guarda en su parte.
  congelacion   el modo run exige task2l_code_hashes.txt y que sha256sum -c pase.

Controles en la semilla 999 (task2l_control.py): brazo "ctrl", redes 210, 219, 228 del conjunto de evaluacion de C11
con los TAG del piloto (4101 estandar, 4102 densas): mismas llamadas, todos los brazos (MT, BOMP en estandar, oraculo).

Partes: task2l_parts/{arm}_seed{s}_{reg}_{cell}_net{j}.json (reanudable)

    OPENBLAS_CORETYPE=Prescott python task2l_run.py control --workers 2
    OPENBLAS_CORETYPE=Prescott python task2l_run.py run --seeds 7171,9393 --workers 2 [--budget]
    python task2l_run.py status
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MPLBACKEND", "Agg")
import argparse
import json
import subprocess
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXPL = os.path.join(HERE, "explore_contrib")
for p in (EXPL, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

import xc  # noqa: E402
import c2_run as C2  # noqa: E402
import c8_lib as C8  # noqa: E402
import c9_lib as C9  # noqa: E402
import c9_run as R9  # noqa: E402
import c10_lib as C10  # noqa: E402
import c10_run as R10  # noqa: E402
import c11_lib as L  # noqa: E402
import c11_run as R11  # noqa: E402
import task2a_confirm as T  # noqa: E402

REQUIRED_KERNEL = "Prescott"
CELLS = T.CELLS
REGIMES = T.REGIMES
SEEDS = (7171, 9393)
TAG_STD, TAG_DENSE = 6101, 6102
KMAX = {"std": 5, "dense": 7}
N_WIN = 5
MT3 = ("MT@S0", "MT@S0+T3", "MT@S2+T3")
BOMP3 = ("BOMP@S0", "BOMP@S0+T3", "BOMP@S2+T3")
ORA_MT, ORA_BOMP = "MT@ORA+T3", "BOMP@ORA+T3"
# brazo: redes del diseno estandar, del denso, de la familia BOMP (estandar) y del oraculo; TAG (estandar, densas)
ARMS = {
    "main": {"std": 160, "dense": 100, "bomp": 60, "ora": 30, "tags": (TAG_STD, TAG_DENSE)},
    "budget": {"std": 120, "dense": 80, "bomp": 60, "ora": 30, "tags": (TAG_STD, TAG_DENSE)},
}
CONTROL_NETS = (210, 219, 228)
CONTROL_TAGS = (4101, 4102)
PARTS = os.path.join(HERE, "task2l_parts")
LOG = os.path.join(HERE, "log2l.txt")
FAIL_FILE = os.path.join(HERE, "task2l_FAIL.txt")
HASHES = os.path.join(HERE, "task2l_code_hashes.txt")
SEED_MODULES = (xc, C9, C10, L, C8)


# ------------------------------------------------------------------------------------------ nucleo
def blas_libs():
    import threadpoolctl
    return [(i.get("architecture"), os.path.basename(i.get("filepath", ""))) for i in threadpoolctl.threadpool_info()
            if i.get("user_api") == "blas"]


def kernel_info():
    return {"OPENBLAS_CORETYPE": os.environ.get("OPENBLAS_CORETYPE"), "blas": blas_libs()}


def check_kernel():
    numpy_arch = [a for a, f in blas_libs() if "scipy" not in f]
    if numpy_arch != [REQUIRED_KERNEL] or os.environ.get("OPENBLAS_CORETYPE") != REQUIRED_KERNEL:
        msg = f"nucleo OpenBLAS {kernel_info()} distinto del requerido {REQUIRED_KERNEL}"
        with open(FAIL_FILE, "a") as fh:
            fh.write(time.strftime("%Y-%m-%d %H:%M:%S ") + msg + "\n")
        raise RuntimeError(msg)


# ------------------------------------------------------------------------------------------ trabajo
def set_seed(seed):
    for mod in SEED_MODULES:
        mod.SEED = seed


def methods_for(kind, j, spec):
    """(familia, terminos, t3, clave) de la red j en el tipo de ventana kind; [] si la red no entra en ese diseno."""
    if spec is None:                                   # control: todos los brazos
        ms = list(MT3) + [ORA_MT] + (list(BOMP3) + [ORA_BOMP] if kind == "std" else [])
    else:
        if j >= spec[kind]:
            return []
        ms = list(MT3)
        if kind == "std" and j < spec["bomp"]:
            ms += list(BOMP3)
        if j < spec["ora"]:
            ms.append(ORA_MT)
            if kind == "std" and j < spec["bomp"]:
                ms.append(ORA_BOMP)
    out = []
    for key in ms:
        fam, rest = key.split("@")
        out.append((fam, rest.replace("+T3", ""), rest.endswith("+T3"), key))
    return out


def part_path(arm, seed, regime, cell, j):
    name = "ctrl" if arm == "ctrl" else "main"          # el brazo budget escribe las mismas partes que main
    return os.path.join(PARTS, f"{name}_seed{seed}_{regime}_{cell}_net{j}.json")


def write_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    os.replace(tmp, path)


def run_job(args):
    arm, seed, regime, cell, j = args
    out = part_path(arm, seed, regime, cell, j)
    if os.path.exists(out):
        return out, "exists", 0.0
    check_kernel()
    if (arm in ARMS) != (seed in SEEDS) or (arm == "ctrl") != (seed == 999):
        raise RuntimeError(f"brazo {arm} con semilla {seed} no permitido")
    spec = None if arm == "ctrl" else ARMS[arm]
    tags = CONTROL_TAGS if arm == "ctrl" else spec["tags"]
    t0 = time.perf_counter()
    set_seed(seed)
    net = R11.Net(cell, j, regime)                      # red, catalogo, S0 = +SC, motores S0 y ORA
    nc = net.nc
    STA_ref, stat_ref = T.network(seed, cell, j)
    assert np.array_equal(nc.STA, STA_ref) and np.array_equal(nc.statics, stat_ref)
    t_setup = time.perf_counter() - t0
    t = time.perf_counter()
    Ecat = [L.energy(D) for D, _s, _t in nc.cat]
    i2 = {}
    net.add("S2", L.s2_terms(nc.STA, nc.cat, nc.corr, Ecat, L.P, info=i2), nc.cat)
    t_s2 = time.perf_counter() - t
    wins = {}
    for kind, tag, draw in (("std", tags[0], xc.draw_sources), ("dense", tags[1], C2.draw_dense)):
        ms = methods_for(kind, j, spec)
        if not ms:
            continue
        st = R9.stream(nc, tag, draw, None)
        rows = []
        for wi in range(N_WIN):
            w = next(st)
            atoms, sec = {}, {}
            for fam, nm, t3, key in ms:
                t = time.perf_counter()
                atoms[key] = R10.pack(net.run(fam, nm, w["MFn"], KMAX[kind], t3))
                sec[key] = time.perf_counter() - t
            rows.append({"w": wi, "k": len(w["truth"]), "truth": [t_.tolist() for t_ in w["truth"]], "taus": w["taus"],
                         "site_logsd": float(np.log(np.asarray(w["site"])).std()), "atoms": atoms, "sec": sec})
        wins[kind] = rows
    assert all(mod.SEED == seed for mod in SEED_MODULES)
    obj = {"task": "2L", "arm": arm, "seed": seed, "regime": regime, "cell": cell, "net": j,
           "tags": {"std": tags[0], "dense": tags[1]}, "kmax": KMAX, "n_win": N_WIN,
           "methods": {kind: [m[3] for m in methods_for(kind, j, spec)] for kind in ("std", "dense")},
           "statics": np.asarray(nc.statics).tolist(), "cat_drawn": int(nc.drawn),
           "terms": {k: np.asarray(v).tolist() for k, v in net.terms.items()},
           "err": {k: L.err(v, nc.statics) for k, v in net.terms.items()},
           "off": {k: L.offset(v, nc.statics) for k, v in net.terms.items()},
           "sigT3": dict(net.sig), "S2": {k: i2[k] for k in ("sig", "n_res", "n_ev")},
           "t_setup": t_setup, "t_s2": t_s2, "t_total": time.perf_counter() - t0, "kernel": kernel_info(),
           "host_time_utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), "windows": wins}
    write_json(out, obj)
    return out, "done", obj["t_total"]


def log(msg):
    line = time.strftime("%Y-%m-%d %H:%M:%S UTC ", time.gmtime()) + msg
    print(line, flush=True)
    with open(LOG, "a") as fh:
        fh.write(line + "\n")


def run_jobs(jobs, workers):
    todo = [a for a in jobs if not os.path.exists(part_path(*a))]
    log(f"{len(todo)} de {len(jobs)} trabajos pendientes; {kernel_info()}")
    if not todo:
        return
    t0 = time.time()
    n = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(run_job, a): a for a in todo}
        for f in as_completed(futs):
            a = futs[f]
            try:
                path, status, tt = f.result()
            except Exception:
                log(f"ERROR {a}: {traceback.format_exc()}")
                with open(FAIL_FILE, "a") as fh:
                    fh.write(time.strftime("%Y-%m-%d %H:%M:%S ") + f"{a}\n{traceback.format_exc()}\n")
                continue
            n += 1
            el = time.time() - t0
            log(f"{status} {os.path.basename(path)} {tt:.1f} s; {n}/{len(todo)}; "
                f"restante estimado {el / n * (len(todo) - n) / 3600:.2f} h")


def jobs_for(seeds, arm):
    """orden: red, semilla, celda, regimen (las dos semillas avanzan a la par)."""
    nets = CONTROL_NETS if arm == "ctrl" else range(max(ARMS[arm]["std"], ARMS[arm]["dense"]))
    return [(arm, seed, regime, cell, j) for j in nets for seed in seeds for cell in CELLS for regime in REGIMES]


def frozen_ok():
    if not os.path.exists(HASHES):
        return False
    r = subprocess.run(["sha256sum", "-c", "--quiet", os.path.basename(HASHES)], cwd=HERE, capture_output=True)
    return r.returncode == 0


def loaded_modules():
    return sorted({os.path.relpath(m.__file__, HERE) for m in list(sys.modules.values())
                   if getattr(m, "__file__", None) and os.path.abspath(m.__file__).startswith(HERE + os.sep)})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("control", "run", "status", "modules"))
    ap.add_argument("--seeds", default="7171,9393")
    ap.add_argument("--budget", action="store_true")
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    os.makedirs(PARTS, exist_ok=True)
    arm = "budget" if a.budget else "main"
    if a.mode == "modules":
        print("\n".join(loaded_modules()))
        return
    if a.mode == "status":
        for seed in SEEDS:
            jobs = jobs_for([seed], arm)
            print(arm, seed, sum(os.path.exists(part_path(*x)) for x in jobs), "/", len(jobs))
        return
    if os.environ.get("OPENBLAS_CORETYPE") != REQUIRED_KERNEL:
        raise SystemExit(f"se corre con OPENBLAS_CORETYPE={REQUIRED_KERNEL}")
    check_kernel()
    if a.mode == "control":
        run_jobs(jobs_for([999], "ctrl"), a.workers)
        return
    seeds = [int(x) for x in a.seeds.split(",")]
    if any(s not in SEEDS for s in seeds):
        raise SystemExit("run solo admite las semillas 7171 y 9393")
    if not frozen_ok():
        raise SystemExit("codigo no congelado: falta task2l_code_hashes.txt o sha256sum -c falla")
    run_jobs(jobs_for(seeds, arm), a.workers)


if __name__ == "__main__":
    main()
