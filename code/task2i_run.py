"""
task2i_run.py - Tarea 2I (PREREGISTRO_1j_1k.txt, seccion "Tarea 2I"): certificado de tensor de momento con ajuste de
Huber escalado por el ruido (HUBER-N) frente a BOMP-PK+SC y a MT-FW+SC con tope 6.  Corredor de las semillas 5151 y
6262 y de los controles en la semilla 999.

Diseno implementado (las definiciones son las del piloto explore_contrib/c8_run.py y c8_lib.py, sin cambios; la semilla
se pasa fijando c8_lib.SEED en el proceso de cada trabajo, unica variable de la que dependen WNet, su catalogo y sus
flujos de ventanas):
  redes         task2a_confirm.network(seed, cell, j); retardos N(0, 0.3 s) por estacion y fase, fijos por red
  mundo         c8_lib.World("iid", sd): gen_realistic por defecto (factores de sitio exp(N(0, sd)) sorteados en cada
                ventana y en cada evento de catalogo); sd 0.5 (primario) y 0.25 / 0.8 (brazos secundarios)
  catalogo, SC  c8_lib.WNet: 20 eventos, flujos [seed, 2003, ci, j, ri, 0..3], en el mismo mundo;
                station_corrections.estimate_station_corrections
  ventanas      c8_lib.WNet.stream: flujos [seed, TAG, ci, j, ri, 41|233|61]
                  TAG 5201  ventanas densas (c2_run.draw_dense, 3-5 eventos), kmax 7
                  TAG 5202  ventanas estandar (task1o_consistent.draw_sources, 0-2 eventos), kmax 5
                5201 y 5202 no aparecen en ningun archivo anterior (grep de todos los TAG y entropias de rng: usados
                2001-2005, 2097, 2098, 26xx, 27xx, 2801, 3001-3008, 3701-3706, 3801-3812, 4201, 4202).
  operadores    xc.Ops(STA, corr, xc.make_grid()) y c8_lib.DOps (rejilla 3 km x 3 profundidades, origenes DEFAULT_TAUS)
  metodos       BOMP-PK+SC      xc.bomp_pk(xops.PD, ...) (= task2d_methods.bomp_pk_record)
                MT-FW+SC        c8_lib.mtfw_general(stat="ls", capv=6.0)
                MT-FW-NOCAP+SC  c8_lib.mtfw_general(stat="ls", capv=None)
                HUBER-N+SC      c8_lib.mtfw_general(stat="huber_n")            (c = c8_lib.C_HUB = 2.7955)
                secundarios     HUBER-N90 (chub 2.5003), HUBER-N99 (chub 3.3682), HUBER-C5 (chub 5.0)
  brazos        prim   j = 0..119, sd 0.5, 5 ventanas densas + 5 estandar por red, los 4 metodos        (primario)
                thr    j = 0..19, las MISMAS ventanas que prim (mismos flujos): HUBER-N (control de identidad) y los
                       tres umbrales secundarios, ventanas densas y estandar                             (secundario)
                s25, s80   j = 0..19, sd de sitio 0.25 / 0.8, los mismos TAG y redes (mismos flujos: mismas fuentes,
                       origenes, mecanismos, ruido y sorteos normales de sitio; solo cambia la sd), catalogo y SC
                       regenerados en ese mundo, 5 densas + 5 estandar, los 4 metodos                    (secundario)
  nucleo        OPENBLAS_CORETYPE=Prescott en todo proceso; cada trabajo comprueba con threadpoolctl que la OpenBLAS de
                numpy cargo Prescott y lo guarda en su parte.
  congelacion   el modo run exige task2i_code_hashes.txt y que sha256sum -c pase.

Controles en la semilla 999 (datos del piloto C8; task2i_control.py):
  cprim  redes 60, 67, 74; TAG 3801 (DENSE-IID, 5 ventanas) y 3802 (STD-IID, 4 ventanas en el piloto)
  cthr   redes 60, 67; TAG 3801 (HSWEEP: N90, N99; HWIDE: C5) y 3802
  cs25, cs80   red 80; TAG 3805 / 3806 (4 ventanas en el piloto), solo densas
  ds25, ds80   redes 60, 67; TAG 3801 / 3802 con sd 0.25 / 0.8: sin referencia en el piloto; solo dan a analyze_2i.py
               (--dry) la misma estructura de brazos que la corrida real para ensayarlo antes de congelar

Partes (una por brazo, semilla, regimen, celda y red; reanudable): task2i_parts/{arm}_seed{s}_{reg}_{cell}_net{j}.json

    OPENBLAS_CORETYPE=Prescott python task2i_run.py control --workers 2
    OPENBLAS_CORETYPE=Prescott python task2i_run.py run --seeds 5151,6262 --arms prim --workers 2
    OPENBLAS_CORETYPE=Prescott python task2i_run.py run --seeds 5151,6262 --arms thr,s25,s80 --workers 2
    python task2i_run.py status
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
import c8_lib as C  # noqa: E402
import task2a_confirm as T  # noqa: E402

REQUIRED_KERNEL = "Prescott"
CELLS = T.CELLS
REGIMES = T.REGIMES
SEEDS = (5151, 6262)
TAG_DENSE, TAG_STD = 5201, 5202
KMAX = {"dense": 7, "std": 5}
N_WIN = 5
PRIMARY = ("BOMP-PK+SC", "MT-FW+SC", "MT-FW-NOCAP+SC", "HUBER-N+SC")
THR = ("HUBER-N+SC", "HUBER-N90+SC", "HUBER-N99+SC", "HUBER-C5+SC")
ENGINE = {
    "MT-FW+SC": dict(stat="ls", capv=6.0),
    "MT-FW-NOCAP+SC": dict(stat="ls", capv=None),
    "HUBER-N+SC": dict(stat="huber_n"),
    "HUBER-N90+SC": dict(stat="huber_n", chub=2.5003),      # raiz del cuantil 0.90 de chi-cuadrado(3)
    "HUBER-N99+SC": dict(stat="huber_n", chub=3.3682),      # raiz del cuantil 0.99
    "HUBER-C5+SC": dict(stat="huber_n", chub=5.0),
}
ARMS = {
    "prim": {"sd": 0.5, "tags": (TAG_DENSE, TAG_STD), "methods": PRIMARY, "nets": tuple(range(0, 120))},
    "thr": {"sd": 0.5, "tags": (TAG_DENSE, TAG_STD), "methods": THR, "nets": tuple(range(0, 20))},
    "s25": {"sd": 0.25, "tags": (TAG_DENSE, TAG_STD), "methods": PRIMARY, "nets": tuple(range(0, 20))},
    "s80": {"sd": 0.8, "tags": (TAG_DENSE, TAG_STD), "methods": PRIMARY, "nets": tuple(range(0, 20))},
}
CONTROL = {
    "cprim": {"sd": 0.5, "tags": (3801, 3802), "methods": PRIMARY, "nets": (60, 67, 74)},
    "cthr": {"sd": 0.5, "tags": (3801, 3802), "methods": THR, "nets": (60, 67)},
    "cs25": {"sd": 0.25, "tags": (3805, None), "methods": PRIMARY, "nets": (80,)},
    "cs80": {"sd": 0.8, "tags": (3806, None), "methods": PRIMARY, "nets": (80,)},
    "ds25": {"sd": 0.25, "tags": (3801, 3802), "methods": PRIMARY, "nets": (60, 67)},
    "ds80": {"sd": 0.8, "tags": (3801, 3802), "methods": PRIMARY, "nets": (60, 67)},
}
PARTS = os.path.join(HERE, "task2i_parts")
LOG = os.path.join(HERE, "log2i.txt")
FAIL_FILE = os.path.join(HERE, "task2i_FAIL.txt")
HASHES = os.path.join(HERE, "task2i_code_hashes.txt")


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
def spec_of(arm):
    return ARMS[arm] if arm in ARMS else CONTROL[arm]


def run_methods(xops, dops, MFn, tr, kmax, methods):
    """las llamadas de c8_run.run / c8_run.methods."""
    rec, loc, sec = {}, {}, {}
    for name in methods:
        t = time.perf_counter()
        if name == "BOMP-PK+SC":
            r, l_ = xc.bomp_pk(xops.PD, MFn, tr, kmax)
        else:
            r, l_ = C.mtfw_general(dops, MFn, tr, kmax, **ENGINE[name])
        sec[name] = time.perf_counter() - t
        rec[name] = [[float(s), int(h)] for s, h in r]
        loc[name] = l_
    return rec, loc, sec


def part_path(arm, seed, regime, cell, j):
    return os.path.join(PARTS, f"{arm}_seed{seed}_{regime}_{cell}_net{j}.json")


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
    spec = spec_of(arm)
    if (arm in ARMS) != (seed in SEEDS) or (arm in CONTROL) != (seed == 999):
        raise RuntimeError(f"brazo {arm} con semilla {seed} no permitido")
    t0 = time.perf_counter()
    C.SEED = seed                                   # la unica dependencia de la semilla en c8_lib.WNet
    world = C.World("iid", spec["sd"])
    net = C.WNet(cell, j, regime, world)
    xops = xc.Ops(net.STA, net.corr, xc.make_grid())
    dops = C.DOps(xops)
    t_setup = time.perf_counter() - t0
    wins = {}
    for kind, tag, draw in (("dense", spec["tags"][0], C2.draw_dense), ("std", spec["tags"][1], xc.draw_sources)):
        if tag is None:
            continue
        st = net.stream(tag, draw)
        rows = []
        for wi in range(N_WIN):
            w = next(st)
            rec, loc, sec = run_methods(xops, dops, w["MFn"], w["truth"], KMAX[kind], spec["methods"])
            rows.append({"w": wi, "k": len(w["truth"]), "truth": [t.tolist() for t in w["truth"]], "taus": w["taus"],
                         "site_logsd": float(np.log(np.asarray(w["site"])).std()),
                         "records": rec, "loc": loc, "sec": sec})
        wins[kind] = rows
    assert C.SEED == seed
    obj = {"task": "2I", "arm": arm, "seed": seed, "regime": regime, "cell": cell, "net": j, "site_sd": spec["sd"],
           "world": world.name(), "tags": {"dense": spec["tags"][0], "std": spec["tags"][1]}, "kmax": KMAX,
           "n_win": N_WIN, "methods": list(spec["methods"]), "corr": np.asarray(net.corr).tolist(),
           "t_setup": t_setup, "t_total": time.perf_counter() - t0, "kernel": kernel_info(),
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


def jobs_for(seeds, arms):
    """orden: brazo, red, semilla, celda, regimen (las dos semillas avanzan a la par)."""
    jobs = []
    for arm in arms:
        for j in spec_of(arm)["nets"]:
            for seed in seeds:
                for cell in CELLS:
                    for regime in REGIMES:
                        jobs.append((arm, seed, regime, cell, j))
    return jobs


def frozen_ok():
    if not os.path.exists(HASHES):
        return False
    r = subprocess.run(["sha256sum", "-c", "--quiet", os.path.basename(HASHES)], cwd=HERE, capture_output=True)
    return r.returncode == 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("control", "run", "status"))
    ap.add_argument("--seeds", default="5151,6262")
    ap.add_argument("--arms", default="prim")
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    os.makedirs(PARTS, exist_ok=True)
    if a.mode == "status":
        for arm in ARMS:
            for seed in SEEDS:
                jobs = jobs_for([seed], [arm])
                print(arm, seed, sum(os.path.exists(part_path(*x)) for x in jobs), "/", len(jobs))
        return
    if os.environ.get("OPENBLAS_CORETYPE") != REQUIRED_KERNEL:
        raise SystemExit(f"se corre con OPENBLAS_CORETYPE={REQUIRED_KERNEL}")
    check_kernel()
    if a.mode == "control":
        run_jobs(jobs_for([999], list(CONTROL)), a.workers)
        return
    seeds = [int(x) for x in a.seeds.split(",")]
    arms = a.arms.split(",")
    if any(s not in SEEDS for s in seeds) or any(x not in ARMS for x in arms):
        raise SystemExit("run solo admite las semillas 5151 y 6262 y los brazos prim, thr, s25, s80")
    if not frozen_ok():
        raise SystemExit("codigo no congelado: falta task2i_code_hashes.txt o sha256sum -c falla")
    run_jobs(jobs_for(seeds, arms), a.workers)


if __name__ == "__main__":
    main()
