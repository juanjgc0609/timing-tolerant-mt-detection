"""
task2f_run.py - Tarea 2F (PREREGISTRO_1j_1k.txt, seccion "Tarea 2F"): el tope por estacion del estadistico de subespacio
del tensor de momento en secuencias densas de eventos solapados.  Corredor de las semillas 4242 y 8888 y de los controles
en la semilla 999.

Diseno implementado (todo como el piloto explore_contrib/c2_ablate.py, con la semilla como parametro):
  redes         task2a_confirm.network(seed, cell, j); retardos N(0, 0.3 s) por estacion y fase, fijos por red
  catalogo      task2a_confirm.catalogue_real, 20 eventos, semillas [seed, 2003, ci, j, ri] (uno por red y regimen, 2A)
  SC            station_corrections.estimate_station_corrections con ese catalogo
  ventanas      gen_realistic.generate_real con los retardos de la red, flujos [seed, TAG, ci, j, ri, 41|233|61]
                (xc.window_stream con la semilla como parametro), filtro adaptado task1o_consistent.matched
                  brazo denso        TAG 4201, c2_run.draw_dense (3-5 eventos), j = 0..95, 10 ventanas, kmax 7
                  brazo estandar     TAG 4202, task1o_consistent.draw_sources (0-2 eventos), j = 96..111, 10 ventanas,
                                     kmax 5
                TAG 4201 y 4202 no aparecen en ningun archivo anterior (grep de todos los TAG y entropias de rng).
  operadores    c6_lib.ModelOps(STA, xc.make_grid(), c6_lib.Homog(), corr) (rejilla 3 km x 3 profundidades, origenes
                DEFAULT_TAUS, modelo homogeneo registrado), c2_ablate.PKDictCap(ops.BD)
  metodos       MT-FW+SC        c6_lib.mtfw_dense(ops, MFn, truth, kmax)
                MT-FW-NOCAP+SC  c6_lib.mtfw_dense(ops, MFn, truth, kmax, cap=False)
                BOMP-PK+SC      c6_lib.bomp_pk_m(ops, MFn, truth, kmax) (= task2d_methods.bomp_pk_record sobre ops.PD)
                BOMP-PK-CAP+SC  task2d_methods.bomp_pk_record(c2_ablate.PKDictCap(ops.BD), MFn, truth, kmax)
                (control, solo semilla 999) BOMP-PK-SEQ+SC  c2_ablate.bomp_seq(ops.PD, MFn, truth, kmax)
  nucleo        OPENBLAS_CORETYPE=Prescott en todo proceso; cada trabajo comprueba con threadpoolctl que la OpenBLAS de
                numpy cargo Prescott y lo guarda en su parte.

Partes (una por semilla, brazo, regimen, celda y red; reanudable): task2f_parts/{arm}_seed{s}_{reg}_{cell}_net{j}.json

    OPENBLAS_CORETYPE=Prescott python task2f_run.py control --workers 2      # semilla 999: reproduccion del piloto
    OPENBLAS_CORETYPE=Prescott python task2f_run.py run --seeds 4242,8888 --workers 2
    python task2f_run.py status
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MPLBACKEND", "Agg")
import argparse
import json
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

import xc  # noqa: E402  (explore_contrib harness; only make_grid and its imports are used)
import c2_run as C2  # noqa: E402
import c2_ablate as CA  # noqa: E402
import c6_lib as L  # noqa: E402
import task2a_confirm as T  # noqa: E402
import task2d_methods as D2  # noqa: E402
from gen_realistic import generate_real  # noqa: E402
from station_corrections import estimate_station_corrections  # noqa: E402
from task1o_consistent import draw_sources, matched  # noqa: E402
from task1y_clutter import CLUT_AMP  # noqa: E402

REQUIRED_KERNEL = "Prescott"
CELLS = T.CELLS
REGIMES = T.REGIMES
SEEDS = (4242, 8888)
N_WIN = 10
ARMS = {
    "dense": {"tag": 4201, "draw": C2.draw_dense, "kmax": 7, "nets": tuple(range(0, 96))},
    "std": {"tag": 4202, "draw": draw_sources, "kmax": 5, "nets": tuple(range(96, 112))},
}
METHODS = ("MT-FW+SC", "BOMP-PK+SC", "BOMP-PK-CAP+SC", "MT-FW-NOCAP+SC")
PARTS = os.path.join(HERE, "task2f_parts")
LOG = os.path.join(HERE, "log2f.txt")
FAIL_FILE = os.path.join(HERE, "task2f_FAIL.txt")

# controles en la semilla 999 (datos del piloto): brazo denso de c2_ablate (redes 10-15, TAG 3002, kmax 7) y brazo
# estandar de C6 (redes 20-21, TAG 3006, draw_sources, kmax 5; metodos HOM)
CONTROL = {
    "cdense": {"tag": C2.TAG, "draw": C2.draw_dense, "kmax": C2.KMAX, "nets": tuple(range(10, 16)), "seq": True},
    "cstd": {"tag": 3006, "draw": draw_sources, "kmax": 5, "nets": (20, 21), "seq": False},
}


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


# ------------------------------------------------------------------------------------------ red y ventanas
def setup_net(seed, cell, j, regime):
    """red, retardos, catalogo de 20 eventos y correcciones de estacion (task2a_confirm; = xc.NetCfg con seed)."""
    ci, _kind, _nst, nclut = T.net_setup(seed, cell)
    ri = REGIMES.index(regime)
    STA, statics = T.network(seed, cell, j)
    cat, drawn, _truth = T.catalogue_real(STA, T.K_CAT, statics, [seed, T.TAG_CAT, ci, j, ri], nclut,
                                          CLUT_AMP[regime])
    corr, _n = estimate_station_corrections(STA, cat)
    return {"ci": ci, "ri": ri, "nclut": nclut, "STA": STA, "statics": statics, "corr": corr, "drawn": drawn}


def window_stream(seed, net, tag, draw):
    """xc.window_stream con la semilla como parametro."""
    b = [seed, tag, net["ci"], net["j"], net["ri"]]
    rng, rng_m, rng_h = (np.random.default_rng(b + [41]), np.random.default_rng(b + [233]),
                         np.random.default_rng(b + [61]))
    while True:
        srcs, taus = draw(rng)
        D, _picks, _tr, _amp = generate_real(net["STA"], srcs, taus, net["nclut"], rng, rng_m, rng_h,
                                             CLUT_AMP[net["regime"]], statics=net["statics"])
        MF, sd = matched(D)
        yield {"MFn": MF / sd, "truth": [np.asarray(s, float) for s in srcs], "taus": [float(t) for t in taus]}


def run_methods(ops, PDc, MFn, tr, kmax, seq=False):
    """las definiciones de c2_ablate.run (mismas llamadas)."""
    rec, loc, sec = {}, {}, {}
    fns = [("MT-FW+SC", lambda: L.mtfw_dense(ops, MFn, tr, kmax)),
           ("MT-FW-NOCAP+SC", lambda: L.mtfw_dense(ops, MFn, tr, kmax, cap=False)),
           ("BOMP-PK+SC", lambda: L.bomp_pk_m(ops, MFn, tr, kmax)),
           ("BOMP-PK-CAP+SC", lambda: ([(float(s), int(h)) for s, h in D2.bomp_pk_record(PDc, MFn, tr, kmax=kmax)],
                                       None))]
    if seq:
        fns.append(("BOMP-PK-SEQ+SC", lambda: (CA.bomp_seq(ops.PD, MFn, tr, kmax), None)))
    for name, fn in fns:
        t = time.perf_counter()
        r, l_ = fn()
        sec[name] = time.perf_counter() - t
        rec[name] = [[float(s), int(h)] for s, h in r]
        if l_ is not None:
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
    spec = ARMS[arm] if arm in ARMS else CONTROL[arm]
    t0 = time.perf_counter()
    net = setup_net(seed, cell, j, regime)
    net.update({"j": j, "regime": regime})
    ops = L.ModelOps(net["STA"], xc.make_grid(), L.Homog(), net["corr"])
    PDc = CA.PKDictCap(ops.BD)
    t_setup = time.perf_counter() - t0
    st = window_stream(seed, net, spec["tag"], spec["draw"])
    wins = []
    for wi in range(N_WIN):
        w = next(st)
        rec, loc, sec = run_methods(ops, PDc, w["MFn"], w["truth"], spec["kmax"], seq=spec.get("seq", False))
        wins.append({"w": wi, "k": len(w["truth"]), "truth": [t.tolist() for t in w["truth"]], "taus": w["taus"],
                     "records": rec, "loc": loc, "sec": sec})
    obj = {"task": "2F", "arm": arm, "seed": seed, "regime": regime, "cell": cell, "net": j, "tag": spec["tag"],
           "kmax": spec["kmax"], "n_win": N_WIN, "cat_drawn": net["drawn"], "corr": np.asarray(net["corr"]).tolist(),
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
    jobs = []
    for seed in seeds:
        for arm in arms:
            spec = ARMS[arm] if arm in ARMS else CONTROL[arm]
            for j in spec["nets"]:
                for cell in CELLS:
                    for regime in REGIMES:
                        jobs.append((arm, seed, regime, cell, j))
    return jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("control", "run", "status"))
    ap.add_argument("--seeds", default="4242,8888")
    ap.add_argument("--arms", default="dense,std")
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    os.makedirs(PARTS, exist_ok=True)
    if a.mode == "status":
        for seed in SEEDS:
            for arm in ARMS:
                jobs = jobs_for([seed], [arm])
                print(seed, arm, sum(os.path.exists(part_path(*x)) for x in jobs), "/", len(jobs))
        return
    if os.environ.get("OPENBLAS_CORETYPE") != REQUIRED_KERNEL:
        raise SystemExit(f"se corre con OPENBLAS_CORETYPE={REQUIRED_KERNEL}")
    check_kernel()
    if a.mode == "control":
        run_jobs(jobs_for([999], ["cdense", "cstd"]), a.workers)
        return
    seeds = [int(x) for x in a.seeds.split(",")]
    if any(s not in SEEDS for s in seeds):
        raise SystemExit("run solo admite las semillas 4242 y 8888")
    run_jobs(jobs_for(seeds, a.arms.split(",")), a.workers)


if __name__ == "__main__":
    main()
