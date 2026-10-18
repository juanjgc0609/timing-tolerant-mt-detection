"""
c19_run.py - exploration C19 (seed 999, exploratory): recall of the Task 2I methods when the sources are not double
couples (c19_lib.FAMILIES).

Design: seed 999 networks task2a_confirm.network(999, cell, j), j = 400.. (never used by an earlier exploration),
cells ladder_4, ladder_6, unilateral_c40, regimes thin and strong; per network a DC catalogue of 20 events and its +SC
corrections (c8_lib.WNet, iid site factors sd 0.5, as in 2I); then, for every family, 5 dense windows (c2_run.draw_dense,
TAG 7101, kmax 7) and 5 standard windows (draw_sources, TAG 7102, kmax 5) from the same streams, so the windows of
different families share sources, origin times, wavelets, coda, clutter and noise and differ only in the mechanism.
TAGs 7101 and 7102 appear in no earlier file.  Methods and engines are those of task2i_run (PRIMARY): BOMP-PK+SC,
MT-FW+SC (cap 6), MT-FW-NOCAP+SC, HUBER-N+SC; hit = epicentre within 5 km.  Per event the fraction of the noise-free
arrival energy inside the 6-D subspace at the true source is stored (clean radiation in the gradient medium, and after
the polarisation jitter and site factors).

    python c19_run.py run --nets 400-419 --workers 4
    python c19_run.py status
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse
import json
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
for p in (HERE, WORK):
    if p not in sys.path:
        sys.path.insert(0, p)

import c19_lib as L  # noqa: E402
import c2_run as C2  # noqa: E402
import c8_lib as C  # noqa: E402
import xc  # noqa: E402
import task2a_confirm as T  # noqa: E402
import task2i_run as R2I  # noqa: E402

TAG_DENSE, TAG_STD = 7101, 7102
KMAX = {"dense": 7, "std": 5}
N_WIN = 5
METHODS = R2I.PRIMARY
PARTS = os.path.join(HERE, "c19_parts")
LOG = os.path.join(HERE, "log_c19.txt")


def part_path(regime, cell, j):
    return os.path.join(PARTS, f"c19_{regime}_{cell}_net{j}.json")


def run_job(args):
    regime, cell, j = args
    out = part_path(regime, cell, j)
    if os.path.exists(out):
        return out, "exists", 0.0
    assert C.SEED == 999
    t0 = time.perf_counter()
    net = L.MNet(cell, j, regime, C.World("iid", 0.5))
    xops = xc.Ops(net.STA, net.corr, xc.make_grid())
    dops = C.DOps(xops)
    fams = {}
    for fam in L.FAMILIES:
        wins = {}
        for kind, tag, draw in (("dense", TAG_DENSE, C2.draw_dense), ("std", TAG_STD, xc.draw_sources)):
            st = net.stream(tag, draw, fam)
            rows = []
            for wi in range(N_WIN):
                w = next(st)
                rec, loc, sec = R2I.run_methods(xops, dops, w["MFn"], w["truth"], KMAX[kind], METHODS)
                frac = [[L.subspace_fraction(net.STA, s, m["clean"]), L.subspace_fraction(net.STA, s, m["pert"])]
                        for s, m in zip(w["truth"], w["meta"])]
                rows.append({"w": wi, "k": len(w["truth"]), "truth": [t.tolist() for t in w["truth"]],
                             "taus": w["taus"], "frac": frac, "records": rec, "loc": loc, "sec": sec})
            wins[kind] = rows
        fams[fam] = wins
    obj = {"exploration": "C19", "seed": 999, "regime": regime, "cell": cell, "net": j,
           "tags": {"dense": TAG_DENSE, "std": TAG_STD}, "kmax": KMAX, "n_win": N_WIN, "methods": list(METHODS),
           "families": list(L.FAMILIES), "corr": np.asarray(net.corr).tolist(),
           "t_total": time.perf_counter() - t0, "numpy": np.__version__,
           "host_time_utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), "windows": fams}
    os.makedirs(PARTS, exist_ok=True)
    tmp = out + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    os.replace(tmp, out)
    return out, "done", obj["t_total"]


def log(msg):
    line = time.strftime("%Y-%m-%d %H:%M:%S UTC ", time.gmtime()) + msg
    print(line, flush=True)
    with open(LOG, "a") as fh:
        fh.write(line + "\n")


def jobs_for(nets):
    """order: network, cell, regime (cells and regimes advance together)."""
    return [(regime, cell, j) for j in nets for cell in T.CELLS for regime in T.REGIMES]


def parse_nets(s):
    a, b = s.split("-")
    return list(range(int(a), int(b) + 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("run", "status"))
    ap.add_argument("--nets", default="400-419")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    os.makedirs(PARTS, exist_ok=True)
    jobs = jobs_for(parse_nets(a.nets))
    todo = [x for x in jobs if not os.path.exists(part_path(*x))]
    if a.mode == "status":
        print(len(jobs) - len(todo), "/", len(jobs))
        return
    log(f"{len(todo)} of {len(jobs)} jobs pending")
    t0, n = time.time(), 0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(run_job, x): x for x in todo}
        for f in as_completed(futs):
            try:
                path, status, tt = f.result()
            except Exception:
                log(f"ERROR {futs[f]}: {traceback.format_exc()}")
                continue
            n += 1
            el = time.time() - t0
            log(f"{status} {os.path.basename(path)} {tt:.1f} s; {n}/{len(todo)}; "
                f"remaining {el / n * (len(todo) - n) / 3600:.2f} h")


if __name__ == "__main__":
    main()
