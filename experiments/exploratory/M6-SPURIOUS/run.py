"""M6-SPURIOUS (calendar M-2): cause of the positive spurious pair gains of M6-CAL step 5. See card.md.
    python run.py [N] [PROCS] > run_out.txt
"""
import importlib.util
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
_argv = sys.argv
sys.argv = ["run.py", "200", "200", "20", "4", "op", "10"]          # near_km 10, the settings of M6-CAL step 5
spec = importlib.util.spec_from_file_location("m6cal", HERE.parent / "M6-CAL" / "run.py")
R = importlib.util.module_from_spec(spec)
sys.modules["m6cal"] = R
spec.loader.exec_module(R)
sys.argv = _argv

from ttmt.stats import exact_mcnemar  # noqa: E402

SEED = 20261127
N = int(sys.argv[1]) if len(sys.argv) > 1 else 600
PROCS = int(sys.argv[2]) if len(sys.argv) > 2 else 4
KINDS = {"C1": 1, "C3": 3}          # one source; target + neighbour + isolated (M6-BIAS.null_scene)


def job(args):
    kind, i, mismatch = args
    R.J.SNR = 12
    key = [SEED, KINDS[kind], i]
    srcs, Y, q = R._scene(key, kind, mismatch)
    a = R.m6.adcg_loop(R.E, Y, q, R.KMAX, R.LAG, thr=-np.inf, free="joint", near_km=R.NEAR_KM)
    sp = [x for x in R.MB.spurious(a, srcs) if x["pair_gain"] is not None]
    if not sp:
        return {"kind": kind, "i": i, "mis": mismatch, "gain": None}
    b = max(sp, key=lambda x: x["pair_gain"])
    d = [float(np.linalg.norm(b["x"][:2] - s[0][:2])) for s in srcs]
    dt = [float(b["tau"] - s[1]) for s in srcs]
    near = int(np.argmin([max(dd / 5.0, abs(tt) / 1.5) for dd, tt in zip(d, dt)]))
    return {"kind": kind, "i": i, "mis": mismatch, "gain": float(b["pair_gain"]), "d": d, "dt": dt, "nearest": near}


def main():
    t0 = time.time()
    jobs = [(k, i, m) for k in KINDS for i in range(N) for m in (True, False)]
    with mp.Pool(PROCS) as p:
        rows = p.map(job, jobs, chunksize=4)
    json.dump(rows, open(HERE / "rows.json", "w"))
    pos = {(r["kind"], r["i"], r["mis"]): (r["gain"] is not None and r["gain"] > 0) for r in rows}
    P = {}
    print(f"M6-SPURIOUS seed {SEED}, {N} scenes per kind, SNR_blk 12, kmax 4, near_km {R.NEAR_KM}, thr -inf")
    for k in KINDS:
        for m in (True, False):
            g = np.array([r["gain"] for r in rows if r["kind"] == k and r["mis"] == m and r["gain"] is not None])
            P[k, m] = np.mean([pos[k, i, m] for i in range(N)])
            q = np.quantile(g, [.5, .9, .95]) if g.size else [np.nan] * 3
            print(f"{k} {'mismatched' if m else 'exact     '}: {g.size} scenes with a spurious near atom; "
                  f"positive gain {int(round(P[k, m] * N))}/{N} = {P[k, m]:.3f}; gain median {q[0]:.1f} q90 {q[1]:.1f} "
                  f"q95 {q[2]:.1f} max {g.max() if g.size else np.nan:.1f}")
    for k in KINDS:
        b = sum(pos[k, i, True] and not pos[k, i, False] for i in range(N))
        c = sum(pos[k, i, False] and not pos[k, i, True] for i in range(N))
        print(f"{k} paired: positive only mismatched {b}, only exact {c}, exact McNemar p {exact_mcnemar(b, c):.2e}")
    posC3 = [r for r in rows if r["kind"] == "C3" and r["mis"] and r["gain"] is not None and r["gain"] > 0]
    share = np.bincount([r["nearest"] for r in posC3], minlength=3) / max(len(posC3), 1)
    print(f"C3 mismatched, positive scenes: nearest source target / neighbour / isolated = "
          f"{share[0]:.2f} / {share[1]:.2f} / {share[2]:.2f} (n {len(posC3)})")
    pm, pe = P["C3", True], P["C3", False]
    b = sum(pos["C3", i, True] and not pos["C3", i, False] for i in range(N))
    c = sum(pos["C3", i, False] and not pos["C3", i, True] for i in range(N))
    h1 = pe <= 0.2 * pm and exact_mcnemar(b, c) < 0.01
    h2 = P["C1", True] <= 0.2 * pm and share[1] >= 0.60
    h3 = pe <= 0.2 * pm and P["C1", False] <= 0.2 * pm and P["C1", True] <= 0.2 * pm
    if pm == 0:
        print("reading (card): no positive spurious gain with neighbours; the tail is not reproduced")
        return
    print(f"reading (card): H1 {h1}, H2 {h2}, H3 {h3} -> "
          + (", ".join(n for n, v in (("H1", h1), ("H2", h2), ("H3", h3)) if v) or "no clear cause"))
    print(f"wall {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
