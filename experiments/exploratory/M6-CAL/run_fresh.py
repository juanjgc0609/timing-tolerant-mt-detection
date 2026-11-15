"""A3: M6-CAL controls on NEW seeds with the step-3 settings frozen (near_km 10, thr 23.8). No calibration here."""
import os, sys, importlib.util, multiprocessing as mp, time
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["OPENBLAS_NUM_THREADS"] = "1"
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.argv = ["run.py", "200", "200", "20", "4", "op", "10"]
spec = importlib.util.spec_from_file_location("m6cal", HERE / "run.py")
R = importlib.util.module_from_spec(spec); sys.modules["m6cal"] = R; spec.loader.exec_module(R)
def job(a):
    return R.job_split(a)
THR, NEW = 23.8, 20261114
if __name__ == "__main__":
    t0 = time.time(); N = int(os.environ.get("N", 300))
    with mp.Pool(4) as pool:
        for kind in ("C1", "C3"):
            out = pool.map(job, [((NEW, 2, len(kind), 1, c), kind, True, THR) for c in range(N)], chunksize=1)
            k = sum(o[0] > 0 for o in out); lo, hi = R.cp(k, N)
            print(kind, "mismatched", N, "scenes split", k, f"[{lo:.3f}, {hi:.3f}]", "greedy0 splits", sum(o[1] for o in out), flush=True)
    print(f"{time.time()-t0:.0f} s")
