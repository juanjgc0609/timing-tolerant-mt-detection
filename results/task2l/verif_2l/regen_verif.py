"""regen_verif.py - regenerate selected 2L parts from scratch in fresh processes (one process per part) into
verif_2l/regen_parts_verif and compare with task2l_parts.  Does not write to any 2L file."""
import os, sys, json, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
OUT = os.path.join(HERE, "regen_parts_verif")

def child(seed, regime, cell, j):
    sys.path.insert(0, WORK)
    import task2l_run as R
    R.PARTS = OUT
    R.LOG = os.path.join(HERE, "regen_log_verif.txt")
    R.FAIL_FILE = os.path.join(HERE, "regen_FAIL_verif.txt")
    print(R.run_job(("main", seed, regime, cell, j)), flush=True)

if __name__ == "__main__":
    if sys.argv[1] == "child":
        child(int(sys.argv[2]), sys.argv[3], sys.argv[4], int(sys.argv[5]))
    else:
        jobs = [l.split() for l in open(os.path.join(HERE, "regen_jobs_verif.txt")) if l.strip()]
        from concurrent.futures import ThreadPoolExecutor
        def go(a):
            return subprocess.run([sys.executable, __file__, "child"] + a, capture_output=True, text=True)
        with ThreadPoolExecutor(2) as ex:
            for a, r in zip(jobs, ex.map(go, jobs)):
                print(a, r.returncode, r.stdout.strip()[-200:], r.stderr.strip()[-500:], flush=True)
