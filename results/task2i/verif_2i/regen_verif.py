"""regen_verif.py - regenerate one 2I part from scratch in a fresh process (writes ONLY into verif_2i/regen_parts_verif)
and compare exactly with the stored part.  usage: regen_verif.py arm seed regime cell j"""
import os, sys, json
os.environ["OPENBLAS_CORETYPE"] = "Prescott"; os.environ["OMP_NUM_THREADS"] = "1"; os.environ["OPENBLAS_NUM_THREADS"] = "1"
HERE = os.path.dirname(os.path.abspath(__file__)); WORK = os.path.dirname(HERE)
sys.path.insert(0, WORK); os.chdir(WORK)
import task2i_run as R
ORIG = R.PARTS
R.PARTS = os.path.join(HERE, "regen_parts_verif"); R.FAIL_FILE = os.path.join(HERE, "regen_FAIL_verif.txt")
arm, seed, reg, cell, j = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4], int(sys.argv[5])
out, status, tt = R.run_job((arm, seed, reg, cell, j))
a = json.load(open(out)); b = json.load(open(os.path.join(ORIG, os.path.basename(out))))
bad = []
for k in a:
    if k in ("t_setup", "t_total", "host_time_utc"): continue
    if k == "windows":
        for kind in a[k]:
            for wa, wb in zip(a[k][kind], b[k][kind]):
                for f in wa:
                    if f == "sec": continue
                    if wa[f] != wb[f]: bad.append((kind, wa["w"], f))
            if len(a[k][kind]) != len(b[k][kind]): bad.append((kind, "len"))
    elif a[k] != b[k]: bad.append(k)
print(os.path.basename(out), status, f"{tt:.1f}s", "stored at", b["host_time_utc"], "IDENTICAL (records, loc, truth, taus, site, corr, meta)" if not bad else f"DIFFERS {bad}", flush=True)
