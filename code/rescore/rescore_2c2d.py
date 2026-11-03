"""rescore_2c2d.py - Tasks 2C and 2D (analyze_2c.py, analyze_2d.py, unchanged) with repeat hits counted as false.
    python rescore/rescore_2c2d.py RULE PARTS_DIR WORKDIR
RULE       legacy  the regenerated records as they are (must reproduce results/task2c/analyze_2c_out.txt and task2d's)
           dedup   froc_strict.drop_duplicates on every record: a repeated event index becomes false.  The 2A, 2C and 2D
                   parts store (score, label) per atom without locations or origin times, so, as for 2A and 2B, the
                   one-to-one reassignment and the origin-time rule LOC+T cannot be applied.
PARTS_DIR  folder with task2a_parts/, task2c_parts/ and task2d_parts/ regenerated with the frozen code
           (posthoc_e52/regen_2a_min.py, then task2c_run.py and task2d_run.py as their queues; error register E52).
Writes WORKDIR/{task2a,task2c,task2d}_parts under RULE, copies the two analyses next to them (they read their parts
from their own folder) and writes WORKDIR/analyze_2c_out.txt and WORKDIR/analyze_2d_out.txt."""
import glob
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)

import froc_strict as FS  # noqa: E402


def dedup_record(r):
    hits = FS.drop_duplicates([int(x[1]) for x in r])
    return [[x[0], h] + list(x[2:]) for x, h in zip(r, hits)]


def main(rule, src, work):
    if rule not in ("legacy", "dedup"):
        raise SystemExit(__doc__)
    os.makedirs(work, exist_ok=True)
    for sub in ("task2a_parts", "task2c_parts", "task2d_parts"):
        os.makedirs(os.path.join(work, sub), exist_ok=True)
        for f in glob.glob(os.path.join(src, sub, "*.json")):
            dst = os.path.join(work, sub, os.path.basename(f))
            if rule == "legacy":
                shutil.copy(f, dst)
                continue
            p = json.load(open(f))
            for w in p.get("windows", []):
                for m, r in w.get("records", {}).items():
                    w["records"][m] = dedup_record(r)
                if "gl_rec_nv" in w:
                    w["gl_rec_nv"] = dedup_record(w["gl_rec_nv"])
            json.dump(p, open(dst, "w"))
    for a in ("analyze_2c.py", "analyze_2d.py"):
        shutil.copy(os.path.join(HERE, a), os.path.join(work, a))
    env = dict(os.environ, PYTHONPATH=HERE)
    for a, o in (("analyze_2c.py", "analyze_2c_out.txt"), ("analyze_2d.py", "analyze_2d_out.txt")):
        subprocess.run([sys.executable, os.path.join(work, a), "--out", os.path.join(work, o)], cwd=work, env=env,
                       check=True, stdout=subprocess.DEVNULL)


if __name__ == "__main__":
    main(*sys.argv[1:4])
