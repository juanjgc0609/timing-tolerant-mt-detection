"""rescore_2a2b.py - Tasks 2A and 2B (analyze_2a.py, analyze_2b.py, unchanged) with repeat hits counted as false.
    python rescore/rescore_2a2b.py RULE WORKDIR
RULE  legacy  stored records as they are (must reproduce results/task2a/analyze_2a_out.txt and task2b's)
      dedup   froc_strict.drop_duplicates on every stored record: a repeated event index becomes false.  The 2A/2B
              parts store (score, label) only, without locations, so the one-to-one reassignment to another unmatched
              event and the origin-time rule LOC+T cannot be applied; the 2A-e operating points are re-frozen on the
              seed-999 calibration under the same rule.  GaMMA and PyOcto store per-window counts per configuration
              and are left as registered.
Copies the stored json files to WORKDIR, rewrites them under RULE and runs the two analyses there (they read their
inputs from the working directory); outputs WORKDIR/analyze_2a_out.txt and WORKDIR/analyze_2b_out.txt."""
import glob
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(HERE, "..", "results")
sys.path.insert(0, HERE)

import froc_strict as FS  # noqa: E402

rule, work = sys.argv[1], os.path.abspath(sys.argv[2])
if rule not in ("legacy", "dedup"):
    raise SystemExit(__doc__)
os.makedirs(work, exist_ok=True)
for f in glob.glob(os.path.join(RES, "task2a", "task2a_*.json")) + glob.glob(os.path.join(RES, "task2b", "task2b_seed*.json")):
    dst = os.path.join(work, os.path.basename(f))
    if rule == "legacy":
        shutil.copy(f, dst)
        continue
    cells = json.load(open(f))
    if isinstance(cells, list):
        for c in cells:
            for m in c.get("methods", {}).values():
                if "records" in m:
                    m["records"] = [[[s, h] for (s, _h), h in zip(r, FS.drop_duplicates([h for _s, h in r]))]
                                    for r in m["records"]]
    json.dump(cells, open(dst, "w"))
env = dict(os.environ, PYTHONPATH=HERE)
run = lambda *a: subprocess.run([sys.executable, *a], cwd=work, env=env, check=True, stdout=subprocess.DEVNULL)  # noqa
if rule == "dedup":
    run(os.path.join(HERE, "analyze_2a.py"), "--freeze")
run(os.path.join(HERE, "analyze_2a.py"), "--out", "analyze_2a_out.txt")
run(os.path.join(HERE, "analyze_2b.py"), "--out", "analyze_2b_out.txt")
