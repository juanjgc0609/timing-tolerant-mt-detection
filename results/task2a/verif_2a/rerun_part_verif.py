import os, sys, json
os.environ.setdefault("TASK2A_SCRATCH", "/tmp/claude-0/-home-claude/d3c612af-4e2f-52d9-b554-c7fa18d867f4/scratchpad/verif")
sys.path.insert(0, "/home/claude/work")
import task2a_confirm as T
seed, regime, cell, j = int(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4])
out = f"/home/claude/work/verif_2a/seed{seed}_{regime}_{cell}_net{j}_verif.json"
r = T.run_job((seed, regime, cell, j, 20, False, out))
orig = json.load(open(f"/home/claude/work/task2a_parts/seed{seed}_{regime}_{cell}_net{j}.json"))
same = {}
for key in ("STA", "statics", "corr", "cat_src", "cat_tau"):
    same[key] = (r[key] if not hasattr(r[key], "tolist") else r[key]) == orig[key]
for w, (a, b) in enumerate(zip(r["windows"], orig["windows"])):
    for k in ("truth", "taus", "records", "per", "picker", "qm_events"):
        if a[k] != b[k]:
            same.setdefault("diff", []).append((w, k))
print(same)
