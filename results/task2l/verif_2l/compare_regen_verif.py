import json, os, glob
SKIP = {"t_setup", "t_s2", "t_total", "host_time_utc", "sec"}
def strip(o):
    if isinstance(o, dict):
        return {k: strip(v) for k, v in o.items() if k not in SKIP}
    if isinstance(o, list):
        return [strip(v) for v in o]
    return o
n_ok = 0
for f in sorted(glob.glob("regen_parts_verif/*.json")):
    a = strip(json.load(open(f)))
    b = strip(json.load(open(os.path.join("..", "task2l_parts", os.path.basename(f)))))
    same = json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)
    n_ok += same
    na = sum(len(x) for k in a["windows"] for w in a["windows"][k] for x in w["atoms"].values())
    print(os.path.basename(f), "IDENTICAL" if same else "DIFFERENT", "atoms compared:", na, "orig written", json.load(open(os.path.join("..", "task2l_parts", os.path.basename(f))))["host_time_utc"])
print("identical:", n_ok)
