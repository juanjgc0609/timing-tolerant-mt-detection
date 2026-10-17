"""assemble.py - builds analysis/ for r0_analyze.py: det_eval_<s>.json = atoms of M2 (this verification, m2/) plus
M1 (run_close) plus M5/M6 when re-run here (m56/); detM3_eval_<s>.json copied from run_close; picks_eval_<s>.json
from picks/.  Only the 'atoms' dict is used by r0_analyze."""
import json, os, shutil
B = "/home/claude/work/r0_verif"
C = "/home/claude/work/r0/run_close"
os.makedirs(f"{B}/analysis", exist_ok=True)
for s in ("N8", "N6", "N4", "U8"):
    out = json.load(open(f"{C}/det_eval_{s}.json"))
    for d, ms in ((f"{B}/m2/det_eval_{s}.json", ["M2"]), (f"{B}/m56/det_eval_{s}.json", ["M5", "M6"])):
        if os.path.exists(d):
            r = json.load(open(d))
            for m in ms:
                if m in r["atoms"]:
                    out["atoms"][m] = r["atoms"][m]
    json.dump(out, open(f"{B}/analysis/det_eval_{s}.json", "w"))
    shutil.copy(f"{C}/detM3_eval_{s}.json", f"{B}/analysis/detM3_eval_{s}.json")
    p = f"{B}/picks/picks_eval_{s}.json"
    if os.path.exists(p):
        shutil.copy(p, f"{B}/analysis/picks_eval_{s}.json")
    print(s, sorted(out["atoms"]), {m: len(v) for m, v in out["atoms"].items()})
