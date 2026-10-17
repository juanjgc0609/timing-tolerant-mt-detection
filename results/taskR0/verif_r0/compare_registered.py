"""compare_registered.py - every point recall and contrast (point, CI, per subnet) of analysis/analysis_eval.json
against the registered results/taskR0/analysis_eval.json; prints mismatches > 1e-9 and a count of matches.
Also diffs the descriptive lines of the two .txt files."""
import json, sys
a = json.load(open("/home/claude/work/repo/results/taskR0/analysis_eval.json"))
b = json.load(open(sys.argv[1]))
nm = nx = 0
miss = []
for k, v in a["point"].items():
    if k not in b["point"]:
        miss.append(k); continue
    if abs(v - b["point"][k]) > 1e-9:
        print("POINT DIFF", k, v, b["point"][k]); nx += 1
    else:
        nm += 1
print("points matched", nm, "differ", nx, "not re-run", len(miss), sorted({k.split('|')[1] for k in miss}))
nm = nx = 0; miss = []
for k, v in a["contrasts"].items():
    if k not in b["contrasts"]:
        miss.append(k); continue
    w = b["contrasts"][k]
    d = [abs(v["point"] - w["point"])] + [abs(x - y) for x, y in zip(v["ci"], w["ci"])] + [abs(x - y) for x, y in zip(v["per"], w["per"])]
    if max(d) > 1e-9:
        print("CONTRAST DIFF", k, v, w); nx += 1
    else:
        nm += 1
print("contrasts matched", nm, "differ", nx, "not computed", len(miss), sorted({k.split('|')[0] for k in miss}))
