"""
thesis_evidence_table.py - builds the evidence table of the thesis from the result files.  Each row of
docs/thesis/evidence_table.tsv names a source file and a regular expression; the script finds the first matching line,
checks that every listed value occurs in it as a whole token (not inside a longer number or word, and "MET" not as
part of "NOT MET"), and prints a Markdown table with the file and line of each value.
A value that is not found stops the script with an error.
    python code/thesis_evidence_table.py            (from the repository root)
"""
import csv
import re
import sys

def token_in(v, line):
    for m in re.finditer(re.escape(v), line):
        a, b = m.start(), m.end()
        before = line[a - 1] if a > 0 else " "
        after = line[b] if b < len(line) else " "
        if before.isalnum() or before == "." or (v[0].isdigit() and before in "+-") or after.isalnum() or \
                (after == "." and b + 1 < len(line) and line[b + 1].isdigit()):
            continue
        if v == "MET" and line[max(0, a - 4):a] == "NOT ":
            continue
        return True
    return False


rows = list(csv.DictReader(open("docs/thesis/evidence_table.tsv"), delimiter="\t"))
out = ["| # | Data | Test | Claim | Value | Source (file:line) |", "|:--|:--|:--|:--|:--|:--|"]
bad = 0
for r in rows:
    lines = open(r["source"]).read().split("\n")
    hit = [(i, s) for i, s in enumerate(lines, 1) if re.search(r["pattern"], s)]
    if not hit:
        print("NO LINE", r["id"], r["pattern"], file=sys.stderr)
        bad += 1
        continue
    n, s = hit[0]
    vals = r["values"].split(";")
    miss = [v for v in vals if not token_in(v, s)]
    if miss:
        print("VALUE NOT IN LINE", r["id"], miss, "|", s.strip(), file=sys.stderr)
        bad += 1
    out.append(f"| {r['id']} | {r['data']} | {r['test']} | {r['claim']} | {'; '.join(vals)} | `{r['source']}:{n}` |")
print("\n".join(out))
sys.exit(1 if bad else 0)
