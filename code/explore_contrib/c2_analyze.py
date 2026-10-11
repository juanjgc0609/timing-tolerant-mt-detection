import sys
import numpy as np
import analyze as A
rows = A.load(sys.argv[1])
ms = ["BOMP-PK", "TS-SFW", "SLIDE", "MT-FW", "MASK6g3", "CEIL-A", "CEIL-B"]
pairs = [("TS-SFW", "BOMP-PK"), ("SLIDE", "BOMP-PK"), ("MT-FW", "BOMP-PK"), ("MASK6g3", "BOMP-PK"),
         ("CEIL-A", "BOMP-PK"), ("CEIL-B", "BOMP-PK")]
txt, *_ = A.table(rows, ms, pairs, title="C2 dense overlaps (3-5 events, kmax 7)")
print(txt)
print({m: round(float(np.mean([r["sec"][m] for r in rows])), 3) for m in rows[0]["sec"]})
# recall without false-alarm constraint (all atoms) by number of events
for k in (3, 4, 5):
    sub = [r for r in rows if r["k"] == k]
    s = f"  k={k} ({len(sub)} windows) any-threshold recall: "
    for m in ("BOMP-PK", "TS-SFW", "SLIDE", "MT-FW", "MASK6g3"):
        h = sum(len({hh for _s, hh in r["rec"][m] if hh >= 0}) for r in sub)
        s += f"{m} {h / sum(r['k'] for r in sub):.3f}  "
    print(s)
