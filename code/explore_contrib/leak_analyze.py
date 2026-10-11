import sys
import numpy as np
import analyze as A
import xc

rows = A.load(sys.argv[1])
base = sys.argv[2] if len(sys.argv) > 2 else "BOMP-PK"
ms = [m for m in rows[0]["rec"] if not m.startswith("CEIL")]
pairs = [(m, base) for m in ms if m != base]
txt, *_ = A.table(rows, ms, pairs, title="leakage suppression")
print(txt)
for rho in (0.4, 0.5, 0.6, 0.7):
    for rad in (15.0, 25.0):
        nm = f"RATIO{rho}r{int(rad)}"
        for r in rows:
            r["rec"][nm] = xc.ratio_filter([tuple(x) for x in r["rec"][base]], r["loc"][base], rho, rad)
        txt, *_ = A.table(rows, [nm, base], [(nm, base)], title=nm, by_cell=False)
        print(txt)
print({m: round(float(np.mean([r["sec"][m] for r in rows])), 3) for m in rows[0]["sec"]})
