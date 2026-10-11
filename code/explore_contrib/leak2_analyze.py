import sys
import numpy as np
import analyze as A
rows = A.load(sys.argv[1])
ms = ["BOMP-PK", "BOMP-PK-MASK", "MT-FW", "MT-FW-MASK", "CEIL-A", "CEIL-B"]
pairs = [("BOMP-PK-MASK", "BOMP-PK"), ("MT-FW-MASK", "MT-FW"), ("MT-FW-MASK", "BOMP-PK"), ("MT-FW", "BOMP-PK"),
         ("CEIL-A", "BOMP-PK"), ("CEIL-B", "BOMP-PK")]
txt, *_ = A.table(rows, ms, pairs, title=sys.argv[2] if len(sys.argv) > 2 else "")
print(txt)
print({m: round(float(np.mean([r["sec"][m] for r in rows])), 3) for m in rows[0]["sec"]})
