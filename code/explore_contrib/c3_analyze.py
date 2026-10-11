import sys
import numpy as np
import analyze as A
allrows = A.load(sys.argv[1])
ms = ["BOMP-PK", "BOMP-PK-8", "BOMP-PK-12", "BOMP-SS", "BOMP-OL", "CEIL-A", "CEIL-B"]
pairs = [("BOMP-PK-8", "BOMP-PK"), ("BOMP-PK-12", "BOMP-PK"), ("BOMP-SS", "BOMP-PK"), ("BOMP-OL", "BOMP-PK"),
         ("BOMP-OL", "BOMP-PK-8"), ("CEIL-A", "BOMP-PK"), ("CEIL-B", "BOMP-PK")]
for c in ["S03sc", "S05nc", "S08nc", "S08sc", "V75sc", "V75nc"]:
    rows = [r for r in allrows if r["cond"] == c]
    if not rows:
        continue
    txt, *_ = A.table(rows, ms, pairs, title=f"C3 condition {c}")
    print(txt)
    print("   median |corr - statics| (s): catalogue %.3f, online %.3f; confident detections learned per network (end) %.1f" % (
        np.median([r["corr_err_cat"] for r in rows]), np.median([r["corr_err_ol"] for r in rows]),
        np.mean([r["n_learn"] for r in rows if r["w"] == max(q["w"] for q in rows)])))
    print("   s/window:", {m: round(float(np.mean([r["sec"][m] for r in rows])), 3) for m in rows[0]["sec"]})
