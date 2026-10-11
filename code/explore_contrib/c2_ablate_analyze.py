"""C2 follow-up analysis: ablations on the 180 dense-overlap windows, plus consistency with the stored c2 records."""
import json
import sys
import numpy as np
import analyze as A

rows = json.load(open(sys.argv[1]))
old = {(r["cell"], r["net"], r["w"]): r for r in json.load(open(sys.argv[2]))}
same_mt = same_bo = 0
for r in rows:
    o = old[(r["cell"], r["net"], r["w"])]
    same_mt += np.allclose(np.array(r["rec"]["MT-FW"]), np.array(o["rec"]["MT-FW"]), atol=1e-9)
    same_bo += np.allclose(np.array(r["rec"]["BOMP-PK"]), np.array(o["rec"]["BOMP-PK"]), atol=1e-9)
print(f"regenerated windows reproduce the stored c2 records: MT-FW {same_mt}/{len(rows)}, BOMP-PK {same_bo}/{len(rows)}")
ms = ["BOMP-PK", "MT-FW", "MT-FW-NOCAP", "BOMP-PK-CAP", "BOMP-PK-SEQ"]
pairs = [("MT-FW", "BOMP-PK"), ("MT-FW-NOCAP", "BOMP-PK"), ("MT-FW", "MT-FW-NOCAP"), ("BOMP-PK-CAP", "BOMP-PK"),
         ("MT-FW", "BOMP-PK-CAP"), ("BOMP-PK-SEQ", "BOMP-PK"), ("MT-FW", "BOMP-PK-SEQ")]
txt, *_ = A.table(rows, ms, pairs, title="C2 ablation (dense overlaps, kmax 7)")
print(txt)
print("s/window:", {m: round(float(np.mean([r["sec"][m] for r in rows])), 3) for m in ms})
