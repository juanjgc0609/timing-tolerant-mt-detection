import json, sys
import numpy as np
import analyze as A

path = sys.argv[1]
rows = A.load(path)
G = ["g3d3", "g3d2", "g6d3", "g6d2", "g9d3", "g9d2"]
NODES = {"g3d3": 675, "g3d2": 450, "g6d3": 192, "g6d2": 128, "g9d3": 75, "g9d2": 50}
out = []
txt, *_ = A.table(rows, ["CEIL-A", "CEIL-B", "g3d3:BOMP-PK"], [("CEIL-A", "g3d3:BOMP-PK"), ("CEIL-B", "g3d3:BOMP-PK")],
                  title="C1 headroom (3 km grid)", by_cell=False)
out.append(txt)
for g in G:
    ms = [f"{g}:BOMP-PK", f"{g}:TS-SFW", f"{g}:MT-FW"]
    if g[:2] != "g3":
        ms += [f"{g}:BOMP-PK-W", f"{g}:TS-SFW-W"]
    pairs = [(f"{g}:TS-SFW", f"{g}:BOMP-PK"), (f"{g}:TS-SFW", "g3d3:BOMP-PK"), (f"{g}:BOMP-PK", "g3d3:BOMP-PK")]
    if g[:2] != "g3":
        pairs += [(f"{g}:TS-SFW-W", f"{g}:BOMP-PK-W"), (f"{g}:TS-SFW-W", "g3d3:BOMP-PK")]
    txt, *_ = A.table(rows, ms + (["g3d3:BOMP-PK"] if g != "g3d3" else []), pairs, title=f"C1 grid {g} ({NODES[g]} nodes)")
    out.append(txt)
    sec = {m: np.mean([r["sec"][m] for r in rows]) for m in ms}
    out.append("  mean s/window: " + ", ".join(f"{m.split(':')[1]} {v:.3f}" for m, v in sec.items()))
print("\n".join(out))
