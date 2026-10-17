"""V6: reproduce (frozen code, read-only use) the catalogue-pool exhaustion of the two failed 3141 networks, and the
catalogue draws of two completed networks as a control."""
import sys
sys.path.insert(0, "/home/claude/work")
import t2m_lib as L2  # noqa: E402
import json

out = []
L2.configure_file("/home/claude/work/t2m_out/cfg_3141.json")
for cell, j in (("ladder_4", 2909), ("ladder_4", 2911), ("ladder_4", 2910), ("ladder_6", 2909)):
    try:
        n = L2.Net2M(cell, j, "REAL")
        d = json.load(open(f"/home/claude/work/t2m_out/parts_3141/REAL_{cell}_{j}.json"))
        out.append(f"{cell} {j}: completed, drawn {n.drawn} (part file: {d['cat_drawn']}), cat slots equal: "
                   f"{n.cat_slots == d['cat_slots']}, pool {len(L2.cat_pool(n.sl))}")
    except RuntimeError as e:
        out.append(f"{cell} {j}: RuntimeError: {e}; pool size {len(L2.cat_pool(L2.slots(j)))}")
L2.configure_file("/home/claude/work/t2m_out/cfg_2718.json")
mx = 0
for cell in L2.CELLS:
    for j in L2.CFG["nets"]:
        mx = max(mx, len(L2.cat_pool(L2.slots(j))))
out.append(f"2718 pool sizes: min {min(len(L2.cat_pool(L2.slots(j))) for j in L2.CFG['nets'])}")
open("/home/claude/work/t2m_verif/v6_failednets_verif.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
