"""C6 per cell and regime: the c6_analyze block (recall, differences, headroom fraction, hit errors) for each cell."""
import sys
import c6_analyze as A
out = []
for path in sys.argv[1:]:
    d = A.load(path)
    reg = d["nets"][0]["regime"]
    for cell in ("ladder_4", "ladder_6", "unilateral_c40"):
        rows = [r for r in d["rows"] if r["cell"] == cell]
        A.block(rows, f"{reg} / {cell}", out, per_cell=False)
print("\n".join(out))
