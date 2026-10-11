"""capfine_analyze_verif.py - CAPV sweep summary on the seed-999 base dense windows (exploratory)."""
import json, os
from stats_verif import froc_point
from mech_analyze_verif import load, OUT, CELLS
W = load("base", ("thin", "strong"))
F = []
for reg in ("thin", "strong"):
    for cell in CELLS:
        for j in range(10, 16):
            F.extend(json.load(open(os.path.join(OUT, f"capfine_{reg}_{cell}_net{j}.json")))["windows"])
assert len(F) == len(W) and all(a["k"] == b["k"] for a, b in zip(W, F))
for w, f in zip(W, F):
    w["rec"].update({k: [tuple(x) for x in v] for k, v in f["rec"].items()})
b = froc_point(W, "BOMP-PK")
print(f"base, thin+strong, 360 windows; BOMP-PK {b[0.1]:.3f}/{b[0.25]:.3f}")
for c in ("2", "3", "4", "5", "6", "8", "12", "24"):
    p = froc_point(W, f"MTFW_cap{c}")
    print(f"  CAPV {c:>2s}: MT-FW {p[0.1]:.3f}/{p[0.25]:.3f}  minus BOMP-PK {p[0.1] - b[0.1]:+.3f}/{p[0.25] - b[0.25]:+.3f}")
p = froc_point(W, "MTFW_nocap")
print(f"  no cap : MT-FW {p[0.1]:.3f}/{p[0.25]:.3f}  minus BOMP-PK {p[0.1] - b[0.1]:+.3f}/{p[0.25] - b[0.25]:+.3f}")
