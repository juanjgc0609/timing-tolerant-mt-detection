"""mech_analyze_verif.py - summary of mech_verif.py (exploratory, seed 999, networks 10-15, 10 dense windows per
network, cell and regime).  Exact pooled FROC (stats_verif.froc_point); CI: stats_verif.hier_boot (cells as strata,
networks, windows; B = 500) for the key contrasts."""
import json
import os

import numpy as np

from stats_verif import froc_point, hier_boot

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "mech_parts")
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
VAR = ("MTFW_cap2", "MTFW_cap3", "MTFW_cap6", "MTFW_cap12", "MTFW_cap24", "MTFW_nocap", "SELcap6_SCOREnocap",
       "SELnocap_SCOREcap6", "BOMP-PK")


def load(cond, regs):
    W = []
    for reg in regs:
        for cell in CELLS:
            for j in range(10, 16):
                p = os.path.join(OUT, f"{cond}_{reg}_{cell}_net{j}.json")
                if not os.path.exists(p):
                    continue
                e = json.load(open(p))
                for w in e["windows"]:
                    W.append({"cell": cell, "net": (reg, j) if len(regs) == 1 else j, "k": w["k"],
                              "rec": {k: [tuple(x) for x in v] for k, v in w["rec"].items() if not k.startswith("_")}})
    return W


checks = []
for reg in ("thin", "strong"):
    for cell in CELLS:
        for j in range(10, 16):
            p = os.path.join(OUT, f"base_{reg}_{cell}_net{j}.json")
            if os.path.exists(p):
                checks += [c for w in json.load(open(p))["windows"] if "_check" in w["rec"] for c in w["rec"]["_check"]]
print(f"re-implementation check vs c6_lib.mtfw_dense (cap 6 and uncapped), base condition: {sum(checks)}/{len(checks)} "
      "identical")
for cond in ("base", "noclut", "noclut_nosite"):
    for regs in (("thin",), ("strong",), ("thin", "strong")):
        W = load(cond, regs)
        if not W:
            continue
        p = {v: froc_point(W, v) for v in VAR}
        print(f"{cond} {'+'.join(regs)}: {len(W)} windows, {sum(w['k'] for w in W)} events")
        print("   " + "  ".join(f"{v} {p[v][0.1]:.3f}/{p[v][0.25]:.3f}" for v in VAR))
        print("   vs BOMP-PK (cap 0.1/0.25): " + "  ".join(
            f"{v} {p[v][0.1] - p['BOMP-PK'][0.1]:+.3f}/{p[v][0.25] - p['BOMP-PK'][0.25]:+.3f}" for v in VAR[:-1]))
        if len(regs) == 2:
            for a, b in (("MTFW_cap6", "MTFW_nocap"), ("MTFW_cap6", "BOMP-PK"), ("SELnocap_SCOREcap6", "MTFW_nocap"),
                         ("SELcap6_SCOREnocap", "MTFW_nocap")):
                d = hier_boot(W, a, b, B=500, seed=7)
                pt = p[a][0.25] - p[b][0.25]
                print(f"   {a} - {b} @0.25: {pt:+.3f} [{np.percentile(d, 2.5):+.3f}, {np.percentile(d, 97.5):+.3f}]")
