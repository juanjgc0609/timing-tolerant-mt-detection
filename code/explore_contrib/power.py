"""Power calculation for the confirmatory designs, from the pilot bootstrap SD of the paired recall difference.
SE(n) = SD_pilot * sqrt(n_pilot / n) (networks per cell, windows per network fixed); one-sided alpha 0.025;
power(n, delta) = Phi(delta / SE(n) - z_{0.975}).  Planning effects: the smallest effect of interest (0.05) and half
of the pilot point estimate (a shrinkage against the winner's curse), whichever is smaller.
    python power.py out/c2_ablate_thin.json[,out/c2_ablate_strong.json] MT-FW BOMP-PK"""
import json
import sys

import numpy as np
from scipy.stats import norm

import xc


def boot_diff(rows, a, b, B=2000, seed=3):
    keys = sorted({(r["cell"], r["net"]) for r in rows})
    idx = {k: [i for i, r in enumerate(rows) if (r["cell"], r["net"]) == k] for k in keys}
    cells = sorted({k[0] for k in keys})
    ks = [r["k"] for r in rows]
    ra = [[tuple(x) for x in r["rec"][a]] for r in rows]
    rb = [[tuple(x) for x in r["rec"][b]] for r in rows]

    def d(ws):
        n = sum(ks[w] for w in ws)
        return np.subtract(xc.recall_caps([ra[w] for w in ws], n)[0], xc.recall_caps([rb[w] for w in ws], n)[0])
    pt = d(list(range(len(rows))))
    rng = np.random.default_rng(seed)
    reps = []
    for _ in range(B):
        ws = []
        for c in cells:
            nl = [k for k in keys if k[0] == c]
            for g in rng.integers(0, len(nl), len(nl)):
                gg = idx[nl[g]]
                ws.extend(np.asarray(gg)[rng.integers(0, len(gg), len(gg))].tolist())
        reps.append(d(ws))
    reps = np.array(reps)
    return pt, reps.std(axis=0), np.percentile(reps, [2.5, 97.5], axis=0), len(keys) // len(cells)


def main():
    rows = []
    for p in sys.argv[1].split(","):
        rr = json.load(open(p))
        rr = rr["rows"] if isinstance(rr, dict) else rr
        tag = p.split("_")[-1].split(".")[0]
        for r in rr:
            r = dict(r)
            r["cell"] = tag + "/" + r["cell"]
            rows.append(r)
    a, b = sys.argv[2], sys.argv[3]
    pt, sd, ci, npc = boot_diff(rows, a, b)
    print(f"{a} - {b}: pilot {len(rows)} windows, {npc} networks per stratum")
    for c, cap in enumerate((0.1, 0.25)):
        print(f"  cap {cap}: point {pt[c]:+.3f}, 95% CI [{ci[0][c]:+.3f}, {ci[1][c]:+.3f}], bootstrap SD {sd[c]:.4f}")
        for delta in sorted({0.05, round(pt[c] / 2, 3)}):
            if delta <= 0:
                continue
            for pw in (0.8, 0.9):
                n = npc * (sd[c] * (norm.ppf(0.975) + norm.ppf(pw)) / delta) ** 2
                print(f"    delta {delta:.3f}, power {pw:.0%}: {int(np.ceil(n))} networks per stratum "
                      f"(x {len(rows) // (npc * len({r['cell'] for r in rows}))} windows)")


if __name__ == "__main__":
    main()
