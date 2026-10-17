"""V5: bootstraps that respect the shared noise hours of the cells (seed 2718, DENSE, LOC, @0.25).
 a) cluster = network index j (the same noise block in all 3 cells): resample the 24 j with replacement, take all three
    cells' networks of each drawn j, then windows with replacement within each (cell, j);
 b) moving-block version: resample contiguous runs of 4 consecutive j (hours are temporally correlated), 6 runs;
 c) registered-type bootstrap (strata = cell) for reference, same seed."""
import sys

import numpy as np

sys.path.insert(0, "/home/claude/work/t2m_verif")
import vlib_verif as V  # noqa: E402

OUT = "/home/claude/work/t2m_verif/v5_hourboot_verif.txt"
B = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
lines = []


def P(s=""):
    print(s, flush=True)
    lines.append(s)
    with open(OUT, "w") as fh:
        fh.write("\n".join(lines) + "\n")


A, H, Bm = "MT@S0", "HUBER-N@S0", "BOMP@S0"
for seed, j0, nj in (("2718", 2800, 24), ("3141", 2900, 17)):
    wins = V.windows(V.load_parts(seed), "DENSE")
    prep = {m: V.prep(wins, m) for m in (A, H, Bm)}
    grp = {}
    for i, w in enumerate(wins):
        grp.setdefault((w["cell"], w["net"]), []).append(i)
    grp = {k: np.array(v) for k, v in grp.items()}
    js = list(range(j0, j0 + nj))
    rng = np.random.default_rng(31415)

    def draw_j(jlist):
        ws = []
        for j in jlist:
            for c in V.CELLS:
                g = grp.get((c, j))
                if g is not None:
                    ws.append(g[rng.integers(0, len(g), len(g))])
        return np.concatenate(ws)

    schemes = {}
    schemes["a) cluster j"] = [draw_j(rng.choice(js, nj, replace=True)) for _ in range(B)]
    L = 4
    starts = list(range(j0, j0 + nj - L + 1))
    nb = int(np.ceil(nj / L))
    schemes["b) moving blocks of 4 j"] = [draw_j([s + k for s in rng.choice(starts, nb) for k in range(L)][:nj])
                                         for _ in range(B)]
    schemes["c) strata = cell (registered type)"] = V.boot_index(wins, B, 31415)
    pt = {m: V.recall_at(*prep[m], np.arange(len(wins)))[1] for m in (A, H, Bm)}
    P(f"seed {seed}: points MT-BOMP {pt[A] - pt[Bm]:+.3f}, HUB-BOMP {pt[H] - pt[Bm]:+.3f}; B = {B}")
    for nm, reps in schemes.items():
        R = {m: np.array([V.recall_at(*prep[m], ws)[1] for ws in reps]) for m in (A, H, Bm)}
        for a in (A, H):
            d = R[a] - R[Bm]
            P(f"   {nm:36s} {a:10s} - BOMP: 95% [{np.percentile(d, 2.5):+.3f}, {np.percentile(d, 97.5):+.3f}] "
              f"sd {d.std():.4f} p(d<=0) {(1 + np.sum(d <= 0)) / (B + 1):.4f}")
