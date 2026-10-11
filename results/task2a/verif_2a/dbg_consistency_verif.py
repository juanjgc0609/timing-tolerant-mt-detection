import os, sys, json
sys.path.insert(0, "/home/claude/work"); sys.path.insert(0, "/home/claude/work/gamma_src")
import numpy as np
import task2a_confirm as T
from real_picker import pick_waveforms
from task1y_clutter import CLUT_AMP, degrade_y
seed, regime, cell, j = 2024, "thin", "ladder_6", 0
part = json.load(open(f"/home/claude/work/task2a_parts/seed{seed}_{regime}_{cell}_net{j}.json"))
ci, kind, nst, nclut = T.net_setup(seed, cell)
STA, statics = T.network(seed, cell, j)
rngs = T.window_rngs(seed, cell, j, regime)
for w in range(3):
    srcs, taus, D, picks, tr, amp = T.generate_window(STA, statics, nclut, CLUT_AMP[regime], rngs)
    n2 = 2 * len(srcs)
    pk, _ = pick_waveforms(D, method="stalta")
    p_dy, _ = degrade_y(picks, amp, len(srcs), rngs[3])
    W0 = part["windows"][w]
    print(W0["truth"], [list(map(float, s)) for s in srcs], W0["taus"], taus)
    print(W0["picker"]["STALTA"], T.pick_counts([pl[:n2] for pl in picks], pk))
