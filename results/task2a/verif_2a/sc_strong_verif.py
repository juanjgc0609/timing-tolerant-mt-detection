"""Exploratory: MT-FW on strong-clutter evaluation windows with no correction, the SC actually used (strong-clutter
catalogue), the SC estimated from the thin-clutter catalogue of the same network, the true statics, and the true
statics plus the median model mismatch."""
import os, sys, json, time
sys.path.insert(0, "/home/claude/work")
import task2a_confirm as T
import numpy as np
import task1x_mtglrt as X
from certified_fw import certified_fw_record
from froc_exact import exact_HF
from task0_common import DEFAULT_TAUS, default_grid
from task1o_consistent import matched
from task1r_fair import recall_at
from task1y_clutter import CLUT_AMP
seed, cells = int(sys.argv[1]), sys.argv[2].split(",")
grid = default_grid(); garr = np.asarray(grid)
out = {}
t0 = time.time()
for cell in cells:
    ci, kind, nst, nclut = T.net_setup(seed, cell)
    recs = {k: [] for k in ("none", "est_strong", "est_thin", "static", "static+mis")}
    ks, rms = [], {"strong": [], "thin": []}
    for j in range(T.N_NETS):
        ps = json.load(open(f"/home/claude/work/task2a_parts/seed{seed}_strong_{cell}_net{j}.json"))
        pt = json.load(open(f"/home/claude/work/task2a_parts/seed{seed}_thin_{cell}_net{j}.json"))
        STA, statics = T.network(seed, cell, j)
        rms["strong"].append(ps["resid_rms"]); rms["thin"].append(pt["resid_rms"])
        SC = {"none": None, "est_strong": np.array(ps["corr"]), "est_thin": np.array(pt["corr"]),
              "static": statics, "static+mis": np.array(ps["expected_corr"])}
        Gs = [X.radiation_operator(g, STA) for g in garr]; Qs = [X.basis(G) for G in Gs]
        rngs = T.window_rngs(seed, cell, j, "strong")
        for w in range(T.N_WIN):
            srcs, taus, D, picks, tr, amp = T.generate_window(STA, statics, nclut, CLUT_AMP["strong"], rngs)
            MF, sd = matched(D); MFn = MF / sd
            ks.append(len(srcs))
            for nm, cc in SC.items():
                recs[nm].append(certified_fw_record(STA, MFn, grid, DEFAULT_TAUS, Gs, Qs, tr, tcorr=cc))
            st = ps["windows"][w]["records"]["MT-FW+SC"]
            assert len(st) == len(recs["est_strong"][-1]) and all(abs(a[0] - b[0]) < 1e-9 * max(1, abs(b[0])) and a[1] == b[1] for a, b in zip(recs["est_strong"][-1], st))
    K = sum(ks)
    out[cell] = {nm: [round(recall_at(*exact_HF(r)[1:], K, c), 3) for c in (0.1, 0.25)] for nm, r in recs.items()}
    out[cell]["resid_rms_strong_cat"] = float(np.sqrt(np.mean(np.square(rms["strong"]))))
    out[cell]["resid_rms_thin_cat"] = float(np.sqrt(np.mean(np.square(rms["thin"]))))
    print(seed, cell, out[cell], f"{time.time() - t0:.0f}s", flush=True)
json.dump(out, open(f"/home/claude/work/verif_2a/sc_strong_{seed}_verif.json", "w"), indent=1)
