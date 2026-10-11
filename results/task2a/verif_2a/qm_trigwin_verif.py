"""Exploratory: QM-3C+SC with the frozen configuration versus the same configuration with the trigger (origin-time)
window widened by 0.5 s on each side, on the regenerated 2A evaluation windows (detect run once per window)."""
import os, sys, json, time, shutil, tempfile
os.environ.setdefault("TASK2A_SCRATCH", "/tmp/claude-0/-home-claude/d3c612af-4e2f-52d9-b554-c7fa18d867f4/scratchpad/verif")
sys.path.insert(0, "/home/claude/work")
import numpy as np
import task2a_confirm as T
import qm_wrapper as QM
from froc_exact import exact_HF
from task0_common import VP, VS
from task1r_fair import recall_at
from task1y_clutter import CLUT_AMP
seed, regime, cells = int(sys.argv[1]), sys.argv[2], sys.argv[3].split(",")
qmc = T.load_qm_cfgs()
cfg = {**QM.QM_CFG, **qmc["QM-3C"]}
tw = cfg["trigger_window"]
wide = {**cfg, "trigger_window": (tw[0] - 0.5, tw[1] + 0.5)}
out = {}
t0 = time.time()
for cell in cells:
    ci, kind, nst, nclut = T.net_setup(seed, cell)
    R = {"orig": [], "wide": [], "k": [], "same_as_stored": 0}
    for j in range(T.N_NETS):
        part = json.load(open(f"/home/claude/work/task2a_parts/seed{seed}_{regime}_{cell}_net{j}.json"))
        STA, statics = T.network(seed, cell, j)
        corr = np.array(part["corr"])
        lut = QM.build_lut(STA, VP, VS, cfg)
        lsc = T.qm_sc_lut(lut, corr, QM.station_names(nst))
        rngs = T.window_rngs(seed, cell, j, regime)
        for w in range(T.N_WIN):
            srcs, taus, D, picks, tr, amp = T.generate_window(STA, statics, nclut, CLUT_AMP[regime], rngs)
            os.makedirs(T.SCRATCH, exist_ok=True)
            work = tempfile.mkdtemp(prefix="qmtw_", dir=T.SCRATCH)
            try:
                QM.write_archive(D, STA, os.path.join(work, "archive"), cfg)
                QM.detect(STA, work, lsc, cfg)
                e1 = QM.run_trigger(work, lsc, cfg, name="orig")
                e2 = QM.run_trigger(work, lsc, wide, name="wide")
            finally:
                shutil.rmtree(work, ignore_errors=True)
            st = part["windows"][w]["qm_events"]["QM-3C+SC"]
            R["same_as_stored"] += int(len(st) == len(e1) and all(np.allclose(a, b, atol=2e-5) for a, b in zip(st, e1)))
            R["orig"].append(QM.qm_records(e1, tr)); R["wide"].append(QM.qm_records(e2, tr)); R["k"].append(len(srcs))
    K = sum(R["k"])
    res = {}
    for key in ("orig", "wide"):
        _t, H, F = exact_HF(R[key])
        res[key] = [round(recall_at(H, F, K, c), 4) for c in (0.1, 0.25)]
    out[cell] = {"recall_orig": res["orig"], "recall_wide": res["wide"], "same_as_stored": R["same_as_stored"], "n": len(R["k"])}
    print(seed, regime, cell, out[cell], f"{time.time() - t0:.0f}s", flush=True)
json.dump(out, open(f"/home/claude/work/verif_2a/qm_trigwin_{seed}_{regime}_verif.json", "w"), indent=1)
