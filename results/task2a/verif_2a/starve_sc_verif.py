"""Verifier experiments on regenerated 2A evaluation windows (deterministic, same RNG streams as run_job):
(a) picker starvation per true event; (b) GaMMA/PyOcto positive controls on true-time picks;
(c) oracle station-correction sanity for MT-FW, QM-3C and GaMMA."""
import os, sys, json, time
os.environ.setdefault("TASK2A_SCRATCH", "/tmp/claude-0/-home-claude/d3c612af-4e2f-52d9-b554-c7fa18d867f4/scratchpad/verif")
sys.path.insert(0, "/home/claude/work"); sys.path.insert(0, "/home/claude/work/gamma_src")
import numpy as np
import task2a_confirm as T
import task1x_mtglrt as X
import qm_wrapper as QM
from certified_fw import certified_fw_record
from froc_exact import exact_HF
from real_picker import match, pick_waveforms
from task0_common import DEFAULT_TAUS, VP, VS, default_grid
from task1o_consistent import matched, score_events
from task1r_baselines import GAMMA_GRID, PYOCTO_GRID, run_gamma
from task1r_fair import recall_at
from task1y_clutter import CLUT_AMP, degrade_y

seed, regime = int(sys.argv[1]), sys.argv[2]
jobs = [tuple(x.split(":")) for x in sys.argv[3].split(",")]          # cell:net
NW = int(sys.argv[4]) if len(sys.argv) > 4 else 20
out_path = f"/home/claude/work/verif_2a/starve_sc_{seed}_{regime}_verif.json"
grid = default_grid(); garr = np.asarray(grid)
qmc = T.load_qm_cfgs()
res = {"per_event": [], "windows": [], "consistency": []}
t0 = time.time()
for cell, j in jobs:
    j = int(j)
    part = json.load(open(f"/home/claude/work/task2a_parts/seed{seed}_{regime}_{cell}_net{j}.json"))
    ci, kind, nst, nclut = T.net_setup(seed, cell)
    STA, statics = T.network(seed, cell, j)
    corr = np.array(part["corr"]); expct = np.array(part["expected_corr"])
    assert np.allclose(STA, np.array(part["STA"])) and np.allclose(statics, np.array(part["statics"]))
    Gs = [X.radiation_operator(g, STA) for g in garr]; Qs = [X.basis(G) for G in Gs]
    lut = QM.build_lut(STA, VP, VS, {**QM.QM_CFG, **qmc["QM-3C"]})
    rngs = T.window_rngs(seed, cell, j, regime); rng_d = rngs[3]
    SCS = {"none": None, "est": corr, "static": statics, "static+mis": expct}
    for w in range(NW):
        srcs, taus, D, picks, tr, amp = T.generate_window(STA, statics, nclut, CLUT_AMP[regime], rngs)
        n2 = 2 * len(srcs)
        true_arr = [pl[:n2] for pl in picks]
        pk, _q = pick_waveforms(D, method="stalta")
        p_dy, _l = degrade_y(picks, amp, len(srcs), rng_d)
        p_dy = [[(t, b) for t, b in pl if t > -100] for pl in p_dy]
        W0 = part["windows"][w]
        same = (W0["truth"] == [list(map(float, s)) for s in srcs]) and (W0["picker"]["STALTA"] == T.pick_counts(true_arr, pk))
        res["consistency"].append(bool(same))
        # (a) per event matched picks
        for e in range(len(srcs)):
            n_lab = n_any = n_1s = 0; errs = []
            for i in range(nst):
                te = true_arr[i][2 * e:2 * e + 2]
                m = match(te, pk[i]); n_lab += len(m); errs += [x[2] for x in m]
                n_any += sum(1 for tt, _b in te if any(abs(tp - tt) <= 0.3 for tp, _bp in pk[i]))
                n_1s += sum(1 for tt, _b in te if any(abs(tp - tt) <= 1.0 for tp, _bp in pk[i]))
            res["per_event"].append({"cell": cell, "net": j, "w": w, "n_arr": 2 * nst, "matched_lab_0.3": n_lab,
                                     "matched_anylabel_0.3": n_any, "matched_anylabel_1.0": n_1s, "errs": errs})
        rec = {"k": len(srcs), "cell": cell, "net": j}
        # (b) GaMMA positive controls
        true_only = [list(pl) for pl in true_arr]
        hit_only = []
        for i in range(nst):
            m = match(true_arr[i], pk[i]); ks = {x[1] for x in m}
            hit_only.append([pk[i][k] for k in sorted(ks)])
        with T.silent():
            rec["G_true"] = [list(score_events(run_gamma(STA, true_only, **c), tr)) for c in GAMMA_GRID]
            rec["G_true_orSC"] = [list(score_events(run_gamma(STA, T.apply_sc(true_only, expct), **c), tr)) for c in GAMMA_GRID]
            rec["G_stalta_hitsonly"] = [list(score_events(run_gamma(STA, hit_only, **c), tr)) for c in GAMMA_GRID]
            rec["P_true"] = [list(score_events(T.pyocto_guarded(STA, true_only, **c), tr)) for c in PYOCTO_GRID]
            for nm, cc in SCS.items():
                P = p_dy if cc is None else T.apply_sc(p_dy, cc)
                rec[f"G_DY_{nm}"] = [list(score_events(run_gamma(STA, P, **c), tr)) for c in GAMMA_GRID]
        # (c) MT-FW and QM-3C with each correction
        MF, sd = matched(D); MFn = MF / sd
        for nm, cc in SCS.items():
            rec[f"MTFW_{nm}"] = certified_fw_record(STA, MFn, grid, DEFAULT_TAUS, Gs, Qs, tr, tcorr=cc)
        if nm == "static+mis":
            pass
        evs = T.run_qm_variants(STA, D, os.path.join(T.SCRATCH, "verif_qm"),
                                {f"QM3C_{nm}": (qmc["QM-3C"], cc) for nm, cc in SCS.items()}, lut=lut)
        for nm in SCS:
            rec[f"QM3C_{nm}"] = QM.qm_records(evs[f"QM3C_{nm}"], tr)
        # consistency with stored records
        def req(a, b):
            return len(a) == len(b) and all(abs(x[0] - y[0]) <= 1e-9 * max(1.0, abs(y[0])) and int(x[1]) == int(y[1]) for x, y in zip(a, b))
        res["consistency"].append(req(rec["MTFW_none"], W0["records"]["MT-FW"]))
        res["consistency"].append(req(rec["MTFW_est"], W0["records"]["MT-FW+SC"]))
        res["consistency"].append([list(x) for x in rec["QM3C_est"]] == W0["records"]["QM-3C+SC"])
        res["consistency"].append(rec["G_DY_est"] == W0["per"]["GaMMA-DY+SC"])
        rec = {k_: ([list(x) for x in v] if isinstance(v, list) else v) for k_, v in rec.items()}
        res["windows"].append(rec)
        print(cell, j, w, f"{time.time() - t0:.0f}s", res["consistency"][-5:], flush=True)
json.dump(res, open(out_path, "w"))
