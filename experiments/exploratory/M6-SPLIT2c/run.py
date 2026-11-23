"""M6-SPLIT2c: M6-SPLIT2 and the director's window-only constraint again, after the audit of 23 November 2026 found
two implementation errors (card.md here, declared before the run). Uses the scene code of M6-SPLIT2b (imported) with
new seeds, the exact (Fourier) delay of the slow S component in the generator (gen_fix.py) and the free one-event
alternative started also from the atom's own times (m6split free_atom_start=True).

    PRECHECK=1 python run.py > precheck_out.txt                  # oracle and no-splitting pre-checks
    ARMS=i,w PART=all python run.py > run_i_w.txt
    ARMS=iii PART=ctrl python run.py > run_iii_ctrl.txt          # each job under the 2 h limit
    ARMS=iii PART=pair python run.py > run_iii_pair.txt
    COMBINE=1 python run.py > run_out.txt                        # reads the saved parts, prints the decision
"""
import importlib.util
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
spec = importlib.util.spec_from_file_location("ms2b", HERE.parent / "M6-SPLIT2b" / "run.py")
B = importlib.util.module_from_spec(spec)
sys.modules["ms2b"] = B
spec.loader.exec_module(B)
import numpy as np  # noqa: E402
import gen_fix  # noqa: E402
from ttmt.v2 import m6, m6split  # noqa: E402

gen_fix.install(B)
B.SEED_C, B.SEED_P, B.SEED_SPLIT, B.SEED_ISO = 20261151, 20261152, 20261153, 20261154   # declared in card.md
SEED_PRE = 20261155
B.SPLIT_ARM.update({"w": True, "orc": True, "pre": False})
B.MISMATCH.update({"w": True, "orc": True, "pre": True})
NCTRL, NPAIR, PROCS, THR = 200, 40, 4, B.THR_REF
ARMS = tuple(os.environ.get("ARMS", "i,w,iii").split(","))
PART = os.environ.get("PART", "all")
R, MB = B.R, B.MB


def loops(arm, Y, q, thr, iso, fs=None):
    if arm == "i":
        return m6.adcg_loop(R.E, Y, q, R.KMAX, R.LAG, thr=thr, free="joint", near_km=R.NEAR_KM)
    if arm == "orc":                                   # pre-check: true splitting held in both hypotheses
        return m6split.adcg_loop(R.E, Y, q, R.KMAX, R.LAG, thr=thr, near_km=R.NEAR_KM, mode="fixed", fixed_split=fs,
                                 free_atom_start=True)
    return m6split.adcg_loop(R.E, Y, q, R.KMAX, R.LAG, thr=thr, near_km=R.NEAR_KM, mode="station",
                             iso_windows=(iso or None) if arm in ("iii", "pre") else None,
                             converge_iso_fits=arm in ("iii", "pre"), free_atom_start=True)


def job_pre(args):
    arm, c = args
    key = (SEED_PRE, 1, c)
    rs = np.random.default_rng([B.SEED_SPLIT] + list(key))
    phi, dt = rs.uniform(0, np.pi, 8), rs.uniform(0.02, 0.2, 8)
    srcs, Y, q, iso = B.scene(key, "C1", arm)
    b = loops(arm, Y, q, -np.inf, iso, (phi, dt / B.K.dt))
    g = [x["pair_gain"] for x in MB.spurious(b, srcs) if x["pair_gain"] is not None]
    return max(g) if g else None


def job_ctrl(args):
    m6split.NONCONV[0] = 0
    arm, c = args
    srcs, Y, q, iso = B.scene((B.SEED_C, 1, c), "C1", arm)
    a = loops(arm, Y, q, THR, iso)
    b = loops(arm, Y, q, -np.inf, iso)
    g = [x["pair_gain"] for x in MB.spurious(b, srcs) if x["pair_gain"] is not None]
    return bool(len(MB.spurious(a, srcs)) > 0), (max(g) if g else None), m6split.NONCONV[0]


def job_pair(args):
    m6split.NONCONV[0] = 0
    arm, dtau, c = args
    srcs, Y, q, iso = B.scene((B.SEED_P, int(100 * dtau), c), "pair", arm, dtau)
    a = loops(arm, Y, q, THR, iso)
    b = loops(arm, Y, q, -np.inf, iso)
    sp = {id(x) for x in MB.spurious(b, srcs)}
    return (bool(len(m6.hits_loctv(a, srcs)) == 2),
            [x["pair_gain"] for x in b if x["pair_gain"] is not None and id(x) not in sp], m6split.NONCONV[0])


def report(arm, ctrl, pairs):
    k = sum(x[0] for x in ctrl)
    g_sp = [x[1] for x in ctrl if x[1] is not None]
    g_t = [g for d in R.DTAUS for o in pairs[str(d)] for g in o[1]]
    rec = {d: float(np.mean([o[0] for o in pairs[str(d)]])) for d in R.DTAUS}
    hits = int(sum(o[0] for d in R.DTAUS for o in pairs[str(d)]))
    lo, hi = R.cp(k, NCTRL)
    nc = sum(x[2] for x in ctrl) + sum(o[2] for d in R.DTAUS for o in pairs[str(d)])
    print(f"   isolated fits that did not converge in arm {arm}: {nc}")
    print(f"arm {arm:4s}: control scenes split at thr {THR}: {k}/{NCTRL} [{lo:.3f}, {hi:.3f}]; near-pair recall "
          + ", ".join(f"{d:.2f}s {v:.3f}" for d, v in rec.items())
          + f", pooled {hits / (len(R.DTAUS) * NPAIR):.3f} ({hits}/{len(R.DTAUS) * NPAIR}); AUC {B.auc(g_t, g_sp):.3f}"
          f" (true gains n {len(g_t)}, median {np.median(g_t):.1f}; largest spurious gain per scene n {len(g_sp)}, "
          f"median {np.median(g_sp):.1f}, max {max(g_sp):.1f})", flush=True)
    return dict(k=k, hits=hits, auc=B.auc(g_t, g_sp))


def decision(res):
    n = len(R.DTAUS) * NPAIR
    for arm in ("w", "iii"):
        if arm not in res or "i" not in res:
            continue
        r, dh = res[arm], res["i"]["hits"] - res[arm]["hits"]
        val = "passes" if r["k"] <= 12 else "FAILS"
        zone = ("passes" if 100 * dh < 5 * n else "useful, owner decides" if 100 * dh <= 15 * n else "COSTS TOO MUCH")
        print(f"DECISION arm {arm} (card.md): validity {r['k']} of {NCTRL} (limit 12): {val}; recall lost {dh} of {n} "
              f"against arm i ({dh / n:+.4f}): {zone}; AUC {r['auc']:.3f}: {'passes' if r['auc'] >= 0.90 else 'FAILS'}")


def main():
    t0 = time.time()
    print(f"M6-SPLIT2c  seeds {B.SEED_C}/{B.SEED_P}/{B.SEED_SPLIT}/{B.SEED_ISO} (pre-check {SEED_PRE}), NISO {B.NISO}, "
          f"exact slow-S delay, free_atom_start, thr {THR}, numpy {np.__version__}", flush=True)
    if os.environ.get("PRECHECK") == "1":
        with mp.Pool(PROCS) as pool:
            for arm, lab in (("orc", "split data, TRUE splitting held (oracle)"), ("pre", "no splitting, model of arm iii")):
                g = [x for x in pool.map(job_pre, [(arm, c) for c in range(20)], chunksize=1) if x is not None]
                med = float(np.median(g))
                print(f"PRE-CHECK {arm} ({lab}): 20 scenes, {len(g)} with a spurious near atom; median largest spurious "
                      f"gain {med:.1f}; > {THR} in {sum(x > THR for x in g)}; values {np.round(sorted(g), 1).tolist()} -> "
                      f"{'ok' if -200 <= med <= -50 else 'STOP: diagnose before running'}", flush=True)
        print(f"{time.time() - t0:.0f} s")
        return
    if os.environ.get("COMBINE") == "1":
        res = {}
        for arm in ("i", "w", "iii"):
            parts = {}
            for f in sorted(HERE.glob(f"part_{arm}_*.json")):
                parts.update(json.loads(f.read_text()))
            if "ctrl" in parts and "pair" in parts:
                res[arm] = report(arm, parts["ctrl"], parts["pair"])
        decision(res)
        return
    with mp.Pool(PROCS) as pool:
        for arm in ARMS:
            out = {}
            if PART in ("all", "ctrl"):
                out["ctrl"] = pool.map(job_ctrl, [(arm, c) for c in range(NCTRL)], chunksize=1)
            if PART in ("all", "pair"):
                out["pair"] = {str(d): pool.map(job_pair, [(arm, d, c) for c in range(NPAIR)], chunksize=1)
                               for d in R.DTAUS}
            (HERE / f"part_{arm}_{PART}.json").write_text(json.dumps(out))
            print(f"arm {arm} part {PART} saved; {time.time() - t0:.0f} s", flush=True)
            if PART == "all":
                report(arm, out["ctrl"], out["pair"])
    print(f"{time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
