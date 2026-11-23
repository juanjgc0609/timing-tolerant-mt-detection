"""M6-SPLIT-AUDIT (card.md here): known-truth recovery, oracle negative control, oracle positive control.

    PART=1 python run.py > part1_out.txt
    PART=2 FIX=1 python run.py > part2_fix_out.txt
    PART=2 FIX=0 python run.py > part2_nofix_out.txt
    PART=3 python run.py > part3_out.txt
"""
import importlib.util
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "M6-SPLIT2c"))
spec = importlib.util.spec_from_file_location("ms2b", HERE.parent / "M6-SPLIT2b" / "run.py")
B = importlib.util.module_from_spec(spec)
sys.modules["ms2b"] = B
spec.loader.exec_module(B)
import numpy as np  # noqa: E402
import gen_fix  # noqa: E402
from ttmt.v2 import m6, m6split  # noqa: E402

gen_fix.install(B)
SEED_REC, SEED_C, SEED_P, SEED_SPLIT = 20261161, 20261162, 20261163, 20261164     # declared in card.md
B.SEED_SPLIT, B.SEED_ISO = SEED_SPLIT, SEED_REC + 1000
B.SPLIT_ARM.update({"orc": True, "rec": True, "recm": True, "i": False})
B.MISMATCH.update({"orc": True, "rec": True, "recm": False, "i": True})
R, MB, K = B.R, B.MB, B.K
PART, FIX = os.environ.get("PART", "1"), os.environ.get("FIX", "1") == "1"
THR, PROCS = B.THR_REF, 4


def truth(key):
    rs = np.random.default_rng([SEED_SPLIT] + list(key))
    return rs.uniform(0, np.pi, 8), rs.uniform(0.02, 0.2, 8)


def phi_err(a, b):
    d = np.abs((a - b) % np.pi)
    return np.degrees(np.minimum(d, np.pi - d))


def job_rec(args):
    arm, c = args
    key = (SEED_REC, 1, c)
    phi, dt = truth(key)
    srcs, Y, q, _iso = B.scene(key, "C1", arm)
    resp, frames, dmax = m6._Resp(K), m6split.Frames(R.E.XN, R.E.STA), m6split.DMAX_S / K.dt
    aw = m6.adcg_loop(R.E, Y, q, 1, R.LAG, thr=None)
    for a in aw:
        a.update(anchor=a["T"].copy(), bound=np.full(a["T"].shape, 2.0), phi=np.zeros(8), dl=np.zeros(8))
    at, _R, _r = m6split.joint_refit(Y, aw, R.E.G, resp, frames, "station", dmax)
    e_d = np.abs(at[0]["dl"] * K.dt - dt)
    e_p = phi_err(at[0]["phi"], phi)
    # descriptive: 10 isolated events per station (arm iii setting), estimated jointly
    iso = []
    for w in range(10):
        B.SPLIT.update(on=True, phi=phi, dt=dt)
        try:
            _s, Yw, qw = R._scene((SEED_REC + 2000, 1, c, w), "C1", B.MISMATCH[arm])
        finally:
            B.SPLIT["on"] = False
        a_ = m6.adcg_loop(R.E, Yw, qw, 1, R.LAG, thr=None)
        for x in a_:
            x.update(anchor=x["T"].copy(), bound=np.full(x["T"].shape, 2.0), phi=np.zeros(8), dl=np.zeros(8))
        iso.append((np.asarray(Yw, float), a_))
    m6split.joint_refit(Y, [], R.E.G, resp, frames, "station", dmax, iso=iso)
    e_d10 = np.abs(iso[0][1][0]["dl"] * K.dt - dt)
    e_p10 = phi_err(iso[0][1][0]["phi"], phi)
    return dt.tolist(), e_d.tolist(), e_p.tolist(), e_d10.tolist(), e_p10.tolist()


def loops(arm, Y, q, thr, fs=None):
    if arm == "i":
        return m6.adcg_loop(R.E, Y, q, R.KMAX, R.LAG, thr=thr, free="joint", near_km=R.NEAR_KM)
    return m6split.adcg_loop(R.E, Y, q, R.KMAX, R.LAG, thr=thr, near_km=R.NEAR_KM, mode="fixed", fixed_split=fs,
                             free_atom_start=FIX)


def job_ctrl(c):
    key = (SEED_C, 1, c)
    phi, dt = truth(key)
    srcs, Y, q, _ = B.scene(key, "C1", "orc")
    a = loops("orc", Y, q, THR, (phi, dt / K.dt))
    b = loops("orc", Y, q, -np.inf, (phi, dt / K.dt))
    g = [x["pair_gain"] for x in MB.spurious(b, srcs) if x["pair_gain"] is not None]
    return bool(len(MB.spurious(a, srcs)) > 0), (max(g) if g else None)


def job_pair(args):
    arm, dtau, c = args
    key = (SEED_P, int(100 * dtau), c)
    phi, dt = truth(key)
    srcs, Y, q, _ = B.scene(key, "pair", arm, dtau)
    a = loops(arm, Y, q, THR, (phi, dt / K.dt))
    return bool(len(m6.hits_loctv(a, srcs)) == 2)


def main():
    t0 = time.time()
    print(f"M6-SPLIT-AUDIT part {PART} (fix {FIX}); seeds {SEED_REC}/{SEED_C}/{SEED_P}/{SEED_SPLIT}; exact slow-S delay; "
          f"thr {THR}; numpy {np.__version__}", flush=True)
    with mp.Pool(PROCS) as pool:
        if PART == "1":
            for arm, lab in (("rec", "mismatched"), ("recm", "matched")):
                out = pool.map(job_rec, [(arm, c) for c in range(100)], chunksize=1)
                dt = np.array([o[0] for o in out]).ravel()
                for j, name in ((1, "one event"), (3, "10 isolated events (descriptive)")):
                    ed, ep = np.array([o[j] for o in out]).ravel(), np.array([o[j + 1] for o in out]).ravel()
                    ok = (ed <= 0.02) & (ep <= 15.0)
                    big = dt >= 0.05
                    print(f"recovery, {lab}, {name}: delta err <= 0.02 s and phi err <= 15 deg at {ok.sum()}/{ok.size} "
                          f"= {ok.mean():.3f} of stations (pass >= 0.80: {'PASSES' if ok.mean() >= 0.80 else 'FAILS'}"
                          f"{'' if j == 1 else ', not a criterion'}); delta >= 0.05 s: {ok[big].mean():.3f}; median "
                          f"delta err {np.median(ed):.3f} s, median phi err {np.median(ep):.1f} deg; delta estimated 0 "
                          f"at {(np.array([o[j] for o in out]).ravel() == dt).sum()} stations", flush=True)
        elif PART == "2":
            out = pool.map(job_ctrl, range(200), chunksize=1)
            k = sum(o[0] for o in out)
            g = [o[1] for o in out if o[1] is not None]
            lo, hi = R.cp(k, 200)
            print(f"negative control, oracle (true splitting held), fix {FIX}: {k}/200 control scenes split at thr {THR}"
                  f" [{lo:.3f}, {hi:.3f}] -> {'<= 12: the model does not invent pairs with the truth' if k <= 12 else '> 12: the model invents pairs even with the truth'}"
                  f"; largest spurious gain median {np.median(g):.1f}, max {max(g):.1f}", flush=True)
        elif PART == "3":
            h = {}
            for arm in ("i", "orc"):
                h[arm] = {d: sum(pool.map(job_pair, [(arm, d, c) for c in range(40)], chunksize=1)) for d in R.DTAUS}
                print(f"pairs arm {arm}: " + ", ".join(f"{d:.2f}s {v}/40" for d, v in h[arm].items())
                      + f"; pooled {sum(h[arm].values())}/160", flush=True)
            dh = sum(h["i"].values()) - sum(h["orc"].values())
            print(f"positive control, oracle: recall lost {dh} of 160 against arm i ({dh / 160:+.4f}) -> "
                  f"{'PASSES' if 100 * dh < 5 * 160 else 'FAILS'} (within 0.05)", flush=True)
    print(f"{time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
