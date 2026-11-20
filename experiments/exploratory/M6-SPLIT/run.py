"""M6-SPLIT: S-wave splitting as a constraint inside the joint fit of M6 (director's request, 20 November 2026;
design in the supervision note of that day, section 3; card.md here, declared before the run). Simulation only.

Arms (same sources and noise; splitting parameters from their own random stream):
  i    no splitting in the data, model without it (ttmt.v2.m6.adcg_loop)
  ii   split S in the data (per station delta ~ U[0.02, 0.2] s, fast axis uniform), model without it (today)
  iii  the same data, model with the constraint (ttmt.v2.m6split, mode "station": phi_i, delta_i <= 0.2 s per
       station, common to all atoms, estimated in both hypotheses)
  iv   the same data, splitting free per event (mode "event"; descriptive, as exploration X1)
Scenes: M6-CAL settings (8-station ring, homogeneous medium, mismatched Ricker pulses, SNR 12 per-block peak,
kmax 4, near_km 10, LOC+T-V). Controls: NCTRL single-source scenes (C1). Pairs: near stratum (1.5-5 km),
dtau 0.1, 0.25, 0.5, 1.0 s, NPAIR per cell.

    python run.py [NCTRL] [NPAIR] [PROCS] > run_out.txt
"""
import importlib.util
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
NCTRL = int(sys.argv[1]) if len(sys.argv) > 1 else 200
NPAIR = int(sys.argv[2]) if len(sys.argv) > 2 else 40
PROCS = int(sys.argv[3]) if len(sys.argv) > 3 else 4
sys.argv = ["run.py", "200", "200", "20", "4", "op", "10"]          # M6-CAL settings (near_km 10)
spec = importlib.util.spec_from_file_location("m6cal", HERE.parent / "M6-CAL" / "run.py")
R = importlib.util.module_from_spec(spec)
sys.modules["m6cal"] = R
spec.loader.exec_module(R)
S, J, MB = R.S, R.J, R.MB

from ttmt.kernel import LEGACY as K  # noqa: E402
from ttmt.radiation import radiation_operator  # noqa: E402
from ttmt.v2 import m6, m6split  # noqa: E402

SEED_C, SEED_P, SEED_SPLIT = 20261121, 20261122, 20261123          # new seeds, never used
THR_REF = 23.8
ARMS = ("i", "ii", "iii", "iv")
_FIELD0 = S.field
SPLIT = {"on": False, "phi": None, "dt": None}


def field(STA, src, tau, m, stat):
    """S.field with the S pulse split per station when SPLIT["on"] (fast axis at phi in the plane normal to the ray,
    slow component delayed by dt); one pulse shape per arrival (the mismatched pulse is drawn once and shifted)."""
    if not SPLIT["on"]:
        return _FIELD0(STA, src, tau, m, stat)
    g = (radiation_operator(src, STA) @ m).reshape(2, len(STA), 3)
    d = np.linalg.norm(STA - src, axis=1)
    tarr = np.stack([d / S.VP, d / S.VS]) + tau + stat
    t = np.arange(S.NT) * S.DT
    D = np.zeros((len(STA), 3, S.NT))
    for i in range(len(STA)):
        D[i] += g[0, i][:, None] * S.ricker(t - tarr[0, i], K.pulse_fc)[None]
        n = (STA[i] - src) / d[i]
        e1 = np.cross([0, 0, 1.0], n)
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(n, e1)
        ph, dt = SPLIT["phi"][i], SPLIT["dt"][i]
        ef = np.cos(ph) * e1 + np.sin(ph) * e2
        es = -np.sin(ph) * e1 + np.cos(ph) * e2
        w = S.ricker(t - tarr[1, i], K.pulse_fc)
        D[i] += (g[1, i] @ ef) * ef[:, None] * w[None] + (g[1, i] @ n) * n[:, None] * w[None]
        D[i] += (g[1, i] @ es) * es[:, None] * np.interp(t - dt, t, w, left=0.0, right=0.0)[None]
    MF = np.array([[np.convolve(D[i, c], K.taps[::-1], mode="same") for c in range(3)] for i in range(len(STA))])
    return MF, g, tarr


S.field = field


def scene(key, kind, split, dtau=None):
    rs = np.random.default_rng([SEED_SPLIT] + list(key))
    SPLIT.update(on=split, phi=rs.uniform(0, np.pi, 8), dt=rs.uniform(0.02, 0.2, 8))
    try:
        return R._scene(key, kind, True, dtau, J.STRATA["near"])
    finally:
        SPLIT["on"] = False


def loops(arm, Y, q, thr):
    if arm in ("i", "ii"):
        return m6.adcg_loop(R.E, Y, q, R.KMAX, R.LAG, thr=thr, free="joint", near_km=R.NEAR_KM)
    return m6split.adcg_loop(R.E, Y, q, R.KMAX, R.LAG, thr=thr, near_km=R.NEAR_KM,
                             mode="station" if arm == "iii" else "event")


def job_ctrl(args):
    arm, c = args
    srcs, Y, q = scene((SEED_C, 1, c), "C1", arm != "i")
    a = loops(arm, Y, q, THR_REF)
    b = loops(arm, Y, q, -np.inf)
    g = [x["pair_gain"] for x in MB.spurious(b, srcs) if x["pair_gain"] is not None]
    return len(MB.spurious(a, srcs)) > 0, (max(g) if g else None)


def job_pair(args):
    arm, dtau, c = args
    srcs, Y, q = scene((SEED_P, int(100 * dtau), c), "pair", arm != "i", dtau)
    a = loops(arm, Y, q, THR_REF)
    b = loops(arm, Y, q, -np.inf)
    sp = {id(x) for x in MB.spurious(b, srcs)}
    return len(m6.hits_loctv(a, srcs)) == 2, [x["pair_gain"] for x in b if x["pair_gain"] is not None
                                              and id(x) not in sp]


def auc(t, s):
    t, s = np.asarray(t, float), np.asarray(s, float)
    if not t.size or not s.size:
        return float("nan")
    return float(np.mean(t[:, None] > s[None, :]) + 0.5 * np.mean(t[:, None] == s[None, :]))


def main():
    t0 = time.time()
    print(f"M6-SPLIT  controls {NCTRL} (one source), pairs {NPAIR} per cell (near), seeds {SEED_C}/{SEED_P}/"
          f"{SEED_SPLIT}, thr_ref {THR_REF}, delta <= {m6split.DMAX_S} s, numpy {np.__version__}", flush=True)
    res = {}
    with mp.Pool(PROCS) as pool:
        for arm in ARMS:
            ctrl = pool.map(job_ctrl, [(arm, c) for c in range(NCTRL)], chunksize=1)
            pairs = {d: pool.map(job_pair, [(arm, d, c) for c in range(NPAIR)], chunksize=1) for d in R.DTAUS}
            k = sum(x[0] for x in ctrl)
            g_sp = [x[1] for x in ctrl if x[1] is not None]
            g_t = [g for d in R.DTAUS for o in pairs[d] for g in o[1]]
            rec = {d: float(np.mean([o[0] for o in pairs[d]])) for d in R.DTAUS}
            hits = int(sum(o[0] for d in R.DTAUS for o in pairs[d]))
            pooled = hits / (len(R.DTAUS) * NPAIR)
            lo, hi = R.cp(k, NCTRL)
            res[arm] = dict(k=k, rec=rec, pooled=pooled, hits=hits, auc=auc(g_t, g_sp))
            print(f"arm {arm:3s}: control scenes split at thr {THR_REF}: {k}/{NCTRL} [{lo:.3f}, {hi:.3f}]; "
                  f"near-pair recall " + ", ".join(f"{d:.2f}s {v:.3f}" for d, v in rec.items())
                  + f", pooled {pooled:.3f}; AUC {res[arm]['auc']:.3f} (true gains n {len(g_t)}, median "
                  f"{np.median(g_t) if g_t else np.nan:.1f}; largest spurious gain per scene n {len(g_sp)}, median "
                  f"{np.median(g_sp) if g_sp else np.nan:.1f}, max {max(g_sp) if g_sp else np.nan:.1f}); "
                  f"{time.time() - t0:.0f} s", flush=True)
    r3, r1 = res["iii"], res["i"]
    npair = len(R.DTAUS) * NPAIR
    dh = r1["hits"] - r3["hits"]                  # exact zone boundaries on counts (no floating-point comparison)
    gap = dh / npair
    zone = "PASS" if 100 * dh < 5 * npair else ("USEFUL, owner decides" if 100 * dh <= 15 * npair else "COSTS TOO MUCH")
    print(f"\nDECISION (card.md): validity, arm iii {r3['k']} of {NCTRL} split (limit 12 of 200): "
          f"{'ok' if r3['k'] <= 12 and NCTRL == 200 else 'FAILS' if r3['k'] > 12 else 'n != 200'}; "
          f"recall gap i - iii {gap:+.3f}: {zone}; AUC iii {r3['auc']:.3f}: {'ok' if r3['auc'] >= 0.90 else 'FAILS'}")
    print(f"{time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
