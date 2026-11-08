"""M5-ACAP, gate G-M5b (no saturation) of the v2 specification (M5; build order step 3, plan T4.13).

Question. On the continuous Ridgecrest development segment (CAL, origins 0-900 s), is the AUC of the ACAP ranking
statistic (relative cap c = 3 x the median block norm of the hypothesis, ttmt.detector.ACAP) on the true hypocentres
against the maxima of windows without an event >= 0.65? The registered cap 6 gave 0.39-0.49 (C18b/C18c,
results/explorations/c18/c18c_*.json).

Design. The oracle of C18c (code/explore_contrib/c18c_oracle.py), unchanged except the statistics: every 3 s window
of the CAL segment, the hard reading of the registered engine (station corrections SC), the score map of each
statistic; null = the maximum over the owned origins (5-8 s) of each window with no catalogued event within 2 s;
event = the maximum within 3 km and 0.5 s of each catalogued M >= 2 origin owned by the window. AUC = P(event >
null) + 0.5 P(tie). Statistics:
  CAP6    registered cap (must reproduce c18c CAP6 exactly: a check of the pipeline)
  NOCAP   unbounded (must reproduce c18c NOCAP exactly)
  ACAP3   the gate statistic (c = 3)
  ACAP2, ACAP5   sensitivity of c (declared, not part of the gate)
Gate (declared before the run): ACAP3 AUC >= 0.65 on each of the four registered subnets N4, N6, N8, U8. Reported with
a bootstrap 95 % interval (events and null maxima resampled independently, B = 2000, seed 20261108), since there are
only 16 events.

    TTMT_DATA_ROOT=... python experiments/exploratory/M5-ACAP/run_b.py [SUBNET ...] > .../run_b_out.txt
Default subnets: N4 (the small first run), then N6 N8 U8 when asked for.
"""
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(ROOT, "code", "r0"))
sys.path.insert(0, os.path.join(ROOT, "src"))
import r0_lib as R  # noqa: E402  (must come first: it sets the registered medium)
import r0_run as RR  # noqa: E402
import c10_lib as C  # noqa: E402
import xc  # noqa: E402

from ttmt.config import dataset_dir  # noqa: E402
from ttmt.detector import ACAP  # noqa: E402

R.DATA = str(dataset_dir("ridgecrest_2019"))
SEED = 20261108
B = 2000
GATE = 0.65
VARIANTS = ("CAP6", "NOCAP", "ACAP3", "ACAP2", "ACAP5")
STORED = os.path.join(ROOT, "results", "explorations", "c18")


def score_map(E, MFn, kind, c):
    tpk, _ = C.choose(E, None, "hard", None, C.SH, R=MFn)
    if kind == "cap":
        vf = C._vectors(MFn, tpk, c).reshape(E.N, E.T, -1)
    else:
        vf = C._vectors(MFn, tpk, None).reshape(E.N, E.T, -1)
    if kind == "acap":
        return ACAP(c).score(vf, E.Q)[0]
    y = np.einsum("gtr,grk->gtk", vf, E.Q)
    return (y * y).sum(-1) / 2.0


def auc(st, nul):
    return float(np.mean([(np.sum(nul < x) + 0.5 * np.sum(nul == x)) / len(nul) for x in st]))


def auc_ci(st, nul, rng):
    bs = [auc(rng.choice(st, len(st)), rng.choice(nul, len(nul))) for _ in range(B)]
    return np.percentile(bs, [2.5, 97.5])


def run(sub):
    t0 = time.time()
    net = RR.Net(sub, 940.0)
    E = net.eng["SC"][0]
    tg = [e for e in net.cat if R.T_CAL[0] <= e["t"] < R.T_CAL[1] and e["M"] >= R.M_TARGET]
    allt = np.array([e["t"] for e in net.cat])
    kinds = {"CAP6": ("cap", 6.0), "NOCAP": ("none", None), "ACAP3": ("acap", 3.0), "ACAP2": ("acap", 2.0),
             "ACAP5": ("acap", 5.0)}
    out = {v: {"null": [], "ev": []} for v in VARIANTS}
    ownt = (E.taus >= R.OWN[0] - 1e-9) & (E.taus < R.OWN[1] - 1e-9)
    for ts in RR.starts("cal"):
        MF, sd = xc.matched(R.window(net.d16, ts))
        MFn = MF / sd
        lo, hi = ts + R.OWN[0], ts + R.OWN[1]
        own = [e for e in tg if lo <= e["t"] < hi]
        quiet = not np.any((allt > lo - 2.0) & (allt < hi + 2.0))
        for v in VARIANTS:
            s = score_map(E, MFn, *kinds[v])
            if quiet:
                out[v]["null"].append(float(np.where(ownt[None, :], s, -np.inf).max()))
            for e in own:
                near = np.hypot(E.XN[:, 0] - e["x"], E.XN[:, 1] - e["y"]) <= 3.0
                tt = np.abs(ts + E.taus - e["t"]) <= 0.5
                out[v]["ev"].append(float(s[np.ix_(near, tt)].max()) if near.any() and tt.any() else np.nan)
    return out, time.time() - t0


def main(subs):
    import json
    print(f"M5-ACAP G-M5b  Ridgecrest CAL oracle (C18c design)  seed {SEED}  B {B}  numpy {np.__version__}")
    print(f"gate: ACAP3 AUC >= {GATE} on every subnet; CAP6 and NOCAP must reproduce c18c exactly\n")
    print(f"{'subnet':6s} {'statistic':9s} {'events':>6s} {'null':>5s} {'AUC':>6s} {'95% CI':>15s} {'c18c':>6s} "
          f"{'>q90':>5s} {'>q99':>5s}")
    verdict = {}
    for sub in subs:
        out, sec = run(sub)
        ref = json.load(open(os.path.join(STORED, f"c18c_{sub}.json")))["summary"]
        rng = np.random.default_rng([SEED, sum(map(ord, sub))])
        for v in VARIANTS:
            st, nul = np.array(out[v]["ev"]), np.array(out[v]["null"])
            st = st[np.isfinite(st)]
            a = auc(st, nul)
            lo, hi = auc_ci(st, nul, rng)
            q90, q99 = np.percentile(nul, [90, 99])
            r = ref.get(v, {}).get("auc")
            rs = f"{r:6.3f}" if r is not None else f"{'-':>6s}"
            print(f"{sub:6s} {v:9s} {len(st):6d} {len(nul):5d} {a:6.3f} [{lo:6.3f},{hi:6.3f}] {rs} "
                  f"{np.mean(st > q90):5.2f} {np.mean(st > q99):5.2f}")
            if r is not None and abs(a - r) > 1e-9:
                print(f"   WARNING: {v} does not reproduce c18c ({a:.6f} vs {r:.6f})")
            if v == "ACAP3":
                verdict[sub] = a >= GATE
        print(f"   {sub}: {sec:.0f} s", flush=True)
    print("\nG-M5b:", "PASS" if all(verdict.values()) else "FAIL", verdict,
          "" if len(verdict) == 4 else "(partial: not all four subnets run)")


if __name__ == "__main__":
    main(sys.argv[1:] or ["N4"])
