"""WHITEN-REVIEW part 2: C18b (Ridgecrest calibration segment only) with Q from diag(1/scale) G. See card.md.

    python run.py check            code checks 2 and 3 (no waveforms beyond the network set-up)
    python run.py SUBNET OUT.json  REG, REG-W, NOCAP, NOCAP-W for one subnet (N4, N6, N8, U8)
    python run.py analyse OUTDIR   check 1 against the frozen C18b JSON, the AUC differences and the bootstrap

The loop is c18b_oracle.main verbatim except for the variants; nothing outside [0, 900) s is scored.
"""
import copy
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
for p in (os.path.join(ROOT, "code"), os.path.join(ROOT, "code", "r0"), os.path.join(ROOT, "code", "explore_contrib")):
    if p not in sys.path:
        sys.path.insert(0, p)

import r0_lib as R  # noqa: E402
import r0_run as RR  # noqa: E402
import c10_lib as C  # noqa: E402
import xc  # noqa: E402
import task1x_mtglrt as X  # noqa: E402
from certified_fw import stack_Q  # noqa: E402
import c18b_oracle as B  # noqa: E402  (score_map, unchanged)

VARIANTS = ("REG", "REG-W", "NOCAP", "NOCAP-W")
SEED = 20261126
NBOOT = 10000
FROZEN = os.path.join(ROOT, "results", "explorations", "c18")


def row_scale(scale):
    """per-row factor of G (rows ordered phase, station, component, as task1x_mtglrt.radiation_operator)."""
    return np.tile(np.repeat(np.asarray(scale, float), 3), 2)


def weighted_Q(Gs, scale):
    w = 1.0 / row_scale(scale)
    return stack_Q([X.basis(w[:, None] * g) for g in Gs])


def engines(net):
    E = net.eng["SC"][0]
    Ew = copy.copy(E)
    Ew.Q = weighted_Q(net.Gs, net.scale)
    return E, Ew


def run_subnet(sub, outp):
    net = RR.Net(sub, 940.0)
    E, Ew = engines(net)
    cat = net.cat
    tg = [e for e in cat if R.T_CAL[0] <= e["t"] < R.T_CAL[1] and e["M"] >= R.M_TARGET]
    allt = np.array([e["t"] for e in cat])
    starts = RR.starts("cal")
    out = {v: {"null": [], "events": []} for v in VARIANTS}
    ownt = (E.taus >= R.OWN[0] - 1e-9) & (E.taus < R.OWN[1] - 1e-9)
    for ts in starts:
        W = R.window(net.d16, ts)
        MF, sd = xc.matched(W)
        inp = {"REG": (E, MF / sd, C.CAPV), "REG-W": (Ew, MF / sd, C.CAPV),
               "NOCAP": (E, MF / sd, None), "NOCAP-W": (Ew, MF / sd, None)}
        lo, hi = ts + R.OWN[0], ts + R.OWN[1]
        own = [e for e in tg if lo <= e["t"] < hi]
        quiet = not np.any((allt > lo - 2.0) & (allt < hi + 2.0))
        for v, (Ev, M_, capv) in inp.items():
            s = B.score_map(Ev, M_, capv, False)
            s_own = np.where(ownt[None, :], s, -np.inf)
            smax = float(s_own.max())
            if quiet:
                out[v]["null"].append(smax)
            for e in own:
                d = np.hypot(E.XN[:, 0] - e["x"], E.XN[:, 1] - e["y"])
                near = d <= 3.0
                tt = np.abs(ts + E.taus - e["t"]) <= 0.5
                st = float(s[np.ix_(near, tt)].max()) if near.any() and tt.any() else float("nan")
                out[v]["events"].append({"t": e["t"], "M": e["M"], "s_true": st, "s_map": smax})
    json.dump({"subnet": sub, "scale": net.scale.tolist(), "detail": out}, open(outp, "w"))
    print(sub, "done", flush=True)


def auc(st, nul):
    st, nul = np.asarray(st), np.sort(np.asarray(nul))
    lo = np.searchsorted(nul, st, "left")
    hi = np.searchsorted(nul, st, "right")
    return float(np.mean((lo + 0.5 * (hi - lo)) / len(nul)))


def check():
    """checks 2 and 3 of the card, on the N8 geometry (set-up reads waveforms only to build the network)."""
    net = RR.Net("N8", 940.0)
    Q1 = weighted_Q(net.Gs, np.ones(len(net.STA)))
    same = np.array_equal(Q1, net.Qst)
    print(f"check 2: scale = 1 gives the registered Q bit for bit: {same}")
    sc = np.asarray(net.scale, float)
    print(f"station_scale max/min {sc.max() / sc.min():.2f}")
    Qw = weighted_Q(net.Gs, sc)
    w = 1.0 / row_scale(sc)
    rng = np.random.default_rng(SEED)
    fw, fu = [], []
    for g in rng.choice(len(net.Gs), 200, replace=False):
        m = rng.normal(size=6)
        v = w * (net.Gs[g] @ m)
        fw.append(float(np.sum((Qw[g].T @ v) ** 2) / np.sum(v ** 2)))
        fu.append(float(np.sum((net.Qst[g].T @ v) ** 2) / np.sum(v ** 2)))
    print(f"check 3: noise-free f, weighted basis min {min(fw):.4f} median {np.median(fw):.4f}; "
          f"unweighted median {np.median(fu):.4f} min {min(fu):.4f}")
    ok = same and min(fw) >= 0.99 and (sc.max() / sc.min() < 3 or np.median(fu) < 0.95)
    print("CHECKS", "PASS" if ok else "FAIL")


def analyse(outdir):
    subs = ("N4", "N6", "N8", "U8")
    D = {s: json.load(open(os.path.join(outdir, f"wr_{s}.json"))) for s in subs}
    print("check 1: REG and NOCAP against the frozen C18b")
    ok1 = True
    for s in subs:
        fz = json.load(open(os.path.join(FROZEN, f"c18b_{s}.json")))
        for v in ("REG", "NOCAP"):
            a_new = auc([e["s_true"] for e in D[s]["detail"][v]["events"]], D[s]["detail"][v]["null"])
            a_old = fz["summary"][v]["auc"]
            nn, no = len(D[s]["detail"][v]["null"]), fz["summary"][v]["n_null"]
            eq = abs(a_new - a_old) < 1e-12 and nn == no
            ok1 &= eq
            print(f"   {s} {v}: AUC {a_new:.4f} frozen {a_old:.4f}, null {nn}/{no} {'same' if eq else 'DIFFERENT'}")
    print("check 1", "PASS" if ok1 else "FAIL")
    rng = np.random.default_rng(SEED)
    ne = len(D["N4"]["detail"]["REG"]["events"])
    nn = len(D["N4"]["detail"]["REG"]["null"])
    S = {(s, v): (np.array([e["s_true"] for e in D[s]["detail"][v]["events"]]), np.array(D[s]["detail"][v]["null"]))
         for s in subs for v in VARIANTS}
    for s in subs:
        sc = np.asarray(D[s]["scale"])
        print(f"{s}: station_scale max/min {sc.max() / sc.min():.2f}")
    for base, w in (("NOCAP", "NOCAP-W"), ("REG", "REG-W")):
        print(f"\n{w} - {base}")
        diffs = []
        for s in subs:
            a0, a1 = auc(*S[(s, base)]), auc(*S[(s, w)])
            q99 = np.percentile(S[(s, w)][1], 99), np.percentile(S[(s, base)][1], 99)
            f99 = np.mean(S[(s, w)][0] > q99[0]), np.mean(S[(s, base)][0] > q99[1])
            med = np.median(S[(s, w)][0] / np.percentile(S[(s, w)][1], 90)), \
                np.median(S[(s, base)][0] / np.percentile(S[(s, base)][1], 90))
            diffs.append(a1 - a0)
            print(f"   {s}: AUC {a0:.3f} -> {a1:.3f} ({a1 - a0:+.3f}); above q99 {f99[1]:.3f} -> {f99[0]:.3f}; "
                  f"median s_true/q90 {med[1]:.3f} -> {med[0]:.3f}")
        bs = np.empty((NBOOT, len(subs)))
        for b in range(NBOOT):
            ie, inn = rng.integers(0, ne, ne), rng.integers(0, nn, nn)
            for k, s in enumerate(subs):
                bs[b, k] = auc(S[(s, w)][0][ie], S[(s, w)][1][inn]) - auc(S[(s, base)][0][ie], S[(s, base)][1][inn])
        pooled = bs.mean(1)
        lo, hi = np.percentile(pooled, [2.5, 97.5])
        per = np.percentile(bs, [2.5, 97.5], axis=0)
        print("   per subnet 95 %: " + ", ".join(f"{s} [{per[0, k]:+.3f}, {per[1, k]:+.3f}]" for k, s in enumerate(subs)))
        m = float(np.mean(diffs))
        zone = ("defect explains part of the failure" if m >= 0.05 and lo > 0 else
                "defect does not explain the failure" if -0.05 < m < 0.05 else
                "unexpected, diagnose" if m <= -0.05 else "between zones (>= +0.05 but interval reaches 0)")
        print(f"   pooled mean {m:+.3f}, 95 % [{lo:+.3f}, {hi:+.3f}] (n events {ne}, null windows {nn}, "
              f"{NBOOT} resamples, seed {SEED}) -> {zone}")


if __name__ == "__main__":
    if sys.argv[1] == "check":
        check()
    elif sys.argv[1] == "analyse":
        analyse(sys.argv[2])
    else:
        run_subnet(sys.argv[1], sys.argv[2])
