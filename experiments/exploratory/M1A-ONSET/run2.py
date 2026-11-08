"""M1A-ONSET, second step: diagnosis of the G-M1a failure on CAL (0.87 < 0.90 for P) and confirmation on EVAL.

Declared before running (card, "Second step"). The registered first result stays as it is (run_out.txt: fails).
1. Diagnosis on CAL: every gated P or S block whose onset residual exceeds 0.25 s is listed with its q value, the
   number of onset maxima in the window, and the residual of the maximum nearest to the prediction.
2. One variant, fixed in advance (no parameter is tuned on these events): "prior" reads the maximum that maximises
   log q - (s - mu)^2 / (2 sigma^2), the per-block local term of M4, with mu = prediction + term and sigma the robust
   sd of the phase. On CAL both come leave-one-event-out from the largest-q readings of the other events.
3. Confirmation on EVAL (840-1960 s, the 72-event part never used to build anything): truth events with M >= 1.0,
   terms and sigma from all CAL events (held out by construction). The gate is read on EVAL, for both readers.
The PhaseNet clone is checked against the commit pinned in code/real_picker.py before any inference.

    python experiments/exploratory/M1A-ONSET/run2.py > experiments/exploratory/M1A-ONSET/run2_out.txt
"""
import subprocess
import sys
from pathlib import Path

import numpy as np

from ttmt import legacy
from ttmt.datasets import HINET as R
from ttmt.r1 import Net, model_times
from ttmt.station_terms import peak_delays
from ttmt.v2.m1 import FLOOR, local_maxima, parabolic_offset

HERE = Path(__file__).resolve().parent
FS = 100.0
GATE = 25.0


def check_phasenet():
    legacy.install()
    import real_picker as RP
    head = subprocess.run(["git", "-C", RP.PN_REPO, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    if head != RP.PN_COMMIT:
        raise SystemExit(f"PhaseNet clone at {head or '?'}, pinned {RP.PN_COMMIT}: run "
                         f"git -C {RP.PN_REPO} checkout {RP.PN_COMMIT} first")
    return head


def phasenet_q(net, name):
    cache = HERE / f"q_{name}.npz"
    if cache.exists():
        return np.load(cache)["P"]
    import r0_picks
    P = r0_picks.phasenet_hour(net.raw).astype(np.float32)
    np.savez_compressed(cache, P=P)
    return P


def maxima(qb, t_pred, half):
    """sub-sample times (s) and q of the local maxima of qb within t_pred +- half."""
    lo, hi = max(1, int(np.floor((t_pred - half) * FS))), min(len(qb) - 2, int(np.ceil((t_pred + half) * FS)))
    if hi <= lo:
        return np.zeros(0), np.zeros(0)
    k = np.arange(lo, hi + 1)
    k = k[local_maxima(qb[None, lo - 1:hi + 2])[0][1:-1]]
    d = parabolic_offset(qb[k - 1], qb[k], qb[k + 1])
    return (k + d) / FS, qb[k]


def blocks(net, P, events):
    """per event, station, phase: prediction (record s), half window, maxima (times, q), energy-peak delay, energy."""
    out = []
    for ev in events:
        if ev["t"] - 6.5 < 0:
            continue
        D = R.window(net.d16, ev["t"] - 6.5)
        src = np.array([ev["x"], ev["y"], ev["z"]], float)
        tl = model_times(src, net.STA)
        d = np.linalg.norm(net.STA - src, axis=1)
        th = np.stack([d / R.vp, d / R.vs], axis=1)
        rows = []
        for i in range(len(net.STA)):
            half = min(2.0, 0.5 * float(tl[1, i] - tl[0, i]))
            a, en = peak_delays(net.STA, D, src, 6.5, R.vp, R.vs, half)
            r = []
            for b in range(2):
                tp = ev["t"] + float(tl[b, i])
                ts, qs = maxima(P[i, :, 1 + b], tp, half)
                r.append({"tp": tp, "ts": ts, "qs": qs, "peak": float(a[i, b] + th[i, b] - tl[b, i]),
                          "en": float(en[i, b])})
            rows.append(r)
        out.append((ev, rows))
    return out


def read_largest(blk, mu=None, sig=None):
    if not len(blk["qs"]):
        return np.nan
    return float(blk["ts"][np.argmax(blk["qs"])] - blk["tp"])


def read_prior(blk, mu, sig):
    if not len(blk["qs"]):
        return np.nan
    d = blk["ts"] - blk["tp"]
    return float(d[np.argmax(np.log(blk["qs"]) - (d - mu) ** 2 / (2 * sig ** 2))])


def rsd(x):
    x = np.asarray(x, float)
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def line(name, res):
    r = np.array(res)
    return (f"   {name:44s} n {len(r):3d}, robust sd {rsd(r):.3f} s, within 0.25 s {np.mean(np.abs(r) <= 0.25):.2f}, "
            f"within 0.5 s {np.mean(np.abs(r) <= 0.5):.2f}")


def main():
    head = check_phasenet()
    cal_net, eval_net = Net("cal"), Net("eval")
    cal = blocks(cal_net, phasenet_q(cal_net, "cal"), cal_net.calev)
    lo, hi = R.segment("eval")
    ev_eval = [e for e in eval_net.cat if lo <= e["t"] < hi and e["M"] >= R.m_terms]
    evl = blocks(eval_net, phasenet_q(eval_net, "eval"), ev_eval)
    ns = len(cal_net.STA)
    L = np.array([[[read_largest(r[i][b]) for b in range(2)] for i in range(ns)] for _, r in cal])   # (e, i, b)
    gate = np.array([[[r[i][b]["en"] >= GATE for b in range(2)] for i in range(ns)] for _, r in cal])
    out = [f"M1A-ONSET step 2  PhaseNet commit {head[:12]}  CAL events {len(cal)}, EVAL events {len(evl)} "
           f"(truth M >= {R.m_terms})", ""]

    # 1. diagnosis on CAL (largest-q reader, leave-one-event-out terms, gated blocks)
    out.append("1. CAL failures of the largest-q reader (gated blocks, |residual| > 0.25 s):")
    out.append("   event    M  station phase  residual  q_read  n_max  residual_nearest_max  q_nearest")
    res_cal = {"largest": [[], []], "prior": [[], []]}
    for e, (ev, rows) in enumerate(cal):
        for i in range(ns):
            for b in range(2):
                if not gate[e, i, b] or not np.isfinite(L[e, i, b]):
                    continue
                oth = [L[k, i, b] for k in range(len(cal)) if k != e and gate[k, i, b] and np.isfinite(L[k, i, b])]
                term = np.median(oth) if oth else 0.0
                allr = [L[k, j, b] - np.median([L[kk, j, b] for kk in range(len(cal)) if kk != k and gate[kk, j, b]
                                                and np.isfinite(L[kk, j, b])] or [0.0])
                        for k in range(len(cal)) if k != e for j in range(ns)
                        if gate[k, j, b] and np.isfinite(L[k, j, b])]
                sig = max(rsd(allr), 0.02)
                res_l = L[e, i, b] - term
                res_p = read_prior(rows[i][b], term, sig) - term
                res_cal["largest"][b].append(res_l)
                res_cal["prior"][b].append(res_p)
                if abs(res_l) > 0.25:
                    blk = rows[i][b]
                    d = blk["ts"] - blk["tp"] - term
                    j = int(np.argmin(np.abs(d)))
                    out.append(f"   {ev['t']:7.1f} {ev['M']:4.1f}  {cal_net.names[i]:6s}   {'PS'[b]}   {res_l:+7.3f}   "
                               f"{blk['qs'].max():5.2f}  {len(d):5d}  {d[j]:+20.3f}  {blk['qs'][j]:9.2f}")
    out += ["", "2. CAL, leave-one-event-out, gated blocks with an onset maximum:"]
    for name in ("largest", "prior"):
        for b in range(2):
            out.append(line(f"{name:8s} {'PS'[b]}", res_cal[name][b]))

    # 3. confirmation on EVAL with CAL terms and sigma
    terms = np.zeros((ns, 2))
    sig = np.zeros(2)
    for b in range(2):
        for i in range(ns):
            v = [L[e, i, b] for e in range(len(cal)) if gate[e, i, b] and np.isfinite(L[e, i, b])]
            terms[i, b] = np.median(v) if v else 0.0
        sig[b] = max(rsd([L[e, i, b] - terms[i, b] for e in range(len(cal)) for i in range(ns)
                          if gate[e, i, b] and np.isfinite(L[e, i, b])]), 0.02)
    out += ["", f"3. EVAL (confirmation), CAL terms, sigma P {sig[0]:.3f} s, S {sig[1]:.3f} s:"]
    for gated in (True, False):
        for name, fn in (("largest", read_largest), ("prior", read_prior)):
            for b in range(2):
                res = []
                for ev, rows in evl:
                    for i in range(ns):
                        blk = rows[i][b]
                        if gated and blk["en"] < GATE:
                            continue
                        s = fn(blk, terms[i, b], sig[b])
                        if np.isfinite(s):
                            res.append(s - terms[i, b])
                tag = "gated" if gated else "every block"
                out.append(line(f"{name:8s} {'PS'[b]} ({tag})", res) if res else f"   {name} {'PS'[b]}: none")
    # the registered energy-peak reading on the same EVAL blocks (terms: median CAL peak delay)
    for b in range(2):
        pt = [np.median([r[i][b]["peak"] for _, r in cal if r[i][b]["en"] >= GATE] or [0.0]) for i in range(ns)]
        res = [rows[i][b]["peak"] - pt[i] for _, rows in evl for i in range(ns) if rows[i][b]["en"] >= GATE]
        out.append(line(f"energy peak (registered) {'PS'[b]} (gated)", res) if res else "   energy peak: none")
    out.append("")
    out.append("G-M1a on EVAL: line 'largest P (gated)' (the reader as registered) within 0.25 s >= 0.90; the 'prior' "
               "line is the declared variant.")
    print("\n".join(out))


if __name__ == "__main__":
    sys.exit(main())
