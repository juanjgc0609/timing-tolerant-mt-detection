"""
timing_cal.py - Task R1, post hoc and descriptive (CAL only): timing residuals of the matched-filter energy-peak
readings of the CAL truth events (magnitude >= 1.0) against the JMA2001 times plus the frozen run's matched-filter
terms, read as in r1_run.mf_terms (search +-min(2 s, half the predicted S-P), energy gate 25).  Leave-one-event-out
terms are also reported, so that the residual scale is not in-sample.
    python timing_cal.py PICKS_ALL.json OUT.txt
"""
import json
import sys

import numpy as np

sys.path.insert(0, "/home/claude/work/r1")
import r1_lib as R  # noqa: E402

R.setup()
import c18_tt as TT  # noqa: E402
import c18j_layered as LJ  # noqa: E402
from station_corrections import peak_delays  # noqa: E402
from task0_common import VP, VS  # noqa: E402


def main(pkp, outp):
    net = R.make_net("cal")
    pt = np.array(json.load(open(pkp))["pick_terms"])
    ds, ens = [], []
    for D, src, tau in net.calwin:
        src = np.asarray(src, float)
        tl = TT.times(LJ.MODEL, LJ.R_VPVS, src, net.STA)
        d = np.linalg.norm(net.STA - src, axis=1)
        th = np.stack([d / VP, d / VS], axis=1)
        dl, en = np.zeros((len(d), 2)), np.zeros((len(d), 2))
        for i in range(len(d)):
            a, b = peak_delays(net.STA, D, src, float(tau), min(2.0, 0.5 * float(tl[1, i] - tl[0, i])))
            dl[i], en[i] = a[i], b[i]
        ds.append(dl + th - tl.T)
        ens.append(en)
    ds, ens = np.array(ds), np.array(ens)          # (events, stations, phases)
    ok = ens >= 25.0
    lines = [f"R1 CAL, {ds.shape[0]} calibration events (truth M >= 1.0), {ds.shape[1]} stations; readings above the gate: "
             f"P {int(ok[:, :, 0].sum())}, S {int(ok[:, :, 1].sum())}"]
    for tag in ("in-sample", "leave-one-event-out"):
        res = {0: [], 1: []}
        for e in range(ds.shape[0]):
            for b in range(2):
                for i in range(ds.shape[1]):
                    if not ok[e, i, b]:
                        continue
                    if tag == "in-sample":
                        term = pt[i, b]
                        oth = [ds[k, i, b] - pt[i, b] for k in range(ds.shape[0]) if ok[k, i, b]]
                    else:
                        oth = [ds[k, i, b] - pt[i, b] for k in range(ds.shape[0]) if k != e and ok[k, i, b]]
                        term = pt[i, b]
                    lag = np.median(oth) if oth else 0.0
                    res[b].append(ds[e, i, b] - (term + lag))
        for b in range(2):
            r = np.array(res[b])
            mad = 1.4826 * np.median(np.abs(r - np.median(r)))
            lines.append(f"   {tag:20s} {'PS'[b]}: n {len(r)}, median {np.median(r):+.3f} s, robust sd {mad:.3f} s, "
                         f"within 0.25 s {np.mean(np.abs(r) <= 0.25):.2f}, within 0.5 s {np.mean(np.abs(r) <= 0.5):.2f}")
    open(outp, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
