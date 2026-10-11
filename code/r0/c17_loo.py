"""C17a (exploratory, CAL only): leave-one-event-out prediction of the timing residuals after the homogeneous model by
(a) a constant term per station-phase (median of the other events), (b) Gaussian-process kriging over source position
with an exponential covariance (nugget, sill, range fitted by least squares to the empirical semivariogram of the other
events), (c) a per station-phase linear model in source depth and epicentral distance (Huber-free least squares).
Reads: station_corrections.peak_delays at the catalogue hypocentre (+-1 s), energy >= 25."""
import sys
import numpy as np
import r0_run as RR
from station_corrections import peak_delays

sub = sys.argv[1] if len(sys.argv) > 1 else "N8"
net = RR.Net(sub, 940.0)
rows = []
for e, (D, src, tau) in zip([e for e in net.calev if e["t"] >= 6.5], net.calwin):
    dl, en = peak_delays(net.STA, D, src, tau)
    for i in range(len(net.STA)):
        for b in range(2):
            if en[i, b] >= 25:
                rows.append((e["t"], src[0], src[1], src[2], i, b, dl[i, b], np.linalg.norm(net.STA[i] - src)))
A = np.array(rows)
ev = np.unique(A[:, 0])


def gp_pred(Xtr, ytr, Xte, nug, sill, rng):
    d = np.linalg.norm(Xtr[:, None] - Xtr[None], axis=-1)
    K = sill * np.exp(-d / rng) + nug * np.eye(len(Xtr))
    k = sill * np.exp(-np.linalg.norm(Xtr - Xte[None], axis=-1) / rng)
    mu = np.median(ytr)
    return mu + k @ np.linalg.solve(K, ytr - mu)


errs = {"const": [], "gp": [], "lin": []}
for e0 in ev:
    tr, te = A[A[:, 0] != e0], A[A[:, 0] == e0]
    # pooled semivariogram of the training events (same station-phase pairs), exponential fit by grid search
    dd, g = [], []
    for i in range(len(tr)):
        for j in range(i + 1, len(tr)):
            if tr[i, 4] == tr[j, 4] and tr[i, 5] == tr[j, 5]:
                dd.append(np.linalg.norm(tr[i, 1:4] - tr[j, 1:4]))
                g.append(0.5 * (tr[i, 6] - tr[j, 6]) ** 2)
    dd, g = np.array(dd), np.array(g)
    best = None
    for nug in np.linspace(0.0, 0.3, 13):
        for sill in np.linspace(0.0, 0.4, 17):
            for rg in (2.0, 4.0, 8.0, 16.0, 32.0):
                f = nug + sill * (1 - np.exp(-dd / rg))
                L = np.median(np.abs(f - g))
                if best is None or L < best[0]:
                    best = (L, nug, sill, rg)
    _L, nug, sill, rg = best
    for r in te:
        m = (tr[:, 4] == r[4]) & (tr[:, 5] == r[5])
        if m.sum() < 3:
            continue
        y = tr[m, 6]
        errs["const"].append(r[6] - np.median(y))
        errs["gp"].append(r[6] - (gp_pred(tr[m, 1:4], y, r[1:4], max(nug, 1e-3), sill, rg) if sill > 0 else np.median(y)))
        X = np.c_[np.ones(m.sum()), tr[m, 3], tr[m, 7]]
        if m.sum() >= 5:
            c, *_ = np.linalg.lstsq(X, y, rcond=None)
            errs["lin"].append(r[6] - np.array([1.0, r[3], r[7]]) @ c)
        else:
            errs["lin"].append(r[6] - np.median(y))
print(sub, "reads", len(A), "events", len(ev), "last fit nugget/sill/range", round(nug, 3), round(sill, 3), rg)
for k, v in errs.items():
    v = np.array(v)
    print(f"  {k:5s} n {len(v)}  robust sd {1.4826 * np.median(np.abs(v - np.median(v))):.3f}  rms {np.sqrt(np.mean(v ** 2)):.3f}"
          f"  |e|<=0.25: {np.mean(np.abs(v) <= 0.25):.2f}")
