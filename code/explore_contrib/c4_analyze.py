"""C4 analysis. Calibration networks CAL (per cell), held-out networks HO. For each method and alpha:
  P  pooled event-aware: threshold from the CAL networks' standard windows of the cell (largest recall with false
     events per window <= alpha); applied to HO networks.
  E  empty-window split conformal (2D-c): threshold = ceil((n+1)(1-alpha))-th smallest per-window max score of the CAL
     networks' empty windows of the cell.
  N  per-network: threshold from the HO network's OWN 20 catalogue records + its empty windows (false events per
     window <= alpha, plug-in), applied to that network's standard windows.
Realised on HO standard windows: recall, false events per window (all windows; windows with events), fraction of HO
networks whose realised rate exceeds 2 alpha. Bootstrap over HO networks (within cell) and windows."""
import json, sys
import numpy as np

ALPHAS = (0.1, 0.25)


def thr_fw(recs, alpha):
    """lowest threshold t such that false events per window <= alpha (records prefix-monotone)."""
    sc = sorted({s for r in recs for s, _h in r}, reverse=True)
    best = np.inf
    n = len(recs)
    allf = [(s, h) for r in recs for s, h in r]
    arr = np.array([s for s, h in allf if h < 0])
    for t in sc:
        if np.sum(arr >= t) / n <= alpha:
            best = t
        else:
            break
    return best


def thr_empty(recs, alpha):
    mx = np.sort([max([s for s, _h in r] or [-np.inf]) for r in recs])
    n = len(mx)
    k = int(np.ceil((n + 1) * (1 - alpha)))
    return float(mx[min(k, n) - 1]) if k <= n else np.inf


def apply(recs, ks, t):
    hits, fal, fal_ev, n_ev = 0, 0, 0, 0
    for r, k in zip(recs, ks):
        hs = {h for s, h in r if s >= t and h >= 0}
        f = sum(1 for s, h in r if s >= t and h < 0)
        hits += len(hs)
        fal += f
        if k:
            fal_ev += f
            n_ev += 1
    return hits, fal, fal_ev, n_ev


def main(path, cal, ho, methods, B=1000):
    rows = json.load(open(path))
    cells = sorted({r["cell"] for r in rows})
    out = []
    rng = np.random.default_rng(99)
    for alpha in ALPHAS:
        out.append(f"=== alpha {alpha} (false events per window) ===")
        for m in methods:
            res = {}
            for strat in ("P", "E", "N"):
                per_net = []           # (cell, net, hits, ntrue, fal, nwin, fal_ev, n_ev)
                for cell in cells:
                    R = [r for r in rows if r["cell"] == cell]
                    if strat == "P":
                        cr = [[tuple(x) for x in r["rec"][m]] for r in R if r["net"] in cal and r["type"] == "eval"]
                        t_cell = thr_fw(cr, alpha)
                    elif strat == "E":
                        cr = [[tuple(x) for x in r["rec"][m]] for r in R if r["net"] in cal and r["type"] == "empty"]
                        t_cell = thr_empty(cr, alpha)
                    for j in ho:
                        ev = [r for r in R if r["net"] == j and r["type"] == "eval"]
                        if not ev:
                            continue
                        if strat == "N":
                            cr = [[tuple(x) for x in r["rec"][m]] for r in R if r["net"] == j and r["type"] in ("cat", "empty")]
                            t = thr_fw(cr, alpha)
                        else:
                            t = t_cell
                        h, f, fe, ne = apply([[tuple(x) for x in r["rec"][m]] for r in ev], [r["k"] for r in ev], t)
                        per_net.append((cell, j, h, sum(r["k"] for r in ev), f, len(ev), fe, ne))
                A = np.array([p[2:] for p in per_net], float)
                rec = A[:, 0].sum() / A[:, 1].sum()
                fr = A[:, 2].sum() / A[:, 3].sum()
                fre = A[:, 4].sum() / max(1, A[:, 5].sum())
                rate_net = A[:, 2] / A[:, 3]
                over = float(np.mean(rate_net > 2 * alpha))
                # bootstrap over networks within cell
                cl = np.array([p[0] for p in per_net])
                bs = []
                for _ in range(B):
                    idx = np.concatenate([rng.choice(np.where(cl == c)[0], (cl == c).sum()) for c in cells])
                    Ab = A[idx]
                    bs.append([Ab[:, 0].sum() / Ab[:, 1].sum(), Ab[:, 2].sum() / Ab[:, 3].sum()])
                bs = np.array(bs)
                lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
                res[strat] = (rec, fr)
                out.append(f"  {m:14s} {strat}: recall {rec:.3f} [{lo[0]:.3f},{hi[0]:.3f}]  false/window {fr:.3f} "
                           f"[{lo[1]:.3f},{hi[1]:.3f}]  (windows with events {fre:.3f})  nets > 2alpha {over:.2f}  "
                           f"per-net rate median {np.median(rate_net):.2f} max {rate_net.max():.2f}")
    return "\n".join(out)


if __name__ == "__main__":
    path = sys.argv[1]
    cal = [int(x) for x in sys.argv[2].split(",")]
    ho = [int(x) for x in sys.argv[3].split(",")]
    methods = sys.argv[4].split(",")
    print(main(path, cal, ho, methods))
