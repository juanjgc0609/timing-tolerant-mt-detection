"""load_verif.py - independent loader / hit labelling / FROC for the 2L parts (no code of analyze_2l.py or c10_lib)."""
import json, os, pickle
import numpy as np
PARTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "task2l_parts")
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
REGIMES = ("thin", "strong")
NNET = {"std": 160, "dense": 100}
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache_verif.pkl")

def load_all():
    if os.path.exists(CACHE):
        return pickle.load(open(CACHE, "rb"))
    W, N = [], {}
    for seed in (7171, 9393):
        for reg in REGIMES:
            for cell in CELLS:
                for j in range(160):
                    e = json.load(open(os.path.join(PARTS, f"main_seed{seed}_{reg}_{cell}_net{j}.json")))
                    assert (e["seed"], e["regime"], e["cell"], e["net"]) == (seed, reg, cell, j)
                    N[(seed, reg, cell, j)] = {k: e[k] for k in ("err", "off", "sigT3", "S2", "terms", "statics", "tags", "kmax", "cat_drawn", "methods", "host_time_utc")}
                    for kind in ("std", "dense"):
                        if j >= NNET[kind]:
                            assert kind not in e["windows"]
                            continue
                        rows = e["windows"][kind]
                        assert len(rows) == 5
                        for w in rows:
                            W.append({"seed": seed, "reg": reg, "cell": cell, "net": j, "kind": kind, "w": w["w"], "k": w["k"],
                                      "truth": np.array(w["truth"], float).reshape(-1, 3), "taus": np.array(w["taus"], float),
                                      "atoms": {m: np.array(a, float).reshape(-1, 7) for m, a in w["atoms"].items()}})
    pickle.dump((W, N), open(CACHE, "wb"))
    return W, N

def label(at, truth, taus, tol_t, tol_km=5.0, shift=0.0):
    """own labelling.  at (n, 7): x y z tau s s0 sref.  tol_t None: nearest epicentre within tol_km.
    else nearest event among those with epicentre <= tol_km and |tau - shift - true tau| <= tol_t."""
    lab = np.full(len(at), -1, int)
    if len(truth) == 0:
        return lab
    for i, a in enumerate(at):
        d = np.hypot(truth[:, 0] - a[0], truth[:, 1] - a[1])
        ok = d <= tol_km
        if tol_t is not None:
            ok &= np.abs(a[3] - shift - taus) <= tol_t + 1e-9
        if ok.any():
            lab[i] = int(np.argmin(np.where(ok, d, np.inf)))
    return lab

def window_rec(at, lab, running=True, key=4):
    """per atom: (score, is_false, is_new_hit) in selection order; score = running minimum of column key."""
    s = at[:, key].copy()
    if running:
        s = np.minimum.accumulate(s) if len(s) else s
    seen, new = set(), np.zeros(len(at), bool)
    for i, h in enumerate(lab):
        if h >= 0 and h not in seen:
            seen.add(h); new[i] = True
    return s, lab < 0, new

def froc(recs, mult, nte, caps=(0.1, 0.25)):
    """exact FROC.  recs: list of (s, false, new) per window; mult: multiplicity per window (bootstrap weights);
    nte: number of true events (weighted).  Recall at realised false/window <= cap, threshold = any score."""
    idx = [i for i in range(len(recs)) if mult[i] > 0 and len(recs[i][0])]
    nw = float(np.sum(mult))
    if not idx or nte <= 0:
        return [0.0 for _ in caps], [np.inf for _ in caps]
    s = np.concatenate([recs[i][0] for i in idx])
    f = np.concatenate([recs[i][1] * mult[i] for i in idx]).astype(float)
    h = np.concatenate([recs[i][2] * mult[i] for i in idx]).astype(float)
    o = np.argsort(-s, kind="stable")
    s, F, H = s[o], np.cumsum(f[o]), np.cumsum(h[o])
    last = np.r_[np.nonzero(np.diff(s) != 0)[0], len(s) - 1]
    F, H, th = F[last], H[last], s[last]
    out, thr = [], []
    for c in caps:
        ok = F / nw <= c + 1e-12
        if ok.any():
            j = int(np.argmax(np.where(ok, H, -1)))
            out.append(float(H[j] / nte)); thr.append(float(th[j]))
        else:
            out.append(0.0); thr.append(np.inf)
    return out, thr
