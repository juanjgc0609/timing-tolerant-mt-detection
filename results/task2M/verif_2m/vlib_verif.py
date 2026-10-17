"""Independent FROC / recall / bootstrap implementation for the 2M verification (no project imports)."""
import glob
import json
import os

import numpy as np

W = "/home/claude/work"
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")


def load_parts(seed, arm="REAL"):
    pdir = os.path.join(W, "t2m_out", f"parts_{seed}")
    nets = []
    for p in sorted(glob.glob(os.path.join(pdir, f"{arm}_*.json"))):
        with open(p) as fh:
            nets.append(json.load(fh))
    return nets


def windows(nets, design):
    out = []
    for n in nets:
        for r in n["rows"]:
            if r["design"] == design:
                out.append({"cell": n["cell"], "net": n["net"], "day": n["day"], "w": r["w"], "k": r["k"],
                            "truth": np.array(r["truth"], float).reshape(-1, 3), "taus": np.array(r["taus"], float),
                            "atoms": r["atoms"], "slot": r["slot"], "kurt": r["kurt"], "maxabs": r["maxabs"]})
    return out


def label_window(atoms, truth, taus, tol_km=5.0, tol_t=None):
    """per atom: index of the truth it hits (nearest qualifying epicentre), -1 if false; effective score =
    running minimum of the scores in emission order."""
    a = np.array(atoms, float).reshape(-1, 7)
    n = len(a)
    lab = np.full(n, -1, int)
    if n and len(truth):
        dd = np.sqrt(((a[:, None, :2] - truth[None, :, :2]) ** 2).sum(-1))
        if tol_t is None:
            j = dd.argmin(1)
            ok = dd[np.arange(n), j] <= tol_km
            lab = np.where(ok, j, -1)
        else:
            q = (dd <= tol_km) & (np.abs(a[:, 3][:, None] - taus[None, :]) <= tol_t + 1e-9)
            dq = np.where(q, dd, np.inf)
            j = dq.argmin(1)
            lab = np.where(q.any(1), j, -1)
    s = np.minimum.accumulate(a[:, 4]) if n else np.zeros(0)
    return s, lab


def prep(wins, method, tol_km=5.0, tol_t=None, kcut=None):
    S, Lb = [], []
    for w in wins:
        at = w["atoms"][method][:kcut] if kcut else w["atoms"][method]
        s, lab = label_window(at, w["truth"], w["taus"], tol_km, tol_t)
        S.append(s)
        Lb.append(lab)
    K = np.array([w["k"] for w in wins])
    return S, Lb, K


def recall_at(S, Lb, K, ws, caps=(0.1, 0.25)):
    """exact FROC over the multiset of windows ws: for each score level t, H(t) = distinct (window draw, truth) hits
    among atoms with score >= t, F(t) = false atoms with score >= t; recall(cap) = max H/E over levels with F/nwin <= cap."""
    nw = len(ws)
    E = K[ws].sum()
    sc = np.concatenate([S[w] for w in ws])
    lab = np.concatenate([Lb[w] for w in ws])
    pos = np.repeat(np.arange(nw), [len(S[w]) for w in ws])
    if E == 0 or not len(sc):
        return np.zeros(len(caps))
    o = np.lexsort((np.arange(len(sc)), -sc))
    sc, lab, pos = sc[o], lab[o], pos[o]
    key = pos * 1000 + lab
    isfalse = lab < 0
    # first occurrence of each hit key in descending-score order
    seen = np.zeros(len(sc), bool)
    hk = key[~isfalse]
    idx = np.nonzero(~isfalse)[0]
    _, first = np.unique(hk, return_index=True)
    seen[idx[first]] = True
    Hc = np.cumsum(seen)
    Fc = np.cumsum(isfalse)
    # evaluate at the end of each block of equal scores
    ends = np.nonzero(np.r_[sc[1:] != sc[:-1], True])[0]
    Hc, Fc = Hc[ends], Fc[ends]
    out = []
    for c in caps:
        ok = Fc <= c * nw + 1e-12
        out.append(Hc[ok].max() / E if ok.any() else 0.0)
    return np.array(out)


def boot_index(wins, B, seed, cells=CELLS):
    """hierarchical: strata = cell; networks with replacement within a cell, then windows with replacement within
    each drawn network."""
    rng = np.random.default_rng(seed)
    groups = {}
    for i, w in enumerate(wins):
        groups.setdefault((w["cell"], w["net"]), []).append(i)
    bycell = {c: [np.array(v) for k, v in sorted(groups.items()) if k[0] == c] for c in cells}
    bycell = {c: v for c, v in bycell.items() if v}
    reps = []
    for _ in range(B):
        ws = []
        for c, gl in bycell.items():
            for gi in rng.integers(0, len(gl), len(gl)):
                g = gl[gi]
                ws.append(g[rng.integers(0, len(g), len(g))])
        reps.append(np.concatenate(ws))
    return reps
