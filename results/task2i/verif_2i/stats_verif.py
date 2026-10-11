"""stats_verif.py - independent re-derivation of the 2I rule estimates from task2i_parts (own loader, own FROC, own
hierarchical bootstrap with its own RNG stream), plus exploratory robustness.  Reads only; writes nothing in 2I files."""
import json, os, sys
import numpy as np
W = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PARTS = os.path.join(W, "task2i_parts")
CELLS = ("ladder_4", "ladder_6", "unilateral_c40"); REGS = ("thin", "strong")
METH = ("BOMP-PK+SC", "MT-FW+SC", "MT-FW-NOCAP+SC", "HUBER-N+SC")
H, BO, MT, NC = "HUBER-N+SC", "BOMP-PK+SC", "MT-FW+SC", "MT-FW-NOCAP+SC"
B = int(os.environ.get("VB", 600))

def load(seed):
    """list of windows: dict(stratum, net, kind, k, truth, taus, det={m: [(score, hit, x, y, z, tau)]}, sec)"""
    wins = []; miss = []; meta = set()
    for ri, reg in enumerate(REGS):
        for ci, cell in enumerate(CELLS):
            for j in range(120):
                p = os.path.join(PARTS, f"prim_seed{seed}_{reg}_{cell}_net{j}.json")
                if not os.path.exists(p): miss.append(p); continue
                e = json.load(open(p))
                assert e["seed"] == seed and e["regime"] == reg and e["cell"] == cell and e["net"] == j and e["arm"] == "prim"
                meta.add((e["tags"]["dense"], e["tags"]["std"], e["site_sd"], e["kmax"]["dense"], e["kmax"]["std"], e["n_win"],
                          tuple(e["methods"]), e["kernel"]["OPENBLAS_CORETYPE"], e["world"]))
                for kind in ("dense", "std"):
                    assert len(e["windows"][kind]) == 5
                    for w in e["windows"][kind]:
                        assert w["k"] == len(w["truth"]) == len(w["taus"])
                        det = {}
                        for m in METH:
                            assert len(w["records"][m]) == len(w["loc"][m]) <= e["kmax"][kind]
                            det[m] = [(r[0], r[1]) + tuple(l) for r, l in zip(w["records"][m], w["loc"][m])]
                        wins.append(dict(st=ri * 3 + ci, reg=reg, cell=cell, net=j, kind=kind, k=w["k"], truth=w["truth"],
                                         taus=w["taus"], det=det, sec=w["sec"]))
    return wins, miss, meta

def relabel(w, m, tol_km=5.0, tol_t=None):
    tr = np.asarray(w["truth"], float).reshape(-1, 3); ta = np.asarray(w["taus"], float); out = []
    for d in w["det"][m]:
        h = -1
        if len(tr):
            dd = np.hypot(tr[:, 0] - d[2], tr[:, 1] - d[3]); ok = dd <= tol_km
            if tol_t is not None: ok &= np.abs(ta - d[5]) <= tol_t + 1e-12
            if ok.any(): h = int(np.argmin(np.where(ok, dd, np.inf)))
        out.append(h)
    return out

class Froc:
    """pooled FROC of one method over a list of windows; labels optional override."""
    def __init__(self, wins, m, labels=None, score=None):
        sc, wi, hit = [], [], []
        for i, w in enumerate(wins):
            lab = [d[1] for d in w["det"][m]] if labels is None else labels[i]
            for a, (d, h) in enumerate(zip(w["det"][m], lab)):
                sc.append(d[0] if score is None else score(w, m, a)); wi.append(i); hit.append(h)
        sc, wi, hit = np.array(sc), np.array(wi, int), np.array(hit, int)
        o = np.lexsort((np.arange(len(sc)), -sc)); sc, wi, hit = sc[o], wi[o], hit[o]
        seen = set(); new = np.zeros(len(sc))
        for i in range(len(sc)):
            if hit[i] >= 0 and (wi[i], hit[i]) not in seen: seen.add((wi[i], hit[i])); new[i] = 1
        self.wi, self.new, self.fa = wi, new, (hit < 0).astype(float)
        self.end = np.r_[np.nonzero(sc[1:] != sc[:-1])[0], len(sc) - 1]
        self.k = np.array([w["k"] for w in wins], float); self.sc = sc
    def recall(self, C, caps=(0.1, 0.25)):
        C = np.atleast_2d(C).astype(float); res = {c: np.zeros(len(C)) for c in caps}
        for a in range(0, len(C), 50):
            c = C[a:a + 50]; wt = c[:, self.wi]
            Hh = np.cumsum(wt * self.new, 1)[:, self.end]; Ff = np.cumsum(wt * self.fa, 1)[:, self.end]
            n = c.sum(1)[:, None]; kk = np.maximum(1, c @ self.k)
            for cap in caps:
                res[cap][a:a + 50] = np.where(Ff <= cap * n + 1e-9, Hh, 0).max(1) / kk
        return res
    def brute(self, cap):
        """explicit threshold loop, point estimate only (cross-check of the vectorised code)."""
        n = len(self.k); best = 0.0
        for t in np.unique(self.sc):
            s = self.sc >= t
            if self.fa[s].sum() <= cap * n + 1e-9: best = max(best, self.new[s].sum())
        return best / self.k.sum()

def boot_counts(wins, rng, B_):
    """hierarchical: strata regime x cell; nets with replacement; within each drawn net its windows with replacement."""
    C = np.zeros((B_, len(wins)))
    groups = {}
    for i, w in enumerate(wins): groups.setdefault(w["st"], {}).setdefault(w["net"], []).append(i)
    for st, nets in groups.items():
        ids = np.array([nets[j] for j in sorted(nets)])            # (nnet, nwin)
        nn, nw = ids.shape
        for b in range(B_):
            pk = rng.integers(0, nn, nn)
            ww = rng.integers(0, nw, (nn, nw))
            np.add.at(C[b], ids[pk[:, None], ww].ravel(), 1.0)
    return C

def contrast(wins, a, b, C, labels=None, score=None):
    fa = Froc(wins, a, None if labels is None else labels[a], score); fb = Froc(wins, b, None if labels is None else labels[b], score)
    one = np.ones((1, len(wins)))
    pa, pb = fa.recall(one), fb.recall(one); ra, rb = fa.recall(C), fb.recall(C)
    out = {}
    for cap in (0.1, 0.25):
        d = ra[cap] - rb[cap]
        out[cap] = dict(pa=pa[cap][0], pb=pb[cap][0], pt=pa[cap][0] - pb[cap][0], lo=np.percentile(d, 2.5), hi=np.percentile(d, 97.5),
                        lo90=np.percentile(d, 5), hi90=np.percentile(d, 95), sd=d.std())
    return out

def fmt(c, cap=0.25):
    x = c[cap]; return f"{x['pa']:.4f} - {x['pb']:.4f} = {x['pt']:+.4f} [95: {x['lo']:+.4f}, {x['hi']:+.4f}] [90: {x['lo90']:+.4f}, {x['hi90']:+.4f}] sd {x['sd']:.4f}"
