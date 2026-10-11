"""C10 analysis: recall at <= 0.1 / <= 0.25 realised false events per window with both hit rules, paired differences,
CEIL-GD headroom, epicentral error of hits, near-miss (5-10 km) share, runtime.  Hierarchical bootstrap (strata regime x
cell, networks then windows; operating point re-chosen per replicate).
    python c10_analyze.py DESIGN [B] [prefix]  ->  out/c10_analysis_DESIGN.txt"""
import json
import sys

import numpy as np

import c10_lib as C
import xc

REGIMES = ("thin", "strong")
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
CAPS = (0.1, 0.25)


def load(design, prefix="c10"):
    rows, nets = [], []
    for reg in REGIMES:
        d = json.load(open(f"out/{prefix}_{design}_{reg}.json"))
        for r in d["rows"]:
            r["regime"] = reg
            rows.append(r)
        nets += d["nets"]
    return rows, nets


def atoms_of(row, name):
    return [{"x": np.array(a[:3]), "tau": a[3], "s": a[4], "s0": a[5], "sref": a[6] if a[6] == a[6] else a[4]}
            for a in row["atoms"][name]]


def build(rows, names, tol_t):
    """records per method and window under a hit rule; '/sref' = the same atoms scored by the refined objective."""
    recs = {}
    for nm in names:
        base, _, key = nm.partition("/")
        out = []
        for r in rows:
            at = atoms_of(r, base)
            out.append(C.record(at, C.label(at, [np.array(t) for t in r["truth"]], r["taus"], tol_t), key or "s"))
        recs[nm] = out
    recs["CEIL-GD"] = [[(float(s), int(h)) for s, h in r["ceil"]] for r in rows]
    return recs


class Boot:
    def __init__(self, rows, recs, sel, B, seed=12345):
        self.methods = list(recs)
        idx = [i for i, r in enumerate(rows) if sel(r)]
        self.idx = idx
        self.ks = np.array([rows[i]["k"] for i in idx])
        g = {}
        for p, i in enumerate(idx):
            r = rows[i]
            g.setdefault((r["regime"], r["cell"], r["net"]), []).append(p)
        self.groups = list(g.values())
        self.strata = [k[:2] for k in g]
        self.S = {m: [np.array([x[0] for x in recs[m][i]], float) for i in idx] for m in self.methods}
        self.H = {m: [np.array([x[1] for x in recs[m][i]], int) for i in idx] for m in self.methods}
        allw = np.arange(len(idx))
        self.pt = {m: self.recall(m, allw) for m in self.methods}
        for m in self.methods[:3]:
            ref = xc.recall_caps([recs[m][i] for i in idx], int(self.ks.sum()), CAPS)[0]
            assert np.allclose(ref, self.pt[m]), (m, ref, self.pt[m])
        rng = np.random.default_rng(seed)
        sidx = {}
        for gi, s in enumerate(self.strata):
            sidx.setdefault(s, []).append(gi)
        self.rep = {m: np.zeros((B, len(CAPS))) for m in self.methods}
        self.ws = []
        for b in range(B):
            ws = []
            for s, gl in sidx.items():
                for gi in np.asarray(gl)[rng.integers(0, len(gl), len(gl))]:
                    gg = np.asarray(self.groups[gi])
                    ws.extend(gg[rng.integers(0, len(gg), len(gg))].tolist())
            ws = np.asarray(ws)
            self.ws.append(ws)
            for m in self.methods:
                self.rep[m][b] = self.recall(m, ws)

    def recall(self, m, ws):
        nte = int(self.ks[ws].sum())
        S, H = self.S[m], self.H[m]
        sc = np.concatenate([S[w] for w in ws])
        if not len(sc) or nte == 0:
            return [0.0] * len(CAPS)
        hit = np.concatenate([H[w] for w in ws])
        pos = np.repeat(np.arange(len(ws)), [len(S[w]) for w in ws])
        key = np.where(hit >= 0, pos * 64 + hit, -1)
        o = np.argsort(-sc, kind="stable")
        sc, hit, key = sc[o], hit[o], key[o]
        _u, first = np.unique(key, return_index=True)
        newhit = np.zeros(len(sc), bool)
        newhit[first] = True
        newhit &= key >= 0
        Hc, Fc = np.cumsum(newhit), np.cumsum(hit < 0)
        last = np.r_[np.where(np.diff(sc) != 0)[0], len(sc) - 1]
        Hc, Fc = Hc[last], Fc[last]
        out = []
        for c in CAPS:
            ok = Fc / len(ws) <= c
            out.append(float(Hc[ok].max() / nte) if ok.any() else 0.0)
        return out

    def ci(self, m):
        return np.percentile(self.rep[m], [2.5, 97.5], axis=0).T

    def diff(self, a, b):
        d = self.rep[a] - self.rep[b]
        return np.subtract(self.pt[a], self.pt[b]), np.percentile(d, [2.5, 97.5], axis=0).T, d.std(axis=0)

    def headroom(self, m, base="BOMP-PK", ceil="CEIL-GD", k=1):
        den = self.pt[ceil][k] - self.pt[base][k]
        r = (self.rep[m][:, k] - self.rep[base][:, k]) / np.maximum(self.rep[ceil][:, k] - self.rep[base][:, k], 1e-9)
        return (self.pt[m][k] - self.pt[base][k]) / den if den > 0 else np.nan, np.percentile(r, [2.5, 97.5])


def geometry(rows, names):
    """per method and window: epicentral errors of the hit atoms (location rule) and, per true event, the class
    0 = hit (an atom within 5 km), 1 = near miss (not hit; nearest FALSE atom at 5-10 km), 2 = other."""
    out = {}
    for nm in names:
        errs, cls = [], []
        for r in rows:
            at = atoms_of(r, nm)
            tr = [np.array(t) for t in r["truth"]]
            lab = C.label(at, tr, r["taus"], None)
            e = [float(np.linalg.norm(a["x"][:2] - tr[h][:2])) for a, h in zip(at, lab) if h >= 0]
            c = []
            for k, t in enumerate(tr):
                if k in lab:
                    c.append(0)
                    continue
                dd = [np.linalg.norm(a["x"][:2] - t[:2]) for a, h in zip(at, lab) if h < 0]
                c.append(1 if dd and 5.0 < min(dd) <= 10.0 else 2)
            errs.append(e)
            cls.append(c)
        out[nm] = (errs, cls)
    return out


def fmt(p, ci):
    return f"{p:+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}]"


def main(design, B=1000, prefix="c10", names=None, out=None):
    rows, nets = load(design, prefix)
    names = names or list(rows[0]["atoms"])
    ext = names + [n + "/sref" for n in names if n.endswith("+T3")]
    lines = []

    def P(s=""):
        print(s, flush=True)
        lines.append(s)
    P(f"C10 ANALYSIS {design}: {len(rows)} windows, {sum(r['k'] for r in rows)} events, "
      f"{len(nets)} network x regime setups; bootstrap B = {B}")
    for reg in REGIMES:
        sg = np.array([n["sig"][0] for n in nets if n["regime"] == reg])
        so = np.array([n["sig_ora"] for n in nets if n["regime"] == reg])
        P(f"  timing level by the catalogue rule, {reg}: median {np.median(sg):.3f} s, quartiles "
          f"{np.percentile(sg, 25):.3f}-{np.percentile(sg, 75):.3f}, range {sg.min():.3f}-{sg.max():.3f}; oracle (true "
          f"model error) P {np.median(so[:, 0]):.3f} S {np.median(so[:, 1]):.3f}")
    cs = np.array([n["cell_sd"] for n in nets])
    cf = np.array([n["cell_sd_fine"] for n in nets])
    P(f"  rms cell term (s): coarse 3 km lattice P {cs[:, 0].mean():.3f} S {cs[:, 1].mean():.3f}; fine 1 km lattice P "
      f"{cf[:, 0].mean():.3f} S {cf[:, 1].mean():.3f}")
    res = {}
    for rule, tol in (("LOC", None), ("LOC+T0.5", 0.5)):
        recs = build(rows, ext, tol)
        for grp in ("both",) + REGIMES:
            bt = Boot(rows, recs, (lambda r: True) if grp == "both" else (lambda r, g=grp: r["regime"] == g),
                      B if grp == "both" else max(200, B // 2))
            res[(rule, grp)] = bt
            P("")
            P(f"== hit rule {rule}, regimes {grp}: {len(bt.idx)} windows, {int(bt.ks.sum())} events")
            P("   recall @0.1 / @0.25 [95%]")
            for m in bt.methods:
                c = bt.ci(m)
                P(f"   {m:16s} {bt.pt[m][0]:.3f} [{c[0][0]:.3f},{c[0][1]:.3f}]  {bt.pt[m][1]:.3f} "
                  f"[{c[1][0]:.3f},{c[1][1]:.3f}]")
            P("   paired differences @0.1 / @0.25 [95%] (bootstrap SD @0.25)" + ("; fraction of the CEIL-GD headroom over "
              "BOMP-PK @0.25" if rule == "LOC" else ""))
            for m in bt.methods:
                if m == "CEIL-GD":
                    continue
                fam = "MT-FW" if m.startswith("MT") else "BOMP-PK"
                oth = "BOMP-PK" if fam == "MT-FW" else "MT-FW"
                s = f"   {m:16s}"
                for b_ in (fam, oth) + (("MT-ORA",) if m.startswith("MT-ORA") and m != "MT-ORA" else ()):
                    if b_ == m or b_ not in bt.pt:
                        continue
                    p, c, sd = bt.diff(m, b_)
                    s += f" - {b_}: {fmt(p[0], c[0])} / {fmt(p[1], c[1])} (sd {sd[1]:.4f});"
                if rule == "LOC":
                    h, hc = bt.headroom(m)
                    s += f" headroom {h:+.2f} [{hc[0]:+.2f},{hc[1]:+.2f}]"
                P(s)
    # ------------------------------------------------------------------ geometry
    P("")
    P("== location diagnostics (all atoms, location rule): median epicentral error of hit atoms (km) [95%, bootstrap over "
      "networks], share of true events hit (atom within 5 km), share NEAR-MISS (not hit, nearest false atom at 5-10 km), "
      "difference of the near-miss share to the family base [95%]")
    geo = geometry(rows, names)
    gk = {}
    for i, r in enumerate(rows):
        gk.setdefault((r["regime"], r["cell"], r["net"]), []).append(i)
    gl = list(gk.values())
    rng = np.random.default_rng(777)
    draws = [rng.integers(0, len(gl), len(gl)) for _ in range(400)]

    def stats(nm, ws):
        e = [x for w in ws for x in geo[nm][0][w]]
        c = np.array([x for w in ws for x in geo[nm][1][w]])
        return (np.median(e) if e else np.nan), np.mean(c == 0), np.mean(c == 1)
    allw = list(range(len(rows)))
    for nm in names:
        pt = stats(nm, allw)
        fam = "MT-FW" if nm.startswith("MT") else "BOMP-PK"
        bs, bd = [], []
        for d in draws:
            ws = [w for g in d for w in gl[g]]
            a = stats(nm, ws)
            bs.append(a)
            bd.append(a[2] - stats(fam, ws)[2])
        bs = np.array(bs)
        ptd = pt[2] - stats(fam, allw)[2]
        P(f"   {nm:14s} error {pt[0]:.2f} [{np.percentile(bs[:, 0], 2.5):.2f},{np.percentile(bs[:, 0], 97.5):.2f}] km   "
          f"hit {pt[1]:.3f}   near-miss {pt[2]:.3f} [{np.percentile(bs[:, 2], 2.5):.3f},{np.percentile(bs[:, 2], 97.5):.3f}]"
          f"   - {fam}: {fmt(ptd, np.percentile(bd, [2.5, 97.5]))}")
    P("")
    P("== runtime, s per window (mean; two processes on two cores): " +
      "; ".join(f"{m} {np.mean([r['sec'][m] for r in rows]):.2f}" for m in list(rows[0]['sec'])))
    P(f"   network setup (operators, catalogue, both engines): {np.mean([n['t_setup'] for n in nets]):.1f} s")
    with open(out or f"out/{prefix}_analysis_{design}.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return res


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 1000, sys.argv[3] if len(sys.argv) > 3 else "c10")
