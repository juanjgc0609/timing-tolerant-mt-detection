"""C9 P4 analysis.  Records are built from the stored atoms; exact FROC recall at <= 0.1 / <= 0.25 realised false events
per window (xc.recall_caps), hierarchical bootstrap (strata regime x cell; networks then windows; operating point
re-chosen per replicate; B = 1000), paired differences to the base detector, fraction of the CEIL-GD headroom captured,
epicentral error of hits, relocation conversion diagnostics.
    python c9_analyze.py DESIGN [...] > out/c9_analysis_DESIGN.txt"""
import json
import os
import sys

import numpy as np

import c9_lib as C
import xc

B = int(os.environ.get("NBOOT", 1000))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
GARR = np.asarray(xc.make_grid())
BASES = ("BOMP-PK", "MT-FW")
RL = [f"{s}_{t}" for t in ("n", "w") for s in ("mt", "pol", "ppol", "amp", "mtos")]
INLOOP = {"MT-FW": ("MT-FW+RLIN-mt", "MT-FW+RLIN-pol", "MT-FW+PEN", "MT-FW+SELPOL", "MT-FW+PENPOL"),
          "BOMP-PK": ("BOMP-PK+RLIN-mt",)}


STRICT_KEEP = ("RL-mt_n", "RL-pol_w", "RL-amp_w")
WORLD_KEEP = ("RR-pol", "RR-pol^2", "RR-mt", "RL-mt_n", "RL-mt_w", "RL-pol_n", "RL-pol_w", "RL-ppol_n", "RL-amp_n",
              "RL-amp_w")
SHORT_KEEP = ("RR-pol", "RR-pol^2", "RR-mt", "RL-mt_n", "RL-pol_w")


def load(design):
    rows = []
    for reg in ("thin", "strong"):
        p = os.path.join(OUT, f"c9_{design}_{reg}.json")
        if os.path.exists(p):
            px = os.path.join(OUT, f"c9x_{design}_{reg}.json")
            ex = {}
            if os.path.exists(px):
                ex = {(r["cell"], r["net"], r["w"]): r for r in json.load(open(px))["rows"]}
            for r in json.load(open(p))["rows"]:
                r["regime"] = reg
                e = ex.get((r["cell"], r["net"], r["w"]))
                if e is not None:
                    assert e["truth"] == r["truth"]
                    b0, b1 = r["atoms"]["MT-FW"], e["atoms"]["MT-FW"]
                    assert [a["n"] for a in b0] == [a["n"] for a in b1] and np.allclose([a["s"] for a in b0],
                                                                                         [a["s"] for a in b1], rtol=1e-9)
                    for m in ("MT-FW+SELPOL", "MT-FW+PENPOL"):
                        r["atoms"][m] = e["atoms"][m]
                rows.append(r)
    return rows


def build(row):
    tr = [np.asarray(t) for t in row["truth"]]
    rec, loc = {}, {}

    def add(name, atoms, key="n", **kw):
        rec[name] = C.record(atoms, GARR, tr, key=key, **kw)
        loc[name] = [GARR[a["n"] if isinstance(key, str) else key[k]] for k, a in enumerate(atoms)]
    for b in BASES:
        at = row["atoms"][b]
        add(b, at)
        add(b + "/raw", at, runmin=False)
        add(b + "/RR-mt", at, runmin=False, score=[a["s"] * a["f"]["rho_mt"] for a in at])
        add(b + "/RR-pol", at, runmin=False, score=[a["s"] * a["f"]["rho_pol"] for a in at])
        for p in (2, 4):                                                    # post hoc exponents
            add(b + f"/RR-pol^{p}", at, runmin=False, score=[a["s"] * a["f"]["rho_pol"] ** p for a in at])
        add(b + "/RR-pol-only", at, runmin=False, score=[a["f"]["rho_pol"] for a in at])
        for k in RL:
            key = [a["f"][k] for a in at]
            add(f"{b}/RL-{k}", at, key=key)
            add(f"{b}/RL-{k}/strict", at, key=key, dup="strict")
        for m in INLOOP[b]:
            if m in row["atoms"]:
                add(m, row["atoms"][m])
    if "ceil" in row:
        rec["CEIL-GD"] = [tuple(x) for x in row["ceil"]]
    return rec, loc


def boot(rows, recs, meths, seed=17):
    keys = sorted({(r["regime"], r["cell"], r["net"]) for r in rows})
    idx = {k: [] for k in keys}
    for i, r in enumerate(rows):
        idx[(r["regime"], r["cell"], r["net"])].append(i)
    strata = sorted({k[:2] for k in keys})
    ks = [r["k"] for r in rows]

    def rc(ws):
        n = sum(ks[w] for w in ws)
        return {m: xc.recall_caps([recs[w][m] for w in ws], n)[0] for m in meths}
    pt = rc(list(range(len(rows))))
    rng = np.random.default_rng(seed)
    reps = {m: [] for m in meths}
    for _ in range(B):
        ws = []
        for s in strata:
            nets = [k for k in keys if k[:2] == s]
            for g in rng.integers(0, len(nets), len(nets)):
                gg = idx[nets[g]]
                ws.extend(np.asarray(gg)[rng.integers(0, len(gg), len(gg))].tolist())
        r_ = rc(ws)
        for m in meths:
            reps[m].append(r_[m])
    return pt, {m: np.array(v) for m, v in reps.items()}


def f3(p, lo, hi, sign=True):
    return f"{p:+.3f} [{lo:+.3f}, {hi:+.3f}]" if sign else f"{p:.3f} [{lo:.3f}, {hi:.3f}]"


def hit_err(rows, recs, locs, m):
    e = []
    for r, rc, lc in zip(rows, recs, locs):
        if m not in lc:
            continue
        for (s, h), x in zip(rc[m], lc[m]):
            if h >= 0:
                e.append(np.linalg.norm(x[:2] - np.asarray(r["truth"][h])[:2]))
    return (float(np.median(e)), len(e)) if e else (np.nan, 0)


def conversion(rows, base, k):
    """atoms of the base detector by distance d0 of the entry node to the nearest true event: hits (d0 <= 5) kept /
    lost by the relocation; near misses (5 < d0 <= 15) converted to hits; far atoms (d0 > 15) converted."""
    c = {"hit": [0, 0], "near": [0, 0], "far": [0, 0]}
    dn0, dn1 = [], []
    for r in rows:
        tr = np.asarray(r["truth"])
        if len(tr) == 0:
            continue
        for a in r["atoms"][base]:
            d0 = np.linalg.norm(tr[:, :2] - GARR[a["n"]][None, :2], axis=1).min()
            d1 = np.linalg.norm(tr[:, :2] - GARR[a["f"][k]][None, :2], axis=1).min()
            cls = "hit" if d0 <= 5 else ("near" if d0 <= 15 else "far")
            c[cls][0] += 1
            c[cls][1] += int(d1 <= 5)
            if cls == "near":
                dn0.append(d0); dn1.append(d1)
    return c, (float(np.median(dn0)) if dn0 else np.nan), (float(np.median(dn1)) if dn1 else np.nan)


def block(rows, title, full=True):
    built = [build(r) for r in rows]
    recs, locs = [b[0] for b in built], [b[1] for b in built]
    meths = [m for m in recs[0] if all(m in rc for rc in recs)]
    meths = [m for m in meths if not m.endswith("/strict") or m.split("/")[1] in STRICT_KEEP]
    if not full:
        keep = WORLD_KEEP if os.environ.get("WORLD") else SHORT_KEEP
        meths = [m for m in meths if "/" not in m or m.split("/", 1)[1] in keep]
    nw, ne = len(rows), sum(r["k"] for r in rows)
    nn = len({(r["regime"], r["cell"], r["net"]) for r in rows})
    print(f"\n{title}: windows {nw}, events {ne}, network-regime units {nn}")
    pt, reps = boot(rows, recs, meths)
    ci = lambda a: np.percentile(a, [2.5, 97.5], axis=0)
    has_ceil = "CEIL-GD" in meths
    print("  recall @0.1 / @0.25 [95%]; paired difference to the base @0.1 / @0.25 [95%]; headroom fraction @0.25; "
          "median epicentral error of hits km (n hits)")
    for m in meths:
        base = "BOMP-PK" if m.startswith("BOMP-PK") else ("MT-FW" if m.startswith("MT-FW") else "BOMP-PK")
        c = ci(reps[m])
        d = reps[m] - reps[base]
        cd = ci(d)
        pd_ = np.subtract(pt[m], pt[base])
        line = (f"    {m:28s} {f3(pt[m][0], c[0][0], c[1][0], False)}  {f3(pt[m][1], c[0][1], c[1][1], False)}")
        if m != base:
            line += f"  | {f3(pd_[0], cd[0][0], cd[1][0])}  {f3(pd_[1], cd[0][1], cd[1][1])}"
        else:
            line += " " * 54
        if has_ceil and m != "CEIL-GD":
            H = reps["CEIL-GD"][:, 1] - reps["BOMP-PK"][:, 1]
            G = reps[m][:, 1] - reps["BOMP-PK"][:, 1]
            with np.errstate(divide="ignore", invalid="ignore"):
                fr = np.where(np.abs(H) > 1e-9, G / H, np.nan)
            v = fr[np.isfinite(fr)]
            Hp = pt["CEIL-GD"][1] - pt["BOMP-PK"][1]
            line += f"  | {(pt[m][1] - pt['BOMP-PK'][1]) / Hp:+.2f} [{np.percentile(v, 2.5):+.2f}, {np.percentile(v, 97.5):+.2f}]"
        he, nh = hit_err(rows, recs, locs, m)
        line += f"  | {he:.2f} ({nh})"
        print(line)
    d = reps["MT-FW"] - reps["BOMP-PK"]
    cd = ci(d)
    pd_ = np.subtract(pt["MT-FW"], pt["BOMP-PK"])
    print(f"    MT-FW - BOMP-PK: @0.1 {f3(pd_[0], cd[0][0], cd[1][0])}  @0.25 {f3(pd_[1], cd[0][1], cd[1][1])}")
    if full:
        print("  relocation diagnostics (post hoc relocations of the base atoms): hits kept / near misses (5-15 km) converted "
              "/ far atoms (> 15 km) converted; median distance of near misses to the nearest event before -> after")
        for b in BASES:
            for k in RL:
                c, d0, d1 = conversion(rows, b, k)
                print(f"    {b:8s} RL-{k:7s} hits kept {c['hit'][1]}/{c['hit'][0]} ({c['hit'][1] / max(1, c['hit'][0]):.2f})  "
                      f"near converted {c['near'][1]}/{c['near'][0]} ({c['near'][1] / max(1, c['near'][0]):.2f})  far converted "
                      f"{c['far'][1]}/{c['far'][0]} ({c['far'][1] / max(1, c['far'][0]):.2f})  near median {d0:.1f} -> {d1:.1f} km")
        print("  consistency fractions of the base atoms (median): hits vs false")
        for b in BASES:
            for key in ("rho_mt", "rho_pol"):
                h, f_ = [], []
                for r, rc in zip(rows, recs):
                    for a, (s, hh) in zip(r["atoms"][b], rc[b]):
                        (h if hh >= 0 else f_).append(a["f"][key])
                if h and f_:
                    # AUC of the fraction as a hit/false discriminator
                    hh_, ff_ = np.array(h), np.array(f_)
                    auc = float(np.mean(hh_[:, None] > ff_[None, :]) + 0.5 * np.mean(hh_[:, None] == ff_[None, :]))
                    print(f"    {b:8s} {key:8s} hits {np.median(h):.3f} (n {len(h)})  false {np.median(f_):.3f} (n {len(f_)})  "
                          f"AUC {auc:.3f}")
    sec = {}
    for r in rows:
        for k, v in r["sec"].items():
            sec.setdefault(k, []).append(v)
    print("  runtime s/window: " + ", ".join(f"{k} {np.mean(v):.3f}" for k, v in sec.items()))


def main():
    for design in sys.argv[1:]:
        rows = load(design)
        if not rows:
            continue
        print("=" * 150)
        print(f"DESIGN {design}")
        world = bool(os.environ.get("WORLD"))
        block(rows, f"{design} pooled (thin + strong)", full=not world)
        for reg in ("thin", "strong"):
            sub = [r for r in rows if r["regime"] == reg]
            if sub and not world:
                block(sub, f"{design} {reg}", full=False)
        for cell in ("ladder_4", "ladder_6", "unilateral_c40"):
            sub = [r for r in rows if r["cell"] == cell]
            if sub and os.environ.get("PERCELL"):
                block(sub, f"{design} {cell} (both regimes)", full=False)


if __name__ == "__main__":
    main()
