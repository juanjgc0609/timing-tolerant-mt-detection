"""
posthoc_defects.py - Task R1, declared post hoc (error register E60 and E61): two defects of the registered analysis,
measured on the stored lists of the frozen run. Nothing is re-detected; the registered numbers stand.

E60. The HYB score is "supporting picks + 0.5 x percentile of the waveform score", and the percentile is taken within
the segment being scored (hyb_dev.pct over the BASE events of the segment, and separately over the P1 atoms of the
segment). The score of one detection therefore depends on every other detection of the same segment, and the same
waveform score maps to a different value on CAL and on EVAL; the CAL threshold carried to EVAL (the deployable
threshold) is a threshold on a scale that changed. Readings: (1) the score scales of the two segments; (2) the
registered rules and the deployable threshold with picks only (W_WAVE = 0: the integer part of the score); (3) the P1
atoms with their percentile taken against the CAL atoms instead of the segment's own atoms (the raw P1 scores are
stored; the raw waveform scores of the BASE events are not, so they keep the registered within-segment percentile).

E61. The 95 per cent intervals of the registered block bootstrap resample 19 blocks of 60 s (18 full blocks and one of
40 s counted as 60 s, verification D6). Readings: (1) the registered intervals; (2) the length of the last block
counted correctly; (3) blocks of 30 s and 120 s; (4) the share of replicates at or below zero, i.e. the atom that puts
the lower bound of HYB - BASE at exactly 0.000 with one discordant target; (5) the number of blocks that carry a
discordant target; (6) a paired-proportion interval that does not resample (Newcombe 1998, method 10) for comparison.

    python code/r1/posthoc_defects.py results/taskR1/run/hyb_eval.json results/taskR1/run/hyb_cal.json \
        results/taskR1/run/p1_eval.json results/taskR1/run/p1_cal.json OUT.txt
"""
import copy
import json
import os
import sys
from math import sqrt

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r1_lib as R  # noqa: E402

R.install()
import r0_analyze as A  # noqa: E402
import r1_analyze as Y  # noqa: E402

Z = 1.959963984540054


def wilson(x, n):
    if n == 0:
        return 0.0, 1.0
    p = x / n
    d = 1 + Z * Z / n
    c = (p + Z * Z / (2 * n)) / d
    h = Z * sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / d
    return c - h, c + h


def newcombe_paired(a, b, c, d):
    """Newcombe (1998) method 10: interval for p1 - p2 of paired proportions, cells a (both), b (1 only), c (2 only),
    d (neither)."""
    n = a + b + c + d
    p1, p2 = (a + b) / n, (a + c) / n
    l1, u1 = wilson(a + b, n)
    l2, u2 = wilson(a + c, n)
    den = (a + b) * (c + d) * (a + c) * (b + d)
    num = a * d - b * c
    phi = 0.0 if den == 0 else (max(num - n / 2, 0) if num > 0 else num) / sqrt(den)
    dd = p1 - p2
    lo = dd - sqrt(max((p1 - l1) ** 2 - 2 * phi * (p1 - l1) * (u2 - p2) + (u2 - p2) ** 2, 0))
    hi = dd + sqrt(max((u1 - p1) ** 2 - 2 * phi * (u1 - p1) * (p2 - l2) + (p2 - l2) ** 2, 0))
    return lo, hi


def segment(name):
    lo, hi = R.T_CAL if name == "cal" else R.T_EVAL
    return lo, hi, (hi - lo) / 3600.0


def targets(cat, lo, hi, block):
    tg = np.array([i for i, e in enumerate(cat) if lo <= e["t"] < hi and e["M"] >= R.M_TARGET])
    return tg, (tg, np.array([int((cat[i]["t"] - lo) // block) for i in tg]))


def with_scores(d, f):
    """copy of a detection file with every HYB score replaced by f(score, detection)."""
    e = copy.deepcopy(d)
    e["HYB"] = [h[:4] + [f(h[4], h)] for h in e["HYB"]]
    return e


def rules(sets, tmask, M, hours, ntar):
    out = {}
    for m in ("HYB", "BASE", "BESTPICK", "P1"):
        out[m] = Y.best_of(sets[m], tmask, M, 30, "LOC+T", hours)[1]
    L = []
    for tag, a, b in (("R1-a", "HYB", "BASE"), ("R1-b", "HYB", "BESTPICK")):
        ga, gb = out[a][1], out[b][1]
        nb_, nc_ = len(ga - gb), len(gb - ga)
        L.append(f"{tag} {a} {len(ga)} vs {b} {len(gb)} of {ntar}, discordant {nb_}/{nc_}, "
                 f"McNemar p = {Y.mcnemar(nb_, nc_):.4f}")
    return out, L


def deployable(dc, de, cat, M):
    clo, chi, ch = segment("cal")
    elo, ehi, eh = segment("eval")
    ctg, _ = targets(cat, clo, chi, A.BLOCK)
    etg, _ = targets(cat, elo, ehi, A.BLOCK)
    cs = A.Set(dc["HYB"], cat, clo, chi)
    es = A.Set(de["HYB"], cat, elo, ehi)
    _r, got_c, thr = Y.operating(cs, (ctg, np.zeros(len(ctg), int)), M, 30, "LOC+T", ch)
    lab = es.lab["LOC+T"]
    sel = es.s > thr
    got = set(int(i) for i in lab[sel & (lab >= 0)]) & set(int(i) for i in etg)
    nf = int((sel & (lab < 0)).sum())
    return thr, len(got_c), len(ctg), len(got), len(etg), nf, nf / eh


def boot(sets, tmask, M, nb, B, seed, hours_last=None):
    """registered resampler (r1_analyze.main); hours_last: true length of the last block in blocks (None: as run)."""
    rng = np.random.default_rng(seed)
    con = (("HYB", "BASE"), ("HYB", "BESTPICK"), ("P1", "BASE"))
    vals = {c: [] for c in con}
    for _b in range(B):
        wb = np.bincount(rng.integers(0, nb, nb), minlength=nb).astype(float)
        if hours_last is not None:
            f = {m: froc_len(sets[m], tmask, M, wb, hours_last) for m in ("HYB", "BASE", "BESTPICK", "P1")}
        else:
            f = {m: A.froc([S for _k, S in sets[m]], tmask, M, nb, wb, 30, "LOC+T")
                 for m in ("HYB", "BASE", "BESTPICK", "P1")}
        for a, b in con:
            vals[(a, b)].append(f[a] - f[b])
    return {c: np.array(v) for c, v in vals.items()}


def froc_len(sets, tmask, M, wb, last_frac):
    """r0_analyze.froc with the hours of the replicate counting the last block at its true length."""
    blk = np.full(len(wb), A.BLOCK)
    blk[-1] = A.BLOCK * last_frac
    hours = float((wb * blk).sum()) / 3600.0
    cap = int(np.floor(30 * hours + 1e-9))
    tw = wb[tmask[1]]
    ntar = tw.sum()
    best = 0.0
    for _k, S in sets:
        if not len(S.s):
            continue
        w = wb[S.blk]
        lab = S.lab["LOC+T"]
        neutral = (lab >= 0) & (M[np.maximum(lab, 0)] < R.M_TARGET)
        fal = lab < 0
        fs, fw = S.s[fal], w[fal]
        o = np.argsort(-fs)
        cum = np.cumsum(fw[o])
        k = np.searchsorted(cum, cap, side="right")
        thr = fs[o][k] if k < len(fs) else -np.inf
        hit = (lab >= 0) & ~neutral & (S.s > thr)
        got = np.zeros(len(M), bool)
        got[lab[hit]] = True
        best = max(best, float((got[tmask[0]] * tw).sum() / max(ntar, 1e-12)))
    return best


def summ(v):
    return (f"median {np.median(v):+.3f}, 95% [{np.percentile(v, 2.5):+.3f}, {np.percentile(v, 97.5):+.3f}], "
            f"share <= 0 {np.mean(v <= 0):.3f}, share = 0 {np.mean(v == 0):.3f}")


def main(hep, hcp, pep, pcp, outp):
    cat = R.catalogue_local(R.frame())
    M = np.array([e["M"] for e in cat])
    de, dc = json.load(open(hep)), json.load(open(hcp))
    elo, ehi, eh = segment("eval")
    clo, chi, ch = segment("cal")
    etg, emask = targets(cat, elo, ehi, A.BLOCK)
    L = ["R1, declared post hoc defects E60 (within-segment percentile score) and E61 (19-block bootstrap)",
         "Inputs: stored lists of the frozen run; truth and catalogue from the MatchLocate2 example (checksums of",
         "data/checksums/hinet_matchlocate2_example.sha256). The registered numbers are not changed.", ""]

    # ------------------------------------------------------------------ E60
    L.append("E60. Within-segment percentile score of HYB")
    p1e = A.dedupe(Y_in(json.load(open(pep))["atoms"]["NOCAP"], elo, ehi))
    p1c = A.dedupe(Y_in(json.load(open(pcp))["atoms"]["NOCAP"], clo, chi))
    ref = np.sort([a[4] for a in p1c])
    for nm, d in (("CAL", dc), ("EVAL", de)):
        fr = np.array([h[4] - np.floor(h[4]) for h in d["HYB"]])
        L.append(f"   {nm}: BASE events {len(d['BASE'])}, P1 atoms {len(d['P1'])}; percentile steps 0.5/{len(d['BASE'])} "
                 f"(BASE) and 0.5/{len(d['P1'])} (P1); HYB fractional parts in [{fr.min():.4f}, {fr.max():.4f}]")
    for q in (0.25, 0.5, 0.75):
        L.append(f"   P1 raw score at within-segment percentile {q:.2f}: CAL {np.quantile(ref, q):.1f}, "
                 f"EVAL {np.quantile([a[4] for a in p1e], q):.1f}")
    assert len(p1e) == len(de["P1"]) and len(p1c) == len(dc["P1"])
    variants = {
        "registered (picks + 0.5 x within-segment percentile)": (de, dc),
        "picks only (W_WAVE = 0)": (with_scores(de, lambda s, h: float(np.floor(s))),
                                     with_scores(dc, lambda s, h: float(np.floor(s)))),
    }
    # P1 atoms with a CAL-referenced percentile; BASE events keep the registered value
    p1key_e = {(round(a[0], 6), round(a[1], 6), round(a[2], 6), round(a[3], 6)): a[4] for a in p1e}
    p1key_c = {(round(a[0], 6), round(a[1], 6), round(a[2], 6), round(a[3], 6)): a[4] for a in p1c}

    def calref(table, d):
        n = 0

        def f(s, h):
            nonlocal n
            k = (round(h[0], 6), round(h[1], 6), round(h[2], 6), round(h[3], 6))
            if k in table:
                n += 1
                return float(np.floor(s)) + 0.5 * float(np.searchsorted(ref, table[k], side="right")) / len(ref)
            return s
        out = with_scores(d, f)
        return out, n

    ve, ne = calref(p1key_e, de)
    vc, nc = calref(p1key_c, dc)
    variants["P1 percentile against the CAL atoms"] = (ve, vc)
    L.append(f"   HYB detections that are P1 atoms (re-scored in the third variant): EVAL {ne}, CAL {nc}")
    for nm, (e_, c_) in variants.items():
        sets = Y.load(e_, cat, elo, ehi)
        out, rl = rules(sets, emask, M, eh, len(etg))
        thr, gc, nct, ge, net_, nf, nfh = deployable(c_, e_, cat, M)
        L.append(f"   {nm}:")
        L.append(f"      in-segment LOC+T 30/h: HYB {len(out['HYB'][1])}, BASE {len(out['BASE'][1])}, "
                 f"BESTPICK {len(out['BESTPICK'][1])} of {len(etg)}")
        L += ["      " + x for x in rl]
        L.append(f"      deployable: CAL threshold {thr:.3f} ({gc}/{nct} on CAL) -> EVAL {ge}/{net_}, {nf} false "
                 f"({nfh:.1f} per hour; budget 30)")
    L.append("")

    # ------------------------------------------------------------------ E61
    L.append("E61. Block bootstrap of the contrasts (LOC+T 30/h, 2000 replicates, seed of the run)")
    sets = Y.load(de, cat, elo, ehi)
    nb = int(np.ceil((ehi - elo) / A.BLOCK))
    last = ((ehi - elo) - (nb - 1) * A.BLOCK) / A.BLOCK
    L.append(f"   registered: {nb} blocks of {A.BLOCK:.0f} s, the last one {last * A.BLOCK:.0f} s long and counted as "
             f"{A.BLOCK:.0f} s")
    v = boot(sets, emask, M, nb, Y.B, Y.SEED)
    for c, x in v.items():
        L.append(f"      {c[0]} - {c[1]}: {summ(x)}")
    L.append("   last block at its true length:")
    v2 = boot(sets, emask, M, nb, Y.B, Y.SEED, hours_last=last)
    for c, x in v2.items():
        L.append(f"      {c[0]} - {c[1]}: {summ(x)}")
    for blen in (30.0, 120.0):
        A.BLOCK = blen
        _tg, m2 = targets(cat, elo, ehi, blen)
        s2 = Y.load(de, cat, elo, ehi)
        nb2 = int(np.ceil((ehi - elo) / blen))
        lst = ((ehi - elo) - (nb2 - 1) * blen) / blen
        L.append(f"   blocks of {blen:.0f} s ({nb2} blocks, last {lst * blen:.0f} s at its true length):")
        v3 = boot(s2, m2, M, nb2, Y.B, Y.SEED, hours_last=lst)
        for c, x in v3.items():
            L.append(f"      {c[0]} - {c[1]}: {summ(x)}")
    A.BLOCK = 60.0
    out, _rl = rules(sets, emask, M, eh, len(etg))
    tb = dict(zip(etg.tolist(), emask[1].tolist()))
    for a, b in (("HYB", "BASE"), ("HYB", "BESTPICK"), ("P1", "BASE")):
        ga, gb = out[a][1], out[b][1]
        disc = sorted((ga - gb) | (gb - ga))
        n = len(etg)
        both, only_a, only_b = len(ga & gb), len(ga - gb), len(gb - ga)
        lo, hi = newcombe_paired(both, only_a, only_b, n - both - only_a - only_b)
        L.append(f"   {a} - {b}: discordant targets {only_a}/{only_b} in {len({tb[i] for i in disc})} of {nb} blocks; "
                 f"difference {(len(ga) - len(gb)) / n:+.3f}, Newcombe paired 95% [{lo:+.3f}, {hi:+.3f}]")
    open(outp, "w").write("\n".join(L) + "\n")
    print("\n".join(L), flush=True)


def Y_in(atoms, lo, hi):
    return [list(a[:5]) for a in atoms if lo <= a[3] < hi]


if __name__ == "__main__":
    main(*sys.argv[1:6])
