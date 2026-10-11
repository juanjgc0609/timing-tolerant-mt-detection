"""Analysis of checks 3 (budget) and 4 (ceiling) of the 2E verification; own FROC code (as froc_rederive_verif.py).
    python3 analyze_verif.py"""
import glob
import json
import os

import numpy as np

HERE = "/home/claude/work/verif_2e"
CAPS = (0.1, 0.25)
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
LAB = {"ladder_4": "4 st.", "ladder_6": "6 st.", "unilateral_c40": "one-sided 8"}
L = []


def P_(s=""):
    print(s, flush=True)
    L.append(s)


def HFmat(recs):
    """per window and threshold (all distinct scores, descending): distinct true events hit and false atoms emitted."""
    allsc = [s for r in recs for s, h in r]
    thr = np.array(sorted(set(allsc), reverse=True)) if allsc else np.array([0.0])
    H = np.zeros((len(recs), len(thr)))
    F = np.zeros((len(recs), len(thr)))
    for w, r in enumerate(recs):
        best = {}
        fs = []
        for s, h in r:
            if h >= 0:
                best[h] = max(best.get(h, -np.inf), s)
            else:
                fs.append(s)
        if best:
            H[w] = (thr[None, :] <= np.array(list(best.values()))[:, None]).sum(0)
        if fs:
            F[w] = (thr[None, :] <= np.array(fs)[:, None]).sum(0)
    return H, F


def recall(HF, C, k):
    H, F = HF
    Hs, Fs = C @ H, C @ F
    n = C.sum(1)[:, None]
    out = {}
    for c in CAPS:
        r = np.where(Fs / n <= c + 1e-12, Hs, 0).max(1)
        out[c] = r / np.maximum(1, C @ k)
    return out


def boot_counts(n, B=2000, seed=0, groups=None):
    rng = np.random.default_rng(seed)
    C = np.zeros((B, n))
    if groups is None:
        idx = rng.integers(0, n, (B, n))
        np.add.at(C, (np.repeat(np.arange(B), n), idx.ravel()), 1.0)
        return C
    ids = sorted(set(groups))
    pos = {g: np.where(np.asarray(groups) == g)[0] for g in ids}
    for b in range(B):
        for gg in rng.integers(0, len(ids), len(ids)):
            p = pos[ids[gg]]
            np.add.at(C[b], p[rng.integers(0, len(p), len(p))], 1.0)
    return C


def compare(recs_a, recs_b, k, groups, seed):
    one = np.ones((1, len(k)))
    Ha, Hb = HFmat(recs_a), HFmat(recs_b)
    pa, pb = recall(Ha, one, k), recall(Hb, one, k)
    C = boot_counts(len(k), seed=seed, groups=groups)
    ba, bb = recall(Ha, C, k), recall(Hb, C, k)
    out = {}
    for c in CAPS:
        d = ba[c] - bb[c]
        out[c] = (pa[c][0], pb[c][0], pa[c][0] - pb[c][0], np.percentile(d, 2.5), np.percentile(d, 97.5))
    return out


# ============================================================================================ check 3
def check3():
    files = sorted(glob.glob(os.path.join(HERE, "check3_parts_verif", "b_*.json")))
    if not files:
        return
    W = []
    for f in files:
        p = json.load(open(f))
        for i, w in enumerate(p["windows"]):
            w["regime"], w["cell"], w["net"] = p["regime"], p["cell"], p["net"]
            W.append(w)
    V_ = ("REG", "X10", "MS", "X10MS")
    P_("## Check 3. Budget handicap (seed 999, 2A calibration networks 0-4, windows 0-%d of each network, %d windows)"
       % (max(1, len(W) // 30) - 1, len(W)))
    same = sum(w["rec"]["REG"] == w["rec"]["stored_REG"] for w in W)
    same_hits = sum([h for s, h in w["rec"]["REG"]] == [h for s, h in w["rec"]["stored_REG"]] for w in W)
    P_(f"- reproducibility: registered TS-SFW+SC recomputed here == stored calibration record bit-for-bit in "
       f"{same}/{len(W)} windows; identical hit pattern in {same_hits}/{len(W)} (the rest differ by floating-point "
       f"path differences amplified by the 50-iteration L-BFGS-B)")
    # diagnostics
    for v in V_:
        sec = np.mean([w["diag"][v]["sec"] for w in W])
        J1 = [w["diag"][v]["J"][0] for w in W if w["diag"][v]["J"]]
        J1r = [w["diag"]["REG"]["J"][0] for w in W if w["diag"][v]["J"]]
        dJ = np.array([(a - b) / abs(b) for a, b in zip(J1, J1r)])
        lit = [x for w in W for it in w["diag"][v]["loc_nit"] for x in it]
        cap = 500 if "X10" in v else 50
        rit = [x for w in W for x in w["diag"][v]["ref_nit"]]
        alt = [x for w in W for x in w["diag"][v]["alt"]]
        acap = 20 if "X10" in v else 5
        fb = np.mean([w["diag"][v]["frac_bound"] for w in W])
        P_(f"- {v:6s}: {sec:5.1f} s/window; J after atom 1 relative to REG: median {np.median(dJ):+.4f}, P10 "
           f"{np.percentile(dJ, 10):+.4f} (negative = lower objective); sliding L-BFGS-B at cap {np.mean(np.array(lit) >= cap):.2f}; "
           f"refinement at cap {np.mean(np.array(rit) >= cap):.2f}; alternations at cap {np.mean(np.array(alt) >= acap):.2f}; "
           f"shifts at the +-0.25 s bound {fb:.3f}")
    k = np.array([w["k"] for w in W], float)
    groups = [f"{w['regime']}_{w['cell']}_{w['net']}" for w in W]
    P_("")
    P_("Recall at <= 0.1 / <= 0.25 false per window (exact FROC, one curve per pool); differences with 95% CI "
       "(bootstrap over networks then windows, 2000 replicates, operating point re-chosen)")
    for pool, sel in (("all", lambda w: True), ("thin", lambda w: w["regime"] == "thin"),
                      ("strong", lambda w: w["regime"] == "strong")):
        idx = [i for i, w in enumerate(W) if sel(w)]
        kk = k[idx]
        gg = [groups[i] for i in idx]
        recs = {m: [W[i]["rec"][m] for i in idx] for m in V_ + ("stored_REG", "BOMP-PK+SC", "MT-FW+SC",
                                                                  "BOMP-PK-SLIDE+SC")}
        one = np.ones((1, len(kk)))
        pts = {m: recall(HFmat(recs[m]), one, kk) for m in recs}
        P_(f"- pool {pool} ({len(idx)} windows, {int(kk.sum())} events): " +
           "; ".join(f"{m} {pts[m][0.1][0]:.3f}/{pts[m][0.25][0]:.3f}" for m in recs))
        for v in ("X10", "MS", "X10MS"):
            c1 = compare(recs[v], recs["REG"], kk, gg, 31)
            c2 = compare(recs[v], recs["BOMP-PK+SC"], kk, gg, 32)
            P_(f"    {v:6s} - REG: " + ", ".join(f"{c1[c][2]:+.3f} [{c1[c][3]:+.3f}, {c1[c][4]:+.3f}]" for c in CAPS) +
               f" | {v} - BOMP-PK+SC: " + ", ".join(f"{c2[c][2]:+.3f} [{c2[c][3]:+.3f}, {c2[c][4]:+.3f}]" for c in CAPS))
        c3 = compare(recs["REG"], recs["BOMP-PK+SC"], kk, gg, 33)
        c4 = compare(recs["REG"], recs["stored_REG"], kk, gg, 34)
        P_(f"    REG    - BOMP-PK+SC: " + ", ".join(f"{c3[c][2]:+.3f} [{c3[c][3]:+.3f}, {c3[c][4]:+.3f}]" for c in CAPS)
           + " | REG(recomputed) - REG(stored): " + ", ".join(f"{c4[c][2]:+.3f}" for c in CAPS))
    # hit changes per atom
    n_ch = {v: [0, 0] for v in ("X10", "MS", "X10MS")}
    for w in W:
        hr = {h for s, h in w["rec"]["REG"] if h >= 0}
        for v in n_ch:
            hv = {h for s, h in w["rec"][v] if h >= 0}
            n_ch[v][0] += len(hv - hr)
            n_ch[v][1] += len(hr - hv)
    P_("- distinct true events detected (any score) gained / lost relative to REG: " +
       "; ".join(f"{v} +{a}/-{b}" for v, (a, b) in n_ch.items()) + f" (of {int(k.sum())} events)")
    P_("")


# ============================================================================================ check 4
def ceiling_records(w, true_key, noise_key):
    r = [(s, i) for i, s in enumerate(w[true_key]) if np.isfinite(s)] + [(s, -1) for s in w[noise_key]]
    r.sort(key=lambda x: -x[0])
    return r[:5]


def realfalse_records(w, true_key, method):
    """true candidates at their oracle score + the method's own false atoms (score as recorded)."""
    r = [(s, i) for i, s in enumerate(w[true_key])] + [(s, h) for s, h in w["rec"][method] if h < 0]
    r.sort(key=lambda x: -x[0])
    return r[:5]


def check4():
    files = sorted(glob.glob(os.path.join(HERE, "check4_parts_verif", "c_*.json")))
    if not files:
        return
    W = []
    for f in files:
        p = json.load(open(f))
        for w in p["windows"]:
            w["regime"], w["cell"], w["net"] = p["regime"], p["cell"], p["net"]
            W.append(w)
    P_("## Check 4. Headroom / ceiling (seed 999, 2A calibration networks 0-4 x 20 windows = %d windows)" % len(W))
    P_("Oracle ('genie') detectors, exact FROC, recall at <= 0.1 / <= 0.25 false per window:")
    P_("  CEIL-G     true events scored by the best BOMP-PK block statistic within 5 km of each source (grid, +-0.25 s "
       "peak reading, full data); false candidates = greedy block maxima on the clutter+noise-only data (event signal "
       "removed exactly). Perfect localisation + perfect deflation: generous ceiling for any grid block detector.")
    P_("  CEIL-GD    as CEIL-G, false candidates from the full data after an ORACLE joint LS deflation of the true "
       "events (truth-start continuous atoms): keeps the false alarms caused by event-model mismatch.")
    P_("  CEIL-C     continuous version: true events scored by the single-atom off-grid energy (TS-SFW refinement from "
       "the truth, 500 it.); false candidates refined the same way on the clutter+noise-only data.")
    P_("  CEIL-G+BF  true events at the CEIL-G score + BOMP-PK+SC's own false atoms (its realised false-alarm process; "
       "can double-count a mislocated event, so it is a lower-side ceiling).")
    names = ("BOMP-PK+SC", "TS-SFW+SC", "MT-FW+SC", "CEIL-G", "CEIL-GD", "CEIL-C", "CEIL-G+BF")

    def recs_of(w):
        return {"BOMP-PK+SC": w["rec"]["BOMP-PK+SC"], "TS-SFW+SC": w["rec"]["TS-SFW+SC f=0.01"],
                "MT-FW+SC": w["rec"]["MT-FW+SC"], "CEIL-G": ceiling_records(w, "g_true", "g_noise"),
                "CEIL-GD": ceiling_records(w, "g_true", "g_resid"), "CEIL-C": ceiling_records(w, "c_true", "c_noise"),
                "CEIL-G+BF": realfalse_records(w, "g_true", "BOMP-PK+SC")}
    P_("")
    P_("| Pool | windows | events | " + " | ".join(names) + " | headroom CEIL-G - BOMP-PK [95% CI] at 0.1 / 0.25 | "
       "CEIL-GD - BOMP-PK at 0.1 / 0.25 |")
    P_("|:--|--:|--:|" + ":-:|" * (len(names) + 2))
    res = {}
    pools = [(f"{LAB[c]} ({r})", (lambda w, c=c, r=r: w["cell"] == c and w["regime"] == r)) for r in ("thin", "strong")
             for c in CELLS] + [("thin (3 cells)", lambda w: w["regime"] == "thin"),
                                ("strong (3 cells)", lambda w: w["regime"] == "strong"), ("all", lambda w: True)]
    for name, sel in pools:
        idx = [i for i, w in enumerate(W) if sel(w)]
        kk = np.array([W[i]["k"] for i in idx], float)
        gg = [f"{W[i]['regime']}_{W[i]['cell']}_{W[i]['net']}" for i in idx]
        R = {m: [recs_of(W[i])[m] for i in idx] for m in names}
        one = np.ones((1, len(kk)))
        pts = {m: recall(HFmat(R[m]), one, kk) for m in names}
        h1 = compare(R["CEIL-G"], R["BOMP-PK+SC"], kk, gg, 41)
        h2 = compare(R["CEIL-GD"], R["BOMP-PK+SC"], kk, gg, 42)
        res[name] = {m: {c: float(pts[m][c][0]) for c in CAPS} for m in names}
        P_(f"| {name} | {len(idx)} | {int(kk.sum())} | " +
           " | ".join(f"{pts[m][0.1][0]:.3f}/{pts[m][0.25][0]:.3f}" for m in names) + " | " +
           " / ".join(f"{h1[c][2]:+.3f} [{h1[c][3]:+.3f}, {h1[c][4]:+.3f}]" for c in CAPS) + " | " +
           " / ".join(f"{h2[c][2]:+.3f} [{h2[c][3]:+.3f}, {h2[c][4]:+.3f}]" for c in CAPS) + " |")
    json.dump(res, open(os.path.join(HERE, "check4_summary_verif.json"), "w"), indent=1)
    # where does BOMP-PK lose events that the ceiling detects? classify missed events at BOMP-PK's 0.25 operating point
    P_("")
    P_("Failure anatomy of BOMP-PK+SC and TS-SFW+SC on these windows (all pools): for every true event, whether the "
       "method ever places an atom within 5 km of it (any score), and whether its best near-truth grid score exceeds "
       "the largest clutter+noise-only score of its window:")
    for m, key in (("BOMP-PK+SC", "BOMP-PK+SC"), ("TS-SFW+SC", "TS-SFW+SC f=0.01")):
        n_ev = n_det = n_strong = n_strong_det = 0
        for w in W:
            hits = {h for s, h in w["rec"][key] if h >= 0}
            nmax = max(w["g_noise"]) if w["g_noise"] else 0.0
            for i, g in enumerate(w["g_true"]):
                n_ev += 1
                n_det += i in hits
                if g > nmax:
                    n_strong += 1
                    n_strong_det += i in hits
        P_(f"- {m}: events {n_ev}; ever hit (any score, <= 5 atoms) {n_det} ({n_det / n_ev:.3f}); events whose "
           f"near-truth grid score beats every clutter+noise-only candidate of the window {n_strong} "
           f"({n_strong / n_ev:.3f}), of which ever hit {n_strong_det} ({n_strong_det / max(1, n_strong):.3f})")
    # continuous refinement from the truth: how far does it drift?
    dd = [x[1][3] for w in W for x in w["c_dbg"]]
    ok = [x[1][2] for w in W for x in w["c_dbg"]]
    P_(f"- truth-start continuous refinement: final 3D distance to the true source median {np.median(dd):.2f} km, "
       f"P90 {np.percentile(dd, 90):.2f} km; stays within 5 km epicentral in {np.mean(ok):.3f} of events (model "
       f"mismatch: 1D gradient medium, source spectra/phase, statics vs the homogeneous detector model)")
    P_("")


def check3_pooled():
    """REG vs X10MS on all budget windows (first batch: windows 0-2; extra batch: windows 3-8 of each network)."""
    W = []
    for pat in ("b_*.json", "bx_*.json"):
        for f in sorted(glob.glob(os.path.join(HERE, "check3_parts_verif", pat))):
            p = json.load(open(f))
            for w in p["windows"]:
                w["regime"], w["cell"], w["net"] = p["regime"], p["cell"], p["net"]
                W.append(w)
    if not W:
        return
    k = np.array([w["k"] for w in W], float)
    groups = [f"{w['regime']}_{w['cell']}_{w['net']}" for w in W]
    P_("### Check 3 (pooled). REG vs X10MS (x10 caps + 9-start refinement) on %d windows (windows 0-8 of each of the "
       "30 seed-999 calibration networks), %d events" % (len(W), int(k.sum())))
    same_hits = sum([h for s, h in w["rec"]["REG"]] == [h for s, h in w["rec"]["stored_REG"]] for w in W)
    P_(f"- registered recompute vs stored: identical hit pattern in {same_hits}/{len(W)} windows")
    for pool, sel in (("all", lambda w: True), ("thin", lambda w: w["regime"] == "thin"),
                      ("strong", lambda w: w["regime"] == "strong")):
        idx = [i for i, w in enumerate(W) if sel(w)]
        kk = k[idx]
        gg = [groups[i] for i in idx]
        recs = {m: [W[i]["rec"][m] for i in idx] for m in ("REG", "X10MS", "BOMP-PK+SC", "BOMP-PK-SLIDE+SC",
                                                            "MT-FW+SC")}
        one = np.ones((1, len(kk)))
        pts = {m: recall(HFmat(recs[m]), one, kk) for m in recs}
        c1 = compare(recs["X10MS"], recs["REG"], kk, gg, 51)
        c2 = compare(recs["X10MS"], recs["BOMP-PK+SC"], kk, gg, 52)
        c3 = compare(recs["REG"], recs["BOMP-PK+SC"], kk, gg, 53)
        P_(f"- {pool} ({len(idx)} windows, {int(kk.sum())} events): " +
           "; ".join(f"{m} {pts[m][0.1][0]:.3f}/{pts[m][0.25][0]:.3f}" for m in recs))
        P_("    X10MS - REG: " + ", ".join(f"{c1[c][2]:+.3f} [{c1[c][3]:+.3f}, {c1[c][4]:+.3f}]" for c in CAPS) +
           " | X10MS - BOMP-PK+SC: " + ", ".join(f"{c2[c][2]:+.3f} [{c2[c][3]:+.3f}, {c2[c][4]:+.3f}]" for c in CAPS) +
           " | REG - BOMP-PK+SC: " + ", ".join(f"{c3[c][2]:+.3f} [{c3[c][3]:+.3f}, {c3[c][4]:+.3f}]" for c in CAPS))
    gain = lost = 0
    ch_w = 0
    for w in W:
        hr = {h for s, h in w["rec"]["REG"] if h >= 0}
        hv = {h for s, h in w["rec"]["X10MS"] if h >= 0}
        gain += len(hv - hr)
        lost += len(hr - hv)
        ch_w += [h for s, h in w["rec"]["REG"]] != [h for s, h in w["rec"]["X10MS"]]
    J1 = np.array([(w["diag"]["X10MS"]["J"][0] - w["diag"]["REG"]["J"][0]) / abs(w["diag"]["REG"]["J"][0])
                   for w in W])
    Jl = np.array([(w["diag"]["X10MS"]["J"][-1] - w["diag"]["REG"]["J"][-1]) / abs(w["diag"]["REG"]["J"][-1])
                   for w in W])
    P_(f"- distinct true events ever hit: X10MS gains {gain}, loses {lost} of {int(k.sum())}; windows whose hit "
       f"pattern changes {ch_w}/{len(W)}; joint objective J after atom 1 (X10MS vs REG, relative): median "
       f"{np.median(J1):+.4f}, P10 {np.percentile(J1, 10):+.4f}, fraction improved by > 1%: {np.mean(J1 < -0.01):.2f}; "
       f"after atom 5: median {np.median(Jl):+.4f}")
    P_("")


if __name__ == "__main__":
    check3()
    check3_pooled()
    check4()
    with open(os.path.join(HERE, "analyze_verif_out.txt"), "w") as fh:
        fh.write("\n".join(L) + "\n")
