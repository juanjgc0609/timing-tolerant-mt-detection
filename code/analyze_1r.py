"""Tarea 1r: recall a tasa de falsos realizada y IC95 por bootstrap pareado de ventanas."""
import json, os
import numpy as np
LAB = {"unilateral_c40": "one-sided, 8 st.", "ladder_4": "4 stations", "ladder_6": "6 stations"}

def rat(H, F, nt, cap):
    rec = H.sum(0) / max(1, nt); fr = F.sum(0) / H.shape[0]; ok = fr <= cap
    return float(rec[ok].max()) if ok.any() else 0.0

def boot(r, a, b, cap, B=2000, rng=np.random.default_rng(0)):
    k = np.array(r["k_eval"]); n = len(k)
    Ha, Fa = np.array(r["methods"][a]["H"]), np.array(r["methods"][a]["F"])
    Hb, Fb = np.array(r["methods"][b]["H"]), np.array(r["methods"][b]["F"])
    d0 = rat(Ha, Fa, k.sum(), cap) - rat(Hb, Fb, k.sum(), cap)
    ds = []
    for _ in range(B):
        i = rng.integers(0, n, n)
        nt = k[i].sum()
        ds.append(rat(Ha[i], Fa[i], nt, cap) - rat(Hb[i], Fb[i], nt, cap))
    lo, hi = np.percentile(ds, [2.5, 97.5])
    return d0, lo, hi

for s in (2024, 7331):
    for r in json.load(open(f"task1r_seed{s}.json")):
        k = sum(r["k_eval"])
        ra = {m: r["methods"][m]["recall_at"] for m in r["methods"]}
        line = f"{s} {r['cond']:<11} {LAB[r['scenario']]:<17} eventos {k:>3} | " + " ".join(
            f"{m}:{ra[m]['0.1']:.2f}/{ra[m]['0.25']:.2f}" for m in ("KIN-LR", "PHYS-LR", "PHYS-SH", "TOL-TB", "WAVE-P"))
        print(line)
        for a, b in (("PHYS-LR", "KIN-LR"), ("PHYS-LR", "PHYS-SH"), ("PHYS-SH", "KIN-LR")):
            d, lo, hi = boot(r, a, b, 0.25)
            print(f"      {a} - {b} a tope 0.25: {d:+.2f} [{lo:+.2f}, {hi:+.2f}]")
    fb = f"task1r_baselines_seed{s}.json"
    if os.path.exists(fb):
        for r in json.load(open(fb)):
            print(f"{s} baselines {LAB[r['scenario']]:<17} " + " ".join(
                f"{m}:{r['methods'][m]['recall_at']['0.1']:.2f}/{r['methods'][m]['recall_at']['0.25']:.2f}" for m in r["methods"]))
