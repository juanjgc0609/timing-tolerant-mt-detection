"""analyze_1s.py - reglas S-0, S-a y S-b de la Tarea 1s.  python analyze_1s.py 200"""
import json
import sys

import numpy as np

REPS = sys.argv[1] if len(sys.argv) > 1 else "200"
CAPS = ("0.1", "0.25")
B = 2000


def boot(res, a, b, cap, rng):
    ca, cb = res["methods"][a]["chosen"][cap], res["methods"][b]["chosen"][cap]
    k = np.asarray(res["k_eval"], float)
    ha = np.asarray(ca[0], float) if ca else np.zeros_like(k)
    hb = np.asarray(cb[0], float) if cb else np.zeros_like(k)
    d0 = ha.sum() / max(1, k.sum()) - hb.sum() / max(1, k.sum())
    n = len(k)
    ds = []
    for _ in range(B):
        i = rng.integers(0, n, n)
        ds.append((ha[i].sum() - hb[i].sum()) / max(1, k[i].sum()))
    return d0, np.percentile(ds, 2.5), np.percentile(ds, 97.5)


out = []
sa_ok, sb_ok, s0_ok = {}, {}, {}
for seed in (2024, 7331):
    try:
        R = json.load(open(f"task1s_seed{seed}_r{REPS}.json"))
    except FileNotFoundError:
        R = json.load(open(f"task1s_seed{seed}.json"))
    rng = np.random.default_rng(seed)
    nsa = nsb = nauc = 0
    ci_excl = False
    base = None
    if REPS == "60":
        base = {r["scenario"]: r for r in json.load(open(f"task1r_baselines_seed{seed}.json"))}
    for r in R:
        ra = {m: r["methods"][m]["recall_at"] for m in r["methods"]}
        line = (f"{seed} {r['scenario']:15s} n={sum(r['k_eval']):4d} ev V/F={r['n_true_ev']}/{r['n_false_ev']} "
                f"AUC pol={r['auc_pol']:.2f} sh={r['auc_sh']:.2f} :: "
                + " | ".join(f"{m} " + "/".join(f"{ra[m][c]:.2f}" for c in ("0.05",) + CAPS) for m in ra))
        out.append(line)
        if base is not None:
            g = base[r["scenario"]]["methods"]["GAMMA"]["recall_at"]
            same = all(abs(g[c] - ra["GAMMA"][c]) < 1e-9 for c in ("0.05",) + CAPS)
            out.append(f"     S-0 GAMMA reproduce 1r: {same} (1r {[round(g[c], 3) for c in ('0.05',) + CAPS]})")
            s0_ok[(seed, r["scenario"])] = same
        nauc += r["auc_pol"] > 0.6 and 0.4 <= r["auc_sh"] <= 0.6
        d_ps = [ra["GAMMA+POL"][c] - ra["GAMMA+SH"][c] for c in CAPS]
        d_pg = [ra["GAMMA+POL"][c] - ra["GAMMA"][c] for c in CAPS]
        nsa += all(x >= 0.05 for x in d_ps + d_pg)
        nsb += ra["GAMMA+SH"]["0.25"] - ra["GAMMA"]["0.25"] < 0.05
        for a, b in (("GAMMA+POL", "GAMMA+SH"), ("GAMMA+POL", "GAMMA"), ("GAMMA+SH", "GAMMA")):
            d, lo, hi = boot(r, a, b, "0.25", rng)
            out.append(f"     {a} - {b} @0.25: {d:+.3f} [{lo:+.3f}, {hi:+.3f}]")
            if a == "GAMMA+POL" and b == "GAMMA+SH" and lo > 0:
                ci_excl = True
    sa_ok[seed] = nsa >= 2 and ci_excl
    sb_ok[seed] = nsb >= 2
    out.append(f"   seed {seed}: celdas S-a {nsa}/3, IC excluye 0: {ci_excl}; celdas S-b {nsb}/3; "
               f"celdas AUC ok {nauc}/3")
out.append(f"S-a cumplida: {all(sa_ok.values())}  S-b cumplida: {all(sb_ok.values())}"
           + (f"  S-0 GAMMA reproduce: {all(s0_ok.values())}" if s0_ok else ""))
txt = "\n".join(out)
print(txt)
open(f"analyze_1s_r{REPS}_out.txt", "w").write(txt + "\n")
