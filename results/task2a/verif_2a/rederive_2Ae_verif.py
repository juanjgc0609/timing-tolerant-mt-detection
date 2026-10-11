import json
import numpy as np
FZ = json.load(open("task2a_frozen_2Ae.json"))
def entries(recs):
    return [(s, w, h) for w, r in enumerate(recs) for s, h in r]
def at_thr(recs, t, K, n):
    e = entries(recs)
    hits = {(w, h) for s, w, h in e if s >= t and h >= 0}
    return len(hits) / K, sum(1 for s, w, h in e if s >= t and h < 0) / n
def freeze(recs, K, n, bud):
    best = None
    for t in sorted({s for s, _w, _h in entries(recs)}, reverse=True):
        r, f = at_thr(recs, t, K, n)
        if f <= bud + 1e-12 and (best is None or (r, -f, t) > (best[1], -best[2], best[0])):
            best = (t, r, f)
    return best
for reg in ("thin", "strong"):
    cal = {c["scenario"]: c for c in json.load(open(f"task2a_seed999_{reg}.json"))}
    for cell in ("ladder_4", "unilateral_c40"):
        c = cal[cell]; K = sum(c["k_eval"]); n = len(c["k_eval"])
        for bud in (0.1, 0.25):
            t, r, f = freeze(c["methods"]["MT-FW+SC"]["records"], K, n, bud)
            fz = FZ[f"{reg}|{cell}|{bud}"]["MT-FW+SC"]
            line = f"{reg} {cell} {bud}: my frozen thr {t:.4f} (cal {r:.3f}/{f:.2f}) vs file {fz[0]:.4f} ({fz[1]:.3f}/{fz[2]:.2f})"
            for s in (2024, 7331):
                e = {x["scenario"]: x for x in json.load(open(f"task2a_seed{s}_{reg}.json"))}[cell]
                re_, fe = at_thr(e["methods"]["MT-FW+SC"]["records"], fz[0], sum(e["k_eval"]), len(e["k_eval"]))
                tq = FZ[f"{reg}|{cell}|{bud}"]["QM-3C+SC"][0]
                rq, fq = at_thr(e["methods"]["QM-3C+SC"]["records"], tq, sum(e["k_eval"]), len(e["k_eval"]))
                line += f" | {s}: MT-FW+SC {re_:.2f}/{fe:.2f}, QM-3C+SC {rq:.2f}/{fq:.2f}"
            print(line)
