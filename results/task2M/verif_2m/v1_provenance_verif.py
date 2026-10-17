"""V1 provenance: parts vs configs, slots vs meta eligibility, slot disjointness, catalogue margins, counts."""
import glob
import json
import os
import sys

import numpy as np

W = "/home/claude/work"
OUT = os.path.join(W, "t2m_verif", "v1_provenance_verif.txt")
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
lines = []


def P(s=""):
    print(s, flush=True)
    lines.append(s)


def own_blocks(cfg):
    blocks = []
    metas = []
    for di, d in enumerate(cfg["noise_dirs"]):
        m = json.load(open(os.path.join(d, "meta.json")))
        metas.append(m)
        el = np.array(m["eligible"])
        n = cfg["nblocks"][di]
        # own array_split: sizes q+1 for the first r blocks
        q, r = divmod(len(el), n)
        a = 0
        for b in range(n):
            sz = q + (1 if b < r else 0)
            blocks.append((di, el[a:a + sz]))
            a += sz
        assert a == len(el)
    return blocks, metas


for seed in ("2718", "3141"):
    cfg = json.load(open(os.path.join(W, "t2m_out", f"cfg_{seed}.json")))
    blocks, metas = own_blocks(cfg)
    nets = list(range(cfg["j0"], cfg["j0"] + len(blocks)))
    pdir = os.path.join(W, "t2m_out", f"parts_{seed}")
    P("=" * 100)
    P(f"SEED {seed}: config nets {nets[0]}-{nets[-1]} ({len(nets)}), days " +
      ", ".join(f"{m['day']} ({len(m['eligible'])} eligible, {len(m['stations'])} stations, {len(m['events_s'])} cat. events)"
                for m in metas))
    present = {}
    nerr = 0
    for cell in CELLS:
        for j in nets:
            fn = os.path.join(pdir, f"REAL_{cell}_{j}.json")
            if not os.path.exists(fn):
                P(f"  MISSING REAL {cell} {j}")
                continue
            d = json.load(open(fn))
            present[(cell, j)] = d
    P(f"  REAL parts present {len(present)}; per cell " +
      ", ".join(f"{c} {sum(1 for k in present if k[0] == c)}" for c in CELLS))
    other = sorted(set(os.path.basename(p) for p in glob.glob(os.path.join(pdir, "REAL_*.json"))) -
                   set(f"REAL_{c}_{j}.json" for c, j in present))
    P(f"  unexpected REAL files: {other}")
    nwin = cfg["nwin"]
    slot_use_by_cell = {c: {} for c in CELLS}
    margin_viol = 0
    for (cell, j), d in sorted(present.items()):
        di, blk = blocks[j - cfg["j0"]]
        m = metas[di]
        errs = []
        if d["seed"] != cfg["seed"]:
            errs.append("seed")
        if d["day"] != m["day"]:
            errs.append(f"day {d['day']} vs {m['day']}")
        if d["arm"] != "REAL" or d["cell"] != cell or d["net"] != j:
            errs.append("header")
        if d["stations"] != [m["stations"][i] for i in d["sub"]]:
            errs.append("stations")
        sta = np.stack([np.array(m["x_km"])[d["sub"]], np.array(m["y_km"])[d["sub"]], np.zeros(len(d["sub"]))], 1)
        if not np.allclose(sta, np.array(d["STA"])):
            errs.append("STA")
        nsta = {"ladder_4": 4, "ladder_6": 6, "unilateral_c40": 8}[cell]
        if len(d["sub"]) != nsta or len(set(d["sub"])) != nsta:
            errs.append("subset size")
        pos = np.arange(len(blk))
        sl = {"cat": blk[pos % 4 == 0].tolist(), "STD": blk[pos % 4 == 1].tolist(),
              "DENSE": blk[pos % 4 == 2].tolist(), "NOISE": blk[pos % 4 == 3].tolist()}
        if sl != d["slots"]:
            errs.append("slots")
        # windows
        rows = d["rows"]
        cnt = {k: sum(r["design"] == k for r in rows) for k in ("STD", "DENSE", "NOISE")}
        if cnt != nwin:
            errs.append(f"window counts {cnt}")
        used = []
        for r in rows:
            if r["slot"] != sl[r["design"]][r["w"]]:
                errs.append(f"row slot {r['design']} {r['w']}")
            used.append(r["slot"])
            if r["k"] != len(r["truth"]) or r["k"] != len(r["taus"]):
                errs.append("k")
            if set(r["atoms"]) != set(cfg["methods"]):
                errs.append("methods")
        allused = used + list(d["cat_slots"])
        if len(set(allused)) != len(allused):
            errs.append(f"slot reused within network ({len(allused) - len(set(allused))})")
        elig = set(m["eligible"])
        if not set(allused) <= elig:
            errs.append("non-eligible slot")
        if not set(allused) <= set(blk.tolist()):
            errs.append("slot outside block")
        # independent catalogue-pool reconstruction
        un = set(sl["NOISE"][:nwin["NOISE"]])
        pool = sl["cat"] + [s for s in sl["NOISE"][::-1] if s not in un] + sl["DENSE"][nwin["DENSE"]:][::-1] + \
            sl["STD"][nwin["STD"]:][::-1]
        if not set(d["cat_slots"]) <= set(pool[:d["cat_drawn"]]):
            errs.append("cat slots not in pool prefix")
        if d["cat_drawn"] > len(pool):
            errs.append("drawn > pool")
        # catalogue-margin check on all used slots (own)
        ev = np.array(m["events_s"])
        for s in allused:
            a, b = s * 40.0, (s + 1) * 40.0
            if a < 60 or b > 86400 - 60 or (len(ev) and np.any((ev > a - 60) & (ev < b + 60))):
                margin_viol += 1
        for s in used:
            slot_use_by_cell[cell].setdefault((m["day"], s), []).append(j)
        if errs:
            nerr += 1
            P(f"  {cell} {j}: ERRORS {errs}")
    P(f"  parts with errors: {nerr}; slots violating the catalogue/day-end margins: {margin_viol}")
    for c in CELLS:
        dup = {s: v for s, v in slot_use_by_cell[c].items() if len(set(v)) > 1}
        P(f"  {c}: window (day, slot) pairs shared between networks: {len(dup)}")
    # pool size per network vs drawn
    dr = [(c, j, d["cat_drawn"]) for (c, j), d in present.items()]
    P(f"  catalogue draws: mean {np.mean([x[2] for x in dr]):.1f}, max {max(x[2] for x in dr)}; "
      f"pool sizes min {min(len(blocks[j - cfg['j0']][1]) // 4 + 0 for j in nets)} (cat list only)")
    for cell in CELLS:
        for j in nets:
            if (cell, j) not in present:
                di, blk = blocks[j - cfg["j0"]]
                pos = np.arange(len(blk))
                sl = {"cat": blk[pos % 4 == 0].tolist(), "STD": blk[pos % 4 == 1].tolist(),
                      "DENSE": blk[pos % 4 == 2].tolist(), "NOISE": blk[pos % 4 == 3].tolist()}
                un = set(sl["NOISE"][:nwin["NOISE"]])
                pool = sl["cat"] + [s for s in sl["NOISE"][::-1] if s not in un] + sl["DENSE"][nwin["DENSE"]:][::-1] + \
                    sl["STD"][nwin["STD"]:][::-1]
                P(f"  missing {cell} {j}: block {len(blk)} slots, catalogue pool size {len(pool)} (log reports drawn = pool)")

with open(OUT, "w") as fh:
    fh.write("\n".join(lines) + "\n")
