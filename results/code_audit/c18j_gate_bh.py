"""c18j_gate_bh.py - independent recomputation of the C18j gate reading (CAL, 30 false/h, LOC+T) from the stored atoms,
with an own labelling/FROC, plus sensitivity of the reading to (a) the owned-origin edges and (b) hits on neutral
(M < 2) events, (c) the de-duplication."""
import json
import sys

import numpy as np

sys.path.insert(0, "/home/claude/work/r0")
import r0_lib as R  # noqa: E402

D = "/home/claude/work/repo/results/explorations/c18/"


def own_label(det, cat, tol):
    out = []
    for d in det:
        best = None
        for j, e in enumerate(cat):
            dd = np.hypot(e["x"] - d[0], e["y"] - d[1])
            if dd <= 5.0 and abs(e["t"] - d[3]) <= tol + 1e-9:
                key = (0 if e["M"] >= 2.0 else 1, dd)
                if best is None or key < best[0]:
                    best = (key, j)
        out.append(-1 if best is None else best[1])
    return out


def dedupe(atoms):
    keep = []
    for r in sorted(atoms, key=lambda r: -r[4]):
        if all(np.hypot(r[0] - k[0], r[1] - k[1]) > 5.0 or abs(r[3] - k[3]) > 1.0 for k in keep):
            keep.append(r)
    return keep


def recall_at(det, cat, lo, hi, rate, tol, dd=True):
    det = [d for d in (dedupe(det) if dd else det) if lo <= d[3] < hi]
    lab = own_label(det, cat, tol)
    tg = {j for j, e in enumerate(cat) if lo <= e["t"] < hi and e["M"] >= 2.0}
    hours = np.ceil((hi - lo) / 60.0) * 60.0 / 3600.0
    cap = int(np.floor(rate * hours + 1e-9))
    order = sorted(range(len(det)), key=lambda k: -det[k][4])
    nf, got, best = 0, set(), 0
    rows = []
    i = 0
    while i < len(order):
        j = i
        s = det[order[i]][4]
        while j < len(order) and det[order[j]][4] == s:
            j += 1
        for k in order[i:j]:
            l = lab[k]
            if l < 0:
                nf += 1
            elif l in tg:
                got.add(l)
            rows.append((s, l, nf, len(got)))
        if nf <= cap:
            best = max(best, len(got))
        i = j
    return best, len(tg), cap, rows


def main():
    cat = R.catalogue_local(R.frame())
    out = []
    for f in ("c18j_detect_f2.0_U8.json", "c18j_detect_f2.0_N8.json", "c18j_detect_f3.0_U8.json",
              "c18j_detect_f3.0_N8.json", "c18i_detect_s2_f2_U8.json", "c18i_detect_s2_f2_N8.json"):
        r = json.load(open(D + f))
        for v, at in r["atoms"].items():
            for tol, rule in ((0.5, "LOC+T"), (3.0, "LOC")):
                b, n, cap, rows = recall_at(at, cat, 0.0, 900.0, 30, tol)
                b2, *_ = recall_at(at, cat, 0.0, 900.0, 30, tol, dd=False)
                out.append(f"{f:32s} {v:6s} {rule:5s} 30/h: {b}/{n} (cap {cap} false); without de-dup {b2}/{n}")
    print("\n".join(out))


if __name__ == "__main__":
    main()
