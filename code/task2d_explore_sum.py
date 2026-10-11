"""task2d_explore_sum.py - resumen de E1 de task2d_explore.py (EXPLORATORIO): distancia de los falsos sobre el umbral
split-conformal congelado a la fuente verdadera mas cercana, en las ventanas de evaluacion de 2A con eventos.

    python task2d_explore_sum.py  -> task2d_explore_e1.txt
"""
import json
import os

import numpy as np

import analyze_2d as A

HERE = os.path.dirname(os.path.abspath(__file__))
BINS = ((5.0, 10.0), (10.0, 20.0), (20.0, np.inf))


def qm_sel(events, min_sep=5.0, max_events=5):
    acc, out = [], []
    for ev in sorted(events, key=lambda e: -e[0]):
        s, p = float(ev[0]), np.asarray(ev[1:3], float)
        if not np.isfinite(s) or any(np.linalg.norm(p - a) <= min_sep for a in acc):
            continue
        acc.append(p)
        out.append((s, p))
        if len(acc) >= max_events:
            break
    return out


def main():
    L = ["E1 (exploratorio): falsos sobre el umbral split-conformal congelado en las ventanas de evaluacion de 2A con "
         "eventos (k >= 1): falsos por ventana con eventos y distancia epicentral del falso a la fuente verdadera mas "
         "cercana (5-10 / 10-20 / >= 20 km). Control: registros de MT-FW+SC recalculados con nodo = registros de 2A."]
    tot_same = tot_n = 0
    for reg in A.REGIMES:
        cal = {cell: A.cal_data(reg, cell) for cell in A.CELLS}
        for s in A.SEEDS:
            for cell in A.CELLS:
                mx = cal[cell][0]
                rows = {}
                for j in range(5):
                    a = A.load(os.path.join(HERE, "task2a_parts", f"seed{s}_{reg}_{cell}_net{j}.json"))
                    e = A.load(os.path.join(HERE, "task2d_explore_parts", f"e1_seed{s}_{reg}_{cell}_net{j}.json"))
                    for wa, we in zip(a["windows"], e["windows"]):
                        if wa["k"] == 0:
                            continue
                        tot_n += 1
                        tot_same += we["same"]
                        tr = np.asarray(wa["truth"], float)
                        rows.setdefault("win", []).append(wa)
                        rows.setdefault("mt", []).append([(sc, h, np.asarray(p)) for sc, h, p in we["rec"]])
                        rows.setdefault("qm", []).append(qm_sel(wa["qm_events"]["QM-3C+SC"]))
                        rows.setdefault("tr", []).append(tr)
                for al in A.ALPHAS:
                    parts = []
                    for m, key in (("MT-FW+SC", "mt"), ("QM-3C+SC", "qm")):
                        q = A.conformal_q(mx[m], al)[0]
                        d = []
                        for ev, tr in zip(rows[key], rows["tr"]):
                            for x in ev:
                                sc, p = x[0], x[-1]
                                if sc <= q:
                                    continue
                                dd = float(np.min(np.linalg.norm(tr[:, :2] - np.asarray(p)[:2], axis=1)))
                                if dd > 5.0:
                                    d.append(dd)
                        d = np.asarray(d)
                        nb = [int(np.sum((d > lo) & (d < hi))) if np.isfinite(hi) else int(np.sum(d >= lo)) for lo, hi in BINS]
                        nb[0] = int(np.sum((d > 5.0) & (d < 10.0)))
                        parts.append(f"{m} {len(d) / len(rows['tr']):.2f} por ventana ({len(d)}: " + "/".join(map(str, nb)) + ")")
                    L.append(f"   {reg:6s} {s} {A.LAB[cell]:12s} alfa {al:<4} ventanas con eventos {len(rows['tr'])}: "
                             + " | ".join(parts))
    L.append(f"   control: registros de MT-FW+SC identicos en {tot_same}/{tot_n} ventanas con eventos")
    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(HERE, "task2d_explore_e1.txt"), "w") as fh:
        fh.write(txt + "\n")


if __name__ == "__main__":
    main()
