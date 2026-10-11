"""task2f_control.py - controles de la Tarea 2F en la semilla 999 (datos del piloto; ninguna semilla 4242/8888).

(1) Reproduccion del piloto denso (explore_contrib/c2_ablate.py; out/c2_ablate_thin.json y c2_ablate_strong.json;
    redes 10-15, TAG 3002, draw_dense, kmax 7): los registros de task2f_run.py (partes cdense_*) deben ser identicos,
    ventana a ventana, a los guardados para MT-FW, MT-FW-NOCAP, BOMP-PK, BOMP-PK-CAP y BOMP-PK-SEQ; el recall agrupado
    por regimen (analyze_2f.Pool) debe igualar al del piloto (xc.recall_caps sobre los registros guardados).
(2) MT-FW-NOCAP+SC y BOMP-PK-SEQ+SC (BOMP-PK con deflacion secuencial) con secuencias de aciertos identicas.
(3) Brazo estandar: registros de MT-FW:HOM y BOMP-PK:HOM de explore_contrib C6 (out/c6_*.json, redes 20-21, TAG 3006,
    draw_sources, kmax 5) reproducidos por task2f_run.py (partes cstd_*).
(4) analyze_2f.Pool frente a froc_exact.exact_HF + analyze_2c.val (la FROC exacta y el recall de 2C-2E) en la
    estimacion puntual y en 200 remuestras jerarquicas.
(5) Estimacion del tiempo del diseno completo a partir de los tiempos de las partes de control.

    python task2f_control.py [--out task2f_control.txt]
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXPL = os.path.join(HERE, "explore_contrib")
for p in (EXPL, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

import analyze_2f as A  # noqa: E402
from analyze_2c import HF, val  # noqa: E402

PARTS = A.PARTS
CELLS, REGIMES = A.CELLS, A.REGIMES
NAME = {"MT-FW+SC": "MT-FW", "MT-FW-NOCAP+SC": "MT-FW-NOCAP", "BOMP-PK+SC": "BOMP-PK", "BOMP-PK-CAP+SC": "BOMP-PK-CAP",
        "BOMP-PK-SEQ+SC": "BOMP-PK-SEQ"}
NAME6 = {"MT-FW+SC": "MT-FW:HOM", "BOMP-PK+SC": "BOMP-PK:HOM"}


def recall_caps_xc(recs, nte):
    import xc
    return xc.recall_caps([[tuple(x) for x in r] for r in recs], nte)[0]


def load_parts(prefix, regime, nets):
    rows = []
    for cell in CELLS:
        for j in nets:
            with open(os.path.join(PARTS, f"{prefix}_seed999_{regime}_{cell}_net{j}.json")) as fh:
                e = json.load(fh)
            for w in e["windows"]:
                rows.append({"cell": cell, "net": j, "w": w["w"], "k": w["k"], "truth": w["truth"], "taus": w["taus"],
                             "rec": w["records"], "sec": w["sec"], "t_setup": e["t_setup"], "kernel": e["kernel"]})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "task2f_control.txt"))
    a = ap.parse_args()
    L = []
    P = L.append
    P("Controles de la Tarea 2F en la semilla 999 (task2f_control.py)")
    allok = True
    secs = {"dense": [], "std": []}
    kernels = set()

    # ------------------------------------------------------------------ (1) + (2)
    P("")
    P("(1) Reproduccion del piloto denso (c2_ablate; redes 10-15, TAG 3002, kmax 7)")
    for reg in REGIMES:
        pil = json.load(open(os.path.join(EXPL, "out", f"c2_ablate_{reg}.json")))
        pil = {(r["cell"], r["net"], r["w"]): r for r in pil}
        mine = load_parts("cdense", reg, range(10, 16))
        same = {m: 0 for m in NAME}
        kok = 0
        seq_eq = 0
        for r in mine:
            p = pil[(r["cell"], r["net"], r["w"])]
            kok += r["k"] == p["k"]
            for m, pm in NAME.items():
                same[m] += r["rec"][m] == p["rec"][pm]
            seq_eq += [h for _s, h in r["rec"]["MT-FW-NOCAP+SC"]] == [h for _s, h in r["rec"]["BOMP-PK-SEQ+SC"]]
            kernels.add(json.dumps(r["kernel"]))
            secs["dense"].append((sum(r["sec"][m] for m in A.METHODS), r["t_setup"] / A.N_WIN))
        n = len(mine)
        ok = n == len(pil) == 180 and kok == n and all(v == n for v in same.values())
        allok &= ok and seq_eq == n
        P(f"  {reg}: {n} ventanas (piloto {len(pil)}); k igual {kok}/{n}; registros identicos (==, bit a bit tras JSON): "
          + ", ".join(f"{m} {v}/{n}" for m, v in same.items()) + f" -> {'OK' if ok else 'DIFERENCIAS'}")
        P(f"  (2) {reg}: MT-FW-NOCAP+SC y BOMP-PK-SEQ+SC con secuencias de aciertos identicas en {seq_eq}/{n} ventanas")
        ks = [r["k"] for r in mine]
        nte = sum(ks)
        P(f"      recall agrupado (3 celdas) tope 0.1 / 0.25: analyze_2f.Pool sobre task2f  |  xc.recall_caps sobre el "
          f"piloto guardado")
        for m, pm in NAME.items():
            pool = A.Pool([r["rec"][m] for r in mine], ks)
            v = pool.recall(np.ones((1, n)))
            prow = [pil[(r["cell"], r["net"], r["w"])]["rec"][pm] for r in mine]
            rp = recall_caps_xc(prow, sum(pil[(r["cell"], r["net"], r["w"])]["k"] for r in mine))
            eq = abs(v[0.1][0] - rp[0]) < 1e-12 and abs(v[0.25][0] - rp[1]) < 1e-12
            allok &= eq
            P(f"      {m:16s} {v[0.1][0]:.3f} / {v[0.25][0]:.3f}  |  {pm:12s} {rp[0]:.3f} / {rp[1]:.3f}  "
              f"{'igual' if eq else 'DISTINTO'}")
    P("  Valores publicados en explore_contrib/out/c2_ablate_analysis.txt (thin) y c2_ablate_strong_analysis.txt:")
    P("    thin   BOMP-PK 0.135/0.230, MT-FW 0.212/0.304, MT-FW-NOCAP 0.131/0.233, BOMP-PK-CAP 0.180/0.291, "
      "BOMP-PK-SEQ 0.131/0.233")
    P("    strong BOMP-PK 0.085/0.190, MT-FW 0.164/0.241, MT-FW-NOCAP 0.085/0.183, BOMP-PK-CAP 0.144/0.215, "
      "BOMP-PK-SEQ 0.085/0.183")

    # ------------------------------------------------------------------ (3)
    P("")
    P("(3) Brazo estandar: C6 HOM (redes 20-21, TAG 3006, draw_sources, kmax 5)")
    for reg in REGIMES:
        c6 = json.load(open(os.path.join(EXPL, "out", f"c6_{reg}.json")))["rows"]
        c6 = {(r["cell"], r["net"], r["w"]): r for r in c6}
        mine = load_parts("cstd", reg, (20, 21))
        same = {m: 0 for m in NAME6}
        kok = 0
        for r in mine:
            p = c6[(r["cell"], r["net"], r["w"])]
            kok += (r["k"] == p["k"]) and (r["truth"] == p["truth"]) and (r["taus"] == p["taus"])
            for m, pm in NAME6.items():
                same[m] += r["rec"][m] == p["rec"][pm]
            kernels.add(json.dumps(r["kernel"]))
            secs["std"].append((sum(r["sec"][m] for m in A.METHODS), r["t_setup"] / A.N_WIN))
        n = len(mine)
        ok = n == 60 and kok == n and all(v == n for v in same.values())
        allok &= ok
        P(f"  {reg}: {n} ventanas; k, verdad y origenes iguales {kok}/{n}; registros identicos: "
          + ", ".join(f"{m} {v}/{n}" for m, v in same.items()) + f" -> {'OK' if ok else 'DIFERENCIAS'}")

    # ------------------------------------------------------------------ (4)
    P("")
    P("(4) analyze_2f.Pool frente a froc_exact.exact_HF + analyze_2c.val")
    worst = 0.0
    for reg in REGIMES:
        mine = load_parts("cdense", reg, range(10, 16))
        ks = np.array([r["k"] for r in mine], float)
        nets = len(range(10, 16))
        CC = A.counts(999, "dense", nets, B_=200)
        C = np.concatenate([CC[(reg, c)] for c in CELLS], axis=1)
        C = np.vstack([np.ones((1, len(mine))), C])
        for m in A.METHODS:
            recs = [r["rec"][m] for r in mine]
            v1 = A.Pool(recs, ks).recall(C)
            v2 = val(HF(recs), C, ks)
            for cap in A.CAPS:
                worst = max(worst, float(np.abs(v1[cap] - v2[str(cap)]).max()))
    allok &= worst < 1e-12
    P(f"  punto y 200 remuestras jerarquicas, 4 metodos, 2 regimenes, 2 topes: diferencia maxima {worst:.1e}")

    # ------------------------------------------------------------------ (5)
    P("")
    P("(5) Estimacion de tiempo del diseno completo (tiempos de las partes de control, 2 procesos en paralelo)")
    tot = 0.0
    for arm, nets in (("dense", 96), ("std", 16)):
        s = np.array(secs[arm])
        per_win = float(s[:, 0].mean() + s[:, 1].mean())
        nwin = nets * 3 * 2 * A.N_WIN * 2
        tot += nwin * per_win
        P(f"  {arm}: {per_win:.2f} s de CPU por ventana (4 metodos {s[:, 0].mean():.2f} + preparacion de la red "
          f"repartida {s[:, 1].mean():.2f}); {nwin} ventanas en las dos semillas -> {nwin * per_win / 3600:.2f} h de CPU")
    P(f"  total {tot / 3600:.2f} h de CPU, unas {tot / 3600 / 2:.2f} h de pared en 2 nucleos")
    P("")
    P(f"Nucleos OpenBLAS registrados en las partes de control: {sorted(kernels)}")
    P(f"RESULTADO GLOBAL DE LOS CONTROLES: {'OK' if allok else 'FALLA'}")
    txt = "\n".join(L) + "\n"
    with open(a.out, "w") as fh:
        fh.write(txt)
    print(txt)


if __name__ == "__main__":
    main()
