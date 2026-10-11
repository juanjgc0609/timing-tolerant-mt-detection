"""task2i_control.py - controles de la Tarea 2I en la semilla 999 (datos del piloto C8; ninguna semilla 5151/6262).

(1) Reproduccion del piloto C8 (explore_contrib/out/c8_DENSE-IID_*.json y c8_STD-IID_*.json; redes 60, 67, 74; TAG 3801,
    5 ventanas, kmax 7; TAG 3802, 4 ventanas en el piloto, kmax 5): los registros y las posiciones de task2i_run.py
    (partes cprim_*) deben ser identicos (==) a los guardados para BOMP-PK, MT-FW, MT-FW-NOCAP y HUBER-N.
(2) Umbrales secundarios (partes cthr_*, redes 60, 67, TAG 3801): HUBER-N, HUBER-N90 y HUBER-N99 identicos a
    c8_DENSE-IID-HSWEEP_*.json y HUBER-C5 a c8_DENSE-IID-HWIDE_*.json.
(3) sd de sitio 0.25 y 0.8 (partes cs25_*, cs80_*; red 80; TAG 3805 / 3806; 4 ventanas en el piloto): los 4 metodos
    identicos a c8_DENSE-IID-S25_*.json / S80.
(4) MT-FW tope 6 (y sin tope) por c8_lib.mtfw_general frente a c6_lib.mtfw_dense (red 60, 3 celdas x 2 regimenes, 5
    ventanas densas y 5 estandar regeneradas): aciertos y posiciones identicos, puntajes con diferencia relativa
    <= 1e-9 (las bases Q de xc.Ops y c6_lib.ModelOps difieren en el ultimo bit; c8_validation.txt, V2-V3).
(5) analyze_2f.Pool (FROC exacta de analyze_2i.py) frente a xc.recall_caps sobre TODO el piloto DENSE-IID y STD-IID, y
    frente a los valores publicados en REPORT_C8.txt.
(6) El bootstrap de analyze_2i.py (estratos regimen x celda) sobre todo el piloto DENSE-IID y STD-IID frente a los
    intervalos y DE publicados en REPORT_C8.txt (otro generador de remuestras: acuerdo aproximado).
(7) Estimacion de tiempo del diseno completo y regla de presupuesto.

    OPENBLAS_CORETYPE=Prescott python task2i_control.py [--out task2i_control.txt]
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

import analyze_2i as A  # noqa: E402
import task2i_run as R  # noqa: E402
import c2_run as C2  # noqa: E402
import c6_lib as L6  # noqa: E402
import c8_lib as C  # noqa: E402
import xc  # noqa: E402

OUT = os.path.join(EXPL, "out")
PN = {"BOMP-PK+SC": "BOMP-PK", "MT-FW+SC": "MT-FW", "MT-FW-NOCAP+SC": "MT-FW-NOCAP", "HUBER-N+SC": "HUBER-N",
      "HUBER-N90+SC": "HUBER-N90", "HUBER-N99+SC": "HUBER-N99", "HUBER-C5+SC": "HUBER-C5"}
CELLS, REGIMES = A.CELLS, A.REGIMES


def pilot(design):
    rows = {}
    for reg in REGIMES:
        d = json.load(open(os.path.join(OUT, f"c8_{design}_{reg}.json")))
        for r in d["rows"]:
            rows[(reg, r["cell"], r["net"], r["w"])] = r
    return rows


def part(arm, reg, cell, j):
    with open(R.part_path(arm, 999, reg, cell, j)) as fh:
        return json.load(fh)


def compare(arm, nets, kind, pil, methods, nwin):
    """numero de ventanas comparadas, y por metodo cuantas tienen registro y posiciones identicos."""
    n = kok = 0
    same = {m: 0 for m in methods}
    for reg in REGIMES:
        for cell in CELLS:
            for j in nets:
                e = part(arm, reg, cell, j)
                for w in e["windows"][kind][:nwin]:
                    p = pil[(reg, cell, j, w["w"])]
                    n += 1
                    kok += int(w["k"] == p["k"] and w["truth"] == p["truth"] and w["taus"] == p["taus"])
                    for m in methods:
                        same[m] += int(w["records"][m] == [list(x) for x in p["rec"][PN[m]]]
                                       and w["loc"][m] == p["loc"][PN[m]])
    return n, kok, same


def pilot_as_data(design):
    """el piloto completo en la estructura de analyze_2i (estratos regimen x celda, ventanas por red y ventana)."""
    pil = pilot(design)
    nets = sorted({k[2] for k in pil})
    nw = max(k[3] for k in pil) + 1
    D = {}
    for reg in REGIMES:
        for cell in CELLS:
            W = {"k": [], "rec": {m: [] for m in A.PRIMARY}}
            for j in nets:
                for w in range(nw):
                    r = pil[(reg, cell, j, w)]
                    W["k"].append(r["k"])
                    for m in A.PRIMARY:
                        W["rec"][m].append([list(x) for x in r["rec"][PN[m]]])
            W["k"] = np.asarray(W["k"], float)
            D[(reg, cell)] = W
    return D, len(nets), nw


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "task2i_control.txt"))
    a = ap.parse_args()
    R.check_kernel()
    out = []
    P = out.append
    P("Controles de la Tarea 2I en la semilla 999 (task2i_control.py)")
    P(f"Nucleo de este proceso: {R.kernel_info()}")
    allok = True

    # ------------------------------------------------------------------ (1)
    P("")
    P("(1) Reproduccion del piloto C8, 4 metodos primarios (redes 60, 67, 74)")
    nets = R.CONTROL["cprim"]["nets"]
    for kind, design, nwin in (("dense", "DENSE-IID", 5), ("std", "STD-IID", 4)):
        n, kok, same = compare("cprim", nets, kind, pilot(design), A.PRIMARY, nwin)
        ok = n == len(nets) * 6 * nwin and kok == n and all(v == n for v in same.values())
        allok &= ok
        P(f"  {design}: {n} ventanas; k, verdad y origenes iguales {kok}/{n}; registro y posiciones identicos (==): " +
          ", ".join(f"{m} {v}/{n}" for m, v in same.items()) + f" -> {'OK' if ok else 'DIFERENCIAS'}")

    # ------------------------------------------------------------------ (2)
    P("")
    P("(2) Umbrales secundarios (redes 60, 67; ventanas densas TAG 3801)")
    nets = R.CONTROL["cthr"]["nets"]
    for design, ms in (("DENSE-IID-HSWEEP", ("HUBER-N+SC", "HUBER-N90+SC", "HUBER-N99+SC")),
                       ("DENSE-IID-HWIDE", ("HUBER-N+SC", "HUBER-C5+SC"))):
        n, kok, same = compare("cthr", nets, "dense", pilot(design), ms, 5)
        ok = n == len(nets) * 30 and kok == n and all(v == n for v in same.values())
        allok &= ok
        P(f"  {design}: {n} ventanas; verdad igual {kok}/{n}; identicos: " +
          ", ".join(f"{m} {v}/{n}" for m, v in same.items()) + f" -> {'OK' if ok else 'DIFERENCIAS'}")
    n, kok, same = compare("cthr", nets, "std", pilot("STD-IID"), ("HUBER-N+SC",), 4)
    ok = kok == n and same["HUBER-N+SC"] == n
    allok &= ok
    P(f"  STD-IID (HUBER-N del brazo cthr, TAG 3802): identicos {same['HUBER-N+SC']}/{n} -> {'OK' if ok else 'DIFERENCIAS'}")

    # ------------------------------------------------------------------ (3)
    P("")
    P("(3) sd de sitio 0.25 y 0.8 (red 80; TAG 3805 / 3806; 4 ventanas del piloto)")
    for arm, design in (("cs25", "DENSE-IID-S25"), ("cs80", "DENSE-IID-S80")):
        n, kok, same = compare(arm, R.CONTROL[arm]["nets"], "dense", pilot(design), A.PRIMARY, 4)
        ok = n == 24 and kok == n and all(v == n for v in same.values())
        allok &= ok
        lsd = np.mean([w["site_logsd"] for reg in REGIMES for cell in CELLS for w in
                       part(arm, reg, cell, 80)["windows"]["dense"]])
        P(f"  {design}: {n} ventanas; verdad igual {kok}/{n}; identicos: " +
          ", ".join(f"{m} {v}/{n}" for m, v in same.items()) + f"; sd media del log de sitio {lsd:.3f} -> "
          f"{'OK' if ok else 'DIFERENCIAS'}")

    # ------------------------------------------------------------------ (4)
    P("")
    P("(4) c8_lib.mtfw_general('ls') frente a c6_lib.mtfw_dense (red 60; ventanas regeneradas)")
    C.SEED = 999
    nwin = eq_cap = eq_nc = eq_part = 0
    mrel = 0.0
    for reg in REGIMES:
        for cell in CELLS:
            net = C.WNet(cell, 60, reg, C.World("iid", 0.5))
            grid = xc.make_grid()
            dops = C.DOps(xc.Ops(net.STA, net.corr, grid))
            mops = L6.ModelOps(net.STA, grid, L6.Homog(), net.corr)
            e = part("cprim", reg, cell, 60)
            for kind, tag, draw, kmax in (("dense", 3801, C2.draw_dense, 7), ("std", 3802, xc.draw_sources, 5)):
                st = net.stream(tag, draw)
                for wi in range(5):
                    w = next(st)
                    nwin += 1
                    for cap, capv, name in ((True, 6.0, "MT-FW+SC"), (False, None, "MT-FW-NOCAP+SC")):
                        r6, l6 = L6.mtfw_dense(mops, w["MFn"], w["truth"], kmax, cap=cap)
                        r8, l8 = C.mtfw_general(dops, w["MFn"], w["truth"], kmax, "ls", capv)
                        same = l6 == l8 and len(r6) == len(r8) and all(h1 == h2 for (_s1, h1), (_s2, h2) in zip(r6, r8))
                        if same and r6:
                            mrel = max(mrel, max(abs(s1 - s2) / abs(s1) for (s1, _h), (s2, _h2) in zip(r6, r8)))
                        if cap:
                            eq_cap += int(same)
                        else:
                            eq_nc += int(same)
                        eq_part += int([[float(s), int(h)] for s, h in r8] == e["windows"][kind][wi]["records"][name])
    ok = eq_cap == nwin and eq_nc == nwin and mrel <= 1e-9 and eq_part == 2 * nwin
    allok &= ok
    P(f"  {nwin} ventanas: aciertos y posiciones (nodo y origen) identicos con tope 6 {eq_cap}/{nwin}, sin tope "
      f"{eq_nc}/{nwin}; diferencia relativa maxima de puntajes {mrel:.1e}; registros regenerados aqui == partes cprim "
      f"{eq_part}/{2 * nwin} -> {'OK' if ok else 'FALLA'}")

    # ------------------------------------------------------------------ (5) + (6)
    P("")
    P("(5) analyze_2f.Pool frente a xc.recall_caps, piloto completo (20 redes por celda y regimen), agrupado")
    pub = {"DENSE-IID": {"BOMP-PK+SC": (0.104, 0.200), "MT-FW+SC": (0.154, 0.242), "MT-FW-NOCAP+SC": (0.104, 0.197),
                         "HUBER-N+SC": (0.172, 0.264)},
           "STD-IID": {"BOMP-PK+SC": (0.397, 0.511), "MT-FW+SC": (0.359, 0.470), "MT-FW-NOCAP+SC": (0.397, 0.511),
                       "HUBER-N+SC": (0.421, 0.497)}}
    pubd = {"DENSE-IID": {("HUBER-N+SC", "BOMP-PK+SC"): "+0.064 [+0.028, +0.105] DE 0.0194",
                          ("HUBER-N+SC", "MT-FW+SC"): "+0.022 [-0.014, +0.055] DE 0.0178",
                          ("MT-FW+SC", "BOMP-PK+SC"): "+0.042 [+0.004, +0.083]"},
            "STD-IID": {("HUBER-N+SC", "BOMP-PK+SC"): "-0.014 [-0.050, +0.042], IC 90 [-0.042, +0.035], DE 0.0233",
                        ("MT-FW+SC", "BOMP-PK+SC"): "-0.041 [-0.086, +0.029], IC 90 [-0.078, +0.020]"}}
    lines6 = []
    for design in ("DENSE-IID", "STD-IID"):
        D, nnet, nw = pilot_as_data(design)
        ks = np.concatenate([D[k]["k"] for k in A.ALL])
        old = A.N_WIN
        A.N_WIN = nw
        CC = A.counts(999, "dense" if design.startswith("DENSE") else "std", nnet)
        e = A.evaluate(D, CC, A.ALL, A.PRIMARY)
        A.N_WIN = old
        for m in A.PRIMARY:
            recs = [[tuple(x) for x in r] for k in A.ALL for r in D[k]["rec"][m]]
            rx = xc.recall_caps(recs, int(ks.sum()))[0]
            eq = abs(e[0][m][0.1] - rx[0]) < 1e-12 and abs(e[0][m][0.25] - rx[1]) < 1e-12
            eqp = abs(round(rx[0], 3) - pub[design][m][0]) < 5e-4 and abs(round(rx[1], 3) - pub[design][m][1]) < 5e-4
            allok &= eq and eqp
            P(f"  {design:9s} {m:15s} Pool {e[0][m][0.1]:.3f} / {e[0][m][0.25]:.3f}; xc.recall_caps {rx[0]:.3f} / "
              f"{rx[1]:.3f} {'igual' if eq else 'DISTINTO'}; publicado {pub[design][m][0]:.3f} / {pub[design][m][1]:.3f} "
              f"{'igual' if eqp else 'DISTINTO'}")
        for (x, y), txt in pubd[design].items():
            d = A.diff(e, x, y, 0.25)
            lines6.append(f"  {design:9s} {x} - {y}: {A.fci(d)}, IC 90 [{d['lo90']:+.3f}, {d['hi90']:+.3f}], DE "
                          f"{d['sd']:.4f}  |  publicado {txt}")
    P("")
    P("(6) Bootstrap de analyze_2i.py (2000 replicas) sobre el piloto completo frente a REPORT_C8.txt (1000 replicas, "
      "otro generador)")
    out.extend(lines6)

    # ------------------------------------------------------------------ (7)
    P("")
    P("(7) Tiempo (partes de control, 2 procesos en paralelo) y regla de presupuesto")
    def tmean(arm):
        nets_ = R.CONTROL[arm]["nets"]
        return float(np.mean([part(arm, reg, cell, j)["t_total"] for reg in REGIMES for cell in CELLS for j in nets_]))
    tp, tt, ts = tmean("cprim"), tmean("cthr"), 0.5 * (tmean("ds25") + tmean("ds80"))
    nprim = len(R.ARMS["prim"]["nets"]) * 6 * 2
    nsec = len(R.ARMS["thr"]["nets"]) * 6 * 2
    cpu_p = nprim * tp
    cpu_s = nsec * tt + 2 * nsec * ts
    P(f"  parte del brazo prim (red + 5 densas + 5 estandar, 4 metodos): {tp:.1f} s; thr: {tt:.1f} s; s25/s80: {ts:.1f} s")
    P(f"  primario: {nprim} partes en las dos semillas -> {cpu_p / 3600:.2f} h de CPU, {cpu_p / 7200:.2f} h en 2 nucleos")
    P(f"  secundarios (thr, s25, s80): {3 * nsec} partes -> {cpu_s / 3600:.2f} h de CPU, {cpu_s / 7200:.2f} h en 2 nucleos")
    tot = (cpu_p + cpu_s) / 7200
    P(f"  total estimado {tot:.2f} h en 2 nucleos (mas unos minutos de analisis). Regla de presupuesto (8 h): "
      f"{'SE ACTIVA: se omiten los brazos secundarios' if tot > 8 else 'NO se activa: se corren los brazos secundarios'}")
    kern = sorted({json.dumps(part(arm, reg, cell, j)["kernel"]) for arm in R.CONTROL for reg in REGIMES
                   for cell in CELLS for j in R.CONTROL[arm]["nets"]})
    P("")
    P(f"Nucleos OpenBLAS registrados en las partes de control: {kern}")
    P(f"RESULTADO GLOBAL DE LOS CONTROLES: {'OK' if allok else 'FALLA'}")
    txt = "\n".join(out) + "\n"
    with open(a.out, "w") as fh:
        fh.write(txt)
    print(txt)


if __name__ == "__main__":
    main()
