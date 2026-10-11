"""task2l_control.py - controles de la Tarea 2L en la semilla 999 (datos del piloto C11; ninguna semilla 7171/9393).

(1) Reproduccion del piloto C11 (explore_contrib/out/c11_eval_{thin,strong}.json; redes 210, 219, 228 del conjunto de
    evaluacion; TAG 4101 estandar y 4102 densas, 5 ventanas): los atomos guardados por task2l_run.py (partes ctrl_*)
    deben ser identicos (==) a los del piloto para MT@S0, MT@S0+T3, MT@S2+T3, MT@ORA+T3 (estandar y densas) y BOMP@S0,
    BOMP@S0+T3, BOMP@S2+T3, BOMP@ORA+T3 (estandar); tambien la verdad, los terminos S0, S2 y ORA y los niveles de T3.
(2) S0 = +SC registrado: (a) en las redes de control, task2a_confirm.catalogue_real +
    station_corrections.estimate_station_corrections (la construccion registrada en 2A) frente al S0 de las partes;
    (b) xc.NetCfg(cell, j, regimen).corr para j = 0..4 frente a las correcciones guardadas por la corrida de
    calibracion de 2A en la semilla 999 (task2a_seed999_{thin,strong}.json).
(3) Mecanismo de la semilla: con task2l_run.set_seed(s) para una semilla de prueba s = 31337 (no registrada ni usada),
    la red y los retardos de xc.NetCfg son task2a_confirm.network(s, ...), distintos de los de 999, el catalogo es
    task2a_confirm.catalogue_real con los flujos [s, 2003, ci, j, ri] y la primera ventana de c9_run.stream es la de
    c8_lib.generate_site con los flujos [s, TAG, ci, j, ri, 41|233|61].
(4) analyze_2l.py (c10_lib.label / record, analyze_2f.Pool, bootstrap por estratos) sobre TODO el conjunto de evaluacion
    de C11 (30 redes por celda y regimen): recall identico al de xc.recall_caps y a los valores publicados en
    REPORT_C11.txt; intervalos frente a los publicados (otro generador de remuestras: acuerdo aproximado).
(5) Estimacion de tiempo del diseno completo y regla de presupuesto (10 h en 2 nucleos).

    OPENBLAS_CORETYPE=Prescott python task2l_control.py [--out task2l_control.txt]
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

import analyze_2l as A  # noqa: E402
import task2l_run as R  # noqa: E402
import c8_lib as C8  # noqa: E402
import c9_run as R9  # noqa: E402
import xc  # noqa: E402
import task2a_confirm as T  # noqa: E402
from station_corrections import estimate_station_corrections  # noqa: E402

OUT = os.path.join(EXPL, "out")
CELLS, REGIMES = A.CELLS, A.REGIMES
TEST_SEED = 31337


def pilot():
    rows, nets = {}, {}
    for reg in REGIMES:
        d = json.load(open(os.path.join(OUT, f"c11_eval_{reg}.json")))
        for r in d["rows"]:
            rows[(reg, r["design"], r["cell"], r["net"], r["w"])] = r
        for n in d["nets"]:
            nets[(reg, n["cell"], n["net"])] = n
    return rows, nets


def part(reg, cell, j):
    with open(R.part_path("ctrl", 999, reg, cell, j)) as fh:
        return json.load(fh)


def same(a, b):
    return json.dumps(a) == json.dumps(b)               # igualdad exacta, con NaN (sref ausente) incluido


def pilot_as_data(rows, design):
    """el piloto completo en la estructura de analyze_2l.load."""
    nets = sorted({k[3] for k in rows})
    Dk = {}
    for reg in REGIMES:
        for cell in CELLS:
            W = {"k": [], "truth": [], "taus": [], "atoms": {}, "sec": {}}
            for j in nets:
                for w in range(A.N_WIN):
                    r = rows[(reg, design, cell, j, w)]
                    W["k"].append(r["k"])
                    W["truth"].append(r["truth"])
                    W["taus"].append(r["taus"])
                    for m, a in r["atoms"].items():
                        W["atoms"].setdefault(m, []).append(a)
            W["k"] = np.asarray(W["k"], float)
            Dk[(reg, cell)] = W
    return Dk, len(nets)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "task2l_control.txt"))
    a = ap.parse_args()
    R.check_kernel()
    out = []
    P = out.append
    P("Controles de la Tarea 2L en la semilla 999 (task2l_control.py)")
    P(f"Nucleo de este proceso: {R.kernel_info()}")
    allok = True
    rows, pnets = pilot()

    # ------------------------------------------------------------------ (1)
    P("")
    P(f"(1) Reproduccion del piloto C11 (redes {R.CONTROL_NETS} x 3 celdas x 2 regimenes = 18 partes, "
      f"{len(R.CONTROL_NETS) * 3} redes distintas)")
    cnt, tot = {}, {}

    def add(k, v):
        cnt[k] = cnt.get(k, 0) + int(v)
        tot[k] = tot.get(k, 0) + 1
    for reg in REGIMES:
        for cell in CELLS:
            for j in R.CONTROL_NETS:
                e = part(reg, cell, j)
                ni = pnets[(reg, cell, j)]
                for k in ("S0", "S2", "ORA"):
                    add("terminos " + k, same(e["terms"][k], ni["terms"][k]))
                    add("nivel T3 " + k, same(e["sigT3"][k], ni["info"]["sigT3"][k]))
                add("error rms S0, S2", all(e["err"][k] == ni["info"]["err"][k] for k in ("S0", "S2")))
                for kind, des in (("std", "STD"), ("dense", "DENSE")):
                    for w in e["windows"][kind]:
                        p_ = rows[(reg, des, cell, j, w["w"])]
                        add(f"{kind} verdad y origenes", w["truth"] == p_["truth"] and w["taus"] == p_["taus"]
                            and w["k"] == p_["k"])
                        for m, at in w["atoms"].items():
                            add(f"{kind} {m}", same(at, p_["atoms"][m]))
    ok = all(cnt[k] == tot[k] for k in cnt)
    allok &= ok
    for k in cnt:
        P(f"  {k:24s} identicos (==) {cnt[k]}/{tot[k]}")
    P(f"  -> {'OK' if ok else 'DIFERENCIAS'}")

    # ------------------------------------------------------------------ (2)
    P("")
    P("(2) S0 = +SC registrado")
    R.set_seed(999)
    n = eq = 0
    for reg in REGIMES:
        for cell in CELLS:
            for j in R.CONTROL_NETS:
                ci, _kind, _nst, nclut = T.net_setup(999, cell)
                STA, statics = T.network(999, cell, j)
                cat, _drawn, _truth = T.catalogue_real(STA, T.K_CAT, statics, [999, T.TAG_CAT, ci, j, REGIMES.index(reg)],
                                                       nclut, T.CLUT_AMP[reg])
                corr, _nacc = estimate_station_corrections(STA, cat)
                n += 1
                eq += int(same(np.asarray(corr).tolist(), part(reg, cell, j)["terms"]["S0"]))
    ok = eq == n
    allok &= ok
    P(f"  (a) task2a_confirm.catalogue_real + estimate_station_corrections frente al S0 de las partes de control: "
      f"identicos (==) {eq}/{n} -> {'OK' if ok else 'DIFERENCIAS'}")
    n = eq = 0
    mx = 0.0
    for reg in REGIMES:
        d2a = json.load(open(os.path.join(HERE, f"task2a_seed999_{reg}.json")))
        for cellres in d2a:
            cell = cellres["scenario"]
            for nw in cellres["networks"]:
                nc = xc.NetCfg(cell, nw["net"], reg)
                n += 1
                eq += int(same(np.asarray(nc.corr).tolist(), nw["corr"]))
                mx = max(mx, float(np.abs(np.asarray(nc.corr) - np.asarray(nw["corr"])).max()))
    ok = eq == n
    allok &= ok
    P(f"  (b) xc.NetCfg.corr (j = 0..4) frente a las correcciones +SC guardadas por 2A en la semilla 999 "
      f"(task2a_seed999_*.json): identicos (==) {eq}/{n}; diferencia maxima {mx:.1e} s -> {'OK' if ok else 'DIFERENCIAS'}")

    # ------------------------------------------------------------------ (3)
    P("")
    P(f"(3) Mecanismo de la semilla (semilla de prueba {TEST_SEED}; ninguna semilla registrada)")
    R.set_seed(TEST_SEED)
    oks = []
    for cell, reg, j in (("ladder_4", "thin", 0), ("unilateral_c40", "strong", 3)):
        nc = xc.NetCfg(cell, j, reg)
        ci, _kind, _nst, nclut = T.net_setup(TEST_SEED, cell)
        ri = REGIMES.index(reg)
        STA, statics = T.network(TEST_SEED, cell, j)
        STA9, stat9 = T.network(999, cell, j)
        cat, _drawn, _truth = T.catalogue_real(STA, T.K_CAT, statics, [TEST_SEED, T.TAG_CAT, ci, j, ri], nclut,
                                               T.CLUT_AMP[reg])
        c1 = np.array_equal(nc.STA, STA) and np.array_equal(nc.statics, statics)
        c2 = not np.array_equal(STA, STA9) and not np.array_equal(statics, stat9)
        c3 = all(np.array_equal(x[0], y[0]) and np.array_equal(x[1], y[1]) and x[2] == y[2] for x, y in zip(nc.cat, cat))
        c4 = True
        for tag, draw in ((R.TAG_STD, xc.draw_sources), (R.TAG_DENSE, R.C2.draw_dense)):
            w = next(R9.stream(nc, tag, draw, None))
            b = [TEST_SEED, tag, ci, j, ri]
            rng, rng_m, rng_h = (np.random.default_rng(b + [41]), np.random.default_rng(b + [233]),
                                 np.random.default_rng(b + [61]))
            srcs, taus = draw(rng)
            D, _pk, _site, _amp = C8.generate_site(STA, srcs, taus, nclut, rng, rng_m, rng_h, T.CLUT_AMP[reg], cfg=None,
                                                   statics=statics)
            c4 &= np.array_equal(D, w["D"]) and [float(t) for t in taus] == w["taus"]
        oks.append(c1 and c2 and c3 and c4)
        P(f"  {cell} {reg} red {j}: red y retardos = task2a_confirm.network({TEST_SEED}) {c1}; distintos de 999 {c2}; "
          f"catalogo = catalogue_real([{TEST_SEED}, 2003, ...]) {c3}; primera ventana estandar y densa = flujos "
          f"[{TEST_SEED}, TAG, ...] {c4}")
    R.set_seed(999)
    ok = all(oks)
    allok &= ok
    P(f"  -> {'OK' if ok else 'FALLA'}")

    # ------------------------------------------------------------------ (4)
    P("")
    P("(4) analyze_2l.py sobre todo el conjunto de evaluacion de C11 (30 redes por celda y regimen) frente a "
      "xc.recall_caps y a REPORT_C11.txt (2000 replicas aqui; el piloto uso otro generador de remuestras)")
    pub_rec = {("std", "LOC+T"): {"MT@S0": (0.256, 0.347), "MT@S0+T3": (0.256, 0.371), "MT@S2+T3": (0.323, 0.441),
                                  "MT@ORA+T3": (0.374, 0.462), "BOMP@S0": (0.238, 0.356), "BOMP@S0+T3": (0.256, 0.398),
                                  "BOMP@S2+T3": (0.284, 0.455), "BOMP@ORA+T3": (0.385, 0.508)},
               ("dense", "LOC+T"): {"MT@S0": (0.055, 0.121), "MT@S0+T3": (0.048, 0.112), "MT@S2+T3": (0.092, 0.176),
                                    "MT@ORA+T3": (0.114, 0.189)},
               ("dense", "LOC"): {"MT@S0": (0.168, 0.272), "MT@S0+T3": (0.173, 0.272), "MT@S2+T3": (0.224, 0.306),
                                  "MT@ORA+T3": (0.211, 0.322)}}
    pub_dif = {("std", "LOC+T", "both", "MT@S2+T3", "MT@S0+T3"): "+0.070 [+0.010, +0.120] DE 0.0282",
               ("std", "LOC+T", "both", "MT@S2+T3", "MT@S0"): "+0.094 [+0.035, +0.151] DE 0.0298",
               ("std", "LOC+T", "both", "MT@ORA+T3", "MT@S0+T3"): "+0.091 [+0.031, +0.141]",
               ("std", "LOC+T", "both", "BOMP@S2+T3", "BOMP@S0+T3"): "+0.057 [-0.003, +0.115] DE 0.0304",
               ("std", "LOC+T", "thin", "MT@S2+T3", "MT@S0+T3"): "+0.073 [-0.021, +0.161]",
               ("std", "LOC+T", "strong", "MT@S2+T3", "MT@S0+T3"): "+0.075 [-0.022, +0.160]",
               ("std", "LOC", "both", "MT@S2+T3", "MT@S0+T3"): "+0.017 [-0.028, +0.068] DE 0.0244",
               ("dense", "LOC+T", "both", "MT@S2+T3", "MT@S0+T3"): "+0.065 [+0.034, +0.094] DE 0.0151",
               ("dense", "LOC+T", "thin", "MT@S2+T3", "MT@S0+T3"): "+0.063 [+0.019, +0.104]",
               ("dense", "LOC+T", "strong", "MT@S2+T3", "MT@S0+T3"): "+0.061 [+0.020, +0.110]",
               ("dense", "LOC", "both", "MT@S2+T3", "MT@S0+T3"): "+0.034 [+0.001, +0.070] DE 0.0179"}
    lines = []
    for kind, design in (("std", "STD"), ("dense", "DENSE")):
        Dk, nnet = pilot_as_data(rows, design)
        CC = A.counts(999, kind, nnet)
        for rule, tol in A.RULES:
            recs = A.build(Dk, tol)
            ms = list(pub_rec.get((kind, rule), {})) or ["MT@S0+T3", "MT@S2+T3"]
            ms = sorted(set(ms) | {x for k_, v in pub_dif.items() if k_[:2] == (kind, rule) for x in k_[3:]})
            for grp in ("both",) + REGIMES:
                if not any(k_[:3] == (kind, rule, grp) for k_ in pub_dif) and grp != "both":
                    continue
                e = A.evaluate(Dk, recs, CC, A.keys_of(grp), ms, nnet)
                if grp == "both":
                    ks = np.concatenate([Dk[k_]["k"] for k_ in A.ALL])
                    for m in ms:
                        rr = [[tuple(x) for x in r] for k_ in A.ALL for r in recs[k_][m]]
                        rx = xc.recall_caps(rr, int(ks.sum()))[0]
                        eq = abs(e[0][m][0.1] - rx[0]) < 1e-12 and abs(e[0][m][0.25] - rx[1]) < 1e-12
                        txt = f"  {kind:5s} {rule:5s} {m:12s} Pool {e[0][m][0.1]:.3f} / {e[0][m][0.25]:.3f}; xc.recall_caps " \
                              f"{rx[0]:.3f} / {rx[1]:.3f} {'igual' if eq else 'DISTINTO'}"
                        allok &= eq
                        if m in pub_rec.get((kind, rule), {}):
                            pb = pub_rec[(kind, rule)][m]
                            eqp = abs(round(rx[0], 3) - pb[0]) < 5e-4 and abs(round(rx[1], 3) - pb[1]) < 5e-4
                            allok &= eqp
                            txt += f"; publicado {pb[0]:.3f} / {pb[1]:.3f} {'igual' if eqp else 'DISTINTO'}"
                        P(txt)
                for k_, pubtxt in pub_dif.items():
                    if k_[:3] != (kind, rule, grp):
                        continue
                    d = A.diff(e, k_[3], k_[4], 0.25)
                    eqp = abs(round(d["point"], 3) - float(pubtxt.split()[0])) < 1.5e-3
                    allok &= eqp
                    lines.append(f"  {kind:5s} {rule:5s} {grp:6s} {k_[3]} - {k_[4]}: {A.fci(d)} DE {d['sd']:.4f}  |  "
                                 f"publicado {pubtxt}  (punto {'igual' if eqp else 'DISTINTO'})")
                if kind == "std" and rule == "LOC+T" and grp == "both":
                    f1 = A.ratio(e, "MT@S2+T3", "MT@S0+T3", "MT@ORA+T3", 0.25)
                    lines.append(f"  std   LOC+T both   fraccion de la ganancia del oraculo (forma T3, MT): {f1[0]:.2f} "
                                 f"[{f1[1]:.2f}, {f1[2]:.2f}]  |  publicado 0.77 [0.18, 1.45]")
    P("  diferencias a tope 0.25 (el punto puede diferir en +-0.001 por el redondeo de los recalls publicados):")
    out.extend(lines)

    # ------------------------------------------------------------------ (5)
    P("")
    P("(5) Tiempo (partes de control, 2 procesos en paralelo) y regla de presupuesto")
    parts = [part(reg, cell, j) for reg in REGIMES for cell in CELLS for j in R.CONTROL_NETS]

    def tsum(e, kind, ms):
        return sum(w["sec"][m] for w in e["windows"][kind] for m in ms)
    over = float(np.mean([e["t_total"] - sum(sum(w["sec"].values()) for k in e["windows"] for w in e["windows"][k])
                          for e in parts]))
    t_std_mt = float(np.mean([tsum(e, "std", R.MT3) for e in parts]))
    t_std_bo = float(np.mean([tsum(e, "std", R.BOMP3) for e in parts]))
    t_std_or = float(np.mean([tsum(e, "std", (R.ORA_MT, R.ORA_BOMP)) for e in parts]))
    t_den_mt = float(np.mean([tsum(e, "dense", R.MT3) for e in parts]))
    t_den_or = float(np.mean([tsum(e, "dense", (R.ORA_MT,)) for e in parts]))
    P(f"  por parte (5 ventanas por diseno; media de 18 partes): preparacion + generacion {over:.1f} s; estandar MT x3 "
      f"{t_std_mt:.1f} s, BOMP x3 {t_std_bo:.1f} s, oraculo (MT y BOMP) {t_std_or:.1f} s; densas MT x3 {t_den_mt:.1f} s, "
      f"oraculo MT {t_den_or:.1f} s")
    for name in ("main", "budget"):
        sp = R.ARMS[name]
        npart = max(sp["std"], sp["dense"]) * 6 * 2
        cpu = 6 * 2 * (max(sp["std"], sp["dense"]) * over + sp["std"] * t_std_mt + sp["bomp"] * t_std_bo +
                       sp["ora"] * (t_std_or + t_den_or) + sp["dense"] * t_den_mt)
        P(f"  diseno {name} (estandar {sp['std']}, densas {sp['dense']}, BOMP {sp['bomp']}, oraculo {sp['ora']} redes): "
          f"{npart} partes en las dos semillas -> {cpu / 3600:.2f} h de CPU, {cpu / 7200:.2f} h en 2 nucleos "
          "(cota superior: la preparacion de las partes sin ventanas densas es menor)")
        if name == "main":
            est = cpu / 7200
    P(f"  Regla de presupuesto (10 h en 2 nucleos): estimacion {est:.2f} h -> "
      f"{'SE ACTIVA: estandar 120 redes y densas 80' if est > 10 else 'NO se activa: estandar 160 redes y densas 100'}")
    kern = sorted({json.dumps(e["kernel"]) for e in parts})
    P("")
    P(f"Nucleos OpenBLAS registrados en las partes de control: {kern}")
    P(f"RESULTADO GLOBAL DE LOS CONTROLES: {'OK' if allok else 'FALLA'}")
    txt = "\n".join(out) + "\n"
    with open(a.out, "w") as fh:
        fh.write(txt)
    print(txt)


if __name__ == "__main__":
    main()
