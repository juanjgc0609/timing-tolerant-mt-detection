"""analyze_2c.py - reglas 2C-a, 2C-b y 2C-c de la Tarea 2C (PREREGISTRO_1j_1k.txt), por regimen de clutter, y tablas.

Datos: MT-FW+SC y QM-3C+SC de las redes 0-4 de 2A (task2a_parts, 20 ventanas por red); BOMP+SC, MTSTACK+SC y GL-FW+SC
de las mismas ventanas (task2c_parts/p1_*, parte (i); se comprueba ventana a ventana que verdad, origenes y k
coinciden); MT-FW+SC, QM-3C+SC, BOMP+SC y MTSTACK+SC de las redes 5-29 (task2c_parts/p2_*, parte (ii), 12 ventanas por
red).
Endpoint como analyze_2a.py: recall de eventos a eventos falsos realizados por ventana <= 0.1 y <= 0.25 sobre la curva
FROC exacta (froc_exact.exact_HF), punto de operacion elegido sobre la propia evaluacion; una diferencia cuenta si es
>= 0.05 - 1e-9.
2C-a (parte (i), 100 ventanas por celda): MT-FW+SC - max(BOMP+SC, MTSTACK+SC) >= 0.05 a tope 0.25 en >= 2 de 3 celdas
  en ambas semillas (estimaciones puntuales). IC 95%: bootstrap pareado por ventanas (2000 replicas,
  default_rng([s, ci, ri]), las remuestras de analyze_2a.py) y por red (5 redes con reemplazo, default_rng([s, ci, 77])
  como verif_2a/cluster_boot_verif.py); el maximo y el punto de operacion se re-eligen en cada replica.
2C-b (parte (i), descriptiva): GL-FW+SC - MT-FW+SC a ambos topes con los mismos IC; iguala si |dif| < 0.05, supera si
  dif >= 0.05, queda por debajo si dif <= -0.05.
2C-c (30 redes por celda: las 5 de 2A con 20 ventanas y las 25 nuevas con 12): la cota inferior del IC 95 del bootstrap
  jerarquico (2000 replicas; en cada una 30 redes con reemplazo y, dentro de cada red elegida, sus ventanas con
  reemplazo; default_rng([s, ci, ri, 2603]); punto de operacion re-elegido en cada replica, con la tasa de falsos
  sobre las ventanas remuestreadas) de MT-FW+SC - QM-3C+SC a tope 0.25 es > 0 en >= 2 de 3 celdas en ambas semillas y
  la estimacion puntual (todas las ventanas de las 30 redes juntas) es >= 0.05. Lo mismo para - BOMP+SC y - MTSTACK+SC
  y a tope 0.1. Holm: p unilateral por bootstrap p = (1 + #{d* <= 0}) / (B + 1), Holm al nivel 0.025 (el de la cota
  inferior de un IC bilateral 95) sobre la familia de las 36 comparaciones de cada regimen (3 comparadores x 2 topes
  x 3 celdas x 2 semillas); la regla 2C-c se lee sin correccion (texto literal) y con Holm.
Exploratorio (declarado): GL-FW+SC con desempate de puntajes iguales por la norma de grupo al entrar; 2C-c solo en
las 25 redes nuevas; fraccion de replicas en que se cumplirian las reglas.

    python analyze_2c.py [--out analyze_2c_out.txt]
"""
import argparse
import json
import os

import numpy as np

from froc_exact import exact_HF

HERE = os.path.dirname(os.path.abspath(__file__))
B = 2000
CAPS = ("0.1", "0.25")
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
LAB = {"ladder_4": "4 st.", "ladder_6": "6 st.", "unilateral_c40": "one-sided 8"}
REGIMES = ("thin", "strong")
SEEDS = (2024, 7331)
EPS = 1e-9
THR = 0.05
ALPHA = 0.025
NETS_A = range(5)
NETS_II = range(5, 30)
M_I = ("MT-FW+SC", "BOMP+SC", "MTSTACK+SC", "GL-FW+SC", "QM-3C+SC")
M_II = ("MT-FW+SC", "QM-3C+SC", "BOMP+SC", "MTSTACK+SC")
COMP_C = ("QM-3C+SC", "BOMP+SC", "MTSTACK+SC")


def ge(d):
    return d >= THR - EPS


def fmt_ci(p, lo, hi):
    return f"{p:+.3f} [{lo:+.3f}, {hi:+.3f}]"


def HF(recs):
    _t, H, F = exact_HF([[tuple(x) for x in r] for r in recs])
    return np.asarray(H, float), np.asarray(F, float)


def val(HFm, C, k):
    """recall a cada tope para cada fila de C (conteos de remuestreo): max sobre umbrales con falsos por ventana
    remuestreada <= tope."""
    H, F = HFm
    Hs, Fs = C @ H, C @ F
    nn = C.sum(1)[:, None]
    kk = np.maximum(1.0, C @ k)
    out = {}
    for cap in CAPS:
        ok = Fs / nn <= float(cap) + 1e-12
        r = np.where(ok, Hs, -1.0).max(axis=1)
        out[cap] = np.where(r < 0, 0.0, r) / kk
    return out


def point(HFm, k):
    v = val(HFm, np.ones((1, len(k))), k)
    return {cap: float(v[cap][0]) for cap in CAPS}


def counts_windows(s, ci, ri, n):
    rng = np.random.default_rng([s, ci, ri])
    idx = rng.integers(0, n, size=(B, n))
    C = np.zeros((B, n))
    np.add.at(C, (np.repeat(np.arange(B), n), idx.ravel()), 1.0)
    return C


def counts_networks(s, ci, net):
    rng = np.random.default_rng([s, ci, 77])
    ids = sorted(set(net))
    C = np.zeros((B, len(net)))
    for b in range(B):
        for j in rng.integers(0, len(ids), len(ids)):
            C[b, net == ids[j]] += 1
    return C


def counts_hier(s, ci, ri, net):
    rng = np.random.default_rng([s, ci, ri, 2603])
    ids = sorted(set(net))
    pos = {j: np.where(net == j)[0] for j in ids}
    C = np.zeros((B, len(net)))
    for b in range(B):
        for jj in rng.integers(0, len(ids), len(ids)):
            p = pos[ids[jj]]
            np.add.at(C[b], p[rng.integers(0, len(p), len(p))], 1.0)
    return C


def holm(pvals, alpha=ALPHA):
    """rechazos de Holm y p ajustados (paso descendente) para un diccionario {clave: p}."""
    keys = sorted(pvals, key=lambda k_: pvals[k_])
    m = len(keys)
    adj, run, rej, stop = {}, 0.0, {}, False
    for i, k_ in enumerate(keys):
        run = max(run, min(1.0, (m - i) * pvals[k_]))
        adj[k_] = run
        if not stop and pvals[k_] <= alpha / (m - i):
            rej[k_] = True
        else:
            stop = True
            rej[k_] = False
    return rej, adj


# ------------------------------------------------------------------------------------------ carga
def load(path):
    with open(path) as fh:
        return json.load(fh)


def tie_break(rec_nv):
    """GL-FW+SC con desempate (exploratorio): dentro de cada grupo de puntajes iguales (mismo lambda de entrada), el
    orden del registro (norma de grupo decreciente) se hace estricto restando 1e-9 lambda por posicion."""
    out, prev, pos = [], None, 0
    for lam, hit, _nv, _n, _t in rec_nv:
        pos = pos + 1 if lam == prev else 0
        prev = lam
        out.append([lam * (1.0 - 1e-9 * pos), hit])
    return out


def cell_data(s, reg, cell):
    """ventanas de la celda: part (i) (redes 0-4) y parte (ii) (5-29). Devuelve dicts de listas."""
    W = {"k": [], "net": [], "rec": {m: [] for m in set(M_I) | set(M_II) | {"GL-FW+SC-TB"}}, "gl": [], "part": []}
    ctl = {"windows_2a": 0, "mismatch": 0, "p1_controls": [], "kernels": set(), "p2_controls": []}
    for j in NETS_A:
        a = load(os.path.join(HERE, "task2a_parts", f"seed{s}_{reg}_{cell}_net{j}.json"))
        c = load(os.path.join(HERE, "task2c_parts", f"p1_seed{s}_{reg}_{cell}_net{j}.json"))
        ctl["p1_controls"].append(c["control"])
        ctl["kernels"].add(tuple(tuple(x) for x in c["kernel"]["blas"]))
        assert c["n_win"] == a["n_win"] == len(c["windows"])
        for wa, wc in zip(a["windows"], c["windows"]):
            same = wa["k"] == wc["k"] and wa["truth"] == wc["truth"] and wa["taus"] == wc["taus"]
            ctl["windows_2a"] += 1
            ctl["mismatch"] += not same
            W["k"].append(wa["k"])
            W["net"].append(j)
            W["part"].append(1)
            for m in ("MT-FW+SC", "QM-3C+SC"):
                W["rec"][m].append(wa["records"][m])
            for m in ("BOMP+SC", "MTSTACK+SC", "GL-FW+SC"):
                W["rec"][m].append(wc["records"][m])
            W["rec"]["GL-FW+SC-TB"].append(tie_break(wc["gl_rec_nv"]))
            W["gl"].append(wc["gl"])
    for j in NETS_II:
        path = os.path.join(HERE, "task2c_parts", f"p2_seed{s}_{reg}_{cell}_net{j}.json")
        if not os.path.exists(path):
            continue
        c = load(path)
        ctl["p2_controls"].append(c["control"])
        ctl["kernels"].add(tuple(tuple(x) for x in c["kernel"]["blas"]))
        for wc in c["windows"]:
            W["k"].append(wc["k"])
            W["net"].append(j)
            W["part"].append(2)
            for m in M_II:
                W["rec"][m].append(wc["records"][m])
    W["k"] = np.asarray(W["k"], float)
    W["net"] = np.asarray(W["net"])
    W["part"] = np.asarray(W["part"])
    return W, ctl


def sub(W, mask, methods):
    idx = np.where(mask)[0]
    return W["k"][idx], W["net"][idx], {m: HF([W["rec"][m][i] for i in idx]) for m in methods}


# ------------------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "analyze_2c_out.txt"))
    a = ap.parse_args()
    L = []
    P = L.append
    verdicts = {}
    near = []
    P("# Tarea 2C - comparadores cercanos, metodo convexo y replicacion por red (analyze_2c.py)")
    P("")
    P("Recall de eventos a eventos falsos realizados por ventana <= 0.1 / <= 0.25 (curva FROC exacta, punto de operacion "
      "elegido sobre la propia evaluacion). Diferencias: punto [IC 95%]. W = bootstrap pareado por ventanas; N = bootstrap "
      "por red; H = bootstrap jerarquico (redes, luego ventanas). 2000 replicas; punto de operacion y maximo re-elegidos "
      "en cada replica.")
    all_data = {}
    for ri, reg in enumerate(REGIMES):
        for s in SEEDS:
            for ci, cell in enumerate(CELLS):
                all_data[(reg, s, cell)] = cell_data(s, reg, cell)
    # ---------------------------------------------------------------------------------------- controles
    P("")
    P("## Controles")
    P("")
    for reg in REGIMES:
        for s in SEEDS:
            for cell in CELLS:
                W, ctl = all_data[(reg, s, cell)]
                p1ok = all(all(v for k_, v in c.items() if k_ not in ("windows_truth_taus", "mtfw_sc_exact"))
                           and c["mtfw_sc_exact"] == 2 for c in ctl["p1_controls"])
                p2 = ctl["p2_controls"]
                P(f"- {reg} {s} {LAB[cell]}: parte (i) {len(ctl['p1_controls'])} redes, control de red y MT-FW+SC exacto "
                  f"(ventanas 0-1) {'correcto' if p1ok else 'FALLA'}; ventanas 2A emparejadas {ctl['windows_2a']} "
                  f"(distintas {ctl['mismatch']}); parte (ii) {len(p2)} redes, MT-FW+SC certificado = denso en "
                  f"{sum(c['mtfw_cert_eq_dense'] for c in p2)} ventanas; nucleos {sorted(ctl['kernels'])}")
                assert ctl["mismatch"] == 0 and p1ok
    # ---------------------------------------------------------------------------------------- parte (i)
    boots = {}
    for ri, reg in enumerate(REGIMES):
        P("")
        P(f"## Regimen {reg}: parte (i), redes 0-4 de 2A (100 ventanas por celda)")
        P("")
        P(f"### Tabla 1 ({reg}). Recall a <= 0.1 / <= 0.25")
        P("")
        P("| Celda (semilla) | eventos | " + " | ".join(M_I) + " | GL-FW+SC-TB (expl.) |")
        P("|:--|--:|" + ":-:|" * (len(M_I) + 1))
        rowsA, rowsB, cntA, detA, clsB = [], [], {s: 0 for s in SEEDS}, {s: [] for s in SEEDS}, []
        repA = {s: [] for s in SEEDS}
        for s in SEEDS:
            for ci, cell in enumerate(CELLS):
                W, _ctl = all_data[(reg, s, cell)]
                k, net, HFs = sub(W, W["part"] == 1, M_I + ("GL-FW+SC-TB",))
                pt = {m: point(HFs[m], k) for m in HFs}
                P(f"| {LAB[cell]} ({s}) | {int(k.sum())} | " + " | ".join(
                    f"{pt[m]['0.1']:.3f}/{pt[m]['0.25']:.3f}" for m in M_I + ("GL-FW+SC-TB",)) + " |")
                Cw = counts_windows(s, ci, ri, len(k))
                Cn = counts_networks(s, ci, net)
                bw = {m: val(HFs[m], Cw, k) for m in HFs}
                bn = {m: val(HFs[m], Cn, k) for m in HFs}
                cellsA = []
                okA = None
                for cap in CAPS:
                    p0 = pt["MT-FW+SC"][cap] - max(pt["BOMP+SC"][cap], pt["MTSTACK+SC"][cap])
                    dw = bw["MT-FW+SC"][cap] - np.maximum(bw["BOMP+SC"][cap], bw["MTSTACK+SC"][cap])
                    dn = bn["MT-FW+SC"][cap] - np.maximum(bn["BOMP+SC"][cap], bn["MTSTACK+SC"][cap])
                    best = "BOMP+SC" if pt["BOMP+SC"][cap] >= pt["MTSTACK+SC"][cap] else "MTSTACK+SC"
                    cellsA.append(f"{best}: {p0:+.3f} W[{np.percentile(dw, 2.5):+.3f}, {np.percentile(dw, 97.5):+.3f}] "
                                  f"N[{np.percentile(dn, 2.5):+.3f}, {np.percentile(dn, 97.5):+.3f}]")
                    if cap == "0.25":
                        okA = ge(p0)
                        repA[s].append(dw >= THR - EPS)
                        if abs(p0 - THR) < 0.005:
                            near.append(f"2C-a {reg} {s} {cell}: {p0:+.4f}")
                        detA[s].append((cell, okA, p0))
                cntA[s] += okA
                sep = {}
                for m in ("BOMP+SC", "MTSTACK+SC"):
                    sep[m] = " / ".join(fmt_ci(pt["MT-FW+SC"][cap] - pt[m][cap],
                                               *np.percentile(bw["MT-FW+SC"][cap] - bw[m][cap], [2.5, 97.5]))
                                        for cap in CAPS)
                rowsA.append(f"| {LAB[cell]} ({s}) | {cellsA[0]} | {cellsA[1]} | {'si' if okA else 'no'} | "
                             f"{sep['BOMP+SC']} | {sep['MTSTACK+SC']} |")
                cb = []
                for cap in CAPS:
                    p0 = pt["GL-FW+SC"][cap] - pt["MT-FW+SC"][cap]
                    dw = bw["GL-FW+SC"][cap] - bw["MT-FW+SC"][cap]
                    dn = bn["GL-FW+SC"][cap] - bn["MT-FW+SC"][cap]
                    cls = "iguala" if abs(p0) < THR - EPS else ("supera" if p0 > 0 else "por debajo")
                    clsB.append((s, cell, cap, cls, p0))
                    p1 = pt["GL-FW+SC-TB"][cap] - pt["MT-FW+SC"][cap]
                    p2 = pt["GL-FW+SC"][cap] - pt["BOMP+SC"][cap]
                    cb.append(f"{p0:+.3f} W[{np.percentile(dw, 2.5):+.3f}, {np.percentile(dw, 97.5):+.3f}] "
                              f"N[{np.percentile(dn, 2.5):+.3f}, {np.percentile(dn, 97.5):+.3f}] {cls} "
                              f"(TB {p1:+.3f}; GL - BOMP {p2:+.3f})")
                rowsB.append(f"| {LAB[cell]} ({s}) | {cb[0]} | {cb[1]} |")
        P("")
        P(f"### Tabla 2 ({reg}). 2C-a: MT-FW+SC - max(BOMP+SC, MTSTACK+SC) (mejor comparador por celda y tope), y contra cada uno (W)")
        P("")
        P("| Celda (semilla) | a 0.1 | a 0.25 | 2C-a (>= 0.05 a 0.25) | - BOMP+SC (0.1 / 0.25) | - MTSTACK+SC (0.1 / 0.25) |")
        P("|:--|:-:|:-:|:-:|:-:|:-:|")
        L.extend(rowsA)
        okA_reg = all(cntA[s] >= 2 for s in SEEDS)
        joint = np.all([np.sum(repA[s], axis=0) >= 2 for s in SEEDS], axis=0)
        P("")
        P(f"- Regla 2C-a ({reg}): " + "; ".join(
            f"semilla {s}: {cntA[s]}/3 celdas (" + ", ".join(f"{LAB[c_]} {'si' if o else 'no'} {x:+.3f}" for c_, o, x in detA[s])
            + ")" for s in SEEDS) + f" -> {'SE CUMPLE' if okA_reg else 'NO SE CUMPLE'} (exploratorio: se cumple en "
            f"{100 * joint.mean():.0f}% de las replicas W)")
        P("")
        P(f"### Tabla 3 ({reg}). 2C-b: GL-FW+SC - MT-FW+SC (descriptiva; iguala si |dif| < 0.05)")
        P("")
        P("| Celda (semilla) | a 0.1 | a 0.25 |")
        P("|:--|:-:|:-:|")
        L.extend(rowsB)
        verdicts[reg] = {"2C-a": okA_reg, "2C-a_counts": cntA, "2C-a_replicate_fraction": float(joint.mean()),
                         "2C-b": [list(x) for x in clsB]}
        # diagnosticos de GL-FW
        gl = [g for s in SEEDS for cell in CELLS for g in all_data[(reg, s, cell)][0]["gl"]]
        pl = [p for g in gl for p in g["per_lambda"]]
        fin = [g["per_lambda"][-1] for g in gl]
        nrec = [len(r) for s in SEEDS for cell in CELLS for r in all_data[(reg, s, cell)][0]["rec"]["GL-FW+SC"][:100]]
        ties = []
        for s in SEEDS:
            for cell in CELLS:
                for r in all_data[(reg, s, cell)][0]["rec"]["GL-FW+SC"][:100]:
                    sc = [x[0] for x in r]
                    ties.append(len(sc) - len(set(sc)))
        P("")
        P(f"Diagnosticos de GL-FW+SC ({reg}, {len(gl)} ventanas, {len(pl)} lambdas resueltos): lambdas resueltos por "
          f"ventana {min(len(g['per_lambda']) for g in gl)}-{max(len(g['per_lambda']) for g in gl)} (mediana "
          f"{np.median([len(g['per_lambda']) for g in gl]):.0f}); parada: registro completo en "
          f"{sum(1 for g in gl if g['stopped'] and g['stopped'][0] == 'record_complete')}, trayectoria agotada en "
          f"{sum(1 for g in gl if not g['stopped'])}; brecha relativa max {max(p[2] for p in pl if np.isfinite(p[2])):.1e} "
          f"(en el ultimo lambda: mediana {np.median([p[2] for p in fin]):.1e}, max {max(p[2] for p in fin):.1e}); "
          f"brecha / (lambda Omega) max {max(p[3] for p in pl):.1e}; KKT fuera del soporte max {max(p[6] for p in pl):.9f}; "
          f"max |c/lambda - 1| en el soporte {max(p[7] for p in pl):.1e}; re-solves {sum(p[8] for p in pl)}; soporte final "
          f"mediana {np.median([p[4] for p in fin]):.0f}, max {max(p[4] for p in fin)}; ventanas con < 5 eventos en el "
          f"registro {sum(n < 5 for n in nrec)}; eventos con puntaje empatado dentro de la ventana: {sum(ties)} de "
          f"{sum(nrec)} (ventanas con algun empate {sum(t > 0 for t in ties)}).")
    # ---------------------------------------------------------------------------------------- parte (ii)
    for ri, reg in enumerate(REGIMES):
        P("")
        P(f"## Regimen {reg}: 2C-c, 30 redes por celda (redes 0-4 de 2A con 20 ventanas + 25 nuevas con 12)")
        P("")
        P(f"### Tabla 4 ({reg}). Recall a <= 0.1 / <= 0.25 (todas las ventanas de las 30 redes)")
        P("")
        P("| Celda (semilla) | redes | ventanas | eventos | " + " | ".join(M_II) + " |")
        P("|:--|--:|--:|--:|" + ":-:|" * len(M_II))
        tests, pts, rows, repC = {}, {}, {}, {s: [] for s in SEEDS}
        rows_new = []
        for s in SEEDS:
            for ci, cell in enumerate(CELLS):
                W, _ctl = all_data[(reg, s, cell)]
                k, net, HFs = sub(W, np.ones(len(W["k"]), bool), M_II)
                nnet = len(set(net))
                pt = {m: point(HFs[m], k) for m in M_II}
                P(f"| {LAB[cell]} ({s}) | {nnet} | {len(k)} | {int(k.sum())} | " + " | ".join(
                    f"{pt[m]['0.1']:.3f}/{pt[m]['0.25']:.3f}" for m in M_II) + " |")
                C = counts_hier(s, ci, ri, net)
                bv = {m: val(HFs[m], C, k) for m in M_II}
                for m in COMP_C:
                    for cap in CAPS:
                        d = bv["MT-FW+SC"][cap] - bv[m][cap]
                        p0 = pt["MT-FW+SC"][cap] - pt[m][cap]
                        lo, hi = np.percentile(d, [2.5, 97.5])
                        pv = (1.0 + np.sum(d <= 0)) / (B + 1.0)
                        tests[(s, cell, m, cap)] = pv
                        pts[(s, cell, m, cap)] = (p0, lo, hi)
                        if m == "QM-3C+SC" and cap == "0.25":
                            repC[s].append((d > 0) & (p0 >= THR - EPS))
                # exploratorio: solo las 25 redes nuevas
                k2, net2, HF2 = sub(W, W["part"] == 2, M_II)
                if len(k2):
                    pt2 = {m: point(HF2[m], k2) for m in M_II}
                    rows_new.append(f"| {LAB[cell]} ({s}) | {len(set(net2))} | {len(k2)} | " + " | ".join(
                        f"{pt2[m]['0.1']:.3f}/{pt2[m]['0.25']:.3f}" for m in M_II) + " | " + " / ".join(
                        f"{pt2['MT-FW+SC'][c_] - pt2['QM-3C+SC'][c_]:+.3f}" for c_ in CAPS) + " |")
        rej, adj = holm(tests)
        P("")
        P(f"### Tabla 5 ({reg}). MT-FW+SC menos cada comparador: punto [IC 95% H]; p unilateral bootstrap; p de Holm "
          f"(familia de {len(tests)} comparaciones del regimen, nivel {ALPHA})")
        P("")
        P("| Celda (semilla) | comparador | a 0.1 | a 0.25 |")
        P("|:--|:--|:-:|:-:|")
        cnt_lit = {s: 0 for s in SEEDS}
        cnt_holm = {s: 0 for s in SEEDS}
        det = {s: [] for s in SEEDS}
        for s in SEEDS:
            for cell in CELLS:
                for m in COMP_C:
                    cc = []
                    for cap in CAPS:
                        p0, lo, hi = pts[(s, cell, m, cap)]
                        key = (s, cell, m, cap)
                        cc.append(f"{fmt_ci(p0, lo, hi)} p={tests[key]:.4f} pHolm={adj[key]:.4f}"
                                  f"{' *' if rej[key] else ''}")
                    P(f"| {LAB[cell]} ({s}) | {m} | {cc[0]} | {cc[1]} |")
                p0, lo, hi = pts[(s, cell, "QM-3C+SC", "0.25")]
                ok_lit = lo > 0 and ge(p0)
                ok_h = rej[(s, cell, "QM-3C+SC", "0.25")] and ge(p0)
                cnt_lit[s] += ok_lit
                cnt_holm[s] += ok_h
                det[s].append((cell, ok_lit, ok_h, p0, lo))
                if abs(p0 - THR) < 0.005:
                    near.append(f"2C-c {reg} {s} {cell}: {p0:+.4f}")
        okC = all(cnt_lit[s] >= 2 for s in SEEDS)
        okCh = all(cnt_holm[s] >= 2 for s in SEEDS)
        jointC = np.all([np.sum(repC[s], axis=0) >= 2 for s in SEEDS], axis=0)
        P("")
        P("(* = rechazo de Holm, MT-FW+SC mayor que el comparador)")
        P("")
        P(f"- Regla 2C-c ({reg}), lectura literal (cota inferior del IC 95 H > 0 y punto >= 0.05, MT-FW+SC - QM-3C+SC a "
          "0.25): " + "; ".join(f"semilla {s}: {cnt_lit[s]}/3 (" + ", ".join(
              f"{LAB[c_]} {'si' if o else 'no'} {p0:+.3f}, cota {lo:+.3f}" for c_, o, _oh, p0, lo in det[s]) + ")"
              for s in SEEDS) + f" -> {'SE CUMPLE' if okC else 'NO SE CUMPLE'}")
        P(f"- Regla 2C-c ({reg}) con Holm (rechazo de Holm y punto >= 0.05): " + "; ".join(
            f"semilla {s}: {cnt_holm[s]}/3" for s in SEEDS) + f" -> {'SE CUMPLE' if okCh else 'NO SE CUMPLE'}")
        P(f"- Exploratorio: la lectura literal se cumple en {100 * np.mean(jointC):.0f}% de las replicas H (d* > 0 y "
          "punto >= 0.05, >= 2 de 3 celdas en ambas semillas).")
        verdicts[reg]["2C-c"] = okC
        verdicts[reg]["2C-c_holm"] = okCh
        verdicts[reg]["2C-c_counts"] = cnt_lit
        verdicts[reg]["2C-c_holm_counts"] = cnt_holm
        verdicts[reg]["holm_rejections"] = {"|".join(map(str, k_)): bool(v) for k_, v in rej.items()}
        P("")
        P(f"### Tabla 6 ({reg}, exploratoria). Solo las 25 redes nuevas (12 ventanas cada una): recall a 0.1 / 0.25 y "
          "MT-FW+SC - QM-3C+SC")
        P("")
        P("| Celda (semilla) | redes | ventanas | " + " | ".join(M_II) + " | MT-FW+SC - QM-3C+SC |")
        P("|:--|--:|--:|" + ":-:|" * (len(M_II) + 1))
        L.extend(rows_new)
    P("")
    P("## Veredictos")
    P("")
    P("| Regla | " + " | ".join(REGIMES) + " |")
    P("|:--|" + ":-:|" * len(REGIMES))
    P("| 2C-a | " + " | ".join("se cumple" if verdicts[g]["2C-a"] else "no se cumple" for g in REGIMES) + " |")
    P("| 2C-c (literal) | " + " | ".join("se cumple" if verdicts[g]["2C-c"] else "no se cumple" for g in REGIMES) + " |")
    P("| 2C-c (Holm) | " + " | ".join("se cumple" if verdicts[g]["2C-c_holm"] else "no se cumple" for g in REGIMES) + " |")
    for g in REGIMES:
        cl = verdicts[g]["2C-b"]
        P(f"- 2C-b {g}: " + "; ".join(f"{s} {LAB[c_]} a {cap}: {cls} ({p0:+.3f})" for s, c_, cap, cls, p0 in cl))
    P("")
    P("Casos a menos de 0.005 del umbral 0.05: " + ("; ".join(near) if near else "ninguno"))
    txt = "\n".join(L)
    print(txt)
    with open(a.out, "w") as fh:
        fh.write(txt + "\n")
    with open(a.out.replace(".txt", ".json"), "w") as fh:
        json.dump(verdicts, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
