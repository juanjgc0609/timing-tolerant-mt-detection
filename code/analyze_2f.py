"""analyze_2f.py - reglas 2F-a (con la lectura estricta de Holm), 2F-b y 2F-c de la Tarea 2F (PREREGISTRO_1j_1k.txt),
diferencias descriptivas por celda, tope 0.1 y tiempos.  Escrito y congelado (task2f_code_hashes.txt) antes de correr
codigo alguno sobre las semillas 4242 y 8888.

Datos: task2f_parts/ (task2f_run.py).  Brazo denso: 96 redes (j = 0..95) x 10 ventanas por celda, regimen y semilla,
kmax 7.  Brazo estandar: 16 redes (j = 96..111) x 10 ventanas, kmax 5.
Endpoint (como 2A-2E): recall de eventos a eventos falsos realizados por ventana <= 0.25 (primario) y <= 0.1
(secundario), curva FROC exacta (todos los puntajes del registro como umbrales; = froc_exact.exact_HF +
analyze_2c.val, comprobado en task2f_control.py), punto de operacion elegido por metodo sobre la propia evaluacion.
Al agrupar celdas (y regimenes) la curva FROC es la del conjunto agrupado: un unico umbral por metodo para todas las
ventanas agrupadas (como el piloto, explore_contrib/power.py y xc.hier_boot).
Bootstrap jerarquico: 2000 replicas; celdas como estratos; en cada celda las redes con reemplazo y, dentro de cada red
elegida, sus ventanas con reemplazo; las mismas remuestras para todos los metodos (diferencias pareadas); punto de
operacion re-elegido en cada replica.  Las redes j se comparten entre regimenes (misma geometria y retardos), por eso la
red elegida se usa en ambos regimenes (las ventanas se remuestrean por separado en cada regimen): en el analisis de un
solo regimen esto es exactamente "celdas como estratos, redes y ventanas", y al agrupar regimenes (2F-b, 2F-c) la red es
la unidad de replicacion.  Generador: default_rng([semilla, indice de celda, 4201 denso | 4202 estandar]).
IC 95: percentiles 2.5 y 97.5.  p unilateral de H0: d <= 0: (1 + #{d* <= 0}) / (B + 1).
Regla 2F-a (primaria), por regimen: MT-FW+SC - BOMP-PK+SC a tope 0.25, agrupando las 3 celdas: cota inferior > 0 y
  estimacion puntual >= 0.03 en ambas semillas.  Lectura estricta: Holm sobre las 4 pruebas (2 regimenes x 2 semillas),
  unilateral 0.025; el regimen cumple si Holm rechaza en ambas semillas y la estimacion puntual es >= 0.03 en ambas.
Regla 2F-b (mecanismo), a tope 0.25 (el endpoint primario; el registro no nombra el tope), agrupando celdas y
  regimenes: MT-FW+SC - MT-FW-NOCAP+SC con cota inferior > 0 en ambas semillas, y |MT-FW-NOCAP+SC - BOMP-PK+SC| <= 0.02
  (estimacion puntual) en ambas semillas.
Regla 2F-c (no dano, brazo estandar): TOST con margen 0.05 de MT-FW+SC - BOMP-PK+SC a tope 0.25, agrupando celdas y
  regimenes, en ambas semillas: p_inf = (1 + #{d* <= -0.05}) / (B + 1) < 0.05 y p_sup = (1 + #{d* >= +0.05}) / (B + 1)
  < 0.05 (equivale a que el IC 90 percentil quede dentro de (-0.05, +0.05)); se reporta tambien el IC 90.
Comparaciones con umbral: >= umbral - 1e-9, <= umbral + 1e-9.
Descriptivo: recall por metodo; diferencias por celda con IC a ambos topes; tope 0.1; BOMP-PK-CAP+SC - BOMP-PK+SC y
  MT-FW+SC - BOMP-PK-CAP+SC; tiempos por ventana; DE bootstrap de la diferencia primaria frente a la DE piloto.

    python analyze_2f.py [--out analyze_2f_out.txt]
"""
import argparse
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PARTS = os.path.join(HERE, "task2f_parts")
B = 2000
CAPS = (0.1, 0.25)
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
LAB = {"ladder_4": "4 st.", "ladder_6": "6 st.", "unilateral_c40": "one-sided 8"}
REGIMES = ("thin", "strong")
SEEDS = (4242, 8888)
METHODS = ("MT-FW+SC", "BOMP-PK+SC", "BOMP-PK-CAP+SC", "MT-FW-NOCAP+SC")
ARMS = {"dense": {"nets": tuple(range(0, 96)), "kmax": 7, "tag": 4201},
        "std": {"nets": tuple(range(96, 112)), "kmax": 5, "tag": 4202}}
N_WIN = 10
EPS = 1e-9
THR_A = 0.03
THR_B = 0.02
MARGIN_C = 0.05
ALPHA_HOLM = 0.025
ALPHA_TOST = 0.05
PILOT_SD = {"thin": 0.0383, "strong": 0.0353}
PAIRS = (("MT-FW+SC", "BOMP-PK+SC"), ("MT-FW+SC", "MT-FW-NOCAP+SC"), ("MT-FW-NOCAP+SC", "BOMP-PK+SC"),
         ("BOMP-PK-CAP+SC", "BOMP-PK+SC"), ("MT-FW+SC", "BOMP-PK-CAP+SC"))


# ============================================================================================ recall (FROC exacta)
class Pool:
    """registros de un conjunto de ventanas.  recall(C) da, para cada fila de C (pesos por ventana: conteos de
    remuestreo; unos = estimacion puntual), el recall a cada tope: max sobre umbrales (todos los puntajes) del numero
    ponderado de eventos verdaderos acertados con falsos ponderados por ventana <= tope; 0 si ninguno."""

    def __init__(self, recs, ks):
        self.ks = np.asarray(ks, float)
        sc, win, hit = [], [], []
        for w, r in enumerate(recs):
            for s_, h in r:
                sc.append(float(s_))
                win.append(w)
                hit.append(int(h))
        sc, win, hit = np.array(sc, float), np.array(win, int), np.array(hit, int)
        o = np.argsort(-sc, kind="stable")
        sc, win, hit = sc[o], win[o], hit[o]
        key = np.where(hit >= 0, win * 64 + hit, -1)
        _u, first = np.unique(key, return_index=True)
        newhit = np.zeros(len(sc), bool)
        newhit[first] = True
        newhit &= key >= 0
        self.win = win
        self.newhit = newhit.astype(float)
        self.false = (hit < 0).astype(float)
        self.ends = np.r_[np.where(np.diff(sc) != 0)[0], len(sc) - 1] if len(sc) else np.zeros(0, int)
        self.nw = len(recs)

    def recall(self, C, chunk=100):
        C = np.atleast_2d(np.asarray(C, float))
        out = {cap: np.zeros(len(C)) for cap in CAPS}
        for a in range(0, len(C), chunk):
            Cc = C[a:a + chunk]
            n = Cc.sum(1)
            kk = np.maximum(1.0, Cc @ self.ks)
            if len(self.win) == 0:
                continue
            Wt = Cc[:, self.win]
            H = np.cumsum(Wt * self.newhit, axis=1)[:, self.ends]
            F = np.cumsum(Wt * self.false, axis=1)[:, self.ends]
            for cap in CAPS:
                ok = F / n[:, None] <= cap + 1e-12
                out[cap][a:a + chunk] = np.where(ok, H, 0.0).max(axis=1) / kk
        return out


# ============================================================================================ datos
def load_arm(seed, arm, parts=PARTS, nets=None, prefix=None):
    """{(regime, cell): {"k", "net", "rec": {m: [...]}, "sec": {m: [...]}, "meta": [...]}}; ventanas ordenadas por red
    (en el orden de nets) y ventana."""
    nets = ARMS[arm]["nets"] if nets is None else nets
    prefix = arm if prefix is None else prefix
    D = {}
    for reg in REGIMES:
        for cell in CELLS:
            W = {"k": [], "net": [], "rec": {m: [] for m in METHODS}, "sec": {m: [] for m in METHODS}, "meta": []}
            for jj, j in enumerate(nets):
                path = os.path.join(parts, f"{prefix}_seed{seed}_{reg}_{cell}_net{j}.json")
                with open(path) as fh:
                    e = json.load(fh)
                assert (e["seed"], e["regime"], e["cell"], e["net"]) == (seed, reg, cell, j)
                assert e["n_win"] == len(e["windows"]) == N_WIN
                W["meta"].append({x: e[x] for x in ("tag", "kmax", "cat_drawn", "t_setup", "t_total", "kernel")})
                for w in e["windows"]:
                    W["k"].append(w["k"])
                    W["net"].append(jj)
                    for m in METHODS:
                        W["rec"][m].append(w["records"][m])
                        W["sec"][m].append(w["sec"][m])
            W["k"] = np.asarray(W["k"], float)
            W["net"] = np.asarray(W["net"], int)
            D[(reg, cell)] = W
    return D


def counts(seed, arm, nnet, nwin=N_WIN, B_=B):
    """{(regime, cell): C (B, nnet * nwin)}: remuestras jerarquicas (celda = estrato; la red elegida se usa en ambos
    regimenes; ventanas remuestreadas por separado en cada regimen)."""
    out = {}
    for ci, cell in enumerate(CELLS):
        rng = np.random.default_rng([seed, ci, ARMS[arm]["tag"]])
        C = {reg: np.zeros((B_, nnet * nwin)) for reg in REGIMES}
        for b in range(B_):
            pick = rng.integers(0, nnet, nnet)
            for reg in REGIMES:
                idx = pick[:, None] * nwin + rng.integers(0, nwin, (nnet, nwin))
                np.add.at(C[reg][b], idx.ravel(), 1.0)
        for reg in REGIMES:
            out[(reg, cell)] = C[reg]
    return out


def pooled(D, CC, keys):
    """Pool por metodo y matriz de conteos concatenada para las (regime, cell) de keys."""
    recs = {m: [r for k_ in keys for r in D[k_]["rec"][m]] for m in METHODS}
    ks = np.concatenate([D[k_]["k"] for k_ in keys])
    C = np.concatenate([CC[k_] for k_ in keys], axis=1)
    return {m: Pool(recs[m], ks) for m in METHODS}, C, int(ks.sum()), len(ks)


def evaluate(pools, C, n):
    """puntos y replicas por metodo."""
    one = np.ones((1, n))
    pt = {m: {cap: float(p.recall(one)[cap][0]) for cap in CAPS} for m, p in pools.items()}
    bs = {m: p.recall(C) for m, p in pools.items()}
    return pt, bs


def diff(pt, bs, a, b, cap):
    d = bs[a][cap] - bs[b][cap]
    p = pt[a][cap] - pt[b][cap]
    lo, hi = np.percentile(d, [2.5, 97.5])
    lo90, hi90 = np.percentile(d, [5.0, 95.0])
    return {"point": p, "lo": float(lo), "hi": float(hi), "lo90": float(lo90), "hi90": float(hi90),
            "p_le0": float((1 + np.sum(d <= 0)) / (len(d) + 1)), "sd": float(d.std()),
            "p_tost_lo": float((1 + np.sum(d <= -MARGIN_C)) / (len(d) + 1)),
            "p_tost_hi": float((1 + np.sum(d >= MARGIN_C)) / (len(d) + 1))}


def holm(pvals, alpha):
    """analyze_2c.holm (paso descendente)."""
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


def fci(d, w90=False):
    if w90:
        return f"{d['point']:+.3f} [{d['lo90']:+.3f}, {d['hi90']:+.3f}]"
    return f"{d['point']:+.3f} [{d['lo']:+.3f}, {d['hi']:+.3f}]"


def yn(x):
    return "SI" if x else "NO"


# ============================================================================================ principal
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "analyze_2f_out.txt"))
    ap.add_argument("--seeds", default="4242,8888")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    L = []
    P = L.append
    P("# Tarea 2F - tope por estacion en secuencias densas (analyze_2f.py)")
    P("")
    P("Recall de eventos a eventos falsos realizados por ventana <= 0.1 / <= 0.25 (FROC exacta, punto de operacion por "
      "metodo sobre la evaluacion; agrupado = un umbral por metodo para las ventanas agrupadas). Diferencias: punto "
      "[IC 95% jerarquico: celdas como estratos, redes y ventanas; 2000 replicas; punto de operacion re-elegido en cada "
      "replica]; p = (1 + #{d* <= 0}) / 2001.")
    res = {}
    ctl = []
    for s in seeds:
        for arm in ARMS:
            nets = ARMS[arm]["nets"]
            D = load_arm(s, arm)
            CC = counts(s, arm, len(nets))
            res[(s, arm)] = (D, CC)
            kern = {tuple(tuple(x) for x in m_["kernel"]["blas"]) for W in D.values() for m_ in W["meta"]}
            env = {m_["kernel"]["OPENBLAS_CORETYPE"] for W in D.values() for m_ in W["meta"]}
            tags = {(m_["tag"], m_["kmax"]) for W in D.values() for m_ in W["meta"]}
            kset = sorted({int(k) for W in D.values() for k in W["k"]})
            nwin = sum(len(W["k"]) for W in D.values())
            nev = int(sum(W["k"].sum() for W in D.values()))
            ctl.append(f"  semilla {s}, brazo {arm}: {len(nets)} redes x 3 celdas x 2 regimenes, {nwin} ventanas, "
                       f"{nev} eventos, k en {kset}; (TAG, kmax) {sorted(tags)}; OPENBLAS_CORETYPE {sorted(env)}; "
                       f"OpenBLAS cargadas {sorted(kern)}")
    P("")
    P("Controles de datos:")
    L.extend(ctl)

    # --------------------------------------------------------------------------- recall por metodo
    P("")
    P("## Recall por metodo, brazo denso (agrupando las 3 celdas) y brazo estandar (agrupando celdas)")
    ev = {}
    for s in seeds:
        for arm in ARMS:
            D, CC = res[(s, arm)]
            for reg in REGIMES:
                keys = [(reg, c) for c in CELLS]
                pools, C, nev, nw = pooled(D, CC, keys)
                ev[(s, arm, reg)] = evaluate(pools, C, nw) + (nev, nw)
            keys = [(r, c) for r in REGIMES for c in CELLS]
            pools, C, nev, nw = pooled(D, CC, keys)
            ev[(s, arm, "both")] = evaluate(pools, C, nw) + (nev, nw)
            for reg in REGIMES:
                for c in CELLS:
                    pools, C, nev, nw = pooled(D, CC, [(reg, c)])
                    ev[(s, arm, reg, c)] = evaluate(pools, C, nw) + (nev, nw)
    for arm in ARMS:
        for grp in REGIMES + ("both",):
            for s in seeds:
                pt, bs, nev, nw = ev[(s, arm, grp)]
                row = f"  {arm:5s} {grp:6s} {s}: {nw} ventanas, {nev} eventos | "
                row += "  ".join(f"{m} {pt[m][0.1]:.3f}/{pt[m][0.25]:.3f}" for m in METHODS)
                P(row)

    # --------------------------------------------------------------------------- 2F-a
    P("")
    P("## Regla 2F-a (primaria): MT-FW+SC - BOMP-PK+SC a tope 0.25, brazo denso, agrupando las 3 celdas")
    pv = {}
    da = {}
    for reg in REGIMES:
        for s in seeds:
            pt, bs, _nev, _nw = ev[(s, "dense", reg)]
            d = diff(pt, bs, "MT-FW+SC", "BOMP-PK+SC", 0.25)
            da[(reg, s)] = d
            pv[(reg, s)] = d["p_le0"]
    rej, adj = holm(pv, ALPHA_HOLM)
    verdict_a, verdict_a_strict = {}, {}
    for reg in REGIMES:
        for s in seeds:
            d = da[(reg, s)]
            P(f"  {reg:6s} {s}: {fci(d)}  p = {d['p_le0']:.4f}  Holm p aj = {adj[(reg, s)]:.4f} "
              f"rechaza = {yn(rej[(reg, s)])}  (DE bootstrap {d['sd']:.4f}; DE piloto {PILOT_SD[reg]:.4f})")
        ok = [da[(reg, s)]["lo"] > 0 and da[(reg, s)]["point"] >= THR_A - EPS for s in seeds]
        ok_s = [rej[(reg, s)] and da[(reg, s)]["point"] >= THR_A - EPS for s in seeds]
        verdict_a[reg] = all(ok) and len(seeds) == 2
        verdict_a_strict[reg] = all(ok_s) and len(seeds) == 2
        P(f"  -> {reg}: cota inferior > 0 y punto >= {THR_A} por semilla {[yn(x) for x in ok]}: 2F-a "
          f"{'SE CUMPLE' if verdict_a[reg] else 'NO SE CUMPLE'}; lectura estricta (Holm, 4 pruebas, 0.025) "
          f"{[yn(x) for x in ok_s]}: {'SE CUMPLE' if verdict_a_strict[reg] else 'NO SE CUMPLE'}")
    P("  Secundario, tope 0.1 (mismo contraste, descriptivo):")
    for reg in REGIMES:
        for s in seeds:
            pt, bs, _nev, _nw = ev[(s, "dense", reg)]
            d = diff(pt, bs, "MT-FW+SC", "BOMP-PK+SC", 0.1)
            P(f"    {reg:6s} {s}: {fci(d)}  p = {d['p_le0']:.4f}")

    # --------------------------------------------------------------------------- 2F-b
    P("")
    P("## Regla 2F-b (mecanismo), brazo denso, agrupando celdas y regimenes, tope 0.25")
    okb = []
    for s in seeds:
        pt, bs, _nev, _nw = ev[(s, "dense", "both")]
        d1 = diff(pt, bs, "MT-FW+SC", "MT-FW-NOCAP+SC", 0.25)
        d2 = diff(pt, bs, "MT-FW-NOCAP+SC", "BOMP-PK+SC", 0.25)
        c1 = d1["lo"] > 0
        c2 = abs(d2["point"]) <= THR_B + EPS
        okb.append(c1 and c2)
        P(f"  {s}: MT-FW+SC - MT-FW-NOCAP+SC {fci(d1)} (cota inferior > 0: {yn(c1)});  "
          f"MT-FW-NOCAP+SC - BOMP-PK+SC {fci(d2)} (|punto| <= {THR_B}: {yn(c2)})")
        d3 = diff(pt, bs, "MT-FW+SC", "BOMP-PK+SC", 0.25)
        P(f"        (referencia, agrupando regimenes: MT-FW+SC - BOMP-PK+SC {fci(d3)})")
    verdict_b = all(okb) and len(seeds) == 2
    P(f"  -> 2F-b {'SE CUMPLE' if verdict_b else 'NO SE CUMPLE'}")
    P("  Tope 0.1 (descriptivo):")
    for s in seeds:
        pt, bs, _nev, _nw = ev[(s, "dense", "both")]
        P(f"    {s}: MT-FW+SC - MT-FW-NOCAP+SC {fci(diff(pt, bs, 'MT-FW+SC', 'MT-FW-NOCAP+SC', 0.1))};  "
          f"MT-FW-NOCAP+SC - BOMP-PK+SC {fci(diff(pt, bs, 'MT-FW-NOCAP+SC', 'BOMP-PK+SC', 0.1))}")

    # --------------------------------------------------------------------------- 2F-c
    P("")
    P(f"## Regla 2F-c (no dano, brazo estandar): TOST, margen {MARGIN_C}, MT-FW+SC - BOMP-PK+SC a tope 0.25, "
      "agrupando celdas y regimenes")
    okc = []
    for s in seeds:
        pt, bs, nev, nw = ev[(s, "std", "both")]
        d = diff(pt, bs, "MT-FW+SC", "BOMP-PK+SC", 0.25)
        c = d["p_tost_lo"] < ALPHA_TOST and d["p_tost_hi"] < ALPHA_TOST
        okc.append(c)
        P(f"  {s}: {nw} ventanas, {nev} eventos; diferencia {d['point']:+.3f}, IC 90 [{d['lo90']:+.3f}, "
          f"{d['hi90']:+.3f}], IC 95 [{d['lo']:+.3f}, {d['hi']:+.3f}]; p(d <= -0.05) = {d['p_tost_lo']:.4f}, "
          f"p(d >= +0.05) = {d['p_tost_hi']:.4f}: equivalente {yn(c)}")
    verdict_c = all(okc) and len(seeds) == 2
    P(f"  -> 2F-c {'SE CUMPLE' if verdict_c else 'NO SE CUMPLE'}")
    P("  Tope 0.1 y por regimen (descriptivo):")
    for s in seeds:
        pt, bs, _nev, _nw = ev[(s, "std", "both")]
        P(f"    {s} ambos, tope 0.1: {fci(diff(pt, bs, 'MT-FW+SC', 'BOMP-PK+SC', 0.1))}")
        for reg in REGIMES:
            pt, bs, _nev, _nw = ev[(s, "std", reg)]
            P(f"    {s} {reg}: tope 0.1 {fci(diff(pt, bs, 'MT-FW+SC', 'BOMP-PK+SC', 0.1))}; tope 0.25 "
              f"{fci(diff(pt, bs, 'MT-FW+SC', 'BOMP-PK+SC', 0.25))}")

    # --------------------------------------------------------------------------- descriptivo
    P("")
    P("## Descriptivo: diferencias agrupando las 3 celdas (brazo denso), ambos topes")
    for reg in REGIMES + ("both",):
        for s in seeds:
            pt, bs, _nev, _nw = ev[(s, "dense", reg)]
            for a_, b_ in PAIRS:
                P(f"  {reg:6s} {s} {a_} - {b_}: tope 0.1 {fci(diff(pt, bs, a_, b_, 0.1))}; tope 0.25 "
                  f"{fci(diff(pt, bs, a_, b_, 0.25))}")
    P("")
    P("## Descriptivo: por celda (brazo denso), recall por metodo y diferencias con IC 95 por celda")
    npos = {0.1: 0, 0.25: 0}
    ncase = 0
    for reg in REGIMES:
        for s in seeds:
            for c in CELLS:
                pt, bs, nev, nw = ev[(s, "dense", reg, c)]
                P(f"  {reg:6s} {s} {LAB[c]:12s} ({nw} ventanas, {nev} eventos): " +
                  "  ".join(f"{m} {pt[m][0.1]:.3f}/{pt[m][0.25]:.3f}" for m in METHODS))
                ncase += 1
                for a_, b_ in PAIRS:
                    d1, d2 = diff(pt, bs, a_, b_, 0.1), diff(pt, bs, a_, b_, 0.25)
                    P(f"      {a_} - {b_}: tope 0.1 {fci(d1)}; tope 0.25 {fci(d2)}")
                    if (a_, b_) == ("MT-FW+SC", "BOMP-PK+SC"):
                        npos[0.1] += d1["point"] > 0
                        npos[0.25] += d2["point"] > 0
    P(f"  MT-FW+SC - BOMP-PK+SC positivo por celda: tope 0.1 {npos[0.1]}/{ncase}, tope 0.25 {npos[0.25]}/{ncase}")
    P("")
    P("## Descriptivo: brazo estandar por celda (MT-FW+SC - BOMP-PK+SC)")
    for reg in REGIMES:
        for s in seeds:
            for c in CELLS:
                pt, bs, nev, nw = ev[(s, "std", reg, c)]
                P(f"  {reg:6s} {s} {LAB[c]:12s} ({nw} ventanas, {nev} eventos): MT-FW+SC "
                  f"{pt['MT-FW+SC'][0.1]:.3f}/{pt['MT-FW+SC'][0.25]:.3f}, BOMP-PK+SC {pt['BOMP-PK+SC'][0.1]:.3f}/"
                  f"{pt['BOMP-PK+SC'][0.25]:.3f}; dif. tope 0.1 {fci(diff(pt, bs, 'MT-FW+SC', 'BOMP-PK+SC', 0.1))}; "
                  f"tope 0.25 {fci(diff(pt, bs, 'MT-FW+SC', 'BOMP-PK+SC', 0.25))}")

    # --------------------------------------------------------------------------- tiempos
    P("")
    P("## Tiempos (s por ventana, media; un proceso por nucleo, 2 procesos en paralelo)")
    for arm in ARMS:
        for s in seeds:
            D, _CC = res[(s, arm)]
            for reg in REGIMES:
                ws = [(W["sec"], W["meta"]) for (r_, _c), W in D.items() if r_ == reg]
                row = "  ".join(f"{m} {np.mean([x for sec, _m in ws for x in sec[m]]):.3f}" for m in METHODS)
                tset = np.mean([mm["t_setup"] for _sec, meta in ws for mm in meta])
                P(f"  {arm:5s} {s} {reg:6s}: {row}; preparacion por red {tset:.1f} s")

    # --------------------------------------------------------------------------- veredictos
    P("")
    P("## Veredictos")
    for reg in REGIMES:
        P(f"  2F-a {reg}: {'SE CUMPLE' if verdict_a[reg] else 'NO SE CUMPLE'}  (estricta Holm: "
          f"{'SE CUMPLE' if verdict_a_strict[reg] else 'NO SE CUMPLE'})")
    P(f"  2F-b: {'SE CUMPLE' if verdict_b else 'NO SE CUMPLE'}")
    P(f"  2F-c: {'SE CUMPLE' if verdict_c else 'NO SE CUMPLE'}")
    any_a = any(verdict_a.values())
    if any_a and verdict_b:
        interp = ("2F-a y 2F-b se cumplen -> el tope por estacion es una contribucion confirmada para secuencias densas "
                  f"(en los regimenes donde 2F-a se cumple: {[r for r in REGIMES if verdict_a[r]]}).")
    elif any_a:
        interp = ("2F-a se cumple y 2F-b no -> la ganancia es real pero no se debe al tope "
                  f"(regimenes con 2F-a: {[r for r in REGIMES if verdict_a[r]]}).")
    else:
        interp = "2F-a falla -> segun el registro, el efecto piloto era un artefacto de seleccion."
    P(f"  Interpretacion registrada: {interp}")
    txt = "\n".join(L) + "\n"
    with open(a.out, "w") as fh:
        fh.write(txt)
    print(txt)


if __name__ == "__main__":
    main()
