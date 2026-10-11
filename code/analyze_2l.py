"""analyze_2l.py - reglas 2L-a, 2L-b, 2L-c y 2L-d de la Tarea 2L (PREREGISTRO_1j_1k.txt) por semilla, el secundario
registrado (familia BOMP) y los puntos descriptivos registrados.  Escrito y congelado (task2l_code_hashes.txt) antes de
correr codigo alguno sobre las semillas 7171 y 9393.

Datos: task2l_parts/ (task2l_run.py).  Diseno estandar: redes j = 0..159 por celda y regimen, 5 ventanas (kmax 5);
diseno denso: j = 0..99, 5 ventanas (kmax 7); familia BOMP: estandar, j = 0..59; brazo oraculo: j = 0..29.  Con --budget
(regla de presupuesto del registro): estandar j = 0..119 y densas j = 0..79.
Reglas de acierto (c10_lib.label, sobre los atomos guardados): LOC+T = epicentro a <= 5 km Y |origen del atomo - origen
verdadero| <= 0.5 s, el evento mas cercano que cumpla ambas (primaria); LOC = evento mas cercano a <= 5 km (la regla
registrada de 2A-2I).  Registro de una ventana: c10_lib.record (minimo corrido del estadistico de seleccion "s").
Endpoint: recall de eventos a eventos falsos realizados por ventana <= 0.25 (primario) y <= 0.1 (secundario), FROC
exacta (analyze_2f.Pool, la de analyze_2i.py: todos los puntajes del registro como umbrales), punto de operacion elegido
por metodo sobre la propia evaluacion; al agrupar regimenes y celdas la FROC es la del conjunto agrupado (un umbral por
metodo); en el analisis por regimen (o por celda) la FROC es la de ese subconjunto.
Bootstrap jerarquico registrado (el de analyze_2i.py): 2000 replicas; estratos regimen x celda (6); en cada estrato las
redes con reemplazo y, dentro de cada red elegida, sus ventanas con reemplazo; las mismas remuestras para todos los
metodos (diferencias pareadas); punto de operacion re-elegido en cada replica.  Generador: default_rng([semilla, TAG
del tipo de ventana (6101 estandar | 6102 densas), indice de regimen, indice de celda, numero de redes]).
Cotas: cota inferior unilateral 97.5 = percentil 2.5 de las diferencias bootstrap; cota inferior unilateral 95 =
percentil 5; IC 95 = percentiles 2.5 y 97.5.  Comparacion de la estimacion puntual: >= umbral - 1e-9.
Contraste de las reglas: MT@S2+T3 - MT@S0+T3 a tope 0.25.
  2L-a superioridad, ventanas estandar, LOC+T, agrupado: cota inferior unilateral 97.5 > 0 y punto >= +0.03.
  2L-b superioridad, ventanas densas, LOC+T, agrupado: cota inferior unilateral 97.5 > 0.
  2L-c no inferioridad, margen 0.02, por regimen (thin, strong), LOC+T, estandar y densas: cota inferior unilateral 95
       > -0.02 en las cuatro combinaciones.
  2L-d no inferioridad, margen 0.02, regla LOC, agrupado, estandar y densas: cota inferior unilateral 95 > -0.02 en las
       dos.
  Cada regla debe cumplirse en AMBAS semillas.
Secundario registrado: BOMP@S2+T3 - BOMP@S0+T3, estandar, LOC+T, 60 redes: cota inferior unilateral 97.5 > 0.
Descriptivo (sin regla): contrastes totales contra el detector registrado (MT@S2+T3 - MT@S0) y MT@S0+T3 - MT@S0; por
  celda; tope 0.1; error rms de los terminos contra los retardos verdaderos (mediana quitada; c11_lib.err) y
  desplazamiento comun; error epicentral y de tiempo de origen de los aciertos (todos los atomos, regla LOC; IC por
  bootstrap de redes dentro de los estratos, 400 replicas); fraccion de la ganancia del oraculo capturada (30 redes;
  forma T3 = (S2+T3 - S0+T3) / (ORA+T3 - S0+T3), forma total = (S2+T3 - S0) / (ORA+T3 - S0); cociente por replica);
  DE bootstrap frente a la planificada; tiempos.

    python analyze_2l.py [--budget] [--out analyze_2l_out.txt]
    python analyze_2l.py --dry --out ...      (ensayo con las partes de control de la semilla 999)
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse
import json
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXPL = os.path.join(HERE, "explore_contrib")
for p_ in (EXPL, HERE):
    if p_ not in sys.path:
        sys.path.insert(0, p_)

from analyze_2f import Pool  # noqa: E402
import c10_lib as C10  # noqa: E402

PARTS = os.path.join(HERE, "task2l_parts")
B = 2000
B_GEO = 400
CAPS = (0.1, 0.25)
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
LAB = {"ladder_4": "4 st.", "ladder_6": "6 st.", "unilateral_c40": "one-sided 8"}
REGIMES = ("thin", "strong")
KINDS = ("std", "dense")
KNAME = {"std": "estandar", "dense": "densas"}
TAGK = {"std": 6101, "dense": 6102}
SEEDS = (7171, 9393)
N_WIN = 5
S0H, S0T, S2T, ORT = "MT@S0", "MT@S0+T3", "MT@S2+T3", "MT@ORA+T3"
B0H, B0T, B2T, BOT = "BOMP@S0", "BOMP@S0+T3", "BOMP@S2+T3", "BOMP@ORA+T3"
MT3 = (S0H, S0T, S2T)
BOMP3 = (B0H, B0T, B2T)
RULES = (("LOC+T", 0.5), ("LOC", None))
CFG = {"std": 160, "dense": 100, "bomp": 60, "ora": 30, "prefix": "main", "nets": None}
CFG_BUDGET = {"std": 120, "dense": 80, "bomp": 60, "ora": 30, "prefix": "main", "nets": None}
CFG_DRY = {"std": 3, "dense": 3, "bomp": 3, "ora": 3, "prefix": "ctrl", "nets": (210, 219, 228)}
EPS = 1e-9
THR_A = 0.03
MARGIN = 0.02
# DE bootstrap del piloto C11 (semilla 999, 30 redes por celda y regimen; REPORT_C11.txt, secciones 6, 7 y 12)
PILOT_SD = {("std", "LOC+T"): 0.0282, ("dense", "LOC+T"): 0.0151, ("std", "LOC"): 0.0244, ("dense", "LOC"): 0.0179,
            "bomp": 0.0304}
PILOT_N = 30
ALL = [(r, c) for r in REGIMES for c in CELLS]


# ============================================================================================ datos
def load(cfg, seed, parts=PARTS):
    """D[kind][(regime, cell)] = {"k", "truth", "taus", "atoms": {m: [...]}, "sec": {m: [...]}} con las ventanas en el
    orden de las redes y de las ventanas (los metodos de un subconjunto de redes ocupan las primeras ventanas), y
    N[(regime, cell)] = lista de metadatos por red.  (None, None) si falta alguna parte."""
    nets = cfg["nets"] if cfg["nets"] is not None else tuple(range(max(cfg["std"], cfg["dense"])))
    D = {kind: {} for kind in KINDS}
    N = {}
    for reg in REGIMES:
        for cell in CELLS:
            N[(reg, cell)] = []
            for kind in KINDS:
                D[kind][(reg, cell)] = {"k": [], "truth": [], "taus": [], "atoms": {}, "sec": {}}
            for ji, j in enumerate(nets):
                path = os.path.join(parts, f"{cfg['prefix']}_seed{seed}_{reg}_{cell}_net{j}.json")
                if not os.path.exists(path):
                    return None, None
                with open(path) as fh:
                    e = json.load(fh)
                assert (e["seed"], e["regime"], e["cell"], e["net"], e["n_win"]) == (seed, reg, cell, j, N_WIN)
                N[(reg, cell)].append({x: e[x] for x in ("tags", "kmax", "err", "off", "sigT3", "S2", "t_setup", "t_s2",
                                                         "t_total", "kernel", "host_time_utc", "cat_drawn")})
                for kind in KINDS:
                    if ji >= cfg[kind]:
                        continue
                    rows = e["windows"][kind]
                    W = D[kind][(reg, cell)]
                    assert len(rows) == N_WIN
                    for w in rows:
                        W["k"].append(w["k"])
                        W["truth"].append(w["truth"])
                        W["taus"].append(w["taus"])
                        for m, a in w["atoms"].items():
                            W["atoms"].setdefault(m, []).append(a)
                            W["sec"].setdefault(m, []).append(w["sec"][m])
    for kind in KINDS:
        for W in D[kind].values():
            W["k"] = np.asarray(W["k"], float)
    return D, N


def atoms_of(packed):
    return [{"x": np.array(a[:3], float), "tau": a[3], "s": a[4]} for a in packed]


def build(Dk, tol_t):
    """registros {clave: {metodo: [registro por ventana]}} con la regla de acierto dada (c10_lib.label / record)."""
    out = {}
    for key, W in Dk.items():
        out[key] = {}
        for m, lst in W["atoms"].items():
            recs = []
            for wi, packed in enumerate(lst):
                at = atoms_of(packed)
                lab = C10.label(at, [np.array(t, float) for t in W["truth"][wi]], W["taus"][wi], tol_t)
                recs.append(C10.record(at, lab, "s"))
            out[key][m] = recs
    return out


# ============================================================================================ bootstrap
def counts(seed, kind, nnet, B_=B):
    """{(regime, cell): C (B, nnet * N_WIN)}: estratos regimen x celda, redes y ventanas remuestreadas de forma
    independiente en cada estrato (analyze_2i.counts)."""
    out = {}
    rows = np.arange(B_)[:, None, None]
    for ci, cell in enumerate(CELLS):
        for ri, reg in enumerate(REGIMES):
            rng = np.random.default_rng([seed, TAGK[kind], ri, ci, nnet])
            pick = rng.integers(0, nnet, (B_, nnet))
            idx = pick[:, :, None] * N_WIN + rng.integers(0, N_WIN, (B_, nnet, N_WIN))
            C = np.zeros((B_, nnet * N_WIN))
            np.add.at(C, (rows, idx), 1.0)
            out[(reg, cell)] = C
    return out


def evaluate(Dk, recs, CC, keys, methods, nnet):
    """puntos y replicas por metodo para los estratos keys agrupados, con las primeras nnet redes de cada estrato."""
    n = nnet * N_WIN
    ks = np.concatenate([Dk[k_]["k"][:n] for k_ in keys])
    C = np.concatenate([CC[k_] for k_ in keys], axis=1)
    assert C.shape[1] == len(ks)
    one = np.ones((1, len(ks)))
    pt, bs = {}, {}
    for m in methods:
        r = [x for k_ in keys for x in recs[k_][m][:n]]
        assert len(r) == len(ks)
        p = Pool(r, ks)
        v = p.recall(one)
        pt[m] = {cap: float(v[cap][0]) for cap in CAPS}
        bs[m] = p.recall(C)
    return pt, bs, int(ks.sum()), len(ks)


def diff(ev, a, b, cap):
    pt, bs = ev[0], ev[1]
    d = bs[a][cap] - bs[b][cap]
    p = pt[a][cap] - pt[b][cap]
    lo, lo95, hi = np.percentile(d, [2.5, 5.0, 97.5])
    n = len(d) + 1
    return {"point": p, "lo": float(lo), "hi": float(hi), "lo95": float(lo95), "sd": float(d.std()),
            "p_le0": float((1 + np.sum(d <= 0)) / n), "p_le_m": float((1 + np.sum(d <= -MARGIN)) / n)}


def ratio(ev, x, base, ora, cap):
    """fraccion de la ganancia del oraculo: punto, IC 95 del cociente por replica, y fraccion de replicas con ganancia
    del oraculo <= 0."""
    pt, bs = ev[0], ev[1]
    den = pt[ora][cap] - pt[base][cap]
    num = pt[x][cap] - pt[base][cap]
    dd = bs[ora][cap] - bs[base][cap]
    r = (bs[x][cap] - bs[base][cap]) / np.where(np.abs(dd) < 1e-9, np.nan, dd)
    lo, hi = np.nanpercentile(r, [2.5, 97.5])
    return (num / den if abs(den) > 1e-9 else float("nan")), float(lo), float(hi), float(np.mean(dd <= 0))


def fci(d):
    return f"{d['point']:+.3f} [{d['lo']:+.3f}, {d['hi']:+.3f}]"


def rci(ev, m, cap):
    lo, hi = np.percentile(ev[1][m][cap], [2.5, 97.5])
    return f"{ev[0][m][cap]:.3f} [{lo:.3f}, {hi:.3f}]"


def yn(x):
    return "SI" if x else "NO"


def cumple(x):
    return "SE CUMPLE" if x else "NO SE CUMPLE"


def keys_of(grp):
    if grp == "both":
        return ALL
    if grp in REGIMES:
        return [(grp, c) for c in CELLS]
    if grp in CELLS:
        return [(r, grp) for r in REGIMES]
    return [grp]


# ============================================================================================ aciertos
def geometry(Dk, m, nnet, seed, kind):
    """atomos acertados (todos los atomos, regla LOC): mediana del error epicentral (km) y de |error de origen| (s),
    fraccion con |error de origen| <= 0.5 s, fraccion de eventos tocados por algun atomo; IC 95 por bootstrap de redes
    dentro de los estratos regimen x celda (B_GEO replicas)."""
    nets = []                                            # por red: (errores (n, 2), eventos tocados, eventos)
    strata = []
    for si, key in enumerate(ALL):
        W = Dk[key]
        for jn in range(nnet):
            e, nh, nk = [], 0, 0
            for wi in range(jn * N_WIN, (jn + 1) * N_WIN):
                at = atoms_of(W["atoms"][m][wi])
                tr = [np.array(t, float) for t in W["truth"][wi]]
                lab = C10.label(at, tr, W["taus"][wi], None)
                e += [(float(np.linalg.norm(a["x"][:2] - tr[h][:2])), abs(a["tau"] - W["taus"][wi][h]))
                      for a, h in zip(at, lab) if h >= 0]
                nh += len(set(h for h in lab if h >= 0))
                nk += len(tr)
            nets.append((np.array(e, float).reshape(-1, 2), nh, nk))
            strata.append(si)
    strata = np.array(strata)

    def st(ids):
        e = np.concatenate([nets[i][0] for i in ids])
        if not len(e):
            return (np.nan, np.nan, np.nan, 0.0)
        return (float(np.median(e[:, 0])), float(np.median(e[:, 1])), float(np.mean(e[:, 1] <= 0.5)),
                sum(nets[i][1] for i in ids) / max(1, sum(nets[i][2] for i in ids)))
    pt = st(range(len(nets)))
    rng = np.random.default_rng([seed, TAGK[kind], 777, nnet])
    bs = []
    for _ in range(B_GEO):
        ids = []
        for si in range(len(ALL)):
            g = np.where(strata == si)[0]
            ids += g[rng.integers(0, len(g), len(g))].tolist()
        bs.append(st(ids))
    lo, hi = np.nanpercentile(np.array(bs), [2.5, 97.5], axis=0)
    return pt, lo, hi


# ============================================================================================ principal
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "analyze_2l_out.txt"))
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--budget", action="store_true")
    ap.add_argument("--boot", type=int, default=B)
    a = ap.parse_args()
    cfg = CFG_DRY if a.dry else (CFG_BUDGET if a.budget else CFG)
    seeds = (999,) if a.dry else SEEDS
    nb = a.boot
    NN = {"std": cfg["std"], "dense": cfg["dense"]}
    L_ = []
    P = L_.append
    P("# Tarea 2L - terminos de estacion por JHD sin picks (S2) en el detector refinado (T3) (analyze_2l.py)" +
      ("   *** ENSAYO, semilla 999, partes de control ***" if a.dry else "") +
      ("   [regla de presupuesto: estandar 120 redes, densas 80]" if a.budget else ""))
    P("")
    P("Recall de eventos a eventos falsos realizados por ventana <= 0.1 / <= 0.25 (FROC exacta, punto de operacion por "
      "metodo sobre la evaluacion; agrupado = un umbral por metodo). Diferencias pareadas: punto [IC 95% jerarquico: "
      f"estratos regimen x celda, redes y ventanas; {nb} replicas; punto de operacion re-elegido en cada replica]. Cota "
      "inferior unilateral 97.5 = percentil 2.5; cota inferior unilateral 95 = percentil 5. Reglas de acierto: LOC+T "
      "(epicentro a 5 km y |error de origen| <= 0.5 s; primaria) y LOC (solo epicentro).")
    P(f"Redes por celda y regimen: estandar {cfg['std']}, densas {cfg['dense']}, familia BOMP (estandar) {cfg['bomp']}, "
      f"oraculo {cfg['ora']}; {N_WIN} ventanas por red.")

    data, meta, recs, CC = {}, {}, {}, {}
    P("")
    P("Controles de datos:")
    for s in seeds:
        D, N = load(cfg, s)
        if D is None:
            P(f"  semilla {s}: PARTES AUSENTES")
            continue
        data[s], meta[s] = D, N
        mt = [x for v in N.values() for x in v]
        kern = sorted({tuple(tuple(x) for x in m_["kernel"]["blas"]) for m_ in mt})
        env = sorted({str(m_["kernel"]["OPENBLAS_CORETYPE"]) for m_ in mt})
        tags = sorted({(m_["tags"]["std"], m_["tags"]["dense"]) for m_ in mt})
        hts = sorted(m_["host_time_utc"] for m_ in mt)
        for kind in KINDS:
            Dk = D[kind]
            nwin = sum(len(W["k"]) for W in Dk.values())
            nev = int(sum(W["k"].sum() for W in Dk.values()))
            kset = sorted({int(k) for W in Dk.values() for k in W["k"]})
            nm = {m: sum(len(W["atoms"].get(m, [])) for W in Dk.values()) for m in
                  sorted({m for W in Dk.values() for m in W["atoms"]})}
            P(f"  semilla {s}, {KNAME[kind]}: {NN[kind]} redes x 3 celdas x 2 regimenes, {nwin} ventanas, {nev} eventos, "
              f"k en {kset}; ventanas por metodo {nm}")
            for rule, tol in RULES:
                recs[(s, kind, rule)] = build(Dk, tol)
            for nnet in sorted({NN[kind], cfg["ora"]} | ({cfg["bomp"]} if kind == "std" else set())):
                CC[(s, kind, nnet)] = counts(s, kind, nnet, nb)
        P(f"      (TAG estandar, TAG densas) {tags}; OPENBLAS_CORETYPE {env}; OpenBLAS cargadas {kern}; {len(mt)} partes "
          f"escritas entre {hts[0]} y {hts[-1]} UTC")
    if any(s not in data for s in seeds):
        P("")
        P("CORRIDA INCOMPLETA: no se evaluan las reglas.")
        open(a.out, "w").write("\n".join(L_) + "\n")
        print("\n".join(L_))
        return

    memo = {}

    def E(s, kind, grp, rule, sub="all"):
        """sub: all (MT, todas las redes) | bomp (primeras redes de la familia BOMP, 6 metodos) | ora (redes del
        oraculo: MT + ORA y, en estandar, BOMP + ORA)."""
        key = (s, kind, grp, rule, sub)
        if key not in memo:
            if sub == "all":
                nnet, ms = NN[kind], MT3
            elif sub == "bomp":
                nnet, ms = cfg["bomp"], MT3 + BOMP3
            else:
                nnet, ms = cfg["ora"], MT3 + (ORT,) + ((BOMP3 + (BOT,)) if kind == "std" else ())
            memo[key] = evaluate(data[s][kind], recs[(s, kind, rule)], CC[(s, kind, nnet)], keys_of(grp), ms, nnet)
        return memo[key]

    # --------------------------------------------------------------------------- recall por metodo
    P("")
    P("## Recall por metodo (familia MT, todas las redes), tope 0.1 / tope 0.25, con IC 95")
    for kind in KINDS:
        for rule, _tol in RULES:
            for grp in ("both",) + REGIMES:
                for s in seeds:
                    e = E(s, kind, grp, rule)
                    P(f"  {KNAME[kind]:8s} {rule:5s} {grp:6s} {s}: {e[3]} ventanas, {e[2]} eventos")
                    for m in MT3:
                        P(f"      {m:11s} {rci(e, m, 0.1)}   {rci(e, m, 0.25)}")

    # --------------------------------------------------------------------------- reglas
    two = len(seeds) == 2
    P("")
    P("## Regla 2L-a (superioridad, ventanas estandar, LOC+T, agrupado): MT@S2+T3 - MT@S0+T3 a tope 0.25")
    ok_a = []
    for s in seeds:
        d = diff(E(s, "std", "both", "LOC+T"), S2T, S0T, 0.25)
        c1, c2 = d["lo"] > 0, d["point"] >= THR_A - EPS
        ok_a.append(c1 and c2)
        P(f"  {s}: {fci(d)}; cota inferior unilateral 97.5 = {d['lo']:+.4f} > 0: {yn(c1)}; punto {d['point']:+.4f} >= "
          f"+{THR_A}: {yn(c2)}; p(d* <= 0) = {d['p_le0']:.4f}; DE bootstrap {d['sd']:.4f} -> {cumple(c1 and c2)}")
    va = all(ok_a) and two
    P(f"  -> 2L-a {cumple(va)} (por semilla: {[yn(x) for x in ok_a]})")
    P("")
    P("## Regla 2L-b (superioridad, ventanas densas, LOC+T, agrupado): MT@S2+T3 - MT@S0+T3 a tope 0.25")
    ok_b = []
    for s in seeds:
        d = diff(E(s, "dense", "both", "LOC+T"), S2T, S0T, 0.25)
        c1 = d["lo"] > 0
        ok_b.append(c1)
        P(f"  {s}: {fci(d)}; cota inferior unilateral 97.5 = {d['lo']:+.4f} > 0: {yn(c1)}; p(d* <= 0) = "
          f"{d['p_le0']:.4f}; DE bootstrap {d['sd']:.4f} -> {cumple(c1)}")
    vb = all(ok_b) and two
    P(f"  -> 2L-b {cumple(vb)} (por semilla: {[yn(x) for x in ok_b]})")
    P("")
    P(f"## Regla 2L-c (no inferioridad, margen {MARGIN}, por regimen, LOC+T, estandar y densas): MT@S2+T3 - MT@S0+T3 a "
      "tope 0.25, cota inferior unilateral 95")
    ok_c = []
    for s in seeds:
        oks = []
        for kind in KINDS:
            for reg in REGIMES:
                d = diff(E(s, kind, reg, "LOC+T"), S2T, S0T, 0.25)
                c1 = d["lo95"] > -MARGIN
                oks.append(c1)
                P(f"  {s} {KNAME[kind]:8s} {reg:6s}: {fci(d)}; cota inferior unilateral 95 = {d['lo95']:+.4f} > -{MARGIN}: "
                  f"{yn(c1)}; p(d* <= -{MARGIN}) = {d['p_le_m']:.4f}; DE bootstrap {d['sd']:.4f}")
        ok_c.append(all(oks))
        P(f"  {s}: -> {cumple(all(oks))}")
    vc = all(ok_c) and two
    P(f"  -> 2L-c {cumple(vc)} (por semilla: {[yn(x) for x in ok_c]})")
    P("")
    P(f"## Regla 2L-d (no inferioridad, margen {MARGIN}, regla LOC, agrupado, estandar y densas): MT@S2+T3 - MT@S0+T3 a "
      "tope 0.25, cota inferior unilateral 95")
    ok_d = []
    for s in seeds:
        oks = []
        for kind in KINDS:
            d = diff(E(s, kind, "both", "LOC"), S2T, S0T, 0.25)
            c1 = d["lo95"] > -MARGIN
            oks.append(c1)
            P(f"  {s} {KNAME[kind]:8s}: {fci(d)}; cota inferior unilateral 95 = {d['lo95']:+.4f} > -{MARGIN}: {yn(c1)}; "
              f"p(d* <= -{MARGIN}) = {d['p_le_m']:.4f}; DE bootstrap {d['sd']:.4f}")
        ok_d.append(all(oks))
        P(f"  {s}: -> {cumple(all(oks))}")
    vd = all(ok_d) and two
    P(f"  -> 2L-d {cumple(vd)} (por semilla: {[yn(x) for x in ok_d]})")

    P("")
    P(f"## Secundario registrado: BOMP@S2+T3 - BOMP@S0+T3, estandar, LOC+T, {cfg['bomp']} redes, tope 0.25 (cota "
      "inferior unilateral 97.5 > 0)")
    ok_s = []
    for s in seeds:
        e = E(s, "std", "both", "LOC+T", "bomp")
        d = diff(e, B2T, B0T, 0.25)
        ok_s.append(d["lo"] > 0)
        P(f"  {s}: {e[3]} ventanas, {e[2]} eventos; recall " + "  ".join(f"{m} {rci(e, m, 0.25)}" for m in BOMP3))
        P(f"      {fci(d)}; cota inferior unilateral 97.5 = {d['lo']:+.4f} > 0: {yn(d['lo'] > 0)}; p(d* <= 0) = "
          f"{d['p_le0']:.4f}; DE bootstrap {d['sd']:.4f}")
    vs = all(ok_s) and two
    P(f"  -> secundario BOMP {cumple(vs)} (por semilla: {[yn(x) for x in ok_s]})")

    P("")
    P("## Veredictos")
    P(f"  2L-a: {cumple(va)}    2L-b: {cumple(vb)}    2L-c: {cumple(vc)}    2L-d: {cumple(vd)}    secundario BOMP: "
      f"{cumple(vs)}")
    if va and vb and vc and vd:
        interp = ("a, b, c y d se cumplen -> contribucion confirmada (calibracion mas refinamiento, para deteccion "
                  "consistente en tiempo de origen).")
    elif not (va and vb):
        interp = ("a o b fallan -> segun el registro, el piloto fue un artefacto de seleccion." +
                  ("" if (vc and vd) else " Ademas falla " + " y ".join(x for x, v in (("c", vc), ("d", vd)) if not v)
                   + "."))
    else:
        interp = ("a y b se cumplen y " + " y ".join(x for x, v in (("c", vc), ("d", vd)) if not v) + " falla -> la "
                  "ganancia existe pero con un costo, y se reporta asi.")
    P(f"  Interpretacion registrada: {interp}")

    # --------------------------------------------------------------------------- descriptivo: contrastes
    P("")
    P("## Descriptivo: contrastes de la familia MT, tope 0.1 y tope 0.25 (IC 95), agrupado y por regimen, ambas reglas")
    P("   terminos = MT@S2+T3 - MT@S0+T3 (el contraste de las reglas); total = MT@S2+T3 - MT@S0 (contra el detector "
      "registrado); T3 = MT@S0+T3 - MT@S0")
    for kind in KINDS:
        for rule, _tol in RULES:
            for grp in ("both",) + REGIMES:
                for s in seeds:
                    e = E(s, kind, grp, rule)
                    for lab, x, y in (("terminos", S2T, S0T), ("total", S2T, S0H), ("T3", S0T, S0H)):
                        d1, d2 = diff(e, x, y, 0.1), diff(e, x, y, 0.25)
                        P(f"  {KNAME[kind]:8s} {rule:5s} {grp:6s} {s} {lab:8s}: tope 0.1 {fci(d1)}; tope 0.25 {fci(d2)}; "
                          f"cota 95 a 0.25 {d2['lo95']:+.3f}")
    P("")
    P("## Descriptivo: por celda (regimenes agrupados) y por celda y regimen, tope 0.25 (IC 95)")
    for kind in KINDS:
        for rule, _tol in RULES:
            for s in seeds:
                for c in CELLS:
                    e = E(s, kind, c, rule)
                    P(f"  {KNAME[kind]:8s} {rule:5s} {s} {LAB[c]:12s} ambos  ({e[3]} v., {e[2]} ev.): " +
                      "  ".join(f"{m} {e[0][m][0.25]:.3f}" for m in MT3) +
                      f" | terminos {fci(diff(e, S2T, S0T, 0.25))}; total {fci(diff(e, S2T, S0H, 0.25))}")
                    for reg in REGIMES:
                        e = E(s, kind, (reg, c), rule)
                        P(f"  {KNAME[kind]:8s} {rule:5s} {s} {LAB[c]:12s} {reg:6s} ({e[3]} v., {e[2]} ev.): " +
                          "  ".join(f"{m} {e[0][m][0.25]:.3f}" for m in MT3) +
                          f" | terminos {fci(diff(e, S2T, S0T, 0.25))}; total {fci(diff(e, S2T, S0H, 0.25))}")
    P("")
    P(f"## Descriptivo: familia BOMP (estandar, {cfg['bomp']} redes) y familia MT en las mismas redes; tope 0.1 / 0.25")
    for rule, _tol in RULES:
        for grp in ("both",) + REGIMES:
            for s in seeds:
                e = E(s, "std", grp, rule, "bomp")
                P(f"  {rule:5s} {grp:6s} {s} ({e[3]} v., {e[2]} ev.): recall 0.1/0.25 " +
                  "  ".join(f"{m} {e[0][m][0.1]:.3f}/{e[0][m][0.25]:.3f}" for m in MT3 + BOMP3))
                for lab, x, y in (("BOMP terminos", B2T, B0T), ("BOMP total", B2T, B0H), ("BOMP T3", B0T, B0H),
                                  ("MT terminos", S2T, S0T)):
                    d1, d2 = diff(e, x, y, 0.1), diff(e, x, y, 0.25)
                    P(f"      {lab:13s}: tope 0.1 {fci(d1)}; tope 0.25 {fci(d2)}")

    # --------------------------------------------------------------------------- oraculo
    P("")
    P(f"## Descriptivo: brazo oraculo ({cfg['ora']} redes por celda y regimen) y fraccion de la ganancia capturada")
    P("   forma T3 = (S2+T3 - S0+T3) / (ORA+T3 - S0+T3); forma total = (S2+T3 - S0) / (ORA+T3 - S0); [IC 95 del cociente "
      "por replica]; (fraccion de replicas con ganancia del oraculo <= 0)")
    for kind in KINDS:
        fams = (("MT", S0H, S0T, S2T, ORT),) + ((("BOMP", B0H, B0T, B2T, BOT),) if kind == "std" else ())
        for rule, _tol in RULES:
            for s in seeds:
                e = E(s, kind, "both", rule, "ora")
                P(f"  {KNAME[kind]:8s} {rule:5s} {s} ({e[3]} v., {e[2]} ev.)")
                for fam, h0, t0, t2, ora in fams:
                    P(f"      {fam}: recall 0.1/0.25 " + "  ".join(f"{m} {e[0][m][0.1]:.3f}/{e[0][m][0.25]:.3f}"
                                                                 for m in (h0, t0, t2, ora)))
                    for cap in CAPS:
                        dS, dO, dG = diff(e, t2, t0, cap), diff(e, ora, t0, cap), diff(e, ora, t2, cap)
                        f1, f2 = ratio(e, t2, t0, ora, cap), ratio(e, t2, h0, ora, cap)
                        P(f"        tope {cap}: S2+T3 - S0+T3 {fci(dS)}; ORA+T3 - S0+T3 {fci(dO)}; ORA+T3 - S2+T3 {fci(dG)}; "
                          f"fraccion T3 {f1[0]:+.2f} [{f1[1]:+.2f}, {f1[2]:+.2f}] ({f1[3]:.3f}); total {f2[0]:+.2f} "
                          f"[{f2[1]:+.2f}, {f2[2]:+.2f}] ({f2[3]:.3f})")

    # --------------------------------------------------------------------------- terminos de estacion
    P("")
    P("## Descriptivo: terminos de estacion; error rms contra los retardos verdaderos, mediana quitada (s): mediana "
      "sobre redes [cuartiles], media; 'mejor' = fraccion de redes con error menor que S0; |desplazamiento| = mediana "
      "de |mediana(termino - retardo)|")
    for s in seeds:
        for reg in REGIMES:
            for tk in ("S0", "S2"):
                row = f"  {s} {reg:6s} {tk}:"
                allv, all0, offs = [], [], []
                for c in CELLS:
                    v = [n_["err"][tk] for n_ in meta[s][(reg, c)]]
                    row += f" {LAB[c]} {np.median(v):.3f}"
                    allv += v
                    all0 += [n_["err"]["S0"] for n_ in meta[s][(reg, c)]]
                    offs += [abs(n_["off"][tk]) for n_ in meta[s][(reg, c)]]
                allv, all0 = np.array(allv), np.array(all0)
                row += (f" | todas ({len(allv)}) {np.median(allv):.3f} [{np.percentile(allv, 25):.3f}, "
                        f"{np.percentile(allv, 75):.3f}] media {allv.mean():.3f}; mejor {np.mean(allv < all0):.2f}; "
                        f"|desplazamiento| {np.median(offs):.3f}")
                P(row)
            ms = [n_ for c in CELLS for n_ in meta[s][(reg, c)]]
            P(f"  {s} {reg:6s} eventos de catalogo relocalizados por S2: media {np.mean([n_['S2']['n_ev'] for n_ in ms]):.1f} "
              f"de 20; escala residual final de S2 (mediana) {np.median([n_['S2']['sig'] for n_ in ms]):.3f} s; nivel de "
              "tiempo de T3 (mediana, s): " + ", ".join(f"{k} {np.median([n_['sigT3'][k][0] for n_ in ms]):.3f}"
                                                        for k in ("S0", "S2", "ORA")) +
              f"; eventos sorteados por catalogo (media) {np.mean([n_['cat_drawn'] for n_ in ms]):.1f}")

    # --------------------------------------------------------------------------- aciertos
    P("")
    P("## Descriptivo: atomos acertados (todos los atomos, regla LOC): mediana del error epicentral (km) y de |error de "
      f"origen| (s) [IC 95, bootstrap de redes por estrato, {B_GEO} replicas], fraccion de aciertos con |error de origen| "
      "<= 0.5 s, fraccion de eventos tocados")
    for kind in KINDS:
        for s in seeds:
            sets = [(MT3, NN[kind])] + ([(BOMP3, cfg["bomp"])] if kind == "std" else []) + \
                   [((S0T, S2T, ORT) + ((B2T, BOT) if kind == "std" else ()), cfg["ora"])]
            for ms, nnet in sets:
                P(f"  {KNAME[kind]:8s} {s}, primeras {nnet} redes por celda y regimen")
                for m in ms:
                    pt, lo, hi = geometry(data[s][kind], m, nnet, s, kind)
                    P(f"      {m:12s} epi {pt[0]:.2f} [{lo[0]:.2f}, {hi[0]:.2f}] km   |dt| {pt[1]:.3f} [{lo[1]:.3f}, "
                      f"{hi[1]:.3f}] s   a 0.5 s {pt[2]:.3f} [{lo[2]:.3f}, {hi[2]:.3f}]   eventos tocados {pt[3]:.3f}")

    # --------------------------------------------------------------------------- DE y tiempos
    P("")
    P("## DE bootstrap de los contrastes de las reglas frente a la planificada (DE piloto x raiz(30 / redes))")
    for s in seeds:
        for kind in KINDS:
            for rule, _tol in RULES:
                d = diff(E(s, kind, "both", rule), S2T, S0T, 0.25)
                sd0 = PILOT_SD[(kind, rule)]
                P(f"  {s} {KNAME[kind]:8s} {rule:5s} agrupado: DE {d['sd']:.4f}; planificada "
                  f"{sd0 * np.sqrt(PILOT_N / NN[kind]):.4f} (piloto {sd0:.4f})")
        d = diff(E(s, "std", "both", "LOC+T", "bomp"), B2T, B0T, 0.25)
        P(f"  {s} BOMP estandar LOC+T: DE {d['sd']:.4f}; planificada "
          f"{PILOT_SD['bomp'] * np.sqrt(PILOT_N / cfg['bomp']):.4f} (piloto {PILOT_SD['bomp']:.4f})")
    P("")
    P("## Tiempos (s por ventana, media; un proceso por nucleo, 2 procesos en paralelo)")
    for s in seeds:
        for kind in KINDS:
            Ws = list(data[s][kind].values())
            ms = sorted({m for W in Ws for m in W["sec"]})
            P(f"  {s} {KNAME[kind]:8s}: " + "  ".join(f"{m} {np.mean([x for W in Ws for x in W['sec'][m]]):.3f}"
                                                     for m in ms))
        mt = [x for v in meta[s].values() for x in v]
        P(f"  {s} preparacion por parte (red, catalogo, S0, motores) {np.mean([x['t_setup'] for x in mt]):.2f} s; S2 "
          f"{np.mean([x['t_s2'] for x in mt]):.2f} s; total por parte {np.mean([x['t_total'] for x in mt]):.1f} s; suma "
          f"{sum(x['t_total'] for x in mt) / 3600:.2f} h de CPU")
    txt = "\n".join(L_) + "\n"
    with open(a.out, "w") as fh:
        fh.write(txt)
    print(txt)


if __name__ == "__main__":
    main()
