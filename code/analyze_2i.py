"""analyze_2i.py - reglas 2I-a, 2I-b y 2I-c de la Tarea 2I (PREREGISTRO_1j_1k.txt) por semilla, y los puntos secundarios
y descriptivos registrados.  Escrito y congelado (task2i_code_hashes.txt) antes de correr codigo alguno sobre las
semillas 5151 y 6262.

Datos: task2i_parts/ (task2i_run.py).  Brazo prim: 120 redes (j = 0..119) por celda y regimen, 5 ventanas densas
(kmax 7) y 5 estandar (kmax 5) por red.  Brazos secundarios (si se corrieron): thr, s25, s80, redes j = 0..19.
Endpoint (como 2A-2F): recall de eventos a eventos falsos realizados por ventana <= 0.25 (primario) y <= 0.1
(secundario), FROC exacta (analyze_2f.Pool: todos los puntajes del registro como umbrales), punto de operacion elegido
por metodo sobre la propia evaluacion; al agrupar regimenes y celdas la FROC es la del conjunto agrupado (un umbral por
metodo), como en el piloto C8 y en 2F.
Bootstrap jerarquico registrado: 2000 replicas; estratos regimen x celda (6 estratos); en cada estrato las redes con
reemplazo y, dentro de cada red elegida, sus ventanas con reemplazo; las mismas remuestras para todos los metodos
(diferencias pareadas); punto de operacion re-elegido en cada replica.  Generador: default_rng([semilla, TAG del tipo
de ventana (5201 | 5202), indice de regimen, indice de celda, numero de redes]).
Cotas: cota inferior unilateral 97.5 = percentil 2.5 de las diferencias bootstrap; intervalo de 90 = percentiles 5 y
95; IC 95 = percentiles 2.5 y 97.5.  Comparaciones con umbral: >= umbral - 1e-9, <= umbral + 1e-9.
Regla 2I-a (superioridad, densas): HUBER-N+SC - BOMP-PK+SC a tope 0.25, agrupando regimenes y celdas: cota inferior
  unilateral 97.5 > 0 y estimacion puntual >= 0.03, en AMBAS semillas.
Regla 2I-b (no inferioridad, densas): HUBER-N+SC - MT-FW+SC: cota inferior unilateral 97.5 > -0.03, en ambas semillas.
Regla 2I-c (equivalencia, estandar): HUBER-N+SC - BOMP-PK+SC: intervalo de 90 dentro de [-0.05, +0.05] (TOST), en
  ambas semillas.  Se informan ademas p(d* <= -0.05) y p(d* >= +0.05) con (1 + #) / (B + 1).
Secundario y descriptivo (sin regla): tope 0.1; por regimen; por celda; MT-FW+SC - BOMP-PK+SC (replica de 2F),
  MT-FW+SC - MT-FW-NOCAP+SC y MT-FW-NOCAP+SC - BOMP-PK+SC; regla de acierto con coincidencia de tiempo de origen
  (epicentro a <= 5 km Y |origen del atomo - origen verdadero| <= 0.5 s, el evento mas cercano que cumpla ambas; la
  regla R5t05 de verif_2f/hitrule_verif.py; puntajes sin cambio); umbrales de Huber 2.5003, 3.3682 y 5 (redes 0..19);
  sd de sitio 0.25 y 0.8 (redes 0..19; referencia sd 0.5 = las mismas redes y flujos del brazo prim); bootstrap de
  sensibilidad en el que la red elegida se usa en ambos regimenes (la red como unidad de replicacion, como en 2F);
  DE bootstrap frente a la planificada; tiempos.

    python analyze_2i.py [--out analyze_2i_out.txt]
    python analyze_2i.py --dry --out ...      (ensayo con las partes de control de la semilla 999)
"""
import argparse
import json
import os

import numpy as np

from analyze_2f import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
PARTS = os.path.join(HERE, "task2i_parts")
B = 2000
CAPS = (0.1, 0.25)
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
LAB = {"ladder_4": "4 st.", "ladder_6": "6 st.", "unilateral_c40": "one-sided 8"}
REGIMES = ("thin", "strong")
KINDS = ("dense", "std")
TAGK = {"dense": 5201, "std": 5202}
SEEDS = (5151, 6262)
N_WIN = 5
PRIMARY = ("BOMP-PK+SC", "MT-FW+SC", "MT-FW-NOCAP+SC", "HUBER-N+SC")
THR_EXTRA = ("HUBER-N90+SC", "HUBER-N99+SC", "HUBER-C5+SC")
THR_C = {"HUBER-N90+SC": 2.5003, "HUBER-N+SC": 2.7955, "HUBER-N99+SC": 3.3682, "HUBER-C5+SC": 5.0}
H, BO, MT, NC = "HUBER-N+SC", "BOMP-PK+SC", "MT-FW+SC", "MT-FW-NOCAP+SC"
CFG = {"prim": ("prim", tuple(range(120))), "thr": ("thr", tuple(range(20))), "s25": ("s25", tuple(range(20))),
       "s80": ("s80", tuple(range(20)))}
CFG_DRY = {"prim": ("cprim", (60, 67, 74)), "thr": ("cthr", (60, 67)), "s25": ("ds25", (60, 67)),
           "s80": ("ds80", (60, 67))}
EPS = 1e-9
THR_A = 0.03
MARGIN_B = 0.03
MARGIN_C = 0.05
TOL_KM = 5.0
TOL_T = 0.5
# DE bootstrap del piloto (semilla 999, 20 redes por celda y regimen; REPORT_C8.txt, secciones 4 y 10)
PILOT_SD = {"a": (0.0194, 20, 5), "b": (0.0178, 20, 5), "c": (0.0233, 20, 4)}


# ============================================================================================ datos
def load(prefix, seed, nets, parts=PARTS):
    """{kind: {(regime, cell): W}}, W = {"k", "rec": {m: [...]}, "loc", "truth", "taus", "sec"}; ventanas ordenadas
    por red (en el orden de nets) y ventana.  None si falta alguna parte."""
    D = {}
    meta = []
    for reg in REGIMES:
        for cell in CELLS:
            for j in nets:
                path = os.path.join(parts, f"{prefix}_seed{seed}_{reg}_{cell}_net{j}.json")
                if not os.path.exists(path):
                    return None, None
                with open(path) as fh:
                    e = json.load(fh)
                assert (e["seed"], e["regime"], e["cell"], e["net"]) == (seed, reg, cell, j)
                meta.append({x: e[x] for x in ("tags", "kmax", "site_sd", "t_setup", "t_total", "kernel", "methods",
                                               "host_time_utc")})
                for kind, rows in e["windows"].items():
                    assert len(rows) == e["n_win"] == N_WIN
                    W = D.setdefault(kind, {}).setdefault(
                        (reg, cell), {"k": [], "rec": {m: [] for m in e["methods"]}, "loc": {m: [] for m in e["methods"]},
                                      "sec": {m: [] for m in e["methods"]}, "truth": [], "taus": [], "logsd": []})
                    for w in rows:
                        W["k"].append(w["k"])
                        W["truth"].append(w["truth"])
                        W["taus"].append(w["taus"])
                        W["logsd"].append(w["site_logsd"])
                        for m in e["methods"]:
                            W["rec"][m].append(w["records"][m])
                            W["loc"][m].append(w["loc"][m])
                            W["sec"][m].append(w["sec"][m])
    for kind in D:
        for W in D[kind].values():
            W["k"] = np.asarray(W["k"], float)
    return D, meta


def head(Dk, nnet):
    """las primeras nnet redes de cada estrato (las partes se cargan en el orden de las redes)."""
    n = nnet * N_WIN
    out = {}
    for key, W in Dk.items():
        out[key] = {"k": W["k"][:n], "truth": W["truth"][:n], "taus": W["taus"][:n], "logsd": W["logsd"][:n],
                    "rec": {m: v[:n] for m, v in W["rec"].items()}, "loc": {m: v[:n] for m, v in W["loc"].items()},
                    "sec": {m: v[:n] for m, v in W["sec"].items()}}
    return out


def relabel(W, methods, tol_t):
    """registros con la regla de acierto recalculada desde las posiciones y origenes guardados: epicentro a <= 5 km
    (y, si tol_t no es None, |origen del atomo - origen verdadero| <= tol_t) del evento mas cercano que lo cumpla.
    Devuelve los registros y el numero de etiquetas distintas de las guardadas."""
    out = {m: [] for m in methods}
    diff_ = 0
    for wi in range(len(W["k"])):
        tr = np.asarray(W["truth"][wi], float).reshape(-1, 3)
        taus = np.asarray(W["taus"][wi], float)
        for m in methods:
            lab = []
            for (s, h0), loc in zip(W["rec"][m][wi], W["loc"][m][wi]):
                h = -1
                if len(tr):
                    d = np.linalg.norm(tr[:, :2] - np.asarray(loc[:2], float), axis=1)
                    ok = d <= TOL_KM
                    if tol_t is not None:
                        ok = ok & (np.abs(taus - loc[3]) <= tol_t + 1e-12)
                    if ok.any():
                        h = int(np.argmin(np.where(ok, d, np.inf)))
                diff_ += int(h != h0)
                lab.append([s, h])
            out[m].append(lab)
    return out, diff_


# ============================================================================================ bootstrap
def counts(seed, kind, nnet, B_=B, shared=False):
    """{(regime, cell): C (B, nnet * N_WIN)}.  Registrado (shared=False): estratos regimen x celda, redes y ventanas
    remuestreadas de forma independiente en cada estrato.  Sensibilidad (shared=True): en cada celda la red elegida
    se usa en ambos regimenes; ventanas remuestreadas por separado."""
    out = {}
    rows = np.arange(B_)[:, None, None]
    for ci, cell in enumerate(CELLS):
        if shared:
            pick_s = np.random.default_rng([seed, TAGK[kind], 7, ci, nnet]).integers(0, nnet, (B_, nnet))
        for ri, reg in enumerate(REGIMES):
            rng = np.random.default_rng([seed, TAGK[kind], ri, ci, nnet])
            pick = rng.integers(0, nnet, (B_, nnet))
            if shared:
                pick = pick_s
            idx = pick[:, :, None] * N_WIN + rng.integers(0, N_WIN, (B_, nnet, N_WIN))
            C = np.zeros((B_, nnet * N_WIN))
            np.add.at(C, (rows, idx), 1.0)
            out[(reg, cell)] = C
    return out


def evaluate(Dk, CC, keys, methods, recs=None):
    """puntos y replicas por metodo para los estratos keys agrupados.  recs: registros alternativos {key: {m: [...]}}."""
    ks = np.concatenate([Dk[k_]["k"] for k_ in keys])
    C = np.concatenate([CC[k_] for k_ in keys], axis=1)
    one = np.ones((1, len(ks)))
    pt, bs = {}, {}
    for m in methods:
        r = [x for k_ in keys for x in (Dk[k_]["rec"][m] if recs is None else recs[k_][m])]
        p = Pool(r, ks)
        v = p.recall(one)
        pt[m] = {cap: float(v[cap][0]) for cap in CAPS}
        bs[m] = p.recall(C)
    return pt, bs, int(ks.sum()), len(ks)


def diff(ev, a, b, cap):
    pt, bs = ev[0], ev[1]
    d = bs[a][cap] - bs[b][cap]
    p = pt[a][cap] - pt[b][cap]
    lo, hi = np.percentile(d, [2.5, 97.5])
    lo90, hi90 = np.percentile(d, [5.0, 95.0])
    n = len(d) + 1
    return {"point": p, "lo": float(lo), "hi": float(hi), "lo90": float(lo90), "hi90": float(hi90),
            "sd": float(d.std()), "p_le0": float((1 + np.sum(d <= 0)) / n),
            "p_le_mb": float((1 + np.sum(d <= -MARGIN_B)) / n),
            "p_tost_lo": float((1 + np.sum(d <= -MARGIN_C)) / n), "p_tost_hi": float((1 + np.sum(d >= MARGIN_C)) / n)}


def fci(d):
    return f"{d['point']:+.3f} [{d['lo']:+.3f}, {d['hi']:+.3f}]"


def fci90(d):
    return f"{d['point']:+.3f} [{d['lo90']:+.3f}, {d['hi90']:+.3f}]"


def rci(ev, m, cap):
    lo, hi = np.percentile(ev[1][m][cap], [2.5, 97.5])
    return f"{ev[0][m][cap]:.3f} [{lo:.3f}, {hi:.3f}]"


def yn(x):
    return "SI" if x else "NO"


def cumple(x):
    return "SE CUMPLE" if x else "NO SE CUMPLE"


ALL = [(r, c) for r in REGIMES for c in CELLS]


def keys_of(grp):
    return ALL if grp == "both" else [(grp, c) for c in CELLS]


# ============================================================================================ principal
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "analyze_2i_out.txt"))
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--boot", type=int, default=B)
    a = ap.parse_args()
    cfg = CFG_DRY if a.dry else CFG
    seeds = (999,) if a.dry else SEEDS
    nb = a.boot
    L = []
    P = L.append
    P("# Tarea 2I - HUBER-N+SC frente a BOMP-PK+SC y MT-FW+SC tope 6 (analyze_2i.py)" +
      ("   *** ENSAYO, semilla 999, partes de control ***" if a.dry else ""))
    P("")
    P("Recall de eventos a eventos falsos realizados por ventana <= 0.1 / <= 0.25 (FROC exacta, punto de operacion por "
      "metodo sobre la evaluacion; agrupado = un umbral por metodo). Diferencias pareadas: punto [IC 95% jerarquico: "
      f"estratos regimen x celda, redes y ventanas; {nb} replicas; punto de operacion re-elegido en cada replica]. Cota "
      "inferior unilateral 97.5 = extremo inferior del IC 95; IC 90 = percentiles 5 y 95.")

    data, CC, CS, metas = {}, {}, {}, {}
    P("")
    P("Controles de datos:")
    for s in seeds:
        for arm, (prefix, nets) in cfg.items():
            D, meta = load(prefix, s, nets)
            if D is None:
                P(f"  semilla {s}, brazo {arm}: NO CORRIDO (partes ausentes)")
                continue
            data[(s, arm)] = D
            metas[(s, arm)] = meta
            kern = sorted({tuple(tuple(x) for x in m_["kernel"]["blas"]) for m_ in meta})
            env = sorted({m_["kernel"]["OPENBLAS_CORETYPE"] for m_ in meta})
            tags = sorted({(m_["tags"]["dense"], m_["tags"]["std"], m_["site_sd"]) for m_ in meta})
            hts = sorted(m_["host_time_utc"] for m_ in meta)
            for kind in KINDS:
                nwin = sum(len(W["k"]) for W in D[kind].values())
                nev = int(sum(W["k"].sum() for W in D[kind].values()))
                kset = sorted({int(k) for W in D[kind].values() for k in W["k"]})
                lsd = float(np.mean([x for W in D[kind].values() for x in W["logsd"]]))
                P(f"  semilla {s}, brazo {arm}, {kind}: {len(nets)} redes x 3 celdas x 2 regimenes, {nwin} ventanas, "
                  f"{nev} eventos, k en {kset}; sd media del log de los factores de sitio por ventana {lsd:.3f}")
            P(f"      (TAG densas, TAG estandar, sd de sitio) {tags}; OPENBLAS_CORETYPE {env}; OpenBLAS cargadas {kern}; "
              f"partes escritas entre {hts[0]} y {hts[-1]} UTC")
        for kind in KINDS:
            for nnet in sorted({len(cfg[arm][1]) for arm in cfg}):
                CC[(s, kind, nnet)] = counts(s, kind, nnet, nb)
            CS[(s, kind)] = counts(s, kind, len(cfg["prim"][1]), nb, shared=True)
    nprim = len(cfg["prim"][1])
    if any((s, "prim") not in data for s in seeds):
        P("")
        P("BRAZO PRIMARIO INCOMPLETO: no se evaluan las reglas.")
        open(a.out, "w").write("\n".join(L) + "\n")
        print("\n".join(L))
        return

    ev = {}
    for s in seeds:
        for kind in KINDS:
            Dk = data[(s, "prim")][kind]
            for grp in REGIMES + ("both",):
                ev[(s, kind, grp)] = evaluate(Dk, CC[(s, kind, nprim)], keys_of(grp), PRIMARY)

    # --------------------------------------------------------------------------- recall por metodo
    P("")
    P("## Recall por metodo (brazo primario), tope 0.1 / tope 0.25, con IC 95")
    for kind in KINDS:
        for grp in ("both",) + REGIMES:
            for s in seeds:
                e = ev[(s, kind, grp)]
                P(f"  {kind:5s} {grp:6s} {s}: {e[3]} ventanas, {e[2]} eventos")
                for m in PRIMARY:
                    P(f"      {m:15s} {rci(e, m, 0.1)}   {rci(e, m, 0.25)}")

    # --------------------------------------------------------------------------- reglas
    P("")
    P("## Regla 2I-a (superioridad, ventanas densas): HUBER-N+SC - BOMP-PK+SC a tope 0.25, agrupando regimenes y celdas")
    ok_a = []
    for s in seeds:
        d = diff(ev[(s, "dense", "both")], H, BO, 0.25)
        c1, c2 = d["lo"] > 0, d["point"] >= THR_A - EPS
        ok_a.append(c1 and c2)
        P(f"  {s}: {fci(d)}; cota inferior unilateral 97.5 = {d['lo']:+.4f} > 0: {yn(c1)}; punto {d['point']:+.4f} >= "
          f"{THR_A}: {yn(c2)}; p(d* <= 0) = {d['p_le0']:.4f}; DE bootstrap {d['sd']:.4f} -> {cumple(c1 and c2)}")
    va = all(ok_a) and len(seeds) == 2
    P(f"  -> 2I-a {cumple(va)} (por semilla: {[yn(x) for x in ok_a]})")
    P("")
    P("## Regla 2I-b (no inferioridad, ventanas densas): HUBER-N+SC - MT-FW+SC a tope 0.25, margen -0.03")
    ok_b = []
    for s in seeds:
        d = diff(ev[(s, "dense", "both")], H, MT, 0.25)
        c1 = d["lo"] > -MARGIN_B
        ok_b.append(c1)
        P(f"  {s}: {fci(d)}; cota inferior unilateral 97.5 = {d['lo']:+.4f} > -{MARGIN_B}: {yn(c1)}; p(d* <= -0.03) = "
          f"{d['p_le_mb']:.4f}; DE bootstrap {d['sd']:.4f} -> {cumple(c1)}")
    vb = all(ok_b) and len(seeds) == 2
    P(f"  -> 2I-b {cumple(vb)} (por semilla: {[yn(x) for x in ok_b]})")
    P("")
    P("## Regla 2I-c (equivalencia, ventanas estandar): HUBER-N+SC - BOMP-PK+SC a tope 0.25, TOST, margen 0.05")
    ok_c = []
    for s in seeds:
        e = ev[(s, "std", "both")]
        d = diff(e, H, BO, 0.25)
        c1 = d["lo90"] >= -MARGIN_C - EPS and d["hi90"] <= MARGIN_C + EPS
        ok_c.append(c1)
        P(f"  {s}: {e[3]} ventanas, {e[2]} eventos; diferencia {d['point']:+.4f}; IC 90 [{d['lo90']:+.4f}, "
          f"{d['hi90']:+.4f}] dentro de [-0.05, +0.05]: {yn(c1)}; IC 95 [{d['lo']:+.3f}, {d['hi']:+.3f}]; p(d* <= -0.05) = "
          f"{d['p_tost_lo']:.4f}, p(d* >= +0.05) = {d['p_tost_hi']:.4f}; DE bootstrap {d['sd']:.4f} -> {cumple(c1)}")
    vc = all(ok_c) and len(seeds) == 2
    P(f"  -> 2I-c {cumple(vc)} (por semilla: {[yn(x) for x in ok_c]})")

    P("")
    P("## Veredictos")
    P(f"  2I-a: {cumple(va)}    2I-b: {cumple(vb)}    2I-c: {cumple(vc)}")
    if va and vb and vc:
        interp = ("a, b y c se cumplen -> la afirmacion 'un ajuste de Huber escalado por el ruido reemplaza al tope "
                  "calibrado' queda confirmada en este diseno sintetico.")
    elif va and not vb:
        interp = ("a se cumple y b no -> HUBER-N ayuda pero el tope es mejor (lectura registrada)." +
                  ("" if vc else " Ademas c no se cumple."))
    elif va:
        interp = ("a y b se cumplen y c no -> la afirmacion de reemplazo NO queda confirmada (exige a, b y c): no se "
                  "demostro la equivalencia en eventos aislados.")
    else:
        interp = "a falla -> segun el registro, el piloto fue un artefacto de seleccion."
    P(f"  Interpretacion registrada: {interp}")

    # --------------------------------------------------------------------------- secundario
    PAIRS_D = ((H, BO), (H, MT), (MT, BO), (MT, NC), (NC, BO), (H, NC))
    P("")
    P("## Secundario / descriptivo: contrastes agrupados y por regimen, tope 0.1 y tope 0.25 (IC 95)")
    P("   (MT-FW+SC - BOMP-PK+SC en ventanas densas por regimen es la replica de 2F-a; MT-FW+SC - MT-FW-NOCAP+SC y "
      "MT-FW-NOCAP+SC - BOMP-PK+SC, la de 2F-b; MT-FW+SC - BOMP-PK+SC en estandar, la de 2F-c)")
    for kind in KINDS:
        for grp in ("both",) + REGIMES:
            for s in seeds:
                e = ev[(s, kind, grp)]
                for x, y in PAIRS_D:
                    d1, d2 = diff(e, x, y, 0.1), diff(e, x, y, 0.25)
                    P(f"  {kind:5s} {grp:6s} {s} {x} - {y}: tope 0.1 {fci(d1)}; tope 0.25 {fci(d2)}; IC 90 a 0.25 "
                      f"[{d2['lo90']:+.3f}, {d2['hi90']:+.3f}]")
    P("  Replica de 2F en estos datos (criterios de 2F aplicados de forma descriptiva; 2F uso 10 ventanas por red):")
    for s in seeds:
        for reg in REGIMES:
            d = diff(ev[(s, "dense", reg)], MT, BO, 0.25)
            P(f"    {s} {reg}: MT-FW+SC - BOMP-PK+SC {fci(d)}: cota inferior > 0 {yn(d['lo'] > 0)}, punto >= 0.03 "
              f"{yn(d['point'] >= 0.03 - EPS)}")
        e = ev[(s, "dense", "both")]
        d1, d2 = diff(e, MT, NC, 0.25), diff(e, NC, BO, 0.25)
        P(f"    {s} ambos: MT-FW+SC - MT-FW-NOCAP+SC {fci(d1)} (cota inferior > 0 {yn(d1['lo'] > 0)}); "
          f"MT-FW-NOCAP+SC - BOMP-PK+SC {fci(d2)} (|punto| <= 0.02 {yn(abs(d2['point']) <= 0.02 + EPS)})")
        d3 = diff(ev[(s, "std", "both")], MT, BO, 0.25)
        P(f"    {s} estandar: MT-FW+SC - BOMP-PK+SC {d3['point']:+.3f}, IC 90 [{d3['lo90']:+.3f}, {d3['hi90']:+.3f}] "
          f"(dentro de +-0.05: {yn(d3['lo90'] >= -MARGIN_C - EPS and d3['hi90'] <= MARGIN_C + EPS)})")

    P("")
    P("## Descriptivo: por celda y regimen (tope 0.25, IC 95)")
    for kind in KINDS:
        for s in seeds:
            Dk = data[(s, "prim")][kind]
            for reg in REGIMES:
                for c in CELLS:
                    e = evaluate(Dk, CC[(s, kind, nprim)], [(reg, c)], PRIMARY)
                    P(f"  {kind:5s} {s} {reg:6s} {LAB[c]:12s} ({e[3]} v., {e[2]} ev.): " +
                      "  ".join(f"{m} {e[0][m][0.25]:.3f}" for m in PRIMARY) +
                      f" | H-N - BOMP {fci(diff(e, H, BO, 0.25))}; H-N - MT-FW {fci(diff(e, H, MT, 0.25))}; MT-FW - BOMP "
                      f"{fci(diff(e, MT, BO, 0.25))}")

    # --------------------------------------------------------------------------- regla de acierto con tiempo de origen
    P("")
    P(f"## Secundario: regla de acierto con coincidencia de tiempo de origen (<= {TOL_KM:g} km Y |error de origen| <= "
      f"{TOL_T} s); puntajes sin cambio; mismo bootstrap")
    for kind in KINDS:
        for s in seeds:
            Dk = data[(s, "prim")][kind]
            r0, mism = {}, 0
            rt = {}
            for k_ in ALL:
                r0[k_], dd = relabel(Dk[k_], PRIMARY, None)
                mism += dd
                rt[k_], _dd = relabel(Dk[k_], PRIMARY, TOL_T)
            P(f"  {kind:5s} {s}: regla registrada recalculada desde las posiciones guardadas: {mism} etiquetas distintas "
              "de las guardadas")
            for grp in ("both",) + REGIMES:
                e = evaluate(Dk, CC[(s, kind, nprim)], keys_of(grp), PRIMARY, recs=rt)
                P(f"      {grp:6s} recall 0.1/0.25: " + "  ".join(f"{m} {e[0][m][0.1]:.3f}/{e[0][m][0.25]:.3f}"
                                                                 for m in PRIMARY))
                for x, y in ((H, BO), (H, MT), (MT, BO)):
                    d1, d2 = diff(e, x, y, 0.1), diff(e, x, y, 0.25)
                    P(f"          {x} - {y}: tope 0.1 {fci(d1)}; tope 0.25 {fci(d2)}; IC 90 a 0.25 "
                      f"[{d2['lo90']:+.3f}, {d2['hi90']:+.3f}]")

    # --------------------------------------------------------------------------- umbrales
    P("")
    P("## Secundario: umbral de Huber (redes del brazo thr; BOMP-PK+SC, MT-FW+SC y HUBER-N+SC de las mismas ventanas "
      "del brazo prim)")
    nthr = len(cfg["thr"][1])
    for s in seeds:
        if (s, "thr") not in data:
            P(f"  {s}: brazo thr no corrido")
            continue
        for kind in KINDS:
            Dp = head(data[(s, "prim")][kind], nthr)
            Dt = data[(s, "thr")][kind]
            same = tot = 0
            for k_ in ALL:
                assert Dp[k_]["truth"] == Dt[k_]["truth"] and Dp[k_]["taus"] == Dt[k_]["taus"]
                for r1, r2, l1, l2 in zip(Dp[k_]["rec"][H], Dt[k_]["rec"][H], Dp[k_]["loc"][H], Dt[k_]["loc"][H]):
                    tot += 1
                    same += int(r1 == r2 and l1 == l2)
                for m in THR_EXTRA:
                    Dp[k_]["rec"][m] = Dt[k_]["rec"][m]
            ms = (BO, MT) + tuple(sorted(THR_C, key=lambda m: THR_C[m]))
            e = evaluate(Dp, CC[(s, kind, nthr)], ALL, ms)
            P(f"  {kind:5s} {s}: {e[3]} ventanas, {e[2]} eventos; HUBER-N+SC regenerado identico al del brazo prim en "
              f"{same}/{tot} ventanas")
            for m in ms:
                c = f"c = {THR_C[m]:.4f}" if m in THR_C else ""
                row = f"      {m:13s} {c:11s} recall {e[0][m][0.1]:.3f}/{e[0][m][0.25]:.3f}"
                if m in THR_C:
                    row += (f" | a 0.25: - BOMP-PK {fci(diff(e, m, BO, 0.25))}; - MT-FW {fci(diff(e, m, MT, 0.25))}" +
                            (f"; - HUBER-N {fci(diff(e, m, H, 0.25))}" if m != H else "") +
                            f" | a 0.1: - BOMP-PK {fci(diff(e, m, BO, 0.1))}")
                else:
                    row += f" | a 0.25: MT-FW - BOMP-PK {fci(diff(e, MT, BO, 0.25))}" if m == MT else ""
                P(row)

    # --------------------------------------------------------------------------- sd de sitio
    P("")
    P("## Secundario: sd de los factores de sitio (redes de los brazos s25 / s80; sd 0.5 = las mismas redes y flujos del "
      "brazo prim)")
    for s in seeds:
        for kind in KINDS:
            for arm, sd in (("s25", 0.25), ("prim", 0.5), ("s80", 0.8)):
                if (s, arm) not in data:
                    P(f"  {kind:5s} {s} sd {sd}: brazo {arm} no corrido")
                    continue
                nn = len(cfg["s25"][1])
                Dk = head(data[(s, arm)][kind], nn)
                e = evaluate(Dk, CC[(s, kind, nn)], ALL, PRIMARY)
                P(f"  {kind:5s} {s} sd {sd}: {e[3]} ventanas, {e[2]} eventos; recall 0.1/0.25 " +
                  "  ".join(f"{m} {e[0][m][0.1]:.3f}/{e[0][m][0.25]:.3f}" for m in PRIMARY))
                for x, y in ((H, BO), (H, MT), (MT, BO), (MT, NC)):
                    d1, d2 = diff(e, x, y, 0.1), diff(e, x, y, 0.25)
                    P(f"          {x} - {y}: tope 0.1 {fci(d1)}; tope 0.25 {fci(d2)}")

    # --------------------------------------------------------------------------- sensibilidad del bootstrap
    P("")
    P("## Sensibilidad (descriptivo): bootstrap con la red como unidad de replicacion compartida por los dos regimenes "
      "(como en 2F); tope 0.25")
    for s in seeds:
        for kind, pairs in (("dense", ((H, BO), (H, MT), (MT, BO))), ("std", ((H, BO), (MT, BO)))):
            e = evaluate(data[(s, "prim")][kind], CS[(s, kind)], ALL, PRIMARY)
            for x, y in pairs:
                d = diff(e, x, y, 0.25)
                P(f"  {kind:5s} {s} {x} - {y}: {fci(d)}; IC 90 [{d['lo90']:+.3f}, {d['hi90']:+.3f}]; DE {d['sd']:.4f}")

    # --------------------------------------------------------------------------- DE y tiempos
    P("")
    P("## DE bootstrap de los contrastes de las reglas frente a la planificada (DE piloto x raiz(redes y ventanas))")
    for s in seeds:
        for rule, kind, x, y in (("a", "dense", H, BO), ("b", "dense", H, MT), ("c", "std", H, BO)):
            sd0, n0, w0 = PILOT_SD[rule]
            plan = sd0 * np.sqrt(n0 / nprim) * (np.sqrt(w0 / N_WIN) if kind == "std" else 1.0)
            d = diff(ev[(s, kind, "both")], x, y, 0.25)
            P(f"  {s} 2I-{rule}: DE {d['sd']:.4f}; planificada {plan:.4f} (piloto {sd0:.4f})")
    P("")
    P("## Tiempos (s por ventana, media; un proceso por nucleo, 2 procesos en paralelo)")
    for s in seeds:
        for arm in cfg:
            if (s, arm) not in data:
                continue
            for kind in KINDS:
                Ws = list(data[(s, arm)][kind].values())
                ms = list(Ws[0]["sec"])
                P(f"  {s} {arm:4s} {kind:5s}: " + "  ".join(f"{m} {np.mean([x for W in Ws for x in W['sec'][m]]):.3f}"
                                                         for m in ms))
            mt = metas[(s, arm)]
            P(f"  {s} {arm:4s} preparacion por parte (red, catalogo, SC, operadores) {np.mean([x['t_setup'] for x in mt]):.1f}"
              f" s; total por parte {np.mean([x['t_total'] for x in mt]):.1f} s; suma {sum(x['t_total'] for x in mt) / 3600:.2f}"
              " h de CPU")
    txt = "\n".join(L) + "\n"
    with open(a.out, "w") as fh:
        fh.write(txt)
    print(txt)


if __name__ == "__main__":
    main()
