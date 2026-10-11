"""analyze_2a.py - reglas 2A-a a 2A-e de la Tarea 2A (PREREGISTRO_1j_1k.txt), por regimen de clutter, y tablas.

Las estimaciones puntuales deciden las reglas (como analyze_1y.py): recall a eventos falsos realizados por ventana
<= 0.1 y <= 0.25 sobre la curva FROC exacta (froc_exact.exact_HF) de cada metodo, con el punto de operacion elegido
sobre la propia evaluacion; GaMMA y PyOcto eligen ademas la configuracion por oraculo (task1y_clutter.run_cell).
IC: bootstrap pareado por ventanas (2000 replicas, las mismas remuestras para todos los metodos de una celda),
re-eligiendo en cada replica el punto de operacion, la configuracion y, para "mejor comparador", el maximo.
2A-e: umbrales y configuraciones congelados en la semilla 999 (task2a_seed999_{regimen}.json) por celda, regimen y
presupuesto, aplicados sin cambios a la evaluacion.
Umbral de las reglas: una diferencia cuenta si es >= 0.05 - 1e-9 (tolerancia de redondeo de coma flotante; las
diferencias de recall son multiplos de 1/eventos); se listan los casos a menos de 0.005 del umbral.

    python analyze_2a.py [--seeds 2024,7331] [--cal 999] [--out analyze_2a_out.txt]
"""
import argparse
import json

import numpy as np

from froc_exact import exact_HF

B = 2000
CAPS = ("0.1", "0.25")
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
LAB = {"ladder_4": "4 st.", "ladder_6": "6 st.", "unilateral_c40": "one-sided 8"}
REGIMES = ("thin", "strong")
SOTA = ("GaMMA-STALTA", "GaMMA-STALTA+SC", "PyOcto-STALTA", "PyOcto-STALTA+SC", "QM", "QM-3C", "QM-3C+SC")
DYS = ("GaMMA-DY", "GaMMA-DY+SC")
ALL = ("MT-FW+SC", "MT-FW", "MT-FW-NODE+SC", "WAVE-PT@20") + SOTA + DYS
EPS = 1e-9
THR = 0.05


def mats(m):
    """(H, F) ventana x umbral (registros, curva exacta) o ventana x configuracion (asociadores)."""
    if "records" in m:
        _t, H, F = exact_HF([[tuple(x) for x in r] for r in m["records"]])
        return np.asarray(H, float), np.asarray(F, float)
    pw = m["per_window"]
    H = np.array([[c[0] for c in w] for w in pw], float)
    F = np.array([[c[1] for c in w] for w in pw], float)
    return H, F


def val_counts(HF, C, kk, n, cap):
    """recall a tope cap para cada fila de C (conteos de remuestreo, B x n): max sobre umbrales o configuraciones
    con tasa de falsos realizada <= cap (0 si ninguno)."""
    H, F = HF
    Hs, Fs = C @ H, C @ F
    ok = Fs / n <= cap + 1e-12
    r = np.where(ok, Hs, -1.0).max(axis=1)
    return np.where(r < 0, 0.0, r) / np.maximum(1.0, kk)


def point(HF, k, cap):
    n = len(k)
    return float(val_counts(HF, np.ones((1, n)), np.array([k.sum()]), n, cap)[0])


def ge(d):
    return d >= THR - EPS


def fmt_ci(p, lo, hi):
    return f"{p:+.2f} [{lo:+.2f}, {hi:+.2f}]"


# ------------------------------------------------------------------------------------------ 2A-e
def freeze_threshold(recs, k, budget):
    """umbral de puntaje de la curva exacta de calibracion con recall maximo y tasa de falsos <= budget
    (empates: menor tasa de falsos, es decir el umbral mas alto). None si ningun umbral cumple."""
    thr, H, F = exact_HF([[tuple(x) for x in r] for r in recs])
    rec = H.sum(0) / max(1, k.sum())
    fr = F.sum(0) / H.shape[0]
    ok = np.where(fr <= budget + 1e-12)[0]
    if len(ok) == 0:
        return None, 0.0, 0.0
    j = min(ok, key=lambda j_: (-rec[j_], fr[j_], -thr[j_]))
    return float(thr[j]), float(rec[j]), float(fr[j])


def apply_threshold(recs, k, t):
    """eventos con puntaje >= t (prefijo de cada lista): recall y falsos por ventana realizados."""
    if t is None:
        return 0.0, 0.0, np.zeros(len(recs)), np.zeros(len(recs))
    h, f = np.zeros(len(recs)), np.zeros(len(recs))
    for w, r in enumerate(recs):
        hits = set()
        for s, hh in r:
            if s >= t:
                if hh >= 0:
                    hits.add(hh)
                else:
                    f[w] += 1
        h[w] = len(hits)
    return float(h.sum() / max(1, k.sum())), float(f.mean()), h, f


def freeze_config(per_window, k, budget):
    H = np.array([[c[0] for c in w] for w in per_window], float)
    F = np.array([[c[1] for c in w] for w in per_window], float)
    rec = H.sum(0) / max(1, k.sum())
    fr = F.sum(0) / H.shape[0]
    ok = [c for c in range(H.shape[1]) if fr[c] <= budget + 1e-12]
    if not ok:
        return None, 0.0, 0.0
    c = min(ok, key=lambda c_: (-rec[c_], fr[c_], c_))
    return int(c), float(rec[c]), float(fr[c])


def apply_config(per_window, k, c):
    if c is None:
        return 0.0, 0.0, np.zeros(len(per_window)), np.zeros(len(per_window))
    h = np.array([w[c][0] for w in per_window], float)
    f = np.array([w[c][1] for w in per_window], float)
    return float(h.sum() / max(1, k.sum())), float(f.mean()), h, f


def lut_clip_fraction(STA, corr, vp=1.73 * 3.0, vs=3.0, fs=50):
    """fraccion de entradas (nodo, estacion, fase) del LUT de QuakeMigrate (x, y en [-21, 22], z en [0, 20] km, 1 km)
    cuyo tiempo corregido rint((t + corr) fs) es negativo: la libreria C de QuakeMigrate las recorta a 0
    (migratelib.c: MAX(0, traveltime)). Devuelve (todas las profundidades, profundidades 5-13 km)."""
    ax = np.arange(-21.0, 22.0 + 1e-9, 1.0)
    zz = np.arange(0.0, 20.0 + 1e-9, 1.0)
    X_, Y_, Z_ = np.meshgrid(ax, ax, zz, indexing="ij")
    n_all = n_dom = c_all = c_dom = 0
    dom = (Z_ >= 5) & (Z_ <= 13)
    for i, sta in enumerate(np.asarray(STA, float)):
        r = np.sqrt((X_ - sta[0]) ** 2 + (Y_ - sta[1]) ** 2 + (Z_ - sta[2]) ** 2)
        for b, v in enumerate((vp, vs)):
            neg = np.rint((r / v + corr[i][b]) * fs) < 0
            c_all += int(neg.sum())
            c_dom += int(neg[dom].sum())
            n_all += neg.size
            n_dom += int(dom.sum())
    return c_all / n_all, c_dom / n_dom


def freeze_all(cal_seed, tag=""):
    """puntos de operacion de 2A-e congelados en la semilla de calibracion, por regimen, celda y presupuesto."""
    out = {}
    for reg in REGIMES:
        CAL = {c["scenario"]: c for c in json.load(open(f"task2a_seed{cal_seed}_{reg}{tag}.json"))}
        for cell in CELLS:
            cal = CAL[cell]
            kc = np.asarray(cal["k_eval"], float)
            for bud in (0.1, 0.25):
                tM, rMc, fMc = freeze_threshold(cal["methods"]["MT-FW+SC"]["records"], kc, bud)
                tQ, rQc, fQc = freeze_threshold(cal["methods"]["QM-3C+SC"]["records"], kc, bud)
                cG, rGc, fGc = freeze_config(cal["methods"]["GaMMA-STALTA+SC"]["per_window"], kc, bud)
                out[f"{reg}|{cell}|{bud}"] = {"MT-FW+SC": [tM, rMc, fMc], "QM-3C+SC": [tQ, rQc, fQc],
                                             "GaMMA-STALTA+SC": [cG, rGc, fGc], "n_cal_windows": len(kc),
                                             "n_cal_events": int(kc.sum())}
    return out


# ------------------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="2024,7331")
    ap.add_argument("--cal", type=int, default=999)
    ap.add_argument("--out", default="analyze_2a_out.txt")
    ap.add_argument("--tag", default="")
    ap.add_argument("--freeze", action="store_true", help="solo congela los puntos de 2A-e (task2a_frozen_2Ae.json)")
    a = ap.parse_args()
    ffile = f"task2a_frozen_2Ae{a.tag}.json"
    if a.freeze:
        import time
        fz = freeze_all(a.cal, a.tag)
        fz["_written"] = time.strftime("%Y-%m-%d %H:%M:%S %Z")
        with open(ffile, "w") as fh:
            json.dump(fz, fh, indent=1)
        for key, v in fz.items():
            print(key, v)
        return
    FZ = json.load(open(ffile))
    chk = freeze_all(a.cal, a.tag)
    assert all(chk[key] == FZ[key] for key in chk), "los puntos congelados no coinciden con la calibracion"
    seeds = [int(s) for s in a.seeds.split(",")]
    L = []
    P = L.append
    verdicts = {}
    P("# Tarea 2A - analisis confirmatorio (analyze_2a.py)")
    P("")
    P(f"Recall de eventos a eventos falsos realizados por ventana <= 0.1 / <= 0.25 (curva FROC exacta; asociadores con "
      f"oraculo de configuracion). IC 95%: bootstrap pareado por ventanas, {B} replicas, re-eligiendo punto de "
      f"operacion, configuracion y maximo en cada replica. Semillas {seeds}; calibracion 2A-e: semilla {a.cal}.")
    near = []
    for reg in REGIMES:
        R = {s: {c["scenario"]: c for c in json.load(open(f"task2a_seed{s}_{reg}{a.tag}.json"))} for s in seeds}
        P("")
        P(f"## Regimen {reg} (clutter {R[seeds[0]][CELLS[0]]['clut_amp']})")
        P("")
        P(f"### Tabla 1 ({reg}). Recall a <= 0.1 / <= 0.25")
        P("")
        P("| Celda (semilla) | eventos | " + " | ".join(ALL) + " |")
        P("|:--|--:|" + ":-:|" * len(ALL))
        rows_diff, rows_best = [], []
        counts = {r_: {s: 0 for s in seeds} for r_ in ("2A-a", "2A-b", "2A-c", "2A-d")}
        detail = {r_: {s: [] for s in seeds} for r_ in ("2A-a", "2A-b", "2A-c", "2A-d")}
        e_rows = []
        e_count = {s: 0 for s in seeds}
        e_count_strict = {s: [0, True] for s in seeds}
        for s in seeds:
            for ci, cell in enumerate(CELLS):
                c = R[s][cell]
                k = np.asarray(c["k_eval"], float)
                n = len(k)
                M = {m: mats(c["methods"][m]) for m in ALL}
                pt = {m: {cap: point(M[m], k, float(cap)) for cap in CAPS} for m in ALL}
                P(f"| {LAB[cell]} ({s}) | {int(k.sum())} | " + " | ".join(f"{pt[m]['0.1']:.2f}/{pt[m]['0.25']:.2f}" for m in ALL) + " |")
                # reglas (estimaciones puntuales)
                best = {cap: max(SOTA, key=lambda m: pt[m][cap]) for cap in CAPS}
                bdy = {cap: max(DYS, key=lambda m: pt[m][cap]) for cap in CAPS}
                da = {cap: pt["MT-FW+SC"][cap] - pt[best[cap]][cap] for cap in CAPS}
                db = {cap: pt["MT-FW+SC"][cap] - pt["MT-FW-NODE+SC"][cap] for cap in CAPS}
                dc = pt["MT-FW+SC"]["0.25"] - pt[bdy["0.25"]]["0.25"]
                dd = pt["MT-FW+SC"]["0.25"] - pt["MT-FW"]["0.25"]
                okA, okB, okC, okD = all(ge(da[x]) for x in CAPS), all(ge(db[x]) for x in CAPS), ge(dc), ge(dd)
                for r_, ok_, v in (("2A-a", okA, da), ("2A-b", okB, db), ("2A-c", okC, {"0.25": dc}), ("2A-d", okD, {"0.25": dd})):
                    counts[r_][s] += ok_
                    detail[r_][s].append((cell, ok_, v))
                    for cap, x in v.items():
                        if abs(x - THR) < 0.005:
                            near.append(f"{reg} {s} {cell} {r_} tope {cap}: diferencia {x:+.4f}")
                # bootstrap pareado
                rng = np.random.default_rng([s, ci, REGIMES.index(reg)])
                idx = rng.integers(0, n, size=(B, n))
                C = np.zeros((B, n))
                np.add.at(C, (np.repeat(np.arange(B), n), idx.ravel()), 1.0)
                kk = C @ k
                bv = {m: {cap: val_counts(M[m], C, kk, n, float(cap)) for cap in CAPS} for m in ALL}
                cells_d = []
                for m in ALL[1:]:
                    cc = []
                    for cap in CAPS:
                        d = bv["MT-FW+SC"][cap] - bv[m][cap]
                        cc.append(fmt_ci(pt["MT-FW+SC"][cap] - pt[m][cap], *np.percentile(d, [2.5, 97.5])))
                    cells_d.append(" / ".join(cc))
                rows_diff.append(f"| {LAB[cell]} ({s}) | " + " | ".join(cells_d) + " |")
                cb = []
                for grp, cand in (("SOTA", SOTA), ("DY", DYS)):
                    for cap in CAPS:
                        mx = np.max(np.stack([bv[m][cap] for m in cand]), axis=0)
                        d = bv["MT-FW+SC"][cap] - mx
                        p0 = pt["MT-FW+SC"][cap] - max(pt[m][cap] for m in cand)
                        cb.append(fmt_ci(p0, *np.percentile(d, [2.5, 97.5])))
                rows_best.append(f"| {LAB[cell]} ({s}) | {best['0.1']} / {best['0.25']} | {cb[0]} | {cb[1]} | "
                                 f"{bdy['0.25']} | {cb[3]} |")
                # 2A-e
                okE = True
                okE_rate = True
                okE_rec = True
                for bud in (0.1, 0.25):
                    fz = FZ[f"{reg}|{cell}|{bud}"]
                    tM, rMc, fMc = fz["MT-FW+SC"]
                    tQ, rQc, fQc = fz["QM-3C+SC"]
                    cG, rGc, fGc = fz["GaMMA-STALTA+SC"]
                    rM, fM, hM, fwM = apply_threshold(c["methods"]["MT-FW+SC"]["records"], k, tM)
                    rQ, fQ, hQ, fwQ = apply_threshold(c["methods"]["QM-3C+SC"]["records"], k, tQ)
                    rG, fG, hG, fwG = apply_config(c["methods"]["GaMMA-STALTA+SC"]["per_window"], k, cG)
                    bestc = "GaMMA-STALTA+SC" if rG >= rQ else "QM-3C+SC"
                    rB = max(rG, rQ)
                    hB = hG if rG >= rQ else hQ
                    d = rM - rB
                    # IC del bootstrap pareado (umbrales congelados, el mejor comparador re-elegido por replica)
                    dm = (C @ hM) / np.maximum(1, kk) - np.maximum(C @ hG, C @ hQ) / np.maximum(1, kk)
                    lo, hi = np.percentile(dm, [2.5, 97.5])
                    rate_ok = fM <= 2 * bud + 1e-12
                    rec_ok = ge(d)
                    okE &= rate_ok and rec_ok
                    okE_rate &= rate_ok
                    okE_rec &= rec_ok
                    if abs(d - THR) < 0.005:
                        near.append(f"{reg} {s} {cell} 2A-e presupuesto {bud}: diferencia {d:+.4f}")
                    e_rows.append(f"| {LAB[cell]} ({s}) | {bud} | {'-' if tM is None else f'{tM:.2f}'} ({rMc:.2f}/{fMc:.2f}) | "
                                  f"{rM:.2f} / {fM:.2f} | {'-' if cG is None else cG} ({rGc:.2f}/{fGc:.2f}) | {rG:.2f} / {fG:.2f} | "
                                  f"{'-' if tQ is None else f'{tQ:.3f}'} ({rQc:.2f}/{fQc:.2f}) | {rQ:.2f} / {fQ:.2f} | "
                                  f"{bestc} | {fmt_ci(d, lo, hi)} | {'si' if rate_ok else 'NO'} |")
                e_count[s] += okE
                e_count_strict[s][0] += okE_rec
                e_count_strict[s][1] &= okE_rate
        P("")
        P(f"### Tabla 2 ({reg}). MT-FW+SC menos cada metodo, a <= 0.1 / <= 0.25: punto [IC 95%]")
        P("")
        P("| Celda (semilla) | " + " | ".join(f"- {m}" for m in ALL[1:]) + " |")
        P("|:--|" + ":-:|" * (len(ALL) - 1))
        L.extend(rows_diff)
        P("")
        P(f"### Tabla 3 ({reg}). MT-FW+SC menos el mejor comparador por celda y tope (maximo re-elegido en cada replica)")
        P("")
        P("| Celda (semilla) | mejor SOTA (0.1 / 0.25) | - max SOTA a 0.1 | - max SOTA a 0.25 | mejor DY (0.25) | - max(GaMMA-DY, GaMMA-DY+SC) a 0.25 |")
        P("|:--|:--|:-:|:-:|:--|:-:|")
        L.extend(rows_best)
        P("")
        P(f"### Tabla 4 ({reg}). 2A-e: puntos de operacion congelados en la semilla {a.cal} (5 redes x 20 ventanas por celda), "
          "aplicados a la evaluacion")
        P("")
        P("Entre parentesis: recall / falsos por ventana en calibracion. Evaluacion: recall / falsos por ventana realizados. "
          "GaMMA: indice de GAMMA_GRID congelado. Diferencia: MT-FW+SC menos el mejor comparador congelado [IC 95%].")
        P("")
        P("| Celda (semilla) | presupuesto | MT-FW+SC umbral (cal) | MT-FW+SC eval | GaMMA-STALTA+SC cfg (cal) | GaMMA eval | "
          "QM-3C+SC umbral (cal) | QM-3C+SC eval | mejor congelado | diferencia | tasa MT-FW+SC <= 2x |")
        P("|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:--|:-:|:-:|")
        L.extend(e_rows)
        P("")
        P(f"### Reglas ({reg})")
        P("")
        v = {}
        for r_ in ("2A-a", "2A-b", "2A-c", "2A-d"):
            ok = all(counts[r_][s] >= 2 for s in seeds)
            v[r_] = ok
            P(f"- {r_}: " + "; ".join(f"semilla {s}: {counts[r_][s]}/3 celdas ("
                                       + ", ".join(f"{LAB[c_]} {'si' if o else 'no'} "
                                                   + "/".join(f"{x:+.2f}" for x in vv.values())
                                                   for c_, o, vv in detail[r_][s]) + ")" for s in seeds)
              + f" -> {'SE CUMPLE' if ok else 'NO SE CUMPLE'}")
        okE = all(e_count[s] >= 2 for s in seeds)
        okE2 = all(e_count_strict[s][0] >= 2 and e_count_strict[s][1] for s in seeds)
        v["2A-e"] = okE
        P(f"- 2A-e (lectura por celda: tasa <= 2x presupuesto y diferencia >= 0.05 a ambos presupuestos en la celda): "
          + "; ".join(f"semilla {s}: {e_count[s]}/3" for s in seeds) + f" -> {'SE CUMPLE' if okE else 'NO SE CUMPLE'}")
        P(f"- 2A-e (lectura estricta: tasa <= 2x en todas las celdas y presupuestos; diferencia en >= 2 de 3 celdas): "
          + "; ".join(f"semilla {s}: diferencia {e_count_strict[s][0]}/3, tasa {'si' if e_count_strict[s][1] else 'no'}"
                      for s in seeds) + f" -> {'SE CUMPLE' if okE2 else 'NO SE CUMPLE'}")
        verdicts[reg] = v
        # picker, SC, tiempos
        P("")
        P(f"### Tabla 5 ({reg}). Picker STA/LTA (congelado) y degrade_y en las ventanas de evaluacion; correcciones de estacion")
        P("")
        P("| Celda (semilla) | STA/LTA: arribos picados a 0.3 s (P / S) | STA/LTA falsos por ventana | DY picados | DY falsos por ventana | "
          "rms(static) | rms(SC - static) | rms(SC - static - desajuste medio) | eventos sorteados por 20 catalogados |")
        P("|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|")
        for s in seeds:
            for cell in CELLS:
                c = R[s][cell]
                p1, p2 = c["picker"]["STALTA"], c["picker"]["DY"]
                st = np.concatenate([np.ravel(nw["statics"]) for nw in c["networks"]])
                rr = np.sqrt(np.mean([nw["resid_rms"] ** 2 for nw in c["networks"]]))
                re = np.sqrt(np.mean([nw["resid_expected_rms"] ** 2 for nw in c["networks"]]))
                P(f"| {LAB[cell]} ({s}) | {p1['frac_true_picked_0.3s']:.3f} ({p1['frac_P_picked']:.3f} / {p1['frac_S_picked']:.3f}) | "
                  f"{p1['false_picks_per_window']:.2f} | {p2['frac_true_picked_0.3s']:.3f} | {p2['false_picks_per_window']:.2f} | "
                  f"{np.sqrt(np.mean(st ** 2)):.3f} | {rr:.3f} | {re:.3f} | "
                  + ", ".join(str(nw["cat_drawn"]) for nw in c["networks"]) + " |")
        P("")
        P(f"QuakeMigrate + SC ({reg}): fraccion de entradas del LUT corregido con tiempo negativo (recortado a 0 por la "
          "libreria C), todas las profundidades / profundidades 5-13 km: " + "; ".join(
              f"{LAB[cell]} ({s}) " + "/".join(f"{x:.1e}" for x in np.mean(
                  [lut_clip_fraction(nw["STA"], nw["corr"]) for nw in R[s][cell]["networks"]], axis=0))
              for s in seeds for cell in CELLS))
    P("")
    P("## Veredictos")
    P("")
    P("| Regla | " + " | ".join(REGIMES) + " |")
    P("|:--|" + ":-:|" * len(REGIMES))
    for r_ in ("2A-a", "2A-b", "2A-c", "2A-d", "2A-e"):
        P(f"| {r_} | " + " | ".join("se cumple" if verdicts[g][r_] else "no se cumple" for g in REGIMES) + " |")
    P("")
    a_ok = [g for g in REGIMES if verdicts[g]["2A-a"]]
    if len(a_ok) == 2:
        P("Interpretacion preregistrada: 2A-a se cumple en ambos regimenes: la ventaja de 1y sobrevive al generador "
          "desajustado, al picker real y a QuakeMigrate.")
    elif len(a_ok) == 1:
        P(f"Interpretacion preregistrada: 2A-a se cumple solo en el regimen {a_ok[0]}: la afirmacion se restringe a ese regimen.")
    else:
        P("Interpretacion preregistrada: 2A-a falla en ambos regimenes: la afirmacion de superioridad sobre el estado del "
          "arte se retira de la propuesta y se reporta como resultado negativo.")
    for g in REGIMES:
        if not verdicts[g]["2A-b"]:
            P(f"2A-b falla en {g}: en ese regimen la ganancia no se atribuye al patron de radiacion.")
    P("")
    P("Casos a menos de 0.005 del umbral 0.05: " + ("; ".join(near) if near else "ninguno"))
    txt = "\n".join(L)
    print(txt)
    with open(a.out, "w") as fh:
        fh.write(txt + "\n")
    with open(a.out.replace(".txt", ".json"), "w") as fh:
        json.dump(verdicts, fh, indent=1)


if __name__ == "__main__":
    main()
