"""
task2a_confirm.py - Tarea 2A (preregistrada en PREREGISTRO_1j_1k.txt): benchmark confirmatorio con un generador
desajustado (gen_realistic.generate_real), un picker operativo (real_picker, STA/LTA congelado en la semilla 999),
correcciones de estacion estimadas de un catalogo (station_corrections) y comparadores reales (GaMMA, PyOcto,
QuakeMigrate).

Diseno implementado (ver task2a_report.txt):
  celdas ladder_4, ladder_6, unilateral_c40 (task1r_fair.SCEN); regimenes thin (clutter 0.35) y strong (1.75)
  por semilla y celda 5 redes fijas, compartidas por los dos regimenes (geometria y retardos); por red y regimen un
  catalogo de 20 eventos y N_WIN ventanas de evaluacion independientes.
  Semillas (todas con numpy SeedSequence de entropia en lista; s = semilla, ci = indice de celda en CELLS,
  j = red, ri = 0 thin / 1 strong):
      red            make_network(kind, nst, seed=[s, 2001, ci, j])
      retardos       default_rng([s, 2002, ci, j]).normal(0, 0.3, (nst, 2))        (P, S; fijos por red)
      catalogo       rng [s, 2003, ci, j, ri, 0], rng_m [.., 1], rng_h [.., 2], errores [.., 3]
      ventanas       rng [s, 2004, ci, j, ri, 41] (fuentes, clutter, ruido), rng_m [.., 233] (mecanismos),
                     rng_h [.., 61] (direccion, sitio, funcion fuente, coda), rng_d [.., 43] (degrade_y)
      nodo NODE      rng_n [s, 2005, ci, j, ri] (una permutacion por ventana, como task1y_clutter.run_cell)
      GaMMA          np.random.seed(42) en cada llamada (task1r_baselines.run_gamma)
  La semilla 999 (calibracion de 2A-e y de QuakeMigrate) usa la misma construccion con s = 999: redes distintas
  de las de 2024 y 7331.

    python task2a_confirm.py selfcheck                      # comprobaciones (generador, QM+SC, identidades)
    python task2a_confirm.py qmcal --workers 2              # recalibracion de QM y QM-3C (semilla 999)
    python task2a_confirm.py smoke --nwin 2                 # prueba de humo y estimacion de tiempo (semilla 999)
    python task2a_confirm.py run --seed 2024 --regime thin --workers 2
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MPLBACKEND", "Agg")
import argparse
import contextlib
import copy
import json
import pickle
import shutil
import sys
import tempfile
import time
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "gamma_src"))
warnings.filterwarnings("ignore")

import qm_wrapper as QM  # noqa: E402
import task1x_mtglrt as X  # noqa: E402
from certified_fw import certified_fw_record, dense_fw_record_c  # noqa: E402
from froc_exact import exact_HF, greedy_record  # noqa: E402
from gen_realistic import generate_real  # noqa: E402
from real_picker import CLASSIC_CFG, match, pick_waveforms  # noqa: E402
from station_corrections import estimate_station_corrections, peak_delays  # noqa: E402
from task0_common import DEFAULT_TAUS, DT, VP, VS, default_grid, make_network  # noqa: E402
from task1o_consistent import draw_sources, matched, score_events  # noqa: E402
from task1r_baselines import GAMMA_GRID, PYOCTO_GRID, run_gamma, run_pyocto  # noqa: E402
from task1r_fair import CAPS, SCEN, recall_at  # noqa: E402
from task1y_clutter import CLUT_AMP, degrade_y  # noqa: E402

CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
REGIMES = ("thin", "strong")
N_NETS = 5
N_WIN = 20                     # ventanas de evaluacion por red (regla de presupuesto: 12 si la estimacion > 10 h)
K_CAT = 20
STATIC_SD_2A = 0.3
CAT_ERR = (1.0, 0.1)           # error de catalogo: km por coordenada, s en el tiempo de origen
TAG_NET, TAG_STAT, TAG_CAT, TAG_WIN, TAG_NODE = 2001, 2002, 2003, 2004, 2005
CAL_SEED = 999
QM_RECAL_FILE = os.path.join(HERE, "qm_recal_2a.json")
QM_PER_NET = 6                 # recalibracion de QM: ventanas 0..5 de cada red (30 por celda, el n_cal de qm_wrapper)
SCRATCH = os.environ.get("TASK2A_SCRATCH", os.path.join(tempfile.gettempdir(), "task2a_scratch"))
PARTS = os.path.join(HERE, "task2a_parts")

GAMMA_SETS = ("GaMMA-STALTA", "GaMMA-STALTA+SC", "GaMMA-DY", "GaMMA-DY+SC")
PYOCTO_SETS = ("PyOcto-STALTA", "PyOcto-STALTA+SC")
REC_METHODS = ("MT-FW", "MT-FW+SC", "MT-FW-NODE+SC", "WAVE-PT@20", "QM", "QM-3C", "QM-3C+SC")


# ------------------------------------------------------------------------------------------ seeds
def net_setup(seed, cell):
    ci = CELLS.index(cell)
    kind, nst, nclut = SCEN[cell]
    return ci, kind, nst, nclut


def network(seed, cell, j):
    ci, kind, nst, _nc = net_setup(seed, cell)
    STA = make_network(kind, nst, seed=[seed, TAG_NET, ci, j])
    statics = np.random.default_rng([seed, TAG_STAT, ci, j]).normal(0.0, STATIC_SD_2A, size=(nst, 2))
    return STA, statics


def window_rngs(seed, cell, j, regime):
    ci = CELLS.index(cell)
    ri = REGIMES.index(regime)
    b = [seed, TAG_WIN, ci, j, ri]
    return (np.random.default_rng(b + [41]), np.random.default_rng(b + [233]), np.random.default_rng(b + [61]),
            np.random.default_rng(b + [43]), np.random.default_rng([seed, TAG_NODE, ci, j, ri]))


# -------------------------------------------------------------------------------------- catalogue
def catalogue_real(STA, K, statics, seeds, n_clutter, clut_amp, gate=25.0, min_arr=3, cat_err=CAT_ERR):
    """station_corrections.catalogue con generate_real y los retardos fijos de la red: un evento (fuente uniforme en
    [-15, 15]^2 x [5, 13] km, origen U[5.5, 7.5] s, doble par al azar) entra al catalogo si al menos min_arr de sus
    arribos estacion/fase superan la compuerta de energia en su pico (busqueda +-1 s alrededor del tiempo homogeneo
    desde la fuente verdadera); la localizacion y el origen del catalogo llevan errores N(0, cat_err)."""
    rng = np.random.default_rng(list(seeds) + [0])
    rng_m = np.random.default_rng(list(seeds) + [1])
    rng_h = np.random.default_rng(list(seeds) + [2])
    rng_e = np.random.default_rng(list(seeds) + [3])
    out, truth, drawn = [], [], 0
    while len(out) < K:
        drawn += 1
        src = np.array([rng.uniform(-15, 15), rng.uniform(-15, 15), rng.uniform(5, 13)])
        tau = float(rng.uniform(5.5, 7.5))
        D, _p, _tr, _a = generate_real(STA, [src], [tau], n_clutter, rng, rng_m, rng_h, clut_amp, statics=statics)
        _dl, en = peak_delays(STA, D, src, tau)
        if (en >= gate).sum() < min_arr:
            continue
        src_c = src + rng_e.normal(0, cat_err[0], 3)
        src_c[2] = float(np.clip(src_c[2], 1.0, 20.0))
        out.append((D, src_c, tau + float(rng_e.normal(0, cat_err[1]))))
        truth.append((src, tau))
    return out, drawn, truth


def mismatch_offsets(STA, cat_true):
    """desajuste medio del modelo homogeneo por estacion y fase sobre las fuentes (verdaderas) del catalogo: la
    parte del retardo que no es el static (referencia para interpretar corr - statics)."""
    from gen_realistic import travel_time_gradient
    off = []
    for src in cat_true:
        d = np.linalg.norm(STA - src, axis=1)
        off.append(np.stack([travel_time_gradient(src, STA, "P") - d / VP, travel_time_gradient(src, STA, "S") - d / VS], 1))
    return np.median(np.array(off), axis=0)


# ------------------------------------------------------------------------------------ QuakeMigrate
def load_qm_cfgs():
    if os.environ.get("TASK2A_QM_FALLBACK") == "1":          # solo para medir tiempos antes de la recalibracion
        return {"QM": dict(QM.QM_CFG), "QM-3C": dict(QM.QM_CFG_3C)}
    with open(QM_RECAL_FILE) as fh:
        d = json.load(fh)
    return {k: {kk: (tuple(v) if isinstance(v, list) else v) for kk, v in d[k]["cfg"].items()} for k in ("QM", "QM-3C")}


def qm_sc_lut(lut, corr, names):
    """copia del LUT de QuakeMigrate con las correcciones sumadas a los tiempos de viaje por estacion y fase
    (lut.traveltimes[estacion]['P' | 'S'], las tablas que serve_traveltimes convierte a muestras)."""
    l2 = copy.deepcopy(lut)
    for i, nm in enumerate(names):
        for b, ph in enumerate(("P", "S")):
            l2.traveltimes[nm][ph] = l2.traveltimes[nm][ph] + float(corr[i, b])
    return l2


def run_qm_variants(STA, D, tmpdir, variants, lut=None):
    """QuakeMigrate (detect + trigger publicos) para varias configuraciones sobre la misma ventana: el archivo
    MiniSEED y el LUT se construyen una vez (ambos dependen solo de la ventana, la red y la rejilla, iguales en
    todas las variantes). variants: {nombre: (cfg, corr o None)}. Devuelve {nombre: eventos}."""
    cfgs = [{**QM.QM_CFG, **c} for c, _ in variants.values()]
    for k in ("grid_x", "grid_y", "grid_z", "node_spacing", "pad_s", "pad_seed"):
        assert all(c[k] == cfgs[0][k] for c in cfgs), k
    os.makedirs(tmpdir, exist_ok=True)
    work = tempfile.mkdtemp(prefix="qm2a_", dir=tmpdir)
    try:
        QM.write_archive(D, STA, os.path.join(work, "archive"), cfgs[0])
        if lut is None:
            lut = QM.build_lut(STA, VP, VS, cfgs[0])
        names = QM.station_names(len(STA))
        out = {}
        for (name, (_c, corr)), cfg in zip(variants.items(), cfgs):
            lu = lut if corr is None else qm_sc_lut(lut, corr, names)
            QM.detect(STA, work, lu, cfg)
            out[name] = QM.run_trigger(work, lu, cfg)
        return out
    finally:
        shutil.rmtree(work, ignore_errors=True)


# -------------------------------------------------------------------------------------- utilities
@contextlib.contextmanager
def silent():
    with open(os.devnull, "w") as dn, contextlib.redirect_stdout(dn), contextlib.redirect_stderr(dn):
        yield


def pyocto_guarded(STA, picks, n_picks, tol):
    """task1r_baselines.run_pyocto, salvo que con menos de n_picks picks en total devuelve [] sin llamar a PyOcto:
    PyOcto 0.2.0 termina con un fallo de segmentacion (codigo C++) cuando la lista de picks esta vacia (hallado en
    la semilla 999: una ventana sin ningun pick STA/LTA). Con menos picks que el minimo por evento n_picks no puede
    haber ningun evento, asi que el resultado es el mismo que el de PyOcto en toda entrada que no falla."""
    if sum(len(pl) for pl in picks) < int(n_picks):
        return []
    return run_pyocto(STA, picks, n_picks=n_picks, tol=tol)


def apply_sc(picks, corr):
    """correcciones restadas a los tiempos de los picks de la estacion i y la fase."""
    return [[(float(t) - float(corr[i, 0 if b == "P" else 1]), b) for t, b in pl] for i, pl in enumerate(picks)]


def nodemap_fn(garr, rng_n):
    m = rng_n.permutation(len(garr))
    for gi in range(len(garr)):
        tries = 0
        while np.linalg.norm(garr[m[gi], :2] - garr[gi, :2]) < 9.0 and tries < 50:
            m[gi] = rng_n.integers(len(garr))
            tries += 1
    return m


def pick_counts(true_arr, picks):
    """por ventana: arribos verdaderos, picados (uno a uno por fase a menos de 0.3 s), por fase, picks, falsos."""
    nt = nh = npk = 0
    ntp = nhp = 0
    for tl, pl in zip(true_arr, picks):
        m = match(tl, pl)
        nt += len(tl)
        nh += len(m)
        npk += len(pl)
        ntp += sum(1 for _t, b in tl if b == "P")
        nhp += sum(1 for j, _k, _e in m if tl[j][1] == "P")
    return {"n_true": nt, "n_hit": nh, "n_true_P": ntp, "n_hit_P": nhp, "n_picks": npk, "n_false": npk - nh}


def generate_window(STA, statics, nclut, clut_amp, rngs):
    rng, rng_m, rng_h = rngs[:3]
    srcs, taus = draw_sources(rng)
    D, picks, tr, amp = generate_real(STA, srcs, taus, nclut, rng, rng_m, rng_h, clut_amp, statics=statics)
    return srcs, taus, D, picks, tr, amp


def stream_2a(seed, regime, cell, per_net):
    """ventanas 0..per_net-1 de cada red del diseno 2A (mismos flujos aleatorios que run_job): [(STA, D, truth, taus)]."""
    out = []
    _ci, _kind, _nst, nclut = net_setup(seed, cell)
    for j in range(N_NETS):
        STA, statics = network(seed, cell, j)
        rngs = window_rngs(seed, cell, j, regime)
        for _w in range(per_net):
            srcs, taus, D, _p, tr, _a = generate_window(STA, statics, nclut, CLUT_AMP[regime], rngs)
            out.append((STA, D, tr, list(taus)))
    return out


# ------------------------------------------------------------------------------------------ one job
def run_job(args):
    """una red de una celda, semilla y regimen: catalogo, correcciones y n_win ventanas con todos los metodos."""
    seed, regime, cell, j, n_win, check, out_path = args
    t_start = time.time()
    ci, kind, nst, nclut = net_setup(seed, cell)
    ri = REGIMES.index(regime)
    clut_amp = CLUT_AMP[regime]
    STA, statics = network(seed, cell, j)
    grid = default_grid()
    garr = np.asarray(grid)
    taus_grid = DEFAULT_TAUS
    qmc = load_qm_cfgs()
    tim = {}

    def tick(key, t0):
        tim.setdefault(key, []).append(time.time() - t0)

    t0 = time.time()
    cat, drawn, cat_truth = catalogue_real(STA, K_CAT, statics, [seed, TAG_CAT, ci, j, ri], nclut, clut_amp)
    corr, nacc = estimate_station_corrections(STA, cat)
    expected = statics + mismatch_offsets(STA, [c[0] for c in cat_truth])
    tick("catalogue+SC", t0)
    t0 = time.time()
    Gs = [X.radiation_operator(g, STA) for g in garr]
    Qs = [X.basis(G) for G in Gs]
    lut = QM.build_lut(STA, VP, VS, {**QM.QM_CFG, **qmc["QM-3C"]})
    tick("operators+LUT", t0)
    rngs = window_rngs(seed, cell, j, regime)
    rng_d, rng_n = rngs[3], rngs[4]
    tmpdir = os.path.join(SCRATCH, f"s{seed}_{regime}_{cell}_{j}")
    W = []
    checks = {"fw_identity": 0, "fw_sc_identity": 0, "node_sc_identity": 0, "qm_identity": 0}
    for w in range(n_win):
        rec = {}
        t0 = time.time()
        srcs, taus, D, picks, tr, amp = generate_window(STA, statics, nclut, clut_amp, rngs)
        n2 = 2 * len(srcs)
        true_arr = [pl[:n2] for pl in picks]
        tick("generate", t0)
        t0 = time.time()
        pk_st, _q = pick_waveforms(D, method="stalta")
        tick("picker", t0)
        t0 = time.time()
        p_dy, _lab = degrade_y(picks, amp, len(srcs), rng_d)
        p_dy = [[(t, b) for t, b in pl if t > -100] for pl in p_dy]
        tick("degrade_y", t0)
        psets = {"STALTA": pk_st, "STALTA+SC": apply_sc(pk_st, corr), "DY": p_dy, "DY+SC": apply_sc(p_dy, corr)}
        per = {}
        t0 = time.time()
        with silent():
            for m in GAMMA_SETS:
                P = psets[m.split("-", 1)[1]]
                per[m] = [list(score_events(run_gamma(STA, P, **cfg), tr)) for cfg in GAMMA_GRID]
        tick("gamma(4 sets x 10)", t0)
        t0 = time.time()
        with silent():
            for m in PYOCTO_SETS:
                P = psets[m.split("-", 1)[1]]
                per[m] = [list(score_events(pyocto_guarded(STA, P, **cfg), tr)) for cfg in PYOCTO_GRID]
        tick("pyocto(2 sets x 15)", t0)
        t0 = time.time()
        MF, sd = matched(D)
        MFn = MF / sd
        rec["MT-FW"] = certified_fw_record(STA, MFn, grid, taus_grid, Gs, Qs, tr)
        rec["MT-FW+SC"] = certified_fw_record(STA, MFn, grid, taus_grid, Gs, Qs, tr, tcorr=corr)
        nm = nodemap_fn(garr, rng_n)
        Gn, Qn = [Gs[k] for k in nm], [Qs[k] for k in nm]
        rec["MT-FW-NODE+SC"] = certified_fw_record(STA, MFn, grid, taus_grid, Gn, Qn, tr, tcorr=corr)
        tick("mtfw(3 variants)", t0)
        t0 = time.time()
        f, _m = X.fields_x(STA, D, grid, taus_grid, nm, Qs)
        rec["WAVE-PT@20"] = greedy_record(grid, f["WAVE-PT"], tr, min_sep_km=20.0)
        tick("fields_x", t0)
        t0 = time.time()
        evs = run_qm_variants(STA, D, tmpdir, {"QM": (qmc["QM"], None), "QM-3C": (qmc["QM-3C"], None),
                                                "QM-3C+SC": (qmc["QM-3C"], corr)}, lut=lut)
        for m in ("QM", "QM-3C", "QM-3C+SC"):
            rec[m] = QM.qm_records(evs[m], tr)
        tick("qm(3 variants)", t0)
        if check and w < 2:
            t0 = time.time()
            r0 = X.fw_record(STA, MFn, grid, taus_grid, Gs, Qs, tr)
            assert r0 == rec["MT-FW"], "certified MT-FW differs from fw_record"
            checks["fw_identity"] += 1
            r1 = dense_fw_record_c(STA, MFn, grid, taus_grid, Gs, Qs, tr, tcorr=corr)
            assert r1 == rec["MT-FW+SC"], "certified MT-FW+SC differs from dense_fw_record_c"
            checks["fw_sc_identity"] += 1
            r2 = dense_fw_record_c(STA, MFn, grid, taus_grid, Gn, Qn, tr, tcorr=corr)
            assert r2 == rec["MT-FW-NODE+SC"], "certified MT-FW-NODE+SC differs from dense_fw_record_c"
            checks["node_sc_identity"] += 1
            if w == 0:
                for m in ("QM", "QM-3C"):
                    ev0 = QM.run_quakemigrate(STA, D, tmpdir, VP, VS, qmc[m])
                    assert ev0 == evs[m], f"shared-archive {m} differs from run_quakemigrate"
                checks["qm_identity"] += 1
            tick("identity checks", t0)
        W.append({"k": len(srcs), "truth": [list(map(float, s)) for s in srcs], "taus": [float(t) for t in taus],
                  "records": {m: [[float(s), int(h)] for s, h in r] for m, r in rec.items()},
                  "per": per,
                  "picker": {"STALTA": pick_counts(true_arr, pk_st), "DY": pick_counts(true_arr, p_dy)},
                  "qm_events": {m: [[float(x) for x in e] for e in evs[m]] for m in evs}})
    shutil.rmtree(tmpdir, ignore_errors=True)
    res = {"seed": seed, "regime": regime, "cell": cell, "net": j, "n_win": n_win,
           "seeds": {"network": [seed, TAG_NET, ci, j], "statics": [seed, TAG_STAT, ci, j],
                     "catalogue": [seed, TAG_CAT, ci, j, ri], "windows": [seed, TAG_WIN, ci, j, ri],
                     "nodemap": [seed, TAG_NODE, ci, j, ri]},
           "STA": STA.tolist(), "statics": statics.tolist(), "corr": corr.tolist(), "nacc": nacc.tolist(),
           "cat_drawn": drawn, "cat_src": [c[1].tolist() for c in cat], "cat_tau": [c[2] for c in cat],
           "resid_rms": float(np.sqrt(np.mean((corr - statics) ** 2))),
           "expected_corr": expected.tolist(), "resid_expected_rms": float(np.sqrt(np.mean((corr - expected) ** 2))),
           "cat_true": [[float(x) for x in c[0]] + [float(c[1])] for c in cat_truth],
           "checks": checks, "windows": W, "timing": tim, "seconds": time.time() - t_start}
    if out_path:
        tmp = out_path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(res, fh)
        os.replace(tmp, out_path)
    return res


# ------------------------------------------------------------------------------------ aggregation
def summarize_records(recs, k):
    _t, H, F = exact_HF([[tuple(x) for x in r] for r in recs])
    return {str(c): recall_at(H, F, int(sum(k)), c) for c in CAPS}


def summarize_per(per_w, k):
    ncfg = len(per_w[0])
    nte = int(sum(k))
    pts = []
    for ci in range(ncfg):
        H = sum(pw[ci][0] for pw in per_w)
        F = sum(pw[ci][1] for pw in per_w)
        pts.append((F / len(per_w), H / max(1, nte)))
    ra = {}
    for c in CAPS:
        ok = [rc for fr, rc in pts if fr <= c]
        ra[str(c)] = max(ok) if ok else 0.0
    return pts, ra


def aggregate(seed, regime, parts):
    """une las redes de cada celda en task2a_seed{s}_{regime}.json (lista por celda, como los JSON de 1y)."""
    out = []
    for cell in CELLS:
        ps = sorted([p for p in parts if p["cell"] == cell], key=lambda p: p["net"])
        W = [w for p in ps for w in p["windows"]]
        k = [w["k"] for w in W]
        methods = {}
        for m in REC_METHODS:
            recs = [w["records"][m] for w in W]
            methods[m] = {"records": recs, "recall_at": summarize_records(recs, k)}
        for m, grid_ in [(g, GAMMA_GRID) for g in GAMMA_SETS] + [(p, PYOCTO_GRID) for p in PYOCTO_SETS]:
            per_w = [w["per"][m] for w in W]
            pts, ra = summarize_per(per_w, k)
            methods[m] = {"configs": grid_, "per_window": per_w, "points": pts, "recall_at": ra}
        pk = {}
        for key in ("STALTA", "DY"):
            c = {f: sum(w["picker"][key][f] for w in W) for f in W[0]["picker"][key]}
            pk[key] = {"frac_true_picked_0.3s": c["n_hit"] / max(1, c["n_true"]),
                       "frac_P_picked": c["n_hit_P"] / max(1, c["n_true_P"]),
                       "frac_S_picked": (c["n_hit"] - c["n_hit_P"]) / max(1, c["n_true"] - c["n_true_P"]),
                       "false_picks_per_window": c["n_false"] / len(W), "picks_per_window": c["n_picks"] / len(W),
                       "n_true_arrivals": c["n_true"]}
        tim = {}
        for p in ps:
            for key, v in p["timing"].items():
                tim.setdefault(key, []).extend(v)
        out.append({"scenario": cell, "seed": seed, "regime": regime, "clut_amp": CLUT_AMP[regime],
                    "n_windows": len(W), "k_eval": k, "net_of_window": [p["net"] for p in ps for _ in p["windows"]],
                    "networks": [{key: p[key] for key in ("net", "seeds", "STA", "statics", "corr", "nacc", "cat_drawn",
                                                          "cat_src", "cat_tau", "cat_true", "resid_rms", "expected_corr",
                                                          "resid_expected_rms", "checks", "seconds")}
                                 for p in ps],
                    "windows_meta": [{"truth": w["truth"], "taus": w["taus"], "qm_events": w["qm_events"]} for w in W],
                    "methods": methods, "picker": pk,
                    "timing_mean_s": {key: float(np.mean(v)) for key, v in tim.items()},
                    "seconds": float(sum(p["seconds"] for p in ps))})
    return out


def run_all(seed, regimes, workers, n_win, tag=""):
    os.makedirs(PARTS, exist_ok=True)
    jobs = []
    for regime in regimes:
        for cell in CELLS:
            for j in range(N_NETS):
                path = os.path.join(PARTS, f"seed{seed}_{regime}_{cell}_net{j}{tag}.json")
                if os.path.exists(path):
                    continue
                jobs.append((seed, regime, cell, j, n_win, True, path))
    # las celdas de 8 estaciones primero (mas lentas), para equilibrar los dos procesos
    jobs.sort(key=lambda a: (-SCEN[a[2]][1], a[1], a[3]))
    print(f"seed {seed}: {len(jobs)} jobs pendientes ({n_win} ventanas por red)", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(workers) as ex:
        for r in ex.map(run_job, jobs):
            ra = {m: summarize_records([w["records"][m] for w in r["windows"]], [w["k"] for w in r["windows"]])
                  for m in ("MT-FW+SC", "QM-3C+SC")}
            print(f"   hecho seed {r['seed']} {r['regime']} {r['cell']} red {r['net']} en {r['seconds']:.0f} s "
                  f"(total {time.time() - t0:.0f} s); resid SC rms {r['resid_rms']:.3f} s; checks {r['checks']} :: "
                  + " | ".join(f"{m} " + "/".join(f"{ra[m][str(c)]:.2f}" for c in CAPS) for m in ra), flush=True)
    for regime in regimes:
        parts = []
        for cell in CELLS:
            for j in range(N_NETS):
                with open(os.path.join(PARTS, f"seed{seed}_{regime}_{cell}_net{j}{tag}.json")) as fh:
                    parts.append(json.load(fh))
        agg = aggregate(seed, regime, parts)
        fn = os.path.join(HERE, f"task2a_seed{seed}_{regime}{tag}.json")
        with open(fn, "w") as fh:
            json.dump(agg, fh)
        for c in agg:
            print(f"   {seed} {regime} {c['scenario']}: " + " | ".join(
                f"{m} " + "/".join(f"{c['methods'][m]['recall_at'][str(x)]:.2f}" for x in CAPS) for m in c["methods"]),
                flush=True)
        print(f"escrito {fn}", flush=True)


# --------------------------------------------------------------------------------- QM recalibration
def qm_recalibrate(workers=2):
    """el procedimiento de qm_wrapper (calibrate_detect + calibrate_trigger: 24 onsets, rejilla de disparo en dos
    fases, criterio recall a <= 0.25 falsos por ventana) sobre ventanas de generate_real de la semilla 999: las
    ventanas 0..5 de cada una de las 5 redes de calibracion de 2A de las celdas de qm_wrapper.CELLS (ladder_6,
    unilateral_c40), regimen thin, 30 por celda (n_cal de qm_wrapper), en lugar de y_stream."""
    out_dir = os.path.join(SCRATCH, "qmcal")
    os.makedirs(out_dir, exist_ok=True)
    summary = {}
    for tag, name, extra in (("default", "QM", None), ("3c", "QM-3C", QM.CHANNELS_3C)):
        pkl = os.path.join(out_dir, f"cal2a_{tag}.pkl")
        if not os.path.exists(pkl):
            cfgs = QM._detect_cfgs(extra)
            jobs, meta = [], {}
            for scen in QM.CELLS.values():
                for k, (S, D, tr, taus) in enumerate(stream_2a(CAL_SEED, "thin", scen, QM_PER_NET)):
                    jobs.append((scen, k, S, D, cfgs, os.path.join(out_dir, "work")))
                    meta[(scen, k)] = (tr, taus)
            res, secs = {}, []
            t0 = time.time()
            with ProcessPoolExecutor(workers) as ex:
                for scen, k, scans, sec in ex.map(QM._cal_worker, jobs):
                    res[(scen, k)] = scans
                    secs.append(sec)
            print(f"   detect {tag}: {len(jobs)} ventanas x {len(cfgs)} onsets en {time.time() - t0:.0f} s", flush=True)
            with open(pkl, "wb") as fh:
                pickle.dump({"cfgs": cfgs, "scans": res, "meta": meta, "sec_per_detect": float(np.mean(secs))}, fh)
        summary[name] = pkl
    return summary


def _trigger_job(args):
    name, pkl, out_json = args
    t0 = time.time()
    rows = QM.calibrate_trigger(pkl, out_json)
    return name, rows, time.time() - t0


def qm_select(pkls, workers=2):
    res = {}
    jobs = [(name, pkl, pkl.replace(".pkl", "_trigger.json")) for name, pkl in pkls.items()]
    with ProcessPoolExecutor(workers) as ex:
        for name, rows, sec in ex.map(_trigger_job, jobs):
            with open(pkls[name], "rb") as fh:
                cfgs = pickle.load(fh)["cfgs"]
            best = max(rows, key=QM._rank)          # la primera fila maxima (orden del archivo), como max() de Python
            dc = cfgs[best["ci"]]
            cfg = {**{k: v for k, v in dc.items()},
                   "trigger_window": (QM.TAU_PRIOR[0] + best["lag"], QM.TAU_PRIOR[1] + best["lag"]),
                   "threshold": best["threshold"], "normalise_coalescence": best["norm"], "smooth_coa": best["smooth"],
                   "marginal_window": best["mw"], "min_event_interval": best["mei"]}
            ties = [r for r in rows if QM._rank(r) == QM._rank(best)]
            pkg = [r for r in rows if tuple(r["bp"]) == (2.0, 16.0) and tuple(r["stalta"]) == (0.2, 1.0)]
            res[name] = {"cfg": {k: (list(v) if isinstance(v, tuple) else v) for k, v in cfg.items()},
                         "best_row": best, "n_rows": len(rows), "n_ties": len(ties),
                         "tied_onsets": sorted({(tuple(r["bp"]), tuple(r["stalta"])) for r in ties}),
                         "package_default_best": (max(pkg, key=QM._rank)["recall_at"] if pkg else None),
                         "trigger_seconds": sec}
            print(f"   {name}: {best['recall_at']} bp {best['bp']} stalta {best['stalta']} lag {best['lag']:.2f} "
                  f"norm {best['norm']} q {best['floor_q']} thr {best['threshold']:.4f} mw/mei {best['mw']}/{best['mei']} "
                  f"smooth {best['smooth']} (empates {len(ties)})", flush=True)
    with open(QM_RECAL_FILE, "w") as fh:
        json.dump(res, fh, indent=1, default=list)
    return res


# ------------------------------------------------------------------------------------------ selfcheck
def selfcheck():
    """(1) generate_real con statics: la secuencia de rng_h no cambia y los tiempos de pick llevan los retardos
    dados; (2) QuakeMigrate + SC: un sintetico homogeneo con retardos conocidos se realinea al sumar esos retardos
    al LUT; (3) los tiempos servidos por el LUT corregido difieren en rint(corr * fs)."""
    from gen_realistic import travel_time_gradient
    from task0_common import generate
    ok = True
    # (1)
    STA, statics = network(CAL_SEED, "ladder_6", 0)
    r1 = [np.random.default_rng(7 + i) for i in range(3)]
    r2 = [np.random.default_rng(7 + i) for i in range(3)]
    srcs, taus = [np.array([2.0, -3.0, 9.0])], [6.2]
    _D1, p1, _t1, _a1 = generate_real(STA, srcs, taus, 35, *r1, 0.35)
    _D2, p2, _t2, _a2 = generate_real(STA, srcs, taus, 35, *r2, 0.35, statics=statics)
    same_seq = all(a.random() == b.random() for a, b in zip(r1, r2))
    e = 0.0
    for i in range(len(STA)):
        for bi, b in enumerate("PS"):
            tt = travel_time_gradient(srcs[0], STA[i], b) + taus[0] + statics[i, bi]
            e = max(e, abs(p2[i][bi][0] - tt))
    clut_same = all(x[2:] == y[2:] for x, y in zip(p1, p2))
    print(f"(1) statics: secuencias rng/rng_m/rng_h iguales tras la llamada: {same_seq}; |t_pick - (t_grad + tau + "
          f"static)| max {e:.2e} s; picks de clutter identicos: {clut_same}")
    ok &= same_seq and e < 1e-9 and clut_same
    # (2) QM + SC: homogeneo (task0_common.generate), una fuente, sin clutter, ruido 1/5, retardos multiplos de 0.5 s
    #     (enteros a 16 Hz y a 50 Hz)
    qmc = load_qm_cfgs() if os.path.exists(QM_RECAL_FILE) else {"QM-3C": QM.QM_CFG_3C}
    cfg = {**QM.QM_CFG, **qmc["QM-3C"]}
    names = QM.station_names(6)
    tmp = os.path.join(SCRATCH, "selfcheck")
    rows = []
    for k in range(4):
        S = make_network("random", 6, seed=[CAL_SEED, 2099, k])
        st = np.random.default_rng([CAL_SEED, 2098, k]).choice([-1.0, -0.5, 0.5, 1.0], size=(6, 2))
        src = np.array([[1.0, 2.0, 9.0], [-6.0, 4.0, 7.0], [5.0, -5.0, 11.0], [-2.0, -8.0, 6.0]][k])
        tau = 6.3
        D0, _p, _tr = generate(S, [src], [tau], n_clutter=0, rng=np.random.default_rng([CAL_SEED, 2097, k]),
                               noise_scale=0.004)
        D1, _p, _tr = generate(S, [src], [tau], n_clutter=0, rng=np.random.default_rng([CAL_SEED, 2097, k]),
                               noise_scale=0.004, statics=st)
        lut = QM.build_lut(S, VP, VS, cfg)
        lsc = qm_sc_lut(lut, st, names)
        d_tt = lsc.serve_traveltimes(50) - lut.serve_traveltimes(50)
        want = np.rint(np.concatenate([st[:, 0], st[:, 1]]) * 50)          # serve order: phase-major, then station
        tt_ok = bool(np.all(d_tt == want[None, None, None, :]))
        res = {}
        for tag, D, lu in (("clean", D0, lut), ("static,noSC", D1, lut), ("static,SC", D1, lsc)):
            os.makedirs(tmp, exist_ok=True)
            work = tempfile.mkdtemp(prefix="qmsc_", dir=tmp)
            QM.write_archive(D, S, os.path.join(work, "archive"), cfg)
            QM.detect(S, work, lu, cfg)
            sc_ = QM.read_scan(work, lu, cfg)
            ev = QM.run_trigger(work, lu, cfg)
            shutil.rmtree(work, ignore_errors=True)
            j = int(np.argmax(sc_["COA"].to_numpy()))
            top = ev[0] if ev else None
            res[tag] = (float(sc_["COA"].iloc[j]), float(np.hypot(sc_["x"].iloc[j] - src[0], sc_["y"].iloc[j] - src[1])),
                        float(sc_["t"].iloc[j]), top, sc_["COA"].to_numpy())
        dtrace = float(np.max(np.abs(res["clean"][4] - res["static,SC"][4])) / np.max(res["clean"][4]))
        rows.append((k, tt_ok, res, dtrace))
        print(f"(2) ventana {k}: LUT corregido = LUT + rint(static*50) en todas las celdas: {tt_ok}; "
              + " | ".join(f"{t}: COA max {v[0]:.3f} a {v[1]:.1f} km, t {v[2]:.2f} s" for t, v in res.items())
              + f" | max |COA_clean - COA_SC| / max COA_clean = {dtrace:.3f}")
        ok &= tt_ok and res["static,SC"][1] <= 2.0 and res["static,SC"][0] >= 0.98 * res["clean"][0] and \
            res["static,SC"][0] > res["static,noSC"][0]
    print("selfcheck", "correcto" if ok else "CON FALLAS")
    return ok, rows


# --------------------------------------------------------------------------------------------- smoke
def smoke(nwin, workers):
    """semilla 999, nwin ventanas por celda y regimen (red 0), todos los metodos; tiempo por componente y total
    extrapolado a las semillas 2024 y 7331 (2 semillas x 2 regimenes x 3 celdas x 5 redes x N_WIN)."""
    jobs = [(CAL_SEED, reg, cell, 0, nwin, True, None) for reg in REGIMES for cell in CELLS]
    t0 = time.time()
    with ProcessPoolExecutor(workers) as ex:
        res = list(ex.map(run_job, jobs))
    wall = time.time() - t0
    lines = [f"prueba de humo: semilla 999, red 0, {nwin} ventanas por celda y regimen, {workers} procesos, "
             f"{wall:.0f} s de pared"]
    per_cell = {}
    for r in res:
        ww = r["timing"]
        comp = {k: float(np.mean(v)) for k, v in ww.items() if k not in ("catalogue+SC", "operators+LUT", "identity checks")}
        per_w = sum(comp.values())
        fixed = sum(float(np.sum(ww[k])) for k in ("catalogue+SC", "operators+LUT"))
        idc = float(np.sum(ww.get("identity checks", [0.0])))
        per_cell[(r["regime"], r["cell"])] = (per_w, fixed, idc)
        lines.append(f"   {r['regime']:6s} {r['cell']:15s} por ventana {per_w:6.2f} s ["
                     + ", ".join(f"{k} {v:.2f}" for k, v in comp.items())
                     + f"]; por red: catalogo+SC y operadores+LUT {fixed:.1f} s, comprobaciones de identidad {idc:.1f} s; "
                     f"checks {r['checks']}")
    for nw in (20, 12):
        tot = 0.0
        for (reg, cell), (pw, fixed, idc) in per_cell.items():
            tot += N_NETS * (nw * pw + fixed + idc)
        tot *= 2                                   # dos semillas
        lines.append(f"   estimacion con {nw} ventanas por red: {tot / 3600:.2f} h de CPU, {tot / 3600 / 2:.2f} h en 2 nucleos")
    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(HERE, "task2a_smoke.txt"), "w") as fh:
        fh.write(txt + "\n")
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["selfcheck", "qmcal", "smoke", "run"])
    ap.add_argument("--seed", type=int, default=CAL_SEED)
    ap.add_argument("--regime", default="both", choices=["thin", "strong", "both"])
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--nwin", type=int, default=N_WIN)
    a = ap.parse_args()
    os.makedirs(SCRATCH, exist_ok=True)
    if a.stage == "selfcheck":
        ok, _r = selfcheck()
        raise SystemExit(0 if ok else 1)
    if a.stage == "qmcal":
        pk = qm_recalibrate(a.workers)
        qm_select(pk, a.workers)
    elif a.stage == "smoke":
        smoke(a.nwin, a.workers)
    else:
        regs = REGIMES if a.regime == "both" else (a.regime,)
        run_all(a.seed, regs, a.workers, a.nwin)


if __name__ == "__main__":
    main()
