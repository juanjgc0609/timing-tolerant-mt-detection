"""
real_picker.py - pickers reales sobre las formas de onda (no el picker parametrico degrade_y).

pick_waveforms(D, dt=DT, cfg=None, method="phasenet") -> (picks, calidad): picks por estacion [(t, 'P'|'S')] en el
formato de task1r_baselines.run_gamma, y por pick la probabilidad de PhaseNet o el maximo de STA/LTA.

1. method="phasenet": PhaseNet (Zhu y Beroza 2019) con el modelo preentrenado ORIGINAL del repositorio
   AI4EPS/PhaseNet (model/190703-214543, checkpoint model_95.ckpt, datos del norte de California), en
   TensorFlow-CPU 2.15.1 (la ultima version con Keras 2, que tf.compat.v1.layers necesita) con el codigo de la red
   del propio repositorio (phasenet/model.py), commit 62005c6195638f88d077f8870ea895971062715a. Sin reentrenar.
   Instalacion fuera de /home/claude/work, para no tocar el entorno compartido (TF 2.15 exige protobuf < 5 y el
   sistema tiene protobuf 7, que usa onnxruntime):
       git clone --depth 1 https://github.com/AI4EPS/PhaseNet.git ~/.cache/phasenet/PhaseNet
       pip install --target ~/.cache/phasenet/tf215 tensorflow-cpu==2.15.1
   tf215 se antepone a sys.path solo mientras se importa TensorFlow y se construye la red; despues pasa al final.
   PHASENET_HOME cambia la raiz (por defecto ~/.cache/phasenet); si falta el repositorio se clona desde GitHub.
   Preproceso como phasenet/predict.py con datos continuos: componentes E, N, Z con E = x, N = y, Z = -z (z es
   profundidad en task0_common); remuestreo de 16 a 100 Hz con resample_poly (25/4); normalize_long de PhaseNet
   (ventanas de 3000 muestras, paso 1500); la traza de 4000 muestras entra entera (red convolucional; da los
   mismos maximos que ventanas de 3000). Posproceso: maximos de la probabilidad P y S con altura >= umbral y
   separacion >= mpd (el detect_peaks de PhaseNet, aqui scipy.signal.find_peaks).
   Comprobaciones: en las muestras etiquetadas del repositorio (dataset/waveform_train) los maximos caen a 0-10
   muestras de itp/its con probabilidad 0.6-0.98; en ventanas de generate_real SIN ruido da 0.8-0.9 en casi
   todos los arribos. Con el ruido del generador (blanco, 0-8 Hz tras remuestrear, por encima de la banda de
   1.4-2.8 Hz de la senal) la red casi no responde: necesita amplitud/ruido de ~8-10 para detectar la mitad de
   los arribos (degrade_y supone 2). Un prefiltro (pasabanda 1-4, 1-6, 0.7-3.5 Hz o pasabajos 4 Hz) no mejora
   el F1 en calibracion; no se usa.
2. method="stalta": picker operacional clasico. Pasabanda Butterworth de orden 4 de fase cero, remuestreo a
   100 Hz, STA/LTA clasico (classic_sta_lta + trigger_onset de ObsPy) sobre la vertical para P y sobre la norma
   horizontal para S; inicio en el disparo, o refinado con AIC de Maeda (aic_simple) o Baer-Kradolfer (pk_baer)
   en [-1.5, +0.5] s del disparo; resolucion opcional de disparos P y S simultaneos por la fraccion de energia
   vertical. El refinamiento (AIC o pk_baer) y la resolucion P/S empeoraron el F1 en calibracion: no se usan.

Calibracion (python real_picker.py --mode calib): semilla 999, ventanas evbase 1000, 20 por celda en las 4 celdas
random-6 / onesided-8 x clut_amp 0.35 / 1.75 con n_clutter 35 (80 ventanas, 840 arribos); el informe usa otras
ventanas (evbase 25000, 40 por celda) y nunca las semillas de evaluacion 2024 y 7331. Criterio: F1 de arribos
verdaderos emparejados uno a uno, por fase, dentro de 0.3 s. Desfase por fase = -(mediana de t_pick - t_verdadero)
a menos de 1 s, refinado una vez con los emparejados: los pickers marcan el INICIO de la fase y el tiempo verdadero
del generador es el CENTRO del pulso de Ricker (PhaseNet pica ~0.35 s antes del centro).
    PhaseNet: umbral P 0.25 y S 0.03 (rejilla 0.02, 0.03, 0.05-0.95), desfase P +0.353 s (solo 5 residuos) y
              S +0.328 s, mpd 0.5 s (valor por defecto de PhaseNet, sin ajustar). F1 de calibracion P 0.025, S 0.135.
    STA/LTA:  banda 1-4 Hz (frente a 1-6 Hz); P: STA 0.3 s, LTA 8 s, disparo 6.0; S: STA 0.6 s, LTA 8 s, disparo
              3.0 (rejilla STA 0.3/0.6, LTA 4/8, disparo 2-8, refinamiento none/aic/baer); fin de disparo 1.0;
              desfase P +0.026 s, S +0.016 s; sin resolucion P/S. F1 de calibracion P 0.126, S 0.193.

    python real_picker.py --mode calib     calibracion (imprime PICK_CFG y CLASSIC_CFG)
    python real_picker.py --mode report    estadisticas del informe, frente a degrade_y en las mismas ventanas
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
import argparse
import importlib.util
import json
import logging
import subprocess
import sys
import time
import warnings
from fractions import Fraction

import numpy as np
from obspy.signal.trigger import aic_simple, classic_sta_lta, pk_baer, trigger_onset
from scipy.interpolate import interp1d
from scipy.signal import butter, find_peaks, resample_poly, sosfiltfilt

from task0_common import DT, T, make_network
from task1o_consistent import NOISE_SD, draw_sources

PN_HOME = os.environ.get("PHASENET_HOME", os.path.expanduser("~/.cache/phasenet"))
PN_REPO = os.path.join(PN_HOME, "PhaseNet")
PN_URL = "https://github.com/AI4EPS/PhaseNet.git"
PN_COMMIT = "62005c6195638f88d077f8870ea895971062715a"
PN_MODEL = os.path.join(PN_REPO, "model", "190703-214543")
PN_TF = os.path.join(PN_HOME, "tf215")
TF_REQUIRED = "2.15.1"
FS = 100.0
WINDOW = 3000
TOL_PICK = 0.3

PICK_CFG = {"p_thr": 0.25, "s_thr": 0.03, "p_shift": 0.353, "s_shift": 0.328, "mpd": 0.5}
CLASSIC_CFG = {"band": (1.0, 4.0), "off": 1.0, "merge": 0.5, "fz": None,
               "p_sta": 0.3, "p_lta": 8.0, "p_on": 6.0, "p_refine": "none", "p_shift": 0.026,
               "s_sta": 0.6, "s_lta": 8.0, "s_on": 3.0, "s_refine": "none", "s_shift": 0.016}

CELLS = [("random", 6, 35, 0.35), ("random", 6, 35, 1.75), ("onesided", 8, 35, 0.35), ("onesided", 8, 35, 1.75)]
SEED, MOFF, EV_CAL, EV_REP, N_CAL, N_REP = 999, 233, 1000, 25000, 20, 40

_PN = None


def _ensure_repo():
    if os.path.isfile(os.path.join(PN_MODEL, "checkpoint")):
        return
    os.makedirs(PN_HOME, exist_ok=True)
    subprocess.run(["git", "clone", "--depth", "1", PN_URL, PN_REPO], check=True)


def load_phasenet():
    """sesion de TensorFlow con la red de PhaseNet y el checkpoint restaurado; se carga una vez por proceso."""
    global _PN
    if _PN is not None:
        return _PN
    _ensure_repo()
    added = os.path.isdir(PN_TF) and PN_TF not in sys.path
    if added:
        sys.path.insert(0, PN_TF)
    try:
        import tensorflow as tf
    except ImportError as e:
        raise RuntimeError(f"falta TensorFlow: pip install --target {PN_TF} tensorflow-cpu=={TF_REQUIRED}") from e
    tf.get_logger().setLevel("ERROR")
    logging.getLogger("tensorflow").setLevel(logging.ERROR)
    tf.compat.v1.disable_eager_execution()
    spec = importlib.util.spec_from_file_location("phasenet_model", os.path.join(PN_REPO, "phasenet", "model.py"))
    mod = importlib.util.module_from_spec(spec)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        spec.loader.exec_module(mod)
        graph = tf.Graph()
        with graph.as_default():
            model = mod.UNet(config=mod.ModelConfig(X_shape=[WINDOW, 1, 3]), mode="pred")
            sess = tf.compat.v1.Session(graph=graph, config=tf.compat.v1.ConfigProto(
                intra_op_parallelism_threads=1, inter_op_parallelism_threads=1))
            tf.compat.v1.train.Saver(tf.compat.v1.global_variables()).restore(sess, tf.train.latest_checkpoint(PN_MODEL))
    if added:
        sys.path.remove(PN_TF)
        sys.path.append(PN_TF)
    _PN = (sess, model, tf.__version__)
    return _PN


def normalize_long(data, window=WINDOW):
    """normalize_long de phasenet/data_reader.py para una traza (nt, 1, nch)."""
    nt, nar, nch = data.shape
    shift = window // 2
    pad = np.pad(data, ((window // 2, window // 2), (0, 0), (0, 0)), mode="reflect")
    t = np.arange(0, nt, shift, dtype="int")
    std = np.zeros([len(t) + 1, nar, nch])
    mean = np.zeros([len(t) + 1, nar, nch])
    for i in range(1, len(std)):
        std[i] = np.std(pad[i * shift:i * shift + window], axis=0)
        mean[i] = np.mean(pad[i * shift:i * shift + window], axis=0)
    t = np.append(t, nt)
    std[-1], mean[-1] = std[-2], mean[-2]
    std[0], mean[0] = std[1], mean[1]
    ti = np.arange(nt, dtype="int")
    si = interp1d(t, std, axis=0, kind="slinear")(ti)
    mi = interp1d(t, mean, axis=0, kind="slinear")(ti)
    tmp = np.sum(si, axis=(0, 1))
    si[si == 0] = 1.0
    out = (data - mi) / si
    nz = np.count_nonzero(tmp)
    if 0 < nz < 3:
        out *= 3.0 / nz
    return out


def phasenet_probs(D, dt=DT):
    """probabilidades de PhaseNet (ns, nt100, 3) [ruido, P, S] a FS Hz para D (ns, 3, nt) en x, y, z."""
    sess, model, _v = load_phasenet()
    fr = Fraction(FS * dt).limit_denominator(1000)
    X = np.stack([D[:, 0], D[:, 1], -D[:, 2]], axis=-1)
    X = resample_poly(X, fr.numerator, fr.denominator, axis=1)
    X = np.stack([normalize_long(x[:, None, :]) for x in X]).astype(np.float32)
    P = sess.run(model.preds, feed_dict={model.X: X, model.drop_rate: 0, model.is_training: False})
    return P[:, :, 0, :]


def probs_to_picks(P, cfg=None):
    """picks por estacion [(t, 'P'|'S')] y sus probabilidades, a partir de phasenet_probs."""
    c = dict(PICK_CFG)
    c.update(cfg or {})
    mpd = max(1, int(round(c["mpd"] * FS)))
    picks, probs = [], []
    for pr in P:
        pl = []
        for k, b in ((1, "P"), (2, "S")):
            idx, prop = find_peaks(pr[:, k], height=c[f"{b.lower()}_thr"], distance=mpd)
            pl += [(float(j / FS + c[f"{b.lower()}_shift"]), b, float(h)) for j, h in zip(idx, prop["peak_heights"])]
        pl.sort()
        picks.append([(t, b) for t, b, _h in pl])
        probs.append([h for _t, _b, h in pl])
    return picks, probs


def classic_traces(D, dt=DT, band=(1.0, 6.0)):
    """pasabanda Butterworth de orden 4 de fase cero y remuestreo a FS: (ns, 3, nt100) en x, y, z."""
    X = sosfiltfilt(butter(4, band, btype="bandpass", fs=1.0 / dt, output="sos"), D, axis=-1)
    fr = Fraction(FS * dt).limit_denominator(1000)
    return resample_poly(X, fr.numerator, fr.denominator, axis=-1)


def _cf(X, phase, sta, lta):
    """STA/LTA clasico: sobre la vertical para P y sobre la norma horizontal para S."""
    x = X[2] if phase == "P" else np.hypot(X[0], X[1])
    return classic_sta_lta(np.ascontiguousarray(x), int(round(sta * FS)), int(round(lta * FS)))


def _refine(X, phase, a, refine, pre=1.5, post=0.5):
    """inicio del arribo cerca del disparo a (muestra): AIC de Maeda (suma sobre las horizontales para S),
    Baer-Kradolfer (pk_baer, en la horizontal de mayor energia para S) o el propio disparo."""
    lo, hi = max(0, a - int(pre * FS)), min(X.shape[1], a + int(post * FS))
    if refine == "none" or hi - lo < 20:
        return a
    comps = [X[2]] if phase == "P" else [X[0], X[1]]
    if refine == "aic":
        s = sum(aic_simple(np.ascontiguousarray(c[lo:hi])) for c in comps)
        return lo + 2 + int(np.argmin(s[2:-2]))
    c = comps[int(np.argmax([np.sum(c[lo:hi] ** 2) for c in comps]))]
    pp, _fm = pk_baer(np.ascontiguousarray(c[lo:hi], dtype=np.float32), int(FS), 20, 30, 7.0, 12.0, 100, 100)
    return lo + int(pp) if pp > 0 else a


def _candidates(X, cf, phase, on, off, refine):
    return [(_refine(X, phase, int(a), refine), float(cf[a:max(a + 1, b)].max())) for a, b in trigger_onset(cf, on, off)]


def classic_picks(D, dt=DT, cfg=None):
    """STA/LTA (classic_sta_lta + trigger_onset) en la vertical (P) y en la horizontal (S), inicio refinado, desfase
    por fase y resolucion de disparos P y S simultaneos (a menos de 'merge' s) por la fraccion de energia vertical
    en los 0.5 s siguientes: P si es >= fz, S si no. Devuelve picks por estacion y el maximo de STA/LTA por pick."""
    c = dict(CLASSIC_CFG)
    c.update(cfg or {})
    X = classic_traces(D, dt, tuple(c["band"]))
    picks, qual = [], []
    for Xi in X:
        cand = {b: _candidates(Xi, _cf(Xi, b, c[f"{b.lower()}_sta"], c[f"{b.lower()}_lta"]), b,
                               c[f"{b.lower()}_on"], c["off"], c[f"{b.lower()}_refine"]) for b in ("P", "S")}
        picks_i, qual_i = _merge_phases(Xi, cand, c)
        picks.append(picks_i)
        qual.append(qual_i)
    return picks, qual


def _merge_phases(Xi, cand, c):
    keep = {b: [True] * len(cand[b]) for b in ("P", "S")}
    if c["fz"] is not None:
        w = int(0.5 * FS)
        for ip, (kp, _q) in enumerate(cand["P"]):
            for is_, (ks, _q2) in enumerate(cand["S"]):
                if abs(kp - ks) <= c["merge"] * FS:
                    k = min(kp, ks)
                    ez = float(np.sum(Xi[2, k:k + w] ** 2))
                    eh = float(np.sum(Xi[:2, k:k + w] ** 2))
                    if ez / max(ez + eh, 1e-30) >= c["fz"]:
                        keep["S"][is_] = False
                    else:
                        keep["P"][ip] = False
    pl = sorted((k / FS + c[f"{b.lower()}_shift"], b, q) for b in ("P", "S")
                for (k, q), ok in zip(cand[b], keep[b]) if ok)
    return [(float(t), b) for t, b, _q in pl], [q for _t, _b, q in pl]


def pick_waveforms(D, dt=DT, cfg=None, method="phasenet"):
    """picker real sobre la ventana D (ns, 3, nt): picks por estacion [(t, 'P'|'S')] y un valor de calidad por
    pick (probabilidad de PhaseNet, o maximo de STA/LTA con method='stalta')."""
    if method == "stalta":
        return classic_picks(D, dt, cfg)
    return probs_to_picks(phasenet_probs(D, dt), cfg)


def match(true_st, picks_st, tol=TOL_PICK):
    """emparejamiento uno a uno por fase dentro de tol, en orden creciente de |dt|: [(j, k, t_pick - t_true)]."""
    cand = [(abs(tp - tt), j, k, tp - tt) for j, (tt, bt) in enumerate(true_st)
            for k, (tp, bp) in enumerate(picks_st) if bp == bt and abs(tp - tt) <= tol]
    cand.sort()
    uj, uk, out = set(), set(), []
    for _a, j, k, e in cand:
        if j not in uj and k not in uk:
            uj.add(j)
            uk.add(k)
            out.append((j, k, e))
    return out


def dataset(kind, nst, nclut, clut_amp, seed, evbase, reps, cfg=None, moff=MOFF):
    """ventanas de generate_real con las semillas de run_cell de task1y (degrade_y sobre los mismos arribos)."""
    from gen_realistic import generate_real
    from task1y_clutter import degrade_y
    rng = np.random.default_rng(seed + 41 + evbase)
    rng_d = np.random.default_rng(seed + 43 + evbase)
    rng_m = np.random.default_rng(seed + moff)
    rng_h = np.random.default_rng(seed + 61 + evbase)
    out = []
    for k in range(reps):
        S = make_network(kind, nst, seed=seed + evbase + k)
        srcs, taus = draw_sources(rng)
        D, p, tr, amp = generate_real(S, srcs, taus, nclut, rng, rng_m, rng_h, clut_amp, cfg)
        p2, _lab = degrade_y(p, amp, len(srcs), rng_d)
        n2 = 2 * len(srcs)
        out.append({"D": D, "true": [pl[:n2] for pl in p], "clut": [pl[n2:] for pl in p],
                    "amp": [a[:n2] for a in amp], "degr": [[(t, b) for t, b in pl if t > -100] for pl in p2]})
    return out


def stats(windows, key):
    """fraccion de arribos picados (por fase y por clase de amplitud), error de tiempo de los picks
    verdaderos, picks falsos por ventana y su origen."""
    ntrue = {"P": 0, "S": 0}
    nhit = {"P": 0, "S": 0}
    err, nfalse, nswap, nclut, npk = [], 0, 0, 0, 0
    edges = np.array([0, 1, 2, 3, 5, 10, np.inf])
    ab_n, ab_h = np.zeros(len(edges) - 1), np.zeros(len(edges) - 1)
    for w in windows:
        for i, (tl, pl) in enumerate(zip(w["true"], w[key])):
            m = match(tl, pl)
            hit = {j for j, _k, _e in m}
            used = {k for _j, k, _e in m}
            err += [e for _j, _k, e in m]
            npk += len(pl)
            for j, (_t, b) in enumerate(tl):
                ntrue[b] += 1
                nhit[b] += j in hit
                c = int(np.clip(np.searchsorted(edges, w["amp"][i][j] / NOISE_SD, side="right") - 1, 0, len(ab_n) - 1))
                ab_n[c] += 1
                ab_h[c] += j in hit
            for k, (tp, bp) in enumerate(pl):
                if k in used:
                    continue
                nfalse += 1
                if any(abs(tp - tt) <= TOL_PICK and bt != bp for tt, bt in tl):
                    nswap += 1
                elif any(abs(tp - tc) <= TOL_PICK for tc, _b in w["clut"][i]):
                    nclut += 1
    nw = len(windows)
    ae = np.abs(np.asarray(err)) if err else np.array([np.nan])
    nt = ntrue["P"] + ntrue["S"]
    rec = (nhit["P"] + nhit["S"]) / max(1, nt)
    prec = (nhit["P"] + nhit["S"]) / max(1, npk)
    return {"n_true": nt, "recall": rec, "recall_P": nhit["P"] / max(1, ntrue["P"]),
            "recall_S": nhit["S"] / max(1, ntrue["S"]), "precision": prec,
            "f1": 2 * rec * prec / max(1e-12, rec + prec), "err_median": float(np.median(ae)),
            "err_p90": float(np.percentile(ae, 90)), "bias": float(np.median(err)) if err else float("nan"),
            "false_per_window": nfalse / nw, "swap_per_window": nswap / nw, "clutter_per_window": nclut / nw,
            "other_per_window": (nfalse - nswap - nclut) / nw, "picks_per_window": npk / nw,
            "recall_by_amp": {f"{edges[c]:g}-{edges[c + 1]:g}": (int(ab_n[c]), float(ab_h[c] / max(1, ab_n[c])))
                              for c in range(len(ab_n))}}


def _cell_windows(evbase, reps):
    out = []
    for kind, nst, nclut, ca in CELLS:
        t0 = time.time()
        ws = dataset(kind, nst, nclut, ca, SEED, evbase, reps)
        for w in ws:
            w["P"] = phasenet_probs(w["D"])
        out.append(((kind, nst, nclut, ca), ws))
        print(f"   {kind}-{nst} clutter {ca}: {reps} ventanas en {time.time() - t0:.0f} s", flush=True)
    return out


def calibrate(cells):
    """PhaseNet: umbral y desfase por fase que maximizan el F1 de arribos verdaderos dentro de 0.3 s."""
    ws = [w for _c, wl in cells for w in wl]
    cfg = dict(PICK_CFG)
    thrs = np.concatenate([[0.02, 0.03], np.round(np.arange(0.05, 0.951, 0.05), 2)])
    for b in ("P", "S"):
        lo = b.lower()
        res = []
        base = dict(cfg, **{f"{lo}_thr": 0.3, f"{lo}_shift": 0.0})
        for w in ws:
            pk, _pr = probs_to_picks(w["P"], base)
            for tl, pl in zip(w["true"], pk):
                for tt, bt in tl:
                    if bt != b:
                        continue
                    d = [tp - tt for tp, bp in pl if bp == b and abs(tp - tt) <= 1.0]
                    if d:
                        res.append(min(d, key=abs))
        shift = -float(np.median(res))
        best = None
        for thr in thrs:
            for sh in (shift,):
                c2 = dict(cfg, **{f"{lo}_thr": float(thr), f"{lo}_shift": sh})
                nt = nh = npk = 0
                errs = []
                for w in ws:
                    pk, _pr = probs_to_picks(w["P"], c2)
                    for tl, pl in zip(w["true"], pk):
                        tl_b = [x for x in tl if x[1] == b]
                        pl_b = [x for x in pl if x[1] == b]
                        m = match(tl_b, pl_b)
                        nt += len(tl_b)
                        nh += len(m)
                        npk += len(pl_b)
                        errs += [e for _j, _k, e in m]
                r, p = nh / max(1, nt), nh / max(1, npk)
                f1 = 2 * r * p / max(1e-12, r + p)
                if best is None or f1 > best[0]:
                    best = (f1, float(thr), sh, r, p, float(np.median(errs)) if errs else 0.0)
        f1, thr, sh, r, p, med = best
        cfg[f"{lo}_thr"], cfg[f"{lo}_shift"] = thr, round(sh - med, 3)
        print(f"   {b}: n residuos {len(res)}, desfase inicial {shift:+.3f} s, umbral {thr:.2f}, F1 {f1:.3f} "
              f"(recall {r:.3f}, precision {p:.3f}), mediana del error firmado {med:+.3f} s -> desfase {cfg[f'{lo}_shift']:+.3f} s",
              flush=True)
    return cfg


def _f1_phase(ws, b, cand_times, shift):
    nt = nh = npk = 0
    errs = []
    for w, ct in zip(ws, cand_times):
        for tl, kl in zip(w["true"], ct):
            tb = [x for x in tl if x[1] == b]
            pb = [(k / FS + shift, b) for k in kl]
            m = match(tb, pb)
            nt, nh, npk = nt + len(tb), nh + len(m), npk + len(pb)
            errs += [e for _j, _k, e in m]
    r, p = nh / max(1, nt), nh / max(1, npk)
    return 2 * r * p / max(1e-12, r + p), r, p, (float(np.median(errs)) if errs else 0.0)


def _shift0(ws, b, cand_times):
    res = []
    for w, ct in zip(ws, cand_times):
        for tl, kl in zip(w["true"], ct):
            for tt, bt in tl:
                d = [k / FS - tt for k in kl if bt == b and abs(k / FS - tt) <= 1.0]
                if d:
                    res.append(min(d, key=abs))
    return -float(np.median(res)) if res else 0.0


def calibrate_classic(cells):
    """rejilla por fase: banda, STA, LTA, umbral de disparo, refinamiento; desfase = mediana de los residuos; se
    elige por F1 dentro de 0.3 s. Despues, el umbral fz de la resolucion P/S por F1 conjunto."""
    ws = [w for _c, wl in cells for w in wl]
    cfg = dict(CLASSIC_CFG)
    best_band = None
    for band in ((1.0, 4.0), (1.0, 6.0)):
        Xs = [classic_traces(w["D"], DT, band) for w in ws]
        found = {}
        for b in ("P", "S"):
            best = None
            for sta, lta in ((0.3, 4.0), (0.3, 8.0), (0.6, 4.0), (0.6, 8.0)):
                cfs = [[_cf(Xi, b, sta, lta) for Xi in X] for X in Xs]
                for on in (2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0):
                    for ref in ("none", "aic", "baer"):
                        ct = [[[k for k, _q in _candidates(Xi, cf, b, on, cfg["off"], ref)] for Xi, cf in zip(X, cw)]
                              for X, cw in zip(Xs, cfs)]
                        sh = _shift0(ws, b, ct)
                        f1, r, p, med = _f1_phase(ws, b, ct, sh)
                        if best is None or f1 > best[0]:
                            best = (f1, r, p, {f"{b.lower()}_sta": sta, f"{b.lower()}_lta": lta, f"{b.lower()}_on": on,
                                               f"{b.lower()}_refine": ref, f"{b.lower()}_shift": round(sh - med, 3)})
            found[b] = best
            print(f"   banda {band} {b}: F1 {best[0]:.3f} (recall {best[1]:.3f}, precision {best[2]:.3f}) {best[3]}", flush=True)
        key = (found["P"][0] + found["S"][0]) / 2
        if best_band is None or key > best_band[0]:
            best_band = (key, band, found)
    _k, band, found = best_band
    cfg["band"] = band
    cfg.update(found["P"][3])
    cfg.update(found["S"][3])
    best = None
    for fz in (None, 0.2, 0.3, 0.4, 0.5, 0.6):
        c2 = dict(cfg, fz=fz)
        nt = nh = npk = 0
        for w in ws:
            pk, _q = classic_picks(w["D"], DT, c2)
            for tl, pl in zip(w["true"], pk):
                nt, nh, npk = nt + len(tl), nh + len(match(tl, pl)), npk + len(pl)
        r, p = nh / max(1, nt), nh / max(1, npk)
        f1 = 2 * r * p / max(1e-12, r + p)
        print(f"   fz {fz}: F1 conjunto {f1:.3f} (recall {r:.3f}, precision {p:.3f})", flush=True)
        if best is None or f1 > best[0]:
            best = (f1, fz)
    cfg["fz"] = best[1]
    return cfg


METHODS = ("PhaseNet", "STA/LTA", "degrade_y")


def report(cells, pn_cfg, cl_cfg):
    out = {}
    for (kind, nst, nclut, ca), ws in cells:
        for w in ws:
            w["PhaseNet"], _pr = probs_to_picks(w["P"], pn_cfg)
            w["STA/LTA"], _q = classic_picks(w["D"], DT, cl_cfg)
        out[f"{kind}-{nst} n_clutter {nclut} clut_amp {ca}"] = _entry(ws)
    return out


def _entry(ws):
    r = {m: stats(ws, "degr" if m == "degrade_y" else m) for m in METHODS}
    r.update(n_windows=len(ws), n_events=int(sum(len(w["true"][0]) // 2 for w in ws)))
    return r


def _print(res):
    for tag, r in res.items():
        print(f"\n=== {tag} | ventanas {r['n_windows']}, eventos {r['n_events']}, arribos {r['PhaseNet']['n_true']}")
        for m in METHODS:
            s = r[m]
            print(f"   {m:<9} picados {s['recall']:.3f} (P {s['recall_P']:.3f}, S {s['recall_S']:.3f}) | error "
                  f"mediana {s['err_median']:.3f} s p90 {s['err_p90']:.3f} s sesgo {s['bias']:+.3f} s | falsos/ventana "
                  f"{s['false_per_window']:.2f} (fase cambiada {s['swap_per_window']:.2f}, clutter "
                  f"{s['clutter_per_window']:.2f}, otros {s['other_per_window']:.2f}) | picks/ventana "
                  f"{s['picks_per_window']:.1f} precision {s['precision']:.3f} F1 {s['f1']:.3f}")
            print("             picados por amplitud/ruido: " + "  ".join(
                f"[{k}] {v[1]:.2f} (n={v[0]})" for k, v in s["recall_by_amp"].items()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="report", choices=["calib", "report"])
    ap.add_argument("--json", default="")
    args = ap.parse_args()
    t0 = time.time()
    _s, _m, v = load_phasenet()
    print(f"PhaseNet {os.path.basename(PN_MODEL)} (commit {PN_COMMIT[:12]}), TensorFlow {v}, carga {time.time() - t0:.1f} s")
    if args.mode == "calib":
        cells = _cell_windows(EV_CAL, N_CAL)
        pn_cfg = calibrate(cells)
        print("PICK_CFG =", json.dumps(pn_cfg), flush=True)
        cl_cfg = calibrate_classic(cells)
        print("CLASSIC_CFG =", json.dumps(cl_cfg), flush=True)
    else:
        cells = _cell_windows(EV_REP, N_REP)
        pn_cfg, cl_cfg = PICK_CFG, CLASSIC_CFG
        print("PICK_CFG =", json.dumps(pn_cfg))
        print("CLASSIC_CFG =", json.dumps(cl_cfg))
    res = report(cells, pn_cfg, cl_cfg)
    res["todas las celdas"] = _entry([w for _c, wl in cells for w in wl])
    _print(res)
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(res, fh, indent=1)
    print(f"\ntotal {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
