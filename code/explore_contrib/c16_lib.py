"""
c16_lib.py - C16: semi-synthetic injection benchmark.  Synthetic events of the project generator (gen_realistic /
c8_lib.generate_site, same source, mechanism, pulse, coda, statics and site draws) placed in the local frame of the
REAL Alto Tiberina stations, with the generator's noise and clutter replaced by REAL continuous noise (c16_data.py).
SEED 999 ONLY.

Arms (same networks, same signals, paired):
  REAL  signal + real noise segment (no clutter, no Gaussian noise).
  GEN   signal + the generator's thin-regime clutter + Gaussian noise drawn from the same rng exactly as
        c8_lib.generate_site does (bit-identical to the standard pipeline; c16_validate.py).
  In REAL the rng draws of the clutter and noise are still consumed (and discarded), so that the sources of every later
  window of the stream and of every catalogue candidate are the same in both arms.

Networks (rules fixed before any detection run): cells of 4 / 6 / 8 stations of the 14 real stations, network index
j (2400-2413), subset rng [999, TAG_SUB, ci, j]:
  ladder_4, ladder_6   uniform random subset without replacement.
  unilateral_c40       one-sided: direction theta ~ U(0, 2 pi); the 8 stations whose azimuth (seen from the frame
                       origin = station centroid) is closest to theta.
  statics N(0, 0.3 s) per station and phase from [999, TAG_STAT, ci, j]; catalogue of 20 events with the registered
  gate and errors (xc.catalogue_cfg logic) from [999, TAG_CAT, ci, j, 0, 0|1|2|3].
Noise segments: the 2037 eligible 40 s slots of the day (c16_meta.json), split into 14 consecutive blocks
(np.array_split of the sorted list); network j uses block j - 2400 in every cell (cells share hours; networks of a cell
do not).  Inside a block, slot position p: p % 4 == 0 catalogue candidates, 1 STD windows, 2 DENSE windows, 3 noise-only
windows, each list in time order; if the catalogue list runs out it continues with the noise-only list from its end.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
for p in (HERE, WORK):
    if p not in sys.path:
        sys.path.insert(0, p)

import c2_run as C2  # noqa: E402
import c8_lib as C8  # noqa: E402
import c10_lib as C  # noqa: E402
import c10_run as R  # noqa: E402
import c11_lib as L  # noqa: E402
import gen_realistic as GR  # noqa: E402
import xc  # noqa: E402
import task2a_confirm as T  # noqa: E402
from station_corrections import estimate_station_corrections, peak_delays  # noqa: E402
from task1y_clutter import CLUT_AMP  # noqa: E402
from task1o_consistent import _KER  # noqa: E402

SD_GEN = GR.SNR * 0.02 * float(np.sqrt(np.sum(_KER ** 2)))   # in-band level of the generator's white noise

SEED = 999
TAG_STD, TAG_DENSE, TAG_CAT, TAG_STAT, TAG_SUB, TAG_NOISE0 = 5161, 5162, 5163, 5164, 5165, 5166
NETS_DEV = (2400, 2401)
NETS_EVAL = tuple(range(2402, 2414))
NBLOCK = 14
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
NSTA = {"ladder_4": 4, "ladder_6": 6, "unilateral_c40": 8}
REGIME = "thin"
RI = 0
NT, DT = C8.NT, C8.DT
NOISE_DIR = os.path.join(HERE, "out", "c16_noise")
METHODS = ["MT@S0", "MT@S0+T3", "MT@S2", "MT@S2+T3", "MT@ORA+T3", "BOMP@S0", "HUBER-N@S0"]

_NOISE = {}


def noise_store():
    if "d" not in _NOISE:
        _NOISE["d"] = np.load(os.path.join(NOISE_DIR, "c16_noise.npy"), mmap_mode="r")
        _NOISE["m"] = json.load(open(os.path.join(NOISE_DIR, "c16_meta.json")))
    return _NOISE["d"], _NOISE["m"]


# ============================================================================================ geometry and slots
def subset(cell, j):
    _d, m = noise_store()
    x, y = np.array(m["x_km"]), np.array(m["y_km"])
    ci = CELLS.index(cell)
    rng = np.random.default_rng([SEED, TAG_SUB, ci, j])
    n = NSTA[cell]
    if cell == "unilateral_c40":
        th = rng.uniform(0, 2 * np.pi)
        az = np.arctan2(y, x)
        dd = np.abs(np.angle(np.exp(1j * (az - th))))
        sub = np.sort(np.argsort(dd, kind="stable")[:n])
    else:
        sub = np.sort(rng.choice(len(x), n, replace=False))
    STA = np.stack([x[sub], y[sub], np.zeros(n)], axis=1)
    return sub, STA


def slots(j):
    _d, m = noise_store()
    el = np.array(m["eligible"])
    blk = np.array_split(el, NBLOCK)[j - NETS_DEV[0]]
    pos = np.arange(len(blk))
    return {"cat": blk[pos % 4 == 0].tolist(), "STD": blk[pos % 4 == 1].tolist(), "DENSE": blk[pos % 4 == 2].tolist(),
            "NOISE": blk[pos % 4 == 3].tolist()}


def real_segment(sub, slot):
    d, _m = noise_store()
    return np.asarray(d[sub, :, slot * NT:(slot + 1) * NT], float)


# ============================================================================================ generator split
def gen_signal(STA, sources, taus, rng_m, rng_h, statics, cfg=None):
    """c8_lib.generate_site lines up to the end of the event loop (verbatim logic): signal + coda, no clutter, no
    noise.  Consumes rng_m and rng_h exactly as generate_site."""
    c = GR.make_cfg(cfg)
    ns = len(STA)
    D = np.zeros((ns, 3, NT))
    picks = [[] for _ in range(ns)]
    st = rng_h.normal(0.0, c["static_sd"], size=(ns, 2))
    if statics is not None:
        st = np.asarray(statics, float).reshape(ns, 2)
    dn = rng_h.normal(0.0, float(np.tan(np.radians(c["dir_deg"]))), size=(ns, 2, 3))
    site = np.exp(rng_h.normal(0.0, c["site_sd"], size=(ns, 2)))
    rx = GR._child(rng_h)
    sos = GR._coda_filter(c["coda_band"])
    TT_ = GR.T
    for src, tau in zip(sources, taus):
        src = np.asarray(src, float)
        d = np.linalg.norm(STA - src, axis=1)
        Mt = GR.rand_dc(rng_m)
        fc = float(rx.uniform(*c["fc_range"]))
        phi = float(np.radians(rx.uniform(-c["phase_max_deg"], c["phase_max_deg"])))
        ew = GR.wavelet_energy(fc)
        tt = {b: GR.travel_time_gradient(src, STA, b, c) for b in ("P", "S")}
        for i in range(ns):
            for bi, (b, A) in enumerate((("P", GR.AP), ("S", GR.AS))):
                ta = float(tt[b][i]) + tau + float(st[i, bi])
                picks[i].append((ta, b))
                vec = GR.radiation_gradient(Mt, src, STA[i], b, c)
                a0 = np.linalg.norm(vec)
                if a0 > 0:
                    w = vec / a0 + dn[i, bi]
                    vec = a0 * site[i, bi] * w / np.linalg.norm(w)
                av = GR.SNR * A * float(np.linalg.norm(vec)) / d[i]
                D[i] += GR.SNR * A * np.outer(vec, GR.source_wavelet(TT_ - ta, fc, phi)) / d[i]
                GR._add_coda(D[i], ta, av, ew, c, sos, rx)
    return D, picks, site


def add_gen_noise(D, STA, sources, n_clutter, rng, clut_amp, cfg=None):
    """c8_lib.generate_site clutter + Gaussian noise (verbatim logic), added IN PLACE to D in the same order."""
    c = GR.make_cfg(cfg)
    ns = len(STA)
    TT_ = GR.T
    ref = np.asarray(sources[0], float) if len(sources) else np.zeros(3)
    scale = float(np.median(np.linalg.norm(STA - ref, axis=1)))
    for _ in range(n_clutter):
        i = int(rng.integers(ns))
        tc = float(rng.uniform(2, TT_[-1] - 4))
        q = rng.normal(size=3)
        q /= np.linalg.norm(q)
        D[i] += GR.SNR * clut_amp * np.outer(q, GR.ricker(TT_ - tc)) / scale
        rng.random()
    D += rng.normal(scale=GR.SNR * c["noise"], size=D.shape)
    return D


def finish(D_sig, STA, sources, n_clutter, rng, clut_amp, arm, sub, slot):
    """total window of an arm; returns D and the noise part N (D = D_sig + N up to the order of the additions)."""
    if arm == "GEN":
        D = D_sig.copy()
        add_gen_noise(D, STA, sources, n_clutter, rng, clut_amp)
        return D, D - D_sig
    add_gen_noise(np.zeros_like(D_sig), STA, sources, n_clutter, rng, clut_amp)       # consume, discard
    N = real_segment(sub, slot)
    return D_sig + N, N


# ============================================================================================ network
class Net16:
    def __init__(self, cell, j, arm, kcat=T.K_CAT):
        self.cell, self.j, self.arm = cell, j, arm
        self.ci = CELLS.index(cell)
        self.sub, self.STA = subset(cell, j)
        self.ns = len(self.STA)
        self.nclut = T.net_setup(SEED, cell)[3]
        self.clut_amp = CLUT_AMP[REGIME]
        self.statics = np.random.default_rng([SEED, TAG_STAT, self.ci, j]).normal(0.0, T.STATIC_SD_2A, (self.ns, 2))
        self.sl = slots(j)
        self.cat, self.drawn, self.cat_truth, self.cat_slots = self.catalogue(kcat)
        self.corr = estimate_station_corrections(self.STA, self.cat)[0]

    def catalogue(self, K, gate=25.0, min_arr=3):
        """xc.catalogue_cfg (= task2a_confirm.catalogue_real) with the signal / noise split."""
        seeds = [SEED, TAG_CAT, self.ci, self.j, RI]
        rng = np.random.default_rng(list(seeds) + [0])
        rng_m = np.random.default_rng(list(seeds) + [1])
        rng_h = np.random.default_rng(list(seeds) + [2])
        rng_e = np.random.default_rng(list(seeds) + [3])
        pool = list(self.sl["cat"]) + list(self.sl["NOISE"])[::-1]
        out, truth, used, drawn = [], [], [], 0
        while len(out) < K:
            src = np.array([rng.uniform(-15, 15), rng.uniform(-15, 15), rng.uniform(5, 13)])
            tau = float(rng.uniform(5.5, 7.5))
            Ds, _p, _s = gen_signal(self.STA, [src], [tau], rng_m, rng_h, self.statics)
            slot = pool[drawn] if self.arm == "REAL" else -1
            drawn += 1
            D, _N = finish(Ds, self.STA, [src], self.nclut, rng, self.clut_amp, self.arm, self.sub, slot)
            _dl, en = peak_delays(self.STA, D, src, tau)
            if (en >= gate).sum() < min_arr:
                continue
            src_c = src + rng_e.normal(0, T.CAT_ERR[0], 3)
            src_c[2] = float(np.clip(src_c[2], 1.0, 20.0))
            out.append((D, src_c, tau + float(rng_e.normal(0, T.CAT_ERR[1]))))
            truth.append((src, tau))
            used.append(int(slot))
        return out, drawn, truth, used

    def setup(self):
        STA, corr = self.STA, self.corr
        self.xops = xc.Ops(STA, corr, R.GRID)
        self.dops = C8.DOps(self.xops)
        sig, _cnt = C.sigma_catalogue(STA, self.cat, corr)
        E = C.Eng(R.GARR, STA, corr, R.TAUS, C.H_COARSE, G=self.dops.Gs, Q=self.dops.Qst, bomp=True)
        assert (E.IDX == self.dops.IDX).all()
        EF = C.Eng(R.FINE, STA, corr, R.TAUS, C.H_FINE, bomp=True)
        self.E0, self.EF0 = E, EF
        self.eng = {"S0": (E, E.prior(sig), EF, EF.prior(sig))}
        self.sig = {"S0": sig.tolist()}
        self.terms = {"S0": corr.copy(), "ORA": self.statics.copy()}
        sig_o = C.sigma_oracle(STA, self.cat_truth)
        EO = C.Eng(R.GARR, STA, self.statics, R.TAUS, C.H_COARSE, G=E.G, Q=E.Q, bomp=True)
        EFO = C.Eng(R.FINE, STA, self.statics, R.TAUS, C.H_FINE, G=EF.G, Q=EF.Q, bomp=True)
        self.eng["ORA"] = (EO, EO.prior(sig_o), EFO, EFO.prior(sig_o))
        self.sig["ORA"] = sig_o.tolist()
        i2 = {}
        c2 = L.s2_terms(STA, self.cat, corr, None, L.P, info=i2)
        sig2, _c = C.sigma_catalogue(STA, self.cat, c2)
        E2, EF2 = L.clone_eng(E, c2), L.clone_eng(EF, c2)
        self.eng["S2"] = (E2, E2.prior(sig2), EF2, EF2.prior(sig2))
        self.sig["S2"] = sig2.tolist()
        self.terms["S2"] = c2
        self.s2info = {k: i2[k] for k in ("sig", "n_res", "n_ev")}

    def run(self, name, MFn, kmax):
        if name == "HUBER-N@S0":
            rec, locs = C8.mtfw_general(self.dops, MFn, [], kmax, stat="huber_n")
            return [{"x": np.asarray(l_[:3], float), "tau": float(l_[3]), "s": float(s), "s0": float(s)}
                    for (s, _h), l_ in zip(rec, locs)]
        fam, rest = name.split("@")
        t3 = rest.endswith("+T3")
        E, pr, EF, prF = self.eng[rest.replace("+T3", "")]
        kw = dict(mode="hard", pr=pr, W=4, EF=EF if t3 else None, prF=prF)
        if fam == "MT":
            return C.mt_loop(E, MFn, kmax, **kw)
        return C.bomp_loop(E, None, MFn, kmax, **kw)

    def stream(self, design):
        tag, draw = {"STD": (TAG_STD, xc.draw_sources), "DENSE": (TAG_DENSE, C2.draw_dense)}[design]
        b = [SEED, tag, self.ci, self.j, RI]
        rng, rng_m, rng_h = (np.random.default_rng(b + [41]), np.random.default_rng(b + [233]),
                             np.random.default_rng(b + [61]))
        k = 0
        while True:
            srcs, taus = draw(rng)
            Ds, picks, site = gen_signal(self.STA, srcs, taus, rng_m, rng_h, self.statics)
            slot = self.sl[design][k] if self.arm == "REAL" else -1
            k += 1
            D, N = finish(Ds, self.STA, srcs, self.nclut, rng, self.clut_amp, self.arm, self.sub, slot)
            MF, sd = xc.matched(D)
            yield {"D": D, "Ds": Ds, "N": N, "MFn": MF / sd, "sd": sd, "truth": [np.asarray(s, float) for s in srcs],
                   "taus": [float(t) for t in taus], "picks": picks, "slot": int(slot)}

    def noise_windows(self):
        rng = np.random.default_rng([SEED, TAG_NOISE0, self.ci, self.j, RI, 41])
        k = 0
        while True:
            Ds = np.zeros((self.ns, 3, NT))
            slot = self.sl["NOISE"][k] if self.arm == "REAL" else -1
            k += 1
            D, N = finish(Ds, self.STA, [], self.nclut, rng, self.clut_amp, self.arm, self.sub, slot)
            MF, sd = xc.matched(D)
            yield {"D": D, "MFn": MF / sd, "sd": sd, "truth": [], "taus": [], "slot": int(slot)}


# ============================================================================================ SNR statistics
def snr_stats(w, ns):
    """per true arrival: peak 3-component |MF(signal)| within +-4 samples of the arrival / window noise level sd (the
    level the detector normalises by); per station: robust in-band level of the noise part / sd_gen."""
    MFs, _ = xc.matched(w["Ds"]) if len(w["truth"]) else (None, None)
    out = []
    if MFs is not None:
        A = np.linalg.norm(MFs, axis=1)
        for i in range(ns):
            for (ta, _b) in w["picks"][i][:2 * len(w["truth"])]:
                k = int(round(ta / DT))
                lo, hi = max(0, k - 4), min(NT, k + 5)
                if hi > lo:
                    out.append(round(float(A[i, lo:hi].max() / w["sd"]), 3))
    MFn_, _ = xc.matched(w["N"])
    lev = (np.median(np.abs(MFn_), axis=(1, 2)) / 0.6745 / SD_GEN).tolist()
    x = MFn_.ravel()
    x = x / (np.median(np.abs(x)) / 0.6745)
    kurt = float(np.mean(x ** 4) / np.mean(x ** 2) ** 2 - 3.0)
    return out, [round(v, 5) for v in lev], round(kurt, 3), round(float(np.max(np.abs(x))), 2)


def pack(atoms):
    return [[float(a["x"][0]), float(a["x"][1]), float(a["x"][2]), float(a["tau"]), float(a["s"]), float(a["s0"]),
             float(a.get("sref", np.nan))] for a in atoms]
