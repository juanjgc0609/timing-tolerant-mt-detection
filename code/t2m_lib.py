"""
t2m_lib.py - Task 2M: semi-synthetic injection benchmark of C16 (explore_contrib/c16_lib.py) with the seed, the noise
day(s), the network indices, the window counts and the method list as parameters.  The generator, the signal / noise
split, the arms, the subset rules, the slot rule, the catalogue, the engines and the methods are those of c16_lib.py
(code copied; c16_*.py are not modified).  configure() must be called before anything else.

Parameters (configure):
  seed        root seed of every random stream ([seed, TAG, ci, j, ...], TAGs as C16: 5161 STD, 5162 DENSE,
              5163 catalogue, 5164 statics, 5165 subsets, 5166 noise-only windows of the GEN arm).
  noise_dirs  list of directories written by t2m_data.py (noise.npy, meta.json), one per noise day, in this order.
              A C16-format directory (c16_noise.npy, c16_meta.json) is also accepted.
  j0, nblocks network j uses noise block j - j0 of the global block list; the global list is, day after day, the
              eligible slots of the day in time order split into nblocks[d] consecutive blocks (np.array_split).
              With noise_dirs = [100723], j0 = 2400, nblocks = [14] this is exactly C16 (block j - 2400).
              apportion(N, counts) gives nblocks for N networks over several days (largest remainder, proportional
              to the number of eligible slots of each day).
  nwin        windows per network and design {"STD": nS, "DENSE": nD, "NOISE": nN}.
  methods     subset of ALL_METHODS (C16: all seven).

Networks (cells share noise hours: network j of every cell uses block j - j0, as in C16):
  The candidate stations of network j are the stations of the day of its block (meta["stations"]; the frame of that
  day: flat earth around 43.315 N, 12.388 E shifted by the centroid of the day's station set).  Subset rng
  [seed, 5165, ci, j]:
    ladder_4, ladder_6  uniform subset of 4 / 6 of the day's stations without replacement (rng.choice(n, k, False)).
    unilateral_c40      theta ~ U(0, 2 pi) (always drawn), and the 8 stations whose azimuth seen from the frame origin
                        is closest to theta (stable argsort of the wrapped angular distance).  With 8 stations this is
                        the WHOLE network for every theta; with 9 it drops the one station farthest in azimuth from
                        theta; with 11 it drops three; with 14 (C16) six.
Slots inside a block (as C16): position p, p % 4 == 0 catalogue candidates, 1 STD, 2 DENSE, 3 noise-only, each in
time order.  Catalogue candidate pool of the REAL arm: the catalogue slots, then the noise-only slots from the end of
the list (C16's rule), skipping the first nN noise-only slots (those used by the noise-only windows), then the unused
DENSE and STD slots from the end.  As long as the catalogue needs fewer than (#cat + #noise - nN) candidates this is
C16's sequence exactly; a pool exhausted raises an error (no silent reuse of a slot).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json
import sys

import numpy as np

WORK = os.path.dirname(os.path.abspath(__file__))
EXPL = os.path.join(WORK, "explore_contrib")
for p in (EXPL, WORK):
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

SD_GEN = GR.SNR * 0.02 * float(np.sqrt(np.sum(_KER ** 2)))

TAG_STD, TAG_DENSE, TAG_CAT, TAG_STAT, TAG_SUB, TAG_NOISE0 = 5161, 5162, 5163, 5164, 5165, 5166
CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
NSTA = {"ladder_4": 4, "ladder_6": 6, "unilateral_c40": 8}
REGIME = "thin"
RI = 0
NT, DT = C8.NT, C8.DT
ALL_METHODS = ["MT@S0", "MT@S0+T3", "MT@S2", "MT@S2+T3", "MT@ORA+T3", "BOMP@S0", "HUBER-N@S0"]
T2M_METHODS = ["MT@S0", "MT@S0+T3", "MT@S2+T3", "BOMP@S0", "HUBER-N@S0"]

CFG = {}
_DAYS = []
_BLOCKS = []


def apportion(n, counts):
    """largest-remainder apportionment of n blocks to days proportional to counts (ties -> earlier day)."""
    counts = np.asarray(counts, float)
    q = n * counts / counts.sum()
    k = np.floor(q).astype(int)
    rest = n - int(k.sum())
    order = sorted(range(len(q)), key=lambda i: (-(q[i] - k[i]), i))
    for i in order[:rest]:
        k[i] += 1
    return k.tolist()


def _load_day(d):
    if os.path.exists(os.path.join(d, "noise.npy")):
        a, m = "noise.npy", "meta.json"
    else:
        a, m = "c16_noise.npy", "c16_meta.json"
    return np.load(os.path.join(d, a), mmap_mode="r"), json.load(open(os.path.join(d, m)))


def configure(seed, noise_dirs, j0, nblocks, nwin, methods=None):
    CFG.clear()
    CFG.update({"seed": int(seed), "noise_dirs": list(noise_dirs), "j0": int(j0), "nblocks": [int(x) for x in nblocks],
                "nwin": dict(nwin), "methods": list(methods or ALL_METHODS)})
    assert all(m in ALL_METHODS for m in CFG["methods"])
    assert len(CFG["nblocks"]) == len(CFG["noise_dirs"])
    _DAYS.clear()
    _BLOCKS.clear()
    for di, d in enumerate(noise_dirs):
        arr, meta = _load_day(d)
        _DAYS.append((arr, meta))
        for blk in np.array_split(np.array(meta["eligible"]), CFG["nblocks"][di]):
            _BLOCKS.append((di, blk))
    CFG["n_networks"] = len(_BLOCKS)
    CFG["nets"] = list(range(CFG["j0"], CFG["j0"] + len(_BLOCKS)))


def configure_file(path):
    c = json.load(open(path))
    configure(c["seed"], c["noise_dirs"], c["j0"], c["nblocks"], c["nwin"], c.get("methods"))
    CFG["file"] = os.path.abspath(path)
    return c


def block(j):
    return _BLOCKS[j - CFG["j0"]]


# ============================================================================================ geometry and slots
def subset(cell, j):
    di, _b = block(j)
    m = _DAYS[di][1]
    x, y = np.array(m["x_km"]), np.array(m["y_km"])
    ci = CELLS.index(cell)
    rng = np.random.default_rng([CFG["seed"], TAG_SUB, ci, j])
    n = NSTA[cell]
    assert len(x) >= n, (cell, len(x))
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
    _di, blk = block(j)
    pos = np.arange(len(blk))
    return {"cat": blk[pos % 4 == 0].tolist(), "STD": blk[pos % 4 == 1].tolist(), "DENSE": blk[pos % 4 == 2].tolist(),
            "NOISE": blk[pos % 4 == 3].tolist()}


def cat_pool(sl):
    nw = CFG["nwin"]
    used_noise = set(sl["NOISE"][:nw["NOISE"]])
    pool = list(sl["cat"]) + [s for s in sl["NOISE"][::-1] if s not in used_noise]
    pool += list(sl["DENSE"][nw["DENSE"]:])[::-1] + list(sl["STD"][nw["STD"]:])[::-1]
    return pool


def real_segment(di, sub, slot):
    d = _DAYS[di][0]
    return np.asarray(d[sub, :, slot * NT:(slot + 1) * NT], float)


# ============================================================================================ generator split (C16)
def gen_signal(STA, sources, taus, rng_m, rng_h, statics, cfg=None):
    """c8_lib.generate_site lines up to the end of the event loop (verbatim logic, = c16_lib.gen_signal)."""
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
    """c8_lib.generate_site clutter + Gaussian noise (verbatim logic, = c16_lib.add_gen_noise), in place."""
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


def finish(D_sig, STA, sources, n_clutter, rng, clut_amp, arm, di, sub, slot):
    if arm == "GEN":
        D = D_sig.copy()
        add_gen_noise(D, STA, sources, n_clutter, rng, clut_amp)
        return D, D - D_sig
    add_gen_noise(np.zeros_like(D_sig), STA, sources, n_clutter, rng, clut_amp)       # consume, discard
    N = real_segment(di, sub, slot)
    return D_sig + N, N


# ============================================================================================ network
class Net2M:
    def __init__(self, cell, j, arm, kcat=T.K_CAT):
        self.seed = CFG["seed"]
        self.cell, self.j, self.arm = cell, j, arm
        self.ci = CELLS.index(cell)
        self.di = block(j)[0]
        self.day = _DAYS[self.di][1].get("day", "100723")
        self.sub, self.STA = subset(cell, j)
        self.stations = [_DAYS[self.di][1]["stations"][i] for i in self.sub]
        self.ns = len(self.STA)
        self.nclut = T.net_setup(self.seed, cell)[3]
        self.clut_amp = CLUT_AMP[REGIME]
        self.statics = np.random.default_rng([self.seed, TAG_STAT, self.ci, j]).normal(0.0, T.STATIC_SD_2A, (self.ns, 2))
        self.sl = slots(j)
        self.cat, self.drawn, self.cat_truth, self.cat_slots = self.catalogue(kcat)
        self.corr = estimate_station_corrections(self.STA, self.cat)[0]

    def catalogue(self, K, gate=25.0, min_arr=3):
        seeds = [self.seed, TAG_CAT, self.ci, self.j, RI]
        rng = np.random.default_rng(list(seeds) + [0])
        rng_m = np.random.default_rng(list(seeds) + [1])
        rng_h = np.random.default_rng(list(seeds) + [2])
        rng_e = np.random.default_rng(list(seeds) + [3])
        pool = cat_pool(self.sl)
        out, truth, used, drawn = [], [], [], 0
        while len(out) < K:
            src = np.array([rng.uniform(-15, 15), rng.uniform(-15, 15), rng.uniform(5, 13)])
            tau = float(rng.uniform(5.5, 7.5))
            Ds, _p, _s = gen_signal(self.STA, [src], [tau], rng_m, rng_h, self.statics)
            if self.arm == "REAL" and drawn >= len(pool):
                raise RuntimeError(f"catalogue slot pool exhausted ({self.cell}, {self.j}, {drawn} drawn)")
            slot = pool[drawn] if self.arm == "REAL" else -1
            drawn += 1
            D, _N = finish(Ds, self.STA, [src], self.nclut, rng, self.clut_amp, self.arm, self.di, self.sub, slot)
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
        b = [self.seed, tag, self.ci, self.j, RI]
        rng, rng_m, rng_h = (np.random.default_rng(b + [41]), np.random.default_rng(b + [233]),
                             np.random.default_rng(b + [61]))
        k = 0
        while True:
            srcs, taus = draw(rng)
            Ds, picks, site = gen_signal(self.STA, srcs, taus, rng_m, rng_h, self.statics)
            slot = self.sl[design][k] if self.arm == "REAL" else -1
            k += 1
            D, N = finish(Ds, self.STA, srcs, self.nclut, rng, self.clut_amp, self.arm, self.di, self.sub, slot)
            MF, sd = xc.matched(D)
            yield {"D": D, "Ds": Ds, "N": N, "MFn": MF / sd, "sd": sd, "truth": [np.asarray(s, float) for s in srcs],
                   "taus": [float(t) for t in taus], "picks": picks, "slot": int(slot)}

    def noise_windows(self):
        rng = np.random.default_rng([self.seed, TAG_NOISE0, self.ci, self.j, RI, 41])
        k = 0
        while True:
            Ds = np.zeros((self.ns, 3, NT))
            slot = self.sl["NOISE"][k] if self.arm == "REAL" else -1
            k += 1
            D, N = finish(Ds, self.STA, [], self.nclut, rng, self.clut_amp, self.arm, self.di, self.sub, slot)
            MF, sd = xc.matched(D)
            yield {"D": D, "MFn": MF / sd, "sd": sd, "truth": [], "taus": [], "slot": int(slot)}


# ============================================================================================ SNR statistics (C16)
def snr_stats(w, ns):
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
