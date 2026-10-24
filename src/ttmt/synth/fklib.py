"""Synthetic event library from the full-wave engine (``ttmt.synth.fk``): events with known moment tensors at the real
station coordinates of a network, for injection into that network's real noise (truth design 2 of ATLAS-0 with real
noise; resolves CR-14: the tensors, hence K_mec/E, are known while noise and station geometry are real).

Each event: a position (inside the network or out to 1.5x its aperture), a mechanism family (DC, DC+CLVD, DC+ISO,
generic, double source, extended source), a magnitude and a Brune source-time function (no Ricker), synthesised in a
true medium that is none of the detector's models (``media.PRESETS``, perturbed per network seed). Effects at the
receiver (CR-09): borehole depth (free-surface ghost, from the engine), shear-wave splitting at a fraction of the
stations, horizontal misorientation at 10-20 per cent of the stations. Arrival times are ray-theory first arrivals in
the true medium (``media.first_arrival``, the engine's own onsets; the validator checks them against the waveforms).
"""
import hashlib
import json
from dataclasses import dataclass, field

import numpy as np

from . import fk
from . import source as SRC
from .media import PRESETS, first_arrival

DEFAULT = {
    "medium": "crust_a", "medium_seed": 1, "v_sd": 0.03, "h_sd": 0.10, "vpvs_sd": 0.02,
    "fs": 100.0, "fmax": 25.0, "f_gauss": 15.0, "gf_len": 30.0, "L_factor": 1.5, "zr": 0.0,
    "depth_range": (2.0, 14.0), "depth_step": 0.5, "out_frac": 0.15, "out_factor": 1.5,
    "families": {"dc": 0.4, "clvd": 0.1, "iso": 0.05, "generic": 0.1, "double": 0.2, "extended": 0.15},
    "mw_range": (0.5, 2.5), "stress_mpa": (1.0, 10.0),
    "split_frac": 0.3, "split_dt": (0.05, 0.2), "misorient_frac": 0.15, "misorient_deg": (10.0, 45.0),
}


def config(cfg=None):
    c = json.loads(json.dumps(DEFAULT))
    c.update(cfg or {})
    return c


def config_hash(c):
    return hashlib.sha256(json.dumps(c, sort_keys=True, default=list).encode()).hexdigest()[:16]


@dataclass
class SynthEvent:
    id: str
    xyz: np.ndarray
    mag: float
    family: str
    subs: list                      # [(dt, tensor 3x3, moment fraction)]
    fc: float
    duration: float
    wave: np.ndarray = None         # (ns, 3, n) velocity, components (E, N, Z up), origin at sample pre * fs
    arr: np.ndarray = None          # (ns, 2) arrival times relative to the origin (s)
    pre: float = 5.0
    meta: dict = field(default_factory=dict)

    @property
    def M(self):
        """total tensor (moment-weighted), for K_mec bookkeeping."""
        return sum(w * T for _dt, T, w in self.subs)


def network_effects(ns, rng, c):
    """per-station receiver effects fixed for a network realisation: splitting (flag, fast azimuth, delay) and
    misorientation angle (rad, 0 for correctly oriented stations)."""
    split = rng.random(ns) < c["split_frac"]
    return {"split": split, "split_az": rng.uniform(0, np.pi, ns),
            "split_dt": np.where(split, rng.uniform(*c["split_dt"], ns), 0.0),
            "misorient": np.where(rng.random(ns) < c["misorient_frac"],
                                  np.radians(rng.uniform(*c["misorient_deg"], ns)) * rng.choice([-1, 1], ns), 0.0)}


def draw_events(STA, n, rng, c, prefix="fk"):
    """n events with positions, mechanisms, magnitudes and source-time functions (no waveforms yet)."""
    STA = np.asarray(STA, float)
    cen = STA[:, :2].mean(0)
    ap = float(np.max(np.linalg.norm(STA[:, :2] - cen, axis=1)))
    fam = list(c["families"])
    pw = np.array([c["families"][k] for k in fam], float)
    out = []
    zgrid = np.arange(c["depth_range"][0], c["depth_range"][1] + 1e-9, c["depth_step"])
    for j in range(n):
        if rng.random() < c["out_frac"]:
            rr = ap * rng.uniform(1.0, c["out_factor"])
        else:
            rr = ap * np.sqrt(rng.random())
        az = rng.uniform(0, 2 * np.pi)
        xy = cen + rr * np.array([np.cos(az), np.sin(az)])
        z = float(rng.choice(zgrid))
        f = fam[int(rng.choice(len(fam), p=pw / pw.sum()))]
        mw = float(rng.uniform(*c["mw_range"]))
        sd = float(np.exp(rng.uniform(*np.log(c["stress_mpa"]))))
        fc = float(SRC.brune_fc(mw, sd))
        base = SRC.mechanism(rng, "dc" if f in ("double", "extended") else f)
        if f == "double":
            subs = [(dt, T, w) for dt, _dx, T, w in SRC.subevents(rng, base, mw, "double")]
        else:
            subs = [(0.0, base, 1.0)]
        dur = float(SRC.rupture_duration(mw, sd)) if f == "extended" else 0.0
        out.append(SynthEvent(f"{prefix}{j:05d}", np.array([xy[0], xy[1], z]), mw, f, subs, fc, dur))
    return out


def true_medium(c):
    m = PRESETS[c["medium"]]
    if c["medium_seed"] is None:
        return m
    return m.perturbed(np.random.default_rng(c["medium_seed"]), c["v_sd"], c["h_sd"], c["vpvs_sd"])


def _apply_receiver(u, fs, t_s, eff, i):
    """splitting of the horizontals from the S arrival on (cosine ramp of 0.1 s), then misorientation (rotation of
    the horizontals by the station's angle). u (3, n) components (E, N, Z up), t_s sample time of S."""
    if eff["split"][i] and np.isfinite(t_s):
        a = eff["split_az"][i]
        ca, sa = np.cos(a), np.sin(a)
        fast = ca * u[0] + sa * u[1]
        slow = -sa * u[0] + ca * u[1]
        n = u.shape[-1]
        f = np.fft.rfftfreq(2 * n, 1 / fs)
        sl = np.fft.irfft(np.fft.rfft(slow, 2 * n) * np.exp(-2j * np.pi * f * eff["split_dt"][i]), 2 * n)[:n]
        t = np.arange(n) / fs
        w = np.clip((t - (t_s - 0.1)) / 0.1, 0.0, 1.0)
        w = 0.5 - 0.5 * np.cos(np.pi * w)
        slow = (1 - w) * slow + w * sl
        u = u.copy()
        u[0], u[1] = ca * fast - sa * slow, sa * fast + ca * slow
    th = eff["misorient"][i]
    if th:
        c, s = np.cos(th), np.sin(th)
        u = u.copy()
        u[0], u[1] = c * u[0] + s * u[1], -s * u[0] + c * u[1]
    return u


def synthesise(events, STA, c, eff, pre=5.0, progress=None, workers=1):
    """fill wave and arr of every event: one engine call per source depth with every (event, station) receiver of
    that depth; then the tensor combination, the source-time function, subevent delays and receiver effects."""
    STA = np.asarray(STA, float)
    ns = len(STA)
    med = true_medium(c)
    nt = int(round(c["gf_len"] * c["fs"]))
    cfg = fk.FKConfig(fs=c["fs"], nt=nt, fmax=c["fmax"], f_gauss=c["f_gauss"], t_shift=1.0, L_factor=c["L_factor"])
    depths = sorted({float(e.xyz[2]) for e in events})
    jobs = []
    for z in depths:
        ev = [e for e in events if float(e.xyz[2]) == z]
        r, phi = [], []
        for e in ev:
            d = STA[:, :2] - e.xyz[:2]
            r.extend(np.hypot(d[:, 0], d[:, 1]).tolist())
            phi.extend(np.arctan2(d[:, 1], d[:, 0]).tolist())
        jobs.append((z, ev, np.maximum(np.array(r), 0.05), np.array(phi)))

    def run(job):
        z, ev, r, phi = job
        return fk.greens(med, z, r, phi, cfg, zr=c["zr"])

    if workers > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(workers) as ex:
            res = list(ex.map(_greens_job, [(med.to_dict(), j[0], j[2], j[3], cfg.as_dict(), c["zr"]) for j in jobs]))
    else:
        res = []
        for jb in jobs:
            res.append(run(jb))
            if progress:
                progress(len(res), len(jobs))
    n_out = int(round((c["gf_len"] + pre) * c["fs"]))
    f = np.fft.rfftfreq(2 * n_out, 1 / c["fs"])
    for (z, ev, _r, _phi), G in zip(jobs, res):
        G = G.reshape(len(ev), ns, 6, 3, nt)
        for k, e in enumerate(ev):
            X = np.hypot(STA[:, 0] - e.xyz[0], STA[:, 1] - e.xyz[1])
            tp = first_arrival(med, z, X, "P")[0]
            ts = first_arrival(med, z, X, "S")[0]
            e.arr = np.stack([tp, ts], axis=1)
            stf = SRC.stf_spectrum(f, "brune_box" if e.duration > 0 else "brune", e.fc, e.duration or None)
            w = np.zeros((ns, 3, n_out))
            for dt, T, wt in e.subs:
                mv = SRC.as_vector(T) * wt
                u = np.einsum("j,sjcn->scn", mv, G[k])                      # (ns, 3, nt) for Gaussian rate at 1 s
                U = np.fft.rfft(u, 2 * n_out, axis=-1)
                # remove the engine's 1 s delay, add pre and the subevent delay; Gaussian anti-alias kept
                U *= stf * np.exp(-2j * np.pi * f * (dt + pre - 1.0))
                w += np.fft.irfft(U, 2 * n_out, axis=-1)[..., :n_out]
            w[:, 2] *= -1.0                                                  # z down -> Z up
            for i in range(ns):
                w[i] = _apply_receiver(w[i], c["fs"], pre + ts[i], eff, i)
            e.wave, e.pre = w, pre
            e.meta = {"medium": med.name, "medium_digest": med.digest()}
    return events


def _greens_job(a):
    from .media import LayeredModel
    md, z, r, phi, cd, zr = a
    cfg = fk.FKConfig()
    cfg.__dict__.update(cd)
    return fk.greens(LayeredModel.from_dict(md), z, r, phi, cfg, zr=zr)
