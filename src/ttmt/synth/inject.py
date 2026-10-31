"""Semi-synthetic records for ATLAS-0: events injected into real continuous noise, with their truth table.

A *realisation* is one noise block of a real record (300 s, from the noise-role blocks of ``realdata.block_roles``)
plus injected targets: isolated events and pairs at controlled origin separations, signal-to-noise ratios and energy
ratios. Sources are either recorded events (``realdata.EGF``; tensors unknown) or full-wave synthetic events with known
tensors (``fklib.SynthEvent``). The signal and the noise are kept apart, so the signal-to-noise ratio of the card is
exact and the contamination variants (CR-24) are recombinations.

Truth (``ttmt.atlas.schema.Truth``): per event id, origin, location, magnitude (catalogue magnitude of the EGF, or
the synthetic Mw, plus log10 of the injection scale in ``meta``), target flag (catalogued events of the noise block are
neutral targets=False), pair id, K_mec/E of the pair (known tensors only), arr_t (true arrivals) and arr_t_obs (with
the optional injected timing error), and the per-station SNR of the card:
    peak |Z| of the event's own injected signal in [tP, tP + 2 s] / std of the noise record's Z in
    [tP - 5 s, tP - 0.5 s], both band-passed 2-20 Hz.
Per target (``meta['targets']``): source id and kind, injection scale, the largest normalised correlation of the
injected waveform with any template or deflation atom of a given library (CR-06), the radiated vectors g_ib measured on
the noise-free injected waveform, and the empirical K_mec/E of the pair from those vectors. ``meta['max_cc']``: the
same correlation per event of the truth, shape (m,), nan for the neutral catalogued events (read by the atlas scope
CC < 0.5).
"""
import hashlib
import json

import numpy as np

from . import realdata as RD
from ..atlas.schema import Truth

SEP_BINS = ((0.0, 0.1), (0.1, 0.5), (0.5, 1.0), (1.0, 3.0))
DEFAULT = {
    "block_s": 300.0, "group_dt": 40.0, "first": 10.0, "p_single": 0.3, "p_self": 0.2,
    "dist_bins": ((0.0, 1.0), (1.0, 5.0), (5.0, 15.0)),
    "sep_bins": SEP_BINS, "snr_range": (3.0, 60.0), "energy_ratio": (1.0, 10.0, 100.0),
    "timing_static_sd": 0.0, "timing_event_sd": 0.0,             # per-station shifts of the injected waveform (s)
    "max_egf_scale": 1.0, "taper_s": 0.5,
}
CONTAMINATION = ((1, 3.0), (1, 10.0), (2, 3.0), (2, 10.0))     # (stations scaled, factor) of CR-24


def config(cfg=None):
    c = json.loads(json.dumps(DEFAULT))
    c.update(cfg or {})
    return c


def config_hash(c):
    return hashlib.sha256(json.dumps(c, sort_keys=True, default=list).encode()).hexdigest()[:16]


# ------------------------------------------------------------------------------------------------ measurements
def card_snr(sig, noise, fs, tP, t0=0.0):
    """per-station SNR of the card: sig, noise (ns, 3, n) from record time t0; tP (ns,) P arrival (nan -> nan)."""
    zs = RD.bandpass(sig[:, 2], fs)
    zn = RD.bandpass(noise[:, 2], fs)
    out = np.full(len(tP), np.nan)
    for i, t in enumerate(tP):
        if not np.isfinite(t):
            continue
        a, b = int(round((t - t0) * fs)), int(round((t + 2.0 - t0) * fs))
        c, d = int(round((t - 5.0 - t0) * fs)), int(round((t - 0.5 - t0) * fs))
        if c < 0 or b > zs.shape[-1]:
            continue
        out[i] = float(np.max(np.abs(zs[i, a:b])) / (np.std(zn[i, c:d]) + 1e-30))
    return out


def radiated_vectors(wave, fs, arr, t0=0.0):
    """g (ns, 2, 3): the band-passed (2-20 Hz) velocity vector at the peak of its norm in [t - 0.05, t + 0.5] s after
    each arrival (P window cut at the S arrival); nan where the arrival is absent."""
    y = RD.bandpass(wave, fs)
    ns = wave.shape[0]
    g = np.full((ns, 2, 3), np.nan)
    for i in range(ns):
        for b in range(2):
            t = arr[i, b]
            if not np.isfinite(t):
                continue
            end = t + 0.5
            if b == 0 and np.isfinite(arr[i, 1]):
                end = min(end, arr[i, 1] - 0.05)
            a, z = int(round((t - 0.05 - t0) * fs)), int(round((end - t0) * fs))
            if z <= a or a < 0 or z > y.shape[-1]:
                continue
            seg = y[i, :, a:z]
            k = int(np.argmax(np.sum(seg * seg, 0)))
            g[i, b] = seg[:, k]
    return g


def kmec_over_e(g1, g2):
    """K_mec / E of a pair (Theory A.6: K_mec = sum_ib |g1 x g2|^2 / |g1 + g2|^2; E the mean event energy
    (sum_ib |g1|^2 + sum_ib |g2|^2) / 2) over the arrivals present in both."""
    ok = np.isfinite(g1).all(-1) & np.isfinite(g2).all(-1)
    if not ok.any():
        return np.nan
    a, b = g1[ok], g2[ok]
    s = np.sum((a + b) ** 2, -1)
    k = np.sum(np.cross(a, b) ** 2, -1) / np.maximum(s, 1e-300)
    E = 0.5 * (np.sum(a * a) + np.sum(b * b))
    return float(np.sum(k) / E) if E > 0 else np.nan


def max_cc(wave, atoms, fs, band=(2.0, 20.0)):
    """largest normalised network cross-correlation (sum over stations and components, lag searched over the full
    overlap) between an injected waveform (ns, 3, n) and any atom of a list of (ns, 3, m) arrays (nan entries are
    skipped stations). Returns (cc, index of the atom)."""
    if not atoms:
        return np.nan, -1
    x = RD.bandpass(np.nan_to_num(wave), fs, band)
    best, arg = -np.inf, -1
    for j, A in enumerate(atoms):
        y = RD.bandpass(np.nan_to_num(A), fs, band)
        n = x.shape[-1] + y.shape[-1] - 1
        nf = 1 << int(np.ceil(np.log2(n)))
        X, Y = np.fft.rfft(x, nf, axis=-1), np.fft.rfft(y, nf, axis=-1)
        cc = np.fft.irfft(X * np.conj(Y), nf, axis=-1).sum(axis=(0, 1))
        nx = np.sqrt(np.sum(x * x))
        ny = np.sqrt(np.sum(y * y))
        v = float(np.max(np.abs(cc)) / (nx * ny + 1e-300))
        if v > best:
            best, arg = v, j
    return best, arg


# ------------------------------------------------------------------------------------------------ sources
class Source:
    """uniform view of an injectable event: wave (ns, 3, n) at fs with the origin at sample pre * fs; arr (ns, 2)
    relative to the origin; xyz, mag, id, kind ('egf' or 'fk'), tensor (3x3 or None)."""

    def __init__(self, id, kind, wave, fs, pre, arr, xyz, mag, tensor=None, family=None):
        self.id, self.kind, self.wave, self.fs, self.pre = id, kind, wave, fs, pre
        self.arr, self.xyz, self.mag, self.tensor, self.family = arr, np.asarray(xyz, float), mag, tensor, family

    @classmethod
    def from_egf(cls, g, fs):
        return cls(g.id, "egf", g.wave, fs, g.t - g.t_cut, g.arr - g.t, g.xyz, g.mag)

    @classmethod
    def from_fk(cls, e, fs):
        return cls(e.id, "fk", e.wave, fs, e.pre, e.arr, e.xyz, e.mag, e.M, e.family)


def _taper(n, fs, ts):
    w = np.ones(n)
    k = int(ts * fs)
    if k > 0:
        r = 0.5 - 0.5 * np.cos(np.pi * np.arange(k) / k)
        w[:k], w[-k:] = r, r[::-1]
    return w


def _shift(x, fs, dt):
    """sub-sample delay of x (..., n) by dt seconds (frequency domain, zero-padded)."""
    n = x.shape[-1]
    f = np.fft.rfftfreq(2 * n, 1 / fs)
    return np.fft.irfft(np.fft.rfft(x, 2 * n, axis=-1) * np.exp(-2j * np.pi * f * dt), 2 * n, axis=-1)[..., :n]


# ------------------------------------------------------------------------------------------------ realisation
def plan(rng, c, n_sources, xyz=None):
    """origin times (record-relative), source indices, pair ids and energy ratios of the targets of one block."""
    groups = []
    t = c["first"]
    while t + c["group_dt"] <= c["block_s"] + 1e-9:
        groups.append(t)
        t += c["group_dt"]
    ev = []
    pid = 0
    for t0 in groups:
        t0 = t0 + rng.uniform(-2.0, 2.0)
        if rng.random() < c["p_single"] or n_sources < 2:
            ev.append({"t": t0, "src": int(rng.integers(n_sources)), "pair": -1, "rel": 1.0})
            continue
        a, b = rng.choice(n_sources, 2, replace=False)
        if xyz is not None:                     # partner from a spatial-separation bin of the spec (0-1, 1-5, 5-15 km)
            lo_d, hi_d = c["dist_bins"][int(rng.integers(len(c["dist_bins"])))]
            d = np.linalg.norm(xyz - xyz[a], axis=1)
            cand = np.flatnonzero((d >= lo_d) & (d < hi_d) & (np.arange(n_sources) != a))
            if len(cand):
                b = int(rng.choice(cand))
        if rng.random() < c["p_self"]:          # the same source twice: co-located, same mechanism (K_mec = 0)
            b = a
        lo, hi = c["sep_bins"][int(rng.integers(len(c["sep_bins"])))]
        dt = rng.uniform(lo, hi) * (1 if rng.random() < 0.5 else -1)
        er = float(rng.choice(c["energy_ratio"]))
        ev.append({"t": t0, "src": int(a), "pair": pid, "rel": 1.0})
        ev.append({"t": t0 + dt, "src": int(b), "pair": pid, "rel": 1.0 / np.sqrt(er)})
        pid += 1
    return ev


def realise(rec, block, sources, seed, c=None, atoms=None, atom_ids=None, label=""):
    """one realisation: (signal, noise, Truth). block: (start, end) record time of the noise block; sources: list of
    Source on the record's stations; seed: a numpy Generator (for ATLAS-0 the anchored stream
    ``ttmt.atlas.seeds.realisation_rng(set_id, i, r)``, ADR-0009 item 5) or an integer; atoms: optional template /
    deflation library for CR-06."""
    c = config(c)
    rng = seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)
    fs = rec.fs
    a, b = block
    noise = rec.cut(a, b).copy()
    ns, _, n = noise.shape
    sig = np.zeros_like(noise)
    ev = plan(rng, c, len(sources), np.array([s.xyz for s in sources]))
    for e in ev:
        e["t"] += a                                     # block-relative -> record time
    st_static = rng.normal(0.0, c["timing_static_sd"], ns) if c["timing_static_sd"] > 0 else np.zeros(ns)
    rows, tmeta = [], []
    for e in ev:
        s = sources[e["src"]]
        snr_t = float(np.exp(rng.uniform(*np.log(c["snr_range"]))))
        # scale for the target network SNR at unit scale, measured against this block's noise
        w = s.wave * _taper(s.wave.shape[-1], fs, c["taper_s"])
        arr_rel = s.arr
        tP = e["t"] + arr_rel[:, 0]
        t_cut = e["t"] - s.pre
        i0 = int(np.floor((t_cut - a) * fs))
        frac = (t_cut - a) * fs - i0
        w = _shift(w, fs, frac / fs)
        tim = rng.normal(0.0, c["timing_event_sd"], ns) if c["timing_event_sd"] > 0 else np.zeros(ns)
        unit = np.zeros_like(noise)
        j0, j1 = max(i0, 0), min(i0 + w.shape[-1], n)
        if j1 <= j0:
            continue
        for i in range(ns):
            wi = _shift(w[i], fs, tim[i] + st_static[i]) if (tim[i] or st_static[i]) else w[i]
            unit[i, :, j0:j1] = wi[:, j0 - i0:j1 - i0]
        snr1 = card_snr(unit, noise, fs, tP + tim + st_static, a)
        med = float(np.nanmedian(snr1)) if np.isfinite(snr1).any() else np.nan
        if not np.isfinite(med) or med <= 0:
            continue
        rows.append((e, s, unit, med, snr_t, tP, tim))
    # pair scales: the first of a pair sets the SNR, the second follows its energy ratio
    scale = {}
    for k, (e, s, unit, med, snr_t, tP, tim) in enumerate(rows):
        if e["pair"] >= 0 and any(r[0]["pair"] == e["pair"] for r in rows[:k]):
            j = next(jj for jj, r in enumerate(rows[:k]) if r[0]["pair"] == e["pair"])
            target = rows[j][4] * e["rel"]
        else:
            target = snr_t
        sc = target / med
        if s.kind == "egf":
            sc = min(sc, c["max_egf_scale"])
        scale[k] = sc
    ids, ts, xyz, mag, arr, arr_obs, snr, pair, gvec = [], [], [], [], [], [], [], [], []
    for k, (e, s, unit, med, snr_t, tP, tim) in enumerate(rows):
        sc = scale[k]
        sig += sc * unit
        ids.append(f"{label}#{k}")
        ts.append(e["t"])
        xyz.append(s.xyz)
        mag.append(s.mag + (np.log10(sc) if s.kind == "fk" else 0.0))
        at = e["t"] + s.arr
        arr.append(at)
        obs = at + (tim + st_static)[:, None]
        arr_obs.append(obs)
        snr.append(card_snr(sc * unit, noise, fs, obs[:, 0], a))
        pair.append(e["pair"])
        gvec.append(radiated_vectors(sc * unit, fs, obs, a))
        keep = [j for j in range(len(atoms or [])) if not atom_ids or atom_ids[j] != s.id]
        seg = sc * unit[:, :, max(0, int((np.nanmin(obs) - 1.0 - a) * fs)):int((np.nanmax(obs) + 6.0 - a) * fs)]
        cc, ia = max_cc(seg, [atoms[j] for j in keep], fs)
        ia = keep[ia] if ia >= 0 else -1
        part = [r_[1] for r_ in rows if r_[0]["pair"] == e["pair"] and r_[1] is not s] if e["pair"] >= 0 else []
        dxy = float(np.linalg.norm(part[0].xyz - s.xyz)) if part else (0.0 if e["pair"] >= 0 else None)
        tmeta.append({"source": s.id, "kind": s.kind, "family": s.family, "scale": float(sc), "pair_dist_km": dxy,
                      "log10_scale": float(np.log10(sc)), "max_cc": cc,
                      "max_cc_atom": (atom_ids[ia] if atom_ids and ia >= 0 else None),
                      "g": gvec[-1].tolist(), "tensor": None if s.tensor is None else np.asarray(s.tensor).tolist()})
    m = len(ts)
    pair = np.array(pair, int)
    kmec = np.full(m, np.nan)
    kmec_emp = np.full(m, np.nan)
    for p in set(pair.tolist()) - {-1}:
        jj = np.where(pair == p)[0]
        if len(jj) != 2:
            continue
        ke = kmec_over_e(gvec[jj[0]], gvec[jj[1]])
        kmec_emp[jj] = ke
        if all(rows[j][1].tensor is not None for j in jj):
            kmec[jj] = ke
    # neutral catalogued events of the noise block
    for e in rec.catalogue:
        if a <= e["t"] < b:
            ids.append(f"cat:{e['t']:.2f}")
            ts.append(e["t"])
            xyz.append([e["x"], e["y"], e["z"]])
            mag.append(e["M"])
            arr.append(np.full((ns, 2), np.nan))
            arr_obs.append(np.full((ns, 2), np.nan))
            snr.append(np.full(ns, np.nan))
    nn = len(ts) - m
    truth = Truth(np.array(ids), np.array(ts), np.array(xyz, float).reshape(-1, 3), np.array(mag, float),
                  np.array(arr, float).reshape(len(ts), ns, 2), np.array(snr, float).reshape(len(ts), ns), a, b,
                  target=np.r_[np.ones(m, bool), np.zeros(nn, bool)], pair=np.r_[pair, -np.ones(nn, int)],
                  kmec=np.r_[kmec, np.full(nn, np.nan)],
                  arr_t_obs=np.array(arr_obs, float).reshape(len(ts), ns, 2),
                  meta={"seed": None if isinstance(seed, np.random.Generator) else int(seed), "config_hash": config_hash(c), "config": c, "dataset": rec.dataset,
                        "block": [a, b], "targets": tmeta, "kmec_empirical": kmec_emp.tolist(),
                        "max_cc": [x["max_cc"] for x in tmeta] + [float("nan")] * nn,
                        "stations": rec.codes, "xyz": rec.xyz.tolist(), "fs": fs,
                        "components": "E, N, Z up", "label": label})
    return sig, noise, truth


def contaminate(noise, stations, factor):
    """noise with the given stations scaled by factor (CR-24)."""
    out = noise.copy()
    out[list(stations)] *= factor
    return out


def contamination_variants(subset, rng):
    """for a 4-station subset: the variants of CR-24, [(name, stations scaled, factor)], stations drawn at random."""
    out = []
    for k, f in CONTAMINATION:
        st = sorted(rng.choice(subset, k, replace=False).tolist())
        out.append((f"x{f:g}_{k}sta", st, f))
    return out


def truth_contaminated(truth, stations, factor):
    """the truth of a contamination variant: SNR of the scaled stations divided by the factor."""
    snr = truth.snr.copy()
    snr[:, list(stations)] /= factor
    t = Truth(truth.id, truth.t, truth.xyz, truth.mag, truth.arr_t, snr, truth.t0, truth.t1, truth.target, truth.pair,
              truth.kmec, truth.arr_t_obs, dict(truth.meta))
    t.meta["contamination"] = {"stations": list(map(int, stations)), "factor": float(factor)}
    return t
