"""
c19_lib.py - exploration C19 (seed 999 only, exploratory): the double-couple bias of the synthetic evidence.

Every synthetic task draws its sources with task1w_mech.rand_dc, a pure double couple, while every detector (BOMP-PK,
MT-FW with or without cap, HUBER-N) scores the full 6-dimensional moment-tensor subspace of the homogeneous
straight-ray operator task1x_mtglrt.radiation_operator.  C19 replaces the mechanism, and only the mechanism, by
families that are non-double-couple but still moment tensors (inside the subspace up to the medium mismatch) and by
families that no point moment tensor produces (outside it):

    DC        task1w_mech.rand_dc (the generator's own draw; the window is bitwise that of c8_lib.generate_site)
    CLVD      R diag(-1, -1, 2) R^T with a random rotation R and random sign
    ISO       +-I (explosion or implosion): P only, no S
    MIX       a random symmetric matrix with i.i.d. N(0, 1) entries (GOE): a generic full moment tensor with random
              isotropic, CLVD and DC parts
    FORCE     a single force f (random unit vector): P = n_r (n_s . f), S = transverse part of f (outside the subspace)
    DC2       two independent double couples at the same hypocentre, the second 0.3 s later, each with half the
              energy (a two-subevent rupture: no single point moment tensor at one origin time)
    SCRAMBLE  the DC's vector amplitudes with a uniformly random direction per station and phase (the radiation keeps
              its amplitude distribution but loses all moment-tensor consistency; a diagnostic bound, not physics)

Energy matching: every family is scaled so that the mean over the sphere of AP^2 |u_P|^2 + AS^2 |u_S|^2 (homogeneous
radiation, the generator's CP and CS normalisation) equals that of a double couple; the comparison is of radiation
patterns at equal radiated energy, not of source strength.  FORCE uses the moment-tensor wavelet (no extra time
derivative), so only its pattern differs.

Pairing: the mechanism consumes rng_m only (the generator's convention), and the extra draws of DC2 / SCRAMBLE come
from a dedicated generator seeded from rng_m; sources, origin times, statics, polarisation jitter, site factors,
wavelets, coda, clutter and noise are identical across families for the same stream.  The catalogue of each network
(and hence its +SC corrections) is always generated with DC sources, so the detectors are the same across families.

    python c19_lib.py        (self-test)
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
for p in (HERE, WORK):
    if p not in sys.path:
        sys.path.insert(0, p)

import gen_realistic as GR  # noqa: E402
import c8_lib as C  # noqa: E402
import task1x_mtglrt as X  # noqa: E402
from task1w_mech import CP, CS, rand_dc  # noqa: E402

FAMILIES = ("DC", "CLVD", "ISO", "MIX", "FORCE", "DC2", "SCRAMBLE")
INSIDE = ("DC", "CLVD", "ISO", "MIX")          # point moment tensors: in the 6-D subspace up to the medium mismatch
DC2_LAG = 0.3


# ============================================================================================ mechanisms
def _fib(n=4000):
    k = np.arange(n) + 0.5
    z = 1.0 - 2.0 * k / n
    r = np.sqrt(1.0 - z * z)
    a = np.pi * (1.0 + 5 ** 0.5) * k
    return np.stack([r * np.cos(a), r * np.sin(a), z], 1)


_N = _fib()


def energy_mt(M):
    p = np.einsum("ki,ij,kj->k", _N, M, _N)
    u = _N @ M.T
    s = u - p[:, None] * _N
    return GR.AP ** 2 * np.mean(p ** 2) / CP + GR.AS ** 2 * np.mean(np.sum(s * s, 1)) / CS


def energy_force(f):
    p = _N @ f
    s = f[None, :] - p[:, None] * _N
    return GR.AP ** 2 * np.mean(p ** 2) / CP + GR.AS ** 2 * np.mean(np.sum(s * s, 1)) / CS


E_DC = energy_mt(rand_dc(np.random.default_rng(0)))


def _rot(r):
    q, rr = np.linalg.qr(r.normal(size=(3, 3)))
    return q * np.sign(np.diag(rr))[None, :]


def draw_mech(family, rng_m):
    """('mt', M) or ('force', f), energy-matched to a double couple; DC consumes rng_m exactly as rand_dc."""
    if family in ("DC", "SCRAMBLE"):
        return "mt", rand_dc(rng_m)
    if family == "DC2":
        return "mt2", (rand_dc(rng_m) / np.sqrt(2.0), rand_dc(rng_m) / np.sqrt(2.0))
    if family == "CLVD":
        R = _rot(rng_m)
        M = R @ np.diag([-1.0, -1.0, 2.0]) @ R.T * (1.0 if rng_m.random() < 0.5 else -1.0)
    elif family == "ISO":
        M = np.eye(3) * (1.0 if rng_m.random() < 0.5 else -1.0)
    elif family == "MIX":
        A = rng_m.normal(size=(3, 3))
        M = (A + A.T) / np.sqrt(2.0)
    elif family == "FORCE":
        f = rng_m.normal(size=3)
        f /= np.linalg.norm(f)
        return "force", f * np.sqrt(E_DC / energy_force(f))
    else:
        raise ValueError(family)
    return "mt", M * np.sqrt(E_DC / energy_mt(M))


def radiation_force(f, src, sta, phase, cfg=None):
    """gen_realistic.radiation_gradient with M n_s replaced by the force f (P: n_r (n_s . f); S: transverse part of
    f at the source, transported to the station in SV / SH)."""
    v0, g = GR.phase_velocity(phase, cfg)
    ns_, nr_, eh = GR._ray(src, sta, v0, g)
    if phase == "P":
        return nr_ * float(ns_ @ f) / np.sqrt(CP)
    ut = f - float(f @ ns_) * ns_
    if eh is None:
        return ut / np.sqrt(CS)
    esh = np.cross(GR.EZ, eh)
    return (float(ut @ np.cross(esh, ns_)) * np.cross(esh, nr_) + float(ut @ esh) * esh) / np.sqrt(CS)


# ============================================================================================ generator
def generate_mech(STA, sources, taus, n_clutter, rng, rng_m, rng_h, clut_amp, family="DC", cfg=None, statics=None,
                  meta=None):
    """c8_lib.generate_site (= gen_realistic.generate_real) with the mechanism of `family`.  Returns D, picks, site,
    amp.  meta (list) receives per event the noise-free arrival vectors (2, ns, 3) at the true times, before and after
    the polarisation jitter and site factors, for the subspace diagnostic."""
    c = GR.make_cfg(cfg)
    ns = len(STA)
    NT = GR.NT
    D = np.zeros((ns, 3, NT))
    picks = [[] for _ in range(ns)]
    amp = [[] for _ in range(ns)]
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
        kind, mech = draw_mech(family, rng_m)
        rs = np.random.default_rng(rng_m.integers(2 ** 63)) if family == "SCRAMBLE" else None
        fc = float(rx.uniform(*c["fc_range"]))
        phi = float(np.radians(rx.uniform(-c["phase_max_deg"], c["phase_max_deg"])))
        ew = GR.wavelet_energy(fc)
        tt = {b: GR.travel_time_gradient(src, STA, b, c) for b in ("P", "S")}
        clean, pert = np.zeros((2, ns, 3)), np.zeros((2, ns, 3))
        for i in range(ns):
            for bi, (b, A) in enumerate((("P", GR.AP), ("S", GR.AS))):
                ta = float(tt[b][i]) + tau + float(st[i, bi])
                picks[i].append((ta, b))
                if kind == "force":
                    vecs = [radiation_force(mech, src, STA[i], b, c)]
                elif kind == "mt2":
                    vecs = [GR.radiation_gradient(m_, src, STA[i], b, c) for m_ in mech]
                else:
                    vecs = [GR.radiation_gradient(mech, src, STA[i], b, c)]
                if family == "SCRAMBLE":
                    q = rs.normal(size=3)
                    vecs = [np.linalg.norm(vecs[0]) * q / np.linalg.norm(q)]
                clean[bi, i] = sum(vecs) * A / d[i]
                out = []
                for vec in vecs:
                    a0 = np.linalg.norm(vec)
                    if a0 > 0:
                        w = vec / a0 + dn[i, bi]
                        vec = a0 * site[i, bi] * w / np.linalg.norm(w)
                    out.append(vec)
                pert[bi, i] = sum(out) * A / d[i]
                av = GR.SNR * A * float(np.linalg.norm(np.concatenate(out))) / d[i]
                amp[i].append(av)
                for k_, vec in enumerate(out):
                    D[i] += GR.SNR * A * np.outer(vec, GR.source_wavelet(TT_ - (ta + k_ * DC2_LAG), fc, phi)) / d[i]
                GR._add_coda(D[i], ta, av, ew, c, sos, rx)
        if meta is not None:
            meta.append({"clean": clean, "pert": pert})
    ref = np.asarray(sources[0], float) if len(sources) else np.zeros(3)
    scale = float(np.median(np.linalg.norm(STA - ref, axis=1)))
    for _ in range(n_clutter):
        i = int(rng.integers(ns))
        tc = float(rng.uniform(2, TT_[-1] - 4))
        q = rng.normal(size=3)
        q /= np.linalg.norm(q)
        D[i] += GR.SNR * clut_amp * np.outer(q, GR.ricker(TT_ - tc)) / scale
        picks[i].append((tc, "P" if rng.random() < 0.5 else "S"))
        amp[i].append(GR.SNR * clut_amp / scale)
    D += rng.normal(scale=GR.SNR * c["noise"], size=D.shape)
    for _ in range(int(c["spikes"])):
        i, k, j = int(rx.integers(ns)), int(rx.integers(3)), int(rx.integers(NT))
        D[i, k, j] += (1.0 if rx.random() < 0.5 else -1.0) * c["spike_amp"] * GR.SNR * c["noise"]
    return D, picks, site, amp


def subspace_fraction(STA, src, v):
    """fraction of the energy of the arrival vectors v (2, ns, 3) inside the detector's 6-D subspace at the true
    source (homogeneous straight-ray operator, task1x_mtglrt)."""
    Q = X.basis(X.radiation_operator(np.asarray(src, float), STA))
    x = np.asarray(v, float).reshape(-1)
    e = float(x @ x)
    return float(np.sum((Q.T @ x) ** 2) / e) if e > 0 else np.nan


class MNet(C.WNet):
    """c8_lib.WNet (iid site factors, sd 0.5, DC catalogue and +SC corrections) with windows of a chosen family."""

    def stream(self, tag, draw, family):
        b = [C.SEED, tag, self.ci, self.j, self.ri]
        rng, rng_m, rng_h = (np.random.default_rng(b + [41]), np.random.default_rng(b + [233]),
                             np.random.default_rng(b + [61]))
        while True:
            srcs, taus = draw(rng)
            meta = []
            D, _picks, site, _amp = generate_mech(self.STA, srcs, taus, self.nclut, rng, rng_m, rng_h, self.clut_amp,
                                                  family=family, cfg=self.world.cfg, statics=self.statics, meta=meta)
            MF, sd = C.xc.matched(D)
            yield {"D": D, "MFn": MF / sd, "truth": [np.asarray(s, float) for s in srcs],
                   "taus": [float(t) for t in taus], "meta": meta, "site": site}


# ============================================================================================ self-test
def _selftest():
    ok = True

    def check(name, val, tol):
        nonlocal ok
        ok &= bool(val < tol)
        print(f"   {name:<66} {val:.3e}  (< {tol:g})  {'ok' if val < tol else 'FAIL'}")

    # 1. DC family reproduces c8_lib.generate_site bitwise, with identical rng states afterwards
    import c2_run as C2
    import task2a_confirm as T
    STA, statics = T.network(999, "ladder_6", 400)
    md = 0.0
    for s in range(4):
        g1 = [np.random.default_rng([s, k]) for k in range(3)]
        g2 = [np.random.default_rng([s, k]) for k in range(3)]
        srcs, taus = C2.draw_dense(g1[0])
        srcs2, taus2 = C2.draw_dense(g2[0])
        D1, _p, s1, _a = C.generate_site(STA, srcs, taus, 6, *g1, 1.75, statics=statics)
        D2, _p, s2, _a = generate_mech(STA, srcs2, taus2, 6, *g2, 1.75, family="DC", statics=statics)
        md = max(md, float(np.max(np.abs(D1 - D2))), float(np.max(np.abs(s1 - s2))),
                 float(abs(g1[0].random() - g2[0].random())), float(abs(g1[1].random() - g2[1].random())))
    check("DC family vs c8_lib.generate_site (max abs diff, rng states)", md, 1e-300)
    # 2. energy matching
    r = np.random.default_rng(11)
    em = 0.0
    for fam in ("CLVD", "ISO", "MIX"):
        for _ in range(20):
            em = max(em, abs(energy_mt(draw_mech(fam, r)[1]) / E_DC - 1))
    for _ in range(20):
        em = max(em, abs(energy_force(draw_mech("FORCE", r)[1]) / E_DC - 1))
    check("energy of CLVD / ISO / MIX / FORCE relative to a DC", em, 1e-12)
    ed = max(abs(energy_mt(rand_dc(r)) / E_DC - 1) for _ in range(50))
    check("energy spread among double couples (Fibonacci quadrature)", ed, 2e-3)
    # 3. ISO radiates no S, force radiation with g = 0 is the analytic one
    c0 = GR.make_cfg({"g": 0.0})
    src, sta = np.array([1.0, 2.0, 8.0]), np.array([20.0, -7.0, 0.0])
    check("ISO: |S| at a station", float(np.linalg.norm(GR.radiation_gradient(np.eye(3), src, sta, "S"))), 1e-12)
    f = np.array([0.3, -0.5, 0.8])
    n = (sta - src) / np.linalg.norm(sta - src)
    ef = max(float(np.linalg.norm(radiation_force(f, src, sta, "P", c0) - n * (n @ f) / np.sqrt(CP))),
             float(np.linalg.norm(radiation_force(f, src, sta, "S", c0) - (f - (f @ n) * n) / np.sqrt(CS))))
    check("FORCE, g = 0: radiation equal to the analytic one", ef, 1e-12)
    # 4. subspace fractions of clean homogeneous radiation: MT families 1, force and scramble below 1
    from task0_common import VS
    c_h = {"g": 0.0, "VS0": VS, "dir_deg": 0.0, "site_sd": 0.0, "static_sd": 0.0}
    rows = {}
    for fam in FAMILIES:
        fr = []
        for s in range(10):
            g = [np.random.default_rng([s, k, 7]) for k in range(3)]
            srcs = [np.array([g[0].uniform(-10, 10), g[0].uniform(-10, 10), g[0].uniform(5, 13)])]
            meta = []
            generate_mech(STA, srcs, [6.0], 0, *g, 1.75, family=fam, cfg=c_h, meta=meta)
            fr.append(subspace_fraction(STA, srcs[0], meta[0]["clean"]))
        rows[fam] = float(np.mean(fr))
    print("   clean homogeneous radiation, mean fraction in the 6-D subspace:",
          ", ".join(f"{k} {v:.3f}" for k, v in rows.items()))
    check("point moment tensors (DC, CLVD, ISO, MIX): 1 - fraction", max(1 - rows[k] for k in INSIDE), 1e-9)
    check("FORCE is outside the subspace (fraction < 0.99)", float(rows["FORCE"] > 0.99), 0.5)
    print("self-test", "passed" if ok else "FAILED")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if _selftest() else 1)
