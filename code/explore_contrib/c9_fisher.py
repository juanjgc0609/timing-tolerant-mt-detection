"""C9 P5: Fisher information / Cramer-Rao bound of the epicentre per channel and cell geometry.  SEED 999 ONLY.
Networks task2a_confirm.network(999, cell, j), j = 130-149 (geometry only); sources and mechanisms rng [999, 3903, ci, j]
(10 per network: x, y uniform +-15 km, z 5-13 km, random double couple).  Homogeneous straight-ray model (VP 5.19, VS 3.0).
Amplitude calibration (noise units per unit of |G m|): 3 x 30 default-generator windows, network 199, streams
[999, 3904, ci, 199, ev, .] (dir_deg 0, site_sd 0).
Parameters: x, y, z, origin time, mechanism (double couple: 3 rotation angles + log M0; or full moment tensor: 6).
Channels and their noise:
  T     P and S arrival times, sd sigma_t (origin time free)
  PP    P polarization unit vector, sd per component  s^2 = tan^2(sigma_pol) + 1 / A_P^2  (A = amplitude in noise sd)
  SP    S polarization unit vector, same form with A_S (depends on the mechanism: nuisance)
  A     log amplitudes of P and S, variance sigma_a^2 + 1 / A^2
  RATIO log(A_S / A_P) per station, variance 2 sigma_a^2 + 1/A_P^2 + 1/A_S^2 (no M0, no inter-station amplitudes)
Bayesian (van Trees) bound with the benchmark's prior: x, y uniform on +-15 km (var 75), z uniform on 5-13 (var 5.33),
origin uniform on 2 s; flat on the mechanism (ridge 1e-8).  Reported: median over sources of sqrt(C_xx + C_yy) in km.
With no data at all this is 12.25 km.
    python c9_fisher.py > out/c9_fisher.txt
"""
import sys

import numpy as np

import c9_lib as C
import c8_lib as C8
import task2a_confirm as T
import xc
from task0_common import VP, VS
from task1w_mech import rand_dc
from task1y_clutter import CLUT_AMP

CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
PRIOR_LOC = np.array([1 / 75.0, 1 / 75.0, 12.0 / 64.0])
PRIOR_TAU = 12.0 / 4.0


def m6(Mt):
    return np.array([Mt[0, 0], Mt[1, 1], Mt[2, 2], Mt[0, 1], Mt[0, 2], Mt[1, 2]])


def rot(w):
    th = np.linalg.norm(w)
    if th < 1e-15:
        return np.eye(3)
    k = w / th
    Kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(th) * Kx + (1 - np.cos(th)) * Kx @ Kx


def forward(STA, x, tau, m):
    """times (2, ns), displacement vectors P and S (2, ns, 3) in |G m| units."""
    G = C.operators_hom(x[None], STA)[0]
    ns = len(STA)
    u = (G @ m).reshape(2, ns, 3)
    d = np.linalg.norm(STA - x, axis=1)
    return np.stack([d / VP + tau, d / VS + tau]), u


def jac(STA, x, Mt, mech):
    """Jacobians of the channel observables w.r.t. theta = (x, y, z, tau, mechanism...)."""
    m0 = m6(Mt)
    npar = 4 + (4 if mech == "dc" else 6)

    def obs(th):
        if mech == "dc":
            R = rot(th[4:7])
            m = m6(R @ Mt @ R.T) * np.exp(th[7])
        else:
            m = m0 + th[4:10]
        t, u = forward(STA, x + th[:3], th[3], m)
        a = np.linalg.norm(u, axis=-1)
        return t, u / a[..., None], np.log(a)
    t0, u0, la0 = obs(np.zeros(npar))
    Jt, Ju, Ja = [], [], []
    for k in range(npar):
        h = 1e-4 * (np.abs(m0).max() if (mech == "mt" and k >= 4) else 1.0)
        e = np.zeros(npar)
        e[k] = h
        tp, up, lp = obs(e)
        tm, um, lm = obs(-e)
        Jt.append((tp - tm) / (2 * h))
        Ju.append((up - um) / (2 * h))
        Ja.append((lp - lm) / (2 * h))
    return np.stack(Jt, -1), np.stack(Ju, -1), np.stack(Ja, -1), np.exp(la0)   # (2,ns,P), (2,ns,3,P), (2,ns,P), (2,ns)


_JC = {}


def bound(STA, x, Mt, kappa, chan, sig_t=0.144, sig_pol=10.0, sig_a=0.5, mech="dc", wave=True, snr=1.0):
    key = (STA.tobytes(), x.tobytes(), Mt.tobytes(), mech)
    if key not in _JC:
        _JC[key] = jac(STA, x, Mt, mech)
    Jt, Ju, Ja, a = _JC[key]
    A = kappa * snr * a                                                     # amplitudes in noise sd
    wn = (1.0 / A ** 2) if wave else np.zeros_like(A)
    P = Jt.shape[-1]
    I = np.zeros((P, P))
    if "T" in chan:
        J = Jt.reshape(-1, P)
        I += J.T @ J / sig_t ** 2
    s2 = np.tan(np.radians(sig_pol)) ** 2 + wn                              # (2, ns)
    if "PP" in chan:
        for i in range(len(STA)):
            I += Ju[0, i].T @ Ju[0, i] / max(s2[0, i], 1e-12)
    if "SP" in chan:
        for i in range(len(STA)):
            I += Ju[1, i].T @ Ju[1, i] / max(s2[1, i], 1e-12)
    if "A" in chan:
        va = sig_a ** 2 + wn
        J = Ja.reshape(-1, P)
        I += J.T @ (J / np.maximum(va.reshape(-1, 1), 1e-12))
    if "RATIO" in chan:
        J = Ja[1] - Ja[0]
        vr = 2 * sig_a ** 2 + wn[0] + wn[1]
        I += J.T @ (J / np.maximum(vr[:, None], 1e-12))
    pr = np.concatenate([PRIOR_LOC, [PRIOR_TAU], np.full(P - 4, 1e-8)])
    # scale-safe inversion
    Itot = I + np.diag(pr)
    sc = 1.0 / np.sqrt(np.diag(Itot))
    Cv = sc[:, None] * np.linalg.pinv(sc[:, None] * Itot * sc[None, :], rcond=1e-13) * sc[None, :]
    return float(np.sqrt(Cv[0, 0] + Cv[1, 1]))


def calibrate():
    """noise units per |G m| unit: median over blocks of |v at the true time| / |G_hom(src) m| (dir 0, site 0)."""
    rat = []
    for cell in CELLS:
        ci, _k, nst, nclut = T.net_setup(C.SEED, cell)
        STA, statics = T.network(C.SEED, cell, 199)
        for ev in range(30):
            sd_ = [C.SEED, 3904, ci, 199, ev]
            rng, rng_m, rng_h = (np.random.default_rng(sd_ + [41]), np.random.default_rng(sd_ + [233]),
                                 np.random.default_rng(sd_ + [61]))
            Mt = rand_dc(np.random.default_rng(sd_ + [233]))
            src = np.array([rng.uniform(-15, 15), rng.uniform(-15, 15), rng.uniform(5, 13)])
            D, picks, _s, _a = C8.generate_site(STA, [src], [6.5], nclut, rng, rng_m, rng_h, CLUT_AMP["thin"],
                                                cfg={"dir_deg": 0.0, "site_sd": 0.0}, statics=statics)
            MF, sd = xc.matched(D)
            tarr = np.array([[picks[i][b][0] for i in range(nst)] for b in range(2)])
            v, _K = C.read_true(MF / sd, tarr)
            _t, u = forward(STA, src, 0.0, m6(Mt))
            au, av = np.linalg.norm(u, axis=-1), np.linalg.norm(v, axis=-1)
            ok = au > np.median(au)
            rat.extend((av[ok] / au[ok]).tolist())
    return float(np.median(rat)), float(np.percentile(rat, 25)), float(np.percentile(rat, 75))


def sample(cell, nsrc=10, nets=range(130, 150)):
    ci = CELLS.index(cell)
    out = []
    for j in nets:
        STA, _st = T.network(C.SEED, cell, j)
        rng = np.random.default_rng([C.SEED, 3903, ci, j])
        for _ in range(nsrc):
            x = np.array([rng.uniform(-15, 15), rng.uniform(-15, 15), rng.uniform(5, 13)])
            out.append((STA, x, rand_dc(rng)))
    return out


def med(S, kappa, chan, **kw):
    v = np.array([bound(STA, x, Mt, kappa, chan, **kw) for STA, x, Mt in S])
    return float(np.median(v)), float(np.mean(v <= 5.0 / np.sqrt(2.0)))


def mc_check(cell, kappa, n=300):
    """Monte Carlo of the maximum-likelihood grid estimators against the bound (no waveform noise): P polarization
    with 10 deg scatter and arrival times with sd 0.144 s; network 130; rms and median epicentral error."""
    from c9_standalone import FINE
    STA, _st = T.network(C.SEED, cell, 130)
    geo = C.Geo(FINE, STA)
    rng = np.random.default_rng([C.SEED, 3905, CELLS.index(cell)])
    eP, eT, bP, bT = [], [], [], []
    for _ in range(n):
        x = np.array([rng.uniform(-15, 15), rng.uniform(-15, 15), rng.uniform(5, 13)])
        Mt = rand_dc(rng)
        nn, d = C.rays_hom(x[None], STA)
        u = nn[0] + rng.normal(0, np.tan(np.radians(10.0)), size=nn[0].shape)
        v = np.stack([100.0 * u / np.linalg.norm(u, axis=1, keepdims=True), np.zeros_like(u)])
        k = int(np.argmax(C.st_ppol_angle(geo, v)))
        eP.append(np.linalg.norm(FINE[k, :2] - x[:2]))
        t = np.stack([d[0] / VP, d[0] / VS]) + 6.0 + rng.normal(0, 0.144, size=(2, len(STA)))
        r = t[None] - np.stack([geo.d / VP, geo.d / VS], axis=1)
        r = r - r.mean(axis=(1, 2), keepdims=True)
        k = int(np.argmin((r * r).sum((1, 2))))
        eT.append(np.linalg.norm(FINE[k, :2] - x[:2]))
        bP.append(bound(STA, x, Mt, kappa, ("PP",), wave=False))
        bT.append(bound(STA, x, Mt, kappa, ("T",)))
    f = lambda a: (float(np.sqrt(np.mean(np.square(a)))), float(np.median(a)))
    return f(eP), f(bP), f(eT), f(bT)


def main():
    kappa, k25, k75 = calibrate()
    print(f"amplitude calibration (network 199, 90 single-event windows, dir 0, site 0): |v| / |G m| median {kappa:.1f} "
          f"(IQR {k25:.1f}-{k75:.1f}) noise sd per unit; analytic value for a matched 2 Hz pulse and noise-only sd: 77")
    S = {c: sample(c) for c in CELLS}
    for c in CELLS:
        aP, aS = [], []
        for STA, x, Mt in S[c]:
            _t, u = forward(STA, x, 0.0, m6(Mt))
            a = kappa * np.linalg.norm(u, axis=-1)
            aP.extend(a[0]); aS.extend(a[1])
        print(f"  {c:15s} amplitude in noise sd (site factor 1): P median {np.median(aP):.1f} (P25 {np.percentile(aP, 25):.1f}, "
              f"P75 {np.percentile(aP, 75):.1f}), S median {np.median(aS):.1f} (P25 {np.percentile(aS, 25):.1f}, P75 "
              f"{np.percentile(aS, 75):.1f}); waveform-noise angular scatter 1/A: P {np.degrees(1 / np.median(aP)):.0f} deg, "
              f"S {np.degrees(1 / np.median(aS)):.0f} deg")
    TTF = ("PP", "SP", "A")
    rows = [
        ("no data (prior only)", (), {}),
        ("T  sigma_t 0.05 s", ("T",), dict(sig_t=0.05)),
        ("T  sigma_t 0.144 s (+-0.25 s tolerance)", ("T",), dict(sig_t=0.144)),
        ("T  sigma_t 0.15 s (measured medium residual)", ("T",), dict(sig_t=0.15)),
        ("T  sigma_t 0.30 s (uncorrected statics)", ("T",), dict(sig_t=0.30)),
        ("PP 0 deg  (waveform noise only)", ("PP",), dict(sig_pol=0.0)),
        ("PP 5 deg", ("PP",), dict(sig_pol=5.0)),
        ("PP 10 deg (nominal)", ("PP",), dict(sig_pol=10.0)),
        ("PP 20 deg", ("PP",), dict(sig_pol=20.0)),
        ("PP 10 deg, no waveform noise", ("PP",), dict(sig_pol=10.0, wave=False)),
        ("PP 5 deg, no waveform noise", ("PP",), dict(sig_pol=5.0, wave=False)),
        ("PP 2 deg, no waveform noise", ("PP",), dict(sig_pol=2.0, wave=False)),
        ("PP+SP 10 deg (DC)", ("PP", "SP"), dict(sig_pol=10.0)),
        ("PP+SP 0 deg (DC)", ("PP", "SP"), dict(sig_pol=0.0)),
        ("A  sd 0.5 (DC)", ("A",), dict(sig_a=0.5)),
        ("A  sd 0.25 (DC)", ("A",), dict(sig_a=0.25)),
        ("A  sd 0 (DC; waveform noise only)", ("A",), dict(sig_a=0.0)),
        ("A  sd 0.5 (full MT)", ("A",), dict(sig_a=0.5, mech="mt")),
        ("A  sd 0 (full MT)", ("A",), dict(sig_a=0.0, mech="mt")),
        ("RATIO S/P sd 0.5 (DC)", ("RATIO",), dict(sig_a=0.5)),
        ("RATIO S/P sd 0 (DC)", ("RATIO",), dict(sig_a=0.0)),
        ("ALL time-free 10 deg, 0.5 (DC)", TTF, dict()),
        ("ALL time-free 10 deg, 0.5 (full MT)", TTF, dict(mech="mt")),
        ("ALL time-free 0 deg, 0 (DC) = oracle world", TTF, dict(sig_pol=0.0, sig_a=0.0)),
        ("ALL time-free 0 deg, 0, SNR x3", TTF, dict(sig_pol=0.0, sig_a=0.0, snr=3.0)),
        ("ALL time-free 0 deg, 0, SNR x10", TTF, dict(sig_pol=0.0, sig_a=0.0, snr=10.0)),
        ("ALL time-free 2 deg, 0.1, SNR x10", TTF, dict(sig_pol=2.0, sig_a=0.1, snr=10.0)),
        ("JOINT T 0.144 + ALL time-free nominal (DC)", ("T",) + TTF, dict()),
        ("JOINT T 0.144 + ALL time-free oracle world", ("T",) + TTF, dict(sig_pol=0.0, sig_a=0.0)),
        ("JOINT T 0.30 + ALL time-free nominal (DC)", ("T",) + TTF, dict(sig_t=0.30)),
    ]
    print("\nBayesian Cramer-Rao bound on the epicentre: median over 200 sources of sqrt(C_xx + C_yy), km")
    print(f"  {'channel':48s} " + " ".join(f"{c:>15s}" for c in CELLS))
    tab = {}
    for name, chan, kw in rows:
        r = [med(S[c], kappa, chan, **kw)[0] for c in CELLS]
        tab[name] = r
        print(f"  {name:48s} " + " ".join(f"{v:15.2f}" for v in r))
    a, b = tab["T  sigma_t 0.144 s (+-0.25 s tolerance)"], tab["JOINT T 0.144 + ALL time-free nominal (DC)"]
    print("  joint / time-only ratio (nominal):               " + " ".join(f"{y / x:15.3f}" for x, y in zip(a, b)))
    b = tab["JOINT T 0.144 + ALL time-free oracle world"]
    print("  joint / time-only ratio (oracle world):          " + " ".join(f"{y / x:15.3f}" for x, y in zip(a, b)))
    a, b = tab["T  sigma_t 0.30 s (uncorrected statics)"], tab["JOINT T 0.30 + ALL time-free nominal (DC)"]
    print("  joint / time-only ratio (sigma_t 0.30, nominal): " + " ".join(f"{y / x:15.3f}" for x, y in zip(a, b)))
    print("\nBreak-even scan: polarization scatter (deg) x SNR multiplier at which PP+SP (DC) reaches the T 0.144 s bound "
          "(median, km)")
    for snr in (1.0, 3.0, 10.0, 30.0):
        print(f"  SNR x{snr:<4g} " + "  ".join(
            f"{sp:g} deg: " + "/".join(f"{med(S[c][::4], kappa, ('PP', 'SP'), sig_pol=sp, snr=snr)[0]:.1f}" for c in CELLS)
            for sp in (0.0, 1.0, 2.0, 5.0, 10.0)))
    print("\nMonte Carlo check of the bound (network 130, 300 sources; grid ML on the 1 km grid; no waveform noise):")
    for c in CELLS:
        eP, bP, eT, bT = mc_check(c, kappa)
        print(f"  {c:15s} PP 10 deg: ML rms {eP[0]:.2f} median {eP[1]:.2f} | bound rms {bP[0]:.2f} median {bP[1]:.2f}     "
              f"T 0.144 s: ML rms {eT[0]:.2f} median {eT[1]:.2f} | bound rms {bT[0]:.2f} median {bT[1]:.2f}")


if __name__ == "__main__":
    main()
