"""Exploration: how much of the identifiability margin of a co-located pair survives when the one-event class
admits shear-wave splitting. Whitened norm (Appendix AB), 2 Hz Ricker, 8-station ring, straight rays, far field.
Classes for the S block: (1) one free pulse (current M6 null); (2) free orthogonal split, any delay (E.4a);
(3) orthogonal split with |delay| <= dmax; (4) CALIBRATED split: axis and delay known per station (no extra freedom):
the pair is two split pulses, the one-event class one split pulse with free time and free vector in the S plane."""
import numpy as np
fc = 2.0
def a(t):
    u = (np.pi*fc*t)**2
    return (1-2*u+u*u/3)*np.exp(-u/2)
EPS = np.linspace(-0.6, 1.6, 881)            # 2.5 ms grid for model pulse times
def R_vec(g1, g2, D):                        # Lemma A.4, vector amplitudes, one free pulse
    f = ((g1[:, None]*a(EPS)[None] + g2[:, None]*a(EPS-D)[None])**2).sum(0)
    return g1@g1 + g2@g2 + 2*g1@g2*a(D) - f.max()
def scal_curve(al, be, D):                   # captured energy of a scalar pulse at each model time
    return (al*a(EPS) + be*a(EPS-D))**2, al*al + be*be + 2*al*be*a(D)
def R_split(g1, g2, e1, e2, D, dmax=None):   # orthogonal split class in the plane (e1, e2); min over axis angle
    best = np.inf
    for ph in np.linspace(0, np.pi/2, 46, endpoint=False):
        ef, es = np.cos(ph)*e1+np.sin(ph)*e2, -np.sin(ph)*e1+np.cos(ph)*e2
        cf, Ef = scal_curve(g1@ef, g2@ef, D); cs, Es = scal_curve(g1@es, g2@es, D)
        if dmax is None:
            cap = cf.max() + cs.max()
        else:
            w = int(round(dmax/(EPS[1]-EPS[0])))
            # max over |i-j|<=w of cf[i]+cs[j]: running max of cs
            from scipy.ndimage import maximum_filter1d
            cap = (cf + maximum_filter1d(cs, 2*w+1, mode="nearest")).max()
        best = min(best, Ef + Es - cap)
    return best
def R_cal(g1, g2, e1, e2, D, ph, dt):        # calibrated split operator: pulse = ef ef^T s(t) + es es^T s(t-dt)
    ef, es = np.cos(ph)*e1+np.sin(ph)*e2, -np.sin(ph)*e1+np.cos(ph)*e2
    # data: sum_k [ (gk.ef) ef s(t-tk) + (gk.es) es s(t-tk-dt) ]; model: same form, one vector h in plane, one time
    # components along ef and es are orthogonal channels; model time shared (slow = fast + dt): joint max over eps
    cf, Ef = scal_curve(g1@ef, g2@ef, D); cs, Es = scal_curve(g1@es, g2@es, D)
    return Ef + Es - (cf + cs).max()         # the same eps in both channels (the dt shift cancels)
def dc(rng):
    s = rng.normal(size=3); s /= np.linalg.norm(s)
    n = rng.normal(size=3); n -= (n@s)*s; n /= np.linalg.norm(n)
    return np.outer(s, n) + np.outer(n, s)
def rays(depth=9.0, dist=15.0, ns=8):
    out = []
    for az in np.arange(ns)*2*np.pi/ns + 0.2:
        g = np.array([dist*np.cos(az), dist*np.sin(az), -depth]); g /= np.linalg.norm(g)
        e1 = np.cross(g, [0, 0, 1.0]); e1 /= np.linalg.norm(e1); e2 = np.cross(g, e1)
        out.append((g, e1, e2))
    return out
def amps(M, ray, vpvs=1.73):
    g, e1, e2 = ray
    Mg = M@g
    return (g@Mg)*g, (Mg - (g@Mg)*g)*vpvs**3
rng = np.random.default_rng(7)
RY = rays(); SNR = 12.0; NP = 40
pairs = []
for _ in range(NP):
    M1, M2 = dc(rng), dc(rng)
    A1 = [amps(M1, r) for r in RY]; A2 = [amps(M2, r) for r in RY]
    for A in (A1, A2):
        E = sum(p@p + s@s for p, s in A); sc = SNR/np.sqrt(E)
        for i in range(len(A)): A[i] = (A[i][0]*sc, A[i][1]*sc)
    cal = [(rng.uniform(0, np.pi), rng.uniform(0.03, 0.15)) for _ in RY]
    pairs.append((A1, A2, cal))
print("co-located pairs with DIFFERENT random double couples, SNR 12 per event; median over", NP, "pairs of delta^2 (whitened)")
print("info limit delta^2 = 10.8; current test 50 % acceptance near delta^2 = 38 (E.3)")
print(f"{'dtau':>6s} {'P':>8s} {'S one pulse':>12s} {'S free split':>13s} {'split<=0.1s':>12s} {'split<=0.2s':>12s} {'calibrated':>11s}")
for D in (0.02, 0.03, 0.05, 0.08, 0.1, 0.15, 0.25, 0.35, 0.5):
    rows = []
    for A1, A2, cal in pairs:
        P = sum(R_vec(A1[i][0], A2[i][0], D) for i in range(8))
        S1 = sum(R_vec(A1[i][1], A2[i][1], D) for i in range(8))
        Sf = sum(R_split(A1[i][1], A2[i][1], RY[i][1], RY[i][2], D) for i in range(8))
        S01 = sum(R_split(A1[i][1], A2[i][1], RY[i][1], RY[i][2], D, 0.1) for i in range(8))
        S02 = sum(R_split(A1[i][1], A2[i][1], RY[i][1], RY[i][2], D, 0.2) for i in range(8))
        Sc = sum(R_cal(A1[i][1], A2[i][1], RY[i][1], RY[i][2], D, *cal[i]) for i in range(8))
        rows.append((P, S1, Sf, S01, S02, Sc))
    m = np.median(np.array(rows), axis=0)
    print(f"{D:6.2f} {m[0]:8.2f} {m[1]:12.1f} {m[2]:13.2f} {m[3]:12.2f} {m[4]:12.2f} {m[5]:11.1f}")
