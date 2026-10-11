"""
c9_lib.py - C9: travel-time-free physics (polarization, amplitudes, S/P geometry) as location evidence.  SEED 999 ONLY.

Geometry helpers (vectorised over nodes), travel-time-free location statistics for a set of per-station P and S peak
vectors v (2, ns, 3) in matched-filter noise units, peak readers without travel times, and detector loops (BOMP-PK,
MT-FW cap 6) that expose the per-atom peak vectors and accept an in-loop relocation rule.

Travel-time-free statistics of a candidate position x (larger = better), all computed from v only:
  mt     ||Q_x^T v||^2                      moment-tensor projection energy (P2: spreading 1/r, radiation pattern,
                                            P and S polarization, all at once; times free)
  pol    sum_i (v_Pi.n_xi)^2 + sum_i (|v_Si|^2 - (v_Si.n_xi)^2)      P along the ray, S transverse to it (P1 + P3)
  ppol   sum_i (v_Pi.n_xi)^2                P polarization only (P1, energy form)
  spol   sum_i (|v_Si|^2 - (v_Si.n_xi)^2)   S transversality only (P3)
  amp    - sum_b var_i (log max(|v_bi|, 1) + log d_xi)               1/r spreading only, no radiation pattern (P2b)
  gfree  explained energy of the fit v_i ~ a_i G_i(x) m with a free gain a_i >= 0 per STATION (P3: S/P amplitude
         ratio and P/S polarization geometry, no inter-station amplitudes); alternating, 3 rounds from the LS fit
n_xi is the straight ray of the detectors' homogeneous model (or the circular ray of the gradient medium).
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

import c8_lib as C8  # noqa: E402
import c6_lib as L  # noqa: E402
import xc  # noqa: E402
import task2c_methods as M  # noqa: E402
import task2d_methods as D2  # noqa: E402
import task1x_mtglrt as X  # noqa: E402
from station_corrections import MIN_PS  # noqa: E402
from task0_common import AP, AS, DT, NT, VP, VS  # noqa: E402
from task1w_mech import CP, CS  # noqa: E402

SEED = 999
EZ = np.array([0.0, 0.0, 1.0])
CAPV = 6.0
R_RELOC = 12.0          # km, epicentral search radius of every relocation (fixed before any run; task range 10-15)
GAP_BLIND = 12
W_WIDE = 24             # samples (1.5 s): half-width of the wide peak re-read
PEN_LAM = 0.5           # weight of the unexplained-energy penalty (fixed before any run)


# ============================================================================================ geometry
def rays_hom(XN, STA):
    """straight rays: unit vectors node -> station (N, ns, 3) and distances (N, ns)."""
    dv = np.asarray(STA, float)[None] - np.asarray(XN, float)[:, None]
    d = np.linalg.norm(dv, axis=2)
    return dv / d[..., None], d


def rays_grad(XN, STA, vs0=2.8, g=0.04):
    """arrival direction at the station of the circular ray in v = v0 + g z (gen_realistic._ray, vectorised).  P and S
    share the ray (same g / v0)."""
    XN = np.asarray(XN, float)
    STA = np.asarray(STA, float)
    dxy = STA[None, :, :2] - XN[:, None, :2]
    Xh = np.maximum(np.linalg.norm(dxy, axis=2), 1e-9)
    eh = np.concatenate([dxy / Xh[..., None], np.zeros(Xh.shape + (1,))], axis=2)
    dz = STA[None, :, 2] - XN[:, None, 2]
    vs = vs0 + g * XN[:, None, 2]
    vr = vs0 + g * STA[None, :, 2]
    gxc = (g * Xh * Xh + dz * (vs + vr)) / (2.0 * Xh)
    gR = np.hypot(gxc, vs)
    nr = (vr[..., None] * eh - (g * Xh - gxc)[..., None] * EZ) / gR[..., None]
    return nr


def operators_hom(XN, STA):
    """task1x_mtglrt.radiation_operator vectorised: G (N, 6 ns, 6), rows (b, i, c)."""
    n, d = rays_hom(XN, STA)
    N, ns = d.shape
    nx, ny, nz = n[..., 0], n[..., 1], n[..., 2]
    c = np.stack([nx * nx, ny * ny, nz * nz, 2 * nx * ny, 2 * nx * nz, 2 * ny * nz], axis=-1)       # (N, ns, 6)
    GP = n[..., :, None] * c[..., None, :] * (AP / np.sqrt(CP)) / d[..., None, None]                # (N, ns, 3, 6)
    z = np.zeros_like(nx)
    B = np.stack([np.stack([nx, z, z, ny, nz, z], -1), np.stack([z, ny, z, nx, z, nz], -1),
                  np.stack([z, z, nz, z, nx, ny], -1)], axis=-2)                                    # (N, ns, 3, 6)
    Pm = np.eye(3) - n[..., :, None] * n[..., None, :]
    GS = np.matmul(Pm, B) * (AS / np.sqrt(CS)) / d[..., None, None]
    return np.concatenate([GP.reshape(N, ns * 3, 6), GS.reshape(N, ns * 3, 6)], axis=1)


def bases(G):
    U, s, _ = np.linalg.svd(G, full_matrices=False)
    keep = s > 1e-8 * s.max(axis=1, keepdims=True)
    return U * keep[:, None, :]


class Geo:
    """time-free operators of a node set for one network."""

    def __init__(self, XN, STA, grad=False, need_g=False):
        self.XN = np.asarray(XN, float)
        self.STA = np.asarray(STA, float)
        self.ns = len(self.STA)
        self.n, self.d = rays_hom(self.XN, self.STA)
        self.ng = rays_grad(self.XN, self.STA) if grad else None
        G = operators_hom(self.XN, self.STA)
        self.Q = bases(G)
        self.G = G if need_g else None


def cap(v, c=CAPV):
    if c is None:
        return v
    nv = np.linalg.norm(v, axis=-1, keepdims=True)
    return v * np.minimum(1.0, c / np.maximum(nv, 1e-12))


# ============================================================================================ statistics
def st_mt(geo, v, sel=None, Q=None):
    Q = geo.Q if Q is None else Q
    Q = Q if sel is None else Q[sel]
    y = np.einsum("r,nrk->nk", v.reshape(-1), Q)
    return (y * y).sum(-1)


def st_pol(geo, v, sel=None, grad=False, parts=False):
    n = geo.ng if grad else geo.n
    n = n if sel is None else n[sel]
    pp = (np.einsum("nic,ic->ni", n, v[0]) ** 2).sum(-1)
    sp = ((v[1] ** 2).sum() - (np.einsum("nic,ic->ni", n, v[1]) ** 2).sum(-1))
    return (pp, sp) if parts else pp + sp


def st_ppol_angle(geo, v, sel=None, grad=False, gate=3.0):
    """P1 grid form: - sum_i w_i sin^2(angle(v_Pi, n_xi)), w_i = min(1, |v_Pi|^2 / gate^2)."""
    n = geo.ng if grad else geo.n
    n = n if sel is None else n[sel]
    a2 = (v[0] ** 2).sum(-1)
    w = np.minimum(1.0, a2 / gate ** 2)
    c2 = np.einsum("nic,ic->ni", n, v[0]) ** 2 / np.maximum(a2, 1e-12)[None]
    return -(w[None] * (1.0 - c2)).sum(-1)


def st_amp(geo, v, sel=None):
    d = geo.d if sel is None else geo.d[sel]
    la = np.log(np.maximum(np.linalg.norm(v, axis=-1), 1.0))                # (2, ns)
    z = la[None] + np.log(d)[:, None, :]                                    # (N, 2, ns)
    return -z.var(axis=-1).sum(-1)


def st_gfree(geo, v, sel=None, rounds=3):
    G = geo.G if sel is None else geo.G[sel]
    N = len(G)
    ns = geo.ns
    Gb = G.reshape(N, 2, ns, 3, 6).transpose(0, 2, 1, 3, 4).reshape(N, ns, 6, 6)       # per station (P3, S3) x 6
    vb = np.concatenate([v[0], v[1]], axis=-1)                                           # (ns, 6)
    GG = np.matmul(Gb.transpose(0, 1, 3, 2), Gb)                                         # (N, ns, 6, 6)
    Gv = np.einsum("nirk,ir->nik", Gb, vb)                                               # (N, ns, 6)
    a = np.ones((N, ns))
    for _ in range(rounds + 1):
        A = np.einsum("ni,nikl->nkl", a * a, GG) + 1e-9 * np.eye(6)
        m = np.linalg.solve(A, np.einsum("ni,nik->nk", a, Gv)[..., None])[..., 0]
        g = np.einsum("nirk,nk->nir", Gb, m)
        a = np.maximum(np.einsum("nir,ir->ni", g, vb), 0.0) / np.maximum((g * g).sum(-1), 1e-12)
    res = ((vb[None] - a[..., None] * g) ** 2).sum((-1, -2))
    return (vb ** 2).sum() - res


def azimuth_intersection(STA, vP, gate=3.0):
    """P1 closed form: horizontal lines through the stations along the horizontal projection of the P vector (sign
    irrelevant for a line), weighted least-squares intersection.  Weights min(1, |v|^2 / gate^2) x horizontal fraction."""
    A = np.zeros((2, 2))
    b = np.zeros(2)
    for s, v in zip(np.asarray(STA, float), vP):
        h = v[:2]
        nh = np.linalg.norm(h)
        if nh < 1e-9:
            continue
        h = h / nh
        w = min(1.0, float(v @ v) / gate ** 2) * (nh ** 2 / max(float(v @ v), 1e-12))
        Pm = np.eye(2) - np.outer(h, h)
        A += w * Pm
        b += w * Pm @ s[:2]
    if np.linalg.cond(A) > 1e8:
        return np.array([np.nan, np.nan])
    return np.linalg.solve(A, b)


# ============================================================================================ peak readers
def read_true(MFn, tarr):
    """vectors at the energy peak within +-2 samples of the true arrival times tarr (2, ns) (oracle read)."""
    ns = MFn.shape[0]
    E = (MFn ** 2).sum(1)
    v = np.zeros((2, ns, 3))
    K = np.zeros((2, ns), int)
    for b in range(2):
        for i in range(ns):
            k0 = int(round(tarr[b, i] / DT))
            ks = np.clip(np.arange(k0 - 2, k0 + 3), 0, NT - 1)
            k = int(ks[np.argmax(E[i, ks])])
            v[b, i], K[b, i] = MFn[i, :, k], k
    return v, K


def read_blind(MFn):
    """per station, the pair of samples kP + GAP_BLIND <= kS maximising E(kP) + E(kS) over the whole window; the
    earlier one is called P.  No travel time, no location.  GAP_BLIND = 12 samples (0.75 s; with the registered GAP_BLIND
    = 6 the pair was often the S peak and its own side lobe in the noise-free validation on network 199)."""
    ns = MFn.shape[0]
    E = (MFn ** 2).sum(1)
    v = np.zeros((2, ns, 3))
    K = np.zeros((2, ns), int)
    for i in range(ns):
        e = E[i]
        pm = np.maximum.accumulate(e)
        pa = np.zeros(NT, int)
        cur = 0
        for k in range(1, NT):
            if e[k] > e[cur]:
                cur = k
            pa[k] = cur
        tot = np.full(NT, -np.inf)
        tot[GAP_BLIND:] = e[GAP_BLIND:] + pm[:NT - GAP_BLIND]
        kS = int(np.argmax(tot))
        kP = int(pa[kS - GAP_BLIND])
        v[0, i], v[1, i] = MFn[i, :, kP], MFn[i, :, kS]
        K[0, i], K[1, i] = kP, kS
    return v, K


def read_wide(MFn, I, W=W_WIDE):
    """wide re-read around an atom's samples I (2, ns): P = energy peak in [I_P - W, min(I_P + W, mid)], S = energy
    peak in [max(I_S - W, mid + 1), I_S + W], mid = (I_P + I_S) // 2."""
    ns = MFn.shape[0]
    E = (MFn ** 2).sum(1)
    v = np.zeros((2, ns, 3))
    K = np.zeros((2, ns), int)
    for i in range(ns):
        kp, ks = int(I[0, i]), int(I[1, i])
        mid = (kp + ks) // 2
        lo, hi = max(0, kp - W), max(min(kp + W, mid), max(0, kp - W))
        k = lo + int(np.argmax(E[i, lo:hi + 1]))
        v[0, i], K[0, i] = MFn[i, :, k], k
        lo, hi = min(max(ks - W, mid + 1), NT - 1), min(NT - 1, ks + W)
        hi = max(hi, lo)
        k = lo + int(np.argmax(E[i, lo:hi + 1]))
        v[1, i], K[1, i] = MFn[i, :, k], k
    return v, K


def time_locate(geo, K, statics=None):
    """arrival-time location from peak samples K (2, ns): homogeneous model, origin time profiled out, LS."""
    t = K * DT - (0.0 if statics is None else np.asarray(statics).T)
    T_ = np.stack([geo.d / VP, geo.d / VS], axis=1)                          # (N, 2, ns)
    r = t[None] - T_
    r = r - r.mean(axis=(1, 2), keepdims=True)
    return -(r * r).sum((1, 2))


# ============================================================================================ atoms: travel-time-free features
STATS = ("mt", "pol", "ppol", "amp")


def near_nodes(garr, src, radius=R_RELOC):
    return np.where(np.linalg.norm(garr[:, :2] - np.asarray(src)[None, :2], axis=1) <= radius)[0]


def stat_on(geo, name, v, sel, Q=None):
    if name == "mt":
        return st_mt(geo, v, sel, Q)
    if name == "pol":
        return st_pol(geo, v, sel)
    if name == "ppol":
        return st_pol(geo, v, sel, parts=True)[0]
    if name == "amp":
        return st_amp(geo, v, sel)
    raise ValueError(name)


def relocate(geo, name, v, gi, Q=None, allowed=None):
    """node maximising the statistic within R_RELOC km (epicentral) of node gi; ties and no-gain keep gi."""
    sel = near_nodes(geo.XN, geo.XN[gi])
    if allowed is not None:
        sel = sel[allowed[sel] | (sel == gi)]
    s = stat_on(geo, name, v, sel, Q)
    k = int(np.argmax(s))
    s0 = s[int(np.where(sel == gi)[0][0])]
    return int(sel[k]) if s[k] > s0 else int(gi)


def atom_features(geo, MFn, gi, I, v_sel, wsite=None):
    """travel-time-free features of an accepted atom.  v_sel: residual vectors at the read samples at selection time
    (uncapped).  Returns consistency fractions at the atom's node and the relocated nodes for each statistic, with the
    narrow (atom's own samples; capped at 6) and the wide (+-1.5 s re-read on the data MFn; capped at 6) vector sets."""
    f = {}
    vn = cap(v_sel)
    vw = cap(read_wide(MFn, I)[0])
    e = float((vn ** 2).sum())
    one = np.array([gi])
    f["rho_mt"] = float(st_mt(geo, vn, one)[0] / max(e, 1e-12))
    f["rho_pol"] = float(st_pol(geo, vn, one)[0] / max(e, 1e-12))
    for tag, v in (("n", vn), ("w", vw)):
        for name in STATS:
            f[f"{name}_{tag}"] = relocate(geo, name, v, gi)
    if wsite is not None:
        f["mtos_n"] = relocate(geo, "mt", vn, gi, Q=wsite)
        f["mtos_w"] = relocate(geo, "mt", vw, gi, Q=wsite)
    return f


# ============================================================================================ detectors
def mtfw_c9(dops, geo, MFn, kmax, reloc=None, pen=0.0, min_sep=5.0, sel=None):
    """c6_lib.mtfw_dense / c8_lib.mtfw_general "ls" cap 6 (identical when reloc is None and pen == 0), returning per
    atom: node, origin index, raw score, read samples, residual vectors at the read samples (uncapped).
    reloc in STATS: after selection the atom is moved to the maximum of that travel-time-free statistic (capped
    vectors, fixed samples) within R_RELOC km among non-excluded nodes; exclusion and deflation use the new node.
    pen > 0: selection statistic  s - pen x (|v|^2 / 2 - s)  (unexplained energy at the read samples penalised).
    sel (post hoc, added after the first look at DENSE thin): "selpol" = selection statistic s x rho_pol, rho_pol the
    fraction of the read energy that is polarization-admissible at the node (P along the ray, S transverse);
    "penpol" = s - PEN_LAM x (polarization-inadmissible energy) / 2."""
    ns = dops.ns
    R = MFn.copy()
    garr, IDX = dops.garr, dops.IDX
    N, T_ = IDX.shape[:2]
    st = np.arange(ns)[None, None, None, :]
    Gs, Q = dops.Gs, dops.Qst
    excluded = np.zeros(N, bool)
    atoms = []
    for _ in range(kmax):
        if excluded.all():
            break
        K, V = C8._leaf(R, CAPV)
        tpk = K[st, IDX]
        vf = V[st, tpk].reshape(N, T_, 6 * ns)
        y = np.einsum("gtr,grk->gtk", vf, Q)
        sc = (y * y).sum(-1) / 2.0
        if pen > 0:
            sc = sc - pen * ((vf * vf).sum(-1) / 2.0 - sc)
        if sel is not None:
            v5 = vf.reshape(N, T_, 2, ns, 3)
            al = np.einsum("ntbic,nic->ntbi", v5, geo.n) ** 2                     # energy along the ray
            etot = (vf * vf).sum(-1)
            epol = al[:, :, 0].sum(-1) + (v5[:, :, 1] ** 2).sum((-1, -2)) - al[:, :, 1].sum(-1)
            if sel == "selpol":
                sc = sc * epol / np.maximum(etot, 1e-12)
            else:
                sc = sc - PEN_LAM * (etot - epol) / 2.0
        sc[excluded] = -np.inf
        best = sc.max(axis=1)
        gi = int(np.argmax(best))
        tj = int(np.argmax(sc[gi]))
        s = float(sc[gi, tj])
        if not np.isfinite(s):
            break
        tp = tpk[gi, tj].copy()
        vcap = vf[gi, tj]
        vraw = R[np.arange(ns)[None, :], :, tp]                                        # (2, ns, 3)
        g0 = gi
        if reloc is not None:
            gi = relocate(geo, reloc, vcap.reshape(2, ns, 3), gi, allowed=~excluded)
        src = garr[gi]
        atoms.append({"n": int(gi), "n0": int(g0), "t": int(tj), "s": s, "I": tp, "v": vraw})
        excluded = excluded | (np.linalg.norm(garr[:, :2] - src[:2], axis=1) <= min_sep)
        m, *_ = np.linalg.lstsq(Gs[gi], vcap, rcond=None)
        pred = (Gs[gi] @ m).reshape(2, ns, 3)
        for b in range(2):
            for i in range(ns):
                t0 = int(tp[b, i])
                lo, hi = max(0, t0 - X.AC_HALF), min(NT, t0 + X.AC_HALF + 1)
                R[i, :, lo:hi] -= np.outer(pred[b, i], X.AC[X.AC_HALF - (t0 - lo):X.AC_HALF + (hi - t0)])
    return atoms


def bomp_c9(PD, geo, MFn, kmax, reloc=None, min_sep=5.0):
    """task2d_methods.bomp_pk_record (identical when reloc is None) returning the atoms with raw scores and vectors.
    reloc: the entering atom (n, t, I) becomes (n', t, I), n' the maximum of the travel-time-free statistic of the
    residual vectors at I (uncapped, as BOMP-PK reads them) within R_RELOC km among non-excluded nodes."""
    garr = PD.BD.garr
    ns = PD.ns
    excluded = np.zeros(PD.N, bool)
    S, atoms = [], []
    R = MFn.copy()
    for _ in range(kmax):
        e, I, _h = PD.energy(R, True)
        e[excluded] = -np.inf
        j = int(np.argmax(e))
        n, t = divmod(j, PD.T)
        s = float(e[n, t]) / 2.0
        if not np.isfinite(s):
            break
        Ia = I[n, t].copy()
        vraw = R[np.arange(ns)[None, :], :, Ia]
        n0 = n
        if reloc is not None:
            n = relocate(geo, reloc, vraw, n, allowed=~excluded)
        src = garr[n]
        atoms.append({"n": int(n), "n0": int(n0), "t": int(t), "s": s, "I": Ia, "v": vraw})
        S.append((n, t, Ia))
        excluded = M.excl_mask(excluded, garr, src, min_sep)
        Kg = PD.gram(S)
        b = PD.corr_atoms(MFn, S).ravel()
        U, *_ = np.linalg.lstsq(Kg, b, rcond=None)
        R = PD.residual(MFn, S, U.reshape(len(S), 6))
    return atoms


def record(atoms, garr, truth, key="n", runmin=True, score=None, dup="free"):
    """record [(score, hit)] in acceptance order.  key: which node labels the atom.  runmin: running minimum of the
    scores (the registered convention).  dup: "free" = the registered metric (an atom within 5 km of a true event is a
    hit, repeated hits are neither hits nor false); "strict" = an atom whose location is within 5 km (epicentral) of an
    earlier atom's location is counted as a FALSE event (relocation may break the 5 km exclusion)."""
    rec, run, prev = [], np.inf, []
    for k, a in enumerate(atoms):
        s = a["s"] if score is None else score[k]
        run = min(run, s)
        src = garr[a[key] if isinstance(key, str) else key[k]]
        h = int(M.hit_of(src, truth))
        if dup == "strict" and any(np.linalg.norm(src[:2] - q[:2]) <= 5.0 for q in prev):
            h = -1
        prev.append(src)
        rec.append((float(run if runmin else s), h))
    return rec
