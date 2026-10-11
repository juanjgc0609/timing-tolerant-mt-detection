"""C7 validation (seed 999, network j = 40 of each cell, window stream TAG 3705 - all new).
A. structural checks of c7_lib against the registered code and against explicit data-domain computations;
B. noise-free invariance: single source at a grid node, origin on the tau grid, homogeneous medium (the detectors'
   medium), no statics, no site/polarisation perturbation, no coda, no clutter, no noise; pulse fc x phi grid (and
   out-of-family pulses).  Reported: first-iteration statistic at the best block within 5 km of the source divided by
   the whitened signal energy (the fraction of the signal energy the detector's block captures).
    python c7_validate.py > out/c7_validation.txt"""
import time

import numpy as np

import c6_lib as L
import c7_lib as C
import xc

T_, D2, M, X = xc.T, xc.D2, xc.M, xc.X
GR = C.GR
grid = xc.make_grid()
garr = np.asarray(grid)
BASES = {"Q1": C.Basis(C.basis_q1(), "Q1"), "PH": C.Basis(C.basis_ph(), "PH"),
         "PHF4": C.Basis(C.basis_phf(4)[0], "PHF4"), "PHF3": C.Basis(C.basis_phf(3)[0], "PHF3")}


def mx(a, b):
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b)))) if np.size(a) else 0.0


def part_a():
    print("A. STRUCTURAL CHECKS")
    B1 = BASES["Q1"]
    print(f"  A1 q=1 kernel == task1x AC: max diff {mx(B1.Xk[0, 0], X.AC):.1e}; PH psi_1 autocorrelation vs AC "
          f"(lags within +-48): {mx(BASES['PH'].Xk[0, 0][64 - 48:64 + 49], X.AC):.1e}")
    for nm in ("PH", "PHF4"):
        Xk = BASES[nm].Xk
        print(f"  A1 {nm}: X(0) = I max dev {mx(Xk[:, :, BASES[nm].XH], np.eye(BASES[nm].q)):.1e}; "
              f"X[j,J](tau) = X[J,j](-tau) max dev {mx(Xk, np.transpose(Xk, (1, 0, 2))[:, :, ::-1]):.1e}")
    cell = "ladder_6"
    nc = xc.NetCfg(cell, 40, "thin")
    ops = xc.Ops(nc.STA, nc.corr, grid)
    mops = L.ModelOps(nc.STA, grid, L.Homog(), nc.corr)
    Q1 = C.QDict(ops.BD, B1, mops.Qst)
    st = xc.window_stream(nc, 3705)
    ws = [next(st) for _ in range(6)]
    rng = np.random.default_rng(5)
    # A2 Gram / correlation / residual for q = 1 vs PKDict
    atoms = []
    for _ in range(4):
        n, t = int(rng.integers(ops.BD.N)), int(rng.integers(ops.BD.T))
        I = ops.BD.idx[n, t] + rng.integers(-4, 5, size=ops.BD.idx[n, t].shape)
        atoms.append((n, t, np.clip(I, 0, xc.NT - 1)))
    K1 = Q1.cross_gram(atoms, atoms)[:, :, 0, :, :, 0]
    K0 = ops.PD.cross_gram(atoms, atoms)
    U = rng.normal(size=(4, 6))
    MFn = ws[0]["MFn"]
    print(f"  A2 q=1 vs PKDict: cross_gram {mx(K1, K0):.1e} (scale {np.abs(K0).max():.1e}); corr_atoms "
          f"{mx(Q1.corr_atoms(MFn[None], atoms)[..., 0], ops.PD.corr_atoms(MFn, atoms)):.1e}; residual "
          f"{mx(Q1.residual(MFn[None], atoms, [u[:, None] for u in U])[0], ops.PD.residual(MFn, atoms, U)):.1e}")
    # A3 mtfw_q (q=1, cap) == c6_lib.mtfw_dense ; A5 NX (bomp_q q=1) vs BOMP-PK
    same_fw, same_nx, dsc = 0, 0, []
    for w in ws:
        r0, l0 = L.mtfw_dense(mops, w["MFn"], w["truth"], 5)
        r1, l1 = C.mtfw_q(Q1, mops.IDX, w["MFn"][None], w["truth"], 5)
        same_fw += ([h for _s, h in r0] == [h for _s, h in r1] and l0 == l1
                    and max(abs(a[0] - b[0]) / abs(a[0]) for a, b in zip(r0, r1)) < 1e-12)
        rb, lb = xc.bomp_pk(ops.PD, w["MFn"], w["truth"], 5)
        rn, ln = C.bomp_q(Q1, w["MFn"][None], w["truth"], 5, "R1")
        same_nx += ([h for _s, h in rb] == [h for _s, h in rn] and [l[:3] for l in lb] == [l[:3] for l in ln])
        dsc.append(abs(rb[0][0] - rn[0][0]) / rb[0][0])
    print(f"  A3 MT-FW-Q with q=1 reproduces c6_lib.mtfw_dense hits, locations and scores (rel. 1e-12): {same_fw}/{len(ws)} windows")
    print(f"  A5 BOMP-PK-NX (q=1, no P/S cross term) vs BOMP-PK: identical hit/location sequence in {same_nx}/{len(ws)}"
          f" windows; first-score relative difference median {np.median(dsc):.1e} (the P/S cross term)")
    # A4 PH / PHF4: filter outputs of explicit data-domain atoms == synth prediction; Gram == explicit inner products
    sd = 3.0
    for nm in ("PH", "PHF4"):
        bs = BASES[nm]
        QD = C.QDict(ops.BD, bs, mops.Qst)
        q, Lh = bs.q, bs.L // 2
        ys, Us = [], []
        for a in atoms[:3]:
            Uc = rng.normal(size=(6, q))
            pred = np.einsum("pick,kj->picj", ops.BD.Gb[a[0]], Uc)
            y = np.zeros((ops.BD.ns, 3, xc.NT))
            for b in range(2):
                for i in range(ops.BD.ns):
                    t0 = int(a[2][b, i])
                    for j in range(q):
                        lo, hi = max(0, t0 - Lh), min(xc.NT, t0 + Lh + 1)
                        y[i, :, lo:hi] += np.outer(pred[b, i, :, j], bs.B[j][Lh - (t0 - lo):Lh + (hi - t0)])
            ys.append(y)
            Us.append(Uc)
        Fy = bs.filt(ys[0] * sd / np.sqrt(C.E_R), sd)
        Fp = -QD.residual(np.zeros_like(Fy), atoms[:1], Us[:1])
        inner = slice(2 * bs.L, xc.NT - 2 * bs.L)
        Kq = QD.cross_gram(atoms[:3], atoms[:3])
        eg = 0.0
        for a in range(3):
            for b in range(3):
                g_exp = float(np.sum(ys[a] * ys[b]))
                g_mod = float(np.einsum("kj,kjlJ,lJ->", Us[a], Kq[a, :, :, b], Us[b]))
                eg = max(eg, abs(g_exp - g_mod) / max(1.0, abs(g_exp)))
        cb = QD.corr_atoms(Fy, atoms[:3])
        ce = max(abs(float(np.sum(cb[b] * Us[b])) - float(np.sum(ys[0] * ys[b]))) / max(1.0, float(np.sum(ys[0] ** 2)))
                 for b in range(3))
        print(f"  A4 {nm}: filt(explicit atom) vs synth prediction (interior) max rel diff "
              f"{mx(Fy[..., inner], Fp[..., inner]) / np.abs(Fp).max():.1e}; Gram vs explicit <y_a, y_b> max rel "
              f"diff {eg:.1e}; corr_atoms(filt(y_a)) . U_b vs <y_a, y_b> {ce:.1e}")
    # A6 rank-1 ALS: monotone objective, rank-1 residual >= full residual, rank-1 <= single-c fit
    bs = BASES["PH"]
    QD = C.QDict(ops.BD, bs, mops.Qst)
    w = ws[1]
    F = bs.filt(w["D"], w["sd"])
    ok = True
    for na in (1, 2, 3):
        at = atoms[:na]
        Uf = C.fit_full(QD, F, at)
        C0 = [np.array([1.0, 0.0])] * na
        U1, _c = C.fit_rank1(QD, F, at, C0)
        U0, _c0 = C.fit_rank1(QD, F, at, C0, iters=1)
        ef = float(np.sum(QD.residual(F, at, Uf) ** 2))
        e1 = float(np.sum(QD.residual(F, at, U1) ** 2))
        e0 = float(np.sum(QD.residual(F, at, U0) ** 2))
        ok &= (ef <= e1 + 1e-9 * ef) and (e1 <= e0 + 1e-9 * e0)
    print(f"  A6 rank-1 ALS on the PH basis: full-fit residual <= rank-1 <= first ALS pass (1-3 atoms): "
          f"{'PASS' if ok else 'FAIL'}")
    # timing
    for nm, (bsn, mode, cap) in {"BOMP-PK": (None, None, None), "NX": ("Q1", "R1", False), "PH": ("PH", "R1", False),
                                 "PH-FULL": ("PH", "FULL", False), "PHF4": ("PHF4", "R1", False),
                                 "PHF4-FULL": ("PHF4", "FULL", False), "PH-CAP": ("PH", "R1", True)}.items():
        t0 = time.time()
        for w in ws:
            if bsn is None:
                xc.bomp_pk(ops.PD, w["MFn"], w["truth"], 5)
            else:
                QD = C.QDict(ops.BD, BASES[bsn], mops.Qst)
                C.bomp_q(QD, BASES[bsn].filt(w["D"], w["sd"], w["MFn"]), w["truth"], 5, mode, cap)
        print(f"  timing {nm:10s} {(time.time() - t0) / len(ws):.3f} s/window (ladder_6, kmax 5, 1 process)")


class FixedPulse:
    def __init__(self, fn):
        self.fn = fn

    def __enter__(self):
        self._sw = GR.source_wavelet
        fn = self.fn
        GR.source_wavelet = lambda t, fc, phi=0.0: fn(t)
        return self

    def __exit__(self, *a):
        GR.source_wavelet = self._sw


def part_b():
    print("\nB. NOISE-FREE INVARIANCE (captured fraction = first-iteration block statistic at the best block within 5 km"
          " / whitened signal energy)")
    cfg = {"VS0": 3.0, "g": 0.0, "static_sd": 0.0, "dir_deg": 0.0, "site_sd": 0.0, "coda_ratio": 0.0, "noise": 0.0,
           "fc_range": (2.0, 2.0), "phase_max_deg": 0.0}
    sd = GR.SNR * 0.02 * np.sqrt(C.E_R)
    nets = [(cell, T_.network(999, cell, 40)[0]) for cell in T_.CELLS]
    dets = {"BOMP-PK": None, "NX": ("Q1", "R1"), "PH": ("PH", "R1"), "PH-FULL": ("PH", "FULL"),
            "PHF4": ("PHF4", "R1"), "PHF4-FULL": ("PHF4", "FULL")}
    pulses = {}
    for fc in (1.4, 1.75, 2.0, 2.4, 2.8):
        for ph in (-45, -22.5, 0, 22.5, 45):
            pulses[("design", fc, ph)] = (lambda fc_, ph_: (lambda t: GR.__dict__["_sw0"](t, fc_, np.radians(ph_))))(fc, ph)
    for fc in (1.4, 2.0, 2.8):
        for ph in (-90, 90):
            pulses[("rot90", fc, ph)] = (lambda fc_, ph_: (lambda t: GR.__dict__["_sw0"](t, fc_, np.radians(ph_))))(fc, ph)
    for fn in (C.gabor, C.berlage):
        for fc in (1.4, 2.0, 2.8):
            for ph in (-90, 0, 90):
                pulses[(fn.__name__, fc, ph)] = (lambda fn_, fc_, ph_: (lambda t: fn_(t, fc_, np.radians(ph_))))(fn, fc, ph)
    GR._sw0 = GR.source_wavelet
    rng = np.random.default_rng([999, 3705, 7])
    cases = []
    for cell, STA in nets:
        ops = xc.Ops(STA, np.zeros((len(STA), 2)), grid)
        mops_Qst = __import__("certified_fw").stack_Q(ops.Qs)
        QDs = {k: C.QDict(ops.BD, BASES[k], mops_Qst) for k in ("Q1", "PH", "PHF4")}
        for _s in range(3):
            n = int(rng.integers(len(garr)))
            while np.any(np.abs(garr[n][:2]) > 15):
                n = int(rng.integers(len(garr)))
            src, tau = garr[n], float(xc.DEFAULT_TAUS[int(rng.integers(5, 26))])
            mseed = int(rng.integers(1 << 30))
            near = np.where(np.linalg.norm(garr[:, :2] - src[None, :2], axis=1) <= 5.0)[0]
            for key, fn in pulses.items():
                with FixedPulse(fn):
                    D, _p, _t, _a = GR.generate_real(STA, [src], [tau], 0, np.random.default_rng(1),
                                                     np.random.default_rng(mseed), np.random.default_rng(2), 0.0,
                                                     cfg=cfg, statics=np.zeros((len(STA), 2)))
                Esig = float(np.sum(D ** 2)) * C.E_R / sd ** 2
                MF, _ = xc.matched(D)
                MFn = MF / sd
                row = {}
                for dn, spec in dets.items():
                    if spec is None:
                        e = ops.PD.energy(MFn, True)[0]
                    else:
                        QD = QDs[spec[0]]
                        e = QD.stat(BASES[spec[0]].filt(D, sd, MFn), spec[1])[0]
                    row[dn] = float(e[near].max()) / Esig
                cases.append((key, cell, row))
    GR.source_wavelet = GR._sw0

    def summ(sel, title):
        print(f"  {title}")
        print("    detector     " + "".join(f"{k:>11s}" for k in ("mean", "min", "max", "sd", "CV")))
        for dn in dets:
            v = np.array([r[dn] for (k, _c, r) in cases if sel(k)])
            print(f"    {dn:12s}" + "".join(f"{x:11.3f}" for x in (v.mean(), v.min(), v.max(), v.std(), v.std() / v.mean())))
    summ(lambda k: k[0] == "design", "design family (fc 1.4-2.8 x phase -45..45; 3 cells x 3 sources x 25 pulses)")
    print("  per phase at fc = 2.0 (mean over 9 source/network cases):")
    for ph in (-45, -22.5, 0, 22.5, 45):
        print(f"    phase {ph:6.1f}: " + "  ".join(f"{dn} {np.mean([r[dn] for (k, _c, r) in cases if k == ('design', 2.0, ph)]):.3f}" for dn in dets))
    print("  per fc at phase 0:")
    for fc in (1.4, 1.75, 2.0, 2.4, 2.8):
        print(f"    fc {fc:4.2f}: " + "  ".join(f"{dn} {np.mean([r[dn] for (k, _c, r) in cases if k == ('design', fc, 0)]):.3f}" for dn in dets))
    summ(lambda k: k[0] == "rot90", "phase +-90 deg (outside the design phase range), fc 1.4/2.0/2.8")
    summ(lambda k: k[0] == "gabor", "Gabor pulse (width 0.6/f0, phase -90/0/90), f0 1.4/2.0/2.8 - out of family")
    summ(lambda k: k[0] == "berlage", "Berlage pulse (causal t^2 exp(-2 pi f0 t) sin, phase -90/0/90) - out of family")


if __name__ == "__main__":
    t0 = time.time()
    part_a()
    part_b()
    print(f"\nvalidation wall time {time.time() - t0:.0f} s")
