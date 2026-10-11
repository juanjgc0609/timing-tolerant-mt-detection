"""C2 follow-up: why does MT-FW beat BOMP-PK in dense overlaps?  Same 180 dense windows as c2_run.py (thin, networks
10-15, TAG 3002, draw_dense, kmax 7), regenerated from the same streams.  Differences between the two detectors:
(1) MT-FW caps each station/phase peak vector at norm CAPV = 6 in its statistic (and in its amplitude fit);
(2) MT-FW scores ||Q^T v||^2 / 2 (no P/S cross-term in the Gram);
(3) MT-FW deflates sequentially (only the new atom's capped fit is subtracted), BOMP-PK refits all atoms jointly.
The running-minimum ("prefix") score is the same in both.  Ablations:
  MT-FW          dense MT-FW (c6_lib, homogeneous model; must reproduce the stored c2 MT-FW records)
  MT-FW-NOCAP    MT-FW without the cap (1)
  BOMP-PK-CAP    BOMP-PK whose selection statistic uses capped peak vectors (1); refit and deflation unchanged
  BOMP-PK-SEQ    BOMP-PK with MT-FW's sequential deflation (3): after each atom only that atom is fitted (uncapped LS
                 on the residual at its peak samples) and subtracted; no joint refit
    python c2_ablate.py thin 10,11,12,13,14,15 10 out/c2_ablate_thin.json"""
import json
import sys
import time

import numpy as np

import c2_run as C2
import c6_lib as L
import xc

M, D2, X = xc.M, xc.D2, xc.X


class PKDictCap(D2.PKDict):
    def energy(self, R, peak):
        BD = self.BD
        I = self.read_idx(R, peak)
        V = self.sample_idx(R, I)                                                  # (N, T, 2, ns, 3)
        nv = np.linalg.norm(V, axis=-1, keepdims=True)
        V = V * np.minimum(1.0, X.CAPV / np.maximum(nv, 1e-12))
        h = np.matmul(V.reshape(self.N, self.T, -1), BD.Gf)
        lag0 = BD.idx[:, :, 1, :] - BD.idx[:, :, 0, :]
        lag = I[:, :, 1, :] - I[:, :, 0, :]
        diff = np.any((M.ac_at(lag) != M.ac_at(lag0)) | (M.ac_at(-lag) != M.ac_at(-lag0)), axis=-1)
        Kinv = BD.Kinv
        if diff.any():
            Kinv = Kinv.copy()
            n_, t_ = np.nonzero(diff)
            Isub = I[n_, t_]
            lg = Isub[:, 1, :] - Isub[:, 0, :]
            wm, wp = M.ac_at(-lg), M.ac_at(lg)
            K = (self.GtG[n_] + np.einsum("mi,mikl->mkl", wm, self.C[n_]) + np.einsum("mi,milk->mkl", wp, self.C[n_]))
            K = 0.5 * (K + np.swapaxes(K, -1, -2))
            Kinv[n_, t_] = self.pinv_sym(K)
        e = np.einsum("ntk,ntk->nt", h, np.matmul(Kinv, h[..., None])[..., 0])
        return e, I, h


def bomp_seq(PD, MFn, truth, kmax, min_sep=5.0):
    garr = PD.BD.garr
    excluded = np.zeros(PD.N, bool)
    rec, run = [], np.inf
    R = MFn.copy()
    for _ in range(kmax):
        e, I, _h = PD.energy(R, True)
        e[excluded] = -np.inf
        j = int(np.argmax(e))
        n, t = divmod(j, PD.T)
        s = float(e[n, t]) / 2.0
        if not np.isfinite(s):
            break
        run = min(run, s)
        src = garr[n]
        rec.append((float(run), int(M.hit_of(src, truth))))
        excluded = M.excl_mask(excluded, garr, src, min_sep)
        a = (n, t, I[n, t].copy())
        K = PD.gram([a])
        b = PD.corr_atoms(R, [a]).ravel()
        u, *_ = np.linalg.lstsq(K, b, rcond=None)
        PD.synth_sub(R, n, a[2], u)
    return rec


def run(regime, nets, nwin, out):
    res = []
    grid = xc.make_grid()
    for cell in C2.CELLS:
        for j in nets:
            nc = xc.NetCfg(cell, j, regime)
            st = xc.window_stream(nc, C2.TAG, draw=C2.draw_dense)
            ops = L.ModelOps(nc.STA, grid, L.Homog(), nc.corr)
            PDc = PKDictCap(ops.BD)
            for wi in range(nwin):
                w = next(st)
                tr = w["truth"]
                row = {"cell": cell, "net": j, "w": wi, "k": len(tr), "rec": {}, "sec": {}}
                for name, fn in (("MT-FW", lambda: L.mtfw_dense(ops, w["MFn"], tr, C2.KMAX)[0]),
                                 ("MT-FW-NOCAP", lambda: L.mtfw_dense(ops, w["MFn"], tr, C2.KMAX, cap=False)[0]),
                                 ("BOMP-PK", lambda: L.bomp_pk_m(ops, w["MFn"], tr, C2.KMAX)[0]),
                                 ("BOMP-PK-CAP", lambda: [(float(s), int(h)) for s, h in
                                                          D2.bomp_pk_record(PDc, w["MFn"], tr, kmax=C2.KMAX)]),
                                 ("BOMP-PK-SEQ", lambda: bomp_seq(ops.PD, w["MFn"], tr, C2.KMAX))):
                    t = time.time()
                    row["rec"][name] = fn()
                    row["sec"][name] = time.time() - t
                res.append(row)
            print(cell, j, "done", time.strftime("%H:%M:%S"), flush=True)
            with open(out, "w") as fh:
                json.dump(res, fh)


if __name__ == "__main__":
    run(sys.argv[1], [int(x) for x in sys.argv[2].split(",")], int(sys.argv[3]), sys.argv[4])
