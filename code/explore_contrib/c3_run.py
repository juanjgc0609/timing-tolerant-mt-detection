"""C3: large timing errors. Conditions: station statics sd 0.3/0.5/0.8 s (fixed per network) with and without
catalogue station corrections, and a velocity-model error (true medium 7.5% slower). Methods:
  BOMP-PK (+-0.25 s), BOMP-PK-8 (+-0.5 s), BOMP-PK-12 (+-0.75 s);
  BOMP-SS: shared shifts within a window (pass 1 wide +-0.75 s, per station/phase median shift of the confident atoms
           (score >= S_CONF), pass 2 narrow +-0.25 s with those shifts added to the corrections);
  BOMP-OL: online station terms: per station/phase median of delays measured on confident detections (score >= S_CONF,
           wide +-0.75 s BOMP, joint P/S peak search +-1 s around the corrected prediction) of all previous windows of the
           network stream, seeded with the catalogue delays when a catalogue exists; detection with +-0.25 s.
Stream: NB burn-in windows (learning only) then NE evaluation windows per network (TAG 3003, same sources for all
conditions)."""
import json, sys, time
import numpy as np
import xc
from task2c_methods import BlockDict
from station_corrections import MIN_PS, peak_delays

CELLS = ("ladder_4", "ladder_6", "unilateral_c40")
TAG = 3003
S_CONF = 60.0          # fixed a priori (thin 2A-e thresholds were 35-43)
VEL = {"VS0": 2.8 / 1.075, "g": 0.04 / 1.075}
CONDS = {"S03sc": (0.3, None, True), "S05sc": (0.5, None, True), "S08sc": (0.8, None, True),
         "S05nc": (0.5, None, False), "S08nc": (0.8, None, False), "V75sc": (0.3, VEL, True),
         "V75nc": (0.3, VEL, False)}


def delays_at(STA, MFn, src, tau, center_corr, search=1.0):
    """absolute delays (peak time - homogeneous prediction) per (station, phase), joint P/S search +-search around
    the corrected prediction; parabolic refinement (station_corrections.peak_delays logic). Returns delays, energies."""
    Et = np.sum(MFn ** 2, axis=1)
    d = np.linalg.norm(STA - src, axis=1)
    W = int(round(search / xc.DT))
    NT = MFn.shape[-1]
    out = np.zeros((len(STA), 2))
    en = np.zeros((len(STA), 2))
    for i in range(len(STA)):
        tp = [d[i] / xc.VP + tau, d[i] / xc.VS + tau]
        win = []
        for b in range(2):
            c = int(round((tp[b] + center_corr[i, b]) / xc.DT))
            win.append(np.arange(max(1, c - W), min(NT - 2, c + W) + 1))
        sc = Et[i, win[0]][:, None] + Et[i, win[1]][None, :]
        sc = np.where(win[1][None, :] >= win[0][:, None] + MIN_PS, sc, -np.inf)
        a, b_ = np.unravel_index(int(np.argmax(sc)), sc.shape)
        ks = (int(win[0][a]), int(win[1][b_]))
        for b in range(2):
            k = ks[b]
            e0, e1, e2 = Et[i, k - 1], Et[i, k], Et[i, k + 1]
            den = e0 - 2 * e1 + e2
            dk = 0.5 * (e0 - e2) / den if den < 0 else 0.0
            out[i, b] = (k + float(np.clip(dk, -0.5, 0.5))) * xc.DT - tp[b]
            en[i, b] = e1
    return out, en


class Online:
    def __init__(self, nc, sc):
        self.nc = nc
        ns = len(nc.STA)
        self.obs = [[[] for _ in range(2)] for _ in range(ns)]
        if sc:
            for D, src, tau in nc.cat:
                dl, en = peak_delays(nc.STA, D, np.asarray(src, float), float(tau))
                for i in range(ns):
                    for b in range(2):
                        if en[i, b] >= 25.0:
                            self.obs[i][b].append(float(dl[i, b]))
        self.n_learn = 0

    def corr(self):
        ns = len(self.nc.STA)
        c = np.zeros((ns, 2))
        for i in range(ns):
            for b in range(2):
                if len(self.obs[i][b]) >= 3:
                    c[i, b] = float(np.median(self.obs[i][b]))
        return c

    def learn(self, STA, MFn, atoms_conf, corr):
        for loc in atoms_conf:
            dl, en = delays_at(STA, MFn, np.asarray(loc[:3]), loc[3], corr)
            for i in range(len(STA)):
                for b in range(2):
                    if en[i, b] >= 25.0:
                        self.obs[i][b].append(float(dl[i, b]))
            self.n_learn += 1


def run(regime, nets, nb, ne, conds, out):
    res = []
    for cname in conds:
        ssd, cfg, sc = CONDS[cname]
        for cell in CELLS:
            for j in nets:
                nc = xc.NetCfg(cell, j, regime, static_sd=ssd, cfg=cfg, sc=sc)
                st = xc.window_stream(nc, TAG)
                ops = xc.Ops(nc.STA, nc.corr, xc.make_grid())
                P8, P12 = xc.PKDictW(ops.BD, 8), xc.PKDictW(ops.BD, 12)
                ol = Online(nc, sc)
                for wi in range(nb + ne):
                    w = next(st)
                    tr = w["truth"]
                    # ---- online: detection with current learned corrections, then learning from this window
                    t0 = time.time()
                    c_ol = ol.corr()
                    BDo = BlockDict(nc.STA, ops.grid, xc.DEFAULT_TAUS, c_ol, Gs=ops.Gs)
                    PDo = xc.D2.PKDict(BDo)
                    ev = wi >= nb
                    if ev:
                        r_ol, _l = xc.bomp_pk(PDo, w["MFn"], tr)
                    PDo12 = xc.PKDictW(BDo, 12)
                    rw, lw = xc.bomp_pk(PDo12, w["MFn"], tr, kmax=3)
                    conf = [l for (s, h), l in zip(rw, lw) if s >= S_CONF]
                    ol.learn(nc.STA, w["MFn"], conf, c_ol)
                    t_ol = time.time() - t0
                    if not ev:
                        continue
                    osc = xc.oracle_scores(ops, w)
                    row = {"cond": cname, "cell": cell, "net": j, "w": wi, "k": len(tr), "truth": [t.tolist() for t in tr],
                           "oracle": osc, "rec": {}, "loc": {}, "sec": {"BOMP-OL": t_ol}, "n_learn": ol.n_learn,
                           "corr_err_ol": float(np.median(np.abs(c_ol - nc.statics))),
                           "corr_err_cat": float(np.median(np.abs(nc.corr - nc.statics)))}
                    row["rec"]["BOMP-OL"] = r_ol
                    for name, PD in (("BOMP-PK", ops.PD), ("BOMP-PK-8", P8), ("BOMP-PK-12", P12)):
                        t = time.time(); r, l = xc.bomp_pk(PD, w["MFn"], tr); row["sec"][name] = time.time() - t
                        row["rec"][name], row["loc"][name] = r, l
                        if name == "BOMP-PK":
                            A, Bc = xc.ceilings(r, osc, tr, l)
                            row["rec"]["CEIL-A"], row["rec"]["CEIL-B"] = A, Bc
                        if name == "BOMP-PK-12":
                            r12, dg12 = r, None
                    # ---- shared shifts within the window
                    t = time.time()
                    dg = {}
                    r1 = xc.D2.bomp_pk_record(P12, w["MFn"], tr, diag=dg)
                    sh = [np.asarray(s_) for (s, _h), s_ in zip(r1, dg["shifts"]) if s >= S_CONF]
                    if sh:
                        shm = np.median(np.stack(sh), axis=0).reshape(2, -1).T * xc.DT     # (ns, 2)
                        BDs = BlockDict(nc.STA, ops.grid, xc.DEFAULT_TAUS, nc.corr + shm, Gs=ops.Gs)
                        r, _l = xc.bomp_pk(xc.D2.PKDict(BDs), w["MFn"], tr)
                    else:
                        r = row["rec"]["BOMP-PK"]
                    row["sec"]["BOMP-SS"] = time.time() - t
                    row["rec"]["BOMP-SS"] = r
                    res.append(row)
                print(cname, cell, j, "done", time.strftime("%H:%M:%S"), "learned", ol.n_learn, flush=True)
                with open(out, "w") as fh:
                    json.dump(res, fh)


if __name__ == "__main__":
    run(sys.argv[1], [int(x) for x in sys.argv[2].split(",")], int(sys.argv[3]), int(sys.argv[4]),
        sys.argv[5].split(","), sys.argv[6])
