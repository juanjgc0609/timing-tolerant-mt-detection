"""
r0_run.py - Task R0: pick-free detectors over the continuous data of one subnetwork.
    python r0_run.py SEGMENT SUBNET METHODS OUT.json        SEGMENT = cal | eval;  METHODS = comma list of M1..M5
Each 40 s window starting at t_start (step 3 s) owns the atoms with origin in [5.0, 8.0) s of window time; every
other atom is explained and deflated but not reported.  Output: per method, the owned atoms of every window
(absolute origin time in s after 19:00:00 UTC), plus the calibration products (corrections, scales, frame).
"""
import json
import sys
import time
import types

import numpy as np

import r0_lib as R                                   # sets VP, VS and the origin range first

import c10_lib as C                                  # noqa: E402
import c11_lib as L11                                # noqa: E402
import c8_lib as C8                                  # noqa: E402
import task1x_mtglrt as X                            # noqa: E402
import xc                                            # noqa: E402
from certified_fw import stack_Q                     # noqa: E402
from station_corrections import estimate_station_corrections  # noqa: E402

METHODS = {
    "M1": "MT-FW cap 6 +SC",
    "M2": "MT-FW cap 6 + S2 + T3",
    "M3": "BOMP-PK +SC",
    "M4": "HUBER-N +SC",
    "M5": "MT-FW cap 6 +SC +T3",
    "M6": "MT-FW cap 6 + S2 + T3, read half-width W = round(2 sigma_S2 / DT) (amendment, secondary)",
}


def grids():
    G = np.asarray(xc.make_grid(depths=R.DEPTHS))
    ax = np.arange(-21.0, 21.0 + 1e-9, C.H_FINE) + 0.37
    F = np.array([[x, y, z] for x in ax for y in ax for z in R.DEPTHS])
    return G, F


class Net:
    def __init__(self, subnet, t_end):
        t0 = time.time()
        self.names = R.SUBNETS[subnet]
        self.centre = R.frame()
        self.cat = R.catalogue_local(self.centre)
        self.STA, sens = R.stations(self.names, self.centre)
        self.raw, d16 = R.load_waveforms(self.names, sens, t_end)
        self.scale = R.station_scale(d16)
        self.d16 = d16 / self.scale[:, None, None]
        # calibration catalogue: CAL events with M >= 2.0, windows with the origin at 6.5 s
        self.calev = [e for e in self.cat if R.T_CAL[0] <= e["t"] < R.T_CAL[1] and e["M"] >= R.M_TARGET]
        self.calwin = []
        for e in self.calev:
            ts = e["t"] - 6.5
            if ts < 0:
                continue
            self.calwin.append((R.window(self.d16, ts), np.array([e["x"], e["y"], e["z"]]), 6.5))
        self.corr_sc, self.nacc = estimate_station_corrections(self.STA, self.calwin)
        self.s2_info = {}
        self.corr_s2 = L11.s2_terms(self.STA, self.calwin, self.corr_sc, None, L11.P, info=self.s2_info)
        G, F = grids()
        self.G, self.F = G, F
        Gs = np.stack([X.radiation_operator(g, self.STA) for g in G])
        Qs = [X.basis(g) for g in Gs]
        self.Qst = stack_Q(Qs)
        self.Gs = Gs
        self.eng = {}
        for nm, corr in (("SC", self.corr_sc), ("S2", self.corr_s2)):
            sig, _ = C.sigma_catalogue(self.STA, self.calwin, corr)
            E = C.Eng(G, self.STA, corr, R.TAUS, C.H_COARSE, G=Gs, Q=self.Qst, bomp=True)
            EF = C.Eng(F, self.STA, corr, R.TAUS, C.H_FINE, bomp=True)
            self.eng[nm] = (E, E.prior(sig), EF, EF.prior(sig), sig)
        E = self.eng["SC"][0]
        self.dops = types.SimpleNamespace(ns=len(self.STA), garr=G, IDX=E.IDX, Gs=Gs, Qst=self.Qst,
                                          Gp=np.linalg.pinv(Gs), taus=np.asarray(R.TAUS, float))
        self.t_setup = time.time() - t0

    def run(self, m, MFn):
        k = R.KMAX
        if m == "M1":
            E = self.eng["SC"][0]
            return [(a["x"], a["tau"], a["s"]) for a in C.mt_loop(E, MFn, k)]
        if m == "M2":
            E, _p, EF, prF, _s = self.eng["S2"]
            return [(a["x"], a["tau"], a["s"]) for a in C.mt_loop(E, MFn, k, EF=EF, prF=prF)]
        if m == "M6":
            E, _p, EF, prF, sig = self.eng["S2"]
            W = int(round(2.0 * float(np.max(sig)) / R.DT))
            return [(a["x"], a["tau"], a["s"]) for a in C.mt_loop(E, MFn, k, W=W, EF=EF, prF=prF)]
        if m == "M5":
            E, _p, EF, prF, _s = self.eng["SC"]
            return [(a["x"], a["tau"], a["s"]) for a in C.mt_loop(E, MFn, k, EF=EF, prF=prF)]
        if m == "M3":
            E = self.eng["SC"][0]
            return [(a["x"], a["tau"], a["s"]) for a in C.bomp_loop(E, None, MFn, k)]
        if m == "M4":
            rec, locs = C8.mtfw_general(self.dops, MFn, [np.array([1e6, 1e6, 0.0])], kmax=k, stat="huber_n")
            return [(np.array(l[:3]), float(l[3]), float(r[0])) for r, l in zip(rec, locs)]
        raise ValueError(m)


def starts(segment):
    lo, hi = R.T_CAL if segment == "cal" else R.T_EVAL
    first = max(0.0, lo - R.OWN[0])
    out = []
    t = first
    while t + R.OWN[0] < hi:
        out.append(round(t, 6))
        t += R.STEP
    return out


def main(segment, subnet, methods, path):
    t_end = (R.T_CAL[1] if segment == "cal" else R.T_EVAL[1]) + 40.0
    net = Net(subnet, t_end)
    res = {"segment": segment, "subnet": subnet, "names": net.names, "centre": net.centre,
           "STA": net.STA.tolist(), "scale": net.scale.tolist(), "corr_sc": net.corr_sc.tolist(),
           "nacc": np.asarray(net.nacc).tolist(), "corr_s2": net.corr_s2.tolist(),
           "s2_info": {k: (v if not isinstance(v, np.ndarray) else v.tolist()) for k, v in net.s2_info.items()
                       if k != "hist"},
           "sig": {k: np.asarray(v[4]).tolist() for k, v in net.eng.items()}, "n_calwin": len(net.calwin),
           "t_setup": net.t_setup, "atoms": {m: [] for m in methods}, "sec": {m: 0.0 for m in methods}}
    ws = starts(segment)
    for wi, ts in enumerate(ws):
        W = R.window(net.d16, ts)
        MF, sd = xc.matched(W)
        MFn = MF / sd
        for m in methods:
            t0 = time.time()
            for x, tau, s in net.run(m, MFn):
                if R.OWN[0] - 1e-9 <= tau < R.OWN[1] - 1e-9:
                    res["atoms"][m].append([float(x[0]), float(x[1]), float(x[2]), float(ts + tau), float(s), wi])
            res["sec"][m] += time.time() - t0
        if wi % 50 == 0:
            print(time.strftime("%H:%M:%S"), subnet, segment, wi, len(ws), {m: round(v, 1) for m, v in
                                                                         res["sec"].items()}, flush=True)
    res["n_windows"] = len(ws)
    with open(path + ".tmp", "w") as fh:
        json.dump(res, fh)
    import os
    os.replace(path + ".tmp", path)
    print("DONE", time.strftime("%H:%M:%S"), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3].split(","), sys.argv[4])
