"""C15 diagnosis of the origin-time error of located atoms (truth used; development networks 1800-1802 only).
For every BASE atom within 5 km of a true event (nearest), with the read samples tp of the detector:
  eB     reported origin minus true origin
  eRT    origin implied by the reads with the TRUE location, TRUE travel times (gradient medium) and TRUE statics
         (weighted median of read time - tt - static, weights = read energy - 3): reading error only
  eRE    origin implied by the reads with the TRUE location and the ENGINE model (homogeneous rays + registered +SC):
         reading error + model error at the true location
  eAE    origin implied by the reads at the ATOM node with the engine model (= LSO): adds the location error
  fsig   fraction of the station-phases with weight > 0 whose read lies within 1 sample of the true arrival
"""
import json
import sys

import numpy as np

import c2_run as C2
import c10_run as R
import c14_lib as L14
import c15_lib as L
import gen_realistic as GR
import xc
from task0_common import DT, VP, VS

out = []
for j in (1800, 1801, 1802):
    for regime in ("strong", "thin"):
        for cell in R.CELLS:
            S = R.Setup(cell, j, regime, oracle=False)
            E = S.eng["cat"][0]
            nc = S.nc
            corr = np.asarray(nc.corr, float)                      # (ns, 2)
            for design, (tag, draw, kmax) in {"STD": (1851, xc.draw_sources, 5),
                                              "DENSE": (1852, C2.draw_dense, 7)}.items():
                st = L14.stream(nc, tag, draw)
                for wi in range(20):
                    w = next(st)
                    T = np.array(w["truth"])
                    if not len(T):
                        continue
                    for a in L.mt_loop15(E, w["MFn"], kmax):
                        d = np.linalg.norm(T[:, :2] - a["x"][None, :2], axis=1)
                        i = int(d.argmin())
                        if d[i] > 5.0:
                            continue
                        src, tau = T[i], w["taus"][i]
                        wt = L.weights(a)                          # (2, ns)
                        if wt.sum() <= 0:
                            continue
                        tread = (a["tp"] + a["sub"]) * DT          # (2, ns)
                        ttT = np.stack([GR.travel_time_gradient(src, nc.STA, "P"),
                                        GR.travel_time_gradient(src, nc.STA, "S")])
                        tarr = tau + ttT + np.asarray(nc.statics, float).T
                        dh = np.linalg.norm(nc.STA - src[None], axis=1)
                        ttE = np.stack([dh / VP, dh / VS]) + corr.T
                        eRT = L._wmedian((tread - ttT - np.asarray(nc.statics, float).T).ravel(), wt.ravel()) - tau
                        eRE = L._wmedian((tread - ttE).ravel(), wt.ravel()) - tau
                        eAE = L.est_lso(E, a, True)[1] - tau
                        m = wt > 0
                        fsig = float(np.mean(np.abs(tread[m] - tarr[m]) <= DT))
                        out.append(dict(design=design, cell=cell, regime=regime, s=a["s"], eB=a["tau"] - tau,
                                        eRT=eRT, eRE=eRE, eAE=eAE, fsig=fsig, dxy=float(d[i]),
                                        dz=float(a["x"][2] - src[2]), nsig=int(m.sum()),
                                        rd=(tread[m] - tarr[m]).tolist()))
json.dump(out, open("out/c15_diag.json", "w"))
for design in ("STD", "DENSE"):
    o = [r for r in out if r["design"] == design]
    print(design, "n", len(o))
    for q in ("all", "top"):
        oo = o if q == "all" else [r for r in o if r["s"] >= np.percentile([x["s"] for x in o], 50)]
        for k in ("eB", "eRT", "eRE", "eAE"):
            e = np.array([r[k] for r in oo])
            print(f"  {q:3s} {k:3s} median {np.median(e):+.3f}  |e| median {np.median(np.abs(e)):.3f}  "
                  f"P(|e|<=0.5) {np.mean(np.abs(e) <= 0.5):.3f}")
        rd = np.concatenate([r["rd"] for r in oo])
        print(f"  {q:3s} read minus true arrival (weighted reads): median {np.median(rd):+.3f}, |.| median "
              f"{np.median(np.abs(rd)):.3f}, within 1 sample {np.mean(np.abs(rd) <= DT):.3f}, within 0.25 s "
              f"{np.mean(np.abs(rd) <= 0.25):.3f}")
    sys.stdout.flush()
