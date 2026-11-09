"""M2-GP: gate G-M2 of the v2 specification (M2, station terms in a Gaussian process; plan T4.8, build order after M6).

In validation by blocks of time, does the GP term (``ttmt.v2.m2``, Matern 3/2 in (x, y, z), hyperparameters by
marginal likelihood on the training blocks only) lower the mean squared error of the held-out calibration readings by
>= 30 % against the constant term, with the coverage of +-2 sigma between 90 % and 98 %?

Data: SCSN analyst P and S picks of the Ridgecrest 2019 sequence, 2019-07-09 to 07-15 (``tools/fetch_scedc_phases.py``;
SHA-256 of the tables checked below), events M >= 1.5 in the box of the fetch tool, stations with >= 300 P and >= 300 S
picks (rule fixed before the run). Local flat-earth frame around 35.80 N, 117.60 W; stations at z = 0 (elevation goes
into the term), events at the catalogue hypocentre. 1-D medium: Hadley and Kanamori (``ttmt.traveltime.HK``),
VP/VS 1.73. Residual = observed travel time - 1-D time; gross errors |residual| > 3 s are removed before anything
(the same rows for every method; their count is reported).

Folds: 5 contiguous blocks of events in origin time (``ttmt.v2.m2.time_blocks``); each block is predicted by terms
fitted on the other four. Methods: none (1-D only), const (mean of the training residuals, sd as tolerance), gp.

    python run.py DATA_DIR [plain|het] [--no-sha] > run_out.txt
Variant "het" (second step, card.md): the GP's white level depends on the SCSN pick weight class (ttmt.v2.m2.GP1 with
groups); "plain" (first step): one white level. --no-sha: the confirmation days (197-203), whose tables have their
own SHA-256 printed and recorded in the card.
"""
import hashlib
import os
import sys
import time

import numpy as np
import pandas as pd

from ttmt.v2.m2 import Const1, GP1, time_blocks
from ttmt.traveltime import HK, R_VPVS, SourceTimes

SHA = {"events.csv": "9e220d4d0dcb4bea33ce00de8fc82c33f7d840bac8f37a06f27ca544ffc8f931",
       "picks.csv": "40781188fdc53338b7480445d98f8ef7d45aaf7a45c950b3f016d66b1e57cdef"}
LAT0, LON0 = 35.80, -117.60
KM_LAT = 111.19
NMIN = 300
GROSS = 3.0
K = 5


def local(lat, lon):
    return (np.asarray(lon) - LON0) * KM_LAT * np.cos(np.radians(LAT0)), (np.asarray(lat) - LAT0) * KM_LAT


def main(d, variant="plain", check=True):
    t0 = time.time()
    print(f"variant {variant}")
    for n, h in SHA.items():
        got = hashlib.sha256(open(os.path.join(d, n), "rb").read()).hexdigest()
        print(f"{n}: sha256 {got if not check else ('OK' if got == h else 'MISMATCH ' + got)}")
    ev = pd.read_csv(os.path.join(d, "events.csv"))
    pk = pd.read_csv(os.path.join(d, "picks.csv"))
    pk["code"] = pk.net + "." + pk.sta
    cnt = pk.groupby(["code", "phase"]).size().unstack().fillna(0)
    codes = sorted(cnt.index[(cnt.P >= NMIN) & (cnt.S >= NMIN)])
    pk = pk[pk.code.isin(codes)].copy()
    st = pk.groupby("code")[["lat", "lon"]].first()
    sx, sy = local(st.lat.values, st.lon.values)
    STA = np.stack([sx, sy, np.zeros(len(st))], 1)
    sidx = {c: i for i, c in enumerate(st.index)}
    ex, ey = local(ev.lat.values, ev.lon.values)
    EX = np.stack([ex, ey, ev.depth_km.values], 1)
    eidx = {e: i for i, e in enumerate(ev.evid)}
    org = pd.to_datetime(ev.origin).astype("int64").values / 1e9
    fold = time_blocks(org, K)
    tt = SourceTimes(HK, R_VPVS)
    T = np.array([tt(EX[i], STA) for i in range(len(EX))])           # (ne, 2, ns)
    e = pk.evid.map(eidx).values
    s = pk.code.map(sidx).values
    b = (pk.phase == "S").values.astype(int)
    pk["res"] = pk.tt_s.values - T[e, b, s]
    pk["e"], pk["s"], pk["b"], pk["fold"] = e, s, b, fold[e]
    gross = np.abs(pk.res.values) > GROSS
    print(f"{len(ev)} events (M >= 1.5), {len(codes)} stations with >= {NMIN} P and S picks, {len(pk)} picks; "
          f"gross |res| > {GROSS} s removed: {gross.sum()} ({100 * gross.mean():.2f} %)")
    pk = pk[~gross]
    print(f"folds (time blocks): {np.bincount(fold).tolist()} events; 1-D medium Hadley-Kanamori, VP/VS {R_VPVS}\n")

    out = []
    for (code, ph), g in pk.groupby(["code", "phase"]):
        X, y, f = EX[g.e.values], g.res.values, g.fold.values
        for k in range(K):
            tr, te = f != k, f == k
            if te.sum() == 0 or tr.sum() < 30:
                continue
            w = g.weight.values if variant == "het" else None
            c = Const1().fit(X[tr], y[tr])
            gp = GP1(seed=k).fit(X[tr], y[tr], None if w is None else w[tr])
            mc, sc = c.predict(X[te])
            mg, sg = gp.predict(X[te], None if w is None else w[te])
            nz = gp._noise(gp.theta, gp._gi(None if w is None else w[te], int(te.sum())))
            out.append(pd.DataFrame({"code": code, "phase": ph, "fold": k, "w": g.weight.values[te], "res": y[te],
                                     "ec": y[te] - mc, "sc": sc, "eg": y[te] - mg, "sg": sg,
                                     "ell": gp.ell, "sf": gp.sf, "sn": gp.sn, "white": nz}))
    R = pd.concat(out, ignore_index=True)

    def summ(D):
        mse0, msec, mseg = (D.res ** 2).mean(), (D.ec ** 2).mean(), (D.eg ** 2).mean()
        return {"n": len(D), "rmse 1-D": np.sqrt(mse0), "rmse const": np.sqrt(msec), "rmse gp": np.sqrt(mseg),
                "MSE drop": 1 - mseg / msec, "MAD const": 1.4826 * np.median(np.abs(D.ec - np.median(D.ec))),
                "MAD gp": 1.4826 * np.median(np.abs(D.eg - np.median(D.eg))),
                "cov2 const": (np.abs(D.ec) <= 2 * D.sc).mean(), "cov2 gp": (np.abs(D.eg) <= 2 * D.sg).mean(),
                "med sigma gp": D.sg.median(), "drop bound": 1 - D.white.mean() / msec}

    pd.set_option("display.width", 200)
    rows = {"all": summ(R), "P": summ(R[R.phase == "P"]), "S": summ(R[R.phase == "S"])}
    for w in sorted(R.w.unique()):
        rows[f"weight {w}"] = summ(R[R.w == w])
    for k in range(K):
        rows[f"fold {k}"] = summ(R[R.fold == k])
    print("held-out readings (every fold predicted by terms fitted on the other four)")
    print(pd.DataFrame(rows).T.round(3).to_string())
    per = pd.DataFrame({key: summ(D) for key, D in R.groupby(["code", "phase"])}).T
    hp = R.groupby(["code", "phase"])[["ell", "sf", "sn"]].median()
    per.index = hp.index
    print("\nper station and phase (median hyperparameters over folds: ell km, sf s, sn s)")
    print(per.join(hp)[["n", "rmse const", "rmse gp", "MSE drop", "drop bound", "cov2 gp", "ell", "sf", "sn"]]
          .round(3).to_string())
    print("\ndrop bound = 1 - (mean white variance of the GP at the held-out readings) / MSE(const): the largest MSE "
          "drop any term can reach if the white part is irreducible reading noise (diagnosis, not a gate)")
    a = rows["all"]
    ok = a["MSE drop"] >= 0.30 and 0.90 <= a["cov2 gp"] <= 0.98
    print(f"\nG-M2: MSE drop {a['MSE drop']:.3f} (>= 0.30), coverage +-2 sigma {a['cov2 gp']:.3f} (0.90-0.98) -> "
          f"{'PASS' if ok else 'FAIL'}")
    print(f"\n{time.time() - t0:.0f} s, numpy {np.__version__}")


if __name__ == "__main__":
    a = [x for x in sys.argv[1:] if not x.startswith("--")]
    sys.exit(main(a[0], a[1] if len(a) > 1 else "plain", "--no-sha" not in sys.argv))
