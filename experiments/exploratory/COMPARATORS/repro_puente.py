"""Reproduction gate G3a, target 1 of the ATLAS-0 card: GaMMA, PyOcto and REAL on the crustal synthetic scenario of
Puente et al. (2025, Seismica 4(2), doi 10.26443/seismica.v4i2.1559) with the parameters of their Supplementary
Tables S3-S5, against the event-level recalls of their Table S10 (the card names S8-S9; in the supplement the crustal
event-level results are Table S10, S8-S9 being the subduction scenario).

Scenario as the paper describes it (Section 2): SCSN station layout, Hadley-Kanamori 1-D model, events uniform in
x [385, 520] km, y [3860, 4040] km (UTM 11N, EPSG 32611), depth 0-30 km, origin times uniform over 24 h; arrival
noise uniform within +-1 % of the travel time; per event a distance cutoff uniform in 70-150 km and 20 % of the
remaining stations dropped; false picks 30, 100 or 300 % of the true picks; 100, 500 and 2000 events. Match
(Figs. S7-S8 "50 % matching threshold"): a predicted event matches the first true event with which it shares at least
half of the larger of the two pick sets (the rule of PyOcto's benchmark, benchmark/eval.py, at 0.5 instead of 0.6).

Deviations, declared: travel times from ttmt.traveltime (flat layered first arrivals) instead of NonLinLoc; the SCSN
layout is the list given with --stations (the paper does not list its stations or date); events with fewer than 10
picks are not generated (no associator can report them at the paper's threshold of 10 picks); GaMMA runs with ncpu = 1
(the paper used 25), which is this project's fixed setting.

    python repro_puente.py fetch-stations stations.csv            # SCEDC FDSN, CI broadband/high-gain, box + 50 km
    python repro_puente.py run --stations stations.csv --out out.json [--events 100 500 2000] [--noise .3 1 3]
                                [--methods gamma pyocto real] [--seed 0]
    python repro_puente.py report out.json                         # recall against Table S10, pass within 5 %
"""
import argparse
import json
import sys
import time

import numpy as np

from ttmt.associate.base import Frame
from ttmt.associate.velocity import Velocity1D
from ttmt.atlas.schema import Picks, Stations

XLIM, YLIM, ZLIM = (385.0, 520.0), (3860.0, 4040.0), (0.0, 30.0)
DAY = 86400.0
PAD_KM = 50.0

# Table S10 (crustal, event level): (n_events, noise fraction) -> (precision, recall, F1); None: not reported
S10 = {
    "gamma": {(100, 0.3): (1.00, 1.00, 1.00), (100, 1.0): (1.00, 0.99, 0.99), (100, 3.0): (0.99, 0.99, 0.99),
              (500, 0.3): (0.99, 0.99, 0.99), (500, 1.0): (1.00, 1.00, 1.00), (500, 3.0): (0.79, 0.99, 0.88),
              (2000, 0.3): (0.98, 0.99, 0.98), (2000, 1.0): (0.57, 0.96, 0.71), (2000, 3.0): None},
    "pyocto": {(100, 0.3): (1.00, 1.00, 1.00), (100, 1.0): (1.00, 0.98, 0.99), (100, 3.0): (1.00, 0.99, 0.99),
               (500, 0.3): (1.00, 0.97, 0.98), (500, 1.0): (1.00, 0.98, 0.99), (500, 3.0): (1.00, 0.98, 0.99),
               (2000, 0.3): (0.99, 0.93, 0.96), (2000, 1.0): (0.99, 0.92, 0.95), (2000, 3.0): (0.73, 0.91, 0.81)},
    "real": {(100, 0.3): (1.00, 0.98, 0.99), (100, 1.0): (1.00, 0.98, 0.99), (100, 3.0): (1.00, 0.99, 0.99),
             (500, 0.3): (1.00, 0.93, 0.96), (500, 1.0): (0.99, 0.93, 0.96), (500, 3.0): (0.98, 0.93, 0.96),
             (2000, 0.3): (0.99, 0.78, 0.87), (2000, 1.0): (0.96, 0.79, 0.87), (2000, 3.0): (0.91, 0.73, 0.81)},
}

# Supplementary Tables S3 (GaMMA 1-D), S4 (PyOcto 1-D), S5 (REAL 0-D), crustal scenario
CFG = {"gamma": dict(as_min_picks=10, as_max_sigma11=2.0, as_method="BGMM"),
       "pyocto": dict(as_min_picks=10, as_tolerance=1.0, as_pick_match_tolerance=0.8),
       "real": dict(as_np0=4, as_nps0=10, as_tint=0.1)}
FIXED = {"gamma": dict(vel={"p": 6.2, "s": 3.4}, dbscan_eps=7.0, dbscan_min_samples=20, oversample_factor=3,
                       max_sigma22=1.0, max_sigma12=1.0, covariance_prior=None, initial_points=None, one_d=True),
         "pyocto": dict(one_d=True, table_delta=1.0, table_xdist=300.0, association_cutoff_distance=200.0,
                        min_node_size=10.0, min_node_size_location=2.5, min_interevent_time=3.0, max_pick_overlap=4,
                        refinement_iterations=3, time_slicing=1200.0, min_pick_fraction=0.0, n_p_picks=5,
                        n_s_picks=5, n_p_and_s_picks=4),
         "real": dict(one_d=False, rx_km=1.0 * 111.19, tdx_km=0.1 * 111.19, tdh_km=8.0, V=(6.2, 3.3, 5.4, 3.3),
                      ns0=4, npsboth0=4, std0=2.0, dtps=0.0, nrt=2.0, drt=2.0, nxd=0.4, rsel=6.0, timeout=4 * 3600)}
TIME_BEFORE = 100.0


def fetch_stations(out):
    from obspy.clients.fdsn import Client
    from pyproj import Transformer
    tr = Transformer.from_crs("EPSG:32611", "EPSG:4326", always_xy=True)
    lon0, lat0 = tr.transform((XLIM[0] - PAD_KM) * 1e3, (YLIM[0] - PAD_KM) * 1e3)
    lon1, lat1 = tr.transform((XLIM[1] + PAD_KM) * 1e3, (YLIM[1] + PAD_KM) * 1e3)
    inv = Client("SCEDC").get_stations(network="CI", channel="HH?,BH?", level="station", minlatitude=lat0,
                                       maxlatitude=lat1, minlongitude=lon0, maxlongitude=lon1,
                                       starttime="2020-01-01", endtime="2020-01-02")
    rows = sorted({(s.code, s.latitude, s.longitude, s.elevation) for n in inv for s in n})
    with open(out, "w") as fh:
        fh.write("code,latitude,longitude,elevation_m\n")
        for r in rows:
            fh.write("{},{:.5f},{:.5f},{:.1f}\n".format(*r))
    src = dict(service="SCEDC FDSN", network="CI", channel="HH?,BH?", box_lonlat=[lon0, lat0, lon1, lat1],
               pad_km=PAD_KM, date="2020-01-01", n_stations=len(rows), sha256=file_sha256(out))
    with open(out + ".source.json", "w") as fh:
        json.dump(src, fh, indent=1)
    print(len(rows), "stations ->", out)


def file_sha256(path):
    import hashlib
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def station_source(path):
    """the provenance written by fetch-stations next to a station file, None if absent or if the file changed since."""
    import os
    side = str(path) + ".source.json"
    if not os.path.exists(side):
        return None
    src = json.load(open(side))
    return src if src.get("sha256") == file_sha256(path) else None


def load_stations(path):
    """Stations in UTM 11N km (z = -elevation km) from a CSV with code, latitude, longitude, elevation_m (or the
    GaMMA test layout: station, longitude, latitude, elevation(m), tab separated)."""
    import pandas as pd
    from pyproj import Transformer
    df = pd.read_csv(path, sep=None, engine="python")
    if "station" in df.columns:
        df = df.rename(columns={"station": "code", "elevation(m)": "elevation_m"})
        df["code"] = df["code"].str.split(".").str[1]
    df = df.drop_duplicates("code").sort_values("code")
    tr = Transformer.from_crs("EPSG:4326", "EPSG:32611", always_xy=True)
    x, y = tr.transform(df["longitude"].to_numpy(), df["latitude"].to_numpy())
    xyz = np.c_[np.asarray(x) / 1e3, np.asarray(y) / 1e3, -df["elevation_m"].to_numpy() / 1e3]
    keep = ((xyz[:, 0] > XLIM[0] - PAD_KM) & (xyz[:, 0] < XLIM[1] + PAD_KM) & (xyz[:, 1] > YLIM[0] - PAD_KM)
            & (xyz[:, 1] < YLIM[1] + PAD_KM))
    return Stations(df["code"].to_numpy()[keep], xyz[keep])


def generate(stations, n_events, noise, rng, vel):
    """(Picks, event of each pick (-1 false), events (n, 4) x y z t)."""
    sta = stations.xyz.copy()
    sta[:, 2] = 0.0
    ev, st, ph, tt, owner = [], [], [], [], []
    while len(ev) < n_events:
        x = np.array([rng.uniform(*XLIM), rng.uniform(*YLIM), rng.uniform(*ZLIM)])
        t0 = rng.uniform(0.0, DAY)
        d = np.hypot(*(sta[:, :2] - x[:2]).T)
        near = np.flatnonzero(d <= rng.uniform(70.0, 150.0))
        keep = near[rng.permutation(len(near))[: int(round(0.8 * len(near)))]]
        if 2 * len(keep) < 10:
            continue
        T = vel.times(x, sta[keep])
        e = len(ev)
        ev.append((*x, t0))
        for p in (0, 1):
            arr = t0 + T[p] * (1.0 + rng.uniform(-0.01, 0.01, len(keep)))
            st += keep.tolist()
            ph += [p] * len(keep)
            tt += arr.tolist()
            owner += [e] * len(keep)
    nf = int(round(noise * len(tt)))
    st += rng.integers(0, len(sta), nf).tolist()
    ph += rng.integers(0, 2, nf).tolist()
    tt += rng.uniform(0.0, DAY, nf).tolist()
    owner += [-1] * nf
    o = np.argsort(tt, kind="stable")
    pk = Picks(np.asarray(st)[o], np.asarray(ph)[o], np.asarray(tt)[o], np.ones(len(tt)))
    return pk, np.asarray(owner)[o], np.asarray(ev)


def match(det_of_pick, owner, n_true, n_pred, threshold=0.5):
    """precision, recall, F1 by the pick-overlap rule (score = shared / max(set sizes) >= threshold; each predicted
    event takes the first true event that qualifies, as PyOcto's benchmark), and the recall when each true event may
    be matched once (n_unique)."""
    pred = {}
    for i, e in enumerate(det_of_pick):
        if e >= 0:
            pred.setdefault(int(e), set()).add(i)
    true = {}
    for i, e in enumerate(owner):
        if e >= 0:
            true.setdefault(int(e), set()).add(i)
    pairs = []
    for pe, ps in sorted(pred.items()):
        cand = {}
        for i in ps:
            if owner[i] >= 0:
                cand[int(owner[i])] = cand.get(int(owner[i]), 0) + 1
        for te in sorted(cand):
            if cand[te] / max(len(ps), len(true[te])) >= threshold:
                pairs.append((pe, te))
                break
    tp = len(pairs)
    prec = tp / n_pred if n_pred else float("nan")
    rec = tp / n_true
    f1 = 2 * prec * rec / (prec + rec) if tp else 0.0
    return dict(precision=prec, recall=rec, f1=f1, tp=tp, n_pred=n_pred, n_true=n_true,
                recall_unique=len({t for _p, t in pairs}) / n_true)


def run(args):
    from ttmt.associate import pipeline as PL
    stations = load_stations(args.stations)
    vel = Velocity1D.hadley_kanamori()
    out = {"stations": len(stations.code), "station_file": args.stations, "station_source": station_source(args.stations),
           "seed": args.seed, "runs": []}
    for n in args.events:
        for noise in args.noise:
            rng = np.random.default_rng([args.seed, n, int(round(noise * 100))])
            pk, owner, ev = generate(stations, n, noise, rng, vel)
            for m in args.methods:
                fr = Frame(stations, XLIM, YLIM, ZLIM, vel, time_before=TIME_BEFORE, workdir=args.workdir,
                           fixed={m: FIXED[m]})
                t0 = time.time()
                try:
                    det, picks = PL.associator(m).associate(fr, pk, CFG[m])
                    r = match(picks.det, owner, len(ev), len(det))
                    r["status"] = "ok"
                except Exception as e:                              # noqa: BLE001
                    r = dict(status="error", error=repr(e)[:300])
                r.update(method=m, n_events=n, noise=noise, n_picks=len(pk), wall_s=round(time.time() - t0, 1))
                out["runs"].append(r)
                print(json.dumps(r), flush=True)
                with open(args.out, "w") as fh:
                    json.dump(out, fh, indent=1)


def report(path, tol=0.05):
    res = json.load(open(path))
    print(f"stations: {res['stations']} ({res['station_file']})")
    print("method   events noise  recall  S10    rel.diff  pass   precision S10")
    passed = {}
    for r in res["runs"]:
        ref = S10[r["method"]].get((r["n_events"], r["noise"]))
        if r["status"] != "ok" or ref is None:
            print(f"{r['method']:8s} {r['n_events']:6d} {r['noise']:4.1f}  {r['status']}  ref {ref}")
            continue
        rel = (r["recall"] - ref[1]) / ref[1]
        ok = abs(rel) <= tol
        passed.setdefault(r["method"], []).append(ok)
        print(f"{r['method']:8s} {r['n_events']:6d} {r['noise']:4.1f}  {r['recall']:.3f}  {ref[1]:.2f}  {rel:+.3f}"
              f"    {'yes' if ok else 'NO ':4s}  {r['precision']:.3f}     {ref[0]:.2f}")
    for m, v in passed.items():
        print(f"{m}: {sum(v)}/{len(v)} conditions within {tol:.0%} of the published recall")
    return passed


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fetch-stations")
    f.add_argument("out")
    r = sub.add_parser("run")
    r.add_argument("--stations", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--events", type=int, nargs="+", default=[100, 500, 2000])
    r.add_argument("--noise", type=float, nargs="+", default=[0.3, 1.0, 3.0])
    r.add_argument("--methods", nargs="+", default=["gamma", "pyocto", "real"])
    r.add_argument("--seed", type=int, default=0)
    r.add_argument("--workdir", default=None)
    p = sub.add_parser("report")
    p.add_argument("path")
    a = ap.parse_args(argv)
    if a.cmd == "fetch-stations":
        fetch_stations(a.out)
    elif a.cmd == "run":
        run(a)
    else:
        report(a.path)


if __name__ == "__main__":
    sys.exit(main())
