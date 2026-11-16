"""B1: is the T-regime merge of frozen PyOcto structural, or travel-time misfit? Oracle picks at the injected times
(OR, as STAGE-LOSS) against oracle picks at the times of the ASSOCIATOR'S OWN 1-D model (ORM: zero travel-time
residual), each with the frozen settings, with min_interevent_time 0.1 s (mit), and with all tolerances 0.1 s (tight).
Per network. Same objective code as STAGE-LOSS (LOC+T-V, 30 false/h)."""
import importlib.util, json, sys
import numpy as np
from pathlib import Path
REPO = str(Path(__file__).resolve().parents[3])
spec = importlib.util.spec_from_file_location("sl", f"{REPO}/experiments/exploratory/STAGE-LOSS/stage_loss.py")
SL = importlib.util.module_from_spec(spec); spec.loader.exec_module(SL)
TU, Picks = SL.TU, SL.Picks
fz = json.load(open(f"{REPO}/experiments/exploratory/COMPARATORS/frozen/phasenet+pyocto.json"))
cfg = dict(fz["params"]); acfg = {k: v for k, v in cfg.items() if k.startswith("as_")}
U = sys.argv[1].rstrip("/") + "/"          # the tuning partition folder (ss_ridgecrest_2019_egf_v1, ss_hinet_2012_egf_v1)
def model_picks(r):
    tr, fr = r.truth, r.data["frame"]
    a = tr.arr_t_obs
    out = np.full(a.shape, np.nan)
    for e in range(tr.n_events):
        tt = fr.velocity.times(tr.xyz[e], fr.stations.xyz)           # (2, ns)
        out[e] = tr.t[e] + tt.T
    out[~np.isfinite(a)] = np.nan
    e, s, ph = np.nonzero(np.isfinite(out)); t = out[e, s, ph]; o = np.lexsort((ph, s, t))
    return Picks(s[o], ph[o], t[o], np.ones(len(o))), np.nanmedian(np.abs(out - a))
def model_times(r):
    tr, fr = r.truth, r.data["frame"]
    out = np.full(tr.arr_t_obs.shape, np.nan)
    for e in range(tr.n_events):
        out[e] = tr.t[e] + fr.velocity.times(tr.xyz[e], fr.stations.xyz).T
    return out
def static_picks(r0, others):
    """injected times minus a per-station, per-phase static term (median residual on the OTHER realisations)."""
    res = np.concatenate([o.truth.arr_t_obs - model_times(o) for o in others])     # (events, ns, 2)
    c = np.nanmedian(res, axis=0)
    a = r0.truth.arr_t_obs - c[None]
    e, s, ph = np.nonzero(np.isfinite(a)); t = a[e, s, ph]; o = np.lexsort((ph, s, t))
    rem = (r0.truth.arr_t_obs - model_times(r0)) - c[None]
    return Picks(s[o], ph[o], t[o], np.ones(len(o))), np.nanmedian(np.abs(rem))
conds = {"PN": ("pn", acfg), "OR": ("or", acfg), "OR-mit": ("or", dict(acfg, as_min_interevent_time=0.1)),
         "OR-tight": ("or", dict(acfg, as_min_interevent_time=0.1, as_tolerance=0.1, as_pick_match_tolerance=0.1)),
         "ORS": ("ors", acfg), "ORS-mit-tol0.3": ("ors", dict(acfg, as_min_interevent_time=0.1, as_tolerance=0.3, as_pick_match_tolerance=0.3)),
         "ORS-tight": ("ors", dict(acfg, as_min_interevent_time=0.1, as_tolerance=0.1, as_pick_match_tolerance=0.1)),
         "ORS-mit-tol0.2": ("ors", dict(acfg, as_min_interevent_time=0.1, as_tolerance=0.2, as_pick_match_tolerance=0.2)),
         "ORM": ("orm", acfg), "ORM-mit": ("orm", dict(acfg, as_min_interevent_time=0.1)),
         "ORM-mit-tol0.3": ("orm", dict(acfg, as_min_interevent_time=0.1, as_tolerance=0.3, as_pick_match_tolerance=0.3)),
         "ORM-tight": ("orm", dict(acfg, as_min_interevent_time=0.1, as_tolerance=0.1, as_pick_match_tolerance=0.1))}
for net in ("ss_ridgecrest_2019_egf_v1", "ss_hinet_2012_egf_v1"):
    recs = TU.load_manifest(U + net + "/tune_phasenet.json", "phasenet+pyocto")
    labels = {}
    for r in recs: labels.setdefault(r.name.split("__")[0], []).append(r)
    print(net, "records", len(recs))
    for name, (kind, ac) in conds.items():
        if not name.startswith("ORS"): continue
        cache, mis = {}, []
        for lab, rs in labels.items():
            r0 = rs[0]
            if kind == "pn": p = SL.phasenet_picks(r0, cfg)
            elif kind == "or": p = SL.oracle_picks(r0.truth)
            elif kind == "ors":
                p, m = static_picks(r0, [x[0] for l2, x in labels.items() if l2 != lab]); mis.append(m)
            else:
                p, m = model_picks(r0); mis.append(m)
            cache[lab] = SL.associator("pyocto").associate(r0.data["frame"], p, ac)[0]
        dets = [cache[r.name.split("__")[0]] for r in recs]
        ev = SL.evaluate(recs, dets)
        T = ev["T"]; iso = ev["iso"]
        print(f"  {name:15s} T {ev['recall_T']:.3f} ({int((ev['got']&T).sum())}/{int(T.sum())})  C {ev['recall_C']:.3f}  iso {ev['recall_iso']:.3f}  n_det {sum(len(d.t) for d in cache.values())}"
              + (f"  median |model - injected| {np.median(mis):.2f} s" if mis else ""), flush=True)
