"""A1: duplicator of frozen PyOcto on the tune partition at the registered operating point (LOC+T-V, 30 false/h).
Reuses STAGE-LOSS/stage_loss.py (same picks, associator, objective). Read-only on the staged tune files."""
import importlib.util, json, os, sys, dataclasses
import numpy as np
from pathlib import Path
REPO = str(Path(__file__).resolve().parents[3])
spec = importlib.util.spec_from_file_location("sl", f"{REPO}/experiments/exploratory/STAGE-LOSS/stage_loss.py")
SL = importlib.util.module_from_spec(spec); spec.loader.exec_module(SL)
TU = SL.TU
fz = json.load(open(f"{REPO}/experiments/exploratory/COMPARATORS/frozen/phasenet+pyocto.json"))
cfg = dict(fz["params"]); acfg = {k: v for k, v in cfg.items() if k.startswith("as_")}
U = sys.argv[1].rstrip("/") + "/"          # the tuning partition folder (ss_ridgecrest_2019_egf_v1, ss_hinet_2012_egf_v1)
recs = [r for m in (U + "ss_ridgecrest_2019_egf_v1/tune_phasenet.json", U + "ss_hinet_2012_egf_v1/tune_phasenet.json")
        for r in TU.load_manifest(m, "phasenet+pyocto")]
print("partition", TU.partition_hash(recs), "frozen", fz["partition"], "records", len(recs))
labels = {}
for r in recs: labels.setdefault(r.name.split("__")[0], []).append(r)
cache = {}
for lab, rs in labels.items():
    p = SL.phasenet_picks(rs[0], cfg)
    cache[lab] = SL.associator("pyocto").associate(rs[0].data["frame"], p, acfg)[0]
dets = [cache[r.name.split("__")[0]] for r in recs]
d0 = dets[0]; print(type(d0), [f.name for f in dataclasses.fields(d0)] if dataclasses.is_dataclass(d0) else dir(d0)[-12:])
import pickle; pickle.dump(dets, open(os.devnull, "wb"))
ev = SL.evaluate(recs, dets)
print(f"comparator      T {ev['recall_T']:.4f}  C {ev['recall_C']:.4f}   (frozen {fz['recall_T']}, {fz['recall_C']})")
def dup(d, h, rng):
    hh = h if np.isscalar(h) else rng.uniform(h[0], h[1], len(d.t))
    kw = {}
    for f in dataclasses.fields(d):
        v = getattr(d, f.name)
        if isinstance(v, np.ndarray) and len(v) == len(d.t):
            kw[f.name] = np.concatenate([v, v])
        else:
            kw[f.name] = v
    kw["t"] = np.concatenate([d.t - hh, d.t + hh])
    return dataclasses.replace(d, **kw)
rng = np.random.default_rng(20261114)
for name, h in [(f"h={x}", x) for x in (0.05, 0.1, 0.15, 0.2, 0.25, 0.35, 0.5)] + [("U[0.05,0.5]", (0.05, 0.5))]:
    e = SL.evaluate(recs, [dup(d, h, rng) for d in dets])
    print(f"duplicator {name:12s} T {e['recall_T']:.4f}  C {e['recall_C']:.4f}  hits T {int((e['got'] & e['T']).sum())} of {int(e['T'].sum())}")
