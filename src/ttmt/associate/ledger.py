"""The tuning ledger of ATLAS-0 (card, "Tuning"): one CSV row per configuration tried, appended, never rewritten.

Path: ``experiments/exploratory/ATLAS-0/tuning_ledger.csv``. The card asks for (method, trial, parameters, objective,
time); the ledger also keeps what an audit of the budget needs: the trial's seed, both regime recalls and their target
counts, the status of the trial (a trial that fails still spends its place in the budget), the code version, the
environment and a hash of the tuning partition, so a reader can tell which data every number came from. The SHA-256
of the file enters the freeze of the atlas.
"""
import csv
import hashlib
import json
import math
import os
import time
from pathlib import Path

from ..config import REPO_ROOT

LEDGER = REPO_ROOT / "experiments" / "exploratory" / "ATLAS-0" / "tuning_ledger.csv"
COLUMNS = ("method", "trial", "seed", "params", "objective", "recall_T", "recall_C", "n_T", "n_C", "hours",
           "status", "error", "wall_s", "time_utc", "code", "env", "partition")


LEDGER_A3 = LEDGER.with_name("tuning_ledger_a3.csv")       # GaMMA under ATLAS-0 amendment 3 (E63)


def require_a3(method, path, action):
    """GaMMA after amendment 3 runs and freezes only in LEDGER_A3: the registered ledger already holds its 60 rows
    (a run would skip them all) and its frozen best is the provisional one of registration 1."""
    if method.endswith("+gamma") and Path(path).resolve() != LEDGER_A3.resolve():
        raise ValueError(f"{action} {method}: GaMMA after ATLAS-0 amendment 3 (E63) uses only --ledger {LEDGER_A3}, "
                         f"not {path}")


def params_json(cfg):
    return json.dumps(cfg, sort_keys=True, separators=(",", ":"))


def read(path=LEDGER):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        r = csv.DictReader(fh)
        if tuple(r.fieldnames or ()) != COLUMNS:
            raise ValueError(f"{path}: columns {r.fieldnames} are not the ledger's {COLUMNS}")
        return list(r)


def append(row, path=LEDGER):
    """append one row (dict with the COLUMNS); refuses a (method, trial) already in the ledger."""
    row = {k: row.get(k, "") for k in COLUMNS}
    rows = read(path)
    if any(r["method"] == row["method"] and r["trial"] == str(row["trial"]) for r in rows):
        raise ValueError(f"trial {row['trial']} of {row['method']} is already in the ledger (append-only)")
    new = not os.path.exists(path) or os.path.getsize(path) == 0
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, lineterminator="\n")
        if new:
            w.writeheader()
        w.writerow({k: _fmt(v) for k, v in row.items()})


def _fmt(v):
    if isinstance(v, float):
        return "nan" if math.isnan(v) else repr(round(v, 6))
    return v


def done(method, path=LEDGER):
    """trial numbers of method already in the ledger (any status)."""
    return {int(r["trial"]) for r in read(path) if r["method"] == method}


def sha256(path=LEDGER):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


def best(method, path=LEDGER):
    """the ledger row of method with the largest objective among the trials that ran (ties: lowest trial)."""
    ok = [r for r in read(path) if r["method"] == method and r["status"] == "ok" and r["objective"] != "nan"]
    if not ok:
        return None
    return max(ok, key=lambda r: (float(r["objective"]), -int(r["trial"])))


def freeze(method, out, path=LEDGER):
    """write the frozen configuration of method (its best trial) with the ledger hash and the budget actually spent;
    this file is what the atlas thread receives."""
    require_a3(method, path, "freeze")
    rows = [r for r in read(path) if r["method"] == method]
    b = best(method, path)
    if b is None:
        raise ValueError(f"no successful trial of {method}")
    rec = dict(method=method, trial=int(b["trial"]), seed=int(b["seed"]), params=json.loads(b["params"]),
               objective=float(b["objective"]), recall_T=float(b["recall_T"]), recall_C=float(b["recall_C"]),
               trials_in_ledger=len(rows), trials_ok=sum(r["status"] == "ok" for r in rows),
               ledger_sha256=sha256(path), code=b["code"], env=b["env"], partition=b["partition"],
               frozen_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    with open(out, "w") as fh:
        json.dump(rec, fh, indent=1, sort_keys=True)
    return rec


def sensitivity(method, path=LEDGER):
    """per parameter: for numeric ones the Spearman correlation with the objective and the objective of the best
    trial in each tercile of the parameter; for categorical ones the best and median objective per option."""
    import numpy as np
    rows = [r for r in read(path) if r["method"] == method and r["status"] == "ok" and r["objective"] != "nan"]
    if not rows:
        return {}
    P = [json.loads(r["params"]) for r in rows]
    y = np.array([float(r["objective"]) for r in rows])
    out = {}
    for k in sorted(P[0]):
        v = [p[k] for p in P]
        if all(isinstance(x, (int, float)) and not isinstance(x, bool) for x in v):
            x = np.array(v, float)
            rx, ry = np.argsort(np.argsort(x)), np.argsort(np.argsort(y))
            rho = float(np.corrcoef(rx, ry)[0, 1]) if len(x) > 2 and x.std() > 0 and y.std() > 0 else float("nan")
            q = np.quantile(x, [1 / 3, 2 / 3])
            terc = np.digitize(x, q)
            out[k] = dict(spearman=rho, best_by_tercile=[float(y[terc == t].max()) if (terc == t).any() else None
                                                         for t in range(3)])
        else:
            keys = sorted({json.dumps(x) for x in v})
            out[k] = {kk: dict(best=float(y[[json.dumps(x) == kk for x in v]].max()),
                               median=float(np.median(y[[json.dumps(x) == kk for x in v]])),
                               n=int(sum(json.dumps(x) == kk for x in v))) for kk in keys}
    return out
