"""The tuning loop of ATLAS-0: the same random search, budget and objective for every method (card, "Tuning";
ADR-0008 item 9; critique items CR-03 and CR-07).

Objective of a configuration: the development recall in T and in C (C includes the isolated events, as in
``ttmt.atlas.froc.gate``) at 30 false detections per hour under LOC+T-V, each weighted as in the endpoint (the
prevalence weights of the records, 1 by default), pooled over the tuning records of the method, with the operating
point chosen on those records (``ttmt.atlas.froc.hits``). The card names the two recalls but not how they are combined
into the one number a search maximises; the default here is their mean, ``objective = (recall_T + recall_C) / 2``, so
neither regime decides alone. Both recalls are kept in the ledger, so the atlas thread can reselect under another
combination without rerunning anything.

Only records of the tuning partition are accepted (``TuningRecord.partition == "tune"``); the split of noise hours
and injection seeds is :func:`tuning_split` (30 / 70 by a hash of the hour, card "Tuning").
"""
import hashlib
import json
import subprocess
import time
import traceback
from dataclasses import dataclass, field

import numpy as np

from ..atlas import froc as F
from . import ledger as L
from . import space as SP

BUDGET = F.BUDGETS[0]
TUNE_FRACTION = 0.30


def tuning_split(key, fraction=TUNE_FRACTION):
    """"tune" or "read" for a noise hour or an injection seed (any str), by the first 8 bytes of its SHA-256."""
    h = int.from_bytes(hashlib.sha256(str(key).encode()).digest()[:8], "big") / 2.0 ** 64
    return "tune" if h < fraction else "read"


@dataclass
class TuningRecord:
    """one record of the tuning partition seen by one station subset. truth: ttmt.atlas.schema.Truth; k and subsets
    as ttmt.atlas.froc.score_record; block_edges: record times (s) that cut the record into blocks; stratum: sequence
    or network id; data: whatever the method's runner needs to read the waveforms or the picks."""
    name: str
    truth: object
    k: int
    subsets: list
    centroid: tuple
    block_edges: np.ndarray
    partition: str
    stratum: int = 0
    weights: np.ndarray = None
    data: object = None
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.partition != "tune":
            raise ValueError(f"record {self.name!r} is in partition {self.partition!r}; tuning reads 'tune' only")
        e = np.asarray(self.block_edges, float)
        if e[0] > self.truth.t0 + 1e-9 or e[-1] < self.truth.t1 - 1e-9 or np.any(np.diff(e) <= 0):
            raise ValueError("block edges must increase and cover the record")
        self.block_edges = e

    def block(self, t):
        return int(np.clip(np.searchsorted(self.block_edges, t, side="right") - 1, 0, len(self.block_edges) - 2))

    @property
    def block_hours(self):
        return np.diff(self.block_edges) / 3600.0


def partition_hash(records):
    """identifies the tuning data of a run (record names and truth sizes), for the ledger."""
    s = json.dumps(sorted((r.name, r.k, int(r.truth.n_events)) for r in records))
    return hashlib.sha256(s.encode()).hexdigest()[:16]


def objective(records, dets, budget=BUDGET):
    """(objective, recall_T, recall_C, n_T, n_C, hours) of one configuration: dets[i] are the Detections of the
    method on records[i]."""
    hours = np.concatenate([r.block_hours for r in records])
    strata = np.concatenate([np.full(len(r.block_hours), r.stratum) for r in records])
    pool = F.Pool(hours, strata)
    off = 0
    for r, d in zip(records, dets):
        ev_block = np.array([r.block(t) for t in r.truth.t], int) + off
        F.score_record(pool, r.truth, r.k, r.subsets, r.centroid, ev_block, lambda t, r=r, o=off: r.block(t) + o,
                       {"m": d}, weights=r.weights)
        off += len(r.block_hours)
    got = F.hits(pool, "m", budget)
    T = pool.ev_regime == "T"
    Cm = (pool.ev_regime == "C") | (pool.ev_regime == "iso")
    rT, rC = F.recall(pool, got, T), F.recall(pool, got, Cm)
    obj = float(np.nanmean([rT, rC])) if np.isfinite([rT, rC]).any() else float("nan")
    return obj, rT, rC, int(T.sum()), int(Cm.sum()), float(hours.sum())


def code_version():
    try:
        sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True,
                             cwd=L.REPO_ROOT).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--", "src"], capture_output=True, text=True,
                               cwd=L.REPO_ROOT).stdout.strip()
        return sha[:12] + ("+dirty" if dirty else "")
    except Exception:                                              # noqa: BLE001
        return "unknown"


def run(method, runner, records, ledger=L.LEDGER, n_trials=SP.N_TRIALS, env="", trials=None, log=print):
    """the random search of one method. runner(config, record) -> Detections. Every trial of space.trials(method)
    not yet in the ledger is run and appended (resumable; a failing trial is appended with status "error" and an
    objective nan, and still counts in the budget). Returns the ledger rows of the method."""
    if not records:
        raise ValueError("no tuning records")
    have = L.done(method, ledger)
    part, code = partition_hash(records), code_version()
    for i, seed, cfg in (trials if trials is not None else SP.trials(method, n_trials)):
        if i in have:
            continue
        t0 = time.time()
        row = dict(method=method, trial=i, seed=seed, params=L.params_json(cfg), code=code, env=env, partition=part)
        try:
            dets = [runner(cfg, r) for r in records]
            obj, rT, rC, nT, nC, hrs = objective(records, dets)
            row.update(objective=obj, recall_T=rT, recall_C=rC, n_T=nT, n_C=nC, hours=hrs, status="ok")
        except Exception as e:                                     # noqa: BLE001
            row.update(objective=float("nan"), status="error",
                       error=(type(e).__name__ + ": " + str(e))[:300].replace("\n", " "))
            log(traceback.format_exc())
        row.update(wall_s=round(time.time() - t0, 2), time_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
        L.append(row, ledger)
        log(f"{method} trial {i}: {row['status']} objective {row.get('objective')} ({row['wall_s']} s)")
    return [r for r in L.read(ledger) if r["method"] == method]


REALISM_LOW, REALISM_HIGH = 0.2, 0.85


def realism_flags(ledger=L.LEDGER):
    """methods whose best development recall in T or in C lies outside [0.2, 0.85]: the design check of H2 must
    cover them (coordinator's brief of 24 October 2026)."""
    out = {}
    for m in sorted({r["method"] for r in L.read(ledger)}):
        b = L.best(m, ledger)
        if b is None:
            out[m] = "no successful trial"
            continue
        bad = {k: float(b[k]) for k in ("recall_T", "recall_C")
               if not REALISM_LOW <= float(b[k]) <= REALISM_HIGH}
        if bad:
            out[m] = bad
    return out


def load_manifest(path, method):
    """tuning records from a JSON manifest written by the data side (generator v2 or the real-data thread):
    {"records": [{"name", "partition", "truth": truth.npz, "stations": {"code": [...], "xyz": [[...]]}, "k",
    "subsets": [[station indices] or null per event], "centroid": [x, y], "block_edges": [...], "stratum",
    "prob": {"original": path, "instance": path, ...}, "frame": {"xlim", "ylim", "zlim", "time_before"},
    "velocity": "hadley_kanamori" | "jma2001" | [[top, vp], ...]}]}. truth.npz holds the fields of
    ttmt.atlas.schema.Truth (id, t, xyz, mag, arr_t, snr, t0, t1, target, pair, kmec, arr_t_obs). Paths are relative
    to the manifest. Only the pick-based pipelines read "prob"."""
    import os

    from ..atlas.schema import Stations, Truth
    from .base import Frame
    from .velocity import Velocity1D
    root = os.path.dirname(os.path.abspath(path))
    man = json.load(open(path))
    out = []
    for r in man["records"]:
        z = np.load(os.path.join(root, r["truth"]), allow_pickle=False)
        f = {k: z[k] for k in z.files}
        truth = Truth(f["id"], f["t"], f["xyz"], f["mag"], f["arr_t"], f["snr"], float(f["t0"]), float(f["t1"]),
                      f.get("target"), f.get("pair"), f.get("kmec"), f.get("arr_t_obs"))
        st = Stations(r["stations"]["code"], r["stations"]["xyz"])
        v = r.get("velocity", "hadley_kanamori")
        vel = getattr(Velocity1D, v)() if isinstance(v, str) else Velocity1D(tuple(map(tuple, v)))
        fr = r.get("frame", {})
        frame = Frame(st, tuple(fr["xlim"]), tuple(fr["ylim"]), tuple(fr.get("zlim", (0.0, 30.0))), vel,
                      time_before=fr.get("time_before", 60.0)) if "xlim" in fr else Frame.around(st, vel)
        prob = {w: os.path.join(root, p) for w, p in r.get("prob", {}).items()}
        subsets = [None if s is None else np.asarray(s, int) for s in r["subsets"]]
        out.append(TuningRecord(r["name"], truth, int(r["k"]), subsets, tuple(r["centroid"]), r["block_edges"],
                                r["partition"], int(r.get("stratum", 0)),
                                data={"frame": frame, "prob": (lambda w, p=prob: p[w or next(iter(p))])}))
    return out
