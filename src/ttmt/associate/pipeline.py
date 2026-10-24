"""Pick-based pipelines (picker + associator, tuned jointly) and the status of every method of the ATLAS-0 card.

A pipeline trial reads the stored probability traces of its picker for the record (pickers.load_probabilities; the
file of the trial's weights when the picker has a weights parameter), extracts the picks with the trial's peak
distance and threshold, and runs the associator with the trial's associator parameters. It returns a
``ttmt.atlas.schema.MethodRun``; :func:`runner` adapts it to ``tune.run``.
"""
import importlib

from ..atlas.schema import MethodRun
from . import pickers

ASSOCIATOR_MODULES = {"gamma": "gamma_assoc", "pyocto": "pyocto_assoc", "real": "real_assoc"}
MIN_PICKS_KEY = {"gamma": "as_min_picks", "pyocto": "as_min_picks", "real": "as_nps0"}

# every method of the card: "ready" (wrapper here, tested), "pending" (what is missing), "not runnable" (why)
STATUS = {
    "phasenet+gamma": "ready", "phasenet+pyocto": "ready", "phasenet+real": "ready",
    "eqtransformer+gamma": "ready", "eqtransformer+pyocto": "ready",
    "phaseno+gamma": "pending: PhaseNO inference to the probability-file layout (notebook code, GPU machine)",
    "phaseno+pyocto": "pending: as phaseno+gamma",
    "eqcorrscan": "pending: wrapper (ttmt-eqc env, Linux or conda-forge Windows build)",
    "quakemigrate": "pending: wrapper and install (not in the comparator environments of PR #8)",
    "malmi": "pending: wrapper and install (not in the comparator environments of PR #8)",
    "eqnet": "not runnable: the 2022 EQNet model has no released weights (environments/DEVIATIONS.md)",
}


def associator(name):
    return importlib.import_module(f".{ASSOCIATOR_MODULES[name]}", __package__)


def run_pipeline(method, cfg, frame, prob_path):
    """one trial of a pick-based pipeline on one record. prob_path: the probability file, or a callable
    weights -> path when the picker's weights are tuned."""
    pk, asc = method.split("+", 1)
    path = prob_path(cfg.get("pk_weights")) if callable(prob_path) else prob_path
    prob, dt, t0, _meta = pickers.load_probabilities(path, frame.codes)
    picks = pickers.picks_from_probabilities(prob, dt, cfg["pk_threshold"], cfg["pk_peak_distance"], t0)
    a = {k: v for k, v in cfg.items() if k.startswith("as_")}
    det, picks = associator(asc).associate(frame, picks, a)
    return MethodRun(method, det, picks, dict(cfg), int(cfg[MIN_PICKS_KEY[asc]]))


def runner(method):
    """tune.run runner: the record's data must be {"frame": Frame, "prob": path or callable weights -> path}."""
    def run(cfg, record):
        return run_pipeline(method, cfg, record.data["frame"], record.data["prob"]).det
    return run
