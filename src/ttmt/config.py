"""Paths from configuration (protocol, Section 5.4).

``TTMT_DATA_ROOT``: folder holding the raw datasets with the sub-folder names of ``data/README.md``
(``skillsbench_eqassoc``, ``MatchLocate2/example``, ``pympa37_altotiberina_100723``, ``pympa37_days``). Default: the
legacy location ``/home/claude/realdata``.
``TTMT_WORK_ROOT``: folder for raw run outputs and regenerated records. Default: ``<repository>/work`` (ignored by git).
"""
import os
from pathlib import Path

LEGACY_DATA_ROOT = "/home/claude/realdata"
REPO_ROOT = Path(__file__).resolve().parents[2]
CODE_DIR = REPO_ROOT / "code"
RESULTS_DIR = REPO_ROOT / "results"

DATASETS = {
    "ridgecrest_2019": "skillsbench_eqassoc",
    "hinet_2012": "MatchLocate2/example",
    "altotiberina_c16": "pympa37_altotiberina_100723",
    "altotiberina_days": "pympa37_days",
}


def data_root():
    return Path(os.environ.get("TTMT_DATA_ROOT", LEGACY_DATA_ROOT))


def dataset_dir(name):
    """folder of one dataset under the data root (name: a key of DATASETS)."""
    return data_root() / DATASETS[name]


def has_dataset(name):
    return dataset_dir(name).is_dir()


def work_root():
    return Path(os.environ.get("TTMT_WORK_ROOT", REPO_ROOT / "work"))
