"""Registry of real datasets and the seal on the reserved confirmatory sequences (plan T3.R3, T3.R4; protocol 2.2).

Every real dataset that enters the project has a manifest ``data/manifests/real/<id>.json`` with its source,
licence, partitions in UTC and, once fetched, the SHA-256 of every file. The reserved confirmatory sequences are
listed in ``data/reserved/reservation.json``; their windows are sealed until the phase-6 registration that uses them.

The seal is enforced in code, not by convention:

- :func:`check_read` takes the network, the station coordinates and the time span a reader is about to read, and
  raises :class:`SealedError` if the read touches a sealed window. It is called before decoding by
  ``ttmt.datasets.Dataset.load_waveforms`` (every registered real dataset, and so ``ttmt.synth.realdata``), by
  ``ttmt.qc.waveform_qc`` with ``samples=True`` and by ``tools/data/seal_listing.py``. A new reader of real data must
  call it too (``tests/test_registry.py`` checks the readers above).
- Two purposes pass a sealed window: ``"fetch"`` (bytes are stored and hashed, never decoded) and ``"qc-metadata"``
  (headers, StationXML and catalogue aggregates only, ``ttmt.qc`` without ``samples``).
- Opening a window for a confirmatory run needs ``TTMT_UNSEAL=<sequence id>:<registration id>``, a registration id
  that the reservation lists for that sequence (added by an amendment commit after the registration is anchored), and
  that id present in ``preregistration/PREREGISTRATION_LOG_EN.txt``. Every unsealed read is appended to
  ``<work root>/unseal_log.jsonl``.

Development partitions of registered datasets are enforced by ``ttmt.partition.Partition`` / ``SegmentReader``; this
module adds the cross-dataset seal and the manifest checks.
"""
import hashlib
import json
import math
import re
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from . import config
from .partition import PartitionError

MANIFEST_DIR = config.REPO_ROOT / "data" / "manifests" / "real"
RESERVATION = config.REPO_ROOT / "data" / "reserved" / "reservation.json"
PREREG_LOG = config.REPO_ROOT / "preregistration" / "PREREGISTRATION_LOG_EN.txt"
PASS_PURPOSES = ("fetch", "qc-metadata")


class SealedError(PartitionError):
    """a read touches a reserved confirmatory window that is still sealed."""


def utc(s):
    """ISO 8601 string (naive = UTC) or epoch seconds -> epoch seconds."""
    if isinstance(s, (int, float)):
        return float(s)
    d = datetime.fromisoformat(str(s).replace("Z", "+00:00"))
    return (d if d.tzinfo else d.replace(tzinfo=timezone.utc)).timestamp()


def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                return h.hexdigest()
            h.update(b)


def canonical_sha256(obj):
    """SHA-256 of the canonical JSON of obj (sorted keys, no spaces)."""
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


SEALED_KEYS = ("rule", "sequences", "excluded")


def sealed_sha256(doc):
    """the hash that seals a reservation: rule, sequences and exclusions. ``unseal`` and ``amendments`` stay outside,
    so a later amendment that opens a sequence for its registration does not change it."""
    return canonical_sha256({k: doc.get(k) for k in SEALED_KEYS})


# ================================================================================================= reservation
@dataclass(frozen=True)
class Reserved:
    id: str
    networks: tuple
    lat: float
    lon: float
    radius_km: float
    start: float
    end: float
    registrations: tuple

    def touches(self, network, lat, lon, t0, t1):
        """True if [t0, t1) overlaps the window and the station lies in the area. Every network is sealed there,
        not only the listed ones: other networks record the same events."""
        if t1 <= self.start or t0 >= self.end:
            return False
        if lat is None or lon is None:
            return True
        d = 111.19 * math.hypot(lat - self.lat, (lon - self.lon) * math.cos(math.radians(self.lat)))
        return d <= self.radius_km


def load_reservation(path=None):
    """(document, [Reserved]); the document's ``sealed_sha256`` must equal :func:`sealed_sha256` of the document."""
    path = Path(path or RESERVATION)
    if not path.exists():
        return None, []
    doc = json.loads(path.read_text())
    got = sealed_sha256(doc)
    if got != doc.get("sealed_sha256"):
        raise SealedError(f"reservation {path} was changed after sealing: sha256 {got} != {doc.get('sealed_sha256')}")
    out = []
    for s in doc["sequences"]:
        w = s["window"]
        out.append(Reserved(s["id"], tuple(s["networks"]), float(w["lat"]), float(w["lon"]), float(w["radius_km"]),
                            utc(w["start"]), utc(w["end"]),
                            tuple(doc.get("unseal", {}).get(s["id"], []))))
    return doc, out


def _unsealed(seq, prereg_log=None):
    tok = os.environ.get("TTMT_UNSEAL", "")
    for part in filter(None, tok.split(",")):
        sid, _, reg = part.partition(":")
        if sid == seq.id and reg and reg in seq.registrations:
            log = Path(prereg_log or PREREG_LOG)
            if log.exists() and re.search(rf"(?<![\w.-]){re.escape(reg)}(?![\w.-])", log.read_text(errors="replace")):
                return reg
    return None


def check_read(network, lat, lon, t0, t1, purpose="analysis", reservation=None, prereg_log=None):
    """raise SealedError if reading [t0, t1) (epoch s or ISO) of a station at (lat, lon) of ``network`` touches a
    sealed window. lat/lon None (unknown station) is treated as inside every reserved area."""
    t0, t1 = utc(t0), utc(t1)
    _, seqs = load_reservation(reservation)
    for s in seqs:
        if not s.touches(network, lat, lon, t0, t1):
            continue
        if purpose in PASS_PURPOSES:
            continue
        reg = _unsealed(s, prereg_log)
        if reg is None:
            raise SealedError(f"{network} ({lat}, {lon}) [{t0:.0f}, {t1:.0f}) touches the sealed confirmatory "
                              f"sequence {s.id!r}; purpose {purpose!r} is not allowed before its registration")
        log = config.work_root() / "unseal_log.jsonl"
        log.parent.mkdir(parents=True, exist_ok=True)
        with open(log, "a") as f:
            f.write(json.dumps({"seq": s.id, "registration": reg, "network": network, "lat": lat, "lon": lon,
                                "t0": t0, "t1": t1, "purpose": purpose,
                                "at": datetime.now(timezone.utc).isoformat()}) + "\n")


# ================================================================================================= manifests
def load_manifest(dataset_id, manifest_dir=MANIFEST_DIR):
    return json.loads((Path(manifest_dir) / f"{dataset_id}.json").read_text())


def verify_files(manifest, root=None):
    """compare the files under ``root`` (default ``$TTMT_DATA_ROOT/<id>``) with the manifest's SHA-256 list.
    Returns (n_ok, missing, mismatched). Hashing reads bytes only, so it is allowed on sealed data."""
    root = Path(root) if root else config.data_root() / manifest["id"]
    ok, missing, bad = 0, [], []
    for rel, sha in sorted(manifest.get("files", {}).items()):
        p = root / rel
        if not p.exists():
            missing.append(rel)
        elif sha256_file(p) != sha:
            bad.append(rel)
        else:
            ok += 1
    return ok, missing, bad
