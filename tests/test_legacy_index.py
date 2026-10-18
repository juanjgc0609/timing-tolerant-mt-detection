"""code/ is frozen evidence (ADR-0002): every file is listed in code/LEGACY_INDEX.sha256 with its SHA-256, and no file
is added, removed or changed."""
import hashlib
from pathlib import Path

CODE = Path(__file__).resolve().parents[1] / "code"
INDEX = CODE / "LEGACY_INDEX.sha256"
# generated locally and ignored by git: bytecode, third-party clones, model files built on first use, run parts
IGNORED_PARTS = {"__pycache__", "gamma_src", "c19_parts", "r0_verif"}
IGNORED_NAMES = {"pyocto_jma2001.model", "LEGACY_INDEX.sha256"}


def _index():
    out = {}
    for line in INDEX.read_text().splitlines():
        h, name = line.split("  ", 1)
        out[name] = h
    return out


def _sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def test_index_lists_every_file_with_its_hash():
    idx = _index()
    assert len(idx) >= 350
    bad = [n for n, h in idx.items() if not (CODE / n).is_file() or _sha(CODE / n) != h]
    assert not bad, f"legacy files changed or missing (check line endings: .gitattributes '* -text'): {bad[:10]}"


def test_no_unlisted_files_in_code():
    idx = _index()
    extra = []
    for p in CODE.rglob("*"):
        rel = p.relative_to(CODE)
        if p.is_dir() or IGNORED_PARTS & set(rel.parts) or p.name in IGNORED_NAMES or p.suffix == ".pyc":
            continue
        if rel.as_posix() not in idx:
            extra.append(rel.as_posix())
    assert not extra, f"files added to code/: {extra[:10]}"


def test_freeze_lists_match_present_files():
    """the freeze lists of R0, R1 and 2M hold for every listed file present in code/ (paths relative to code/)."""
    lists = [CODE / "r1" / "FREEZE_R1.sha256", CODE / "r0" / "R0_FREEZE.txt",
             CODE.parent / "results" / "task2M" / "T2M_FREEZE.sha256"]
    checked = 0
    for fl in lists:
        for line in fl.read_text().splitlines():
            parts = line.split()
            if len(parts) < 2 or len(parts[0]) != 64:
                continue
            p = CODE / parts[1].removeprefix("./")
            if p.is_file():
                assert _sha(p) == parts[0], (fl.name, parts[1])
                checked += 1
    assert checked >= 300
