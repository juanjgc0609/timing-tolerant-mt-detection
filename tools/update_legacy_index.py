"""Append the files new to code/ to code/LEGACY_INDEX.sha256, keeping the list sorted.

The index can grow but never re-hash frozen code: if any file already listed is missing or has changed bytes, the
script changes nothing and exits with status 1. Work that adds files to code/ runs it and commits the index.
    python tools/update_legacy_index.py [--check]       (--check: report only, exit 1 if the index is not current)
"""
import hashlib
import subprocess
import sys
from pathlib import Path

CODE = Path(__file__).resolve().parents[1] / "code"
INDEX = CODE / "LEGACY_INDEX.sha256"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def tracked():
    """files under code/ that git tracks or would track (untracked and not ignored), the index excluded."""
    out = subprocess.run(["git", "ls-files", "--cached", "--others", "--exclude-standard", "--", "."], cwd=CODE,
                         capture_output=True, text=True, check=True).stdout.splitlines()
    return sorted({f for f in out if f and f != INDEX.name and (CODE / f).is_file()})


def main(argv):
    check = "--check" in argv
    idx = {}
    for line in INDEX.read_text().splitlines():
        h, name = line.split("  ", 1)
        idx[name] = h
    bad = [n for n, h in idx.items() if not (CODE / n).is_file() or sha(CODE / n) != h]
    if bad:
        print("REFUSED: indexed legacy files are missing or changed (the index is append-only):", *bad, sep="\n  ")
        return 1
    new = [f for f in tracked() if f not in idx]
    for f in new:
        print("new:", f)
    if not new:
        print(f"index current ({len(idx)} files)")
        return 0
    if check:
        return 1
    idx.update({f: sha(CODE / f) for f in new})
    with open(INDEX, "w", newline="\n") as fh:
        for name in sorted(idx, key=lambda s: s.encode()):
            fh.write(f"{idx[name]}  {name}\n")
    print(f"appended {len(new)}; index now {len(idx)} files")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
