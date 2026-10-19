"""Run the frozen legacy scripts of ``code/`` unchanged, with their hard-coded paths redirected at import time.

The legacy bytes are never edited (``code/LEGACY_INDEX.sha256``). The scripts expect the flat working folder
``/home/claude/work`` (= ``code/``) on ``sys.path`` and read their data from module constants (``r0_lib.DATA``,
``r1_lib.DATA``). This harness puts ``code/`` and its sub-folders on ``sys.path`` and installs an import hook that sets
those constants right after the module body runs, before any other module can read them.

    python -m ttmt.legacy SCRIPT.py [ARGS ...]          e.g.  python -m ttmt.legacy code/r1/r1_analyze.py eval ...
"""
import importlib.abc
import importlib.machinery
import os
import runpy
import sys

from . import config
from .runtime import pin_blas

LEGACY_DIRS = ("", "explore_contrib", "r0", "r1", "hyb")


def default_patches():
    """module name -> {constant: value} applied when the module is first imported."""
    return {
        "r0_lib": {"DATA": str(config.dataset_dir("ridgecrest_2019"))},
        "r1_lib": {"DATA": str(config.dataset_dir("hinet_2012"))},
    }


class _PatchLoader(importlib.abc.Loader):
    def __init__(self, inner, values):
        self.inner, self.values = inner, values

    def create_module(self, spec):
        return self.inner.create_module(spec)

    def exec_module(self, module):
        self.inner.exec_module(module)
        for k, v in self.values.items():
            setattr(module, k, v)


class _PatchFinder(importlib.abc.MetaPathFinder):
    def __init__(self, patches):
        self.patches = patches

    def find_spec(self, name, path=None, target=None):
        if name not in self.patches:
            return None
        spec = importlib.machinery.PathFinder.find_spec(name, path)
        if spec is None or spec.loader is None:
            return None
        spec.loader = _PatchLoader(spec.loader, self.patches[name])
        return spec


def install(patches=None, code_dir=None):
    """put the legacy folders on sys.path and register the constant patches (idempotent)."""
    code_dir = str(code_dir or config.CODE_DIR)
    sys.dont_write_bytecode = True       # as PYTHONDONTWRITEBYTECODE=1 of code/r1/env.sh: code/ stays untouched
    for sub in reversed(LEGACY_DIRS):
        p = os.path.join(code_dir, sub) if sub else code_dir
        if p not in sys.path:
            sys.path.insert(0, p)
    finder = next((f for f in sys.meta_path if isinstance(f, _PatchFinder)), None)
    if finder is None:
        finder = _PatchFinder({})
        sys.meta_path.insert(0, finder)
    # defaults never override a constant already registered; explicit patches do, constant by constant
    for mod, vals in default_patches().items():
        for k, v in vals.items():
            finder.patches.setdefault(mod, {}).setdefault(k, v)
    for mod, vals in (patches or {}).items():
        finder.patches.setdefault(mod, {}).update(vals)
    return finder


def run_script(script, argv, patches=None):
    """run one legacy script as __main__ with argv (list, without the script name)."""
    install(patches)
    script = os.path.abspath(script)
    d = os.path.dirname(script)
    if d not in sys.path:
        sys.path.insert(0, d)
    old = sys.argv
    sys.argv = [script] + list(argv)
    try:
        runpy.run_path(script, run_name="__main__")
    finally:
        sys.argv = old


def main():
    pin_blas()
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    run_script(sys.argv[1], sys.argv[2:])


if __name__ == "__main__":
    main()
