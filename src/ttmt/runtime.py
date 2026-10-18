"""Process-level numerical settings shared by the tests, the runners and the legacy harness.

The frozen runs pinned one BLAS thread and, from Task 2B on, the OpenBLAS kernel (error E31: numpy's OpenBLAS picks
its compute kernel when it loads, and the SkylakeX and Prescott kernels round matrix products differently). The
registered outputs were computed with ``OPENBLAS_CORETYPE=Prescott`` (the default kernel of the original machine) and
``OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=1`` (``code/r1/env.sh``, ``code/r0_lib.py``).

:func:`pin_blas` must run before numpy is first imported: OpenBLAS reads these variables once, when it loads. Every
entry point (``tests/conftest.py``, ``ttmt.legacy``, the regression runners) calls it first.
"""
import os
import sys

LEGACY_CORETYPE = "Prescott"
LEGACY_THREADS = 1
_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")


def pin_blas(coretype=LEGACY_CORETYPE, threads=LEGACY_THREADS, strict=True):
    """Set the BLAS kernel and thread count of this process the way the legacy runs did.

    Returns the environment values that were set. With ``strict`` (default) it raises if numpy is already imported
    with different values, because the setting would then be ignored silently."""
    want = {"OPENBLAS_CORETYPE": str(coretype), **{v: str(int(threads)) for v in _VARS}}
    if "numpy" in sys.modules:
        bad = {k: os.environ.get(k) for k, v in want.items() if os.environ.get(k) != v}
        if bad and strict:
            raise RuntimeError(f"numpy was imported before pin_blas(); these BLAS settings would be ignored: {bad}")
    for k, v in want.items():
        os.environ[k] = v
    return want


def blas_info():
    """The BLAS libraries actually loaded (threadpoolctl), for the record of a run: kernel, threads, version."""
    from threadpoolctl import threadpool_info
    out = []
    for d in threadpool_info():
        if d.get("user_api") == "blas":
            out.append({k: d.get(k) for k in ("internal_api", "version", "architecture", "num_threads", "filepath")})
    return out


def is_exact_platform():
    """True where the regressions demand bit-for-bit equality (Linux, as the frozen runs). Elsewhere the regression
    tests apply the declared tolerance (``ttmt.testing.Tolerance``)."""
    return sys.platform.startswith("linux")
