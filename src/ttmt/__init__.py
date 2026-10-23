"""ttmt: timing-tolerant moment-tensor detection (package of Phase 1, ADR-0002).

Every function a confirmatory test uses lives here, ported from the frozen scripts of ``code/`` and numerically
identical to them under the legacy conventions. Conventions that change numbers are explicit parameters
(``ttmt.conventions``). Importing the package does not import numpy, so ``ttmt.runtime.pin_blas()`` can run first.
"""
__version__ = "0.1.0"
