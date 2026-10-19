"""
run_nearest.py - post hoc (E52): run Task 2C or Task 2D with the exact-time readers moved from the floor sample to the
nearest sample, everything else as frozen.

The frozen engines compute the predicted sample as int(t / DT) (certified_fw._idx_node), which reads on average half a
sample early. Only the readers that read AT the predicted sample are handicapped by it; readers with the +-4-sample
peak search only move their window by at most one sample. This wrapper changes the convention for the exact-time
readers only, in memory, without editing the frozen files (their hashes still check):

  2C  task2c_methods._idx_node -> nearest. It is used only by BlockDict, i.e. by BOMP+SC, MTSTACK+SC and GL-FW+SC.
      MT-FW+SC and QM-3C+SC do not go through it (certified_fw keeps its own name) and are recomputed unchanged.
  2D  task2d_methods.node_vectors_v with peak=False (the MT-FW-EX+SC ablation) -> nearest. The peak-search methods
      (BOMP-PK+SC, MTSTACK-PK+SC, GL-FW-PK+SC, MT-FW-NOCAP+SC) keep the frozen centre, so the comparison with them is
      the one E52 says is handicapped.

Usage (in a copy of code/ whose task2a_parts/ holds the regenerated 2A parts):
    OPENBLAS_CORETYPE=Prescott python posthoc_e52/run_nearest.py 2c run --part 2 --workers 4
    OPENBLAS_CORETYPE=Prescott python posthoc_e52/run_nearest.py 2d run --part 1 --workers 4
"""
import os
import sys

import numpy as np

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if CODE not in sys.path:
    sys.path.insert(0, CODE)

import certified_fw as CF  # noqa: E402
import task1x_mtglrt as X  # noqa: E402
from task0_common import DT, NT, VP, VS  # noqa: E402


def idx_nearest(STA, src, taus, tcorr):
    """certified_fw._idx_node with np.rint instead of truncation."""
    d = np.linalg.norm(STA - src, axis=1)
    tt = np.stack([d / VP, d / VS])[None] + taus[:, None, None]
    if tcorr is not None:
        tt = tt + tcorr.T[None]
    return np.clip(np.rint(tt / DT).astype(int), 0, NT - 1)


def patch_2c():
    import task2c_methods as M
    assert M._idx_node is CF._idx_node
    M._idx_node = idx_nearest


def patch_2d():
    import task2d_methods as D
    frozen = D.node_vectors_v

    def node_vectors_v(MFn, STA, src, taus, tcorr, peak=True, cap=True):
        if peak:
            return frozen(MFn, STA, src, taus, tcorr, peak, cap)
        st = np.arange(len(STA))
        idx = idx_nearest(STA, src, taus, tcorr)
        v = MFn[st[None, None, :], :, idx]
        if cap:
            nv = np.linalg.norm(v, axis=-1, keepdims=True)
            v = v * np.minimum(1.0, X.CAPV / np.maximum(nv, 1e-12))
        return v, idx

    D.node_vectors_v = node_vectors_v


def main():
    task = sys.argv[1]
    sys.argv = [sys.argv[0]] + sys.argv[2:]
    if task == "2c":
        patch_2c()
        import task2c_run as run
    elif task == "2d":
        patch_2d()
        import task2d_run as run
    else:
        raise SystemExit("task must be 2c or 2d")
    run.main()


if __name__ == "__main__":
    main()
