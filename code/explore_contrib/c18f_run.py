"""
c18f_run.py - exploration C18f (seed 999; the same networks 3100-3107, TAGs 5181-5186 and windows as C18e, so the
methods can be paired window by window with c18e_run.py).  Question raised by the fifth independent audit: is the
advantage of the bounded statistics in Task 2M an advantage over BOMP-PK as such, or over a comparator without any
per-station robustness?  Two robust versions of the comparator:
  BOMP-CAP6@S0  BOMP-PK whose block energies are computed from peak vectors capped at 6 noise deviations (the pursuit,
                the Gram refit and the deflation are unchanged)
  BOMP-STA@S0   BOMP-PK on matched-filter output normalised per station by its own robust level in the window
    python c18f_run.py ARM PARTSDIR [j0,j1,...]
"""
import os
import sys

import numpy as np

W = "/home/claude/work"
sys.path[:0] = [os.path.join(W, "explore_contrib"), W]
import t2m_lib as L2  # noqa: E402
import t2m_run as TR  # noqa: E402

import c10_lib as C  # noqa: E402

L2.TAG_STD, L2.TAG_DENSE, L2.TAG_CAT, L2.TAG_STAT, L2.TAG_SUB, L2.TAG_NOISE0 = 5181, 5182, 5183, 5184, 5185, 5186
_bomp_energy_orig = C.bomp_energy
_vectors_orig = C._vectors


def _bomp_energy_cap(E, R, tpk, sel=None):
    C._vectors = lambda Rr, t, capv: _vectors_orig(Rr, t, 6.0)
    try:
        return _bomp_energy_orig(E, R, tpk, sel)
    finally:
        C._vectors = _vectors_orig


class Net18f(L2.Net2M):
    def run(self, name, MFn, kmax):
        E, pr, _EF, _prF = self.eng["S0"]
        if name == "BOMP-CAP6@S0":
            C.bomp_energy = _bomp_energy_cap
            try:
                return C.bomp_loop(E, None, MFn, kmax, mode="hard", pr=pr, W=4)
            finally:
                C.bomp_energy = _bomp_energy_orig
        if name == "BOMP-STA@S0":
            sds = np.median(np.abs(MFn), axis=(1, 2)) / 0.6745
            return C.bomp_loop(E, None, MFn / np.maximum(sds, 1e-12)[:, None, None], kmax, mode="hard", pr=pr, W=4)
        return super().run(name, MFn, kmax)


L2.Net2M = Net18f
CFG = dict(seed=999, noise_dirs=[os.path.join(W, "t2m_out/noise/100723_qc1")], j0=3100, nblocks=[12],
           nwin={"STD": 6, "DENSE": 16, "NOISE": 4})

if __name__ == "__main__":
    arm, pdir = sys.argv[1], sys.argv[2]
    L2.configure(CFG["seed"], CFG["noise_dirs"], CFG["j0"], CFG["nblocks"], CFG["nwin"], ["MT@S0"])
    L2.CFG["methods"] = ["BOMP-CAP6@S0", "BOMP-STA@S0"]
    nets = L2.CFG["nets"] if len(sys.argv) <= 3 else [int(x) for x in sys.argv[3].split(",")]
    os.makedirs(pdir, exist_ok=True)
    for j in nets:
        for cell in L2.CELLS:
            try:
                TR.run_net(arm, cell, j, pdir)
            except RuntimeError as e:
                print("FAILED", cell, j, e, flush=True)
