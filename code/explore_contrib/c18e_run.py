"""
c18e_run.py - exploration C18e (seed 999 only; Alto Tiberina day 2010-07-23 with quality control qc1, the noise of
t2m_out/noise/100723_qc1; networks j = 3100-3111, new; TAGs 5181-5186, new, checked by grep to appear in no earlier
file).  Question after C18c: do statistics that bound a station's influence less than the cap 6 (cap 24, logarithmic
compression with c = 6) keep the robustness against real-noise transients that Task 2M confirmed for the cap and
HUBER-N against BOMP-PK?  Same generator, windows, endpoints and code paths as Task 2M (t2m_lib.Net2M, t2m_run.run_net);
only the TAGs, the networks and the method list change.
    python c18e_run.py ARM PARTSDIR [j0,j1,...]
Methods: MT@S0 (cap 6), MT-CAP24@S0, MT-LOG6@S0, BOMP@S0, HUBER-N@S0.
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
_vectors_orig = C._vectors


def _vectors_var(Rr, tpk, capv):
    if capv is not None and capv < 0:
        v = _vectors_orig(Rr, tpk, None)
        c = -capv
        e = (v * v).sum(-1, keepdims=True)
        return v * np.sqrt(c * c * np.log1p(e / (c * c)) / np.maximum(e, 1e-24))
    return _vectors_orig(Rr, tpk, capv)


C._vectors = _vectors_var
VAR = {"MT-CAP24@S0": 24.0, "MT-LOG6@S0": -6.0}


class Net18(L2.Net2M):
    def run(self, name, MFn, kmax):
        if name in VAR:
            E, pr, _EF, _prF = self.eng["S0"]
            return C.mt_loop(E, MFn, kmax, mode="hard", pr=pr, W=4, capv=VAR[name])
        return super().run(name, MFn, kmax)


L2.Net2M = Net18
CFG = dict(seed=999, noise_dirs=[os.path.join(W, "t2m_out/noise/100723_qc1")], j0=3100, nblocks=[12],
           nwin={"STD": 6, "DENSE": 16, "NOISE": 4})

if __name__ == "__main__":
    arm, pdir = sys.argv[1], sys.argv[2]
    L2.configure(CFG["seed"], CFG["noise_dirs"], CFG["j0"], CFG["nblocks"], CFG["nwin"], ["MT@S0"])
    L2.CFG["methods"] = ["MT@S0", "MT-CAP24@S0", "MT-LOG6@S0", "BOMP@S0", "HUBER-N@S0"]
    nets = L2.CFG["nets"] if len(sys.argv) <= 3 else [int(x) for x in sys.argv[3].split(",")]
    os.makedirs(pdir, exist_ok=True)
    for j in nets:
        for cell in L2.CELLS:
            try:
                TR.run_net(arm, cell, j, pdir)
            except RuntimeError as e:
                print("FAILED", cell, j, e, flush=True)
