"""M5-ACAP, gate G-M5a (robustness) of the v2 specification (M5; build order step 3, plan T4.13).

Question. On the 2M transient day (Alto Tiberina 2010-12-31, real noise with transients; Task 2M primary seed 2718),
is the recall of the ACAP ranking (relative cap c = 3 x median block norm) at most 0.03 below that of the registered
cap 6, with the LOC+T hit rule?

Design. The registered Task 2M pipeline (code/t2m_lib.py, code/t2m_run.py, configuration
results/task2M/cfg_2718.json: seed, networks, cells, windows), REAL arm, with three methods that differ ONLY in the
block statistic of the MT loop (station terms S0, hard reading, W = 4, no refinement, kmax as registered):
  MT@S0      registered cap 6 (identical to the 2M method of that name: a check of the pipeline)
  ACAP@S0    relative cap, c = 3 (ttmt.detector.ACAP's map, applied to the read vectors of c10_lib.mt_loop)
  NOCAP@S0   unbounded (the control that the transients hurt)
  ACAP5@S0   c = 5, sensitivity (added after the small run, outside the gate)
Endpoint: DENSE windows, LOC+T (5 km, 0.5 s), recall at <= 0.25 realised false events per window (the operating point
of 2M, re-chosen in every bootstrap replicate), hierarchical bootstrap of t2m_analyze (strata = cell, networks then
windows), B = 2000. Gate (declared before the run): ACAP@S0 - MT@S0 has its 2.5th percentile > -0.03 (the
non-inferiority reading of 2M-d). Also reported: LOC, STD, and per cell.

    python run_a.py run NOISE_DIR PARTSDIR [j,j,...|all] [cell,cell,...]
    python run_a.py analyze PARTSDIR [B] > run_a_out.txt
NOISE_DIR: the 101231 noise directory written by code/t2m_data.py (code/t2m_prep.sh, step "new"); its noise.npy is
checked against results/task2M/noise_arrays.sha256 (a mismatch is reported, not fatal: last-bit differences of a
different NumPy/OpenBLAS). Small first run: nets 2800 (one network per cell, 78 windows).
"""
import hashlib
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(ROOT, "code"))
import numpy as np  # noqa: E402

import t2m_lib as L2  # noqa: E402
import c10_lib as C  # noqa: E402

CFG = os.path.join(ROOT, "results", "task2M", "cfg_2718.json")
SHA = os.path.join(ROOT, "results", "task2M", "noise_arrays.sha256")
METHODS = ["MT@S0", "ACAP@S0", "NOCAP@S0", "ACAP5@S0"]
ACAP_C = {"ACAP@S0": 3.0, "ACAP5@S0": 5.0}         # ACAP5: sensitivity of c (declared after the small run, not gated)
REG_PARTS = os.path.join(ROOT, "results", "task2M", "parts_2718")


class AcapSpec(float):
    """marker passed as capv to c10_lib.mt_loop: value inf (no absolute cap), map = relative cap c x median."""
    def __new__(cls, c):
        o = float.__new__(cls, np.inf)
        o.c = float(c)
        return o


_vectors0 = C._vectors


def _vectors(R, tpk, capv):
    if not isinstance(capv, AcapSpec):
        return _vectors0(R, tpk, capv)
    v = _vectors0(R, tpk, None)                                   # (..., 2, ns, 3)
    r = np.linalg.norm(v, axis=-1)                                # (..., 2, ns)
    med = np.median(r.reshape(r.shape[:-2] + (-1,)), axis=-1)[..., None, None]
    return v * np.minimum(1.0, capv.c * med / np.maximum(r, 1e-12))[..., None]


C._vectors = _vectors
_run0 = L2.Net2M.run


def _run(self, name, MFn, kmax):
    if name in ("ACAP@S0", "ACAP5@S0", "NOCAP@S0"):
        E, pr, _EF, prF = self.eng["S0"]
        capv = AcapSpec(ACAP_C[name]) if name in ACAP_C else None
        return C.mt_loop(E, MFn, kmax, mode="hard", pr=pr, W=4, EF=None, prF=prF, capv=capv)
    return _run0(self, name, MFn, kmax)


L2.Net2M.run = _run
L2.ALL_METHODS.extend(m for m in METHODS if m not in L2.ALL_METHODS)


def check_noise(ndir):
    ref = {ln.split()[1]: ln.split()[0] for ln in open(SHA) if ln.strip()}
    h = hashlib.sha256(open(os.path.join(ndir, "noise.npy"), "rb").read()).hexdigest()
    ok = h == ref["t2m_out/noise/101231/noise.npy"]
    print(f"noise {ndir}: sha256 {h[:16]}... {'matches the registered array' if ok else 'DIFFERS from the registered array'}",
          flush=True)
    return ok


def run(ndir, pdir, nets, cells):
    import t2m_run
    c = json.load(open(CFG))
    check_noise(ndir)
    L2.configure(c["seed"], [ndir], c["j0"], c["nblocks"], c["nwin"], METHODS)
    nets = L2.CFG["nets"] if nets == "all" else [int(x) for x in nets.split(",")]
    os.makedirs(pdir, exist_ok=True)
    for j in nets:
        for cell in cells:
            t2m_run.run_net("REAL", cell, j, pdir)


def check_registered(nets):
    """MT@S0 of every part against the registered Task 2M part of the same (cell, network): same packed atoms."""
    n_same = n_cmp = 0
    for n in nets:
        fn = os.path.join(REG_PARTS, f"REAL_{n['cell']}_{n['net']}.json")
        if not os.path.exists(fn):
            continue
        with open(fn) as fh:
            ref = json.load(fh)
        n_cmp += 1
        same = len(ref["rows"]) == len(n["rows"])
        for r0, r1 in zip(ref["rows"], n["rows"]):
            a0, a1 = r0["atoms"]["MT@S0"], r1["atoms"]["MT@S0"]
            same &= len(a0) == len(a1) and all(np.allclose(np.asarray(u, float), np.asarray(v, float), rtol=1e-6,
                                                           atol=1e-6, equal_nan=True) for u, v in zip(a0, a1))
        n_same += bool(same)
    print(f"pipeline check: MT@S0 identical to the registered part in {n_same} of {n_cmp} parts compared")


def analyze(pdir, B):
    import c10_analyze as A
    import c16_analyze as Z
    import t2m_analyze as T
    nets, rows = Z.load(pdir, "REAL")
    print(f"M5-ACAP G-M5a  Task 2M seed 2718 REAL arm  {len(nets)} parts  B {B}  numpy {np.__version__}")
    print("gate: ACAP@S0 - MT@S0, DENSE LOC+T @0.25, 2.5th percentile > -0.03\n")
    check_registered(nets)
    names = [m for m in METHODS if m in rows[0]["atoms"]]
    verdict = None
    for des, rule, tol in (("DENSE", "LOC+T", 0.5), ("DENSE", "LOC", None), ("STD", "LOC+T", 0.5)):
        sub = [r for r in rows if r["design"] == des]
        recs = Z.records(sub, names, tol)
        for grp in ("all",) + L2.CELLS:
            sel = (lambda r: True) if grp == "all" else (lambda r, g=grp: r["cell"] == g)
            bt = A.Boot(sub, recs, sel, B if grp == "all" else max(500, B // 2))
            print(f"== {des} {rule} {grp}: {len(bt.idx)} windows, {int(bt.ks.sum())} events")
            for m in names:
                ci = bt.ci(m)
                print(f"   {m:9s} @0.1 {bt.pt[m][0]:.3f} [{ci[0][0]:.3f},{ci[0][1]:.3f}]  "
                      f"@0.25 {bt.pt[m][1]:.3f} [{ci[1][0]:.3f},{ci[1][1]:.3f}]")
            for a in ("ACAP@S0", "ACAP5@S0", "NOCAP@S0"):
                if a in names:
                    print(f"   {a} - MT@S0 @0.25 {T.f3(T.contrast(bt, a, 'MT@S0', 1))}")
            if (des, rule, grp) == ("DENSE", "LOC+T", "all"):
                verdict = T.contrast(bt, "ACAP@S0", "MT@S0", 1)["lo"] > -0.03
        print()
    print("G-M5a:", "PASS" if verdict else "FAIL", "(partial run)" if len(nets) < 72 else "")


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else "2800",
            sys.argv[5].split(",") if len(sys.argv) > 5 else list(L2.CELLS))
    elif sys.argv[1] == "analyze":
        analyze(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 2000)
    else:
        raise SystemExit(__doc__)
