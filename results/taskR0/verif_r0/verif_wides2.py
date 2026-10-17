"""
verif_wides2.py - POST HOC sensitivity (not registered): frozen M2 with the first catalogue read of BOTH the +SC
estimator and the S2 joint estimation (c11_lib.P['search0']) widened from +-1.0 s to +-2.0 s (CAL catalogue only;
diag_corrections.py shows S delays truncated at the +-1 s edge).  Everything else as r0_run.Net; in-memory changes only.
    python verif_wides2.py SUBNET OUTDIR [CHUNK]
"""
import os
import sys

sys.path.insert(0, "/home/claude/work/r0")
os.chdir("/home/claude/work/r0")
import r0_lib as R  # noqa: E402,F401
import r0_run as RR  # noqa: E402
import station_corrections as SCm  # noqa: E402

RR.L11.P["search0"] = 2.0
_orig = SCm.estimate_station_corrections
RR.estimate_station_corrections = lambda STA, cal, search=1.0, gate=25.0, joint=True: _orig(STA, cal, 2.0, gate, joint)
import r0_run_chunks as RC  # noqa: E402

if __name__ == "__main__":
    RC.run("eval", sys.argv[1], ["M2"], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 50)
