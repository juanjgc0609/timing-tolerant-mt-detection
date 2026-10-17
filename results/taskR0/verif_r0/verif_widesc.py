"""
verif_widesc.py - POST HOC sensitivity (not part of the registered test): the frozen M2 run through r0_run_chunks.run
with ONE change, the +SC station corrections estimated on the CAL catalogue with a +-2.0 s peak search instead of
the frozen +-1.0 s (diag_corrections.py shows the S delays truncated at the +-1 s edge).  S2 terms, sigma and engines
are then built from those corrections exactly as in r0_run.Net.  Frozen files are not modified (in-memory wrap).
    python verif_widesc.py SUBNET OUTDIR [CHUNK]
"""
import os
import sys

sys.path.insert(0, "/home/claude/work/r0")
os.chdir("/home/claude/work/r0")
import r0_lib as R  # noqa: E402,F401
import r0_run as RR  # noqa: E402
import station_corrections as SCm  # noqa: E402

_orig = SCm.estimate_station_corrections


def wide(STA, cal, search=1.0, gate=25.0, joint=True):
    return _orig(STA, cal, search=2.0, gate=gate, joint=joint)


RR.estimate_station_corrections = wide
import r0_run_chunks as RC  # noqa: E402

if __name__ == "__main__":
    RC.run("eval", sys.argv[1], ["M2"], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 50)
