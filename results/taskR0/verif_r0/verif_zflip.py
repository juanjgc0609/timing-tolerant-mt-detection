"""
verif_zflip.py - POST HOC sensitivity (not part of the registered test): the frozen M2 run through r0_run_chunks.run
with ONE change, the vertical component fed to the detector as +Z (up) instead of the registered -Z.  The frozen
files are not modified: r0_lib.load_waveforms is wrapped in memory before r0_run.Net is built.
    python verif_zflip.py SUBNET OUTDIR [CHUNK]     then r0_run_chunks.py merge ... as usual
"""
import os
import sys

sys.path.insert(0, "/home/claude/work/r0")
os.chdir("/home/claude/work/r0")
import r0_lib as R  # noqa: E402

_orig = R.load_waveforms


def load_flipped(names, sens, t_end):
    raw, d16 = _orig(names, sens, t_end)
    d16 = d16.copy()
    d16[:, 2] *= -1.0                                   # undo the registered -Z: detector sees +Z (up)
    return raw, d16


R.load_waveforms = load_flipped
import r0_run_chunks as RC  # noqa: E402

if __name__ == "__main__":
    RC.run("eval", sys.argv[1], ["M2"], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 50)
