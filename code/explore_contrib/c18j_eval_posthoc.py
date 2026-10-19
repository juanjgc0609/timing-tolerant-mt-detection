"""c18j_eval_posthoc.py - declared post hoc (error register E59): the C18j configuration read out of sample.

C18j estimated its layered station terms on the 16 Ridgecrest CAL events and counted the same 16 events as detected
(E53). This script runs c18j_layered.main unchanged, with the terms still estimated on CAL only, on the windows of the
R0 evaluation segment [900, 3560) s, whose events never enter the terms. The R0 evaluation segment has already been
used by the registered R0 test, so this is a descriptive reading of used data, not a test.

Two changes, made in memory: r0_run.starts returns the EVAL windows, and the waveforms are read to the end of the hour
(3600 s, as the R0 evaluation run) instead of 940 s.
    python c18j_eval_posthoc.py SUBNET VARIANTS OUT.json          (fc 2 Hz, as the C18j gate configuration)
    python c18j_eval_posthoc.py SUBNET VARIANTS OUT.json cal      (control: the CAL windows with t_end 3600)
"""
import sys

import c18j_layered as J

RR = J.RR
_starts = RR.starts
_Net = RR.Net


def main(sub, variants, outp, seg="eval"):
    RR.starts = lambda _segment: _starts(seg)
    RR.Net = lambda s, _t_end: _Net(s, 3600.0)
    J.main("detect", 2.0, sub, [variants, outp])


if __name__ == "__main__":
    main(*sys.argv[1:5])
