"""
c18i_widesc.py - exploration C18i (Ridgecrest CAL only), after the independent verification of R0 found that the
station-term estimator searches only +-1.0 s around the homogeneous prediction while the real S delays reach 1.7 to
1.9 s.  The catalogue station corrections (+SC) of r0_run.Net are estimated here with a +-SEARCH s peak search
(everything else unchanged), and then
  oracle  the C18c oracle with the cap 6, LOG6 and no cap                      python c18i_widesc.py oracle SEARCH SUBNET OUT.json
  detect  the C18d detector with the variants given, optionally with the matched-filter kernel at FC Hz
                                                    python c18i_widesc.py detect SEARCH SUBNET VARIANTS FC OUT.json
"""
import functools
import sys

sys.path.insert(0, "/home/claude/work/r0")
import r0_lib as R  # noqa: E402,F401
import r0_run as RR  # noqa: E402
import station_corrections as SC  # noqa: E402

if __name__ == "__main__":
    mode, search = sys.argv[1], float(sys.argv[2])
    RR.estimate_station_corrections = functools.partial(SC.estimate_station_corrections, search=search)
    if mode == "oracle":
        import c18c_oracle as OC  # noqa: E402
        OC.VARIANTS = ("CAP6", "LOG6", "NOCAP")
        OC.main(sys.argv[3], sys.argv[4])
    else:
        fc = float(sys.argv[5])
        if fc != 2.0:
            import c18h_kernel as K  # noqa: E402
            K.set_kernel(fc)
        import c18d_detect as D  # noqa: E402
        D.main(sys.argv[3], sys.argv[4].split(","), sys.argv[6])
