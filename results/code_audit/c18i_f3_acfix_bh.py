"""c18i_f3_acfix_bh.py - re-run of the C18i detector at FC = 3 Hz with +-2 s station terms (c18i_widesc.py detect 2.0
SUBNET LOG6 3.0) but with the import order of c18h_detect.py (set_kernel BEFORE task1x_mtglrt is imported), so that the
deflation kernel AC is the one C18h used.  Only the import order differs from the stored C18i run.
    python c18i_f3_acfix_bh.py SUBNET OUT.json"""
import functools
import sys

sys.path[:0] = ["/home/claude/work/r0", "/home/claude/work/explore_contrib", "/home/claude/work"]
import c18h_kernel as K  # noqa: E402  (imports r0_lib and task1o_consistent only)

K.set_kernel(3.0)
import r0_run as RR  # noqa: E402
import station_corrections as SC  # noqa: E402
import task1x_mtglrt as X  # noqa: E402

RR.estimate_station_corrections = functools.partial(SC.estimate_station_corrections, search=2.0)
import c18d_detect as D  # noqa: E402

if __name__ == "__main__":
    c = len(X.AC) // 2
    print("AC(+-2) %.3f AC(+-4) %.3f" % (X.AC[c + 2], X.AC[c + 4]), flush=True)
    D.main(sys.argv[1], ["LOG6"], sys.argv[2])
