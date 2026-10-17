"""c18h_detect.py - exploration C18h, detector level on Ridgecrest CAL: the full loop of c18d_detect.py with the
matched-filter kernel at centre frequency FC (c18h_kernel.set_kernel) and the statistics given.
    python c18h_detect.py FC SUBNET VARIANTS OUT.json"""
import sys

import c18h_kernel as K

if __name__ == "__main__":
    K.set_kernel(float(sys.argv[1]))
    import c18d_detect as D  # noqa: E402
    D.main(sys.argv[2], sys.argv[3].split(","), sys.argv[4])
