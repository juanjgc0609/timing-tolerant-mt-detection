"""ac_order_bh.py - which deflation / Gram kernel (task1x_mtglrt.AC) is in force when the matched-filter kernel is
switched to FC = 3 Hz, for the import orders of c18h_detect.py, c18i_widesc.py (detect) and c18j_layered.py.
Each case runs in a fresh interpreter and reproduces only the import statements and the set_kernel call of the script.
Also: the matched-filter response of a 2 Hz Ricker arrival for each kernel, and the 3 Hz kernel response to a 3 Hz Ricker,
to show what shape each AC assumes."""
import subprocess
import sys

PY = sys.executable
PRE = ("import sys, numpy as np; sys.path[:0]=['/home/claude/work/r0','/home/claude/work/explore_contrib','/home/claude/work'];")
CASES = {
    "c18h_detect (set_kernel, then import c18d_detect)":
        PRE + "import c18h_kernel as K; K.set_kernel(3.0); import c18d_detect as D;",
    "c18i_widesc detect FC=3 (import r0_run, then set_kernel)":
        PRE + "import r0_lib as R, r0_run as RR, station_corrections as SC; import c18h_kernel as K; K.set_kernel(3.0); "
              "import c18d_detect as D;",
    "c18j_layered FC=3 (import r0_run at top, set_kernel in main)":
        PRE + "import r0_lib as R, r0_run as RR, c18_tt as TT; import c18h_kernel as K; K.set_kernel(3.0); "
              "import c10_lib as C, xc;",
    "reference 2 Hz (no set_kernel)":
        PRE + "import r0_lib as R, r0_run as RR;",
}
TAIL = ("import task1x_mtglrt as X, task1o_consistent as T1O, task2c_methods as M; from task0_common import DT;"
        "ac=X.AC; c=len(ac)//2; half=lambda a: int(np.sum(a>=0.5));"
        "print('  X.AC: samples >= 0.5 of peak', half(ac), ' AC(+-2) %.3f AC(+-4) %.3f AC(+-6) %.3f' % (ac[c+2], ac[c+4], ac[c+6]),"
        "'| M.AC is X.AC:', M.AC is X.AC, '| kernel in force for matched(): fc=%s' % ('3' if abs(T1O._KER[len(T1O._KER)//2+2] - (1-2*(np.pi*3*2*DT)**2)*np.exp(-(np.pi*3*2*DT)**2))<1e-12 else '2'))")


def main():
    lines = []
    for name, code in CASES.items():
        r = subprocess.run([PY, "-c", code + TAIL], capture_output=True, text=True)
        lines.append(name)
        lines.append(r.stdout.rstrip() or r.stderr[-800:])
    code = "\n".join(x.strip() for x in PRE.split(";")) + """
from task0_common import DT
rk = lambda t, f: (1 - 2 * (np.pi * f * t) ** 2) * np.exp(-(np.pi * f * t) ** 2)
t = np.arange(-3, 3 + 1e-9, DT)
tk = np.arange(-1, 1 + 1e-9, DT)
k3, k2 = rk(tk, 3.0), rk(tk, 2.0)
def ac(w, k):
    a = np.convolve(w, k[::-1], mode="same")
    return a / a.max()
for nm, a in (("2Hz wavelet * 2Hz kernel", ac(rk(t, 2), k2)), ("2Hz wavelet * 3Hz kernel", ac(rk(t, 2), k3)),
              ("3Hz wavelet * 3Hz kernel", ac(rk(t, 3), k3))):
    c = len(a) // 2
    print("  %-26s width(>=0.5) %d samples  a(+-1) %.3f a(+-2) %.3f a(+-4) %.3f a(+-6) %.3f min %.3f" % (
        nm, np.sum(a >= 0.5), a[c + 1], a[c + 2], a[c + 4], a[c + 6], a.min()))
"""
    r = subprocess.run([PY, "-c", code], capture_output=True, text=True)
    lines.append("shapes:")
    lines.append(r.stdout.rstrip() or r.stderr[-800:])
    print("\n".join(lines))


if __name__ == "__main__":
    main()
