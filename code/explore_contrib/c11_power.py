"""C11 power calculation for the draft preregistration: normal approximation with the evaluation-set bootstrap SD of a
paired difference, scaled by sqrt(n_pilot / n) networks per cell and regime.
    python c11_power.py"""
from math import erf, sqrt


def Phi(x):
    return 0.5 * (1 + erf(x / sqrt(2)))


def power_sup(delta, sd, z=1.96):
    """P(lower one-sided 97.5% bound > 0)"""
    return Phi(delta / sd - z)


def power_point(delta, sd, thr):
    return Phi((delta - thr) / sd)


def power_ni(delta, sd, margin, z=1.645):
    """P(lower one-sided 95% bound > -margin)"""
    return Phi((delta + margin) / sd - z)


def table(name, sd_pilot, n_pilot, deltas, ns, kind, margin=0.02, thr=0.03):
    print(name, f"(pilot SD {sd_pilot:.4f} at {n_pilot} networks per cell and regime)")
    for n in ns:
        sd = sd_pilot * sqrt(n_pilot / n)
        if kind == "sup":
            s = "; ".join(f"true {d:+.3f}: bound {power_sup(d, sd):.2f}, point>={thr} {power_point(d, sd, thr):.2f}" for d in deltas)
        else:
            s = "; ".join(f"true {d:+.3f}: {power_ni(d, sd, margin):.2f}" for d in deltas)
        print(f"   n = {n:4d}  SD {sd:.4f}  {s}")


if __name__ == "__main__":
    import sys
    args = [float(x) for x in sys.argv[1:]]
    # usage: sd_primary_std sd_ni_std_loc sd_ni_dense_loct sd_ni_dense_loc
    n0 = 30
    table("superiority, STD LOC+T @0.25", args[0], n0, (0.03, 0.04, 0.05, 0.06), (60, 100, 160, 240), "sup")
    for nm, sd in zip(("STD LOC", "DENSE LOC+T", "DENSE LOC"), args[1:]):
        table(f"non-inferiority margin 0.02, {nm} @0.25", sd, n0, (-0.01, 0.0, 0.01, 0.02), (60, 100, 160, 240), "ni")
