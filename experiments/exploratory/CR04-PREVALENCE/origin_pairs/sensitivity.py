"""Sensitivity of the catalogue bound (ridgecrest_origin_pairs.py) to the rate window, the magnitude cut and the start
day. Same input and box. Usage: python -I sensitivity.py 2019.catalog"""
import datetime as dt
import sys

import numpy as np

t, la, lo, M = [], [], [], []
for line in open(sys.argv[1]):
    if line.startswith('#') or len(line) < 60:
        continue
    f = line.split()
    if f[2] != 'eq' or '.' not in f[1]:
        continue
    d = dt.datetime.strptime(f[0] + ' ' + f[1][:11], '%Y/%m/%d %H:%M:%S.%f').replace(tzinfo=dt.timezone.utc)
    t.append(d.timestamp()); la.append(float(f[6])); lo.append(float(f[7])); M.append(float(f[4]))
t, la, lo, M = map(np.array, (t, la, lo, M))
a = dt.datetime(2019, 7, 4, tzinfo=dt.timezone.utc).timestamp()
box = (t < a + 30 * 86400) & (abs(la - 35.75) < 0.5) & (abs(lo + 117.6) < 0.5)


def shares(T, half_window):
    T = np.sort(T)
    d = np.diff(T)
    obs = np.zeros(T.size, bool); obs[1:] |= d < 1; obs[:-1] |= d < 1
    lam = (np.searchsorted(T, T + half_window) - np.searchsorted(T, T - half_window) - 1) / (2 * half_window)
    return 100 * obs.mean(), 100 * (1 - np.exp(-2 * lam)).mean()


print("variant\tn\tobserved %\tPoisson %\tratio")
for name, start, mmin, hw in ([("window 60 s", a, -9, 30), ("window 300 s", a, -9, 150),
                               ("window 600 s (main)", a, -9, 300), ("window 1800 s", a, -9, 900),
                               ("window 3600 s", a, -9, 1800)]
                              + [(f"M >= {m}", a, m, 300) for m in (1.5, 2.0, 2.5)]
                              + [("from 9 July", a + 5 * 86400, -9, 300)]):
    T = t[box & (t >= start) & (M >= mmin)]
    o, p = shares(T, hw)
    print(f"{name}\t{T.size}\t{o:.2f}\t{p:.2f}\t{o / p:.2f}")
