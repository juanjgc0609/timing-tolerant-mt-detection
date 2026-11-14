"""Close pairs in the SCSN catalogue of Ridgecrest, 2019-07-04 + 30 days, box 35.25-36.25 N, 118.1-117.1 W.
Catalogue fields only. Input: scedc-pds earthquake_catalogs/SCEC_DC/2019.catalog. Usage: python -I this.py 2019.catalog
Origin-time differences, not arrival separations (the CR-04 estimand uses predicted arrivals)."""
import datetime as dt
import sys

import numpy as np

t, la, lo = [], [], []
for line in open(sys.argv[1]):
    if line.startswith('#') or len(line) < 60:
        continue
    f = line.split()
    if f[2] != 'eq' or '.' not in f[1]:
        continue
    d = dt.datetime.strptime(f[0] + ' ' + f[1][:11], '%Y/%m/%d %H:%M:%S.%f').replace(tzinfo=dt.timezone.utc)
    t.append(d.timestamp()); la.append(float(f[6])); lo.append(float(f[7]))
t, la, lo = map(np.array, (t, la, lo))
a = dt.datetime(2019, 7, 4, tzinfo=dt.timezone.utc).timestamp()
s = (t >= a) & (t < a + 30 * 86400) & (abs(la - 35.75) < 0.5) & (abs(lo + 117.6) < 0.5)
o = np.argsort(t[s]); T = t[s][o]
X = (lo[s][o] + 117.6) * 111.19 * np.cos(np.radians(35.75)); Y = (la[s][o] - 35.75) * 111.19
d, r = np.diff(T), np.hypot(np.diff(X), np.diff(Y))
print('events', s.sum())
for th in (0.1, 0.5, 1.0, 2.0):
    k = int(((d < th) & (r < 10)).sum())
    print(f'next event within {th} s and 10 km: {k} ({100 * k / s.sum():.2f} %)')
obs = np.zeros(T.size, bool); obs[1:] |= d < 1; obs[:-1] |= d < 1
lam = (np.searchsorted(T, T + 300) - np.searchsorted(T, T - 300) - 1) / 600.0
print(f'share with another event within 1 s, any distance: observed {100 * obs.mean():.2f} %, '
      f'Poisson at the local 10-min rate {100 * (1 - np.exp(-2 * lam)).mean():.1f} %')
