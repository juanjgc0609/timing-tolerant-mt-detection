# development (CAL only): matched-filter station terms with the reading search bounded by half the predicted S-P time
import sys, numpy as np, json
sys.path.insert(0, '/home/claude/work/r1')
import r1_lib as R; R.setup()
import c18_tt as TT, c18j_layered as LJ
from station_corrections import peak_delays
from task0_common import VP, VS
net = R.make_net('cal')
def corr(frac, cap=2.0, gate=25.0):
    ds, ens = [], []
    for D, src, tau in net.calwin:
        src = np.asarray(src, float)
        tl = TT.times(TT.JMA2001, 1.73, src, net.STA)            # (2, ns)
        d = np.linalg.norm(net.STA - src, axis=1)
        th = np.stack([d / VP, d / VS], axis=1)
        dl, en = np.zeros((len(d), 2)), np.zeros((len(d), 2))
        for i in range(len(d)):
            s = min(cap, frac * (tl[1, i] - tl[0, i])) if frac else cap
            a, b = peak_delays(net.STA, D, src, float(tau), s)
            dl[i], en[i] = a[i], b[i]
        ds.append(dl + th - tl.T); ens.append(en)
    ds, ens = np.array(ds), np.array(ens)
    ok = ens >= gate
    c = np.zeros(ds.shape[1:]); sd = np.zeros(ds.shape[1:])
    for i in range(ds.shape[1]):
        for b in range(2):
            if ok[:, i, b].any():
                c[i, b] = np.median(ds[ok[:, i, b], i, b]); sd[i, b] = np.median(np.abs(ds[ok[:, i, b], i, b] - c[i, b]))
    return c, sd
pk = np.array(json.load(open('/home/claude/work/r1/dev/picks_cal.json'))['pick_terms'])
for frac in (None, 0.5):
    c, sd = corr(frac)
    print('frac', frac)
    print('  MF  P', np.round(c[:, 0], 2), ' S', np.round(c[:, 1], 2))
    print('  MAD P', np.round(sd[:, 0], 2), ' S', np.round(sd[:, 1], 2))
    print('  MF-pick P', np.round(c[:, 0] - pk[:, 0], 2), ' S', np.round(c[:, 1] - pk[:, 1], 2))
