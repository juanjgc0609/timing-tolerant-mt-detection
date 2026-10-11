import sys, time
import r0_run as RR, r0_lib as R, numpy as np, xc
net = RR.Net('N8', 940.0)
ws = RR.starts('cal')[100:103]
for m in sys.argv[1].split(','):
    t = time.time()
    n = 0
    for ts in ws:
        W = R.window(net.d16, ts); MF, sd = xc.matched(W)
        n += len(net.run(m, MF / sd))
    print(m, 'sec/window', round((time.time() - t) / len(ws), 2), 'atoms', n, flush=True)
