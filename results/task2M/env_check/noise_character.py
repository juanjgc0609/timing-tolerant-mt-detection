"""Noise character per day (noise only, no injection, no detector): for every eligible 40 s slot and all stations of
the day, the matched-filter statistics of t2m_lib.snr_stats (excess kurtosis and max |MF| / robust sd)."""
import json, os, sys
import numpy as np
W = "/home/claude/work"
sys.path[:0] = [os.path.join(W, "explore_contrib"), W]
import xc
from task0_common import NT
for d in sys.argv[1:]:
    a = np.load(f"{W}/t2m_out/noise/{d}/noise.npy", mmap_mode="r")
    m = json.load(open(f"{W}/t2m_out/noise/{d}/meta.json"))
    ku, mx = [], []
    for s in m["eligible"]:
        i0 = int(round(s * 40 * 16))
        N = np.asarray(a[:, :, i0:i0 + NT], float)
        MF, _ = xc.matched(N)
        x = MF.ravel(); x = x / (np.median(np.abs(x)) / 0.6745)
        ku.append(np.mean(x ** 4) / np.mean(x ** 2) ** 2 - 3.0); mx.append(np.max(np.abs(x)))
    ku, mx = np.array(ku), np.array(mx)
    print(f"{d}: slots {len(ku)}; excess kurtosis median {np.median(ku):.2f} p90 {np.percentile(ku, 90):.1f}; "
          f"max|MF|/sd median {np.median(mx):.1f}; fraction of slots with max > 10 {np.mean(mx > 10):.2f}")
