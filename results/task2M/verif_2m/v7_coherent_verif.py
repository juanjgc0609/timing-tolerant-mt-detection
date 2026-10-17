"""V7: (a) network-coherent transients (possible uncatalogued real earthquakes) among the slots used by the 2718 DENSE
windows, and the contrast without those windows; (b) BOMP score scale vs transient amplitude (numerical sanity)."""
import sys
import numpy as np
sys.path.insert(0, "/home/claude/work/t2m_verif")
import vlib_verif as V  # noqa: E402
z = np.load("/home/claude/work/t2m_verif/v3_station_slotstats_verif.npz")
mx, rs, el = z["mx"], z["rs"], z["el"]
r = mx / np.median(rs, axis=1)[:, None]
coh = set(el[(r > 10).sum(0) >= 5].tolist())
out = [f"(a) eligible slots of 101231 with >= 5 of 10 stations above 10x their robust level: {len(coh)}"]
wins = V.windows(V.load_parts("2718"), "DENSE")
bad = np.array([w["slot"] in coh for w in wins])
out.append(f"    DENSE windows on those slots: {bad.sum()} of {len(wins)}")
ws = np.nonzero(~bad)[0]
p = {m: V.recall_at(*V.prep(wins, m), ws)[1] for m in ("MT@S0", "HUBER-N@S0", "BOMP@S0")}
out.append(f"    without them: MT-BOMP {p['MT@S0'] - p['BOMP@S0']:+.3f}, HUB-BOMP {p['HUBER-N@S0'] - p['BOMP@S0']:+.3f}")
# (b) BOMP first-atom score vs maxabs^2 of the window noise
s1 = np.array([w["atoms"]["BOMP@S0"][0][4] for w in wins])
m1 = np.array([w["atoms"]["MT@S0"][0][4] for w in wins])
mxw = np.array([w["maxabs"] for w in wins])
hi = mxw > 30
out.append(f"(b) windows with noise max|MF|/sd > 30: {hi.sum()}; BOMP first-atom score / max^2: median "
           f"{np.median(s1[hi] / mxw[hi] ** 2):.3f}, max {np.max(s1[hi] / mxw[hi] ** 2):.3f}; "
           f"Spearman-type rank corr(score, max) {np.corrcoef(np.argsort(np.argsort(s1[hi])), np.argsort(np.argsort(mxw[hi])))[0, 1]:.2f}")
out.append(f"    MT@S0 first-atom score in those windows: median {np.median(m1[hi]):.0f}, max {m1[hi].max():.0f} "
           f"(all windows median {np.median(m1):.0f}); BOMP median {np.median(s1[hi]):.0f}, max {s1[hi].max():.0f} "
           f"(all windows median {np.median(s1):.0f})")
i = int(np.argmax(s1))
out.append(f"    largest BOMP score {s1[i]:.0f} in window with max|MF|/sd {mxw[i]:.0f} (ratio {s1[i] / mxw[i] ** 2:.2f})")
open("/home/claude/work/t2m_verif/v7_coherent_verif.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
