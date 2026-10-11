"""symmetry_verif.py - integrity / symmetry checks on the stored parts."""
import numpy as np, json, os, glob, collections
from load_verif import *
W, N = load_all()
print("# parts metadata")
print("tags:", collections.Counter((v["tags"]["std"], v["tags"]["dense"]) for v in N.values()))
print("kmax:", collections.Counter((v["kmax"]["std"], v["kmax"]["dense"]) for v in N.values()))
print("methods std by net range:", collections.Counter((min(k[3] // 30, 2) if k[3] < 60 else 2 + (k[3] >= 100), tuple(v["methods"]["std"])) for k, v in N.items()).most_common())
print("methods dense:", collections.Counter((k[3] < 30, k[3] < 100, tuple(v["methods"]["dense"])) for k, v in N.items()).most_common())
ker = collections.Counter()
for f in glob.glob(os.path.join(PARTS, "main_*.json")):
    e = json.load(open(f)); ker[(e["kernel"]["OPENBLAS_CORETYPE"], e["kernel"]["blas"][0][0], e["arm"], e["task"])] += 1
print("kernel/arm:", ker)
print("\n# per arm: atoms per window, origin range, min pairwise epicentral distance, grid membership")
for kind in ("std", "dense"):
    for m in ("MT@S0", "MT@S0+T3", "MT@S2+T3", "BOMP@S0+T3", "BOMP@S2+T3"):
        n, tmin, tmax, dmin, zs, fx = collections.Counter(), 9, 0, 99, set(), 0
        for w in W:
            if w["kind"] != kind or m not in w["atoms"]:
                continue
            a = w["atoms"][m]; n[len(a)] += 1
            tmin, tmax = min(tmin, a[:, 3].min()), max(tmax, a[:, 3].max()); zs |= set(a[:, 2].tolist())
            d = np.hypot(a[:, None, 0] - a[None, :, 0], a[:, None, 1] - a[None, :, 1]) + 99 * np.eye(len(a)); dmin = min(dmin, d.min())
            fx = max(fx, np.abs((a[:, :2] - 0.37) - np.round(a[:, :2] - 0.37)).max())
        if n:
            print(f"{kind:5s} {m:11s}: atoms/window {dict(n)}; origins [{tmin:.2f}, {tmax:.2f}]; depths {sorted(zs)}; min pairwise epicentral distance {dmin:.2f} km; max deviation from the 1 km lattice {fx:.1e}")
print("\n# T3 definition: first atom of MT@S0+T3 vs first atom of MT@S0 (same hard selection on the same engine and data)")
dd, dtm, same_s = [], [], 0
for w in W:
    a, b = w["atoms"]["MT@S0"], w["atoms"]["MT@S0+T3"]
    dd.append(np.hypot(*(a[0, :2] - b[0, :2]))); dtm.append(abs(a[0, 3] - b[0, 3])); same_s += a[0, 4] == b[0, 4]
print(f"max epicentral move {max(dd):.2f} km (limit 10), max origin move {max(dtm):.2f} s (limit 1), first-atom score identical in {same_s}/{len(W)} windows, moved in {np.mean(np.array(dd) > 0):.3f}")
print("\n# T3 timing level (sigma_catalogue recomputed with each term set): S0 vs S2")
for cell in CELLS:
    for reg in REGIMES:
        s0 = np.array([v["sigT3"]["S0"][0] for k, v in N.items() if k[1] == reg and k[2] == cell]); s2 = np.array([v["sigT3"]["S2"][0] for k, v in N.items() if k[1] == reg and k[2] == cell])
        print(f"{cell:15s} {reg:6s}: sig S0 median {np.median(s0):.3f} [q10 {np.quantile(s0,.1):.3f}, q90 {np.quantile(s0,.9):.3f}]  sig S2 median {np.median(s2):.3f} [q10 {np.quantile(s2,.1):.3f}, q90 {np.quantile(s2,.9):.3f}]  S2 < S0 in {np.mean(s2 < s0):.2f}")
print("\n# host_time_utc of main parts: earliest", min(v["host_time_utc"] for v in N.values()), "latest", max(v["host_time_utc"] for v in N.values()))
