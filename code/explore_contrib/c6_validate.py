"""C6 validation: the model-generic operators reduce to the registered ones for the homogeneous model, the gradient
operators reproduce the generator, the dense MT-FW reproduces certified_fw_record, and the JHD calibration runs.
Seed 999, networks 10-11 (pilot networks already used by C1-C3, not the C6 evaluation networks)."""
import sys
import time

import numpy as np

import c6_lib as L
import xc
import gen_realistic as GR
from station_corrections import estimate_station_corrections
from task1w_mech import rand_dc

ok = True


def chk(name, val, tol):
    global ok
    good = bool(val <= tol)
    ok &= good
    print(f"  {name:<70s} {val:.3e} (<= {tol:g}) {'ok' if good else 'FAIL'}", flush=True)


rng = np.random.default_rng(5)
STA = np.array([[20.0, 5.0, 0.0], [-30.0, 12.0, 0.0], [5.0, -40.0, 0.0], [44.0, 30.0, 0.0]])
srcs = [np.array([rng.uniform(-15, 15), rng.uniform(-15, 15), rng.uniform(5, 13)]) for _ in range(20)]
# 1. homogeneous model == registered
H = L.Homog()
e1 = max(np.abs(H.tt(s[None], STA)[0] - np.stack([np.linalg.norm(STA - s, axis=1) / xc.VP,
                                                    np.linalg.norm(STA - s, axis=1) / xc.VS])).max() for s in srcs)
chk("Homog tt vs d/VP, d/VS", e1, 1e-12)
e2 = max(np.abs(L.radiation_operator_m(s, STA, H) - xc.X.radiation_operator(s, STA)).max() for s in srcs)
chk("Homog radiation operator vs task1x radiation_operator", e2, 1e-12)
# 2. gradient model == generator
Gm = L.Gradient()
e3 = max(max(np.abs(Gm.tt(s[None], STA)[0, b] - GR.travel_time_gradient(s, STA, "PS"[b])).max() for b in range(2))
         for s in srcs)
chk("Gradient tt vs gen_realistic.travel_time_gradient", e3, 1e-12)
e4 = 0.0
for s in srcs:
    Mt = rand_dc(rng)
    m = np.array([Mt[0, 0], Mt[1, 1], Mt[2, 2], Mt[0, 1], Mt[0, 2], Mt[1, 2]])
    G = L.radiation_operator_m(s, STA, Gm).reshape(2, len(STA), 3, 6)
    for b, A in ((0, xc.X.AP), (1, xc.X.AS)):
        for i in range(len(STA)):
            ref = GR.radiation_gradient(Mt, s, STA[i], "PS"[b]) * A / np.linalg.norm(STA[i] - s)
            e4 = max(e4, float(np.abs(G[b, i] @ m - ref).max()))
chk("Gradient radiation operator . m vs gen_realistic.radiation_gradient", e4, 1e-12)
# 3. two-layer: v1 = v2 -> straight; Snell consistency dT/dX = sin(i_r)/v1
TL = L.TwoLayer(3.0, 3.0, 5.0, 1.73)
e5 = max(np.abs(TL.tt(s[None], STA) - H.tt(s[None], STA)).max() for s in srcs)
chk("TwoLayer(v1 = v2) tt vs homogeneous", e5, 1e-9)
TL = L.TwoLayer(2.9, 3.5, 4.0, 1.73)
e6 = 0.0
for s in srcs:
    Xh = np.array([10.0, 25.0, 40.0])
    for b in range(2):
        t, _ss, sr = TL._solve(Xh, np.full(3, s[2]), b)
        tp, _a, _b = TL._solve(Xh + 1e-5, np.full(3, s[2]), b)
        v1 = TL.v1 * (TL.vpvs if b == 0 else 1.0)
        e6 = max(e6, float(np.abs((tp - t) / 1e-5 - sr / v1).max()))
chk("TwoLayer dT/dX vs ray parameter sin(i_r)/v1", e6, 1e-5)

# 4. operators on a real network: ModelOps(Homog) == xc.Ops; dense MT-FW == certified_fw_record; corr_model == SC
t0 = time.time()
nc = xc.NetCfg("ladder_6", 10, "thin")
print(f"  NetCfg built in {time.time() - t0:.1f} s", flush=True)
grid = xc.make_grid()
t0 = time.time()
ops0 = xc.Ops(nc.STA, nc.corr, grid)
t_ops0 = time.time() - t0
Ec = L.cat_energy(nc.cat)
cH = L.corr_model(nc.STA, nc.cat, H, Ec)
chk("corr_model(Homog) vs estimate_station_corrections", float(np.abs(cH - nc.corr).max()), 1e-12)
t0 = time.time()
opsH = L.ModelOps(nc.STA, grid, H, nc.corr)
print(f"  ModelOps {time.time() - t0:.2f} s (xc.Ops {t_ops0:.2f} s)", flush=True)
chk("BlockDict idx (ModelOps Homog vs registered)", float(np.abs(opsH.BD.idx - ops0.BD.idx).max()), 0)
chk("BlockDict Kinv, max abs diff / max |Kinv| (operators differ by ~1e-17)", float(np.abs(opsH.BD.Kinv - ops0.BD.Kinv).max() / np.abs(ops0.BD.Kinv).max()), 1e-6)
st = xc.window_stream(nc, 3906)
dmt = dbo = 0.0
tb = tm = tc = 0.0
for _ in range(8):
    w = next(st)
    t = time.time(); r0, l0 = xc.bomp_pk(ops0.PD, w["MFn"], w["truth"]); tb -= time.time() - t
    t = time.time(); r1, l1 = L.bomp_pk_m(opsH, w["MFn"], w["truth"]); tb += time.time() - t
    dbo = max(dbo, max(abs(a[0] - b[0]) + abs(a[1] - b[1]) for a, b in zip(r0, r1)), float(len(r0) != len(r1)))
    t = time.time(); r2 = xc.mtfw(ops0, w["MFn"], w["truth"]); tc += time.time() - t
    t = time.time(); r3, l3 = L.mtfw_dense(opsH, w["MFn"], w["truth"]); tm += time.time() - t
    dmt = max(dmt, max(abs(a[0] - b[0]) + abs(a[1] - b[1]) for a, b in zip(r2, r3)), float(len(r2) != len(r3)))
chk("BOMP-PK records (ModelOps Homog vs registered), 8 windows", dbo, 1e-9)
chk("MT-FW records (dense vs certified_fw_record), 8 windows", dmt, 1e-9)
print(f"  s/window: dense MT-FW {tm / 8:.3f}, certified MT-FW {tc / 8:.3f}", flush=True)
# 5. JHD calibration on this network
for fam in ("homog", "gradient", "twolayer"):
    t = time.time()
    mdl, info = L.jhd_fit(nc.STA, nc.cat, fam, Ec, return_info=True)
    print(f"  JHD {fam}: {info}  ({time.time() - t:.1f} s)", flush=True)
print("VALIDATION", "PASS" if ok else "FAIL")
