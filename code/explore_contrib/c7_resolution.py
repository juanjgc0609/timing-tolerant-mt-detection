"""C7 mechanism diagnostic (seed 999, networks j = 40, 41 of each cell; sources from rng [999, 3706, ci, j]; new).
Noise-free, clutter-free single-event windows; first-iteration block statistic of each detector over the whole grid.
  EXACT  homogeneous medium (the detectors' own), no statics, no site/polarisation perturbation, no coda; design pulse
         family (fc U[1.4, 2.8], phase U[-45, 45]); zero corrections.
  REAL   default gen_realistic (gradient medium, site factors, polarisation deviation, coda) with the network's statics
         and the catalogue +SC corrections (xc.NetCfg), design pulse family; no noise, no clutter.
Reported per detector: fraction of windows whose global maximum block is > 5 km (epicentral) from the source (the
first atom would be a false atom plus a missed event), median of the "sidelobe ratio" max_{5-10 km} e / max_{<=5 km} e,
and the captured fraction (max_{<=5 km} e / whitened signal energy).
    python c7_resolution.py > out/c7_resolution.txt"""
import numpy as np

import c7_lib as C
import xc

GR = C.GR
T_ = xc.T
grid = xc.make_grid()
garr = np.asarray(grid)
BASES = {"Q1": C.Basis(C.basis_q1(), "Q1"), "PH": C.Basis(C.basis_ph(), "PH"), "PHF": C.Basis(C.basis_phf(4)[0], "PHF4")}
DETS = {"BOMP-PK": None, "PH": ("PH", "R1", False), "PH-FULL": ("PH", "FULL", False), "PHF": ("PHF", "R1", False),
        "CAP": ("Q1", "R1", True), "PH-CAP": ("PH", "R1", True)}
NSRC = 25


def run(setting):
    res = {d: {"far": [], "side": [], "capt": []} for d in DETS}
    for cell in T_.CELLS:
        for j in (40, 41):
            if setting == "EXACT":
                STA = T_.network(999, cell, j)[0]
                corr = np.zeros((len(STA), 2))
                statics = np.zeros((len(STA), 2))
                cfg = {"VS0": 3.0, "g": 0.0, "dir_deg": 0.0, "site_sd": 0.0, "coda_ratio": 0.0, "noise": 0.0}
            else:
                nc = xc.NetCfg(cell, j, "thin")
                STA, corr, statics = nc.STA, nc.corr, nc.statics
                cfg = {"noise": 0.0}
            ops = xc.Ops(STA, corr, grid)
            from certified_fw import stack_Q
            Qst = stack_Q(ops.Qs)
            QDs = {k: C.QDict(ops.BD, b, Qst) for k, b in BASES.items()}
            ci = T_.CELLS.index(cell)
            rng = np.random.default_rng([999, 3706, ci, j, 0 if setting == "EXACT" else 1])
            sd = GR.SNR * 0.02 * np.sqrt(C.E_R)
            for _ in range(NSRC):
                src = np.array([rng.uniform(-15, 15), rng.uniform(-15, 15), rng.uniform(5, 13)])
                tau = float(rng.uniform(5.5, 7.5))
                D, _p, _t, _a = GR.generate_real(STA, [src], [tau], 0, np.random.default_rng(rng.integers(1 << 30)),
                                                 np.random.default_rng(rng.integers(1 << 30)),
                                                 np.random.default_rng(rng.integers(1 << 30)), 0.0, cfg=cfg,
                                                 statics=statics)
                Esig = float(np.sum(D ** 2)) * C.E_R / sd ** 2
                MF, _ = xc.matched(D)
                MFn = MF / sd
                dist = np.linalg.norm(garr[:, :2] - src[None, :2], axis=1)
                for dn, spec in DETS.items():
                    if spec is None:
                        e = ops.PD.energy(MFn, True)[0]
                    else:
                        e = QDs[spec[0]].stat(BASES[spec[0]].filt(D, sd, MFn), spec[1], spec[2])[0]
                    en = e.max(axis=1)
                    near, ring = en[dist <= 5].max(), en[(dist > 5) & (dist <= 10)].max()
                    res[dn]["far"].append(float(dist[int(np.argmax(en))] > 5))
                    res[dn]["side"].append(float(ring / near))
                    res[dn]["capt"].append(float(near / Esig))
    print(f"{setting}: {len(res['BOMP-PK']['far'])} noise-free single-event windows (3 cells x 2 networks x {NSRC})")
    print("   detector     P(global max > 5 km)   sidelobe ratio median [P25, P75]   captured fraction median")
    for dn in DETS:
        r = res[dn]
        s = np.array(r["side"])
        print(f"   {dn:10s}   {np.mean(r['far']):8.3f}              {np.median(s):.3f} [{np.percentile(s, 25):.3f}, "
              f"{np.percentile(s, 75):.3f}]            {np.median(r['capt']):.3f}")
    b = np.array(res["BOMP-PK"]["side"])
    for dn in DETS:
        if dn != "BOMP-PK":
            print(f"   paired sidelobe ratio {dn} / BOMP-PK: median {np.median(np.array(res[dn]['side']) / b):.3f}, "
                  f"fraction > 1: {np.mean(np.array(res[dn]['side']) > b):.3f}")


if __name__ == "__main__":
    run("EXACT")
    run("REAL")
