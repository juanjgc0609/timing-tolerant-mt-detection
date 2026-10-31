"""Validation of the synthetic generator v2 and of the semi-synthetic sets (report: results/synth_v2/).

    TTMT_DATA_ROOT=... TTMT_WORK_ROOT=... python tools/validate_synth.py [--fk-lib PKL ...] > results/synth_v2/validation_out.txt

1. engine against the closed-form whole-space solution (Aki and Richards 4.29), six tensor components, three distances,
   and its convergence with the period L of the wavenumber sum;
2. free surface: vertical-incidence P doubles at the surface;
3. borehole receiver: the free-surface ghost at 2 zr / vp with the expected sign and spreading (CR-09);
4. layered medium: engine onsets against the own ray tracer and against the detector's (independent) tracer;
5. inverse-crime distance: the fraction of the noise-free arrival energy of v2 events inside the detector's 6-D
   moment-tensor subspace at the true source, seen through the detector's own processing (16 Hz, 0.5-6 Hz, matched
   filter at the true arrival times), against C19's 0.81 for the v1 generator;
6. EGF libraries and sets: counts, pick quality, SNR, separation, CC to atoms (CR-06), K_mec/E (CR-14), partitions.
"""
import argparse
import glob
import json
import pickle
import time

from ttmt.runtime import pin_blas

pin_blas()

import numpy as np  # noqa: E402

from ttmt.synth import analytic, fk, media, sets  # noqa: E402


def line(s=""):
    print(s, flush=True)


def check(name, val, tol, ok_if="lt"):
    ok = val < tol if ok_if == "lt" else val > tol
    line(f"   {name:<72} {val:10.3e}  ({'<' if ok_if == 'lt' else '>'} {tol:g})  {'ok' if ok else 'FAIL'}")
    return ok


def s1_whole_space():
    line("1. Engine against the whole-space solution (vp 6, vs 3.5, rho 2.7; source 5 km deep; Gaussian moment rate, "
         "f_gauss 4 Hz)")
    m = media.whole_space(6.0, 3.5, 2.7)
    r, phi = np.array([3.0, 8.0, 15.0]), np.array([0.3, 2.0, -1.1])
    t = np.arange(800) / 25.0
    A = [analytic.whole_space_velocity([0, 0, 5.0], [r[i] * np.cos(phi[i]), r[i] * np.sin(phi[i]), 0.0],
                                       6.0, 3.5, 2.7, t, 4.0, 1.5) for i in range(3)]
    ok = True
    for Lf in (1.5, 3.0, 6.0):
        cfg = fk.FKConfig(fs=25, nt=800, f_gauss=4.0, t_shift=1.5, L_factor=Lf)
        t0 = time.time()
        x = fk.greens(m, 5.0, r, phi, cfg, free_surface=False)
        e = np.array([[np.max(np.abs(x[i, j] - A[i][j])) / np.max(np.abs(A[i][j])) for j in range(6)]
                      for i in range(3)])
        line(f"   L factor {Lf}: max relative error per distance {np.round(e.max(1), 5).tolist()} "
             f"(worst component {fk.MT_ORDER[int(e.max(0).argmax())]}), {time.time() - t0:.0f} s")
        if Lf == 3.0:
            ok &= check("default L factor 3: worst relative error, 6 components x 3 distances", float(e.max()), 0.01)
    return ok


def _peaks(x, fs, windows):
    return [float(np.max(np.abs(x[int(a * fs):int(b * fs)]))) for a, b in windows]


def s2_free_surface_and_ghost():
    line("2. Free surface: vertical-incidence P on a half-space (source Mzz below the receiver). The ratio tends to 2 as")
    line("   the near- and intermediate-field terms (relative size ~ wavelength / depth) vanish")
    m = media.whole_space(6.0, 3.5, 2.7)
    ok = True
    for zs in (10.0, 30.0):
        cfg = fk.FKConfig(fs=50, nt=1500, f_gauss=6.0, t_shift=1.0, L_factor=3.0)
        w_fs = fk.greens(m, zs, [0.01], [0.0], cfg, free_surface=True)[0, 2, 2]
        w_ws = fk.greens(m, zs, [0.01], [0.0], cfg, free_surface=False)[0, 2, 2]
        tp = 1.0 + zs / 6.0
        sl = slice(int((tp - 0.3) * 50), int((tp + 0.3) * 50))
        a, b = np.linalg.norm(w_fs[sl]), np.linalg.norm(w_ws[sl])
        line(f"   depth {zs:4.0f} km: surface / whole-space P pulse (L2 over +-0.3 s) = {a / b:.4f}")
    ok &= check("|ratio - 2| at 30 km", abs(a / b - 2.0), 0.005)
    line("3. Borehole receiver 2 km deep (CR-09): direct P and the free-surface ghost (source below the receiver)")
    zr = 2.0
    for zs in (20.0, 40.0):
        cfg = fk.FKConfig(fs=50, nt=2000, f_gauss=6.0, t_shift=1.0, L_factor=3.0)
        w_bh = fk.greens(m, zs, [0.01], [0.0], cfg, zr=zr)[0, 2, 2]
        t1, t2 = 1.0 + (zs - zr) / 6.0, 1.0 + (zs + zr) / 6.0
        i1 = int((t1 - 0.2) * 50) + int(np.argmax(np.abs(w_bh[int((t1 - 0.2) * 50):int((t1 + 0.2) * 50)])))
        i2 = int((t2 - 0.2) * 50) + int(np.argmax(np.abs(w_bh[int((t2 - 0.2) * 50):int((t2 + 0.2) * 50)])))
        exp = (zs - zr) / (zs + zr)
        # pulse amplitudes as signed L2 norms over +-0.15 s (sampled peaks of a 0.04 s pulse carry a few % error)
        n1 = w_bh[int((t1 - 0.15) * 50):int((t1 + 0.15) * 50)]
        n2 = w_bh[int((t2 - 0.15) * 50):int((t2 + 0.15) * 50)]
        rat = np.sign(w_bh[i2] * w_bh[i1]) * np.linalg.norm(n2) / np.linalg.norm(n1)
        line(f"   source {zs:3.0f} km: direct peak {i1 / 50:.3f} s (ray {t1:.3f}), ghost {i2 / 50:.3f} s (ray {t2:.3f}); "
             f"ghost / direct {rat:+.4f} (far-field image: +{exp:.4f}); relative difference {rat / exp - 1:+.4f}")
    ok &= check("|ghost delay - 2 zr / vp| (s)", abs((i2 - i1) / 50 - 2 * zr / 6.0), 0.021)
    ok &= check("|ghost / direct / far-field ratio - 1| at 40 km (tends to 0 with depth)", abs(rat / exp - 1), 0.005)
    return ok


def s4_layered():
    line("4. Layered medium crust_a: P arrival of the engine against ray theory (own tracer) and the detector's tracer.")
    line("   Engine arrival: time of the peak of the displacement-vector norm of the P pulse (the far-field displacement")
    line("   of a Gaussian moment rate is a Gaussian centred on the arrival); elastic copy of the model (Q = 1e9) and")
    line("   the model with its Q (constant-Q dispersion makes the 10 Hz band arrive slightly earlier than the 1 Hz")
    line("   reference of the ray times)")
    from ttmt import traveltime as TT
    m = media.PRESETS["crust_a"]
    me = media.LayeredModel(m.thickness, m.vp, m.vs, m.rho, name="crust_a_elastic")
    X = np.array([5.0, 15.0, 30.0, 45.0, 60.0])
    zs = 8.0
    tp, _p, kind = media.first_arrival(m, zs, X, "P")
    ref = TT.first_arrival([(float(a), float(v)) for a, v in zip(m.tops, m.vp)], zs, X, 1.0)
    ok = check("own tracer against ttmt.traveltime, max |dt| (s)", float(np.max(np.abs(tp - ref))), 1e-6)
    cfg = fk.FKConfig(fs=50, nt=2000, fmax=20.0, f_gauss=10.0, t_shift=1.0, L_factor=1.5)
    worst = 0.0
    for mm in (me, m):
        x = fk.greens(mm, zs, X, np.zeros(len(X)), cfg, quantity="displacement")
        out = []
        for i in range(len(X)):
            u = x[i, 2] + x[i, 0]
            env = np.sqrt(np.sum(u * u, 0))
            k0, k1 = int((1.0 + tp[i] - 0.1) * 50), int((1.0 + tp[i] + 0.1) * 50) + 1
            k = k0 + int(np.argmax(env[k0:k1]))
            # parabolic refinement of the peak
            y0, y1, y2 = env[k - 1], env[k], env[k + 1]
            d = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
            out.append((k + d) / 50 - 1.0 - tp[i])
        line(f"   {mm.name:<16} engine - ray (s) at X {X.tolist()} km: {np.round(out, 3).tolist()} (rays: {kind})")
        if mm is me:
            worst = float(np.max(np.abs(out[:4])))
    line("   at 60 km the direct P grazes the 4-12 km layer and is weaker than the head wave 0.04 s behind it and the")
    line("   wide-angle reflections, so the peak of the +-0.1 s window is not the first arrival (not an engine error)")
    ok &= check("elastic model, X <= 45 km: max |engine P - ray P| (s)", worst, 0.01)
    return ok


def _detector_view(wave, fs, arr, pre):
    """noise-free arrival vectors as the detector sees them: 0.5-6 Hz, 16 Hz, components (E, N, down), matched filter
    (2 Hz Ricker), sampled at the nearest sample of each true arrival. (2, ns, 3)."""
    from scipy.signal import butter, resample_poly, sosfiltfilt
    from ttmt.kernel import LEGACY
    x = wave.copy()
    x[:, 2] *= -1.0
    sos = butter(4, [0.5, 6.0], btype="band", fs=fs, output="sos")
    x = sosfiltfilt(sos, x, axis=-1)
    x = resample_poly(x, 4, 25, axis=-1) if fs == 100.0 else resample_poly(x, 16, int(fs), axis=-1)
    MF, _sd = LEGACY.matched(x)
    v = np.zeros((2, wave.shape[0], 3))
    for b in range(2):
        for i in range(wave.shape[0]):
            k = int(round((pre + arr[i, b]) * 16.0))
            v[b, i] = MF[i, :, min(max(k, 0), MF.shape[-1] - 1)]
    return v


def s5_subspace(libs):
    line("5. Inverse-crime distance: fraction of the noise-free arrival energy inside the detector's 6-D subspace at "
         "the true source (C19, v1 generator: 0.81 mean, 0.70 10th percentile)")
    from ttmt import radiation as RAD
    ok = True
    for path in libs:
        d = pickle.load(open(path, "rb"))
        STA = np.asarray(d["xyz"], float)
        fr = {}
        for e in d["events"]:
            v = _detector_view(e.wave, d["cfg"]["fs"], e.arr, e.pre)
            Q = RAD.basis(RAD.radiation_operator(np.asarray(e.xyz, float), STA))
            xv = v.reshape(-1)
            fr.setdefault(e.family, []).append(float(np.sum((Q.T @ xv) ** 2) / (xv @ xv)))
        line(f"   {path.split('/')[-1]} ({d['cfg']['medium']}, zr {d['cfg']['zr']} km, {len(d['events'])} events):")
        for f, a in sorted(fr.items()):
            a = np.array(a)
            line(f"      {f:<9} n {len(a):3d}  mean {a.mean():.3f}  p10 {np.percentile(a, 10):.3f}  p90 "
                 f"{np.percentile(a, 90):.3f}")
        allv = np.concatenate([np.array(a) for a in fr.values()])
        ok &= check("mean fraction below the v1 generator's 0.81 (further from the detector)", float(allv.mean()), 0.81)
    return ok


def s6_sets():
    line("6. Semi-synthetic sets (data/manifests/synth)")
    from ttmt.atlas import cells
    for mf in sorted(glob.glob(str(sets.MANIFEST_DIR / "*.json"))):
        man = json.loads(open(mf).read())
        sid = man["set_id"]
        line(f"   {sid}: {len(man['realisations'])} realisations ({len(man['partitions']['tune'])} tune, "
             f"{len(man['partitions']['read'])} read), source {man['source_meta'].get('kind')} "
             f"(n {man['source_meta'].get('n')}), {len(man['stations'])} stations, anchor {man['anchor']!r}")
        snr, sep, cc, km, kme, nt, npair, sc = [], [], [], [], [], 0, 0, []
        for part in ("tune", "read"):
            rd = sets.SetReader(sid, part)
            for lab in rd:
                t = rd.truth(lab)
                tg = t.target
                nt += int(tg.sum())
                npair += len(set(t.pair[tg]) - {-1})
                snr.extend(np.nanmedian(t.snr[tg], 1).tolist())
                sep.extend(cells.separation(t.arr_t_obs[tg]).tolist())
                cc.extend([x["max_cc"] for x in t.meta["targets"]])
                sc.extend([x["scale"] for x in t.meta["targets"]])
                km.extend(t.kmec[tg][np.isfinite(t.kmec[tg])].tolist())
                kme.extend([x for x in t.meta["kmec_empirical"] if x is not None and np.isfinite(x)])
        snr, sep, cc = np.array(snr, float), np.array(sep, float), np.array(cc, float)
        line(f"      targets {nt}, pairs {npair}; network SNR quantiles 10/50/90 "
             f"{np.round(np.nanpercentile(snr, [10, 50, 90]), 1).tolist()}; SNR bins <4/4-8/8-20/>=20: "
             f"{np.bincount(cells.bin_of(snr, cells.SNR_EDGES)[np.isfinite(snr)], minlength=4).tolist()}")
        sb = cells.bin_of(sep, cells.SEP_EDGES)
        line(f"      separation bins {cells.SEP_LABELS}: {np.bincount(sb[sb >= 0], minlength=len(cells.SEP_LABELS)).tolist()}"
             f"; regime T [0.1, 1) s: {int(np.sum((sep >= 0.1) & (sep < 1.0)))}")
        line(f"      max CC to the atom library: median {np.nanmedian(cc):.3f}, p90 {np.nanpercentile(cc, 90):.3f}, "
             f"fraction < 0.5 {np.mean(cc < 0.5):.3f} (CR-06); injection scale median {np.median(sc):.3g}")
        if km:
            line(f"      K_mec/E (known tensors) quantiles 10/50/90 {np.round(np.percentile(km, [10, 50, 90]), 4).tolist()}")
        if kme:
            line(f"      K_mec/E empirical (measured vectors) 10/50/90 {np.round(np.percentile(kme, [10, 50, 90]), 4).tolist()}")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fk-lib", nargs="*", default=[])
    ap.add_argument("--skip-engine", action="store_true")
    a = ap.parse_args()
    ok = True
    if not a.skip_engine:
        ok &= s1_whole_space()
        ok &= s2_free_surface_and_ghost()
        ok &= s4_layered()
    if a.fk_lib:
        ok &= s5_subspace(a.fk_lib)
    s6_sets()
    line("validation " + ("passed" if ok else "WITH FAILURES"))


if __name__ == "__main__":
    main()
