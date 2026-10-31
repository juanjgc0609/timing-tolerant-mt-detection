# Generator v2: validation and first semi-synthetic batch

Validator output: `validation_out.txt` (`python tools/validate_synth.py`, run 2026-10-30 on this branch). Sets built
2026-10-25 with the code of b66bbbd (same generator and configuration as this head; the later commits change
docstrings, the reader's pre-anchor guard and the data card only, so `inject_config_hash` a2af3a6066f1f2a5 is
unchanged). Data card: `data/cards/synth_v2.md`.

## Engine (sections 1-4): passed

| Check | Result | Tolerance |
|---|---|---|
| Whole space (Aki-Richards 4.29), 6 components x 3 distances, default L factor 3 | worst relative error 6.1e-4 | < 0.01 |
| Same, L factor 1.5 / 6 | 6.4e-3 / 6.1e-4 (converged at 3) | |
| Free surface, vertical-incidence P, source 30 km | surface / whole-space = 1.9998 | \|r - 2\| < 0.005 |
| Borehole receiver 2 km deep (CR-09): ghost delay | 6.7e-3 s off 2 zr / vp | < 0.021 s |
| Borehole: ghost / direct against far-field image, source 40 km | 3.9e-4 | < 0.005 |
| Own ray tracer against `ttmt.traveltime` | 8e-13 s | < 1e-6 s |
| Layered crust_a, elastic: engine P peak against ray P, X <= 45 km | 2.4e-3 s | < 0.01 s |

With Q the P peak arrives 3-11 ms early, as constant-Q dispersion predicts for the 10 Hz band against the 1 Hz
reference. At 60 km the direct P grazes a layer and the window peak is the head wave behind it, so the check stops at
45 km (not an engine error).

## First batch (section 6): EGF injections into real noise

All realisations seeded from the ATLAS-0 anchor (ADR-0009 item 5, `realisation_rng(set_id, i, r)`, master
8e1a8d5630dc0d52); partition by `realdata.tune_read_split` (ADR-0009 item 4).

| Set | Realisations (tune / read) | Sources | Targets | Pairs | In regime T [0.1, 1) s | Network SNR 10/50/90 | Max CC to atoms: median, fraction < 0.5 |
|---|---|---|---|---|---|---|---|
| `ss_ridgecrest_2019_egf_v1` | 70 (10 / 60) | 12 EGFs, 50 atoms | 830 | 340 | 379 | 0.7 / 4.5 / 16.6 | 0.116, 0.831 |
| `ss_hinet_2012_egf_v1` | 40 (10 / 30) | 8 EGFs, 79 atoms | 432 | 197 | 146 | 0.0 / 1.6 / 9.3 | 0.645, 0.000 |

Tuning blocks: Ridgecrest 2100 s; Hi-net 0 s (from the fallback, no Hi-net noise block hashes below 0.3).

## Findings to carry

1. **CR-06, Hi-net.** Every Hi-net injection correlates above 0.5 with some catalogued event of the same swarm
   (median 0.645), so Hi-net has no CC < 0.5 stratum: there a template or deflation method sees near-copies of its
   targets. The CC < 0.5 stratum comes from Ridgecrest (83 per cent of its targets) and, once built, from the
   full-wave sets. The atlas should report Hi-net as the high-CC stratum, not pool it.
2. **SNR.** EGFs are never amplified (scale <= 1), so the network SNR is at most the recorded one: low for Hi-net
   (median 1.6, 309 of 432 targets below 4). The SNR 6-120 range of spec §13 is reached only with the full-wave sets.
3. **K_mec/E.** EGF tensors are unknown, so `kmec` is nan; `meta.kmec_empirical` (measured radiated vectors) has
   median 0.03 in both sets. Known-tensor K_mec/E (CR-14) needs the full-wave sets.
4. **CR-05 domain gate.** The generator stores what the gate needs (per-target SNR, separation, CC, K_mec); the gate
   itself is run by the atlas thread.
5. **Alto Tiberina** is not used: its registered copy is decimated to 25 Hz, below the 2-20 Hz band of the card.

## Not yet in this batch

The full-wave sets with known tensors (`ss_*_fk_v1`, libraries `library_rng(0)` Ridgecrest and `library_rng(1)`
Hi-net). Their libraries were lost when the cloud container restarted and are being rebuilt; they will come with
their manifests and section 5 of the validator (subspace fraction).

## Rebuilding the files

The waveform files (318 MB) are not in git. On a machine with the registered data:

    TTMT_DATA_ROOT=... TTMT_WORK_ROOT=... python tools/build_semisynth.py --dataset ridgecrest_2019 --source egf --per-block 10
    TTMT_DATA_ROOT=... TTMT_WORK_ROOT=... python tools/build_semisynth.py --dataset hinet_2012 --source egf --per-block 10

then check against the manifests in git with `SetReader(set_id, partition, verify="fingerprint")` (the bytes can
differ across platforms, the float fingerprints to 1e-6).
