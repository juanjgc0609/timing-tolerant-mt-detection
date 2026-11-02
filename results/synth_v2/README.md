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

## Known-tensor sets (sections 5 and 6, run 2026-11-01)

Full-wave libraries from `library_rng(0)` (Ridgecrest, crust_a perturbed, 120 events) and `library_rng(1)` (Hi-net,
crust_b perturbed, borehole 0.1 km, 80 events), injected into the same noise blocks and partitions as the EGF sets.

| Set | Realisations (tune / read) | Targets | Pairs | In regime T | Network SNR 10/50/90 | Max CC to atoms: median, fraction < 0.5 | K_mec/E 10/50/90 |
|---|---|---|---|---|---|---|---|
| `ss_ridgecrest_2019_fk_v1` | 70 (10 / 60) | 834 | 344 | 259 | 1.5 / 8.7 / 39.4 | 0.072, 1.000 | 0.000 / 0.013 / 0.199 |
| `ss_hinet_2012_fk_v1` | 40 (10 / 30) | 477 | 197 | 195 | 1.7 / 8.1 / 36.0 | 0.100, 1.000 | 0.001 / 0.024 / 0.184 |

- Inverse-crime distance (section 5): the mean fraction of the noise-free arrival energy inside the detector's 6-D
  subspace at the true source is 0.41 (Ridgecrest) and 0.40 (Hi-net), against 0.81 for the v1 generator (C19): the v2
  waveforms are far from the detector's own forward model, for every mechanism family.
- These sets give what the EGF sets lack: a CC < 0.5 stratum on both networks (CR-06), K_mec/E from known tensors
  (CR-14; computed on the noise-free radiated vectors, so `kmec` equals `meta.kmec_empirical`), and SNR up to ~40
  (half the targets above 8).
- **Magnitude (corrected after #26).** The first build of these sets wrote `mag` = library Mw + log10 scale, which
  is wrong in two ways: scaling a waveform by s scales the moment by s, so Mw moves by (2/3) log10 s, not log10 s
  (Ridgecrest came out about 0.7 units too low: scale median 0.0099); and on Hi-net the scale also carries the unknown
  counts-to-m/s factor (median 2.6e6). Now `mag` = library Mw + (2/3) log10 scale on Ridgecrest (response-scaled
  m/s), the moment magnitude of the injected source with the library event's source-time function; and nan on Hi-net
  (counts). Every target keeps `mw_library` and `scale` in `meta.targets`. On Hi-net use the SNR. EGF targets keep
  the catalogue magnitude of the recorded event (scale <= 1, median 0.99).

## Rebuilding the files

The waveform files are not in git. On a machine with the registered data:

    TTMT_DATA_ROOT=... TTMT_WORK_ROOT=... python tools/build_semisynth.py --dataset ridgecrest_2019 --source egf --per-block 10
    TTMT_DATA_ROOT=... TTMT_WORK_ROOT=... python tools/build_semisynth.py --dataset hinet_2012 --source egf --per-block 10
    TTMT_DATA_ROOT=... python tools/build_fk_library.py ridgecrest_2019 120 0 0.0 fk_ridgecrest_2019_0.pkl 4
    TTMT_DATA_ROOT=... python tools/build_fk_library.py hinet_2012 80 1 0.1 fk_hinet_2012_1.pkl 4
    TTMT_DATA_ROOT=... TTMT_WORK_ROOT=... python tools/build_semisynth.py --dataset ridgecrest_2019 --source fk --fk-lib fk_ridgecrest_2019_0.pkl --per-block 10
    TTMT_DATA_ROOT=... TTMT_WORK_ROOT=... python tools/build_semisynth.py --dataset hinet_2012 --source fk --fk-lib fk_hinet_2012_1.pkl --per-block 10

(the libraries take about 2 h on 4 cores and resume after an interruption), then check against the manifests in git with `SetReader(set_id, partition, verify="fingerprint")` (the bytes can
differ across platforms, the float fingerprints to 1e-6).
