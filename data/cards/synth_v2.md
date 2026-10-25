# Data card: semi-synthetic sets of the generator v2 (ATLAS-0)

- What: real continuous noise of a registered network with injected events, and the truth table of every injected
  event (`ttmt.atlas.schema.Truth`). Code: `src/ttmt/synth/` (generator v2). Build: `tools/build_semisynth.py`;
  validation: `tools/validate_synth.py`, report `results/synth_v2/README.md`.
- Noise: Ridgecrest 2019 (one hour, 11 HH stations, 100 Hz, response-scaled velocity) and Hi-net 2012 (2000 s, 7
  borehole stations, 100 Hz, counts), read through `ttmt.datasets` (data cards `ridgecrest_2019.md`,
  `hinet_swarm_2012.md`; both development only). Alto Tiberina 2010 is not used: its registered copy is decimated to
  25 Hz (Nyquist 12.5 Hz), below the 2-20 Hz band of the ATLAS-0 card and the 100 Hz of the pickers.
- Time blocks: each record is cut into 300 s blocks; a seeded hash gives every block one role, EGF or noise
  (`realdata.block_roles`). Recorded events used as injections (EGFs) come only from EGF blocks; injections go only
  into noise blocks. The catalogued events of a noise block stay in the noise as neutral events (target False).
- Partitions: the noise blocks are split 30/70 into `tune` and `read` by a hash of the block start, with at least one
  block in each (`realdata.partition_blocks`); every realisation inherits the partition of its block, so tuning and
  reading never share a noise sample or an injection seed. `ttmt.synth.sets.SetReader` serves one partition and
  raises `PartitionError` on any other.
- Sources:
  - `egf`: recorded events of the EGF blocks with no comparable event within [-8, +30] s (Hi-net: [-4, +15] s, since
    the swarm is dense), cut [-5, +35] s (Hi-net [-4, +20] s). Arrival times are AIC onsets measured on the EGF
    itself around a 1-D prediction (Hadley-Kanamori; JMA2001); an onset more than 0.5 s (P) or 0.8 s (S) from the
    prediction, or with pick SNR below 4, is dropped (nan). Their precision is that of an AIC pick on a clean
    record, a few hundredths of a second, not exact. Tensors unknown: `kmec` is nan; the measured-vector value is in
    `meta.kmec_empirical`.
  - `fk`: full-wave synthetic events with known tensors (`ttmt.synth.fklib`): layered medium with a free surface and
    constant Q that is none of the detector's models, Brune or extended source-time functions (no Ricker), DC, DC +
    CLVD, DC + ISO, generic, double-source and extended families, positions out to 1.5 times the aperture, half of
    them clustered near another event; receiver effects: borehole depth (Hi-net, 0.1 km assumed), shear-wave
    splitting at 30 per cent of the stations, misorientation of 10-45 degrees at 15 per cent. Exact arrival times
    (ray theory in the true medium, checked against the engine to 0.01 s).
- Realisation: one noise block with groups every 40 s: an isolated event (30 per cent) or a pair whose origin
  separation is drawn from the bins < 0.1, 0.1-0.5, 0.5-1, 1-3 s, whose partner is taken from a spatial-separation bin
  (0-1, 1-5, 5-15 km) when one is available, or is the same source again (20 per cent: co-located, same mechanism),
  with energy ratio 1, 10 or 100. The first event of a group gets a network SNR drawn log-uniform in [3, 60]; EGFs
  are never amplified (scale <= 1, so their SNR is at most their recorded one).
- Truth per target: arr_t (true arrivals) and arr_t_obs (with the optional injected timing error, zero by default),
  per-station SNR as the ATLAS-0 card defines it (exact: signal and noise are stored apart), pair id, kmec when the
  tensors are known, and in `meta`: seed, configuration and its hash, source id and kind, injection scale, spatial
  separation of the pair, the radiated vectors g_ib measured on the noise-free injected waveform, the largest
  normalised network correlation with the default atom library (CR-06), and the 4-station contamination variants
  (CR-24).
- CR-06 atom library: by default every catalogued event of the noise blocks, cut [-5, +35] s
  (`realdata.template_candidates`), that is, the events a template or deflation method would take from the record
  the injections go into. A method whose dictionary differs recomputes the covariate with `inject.max_cc`.
- Seeds: every realisation is one record of the anchored streams of ATLAS-0 (`ttmt.atlas.seeds`, master
  8e1a8d5630dc0d52): injections, timing errors and amplitudes from `rng("inject", record)`, contamination stations from
  `rng("subset", record)`. Record indices: Ridgecrest EGF 0-69, Hi-net EGF 1000-1039, known-tensor sets 2000+ and
  3000+. The full-wave event libraries draw from `rng("inject", 1_000_000 + j)`, a range no record uses.
- Partitions deviate from the card's "hash of the hour": the records are one hour (Ridgecrest) and 2000 s (Hi-net), so
  the hash is of the 300 s noise block.
- Files: under `$TTMT_WORK_ROOT/synth/<set id>/` (not in git); manifests with SHA-256 and float fingerprints in
  `data/manifests/synth/`. A rebuild on another platform is checked with `SetReader(..., verify="fingerprint")`.
