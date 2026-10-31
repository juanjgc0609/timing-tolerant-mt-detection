# ATLAS-0 run log

Dated notes of the atlas run on the reading partition (`ttmt.atlas.run`). Results are not recorded here before the
reproduction gate passes (`gate.json`).

## 2026-10-30: bootstrap unit and noise blocks

- The registered bootstrap unit stays the injected realisation (card, "Regimes and strata"). In the first
  development batch (PR #20) each 300 s noise block carries 10 realisations, and the reading partition has 6 distinct
  noise blocks for Ridgecrest and 3 for Hi-net. Realisations on one noise block share every noise sample, so their
  noise-driven false detections are not independent and the registered intervals can look more precise than they are.
- Every scope therefore also reports `noise_blocks` and a sensitivity row, `sensitivity_noise_block`, whose bootstrap
  draws whole noise blocks within the set strata (streams ("bootstrap", 2) and ("bootstrap", 3)). H-ATLAS is read
  against both rows (`h_atlas_reading`); a reading that holds only under the registered unit is reported as such.
- With 3 noise blocks, Hi-net is effectively a single stratum: its own scope is descriptive only.
