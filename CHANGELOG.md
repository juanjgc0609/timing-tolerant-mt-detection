# Changelog

Dated record of changes to the project's process, code and documents. Results are recorded in the sealed log
(`preregistration/`), errors in `docs/artifact_manifest_and_error_log.md`.

## 2026-11-08

- M5 decision step (exploratory, ADR-0010 pivot): `ttmt.v2.m5` with the mixture e-value S_mix (theory D.2) and
  Student t, gaussian and radial (kernel bulk + Pareto tail) backgrounds; tests `tests/test_m5.py` (in CI).
  Exploration M5-EVALUE, gate G-M5c, on Task 2M real noise (owner's laptop, CPU): FAILS on both seeds at p = 1e-3
  (radial background 5.8e-3 on 101231, 1.2e-3 on 100210-11; passes at 1e-2); the gaussian background fails by two
  orders of magnitude on the transient day, as D.1(a) predicts. The declared second step found the cause in
  dependence between blocks (simultaneous tails across stations, 1.7x the independent rate), not in the marginal,
  and pivoted to M7 by its rule; the held-out confirmation day stays unused. A third attempt, the average of
  per-station e-values (valid under any dependence across stations), already misses the level on a fake i.i.d. null
  (1.4e-3 at 1e-3), so it was not run on real data. Negative result recorded; the decision rests on M7.
- Thesis: Section 8 reports the amendment-4 GaMMA tuning (60 of 60 valid, best 0.2356 against PyOcto 0.2589); with one
  search per method and no intervals neither is ordered above the other.
- M6, first version (plan T4.19, exploratory): `ttmt.v2.m6.adcg_loop`, joint refit of all accepted atoms (closed-form
  tensors, bounded block times, re-assignment among onset maxima) with the one-against-two resolvability test in
  place of the 5 km exclusion; tests `tests/test_m6.py` (in CI). Exploration M6-JOINT (simulation, SNR 12): pair
  recall of near pairs (1.5-5 km) at 0.1-1 s rises from 0.35-0.85 (greedy, 5 km exclusion) to 0.75-1.00; the greedy
  loop without the exclusion recovers the same pairs (0.70-1.00) but splits 40 of 40 single events, M6 0 of 40: the
  gain is the exclusion removed, made admissible by the resolvability test. Far pairs unchanged. Real G-M6a
  (semi-synthetic, real noise) not measured yet.
- Model v2, first components (exploratory, plan T4.1, T4.13, T4.15): subpackage `ttmt.v2` with M0 (`m0`: Tyler noise
  covariance per station, whitening of traces and of radiation operators, orientation audit), M1 (`m1`: onset
  evidence, isotonic calibration, `OnsetEngine` on OnsetDelay terms, `OnsetReader`, sub-sample vector reading,
  onset-to-peak lag), M4 (`m4.joint_read`) and the v2 loop (`loop.v2_loop`: onset or joint reading, ACAP ranking,
  deflation at the waveform-refined time with the tensor refitted). `ttmt.detector`: `Engine.terms_type` and
  `tpred`, `greedy_loop(return_residual=)`, `subsample_delay(half=)`, deflation convention `"refit"`; defaults and
  registered numbers unchanged. Tests `tests/test_v2.py` (in CI).
- Exploration M6C-SUBSAMPLE (gate G-M6c): passes for isolated events (median peak residual 0.035 v2, 0.043 E62
  fractional, 0.197 grid), fails with two overlapping sources (0.06-0.12) for every deflation that keeps the other
  event; the fix is the joint refit of M6.
- Exploration M1A-ONSET (gate G-M1a, owner's laptop, CPU): the first step FAILED on CAL (onset reading, P 0.87
  within 0.25 s < 0.90; energy peak 0.55). A second step declared in the card before running (diagnosis, one fixed
  "prior" variant, confirmation on the held-out EVAL events) PASSED on EVAL: P 0.92 for the reader as registered,
  0.97 for the prior variant, energy peak 0.52 on the same gated blocks. The 0.92 is over the 61 of 103 gated P
  blocks that have an onset maximum; counting the 42 unread blocks as misses it is 0.54 against 0.52: the gain is
  in precision, not in coverage (coverage is the C-regime picker loss, step G-M1b). The prior reading is now the
  default reader of the v2 loops (`ttmt.v2.loop.default_reader`, `OnsetReader(sigma=)`).
- GaMMA re-tuned under ATLAS-0 amendment 4 (DBSCAN in the search space): 60 of 60 trials valid, best objective
  0.2356 (trial 53), frozen in `frozen/phasenet+gamma_a4.json`; it replaces the registered and amendment-3 GaMMA
  configurations. The stopped amendment-3 ledger (3 trials) is kept and not used.
- Thesis: Section 8 status paragraph updated for gate G3a passed (2 of 2: PyOcto on Puente 9/9, QuakeMigrate on
  Rutford exact; GaMMA 0/8 and ETHZ not passed; four declared deviations of the QuakeMigrate run) and amendment 4
  (449d39b); GaMMA stays provisional until the a4 run ends; no claim over GaMMA.

## 2026-11-04

- ATLAS-0 amendment 4 (draft): GaMMA's DBSCAN (`dbscan_eps` [3, 20] s, `dbscan_min_samples` {3..20}) joins its search space, four associator parameters as PyOcto; 60 trials into `tuning_ledger_a4.csv` (`ledger.require_a4`); the amendment-3 run was stopped on the owner's decision and is not used.
- GaMMA's gate failure (target 1, 0/8) diagnosed: the DBSCAN pre-clustering of Table S3 (min_samples 20) against the station density, not E63 and not the travel-time grid depth (`experiments/exploratory/COMPARATORS/README.md`); `gamma_assoc` gains `eikonal_zmax` (default None, registered behaviour), `repro_puente.py` gains `fetch-stations --all-networks` and `run --post-hoc`, which gate.py never counts.

## 2026-11-03

- Thesis: the [pending #33] marks and the draft note removed after PR #33 merged.
- Thesis (pending PR #33, decision D10): E52 nearest-sample re-run of 2C and 2D (no verdict changes; 2C-a and 2D-b
  strong 2 of 3 cells on seed 7331, still met); 2C and 2D rescored with repeats as false and 2L one to one (LOC+T and
  LOC), every verdict stands; 2E not regenerated; E64 (deployable R1 threshold, 41.8 false/h against 30) cited in
  Sections 6.6 and 7.4; error count 64.

## 2026-11-02

- Thesis: Section 8 gains a provisional status paragraph of the doctoral line (gate G3a: PyOcto 9 of 9 and GaMMA 0 of
  8 on Puente, ETHZ not reproduced, QuakeMigrate running; tuning: phasenet+pyocto 0.2589 frozen, phasenet+gamma
  provisional by E63; amendment 3 anchored at 83c66e4), with no claim over GaMMA; error count 63.

## 2026-11-01

- E63: GaMMA 80394dd cannot run method GMM (`gamma/_base.py:140` reads `covariance_prior`, which only the Bayesian mixture has), the cause of the 37 failed GMM tuning trials. Repair prepared as ATLAS-0 amendment 3 (draft, not anchored): method pinned to BGMM keeping the registered draws, the GaMMA wrapper refuses GMM, tests; 60 BGMM trials to run on the device into `tuning_ledger_a3.csv`.
- Atlas noise-block sensitivity: the plan carries the set's `dataset`, and the sensitivity bootstrap keys the noise
  block by dataset and draws together the set strata that share one (the EGF and known-tensor sets of a network use
  the same noise, PR #26). Noted in `experiments/exploratory/ATLAS-0/run_log.md`.
- CR-04 prevalence rule (`docs/protocol/CR04_PREVALENCE.md`): ATLAS-0, G4 and tuning keep weights 1, as the anchored
  card states; R2 and R3 take the measured stratum shares (separation, station subset, distance) of the densest
  reviewed catalogues, with a block-bootstrap interval, reported as a lower bound with a Poisson reference. The
  measurement runs on the author's machine; the R2 and R3 drafts point their weight tables to the rule.
- `ttmt.atlas.prevalence` and `experiments/exploratory/CR04-PREVALENCE/run.py` compute that table from a catalogue
  and a station list (windowed `cells.separation`, 1-hour block bootstrap, Poisson reference, dead time);
  `tests/test_prevalence.py` joins the CI test list.

## 2026-11-04

- Reproduction gate G3a passes, 2 of 2 required targets: PyOcto on Puente et al. (2025), 9/9 conditions within 5 %
  of Table S10, and QuakeMigrate 1.2.2 on the Rutford example of Winder et al. (2026), 49,249 / 47,354 / 20,015
  events, equal to Section 4.1. Not passed: GaMMA on Puente (0/8) and PhaseNet on ETHZ (better than published).
  Evidence in `experiments/exploratory/COMPARATORS/g3a/`; deviations of the QuakeMigrate run (newer dependencies,
  one waveform file re-fetched, `locate` resumed at event 48,454 with the same settings) in the comparators README.

## 2026-10-31

- Registered ATLAS-0 tuning of phasenet+pyocto (best objective 0.2589, 60 of 60 trials valid) and phasenet+gamma (0.2338, 23 of 60 valid: the 37 GMM trials fail in GaMMA 80394dd, a declared budget inequality, CR-07), with the frozen configurations; the owner kept these rows as registered (31 October); the comparison is provisional, and repairing GaMMA and repeating 60 BGMM trials in a new ledger is to be decided before any claim that a method beats GaMMA (`experiments/exploratory/COMPARATORS/README.md`).
- Reproduction gate (G3a, Puente et al. 2025): the 100-station floor of `gate.py` (added in PR #22, unsourced, not in
  the card) is replaced by the fixed station query of `repro_puente.py fetch-stations` (SCEDC, CI, HH?/BH?,
  2020-01-01, scenario box plus 50 km), checked from the provenance file and checksum of the list. The query, the
  5 % tolerance and the targets are unchanged; the SCEDC list has 54 stations. When the floor was changed, no
  Puente run on that list had finished.
- Reproduction gate, target 2 (ETHZ): the device run gives every P-pick metric better than published (MAE 0.086 s
  against 0.12, RMS 0.290 s against 0.42), outside the two-sided 5 % tolerance; not passed, tolerance unchanged. The
  ETHZ report cites the commit that verified the published values (e5be9cf).
- Reproduction gate, target 3: `repro_quakemigrate.py`, QuakeMigrate 1.2.2 on the Rutford icequake example of Winder
  et al. (2026), against the counts of Section 4.1 (49,249 triggered, 47,354 located, 20,015 filtered), fixed in
  `published_quakemigrate.json` before the run. Chosen by the owner as the second target after ETHZ failed.

## 2026-10-30

- Atlas run harness: `ttmt.atlas.plan` (the run plan of a set manifest: station subsets k = 4, 6, 8, 12 by
  `cells.station_subset`, contamination variants on the 4-station subset, CR-24) and `ttmt.atlas.run` (pools the
  reading runs by realisation, reports regimes, H-ATLAS loss and the duplicator control per scope: all sets, each set,
  correlation < 0.5 (CR-06) and each contamination variant). The run refuses the tuning partition and needs
  `experiments/exploratory/ATLAS-0/gate.json` with at least two passed reproduction targets.
- Atlas sensitivity to the bootstrap unit (review of PR #21): every scope also reports the number of distinct noise
  blocks and a sensitivity row that resamples whole noise blocks (`froc.describe` / `pooled_loss` with `cluster`),
  with H-ATLAS read against both; the registered unit stays the realisation. Noted in
  `experiments/exploratory/ATLAS-0/run_log.md`.

## 2026-10-25

- Generator v2 and semi-synthetic sets for ATLAS-0 (CR-05, CR-06, CR-09, CR-14, CR-24): `ttmt.synth` (full-wave
  layered engine with free surface, Q and borehole receivers; no operator shared with the detector, checked by
  `tests/test_synth.py`), injection of recorded events (EGF) and of full-wave events with known tensors into the real
  noise of Ridgecrest and Hi-net, with the `ttmt.atlas.schema.Truth` of every record, the card SNR, pair ids, K_mec/E,
  the largest correlation with the atom library and the 4-station contamination variants. Seeds `ttmt.atlas.seeds.realisation_rng`
  and `library_rng` (ADR-0009 item 5; helpers ported from #19). Tuning/reading partition by a hash of the 300 s noise block (the records are one hour
  or less, so a hash of the hour cannot split them). Data card `data/cards/synth_v2.md`; validation
  `tools/validate_synth.py`, report `results/synth_v2/README.md`. `SetReader` refuses a set not seeded from the anchor
  (`pre_anchor`) unless `allow_pre_anchor=True`.

## 2026-10-24

- ATLAS-0 amendment 1 (ADR-0009), before any tuning trial: REAL `tint` and PyOcto `min_interevent_time` log-uniform
  in [0.1, 5] s; tuning objective is the mean of the development recalls in T and C (isolated events in C);
  PhaseNet `instance` overlap with Alto Tiberina 2010 declared and reported separately (CR-22); tuning/reading split
  by 300 s noise block (`realdata.tune_read`, fallback both ways); one seed rule for realisations and the event
  library from the anchored master (`seeds.realisation_rng`, `seeds.library_rng`). Anchor
  `reg/ATLAS-0/amendment-1` is the owner's step.
- ATLAS-0 anchored (tag `reg/ATLAS-0/registration-1` on 7bafb7a, Release `created_at` 2026-10-24T22:37:39Z). Seeds
  derived after the anchor (ADR-0004 item 4) by `ttmt.atlas.seeds`: master `8e1a8d5630dc0d52`, named child streams
  (tuning, inject, subset, bootstrap, signflip) and the 60-trial tuning seed list, stored in
  `experiments/exploratory/ATLAS-0/seeds.json` (`python -m ttmt.atlas.seeds --check`). `froc.describe`, `froc.pooled_loss`
  and the Monte Carlo signs of `froc.gate` draw from the bootstrap and signflip streams by default.
- Atlas follow-up to PR #13: the ATLAS-0 card quotes the amended item 7 of ADR-0008 (H2 by the clustered score
  test at the level of the design check, at least 30 blocks; Holm threshold of non-inferiority at 0.03 scaled by that
  level); real-data pairs closer than 9 sigma_t are excluded from T (`ttmt.atlas.froc.truth_precision_targets`,
  `score_record(..., sigma_t=...)`).
- Atlas (phase 3; K1, K6, K9, CR-02, CR-03, CR-05, CR-14, CR-24): ATLAS-0 card final draft, ready to be anchored. It
  copies ADR-0008 (regimes T and C, LOC+T-V with the duplicator control, H1 then H2, G4 primary on semi-synthetic
  data), keeps the 240 cells as descriptive strata, restates H-ATLAS on the pooled regime T, fixes the 60-trial tuning
  budget with its ledger, names `prototype-v1` as the internal baseline, and adds the fine-tuned-picker comparator, the
  domain gate per generator, the K_mec/E scope and the 4-station contamination stratum. Harness `ttmt.atlas` (schema,
  cells, scoring on `ttmt.stats`, pick and event diagnostics) with `tests/test_atlas.py`, run in CI.
- Statistics (reviews of PR #11 and #13): H2 of ADR-0008 runs a clustered score test (`ttmt.stats.ni_score_pvalue`)
  instead of the margin sign-flip test: variance under the restricted null (cluster-robust, floored at Tango's 1998
  multinomial variance, so no discordance on few targets is not non-inferior) and a Student t reference with
  Satterthwaite degrees of freedom. Its level is fixed by a design check (`sample_size.py calibrate`) over a range of
  recall, correlations and targets per block taken from development data (`sample_size.py estimate`); the earlier
  forms had size up to 0.033 (sign-flip) and 0.06 to 0.11 (design-effect variance, large unequal blocks) at the least
  favourable configuration of the intersection-union test. H2 needs at least 30 blocks, as H1. `froc.gate` takes the
  H2 level. `tau` in the sample-size model no longer moves the marginal recall. R2 tests the target regime only for
  pairs at least nine times the truth's origin-time error apart (SAP 3.4).

- Statistics (K1, K2, K9, CR-02, CR-07, CR-08): ADR-0008 (accepted by the author on 24 October), one rule for "the model beats the comparators"
  (pooled regimes instead of counted cells, LOC+T-V hit rule with a duplicator control, block sign-flip tests,
  intersection-union over comparators, fixed sequence H1 -> H2, equal tuning budget, G4 on semi-synthetic data first).
  Analysis plan `docs/protocol/PHASE6_SAP.md`, conformance checklist, drafts of R2 and R3, `ttmt.stats` with its tests
  and the sample-size tool `experiments/exploratory/P6-SAP/`.

## 2026-10-23

- Thesis: Sections 5.2, 6.2, 7.3, 7.4 and Table 2 updated for the 2I, 2A, 2B, R0 and R1 rescoring (rule 2I-a not met
  under one-to-one LOC+T for seed 6262, +0.019 against 0.03; the other verdicts unchanged); the coverage gap of 2C, 2D,
  2E and 2L (records not stored, not rescored) stated; Table 1 gains row E6b.
- Executive summary: Task 2I gains the post hoc LOC+T reading (2I-a not met on seed 6262); thesis wording on 2A-e
  (differences move by at most 0.08) and on the cause of the 2I-a failure (the origin-time test) corrected per review.
- Scoring (E58, CR-21): the 2I primary arm regenerated with the frozen code and rescored under the registered and
  one-to-one rules, with and without origin time; 2A and 2B rescored with repeats as false (`rescore/rescore_2a2b.py`);
  R0 and R1 rescored with `--match one_to_one`. Regenerated 2F and 2I parts archived with checksums in
  `results/task2f/records/` and `results/task2i/records/`; regression tests (d) and (e) reproduce the registered 2F and
  2I analyses from them. Coverage table in `results/rescore_duplicates/README.md`.

## 2026-10-19 (thesis integration)

- Thesis: title page with the director; Section 4 corrected (AA.9.5(b) with Lipschitz constant 1: bound 14.59 instead of 28.50) and extended with a
  summary of Appendix AB; Sections 5.2, 6.2, 6.3, 7.1, 7.3, 7.4 rewritten for the one-to-one LOC+T rescoring of 2F
  and 2M (E58); C18j out of sample and the R1 analysis defects (E59 to E61) and the deflation on the sample grid (E62) added; exploration C19 added to the threats to validity; scope stated per ADR-0006; Table 1 gains rows
  E5b, E5c, E8b, E8c and a corrected E16.
- Appendix AB (`docs/proposal/Supplementary_Evidence_Appendix_AB.md`): identifiability of overlapping events, greedy
  recovery with timing error, the adaptive robust statistic, theory for the components of the new model; scripts and
  outputs in `experiments/exploratory/THEORY-5/`.
- AA.9.5 corrected with an erratum; Doctoral Proposal (8b) and (20) given their missing assumptions.
- Error register: the E57 entry, lost when E58 was merged, restored; entries in order E57, E58, E59 to E61, E62; count 62.

## 2026-10-18 (phase 1: repository and core package)

- Line endings: `.gitattributes` `* -text`, so that every file is checked out byte for byte. With `core.autocrlf=true`
  (the Windows default of many installations), git had rewritten 1014 text files to CRLF, and every freeze list
  (`code/r1/FREEZE_R1.sha256`, `data/checksums/*`) failed to verify. After a re-checkout the freeze lists of R0, R1
  and 2M match every file present.
- `code/` declared frozen and hash-indexed (`code/LEGACY_INDEX.sha256`; checked by `tests/test_legacy_index.py`). Its
  bytes are unchanged. It is append-only: `python tools/update_legacy_index.py` appends new files (the five of PR #4
  were added that way; 355 files). It refuses when an indexed file has changed or is missing. Unindexed new files make
  the test skip with their list, and CI prints them.
- Package `src/ttmt/` (`pyproject.toml`, `pip install -e .`): detector core (engine, energy-argmax reader, radial maps
  and block statistics, MT-FW, BOMP-PK, HUBER-N, deflation), radiation operator and subspace, homogeneous and layered
  travel times, matched filter, typed station terms, hit rules and exact FROC, partition guard, dataset readers of R0
  and R1, the R1 pipeline and the analyses of 2M, R0 and R1. Known defects are explicit, versioned conventions with
  the legacy value as default: sample rounding (E52), `Kernel` (E51), `PeakDelay` / `OnsetDelay` (E54),
  `SegmentReader` (E56), and hit rule `match="legacy" | "one_to_one"` (ADR-0007). `ttmt.runtime.pin_blas()` (E31).
  `ttmt.legacy` runs the frozen scripts unchanged, with their data paths redirected at import time.
- E62: the greedy deflation subtracts each atom on the integer sample grid (5.6 % of an arrival's energy left on
  average, 16.3 % at half a sample, for the 2 Hz filter). Convention `deflation="grid"` (legacy, default) or
  `"fractional"` (`Kernel.response(delta)`). Exploratory block statistic `ACAP` (cap relative to the median block
  norm of a hypothesis), not used by any registered run.
- Environments: `environments/core.lock.txt`, `core.yml` and `DEVIATIONS.md` (Python 3.12; the frozen pins exactly).
  Data fetch at the registered commits with checksum verification: `tools/fetch_data.sh`.
- Tests: unit tests, bit-for-bit equivalence with the legacy functions, and regressions against the frozen outputs.
  CI on ubuntu-latest (`.github/workflows/ci.yml`).

## 2026-10-18 (review of the plan)

- Review of the phase 0 to 7 plan (`docs/reviews/2026-10-18_plan_review.md`).
- ADR-0004 revised: a git commit date is not an external timestamp; registrations are anchored by an annotated tag and
  a GitHub Release, OpenTimestamps and, for confirmatory tests on real data, OSF Registries; synthetic seeds derived
  from the anchor. Protocol Section 6 updated.
- ADR-0006: the undergraduate thesis and the doctoral line are separate products (amends ADR-0001). ADR-0007: primary
  hit rule LOC+T with one-to-one matching for future tests.
- Roadmap: phase 2 runs before phase 1; `code/` stays in place; environment lock decided first; phase 2 extended;
  atlas regimes and counting rule fixed beforehand (`experiments/exploratory/ATLAS-0/card.md`); core comparator set
  with EQNet; one primary theory line; phase 6 sized by `experiments/exploratory/P6-power/`.
- ADR-0003 and the traceability table corrected on the evidence for PhaseNO, the closest prior work and the timing
  error of R0.
- Data cards of the three held datasets (`data/cards/`) and the `experiments/` layout.

## 2026-10-18

- Scoring (E58): one-to-one hit rule `code/froc_strict.py`; post hoc rescoring of 2F, 2I and 2M under the
  registered and corrected rules (`code/rescore/`, `results/rescore_duplicates/`).
- Governance (phase 0): research protocol and templates (`docs/protocol/`), traceability of the original plan
  (`docs/TRACEABILITY.md`), decision records ADR-0001 to ADR-0005 (`docs/decisions/`), data access
  (`docs/data_access.md`), roadmap (`docs/ROADMAP.md`).
- Target level fixed at doctoral (ADR-0001); baseline 53.3/100 on the doctoral rubric.
- Thesis Sections 6.2, 6.4, 7.2, 7.4 and 8: synthetic generator described as it is coded (gradient medium, curved
  rays, variable wavelet, coda) and its timing error stated in full; error E57.

## 2026-10-17 to 2026-10-18

- Task R1 registered, amended, run once and verified (negative); errors E51 to E56 logged; thesis document
  `docs/thesis/Thesis.md` with a generated evidence table and three figures; audits 6 and 7.
