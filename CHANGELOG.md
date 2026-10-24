# Changelog

Dated record of changes to the project's process, code and documents. Results are recorded in the sealed log
(`preregistration/`), errors in `docs/artifact_manifest_and_error_log.md`.

## 2026-10-24

- Atlas (phase 3; K1, K6, K9, CR-02, CR-03, CR-05, CR-14, CR-24): ATLAS-0 card final draft, ready to be anchored. It
  copies ADR-0008 (regimes T and C, LOC+T-V with the duplicator control, H1 then H2, G4 primary on semi-synthetic
  data), keeps the 240 cells as descriptive strata, restates H-ATLAS on the pooled regime T, fixes the 60-trial tuning
  budget with its ledger, names `prototype-v1` as the internal baseline, and adds the fine-tuned-picker comparator, the
  domain gate per generator, the K_mec/E scope and the 4-station contamination stratum. Harness `ttmt.atlas` (schema,
  cells, scoring on `ttmt.stats`, pick and event diagnostics) with `tests/test_atlas.py`, run in CI.
- Statistics (review of PR #11): H2 of ADR-0008 runs the studentised sign-flip test at a level calibrated by
  simulation at the least favourable configuration of the intersection-union test (one comparator at the margin),
  with at least 50 blocks (`sample_size.py calibrate`; 0.015 at 50 blocks); the nominal 0.025 gave size 0.026 to
  0.033. `tau` in the sample-size model no longer moves the marginal recall. R2 tests the target regime only for pairs
  at least nine times the truth's origin-time error apart (SAP 3.4).

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
