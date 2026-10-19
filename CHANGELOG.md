# Changelog

Dated record of changes to the project's process, code and documents. Results are recorded in the sealed log
(`preregistration/`), errors in `docs/artifact_manifest_and_error_log.md`.

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
