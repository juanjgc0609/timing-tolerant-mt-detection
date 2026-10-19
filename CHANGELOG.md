# Changelog

Dated record of changes to the project's process, code and documents. Results are recorded in the sealed log
(`preregistration/`), errors in `docs/artifact_manifest_and_error_log.md`.

## 2026-10-19 (thesis integration)

- Thesis: title page with the director; draft-status note and **[pending #N]** markers for results of unmerged pull
  requests; Section 4 corrected (AA.9.5(b) with Lipschitz constant 1: bound 14.59 instead of 28.50) and extended with a
  summary of Appendix AB; Sections 5.2, 6.2, 6.3, 7.1, 7.3, 7.4 rewritten for the one-to-one LOC+T rescoring of 2F
  and 2M (E58); C18j out of sample and the R1 analysis defects added as pending (E59 to E61); deflation on the sample
  grid (E62, pending); exploration C19 added to the threats to validity; scope stated per ADR-0006; Table 1 gains rows
  E5b, E5c, E8b, E8c and a corrected E16.
- Appendix AB (`docs/proposal/Supplementary_Evidence_Appendix_AB.md`): identifiability of overlapping events, greedy
  recovery with timing error, the adaptive robust statistic, theory for the components of the new model; scripts and
  outputs in `experiments/exploratory/THEORY-5/`.
- AA.9.5 corrected with an erratum; Doctoral Proposal (8b) and (20) given their missing assumptions.
- Error register: the E57 entry, lost when E58 was merged, restored; count 58.

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
