# Roadmap

Target: a doctoral-level thesis (ADR-0001). Every phase ends at a gate written here before the phase starts; a gate that
is not met is reported and the plan is revised in a decision record, not silently. Baseline of 18 October 2026: the
prototype detector (MT-FW and BOMP-PK with bounded influence), eleven registered tests, two negative real-data tests,
audit score 53.3/100 at doctoral level.

## Phase 0. Governance (done 18 October 2026)

Research protocol and templates (`docs/protocol/`), traceability (`docs/TRACEABILITY.md`), decision records
(`docs/decisions/`), data access (`docs/data_access.md`), this roadmap, changelog.
Decisions of the author (18 October 2026): code under MIT and documents under CC BY 4.0; data access by allowing the
domains of `docs/data_access.md` (ADR-0005, option a); learned components trained on the author's local GPU; the
thesis is presented at Universidad Icesi (programme and director to be added to the front matter).

## Phase 1. Repository and core package

Deliverables: `src/ttmt/` package with configuration objects; `archive/legacy/` with the current code byte-identical
and hash-indexed; `experiments/registered/` and `experiments/exploratory/` with one folder and card per task;
environment lock files per scenario; unit tests of the core functions (sample index convention, subspace projection,
radial maps, deflation, hit rule, exact FROC, partition boundary); regression tests against frozen outputs; continuous
integration with a smoke test; analyses of archived results runnable without raw data.
Gate: the regression suite reproduces the registered numbers of 2F, 2M, R0 and R1 from the package.

## Phase 2. Integrity of the existing results (declared post hoc)

1. Duplicate detections counted separately in 2F, 2I, 2M and R0; the effect of counting them as false on G1-R.
2. Tasks 2C and 2D re-read at the nearest sample (E52).
3. Whether the synthetic generator shares the wavelet and the radiation operator with the detector.
4. The within-segment percentile score of the R1 hybrid and the 19-block bootstrap, declared as defects.
Gate: each item reported in the log and in the evidence volume, with its effect on the conclusions.

## Phase 3. Failure atlas of the published methods

Deliverables: one benchmark harness with common hit rules, exact FROC and stratification by inter-event time,
magnitude, signal-to-noise ratio, station count, distance and overlap; runs of PhaseNet, PhaseNO, EQTransformer,
GaMMA, PyOcto, REAL, QuakeMigrate, EQcorrscan, Match and Locate, MALMI, BOMP-PK and MT-FW at their best; internal
diagnostics (merged and missed picks against arrival separation, split and merged events, sensitivity to the velocity
model); a second synthetic generator independent of the detector's forward model; development-only real data.
Gate: the harness reproduces published numbers of at least two comparators on their own example data, and the atlas
identifies regimes where every published method loses a measurable fraction of events.

## Phase 4. The new model, component by component

Order and gates (each measured on development data with an ablation that isolates it; see `docs/TRACEABILITY.md`,
Section 3):
1. Learned or eikonal 3-D travel times with source-dependent terms.
2. Evidence from picker probability fields; multi-station coupling constrained by travel times and moment-tensor
   consistency.
3. Adaptive bounded-influence statistic.
4. Joint multi-event refinement; deflation by empirical Green's functions.
5. False-alarm calibration with a deployable threshold.
6. Full-wave modelling, only if its gate in ADR-0003 is met.
Gate of the phase: on development data the model beats the best comparator of the atlas in the majority of its
regimes at the registered false-detection budgets, with a deployable threshold.

## Phase 5. Theory (in parallel with phases 3 and 4)

Identifiability of overlapping events; recovery guarantees for block-greedy pursuit with a perturbed dictionary and a
timing tolerance; an adaptive bounded-influence statistic that does not saturate; false-alarm control of the full
detector (decoys or conformal). Gate: each result stated with complete proofs and checked numerically.

## Phase 6. Confirmatory tests

At least two registered tests on unseen real sequences from networks not used in development, anchored on GitHub
before the run (ADR-0004), with comparators at their best, power computed before the run, deployable thresholds, and a
verification that includes a number re-derived by a person.

## Phase 7. Thesis

The thesis is rewritten around the atlas, the new model, the theory and the confirmatory tests; the present
benchmark, G1-R, R0 and R1 become the chapter on the prototype. The thesis is concluded only after phase 6.

## Risks

Data access (ADR-0005); compute for learned components (two cores in the working environment); a new model that does
not win in most regimes, in which case the thesis reports where and why each method fails, which remains a
contribution; scope, controlled by the gates.
