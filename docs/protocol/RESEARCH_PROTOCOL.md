# Research protocol

This protocol governs every experiment of the project from 18 October 2026 on. It formalises the practice of the
registered tasks 2E to R1 and closes the gaps that the audits and the error register found in them. The target level
of the work is a doctoral thesis, and every audit is scored on the doctoral rubric (Section 9).

## 1. Stages

Every piece of work belongs to exactly one stage, and moves forward only through the gate of that stage.

| Stage | Purpose | Data allowed | Output | Gate to leave the stage |
|:--|:--|:--|:--|:--|
| S0 Idea | State a question and the mechanism it tests | none | Experiment card (template `templates/experiment_card.md`) | The card names the hypothesis, the comparator, the endpoint and what result would refute it |
| S1 Exploration | Build and debug; measure on development data | development partitions and calibration seeds only | Exploration report in `experiments/exploratory/<id>/` with the commands and outputs | A gate written in the card before the exploration began is met, or the line is closed and reported |
| S2 Registration | Fix the confirmatory test before any evaluation sample is read | none new | Registration (template `templates/registration.md`) appended to the sealed log, pushed to the public repository | The registration commit is on GitHub (external timestamp) and its hash was sent to the author |
| S3 Freeze | Fix the code, data and environment | development only | Freeze list (SHA-256 of code, data, environment lock, model weights) and, if needed, an amendment (template `templates/amendment.md`) appended and pushed before the run | Freeze list committed and pushed before the run |
| S4 Run | Execute once | evaluation partition | Raw outputs archived in `results/` | Single run completed; any interruption logged with its cause |
| S5 Verification | Re-derive the numbers independently | the archived outputs and the data | Verification report (template `templates/verification.md`) | Every registered number reproduced or every discrepancy explained; at least one number re-derived by a person for confirmatory tests on real data |
| S6 Reporting | Record and publish | none | Result entry appended to the log (template `templates/result.md`), documents updated from the result files | Every number in a document traced to a file and line |

## 2. Data partitions

1. Every dataset is split before any use into development and evaluation partitions, written in its data card
   (`data/cards/<dataset>.md`). A dataset used for evaluation once is development data from then on.
2. No read of development may touch an evaluation sample, including the context margins of analysis windows. The
   reader enforces the boundary: a request that crosses it fails (error E56 is the reason for this rule).
3. Station terms, thresholds, hyperparameters and model weights are estimated on development data only. Terms used to
   score an event are estimated without that event (error E53).
4. Synthetic seeds are drawn once per registered test and never reused; follow-up analyses on existing windows are
   labelled as such and are not confirmatory.

## 3. Comparators

1. Every comparator is run at its best: the strongest available configuration of each published method, its
   hyperparameters tuned on development data with the same budget as the proposed method, and the same inputs where
   the methods allow it (the same velocity model, the same station terms of the matching type).
2. Station terms are typed: delays of matched-filter energy peaks are used only by readers of energy peaks; delays of
   onsets only by pick-based methods (error E54).
3. The comparator set of a test is fixed in its registration and includes, unless the registration states why not:
   the best neural picker available with the best associator, template matching where templates exist, and the
   previous version of the proposed method.

## 4. Endpoints and statistics

1. The primary endpoint is recall of catalogued events at a false-detection rate fixed in the registration, by exact
   FROC, with the hit rule stated in every table (epicentre rule or origin time required).
2. A deployable operating point fixed on development data is reported next to the oracle operating point chosen on
   evaluation, for every method.
3. Decision rules are checked by simulation before registration: false-positive rate at the null and power at the
   effect of interest, for the actual number of targets and the dependence structure of the data.
4. Duplicate detections of an event that is already hit are reported separately from false and neutral detections.
5. Intervals respect the dependence of the data (blocks of time, networks, shared noise hours); the number of blocks
   is reported and a design with fewer than 30 blocks reports the interval as indicative.

## 5. Code

1. Confirmatory runs use the package `src/ttmt/` at a tagged version. Exploratory code may live in
   `experiments/exploratory/` but must import the package for every function that a later confirmatory test will use.
2. Conventions that change numbers (sample rounding, kernels, normalisations) are explicit parameters with a named
   default, never module-level state changed at run time (errors E51, E52).
3. Every change to the package is covered by a unit test or a regression test against frozen outputs; the regression
   suite reproduces the registered results of 2F, 2M, R0 and R1 before any release.
4. Paths come from configuration; analyses of archived outputs run without the raw data.

## 6. Registration and anchoring

1. The sealed log (`preregistration/original_sealed_es/`) stays append-only, in Spanish, with the line-aligned English
   translation and the hash chain, as before.
2. In addition, every registration, amendment and freeze list is committed and pushed to the public repository before
   the run. The GitHub commit time is the external timestamp; the run must not start before the push is confirmed.
3. The log entry of a result states the commit of the registration and of the freeze.

## 7. Verification

1. Verification of a confirmatory test is done by a process that did not write the code (a separate AI session) and,
   for tests on real data, includes at least one number re-derived by a person.
2. Verifications are described as what they are: re-derivations by a separate AI session are not human peer review.

## 8. Errors and deviations

Every error found, at any stage, enters the register (`docs/artifact_manifest_and_error_log.md`) with its cause,
detection, scope and correction, whether or not it changes a number. A deviation from a registration is logged before
the result is read whenever it is known before.

## 9. Assessment

Audits use the fixed rubric of the earlier audits (novelty and positioning 20 %, theoretical rigor 15 %, experimental
design 20 %, empirical evidence 15 %, internal coherence and writing 10 %, feasibility 10 %, reproducibility 10 %)
scored at doctoral level. The baseline of 18 October 2026 is 53.3/100 (audit 7, doctoral column). Audits by AI
sessions are recorded as such.
