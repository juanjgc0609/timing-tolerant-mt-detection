# Statistical analysis plan, phase 6 (and the phase-4 gate)

Version 0.1, 23 October 2026, statistical design thread. Draft: it becomes binding for a test when the text is copied
into that test's registration and anchored (ADR-0004). It applies the single beat rule of ADR-0008 to the two
confirmatory tests of phase 6 (registrations R2 and R3, drafts in `preregistration/drafts/`) and to the phase-4 gate
G4. Code: `src/ttmt/stats.py` (every decision rule) and `experiments/exploratory/P6-SAP/sample_size.py` (sample size).
The reviewer checks each analysis against `docs/protocol/SAP_CONFORMANCE_CHECKLIST.md`.

## 1. Tests covered

| Test | Data | Truth | Regimes tested | Status |
|:--|:--|:--|:--|:--|
| G4 (phase 4) | development: semi-synthetic (primary), real development (secondary) | injection; T3.V truth for real | T and C (ADR-0008) | gate, not confirmatory |
| R2 | real sequences reserved under T3.R4, from networks not used in development | dense or nodal network decimated to 4-12 station subsets; reviewed catalogue | T only where the truth resolves pairs at 0.1-1 s, C always | confirmatory |
| R3 | semi-synthetic: recorded events of the reserved networks injected into their own continuous noise | injection | T and C | confirmatory |

The doctoral claim needs H1 and H2 met in R3 and H2 met in R2 (H1 in R2 too when R2 can test T). The two
registrations answer different questions, each at its own level 0.025; the claim is their conjunction, so no
multiplicity adjustment between them is needed.

## 2. Hypotheses

As in ADR-0008, items 1 to 7: H1 superiority in the target regime T (arrival separation 0.1 to 1 s), then H2
non-inferiority by 0.05 in the complement C, each against every comparator (intersection-union), fixed sequence,
one-sided 0.025. Secondary Holm family at 0.025 after H1 and H2: superiority in [0.1, 0.5) s and in [0.5, 1) s;
non-inferiority in C by 0.03; the same rule at 10 and 100 false detections per hour; the same rule at the deployable
threshold (Section 7.6). Descriptive: every cell, the epicentre-only rule, duplicates, origin and location errors,
per-sequence estimates.

Estimand: the difference in recall, model minus comparator, over the targets of the regime, weighted to the real
prevalence of its cells (ADR-0008 item 4), at the operating point where each method has at most 30 false detections
per hour of evaluated record.

## 3. Data, partitions and eligibility

3.1 **Reservation.** The R2 and R3 sequences are chosen and sealed (T3.R4) before any of their waveforms is read: the
data card lists the files and their SHA-256 and is anchored. Only the catalogue and station metadata may be read
before the run, to compute eligibility, regimes and the expected number of targets.

3.2 **Causal split (CR-15).** Anything a method estimates from events of the test network (station terms, the
empirical Green's function dictionary, the moment-tensor prior, background statistics) uses only records that end
before the evaluation period starts, with a gap of at least one hour. Validation of such estimates at G4 leaves out
whole spatial clusters (events within 5 km) as well as time blocks, and reports coverage against the distance to the
nearest calibration event.

3.3 **Eligibility, decided before the run by a frozen script.** These rules exist because real data come incomplete:
- a station-hour is unusable when more than 5 % of its samples are gaps, flat fill (runs of identical samples of at
  least 1 s) or clipped values, or when the channel metadata (response, orientation) are missing; flat fill is
  replaced by gaps, never by zeros;
- a target is an event of the reviewed catalogue with magnitude at or above the completeness magnitude Mc of the
  catalogue, estimated by maximum curvature plus 0.2 on the catalogue alone, raised after large shocks by the
  time-dependent Mc(t) of Helmstetter et al. (2006), and with at least 4 usable stations over
  [origin - 10 s, origin + 60 s];
- events below Mc(t) or with fewer than 4 usable stations are excluded events: a detection matched to them by
  LOC+T-V is neutral (neither hit nor false);
- a time block enters only if at least 50 % of its station-hours are usable; excluded blocks are listed with their
  cause before the run.
A method that fails on a block (crash, timeout) scores no detection there (counted as misses); a sensitivity analysis
removes such blocks for every method.

3.4 **Imperfect truth (T3.V3).** A detection that matches no catalogued event is false in the primary analysis. In a
sensitivity analysis it is neutral when it matches, by LOC+T, an event of an independent enhanced catalogue produced
by none of the compared methods (named in the registration). A catalogue produced by template matching is never used
to score template matching.

3.5 **Spent data (CR-23).** A ledger (`data/cards/USE_LEDGER.md`, to be created by the data thread) counts the
configurations evaluated on each segment. A gate is evaluated only on a segment with no earlier configuration count;
the Ridgecrest development hour is spent (about 20 configurations, E53), so G-M5b and G4 use fresh segments.

3.6 **Pretrained weights (CR-22).** The registration lists, for every pretrained model, its training sets with region
and period; a model whose training data overlap the test network and period is replaced or its overlap declared.

## 4. Methods and tuning

Comparators and the equal tuning budget follow ADR-0008, items 6 and 9, and the protocol, Section 3. The tuning ledger
and the frozen configuration of every method are part of the freeze list. Every method receives the same waveforms,
station list, velocity model and typed station terms.

## 5. Hit rule

LOC+T-V (ADR-0008 item 2, `ttmt.stats.label_voronoi`), with the duplicator negative control of item 3 run on every
analysis. Labels are computed once on each method's full detection list (prefix property).

## 6. Regime assignment

Separation, signal-to-noise ratio and station subset are computed from the catalogue (R2) or the injection design
(R3, G4) and the frozen eligibility script, never from method outputs. Signal-to-noise ratio uses the P window of the
recorded event before injection (R3) or of the record (R2).

## 7. Analysis

7.1 **Operating point.** Per method and budget, the lowest score threshold with at most floor(budget x hours) false
detections over the evaluated record (`operating_point`); tied scores enter together.

7.2 **Paired block sums.** Blocks: in R3 and G4 one injected network realisation (one noise block of one network with
one station subset); in R2 consecutive time blocks of 1 hour, or of T/30 when a sequence holds fewer than 30 hours,
within each sequence. d_j and n_j as in ADR-0008 item 5 (`block_sums` with the prevalence weights).

7.3 **Test.** Block sign-flip p-value, exact when at most 16 blocks are non-zero, otherwise 100,000 random sign
vectors from seed 20261023 (`signflip_pvalue`); the clustered McNemar z of Durkalski et al. is reported next to it.

7.4 **Multiplicity.** Intersection-union over comparators (`iut`), fixed sequence H1 then H2 (`fixed_sequence`), Holm
over the secondary family (`holm`).

7.5 **Intervals.** Stratified block bootstrap, 10,000 replicates, seed 20261023, strata = sequence (R2) or network
(R3), operating points and the best comparator recomputed in each replicate (`block_bootstrap`); percentile 95 %
interval; indicative with fewer than 30 blocks or an atom above 5 % at a bound. The number of blocks is reported.

7.6 **False-alarm budget check.** At the deployable threshold (fixed on development), the realised rate of each method
passes when the quasi-Poisson score test rejects "rate at least 1.2 x budget" at one-sided 0.025 (`budget_check`;
the Wald bound under-covers for overdispersed counts, measured in `tests/test_stats.py`). The deployable secondary
hypothesis is claimed only for a method that passes. The check is stratified (CR-01): blocks without catalogued
events within 15 km and 4 s before the window, and windows with a parent event (M >= 2 within 3 to 15 km and 1 to
4 s), reported separately; for the model's e-value rule, R3 also reports the rate of S_mix >= log(1/p) on hypotheses
placed at 3 to 15 km and 1 to 4 s from injected events, and the guarantee of M7 is declared refuted on those
hypotheses when that rate exceeds 2p.

7.7 **Software.** The analysis is one script, frozen with the code, that reads the archived detection lists and
writes every registered number with its source line.

## 8. Stopping and futility

Each test is a single run on the full reserved data; no data are added after any result is read. One interim look is
allowed in R2 only, after the first reserved sequence, and only for futility, which is non-binding (it cannot raise the
error rate): if the conditional power of H1 (or of H2 when R2 cannot test T) at the registered effect, given the
observed blocks, is below 0.10, the author may stop and report the test as stopped for futility. Before the run, if the
eligible targets after Section 3.3 fall below the number that gives power 0.50, the run does not start and an
amendment is required. If accrual stops for reasons outside the data (access lost), the test is reported with the
blocks run and its achieved power, as underpowered.

## 9. Sample size

With `sample_size.py solve`, using the recall, method correlation rho and intra-block correlation icc estimated on
development hits (`sample_size.py estimate`), the targets per block and the share of targets in T of the reserved
catalogue. Planning effects: +0.05 in T for R3 and G4, +0.10 for R2; true difference 0 in C. Reference numbers are in
`experiments/exploratory/P6-SAP/sample_size_out.txt`; the registration quotes the run made with its own parameters.

## 10. Deviations, anchoring and verification

A deviation known before the result is logged before the result is read (protocol, Section 8). Registration,
amendments and freeze are anchored by tag, Release and OpenTimestamps, and R2 and R3 on OSF Registries (ADR-0004); the
author holds the OSF account and submits. Verification is a reimplementation by a separate AI session, so described
(CR-19), plus a number re-derived by a person: for R2, the hit or miss of ten discordant targets chosen by the seed
of ADR-0004 item 4.
