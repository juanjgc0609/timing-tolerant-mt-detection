# Statistical analysis plan, phase 6 (and the phase-4 gate)

Version 0.3, 21 November 2026 (v0.1, 23 October 2026, statistical design thread; v0.2, 14 November 2026, drops the e-value guarantee of 7.6 after G-M5c; v0.3 changes are listed in Section 11, with the ones that need the owner's approval marked). Draft: it becomes binding for a test when the text is copied
into that test's registration and anchored (ADR-0004). It applies the single beat rule of ADR-0008 to the two
confirmatory tests of phase 6 (registrations R2 and R3, drafts in `preregistration/drafts/`) and to the phase-4 gate
G4. Code: `src/ttmt/stats.py` (every decision rule) and `experiments/exploratory/P6-SAP/sample_size.py` (sample size).
The reviewer checks each analysis against `docs/protocol/SAP_CONFORMANCE_CHECKLIST.md`.

## 1. Tests covered

| Test | Data | Truth | Regimes tested | Status |
|:--|:--|:--|:--|:--|
| G4 (phase 4) | development: semi-synthetic (primary), real development (secondary) | injection; T3.V truth for real | T and C (ADR-0008) | gate, not confirmatory |
| R2 | real sequences reserved under T3.R4, from networks not used in development | dense or nodal network decimated to 4-12 station subsets; reviewed catalogue | T only where the truth resolves pairs at 0.1-1 s with origin-time error below a ninth of their separation (3.4), C always | confirmatory |
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

Ambition and sufficiency (ADR-0013, PR #72) [v0.3]. The two levels of ADR-0013 (ambition +0.05, sufficient more
than +0.02 over the best tuned comparator) describe the point estimate of the H1 difference with its interval; they
are not tests. The confirmatory claim is H1 (superiority, difference > 0, one-sided 0.025). "Sufficient" or "ambition
met" is never written as confirmed unless H1 is rejected; the report gives the zone of the point estimate next to the
H1 result. A +0.02 effect is not detectable at any feasible size: by scaling item 4 of the reference output by
(0.05 / 0.02)², about 6 200 targets in T, that is hundreds of thousands of catalogued events at the prevalence of
Section 9 (an approximation, not a solved design).

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

3.4 **Imperfect truth (T3.V3).** Truth precision in T: LOC+T-V counts a detection for one event of a co-located pair
only within a third of their separation of it (0.033 to 0.33 s in T). If the truth's own origin-time error is
comparable, the truth decides hits in T for every method and pulls the paired difference toward zero. Before the run,
the origin-time standard deviation sigma_t of the truth catalogue is measured on T3.V (against a higher-resolution
reference or by the scatter of repeated relocations), and a pair enters T in R2 only if its separation is at least
9 sigma_t; pairs closer than that are excluded events. If no pair qualifies, R2 does not test T.
A detection that matches no catalogued event is false in the primary analysis. In a
sensitivity analysis it is neutral when it matches, by LOC+T, an event of an independent enhanced catalogue produced
by none of the compared methods (named in the registration). A catalogue produced by template matching is never used
to score template matching.

3.5 **Spent data (CR-23).** A ledger (`docs/DATA_USE_LEDGER.md` [v0.3]) counts the
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
Serial dependence in R2 [v0.3]: the lag-1 autocorrelation of d_j over consecutive blocks of each sequence is reported.
If it exceeds 0.2 (pooled over sequences), the primary analysis uses blocks of 2 hours (or 2T/30), fixed by this rule
before any result is read; the 1-hour analysis is then reported as a sensitivity analysis.

7.3 **Test.** Block sign-flip p-value, exact when at most 16 blocks are non-zero, otherwise 100,000 random sign
vectors from seed 20261023 (`signflip_pvalue`); the clustered McNemar z of Durkalski et al. is reported next to it.
H1 uses the raw sum (exact). Non-inferiority tests (H2 and the secondary margin 0.03) use the clustered score test
`ni_score_pvalue` from the per-block discordant counts (`block_counts`): variance under the restricted null
(cluster-robust, floored at Tango's 1998 multinomial variance) and a Student t reference with Satterthwaite degrees of
freedom. They run at 0.025 unless the design check of ADR-0008 item 7 lowers the level: `sample_size.py calibrate`
over the registered range of recall, method correlation, intra-block correlation and mean and dispersion of targets
per block, whose ends come from development data (`sample_size.py estimate`); range and level are quoted in the
registration. The tests are conditional on the operating points of 7.1, which are chosen on all evaluated blocks; the
bootstrap of 7.5 recomputes them in each replicate.

7.4 **Multiplicity.** Intersection-union over comparators (`iut`), fixed sequence H1 then H2 (H2 at the level of the design
check), Holm over the secondary family (`holm`). Sizes are checked at the least favourable
configuration of the intersection-union test, one comparator at the boundary and the others clearly beaten, which
equals the size of one pairwise test.

7.5 **Intervals.** Stratified block bootstrap, 10,000 replicates, seed 20261023, strata = sequence (R2) or network
(R3), operating points and the best comparator recomputed in each replicate (`block_bootstrap`); percentile 95 %
interval; indicative with fewer than 30 blocks or an atom above 5 % at a bound. The number of blocks is reported.

7.6 **False-alarm budget check.** At the deployable threshold (fixed on development), the realised rate of each method
passes when the quasi-Poisson score test rejects "rate at least 1.2 x budget" at one-sided 0.025 (`budget_check`;
the Wald bound under-covers for overdispersed counts, measured in `tests/test_stats.py`). The deployable secondary
hypothesis is claimed only for a method that passes. [v0.3] The registration states the hours of evaluated record;
with fewer than 37 hours (the length that gives the check power 0.80 at overdispersion 5, `sample_size_out.txt`
item 5) the deployable hypothesis is not tested in that registration. The check is stratified (CR-01): blocks without catalogued
events within 15 km and 4 s before the window, and windows with a parent event (M >= 2 within 3 to 15 km and 1 to
4 s), reported separately; for the model's detection statistic S_mix, R3 also reports the rate of S_mix >= log(1/p)
on hypotheses placed at 3 to 15 km and 1 to 4 s from injected events as a measured false-detection rate. No
calibration guarantee is claimed (v0.2): S_mix was found not calibrated on development (G-M5c, PR #45) and the cause,
marginal or dependence, has not been separated (T4.27); false detections are reported as FROC and measured rate only.

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
catalogue. Planning effects: +0.05 in T for R3 and G4, +0.10 for R2; true difference 0 in C.
Feasibility [v0.3]: +0.05 in T needs about 1 000 targets in T at power 0.80 (item 4 of the reference output: 992 at
rho 0.7, icc 0.1, f 0.3). With the close-pair prevalence measured on a development catalogue (CR-04, 1.42-1.67 %, a
strong lower bound because of catalogue dead time; 16.8 % under the Poisson reference) that is between about 6 000 and
70 000 catalogued events, before the truth-precision and station rules of 3.3-3.4 remove more. So +0.05 in T is
planned for R3 (and the simulated gate G4), not R2; R2 plans +0.10 and, by Section 8, does not test T if its eligible targets give power below 0.50.
The comparator recall assumed here (0.40) is above the measured one (0.215-0.233 in T with station terms, tune
partition); `sample_size.py solve` is re-run with the recall, rho and icc estimated after pilot P1, before any
registration is anchored (rule 30). Reference numbers are in
`experiments/exploratory/P6-SAP/sample_size_out.txt`; the registration quotes the run made with its own parameters.

## 10. Deviations, anchoring and verification

A deviation known before the result is logged before the result is read (protocol, Section 8). Registration,
amendments and freeze are anchored by tag, Release and OpenTimestamps, and R2 and R3 on OSF Registries (ADR-0004); the
author holds the OSF account and submits. Verification is a reimplementation by a separate AI session, so described
(CR-19), plus a number re-derived by a person: for R2, the hit or miss of ten discordant targets chosen by the seed
of ADR-0004 item 4.

## 11. Changes from v0.2 (21 November 2026)

| # | Section | Change | Kind | Owner's approval (rule 14) |
|:--|:--|:--|:--|:--|
| 1 | 9 | Feasibility paragraph: +0.05 in T planned for R3 only; R2 plans +0.10; recall assumption to be re-estimated after P1 | states what the numbers already imply; no threshold changes | not needed (wording); owner informed |
| 2 | 2 | ADR-0013 levels are descriptive zones of the point estimate, not tests; "sufficient" never claimed as confirmed without H1 | interpretation of a decision already approved (ADR-0013) | **needed** |
| 3 | 9 | Re-run of the sample size with measured recall after P1 | already required by rule 30 | not needed |
| 4 | 7.2 | R2 serial-dependence rule: lag-1 autocorrelation of d_j > 0.2 -> 2-hour blocks primary | new analysis rule, fixed before any R2 data | **needed** |
| 5 | 3.5 | Ledger path `docs/DATA_USE_LEDGER.md` | cohesion | not needed |
| 6 | 7.6 | Deployable hypothesis not tested with fewer than 37 hours of record | new condition, fixed before any R2 data | **needed** |

No data of R2 or R3 has been read; every change is made before any result it could depend on.
