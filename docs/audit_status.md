# Audits of the work and plan to strengthen it

## Audits 6 and 7 (17 and 18 October 2026): the work judged as an undergraduate thesis

After the author stated that the work is an undergraduate thesis (trabajo de grado), two further audits applied the
same fixed rubric calibrated for that level, each also giving the doctoral-level score for comparison. Like the
earlier ones, they were made by separate AI sessions that had not produced the work and did not read the previous
audits; they are not human committees. Audit 6 (`audits/independent_audit_6.md`) gave 67.5/100 (53.5 at doctoral
level) and required a statement of how AI tools were used, a task-by-task account of the protocol, and corrections of
claims on seeds, on the theory and on R1. After those corrections, audit 7 (`audits/independent_audit_7.md`) gave
70.3/100 (53.3 at doctoral level), approval with text corrections, and no distinction; its corrections were applied
in commit 3b54c75. The rise between the two audits comes from corrections of the text and from disclosure, not from
new positive evidence.

| Criterion | Weight | Audit 6 | Audit 7 |
|:--|:-:|:-:|:-:|
| Novelty and positioning | 20 % | 6.5 | 6.0 |
| Theoretical rigor | 15 % | 6.5 | 7.0 |
| Experimental design | 20 % | 7.5 | 8.0 |
| Empirical evidence | 15 % | 6.5 | 6.5 |
| Internal coherence and writing | 10 % | 6.0 | 7.0 |
| Feasibility | 10 % | 6.5 | 7.0 |
| Reproducibility | 10 % | 7.5 | 8.0 |
| Total (undergraduate) | | 67.5 | 70.3 |
| Total (doctoral comparison) | | 53.5 | 53.3 |

## Independent audits

Between 26 September and 2 October 2026 the proposal underwent four independent audits, and on 16 October a fifth, carried out by
reviewers who did not take part in producing the work and who did not read the previous audits or this
document. All four use the same fixed rubric (novelty and positioning 20 %, theoretical rigor
15 %, experimental design 20 %, empirical evidence 15 %, internal coherence and writing 10 %,
feasibility 10 %, reproducibility 10 %; total equal to the weighted sum times 10) and checked figures
against the results files. The full texts are in `audits/independent_audit_1.md`,
`audits/independent_audit_2.md`, `audits/independent_audit_3.md` and
`audits/independent_audit_4.md` and `audits/independent_audit_5.md`. They are, to date,
the reference evaluation of the proposal; the scores in the following sections (60, 63 and
68 out of 100) are internal, were assigned by the same process that produced the work, and are kept only as
history.

**Table A. Score per criterion in the four independent audits (0 to 10).**

| Criterion | Weight | Audit 1 | Audit 2 | Audit 3 | Audit 4 | Audit 5 |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| Novelty and positioning | 20 % | 5.0 | 5.0 | 4.0 | 5.5 | 4.0 |
| Theoretical rigor | 15 % | 6.0 | 6.0 | 5.5 | 6.0 | 5.0 |
| Experimental design | 20 % | 6.0 | 5.0 | 6.0 | 7.0 | 6.5 |
| Empirical evidence | 15 % | 5.0 | 4.0 | 4.5 | 4.5 | 4.0 |
| Internal coherence and writing | 10 % | 4.0 | 4.0 | 4.0 | 5.5 | 4.0 |
| Feasibility | 10 % | 5.0 | 5.0 | 5.0 | 5.0 | 5.0 |
| Reproducibility | 10 % | 8.0 | 7.0 | 7.0 | 6.5 | 7.0 |
| **Rubric total** | | **55.5, rounded to 56/100** | **51/100** | **51.1, rounded to 51/100** | **57.75, rounded to 57.8/100** | **50.5/100** |
| Simulated decision | | major revision | major revision | major revision | major revision | major revision |
| Estimate by a reasonable committee | | not estimated | about 60/100 | about 62/100 | about 70/100 | about 53/100 (50 to 58) |

**Audit 5 (50.5/100, 16 October 2026).** It evaluated the version with R0 (negative), Task 2M and exploration C18. It
checked 28 numbers (24 matched exactly; the four mismatches were attribution, staleness and inventory errors, corrected
the same day) and the full hash chain (28 hashes, all byte prefixes). Its lower score than Audit 4 reflects the first
real-data evidence: the primary hypothesis G1 failed on Ridgecrest and the text had not been restructured around it;
little method novelty survives the project's own ablations (BOMP-PK with the same peak search matches MT-FW, the
convex method was refuted, and clipping is known practice), and the theory covered the refuted convex method rather
than the deployed statistic; the 2M comparator, BOMP-PK without any per-station robustness, is weak, and the primary
day was chosen knowing its noise had transients; the bounded statistic confirmed in 2M is the one that saturates on
real earthquakes; and the documents are too long. Its top actions: restructure the proposal around the evidence,
write the theory of the deployed statistic (begun in Appendix AA.9), re-plan the first year around the
arrival-reading gate (done in Section 19), test G1-R against a comparator with per-station robustness, and pass the
gate and run the registered real test on the Hi-net swarm.

**Audit 1 (56/100).** It evaluated the version with Tasks 2A and 2B already executed. It found all the
checked figures in agreement with the files and valued the error log, the preregistration and the
bit-for-bit controls, but pointed out that G1 rested only on synthetics from a generator of the same code,
that the closest prior work had not been run (the block orthogonal matching pursuit of Vera Rodriguez et al.
[129] and mechanism-corrected stacking [122], [133]), that the method defended, the group BLASSO,
was not implemented and all the evidence was from MT-FW, that the inference did not use the network as the unit of
replication (five networks per cell), that the document contradicted itself on the primary question and the
primary endpoint, that ablation 2B did not isolate the sign of the radiation pattern, that the pick-based
associators were left without picks by an out-of-domain STA/LTA, that there was no deployable false-alarm
control and that the scope was disproportionate.

**Audit 2 (51/100; committee, about 60).** It evaluated the version with Task 2C. Its main weaknesses
were the exclusively synthetic evidence outside the declared regime; the comparators at a
structural disadvantage, in particular the block pursuit and the stacking read at the exact time
of a coarse grid, without the timing tolerance of MT-FW; the operating points chosen by oracle
on the evaluation; a convex framework that gave the thesis its title and performed below the deployed greedy
detector; the preregistration time of 2B (14:30 UTC) later than its first files; and the misalignment
between the title, Sections 11 and 13 and the budget, on the one hand, and G1, on the other.

**Audit 3 (51/100; committee, about 62).** It evaluated the version with Task 2D and the separate
supplementary volume. Its critical weakness is one of positioning: with the same peak search, BOMP-PK
matches MT-FW, so the announced methodological contribution is not distinctive, and the distinctive one, the group
BLASSO with off-grid station time shifts, is not implemented and has
uncited prior work (the Cut-and-Paste method of Zhu and Helmberger, 1996, and the joint inversion of
moment tensor and station shifts of Pham, 2024); even so, the title ("with a Convex
Certificate") and G1 described the certificate of the convex problem. It added the absence of a deployable
operating point (1.2 to 2.2 false events per window with events with conformal thresholds), the weak
statistical support (no 2D rule survives Holm with five networks), the integrity of the preregistration
verifiable only for 2D, with an undeclared edit of `analyze_2b.py` 18 s after its first
output, the risk of transfer to the RSNC (reliable polarization at 1 of 5 stations), the gap between
the theory of Appendix AA and the deployed detector, and 21 inconsistencies with line numbers, among them
a misrepresented citation ([16]).

**What was done in response.** Tasks 2A and 2B had been executed before the first independent
audit, which evaluated them. In response to the first, Task 2C was preregistered and executed: the block
pursuit and the mechanism-corrected stacking on the same 1200 windows, the on-grid group BLASSO
solved with a certified duality gap, and the comparison with QuakeMigrate replicated over 30
networks per cell with hierarchical bootstrap and Holm; the result was that the on-grid convex problem
falls below MT-FW in 23 of 24 cases (E33), and the document was restructured around a single
primary hypothesis and a single endpoint, with an Abstract of at most 300 words, a plan and publication
units that contain the detector, and Appendices G to AA in a supplementary volume. In response to
the second, Task 2D was preregistered, the first with SHA-256 hashes of the registration and of the code before
any result: it gave the prior methods the same peak search, ablated the peak search and the
cap of MT-FW and froze conformal thresholds; rule 2D-a failed (BOMP-PK matches MT-FW, E35), the
advantage was attributed to the peak search, and the erratum on the time of 2B was recorded (E34). In
response to the third, without new experiments, the proposal was repositioned: the title became
*Timing-Tolerant Moment-Tensor Detection of Overlapping Earthquakes: A Group BLASSO with Station Time
Shifts*; G1 concerns the family of timing-error-tolerant moment-tensor detectors with
deflation (MT-FW and BOMP-PK) against published associators and migration-based locators; a
methodological hypothesis, G1-M, untested, requires that the group BLASSO with station shifts
solved off the grid outperform BOMP-PK with the same tolerance; Cut-and-Paste and
Pham (2024) are cited as prior work ([143], [144]); citation [16] and the 21 inconsistencies were corrected; Tasks
2B and 2C were reclassified as registered in a local append-only file with a date not
externally verifiable, and the edit of `analyze_2b.py` was declared and it was checked, by reconstructing the
previous script, that it changed no rule or figure (E36).

**What was done between the third and the fourth audit.** The methodological hypothesis G1-M was preregistered as
rule 2E-a, with the hash of the registration communicated before writing code, and Task 2E refuted it: the group
BLASSO with station and phase shifts, solved off the grid, gained between +0.000 and +0.039
over BOMP-PK with all intervals containing zero, and an independent verification found no
implementation error or budget disadvantage (with ten times more iterations the difference was −0.014). An
exploration on the calibration seed (candidates C1 to C8, declared exploratory and with its multiplicity)
searched for where margin remained, and the candidates for the next two tasks came out of it. Task 2F, on two new seeds and with 96 networks per
cell, confirmed that the per-station cap of MT-FW raises recall in dense sequences of three to five
events (+0.043 to +0.074 over BOMP-PK), and its verification narrowed the result: +0.024 to +0.045 if
the origin time is required, a gain tied to the level of the cap and a mechanism of score calibration, not of
event separation. Task 2I, on another two new seeds and with 120 networks per cell and regime,
confirmed the same with a Huber fit whose level is set by rule from the noise (+0.058 and +0.053,
non-inferior to the cap and equivalent to BOMP-PK in standard windows when pooling regimes), and its verification
added that with origin time the gain is +0.031 and +0.019, that equivalence does not hold with
weak clutter on one seed, that the fit costs about nine times more and that it does not find more events.
The secondary arms of 2I, of level and of dispersion, were still running at the time of the fourth
audit and nothing is reported about them.

**Audit 4 (57.8/100; committee, about 70).** It evaluated the version with Tasks 2E, 2F and 2I. It reproduced
the twelve hashes of the chain and checked fourteen figures against the results files, all fourteen
matching. Its verdict is that the main defect is one of substance and not of honesty. All the evidence
for the primary hypothesis is synthetic, from a single generator and outside the declared regime; G1, as
formulated, does not hold at 0.1 against GaMMA with optimistic picks, and there is no deployable
operating point; the advantage over the state of the art rests on a heuristic peak search that the
prior work [129] matches as soon as it is granted it, and the comparison with the associators uses a picker that
leaves them without picks; the title and G1-R called "certificate" a statistic that the text itself declares
is not one; the effect of G1-R is small and its equivalence margin (0.05) is larger than its
superiority threshold (0.03); the only real datum, reliable polarization at 1 of 5 stations of the RSNC, goes against
the assumption of the radiation operator and the plan did not address it for the primary detector; the
scope is disproportionate and Sections 5 to 12 play no part in any confirmatory test; the theory of
Appendix AA covers the convex problem and not the deployed statistic; the preregistration is local, with the
hashes stored next to the registration and without a third-party timestamp; and the length, about 130 000 words across
the two volumes, prevents evaluation. It listed fourteen inconsistencies with line numbers and no
numerical discrepancy between the text and the files.

**What was done in response to the fourth.** Text changes only, without experiments. The title became
*Timing-Tolerant, Bounded-Influence Moment-Tensor Detection of Overlapping Earthquakes: A Preregistered
Benchmark and a Real-Data Test Plan*, and "certificate" is reserved for the dual certificate of the
convex problem; the deployed object is called the bounded-influence moment-tensor subspace statistic
(E40). The one-sentence thesis describes what passed the tests, a greedy search with timing tolerance and
bounded influence, and leaves the convex problem as the framework. The sentence on the calibrated constant is restricted
to the registered primary arm (E41). Section 16 declares that G1' depends on amplitudes and polarities and
registers an arm with only the stations that pass the orientation and polarization audit and a sign-blind
fallback arm (E42). The following are declared: the difference in thresholds between G1 and G1-R and its reason, the hit
rule with origin time (origin error of at most 0.5 s) as primary in every future task
with equivalence margins of 0.02 per regime, a single multiplicity rule per family (Holm for the
primary endpoint, Benjamini-Hochberg for the secondary ones), and the local nature of the preregistration with a
public repository and timestamped registration as a deliverable of the first semester; and the
inconsistencies of date, of the count of abandoned claims and of tasks, of family E, of package
E4 and of the freezing of margins were corrected. Reference [9] was verified on the journal's page and is kept.
Neither the length nor the scope was cut, and the Abstract remains at 297 words.

**After the fourth audit: Task 2L (5 October 2026).** The fourth audit (57.8/100) evaluated the
version with Tasks 2E, 2F and 2I and predates exploration C11 and Task 2L, which confirmed on two
never-used seeds the calibration hypothesis G1-C (station terms estimated without picks inside the
refined detector; Appendix Y.28), with independent verification and a narrow scope: gain
concentrated in the one-sided cell of 8 stations, nothing demonstrated with 4 stations, conditional on the
refinement and on catalog errors of known size. No new scoring with the rubric has been done:
57.8 remains the latest audited figure and does not include 2L. The evidence of 2L is also synthetic and from the
same generator.

**After the fourth audit: explorations C12, C13 and C14 (5 to 9 October 2026).** After 2L,
three explorations were carried out on the calibration seed, without a registered hypothesis and with the
registered detector of the previous tasks (Appendix Y.29): the restriction of the subspace statistic to double
couples, whose null law was derived with the tube formula and which does not change recall; a calibration of the threshold
by decoys with the station moveout permuted, which reproduces without labels the operating point
under the epicenter-only rule and not under the origin-time rule; and remedies for self-interference after
deflation, negative for recall, whose diagnosis places the largest class of false events with origin time in
well-located detections of real events with the origin shifted. None changes a registered rule or
meets any of the six pending requirements for 75/100; in particular, the decoy calibration is
not yet the deployable operating point that is required, because it is from a single seed, has no guarantee and
does not hold under the primary hit rule. No new scoring with the rubric has been done for these
explorations either: 57.8 remains the latest audited figure.

**After the fourth audit: C15, C16 and the first test with real data (9 and 10 October 2026).**
C15 (exploratory) showed that the origin error of well-located detections is a location error and that
no second-stage estimate nor cyclic re-estimation corrects it in this generator. C16 (exploratory,
real noise of 23 July 2010 from the Alto Tiberina network with injected synthetic events) strengthened G1-R:
with real noise BOMP-PK collapses and the cap and the Huber estimator do not, although part of the transients turned out to be
telemetry padding. This advances the second requirement (second noise model), without meeting it: the
source generator is the same. Task R0, registered before touching the data, was the first test with real data (one
hour of Ridgecrest 2019) and gave a clear negative result: −0.383 [−0.481, −0.310] against the best associator with
PhaseNet, due to timing error of the homogeneous model (0.35 to 0.41 s) and to saturation of the score in a
continuous sequence (Appendix Y.30). It meets in spirit the first requirement (a real result, whatever
its sign), although not on the RSNC. No new scoring with the rubric has been done; an auditor would now have
real evidence that the synthetic advantage does not transfer without a better travel-time model, which lowers the strength of
G1 and raises the credibility of the program.

**What is missing for 75/100.** The target of 75/100 is still not reached on the rubric: the fourth audit
gives 57.8, estimates that the text and timestamping changes are worth about 5 points (up to about 63) and that with
experimental work about 73 would be reached. What is missing is, in this order, the first result of G1' on the
decimated RSNC, with frozen threshold and false events per day, whatever its sign; a second independent
generator (another synthesis code, real noise, sources outside the double couple) for G1 and G1-R; a deployable
operating point, calibrated on distinct networks over windows with events; a theoretical result
on the deployed statistic, not on the convex problem; the reduction of the scope and the length to
a three-year core with at most six hypotheses and a body of about 12 000 words; and a
public registration timestamped by a third party. None of the six is done. Of the requirements set by the third
audit, the first was met with a negative result (Task 2E) and the others remain pending; the text
that follows keeps them as they were written then. The third audit estimates that the
rewriting and the repositioning, without new results, are worth about 5 points (up to about 56), and that with
new results about 66 would be reached. Reaching 75 requires evidence that does not yet exist. First,
implement the group BLASSO with station and phase shifts, solved off the grid by sliding
Frank-Wolfe with debiasing refit, extend the theory of Appendix AA to that model and show
that it outperforms BOMP-PK with the same tolerance on the windows of 2A and over 30 networks per cell, under a
preregistered rule with an external timestamp (signed commit in a public repository, OSF or RFC 3161 timestamp),
reporting the null if it does not outperform it. Second, execute Tasks 1p (real noise) and 1v (3D Green's
functions) and add sources outside the double couple. Third, obtain a deployable operating point: a threshold
calibrated on distinct networks, with windows that contain events, that meets a budget of false events per
window counted over all windows. Fourth, run a real pilot of G1' on a subnetwork of the RSNC
with audited orientation, threshold frozen on noise windows and false events per day, with MT-FW and BOMP-PK
against QuakeMigrate and MALMI. In addition, the three audits ask for a much shorter main document
(the third proposes 15 000 to 20 000 words; today it has about 44 000) and for moving families F, C-II and C-V,
the learned residual and optimal transport to future work, and that has not been done either.

---

## Previous internal audit (history)

The following sections are the internal audit of 22 and 23 September 2026 and its addenda, carried out
within the same process that produced the work. Their scores are internal and were superseded by the
independent audits above.

**Audited document:** `Doctoral_Proposal.md`
(version of 22 September 2026, with Appendices T.1 and T.2), together with `tasks_and_real_data_plan.md`,
`executive_summary.md` and `task1_step_by_step.md`.

**Interpretation of "Task 1 has not yet been executed".** The registered Task 1, which called for several
stations, two regimes (background and active sequence) and the band of the proposal, has been executed
only in its minimal part: one valid station (ZAR), background regime and 1-10 Hz band. This audit
treats everything that goes beyond T.1 as not executed (1c, 1d, 1f, 1g, 1h, 1i, the active regime and
the second domain) and does not presuppose any result of those tasks.

## Addendum of 23 September 2026: corrections applied and results of the controls

This addendum summarizes what was done after the audit. The body of the document, from Section 1 on,
is the audit of the previous version and is kept as is.

**Corrections applied to the proposal.**

- The twelve contradictions of Section 6 were resolved: the count of pre-studies, the
  evidence taxonomy, the references to equations, the status "Proposition at leading order",
  the broken sentence of (8c), the threshold of B3, the corollaries, the claims of 13.0a and of
  Section 15, the replication across seeds and the curve metric.
- The Abstract was rewritten in about 380 words, without selective claims.
- The literature was expanded with 20 verified references ([88] to [107]) and a new Section 3.2b. The
  positioning with respect to NET-VISA was made explicit: it already uses per-station azimuth likelihoods
  and a model of false arrivals.
- C-IV is formulated as a per-station likelihood ratio, equations (9d)-(9e), with the reliability
  $\hat r_s$ estimated with the same algebra as $\eta$.
- The following were registered: the real pilot by decimation, the calibration sequence of C-IV (orientation
  audit, incidence control, band, catalogs and real PhaseNet residuals), four
  new risks and the baseline tuning protocol.

**Controls executed (Appendix Y).** The decision rules were written before each run,
in `PREREGISTRO_1j_1k.txt`.

| Control | Result |
|:--|:--|
| 1j, ties (K1 of this audit) | Confirmed. Three ablations that break ties without polarization match the prototype in one-sided networks, on both seeds, and a continuous residual kernel outperforms it with 3 and 4 stations |
| 1j-b, pick and velocity error (exploratory) | Polarization added to the kinematic kernel contributes in 1 of 9 cells (4 stations, 0.2 s): 3 versus 19 false events on one seed, 3 versus 6 on the other, at a cost in recall. Rule not met |
| 1k, empty windows and absolute threshold (K2) | The absolute threshold solves the problem of empty windows for the well-calibrated methods (at most one emission in 20-32 empty windows). The best kinematic score is never worse than those that use polarization, so G1' does not hold |
| Baselines tuned on calibration (K3) | GaMMA competes in the one-sided and spread cells; neither GaMMA nor PyOcto reaches the target recall in the sparse ones |

**Same-day update: the negative was also audited, and it did not hold up.** Faced with the question of whether
the negative was an implementation issue, the generator was audited with the same suspicion as the
positives before. It detected all arrivals with exact time even at signal-to-noise of 1.1 to 2.5, whereas the
polarization at that signal-to-noise separates with AUC 0.78; moreover, the combination was not a likelihood.
With a consistent generator and the likelihood (9d) (Task 1o, Appendix Y.5), the first reading,
at the nominal budget, gave 0.18 to 0.25 of recall over the kinematic likelihood and an advantage over
GaMMA and PyOcto. That reading did not hold up under audit either: budget overrun, no negative
control, weak kinematic reference and baselines with another domain. Repeated with a shuffled-polarization
control, the realized false-event rate, an exact FROC curve, GaMMA with seed and 200
windows per cell (Appendices Y.6 to Y.9), polarization seemed to add detections within the
search and on the waveform. A second independent verification found that the generator
radiated every S wave as pure SV, the direction that those methods predict, and that the shuffled control was not
null. With double-couple sources and controls that take the geometry of another hypothesis (Appendix
Y.10), the conclusion narrows again: the per-pick likelihood does not outperform its control; the
migration with the constraints that hold for any source (P along the ray, S transverse)
does, by +0.05 to +0.19, and it outperforms direction-free migration by +0.07 to +0.21 and GaMMA by +0.12 to
+0.35 at 0.25 false events per window (a margin later found to be inflated by E26), on two seeds and in a replicate with new windows; as a filter
on the GaMMA catalog, polarization adds nothing. The certified pruning of the hypothesis space (Appendix Z)
gives catalogs identical to the dense search with 2.3 to 3 times fewer evaluations. The following paragraphs
describe the state before this update.

**What this meant for the thesis before Task 1o.**

- The central advantage that the previous version reported was an artifact of ties. In the
  synthetic generators, polarization widens the discrimination margin but has not been shown to
  reduce false events below what the best kinematic score achieves.
- G1 was reformulated as a hypothesis on real data. The generators give the kinematic score more than
  real data will give: exact picks or picks with Gaussian error, exact model, random clutter.
- The priority becomes measuring the real errors of PhaseNet on the RSNC (Task 1n) and the pilot
  by decimation.
- What does not depend on that question remains intact: the observability certificate and the
  per-station reliability methodology. The ladder is kept as an organizing idea, indexed
  by effective kinematic information, but its intermediate rung has to be measured again against the
  best kinematic score.

**Internal score of the corrected version: 63/100, major revisions**, for reasons different from those of
before:

- Internal coherence goes from 4.5 to 7.5.
- Experimental design goes from 5.0 to 6.5, due to the controls, the absolute threshold and the tuned
  baselines.
- Reproducibility goes from 8.0 to 8.5.
- Novelty drops from 7.0 to 6.0, because NET-VISA narrows the gap.
- The empirical evidence of the core drops from 5.0 to 4.0, because the main advantage turned out to be an
  artifact.
- Theoretical rigor (7.0) and feasibility (5.5) do not change.

The proposal is now honest about its central question, and that question remains open. A third
verification found that since Task 1o GaMMA had been receiving undetectable clutter picks (E26), so
its earlier margins were inflated. With consistent clutter (Task 1y), the reformulated method
(sparse recovery of a measure with moment-tensor amplitudes, solved by Frank-Wolfe)
outperforms GaMMA by +0.11 to +0.28 with weak clutter and by +0.17 to +0.43 with strong clutter, and also its
control with the geometry of another hypothesis. Score changes:

- Empirical evidence goes from 4.0 to 5.0: a positive with controls, consistent clutter and new
  windows, fragile to delays larger than the window and still without the published migration-based
  locators.
- Experimental design goes from 6.5 to 7.0: controls with the same degrees of freedom, exact FROC
  curve, deterministic baseline and registered rules.
- Theoretical rigor goes from 7.0 to 7.5: the method is now a convex problem over measures with
  a dual certificate, certified oracle and random-field threshold, all with established theory.
- Novelty goes from 6.0 to 6.5, as a combination. Moment-tensor migration exists in
  dense arrays and at long periods (Chambers, Anikiev, GRiD MT). There is no certified search nor
  threshold for the maximum over the hypothesis space in seismology, as far as was found.
- Feasibility goes from 5.5 to 6.0: the core is smaller than the optimal transport and beams plan.

The internal score stands at **68/100** (6.3 + 0.15 x 1.0 + 0.20 x 0.5 + 0.15 x 0.5 + 0.20 x 0.5 + 0.10 x
0.5), still in major revisions until the method is tested with a real picker, against the published
locators and on real data.

---

**What was verified in addition to reading.** The figures of T.2 were recomputed from
`results/task1e/*.json`. Those of X.2 and X.3 were checked against the original logs
`v6_onesided.txt` and `v7_onesided.txt`. The catalog emission logic and the relative threshold logic were read in the version of
`task0_common.py` used by Task 0 and T.2 (`catalog_from_score_field(..., relative=True)`). The verifications of the
estimator $\eta$ were read in `results/task1e/verification.txt`. The calculations of Tables 3 and 5 and
of equation (25) were redone by hand and are correct.

---

## 1. Editorial verdict

**Internal score: 60/100. Major revisions before the proposal defence.**

The proposal has virtues that a committee will recognise. It cultivates preregistration and preserves
negative results. Appendix W audits the statistics of the pre-studies a posteriori. Appendix X.1
publishes the two mistaken protocols in addition to the correct one. In T.1 the chance correction of
the coherence is rigorous, and both the per-station validity and the protocol deviations are declared
honestly. On the theoretical side, the observability certificate (8a)-(8d) is correct and verifiable
algebra, and the Lipschitz bound (9a) is independent of the medium. All of this is above the average
of a doctoral proposal.

The diagnosis is different when one reads the data that support the central hypothesis. G1 appears in
the master table as "supported on both endpoints", but:

1. Against the strongest published associator, the advantage in false-event rate does not exist in
   the cells where G1 is declared. In the one-sided cells GaMMA obtains recall 30/30 and 0/30 false
   events in both seeds, the same as the constrained prototype. McNemar gives $p=1.00$ in all four cases.
2. The advantage over the ABL ablation may be due largely to the ties of a discrete score and not
   to the physics (Defect K1).
3. Protocol P3, which normalises the threshold by the maximum of each realisation, cannot operate on
   windows without an event. It is therefore not transferable to real continuous data.

None of this invalidates the programme, but it forces a reformulation of G1 before investing in real
data, because the real pilot (1g) would inherit all three problems.

Added to this is a credibility problem of lesser scientific weight and high cost before a committee:
the document contradicts itself in at least twelve verifiable points (Section 6 of this audit), some
of them in central claims.

---

## 2. Evaluation by criteria

| Criterion | Weight | Score (0-10) | Justification |
|:--|:--|:--|:--|
| Novelty and positioning | 20% | 7.0 | The reformulation "physics buys rejection, not resolution" is original and is supported by the sufficiency ladder. SIGVISA, the closest antecedent, is named, but the comparison with it has not been executed |
| Theoretical rigour | 15% | 7.0 | Corollaries 1-3 correct. The Proposition is declared "proved" with a derivation of the leading term, without bounding the remainder $O(V/(f_{\max}d_{\min}))$ |
| Experimental design | 20% | 5.0 | The relative threshold, the one-source-per-window generator and the ABL ties make the catalogue endpoint measure something different from what is declared (K1, K2) |
| Empirical evidence | 15% | 5.0 | A single generator and a prototype that is not the model. The real evidence is one station in a band that is not the band of the proposal |
| Internal coherence and writing | 10% | 4.5 | Violated evidence taxonomy, contradictory counts, wrong equation references and claims contradicted by the appendices themselves |
| Feasibility | 10% | 5.5 | 48 months with five families, optimal transport, Gaussian beams, learned residual and SO(3) equivariance. Section 23 declares a core much smaller than the plan |
| Reproducibility | 10% | 8.0 | Scripts, seeds, cross-checks and an errata manifest. Exception: X cannot be regenerated realisation by realisation (T.2 acknowledges this) |
| **Weighted total (internal)** | | **6.0** | **60/100** |

---

## 3. Critical defects

### K1. The advantage of G1 over ABL may be a tie effect and not an effect of the physics

**What the document says.** X.3: in the two one-sided cells the constrained prototype produces zero
false events in sixty realisations, against 0.30-0.43 for ABL, with $p\le0.0039$. This supports G1
in the master table (line 570).

**What the code shows.** `catalog_from_score_field` emits, in score order, with 5 km of separation
and up to five events, all candidates that exceed the threshold. In those cells the calibrated
threshold is $1.00$ of the maximum for ABL in all four cell-seed combinations, and for PHYS in three
of the four ($0.95$ in the remaining one; registrations v6 and v7, and the ideal condition of T.2).
This produces two different behaviours:

1. ABL scores a fraction of hits, a rational number whose maximum is exactly $1.000$ in each
   realisation. Several grid nodes reach it at the same time along the trade-off curve left by a
   one-sided network. All those nodes are emitted, and one more than 5 km away is enough to count a
   false event. In those cells ABL has recall 30/30 in both seeds: it never misses the true event,
   and its false events are additional tied events.
2. PHYS multiplies that score by a continuous term, so its maximum is almost surely unique and at
   threshold $1.00$ it emits exactly one event. In those cells its "false-event rate" is the error
   rate of the argmax.

The comparison therefore sets a set of tied maxima against a single argmax. Part of the difference
may come from breaking ties with any continuous term, not from the polarisation. There is a strong
indication in favour of this reading: GaMMA fits travel-time residuals without using polarisation and
also obtains 0/30 false events in those cells.

**Why it is critical.** If an ABL with ties broken by time residual, without polarisation, also gives
zero false events in the one-sided cells, the evidence for G1 in its current form disappears. The
contribution of polarisation would be restricted to the cells where ties do not dominate: the
four-station ladder, where PHYS gives 10/30 against 25/30, and the resistance to bad stations of
T.2.

**Correction.** A new synthetic control, "ABL-cont", with the existing code:
$s=s_{\mathrm{ABL}}-\epsilon\,\mathrm{RMS}(\Delta t)$, or a Gaussian kernel of residuals instead of
the hard tolerance. It is run with protocol P3, the same seeds and paired McNemar against PHYS. The
decision rule must be preregistered: if the difference in false events between ABL-cont and PHYS in
the one-sided cells is less than 0.05, G1 is reformulated as indicated in Section 5.

### K2. The catalogue endpoint, as implemented, does not measure what it claims to measure and is not transferable to real data

**What the document says.** Section 13.0b distinguishes $E_{\mathrm{loc}}$ (argmax) from
$E_{\mathrm{cat}}$ (emitted catalogue) and presents X as a measurement of $E_{\mathrm{cat}}$.

**What the code shows.** Each realisation has exactly one true source, and the threshold is a
fraction of the maximum of that same realisation. There are three consequences:

1. With threshold $1.00$ and a continuous score, $E_{\mathrm{cat}}$ reduces to $E_{\mathrm{loc}}$: a
   false event implies a missed event. The distinction that Section 13.0b presents as a finding is not
   exercised in three of the four cell-seed combinations where G1 is declared.
2. The unknown cardinality, which Section 7.1 places at the centre of the inverse problem, is never
   put to the test: no window has zero, two or three events.
3. **A window without an event always emits at least one**, because its maximum is by definition 100%
   of itself. In continuous data most windows are empty. P3 cannot be run on real data without
   changing the threshold, so the real pilot of G1 (Task 1g) cannot inherit the protocol that
   supports G1. The argument of X.1 ("thresholding a normalized statistic is standard") does not
   apply: the standard normalisation in detection (CFAR) divides by the background level, not by a
   maximum that contains the signal.

**Correction.** A generator with a number of events per window $K\sim\{0,1,2,3\}$ that includes
clutter-only windows. An absolute threshold calibrated on null windows for a target false-alarm rate
per hour (Neyman-Pearson criterion). Primary endpoint: the curve of recall against false alarms per
unit time, with the operating point of Section 14 read on that curve. It is the same endpoint that the
real pilot will have, and this links the synthetic programme with the real one.

### K3. G1 is not demonstrated against the state of the art, and the text claims the opposite

**What the document says.** Section 13.0a, lines 1858-1860: the constrained model "beats both
published associators decisively in well-covered networks under heavy clutter and against a distant
coherent decoy, and loses to GaMMA in one sparse-clutter cell". The Summary presents the three-station
recall (29/30 against 7-9) and the elimination of false events in one-sided networks as main
results.

**What X.2 and the registrations show.** In the two one-sided cells, where G1 is declared, GaMMA
has recall 1.00 and 0.00 false events in both seeds (registrations v6 and v7), the same as PHYS. In
"spread, clutter 40", "decoy 4 km" and "decoy 12 km", GaMMA has 0.00 false events in both seeds: the
sentence about the distant coherent decoy is false with respect to both baselines. PHYS does beat
GaMMA on the ladder: at six stations in both seeds (0.03 against 0.40 and 0.37; $p=0.001$ and
$0.006$), at four stations in one (0.43 against 0.73, $p=0.049$; in the other 0.30 against 0.60,
$p=0.078$) and at twelve stations in one seed, where GaMMA emits false events in 29 of 30
realisations, and in the other in none. Those advantages appear precisely where GaMMA fails at tasks
that the trivial ablation solves: at six stations ABL has recall 1.00 and GaMMA 0.67-0.77. This points
to a configuration of GaMMA unsuited to this generator, rather than to a physical advantage. At three
stations all methods emit false events in 60-100% of the realisations (PHYS 0.80 and 0.83): the recall
of 29/30 comes with a catalogue that almost always contains a spurious event, and the Summary does not
say so.

Moreover, Section 15 states "Travel-time baselines, none of which has yet been executed", which
contradicts X. And the baselines ran with default hyperparameters, except for the velocity model. The
oscillations between seeds in the same cell (PyOcto 1.00 against 0.10 and 0.97 against 0.00; GaMMA
0.97 against 0.00 at twelve stations) indicate an unsuitable configuration, not a property of the
methods: **WEAK BASELINE**.

**Why it is critical.** An external reviewer will read X.2 and the registrations and will find that
the sentence of Section 13.0a is partly false, that the remaining advantages over GaMMA rest on a
badly configured baseline and that the Summary reports selectively. This compromises the credibility
of everything else.

**Correction.** Rewrite Section 13.0a and the Summary with the data of X.2. Tune GaMMA and PyOcto
with the same calibration budget as PHYS: hyperparameter search on the calibration realisations and
freeze before evaluation. Always report recall and false events together. Run SIGVISA or justify in
writing why it cannot be run.

### K4. The validity of T.1 is narrower than its verdict

**The labelling depends on an incomplete catalogue.** The reviewed ISC catalogue is used, whose
completeness in Colombia is lower than that of the SGC catalogue. A pick from an uncatalogued
microearthquake is labelled as clutter. Moreover, $\eta$ measures coherence with the reference event
(the most recent catalogued one within 200 km), not with any physically realisable source. The
danger recorded by Appendix T is precisely a realisable impostor, such as an aftershock or a
neighbouring event. If that impostor is not in the catalogue, $\eta$ cannot detect it. The correct
verdict is narrower than the written one: at ZAR, over 90 days of background, the false picks are not
coherent with the preceding catalogued event. The sentence of W and T.1 "real false picks behave as
the randomly polarized clutter the pre-studies assume" is *Unproven speculative assertiveness* in its
general scope.

**There is a physical explanation of CHI that was not evaluated.** T.1 reports, for CHI, that the
near-vertical arrivals are the worst (0.136 with incidence below 20°) and that the S picks hit 0.677,
a figure it declares unexplained. But it does not connect the two facts, does not give the breakdown
by depth and does not report the control restricted to larger incidences. With near-vertical
incidence, the horizontal component of the P is small and dominated by noise, whereas the SV arrives
almost horizontally and in the radial plane. That would give exactly random P and S aligned with the
back-azimuth. If CHI records mostly intermediate-depth events, such as those of the Bucaramanga Nest,
the incidence geometry would suffice to explain it without invoking the installation. The text states
that "what remains is a property of the station, its installation or its local structure"; before
sustaining that conclusion one must report the control restricted to incidences of $30^\circ$ or
more and broken down by depth. If the hypothesis is confirmed, the reliability becomes a function of the station and of the
source-station geometry, not of the station alone. Moreover, a radial S channel could recover the
direction where the P fails.

**There was an uncorrected interim look and the estimator is verified with a single replicate.** ZAR
was looked at twice and the sample was enlarged because the first interval ended at the threshold:
this is optional stopping, which T.1 declares, without spending $\alpha$. `verification.txt` recovers
proportions of 0, 0.3 and 0.7 as 0.098, 0.223 and 0.776, with a single replicate each. Its section 4
repeats the case $\eta=0$ 20 times and measures how many times the verdict LOW comes out (5/20 with
50 picks, 20/20 with 300 or more), which bounds the verdict error. "$\eta$ recovers planted shares" is
correct within the interval, but neither the bias nor the coverage of the CI95 is characterised.

**Correction.** Relabel with the SGC catalogue and report the sensitivity of $\eta$ to the catalogue.
Search by template correlation for uncatalogued events that explain clutter picks. Stratify the
control by incidence and depth at CHI, OCA and ZAR. Calibrate $\eta$ by Monte Carlo (200 replicates
per proportion: bias and coverage of the CI95). Preregistered sequential design for 1c.

### K5. The gap between the prototype and the model, and between the synthetic and the real, grows instead of closing

All the empirical support comes from the same generator and from the discrete prototype of C-III and
C-IV (Section 13.0b acknowledges this). The proposed model, with transport misfit, beams, interior
point and five families, has no implementation at all. Three facts from the proposal itself narrow
the transfer:

1. The "waveform" rung of the ladder, at three stations, is an oracle matched filter: exact velocity
   model and correct wavelet (Appendix O).
2. Section 16 establishes that the public velocity models of Colombia do not support deterministic
   waveform fitting at 2 Hz. Rung 3 is therefore unverifiable in the main real domain.
3. C-IV, the primary discriminant, works at one of five stations and in a band of 1-10 Hz, whereas
   the band of the proposal is 0.1-2 Hz.

Appendix W requires a second generator and the document does not plan it in the schedule.

**Consequence.** The only branch that could test G1 on real data within the thesis time frame is
kinematics plus polarisation (C-III + C-IV). The plan must say so and be organised around it.
Section 23 says "the core is frozen", but the 48-month plan is still sized for the full programme.

---

## 4. Prioritised roadmap

The sequence follows the order in which the results can change the hypothesis. First come the cheap
synthetic controls that decide whether G1 survives in its current form (Phase 1). Next, the real data
already downloaded (Phase 2). Only then come the new acquisition and the real pilot. The identifiers
1c-1i are those of `tasks_and_real_data_plan.md`; the identifiers 1j-1n are new and are proposed here.

| Phase | Task | What it is | Effort | Data | Deciding criterion (preregister beforehand) | If it fails |
|:--|:--|:--|:--|:--|:--|:--|
| 0 | Corrections | Those of Section 6 | 2-3 days | none | - | - |
| 1 | **1j** ABL-cont | Tie control of K1 | 2-3 days | synthetic | Difference in false events ABL-cont against PHYS $<0.05$ in one-sided cells | G1 is reformulated without the one-sided cell (Section 5) |
| 1 | **1k** Empty windows | Generator with $K\in\{0,1,2,3\}$, absolute threshold calibrated on null windows, recall against false alarms per hour | 1-2 weeks | synthetic | PHYS with $\ge0.05$ fewer false alarms at non-inferior recall (margin 0.02) in $\ge$ half of the cells | The advantage was an artefact of the relative threshold; it is reported and X is rewritten |
| 1 | **X-bis** Tuned baselines | GaMMA and PyOcto with hyperparameter search on the calibration | 1 week | synthetic | The curve is reported; there is no success threshold | - |
| 1 | **1n** Calibration of $\eta$ | Bias and coverage by Monte Carlo | 1 day | synthetic | Coverage of the CI95 $\ge0.93$ | Correct the interval (BCa or profile) before 1c |
| 2 | **1f** Robustness of ZAR | Existing script, 8 subpopulations | 1 day | already downloaded | No subpopulation with the entire CI above 0.20 | Report the subpopulation; T.1 is qualified |
| 2 | **1l** Incidence at CHI | Control stratified by incidence and depth; radial S channel | 3-4 days | already downloaded | P control $\ge0.60$ at incidence $\ge30^\circ$ | The station explanation stands, now supported |
| 2 | **1m** Catalogue | Relabel with SGC; templates on clutter | 1-2 weeks | SGC and already downloaded | $\eta$ with SGC within the CI with ISC | $\eta$ depends on the catalogue; T.1 is reformulated |
| 2 | **1d** Band 0.1-2 Hz | Control and $\eta$ in the band of the proposal | 3-4 days | already downloaded | Control $\ge0.60$ at ZAR in the band | C-IV operates outside the band of the proposal; declare the change of band |
| 3 | **1c** Network map | Control at 20-40 stations; distribution of $r_s$ | 3-6 weeks | new, event windows | The distribution is reported; there is no success threshold | - |
| 3 | **1h** Per-station likelihood | Calibrated likelihood ratio (Section 5); T.2 redone with continuous $r_s$ taken from 1c, new seeds, GaMMA included | 2-3 weeks | synthetic | PHYS-LR against tuned GaMMA: $\ge0.05$ fewer false alarms at non-inferior recall with the measured $S_{\mathrm{eff}}$ | The advantage does not survive the real network composition: central negative result, publishable |
| 4 | **1g** Real pilot by decimation | See Section 5.4 | 2-3 months | new | Same criterion as 1h, on real data | G1 remains restricted to the synthetic; the thesis becomes a characterisation thesis |
| 4 | **1i** Active sequence | $\eta$ during aftershocks or in the Nest | 1 month | new | $\eta<0.20$ | Realisable impostors in the real data: Appendix T becomes central |
| 5 | Second domain and second generator | Requirement of Appendix W | E1-E4 | - | Ladder and hierarchy replicated | The ladder is a property of the generator |

Phase 1 costs less than a month, requires no data and can change the formulation of G1. Running 1c
or 1g before it means investing in data to test a hypothesis that may have to be rewritten.

---

## 5. Strategic reformulation

### 5.1 G1, rewritten against the state of the art and with the operational endpoint

> **G1'.** In continuous windows with an unknown number of events, including empty windows, a
> score that adds to the kinematic admissibility C-III a polarisation likelihood ratio calibrated per
> station reduces the false alarms per unit time by at least 0.05, at non-inferior recall (margin
> 0.02), against the best published associator tuned with the same budget. The reduction occurs in
> the regime $S_{\mathrm{eff}}\in[S_1,S_2]$ and vanishes outside it.

This form has three advantages over the current one. It compares with GaMMA and not with an in-house
ablation. It uses the endpoint that exists in real data. And it turns the result of T.1 into a
parameter of the model, instead of a caveat.

### 5.2 A single quantity for T.1 and T.2: per-station reliability as a mixture weight

Let $\delta$ be the folded angle between the measured polarisation and the predicted direction. Under
hypothesis $H_1$ ("the arrival comes from the candidate"), a station with reliability $r_s$ produces

$$p_{1,s}(\delta)=r_s\,p_{\mathrm{vM}}(\delta;\kappa_s)+(1-r_s)\,p_{0,s}(\delta),$$

and under $H_0$ ("the arrival is unrelated to the candidate") it produces the null geometric density
$p_{0,s}(\delta)$. That density is the same one whose mean gives the $f_0$ of T.1. The control of T.1
estimates $r_s$ with the same algebra as $\eta$:

$$f_{\mathrm{ctrl},s}=r_s\,q_s+(1-r_s)\,f_{0,s}\quad\Longrightarrow\quad
\hat r_s=\frac{f_{\mathrm{ctrl},s}-f_{0,s}}{q_s-f_{0,s}},\qquad q_s=\int_0^{20^\circ}p_{\mathrm{vM}}(\delta;\kappa_s)\,d\delta .$$

The polarisation contribution of a candidate $h$ is then

$$\ell_{\mathrm{pol}}(h)=\sum_{i\in\mathcal M(h)}\log\frac{p_{1,s_i}(\delta_i)}{p_{0,s_i}(\delta_i)} .$$

This form has four properties that today are obtained with separate patches:

1. With $r_s=0$ the contribution is exactly zero. The neutrality of PHYS-N follows from the
   formulation and ceases to be a rule designed after seeing PHYS-G.
2. It is a sum, so it does not reward candidates that explain less. It absorbs correction (20c).
3. It is continuous, so it breaks ties with evidence and not with artefacts, which answers K1.
4. If 1i or 1m find $\eta>0$, the clutter itself is a mixture and
   $p_{0,s}\leftarrow\eta\,p_{\mathrm{vM}}+(1-\eta)\,p_{0,s}$.

The two real measurements of the programme, $r_s$ and $\eta$, thus become the two parameters of
constraint C-IV. This gives T.1 a structural role in the thesis and not only a validation role.

### 5.3 The ladder as a function of reliability

The number of stations $S$ is replaced by $S_{\mathrm{eff}}=\sum_s \hat r_s$ in the sufficiency
ladder. Hypothesis C1 is reformulated as follows: the crossover $S^\ast$ depends on $S_{\mathrm{eff}}$
and not on $S$. It is falsifiable in synthetic data, with T.2 redone with continuous $r_s$ (1h), and in
real data, with the decimation of 5.4. If it is confirmed, 1c ceases to be an inventory and becomes a
map of where the Colombian network needs physics and where arrival times suffice. It is the product
of greatest applied impact of the thesis.

### 5.4 Design of the real pilot (1g): decimation of a dense network

A zone of the RSNC with dense coverage and an SGC catalogue complete above $M_c$ is chosen. The
catalogue of the full network serves as ground truth. Subnetworks of 3, 4 and 6 stations, including
one-sided subnetworks, are built by station selection. All methods (GaMMA, PyOcto, ABL-cont,
C-III + C-IV-LR) receive the same PhaseNet picks. The endpoint is the curve of recall against false
alarms per day for events with $M\ge M_c$. A detection without a match is reviewed by template
correlation against the full network before being counted as false, as required by the governance
of Section 16. The design requires no velocity model for waveforms (it uses a 1-D model for times
and the bound $V_{\min}$ for C-III). It reproduces in real data the same cells as the synthetic
(ladder and one-sided network), so it is the direct transfer test that is missing today. The limit
must be declared: it tests the kinematic and polarisation prototype, not the full waveform model.

### 5.5 Scope reduction

Section 23 defends three claims. The learned residual, the SO(3) equivariance and the chain L0-L3
are not necessary for any of them. It is recommended to move them to an annex of "conditional
extensions" with zero budget until G1' has been decided. The 600-1200 GPU-h freed can cover the
second generator and the second domain that W requires. Optimal transport is retained: it is
necessary for the waveform branch (F2) if that branch survives in the external domain.

---

## 6. Writing and internal-coherence observations

Each point is verifiable in the text and its correction does not change any result.

| Location | Problem | Correction |
|:--|:--|:--|
| L36, L3908 against L1906; `executive_summary` L193 and L215, `task1_step_by_step` L41, `tasks_and_real_data_plan` L3 and L828 | "sixteen" against "seventeen" in the proposal; "diecisiete" and "quince" in the Spanish documents | Count appendices G-V and fix a single figure in all five documents |
| L521-533 against L550, L564, L570, L590, L652, L683, L710, L736, L809, L3053, L3059, L3125 | The taxonomy admits five labels, but the text uses "met", "met decisively", "supported" and "confirmed" | Replace with Numerical evidence or Preliminary |
| L570 (G1) and L550 (A2) | G1 "supported on both endpoints" without the qualifier from T.1, T.2 or the parity with GaMMA; A2 does not reflect the per-station dependence | "Preliminary: against ABL in one-sided cells; not against GaMMA (X.2); conditional on per-station reliability (T.2)" |
| L1858-1860 | Claim (ii) contradicted by X.2 (K3) | Rewrite with the data of X.2 and X.4 |
| L1896 against L3958 | "(20c) ... has not been run", but protocol P2 used it | "Used in P2; its effect was secondary (X.1)" |
| L2035 | "none of which has yet been executed", when GaMMA and PyOcto were in fact executed | "GaMMA and PyOcto executed on synthetic data (X); the others pending" |
| L3972-3973 | "the two seeds agree cell by cell": false for PyOcto (1.00 against 0.10; 0.97 against 0.00) | Restrict the sentence to ABL and PHYS |
| L4034-4036 | "the single cell where the two seeds disagree" in the curve metric: false, because the ladder cells also disagree (3 stations: not reached against 0.900; 4: not reached against 0.533; 6: 0.167 against 0.033) | Report the discrepancy in all cells |
| L1086 | "Status: proved as above": the text derives the leading term but does not bound the remainder $O(V/(f_{\max}d_{\min}))$ | "Proposition, leading-order derivation; remainder bound part of T1", or add the bound in an annex |
| L1057-1063 | Broken sentence ("the crude but explicit bound The following separable lower bound...") and a stray fragment ("precisely the number of captured arrivals...") | Rewrite the paragraph of (8c) |
| L1093 against L2453 and L2479 | B3 "$\lvert\rho\rvert\ge0.6$ met" against the registered Spearman threshold $\ge0.5$ | Use 0.5 and cite which test (App. I or App. S) meets it |
| L2229 against L2165 and Section 6.1 | "three corollaries" against "Corollaries 1-4" (the text has four) | "four corollaries" |
| L2058 | "remove the analytic factorization (18)": the factorisation is (22); (18) is the information decomposition | Change to (22) |
| L2236 | "Rank result (17)": (17) is the Lagrangian; the rank result is (19)-(20) | Change to (19)-(20) |
| Sections 12-14 | Numbering (23), then (20a)-(20c), then (25); (24) is missing | Renumber in order of appearance |
| L2300 | "The core is frozen. The exploratory programme is closed", with T.2 added and tasks 1c-1i registered | "The confirmatory core is frozen; real-data calibration of C-IV is open (T.1, T.2)" |
| Summary (about 520 words) | Three selective claims: the three-station recall without mentioning false events of 60-100%; the elimination of false events in one-sided networks without mentioning that GaMMA also obtains zero; "the synthetic advantage survives the measured network composition", when one of four reliable stations was tested, with PHYS-N designed a posteriori and one replicate seed (*Unproven speculative assertiveness*) | Reduce to about 300 words and report recall and false events together |
| T.1, "η recovers planted shares" | One replicate per proportion; 0.098 for 0 and 0.223 for 0.3 | "within its interval in single replicates; bias and coverage under calibration (1n)" |
| Table 8, months 1-6 | Mixes progress already made ("already established") with plan | Move the progress to Section 4 and leave in the table only the pending work |
| Title | Title of 9 words plus subtitle of 12 | Acceptable for a proposal; for papers, keep only the title |

---

## 7. Final recommendation

**Major revisions. The proposal is defensible after Phase 0 and Phase 1, not before.**

The programme has a solid and original core: physics as a constraint that rejects false events, the
sufficiency ladder and an observability certificate usable for network design. It also has an
uncommon culture of evidence. The central hypothesis, in its current wording, is supported by a
comparison against an in-house ablation, with a threshold that does not exist in continuous data and
with possible ties in the score of the comparator. Against GaMMA, in the cells where it is declared,
the difference is nil.

The strategic recommendation is not to produce more pre-studies of the G-W type. The three synthetic
controls that decide the form of G1 (1j, 1k, X-bis) must be run first. Next, reformulate C-IV as a
per-station likelihood ratio, so that T.1 ceases to be a caveat and becomes the procedure that
estimates the parameters of a constraint. Finally, take G1' to real data through the decimation of a
dense network. If G1' survives that sequence, the thesis will have a claim proven against the state
of the art and on real data. If it does not survive, the characterisation of where the Colombian
network needs physics and where it does not remains a complete and publishable thesis, and the
preceding plan produces it in both cases.
