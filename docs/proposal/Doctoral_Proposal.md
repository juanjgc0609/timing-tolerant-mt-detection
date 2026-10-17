---
title: "Timing-Tolerant, Bounded-Influence Moment-Tensor Detection of Overlapping Earthquakes: A Preregistered Benchmark and a Real-Data Test Plan"
subtitle: "Event Recall at a Fixed False-Event Rate in Sparse Seismic Networks, with Constraints Generated from Elastodynamics"
author: "Juan Jose Gordillo Cordoba"
date: "October 2026"
---

**DOCTORAL RESEARCH PROPOSAL**

# Timing-Tolerant, Bounded-Influence Moment-Tensor Detection of Overlapping Earthquakes: A Preregistered Benchmark and a Real-Data Test Plan
## Event Recall at a Fixed False-Event Rate in Sparse Seismic Networks

*Timing-tolerant moment-tensor detectors with deflation, a moment-tensor subspace statistic that
bounds the contribution of each station, their preregistered synthetic benchmark against published
associators and migration locators, and a test on a decimated national network*

*The title names detection rather than location because the measured benefit of physical evidence is
rejection of false events, not spatial resolution (Appendices H and N). It names bounded influence
because that is the element of the method whose effect two preregistered tests confirmed, in
dense sequences and with a modest size (Appendices Y.24, Y.26), the real-data test plan because all
evidence to date is synthetic, and the benchmark because its negative results
are part of the contribution: the group BLASSO with station time shifts solved off the grid, registered
as the methodological hypothesis G1-M, did not beat block pursuit with equal tolerance (Appendix Y.22)
and remains the framework of the theory rather than the detector. The primary hypothesis, G1, the
refuted hypothesis G1-M, the bounded-influence hypothesis G1-R and the calibration hypothesis G1-C,
which one further preregistered test confirmed within a narrow scope (Appendix Y.28), are stated in
Section 4.3; the design
of G1 is in Section 13.0c.*

**Student**

Juan Jose Gordillo Cordoba

**Proposed Supervisor**

Uram Anibal Sosa Aguirre

**Research areas**

Applied Mathematics; Constrained Optimization; Inverse Problems; Numerical PDEs;
Computational Science; Computational Seismology; Machine Learning

**Document conventions.** A4 page size, IEEE-aligned numeric citation conventions. The governing
physical model is three-dimensional heterogeneous linear elastodynamics. Learned components are
computational representations of the physical state and never substitute for physical inference.
The target operating regime is declared in Section 1.5, which also states that the synthetic evidence
gathered so far lies outside it; quantitative margins are registered in Appendix D and the
computational budget is itemized in Section 12. The programme is informed by
sixteen numerical pre-studies reported in Appendices G through V, a cross-study statistical audit in
Appendix W, a comparison with published associators in Appendix X, controls for the catalogue
endpoint in Appendix Y, and a first real-data measurement in Appendix T.1. Negative results are retained because they delimit the scientific claims and the scope of
the confirmatory programme. Appendices A to F are in this document; Appendices G to AA, which hold the
pre-studies, the preregistered benchmarks and the theory, are in the supplementary evidence volume
`Supplementary_Evidence_Appendices_G_to_AA.md`, with the same lettering and numbering.

---

# Abstract

This thesis addresses detection of overlapping earthquakes in sparse networks where picking discards
kinematic information. Its primary hypothesis, G1, is that a timing-tolerant moment-tensor detector with
deflation, applied to three-component waveforms before picking, detects more events than published
associators and migration locators at equal realized false-event rates.

On synthetics, in the first of eight preregistered tasks (2A), such detectors exceeded GaMMA, PyOcto and
QuakeMigrate by 0.09 to 0.49 in every cell, with operating points chosen on evaluation windows; independent
verifications exist for Tasks 2A, 2E, 2F, 2I and 2L, and for the real-noise Task 2M. Convex formulations
did not beat block orthogonal matching pursuit with the same peak search (G1-M refuted). A
bounded-influence statistic did in dense sequences, in two verified tests (G1-R: 0.04 to 0.07 with a
per-station cap, 0.058 and 0.053 with a noise-scaled Huber fit), and station terms estimated without
picks raised origin-time-matched recall by 0.057 and 0.062 (G1-C), mainly in one-sided eight-station
networks. With synthetic events injected into real noise of days no earlier analysis used, a registered
test confirmed the bounded-influence ordering (G1-R, +0.228 and +0.241 over the uncapped statistic on a day
with transients, +0.057 and +0.120 on near-Gaussian days); the gain of the cap is concentrated in the hours
with transients.

The first registered real-data test, one hour of the 2019 Ridgecrest aftershock sequence on four
subnetworks, was negative. At 30 false detections per hour with origin times
matched, the thesis detector recalled at most 0.022 of the 45 catalogued events of magnitude 2 or more,
against up to 0.667 for PhaseNet with PyOcto (pooled difference -0.383 [-0.481, -0.310]). Pick-free
detections were mislocated by 8 to 18 km in median, because a homogeneous medium with constant station terms left
timing errors of 0.35 to 0.41 s that its tolerance did not absorb, and its scores did not separate real
events in continuous seismicity. A follow-up exploration found that a layered or fitted medium does not make the arrival readings more
predictable, and that a logarithmic score restores the ranking of real events while keeping most of the robustness
to transients on one exploratory seed, but does not locate them.

The contribution is a preregistered benchmark with its positive and negative results and the measured
conditions under which pick-free physics helps. Transfer to real data requires a reading of arrivals that
survives dense seismicity, the remaining work before the test on the decimated Colombian network.

---

## Contents

- [1. Problem, Finding, and Thesis](#1-problem-finding-and-thesis)
- [2. Scope and Significance](#2-scope-and-significance)
- [3. State of the Art and Research Gap](#3-state-of-the-art-and-research-gap)
- [4. Questions, Hypotheses, Objectives](#4-questions-hypotheses-objectives)
- [5. Physical Foundations](#5-physical-foundations)
- [6. Five Sources of Physical Constraint](#6-five-sources-of-physical-constraint)
- [7. The Inverse Problem](#7-the-inverse-problem)
- [8. Constrained Optimization and Sensitivities](#8-constrained-optimization-and-sensitivities)
- [9. Identifiability](#9-identifiability)
- [10. Differentiable Wavefield Representation](#10-differentiable-wavefield-representation)
- [11. Computational Architecture](#11-computational-architecture)
- [12. Computational Budget](#12-computational-budget)
- [13. Experimental Design](#13-experimental-design)
- [14. Metrics and Statistical Protocol](#14-metrics-and-statistical-protocol)
- [15. Baselines and Ablations](#15-baselines-and-ablations)
- [16. Test Domains, Data, and Governance](#16-test-domains-data-and-governance)
- [17. Risks and Go/No-Go Criteria](#17-risks-and-gono-go-criteria)
- [18. Reproducibility and Software](#18-reproducibility-and-software)
- [19. Forty-Eight-Month Two-Track Plan](#19-forty-eight-month-two-track-plan)
- [20. Contributions and Publication Strategy](#20-contributions-and-publication-strategy)
- [21. Alignment with the Supervisor's Research Profile](#21-alignment-with-the-supervisors-research-profile)
- [22. Scientific Integrity](#22-scientific-integrity)
- [23. Conclusion](#23-conclusion)
- [Appendix A. Notation](#appendix-a-notation)
- [Appendix B. Acceptance Tests](#appendix-b-acceptance-tests)
- [Appendix C. Experimental Matrix](#appendix-c-experimental-matrix)
- [Appendix D. Pre-Registered Thresholds](#appendix-d-pre-registered-thresholds)
- [Appendix E. Literature-Search Protocol](#appendix-e-literature-search-protocol)
- [Appendix F. References](#appendix-f-references)
- [Appendices G to AA. Supplementary Evidence Volume](#appendices-g-to-aa-supplementary-evidence-volume)

---

# 1. Problem, Finding, and Thesis

## 1.1 Problem statement

Given three-component recordings from a seismic network over an association window, determine how
many events occurred, and for each one its epicenter, depth, origin time and approximate source
mechanism, in the presence of structured false arrivals, heterogeneous and imperfectly known media,
irregular station geometry and overlapping signals. The unknown is a set of unknown cardinality,
so the problem is simultaneously one of model selection, continuous parameter estimation and
constrained optimization.

The standard pipeline reduces each record to picked arrival times and then clusters those times.
That reduction is the origin of most of the difficulty. An arrival time is a scalar coincidence:
any impulsive noise burst within the tolerance window is admissible evidence for any hypothesis
that predicts an arrival there. The association literature reports the consequence directly, that
high event rates and overlapping events remain difficult and that in large-aperture networks
travel-time information alone cannot discriminate spurious alignments of real moveout with false
picks [6]. That published statement motivates the present work; it is not evidence produced by it.
The evidence produced here is stated in Section 1.6 and labelled by the taxonomy of Section 4.3.

## 1.2 The finding that determines the approach

Rather than assume that using the full wavefield helps, and assume where it helps, two numerical
pre-studies were run before the research programme was finalized. Both compare, under matched geometry, matched
band and equal signal energy, a full three-component physical prediction (polarization, geometrical
spreading, P and S branches) against an arrival-time-only prediction carrying the same $t_P,t_S$
with no polarization and no amplitude. Details and code are in Appendix H.

**Resolution between two nearby true sources improves only marginally.** The separation at which
two hypotheses become distinguishable improves by six to fourteen per cent, and the improvement is
flat from one station to twelve. Physical evidence does not sharpen the resolution limit in any
practically decisive way.

**Rejection of false hypotheses improves by a factor of four.** With one true source embedded in an
increasing number of randomly polarized false arrivals, the discrimination margin between the true
hypothesis and the best false hypothesis more than five kilometres away fell from $3.93\sigma$ to
$1.33\sigma$ for the arrival-time representation, while the physical representation held between
$16.4\sigma$ and $11.6\sigma$. At the highest clutter level a false hypothesis reached $83\%$ of
the true arrival-time score but only $30\%$ of the true physical score. That is a statement about
the **margin**. Whether a larger margin becomes fewer false events in an emitted catalogue, against
the best kinematic score rather than a coincidence count, is a separate question; Appendix Y finds
that the answer depends on how much the picker discards, and Section 1.6 states the evidence.

The mechanism is not subtle. Clutter satisfies a time window because a window is one scalar test.
Clutter cannot fake a joint pattern of polarization directions, amplitude decay with distance and
two branches at the correct speed ratio, because those are coupled constraints across components
and across stations.

## 1.3 The thesis, in one sentence

> **Overlapping earthquakes in sparse networks are detected, before picking, by a greedy pursuit over
> moment-tensor subspaces that tolerates what the model gets wrong in timing at each station and bounds
> the influence of what it gets wrong in amplitude; constrained optimization with constraints generated
> from elastodynamics is the framework in which that detector is derived and analysed, not the method
> that passed the tests.**

The first clause states what passed synthetic tests and the second names the framework. The framework
reads each additional piece of physics as an additional constraint on the feasible set and clutter as
what falls outside it, and what the pre-study measured under the declared conditions is its margin
form. At the catalogue endpoint, what has passed synthetic
tests is timing-tolerant moment-tensor detection with deflation, of which the greedy detector MT-FW and
block orthogonal matching pursuit [129] given the same peak search are equivalent members in standard
windows, not the certified convex solution itself. On the grid, with or without that peak search, and
off the grid with a time shift per station and phase entered explicitly as a nuisance parameter, the
convex problem did not beat block pursuit with equal tolerance, which refutes the methodological
hypothesis G1-M in this synthetic design (Section 1.6, Appendices Y.20 to Y.22). The optimization
statement therefore supplies the framework, the dual certificate of the convex problem and the theory.
The element of the deployed statistic that preregistered tests confirmed is a robust one: bounding the influence of any single
station and phase block whose amplitude the model gets wrong, which in dense sequences of overlapping
events raises recall over block pursuit with equal timing tolerance (hypothesis G1-R). Two forms were
tested, each on its own pair of fresh seeds: clipping the contribution of each block at a calibrated
level (Task 2F, Appendix Y.24) and a Huber fit whose level is fixed by rule from the noise distribution
(Task 2I, Appendix Y.26). Modelling those amplitudes as station terms inside the statistic did not
replace the clipping in an exploration (Appendix Y.25). A second element was confirmed in one further
preregistered test, and it concerns calibration and not the statistic: station terms estimated without
picks from the catalogued events raise the origin-time-matched recall of the detector with a local
timing refinement, a gain concentrated in one network geometry (hypothesis G1-C, Task 2L, Appendix
Y.28).
Whether the result holds on real data is the question the confirmatory programme is organized around.
The thesis statement also fixes the priority order of the whole programme: the one primary endpoint
is event recall at a fixed realized false-event rate (hypothesis G1, Section 13.0b), and localization
accuracy is secondary.

## 1.4 Why this belongs to applied mathematics

The research object is a constrained, non-convex inverse problem with continuous and discrete
latent variables, a PDE-governed forward operator and inequality constraints of physical origin.
Four specific mathematical challenges are named rather than asserted generically.

1. **Constraint generation from a PDE.** Turning statements about wave propagation - observability,
   kinematic admissibility, causality, and non-detection under a declared sensing model, which unique
   continuation motivates but does not imply - into a feasible set that a primal-dual solver can
   handle, and characterizing that set.
2. **Well-posedness of the reduced problem.** Establishing when the reduced Hessian is positive
   definite, which is where the geometric observability condition enters as a coercivity statement.
3. **Non-convexity by phase sensitivity.** A time-domain waveform misfit is non-convex in source
   position and origin time, the cycle-skipping obstruction of full-waveform inversion [29], [30].
4. **Sensitivity reliability.** Gradients driving a constrained solver must be validated against an
   independent adjoint or analytic reference, never accepted because automatic differentiation
   returned a number [76].

Machine learning contributes differentiability, amortization and continuity in source coordinates.
It does not define the scientific question.

## 1.5 Declared operating regime

| Quantity | Primary regime | Stress-test regime |
|:--|:--|:--|
| Model domain | $120\times120\times40$ km | reduced 2-D / 2.5-D control |
| Frequency band | $0.1-2$ Hz | to $5$ Hz, reduced dimension only |
| Minimum shear speed | $2.0$ km/s | $1.2$ km/s |
| Window duration | $\le120$ s | $\le120$ s |
| Stations | $1-40$ three-component | $1-60$ |
| Events per window | $1-8$ | to $30$ |
| False-arrival rate | $0-60$ per window | to $200$ |
| Approximate magnitudes | $M\gtrsim2.5$ regional | $M\gtrsim1.5$ local dense |

Explicit-solver cost scales as $f_{\max}^4$ (Section 12), so $f_{\max}$ decides executability and
is declared in advance. The band limit restricts the accessible magnitude range and is therefore a
scientific limitation stated here rather than discovered later.

The evidence gathered so far does not lie inside this regime, and no result in this document shows
that it transfers to it. The synthetic tests on which G1 rests (Appendices Y.10 to Y.21) use Ricker
wavelets centred at 2 Hz, or at 1.4 to 2.8 Hz from Y.14 on, data sampled at 16 Hz, sources in a $30\times30$ km
epicentral square at 5 to 13 km depth, a homogeneous or one-dimensional gradient medium and networks
of four to eight stations. The band overlaps the primary band only at its upper edge and the source domain is a
quarter of its width. What holds at 0.1 to 2 Hz over $120\times120\times40$ km is therefore open, and
the real-data decimation test of Section 16 runs at the band that the RSNC data and the selected
velocity model support, determined per domain by the velocity-model gate.

## 1.6 The evidence to date, and what the thesis defends

The problem is posed as constrained inversion, from the three-component wavefield, for a sparse
source measure of unknown cardinality, with constraints generated from elastodynamics (Section 6).

Numerical pre-studies on a synthetic generator found that physical evidence barely sharpens
resolution between nearby sources (six to fourteen per cent) but enlarges fourfold the margin by
which clutter-supported hypotheses are rejected. Controls then tested whether that margin becomes
fewer false events in an emitted catalogue. The earlier advantage over a coincidence-count ablation
was a tie effect of that ablation's discrete score. In a generator that picks every arrival exactly,
a continuous kinematic score matched the polarization prototype; but that generator gives timing
information no picker has. In a physically consistent generator, where detection, pick timing and
polarization degrade together with the signal-to-noise ratio, more than half of the arrivals are
never picked, sources have random double-couple mechanisms, clutter follows the same detection rule
as signal, and station delays, polarization deviations and site amplitudes are unknown to the
method, the detector that projects the three-component records of all stations onto the
six-dimensional moment-tensor subspace and explains detected events away by a greedy matching
pursuit whose selection step is the Frank-Wolfe oracle exceeded seeded, tuned GaMMA by 0.11 to 0.28
in event recall at 0.25 realized false events per window with weak clutter and by 0.17 to 0.43 with
clutter made of visible transients, in two seeds on fresh windows, while the same iteration with the
geometry of another hypothesis lost 0.13 to 0.41. A preregistered confirmatory benchmark then removed
the remaining shortcuts together: a medium, wavelet, coda and fixed station delays the detector does
not know, an operational STA/LTA picker, station corrections estimated from a catalogue with location
errors, and QuakeMigrate among the comparators. The detector with corrections exceeded the best of
GaMMA, PyOcto and QuakeMigrate by 0.09 to 0.49 in event recall in every cell of both seeds and both
clutter regimes. QuakeMigrate was the effective competitor, because the picker starves the
pick-based associators; the intervals are per window, over five networks per cell; and the sources
remain double couples inside the detector's hypothesis class. The operating point of every method was
chosen on the evaluation windows; the deployable rule, with thresholds frozen on other networks,
fails under its strict reading because the frozen threshold realized up to 2.3 times its budget in
one cell; and with optimistic synthetic picks and weak clutter GaMMA was ahead of the detector at 0.1
false events per window in four of the six cells, so that in that arm G1, which requires a gain of
0.05 at both 0.1 and 0.25, is not met. A preregistered ablation on the same windows, which
replaces the radiation operator by the identity inside the same solver, showed that the operator adds
0.06 to 0.35 in recall with strong clutter, where without it the detector only matched QuakeMigrate,
but with weak clutter its rule held only at the minimum of two of three cells, and there much of the
advantage over QuakeMigrate did not need the operator; because the identity removes sign,
polarization constraint and amplitude geometry together, it does not isolate the sign of the
radiation pattern. A third preregistered test ran the closest antecedents on the same windows and
replicated the comparison with QuakeMigrate over thirty networks per cell. In the forms implemented,
which read the matched filter at the exact predicted times of a grid of 3 km and 0.1 s, block
orthogonal matching pursuit [129] and mechanism-corrected stacking [122] fell below the detector by
0.08 to 0.52 at 0.25 false events per window in every cell, although in the two strong-clutter
four-station cells the advantage over block matching pursuit is not significant after Holm's
correction; against QuakeMigrate, hierarchical
intervals over thirty networks excluded zero in all twelve cells, and after Holm's correction every
cell but the strong-clutter four-station cell of one seed remained significant; and the group BLASSO
on the same grid, solved with certified duality gaps, fell 0.06 to 0.61 below the detector in 23 of
24 cases and never exceeded it. An exploratory diagnosis found that with sources on the grid nodes
block matching pursuit is the best of these methods (0.97 against 0.76 at 0.25 false events per
window with four stations), which pointed to the peak search of the detector, a search within
$\pm0.25$ s of the predicted time at each station and phase, as what absorbs grid and timing errors. A
fourth preregistered test, prompted by a second independent audit, gave the antecedents that same peak
search and ablated the peak search and the per-station cap of the detector on the same windows
(Appendix Y.21). Its central result is negative for the distinctiveness of the detector. With equal
tolerance, block orthogonal matching pursuit matched it: the differences at 0.25 false events per
window ran from -0.155 to +0.091 and every hierarchical interval contained zero, so the registered
rule failed in both clutter regimes and, as registered, the advantage of the detector over [129] and
[122] is attributed to the peak search alone. Removing the peak search from the detector cost 0.06 to
0.41 in all twelve cases, removing the cap changed recall by amounts whose intervals all contained
zero, the group BLASSO on the grid given the same peak search stayed 0.11 to 0.30 below the detector in
11 of 12 cases, and block matching pursuit with the peak search exceeded QuakeMigrate by 0.17 to 0.49.
The peak search is a trade-off rather than a free gain: with exact timing and sources on the nodes it
cost 0.10 to 0.14 on six and eight stations (exploratory). The same test froze split-conformal
thresholds on noise-only windows of other networks before evaluation; its rule held cell by cell, but
those thresholds control false alarms only in windows without events, and in weak-clutter windows
containing events the detector emitted 1.2 to 2.2 false events per window, and over all windows three
to eight times as many false events as QuakeMigrate, so they do not give a deployable operating point
for a false-event budget. With five networks per cell none of its three rules survives Holm's correction.
A fifth preregistered test, whose registry hash was communicated outside the computing environment
before any of its code was written, ran the principled form of that tolerance: the group BLASSO with a
time shift per station and phase, solved off the grid by sliding Frank-Wolfe, against block orthogonal
matching pursuit with the peak search, on thirty networks per cell (Appendix Y.22). It gained 0.000 to
0.039 in recall at 0.25 false events per window, every hierarchical interval contained zero, no cell
passed its rule in either seed or regime, and the group penalty was inert at its calibrated value; as
registered, the methodological hypothesis G1-M is refuted in this synthetic design. An independent
verification found no implementation error, since the continuous atom reproduces the block of the grid
method to machine precision, and no budget handicap, since ten times the iterations with a nine-start
refinement changed recall by -0.014. It also found that the benchmark is not saturated: oracle ceilings
that keep the false events caused by model mismatch lie 0.08 to 0.24 above block pursuit at 0.25, and the
missed events are absorbed by atoms mislocated 5 to 10 km away, beyond the reach of local sliding. An
exploration on the calibration seed then examined seven candidate sources of gain (Appendix Y.23):
coarser grids, larger timing errors, a physical one-dimensional velocity model inside the detectors,
oracle or calibrated, and a source-pulse-invariant subspace gave none, and an attribution run showed
that the headroom persists when the medium of the generator is exactly the homogeneous medium of the
detectors; the one candidate that met the exploratory criterion was the per-station cap of MT-FW in
windows with three to five overlapping events. Selected among those candidates, it was registered as a
sixth test on two seeds never used before, with its hash communicated before the run and its code frozen
by hashes (Appendix Y.24). In dense windows MT-FW with the cap exceeded block pursuit by 0.072 and 0.074
with weak clutter and by 0.051 and 0.043 with strong clutter at 0.25 false events per window, with 96
networks per cell, intervals excluding zero and all four tests surviving Holm's correction; without the
cap MT-FW equalled block pursuit within 0.002, so the cap is the cause; and in standard windows the two
were equivalent within the margin of 0.05, narrowly in one seed, with negative point estimates in strong
clutter that leave a cost of about 0.03 on isolated events unexcluded. An adversarial verification
confirmed the result and narrowed it. With the origin time also required to match within 0.5 s the gain
is 0.024 to 0.045. It belongs to the cap level fixed in an earlier task on other data, where it peaks,
and it reverses at low levels and vanishes at high ones. Its mechanism is bounded influence: the cap
limits what any single station and phase block contributes when the radiation model gets its amplitude
wrong. The amplitude factors of the generator are redrawn in every window, so they are perturbations of
each window and not station terms; removing them lowers the gain to 0.023 with an interval containing
zero, yet at a small spread of those factors the cap still gains where a statistic given the true
factors does not (Appendix Y.25, exploratory). The gain also comes from a better calibration of the
score across windows, not from better separation of events within a window. Even with the cap, recall at 0.25 false
events per window in dense windows is only 0.22 to 0.29.
A seventh preregistered test replaced the calibrated cap by a level that no detection result selects
(Appendix Y.26). It came from an exploration on the calibration seed in which this variant was the best
of 11 candidates, a selection that its registration declares and for which it planned on a smaller
effect. The statistic applies the Huber function to the norm of the residual of each station and phase
3-vector, with the level fixed by rule at the square root of the 0.95 quantile of a chi-square with
three degrees of freedom, 2.7955 noise standard deviations, four reweighting passes from the
least-squares fit, and deflation with the robust fit. The registry hash was sent to the user before any
run and the code was frozen by hashes; the test used two further seeds never used before, with 120
networks per cell and regime and five dense and five standard windows per network; and its run was
interrupted twice and resumed from its stored parts with unchanged code. Its three rules passed in both
seeds. In dense windows the Huber statistic exceeded block pursuit by 0.058 and 0.053 at 0.25 false
events per window, with intervals from 0.039 to 0.076 and from 0.037 to 0.070; it was not inferior to
MT-FW with the cap, the differences being -0.008 and -0.002 with lower bounds of -0.023 and -0.018
against a margin of -0.03; and in standard windows, pooled over regimes, it was equivalent to block
pursuit within 0.05, with differences of +0.010 and +0.007. On the same windows the cap repeated the
gain of the sixth test, 0.066 and 0.055. An independent verification reproduced the provenance, the
symmetry of the comparison and the statistics, regenerated six parts from the three sessions of the
interrupted run identically, and narrowed the scope. With the origin time required to match within 0.5 s
the gain is 0.031 and 0.019, and the point criterion of 0.03 would fail in one seed. The equivalence on
isolated events holds pooled and not regime by regime: with weak clutter in one seed the Huber
statistic is 0.032 below block pursuit, with an interval excluding zero. Recall with every atom
accepted does not rise, so the gain is a better ranking and selection of atoms, about half of it the
reliability of the score across windows and half the choice of atoms, and not more events found. At 0.1
false events per window the Huber statistic is slightly below the cap in one seed, by 0.018. It costs
about nine times as much as MT-FW with the cap. Absolute recall in dense windows stays low, 0.235 and
0.251 at 0.25 false events per window. The secondary arms of that test, on 20 networks per cell and
regime, are descriptive and were not independently verified. Huber levels of 2.5 and 3.4 noise standard
deviations are indistinguishable from the registered 2.7955, within 0.006 in dense windows, and a level
of 5 is weaker. At an amplitude spread of 0.25 the cap gives nothing over block pursuit, 0.005 and 0.009,
and the Huber fit little, 0.036 and 0.002. At a spread of 0.8 the Huber fit beats block pursuit by 0.100
and 0.072 but falls behind the cap by 0.055 and 0.069, with an interval that excludes zero in one seed,
so its non-inferiority to the cap does not extend to that spread.
Four earlier positives did not survive auditing: a tie effect, a reading at a
nominal budget, a per-pick polarization likelihood that gained only while synthetic S waves had the
polarization it predicts, and margins over GaMMA inflated by clutter picks no picker would produce.
The method is fragile to unmodelled station delays beyond its search window; corrections estimated
from twenty catalogued events repair that when clutter is weak, but not when strong clutter
contaminates the catalogue. On the Colombian network, single-station P polarization passed a
reliability gate at one of five stations, and there real false picks were not coherent with the
preceding catalogued event.

The thesis defends one primary hypothesis inside four contributions (Section 20.0) and reports the
outcome of two methodological hypotheses. The primary hypothesis, G1, concerns the family of
timing-tolerant moment-tensor detectors with deflation, which includes the greedy MT-FW and block
orthogonal matching pursuit [129] given the same tolerance (BOMP-PK): that such a detector, applied to
the three-component waveform before picking and thresholded at a false-alarm rate calibrated on
noise-only windows, detects more events than the published associators and migration locators at the
same realized false-event rate. The empirical advantage over pick-based associators and QuakeMigrate on
these synthetics belongs to that family, and the ingredient that matters for it is tolerance to timing
error from grid discretization, residual station delays and velocity mismatch, supplied by a heuristic
peak search within $\pm0.25$ s. In standard windows MT-FW is not a better detector than [129] with equal
tolerance; in dense sequences its per-station cap makes it one. For that family G1 has numerical support
on synthetics, with geometry-breaking controls and in preregistered benchmarks, against the published
associators, QuakeMigrate and the exact-time antecedents run with the operational picker (Appendices
Y.18 to Y.21); as formulated it is not met at 0.1 false events per window against GaMMA fed optimistic
picks with weak clutter. The group BLASSO over moment-tensor-valued measures is the framework in which
the detection statistic is an optimality certificate, and the object of the theory, but not the best
detector: on a coarse grid its certified solution falls below the family with or without the peak
search, and off the grid with a time shift per station and phase it does not beat BOMP-PK at equal
tolerance, so the methodological hypothesis G1-M is refuted in this design (Appendices Y.20 to Y.22);
and the statistic of the family, with its peak search and cap, is not that certificate (Section 7.5,
Appendix AA.8). Station time shifts estimated jointly with the moment tensor of one known event are
prior art [143], [144], and their use as nuisance parameters of a detection certificate, the element
that G1-M added, is reported as a negative result. The methodological hypothesis that preregistered
tests did confirm, G1-R, is narrower: a bounded-influence moment-tensor subspace statistic, which
limits what the peak vector of any single station and phase block can contribute to the projection on
the moment-tensor subspace, improves event recall at a fixed false-event rate in dense sequences of
overlapping events over block pursuit with equal timing tolerance, without harming isolated events when
regimes are pooled (Appendices Y.24 to Y.26). It is supported by two preregistered synthetic tests, each
on its own fresh seeds, each with the hash of its registration held outside the computing environment
before the run, and each independently verified: Task 2F with a cap at a calibrated level, and Task 2I
with a Huber level fixed by rule from the noise distribution, which replaces the calibrated cap by a
level fixed by rule from the noise distribution in the registered primary arm; its descriptive arms find
the level immaterial between 2.5 and 3.4 noise standard deviations and the Huber fit behind the cap at an
amplitude spread of 0.8. The preregistration records of all eight tasks are local: a registry file with a
chain of SHA-256 hashes kept in the computing environment of the runs, the hashes of Tasks 2E, 2F, 2I and 2L
having also been sent in the working chat before those runs. They are not public and not sealed by a
third party, so they show consistency of the record rather than priority, and a public repository with
a sealed registry is a deliverable of the first semester.
The methodological contribution is stated accordingly, precisely and modestly: that measured design
result, with its mechanism, together with the preregistered benchmark and its negative results, G1-M
among them. Its limits belong to the statement. The evidence is synthetic and comes from one generator
family, of which the seeds are independent draws. The effect is modest, 0.04 to 0.07 under the
registered hit rule and 0.02 to 0.045 when origin times must match. The Huber form costs about nine
times the capped one, and it ranks and selects atoms better without finding more events. Clipping
per-station contributions is known practice in robust stacking and Huber fitting is standard robust
statistics, so what is claimed is the measured and explained effect of bounded influence inside a
moment-tensor detector with deflation, not an algorithm. The exploration that preceded Task 2I does not
support modelling amplitude or site terms per station inside the statistic: the generator contains no
persistent station term, and where one was added neither catalogue estimates nor the true factors took
the place of the cap (Appendix Y.25). The measured headroom, 0.08 to 0.24 below a realistic oracle
ceiling, is not explained by the velocity model or the source pulse, about three quarters of it remains
unexplained by any treatment of amplitudes, and it remains the open problem. Two exploratory studies on
the calibration seed, with oracles and with post hoc analyses declared as such, narrowed what that
problem is (Appendix Y.27). Evidence that does not use travel times does not locate in this benchmark.
By a Cramer-Rao bound the arrival times fix the epicentre to 0.4 to 0.7 km, and polarization, amplitudes
and the S to P ratio together to 4.4 to 5.6 km; a joint estimator is bounded at an improvement of 4 to 5
per cent; relocating detections by such evidence loses 0.11 to 0.48 of recall, even with the true
amplitude factors; and the limit is the waveform noise, the P arrival being 1.5 to 1.8 noise standard
deviations, which is 32 to 39 degrees of angular noise. Weighting by polarization admissibility gave a
small gain, the best of about 50 candidates and chosen post hoc, that vanishes at a polarization scatter
of 20 degrees. A Gaussian timing prior inside the detector did not help either, because station
corrections from twenty catalogued events are wrong by 0.12 to 0.24 s and a cell of the 3 km grid adds
0.18 s for P and 0.30 s for S, so that the timing budget of a true hypothesis, 0.2 to 0.4 s, is as wide
as the tolerance. With the true station delays the picture changes. They alone add 0.047 in dense and
0.098 in standard windows, and with a timing prior at the level of the model error 0.081 and 0.146,
which for isolated events reaches the realistic ceiling, 0.592 against 0.603; and with timing made
perfect 37 to 39 per cent of the isolated events still lose to a false node at the first iteration, for
reasons of amplitude. A local refinement of each detection by the timing likelihood improved location
and origin time of isolated events, with a gain in recall only when a hit must also match the origin
time. On this evidence the headroom divides into a part that is the accuracy of the station terms and a
part that is amplitude, weak arrivals in clutter and coda, and the open problem is calibration and
amplitude evidence, not a further location principle. These are pilot estimates on one seed and bounds
from oracles that a real network does not have. One preregistered test has since confirmed the calibration part within a narrow scope,
hypothesis G1-C (Task 2L, Appendix Y.28). After an exploration on the calibration seed whose plan was
frozen with a hash before any detection result, station terms estimated without picks, by a joint
determination of the hypocentres of the twenty catalogued events and of the terms from matched-filter
energy peaks, replaced the registered corrections inside the capped moment-tensor detector with the
local timing refinement, on two seeds never used before, with the registry hash sent before any run and
the code frozen. With a hit required to match the origin time within 0.5 s, recall at 0.25 false events
per window rose by $+0.057$ [$+0.034$, $+0.081$] and $+0.062$ [$+0.036$, $+0.083$] in standard windows
and by $+0.038$ [$+0.024$, $+0.054$] and $+0.030$ [$+0.015$, $+0.044$] in dense windows, with no loss in
either clutter regime or under the location-only rule within a margin of 0.02, and with the same
direction in block pursuit, $+0.068$ and $+0.071$. An independent verification found no leakage and
reproduced the rules with its own code. Its scope findings belong to the statement: the pooled gain is
essentially that of the one-sided eight-station cell, $+0.108$ and $+0.113$, and with four stations no
gain is demonstrated; the estimated terms are worse than the registered ones in about 22 per cent of the
networks; the gain is conditional on the refinement, the exploration having found nothing measurable
without it; the estimator assumes catalogue errors of known size, 1 km and 0.1 s; it captures 0.28 to
0.61 of the gain of the true delays; the evidence is synthetic and from one generator family; and joint
hypocentre determination with station terms is classical in location from picks, so that what is
measured is their estimation without picks and their effect inside a pick-free detector, not an
algorithm. The
benchmark and the real-data test G1' are contributions whatever the outcome of G1 on real data.
The mathematical contribution is the theory of that convex problem and its certificate
(Appendix AA states what is proved about it and what remains conjecture) together with a
leading-order observability certificate, computable from geometry alone, that predicts the
conditioning of the source problem and maps in advance where false events are likely. The scientific
contribution is a sufficiency ladder in which physical evidence buys rejection rather than
resolution and the physics required is a function of the effective kinematic information, set by
station count, picker miss rate and timing precision. The structural contribution is a division of
labour among constraints. G1 is tested on real data by decimating a dense part of the Colombian
network, after an orientation audit and a measurement of real pick errors, against pick-based
associators and against waveform-migration locators. For the kinematic likelihood, a search over the
hypothesis space pruned with Lipschitz bounds returns catalogues identical to the dense search with
two to three times fewer evaluations and 1.7 to 2.0 times less wall-clock time; for the moment-tensor
statistic the certified search is exact but prunes almost nothing, so its value is exactness, not
speed (Appendix Z). Inference uses primal-dual constrained optimization; a learned representation is
an acceleration layer behind a gate, not a contribution of the thesis.

The first registered contact with real data qualifies everything above (Appendix Y.30). On one hour of
the 2019 Ridgecrest aftershock sequence, with four subnetworks of four to eight broadband stations and the
analyst catalogue as truth, the refined detector of G1-C recalled at most 0.022 of the 45 events of
magnitude 2 or more at 30 false detections per hour with origin times matched, against up to 0.667 for
PhaseNet with PyOcto; the pooled difference from the best pick-based comparator was -0.383 [-0.481,
-0.310], and a reading tolerance widened by rule to the measured timing level did not change it. Two
mechanisms were measured. The homogeneous medium with constant station terms leaves timing errors of
0.35 to 0.41 s that grow with the separation between sources, beyond the 0.25 s reading that the
synthetic generator never challenged, so detections land 8 to 18 km from their events in median; and a
sequence that fills every window with energy saturates the capped statistic, so the scores of real events
do not stand out. Pick-based association with tolerances of 0.5 to 2 s absorbs that timing error. An
exploration with real noise and real geometry around injected synthetic events pointed the other way for
G1-R (Appendix Y.29c), and a registered test on noise of days that no earlier analysis used confirmed it
(Task 2M, Appendix Y.31): on a day with transients the cap and the Huber fit gained +0.228 [+0.148,
+0.298] and +0.241 [+0.179, +0.293] over BOMP-PK, on near-Gaussian days +0.057 and +0.120, and the
advantage of the cap is concentrated in the hours with transients. An exploration on the development
segment of the Ridgecrest hour then separated the two R0 mechanisms (Appendix Y.32). A layered or fitted
medium does not make the energy peaks more predictable: a quarter to a half of the readings of the larger
events fall outside 0.25 s whatever the medium. The cap does make real events indistinguishable from the
background even at their true hypocentres, and a logarithmic compression of the arrival energy restores
most of that ranking while keeping most of the robustness to transients on one exploratory seed, but neither change
locates the events.
The thesis therefore defends the benchmark, its positive and negative results and the measured
conditions under which pick-free physics helps; transfer to real sequences is conditional on a reading of
arrivals that survives dense seismicity, and the next registered real test waits until a candidate passes
a gate fixed in advance on development data, before the Colombian test G1'.

---

# 2. Scope and Significance

The target is epicenter, depth, origin time, event count and approximate mechanism. Earthquake
prediction, dynamic rupture simulation, tomography as an independent goal and unrestricted
full-waveform inversion of Earth structure are outside the claim.

The physical model is 3-D isotropic linear elasticity with spatially varying density and Lame
parameters. Anisotropy and nonlinear constitutive behaviour are extensions. Anelastic attenuation
is treated separately: because predicted amplitude governs the detection model of Section 6.2, an
elastic model cannot calibrate detectability, so that model is empirical (Section 6.2).

**Where the result matters, globally.** The regimes in which the measured advantage is largest -
sparse networks and high false-arrival rates - are not a local peculiarity.

- **Treaty monitoring.** The global monitoring network is sparse, inter-station coherence is weak,
  and pipelines commonly produce events by combining detections from unrelated small events;
  analysts reject them using precisely the non-detecting stations that automatic association
  ignores [23].
- **Aftershock sequences.** The first hours after a mainshock are when overlap is maximal, catalogue
  quality is worst and operational demand is highest.
- **Induced and volcanic seismicity.** High event rates, dense clustering, small magnitudes.
- **Ocean-bottom and temporary deployments.** Few stations, high noise, no dense backbone.
- **Planetary seismology.** The InSight mission operated a **single** seismometer on Mars. From it,
  174 marsquakes were identified in the first ten months of data, in two distinct populations: 150
  small-magnitude high-frequency events propagating at crustal depths and 24 low-frequency
  subcrustal events of magnitude $M_w$ 3 to 4 [85]. Locations were obtained from the P-S
  differential time and the polarization back-azimuth. The single-station problem is canonical
  there, not a workaround, and the two-population structure is precisely a cardinality-and-clutter
  problem under one instrument.

**Five forms of significance.** Statistically, how robust likelihoods and explicit background models
affect structured false-positive attribution. Mathematically, constrained non-convex inverse
optimization with a PDE state and physically generated feasible sets. Geophysically, which physical
constraints carry attribution information and in what regime. Computationally, whether a
differentiable representation reduces repeated forward and storage cost without changing
conclusions. Operationally, whether epicenters can be determined where current pipelines require
three or more detecting stations.

---

# 3. State of the Art and Research Gap

## 3.1 Travel-time association

GaMMA established a Bayesian mixture formulation [3]; Neuma added an explicit noise class and a
physics-informed Eikonal forward model [4]; REAL made grid search with pick counting the
lightweight standard [88]; PhaseLink and graph-neural association introduced learned relational
strategies [5], [6]; PyOcto is the high-throughput reference [7]; EikoNet made differentiable
travel-time fields practical [8]; HARPA combined travel-time neural fields with an unbalanced
optimal-transport comparison of observed and predicted arrival sets [9]; VORA uses Voronoi station
neighbourhoods with optional sub-clustering for overlapping events [19]. All consume arrival times
or arrival distributions. The reported open problems are high event rates, overlap, and the
inability of travel-time information alone to reject spurious alignment with false picks [6].

Two recent comparisons fix what is and is not known about false events. A synthetic benchmark of
PhaseLink, REAL, GaMMA, GENIE and PyOcto with injected noise picks and high event density found
GENIE and PyOcto the most robust, with F1 scores above 0.8 even in the hardest scenarios; PhaseLink
degraded most, to near zero in the most complex subduction cases, while GaMMA outperformed PhaseLink but
struggled with accuracy and scalability at high noise and event density, and REAL lost recall under
extreme conditions [16]. A
comparison on dense nodal arrays found small differences between methods below about a hundred
stations and false events generated by errant picks in very large networks [89]. Neither study
addresses the regime this proposal targets, **sparse or one-sided networks with a high false-pick
rate**, and neither uses three-component polarization. In operation, a PhaseNet plus GaMMA pipeline
in Southern California reported a false-event rate of about five per cent [101], and regional
machine-learning catalogues are documented to add poorly constrained small events and
non-earthquake signals [100]; an automated catalogue of the entire Chilean margin, built on a
trench-parallel and therefore largely one-sided network, is the closest regional analogue to the
geometry class of Appendix M [102]. For Colombia, SeisMonitor applied machine-learning picking and
association to three national networks over 2016-2022 with detection comparable to the standard
catalogues [90]; it is the local reference against which any real-data result here must be
compared.

## 3.2 Waveform-level monitoring: the closest prior art

Two systems must be confronted directly because they already perform Bayesian inference over an
unknown number of events. NETVISA replaced heuristic network processing with Bayesian inference in
a model including probabilities of false and missing detections, with a large reduction in missed
events relative to the previous production system [17]. Its generative model already scores more
than arrival time: arrival azimuth and slowness enter through station- and phase-specific residual
distributions, false arrivals are a Poisson process per station with attributes drawn independently
of any event, and coda arrivals are a separate class correlated with their parent arrival [17]. Its
operational use at the International Data Centre still produces substantially more automatic events
than survive analyst review [106]. **Per-station azimuth likelihoods and an explicit false-arrival
model are therefore established prior art**, and the polarization channel of this proposal is
positioned against them in Section 6.4. SIGVISA extended this by modelling
waveforms rather than discrete detections, with interpretable latent variables for arrival times
and amplitudes, waveform structure represented by Gaussian processes over wavelet coefficients,
and reversible-jump MCMC over event histories, recovering substantially more events and reducing
mean location error by a large factor [18].

**Waveform-level Bayesian monitoring with unknown cardinality therefore already exists and is treated as established prior art.** The boundary is the nature of the forward model. SIGVISA's
waveform prediction is statistical: it learns waveform detail from historical events at the same
location and degrades to parametric envelopes where no historical seismicity exists, which makes it
effectively a learned template prior that excels at repeating sources. The forward model here is
physical: an elastodynamic PDE with a moment-tensor source, which needs no historical seismicity at
the candidate location and predicts rather than recalls. The registered prediction is that the two
dominate in different regimes, and locating that crossover is itself a contribution.

The detector of Section 7.5 has closer antecedents than these Bayesian systems, and they lie in the
microseismic and source-inversion literature. Vera Rodriguez, Sacchi and Gu [129] posed the
simultaneous recovery of origin time, hypocentre and moment tensor as a sparse representation over a
dictionary of Green's functions indexed by hypocentre and origin time, with one six-component
moment-tensor block per atom, solved it by a modified block orthogonal matching pursuit, recovered
several simultaneous sources and proposed the scheme for continuous monitoring; the same authors
placed it in a compressive-sensing framework [130]. Sparsity-aware localization of multiple
microseismic events blind to the source time function [131] and fast sparsity-promoting microseismic
source estimation [141] belong to the same line. The moment tensor as a linear block unknown, group
sparsity as the relaxation of the event count, and greedy deflation over several simultaneous sources
are therefore prior art. What these works do not contain, as far as the targeted searches found, is a
formulation over measures off the grid, a dual certificate used as the detection statistic, a
false-alarm threshold calibrated on noise-only windows, or a certified search of the hypothesis space;
they also invert full waveforms on a grid, whereas the reduced model of Section 7.5 uses P and S peak
amplitudes.

Station time shifts estimated jointly with the moment tensor are prior art as well, in source inversion
for one known event. The Cut-and-Paste method [143] fits body-wave and surface-wave windows separately
and lets each window shift in time to absorb velocity-model error while it searches the mechanism, and
a gradient-based joint inversion of a point-source moment tensor, its centroid depth and
station-specific time shifts, with the misfit made differentiable through the spectral domain, treats
those shifts as proxies for unmodelled three-dimensional structure [144]. Both estimate the source of
an event already detected and located approximately; neither uses the shifts inside a detection
statistic over an unknown number of overlapping events, calibrates a false-alarm rate or certifies the
search over shifts, which is the use that hypothesis G1-M made of them and that the preregistered Task 2E
tested (Appendix Y.22).

Continuous grid-based moment-tensor monitoring is operational. GRiD MT, proposed by Kawakatsu [132]
and operational since 2003 [123], fits a moment tensor at every node of a spatial grid on long-period
records and declares an event from the variance reduction of that fit. Up to normalisation, the
variance reduction is the matched-subspace generalized likelihood ratio over the six-dimensional
moment-tensor subspace [118], [137], which is the certificate of Section 7.5 at zero; the identity of
Proposition AA.2 gives an optimisation meaning to an operational statistic rather than introducing a
new statistic. Joint location and mechanism by mechanism-corrected stacking is mature on dense
microseismic arrays, through moment-tensor migration imaging [121], joint location and mechanism
inversion [122], [133], [134] and back-projection stacking of P and S waves [112], with a systematic
comparison of migration-based location and detection methods [135] and a review of waveform-based
location across scales [136]. Master-event waveform stacking already estimates station and
travel-time corrections from reference events [138], which anticipates the catalogue corrections of
Appendix Y.16. Among open migration codes QuakeMigrate [109] is the maintained reference and is the
one used as a comparator here; on the association side GENIE [6], a graph-neural associator, and HARPA
[9], which compares arrival sets by unbalanced optimal transport over travel-time neural fields, are
the nearest learned neighbours, and both consume picks.

Pick-free waveform location by regression also exists - ArrayConvNet, MSLOC, SUGAR, Loc-Flow,
QNO [10], [21] - and network-wide multi-station picking is established [13]. These are regression
or classification, without a constrained inverse formulation or explicit cardinality; the
multi-task graph approach PLAN reports that with several events in a window it returns only the
first [20]. Overlapping-signal separation in seismology is done by blind, audio-style source
separation with no physics and no attribution [22].

## 3.2b Single-station polarization, and why it fails at some stations

Single-station back-azimuth from three-component data is established, classically and with deep
learning, including estimators that return a confidence [97], [98], a receiver-function rotation
that resolves the 180-degree ambiguity [105], and single-station location with neural priors inside
an eikonal mixture likelihood [104]. Its failure modes are also documented. Horizontal components
of permanent stations are misoriented more often than assumed: in a national network of Iran 36 of
121 sites changed orientation over time [91], in Thailand eight of 45 surface stations were off by
more than fifteen degrees and most borehole sensors were misoriented [92], and a network-wide
catalogue for Australia found orientation changes after equipment swaps [93]; P-wave particle motion
[91] and Rayleigh-wave arrival angles [94] are the standard independent checks. Where orientation is
correct, P polarization is still frequency- and back-azimuth-dependent through anisotropy and
lateral heterogeneity [96] and is dominated by the shallow structure beneath the station [95]. No
published orientation study of the Colombian network was found by the search of Appendix E. These
results bear directly on Appendix T.1, where single-station P polarization passed its reliability
gate at one of five Colombian stations: misorientation, near-surface structure and incidence
geometry are the three candidate causes, and Section 16 registers the tests that separate them.

## 3.3 Non-detection evidence

Non-detection is established, not new. NETVISA incorporates missing-detection probabilities [17],
older network-processing work proposed association-strength measures comparing observed and
expected amplitudes at detecting and non-detecting stations, and recent work notes that automatic
pipelines still ignore non-detecting stations as if the missingness were uninformative while human
analysts use exactly those to reject false events [23]. What is new here is the change of status:
Section 6.2 replaces the calibrated Bernoulli treatment with a **certified exclusion set under a
declared sensing and noise model**, motivated by unique continuation but not derived from it, and
retains the Bernoulli treatment as the softest rung of its relaxation ladder.

## 3.4 Differentiable wavefield representation

Foundational PINNs [55]; wave propagation and FWI with PINNs [56]; semi-infinite domains and
absorbing boundaries [57]. The documented limitations are treated here as design constraints, not
as risks to be discovered: spectral bias toward low frequencies [62], loss-imbalance pathologies
analyzable through the neural tangent kernel [61], and characteristic failure modes on
high-frequency problems [60], with Fourier feature embeddings as the standard mitigation [63].
Frequency-domain and operator approaches are stronger: a U-shaped neural operator for the 3-D
elastic wave equation in the frequency domain reduced memory by a factor of forty relative to a
time-domain formulation, ran about a hundred times faster than a spectral-element solver, and
reduced FWI time by a factor of about three hundred and fifty [58]; meta-learned low-rank PINNs
link frequency content to the singular values of layer weights [65]; physics-informed neural
operators for wavefields remain early-stage with scalability and labelled-data availability as the
open constraints [66]. The documented search of Appendix E did not identify a published result
meeting all of the following simultaneously: time domain, three dimensions, heterogeneous elastic
medium, seismological frequency band, and waveform-level accuracy sufficient for the likelihood used
here. This proposal therefore does not build on one
(Section 10).

**Two novelty claims that are explicitly *not* made.** Multi-source PINN forward modelling exists:
an encoded multi-source physics-informed network for wave propagation with sources at arbitrary
positions and excitation times was published in 2026 [87], so multi-source PINN propagation is not
available as a contribution. Single-station earthquake location from waveform-derived epicentral
distance and back-azimuth also exists, including recent machine-learning implementations trained on
global catalogues and fine-tuned regionally [86], so single-station location is not available as a
contribution either. The remaining contribution investigated here is the
**multi-event** single-station case (C3) and the constraint-generation framework of Section 6.

## 3.5 Structure that this proposal imports

Recovery theory for sparse measures was developed for band-limited, essentially
translation-invariant operators [46], [47], with the exact-support result of [47] proved in one
dimension and the multidimensional case, in which the relative position of the spikes governs
stability, analysed in [140]; all of it concerns scalar amplitudes. The seismic operator is neither
translation-invariant nor scalar-valued, so those theorems do not transfer: Section 7.5 and
Appendix AA borrow the formulation and the form of the non-degenerate source condition, and leave
recovery over the continuum as a conjecture. Coherence-based and convex-geometric
demixing guarantees [48], [49] and phase-transition analyses [50], [51] do not require translation
invariance. Multi-object estimation with missed detections and clutter is mature [52], [53], with
GOSPA as the correct set metric [54]. Joint parameter estimation over overlapping *physical*
waveforms is established in gravitational-wave analysis through hierarchical subtraction, joint
estimation and normalizing flows [70], [71], and has not been imported into seismology.

## 3.6 The gap

Within the databases, query families, date range and inclusion criteria documented in Appendix E,
no indexed study was found that poses epicenter determination as a constrained inverse problem whose feasible set is
generated from elastodynamic structure - observability, medium-independent kinematic admissibility,
polarization rank and source-spectrum bounds, with a non-detection exclusion that unique continuation
motivates but does not imply - or that measures which of those constraints carries the attribution
information. That is the gap, and Section 6 is the response. The gap is
narrower than it first appears, and is stated narrowly: per-station azimuth likelihoods and
explicit false-arrival models exist [17], physics-informed and transport-based associators exist
[4], [9], and single-station back-azimuth estimation exists [97], [104]. What is not found is (i) a
hard kinematic admissibility test valid in any medium with a velocity lower bound, combined with
(ii) a three-component polarization likelihood whose per-station reliability is measured on
catalogued events with a chance correction, and (iii) an evaluation of both on the catalogue
endpoint in sparse and one-sided networks against tuned published associators.

For the waveform detector of Section 7.5 the gap is narrower still, because most of its ingredients
are anticipated by the work reviewed in Section 3.2: block-sparse moment-tensor inversion [129], [130],
grid moment-tensor monitoring [132], [123], mechanism-corrected stacking [121], [122], [133], [134],
[135], station corrections from reference events [138] and extreme-value thresholds [125]. The targeted
searches, which are not an exhaustive review, found no prior for four elements: the off-the-grid
group BLASSO over moment-tensor-valued measures posed as earthquake detection and location; the
explicit identity between its certificate at zero and the matched-subspace generalized likelihood
ratio over the six-dimensional moment-tensor subspace of P and S peak amplitudes, a modest unifying
observation that follows from the optimality conditions of the group lasso; an exact search of the
certificate with a certified cell bound; and split-conformal calibration of the network false-alarm
rate, which transfers conformal outlier testing [139] to network detection. They also found no paper
from 2023 to 2026 reporting a controlled, preregistered benchmark of a mechanism-aware pick-free
method against GaMMA, PyOcto and QuakeMigrate on overlapping events; the protocol of [16] is the
template for one. Absence in these searches is the claim, not priority.

Methods that locate without picked travel times were mapped separately, because a channel independent
of the arrival times is the natural candidate for the events that the detector places 5 to 10 km from
the truth (Appendix Y.27b). Of eight families, four carry information that is independent of the arrival
times on an existing three-component network: the direction of P polarization, as in single-station
back-azimuth estimation by polarization analysis [145]; the amplitude decay between stations, as in
amplitude source location [146]; the amplitude ratio between phases at one station; and high-frequency
envelope energy [147]. Rotation and wavefield gradients are independent as well and need instruments
that a sparse network lacks [148]. The other families hide the kinematics: time-reverse and
wave-equation imaging in the propagator, with requirements on aperture and station density that four to
eight stations do not meet [149], matched field processing and back-projection in the phase of the
Green's function [150], and end-to-end learned location in its weights [97]. The verdict for this thesis
is negative and measured. In the benchmark the independent channels are limited by information, 4.4 to
5.6 km against 0.4 to 0.7 km for the arrival times by a Cramer-Rao bound, and relocating detections with
them lost recall (Appendix Y.27, exploratory). The map is a scan and not a systematic review; more than
a third of its entries were read at the level of bibliographic metadata, and the protocol of Appendix E
governs any claim of absence drawn from it.

That benchmark must contain the methods that anticipate the formulation. Besides GaMMA [3], PyOcto
[7] and QuakeMigrate [109], the preregistered Task 2C ran block orthogonal matching pursuit in the form
of [129] and mechanism-corrected stacking in the form of [122] on the same windows, hypothesis grid,
station corrections and hit rule (Section 1.6, Appendix Y.20). Both read the matched filter at the
exact predicted times of a shared grid of 3 km and 0.1 s, and the preregistered Task 2D gave both the
peak search of the detector within $\pm0.25$ s (Appendix Y.21). With equal tolerance block orthogonal
matching pursuit matched the detector in every cell, so the benchmark supports an advantage of
timing-tolerant moment-tensor detection with deflation, a family that includes [129] so equipped, over
pick-based associators and QuakeMigrate, and no advantage of the detector over sparse moment-tensor
inversion; a finer grid, the full-waveform dictionaries of [129], the stacking of [133] and [112] and
MALMI [110] have not been run. For the detector, the moment-tensor pursuit with deflation is prior art,
and the principled treatment of the timing tolerance with the time shifts as nuisance parameters of the
convex problem, the element that the targeted searches did not find beyond the one-event inversions of
[143] and [144], did not beat block pursuit with equal tolerance in the preregistered Task 2E (Appendix
Y.22). What remains is narrower: the measured and explained effect of a bounded-influence
moment-tensor subspace statistic in dense sequences of overlapping events, shown with a per-station cap at a
calibrated level and with a Huber fit at a level fixed by rule from the noise (Appendices Y.24, Y.26).
Clipping or normalising the contribution of each station is established practice in robust stacking and
back-projection, the family reviewed in [136], and Huber fitting is standard robust statistics, so the
gap that result fills is empirical rather than algorithmic, and the search protocol of Appendix E is to be executed before any claim of novelty is made
for it. Because both antecedents also exploit the radiation pattern, the
contribution of the radiation operator is isolated only by the ablation Task 2B, which replaces the
operator by the identity inside the same solver and so does not separate the sign of the pattern from
the polarization constraint and the amplitude geometry (Appendix Y.19).

**Table 1. Research boundary.**

| Line | Established capability | Role here |
|:--|:--|:--|
| Travel-time association | GaMMA, PyOcto, HARPA, VORA | Counterfactual baseline |
| Waveform coherence | Source scanning, matched filter | Cheap waveform baseline |
| Waveform migration | QuakeMigrate | Comparator in every cell |
| Sparse moment-tensor inversion | Block OMP over Green's-function dictionaries [129], [130]; sparsity-promoting multi-event estimation [131], [141] | Closest prior art for the detector; run as block OMP on the shared grid with exact times (App. Y.20) and with the detector's peak search, where it matches the detector (App. Y.21) |
| Moment-tensor inversion with station time shifts | Cut-and-Paste [143]; gradient-based joint inversion of moment tensor and station-specific shifts [144] | Prior art for the time shifts of G1-M; one known event, no detection statistic; as nuisance parameters of the detection certificate they did not beat block pursuit (App. Y.22) |
| Robust stacking and back-projection | Per-station normalisation, clipping and sign stacking [136] | Prior art for the per-station cap whose effect in dense sequences G1-R measures (App. Y.24); Huber and other M-estimators, standard in robust statistics, are prior art for the fit of Task 2I (App. Y.26) |
| Grid moment-tensor monitoring and mechanism-corrected stacking | GRiD MT; moment-tensor migration; joint location and mechanism | Prior art for the statistic and for mechanism awareness; stacking in the form of [122] run on the shared grid, with and without the peak search (App. Y.20, Y.21) |
| Bayesian monitoring | NETVISA, SIGVISA | Closest prior art; per-station azimuth likelihood and false-arrival model (NETVISA); statistical waveform model (SIGVISA) |
| Single-station direction | P polarization, BAZ networks, orientation checks | Reliability of C-IV per station; orientation audit before gating |
| Pick-free location | ArrayConvNet, MSLOC, QNO | Regression reference |
| PINNs / neural operators | Frequency-domain 3-D elastic operators | Differentiable representation |
| Sparse recovery, RFS | BLASSO, PHD/PMBM, GOSPA | Imported machinery |
| PDE control theory | Observability; unique continuation | **Constraint generation from observability, no prior found (Appendix E); unique continuation motivates C-II only** |
| Constrained geophysical inversion | Sosa and co-authors [72], [73] | Optimization foundation |
| **Proposed integration** | Physically constrained epicenter determination; a bounded-influence moment-tensor subspace statistic, with per-station clipping at a calibrated level or a Huber fit at a level fixed by rule from the noise, certified search and calibrated threshold | Contribution sought |

---

# 4. Questions, Hypotheses, Objectives

## 4.1 Primary question

One primary question, one methodological question and one secondary question, in the order in which
the evidence places them.

**Primary.** In the regime where the picker discards kinematic information (few stations, missed
arrivals, pick timing error), does a timing-tolerant moment-tensor detector with deflation, applied to
the three-component waveform before picking, detect more events than published associators and
waveform-migration locators at the same realized false-event rate? The family includes the greedy
detector MT-FW of Section 7.5 and block orthogonal matching pursuit [129] given the same peak search
(BOMP-PK), which match each other on the synthetics (Appendix Y.21). This is hypothesis G1 of Section
4.3. Its endpoint,
event recall at a fixed realized false-event rate per window read on the exact FROC curve at 0.1 and
0.25 false events per window, and on continuous real data per day, is the one primary endpoint of the
thesis (Section 13.0b).

**Methodological.** Does the group BLASSO of Section 7.5 with a time shift per station and phase as a
nuisance parameter, solved off the grid by sliding Frank-Wolfe, detect more events than BOMP-PK at the
same timing tolerance and realized false-event rate? This is hypothesis G1-M of Section 4.3, and the
preregistered Task 2E answered no in this synthetic design (Appendix Y.22). The narrower methodological
question that two preregistered tests answered yes is G1-R: whether bounding the contribution of each
station and phase in the moment-tensor subspace statistic raises recall in dense sequences of overlapping
events (Appendices Y.24, Y.26). A third question, on calibration, one preregistered test answered yes
within a narrow scope, G1-C: whether station terms estimated without picks from the catalogued events
raise the origin-time-matched recall of the detector with a local timing refinement (Appendix Y.28).
These questions decide what methodological advance the thesis contains beyond
the benchmark, not whether G1 holds.

**Secondary.** What is the *minimum* physical information sufficient to reject false hypotheses
reliably, as a function of the effective observation count? Concretely: which channels of
elastodynamic structure carry the rejection, is there a proper subset of them that is statistically
indistinguishable from the full wave-equation treatment, and does that subset degrade more slowly
under velocity-model error than methods that require a velocity model? Its endpoints, the
false-association rate at non-inferior recall and the discrimination margin, are secondary: they
decompose which physics earns G1 and map where it holds, and they cannot by themselves make the
thesis succeed or fail.

The pre-studies of Appendices G to S repeatedly answered the secondary question at the level of the
margin, and the answer they converge on is *less than a full wavefield*. Whether a larger margin
becomes more detections at equal false alarms is a catalogue-level question, and the catalogue-level
programme of Appendices X and Y converged on G1 as the claim that must survive; that is why G1 is the
primary question and sufficiency the secondary one.

## 4.2 Secondary questions

1. Which of the five constraint families of Section 6 carries the attribution information, and in
   what regime does each become active or inert?
2. How does the advantage scale with false-arrival rate, station count and network aperture?
3. Under what conditions can an event be formed and characterized from a single three-component
   station in the presence of other, overlapping events?
4. Does an optimal-transport misfit with multiscale continuation remove the cycle-skipping
   obstruction that a time-domain quadratic misfit exhibits under realistic velocity-model error?
5. How accurate must the differentiable representation and its source gradients be before
   approximation error changes conclusions?
6. Can posterior or approximate-posterior uncertainty remain calibrated under separate statistical
   and physical mismatch, given that parametric waveform likelihoods are documented to be
   overconfident under theory error [27], [28]?
7. When is a learned representation preferable to analytic beam synthesis or a stored
   strain-Green's-tensor database, once training, storage and interpolation error are accounted for?

## 4.3 Falsifiable hypotheses

Margins are registered defaults, and the complete table of thresholds is Appendix D. The margins of
the constraint-attribution hypotheses are frozen after the Tier-1 pilot, which supplies their paired
standard deviations. That pilot does not test G1 (Section 13.0): the thresholds of G1, G1-M, G1-R and G1-C
were fixed in the registration of each synthetic task before its run, and that of G1' is fixed before
any RSNC evaluation, with its sample size taken from the standard deviation measured at three, four
and six stations. Several hypotheses below are stated as **registered nulls**, because the pre-studies of
Appendices G to S measured against them; reporting them as nulls rather than omitting them is what
makes the remaining positive claims interpretable.

### The structure of the thesis in one figure

```
                    SEISMIC OBSERVATIONS
                            |
                            v
                 How many EFFECTIVE
                  observations do I have?
                            |
              +-------------+-------------+
              v             v             v
            many          few          extreme
              |             |             |
              v             v             v
        travel-time    polarization    waveform
          sufficient      needed        needed
              |             |             |
              +-------------+-------------+
                            v
                  PHYSICAL FEASIBILITY
                            |
          +-----------------+----------------+
          v                 v                v
      kinematics       polarization     observability
      model-robust     clutter reject     geometry
          |                 |                |
          +-----------------+----------------+
                            v
                  FALSE-EVENT REJECTION
                            |
                            v
                    SOURCE ATTRIBUTION
```

The upper half is the **sufficiency ladder**: how much physics is required is a function of the
effective observation count $N_{\mathrm{eff}}$, which is reduced both by sparse deployment and by
picker failure. The lower half is the **division of labour among constraints**: each one does a
different mathematical job in identifying the feasible set. The output is **false-event rejection**,
which is what the physics buys, rather than spatial resolution, which it does not.

### Three organizing ideas, subordinate to G1

**1. Physics buys rejection, not resolution.** Resolution between nearby true sources improves by
$6$-$14\%$ while the discrimination margin improves fourfold (Appendix H). Whether that margin
becomes fewer false events in an emitted catalogue, and against which comparator, is what
Appendices X and Y measure; the answer is regime-dependent and is stated there.

**2. There is no universal minimal set; there is a ladder.** The correct question is not "which
constraints are minimal" but "what is the minimum physical information as a function of
$N_{\mathrm{eff}}$". Station count, picker miss rate and pick timing precision all enter
$N_{\mathrm{eff}}$, which is why exact picks and realistic picks give different answers (Appendix Y,
D1).

**3. The constraints do different jobs.** Polarization supplies discrimination margin but not model
independence; kinematic admissibility supplies model independence but not discrimination;
observability supplies identifiability and design, and does **not** classify hypotheses. That
division of labour, not a list of five modules, is where the applied mathematics lives.

### The one claim that must survive

$$\boxed{\begin{array}{c}\text{This dissertation tests whether timing-tolerant moment-tensor detection with deflation}\\
\text{detects more events than published associators and migration locators at the same realized}\\
\text{false-event rate where the picker discards kinematic information (G1), reports that a group}\\
\text{BLASSO with station time shifts off the grid does not beat block pursuit (G1-M) and that a}\\
\text{bounded-influence statistic does in dense sequences (G1-R) and station terms estimated}\\
\text{without picks do in the refined detector (G1-C), and characterizes}\\
\text{where this holds as a function of the effective observation count.}\end{array}}$$

Everything else is subordinate. If G1 fails on real data, the thesis fails as a claim about
detection, and what remains is a characterization of where physical evidence does not pay; if it
holds, the rest is elaboration. G1-M failed on synthetics, and that outcome is reported (Appendix Y.22);
the methodological contribution is the bounded-influence statistic of G1-R with its measured effect
and mechanism, supported by Tasks 2F and 2I, the calibration result of G1-C, supported by Task 2L, the benchmark with its negative results and the real-data test, and the theory of Appendix AA is the
mathematical contribution.

### Component hierarchy

Stated explicitly so that no component is mistaken for the thesis.

| Tier | Components |
|:--|:--|
| **Core** | detection at a fixed realized false-event rate by timing-tolerant moment-tensor detection with deflation (G1, primary); the bounded-influence moment-tensor subspace statistic in dense sequences, with a per-station cap at a calibrated level or a Huber fit at a level fixed by rule from the noise (G1-R, methodological; Tasks 2F and 2I); station terms estimated without picks inside the detector with a local timing refinement (G1-C, methodological; Task 2L); sufficiency ladder; division of labour among constraints |
| **Mathematical support** | observability to coercivity and local identifiability; model-independent kinematic feasibility |
| **Demonstrative cases** | one-sided geometric failure class; single-station multi-event |
| **Extensions** | certified exclusion from silence; source-time admissible cone |
| **Conditional, behind gates** | optimal transport misfit; Gaussian beams with learned residual |
| **Abandoned as claims** | resolution gain (E1); $V_P/V_S$ as the controlling parameter (E2); second-order rank as a cardinality estimator (E5); velocity-model-free *location*; equal importance of all five constraint families; the off-grid group BLASSO with station time shifts as a better detector than block pursuit (G1-M) |

The abandoned row is not a list of failures. Each entry is a registered null measured in the
exploratory programme, and each one narrows the thesis to what the evidence supports.

**A transversal methodological result.** Across studies the informative endpoint changed with the
regime: where resolution saturated, the false-event rate discriminated; where the false-event rate
saturated at zero, the discrimination margin did. The principle that follows is registered as a
contribution in its own right: **endpoint selection must follow the information regime, not
convenience**, and which endpoint saturated must be reported rather than treated as a null.

### Evidence taxonomy

Five labels state the status of every hypothesis in the master table below and no others are used
there, so that the kind of claim each result supports is explicit.

| Label | Meaning |
|:--|:--|
| **Theorem** | Proved, with its hypotheses stated |
| **Proposition** | Proved for a declared reduced model; extension is conjecture |
| **Numerical evidence** | Measured under one generator; replicated where stated |
| **Preliminary** | Measured, but below the replication or power required to establish it |
| **Confirmatory target** | Registered, not yet run |

**Exploratory versus confirmatory.** Every pre-study in Appendices G to X, and all three protocols
of Appendix X.1, are **exploratory protocol development**. None is confirmatory, including the
final protocol P3: it was adopted after observing that P2 produced an implausible result, and
although the argument for a relative threshold does not depend on any outcome, the order of events
was that. The confirmatory test of G1 on real data is the network-decimation test of Section 16, with
margins frozen in advance, and it has not been run; the Tier-1 pilot of Section 13.0 calibrates the
pipeline and does not test G1. The preregistered results in this document are seven synthetic tasks:
Task 2A (Appendix Y.18), whose design, comparators and decision rules were registered before any
evaluation run; its ablation Task 2B (Appendix Y.19); its extension Task 2C to the closest comparators,
the grid solver and thirty networks per cell (Appendix Y.20); Task 2D, which gave those comparators the
timing tolerance of MT-FW, ablated its peak search and its cap and froze split-conformal thresholds
before evaluation (Appendix Y.21); Task 2E, which ran the off-grid group BLASSO with time shifts against
block pursuit with equal tolerance on the windows of Tasks 2A and 2C (Appendix Y.22); Task 2F,
which tested the per-station cap in dense windows on two seeds never used before, after an exploration
on the calibration seed that is reported as exploratory (Appendices Y.23, Y.24); and Task 2I, which
tested a Huber fit with a level fixed by rule from the noise on two further seeds never used before,
after a second exploration reported as exploratory (Appendices Y.25, Y.26), and whose secondary arms
are descriptive. Two further explorations on the calibration seed, on evidence that does not use travel
times and on the timing budget, are reported as exploratory (Appendix Y.27), and Task 2L, on two more
seeds never used before, tested station terms estimated without picks inside the refined detector after
a third exploration, also reported as exploratory (Appendix Y.28). Three explorations that followed
Task 2L on the calibration seed, on a double-couple statistic, on label-free false discovery control by
decoys and on self-interference after deflation, are reported as exploratory and test no registered
hypothesis (Appendix Y.29). The registrations of Tasks 2A
to 2C were appended to a local registry file, `PREREGISTRO_1j_1k.txt`, with times typed by hand rather
than read from a clock, and no hash or external record fixes their content at the moment of writing.
For Task 2A the file-system record is consistent with registration before evaluation: the registry
file was created at 20:14:09 UTC on 25 September, the first file of the task at 23:33:22 UTC, and the
operating points of rule 2A-e were frozen at 01:12:25 UTC on 26 September, before the first evaluation
output at 01:12:53 UTC; the typed time of the entry, 20:15 UTC, is not a measurement. For Tasks 2B and
2C the file system cannot date the text of the entries, because it records only the last
modification of the registry. The local transcript of the working session records the two appends at
14:14:04 and 15:04:21 UTC on 26 September, before the first code file of each task (14:16:49 UTC and
about 15:12 UTC), but that transcript is also local and is not an external timestamp; the typed time of
the 2B entry, 14:30 UTC, is wrong, and the erratum is recorded in the registry. The analysis script of
Task 2B, `analyze_2b.py`, was also edited at 14:26:33 UTC, eighteen seconds after its first output, to
add one exploratory line per regime with the fraction of bootstrap replicates in which rule 2B holds;
the rule, the tables and every other line of the output are identical before and after that edit, as
a reconstruction of the earlier script confirms, and the edit is declared in Appendix Y.19 and in the
error register (E36). Tasks 2B and 2C are therefore classified as registered in an append-only local
file with timing that is not externally verifiable, and Task 2A as the same with a consistent
file-system record. Task 2D is the first task with SHA-256 records: its registry entry, a budget note
and its result each carry a clock time and a hash of the registry file in `PREREGISTRO_hashes.txt`, and
the hashes of its code were recorded before any result existed. None of the registrations of Tasks 2A
to 2D was sealed externally. Those of Tasks 2E, 2F, 2I and 2L carry, besides their clock times and hashes, a
record of their content held outside the container: the hash of the registry was communicated in the
working session at the moment of registration, before any code of Task 2E was written and before any
run of Tasks 2F, 2I and 2L, whose code was frozen by hashes before their seeds were generated. That is the
sense in which those registrations are externally timestamped, and its limit: the record is neither
public nor issued by a third party, and a signed commit to a public repository, an OSF registration or
an RFC 3161 token remains the standard for the confirmatory real-data test. These tasks confirm claims about synthetic
data only, and no other result should be read as confirmatory. A preregistered synthetic test does not
create a sixth label: its result is numerical evidence, marked as preregistered. In Tasks 2A to 2C,
in part (i) of Task 2D and in Tasks 2E, 2F, 2I and 2L, the operating point of every method was chosen on the
evaluation windows;
the rule that froze thresholds on other networks (2A-e) fails under its strict reading (Appendix
Y.18), and the split-conformal thresholds that Task 2D froze on empty windows pass cell by cell but
control false alarms only in windows without events (Appendix Y.21); a confirmatory result for G1 at
a deployable operating point exists neither on synthetics nor on real data.

**One primary hypothesis, one primary endpoint.** G1 is the only hypothesis whose outcome can make the
thesis fail; G1-M, G1-R and G1-C decide what methodological advance the thesis contains beyond the benchmark,
and the failure of G1-M is reported rather than fatal. The endpoint of G1, event recall at a fixed realized false-event rate on the exact FROC
curve read at 0.1 and 0.25 false events per window, is the only primary endpoint. Every other
endpoint in the table below (false-event or false-association rate at non-inferior recall, margin,
coherence, coverage, cost) is secondary or exploratory; each is stated where the hypothesis is
defined.

### Master table

Every other section of this document refers to hypotheses by the identifiers in this table and by no
others.

| ID | Statement | Endpoint | Threshold | Status |
|:--|:--|:--|:--|:--|
| **A1** | ordering of channels by contribution is stable | rank correlation across cells | $\ge0.8$ over $\ge20$ cells | **Preliminary**; eight cells only (App. M) |
| **A2** | the registered robust subset C-III + C-IV matches the full stack | false-event rate | within $0.02$, $\ge80\%$ of confirmatory cells | **Preliminary**, in the discriminating one-sided cells only; on real data C-IV is station-dependent (App. M, O, W, T.1, T.2) |
| **A3** | model independence is achievable by an estimator using only a velocity bound and measured directions | degradation trend under model error; absolute accuracy cost | no systematic degradation trend; cost reported | **Preliminary**; instance found, with a cost of about an order of magnitude in localization accuracy. The claim is model-independent *feasibility filtering*, never model-free *location* (App. U) |
| **A4** | minimal set is cost-effective | benefit per unit cost | $\ge90\%$ at $<1\%$ | **Preliminary** (App. N) |
| **B1** | certificate survives interfaces | $\rho$ vs measured $\gamma$ | $\ge0.95$ | **Numerical evidence**, narrowly, with quantitative sharpness degraded under interfaces (App. Q) |
| **B2** | observability design beats azimuthal gap | false-event rate | $\ge20\%$ improvement | **Confirmatory target**, false-event-rate endpoint (App. N) |
| **B3** | $\gamma$ predicts where discrimination is hard | Spearman, local margin | $\ge0.5$ | **Preliminary**, one network, $\rho=0.573$ (App. S) |
| **C1** | advantage decreases with station count | crossover $S^\ast$ | report with CI | **Numerical evidence against the coincidence-count ablation only**, $S^\ast\approx5$-$6$ (App. O); against a continuous kinematic score the four-station rung disappears with exact picks (Y.1), so $S^\ast$ must be re-measured. At three stations the prototype recovered $29/30$ against $7$-$9/30$ for GaMMA and PyOcto in X.4 ($25$ and $27/30$ in the rerun of Y.1), but every method emits false events in $60$-$100\%$ of realizations there, and the kinematic kernel emits fewer than the prototype |
| **C2** | a geometric failure class exists in real networks | false-event rate | $\ge3$ real networks | **Confirmatory target** |
| **C3** | single station, multiple events | count and epicentre | see entry | **Preliminary** for $\le2$ sources, negative for three (App. P) |
| **D1** | pick-free value grows with picker miss rate | recovery rate | report threshold | **Numerical evidence**; threshold in $(0.3,0.5]$ (App. R, W) |
| **E1** | resolution gain $\le15\%$ | separation at fixed coherence | registered null | **Numerical evidence** for the null in the tested setting; replication required (App. H, W) |
| **E2** | $\kappa$ independent of $V_P/V_S$ | resolution constant | registered null | **Numerical evidence** for the null in the tested setting; replication required (App. G) |
| **E3** | overlap alone does not defeat association | ambiguity set | registered null | **Preliminary**, null in the tested setting (App. K) |
| **E4** | no graceful degradation in location | slope, km per 10% error | registered null | **Numerical evidence** for the null in the tested setting; replication required (App. L, W) |
| **E5** | second-order rank does not count sources | eigenvalue gap, order selection | registered null | **Numerical evidence** for the null (App. V.1) |
| **F1** | waveform gate | false-event rate beyond minimal set | $\ge0.03$ | **Numerical evidence** with an oracle matched filter, scoped to $S\le3$; not testable on RSNC data (App. O, Section 16) |
| **F2** | transport misfit beats time-domain | basin convergence | $\ge20$ points | **Confirmatory target**, conditional on F1 |
| **F3** | representation fidelity at level L1 | field, receiver, gradient | see Appendix D | **Confirmatory target**, conditional on F1 |
| **F4** | symmetries exact by architecture | machine precision, budget | see Appendix D | **Confirmatory target**, conditional on F1 |
| **F5** | compute or storage break-even | $N_{\mathrm{break}}$ or ratio | $\le10^4$ or $\ge10\times$ | **Confirmatory target**, conditional on F1 |
| **G1** | a timing-tolerant moment-tensor detector with deflation (MT-FW or BOMP-PK, Section 7.5), applied to the waveform before picking, adds detections over the best published associator and migration locator, **thesis hypothesis** (restated, Family G) | event recall at a realized false-event rate per window, exact FROC curve | $\ge0.05$ | **Refuted on the first real test, supported on synthetics.** Real data: Task R0 (Y.30), one hour of the 2019 Ridgecrest sequence, $-0.383$ [$-0.481$, $-0.310$] against PhaseNet with PyOcto at 30 false detections per hour with origin times matched, so G1 as stated does not hold on a real dense sequence; exploration C18 (Y.32) places the failure in the reading of arrivals, and a registered gate conditions any further real test. Synthetic: geometry-breaking controls (Y.10 to Y.12) and the four preregistered tasks that bear on G1, Tasks 2A to 2D (Y.18 to Y.21): at least $0.05$ at both caps over GaMMA, PyOcto and QuakeMigrate run with the operational picker in all 12 cells, and at $0.25$ over block orthogonal matching pursuit and mechanism-corrected stacking read at exact times, not significant after Holm's correction against the former in the strong-clutter 4-station cells; given the same peak search, block orthogonal matching pursuit matches MT-FW (rule 2D-a fails), so the evidence supports timing-tolerant moment-tensor detection with deflation as a family, not MT-FW specifically; not met at $0.1$ against GaMMA fed optimistic picks with weak clutter, ahead in four of six cells; operating points chosen on the evaluation windows, and no deployable frozen operating point; Tasks 2E, 2F, 2I and 2L concern G1-M, G1-R and G1-C, not the comparison with published methods. Full statement in Section 1.6 |
| **G1-M** | the group BLASSO with a time shift per station and phase, solved off the grid by sliding Frank-Wolfe, adds detections over BOMP-PK at equal timing tolerance, **methodological hypothesis** (Family G) | event recall at a realized false-event rate per window, exact FROC curve, on the windows of Task 2A and thirty networks per cell | $\ge0.05$ at $0.25$ | **Numerical evidence** against it, preregistered (Task 2E, Y.22): $+0.000$ to $+0.039$ over BOMP-PK at $0.25$, every hierarchical interval containing zero, rule 2E-a failing in both regimes; refuted in this synthetic design, with no implementation error or budget handicap found by an independent verification |
| **G1-R** | a bounded-influence moment-tensor subspace statistic, which clips each station and phase vector at a fixed level (MT-FW, Section 7.5) or fits it with a Huber function at a level fixed by rule from the noise (HUBER-N), adds detections over BOMP-PK at equal timing tolerance in windows with three to five overlapping events without harming isolated events, **bounded-influence hypothesis** (Family G) | event recall at $0.25$ realized false events per window, exact FROC curve, dense windows pooled over cells, two fresh seeds per task | $\ge0.03$ with lower bound $>0$ | **Numerical evidence**, preregistered twice. Task 2F (Y.24): $+0.043$ to $+0.074$ over BOMP-PK in both regimes and seeds with the cap, surviving Holm's correction; removing the cap removes the gain; equivalent within $0.05$ in standard windows; $+0.024$ to $+0.045$ with origin times matched; tied to the calibrated cap level. Task 2I (Y.26): $+0.058$ and $+0.053$ with the Huber fit, not inferior to the cap within $0.03$ and equivalent to BOMP-PK within $0.05$ in standard windows when pooled, though not with weak clutter in one seed; $+0.031$ and $+0.019$ with origin times matched; about nine times the cost of the cap; no more events found with all atoms accepted; in its descriptive arms the level is immaterial between 2.5 and 3.4 noise standard deviations and the Huber fit is behind the cap at an amplitude spread of 0.8. Mechanism, bounded influence of each station and phase block under amplitude model error (exploratory, Y.24 to Y.26). Task 2M (Y.31), synthetic events in real noise of independent days of the Alto Tiberina network: $+0.228$ [$+0.148$, $+0.298$] with the cap and $+0.241$ [$+0.179$, $+0.293$] with HUBER-N on a day with transients (primary, one day), $+0.057$ and $+0.120$ on near-Gaussian days (secondary); the cap's gain is concentrated in the hours with transients, HUBER-N's persists in quiet noise; independently verified; in an exploration on seed 999, capping the energies of BOMP-PK recovered almost all of the gap (C18f, Y.32: $+0.026$ [$-0.009$, $+0.081$] left), which suggests that the effect belongs to the per-station bound rather than to the pursuit. Synthetic signals, one real network |
| **G1-C** | station terms estimated without picks, by a joint determination of the catalogued hypocentres and the terms from matched-filter energy peaks, in place of the registered catalogue corrections, add detections to the capped moment-tensor detector with a local timing refinement when a hit must also match the origin time, in standard and in dense windows, without loss per clutter regime or under the location-only rule, **calibration hypothesis** (Family G) | event recall at $0.25$ realized false events per window, exact FROC curve, hits with the origin time within 0.5 s, pooled over cells, two fresh seeds | $\ge0.03$ in standard windows, lower bound $>0$ in standard and dense windows; non-inferiority margin $0.02$ | **Numerical evidence**, preregistered once (Task 2L, Y.28), after an exploration on the calibration seed: $+0.057$ and $+0.062$ in standard and $+0.038$ and $+0.030$ in dense windows, every rule passing in both seeds; $+0.068$ and $+0.071$ in block pursuit; independently verified, no leakage found. Scope: the pooled gain is essentially that of the one-sided eight-station cell, $+0.108$ and $+0.113$, with no gain demonstrated at four stations; terms worse than the registered ones in about 22 per cent of the networks; conditional on the refinement; catalogue errors of known size assumed; 0.28 to 0.61 of the gain of the true delays. Synthetic, one generator family |
| **G2** | clutter scaling | margin degradation ratio | $\le0.5$ | **Numerical evidence** (App. H) |
| **G3** | calibrated uncertainty | coverage error at 90% | $\le0.05$ | **Confirmatory target** |
| **G4** | silence channel adds value | false-source rate | $\ge0.03$ | **Preliminary**, conditional on distance spread (App. J) |

**The dissertation hypothesis is G1**, and everything else is subordinate to it: its endpoint is the
primary endpoint of the thesis (Section 13.0b), and the endpoints of every other row are secondary.
G1-M and G1-R use the same endpoint to decide the methodological questions, and G1-C uses it with hits
that also match the origin time. A1
to A4 decompose *which* physics earns it, B1 to B3 establish *when the problem is well posed at all*,
C1 to C3 and D1 establish *where* the effect lives, E1 to E5 record what was ruled out, F1 to F5 gate
the computational branch, and G2 to G4 are secondary endpoints of the same experiments.

### Family A. Sufficiency: how little physics is enough

This family addresses the secondary question of Section 4.1, how much physical information is
actually required for reliable attribution, and identifies the regimes in which each information
channel becomes active. Its endpoint, the false-association rate at non-inferior recall, is secondary
to that of G1.

**A1 (stable ordering).** The physical information channels admit a strict ordering by marginal
contribution to false-association control, and that ordering is stable across network geometry and
clutter level. Registered: the ordering measured in Appendix J - polarization, then moment consistency, then kinematic
admissibility, then silence - reproduces with rank correlation
$\ge0.8$ across at least twenty (geometry, clutter) cells. The null is that the ordering is
geometry-dependent and therefore not a property of the physics. **Preliminary evidence from
Appendix M: the ordering polarization $>$ moment consistency $>$ kinematic admissibility $>$
silence reproduced identically in all eight cells tested, mean pairwise rank correlation
$+1.00$, below the twenty cells the hypothesis registers.**

**A2 (minimal sufficient set).** A proper subset of the channels is statistically indistinguishable
from the full stack. Registered: the C-III and C-IV pair alone attains a false-association rate
within $\delta=0.02$ of the full stack at non-inferior recall, in at least eighty per cent of
confirmatory cells. **Preliminary evidence from Appendix M: three-component polarization alone matched the full
stack in all eight cells tested, and adding the remaining constraints sometimes made results worse.
Appendix W shows that only the one-sided cells discriminate; in the others every channel succeeds,
so equality there is not evidence of sufficiency.** On real data the pair is further conditioned
by the per-station reliability of C-IV (Appendix T.1, T.2). Those
cells used a single event, a homogeneous medium with an exact velocity model, a perfect picker and
six to eight realizations each, so the correct reading is that polarization alone sufficed **in the
preliminary cells tested**. Whether it remains sufficient under heterogeneous media, model error,
realistic detection and multiple events is the confirmatory question, and it is registered as such
rather than settled. The operational statement being tested is that false-event rejection may not
require solving the wave equation.

Read with Appendix L, which found polarization model-dependent through ray geometry, this sharpens
the minimal set rather than contradicting it: **polarization supplies clutter immunity but not model
independence, kinematic admissibility supplies model independence but not discrimination**, and
neither alone suffices where Earth models are poor, which is the registered rationale for keeping
both.

**The minimal sufficient set is regime-dependent, and that is the sharpened form of A2.**
Appendix O measures the dependence against the coincidence-count ablation (C1, F1 below), and
Appendix P extends it to one station (C3). The primary question of this proposal therefore has a
graded answer rather than a single one - **the amount of physics required is set by how sparse the
network is** - and that ladder, not a fixed minimal set, is the object the confirmatory programme must
establish.

**A3 (model-independent estimation).** Under velocity-model error, prediction-based channels exhibit
indistinguishable degradation slopes in the tested gradient-medium experiment, and polarization inherits
model error through the ray geometry (Appendix L). Appendix U provides a positive instance of
model-independent estimation: an estimator using measured directions plus only a velocity bound shows
no systematic degradation trend across the tested model-error range, at a substantial cost in absolute
location accuracy. What is registered is sharper: **C-III is the only physical feasibility test in the
present constraint stack whose validity does not require a detailed velocity model**, because it
invokes only a lower velocity bound and the eikonal Lipschitz property, and its value must be measured
on a rejection endpoint, not a location endpoint, because a feasibility filter can remove hypotheses
without moving the optimizer. The robustness of this rejection channel under model error remains a
confirmatory endpoint.

**A4 (cost-effectiveness).** The registered robust pair C-III+C-IV is tested against the full constraint
stack on the false-association endpoint. Appendix N reports a cost of $0.46\%$ of the beam-plus-transport
waveform condition and $0.00003\%$ of a direct elastodynamic solve under the measured implementation.

### Family B. Observability, now that the certificate is established at leading order for the reduced model

Section 6.1 establishes the certificate for the reduced model and Appendix I verifies it at
$\rho=+0.9998$ in homogeneous and smooth-gradient media, so what follows are the parts that remain
open.

**B1 (interfaces): numerical evidence, narrowly, with the branch sum shown to be necessary.** Appendix Q tests a
two-layer model in which each phase arrives twice, as a direct ray and as a critically refracted
head wave. The eikonal identity $|\nabla_{\mathbf{x}_s}t|=1/V_1$ holds exactly for both branches,
including the head wave whose take-off angle is the critical angle regardless of distance. With the
branch sum the closed-form predictor attains $\rho=+0.963$ against the numerically computed
coercivity, just above the registered $0.95$; dropping the head wave and keeping only the direct
branch drops it to $+0.871$, so multipathing must be summed rather than ignored. The head wave
contributes a factor of $2.73$ to the median leading-order constant.

The predictor is nonetheless **markedly looser with interfaces**: ninety-five per cent of geometries
fall within $59\%$ of the prediction, against $11\%$ in the homogeneous and smooth-gradient cases,
partly because of the crude head-wave amplitude law of the test, which T1 separates from the theory. A
Spearman correlation of $0.963$ is numerical support for the branch-summed predictor in the tested
two-layer configuration, **not** a proof that the Proposition extends to interfaces.

**B2 (design beats the heuristic), false-event-rate endpoint.** The hypothesis
registered *median epicentral error* as its endpoint. Appendix N shows that endpoint cannot decide
it: with continuous optimization the location error saturates at zero for every design in the
regimes tested, so no design can be distinguished from another on it. The endpoint is therefore
changed to the **false-event rate**, where a signal does appear - at four stations, the
observability-optimal network produced no false events while the azimuthal-gap-optimal network
produced them in a fifth of realizations, with the observability constant differing between the two
designs by two orders of magnitude. Ten realizations per design is not enough to register that as a
result, so B2 is retained as open with the false-event-rate endpoint and a power-derived
replication count. This is an instance of the endpoint rule stated with the component hierarchy. Appendix I gives this favourable
prior odds and shows the mechanism: azimuthal gap misranks linear arrays and rings.

**B3 (observability predicts where clutter wins): preliminary evidence, one network.** Appendix J
found the observability constant is not a discriminating constraint when used as a feasibility
mask. Its proper role is different, and Appendix S supports it: moving a true source over twenty-five locations under one
fixed asymmetric network, across which the constant varies by a factor of eighty-six, the Spearman
correlation between the local constant and the local false-hypothesis margin is $+0.573$, above the
registered $0.5$. The lowest-constant quartile has a median margin of $0.93$ against $4.24$ in the
highest, a factor of $4.6$.

The operational statement is that **the parts of a monitoring region that will be vulnerable to
false events can be computed in advance from geometry alone**, before any data are recorded. The
role is therefore

$$\text{observability}\ \longrightarrow\ \text{identifiability and design map},\qquad\text{not}\qquad \text{observability}\ \longrightarrow\ \text{classifier},$$

which is the reading the measurements support and the one carried forward. That
completes the rehabilitation of C-I: it does not judge hypotheses, it maps where judging them will
be hard.

The correlation with the local false-event *rate* was only $+0.147$, but that endpoint carried
almost no variance here - the rate was zero in three quartiles and $0.06$ in the worst - which is the
mirror image of Appendix N and an instance of the endpoint rule stated with the component hierarchy.

### Family C. Where the advantage lives

**C1 (station-count scaling).** The advantage of the physical constraint stack over travel-time
association is a decreasing function of station count and crosses zero at some $S^\ast$. The null
is that the advantage is independent of station count; the registered deliverable is $S^\ast$ with
a confidence interval, not a claim that the advantage holds everywhere. **Numerical evidence:**
Appendix O finds the advantage over the coincidence-count ablation is $+0.88$ in false-event rate
at four stations, $+0.25$ at three and zero at six and twelve, so $S^\ast\approx5$-$6$ against that
comparator; five stations were not tested. Against a continuous kinematic score with exact picks
the four-station advantage disappears (Appendix Y.1), so $S^\ast$ is registered for re-measurement
against that score and with realistic pick errors.
At three stations every method, including the constrained prototype and both published
associators, emits a false event in 60 to 100 per cent of realizations (Appendix X.4); the
prototype's advantage there is recall, not a usable catalogue. With per-station reliability the
axis becomes the effective number of reliable stations $S_{\mathrm{eff}}$ of Section 6.4.

**C2 (geometric failure of travel-time association).** There is a class of network geometries -
those with low variance in epicentral distance - on which travel-time association has essentially
zero discrimination margin at any clutter level, and that class contains real deployed networks.
Appendix J measured a margin of exactly $0.00$ for a ring of eight equidistant stations while the
same score reached $2.60$ with distances spread from 6 to 62 km. Registered: identify at least
three real networks from public station metadata whose distance variance for a specified target
region places them in the failure class, and reproduce the collapse on their real geometries.
**Strengthened by Appendix M: the failure class is broader than equidistant rings. In a one-sided
network - stations confined to a half-plane, which describes coastal, border and single-flank
ocean-bottom deployments - the travel-time ablation produced a false event that tied or beat the
true source in every single realization, at both clutter levels and at every level of adversarial
decoy support, while polarization held the false-event rate at $0.00$ to $0.12$.** Two
qualifications follow from later measurements. GaMMA, which fits travel-time residuals rather than
counting coincidences, does not fail in that geometry (Appendix X.2), and part of the ablation's
failure is a tie effect of its discrete score (Appendix Y). The failure class is therefore a
property of coincidence counting under one-sided coverage, not of travel-time information as such.

**C3 (single-station multi-event): preliminary evidence for up to two sources, negative for
three.** Appendix P tests it directly. With one three-component station, hierarchical subtraction and a matched-filter
score, the two-branch representation recovers a single source to $1.6$ km and two overlapping
sources to $1.6$ and $1.0$ km, while the single-branch representation gives $22.0$ km on one source
and fails on one of two. At three sources both representations fail, with errors of $6$ to $22$ km
and four spurious detections. The registered capability boundary is therefore **two overlapping
events from one station**, with depth fixed and only epicentral position, origin time and source
count claimed. The registered hypothesis is that, from a single three-component station, in a
window with an unknown number of overlapping sources and structured clutter, the two-branch plus
polarization representation recovers source count and epicentral parameters at a rate the
single-branch representation cannot attain at any SNR. Depth is ill-conditioned by construction and
is reported with its Cramer-Rao bound only. On real data the hypothesis applies only to a station
that passes the polarization gate of Appendix T.1.

### Family D. Detection and association are one problem

**D1 (coupling): numerical evidence, with the threshold measured.** Appendix R parameterizes the picker
miss rate directly rather than simulating a detector, which makes the registered threshold a
function of a quantity real pickers report. The pick-free waveform channel recovers the source in
every realization at every miss rate, being by construction insensitive to missing picks, while
pick-based channels degrade from a miss rate of about $0.2$ and collapse beyond $0.6$: at a miss
rate of $0.8$, travel-time association recovers nothing while the pick-free channel recovers
everything. Polarization partially compensates, recovering a third of cases where travel time
recovers none, because it extracts more information from each surviving pick.

**This unifies with the sufficiency ladder.** What governs how much physics is required is not the
deployed station count but the **effective number of usable observations**, and that can be reduced
either by sparse deployment or by picker failure. Six stations at a miss rate of $0.8$ leave about
one station's worth of picks, which is precisely the regime in which Appendices O and P found the
waveform condition mandatory. The two axes are the same axis.

The registered form of D1 is the disciplined replacement for the inconclusive overlap study. Under a
realistic picker, the marginal value of pick-free waveform evidence is an increasing function of
the fraction of arrivals the picker misses, and exceeds the value of any association-stage
improvement once the miss rate passes a threshold. The registered deliverable is that threshold
with a confidence interval. The pre-study of Appendix K failed to settle this because its detector
reached only sixty-two per cent on an isolated event, so between-realization variance swamped the
coda effect; the confirmatory version requires a calibrated deep-learning picker, many
realizations, and a coda model fitted to real data.

### Family E. Registered nulls

These are stated as nulls because the evidence points against the original claim. Reporting them
as nulls is the alternative to quietly deleting hypotheses that did not survive.

**E1 (resolution).** Physical evidence does **not** improve resolution between nearby true sources
by more than fifteen per cent relative to arrival-time evidence at equal energy. Appendix H
measured six to fourteen per cent, flat in station count, so the fifteen per cent level is stated
as a null rather than as a target.

**E2 (velocity contrast).** The noise-limited resolution constant does not depend on $V_P/V_S$.
Appendix G measured it flat to $0.8\%$ over $[1.50,1.95]$.

**E4 (no graceful degradation in location).** No physical channel achieves materially slower
degradation of epicentral accuracy than travel-time association under structural velocity-model
error. Appendix L measured slopes within twenty-five per cent of one another over $0$ to $43\%$
error. The corresponding positive formulation of A3 uses a rejection endpoint rather than a location endpoint.

**E5 (second-order rank).** The number of sources is not recoverable from the rank of a
second-order matrix, in either the component or the station domain. Appendix V.1 measures this: the
inter-station cross-spectral matrix recovers $K=1$ with an eigenvalue gap of $2.1$ but returns
$1$ for $K=2$ and $1$ to $2$ for $K=3$ with gaps of $1.2$ to $1.6$, invariant to the temporal
separation of the events. This closes the signal-subspace family for the problem.

**E3 (overlap).** Overlapping events do not by themselves defeat association scoring in a
well-distributed network with complete picks. Appendix K found an empty ambiguity set across
separations from 3 to 24 km and offsets from 0.3 to 6 s. The corresponding positive claim is D1.

### Family F. The waveform branch, behind a gate

The following are retained but no longer carry the programme, because the evidence that motivated
them has weakened. They execute only if the Family A pilot shows that the minimal sufficient set
leaves a measurable gap.

**F1 (gate): numerical evidence, and the regime is identified.** Appendix O finds that a matched-filter waveform
condition reduces the false-event rate from $0.75$ to $0.00$ at three stations, far beyond the
registered $\delta=0.03$, while adding nothing at four stations or more, where the minimal set
already reaches zero. The waveform branch is therefore **not** closed; it is scoped. Its registered
regime is the sparsest networks, at or near the identifiability limit, and its budget is justified
only there. The matched filter of Appendix O is an oracle: exact velocity model, exact wavelet and
exact polarization. Where no velocity model supports waveform fitting at the working band, which is
the Colombian case of Section 16, this rung of the ladder cannot be tested on real data and is
claimed for synthetics and for the external domain only.

Sampling adequacy is decisive for this channel and is registered as a condition of the test: a
matched-filter correlation peak is about half a period wide, so origin time must be sampled at a
tenth of a period. Quarter-period sampling undersamples the peak and reverses the conclusion, as
Appendix O documents.

**F2 (misfit design).** Conditional on F1. In a $2\times2$ design crossing misfit type with
schedule, the main effect of optimal transport on the fraction of runs reaching the correct basin
is at least twenty percentage points under velocity-model error exceeding a quarter period.

**F3 (representation fidelity).** Conditional on F1. At viability level L1 the beam-plus-residual
representation meets the field, receiver and gradient tolerances of Appendix D and degrades
end-to-end F1 by at most $0.01$.

**F4 (symmetry-exact architecture).** Conditional on F1. Reciprocity and SO(3) covariance hold to
machine precision at random untrained weights, and the F3 tolerances are reached with at most half
the budget of a penalty-based baseline.

**F5 (compute break-even).** Conditional on F1. Break-even at $N_{\mathrm{break}}\le10^4$ source
evaluations, or a tenfold storage reduction at equal accuracy.

### Family G. The primary hypothesis, the methodological hypotheses and three secondary endpoints

**G1 (detection at a fixed realized false-event rate; the primary hypothesis and endpoint), restated.** In the regime where the picker
discards kinematic information (few stations, missed arrivals, pick timing error), a timing-tolerant
moment-tensor detector with deflation, which projects the migrated three-component waveform of all
stations onto the moment-tensor subspace allowed by the physics at the level of the fidelity ladder each
station supports, reads each station and phase within a tolerance around its predicted time and
explains detected events away, raises event recall by at least $0.05$ at a realized rate of $0.1$ and
$0.25$ false events per window, relative to its control with the geometry of another hypothesis, to the
best-tuned published associator and to the published migration locators, read on the exact FROC curve
with the threshold of Section 7.5. The members tested are the greedy group matching pursuit MT-FW,
whose selection step is the exact certified Frank-Wolfe oracle, and block orthogonal matching pursuit
[129] given the same peak search, BOMP-PK, which matches it (Appendix Y.21). The group BLASSO (16a)
solved on the grid falls below both, with or without the peak search (Appendices Y.20, Y.21), and the
statistic of the family, with its peak search and cap, is not the certificate of (16a) (Section 7.5,
Appendix AA.8). The original comparator, the coincidence-count ablation ABL, is retired as the
reference for G1 because its discrete score ties (Appendix Y.1).

**Status: numerical evidence on synthetics only.** The evidence, its qualifications and the four
positives withdrawn on the way are stated once, in Section 1.6, and in full in Appendices Y.1 to Y.21.
The effect is carried by the pick-free coherent detector with explaining-away, not by per-pick
polarization with realistic mechanisms (Y.10), by polarization as a filter on an associator's
catalogue (Y.7), or by a kinematic advantage that survives exact picks (Y.2, Y.3). Of the eight
preregistered synthetic tasks, four bear on G1 (Tasks 2A to 2D), one on G1-M (Task 2E), two on G1-R
(Tasks 2F and 2I) and one on G1-C (Task 2L). G1 holds in those four against the published comparators run with an operational
picker (Y.18 to Y.21), within the qualifications stated there, for the family of timing-tolerant moment-tensor
detectors with deflation, of which MT-FW is one member and not a better one than block orthogonal
matching pursuit given the same tolerance (Y.21); and as formulated it is not met at $0.1$ against GaMMA fed
optimistic picks with weak clutter, which is ahead in four of six cells (Y.18). The designs of Y.5 to
Y.10 followed Y.1 to Y.3, several rules were written after parts of their data had been seen, and the
exact-curve endpoint was adopted after a registered rule had failed on the coarse sweep; all of this is
declared in Appendix Y. Whether real pick errors, real media and real clutter place the RSNC in that
regime is the measurement the programme turns on. Appendix H supplies the mechanism at the margin
level; A1 and A2 decompose which physics earns it.

**G1-M (the methodological hypothesis).** The group BLASSO (16a) with a time shift per station and
phase as a nuisance parameter, bounded by the tolerance of the peak search, solved off the grid by
sliding Frank-Wolfe with a debiasing refit, raises event recall over BOMP-PK by at least $0.05$ at
$0.25$ false events per window, at equal tolerance and with the same station corrections, on the
windows of Task 2A and on thirty networks per cell. **Status: numerical evidence against it,
preregistered.** Task 2E registered it as rule 2E-a, with the registry hash communicated outside the
computing environment before any code was written, and ran it on the windows of Task 2A and thirty
networks per cell: the gain over BOMP-PK was $+0.000$ to $+0.039$ at $0.25$, every hierarchical interval
contained zero and the rule failed in both regimes (Appendix Y.22). An independent verification found no
implementation error and no budget handicap. G1-M is refuted in this synthetic design. Joint estimation
of station time shifts and a moment tensor for one known event is prior art [143], [144]; G1-M concerned
the shifts as nuisance parameters of a detection certificate over an unknown number of events. The null
is reported; G1 is unaffected, because it concerns the family, and the methodological contribution rests
on G1-R.

**G1-R (the bounded-influence hypothesis).** In windows with three to five overlapping events, a
moment-tensor subspace statistic that bounds the influence of the peak vector of each station and phase
before or while projecting it on the moment-tensor subspace, and deflates sequentially, raises event
recall over BOMP-PK by at least $0.03$ at $0.25$ false events per window with an interval excluding zero,
in two fresh seeds, and the gain is due to the bounding. Its first form is the statistic of MT-FW, which
clips each peak vector at a fixed norm; its second replaces the clipping by a Huber fit whose level is
fixed by rule from the noise distribution. **Status: numerical evidence, preregistered twice (Task 2F,
Appendix Y.24; Task 2I, Appendix Y.26).** In Task 2F the rule held in both clutter regimes and both seeds, $+0.043$ to $+0.074$,
also under Holm's correction; removing the cap brought MT-FW to BOMP-PK within $0.002$; and in standard
windows the two were equivalent within $0.05$. An adversarial verification narrowed its scope: the gain
is $+0.024$ to $+0.045$ when origin times must also match, it belongs to the cap level fixed in Task 1x
and reverses at low levels, and its mechanism, exploratory, is bounded influence of each station and
phase block under amplitude model error and a better calibration of the score across windows rather
than a better separation of events. The amplitude factors of the generator are redrawn in every window
and are not station terms; removing them reduces the gain, and at a small spread of them the cap still
helps where a statistic given the true factors does not (Appendix Y.25). The hypothesis was selected among exploratory candidates on the calibration seed (Appendix Y.23),
its seeds were untouched by that selection, and per-station clipping is known in robust stacking. An
exploration on the calibration seed does not support amplitude or site terms per station as its
generalization (Appendix Y.25). Task 2I tested the form that the exploration did support, HUBER-N, the
best of 11 candidates there, on two further fresh seeds with 120 networks per cell and regime, with its
registry hash sent before any run and its code frozen. Its three rules held in both seeds: $+0.058$
[$+0.039$, $+0.076$] and $+0.053$ [$+0.037$, $+0.070$] over BOMP-PK in dense windows; $-0.008$ and
$-0.002$ against MT-FW with the cap, with lower bounds above the margin of $-0.03$; and equivalence with
BOMP-PK within $0.05$ in standard windows. Its verification confirmed the result and stated its scope:
$+0.031$ and $+0.019$ with origin times matched, equivalence on isolated events not met with weak
clutter in one seed, no more events found with every atom accepted, slightly below the cap at $0.1$ in
one seed, and about nine times the cost of MT-FW. The secondary arms of Task 2I, descriptive and on 20
networks per cell and regime, add two things. The result is robust to the level: Huber levels of 2.5 and
3.4 noise standard deviations give the recall of the registered one, within $0.006$ in dense windows,
and a level of 5 gives less. It is not robust to a large amplitude spread: at a spread of 0.8 the Huber
fit still exceeds BOMP-PK, by $+0.100$ and $+0.072$, but is behind the cap by $0.055$ and $0.069$, with
an interval that excludes zero in one seed, while at a spread of 0.25 neither form gains reliably. G1-R
therefore rests on two preregistered, independently verified synthetic tests on one generator family,
the second of which replaces the calibrated cap by a level fixed by rule from the noise distribution at
the nominal amplitude spread; where the amplitude error of the model is large the calibrated cap is the
better form, and the size of the effect is modest (Section 1.6). A note on what G1-R does not address:
two explorations place the remaining headroom in the accuracy of the station terms and in the amplitude
of weak arrivals, not in the statistic that G1-R varies (Appendix Y.27); the first of these is the
subject of G1-C below. G1-R is a weaker
criterion than G1, and the difference is deliberate. G1 asks for at least $0.05$ at both $0.1$ and
$0.25$ false events per window against methods of another kind, where exploratory differences were
several times that size. G1-R asks for at least $0.03$ at $0.25$ only, because it compares two members
of one family that differ in the bound alone, its pilot effect was about $0.07$ and its sample size was
planned on half of that, and the reading at $0.1$ was registered as descriptive. Its equivalence margin
in standard windows, $0.05$, is larger than its own superiority threshold, so a cost on isolated events
of the size of the gain in dense ones is not excluded by the registered rule; tasks registered from
here on use a margin of $0.02$ per regime (Section 14).

**G1-C (the calibration hypothesis).** With the detector held fixed, the capped moment-tensor statistic
of MT-FW with a local refinement of each detection by the timing likelihood, station terms estimated
without picks from the twenty catalogued events raise event recall over the registered catalogue
corrections at $0.25$ false events per window when a hit must also match the origin time within 0.5 s:
by at least $0.03$ with an interval excluding zero in standard windows, with an interval excluding zero
in dense windows, and without loss, within a margin of $0.02$, in either clutter regime or under the
location-only hit rule, in two fresh seeds. The estimator is a joint determination of the hypocentres of
the catalogued events and of the terms, with the catalogue as a Gaussian prior, whose data are
matched-filter energy peaks read again as hypocentres and terms move. **Status: numerical evidence,
preregistered once (Task 2L, Appendix Y.28).** The candidate came from an exploration on the calibration
seed whose plan, with its primary contrast, was frozen with a hash before any detection result; the
pilot estimate was $+0.070$ [$+0.010$, $+0.120$] and the registration planned on $+0.05$. Task 2L ran on
two seeds never used before, with 160 networks per cell and regime in the standard design and 100 in the
dense one, its registry hash sent before any run and its code frozen. Its four rules held in both seeds:
$+0.057$ [$+0.034$, $+0.081$] and $+0.062$ [$+0.036$, $+0.083$] in standard windows; $+0.038$ [$+0.024$,
$+0.054$] and $+0.030$ [$+0.015$, $+0.044$] in dense windows; lower bounds per regime between $+0.007$
and $+0.052$; and $+0.036$ and $+0.050$ in standard and $+0.020$ and $+0.035$ in dense windows under the
location-only rule, for which only non-inferiority was registered. The registered secondary, the same
contrast in block pursuit, gave $+0.068$ [$+0.026$, $+0.107$] and $+0.071$ [$+0.029$, $+0.114$]. An
independent verification found no leakage, the terms recomputed without access to the true delays or the
true catalogue being identical in 300 of 300 units, reproduced the rules with its own code, found the
gain positive at every origin-time window tested from 0.2 to 1.0 s, and attributed it to the terms and
not to the timing level of the refinement. It also narrowed the scope. The pooled gain is essentially
that of the one-sided eight-station cell, $+0.108$ and $+0.113$ in standard windows; with six stations
it is $+0.048$ and $+0.027$, the second with an interval containing zero, and with four stations no gain
is demonstrated in either seed. The estimated terms are worse than the registered ones in about 22 per
cent of the networks, and in six-station networks of that kind recall falls in standard windows. The
gain is conditional on the refinement: Task 2L has no arm with the new terms and the hard window, and in
the exploration the terms gave nothing measurable without it. The estimator assumes a catalogue of
twenty events with errors of known size, 1 km and 0.1 s, and in the exploration a prior of 2 km returned
the term error to that of the registered corrections. It captures 0.28 to 0.61 of the gain of the true
delays, which remain ahead. The evidence is synthetic and from one generator family. Joint hypocentre
determination with station terms is classical in location from picks; what G1-C measures is the
estimation of the terms without picks and their effect inside a pick-free detector at a fixed
false-event rate, and it is a calibration result, not a detector. G1-C is the first hypothesis registered
with the origin-time-matched hit rule as primary and with non-inferiority per regime at a margin of
$0.02$ (Section 14); it has one preregistered test where G1-R has two, and the two have not been run
together.

G2 to G4 are secondary endpoints: they characterize the margin, the silence channel and uncertainty
calibration, and none of them can make the thesis succeed or fail.

**G2 (clutter scaling).** As the false-arrival rate rises, the discrimination margin of the
constrained model degrades at less than half the rate of the travel-time model, and no false
hypothesis exceeds fifty per cent of the true score.

**G4 (silence channel, secondary).** Enabling the certified exclusion channel of Section 6.2 reduces the
false-source rate by at least $0.03$ at recall non-inferior within $\delta_{\mathrm{rec}}=0.02$,
and never removes a true source outside the declared erosion tolerance. Appendix J found this
**conditional on station distance spread**: in a network whose stations sit at similar distances
the channel is inert or harmful, so distance spread is a registered factor rather than an
assumption.

**G3 (uncertainty).** Coverage error at the nominal $90\%$ level is $\le0.05$ and strictly smaller
than the travel-time model at comparable sharpness. The parametric likelihood and a
simulation-based comparator are both reported.

## 4.4 General objective

To test whether timing-tolerant moment-tensor detection with deflation, applied to the
three-component waveform before picking and thresholded at a false-alarm rate calibrated on noise-only
windows, detects more events at a fixed realized false-event rate than published associators and
migration locators where the picker discards kinematic information (G1), on synthetics and on a
decimated national network; to establish which element of the moment-tensor detection statistic carries a
methodological advance over block pursuit with equal tolerance, where the off-grid group BLASSO with
time shifts did not (G1-M) and a bounded-influence statistic did in dense sequences (G1-R), with a
calibrated clipping level and with a Huber level fixed by rule from the noise distribution, and to
establish how far station terms estimated without picks raise the recall of the refined detector
(G1-C); and, secondarily, to formulate, analyze and implement the constrained inverse framework in which the
detection statistic is an optimality certificate, with the feasible set
generated from elastodynamic structure and with identifiability, optimization behaviour,
computational cost and uncertainty calibration treated as scientific quantities.

## 4.5 Specific objectives

Objective 1 carries the primary hypothesis; objectives 2 to 12 support it, and objectives 6 and 8
are conditional on gate F1 (Section 4.3, Family F).

1. Build the bounded-influence moment-tensor subspace statistic as a method: the per-station clipped
   statistic of MT-FW, whose effect in dense sequences was confirmed in Task 2F (G1-R, Appendix Y.24), and
   the Huber fit on the residual of each station and phase block with its level fixed by rule from the
   noise distribution, confirmed in Task 2I as not inferior to the capped statistic (Appendix Y.26);
   extend that evidence beyond its limits, to a second generator family, a held-out operating point, hits
   with the origin time matched and a lower cost of the Huber fit; keep as references BOMP-PK
   with and without the cap and the grid and off-grid convex solvers of Tasks 2C to
   2E, which do not beat BOMP-PK at equal tolerance (Appendices Y.20 to Y.22); study the measured headroom,
   about three quarters of which remains unexplained by any treatment of amplitudes (Appendix Y.25),
   starting from the lever that two explorations pointed to and that Task 2L confirmed within a narrow
   scope, station terms estimated without picks (G1-C, Appendices Y.27, Y.28), and extending that result
   to catalogue errors of unknown size, to the geometries in which it was not demonstrated and to its
   combination with the bounded-influence statistic, read against the oracle of the true delays and with
   hits that also match the origin time as the primary hit rule; and test G1 with the best
   timing-tolerant detector on synthetics and on the RSNC by network decimation (Section 16), against published associators, migration locators,
   block-sparse moment-tensor inversion [129] with and without that tolerance, and mechanism-corrected
   stacking [122], [133].
2. Derive and implement the five constraint families of Section 6, each as an explicit component of
   the feasible set, and characterize when each is active.
3. Establish the observability-to-coercivity translation and use it as a well-posedness certificate
   and as a network-design criterion.
4. Formulate silence as a certified exclusion set under an explicit sensing and noise model, motivated
   by unique continuation but not derived from it (Sections 3.3, 6.2), and quantify its tolerance under
   model error.
5. Formulate the additive multi-source observation model over a sparse source measure with unknown
   cardinality.
6. Conditional on F1: design and validate the misfit stack, optimal transport with multiscale
   frequency continuation, with a whitened Student-$t$ final stage and a simulation-based comparator.
7. Adopt and validate an established elastodynamic solver as physical reference and build a
   reciprocity-based strain-Green's-tensor reference.
8. Conditional on F1: build the differentiable representation as analytic Gaussian beams plus a
   symmetry-exact learned residual, with pre-registered field and gradient tolerances and hard abort
   triggers; it is an option, not a contribution of the thesis.
9. Design constrained optimization and sensitivity algorithms, including adjoint and analytic
   gradients with independent cross-validation.
10. Quantify per-constraint contribution, representation fidelity and compute/storage break-even.
11. Develop and evaluate single-station multi-event epicenter determination.
12. Establish leakage-controlled transfer protocols on at least two contrasting real test domains.

---

# 5. Physical Foundations

## 5.1 Governing equations

With $\mathbf{u}$ particle displacement, $\rho$ density, $\boldsymbol{\sigma}$ Cauchy stress and
$\mathbf{f}$ the body-force source,

$$\rho\,\partial_t^2\mathbf{u}=\nabla\cdot\boldsymbol{\sigma}+\mathbf{f},\qquad
\boldsymbol{\sigma}=\lambda(\nabla\cdot\mathbf{u})\mathbf{I}+2\mu\,\boldsymbol{\epsilon},\qquad
\boldsymbol{\epsilon}=\tfrac12\big(\nabla\mathbf{u}+(\nabla\mathbf{u})^{\mathsf T}\big), \tag{1}$$

giving the heterogeneous elastodynamic system

$$\rho\,\mathbf{u}_{tt}=\nabla[\lambda(\nabla\cdot\mathbf{u})]
+\nabla\cdot\{\mu[\nabla\mathbf{u}+(\nabla\mathbf{u})^{\mathsf T}]\}+\mathbf{f},
\qquad V_P=\sqrt{\tfrac{\lambda+2\mu}{\rho}},\quad V_S=\sqrt{\tfrac{\mu}{\rho}}. \tag{2}$$

P and S arise from the elastic system; they are not imposed as independent scalar fields [1], [2].

## 5.2 Source, initial and boundary conditions

$$f_i(\mathbf{x},t)=-\sum_j M_{ij}s(t-\tau_s)\,\partial_{x_j}\delta(\mathbf{x}-\mathbf{x}_s), \tag{3}$$

with $s_k$ a relative-time, normalized source-time function satisfying $\int s_k(t)\,dt=1$; the physical moment-rate history is shifted by $\tau_k$ and constrained as in Section 6.5. The normalization removes scale ambiguity between source strength and source-time shape. For numerics the
delta is mollified with width $h_s$; convergence as $h_s\to0$ is an acceptance test. Initial
conditions are quiescent; the top boundary is free surface and truncated boundaries use PML.

## 5.3 Observation operator

$$\mathbf{d}_s(t)=\mathcal{O}_s[\mathbf{u}](t)+\boldsymbol{\varepsilon}_s(t), \tag{4}$$

with $\mathcal{O}_s$ carrying receiver location, component projection, sampling and a declared
instrument-response approximation. In multi-source windows the prediction is generated from the
total field, never from isolated candidate windows, which is what prevents overlapping energy from
being credited independently to several sources. A first-arrival operator is defined as a derived
diagnostic map only.

## 5.4 Linearity, the Green's tensor, and what must be computed rather than learned

For the linear system (2),

$$u_n(\mathbf{r},t)=\sum_{ij}M_{ij}\,H_{nij}(\mathbf{r},\mathbf{x}_s;\,t-\tau_s),
\qquad H_{nij}=\partial_{x_j}G_{ni}*s, \tag{5}$$

with $G$ the elastodynamic Green's tensor. Three consequences are used throughout. Origin time is
an exact temporal shift. The moment tensor enters linearly, so $\partial u/\partial\mathbf{M}$ is
exact and estimating $\mathbf{M}$ at fixed $(\mathbf{x},\tau)$ is a linear subproblem. And by
source-receiver reciprocity the strain Green's tensor for all source positions follows from $3S$
simulations with forces at the $S$ receivers [37], [38], [39], independently of the number of
candidates. A learned representation is therefore never asked to model the $\tau$ or $\mathbf{M}$
dependence, and the reciprocity reference is cheap in compute but expensive in storage, which is
where a continuous representation has a measurable value proposition (Section 10.6).

## 5.5 Helmholtz decomposition, with its scope declared

$$\mathbf{u}=\nabla\phi+\nabla\times\boldsymbol{\psi},\qquad \nabla\cdot\boldsymbol{\psi}=0, \tag{6}$$

which for a locally homogeneous region separates (2) into scalar problems at $V_P$ and $V_S$. In a
smoothly heterogeneous medium the potentials do **not** decouple; the correct statement is a
coupled system

$$\ddot\phi=V_P^2\nabla^2\phi+f_\phi+\mathcal{K}_\phi[\phi,\boldsymbol{\psi};\nabla\lambda,\nabla\mu],
\qquad
\ddot{\boldsymbol{\psi}}=V_S^2\nabla^2\boldsymbol{\psi}+\mathbf{f}_\psi+\mathcal{K}_\psi[\phi,\boldsymbol{\psi};\nabla\mu], \tag{7}$$

with $\mathcal{K}$ vanishing only where the elastic parameters are locally constant. Two further
caveats are declared: on a bounded domain with a free surface the decomposition is not unique
without decay conditions, and Rayleigh waves are precisely the coupled P-SV motion that no such
split isolates. The registered scope is therefore: the validated core uses piecewise-homogeneous
media where the coupling is confined to interfaces and appears as mode conversion in the boundary
conditions; for smoothly heterogeneous media $\mathcal{K}$ is treated as a perturbation whose
magnitude is **measured** on the synthetic generator and reported; and the residual is always
evaluated on the reconstructed vector field, never on the potentials alone, so that any coupling
neglected in synthesis reappears in the misfit rather than being silently discarded. Numerical and
learned elastic separation based on (6) is established [42] and is used here as a modelling device,
not a post-processing filter.

Within that scope the decomposition earns its place for one reason above all: a single event
produces two coherent arrivals whose separation $\Delta t_{PS}=d(1/V_S-1/V_P)$ depends on distance
alone and is independent of origin time. This is what separates origin time from location and is
the structural fact behind the single-station result of Section 9.

---

# 6. Five Sources of Physical Constraint

This section develops the constraint framework behind the secondary hypotheses; G1, G1-M, G1-R and G1-C
rest on the detector and the framework of Section 7.5. Each subsection takes a statement about elastic wave propagation
and turns it into an explicit component of the feasible set of the optimization problem of
Section 7. Each is independently ablatable, so that A1 and A2 - which constraints actually carry the
attribution information - is answerable rather than rhetorical.

Throughout, $\mathcal{E}=\{(\mathbf{x}_k,\tau_k,\mathbf M_k,s_k)\}_{k=1}^{K}$ denotes the candidate set, where $\mathbf M_k$ is the physical moment-tensor coefficient and $s_k$ is a normalized source-time function. The derived activity is $a_k=\|\mathbf M_k\|_F/M_{\rm ref}\ge0$, with $M_{\rm ref}$ a declared reference scale; activity thresholding therefore defines existence without introducing an independent activation variable. $\mathcal{C}$ is an *a priori* fixed partition of station-time cells, constructed from network geometry and the window alone, independent of the data and of $\mathcal{E}$.

## 6.1 Constraint I: observability certificates and coercivity

**Physical statement.** Whether a network can determine sources in a region is the observability
question for the wave equation. Control theory answers it geometrically: the wave equation is
observable from a region if and only if every ray of geometrical optics reaches that region within
the observation time - the geometric control condition of Bardos, Lebeau and Rauch [43], with the
broader control framework of [45] and inverse-problem context in [44]. The condition holds in
heterogeneous media and requires no translation invariance, which is exactly why it survives where
super-resolution theory does not.

**Translation into optimization.** Observability is an inequality
$\|\mathbf{u}_0\|^2\le C\int_0^T\!\!\int_{\Gamma}|\mathcal{O}[\mathbf{u}]|^2$, i.e. a lower bound on
how much of the state the measurements see. Linearizing the source-to-data map $F$ about a candidate
and writing $J$ for its Jacobian, the same inequality delivers a lower bound on the smallest
singular value of $J$ restricted to the admissible source manifold, hence a **coercivity constant**
$\gamma>0$ for the reduced Hessian $J^{\mathsf T}\Sigma^{-1}J$ on that manifold. Three deliverables
follow:

$$\gamma(\mathcal{S},T,\Omega)\;=\;\sigma_{\min}^2\!\left(\Sigma^{-1/2}J\big|_{\mathcal{T}}\right)
\;\ge\;C^{-1}(\text{GCC data}), \tag{8}$$

1. **A local identifiability certificate, and only that.** Where $\gamma>0$ the Gauss-Newton
   curvature $J^{\mathsf T}\Sigma^{-1}J$ is positive definite on the admissible tangent space, which
   certifies local identifiability and bounds the linearized estimator variance. Strong convexity of
   the fully nonlinear constrained objective, existence and uniqueness of the constrained minimizer,
   and convergence rates for the primal-dual iteration are **separate results requiring additional
   hypotheses** - second-order regularity, control of the residual Hessian term, a constraint
   qualification, and conditions on the admissible set - and are registered as work rather than
   asserted here [74], [75]. Where $\gamma\to0$ the proposal reports non-identifiability rather than
   a number with a fabricated error bar.
2. **A network-design criterion.** Maximizing $\gamma$ over candidate station placements is a
   design problem with a geometric interpretation: place stations so that rays from the target
   volume are captured. This converts "where should the next station go" from heuristic to
   optimization, and it is a deliverable of independent operational value.
3. **A constraint in the solver.** The trust region and the merit-function parameters are set from
   $\gamma$ rather than tuned, and the source volume is restricted to the sub-region where
   $\gamma$ exceeds a registered floor.

**Proposition (homogeneous half-space, high-frequency limit, known mechanism).** Let the source
parameters be $\theta=(\mathbf{x},\tau)\in\mathbb{R}^4$, let $g$ be the band-limited source pulse,
and for station $s$ and branch $b\in\{P,S\}$ write $t^s_b=\tau+d_s/V_b$ and

$$\mathbf{v}^s_b\;=\;\partial_\theta t^s_b\;=\;\big(-\mathbf{n}_s/V_b,\;1\big)\in\mathbb{R}^4,
\qquad w^s_b=\frac{A_b^2}{d_s^2}. \tag{8a}$$

Let $\mathcal{A}(T)=\{(s,b):\ \mathrm{supp}\,g(\cdot-t^s_b)\subset[0,T]\}$ be the set of arrivals
captured by the recording window. Because the P polarization is longitudinal and the SV
polarization orthogonal to it, the two branch contributions are orthogonal in the three-component
space, and the whitened Fisher matrix of the reduced problem satisfies

$$\gamma\;=\;\frac{\|\dot g\|^2_{[0,T]}}{\sigma^2}\;
\lambda_{\min}\!\Big(\!\!\sum_{(s,b)\in\mathcal{A}(T)}\!\! w^s_b\,\mathbf{v}^s_b(\mathbf{v}^s_b)^{\mathsf T}\Big)
\;+\;O\!\big(V/(f_{\max}d_{\min})\big). \tag{8b}$$

The leading term is computable in closed form from station directions, distances and velocities
alone: no wave simulation, no Jacobian, no data. Three corollaries follow, and each is one half or
one ingredient of the geometric control condition.

**Corollary 1 (time).** If $\mathcal{A}(T)=\emptyset$ the leading term vanishes identically and
$\gamma$ is of subleading order: the problem is not observable until at least one ray has reached a
receiver within $T$. Moreover the leading term is non-decreasing in $T$ and constant once
$T\ge\max_s d_s/V_S+\tau$, so the useful window length is predictable from geometry alone.

**Corollary 2 (single-branch cone degeneracy).** With a single branch, $\mathbf{u}=(\mathbf{a},c)$
annihilates the leading term iff $\mathbf{a}\cdot\mathbf{n}_s=cV$ for every captured $s$, i.e. iff
all station directions lie on a common cone about $\mathbf{a}$. The great-circle case $c=0$ is the
planar or linear array; the equal-take-off-angle case is a ring of stations at equal epicentral
distance. **The reduced single-branch model therefore becomes rank-deficient for planar and equal-distance
idealized geometries**, and what fails is the separation of depth from origin time. Whether such
degeneracies are prevalent in deployed networks is an empirical question, evaluated on real network
geometries under C2 rather than asserted here.

**Corollary 3 (branch diversity).** With both branches at a station,
$\mathbf{v}^s_P-\mathbf{v}^s_S=(\mathbf{n}_s\kappa,0)$ with $\kappa=1/V_S-1/V_P>0$, so the span
contains pure-space directions and the origin time is eliminated: this is
$\Delta t_{PS}$ independence of $\tau$ in algebraic form. The leading term is then positive iff
$\{\mathbf{n}_s\}_{s\in\mathcal{A}}$ spans $\mathbb{R}^3$. The following separable lower bound is
used as a **design proxy**; its derivation beyond the tested reduced setting is part of T1.

$$\gamma\;\ge\;\frac{\|\dot g\|^2}{\sigma^2}\;S_{\mathcal A}\;\mu_{\min}(\mathcal A)\;\min_s\lambda_{\min}(Q_s), \tag{8c}$$

where $Q_s$ is the $2\times2$ weighted Gram matrix of the two branch vectors after restriction to their two-dimensional branch span, and $\mu_{\min}(\mathcal A)$ is the smallest eigenvalue of the direction covariance $\frac{1}{S_\mathcal{A}}\sum_s\mathbf{n}_s\mathbf{n}_s^{\mathsf T}$. Because (8c) is a separable design proxy rather than the exact coercivity formula, it is used to
rank geometries, and the exact leading-order quantity in (8b) remains the primary predictor. Its
three factors are the number of captured arrivals, the directional coverage and the branch
diversity.

Note that the velocity contrast enters (8c) through $\kappa^2$, improving the worst-case
conditioning. This does **not** contradict the withdrawn conjecture of Appendix G, which concerned
the resolution constant at equal signal energy; (8c) is a statement about rank and worst-case
conditioning, not about resolution.

**Corollary 4 (heterogeneous media).** The Proposition is not restricted to a homogeneous medium.
By Fermat's principle and the eikonal equation, the derivative of the two-point travel time with
respect to the source position is minus the slowness vector at the **source end** of the connecting
ray, with $|\mathbf{p}|=1/V(\mathbf{x}_s)$ and direction along the ray tangent there, independently
of the path taken. Hence (8a) generalizes verbatim as

$$\mathbf{v}^s_b=\big(-\mathbf{p}^s_b,\;1\big),\qquad \big|\mathbf{p}^s_b\big|=\frac{1}{V_b(\mathbf{x}_s)}, \tag{8d}$$

with $\mathbf{p}^s_b$ the source-side slowness vector of the ray of branch $b$ reaching station $s$,
and everything above carries over. Two consequences are not cosmetic. The velocity entering the
certificate is the **local** velocity at the source, not a path average. And the direction is the
ray **take-off** direction, which in a medium with a velocity gradient differs substantially from
the straight line to the station. The geometric coverage quantity must therefore be computed on
take-off directions obtained by ray tracing, which is what distinguishes this certificate from a
geometric-dilution-of-precision heuristic rather than merely renaming one.

**Status: Proposition at leading order.** The leading term of (8b) and Corollaries 1-4 follow
from the derivation above; the remainder $O(V/(f_{\max}d_{\min}))$ is stated, not bounded, and its
explicit bound is part of T1. The leading-order statement is verified numerically in both the
homogeneous and the smooth-gradient heterogeneous case. A pre-study reported in Appendix I tested
the full chain on 96 network geometries spanning six structural families. A purely geometric
ray-coverage quantity, computed from station directions alone without ever forming a Jacobian or
touching a waveform, predicts the coercivity constant of the full physical operator with Spearman
$\rho=+0.82$, and predicts the Monte-Carlo epicentral error directly with $\rho=-0.91$; the
classical azimuthal-gap heuristic reaches only $+0.45$ and $-0.66$ on the same data, and station
count $+0.39$ and $-0.48$. These correlations are larger in magnitude than $0.5$, the level that B3
registers for a different endpoint, the local false-hypothesis margin tested in Appendix S, so they are
not a test of B3; what they show is that the control-theoretic quantity outperforms the standard
seismological heuristic on the same data.

The time half of the condition was tested separately and behaves as the theory requires: truncating
the recording window before any ray arrives leaves the Fisher matrix numerically singular, the
first captured arrival raises $\gamma$ by nine orders of magnitude, and the constant saturates
exactly once the slowest arrival of the farthest station is inside the window, with no further gain
afterwards. No purely spatial heuristic can see this behaviour.

The closed-form leading term of (8b) reproduces the numerically computed coercivity with Spearman
$\rho=+0.9998$ over 120 non-degenerate geometries, and after removing a single geometry-independent
discretization constant, $95\%$ of geometries fall within $11\%$ of the closed-form prediction with
median ratio $1.0000$. Corollaries 2 and 3 hold to machine precision: for a ring of equidistant
stations the single-branch leading term is $O(10^{-16})$ while the two-branch term is $5\times
10^{-1}$, and for a linear array every branch set gives $O(10^{-14})$, so that array is observable
only through the subleading polarization term. The bound (8c) held in $87$ of $87$ cases where it is
positive, but is loose by a median factor of $36$: it certifies feasibility, it does not predict
conditioning. For prediction, the exact leading term of (8b) is the object to use.

In a medium with a constant velocity gradient, for which the two-point travel time is known in
closed form, the take-off-slowness form (8d) reproduces the numerically computed coercivity with
$\rho\approx0.999$ and median deviation below $1\%$, while the naive straight-line version - the
one a dilution-of-precision heuristic would use - degrades to between $40\%$ and $90\%$ median
deviation as soon as the gradient is non-zero. The eikonal identity
$|\nabla_{\mathbf{x}_s}T|=1/V(\mathbf{x}_s)$ was verified to a relative error of $10^{-11}$,
confirming the closed-form travel times and their derivatives.

**Remaining scope.** The pre-studies use a homogeneous half-space or a smooth constant
gradient, direct arrivals only, a single source, no clutter, and local estimation started near the
truth, so it verifies local
well-posedness rather than global attribution. The Proposition is asymptotic in $f_{\max}d_{\min}/V$; near the degeneracy manifold the leading
term vanishes and the subleading polarization contribution dominates, which is exactly the regime
where the closed-form prediction loses accuracy and where C3's statement about ill-conditioned
depth applies. What remains for T1 is a medium with interfaces, where the two-point ray is not
unique and multipathing, triplications and head waves must be handled; the certificate then becomes
a sum over ray branches and its degeneracy analysis must account for them. In a heterogeneous medium the coverage quantity
requires ray tracing, which is the genuine test of the condition and is work package T1. The
falsifiable content of B1 and B3 therefore stands for the heterogeneous case; what Appendix I removes is
the possibility that the certificate was decorative from the start.

## 6.2 Constraint II: certified non-detection constraints under an explicit sensing model

*An extension, not a load-bearing component: its value is conditional on station distance spread
(Appendix J), and the classical unique-continuation theorem does not apply to a point sensor. What
is constructed is an exclusion set certified against a declared sensing and noise model.*

**Physical statement, with its logical status declared.** Solutions of the wave equation obey
unique continuation: vanishing on an **open space-time set** forces vanishing in the corresponding
domain of dependence, under conditions on the coefficients [44]. That theorem is the *motivation*
for this constraint and **not its derivation**, and the distinction matters enough to state before
anything is built on it. A seismic station is a point sensor recording noise-limited data over a
finite interval; it does not deliver vanishing on an open space-time set, and no theorem of unique
continuation applies to it directly. The literature is explicit that observation on space-time
regions, boundary observation and point observation are different problems, and there is dedicated
work on unique continuation for the discretized wave equation with noisy data precisely because the
passage from continuum theory to discrete noisy sensors cannot be made automatically.

What can be constructed instead is a **certified exclusion set under explicit observation
assumptions**, which is what this section does. The assumptions are declared, not implied:

1. a declared observation operator and instrument response;
2. a declared temporal bandwidth and window;
3. a calibrated detector threshold with a stated false-negative rate;
4. the station geometry;
5. a finite, characterized noise floor;
6. an admissible coefficient class for the medium, namely the piecewise-smooth class of
   Section 5.5;
7. a minimum observation time consistent with the geometric control condition of Section 6.1.

Under those seven conditions, and only under them, a station-time cell classified as silent by the pre-registered detector
over the declared interval is incompatible with any source whose predicted radiation would exceed
the floor there. The exclusion is therefore certified relative to a declared sensing model, not
derived from a PDE theorem.

**Translation into optimization.** Let $\mathcal{C}_0(\mathbf d)\subset\mathcal{C}$ be the cells classified as silent by a pre-registered detector applied to the observations. The partition $\mathcal{C}$ is data-independent; only the silent-cell label is data-derived. Each such cell defines an **exclusion set**

$$\mathcal{X}_c=\Big\{(\mathbf{x},\tau,\mathbf M,s)\ :\ \big\|\mathcal{O}_{s_c}\big[\mathbf{u}(\mathbf{x},\tau,\mathbf M,s)\big]\big|_{\mathcal{T}_c}
\;>\;\theta_c\Big\}, \qquad
(\mathbf{x}_k,\tau_k,\mathbf M_k,s_k)\notin\bigcup_{c\in\mathcal{C}_0(\mathbf d)}\mathcal{X}_c\ \ \text{whenever}\ a_k>0, \tag{9}$$

with $\theta_c$ a declared amplitude floor. Under the seven assumptions above this is a feasibility
constraint on the source measure rather than a likelihood factor, which is a change of status
relative to existing monitoring systems that treat silence probabilistically [17] or ignore it
[23]. It is not, and must not be presented as, a consequence of unique continuation alone.

**The honesty conditions, stated in advance.** Unique continuation is a statement about exact
solutions of the exact equation. Three things break it in practice, and each is handled explicitly.
Amplitude prediction requires attenuation and site response that the elastic model does not
contain, so $\theta_c$ is calibrated **empirically** per station from continuous data with
template-confirmed events, in the spirit of probability-based network completeness [24], and not
derived from the elastic forward model. Model error blurs the boundary of $\mathcal{X}_c$, so the
exclusion set is **eroded** by a registered tolerance $\eta$ before use, and the erosion is
calibrated so that no true source in the pilot is excluded. And the coefficient regularity
conditions of unique continuation are not satisfied by arbitrary media, so the theoretical
statement is claimed only for the piecewise-smooth class declared in Section 5.5, with the general
case handled by the eroded empirical set. G4 tests the value of the resulting constraint and its
failure mode is measured, not assumed away.

**A relaxation ladder.** The constraint is implemented in three registered variants so that its
aggressiveness is a controlled factor: hard exclusion, soft exclusion with slack and penalty, and
the calibrated-Bernoulli likelihood of the previous formulation, which becomes the softest rung and
the ablation baseline.

## 6.3 Constraint III: medium-independent kinematic admissibility

**Physical statement.** Betti-type reciprocity identities among source-receiver pairs are not the
right instrument here: with a single event there is no second elastodynamic state, so those
identities are not testable from the data this problem has. Reciprocity enters the proposal in its
proper place, as an exact symmetry of the learned representation (Section 10.4). The medium-independent content that *is* available comes from the eikonal
equation. Since $|\nabla T|=1/V\le1/V_{\min}$ wherever the velocity is bounded below, the travel-time
field is Lipschitz with constant $1/V_{\min}$ **in any medium**, so two arrivals of the same phase
belonging to the same event must satisfy

$$\big|t_i-t_j\big|\;\le\;\frac{\|\mathbf{r}_i-\mathbf{r}_j\|}{V_{\min}}, \tag{9a}$$

and the two branches at one station must satisfy $t_S>t_P$ whenever both are paired to the same
source and branch family under the declared model. These require **no velocity model, no source
location and no ray tracing** - only a lower bound on velocity - which is exactly the property a
test must have to transfer between regions whose Earth models differ in quality.

The monotonicity of $\Delta t_{PS}=t_S-t_P$ with $t_P$ is deliberately **not** part of C-III. It
holds in a homogeneous or smoothly graded medium with a single
branch family, and fails in exactly the setting Appendix Q constructs, where interfaces, head waves
and multiple ray families break the pairing. It is retained only as a model-specific diagnostic,
declared as such, and $\Delta t_{PS}$ is used as a locating observable only where the branch pairing
has been identified.

**Translation into optimization.** (9a) and the ordering condition enter the feasible set as the
**inequality family** $g_{\mathrm{kin}}(\mathcal{E},\mathbf{d})\le0$ of (15), not as equalities:
they are **feasibility constraints on pick subsets**, and they act at a different level from the other families: they prune
which *combinations* of arrivals could be one event, before any location hypothesis is posited,
whereas C-IV and C-V score a hypothesis once posited. Reciprocity survives in the proposal in its
proper place, as an exact symmetry imposed architecturally on the learned representation
(Section 10.4), rather than as a data-level identity it cannot support.

**Measured.** In the pre-study of Appendix K these constraints reject $0.0\%$ of uncontaminated synthetic truths, as they must under the declared assumptions, and reject $91$-$100\%$ of pick sets assembled from random
arrivals. In the realistic case of a true event with some of its picks replaced by clutter, the
rejection rate is $43\%$ for one substituted pick, $76\%$ for two and $97\%$ for four. C-III is
therefore the only family that rejects no uncontaminated truth under its declared assumptions, and
the cheapest to evaluate.

## 6.4 Constraint IV: polarization rank

**Physical statement, correctly bounded.** Over a short window the three-component particle motion
at a station has a $3\times3$ covariance whose spectrum is interpretable: a dominant single
eigenvalue indicates rectilinear motion, two comparable eigenvalues indicate elliptical motion from
a surface wave or from two superposed arrivals, and three comparable eigenvalues indicate noise.

It is tempting, and wrong, to read this as a source counter. Two sources whose polarization vectors
at a station are collinear produce an effective rank of one; three sources spanning a plane produce
two; and noise raises the effective rank to three regardless of the source count. The correct
statement is therefore a **bound**, not an estimate:

$$\operatorname{rank}_{\mathrm{eff}}\Sigma_s(t)\ \text{bounds the dimension of the observed
polarization subspace, and nothing more.} \tag{9b}$$

It does **not** bound the number of sources, in either direction, and Appendix V.1 shows that
moving to the inter-station cross-spectral matrix does not repair this: cardinality is not
recoverable from the rank of a second-order matrix in either domain. Five sources whose station
polarizations happen to lie near a common plane produce an effective rank of two, so
$K>\operatorname{rank}_{\mathrm{eff}}\Sigma$ is entirely possible; and noise inflates the rank
independently of $K$. Deriving any statement about $K$ from the rank requires an additional,
explicitly verified incoherence condition on the source polarizations, which this proposal does not
assume and does not use.

The effective rank is computed from an eigenvalue-ratio threshold declared in advance rather than
from an exact rank test, since noise makes exact rank meaningless. What the channel supplies
without picks, without a network and without travel times is therefore a dimension bound on the
polarization subspace, and it is used as such.

**Translation into optimization.** Let $\Sigma_{s}(t)$ be the observed local covariance and
$\widehat\Sigma_s(t;\mathcal{E})$ the covariance predicted by the active set. A single constraint
form is used, and it is a **subspace-consistency** condition, not a cardinality bound:

$$\operatorname{rank}(P_{\mathrm{pred}})\le r_{\mathrm{obs}},\qquad
\left\|(I-P_{\mathrm{obs}})P_{\mathrm{pred}}\right\|_F\le\epsilon, \tag{9c}$$

where $P_{\mathrm{pred}}$ and $P_{\mathrm{obs}}$ are orthogonal projectors onto the dominant predicted
and observed polarization subspaces, and $r_{\mathrm{obs}}$ is the declared effective rank. This
projector condition remains well defined when the two subspaces have different dimensions and encodes
that the predicted subspace must fit inside, and be aligned with, the observed one; **no inference
about the number of sources is drawn from the rank.** This is the polarization channel that the pre-study of Appendix H
identified as the discriminant clutter cannot imitate, since randomly polarized clutter has no
consistent dominant direction across stations. Cardinality is handled where it belongs, by the group/measure sparsity term and the derived activity thresholding of
Section 7, and Appendix V.1 records that the obvious repair - moving the rank test to the
inter-station cross-spectral matrix - does not work either.

**Role in the present formulation.** Polarization is a primary channel because the pre-studies identify it as
the dominant discriminator under the tested clutter regimes. Its effective-rank representation makes it an
explicit constraint rather than a misfit term.
The associated cost is declared: it requires accurate inter-component relative amplitudes and
correct instrument orientation, both measured and reported before the constraint is enabled.

**Measured on real data, the cost is per station.** Appendix T.1 runs this check on five stations
of the Colombian network: single-station P polarization recovers the back-azimuth of catalogued
events within $20^\circ$ in $0.79$ of cases at one station and falls short of the gate at the
others: at the chance rate at CHI, where orientation errors, channel reversals, mis-association and
mid-period instrument changes were excluded, not distinguishable from chance at OCA, and
undecided on too few picks at SJC. C-IV is therefore enabled station by station, behind a
reliability gate measured on catalogued events, and its weight in the stack is a measured
property of the station rather than a constant. How a failing station is treated matters:
Appendix T.2 finds that switching C-IV off there does not repair it, because a hypothesis supported
only by failing stations then escapes the polarization test, whereas giving their arrivals a
neutral value measured on catalogued events recovers, with one station in four reliable, the
advantage over the coincidence-count ablation. Appendix Y shows that advantage is itself partly a
tie effect of the ablation's score, so what T.2 establishes robustly is the ordering of the
treatments of a failing station, not the size of the advantage.

**Per-station likelihood, registered for the confirmatory programme.** The neutral value of T.2 is a
rule; the formulation that makes it exact is a likelihood ratio whose parameters are the two
real-data measurements of Appendix T.1. Let $\delta\in[0^\circ,90^\circ]$ be the folded angle between
the measured polarization of an arrival and the direction predicted by a candidate. Under $H_1$,
the arrival belongs to the candidate, and at station $s$

$$p_{1,s}(\delta)=r_s\,p_{\mathrm{vM}}(\delta;\kappa_s)+(1-r_s)\,p_{0,s}(\delta), \tag{9d}$$

where $p_{\mathrm{vM}}$ is a folded von Mises density with concentration $\kappa_s$ and $p_{0,s}$ is
the geometric null density of the station, the density whose mass below $20^\circ$ is the chance
rate $f_{0,s}$ of T.1. Under $H_0$, the arrival is unrelated to the candidate and has density
$p_{0,s}$, replaced by $\eta\,p_{\mathrm{vM}}+(1-\eta)\,p_{0,s}$ if clutter is found to be partly
coherent. The reliability $r_s$ is a mixture weight and is estimated from the gate exactly as $\eta$
is estimated from clutter,

$$\hat r_s=\frac{f_{\mathrm{ctrl},s}-f_{0,s}}{q_s-f_{0,s}},\qquad
q_s=\int_0^{20^\circ}p_{\mathrm{vM}}(\delta;\kappa_s)\,d\delta, \tag{9e}$$

and the polarization contribution of a candidate $h$ with matched arrivals $\mathcal M(h)$ is
$\ell_{\mathrm{pol}}(h)=\sum_{i\in\mathcal M(h)}\log\big[p_{1,s_i}(\delta_i)/p_{0,s_i}(\delta_i)\big]$.
Four properties follow without further rules: a station with $r_s=0$ contributes exactly zero, so
the neutrality of T.2 is a consequence rather than a design; the score is a sum, so it cannot reward
candidates that explain less (the defect of Section 13.0b); it is continuous, so it breaks ties with
evidence (Appendix Y); and the binary gate becomes a continuous weight, with
$S_{\mathrm{eff}}=\sum_s\hat r_s$ as the effective number of reliable stations on which the
sufficiency ladder is re-indexed. **Relation to prior art.** A per-station, per-phase azimuth
residual distribution is already part of NETVISA [17]; what is added here is the estimation of its
reliability weight from catalogued events with an exact chance correction, its use as a hard
subspace test alongside C-III, and its evaluation on the catalogue endpoint in sparse networks. The
confirmatory programme uses (9d)-(9e) in place of the multiplicative mean-cosine score of the
prototype.

## 6.5 Constraint V: source-spectrum and source-time-function bounds

*A regularizer that shapes the feasible set, not an independent scientific claim: the
isotropic-fraction test on the moment tensor did not discriminate (Appendix V.3), and the useful
content is the admissible cone of the source-time function.*

**Physical statement, scoped to the declared source parameterization.** These bounds are properties
of the source model adopted in Section 5.2, namely the separable form
$M_{ij,k}(t)=M_{ij,k}\,s_k(t-\tau_k)$ with a scalar source-time function $s_k$ defined in relative time and a fixed mechanism, and they are
**not** presented as universal laws holding for an arbitrary time-dependent moment tensor, for which
componentwise non-negativity is not even well defined. Within the declared class, each $s_k$ is non-negative, causal, normalized by $\int s_k(t)\,dt=1$, and of finite duration bounded by the source-moment scale, with a low-frequency
plateau, a corner frequency tied to source size and high-frequency decay. Causality imposes analytic constraints on the complex spectrum; Hilbert-transform amplitude--phase relations require additional minimum-phase assumptions and are therefore not imposed as a universal hard constraint here.
Extending the class - to a time-varying mechanism, for instance - requires re-deriving the
admissible cone rather than reusing this one.

**Translation into optimization.** These become **bound and cone constraints** on the source
parameters. Define $M_{0,k}:=\|\mathbf M_k\|_F$ as the source-moment scale used by the declared source parameterization.

$$s_k(t)\ge0,\qquad \int s_k(t)\,dt=1,\qquad \mathrm{supp}\,s_k\subset[0,T_{\max}(M_{0,k})],\qquad
\big|\widehat{s}_k(\omega)\big|\le \Phi(\omega;M_{0,k},f_{c,k}), \tag{10}$$

together with the declared mechanism restrictions on $\mathbf M_k$. These restrictions are not treated as a primary clutter discriminator: Appendix V.3 finds that the isotropic-fraction test on recovered tensors is not selective. The C-V contribution instead acts through the admissible source-time-function cone; under the declared source parameterization, fitting a false hypothesis can force $s_k$ outside this cone. This is the most direct point of contact with the supervisor's published line on physically bounded inverse formulations [72], [73].

## 6.6 Why constraints enlarge the discrimination margin

The pre-study in Appendix J evaluates the families individually against a true source embedded in
clutter. The results show that the five families are **not co-equal, and one of them is not a
primary discrimination constraint.**

| Family | Effect on discrimination margin | Status |
|:--|:--|:--|
| C-IV polarization | large and consistent on the margin; removing it can invert the margin | **primary discriminator on the margin** where the station passes the gate (T.1); at the catalogue endpoint it pays only inside the coherent moment-tensor detector applied to the waveform before picking (Section 1.6) |
| C-V source bounds | small, consistent in sign across clutter levels and geometries | secondary discriminator |
| C-II silence | positive with distance spread, inert or harmful in a uniform-distance ring | conditional on geometry |
| C-I observability | inert to slightly harmful for discrimination | **well-posedness certificate, not a discriminator** |
| C-III kinematic admissibility | zero rejection of uncontaminated truths; rejects 43-97% of clutter-contaminated pick sets | **combinatorial pre-filter, medium-independent** |

The role of C-I is *conditioning*: it indicates whether the problem is locally solvable and supports
network design, rather than directly deciding whether a particular hypothesis is false. The C-II result is also informative: in a network whose stations are
all at similar distance, every hypothesis predicts similar amplitudes everywhere and silence
carries almost no information; the constraint acquires value only when distances vary, which is a
condition on deployment that the proposal must state rather than assume.

Taken together, the families provide complementary physical information at different stages of inference. A travel-time associator
imposes one scalar coincidence test per predicted arrival. The constrained formulation imposes,
for the same hypothesis, a coercivity requirement on the geometry, an exclusion from every silent
cell, medium-independent kinematic inequalities, a rank and subspace condition on particle motion, and
bound constraints on the source time function. Clutter passes the first easily. Passing all of them
simultaneously is unlikely for clutter that is not itself a physically realizable source, which is the
hardest case for every channel (Appendix T).

That is the qualitative explanation of the quantitative result measured in Appendix H: a fourfold
larger discrimination margin, which lost 29 per cent of its value over the clutter range against 66 per
cent for the arrival-time margin, a relative degradation about 2.3 times slower. Appendix J then
locates *which* part of the physics does the work, and the answer is mostly polarization, which is
consistent with Appendix H rather than an additional assumption. On real data that answer holds
where the station passes the polarization gate of Appendix T.1, which was one of five stations
tested; the hierarchy is a statement about usable polarization, and how much of a real network
supplies it is an open measurement (Section 17). And the hierarchy was measured against a
coincidence count; against a continuous kinematic score its catalogue-level value is the one stated
in Section 1.6. It is also why the
experimental design of Section 13 is built around leave-one-out ablation of constraint families
rather than around a single global comparison: a total improvement without per-constraint
attribution would not be evidence of scientific understanding.

---

# 7. The Inverse Problem

## 7.1 Data and unknowns

An association window contains $N$ records $\mathbf{y}_i=(\mathbf{w}_i,\mathbf{r}_i,s_i)$ with
$\mathbf{w}_i$ a three-component waveform window, $\mathbf{r}_i$ receiver location and $s_i$ station
identity. Picker timestamps, where available, are used only for indexing and windowing, never as the
governing observable. The unknown is a tensor-valued discrete measure over source space,

$$\mu=\sum_{k=1}^{K}\delta_{(\mathbf{x}_k,\tau_k)}\otimes\mathbf M_k, \tag{11}$$

with $K$ an over-complete candidate budget and source-time functions $s_k$ satisfying the normalized
admissible class of C-V. Define the derived activity $a_k=\|\mathbf M_k\|_F/M_{\rm ref}\ge0$;
there is no separate multiplicative activation variable, so the inverse problem has no activation--tensor
scale degeneracy. Binary existence is recovered by a calibrated threshold on $a_k$ after convergence.

The merge metric is the same normalized source-space metric used for set evaluation,
$d_{\mathrm{src}}^2=\|\Delta\mathbf x\|^2/\sigma_x^2+\Delta\tau^2/\sigma_\tau^2$. The merge
threshold is calibrated once in the pilot and frozen before confirmatory testing. After convergence, candidates
within that threshold are clustered for set counting; their **moment-rate histories** are summed for reporting,
rather than forcing two distinct source-time shapes into a single separable parameterization. The resulting
cluster is counted as one recovered source.

## 7.2 Additive observation model

Because the constitutive model is linear in the source moment tensor, the record is the **sum** of
source contributions plus background, not a mixture over exclusive assignments:

$$\mathbf{w}_i=\Phi\Big[\mathcal{O}_{s_i}\Big(\sum_k \mathbf{u}_k(\mathbf{x}_k,\tau_k,\mathbf M_k,s_k)\Big)\Big]+\mathbf{b}_i+\boldsymbol{\varepsilon}_i,
\qquad
\mathbf{r}_i(\mathcal{E})=W_i\big(\mathbf{w}_i-\widehat{\mathbf{w}}_i(\mathcal{E})\big), \tag{12}$$

where $\mathbf u_k$ is the response to moment tensor $\mathbf M_k$ and normalized source-time
function $s_k$. $W_i$ is a whitening operator estimated only on calibration data. There is therefore
no latent exclusive assignment and no responsibility variable in the primary model; a normalized energy
share $\gamma_{ik}$ is computed post hoc as a diagnostic so that pick-level metrics remain comparable
with mixture-based baselines. The exclusive-assignment mixture is retained solely as the
"no-superposition" ablation and as the GaMMA-comparable baseline; the two models are never mixed in
one objective.

## 7.3 Misfit design and cycle skipping

A quadratic-type misfit in whitened residual space is non-convex in $(\mathbf{x}_k,\tau_k)$ through
phase: a half-period displacement gives a residual larger than predicting zero, which at $2$ Hz is
$0.25$ s, about $0.9$ km of position error or $2\%$ velocity error at regional distance. Since the
medium is approximate by construction, a naive time-domain waveform misfit measures velocity-model
error rather than attribution evidence [29], [30]. The misfit is therefore designed:

*The waveform channel has measured value at three stations (Appendix O); that optimal transport is
the reason for that value is **not** established and is registered as F2. Transport is a candidate
mechanism for optimizing the channel, not the explanation of why the channel works.*

1. **Primary: optimal transport, with its regime declared.** The comparison is a **Kantorovich**
   problem. Because seismograms are signed, each trace is lifted to two non-negative channels by the
   Jordan decomposition $u=u_+-u_-$, with $u_+=\max(u,0)$ and $u_-=\max(-u,0)$. Predicted and observed
   lifted traces define non-negative measures $\alpha$ and $\beta$ on the declared space-time--component
   grid. The formulation is unbalanced: the transport plan $\pi$ may create or destroy mass at a declared
   KL cost rather than being forced to satisfy exact marginal constraints. No deterministic Monge map is
   claimed or required. Graph-space transport or the lifted construction is used according to the
   pre-registered stage [31], [32], [33]. Following the entropic regularization principle of Cuturi [34],
   the unbalanced problem solved is

   $$\pi^\star_\varepsilon=\arg\min_{\pi\ge0}\left\{\langle\pi,C\rangle-\varepsilon H(\pi)\right.
   \left.+\lambda_u D_{\mathrm{KL}}(\pi\mathbf{1}\Vert\alpha)+\lambda_u D_{\mathrm{KL}}(\pi^{\mathsf T}\mathbf{1}\Vert\beta)\right\}, \tag{13}$$

   with $C$ the declared ground cost on the space-time--component grid, $H$ the entropy and $\lambda_u$
   the registered marginal-relaxation weight. The same unbalanced transport functional is used in the
   three terms of the debiased divergence. Regularization is
   registered rather than tuned: $\varepsilon$ is set to $10^{-2}$ of the median squared ground cost
   at each continuation stage, so that it scales with the band, and annealed by a factor of two per
   stage. Two consequences are declared. The plan $\pi^\star_\varepsilon$ carries an **entropic
   bias** relative to exact optimal transport, which pushes the misfit toward blurred couplings, so
   the reported quantity is the **debiased Sinkhorn divergence**, defined from the same unbalanced
   transport functional as $S_\varepsilon(\alpha,\beta)=W_\varepsilon(\alpha,\beta)-\tfrac12W_\varepsilon(\alpha,\alpha)-\tfrac12W_\varepsilon(\beta,\beta)$,
   which vanishes at $\alpha=\beta$ whereas the raw entropic cost does not. And sensitivity to
   $\varepsilon$ is an explicit ablation over $\{10^{-3},10^{-2},10^{-1}\}$ of the median cost,
   because an $\varepsilon$ chosen to make a method look good is a known failure mode. The stage misfit $\mathcal J_{\mathrm{OT}}$ is the sum of these debiased station/component divergences under the declared normalization and weighting.
   The Sinkhorn formulation follows the computational construction of [34], while implementation uses the stabilized log-domain scaling of [35]; the Gibbs kernel
   $\exp(-C/\varepsilon)$ is never exponentiated directly. Sliced or quantile transport,
   $O(n\log n)$, is the registered fallback [36].
2. **Multiscale continuation.** A registered band ladder $\mathcal{F}_1=[0.1,0.4]$,
   $\mathcal{F}_2=[0.1,0.8]$, $\mathcal{F}_3=[0.1,1.4]$, $\mathcal{F}_4=[0.1,2.0]$ Hz, each stage
   initialized from the previous [30]. Ladder membership is a controlled factor.
3. **Final stage and probabilistic model.** A whitened multivariate Student-$t$ on
   $\mathbf{r}_i(\mathcal{E})$ at the last stage supplies the likelihood used for calibrated
   uncertainty; running it from a cold start is the ablation testing F2.
4. **Branch weights, fitted.** Appendix G measured that at equal energy and $S\ge2$ the joint P+S
   operator is $10$-$19\%$ worse in resolution than S alone, because energy in the faster branch
   buys less positional sensitivity. Branch weights are therefore fitted on calibration data, not
   set to unity, and the equal-weight variant is ablated.

**Probabilistic status.** $\mathcal{J}_{\mathrm{OT}}$ is a misfit, not a log-likelihood. The
registered default is that all calibration statements come from the final-stage Student-$t$ model,
with the transport stages treated as an initialization mechanism; a Gibbs posterior
$\pi_\beta\propto\exp(-\beta\mathcal{J}_{\mathrm{OT}})p$ with temperature calibrated on the
calibration partition is reported as a secondary analysis. Given documented evidence that
parametric waveform likelihoods are overconfident under real noise and theory error [27], [28], a
simulation-based likelihood is a **registered comparator**, and model error is represented as a
separate covariance contribution estimated from the difference between generating and inference
physics rather than absorbed into observation noise.

## 7.4 The constrained problem

At continuation stage $m$,

$$\min_{\mathcal{E},\ \xi\ge0}\quad
\mathcal{J}^{(m)}_{\mathrm{mis}}(\mathcal{E})\;+\;\lambda_M\sum_k \|\mathbf M_k\|_F/M_{\rm ref}\;+\;\lambda_\Theta R_\Theta\;+\;J_{\mathrm{pen}}(\xi) \tag{14}$$

subject to, for all $k$ with $a_k>0$ and all cells $c$:

$$\underbrace{\gamma(\mathcal{S},T,\Omega_k)\ge\gamma_{\min}}_{\text{C-I observability}}\quad
\underbrace{(\mathbf{x}_k,\tau_k)\notin\mathcal{X}_c^{\,\eta}}_{\text{C-II exclusion}}\quad
\underbrace{g_{\mathrm{kin}}(\mathcal{E},\mathbf{d})\le0}_{\text{C-III kinematic admissibility}}\quad
\underbrace{\left\|(I-P_{\mathrm{obs}})P_{\mathrm{pred}}\right\|_F\le\epsilon+\xi_c}_{\text{C-IV polarization subspace}}\quad
\underbrace{(\mathbf M_k,s_k)\in\mathcal{S}_{\mathrm{adm}}}_{\text{C-V source bounds}} \tag{15}$$

together with the derived activity $a_k=\|\mathbf M_k\|_F/M_{\rm ref}$ and an energy-admissibility budget per cell,

$$\sum_k\big\|\mathcal{O}_{s_c}(\mathbf{u}_k)\big\|^2_{W,\mathcal{T}_c}\le(1+\varrho)\big\|\mathbf{w}_{s_c,\mathcal{T}_c}\big\|^2_W+\xi_c, \tag{16}$$

which is convex in the moment-tensor coefficients $\{\mathbf M_k\}$ for fixed positions, origin times and source-time functions because the forward map is linear in $\mathbf M_k$. Feasibility of the truth is not assumed: under
superposition the cross terms can be negative, so (15) can be violated by the true solution under
destructive interference. The constraint is therefore **soft by default**, with $\xi_c$ always
admissible and penalized, and $\varrho$ set from the measured distribution of the coherence
statistic so that the truth is feasible in at least $95\%$ of cells; the hard variant is retained
only as an ablation, precisely so that the recall cost of hardening is measured rather than assumed
to be zero. The sparsity term is the convex relaxation of cardinality and the discrete counterpart
of the total-variation norm of the measure (11); the source-set prior is an inhomogeneous Poisson
point process with declared seismicity bounds [25].

## 7.5 Operational form: the group BLASSO, its dual certificate, and a certified Frank-Wolfe

The measure (11) and the sparsity term of (14) define, once the misfit is the whitened least-squares
term of the final stage, a convex problem over tensor-valued measures,

$$\hat\mu\in\arg\min_{\mu}\ \tfrac12\|\mathbf d-\Phi\mu\|^2_{\Sigma^{-1}}+\lambda\,|\mu|_{\mathrm{TV},G},\qquad
|\mu|_{\mathrm{TV},G}=\sum_k\big\|\Sigma^{-1/2}G(\theta_k)\mathbf m_k\big\|_2, \tag{16a}$$

the BLASSO of [113], [114], [46] with amplitudes in $\mathbb R^6$ and each atom weighted by the energy its
moment tensor produces in the whitened network. Its optimality condition is a certificate: with
$\mathbf r=\mathbf d-\Phi\hat\mu$ and $P_\theta$ the orthogonal projector onto the range of
$\Sigma^{-1/2}G(\theta)$, the function $\eta(\theta)=\|P_\theta\Sigma^{-1/2}\mathbf r\|/\lambda$ satisfies
$\eta\le1$ on the hypothesis space and $\eta=1$ on the support of $\hat\mu$ (Proposition AA.1). These
identities hold for minimizers of (16a); they need not hold for the output of the greedy iteration
piloted in Appendix Y.11. At $\mu=0$, $\lambda^2\eta^2$ is
the generalized likelihood ratio with the moment tensor as a linear nuisance, the matched-subspace
detector of [118], [137] applied jointly to all stations and to P and S, and up to normalisation the
variance reduction on which grid-based moment-tensor monitoring detects [132], [123]. Four consequences organize the
method.

1. **Mechanism invariance without a mechanism search.** The unknown moment tensor enters linearly and
   is solved in closed form at every hypothesis; the detector demands consistency of signs and
   relative amplitudes across stations, which no scalar characteristic function retains.
2. **Counting, locating and separating in one iteration.** The total-variation term is the convex
   relaxation of the number of events [115]. The sliding Frank-Wolfe algorithm
   [116], [117] adds at each step the global maximizer of $\eta$, refits the amplitudes and moves the
   positions; subtracting the fitted contribution explains an event away, so sidelobes and
   overlapping events are handled by the fit rather than by suppressing neighbours. The pilot of
   Appendix Y.11, labelled MT-FW, is a simpler relative: a greedy group matching pursuit, close to the
   block orthogonal matching pursuit of [129], whose selection step is the Frank-Wolfe oracle, which fits the new atom by unpenalized least squares
   (its amplitude is not shrunk), never refits earlier atoms or moves positions, caps the
   per-station vectors, and excludes from reselection every node within 5 km epicentral distance of
   an accepted atom. That exclusion was fixed before evaluation and is much smaller than the 20 km
   suppression that the single-pass WAVE-MT needs to reach a comparable recall. With the cap active,
   a deflation step does not fully remove an event (Proposition AA.4(g)), a plausible but untested
   reason why some exclusion is needed.
3. **A certified oracle.** On the grid, Frank-Wolfe with an exact linear oracle brings the objective of
   (16a) to within $2\|\mathbf d\|^4_{\Sigma^{-1}}/(\lambda^2(k+2))$ of its minimum after $k\ge1$ steps,
   the $O(1/k)$ rate of [124], and an oracle exact to within $\varepsilon'$ in $\lambda\eta$ adds a term
   proportional to $\varepsilon'$ (Proposition AA.4(a)-(c)). The branch-and-bound of Appendix Z
   supplies the exact oracle: no projector increases a norm, so the energy of the three-component
   records in a time window widened by $h/V_{\min}+\delta\tau$ bounds $\eta$ over a whole cell, and a
   cell is discarded only when that bound is below the current maximum (Proposition AA.3, Corollary
   AA.1). The bound is valid but loose: measured on the pilot it avoids about 4 per cent of the exact
   evaluations, so the certified oracle is implemented as a batched evaluation with a floating-point
   enclosure, and its value is exactness rather than speed (Appendix Z). The rate holds for Frank-Wolfe and its fully corrective and sliding variants, not for the
   greedy deflation of the pilot, for which no rate is available and, with the cap active, not even a
   monotone decrease of the residual (Proposition AA.4(e)-(h)); the duality gap of Proposition
   AA.4(d), evaluated after a refit on the grid dictionary, is the registered measure of how far its
   output is from a minimizer.
4. **A threshold with a prescribed false-alarm rate.** In a window without events the detector
   declares an event exactly when the maximum $M$ of its statistic over the hypothesis space exceeds
   the threshold, so the per-window false-alarm rate is $P(M>u)$ (Proposition AA.2(c) for the BLASSO;
   for the pilot, because its first atom is the global maximum). The primary method sets $u$ from
   noise-only windows recorded by the same network in the same noise regime: split-conformal
   calibration, transferred from conformal outlier testing [139], gives $P(M>u)\le\alpha$ in finite samples whenever calibration and operational windows
   are exchangeable (Proposition AA.5(c)), and a generalized extreme-value fit to the calibration
   maxima [125] extrapolates below $1/n_{\mathrm{cal}}$. Random-field theory gives only conservative
   prior thresholds. Under whitened Gaussian noise the idealized certificate $\lambda^2\eta^2$ at
   $\mu=0$ is $\chi^2_6$ at each hypothesis, and it is a $\chi^2_6$ random field, whose maximum has a
   tail approximated by the expected Euler characteristic of its excursion sets [119], [120], [127],
   only under a symmetry assumption on the moment-tensor subspaces that is not expected to hold for
   them (assumption (G) of Proposition AA.5(b)); the implemented statistic, with its peak search
   within $\pm0.25$ s and its per-station cap, is not a chi-square field, nor $\chi^2_6$ at a single
   hypothesis. Measured on 300 held-out Gaussian-noise windows in each of six settings (Appendix
   Y.17), the realized false-alarm rates at $\alpha=0.01$, $0.05$, $0.10$ and $0.25$ were $0.007$,
   $0.037$, $0.095$ and $0.239$ with split-conformal thresholds, and the expected Euler characteristic
   overestimates $P(M>u)$ by a factor of $20$ to $60$ because the field decorrelates within one grid
   step. Strong clutter made of visible transients breaks every threshold calibrated on Gaussian noise
   (realized $0.60$, $0.78$, $0.90$ and $0.96$ with split-conformal), and recalibration on clutter
   windows with the same clutter law restores calibration; the extreme-value and Bonferroni figures are
   in Appendix Y.17. The calibration is therefore repeated whenever the noise or clutter regime
   changes, and the expected-Euler-characteristic and Bonferroni thresholds serve only as conservative
   priors before noise-only windows are available.

**The convex problem is the framework, not the best detector on a coarse grid.** Task 2C (Appendix
Y.20) solved (16a) restricted to the grid of Appendix Y.18, 675 nodes at 3 km spacing and 31 origin
times at 0.1 s steps, on the matched-filter dictionary with exact predicted times and the
block-whitened group norm of (16a), by fully corrective Frank-Wolfe, with relative duality gaps of at
most $1.4\times10^{-8}$ and the optimality conditions of Proposition AA.1 verified at every solved
point of the path, scoring each block by the largest $\lambda$ of a 30-point path at which it enters
the support. On the 1200 windows of Task 2A its recall stayed below that of MT-FW by $0.06$ to $0.61$
in 23 of 24 combinations of cell, seed and false-event cap, and equal to it in the last, and at or
slightly below that of block orthogonal matching pursuit on the same dictionary, whose selection
statistic is its certificate. The certified minimizer therefore does not deliver the performance of
MT-FW. What does are ingredients that (16a) on a grid does not contain, above all the peak search
within $\pm0.25$ s, which an exploratory diagnosis identifies as what absorbs the timing errors of a
3 km, 0.1 s grid and of the residual station corrections: with sources on the nodes the exact-time
block pursuit is the best method. The group BLASSO remains the framework in which the detection
statistic is an optimality certificate, the source of that certificate and of its likelihood ratio at
zero (Proposition AA.2), and the object of Appendix AA; on a coarse grid with timing errors the
detector to deploy is a timing-tolerant greedy pursuit, MT-FW or, equivalently on these synthetics,
block orthogonal matching pursuit with the same peak search.

**The peak search, not MT-FW, carries the advantage.** Task 2D (Appendix Y.21) gave the antecedents
and the grid problem the peak search of MT-FW. Block orthogonal matching pursuit in the form of [129],
reading each block at the energy peaks of the current residual within $\pm0.25$ s, matched MT-FW in
every cell, with differences of -0.155 to +0.091 at 0.25 false events per window and every
hierarchical interval containing zero; removing the peak search from MT-FW cost 0.06 to 0.41, removing
its cap changed little, and without the cap MT-FW is essentially that pursuit, the two differing only
by the P and S terms of the Gram matrix and by deflation in place of an orthogonal refit. The grid
problem on a dictionary grown from atoms read at those peaks, GL-FW-PK, stayed 0.11 to 0.30 below
MT-FW in 11 of 12 cases, and its certificate covers only the dictionary that the heuristic peak
oracle builds. The element that matters is therefore tolerance to timing error, and at present it is
a heuristic: an argmax per station and phase within a fixed window, with no model of the error, no
penalty on the shifts and no certificate over them. Its principled form makes the shifts part of the
problem: each atom carries a time shift $\delta_{i,b}$ per station $i$ and phase $b$, bounded by a
declared tolerance or penalized, the certificate becomes the maximum of $\eta$ over positions, origin
times and shifts, and the sliding step of consequence 2 moves all three continuously, which removes
the discretization error of the grid and models the residual timing error of the corrections and of
the velocity model, with a debiasing refit on the selected support against the shrinkage of the
penalized fit. Appendix AA.8 notes that within the search window the implemented statistic is a lower
bound on a shift-maximized likelihood ratio. Task 2E made that shift-maximized problem the object of
the optimization, with analytic gradients, a convex step solved to a certified duality gap and a
calibrated penalty, and ran it against BOMP-PK on the windows of Task 2A and thirty networks per cell
(Appendix Y.22): it gained $+0.000$ to $+0.039$ at 0.25 false events per window, every hierarchical
interval containing zero, so the methodological hypothesis G1-M is refuted in this design. Started from
a grid node, the local refinement cycle-skips by about one wavelet period at some stations; at the
calibrated value the group penalty is inert; and the events that block pursuit misses are absorbed by
atoms mislocated 5 to 10 km away, which local sliding does not move. Station time shifts estimated
jointly with the moment tensor of one known event are prior art in source inversion [143], [144].

**A bounded-influence statistic: what the cap does.** Task 2D found the per-station cap immaterial in standard
windows, and Task 2F agrees there; in windows with three to five overlapping events it is the whole
difference between MT-FW and BOMP-PK (Appendix Y.24). The capped statistic replaces each peak vector
$\mathbf v_{i,b}$ by $\mathbf v_{i,b}\min(1,\kappa/\|\mathbf v_{i,b}\|)$ before the projection onto the
moment-tensor subspace, with $\kappa=6$ in units of the noise standard deviation, fixed in Task 1x; it is
a robust statistic of Huber type, not the dual certificate of (16a), which bounds the influence of any one station and phase. In dense
windows it raised recall at 0.25 false events per window by $0.043$ to $0.074$ over BOMP-PK in two fresh
seeds, and without it MT-FW equals BOMP-PK. The verification of that task located the effect: the cap
bounds the influence of any single station and phase block whose amplitude the radiation operator gets
wrong, it calibrates the score better across windows under one global threshold, its gain depends on
the level $\kappa$, and it does not separate overlapping events better within a window. The amplitude
factors of the generator are drawn afresh in every window, so they are perturbations of each window and
not station terms; removing them reduces the gain of the cap, but at a small spread of them the cap still
gains where a statistic given the true factors does not (Appendix Y.25, exploratory). The mechanism is
therefore bounded block influence under amplitude model error, not the correction of a site term. A gain
$g_{i,b}>0$ per station and phase in the operator, $\Sigma^{-1/2}\mathrm{diag}(g)G(\theta)$, estimated
from catalogued events as the station time corrections are, was examined in that exploration with
persistent factors added to the generator, and neither the estimated gains nor the true ones took the
place of the cap. The form that the exploration supports keeps the bounded influence and removes the
calibrated level: a Huber M-estimate of the moment tensor on the residual norms of the blocks,
$T=\sum_b\rho_c(\|\mathbf v_b\|)-\min_a\sum_b\rho_c(\|\mathbf v_b-Q_b a\|)$, with $b$ running over station
and phase blocks, $Q$ an orthonormal basis of the radiation subspace of the candidate, $\rho_c$ the Huber
function and the threshold $c$ fixed from the noise distribution as the square root of the 0.95 quantile
of a chi-square with three degrees of freedom, $c=2.7955$, with four passes of iteratively reweighted
least squares from the least-squares fit and deflation with the robust fit. In the pilot, as the best of
11 candidates, it matched $\kappa=6$ at amplitude spreads of 0.5 or less, was worse at 0.8 and cost nine
times as much. The preregistered Task 2I confirmed it on two fresh seeds (Appendix Y.26): $+0.058$ and
$+0.053$ over BOMP-PK in dense windows at 0.25 false events per window, not inferior to $\kappa=6$
within $0.03$, and equivalent to BOMP-PK within $0.05$ in standard windows when pooled. Its descriptive
arms found the level immaterial between 2.5 and 3.4 noise standard deviations and the Huber fit behind
$\kappa=6$ by $0.055$ and $0.069$ at an amplitude spread of 0.8. The Huber form
is not a better detector than the cap: it costs about nine times as much, it is slightly below the cap
at 0.1 false events per window in one seed, and with every atom accepted it finds no more events, so
that its gain, like that of the cap, is in the ranking and selection of atoms. What it changes is the
origin of the level, a rule on the noise distribution in place of a calibration. Neither form closes the
headroom of 0.08 to 0.24 below a realistic oracle ceiling, which is not explained by the velocity model
or the source pulse, and of which no treatment of amplitudes, the true factors included, captured more
than about a quarter, so that about three quarters of it remains unexplained (Appendices Y.22, Y.23,
Y.25).

**The timing budget.** Two exploratory studies on the calibration seed asked whether that headroom could
be reached by a soft timing prior in place of the hard window, or by evidence that does not use travel
times, and both answers are negative with a mechanism (Appendix Y.27). The shift of a true arrival from
the prediction of its nearest node has three terms: the model error after a perfect station term, 0.05
to 0.07 s; the error of station corrections estimated from twenty catalogued events, 0.12 to 0.24 s rms;
and the position of the source inside a cell of the 3 km grid, 0.18 s for P and 0.30 s for S. The budget
is therefore 0.2 to 0.4 s per arrival, a Gaussian prior of that honest width is nearly flat inside
$\pm0.25$ s, and charging the shifts gained nothing in either detector, while widening the window to
$\pm0.5$ s lost recall. With the true station delays the same prior at the level of the model error
gains 0.081 in dense and 0.146 in standard windows, where it reaches the realistic ceiling, so the form
of the prior is right and its inputs are not accurate enough. Travel-time-free evidence fails for another
reason, information: with a P arrival of 1.5 to 1.8 noise standard deviations the polarization and
amplitude channels together bound the epicentre at 4.4 to 5.6 km against 0.4 to 0.7 km for the times,
and adding them to the times improves the bound by 4 to 5 per cent. What these studies leave as the
lever is the accuracy of the station time corrections, and after it the amplitude evidence for weak
arrivals. The preregistered Task 2L tested that lever (hypothesis G1-C, Appendix Y.28): station terms
estimated without picks from the same twenty catalogued events, by a joint determination of hypocentres
and terms on matched-filter energy peaks, raised the origin-time-matched recall of the capped detector
with the local refinement by $+0.057$ and $+0.062$ in standard and by $+0.038$ and $+0.030$ in dense
windows, a gain that is essentially that of the one-sided eight-station cell, is conditional on the
refinement and reaches 0.28 to 0.61 of that of the true delays.

The non-degenerate source condition under which the BLASSO recovers the support exactly, proved in
[47] for scalar amplitudes on a one-dimensional domain and analysed for scalar amplitudes in several
dimensions in [140], has a global part, that no other hypothesis
explains the noise-free signal, and a local part, the curvature of the certificate at each source. For
one source the local part is the efficient Fisher information for position and origin time with the
moment tensor as a nuisance parameter, which is bounded below by the leading term of the observability
certificate of Section 6.1 with each station and branch weighted by the energy that the actual mechanism
radiates to it, so that a station on a nodal plane of a branch contributes nothing through that branch
(Proposition AA.6). The observability certificate therefore predicts where the local condition holds;
that the group BLASSO then recovers the support exactly is Conjecture AA.1, and whether two events merge
is decided by a computable pre-certificate rather than by a closed form. The level of physics in $G$ follows a ladder chosen per station by
measured fidelity: travel times only; ray geometry (P along the ray, S across it), which ignores
amplitudes; the moment-tensor radiation subspace; and full Green's functions in a 3D medium, computed
with neural travel times [8] and source-dependent neural operators [126]. Station time corrections from
catalogued events and a theory-error covariance enter as nuisance
parameters, and amplitude errors are handled by bounding the influence of each block (Appendix Y.25). The learned components (Green's function surrogates, a heavy-tailed noise model in place
of the per-station cap, signed characteristic functions and amortized proposals for the oracle) only
accelerate or robustify the computation; the certified oracle guarantees that a proposal network
cannot make the search miss a hypothesis the exact search would find.

What is not new, and not claimed, is set out in Sections 3.2 and 3.6. The block-sparse moment-tensor
dictionary with greedy multi-source deflation is that of [129], [130], related to [131] and [141], and
MT-FW is in its algorithmic core a grid relative of that pursuit; the statistic at zero is, up to
normalisation, the variance reduction of operational grid moment-tensor monitoring [132], [123], so
Proposition AA.2 gives an optimisation meaning to an operational statistic; mechanism-corrected
stacking [121], [122], [133], [134], [112], [135], [136], subspace detection [118], [137], master-event
corrections [138], station time shifts estimated jointly with a moment tensor [143], [144],
extreme-value thresholds [125], conformal outlier testing [139] and the BLASSO
machinery [113], [114], [46], [47], [116], [117], [124], [140] are prior art, and QuakeMigrate [109],
GENIE [6] and HARPA [9] are comparators or reference points, not components. Task 2D sharpened this:
given the same peak search, block orthogonal matching pursuit [129] matched MT-FW in standard windows,
and the only element of MT-FW beyond [129] so equipped, its per-station cap, had a small effect there
that was not significant, so neither the greedy pursuit, nor its deflation, nor the cap, nor the peak
search is claimed as new, and neither are station time shifts as such. The principled treatment of the
timing tolerance described above, time shifts of the kind that [143] and [144] estimate for one known
event made nuisance parameters of a detection certificate over an unknown number of events, was tested
and did not beat BOMP-PK (G1-M, Appendix Y.22). Per-station clipping is established practice in robust
stacking and back-projection, the family reviewed in [136]. What the method adds is empirical and
specific: the measured and explained effect of that clipping inside a moment-tensor detector with
deflation in dense sequences (G1-R, Appendix Y.24), and of its bounded-influence form with a Huber fit,
itself standard robust statistics, at a level fixed by rule from the noise (Task 2I, Appendix Y.26),
together with the measured effect of station terms estimated without picks inside the detector with a
local timing refinement, joint hypocentre determination with station terms being classical in location
from picks (G1-C, Task 2L, Appendix Y.28). What remains, stated as absence in targeted searches rather than as priority,
is listed at the end of Section 3.6.

---

# 8. Constrained Optimization and Sensitivities

## 8.1 Stage and block structure

The coupled problem is non-convex because activations, source parameters and learned components
interact. Updates are staged rather than monolithic. For continuation stage $m$ and inner iteration
$q$: the moment tensors are updated first as a linear or robust-linear subproblem, which removes the moment-amplitude search from the non-convex positional search because the forward map is linear in $\mathbf M_k$ at fixed source coordinates and $s_k$; then positions and origin times subject to
the feasible set (15); then moment amplitudes and slacks, a convex or group-sparse subproblem given $\Theta$; then, in
designated experiments only, the parameters of the learned representation. Medium updates are
restricted to designated experiments after the fixed-medium problem passes its gates.

## 8.2 KKT conditions and diagnostics

With inequality constraints $g_j\le0$ collecting (15)-(16), including the kinematic admissibility family C-III,

$$\mathcal{L}(x,\lambda)=f(x)+\sum_j\lambda_j g_j(x),\qquad \lambda_j\ge0, \tag{17}$$

Stationarity, primal and dual feasibility and complementarity are solver diagnostics logged every
iteration. Because the problem is non-convex, KKT is a necessary first-order condition under
constraint qualification and is never read as a certificate of global optimality. The observability
constant $\gamma$ of C-I is what supplies the local second-order statement where it is positive.

## 8.3 Sensitivities: analytic, adjoint, cross-validated

Three routes are used and checked against each other: analytic in $\mathbf{M}$ and $\tau$, exact by
(5) with no adjoint required; adjoint in $\mathbf{x}_s$ and in the medium, with
$[R_\mathbf{u}]^{\mathsf T}\boldsymbol{\lambda}=-J_\mathbf{u}$ and reduced derivative
$dJ/d\Theta=J_\Theta+\langle\boldsymbol{\lambda},R_\Theta\rangle$ following standard adjoint
methodology [40]; and automatic differentiation through the learned representation [76], accepted
only after agreement with the first two within the F3 tolerances. The learned model never defines
its own correctness criterion.

## 8.4 Methods and basin analysis

Primary solver family: primal-dual interior point, with SQP as comparator [74], [75]. Basin
dependence is studied through matched multi-start runs with and without continuation, which is the
$2\times2$ design of F2. Logged per run: primal and dual residuals, complementarity, feasibility
violation, objective, iterations, forward evaluations, wall-clock, and the active-constraint set,
so that which constraints bind is data rather than speculation.

---

# 9. Identifiability

## 9.1 Local identifiability from observability

Where the coercivity constant (8) is positive, the Gauss--Newton reduced Hessian is positive definite on
the admissible tangent space, providing a local curvature and identifiability certificate. Strong
convexity and local uniqueness of the fully nonlinear constrained objective require additional
second-order regularity, control of residual-Hessian terms, and constraint qualifications, and are not
asserted by (8) alone. Cramer--Rao bounds derived from $\mathcal{I}(\Theta)=\mathbb{E}[\nabla_\Theta\log p\,\nabla_\Theta\log p^{\mathsf T}]$ are reported as a statistical resolution benchmark when the
likelihood assumptions are satisfied [26]. Known degeneracies are quantified
rather than ignored: the depth-origin-time trade-off at limited aperture, and the null space of
$\mathbf{M}$ under restricted azimuthal coverage, in particular the weak resolution of the isotropic
component. Moment-tensor components are reported only where the conditioning analysis indicates
identifiability.

## 9.2 The rank budget at one station

The minimum-station count for the two-branch operator is derived here rather than asserted, because
the four-parameter problem and the reduced model of Appendix G give different answers and the
difference is exactly where a reader should press.

**Split the information.** Write

$$\mathcal{I}_{\mathrm{full}}=\mathcal{I}_{\mathrm{arrival}}+\mathcal{I}_{\mathrm{pol}}+\mathcal{I}_{\mathrm{amp}}+\cdots \tag{18}$$

The leading term of the Proposition of Section 6.1 is $\mathcal{I}_{\mathrm{arrival}}$ only. Let $\mathcal S_2$ denote the stations for which both P and S branches are captured in the declared window.

**Arrival information is rank-deficient at one station.** Both branch vectors at a station,
$\mathbf{v}_P=(-\mathbf{n}/V_P,1)$ and $\mathbf{v}_S=(-\mathbf{n}/V_S,1)$, share the same direction
$\mathbf{n}$, so they span at most two dimensions of $\mathbb{R}^4$:

$$\operatorname{rank}\mathcal{I}_{\mathrm{arrival}}\big|_{S=1}\le 2<4. \tag{19}$$

In general $(\mathbf{a},c)$ annihilates the arrival term iff
$\mathbf{a}\cdot\mathbf{n}_s=cV_b$ for every captured station and both branches; the two branches at
a single station force $c=0$, and then $\mathbf{a}\perp\mathbf{n}_s$ for all $s$. Therefore

$$\operatorname{rank}\mathcal{I}_{\mathrm{arrival}}(\mathbf{x},\tau)=4
\iff\{\mathbf{n}_s\}_{s\in\mathcal{A}}\ \text{spans}\ \mathbb{R}^3
\implies S^\ast_{2,\mathrm{full}}=3, \tag{20}$$

For the full four-parameter problem, $S^\ast_{2,\mathrm{full}}=3$ under the arrival-rank criterion. The reduced Appendix G result $S^\ast_{2,\mathrm{reduced}}=1$ is valid only for its depth-fixed, one-dimensional source-coordinate model, where the unknown vector is two-dimensional. **It does not transfer to the full four-parameter problem**, and no claim in this document rests on its doing so.

**What two branches do supply at one station** is the separation of distance from origin time,
since $t_S-t_P$ is independent of $\tau$ while a single branch confounds the two irrecoverably.
The remaining directions must come from the subleading terms of (18). Polarization supplies the
direction of $\mathbf{n}$, in principle two angles, but the incidence component is weakly
conditioned at a free surface where the apparent incidence angle is distorted by the surface
conversion. The following is an interpretive rank budget, not a theorem about the full nonlinear inverse problem: the practical information at one three-component station is therefore

$$\underbrace{2}_{\text{distance and }\tau\ \text{from arrivals}}
+\underbrace{1}_{\text{back-azimuth from polarization}}
+\underbrace{\varepsilon}_{\text{incidence, weakly conditioned}}, \tag{21}$$

three well-conditioned directions against four unknowns, which is precisely why depth is
ill-conditioned from one station. This accounts for the measurements of Appendices G and P without
invoking a rank property the operator does not have: Appendix P fixes depth, leaving
$(x,y,\tau)$, and recovers two overlapping sources to about a kilometre.

**Full four-parameter single-station identifiability is a conjecture outside the reduced setting
and is not claimed anywhere in this document.** What is claimed is the operational statement, with
its condition attached: from one three-component station, with depth fixed or independently
constrained, epicentral position, origin time and source count are recoverable for up to two
overlapping events. Single-station location from waveform-derived distance and back-azimuth is
itself established in the literature, including recent machine-learning implementations [86]; the
narrower claim here is the **multi-event** case.

## 9.3 Separability of overlapping sources

Guarantees for recovering a sparse measure from low-pass measurements require a band-limited,
essentially translation-invariant operator [46], [47]. The seismic operator - heterogeneous,
finite irregular aperture, vector-valued, unknown mechanism - is not one, so that theory is not
imported here. The route retained is coherence-based and
convex-geometric demixing [48], [49], whose conditions are stated through the mutual coherence of
the components rather than a Rayleigh distance and therefore require no translation invariance. The
target statement is exact support recovery whenever $\mu(k,k')\,s_{\max}<\theta$ for all active
pairs. The known weakness is declared: plain mutual coherence bounds suffer the square-root
bottleneck and go vacuous exactly in the near-degenerate regime of interest, so the programme
commits to the sharper route, a convex-geometric statistical-dimension condition or a restricted
eigenvalue over the active set, and reports the plain bound only as a first, weak certificate.

## 9.4 What has already been refuted

Whether the noise-limited resolution constant improves with the velocity contrast $V_P/V_S$ was
tested and refuted (registered null E2, Appendix G); at equal energy with $S\ge2$ the joint P+S
operator is $10$-$19\%$ worse than the S branch alone. The measured constant
$\Delta_{\mathrm{res}}=\kappa V_S/f_{\max}$ with $\kappa\approx0.10$-$0.13$ at SNR 20 is retained, and
the noiseless vanishing-derivative precertificate held at every separation tested down to $20$ m, so
the binding limit is noise and coherence, not separation.

---

# 10. Differentiable Wavefield Representation

## 10.1 Why not a time-domain PINN as the primary object

The documented search of Appendix E did not identify a published time-domain PINN reproducing 3-D
heterogeneous elastic wavefields at seismological frequencies with the accuracy a waveform-level
likelihood requires.
The documented obstacles are summarized in the comprehensive review [59], including spectral bias [62], loss-imbalance pathologies [61] and characteristic
failure modes on oscillatory problems [60]. The successes are elsewhere: frequency-domain operators
for the 3-D elastic equation, with a fortyfold memory reduction, roughly hundredfold speedup over a
spectral-element solver and a factor of about three hundred and fifty in FWI time [58]. The primary
representation here is therefore not a vanilla PINN.

## 10.2 Analytic beams as the primary synthesis

The forward map is synthesized by Gaussian beam summation: each beam follows a ray with complex
width evolving by a Riccati equation, in closed form and without ray-theoretical singularities at
caustics [41]. This permits **analytic differentiation of the beam approximation** with respect to source position,
with no spectral bias, no training and no abort risk; the approximation itself remains asymptotic
and its tolerance against the reference is measured as a function of frequency and heterogeneity. The design principle is deliberate: **what can be computed is not learned.**

## 10.3 Learned residual, factorized

Only the heterogeneity residual is learned, against the analytic kernel; $\mathcal R_\phi$ is dimensionless and has the same tensor shape as the kernel, and $\odot$ denotes componentwise multiplication.

$$\mathbf{H}_\phi(\mathbf{x},\omega;s)=\mathbf{H}^{\mathrm{beam}}(\mathbf{r}_s,\mathbf{x};\omega)\odot\big[\mathbf{1}+\mathcal{R}_\phi(\mathbf{x},\omega;s)\big], \tag{22}$$

so that the singularity, the geometrical spreading and the leading wavefront kinematics are carried
in closed form and $\mathcal{R}_\phi$ is bounded and smooth away from strong interfaces. The
learning problem changes from reproducing an oscillatory 3-D elastic field to correcting an
asymptotic approximation, which is where networks succeed; latent representations provide a relevant precedent for reducing the learned state before applying physics-informed constraints [64]. Training is performed in the
**frequency domain**, which removes the time dimension and, more importantly, matches the
continuation ladder of Section 7.3 stage by stage: at stage $m$ the representation is only required
to be accurate within $\mathcal{F}_m$ [65]. Removing (22) and training the unfactorized network
under identical budget is a mandatory ablation.

## 10.4 Symmetry-exact architecture

Two exact symmetries of elastodynamics are imposed by architecture rather than by loss penalty, so
that they hold identically for every set of weights.

**Reciprocity.** Symmetrizing over source-receiver exchange,
$\mathbf{H}^{\mathrm{sym}}_\phi=\tfrac12[\mathbf{H}_\phi(\mathbf{r},\mathbf{x})+\mathbf{H}_\phi^{\mathsf T}(\mathbf{x},\mathbf{r})]$,
satisfies $G_{ni}(\mathbf{r},\mathbf{x})=G_{in}(\mathbf{x},\mathbf{r})$ identically. This halves
the effective function class and adds a reciprocal consistency relation for every source--receiver
pair; this is a structural constraint, not a claim that the resulting samples are statistically independent.

**Rotational covariance.** Three-component data are vector-valued and the physics is covariant under
rotation of the sensor frame. Standard architectures treat the three components as independent
channels and must learn that covariance from data, which wastes capacity and generalizes poorly
across stations whose orientation is uncertain. SO(3)-equivariant designs obtain it by construction
[67], [68], [69]. Seismology handles sensor orientation as a hardware calibration problem rather
than an architectural one, so this is an open design space, and it is the natural architecture for
a formulation in which polarization is a primary constraint (Section 6.4).

The resulting object is called a **structure-preserving differentiable representation**, and the
term is defined narrowly so that it cannot be read as more than it is: linearity in $\mathbf{M}$, time-shift equivariance in $\tau$,
source-receiver reciprocity and rotational covariance hold **identically for every set of weights**,
which is an algebraic property of the architecture and is demonstrable. It does **not** mean that
the learned representation solves the elastodynamic equation exactly, or that it is numerically
exact in any sense; only the four named symmetries are exact, and the heterogeneity residual
remains an approximation with the error budget of Section 10.5. F4 tests both the exactness, with
random untrained weights, and the training-efficiency claim against a penalty-based baseline.

## 10.5 Viability ladder and abort triggers

L1 is the registered target level of this conditional line; L0 is a floor that must succeed for the
line to continue; L2 and L3 are stretch.

**Table 2. Viability ladder.**

| Level | Physics | Band | Budget cap | Abort trigger |
|:--|:--|:--|:--|:--|
| L0 floor | 2-D heterogeneous elastic | $\mathcal{F}_1$ | 20 GPU-h | $L^2$ error $>10\%$ at 50% of cap |
| **L1 primary** | **2.5-D / layered 3-D** | **to $\mathcal{F}_2$** | **150 GPU-h per medium** | **$>6\%$ at 50% of cap** |
| L2 stretch | full 3-D | to $\mathcal{F}_3$ | 300 GPU-h per medium | as L1 |
| L3 stretch | full 3-D | to $\mathcal{F}_4$ | 300 GPU-h per medium | as L1 |

The association study runs on the beam and strain-Green's-tensor references and **does not depend
on any level succeeding**; the representation affects F3, F4 and F5 only. Budget caps are hard:
exceeding a cap is itself the abort condition, because discretionary extension is how a failing
component consumes a thesis.

## 10.6 Comparison of representations

| Representation | Cost per source evaluation | Storage | $\partial/\partial\mathbf{x}_s$ | Error |
|:--|:--|:--|:--|:--|
| Direct solver per candidate | $O(f_{\max}^4)$ | negligible | adjoint | discretization |
| Reciprocity / SGT database | table lookup + contraction | $\sim10^{12}$ B | grid interpolation | interpolation |
| **Analytic beams** | closed form, $\sim$ms | negligible | exact | asymptotic |
| Beams + learned residual | network evaluation | $\sim10^{8}$-$10^{9}$ B | exact by autodiff | learned |
| Neural operator | batched | intermediate | continuous | learned |

F5 is therefore precise: does the learned residual deliver tenfold storage reduction and continuous
source-position derivatives at accuracy within the F3 tolerances, relative to the SGT reference, and
does it beat analytic beams by enough to justify existing at all? A negative answer is publishable.

---

# 11. Computational Architecture

Four roles are separated: candidate generation, physical forward modelling, constrained attribution,
uncertainty evaluation. Candidate generators - heuristic clustering, a lightweight GNN, network-wide
neural pickers [13], or transport-informed proposals - never determine final attribution.

Components: three-component windows and metadata; candidate proposals; an adopted elastodynamic
solver as reference; receiver-side reciprocity and strain-Green's-tensor reference; analytic beam
synthesis; the symmetry-exact learned residual; declared observation and feature operators; the
five constraint modules of Section 6, each independently switchable; the misfit stack with
continuation; primal-dual and SQP optimization with cross-validated sensitivities; posterior and
simulation-based calibration; and the transfer and reproducibility package.

The hypothesis-space search is certified to return what the dense search returns. For the kinematic
likelihood, branch-and-bound over an octree with bounds derived from the eikonal Lipschitz property
returns identical catalogues with two to three times fewer evaluations (the bound for the
polarization term is open). For the moment-tensor statistic of MT-FW, the corresponding energy bound
is valid but avoids only 4 per cent of the exact evaluations, and the certified oracle is instead a
batched evaluation of every hypothesis with a floating-point enclosure and an exact recheck of the
winner, which reproduces the dense Frank-Wolfe records bit for bit and is 4.3 times faster through
vectorization, not pruning (Appendix Z). p4est [108] is reserved for distributing that search over a
cluster and for storing hypothesis spaces and 3-D travel-time tables that exceed the memory of one
machine at national scale.

**The primary comparison is the test of G1** (Section 13.0c). It runs MT-FW with catalogue station
corrections, its geometry-breaking control and its radiation-operator ablation against the published
associators (GaMMA and PyOcto, with REAL and GENIE registered to join them, Section 15), the published
migration locators (QuakeMigrate, then MALMI [110] and the stacking of [112]), the closest
antecedents, block orthogonal matching pursuit [129] and mechanism-corrected stacking [122], and the
grid and sliding solvers of the group BLASSO, all on the same windows, hypothesis grid, station
corrections, hit rule and false-event cap. Its pipeline is short: three-component windows, a declared
picker for the pick-based methods, matched filtering for the pick-free ones, the certified oracle of
Appendix Z, a threshold calibrated on noise-only windows, and the scorer of Section 13.0b. It needs
neither the elastodynamic solver nor the learned representation, and on real data it runs as the
network-decimation test of Section 16.

**The secondary comparisons are staged** so that scientific value is separated from representation
error. Stage A compares the constrained physical model against the ablation condition (ABL), which
measures what the physical channels add to the same search and is never a proxy for a published
method (Section 13.0a), against the waveform-coherence baseline and against the signal-based Bayesian
reference [18], with identical candidates, continuation schedule and optimization budget. Stage B
ablates constraint families one at a time, which is what makes A1 and A2 answerable. Stage C replaces
the reference synthesis with the learned representation, everything else frozen, and runs only if gate
F1 passes.

---

# 12. Computational Budget

## 12.1 Scaling

$$h=\frac{V_{S,\min}}{p\,f_{\max}},\qquad N_{\mathrm{dof}}=\frac{|\Omega|}{h^3},\qquad
\Delta t=\frac{C h}{V_{P,\max}},\qquad \text{Work}\propto f_{\max}^4. \tag{23}$$

**Table 3. Reference solver cost per simulation** ($|\Omega|=120\times120\times40$ km,
$V_{S,\min}=2$, $V_{P,\max}=8$ km/s, $p=8$, $C=0.4$, $T=120$ s, $2$ TFLOP/s sustained).

| $f_{\max}$ | $h$ | $N_{\mathrm{dof}}$ | steps | work | wall-clock |
|:--|:--|:--|:--|:--|:--|
| 1 Hz | 250 m | $3.7\times10^7$ | $9.6\times10^3$ | 0.07 PFLOP | $\approx35$ s |
| 2 Hz | 125 m | $2.9\times10^8$ | $1.9\times10^4$ | 1.1 PFLOP | $\approx10$ min |
| 5 Hz | 50 m | $4.6\times10^9$ | $4.8\times10^4$ | 44 PFLOP | $\approx6$ h |

## 12.2 Cost per objective evaluation

**Table 4. Decomposition ($K=20$).**

| Component | Unit | Per objective evaluation |
|:--|:--|:--|
| Analytic beam synthesis per candidate | $\approx0.5$ ms | $\approx10$ ms |
| SGT lookup and contraction per candidate | $\approx1$ ms | $\approx20$ ms |
| Learned residual per candidate | $\approx2$ ms | $\approx40$ ms |
| Transport misfit and gradient, 90 traces | $\approx25$ ms | $\approx25$ ms |
| Constraint evaluation C-I to C-V | $<3$ ms | $<3$ ms |
| **Total, beams** | - | $\approx38$ ms |
| **Total, SGT** | - | $\approx48$ ms |

## 12.3 Itemized budget

Constraint-attribution design (secondary): 100 realizations $\times$ 20 conditions $\times$ 200 iterations $\times$ 4
stages $\approx1.6\times10^6$ objective evaluations; the F2 study adds $\approx6.4\times10^6$;
constraint ablation multiplies the confirmatory block by six (all-on plus five leave-one-out).

**Table 5. Budget.**

| Item | Basis | GPU-hours |
|:--|:--|:--|
| Solver validation, mollifier convergence | $\approx200$ simulations | 35 |
| SGT reference, 12 synthetic media | $3S=90$ sims $\times$ 0.16 h $\times$ 12 | 170 |
| SGT reference, two real transfer domains | 90 sims each | 60 |
| Confirmatory runs with constraint ablation | $9.6\times10^6\times38$ ms | 100 |
| F2 misfit and continuation study | $6.4\times10^6\times38$ ms | 68 |
| G1 synthetic arm: sliding Frank-Wolfe, finer-grid block pursuit, MALMI and [112] on the Task 2A windows; second generator | planning allowance, CPU-dominated | 60 |
| Stress tier and regime-map surrogate | LHS plus boundary runs | 60 |
| Uncertainty calibration, SBI training and resampling | - | 140 |
| Representation L0 floor | cap | $\le20$ |
| Representation L1, 4 training media | cap 150 h each, provisional | $\le600$ |
| Representation L2/L3, conditional | cap 300 h, $\le4$ media | $\le1200$ |
| OOD evaluation on 8 held-out media | inference only | 20 |
| G1' decimation test on the RSNC: picking, MT-FW, QuakeMigrate, MALMI, noise-window calibration | planning allowance | 80 |
| Real-data transfer to the second domain | - | 40 |
| Contingency 15% | - | 215 |
| **Total without stretch levels** | - | **$\approx1{,}670$** |
| **Total with stretch levels** | - | **$\approx2{,}870$** |

Storage: $\approx15$ TB scratch for SGT databases, plus $\approx3$ TB for the transfer domains.

**Allocation between G1 and the conditional lines.** The test of G1 is cheap in compute: it needs no
elastodynamic solve and no learned representation, and the seed-999 smoke test of Task 2C estimated
0.79 hours on two cores for its network-replication part (Appendix Y.20). Of the $\approx1{,}670$
GPU-hours without stretch levels, 140 fund G1 directly (the synthetic arm, 60, and the RSNC decimation
test, 80); 605 fund the secondary programme (solver validation, SGT references, constraint-attribution
runs, stress tier, uncertainty calibration and transfer to the second domain); 708 are conditional on
gate F1 (the F2 misfit study, the representation levels L0 and L1 and the out-of-distribution
evaluation); and 215 are contingency. The stretch levels add 1,200 GPU-hours, all conditional. Most of
the budget thus funds lines that cannot decide the thesis because G1 is inexpensive, and the reduction
order below never reaches the two G1 lines.

**Status of the whole budget.** Tables 3 to 5 constitute a **planning model**, not a measurement.
The discretization rule $h=V_{S,\min}/(p f_{\max})$ is standard, but the wall-clock figures depend
on stencil, implementation, parallel efficiency, memory hierarchy, checkpointing, output volume and
boundary treatment, none of which are fixed until the solver is chosen in WP E2. Every figure is
therefore to be re-measured in E1 and E2 and the table re-issued; the ordering between architectures,
which is what the argument rests on, is expected to be less sensitive to these factors than the
absolute numbers, and that expectation is itself checked by the re-measurement.

**Status of the representation figures.** L1 and L2/L3 are *caps*, not measurements; no published
result establishes the training cost of a factorized elastic residual at these bands. They are
provisional planning figures to be replaced by the L0 and early-L1 pilots and are enforced as hard
abort triggers. Every other line except the two G1 allowances derives from the cost model of Tables 3
and 4, which is itself part of the planning model and is re-measured in E1 and E2. The
representation is trained on four media only; the remaining eight are held out for
out-of-distribution evaluation, since training on all twelve would inflate the budget and destroy
the OOD test that F3 depends on.

**Budget verification is a first-quarter deliverable.** If the allocation is below
$\approx1{,}670$ GPU-hours the reduction order is fixed in advance: stretch levels, then number of
media, then number of constraint-attribution conditions, and last $f_{\max}$, because $f_{\max}$
changes scientific scope rather than sample size. The two G1 lines are not reduced.

---

# 13. Experimental Design

The primary experiment is the test of G1, whose design is Section 13.0c: a synthetic benchmark in the
design of Tasks 2A to 2D (Appendices Y.18 to Y.21) and the real-data network-decimation test of Section
16. The tiers of Sections 13.1 to 13.5 serve the secondary endpoints; they report $E_{\mathrm{det}}$ for
every condition as a secondary reading and cannot by themselves decide G1.

## 13.0a The ablation condition, the published baselines, and what each establishes

Two different comparisons appear in this document and conflating them was the largest defect of
earlier drafts.

**The ablation condition, ABL.** Every pre-study in Appendices G to W compares the constrained
model against a condition implemented here: a candidate scores the fraction of its predicted
$t_P$ and $t_S$ arrivals falling within a declared tolerance of some pick,

$$s_{\mathrm{ABL}}(\mathbf{x},\tau)=\frac{1}{2N}\,\#\{(i,b):\ \exists\,p\in\mathcal{P}_i,\ |t^i_b(\mathbf{x},\tau)-p|<\delta\}, \tag{24a}$$

and the constrained score is the same object with the physical channels switched on,
$s_{\mathrm{PHYS}}=s_{\mathrm{ABL}}\cdot(\frac1M\sum_i|\hat{\mathbf{p}}_i^{\mathsf T}\mathbf{p}_i|)^2$
under the C-III admissibility filter. ABL is therefore an **ablation**, not a competitor: it
answers "what do C-III and C-IV add to the same search without them?" It is **not** a proxy for
published associators and is never reported as one.

**The published baselines.** GaMMA [3] fits a Bayesian Gaussian mixture over event components with
amplitude information; PyOcto [7] partitions four-dimensional space-time with a dedicated pruning
strategy. Both were executed on the same synthetic pick sets under an operating-point sweep
(Appendix X), with positive controls and paired exact McNemar tests. On the pick-free side,
QuakeMigrate 1.2.2 [109], the published waveform-migration locator, was run as released, with its own
STA/LTA onsets and trigger calibrated on a separate seed (Appendix Y.15).

$$\text{(i)}\quad \mathcal{J}_{\mathrm{phys}}\succ\mathcal{J}_{\mathrm{ABL}}
\qquad\qquad
\text{(ii)}\quad \mathcal{J}_{\mathrm{phys}}\succ\mathcal{J}_{\mathrm{GaMMA}},\ \mathcal{J}_{\mathrm{PyOcto}} \tag{24b}$$

**Statement (i) holds against ABL and does not survive a fair ablation.** Under the final protocol
of Appendix X the constrained prototype produces no false event in sixty one-sided realizations
against $0.30$-$0.43$ for ABL, $p\le0.0039$, but every ablation that breaks the ties of ABL's discrete
score with a continuous kinematic term, without polarization, does as well, and at three and four
stations better (Appendix Y.1). Statement (i) is therefore replaced by the question whether
polarization adds to the **best** kinematic score: with every arrival picked exactly it does not (Y.2,
Y.3); with degraded picks and random double-couple mechanisms the coherent moment-tensor detector with
explaining-away does, by $0.13$ to $0.41$ in event recall over its geometry-breaking control (Y.12),
and the per-pick likelihood (9d) does not (Y.10).

**Statement (ii) is G1 in its restated form (Section 4.3), and its evidence is stated once, in
Section 1.6.** In the consistent generator and in the preregistered Tasks 2A to 2D (Appendices Y.12,
Y.18 to Y.21) the coherent moment-tensor detector exceeds tuned GaMMA, PyOcto, QuakeMigrate and the
exact-time forms of block orthogonal matching pursuit and mechanism-corrected stacking, and block
orthogonal matching pursuit given the same peak search matches it, with the qualifications stated
there; as
formulated, G1 is not met at 0.1 against GaMMA fed optimistic picks with weak clutter; MALMI [110] and
[112] have not been run; and the result is open on real data. In the exact-pick generator of Appendix
X, GaMMA produces no false event under library defaults (X.2) and two in thirty when tuned (Y.1), and
the cells in which the prototype beats it are mostly those in which GaMMA does not reach the
calibration recall target (Y.1). Earlier margins over GaMMA (Y.5 to Y.10) were inflated by E26, and
earlier protocols reported other results; Appendix X.1 documents why.

## 13.0b The endpoints, named, and the primary one

The primary endpoint of the thesis is that of G1,

$$E_{\mathrm{det}}:\ \text{event recall at a fixed realized false-event rate, read on the exact FROC curve at }0.1\text{ and }0.25\text{ false events per window.}$$

A hit is a true source within 5 km epicentral distance of an emitted event, scored by the same rule
for every method, and an emitted event that matches no true source is a false event. On synthetics the
operating point is read on the evaluation windows for every method, as in Tasks 2A to 2C and in part
(i) of Task 2D, whose part (ii) froze split-conformal thresholds on empty windows (Appendix Y.21); on
continuous real data (G1', Section 16 and Appendix D) the rate is counted per day and the threshold is
calibrated on noise-only windows of the same network. Every other endpoint in this document, including
the two below and the false-association rate at non-inferior recall of Section 13.4, is secondary to
$E_{\mathrm{det}}$.

Two further endpoints appear in the pre-studies and are named here so that no passage can be read as
the other.

$$E_{\mathrm{loc}}:\ \text{does the argmax of the score field lie within }5\text{ km of a true source?}$$
$$E_{\mathrm{cat}}:\ \text{does the emitted catalogue contain any event beyond }5\text{ km of every true source?}$$

Appendices G to W measure $E_{\mathrm{loc}}$. Appendix X measures $E_{\mathrm{cat}}$, but with
exactly one source per realization and, for a continuous score at a threshold equal to the
realization maximum, a single emitted event, so there it largely reduces to $E_{\mathrm{loc}}$;
Appendix Y.3 measures $E_{\mathrm{cat}}$ with zero to two sources per window and an absolute
threshold. A method can
be strong on one and weak on the other, and the constrained prototype was, which is the subject of
the rest of this subsection.

Three causes of the difference were separated by running three protocols (Appendix X.1): in order of
size, the threshold scale, the geometry sampling, and a defect in the polarization estimator, which
multiplied the score by the squared **mean** cosine over matched arrivals and so rewarded candidates
that explain less, a defect invisible on an argmax endpoint and damaging on a catalogue endpoint. The
summed form, registered after P1, used in P2 and evaluated under P3 in Appendix Y as PHYS-S,

$$s_{\mathrm{pol}}=\frac{1}{2N}\sum_{i=1}^{M}\big|\hat{\mathbf{p}}_i^{\mathsf T}\mathbf{p}_i\big|^2, \tag{24c}$$

cannot be inflated by a small matched set. The per-station likelihood of Section 6.4 supersedes both
forms in the confirmatory programme, because it is also a sum and additionally makes a failing station
exactly neutral.

**A symmetry that must be stated.** The sixteen pre-studies of Appendices G to V, the baseline study
of Appendix X and the controls of Appendix Y use the same diagnostic implementation of C-III and
C-IV over a discrete search: not the elastodynamic forward model, not
Gaussian beams, not the learned residual, not continuous optimization. If that prototype does not
represent the proposed model when it loses, it does not represent it when it wins either. The
results of Appendices G to Y are therefore claims about **the constraint prototype**, and the
confirmatory programme is what tests the model.

## 13.0c The test of G1

**Synthetic arm.** The design is that of Tasks 2A to 2D, fixed in their registrations (Appendices Y.18
to Y.21). A generator the detectors do not know (Appendix Y.14, with station delays fixed per
network) produces windows with zero, one or two double-couple sources, in
three cells (random networks of four and six stations, and one-sided networks of eight) and two
clutter regimes, in two evaluation seeds never used for calibration; every calibration uses a third
seed on distinct networks. Each cell has five networks of twenty windows, 1200 windows in all, and
the comparison with the strongest comparator is replicated over thirty networks and 400 windows per
cell. The methods are MT-FW with station corrections estimated from twenty catalogued events carrying
location errors (MT-FW+SC); its geometry-breaking control, which uses the moment-tensor subspace of
another node; its radiation-operator ablation, which replaces the operator by the identity inside the
same solver; GaMMA and PyOcto fed a calibrated STA/LTA picker and, as an optimistic bound, the
parametric picks; QuakeMigrate in its three-component configuration with the same corrections; block
orthogonal matching pursuit in the form of [129] and mechanism-corrected stacking in the form of
[122], read at the exact predicted times and with the peak search of MT-FW; MT-FW without its peak
search and without its cap; and the group BLASSO (16a) solved on the grid to certified optimality, on
the exact-time dictionary and on a dictionary grown at the peaks. All share the station
corrections and the 5 km hit rule. The endpoint is $E_{\mathrm{det}}$ at 0.1 and 0.25 realized false
events per window; the network-level decisions use a hierarchical bootstrap over networks and then
windows, with Holm's correction [142] within each clutter regime. Task 2D (Appendix Y.21) reports what
equal tolerance and a frozen threshold change. With equal tolerance block orthogonal matching pursuit
matched MT-FW (-0.155 to +0.091 at 0.25, every hierarchical interval containing zero), the peak search
mattered (0.06 to 0.41 in 12 of 12 cases) and the cap did not; and split-conformal thresholds frozen
on noise-only windows of other networks met their rule cell by cell but, being calibrated on empty
windows, left 1.2 to 2.2 false events per weak-clutter window containing events. The synthetic arm
therefore tests the family of timing-tolerant moment-tensor detectors, carries BOMP-PK as the
comparator that any proposed method must beat, and still lacks a deployable operating point
calibrated for a false-event budget over all windows; with five networks per cell none of the rules
of Task 2D survives Holm's correction. Two further preregistered tasks used the same machinery for the
methodological hypotheses rather than for G1. Task 2E ran the off-grid group BLASSO with time shifts
against BOMP-PK on the windows of Task 2A and thirty networks per cell and refuted G1-M (Appendix Y.22).
Task 2F, on two fresh seeds with 96 networks per cell, added windows with three to five overlapping
events within 9 km and confirmed that the per-station cap raises recall there (G1-R, Appendix Y.24).
Task 2I, on two further fresh seeds with 120 networks per cell and regime and five dense and five
standard windows per network, confirmed the same for a Huber fit with a level fixed by rule from the
noise, not inferior to the cap (Appendix Y.26). Task 2L, on two more fresh seeds with 160 networks per
cell and regime in standard and 100 in dense windows, confirmed that station terms estimated without
picks raise the origin-time-matched recall of the capped detector with the local timing refinement, a
gain concentrated in the one-sided eight-station cell (G1-C, Appendix Y.28). The dense design joins the synthetic arm for any method
proposed next, with BOMP-PK, the capped statistic and the Huber statistic as the references to beat.
In Task 2L and in every task registered from here on the primary hit rule is the origin-time-matched one, which
requires an origin-time error of at most 0.5 s besides the 5 km distance and which both verifications
showed to matter in dense windows; the location-only rule is kept as a second reading, and equivalence
is decided per clutter regime with a margin of 0.02.

**Real-data arm.** The decisive test is the decimation of a dense part of the RSNC registered in
Section 16 as G1', with its comparators, its two declared pickers, a threshold calibrated on noise-only
windows of the same network and recall against false events per day above $M_c$. Its sample size is
fixed from the paired standard deviation measured at those station counts (Section 13.0), and its
first result is a month-6 gate (Section 19) whatever its sign.

## 13.0 The Tier-1 pilot, defined

Several registered quantities depend on "the pilot", so it is specified. Twenty independent
realizations in each of four anchor conditions spanning the corners of the clutter and SNR factors,
at $f_{\max}=0.8$ Hz, thirty stations, one medium. Deliverables: the paired standard deviation
$\sigma_d$ per registered endpoint, feeding the sample-size formula (25); measured wall-clock per
objective evaluation against Table 4; the empirical distribution of the cell coherence statistic,
fixing $\varrho$; the erosion tolerance $\eta$ for the exclusion sets of C-II, calibrated so that no
true source is excluded; and calibration of the whitening operator and $\nu_w$. Run once, frozen,
never reused as confirmatory data. The pilot additionally includes a **design-validity step**:
each constraint family must be demonstrated to bind on a case constructed so that it must
bind, before any confirmatory condition is run. Appendix L explains why this step is not
optional.

Thirty stations lie where C1 predicts no advantage of physical evidence over travel time: Appendix O
found the advantage over the coincidence-count ablation to be zero at six and twelve stations, with
$S^\ast\approx5$-$6$. The pilot therefore does not test G1, and a null there would not count against
it. What it tests is the pipeline on a dense configuration, and what it calibrates is listed above.
The paired standard deviation it supplies serves the constraint-attribution conditions; for G1 it is
re-estimated at the station counts where G1 is tested, the decimated sub-networks of three, four and
six stations of Section 16, before any sample size for G1 is fixed.

## 13.1 Tiers

**Table 6. Tiers.**

| Tier | Purpose | Design |
|:--|:--|:--|
| 0 Primary | G1 on synthetics; G1' on the RSNC | Section 13.0c: benchmark in the design of Tasks 2A to 2D with a threshold frozen before evaluation; network decimation of the RSNC (Section 16) |
| 1 Reproducibility | Verify solver, SGT, beams, baselines | Low and medium event rates, controlled media, all mandatory baselines |
| 2 Constraint attribution | G2, G4, A1, A2, E1, and per-constraint effects; $E_{\mathrm{det}}$ as a secondary reading | All-on plus five leave-one-out conditions, crossed with clutter rate and station count |
| 2b Design and mechanism | B2, B3, C3, F2, F4 | Observability sweep over geometries; single-station study; $2\times2$ misfit design; symmetry ablation |
| 3 Stress and failure map | Boundaries | High overlap and multiplicity, very high clutter, sparse networks, strong heterogeneity, 5 Hz reduced-dimension control |

## 13.2 Generator and factors

Each realization specifies event count and geometry, origin times, moment tensors, heterogeneous
$\rho,\lambda,\mu$, station geometry and count, bandwidth, SNR, structured false arrivals
(coda-associated, clustered impulsive, near-event distractors, phase confusion), overlap ratio, and
optional model mismatch. Data-generating and inference physics are separated so that misspecification
is studied rather than assumed away. Factor families: event count and separation; **false-arrival
rate**, which the pre-study identifies as the decisive factor; **station count and aperture**;
SNR and bandwidth; heterogeneity and interface complexity; **station distance spread**, which Appendix J
identifies as the factor that switches the silence constraint on; velocity-model error, which activates
cycle skipping; overlap ratio; mechanism diversity; misfit type and continuation schedule;
constraint family on/off; representation; initialization and budget; slack margins.

A full Cartesian product is not executable and not required; the choice of screening and subsequent confirmation follows the experimental-design logic of [80]. Tier 2 uses a definitive screening
design [77] to identify active main effects and two-factor interactions, followed by a full
factorial on the factors that screen active. Tier 3 uses a Latin-hypercube sweep with a
Gaussian-process surrogate, with confirmatory runs only at surrogate-identified boundaries.

## 13.3 Reference validation

Analytic and benchmark solutions; moment-tensor radiation patterns; energy and receiver
consistency; mesh and time refinement; mollifier convergence; agreement between direct multi-source
simulation and superposition; SGT-assembled seismograms against direct simulation on held-out
source positions; and beam synthesis against the reference within a declared asymptotic tolerance,
with the tolerance itself reported as a function of frequency and heterogeneity.

## 13.4 Constraint-attribution experiment

The experiment that decomposes which constraint earns the effect of G1; it does not test G1, which is
tested in Section 13.0c. Matched realizations,
identical candidates, noise, schedule and budget. Conditions: all constraints on; five leave-one-out
conditions; ablation condition (ABL); waveform-coherence baseline; signal-based Bayesian baseline.
Its endpoint, the false-association rate at non-inferior recall, is secondary to the thesis endpoint
$E_{\mathrm{det}}$ of Section 13.0b, which is reported for every condition; further secondary
endpoints are GOSPA, event-count error, under- and over-segmentation and localization error. Crossed
with clutter rate at four levels and station count at the five levels of Appendix D, which are the two
axes on which the pre-study predicts the effect lives.

## 13.5 Single-station and design experiments

The single-station study places one three-component station in windows containing one to four
overlapping sources plus clutter, and reports source-count accuracy, epicentral error, origin-time
error and depth Cramer-Rao bound, against the single-branch representation and against the
ablation condition (ABL). The observability study computes $\gamma$ over a family of network
geometries and correlates it with empirical estimator variance, which is B3, and additionally
evaluates the placement optimization of Section 6.1 against uniform and expert-chosen placements.

---

# 14. Metrics and Statistical Protocol

**Primary endpoint.** Event recall at a fixed realized false-event rate on the exact FROC curve, read
at 0.1 and 0.25 false events per window on synthetics and per day on continuous real data
($E_{\mathrm{det}}$, Section 13.0b; hypothesis G1). Every metric below is secondary to it. For
$E_{\mathrm{det}}$ the replication unit is the network, not the window. Tasks 2A and 2B decided their
rules on intervals over windows; Task 2C applied the declared unit to the comparison with QuakeMigrate,
with thirty networks per cell, a hierarchical bootstrap over networks and then windows, and Holm's
correction [142] (Appendix Y.20). Task 2D used the same bootstrap over the five networks of each cell
of Task 2A, with Holm's correction over 36 comparisons per regime; at that size none of its rules
survives the correction, although the effect of the peak search was positive in 12 of 12 cases
(Appendix Y.21). Task 2E used the thirty networks per cell of Task 2C with Holm's correction over the six
tests of cell and seed of each regime, and Task 2F used 96 networks per cell and seed, with the sample
size fixed in advance from pilot standard deviations at half the pilot effect and Holm's correction over
four tests; its rules passed with every $p$-value at the bootstrap floor (Appendices Y.22, Y.24). Task
2I used 120 networks per cell, regime and seed, sized on an effect below the pilot estimate because its
candidate had been selected among 11, with a superiority rule, a non-inferiority rule with margin $0.03$
and an equivalence rule with margin $0.05$, each required in both seeds and pooled over regimes and
cells; all three passed (Appendix Y.26). Task 2L used 160 networks per cell, regime and seed in standard
and 100 in dense windows, sized on an effect of $0.05$ against a pilot estimate of $0.070$, with two
superiority rules under the origin-time-matched hit rule and two non-inferiority rules with margin
$0.02$, per regime and under the location-only rule, each required in both seeds; all four passed
(Appendix Y.28). In
dense windows a hit rule that ignores the origin time is loose: the gain of Task 2F falls from $0.043$ to
$0.074$ to $0.024$ to $0.045$, and that of Task 2I from $0.058$ and $0.053$ to $0.031$ and $0.019$, when
the origin must match within 0.5 s. The origin-time-matched rule, an origin-time error of at most 0.5 s
together with the distance criterion, is therefore the primary hit rule of every task registered after
Task 2I, of which Task 2L is the first, with the location-only rule as a second reading, and equivalence
or non-inferiority is decided regime by regime
with a margin of 0.02 instead of pooled within 0.05; in Task 2I one per-regime value fell outside the
pooled margin. The other comparisons of Tasks 2A and 2B have not been repeated at that scale. One
multiplicity rule applies per family. Every decision on the primary endpoint $E_{\mathrm{det}}$, in
the synthetic tasks and in G1' on real data, uses Holm's correction [142] within its declared family of
tests, which controls the family-wise error rate; a rule required to hold in each of two seeds, as in
Tasks 2I and 2L, needs no correction. The secondary endpoints of the constraint-attribution experiments use the
Benjamini-Hochberg procedure at a false-discovery rate of 0.05.

**Attribution.** GOSPA as the set metric [54], with base distance
$d^2=\|\Delta\mathbf{x}\|^2/\sigma_x^2+\Delta\tau^2/\sigma_\tau^2$, $\sigma_x=2$ km,
$\sigma_\tau=1$ s, cutoff $c=2.5$, order $p=2$, $\alpha=2$, decomposed into localization, missed and
false components. Also event-level and pick-level precision, recall and F1 for comparability;
false-association fraction at non-inferior recall; event-count error and segmentation errors;
horizontal, depth and origin-time errors reported against the Cramer-Rao floor; moment-tensor error
decomposed into deviatoric and isotropic parts, reported only where conditioning permits.

**Representation.** Relative $L^2$ and RMS field error against the reference, stratified by band;
receiver error on held-out sources and media; PDE residual and boundary error; source-gradient
cosine similarity and relative magnitude against adjoint or analytic reference; OOD error across
held-out velocity families; scaling with collocation density, domain size, duration and storage.

**Uncertainty.** Brier score, reliability diagrams and expected calibration error for attribution
probabilities; empirical coverage, sharpness, CRPS and PIT for continuous parameters. Coverage is
never reported without sharpness.

**Statistical unit and power.** For $E_{\mathrm{det}}$ the replication unit is the network, as stated
above, with realizations nested within networks. For the secondary endpoints of the
constraint-attribution experiments, whose realizations each draw their own network, the independent
replication unit is the complete realization. In both cases picks within a realization are not
independent. Required $n$ per condition:

$$n\ \ge\ \frac{(z_{1-\alpha/2}+z_{1-\beta})^2\sigma_d^2}{\Delta^2}\ \approx\ \frac{7.85\,\sigma_d^2}{\Delta^2}, \tag{25}$$

with $\Delta$ the registered margin and $\sigma_d$ the paired standard deviation from the pilot; when the pilot
rejects normality of the paired differences the sample size is instead recomputed by a
simulation-based paired-permutation power calculation on the pilot data rather than by a fixed
inflation factor, and non-inferiority arms use the one-sided quantile. The working target $n=100$ corresponds to
$\sigma_d\approx0.178$; a larger pilot $\sigma_d$ either raises $n$ or drops the condition from the
confirmatory family, and that decision is registered before any confirmatory run.

**Comparisons.** Paired effect estimates with BCa bootstrap, $B\ge10^4$; paired permutation or
Wilcoxon signed-rank as pre-specified; Friedman with Nemenyi post-hoc for multi-method comparisons
[78]; every $p$-value accompanied by an effect size, $\hat{A}_{12}$ or Cliff's $\delta$ [79];
Holm's correction for decisions on $E_{\mathrm{det}}$ and Benjamini-Hochberg FDR at $0.05$ for the
family of secondary endpoints, as stated above; all thresholds frozen before held-out
testing.

**Error analysis.** Every confirmatory condition reports not only its rate but a characterization of
its failures: for false events, the geometric relation between the spurious solution and the truth
and which channel admitted it; for missed events, which constraint excluded them and whether the
exclusion was within its declared tolerance. The pre-studies of Appendices G to W report rates and,
in several cases, the mechanism of failure in prose; a systematic failure taxonomy per condition is
registered work and is currently absent; this remains an explicit methodological gap.

**Science per compute.** Training, reference construction, inference, memory, storage and
deployment counts reported separately, never aggregated into a single speed claim.

---

# 15. Baselines and Ablations

**Travel-time baselines.** GaMMA [3] and PyOcto [7] have been executed on the synthetic generator
(Appendices X and Y); REAL [88], GENIE [6], Neuma [4], HARPA [9], PhaseLink [5] and VORA [19] have
not. The pre-studies of Appendices G to W used the ablation condition (ABL) in their place; see
Section 13.0a. **Tuning protocol, registered.** Every baseline is tuned on the calibration
realizations only, over a declared grid of its principal hyperparameters (for GaMMA at least the minimum
pick count, the time-residual scale and the clustering radius, the last fixed at 10 in the grid of
Appendix Y for computational cost; for PyOcto the minimum pick count and the travel-time
tolerance), and frozen before evaluation, with the same calibration budget as the
constrained model. Library defaults are not an acceptable baseline configuration: Appendix X used
them and its baselines oscillated between seeds in a way that indicates misconfiguration rather
than method behaviour. Running REAL and GENIE, which the benchmark of [16] identifies as a light
and a robust associator respectively, is the first confirmatory experiment. All pick-consuming methods
receive an identical pick set from a fixed declared picker [11], [12], [13], chosen with reference
to the quantitative picker comparison of [103], so that differences are
attributable to the association model and not to detection. On synthetic waveforms that picker is a
calibrated STA/LTA, because PhaseNet is out of its domain on isolated pulses in white noise
(Appendix Y.14). In Task 2A that STA/LTA picked 10 to 17 per cent of the arrivals within 0.3 s, too
few for GaMMA and PyOcto to form most events, so the comparison against them holds for that picker
only (Appendix Y.18).

**Waveform baselines, mandatory.** Source-scanning and back-projection [14]; waveform migration and
stacking of onset functions, in QuakeMigrate [109] and in MALMI, which migrates learned phase
probabilities [110]; array-based waveform correlation, i.e. matched filter [15]; pick-free regression
location [21]; continuous-waveform neural-operator monitoring [10]. Appendix Y.10 makes the
migration locators the decisive comparators. QuakeMigrate 1.2.2 [109] has been run as released, with
its documented three-component channel mapping (Appendix Y.15), and in Tasks 2A and 2C it was the
strongest comparator (Appendices Y.18, Y.20); MALMI [110] has not been run. The closest antecedents,
block orthogonal matching pursuit [129] and mechanism-corrected stacking [122], were run in Task 2C
with exact predicted times on a shared grid and in Task 2D with the peak search of MT-FW, where block
orthogonal matching pursuit matched MT-FW (Appendix Y.21); that variant, BOMP-PK, is a mandatory
comparator for any proposed method, and a finer grid and the full-waveform dictionaries of [129]
remain to be run (Section 13.0c). Polarization already enters migration-based location in other forms, for example through
a polarization norm of the horizontal components in the P and S back-projection stacking of [112];
the claim of Appendix Y.10 is the marginal value of the ray-geometry constraints under controls, not
the first use of polarization in migration, and [112] is a required comparator.

**Azimuth-likelihood baseline, mandatory.** A travel-time associator augmented with a per-station
Laplacian back-azimuth residual likelihood and a Poisson false-arrival model, in the form of NETVISA
[17], with the back-azimuth taken from single-station P polarization. This is the fair test of
whether the constrained formulation adds anything beyond an established probabilistic treatment of
the same observable.

**Signal-based Bayesian baseline, mandatory.** SIGVISA [18], as the closest prior art and the
statistical-forward-model counterpart. Reported separately for source locations with and without
historical seismicity available for training, since that is the axis on which the two approaches
are predicted to differ.

**Internal counterfactuals.** Travel-time mode of the same constrained framework; mixture rather
than additive observation model; SGT reference, analytic beams, learned residual; each constraint
family on and off.

**Mandatory ablations.** Remove each of C-I to C-V individually and all together; strict versus soft
versus removed energy admissibility, with margin sweep; replace transport misfit by time-domain,
and remove continuation; replace additive superposition by exclusive assignment; remove waveform
evidence, retaining only the derived arrival diagnostic; remove the analytic factorization (22);
remove reciprocity symmetrization and SO(3) equivariance, replacing them by loss penalties; restrict
to a single station; remove moment-tensor information; remove learned initialization; remove
posterior calibration; replace the parametric likelihood by the simulation-based comparator.

A global improvement without per-constraint attribution is not treated as evidence of scientific
contribution.

---

# 16. Test Domains, Data, and Governance

Primary quantitative claims are made on synthetics, where the medium is known exactly and truth is
complete. This is the declared design, not a fallback.

**Real test domains.** At least two contrasting domains are required, chosen so that conclusions are
not an artefact of one network: a dense, well-instrumented, well-characterized region with an
independent high-quality catalogue, which is the hard case for the thesis because the pre-study
predicts the advantage is smallest there; and a sparse or irregular region with high false-arrival
rates, which is where the advantage is predicted to be largest. Public FDSN-compatible archives
supply both. The **primary real-data application domain is the Colombian RSNC** [81], [82], selected because it provides a nationally relevant and geographically heterogeneous test bed. A second contrasting external domain is required only for transfer validation, so that Colombian results are not treated as universally representative.

**Velocity-model gate, blocking.** Waveform-level inference requires a documented 3-D velocity model, and the achievable band is bounded by its resolution. For Colombia, define $f_{\mathrm{CO,max}}$ as the highest frequency at which the selected RSNC velocity model passes the pre-registered waveform-support test; all deterministic waveform claims on RSNC data are restricted to $f\le f_{\mathrm{CO,max}}$. For any candidate domain, the first work
package must identify the model, document its resolution and determine empirically the maximum
frequency at which it supports waveform fitting, by comparing synthetic and observed records for
well-located catalogue events. For the Colombian case the public candidates are a local-earthquake
tomography model [84] and ambient-noise crustal models on grids of order half a degree [83]; a
lateral resolution of order fifty kilometres with sensitivity from surface waves of six to
thirty-eight seconds cannot support deterministic waveform fitting at 2 Hz, so real-data waveform
claims there would run at long period and correspondingly larger magnitudes. That determination is
made per domain and reported, and under no circumstance is transfer executed at a band the model
does not support; the fallback is envelope- and polarization-level transfer with waveform-level
claims restricted to synthetics. Polarization-level transfer is itself conditional: Appendix T.1
finds it usable at one of five RSNC stations in the 1-10 Hz band, so the fallback applies to the
stations that pass the polarization gate, not to the network as a whole.

**Real-data calibration of C-IV, registered.** Before $\hat r_s$ of (9e) is used anywhere, each
station goes through the following sequence, in this order, because each step removes one
explanation of a failing gate. (1) Orientation audit by P-wave particle motion [91] and Rayleigh-wave
arrival angles [94], including changes over time [93]; a misoriented station is corrected, not
discarded. (2) Gate stratified by incidence angle and event depth, because near-vertical P has
little horizontal motion and SV at steep incidence is radial, which is what the CHI measurements of
T.1 show. (3) Gate in the band of Section 1.5, $0.1$-$2$ Hz, as well as in 1-10 Hz. (4) Labels built
from the SGC catalogue and the SeisMonitor catalogue [90] as well as the ISC catalogue, with clutter
picks screened by template matching for uncatalogued events. (5) Monte Carlo calibration of the
bias and interval coverage of the estimators of $\eta$ and $\hat r_s$. (6) A network-wide
reliability survey giving the distribution of $\hat r_s$ over the RSNC. (7) The clutter measurement
of T.1 repeated during an active sequence. These are tasks 1c to 1n of the companion task plan.

**Real-data test of G1', registered: network decimation.** A dense part of the RSNC with a complete
SGC catalogue above its completeness magnitude $M_c$ supplies the truth. Sub-networks of three,
four and six stations, including one-sided subsets, are formed by selecting stations, and every
method (GaMMA, PyOcto, REAL, the azimuth-likelihood baseline of Section 15, the continuous kinematic
likelihood with C-III, and the same with the polarization likelihood (9d) added) receives the same
PhaseNet picks, and the same likelihood with the directions of another hypothesis is run as its
negative control (Appendix Y.10). A pick-free arm runs on the same decimated waveforms: the two
timing-tolerant moment-tensor detectors with explaining-away of G1, MT-FW (Section 7.5) and block
orthogonal matching pursuit with the same peak search (BOMP-PK, Appendix Y.21), with station time
corrections calibrated on catalogued events, each against its geometry-breaking control, against the
mechanism-independent migration, against the polarity-corrected stacking of [112] and against
QuakeMigrate [109] and MALMI [110] as published, with the velocity model the pick-based methods use.
MT-FW against BOMP-PK is reported as a secondary comparison, which on real sequences with overlapping
events is the real-data reading of G1-R; the time-shift solver of G1-M does not join the arm, because it
did not beat BOMP-PK on synthetics (Appendix Y.22).
Station corrections are estimated from catalogued events by the joint P-before-S search of Appendix
Y.16, which on synthetics recovered the recall of the true statics from five to twenty events, and
their residuals are reported per station. The detection threshold is set by split-conformal or GEV
calibration on noise-only windows of the same network and noise regime, and recalibrated when that
regime changes, because strong transient clutter raised the realized false-alarm rate of every
threshold calibrated on Gaussian noise far above nominal, on average 19 to 22 times for the empirical
ones (Appendix Y.17). Because a threshold calibrated on noise-only windows does not bound false events
in windows that contain events, where on synthetics the detector emitted 1.2 to 2.2 per window with
weak clutter (Appendix Y.21), the false events of G1' are counted per day over all continuous data,
events included, and the realized rate is reported. The PhaseNet picks of the pick-based arm are within PhaseNet's domain on real data;
its failure on isolated synthetic pulses in white noise (Appendix Y.14) does not transfer, and a
calibrated STA/LTA runs as a second declared picker so that the comparison does not rest on one
picker. The kinematic likelihood alone is the reference against which the value of (9d) is
measured, because Appendix Y showed that a coincidence count is not a fair one. The endpoint is recall against false events per
day for events above $M_c$; an unmatched detection is screened by template matching against the full
network before it is counted as false, as the governance paragraph below requires, and the
false-positive calibration of template matching itself follows [99]. The design reproduces on real
data the cells in which the synthetic results were obtained and tests the prototype, not the full
waveform model. It needs no three-dimensional velocity model of waveform resolution, because nothing is
fitted in phase beyond the travel times of the model that the pick-based methods use. It is not
independent of the model, however: the radiation operator of the pick-free arm predicts the relative
amplitudes and polarities of P and S on three components from take-off directions in that model, and it
assumes P motion along the ray and S motion transverse to it at every station, whereas single-station P
polarization passed the reliability gate at one of five RSNC stations (Appendix T.1) and Task 2B
indicates that with strong clutter the advantage depends on that operator. Two further arms are
therefore registered with G1'. In the audited arm the moment-tensor subspace is built only from the
stations that pass the orientation audit and the polarization gate of the calibration sequence above,
the remaining stations contributing the energy of their peak vectors without sign or direction; the arm
is not run on a sub-network in which fewer than two stations pass, and that count is reported. In the
sign-blind arm, the fallback, every station contributes only the norm of the three-component peak
vector of each phase, which discards polarity and component direction and makes the statistic a
mechanism-independent energy stack with the same peak search, bound on influence and deflation. Both
arms run for MT-FW and BOMP-PK beside the full operator, each with its geometry-breaking control, and
the fraction of stations that pass the audit is reported per sub-network. G1' is decided on the
full-operator arm as registered in Appendix D; the two arms measure how much of its result survives
where the polarization assumption is not verified.

**Governance.** Held-out stations, dates and intervals are never used for tuning, and access
timestamps are stored so leakage can be audited rather than promised. Newly detected real events are
not counted as quantitative true positives without independent corroboration. Public data use
follows provider terms; only selections, configurations and derived products are redistributed.

---

# 17. Risks and Go/No-Go Criteria

**Table 7. Risks and gates.**

| Risk | Consequence | Response |
|:--|:--|:--|
| Cycle skipping under velocity error | Misfit measures model error, not evidence | Transport misfit with continuation is default; if F2 fails, restrict waveform claims to the band where it holds and report the boundary |
| Unique continuation invalid for the medium class | C-II unsound | Theoretical claim restricted to the declared piecewise-smooth class; general case uses the eroded empirical exclusion set; relaxation ladder falls back to the Bernoulli channel |
| Exclusion sets remove true sources | Recall loss | Erosion tolerance calibrated on the pilot so no true source is excluded; soft variant with slack; failure rate reported |
| Observability constant not predictive | C-I decorative | Report B3 as negative; certificate retained only as a feasibility restriction |
| Amplitude floors mis-set, catalogue-truncation bias | Biased detection thresholds | Fit on continuous data with template-confirmed events; report sensitivity to a declared perturbation |
| Polarization degraded by orientation or site effects | C-IV unreliable | **Observed**: the polarization gate passes at one of five RSNC stations (Appendix T.1). Per-station gate on catalogued events, with fixed-rotation, reflection and temporal-stability checks, before enabling; arrivals at failing stations enter C-IV as neutral evidence, never as absent (Appendix T.2); report the network fraction that passes |
| Catalogue-endpoint advantage is an artefact of scoring or thresholding | G1 unsupported | **Observed risk**: a discrete score ties and a threshold relative to the realization maximum cannot run on empty windows. Tie-breaking ablations and an absolute threshold calibrated on noise-only windows are part of the protocol (Appendix Y); a threshold frozen on empty windows does not bound false events in windows with events, 1.2 to 2.2 per window with weak clutter (Appendix Y.21); an exploratory decoy pass with permuted station moveouts reproduced the labelled operating point without labels under the location-only rule but not under the origin-time rule (Appendix Y.29) |
| Published baselines misconfigured | Advantage over the state of the art overstated | Tuning on calibration data over a declared grid, frozen before evaluation (Section 15); oscillation between seeds reported as a configuration symptom |
| Unmodelled station delays beyond the migration window | MT-FW recall falls (0.37/0.42 with delays of standard deviation 0.4 s, Appendix Y.11) | **Observed on synthetics**: corrections from five to twenty catalogued events recover the recall of the true statics in the 4-station and one-sided cells, but in Task 2A they helped with weak clutter only, because strong clutter contaminated the catalogue windows (Appendices Y.16, Y.18). Residual statics reported per station |
| Convex solver not competitive, on the grid or off it | The certified group BLASSO is not the best detector, and the theory covers a method that is not deployed | **Observed on synthetics** (Appendices Y.20 to Y.22): on the grid of Task 2A the certified solution stays 0.06 to 0.61 below MT-FW in 23 of 24 cases, and 0.11 to 0.30 below in 11 of 12 when given the same peak search; off the grid, sliding Frank-Wolfe with per-station time shifts gained $+0.000$ to $+0.039$ over BOMP-PK with every interval containing zero (Task 2E), and that negative is reported. A timing-tolerant greedy pursuit remains the deployed detector, and the theory of Appendix AA is the theory of the framework, extended next to the bounded-influence statistic (Section 19, T2) |
| The advantage belongs to the timing tolerance, not to the method | MT-FW is not distinctive: its advantage over the published comparators is shared with block orthogonal matching pursuit [129] given the same tolerance | **Observed on synthetics** (Appendices Y.21, Y.22, Y.24): with the same peak search block orthogonal matching pursuit matches MT-FW in every standard-window cell (-0.155 to +0.091, every interval containing zero), and the off-grid group BLASSO with time shifts does not beat it (Task 2E). The one element of MT-FW with a confirmed effect is its per-station cap, in dense sequences only and modest in size (Task 2F), and a Huber fit with a level fixed by rule from the noise gives the same gain at about nine times the cost (Task 2I, Appendix Y.26); the methodological contribution is stated at that size, as a measured design result with the benchmark and its negative results |
| A bounded-influence statistic tied to a tuning constant and to the amplitude model of the generator | The confirmed gain of the cap does not transfer: it peaks at the fixed cap level and depends on the lognormal amplitude factors of the generator, which are redrawn in every window | **Observed on synthetics** (Appendices Y.24, Y.25, exploratory): the gain reverses at a cap of 3 or less, vanishes at 24, and falls to $+0.023$ with an interval containing zero without the amplitude factors, although at a spread of 0.25 the cap still gains $+0.043$, with an interval containing zero, where the true factors gain nothing; catalogue amplitude terms do not replace the cap. The level of the cap stays at the value fixed in Task 1x in every confirmatory run. Task 2I replaced the calibrated cap by a level fixed by rule from the noise distribution in its registered primary arm, at the nominal spread: the Huber fit at a level fixed by rule from the noise was not inferior to the cap in two fresh seeds (Appendix Y.26). The dependence on the generator stays: both tasks use one generator family, the descriptive arms of Task 2I show a level that is immaterial between 2.5 and 3.4 noise standard deviations and a Huber fit behind the cap by $0.055$ and $0.069$ at a spread of 0.8, as in the pilot, with neither form gaining reliably at a spread of 0.25, and the gain shrinks to $+0.02$ to $+0.045$ with origin times matched. A second generator family is run before real data, where the capped and the Huber statistics are both reported |
| Recall limited by the accuracy of the station terms and by weak amplitudes, not by the detection statistic | Further refinement of the statistic, soft timing priors and travel-time-free evidence yield little, and the headroom stays open | **Observed on synthetics, exploratory** (Appendix Y.27): with catalogue corrections wrong by 0.12 to 0.24 s no timing prior gained recall, relocation by polarization or amplitude lost 0.11 to 0.48, and only the true station delays moved recall, by 0.05 to 0.15. Task 2L then confirmed that station terms estimated without picks recover part of that gain in the refined detector, 0.28 to 0.61 of the oracle's, mainly in one-sided eight-station networks and not at four stations, with terms worse than the registered ones in about 22 per cent of the networks and a catalogue prior of known width (G1-C, Appendix Y.28); self-calibration from the detections added nothing measurable in the exploration that preceded it. An exploration that followed found that, at the operating point, the largest class of false atoms under the origin-time rule is the detection of a real event within about 3 km with its origin time off by a median of 0.8 s, which remedies of the deflation do not reduce (Appendix Y.29). The sensitivity to catalogue errors of unknown size is tested on synthetics before real data; on real data, where no oracle exists, residual station terms are reported per station and hits with the origin time matched are a registered endpoint; if calibration does not move recall, the headroom is reported as limited by amplitude. **Observed on real data** (Task R0, Appendix Y.30): after constant station terms the timing residuals of a real aftershock sequence kept a scale of 0.35 to 0.41 s, growing with source separation, the detector missed almost every event and the test was negative. An exploration on the development segment (Appendix Y.32) found that neither a layered nor a fitted medium makes the energy peaks more predictable (a quarter to a half of the readings of the larger events stay beyond 0.25 s), so the planned response before G1' is a reading of arrivals that survives a dense sequence, not a richer medium |
| Radiation operator invalid at stations that fail the polarization audit | The primary detector imposes P motion along the ray and S motion transverse to it at every station, and its advantage with strong clutter depends on that operator (Task 2B) | **Observed on real data for single-station polarization**: the gate passes at one of five RSNC stations (Appendix T.1). G1' carries an arm restricted to stations that pass the orientation and polarization audit and a sign-blind, amplitude-only fallback arm (Section 16); the fraction of stations passing is reported per sub-network, and the audited arm is not run where fewer than two stations pass |
| Synthetic pick model more favourable than a real picker | Comparisons with pick-based baselines biased | **Observed**: the parametric pick model is optimistic for GaMMA in recall and timing, while a calibrated STA/LTA starves GaMMA and PyOcto, so claims against them are stated per picker; with the parametric picks GaMMA is ahead of MT-FW at 0.1 in four of six weak-clutter cells (Appendices Y.14, Y.18). Real-data comparison with PhaseNet picks in their domain (Section 16) |
| False-alarm threshold miscalibrated under non-Gaussian noise | Realized false-alarm rate far above nominal | **Observed on synthetics**: strong clutter raises the realized false-alarm rate of every threshold calibrated on Gaussian noise far above nominal, on average 19 to 22 times for the empirical predictors. Primary threshold by split-conformal or GEV on noise-only windows of the same network and regime, recalibrated when the regime changes; EEC or Bonferroni only as a conservative prior (Appendix Y.17) |
| Real clutter coherent with uncatalogued sources | $\eta$ underestimates coherence | Several catalogues and template-matching screening of clutter picks (Section 16) |
| Horizontal misorientation mistaken for physics | Stations wrongly excluded from C-IV | Orientation audit before the gate (Section 16) |
| Coherence bounds vacuous in the regime of interest | Weak separability theory | Commit to the convex-geometric route; report the plain bound as a weak first certificate only |
| Representation ladder stalls at L0 or L1 | No learning component | Hard caps with abort triggers; association line runs on beams and SGT and is unaffected; negative result reported with diagnostics |
| Beam asymptotics inadequate at the working band | Synthesis error | Tolerance measured against the reference as a function of frequency and heterogeneity; fall back to SGT synthesis |
| Transport misfit cost underestimated | Budget overrun | Measured in the pilot against Table 4; fall back to sliced transport before reducing conditions |
| Parametric likelihood overconfident under theory error | G3 invalid | Explicit model-error covariance plus registered simulation-based comparator; report both |
| No validated velocity model at the working band in a domain | Transfer invalid | Run at the supported band or restrict to envelope and polarization level |
| Scope too large for the schedule | Infeasible | Two-track plan with independent gates (Section 19); theoretical track survives experimental failure and vice versa |

The programme follows graceful degradation: failure of the learned representation, of any single
constraint family, or of any optional physics extension must not invalidate the test of G1.

---

# 18. Reproducibility and Software

Every experiment records dataset version, preprocessing, normalization, source parameterization,
medium generation, architecture, initialization, seeds, collocation rules, optimizer settings,
active constraint sets, continuation schedule, stopping criteria, hardware, software versions,
GPU-hours and partitions. Seed fixing alone is acknowledged insufficient, so library versions,
hardware and determinism settings are recorded too. Fixed seeds must reproduce benchmark statistics
within a declared tolerance.

**Declared stack**, pinned by version and revisable only with documented justification: an adopted
elastodynamic solver, SPECFEM3D_Cartesian or SeisSol, chosen in the second work package after the
Appendix B tests; solver-native SGT output with an in-house contraction layer; an in-house Gaussian
beam synthesis module validated against the reference; ObsPy and FDSN services for waveform I/O;
log-domain stabilized Sinkhorn [35] validated against a public transport library, with sliced
transport as fallback [36]; JAX or PyTorch with double-precision residual evaluation for the learned
component; an in-house primal-dual interior-point solver with IPOPT and an SQP solver as
comparators; neural density estimation for the simulation-based comparator [27]; and versioned
statistical analysis scripts stored with results.

Synthetic realizations are regenerated from seed and configuration rather than archived in bulk;
SGT databases are scratch artefacts with documented regeneration recipes. Negative results remain
part of the record: OOD failures, invalid constraints, failed calibration gates, unsuccessful
initializers, and regimes where a constraint family or the learned representation is not worth its
cost.

---

# 19. Forty-Eight-Month Two-Track Plan

The programme runs two tracks with independent gates. The theoretical track consumes no GPU time
and starts immediately; the experimental track consumes the budget of Section 12. Neither can be
killed by the failure of the other, which is the structural reason for the two-track design rather
than a longer single sequence.

**Table 8. Plan.**

| Months | Theoretical track | Experimental track |
|:--|:--|:--|
| 1-6 | T1: bound the remainder of the Proposition of Section 6.1; extend the certificate to several interfaces, triplications, converted phases and a head-wave amplitude law; exclusion sets and their tolerance (C-II); medium-independent kinematic admissibility (C-III) | E1: after the negative real test R0 (Appendix Y.30) and exploration C18 (Appendix Y.32), the reading of arrivals in a dense real sequence (an envelope reading or a higher-frequency matched filter, since down-weighting the S readings fails, joint relocation of the readings with a coherence check, the logarithmic statistic LOG6), developed only on the development segment of the Ridgecrest hour until a candidate passes the registered gate (at least 5 of the 16 development events with origin times matched at 30 false detections per hour on N8 or U8) or the line is closed and reported; if the gate is passed, a registered test on the unread Hi-net swarm against PhaseNet with PyOcto with an operating point fixed in advance; the first real-data decimation test of G1 on the RSNC only after a positive real test; a public repository with the code, the pinned environment and a registry sealed by a third party; real-data calibration of C-IV on the RSNC (tasks 1c-1n of Section 16); second synthetic generator; **verify compute allocation** |
| 1-12 | T2, started in month 1 as the first mathematical deliverable and the key research question of the mathematical track: the bounded-influence moment-tensor subspace statistic, in its capped form confirmed in Task 2F and in its Huber form with a level fixed by rule from the noise distribution confirmed in Task 2I (Appendices Y.24 to Y.26), with the convex problem of which it would be the dual certificate, a statement on its false-alarm rate under noise (Appendix AA.9 bounds the influence of one station on the capped and logarithmic statistics, not on HUBER-N, and shows that a capped score saturates when the background reaches the cap at every station, the trade-off measured in Task 2M and exploration C18; the false-alarm statement remains open), its reading as an M-estimator, the extension of Appendix AA to it and a cheaper solver for the Huber fit, and the open problem of the headroom, about three quarters of which remains unexplained by any treatment of amplitudes, approached first through station terms, the lever that two explorations identified and that Task 2L confirmed within a narrow scope with terms estimated without picks (G1-C, Appendices Y.27, Y.28), and next through the sensitivity of that estimator to the catalogue, its dependence on geometry and its combination with the bounded-influence statistic; the baselines are the capped statistic (MT-FW) and BOMP-PK, on the dense designs of Tasks 2F and 2I with origin-matched hits as the primary hit rule; amplitude or site terms per station are not pursued, an exploration having found that they do not replace the cap (Appendix Y.25); the off-grid group BLASSO with time shifts did not beat BOMP-PK in Task 2E and is closed as a detector (Appendix Y.22), and the convex method on the grid falls below MT-FW with or without the peak search (Appendices Y.20, Y.21); from month 7, per-station polarization likelihood and subspace conditions (C-IV), source-time-function admissible cone (C-V) and feasibility geometry of the combined set | E2 (months 7-12): G1' by network decimation at full scale (Section 16: both pickers, all comparators including BOMP-PK, MALMI [110] and [112]); the secondary arms of Task 2I, other Huber levels and amplitude spreads, reported as registered, and the bounded-influence statistic, capped and Huber, against BOMP-PK on a second generator family with a held-out operating point and origin-matched hits, registered before it is run; the headroom of Appendix Y.22 studied through the global placement of atoms, with block orthogonal matching pursuit [129] on a finer grid as a further reference; adopt and validate solver; build SGT reference; validate beam synthesis |
| 13-20 | T3: identifiability - Fisher and Cramer-Rao analysis, rank results (19)-(20), coherence and convex-geometric separability | E3: conditional on F1, misfit stack, transport plus continuation, $2\times2$ cycle-skipping study (F2) |
| 18-28 | T4: well-posedness of the constrained problem, convergence of the primal-dual iteration under (15) | E4: constrained attribution, constraint-family ablation, clutter and station-count sweeps (G2, G4, A1, A2; the experiment attributes the effect of G1 to constraints and does not test G1, which is tested in E1 and E2) |
| 26-34 | T5: network design by observability maximization; theory paper consolidation | E5: single-station multi-event study (C3); observability-design experiment (B2, B3) |
| 30-40 | T6: theoretical analysis of the symmetry-exact representation | E6: conditional on F1, representation ladder L0 then L1; symmetry ablation; break-even (F3, F4, F5) |
| 36-44 | T7: uncertainty theory under model error | E7: calibration, simulation-based comparator, transfer on two real domains (G3) |
| 42-48 | T8: synthesis, thesis chapters | E8: final benchmark, software release, reproducibility package |

Gates: allocation confirmed or scope reduced in the declared order by month 6; the arrival-reading gate passed, or the line closed and the scope of the thesis stated as real noise with synthetic signals, by month 6, whatever the outcome; a registered real test on the Hi-net swarm reported by month 9 if the gate is passed; the public repository with a sealed registry in place by month 6; the replication of G1-R on a second generator family, capped and Huber statistics against BOMP-PK with a held-out operating point, reported by month 9, whatever its sign; Appendix B tests 1-8
passed by month 12; L0 by month 34 or the representation line is closed and its budget redirected
to T5 and E5.

T2 runs from month 1 in parallel with T1 because Tasks 2E, 2F and 2I (Appendices Y.22, Y.24, Y.26) placed
the methodological question on it: the time-shift solver did not beat BOMP-PK, the per-station cap did in
dense sequences, an exploration found that amplitude or site terms per station do not replace it
(Appendix Y.25), and a bounded-influence fit with a level fixed by rule from the noise distribution
matched the calibrated cap in its registered test. That result is synthetic, from one generator family
and modest in size, so the track does not treat it as settled: if the replication on a second generator
family fails, the methodological contribution is stated as specific to the generator of Tasks 2A to 2F, 2I and 2L,
with the benchmark and its negative results, and that outcome is reported as such. The headroom, about
three quarters of which remains unexplained by any treatment of amplitudes, stays the open problem of
the track in either case. Two explorations on the calibration seed give it a first direction (Appendix
Y.27): with the true station delays the detector reaches the realistic ceiling for isolated events, while
soft timing priors at catalogue-grade corrections and travel-time-free evidence do not help. Task 2L
took that step with station terms estimated without picks from the catalogue alone, registered with the
true-delay oracle beside it and with origin-time-matched hits as the primary hit rule, and confirmed a
gain that is concentrated in one-sided eight-station networks and reaches 0.28 to 0.61 of the oracle's
(G1-C, Appendix Y.28); self-calibration from the detections added nothing measurable in the exploration
that preceded it. The next synthetic steps are the sensitivity of that estimator to catalogue errors of
unknown size, its dependence on geometry and its combination with the bounded-influence statistic; the
amplitude evidence for weak arrivals comes after them.

Three explorations on the calibration seed that followed Task 2L, with no registered hypothesis and with
the detector of the earlier tasks, bear on those steps (Appendix Y.29, exploratory). Restricting the
moment-tensor statistic to double couples lowers the threshold at a fixed false-alarm rate by about 10
per cent, as a tube formula with 4.4 to 4.5 effective degrees of freedom predicts and simulation
confirms, but recall did not change, $+0.005$ [$-0.021$, $+0.033$] in standard windows, because the hits
lose as much energy to the restriction as the threshold drops; the null law is kept as the calibration of
a test of mechanism, not as a detector. A target-decoy pass with permuted station moveouts reproduced the
operating point of the benchmark without labels, within 0.03 of its recall, and kept the false discovery
rate within a nominal 0.2 under the location-only hit rule, but not under the origin-time rule and not
reliably in the one-sided geometry at lower levels, so a deployable operating point is still not
established. Remedies for the self-interference of the greedy deflation halved the false atoms that the
events induce without moving recall at the operating point under the origin-time rule, and the most
aggressive of them erased about a quarter of the overlapping events found at any score. The diagnosis
behind that negative orders the next steps. At the operating point the largest class of false events
under the origin-time rule is the detection of a real event placed within about 3 km of it with the
origin time off by a median of 0.8 s, mostly early and usually the only detection of that event: 55 to 58
per cent of the false atoms in standard and 80 to 95 per cent in dense windows on the development
networks, and about a third and two thirds on unseen ones. No decoy imitates that class and no remedy of
the deflation reduces it, so the origin time of well-placed detections, with the terms of Task 2L and a
better timing model, is the next target of the calibration work, to be explored before anything is
registered.

Forty-eight months rather than thirty-six is a deliberate consequence of adding a theoretical track
with real content. Compressing it would mean either dropping the mathematical track, which carries
the bounded-influence statistic and the theory of Appendix AA, or claiming a schedule the work does not fit.

---

# 20. Contributions and Publication Strategy

## 20.0 Four contributions, one primary hypothesis

The programme is stated as one primary hypothesis, G1, and three methodological hypotheses, G1-M, refuted,
G1-R, confirmed on synthetics in two preregistered tests, and G1-C, confirmed on synthetics in one, inside four contributions rather than as an enumeration of components, because a list of ten invites the reading that ten separate things must
all succeed. They are the four stated in Section 1.6 and Section 23; the Abstract reports only the methodological one and its evidence.

**Methodological, carrying the primary hypothesis G1, the bounded-influence hypothesis G1-R and the
calibration hypothesis G1-C.**
Timing-tolerant moment-tensor detection with deflation for overlapping earthquakes, applied to the
three-component waveform before picking, thresholded at a false-alarm rate calibrated on noise-only
windows and tested against published associators and migration locators (G1). On these synthetics the
greedy MT-FW and block orthogonal matching pursuit [129] given the same peak search (BOMP-PK) are
equivalent members of that family in standard windows, and the empirical advantage over pick-based
associators and QuakeMigrate belongs to the family (Appendix Y.21). The group BLASSO of Section 7.5 is
the framework, not the best detector: its certified solution on the grid falls below the family with or
without the peak search, and its off-grid form with a time shift per station and phase does not beat
BOMP-PK at equal tolerance, which refutes G1-M in this design (Appendices Y.20 to Y.22). The
methodological deliverable is therefore stated precisely and modestly. G1-R is supported by two
preregistered, independently verified synthetic tests whose registration hashes were sent in the
working session before each run: Task 2F, with a per-station cap at a calibrated level, and Task 2I,
with a Huber fit whose level is fixed by rule from the noise distribution, which replaces the calibrated
cap by a level fixed by rule in the registered primary arm (Appendices Y.24, Y.26). What
they establish is that a bounded-influence moment-tensor subspace statistic improves event recall at a
fixed false-event rate in dense overlapping sequences of three to five events over block orthogonal
matching pursuit with equal timing tolerance, by 0.04 to 0.07 with the cap and by 0.058 and 0.053 with
the Huber fit, without harming isolated events when regimes are pooled. The limits are part of the
deliverable. The evidence is synthetic and from one generator family. The size is modest, and it
shrinks under a strict hit rule that also requires the origin time, to 0.024 to 0.045 for the cap and to
0.031 and 0.019 for the Huber fit. The Huber fit costs about nine times as much as the cap, is slightly
below it at 0.1 false events per window in one seed, and finds no more events when every atom is
accepted, so the gain is in ranking and selecting atoms; its equivalence on isolated events is not met
with weak clutter in one seed; and in its descriptive arms the level is immaterial between 2.5 and 3.4
noise standard deviations, while at an amplitude spread of 0.8 the fit falls behind the cap. Clipping per-station contributions is known practice in robust stacking
[136] and Huber fitting is standard robust statistics, so the novelty is the measured design result
inside a moment-tensor detector with deflation, with its mechanism of bounded influence of each station
and phase block under amplitude model error (Appendix Y.25), and not an algorithm; the same exploration
found that amplitude or site terms per station, estimated or known, do not replace the cap. The other deliverables are the preregistered benchmark, with its negative results as
contributions in their own right: the grid and off-grid convex solvers do not beat block pursuit with
equal tolerance, joint refitting and the cross term of P and S do not matter, and neither the velocity
model nor the source pulse explains the headroom of 0.08 to 0.24 below a realistic oracle ceiling
(Appendices Y.20 to Y.23), of which no treatment of amplitudes captures more than about a quarter
(Appendix Y.25, exploratory); and the test of G1 against published associators, migration locators,
block-sparse moment-tensor inversion [129] with and without the tolerance, and mechanism-corrected
stacking [122], [133], on synthetics and on the RSNC by network decimation (Section 13.0c); the benchmark
and G1' stand as contributions whatever the fate of G1-R on a second generator or on real data. About
three quarters of the measured headroom remains unexplained, and that is the open problem that the
methodological track keeps in view. Two explorations locate it in the accuracy of the station terms and
in the amplitude evidence for weak arrivals, and not in travel-time-free physics or a softer timing
tolerance (Appendix Y.27). G1-C is the first preregistered result on that part: in Task 2L,
independently verified, station terms estimated without picks from the catalogued events raised the
origin-time-matched recall of the capped detector with a local timing refinement by 0.057 and 0.062 in
standard and by 0.038 and 0.030 in dense windows, with no loss in either clutter regime or under the
location-only rule (Appendix Y.28). Its limits are part of the deliverable as well: the gain is
essentially that of one-sided eight-station networks and is not demonstrated with four stations, the
estimated terms are worse than the registered ones in about 22 per cent of the networks, the result is
conditional on the refinement and on a catalogue with errors of known size, it reaches 0.28 to 0.61 of
the gain of the true delays, and joint hypocentre determination with station terms is classical, so
that it is a measured calibration result and not an algorithm. All present evidence
is synthetic (Section 1.6). The comparison with pick-based associators and QuakeMigrate is the part for
which the targeted searches found no 2023 to 2026 precedent. The scheme is evaluated inside a declared
computational budget, with the learned representation as an acceleration layer behind a gate rather
than as a contribution of the thesis.

**Mathematical.** A translation of structural information in elastodynamics into explicit
constraints on a non-convex inverse problem, with the observability certificate established at leading order for a
reduced model and verified numerically in three classes of medium; and the theory of the group BLASSO
detector of Section 7.5 (Appendix AA). Its on-grid counterpart, block-sparse moment-tensor inversion
[129], and its detection statistic, the variance reduction of grid-based moment-tensor monitoring
[132], [123], are prior art; targeted searches found no prior for the off-the-grid formulation over
moment-tensor-valued measures. The theory comprises existence and optimality conditions; the identity
between the certificate at zero and the matched-subspace generalized likelihood ratio, a modest
observation that follows from group-lasso optimality; validity of the pruning bound, which makes the
search of the certificate exact; the Frank-Wolfe rate with an exact oracle, a standard result [124]
applied to this problem; finite-sample validity of the split-conformal threshold, a transfer of
conformal outlier testing [139]; and a lower bound on the curvature of the certificate by the
energy-weighted observability, with exact recovery over the continuum left as a conjecture. The
theory concerns the convex problem and its certificate; since neither its grid solution nor its
off-grid solution with time shifts is the best detector (Appendices Y.20 to Y.22), it is the theory of
the framework, and extending it to the bounded-influence statistic (T2), whose clipped
form is the one element with a confirmed effect, is what would bring it closer to the deployed
detector; Appendix AA.8 bounds the implemented statistic, within its search window, by a
shift-maximized likelihood ratio.

**Scientific.** Physical evidence buys discrimination margin rather than resolution, and its usefulness
is **not uniform**: it depends on the effective kinematic information, set by station count, picker
miss rate and pick timing precision, and on the type of clutter. The object is the boundary

$$\text{arrival-time sufficient}\ \leftrightarrow\ \text{polarization necessary}\ \leftrightarrow\ \text{waveform necessary}$$

as a function of the effective observation count, measured against published software and
replicated across independent seeds. The measured boundary so far is that a continuous kinematic
score suffices when every arrival is picked exactly, and that polarization pays once the picker
misses arrivals, provided the physics is applied to the waveform as a coherent moment-tensor
detector with explaining-away rather than to a kinematic catalogue (Appendix Y.12, Section 1.6);
with weak clutter much of the advantage over QuakeMigrate does not need the radiation operator
(Appendix Y.19). Whether the RSNC lies on that side is the confirmatory question.

**Structural.** The constraints do different mathematical jobs: polarization supplies discrimination
margin but not model independence, kinematic admissibility supplies model independence but not
discrimination, and observability supplies identifiability and a map of where discrimination will be
hard, and does not classify hypotheses (Section 4.3; Appendices J, L, M, S, U). The evidence for this
division is preliminary: eight cells for the ordering of channels, one network for the observability
map, and a model-independent estimator whose localization accuracy is about an order of magnitude
worse.

The component-level claims below expand these four contributions and are subordinate to them and to
G1.

**Component-level claims.**

- **P1. Constraint generation from elastodynamics for source attribution.** Five families -
  observability certificates, certified exclusion from silence, medium-independent kinematic
  admissibility, polarization rank, and source-spectrum bounds - each stated from elastodynamic
  structure (the exclusion from silence is motivated by unique continuation but not derived from it,
  Section 3.3), implemented as an explicit component of a feasible set, and individually ablated so
  that their contributions are attributable rather than aggregate.
- **P2. Silence as a candidate hard-exclusion mechanism under an explicit sensing and noise model.**
  Promotion of non-detection from calibrated likelihood, which is established [17], to a certified
  exclusion set under the seven observation assumptions of Section 6.2, with its erosion tolerance
  under model error measured and a probabilistic fallback retained. The exclusion is hard relative
  to a declared sensing model, not as a consequence of unique continuation alone, and the
  relaxation ladder means the channel degrades to the established likelihood treatment rather than
  failing.
- **P3. Observability as a well-posedness certificate and a design criterion.** Translation of the
  geometric control condition into a coercivity constant for the reduced problem, used both to
  certify identifiability and to place stations, delivered as a Proposition at leading order with four
  corollaries, Eqs. (8a)-(8d), including the rank deficiency of the reduced single-branch model on
  planar and equal-distance idealized geometries and its algebraic removal by two-branch sampling;
  verified numerically at $\rho=+0.9998$, better than the azimuthal-gap heuristic, and different from
  a dilution-of-precision heuristic by 40-90% in any medium with a velocity gradient (Appendix I).
  Whether such geometries are prevalent in deployed networks is registered as C2.
- **P4. Single-station multi-event epicenter determination.** Rank results (19)-(20) and the
  associated method, with the depth ill-conditioning declared rather than hidden.
- **P5, conditional. Structure-preserving differentiable representation.** Analytic beams with a
  learned residual in which linearity, time-shift equivariance, reciprocity and rotational covariance
  hold by architecture rather than by penalty. It is pursued only if gate F1 passes and reported only
  if level L1 passes; it is an option, not a contribution of the thesis.

**Publication units, in priority order; the first four depend neither on the learned representation
nor on the transport misfit:**

1. *Timing-tolerant moment-tensor detection of overlapping earthquakes: a preregistered synthetic
   benchmark against associators, migration locators, block-sparse inversion with and without timing
   tolerance and mechanism-corrected stacking, the negative result for grid and off-grid convex
   solvers, and a bounded-influence statistic for dense sequences, with a calibrated cap and with a
   Huber level fixed by rule from the noise* (Tasks 2A to 2F and 2I, T2, E2, month 14); the secondary
   arms of Task 2I, the calibration result of Task 2L and the replication on a second generator join the unit with their registered
   outcomes, whatever their sign.
2. *Detection at a fixed false-event rate on a decimated national network* (G1' on the RSNC, E1-E2,
   month 20).
3. *Physically constrained epicenter determination: constraint generation and false-event rejection*
   (T1-T2 with E4, month 30).
4. *Observability certificates for seismic source identifiability and network design* (T1, T5, month 34).
5. *Single-station multi-event epicenter determination* (T3, E5, month 38).
6. *Transport misfit and multiscale continuation for waveform-level source attribution* (E3),
   conditional on F1.
7. *Structure-preserving neural representations of elastodynamic propagators* (T6, E6, month 44),
   conditional on F1 and submitted only if L1 passes, otherwise folded into 3 as a negative-result
   section.

Realistic output over 48 months is four submitted papers, units 1 to 4, with unit 5 in preparation;
units 6 and 7 are conditional options, not planned output.

---

# 21. Alignment with the Supervisor's Research Profile

The supervisor is a mathematician with a master's degree in mathematics from Universidad del Valle
and a doctorate in computational science from the University of Texas at El Paso on constrained
optimization schemes for geophysical inversion of seismic data, with published work on a constrained
optimization framework for joint inversion of geophysical data sets, including physically bounded
inverse formulations [72], [73].

The alignment is structural rather than thematic. This proposal keeps the same mathematical
backbone - a constrained inverse problem governed by a physical forward operator, solved by
primal-dual methods with physically motivated bounds - and changes the inverse target from Earth
structure to the source set. Four specific points of contact: Section 6.5 is literally a physically
bounded inverse formulation of the source-time function, which is the supervisor's published line;
Section 6.1 delivers a coercivity constant, which is a well-posedness statement about a constrained
optimization problem; Section 6.2 delivers a feasible-set restriction; and Section 8 is
interior-point optimization with KKT diagnostics and cross-validated sensitivities. The PDE-analytic
content of C-I and C-II is deliberately consumed as constraint generation rather than pursued as
analysis for its own sake, precisely so that the core of the thesis sits inside the supervisor's
competence rather than outside it. Machine learning is one computational layer inside an
applied-mathematics problem, and the project avoids a generic "AI for seismology" identity.

---

# 22. Scientific Integrity

No human subjects or personal data. The integrity risks are scientific: leakage, hidden reuse of
held-out domains, selective reporting, post-hoc threshold tuning and overstated uncertainty.

Freeze partitions before confirmatory testing. Freeze the margins of the constraint-attribution family after the pilot and
those of G1' before any RSNC evaluation, in both cases before confirmatory runs. Pre-register primary metrics, matching rules, non-inferiority margins and
model-comparison rules. Report effect sizes with every significance statement. Preserve negative
results and failed gates, including the six claims abandoned in Section 4.3, two of them refuted by their own pre-studies in
Appendices G and H, the four positives withdrawn after auditing (Appendix Y), the failed rule 2D-a,
which removed the distinctiveness of MT-FW over [129] in standard windows (Appendix Y.21), the failed
rule 2E-a, which refuted G1-M (Appendix Y.22), the null and negative exploratory candidates of Appendix
Y.23, reported with the multiplicity from which G1-R was selected, the selection among 11 candidates
that preceded Task 2I, the null and negative explorations that followed Task 2L (Appendix Y.29), the
negative real-data test R0 on the 2019 Ridgecrest sequence (Appendix Y.30), the null timing and detection
remedies of exploration C18 (Appendix Y.32), and the narrower scope that the independent verifications of
Tasks 2F, 2I, 2L and 2M placed on their own positives (Appendices Y.24, Y.26, Y.28, Y.31). Report a run that was interrupted and resumed as
such, with the times of each resumption and a check that the code did not change. Report arms that are
still running as unfinished and say nothing about their content. State the mechanism of a
confirmed effect only at the scope that an independent verification supports. Use
independent adjoint or analytic references for all learned gradients. Do not infer posterior
uncertainty from optimizer curvature alone. Report training, storage and deployment costs
separately. Re-execute the search protocol of Appendix E before each submission and revise novelty
claims accordingly. Document every deviation from the registered protocol.

---

# 23. Conclusion

**The confirmatory core is defined.** What the thesis defends is one primary hypothesis inside four
contributions (Section 20.0). The primary hypothesis, G1, is that timing-tolerant moment-tensor detection
with deflation, of which MT-FW and block orthogonal matching pursuit given the same peak search are
members, detects more events at a fixed realized false-event rate than published associators and
migration locators where the picker discards kinematic information. The contributions are that detector
family with its bounded-influence statistic and its benchmark; its theory together with a leading-order
observability certificate; a sufficiency ladder in which physics buys discrimination margin rather
than resolution and the physics required is a function of the effective kinematic information; and a
division of labour among constraints, supported by a model-independent feasibility filter. The
evidence for G1 is synthetic and is stated once, with its qualifications, in Section 1.6. Of the eight
preregistered synthetic tasks, four bear on G1 (Tasks 2A to 2D), one on G1-M (Task 2E), two on G1-R
(Tasks 2F and 2I) and one on G1-C (Task 2L). In Tasks 2A to 2C the greedy detector MT-FW exceeded, in point estimate and at 0.25 false
events per window, every comparator run with an operational picker in every cell, although after
Holm's correction the advantage in the strong-clutter four-station cells is not significant over block
orthogonal matching pursuit read at exact times in either seed, nor, in one seed, over QuakeMigrate and
mechanism-corrected stacking (Appendix Y.20). The operating points were chosen on the evaluation
windows, the frozen-threshold rule failed under its strict reading, MT-FW was behind GaMMA fed
optimistic picks at 0.1 in four of six weak-clutter cells, and the band and source domain lie outside the declared
regime of Section 1.5. Task 2D showed that the lead is not MT-FW's own:
given the same peak search, block orthogonal matching pursuit [129] matched it in every cell, the peak
search carried the advantage and the cap did not, and split-conformal thresholds frozen on empty
windows left 1.2 to 2.2 false events per weak-clutter window containing events. Task 2E refuted the
methodological hypothesis G1-M: the group BLASSO with time shifts per station and phase, solved off the
grid by sliding Frank-Wolfe, gained 0.000 to 0.039 over block pursuit with equal tolerance, every
interval containing zero, and an independent verification found no implementation error or budget
handicap. Task 2F, on fresh seeds after an exploration on the calibration seed, confirmed a smaller
methodological claim, G1-R: clipping the contribution of each station in the moment-tensor subspace statistic
raises recall in dense sequences by 0.04 to 0.07, and by 0.024 to 0.045 with origin times matched,
because it bounds the influence of each station and phase block under amplitude model error, not
because it separates events better. Task 2I, on two further fresh seeds after a second exploration,
confirmed the same claim with a Huber fit whose level is fixed by rule from the noise distribution: it
exceeded block pursuit by 0.058 and 0.053 in dense sequences, was not inferior to the calibrated cap
and, pooled, was equivalent to block pursuit on isolated events; with origin times matched the gain is
0.031 and 0.019, the fit costs about nine times as much as the cap, and it ranks and selects atoms
better without finding more events (Appendix Y.26). The empirical advantage over pick-based associators and QuakeMigrate belongs to timing-tolerant
moment-tensor detection with deflation as a family; the convex problem is the framework, the certificate
and the theory rather than the best detector; and the methodological contribution of the thesis is a
measured design result, supported by two preregistered and independently verified synthetic tests on one
generator family: a bounded-influence moment-tensor subspace statistic improves recall at a fixed
false-event rate in dense overlapping sequences over block pursuit with equal timing tolerance, without
harming isolated events when pooled, and in the registered primary arm its level was fixed by rule from
the noise distribution instead of calibrated; its descriptive arms find that level immaterial between
2.5 and 3.4 noise standard deviations and the Huber fit behind the cap at a large amplitude spread. Its size is modest, Huber fitting and per-station clipping are standard, and the benchmark with
its negative results stands beside it; an exploration found that amplitude or site terms per station do
not replace the cap (Appendix Y.25). The headroom of 0.08 to 0.24 below a realistic oracle ceiling, which
neither the velocity model nor the source pulse explains and of which about three quarters remains
unexplained by any treatment of amplitudes, is the open problem. Two explorations indicate that it is a
problem of calibration and of amplitude evidence: the true station delays bring the detector to the
ceiling for isolated events, and travel-time-free evidence is limited by information in this benchmark
(Appendix Y.27). Task 2L, on two more fresh seeds, confirmed a first result on the calibration part,
G1-C: station terms estimated without picks raised the origin-time-matched recall of the refined
detector by 0.057 and 0.062 in standard and by 0.038 and 0.030 in dense windows, a gain concentrated in
one-sided eight-station networks, not demonstrated with four stations, conditional on the refinement and
on catalogue errors of known size, and 0.28 to 0.61 of that of the true delays (Appendix Y.28). Whether the result holds on real data, and against the published
migration locators there, is the question of Section 16. An exploration that followed Task 2L found that,
at the operating point, the largest class of false events under the origin-time rule is the detection of a
real event within a few kilometres with its origin time off by a median of 0.8 s, which points the
calibration work at the origin time of well-placed detections (Appendix Y.29, exploratory).
Everything else is an extension, a gated option, or a registered null. The synthetic programme is not
closed: MALMI [110], the stacking of [112], real noise (Task 1p), three-dimensional Green's functions
(Task 1v) and sources outside the double-couple class remain untested. The real-data decimation test of G1 and the
real-data calibration of the polarization channel (Section 16) are the open work on which the thesis
turns.

The proposal begins from a measurement rather than from a preference: under matched conditions the
wavefield barely sharpens the resolution between nearby true sources and markedly enlarges the margin
by which clutter-supported hypotheses are rejected, and Appendix Y shows that the margin matters for
the catalogue when the picker misses arrivals, not when every arrival is picked exactly. That fixes
the programme: detection at a fixed realized false-event rate is the primary endpoint, the continuous
kinematic likelihood is the reference against which every physical channel must earn its place, and
sparse networks with degraded kinematic information are the target regime.

The organizing idea of the framework is that every piece of physics is a constraint, and clutter is
what falls outside the feasible set. What passed the synthetic tests is not the solution of that
constrained problem but a greedy pursuit with timing tolerance and bounded influence.

The organizing idea is about the feasible set; whether a given constraint removes false events that the
kinematic likelihood does not is an empirical question, answered positively on the consistent
synthetic generator for the coherent moment-tensor detector applied to the waveform, negatively for
polarization used per pick with realistic mechanisms or as a filter on a catalogue, and open on real
data. Against a coherent nearby decoy - a source that
is itself physically realizable - the polarization channel's spurious rate rises to
$0.40\ [0.23,0.59]$ (Appendix T). **Physically realizable impostors are the hardest case for every channel**, and mapping the difficulty from random clutter through
structured clutter to adversarial realizable sources is a primary line of the confirmatory
programme rather than a caveat.

The first measurement on real data (Appendix T.1) points the same way for clutter and adds a
limitation of a different kind. At the one RSNC station where polarization is measurable, real false
picks are not coherent with the preceding catalogued event, $\eta=0.054\ [-0.049,0.164]$; coherence
with uncatalogued sources is not measured by that quantity and is registered work. At four of the five
stations tested, however, single-station polarization falls short of the reliability gate, so the
binding real-data limitation of the polarization channel is that it is available only at some sites,
which makes per-station reliability an input to the constrained problem rather than an assumption of
it; treating unreliable stations as neutral evidence is what keeps the constraint from becoming
harmful (Appendix T.2).

The thesis therefore builds that feasible set explicitly, from
observability, a non-detection exclusion that unique continuation motivates but does not imply,
medium-independent kinematic admissibility, polarization rank and source-spectrum admissibility, and
asks which of those constraints actually carries the information. Positive results would give a
constrained formulation of epicenter determination with attributable components, a well-posedness
certificate that doubles as a network-design criterion, and a route to forming events where current
pipelines cannot. Negative results would delimit, with quantified margins, which constraints are
inert and where the added physics does not justify its cost. Two hypotheses have already been
refuted by their own pre-studies (E1 and E2), three further claims were abandoned, which with G1-M makes the six entries of the abandoned row of Section 4.3,
four positives were withdrawn after auditing, the distinctiveness of MT-FW over [129] in standard
windows did not survive its own preregistered test, and the methodological hypothesis G1-M was refuted
by its own (Appendix Y); all are reported as such, which is the
standard the rest of the programme is expected to meet.

---

# Appendix A. Notation

| Symbol | Meaning |
|:--|:--|
| $N,\ K$ | Records in a window; over-complete candidate budget |
| $\mathbf{w}_i,\mathbf{r}_i,s_i$ | Waveform window, receiver location, station identity |
| $\mathcal{E}$ | Candidate set $\{(\mathbf{x}_k,\tau_k,\mathbf M_k,s_k)\}$ |
| $\mu$ | Source measure, Eq. (11) |
| $a_k$ | Derived activity $\|\mathbf M_k\|_F/M_{\rm ref}$; threshold defines existence |
| $V=(\rho,\lambda,\mu)$ | Medium parameters |
| $M_{\mathrm{ref}},\mathbf M_k$ | Reference moment scale and physical moment-tensor coefficient |
| $G,\ H$ | Green's tensor and receiver-contracted spatial derivative, Eq. (5) |
| $H^{\mathrm{beam}}$ | Analytic Gaussian-beam kernel, Eq. (22) |
| $\mathcal{R}_\phi$ | Learned heterogeneity residual |
| $\phi,\ \boldsymbol{\psi}$ | Helmholtz P and S potentials, Eq. (6) |
| $\mathcal{K}_\phi,\mathcal{K}_\psi$ | Potential coupling terms in heterogeneous media, Eq. (7) |
| $\mathcal{O}_s,\ \Phi,\ W_i$ | Station observation, feature and whitening operators |
| $\mathcal{C},\ \mathcal{C}_0(\mathbf d)$ | Fixed station-time partition; data-derived silent cells under a pre-registered detector |
| $\mathcal{X}_c,\ \mathcal{X}_c^{\eta}$ | Exclusion set of cell $c$ and its erosion by tolerance $\eta$, Eq. (9) |
| $\gamma,\ \gamma_{\min}$ | Observability-derived coercivity constant and its floor, Eq. (8) |
| $\Sigma_s(t)$ | Instantaneous three-component covariance; its rank bounds coherent wavefronts |
| $P_{\mathrm{pred}},P_{\mathrm{obs}}$ | Orthogonal projectors onto predicted and observed polarization subspaces |
| $s_k,\ \mathcal{S}_{\mathrm{adm}}$ | Normalized source-time function and its admissible set, Eq. (10) |
| $\varrho,\ \xi_c$ | Energy-admissibility margin and cell slack, Eq. (16) |
| $\mathcal{J}_{\mathrm{OT}},\ \mathcal{F}_m$ | Transport misfit; frequency band of continuation stage $m$ |
| $\mathcal{I}(\Theta),\ \Delta_{\mathrm{res}},\ \kappa$ | Fisher information; resolution limit and its constant |
| $S^\ast_{b,\mathrm{reduced}},\ S^\ast_{b,\mathrm{full}}$ | Reduced and full-model station-count thresholds for identifiability |
| $\mu(k,k')$ | Mutual coherence of two source components |
| L0-L3 | Viability levels of the representation ladder |

---

# Appendix B. Acceptance Tests

1. Solver reproduces analytic and benchmark solutions; converges under mesh and time refinement.
2. Moment-tensor parameterization reproduces expected radiation patterns, polarity and symmetry.
3. Station records converge as source mollifier width tends to zero.
4. Free-surface and PML treatments agree with the reference within tolerance.
5. Direct multi-source simulation agrees with superposition of single-source solutions.
6. SGT-assembled seismograms agree with direct simulation on held-out source positions.
7. Beam synthesis agrees with the reference within a declared asymptotic tolerance, reported as a
   function of frequency and heterogeneity.
8. Helmholtz coupling magnitude $\|\mathcal{K}\|$ measured on the generator and reported; residual
   always evaluated on the reconstructed vector field.
9. Observability constant $\gamma$ computed for every geometry used, and reported with the results.
10. Exclusion sets: erosion tolerance calibrated so no pilot true source is excluded; exclusion
    rate of true sources reported for every confirmatory condition.
11. Kinematic admissibility: the Lipschitz and ordering tests accept every uncontaminated synthetic truth and reject a declared fraction of contaminated pick sets.
12. Polarization: instrument-orientation error and inter-component amplitude fidelity measured
    before C-IV is enabled.
13. Source-cone constraints verified to reject deliberately non-physical source-time functions.
14. Transport misfit monotone in a controlled pure-shift experiment where the time-domain misfit is
    verified non-monotone beyond a half period.
15. Measured wall-clock per objective evaluation within a factor of two of Table 4.
16. Continuation: stage-$m$ solutions remain in the correct basin when passed to stage $m+1$.
17. Representation field and gradient errors within the F3 tolerances at the declared level.
18. Symmetry exactness: reciprocity and SO(3) covariance verified with random untrained weights.
19. Factorization ablation reported: unfactorized network trained under identical budget.
20. Association equivalence: replacing reference synthesis by the learned representation does not
    change conclusions beyond the F3 end-to-end tolerance.
21. Energy admissibility: duplicated near-identical candidates violate (15) as expected; the truth
    remains feasible at the calibrated margin in at least 95% of pilot cells.
22. Leakage: no held-out station, date or interval used for tuning; access timestamps stored.
23. Budget: measured GPU-hours and storage within the Section 12 declaration; overruns reported.
24. Reproducibility: fixed seeds and declared environments reproduce benchmark statistics within
    tolerance.
25. Negative-result gate: failure of the representation, of a constraint family, or of the
    single-station line never triggers silent substitution by an unvalidated model.

---

# Appendix C. Experimental Matrix

| Code | Evidence and constraints | Synthesis | Primary question |
|:--|:--|:--|:--|
| TT | Arrival times only | travel-time model | Counterfactual |
| COH | Waveform coherence | back-projection, matched filter | Is constrained inference necessary? |
| SIG | Statistical waveform model | SIGVISA | Physics versus historical statistics |
| PH-ALL | All five constraints | beams or SGT | Full model |
| PH-noI | C-II to C-V | beams or SGT | Value of observability restriction |
| PH-noII | C-I, C-III to C-V | beams or SGT | Value of silence as exclusion |
| PH-noIII | C-I, C-II, C-IV, C-V | beams or SGT | Value of kinematic admissibility |
| PH-noIV | C-I to C-III, C-V | beams or SGT | Value of polarization rank |
| PH-noV | C-I to C-IV | beams or SGT | Value of source-spectrum bounds |
| PH-MIX | All, exclusive assignment | beams or SGT | Value of additivity |
| PH-LR | All | beams + learned residual | Does the representation preserve conclusions? |
| PH-1ST | All, one station | beams | Single-station multi-event feasibility |
| PH-OT | All, transport misfit ladder | beams | Cycle-skipping control |

---

# Appendix D. Pre-Registered Thresholds

Registered defaults. Those of the constraint-attribution family are frozen after the pilot, which does
not test G1; those of G1, G1-M, G1-R and G1-C were fixed in the registration of each synthetic task, and that
of G1' is fixed before any RSNC evaluation. All precede any confirmatory run. Revisions are
documented with date, justification, and confirmation that they preceded confirmatory testing.

| Quantity | Value |
|:--|:--|
| Recall non-inferiority margin | 0.02 absolute |
| G1, the primary endpoint of the thesis | $\ge0.05$ absolute gain in event recall at realized rates of $0.1$ and $0.25$ false events per window on the exact FROC curve (synthetic form registered as rule 2A-a, Appendix Y.18, and extended to the closest antecedents and to thirty networks per cell with a hierarchical interval as rules 2C-a and 2C-c, Appendix Y.20; rule 2D-a, which compared MT-FW with the other timing-tolerant members of the family and failed, Appendix Y.21, makes BOMP-PK the baseline of G1-M and does not test G1; real-data form G1' below). False-association reduction at non-inferior recall, the endpoint of A2, F1, G4 and Section 13.4, is secondary |
| G1-M, the methodological hypothesis | $\ge0.05$ higher event recall than BOMP-PK at $0.25$ realized false events per window on the exact FROC curve, at equal timing tolerance and with the same station corrections, on the windows of Task 2A and thirty networks per cell, decided on the lower bound of the hierarchical interval; registered as rule 2E-a, which failed in both regimes (Appendix Y.22) |
| G1-R, the bounded-influence hypothesis | MT-FW+SC minus BOMP-PK+SC $\ge0.03$ at $0.25$ realized false events per window on the exact FROC curve, with the lower bound of the hierarchical interval $>0$, in dense windows of three to five events pooled over cells, per clutter regime, in two fresh seeds; the gain attributed to the cap if MT-FW without it is within $0.02$ of BOMP-PK; equivalence within $0.05$ in standard windows (rules 2F-a to 2F-c, Appendix Y.24). For the Huber statistic with the level $2.7955$, in two further fresh seeds and pooled over regimes and cells: HUBER-N+SC minus BOMP-PK+SC $\ge0.03$ in dense windows with a one-sided 97.5 per cent lower bound $>0$; HUBER-N+SC minus MT-FW+SC with a lower bound $>-0.03$ in dense windows; and the 90 per cent interval of HUBER-N+SC minus BOMP-PK+SC inside $[-0.05,+0.05]$ in standard windows (rules 2I-a to 2I-c, Appendix Y.26), with the hit rule that also matches the origin time as a secondary reading; all six rules passed. The threshold is lower than that of G1 and is read at $0.25$ only, because G1-R compares two members of one family that differ in the bound alone and was planned on half a pilot effect of about $0.07$; its equivalence margin of $0.05$ exceeds that threshold, and tasks registered after Task 2I use the origin-time-matched hit rule as primary and a margin of $0.02$ per regime |
| G1-C, the calibration hypothesis | MT-FW with the cap of 6 and the local timing refinement, with station terms from the pick-free joint determination, minus the same detector with the registered catalogue corrections, at $0.25$ realized false events per window on the exact FROC curve, pooled over cells, in two fresh seeds: with the origin-time-matched hit rule, $\ge0.03$ in standard windows with a one-sided 97.5 per cent lower bound $>0$, and a one-sided 97.5 per cent lower bound $>0$ in dense windows; a one-sided 95 per cent lower bound $>-0.02$ per clutter regime with that hit rule, and pooled under the location-only rule, in standard and in dense windows (rules 2L-a to 2L-d, Appendix Y.28); all four rules passed. Superiority under the location-only rule is not registered |
| G2 margin degradation ratio; false-score ceiling | $\le0.5$ of the travel-time rate; $\le0.5$ of true score |
| A1 per-constraint effect size | $|\hat{A}_{12}-0.5|\ge0.075$ at FDR 0.05 |
| G4 false-source reduction; true-source exclusion | $\ge0.03$; zero outside tolerance $\eta$ |
| B3 certificate correlation | Spearman $|\rho|\ge0.5$ |
| E1 median horizontal improvement (registered null) | $\ge15\%$ |
| C3 endpoints | epicentral position, origin time, source count; depth with CRB only |
| F2 basin advantage | $\ge20$ percentage points, $\ge200$ multi-starts per cell |
| F3 field / receiver / gradient / end-to-end | $\le3\%$ / $\le5\%$ / $\ge0.99$ and $\le5\%$ / $\le0.01$ F1 |
| F4 symmetry exactness; budget | machine precision at random weights; $\le50\%$ of penalty baseline |
| G3 coverage error at 90% | $\le0.05$ |
| F5 break-even / storage | $N_{\mathrm{break}}\le10^4$ or $\ge10\times$ storage reduction |
| Resolution constant $\kappa$ | measured 0.105-0.133 at SNR 20; re-estimated per geometry |
| Minimum station count | $S^\ast_{2,\mathrm{reduced}}=1$, $S^\ast_{1,\mathrm{reduced}}=2$ in Appendix G; for the full four-parameter arrival model $S^\ast_{2,\mathrm{full}}=3$ under Section 9.2; both re-tested in 3-D |
| Clutter levels | 0, 15, 30, 60 false arrivals per window; stress to 200 |
| Station counts | 1, 2, 4, 12, 30 in the Tier-2 sweep of Section 13.4; 30 in the Tier-1 pilot, which does not test G1 (Section 13.0); 3, 4 and 6 in the G1' decimation test (Section 16) |
| Velocity-model error levels | 0, 1, 2, 4, 8 per cent RMS |
| Energy margin $\varrho$ | default 0.5; ablated over 0.2, 0.5, 1.0, unconstrained |
| Exclusion erosion $\eta$ | calibrated on the pilot; no true source excluded |
| Reference moment scale $M_{\mathrm{ref}}$ | fixed before confirmatory testing; used only to nondimensionalize activity and sparsity |
| Activity threshold $a_{\min}$; merge threshold | calibrated in the pilot, then frozen before confirmatory testing |
| Observability floor $\gamma_{\min}$ | set from the pilot; source volume restricted accordingly |
| GOSPA | $\sigma_x=2$ km, $\sigma_\tau=1$ s, $c=2.5$, $p=2$, $\alpha=2$ |
| Multiplicity; bootstrap | Holm's correction within each declared family of tests on $E_{\mathrm{det}}$ (G1, G1-M, G1-R, G1-C, G1'); Benjamini-Hochberg FDR 0.05 for the secondary endpoints; BCa with $B\ge10^4$ |
| A1 ordering stability | rank correlation $\ge0.8$ over $\ge20$ cells; preliminary eight-cell result reported separately |
| A2 minimal sufficient set | registered pair C-III+C-IV within $\delta=0.02$ of full stack in $\ge80\%$ of confirmatory cells |
| A3 model-free construction | report rejection robustness under model error and absolute location-accuracy cost; location improvement is not a requirement |
| A4 cost-effectiveness | $\ge90\%$ of benefit at $<1\%$ of cost |
| B1 interfaces | $\rho\ge0.95$ in the tested branch-summed two-layer interface model; triplications remain confirmatory |
| B2 design advantage | $\ge0.20$ reduction in false-event rate at equal station count (false-event-rate endpoint from median epicentral error, which saturates; Appendix N) |
| B3 local prediction | Spearman $\ge0.5$ between local $\gamma$ and local margin |
| C1 crossover | report $S^\ast$ with confidence interval |
| C2 failure class | $\ge3$ real networks identified and collapse reproduced |
| D1 coupling threshold | report picker miss rate at crossover, with confidence interval |
| F1 gate | $\ge0.03$ false-association reduction beyond the minimal set, else branch closed |
| G1' (continuous-data form of G1) | $\ge0.05$ higher recall of events above $M_c$ for the timing-tolerant detectors of the pick-free arm, MT-FW and BOMP-PK, each tested with Holm's correction over the two and G1' holding for the family if it holds for either, than the best comparator (tuned GaMMA, PyOcto and REAL fed the same PhaseNet picks, and with the calibrated STA/LTA as second picker; the azimuth-likelihood baseline; the continuous kinematic likelihood with and without (9d); the mechanism-independent migration; the stacking of [112]; QuakeMigrate and MALMI), at realized false-event rates per day counted over all continuous data, events included, at the two rates per day fixed before any RSNC evaluation, after template-matching screening of unmatched detections, in at least two of the three decimation levels of three, four and six stations, with a hierarchical bootstrap interval over sub-networks and then days excluding zero; the sample size follows from (25) with $\sigma_d$ measured at those station counts (Sections 13.0, 16), and each detector is also reported against its geometry-breaking control. Two further arms are reported with it, one whose moment-tensor subspace uses only stations that pass the orientation and polarization audit and a sign-blind, amplitude-only fallback (Section 16); the decision is taken on the full-operator arm. (Appendix Y.3 used an earlier form: fewer false events at a fixed recall target) |
| Tie control (Appendix Y) | if the best tie-breaking ablation without polarization is within $0.05$ of PHYS in false-event rate in both one-sided cells and both seeds, the one-sided advantage is attributed to tie-breaking |
| Polarization reliability $\hat r_s$ | estimated by (9e) after the orientation audit; stations enter with weight $\hat r_s$, no binary exclusion |
| Confirmatory $n$ | from (25) with pilot $\sigma_d$; provisional target 100 |
| Representation ladder | L1 primary; L0 abort $>10\%$ at 50% of 20 GPU-h; L1 abort $>6\%$ at 50% of 150 GPU-h per medium |
| Frame sparsity and cost targets | planning figures only; replaced by the first-track measurement |

---

# Appendix E. Literature-Search Protocol

Databases: Scopus, Web of Science, IEEE Xplore, AGU/Wiley, Oxford/GJI, arXiv (physics.geo-ph,
cs.LG, math.OC, math.AP), with citation chaining. Date range 2016 to search date, classical
references without date limit. Query families: (i) *phase association* AND (*Bayesian* OR *mixture*
OR *optimal transport* OR *graph* OR *waveform*); (ii) (*physics-informed neural network* OR *neural
operator*) AND (*elastic wave* OR *elastodynamic*) AND (*multi-source* OR *source encoding*);
(iii) (*observability* OR *unique continuation* OR *geometric control condition*) AND (*wave
equation*) AND (*inverse source*); (iv) (*multi-object* OR *random finite set*) AND (*seismic*);
(v) (*sparse spike* OR *super-resolution* OR *demixing*) AND (*source localization*);
(vi) *non-detection* OR *missed detection* AND *seismic monitoring*; (vii) (*equivariant* OR
*SO(3)*) AND (*seismic* OR *three-component*); (viii) (*moment tensor*) AND (*time shift* OR *station
correction* OR *static*) AND (*detection* OR *inversion* OR *sparse*). Inclusion: peer-reviewed or preprint work presenting
a method, formulation or benchmark relevant to at least two of waveform-level association,
multi-source inversion, constraint generation from wave physics, differentiable wave representation,
detection modelling. Exclusion: observational catalogues without method, non-transferable
applications, duplicate versions. Record counts at identification, screening, eligibility and
inclusion are exported to the repository. The protocol has not yet been executed with recorded counts:
the absences reported in Sections 3.2 and 3.6 come from targeted searches, and its first full execution,
including family (viii), is a deliverable of the first quarter. Re-executed before each submission and before any novelty
claim is committed to writing.

---

# Appendix F. References

References are presented in numeric citation order. Preprints are identified as such, and author names
are reported only to the extent supported by the cited record.

[1] S. Stein and M. Wysession, *An Introduction to Seismology, Earthquakes, and Earth Structure*. Malden, MA: Blackwell, 2003.

[2] K. Aki and P. G. Richards, *Quantitative Seismology*, 2nd ed. Sausalito, CA: University Science Books, 2002.

[3] W. Zhu, I. W. McBrearty, S. M. Mousavi, W. L. Ellsworth, and G. C. Beroza, "Earthquake phase association using a Bayesian Gaussian mixture model," *J. Geophys. Res. Solid Earth*, vol. 127, e2021JB023249, 2022, doi: 10.1029/2021JB023249.

[4] Z. E. Ross, W. Zhu, and K. Azizzadenesheli, "Neural mixture model association of seismic phases," arXiv:2301.02597, 2023.

[5] Z. E. Ross, Y. Yue, M.-A. Meier, E. Hauksson, and T. H. Heaton, "PhaseLink: A deep learning approach to seismic phase association," *J. Geophys. Res. Solid Earth*, vol. 124, no. 1, pp. 856-869, 2019, doi: 10.1029/2018JB016674.

[6] I. W. McBrearty and G. C. Beroza, "Earthquake phase association with graph neural networks," *Bull. Seismol. Soc. Am.*, vol. 113, no. 2, pp. 524-547, 2023, doi: 10.1785/0120220182.

[7] J. Munchmeyer, "PyOcto: A high-throughput seismic phase associator," *Seismica*, vol. 3, no. 1, 2024, doi: 10.26443/seismica.v3i1.1130.

[8] J. D. Smith, K. Azizzadenesheli, and Z. E. Ross, "EikoNet: Solving the Eikonal equation with deep neural networks," *IEEE Trans. Geosci. Remote Sens.*, vol. 59, no. 12, pp. 10685-10696, 2021, doi: 10.1109/TGRS.2020.3039165.

[9] C. Shi, G. Poggiali, C. Marone, M. V. de Hoop, and I. Dokmanic, "High-rate phase association with travel time neural fields," *Nature Communications*, vol. 17, art. 8140, 2026, doi: 10.1038/s41467-026-74092-y; preprint arXiv:2307.07572.

[10] H. Sun, "Microseismic monitoring with the quake neural operator," *Nature Communications*, vol. 17, art. 5677, 2026, doi: 10.1038/s41467-026-73965-6.

[11] W. Zhu and G. C. Beroza, "PhaseNet: A deep-neural-network-based seismic arrival-time picking method," *Geophys. J. Int.*, vol. 216, no. 1, pp. 261-273, 2019.

[12] S. M. Mousavi, W. L. Ellsworth, W. Zhu, L. Y. Chuang, and G. C. Beroza, "Earthquake transformer," *Nature Communications*, vol. 11, art. 3952, 2020.

[13] H. Sun, Z. E. Ross, W. Zhu, and K. Azizzadenesheli, "Phase neural operator for multi-station picking of seismic arrivals," *Geophys. Res. Lett.*, vol. 50, e2023GL106434, 2023.

[14] H. Kao and S.-J. Shan, "The source-scanning algorithm," *Geophys. J. Int.*, vol. 157, no. 2, pp. 589-594, 2004.

[15] S. J. Gibbons and F. Ringdal, "The detection of low magnitude seismic events using array-based waveform correlation," *Geophys. J. Int.*, vol. 165, no. 1, pp. 149-166, 2006.

[16] J. Puente Huerta, C. Sippl, J. Munchmeyer, and I. W. McBrearty, "Benchmarking seismic phase associators: Insights from synthetic scenarios," *Seismica*, vol. 4, no. 2, 2025, doi: 10.26443/seismica.v4i2.1559.

[17] N. S. Arora, S. Russell, and E. B. Sudderth, "NET-VISA: Network processing vertically integrated seismic analysis," *Bull. Seismol. Soc. Am.*, vol. 103, no. 2A, pp. 709-729, 2013, doi: 10.1785/0120120107.

[18] D. A. Moore and S. J. Russell, "Signal-based Bayesian seismic monitoring," in *Proc. AISTATS*, PMLR vol. 54, pp. 1293-1301, 2017; see also D. A. Moore, Ph.D. dissertation, Univ. of California, Berkeley, 2016.

[19] J. Song et al., "VORA: Rapid association of earthquake phases from local to global," preprint, 2026, arXiv:2607.10450.

[20] X. Si, X. Wu, Z. Li, S. Wang, and J. Zhu, "An all-in-one seismic phase picking, location, and association network for multi-task multi-station earthquake monitoring," *Communications Earth & Environment*, vol. 5, no. 1, art. 22, 2024, doi: 10.1038/s43247-023-01188-4.

[21] J. Zhang, H. Zhu, Z. Li, X. Wu, and J. Zhang, "Multistation seismic location via machine learning," *Geophys. J. Int.*, vol. 241, no. 3, pp. 1852-1866, 2025.

[22] A. Novoselov, P. Balazs, and G. Bokelmann, "SEDENOSS: Separating and denoising seismic signals with dual-path recurrent neural network architecture," *J. Geophys. Res. Solid Earth*, vol. 127, no. 3, e2021JB023183, 2022.

[23] S. Cohen, D. M. Steinberg, Y. Radzyner, and Y. Ben Horin, "Expert-guided class-conditional goodness-of-fit scores for interpretable classification with informative missingness: an application to seismic monitoring," preprint, 2026, arXiv:2604.14809.

[24] D. Schorlemmer and J. Woessner, "Probability of detecting an earthquake," *Bull. Seismol. Soc. Am.*, vol. 98, no. 5, pp. 2103-2117, 2008.

[25] Y. Ogata, "Statistical models for earthquake occurrences and residual analysis for point processes," *J. Amer. Statist. Assoc.*, vol. 83, no. 401, pp. 9-27, 1988.

[26] Z. Duputel, L. Rivera, Y. Fukahata, and H. Kanamori, "Uncertainty estimations for seismic source inversions," *Geophys. J. Int.*, vol. 190, no. 2, pp. 1243-1256, 2012, doi: 10.1111/j.1365-246X.2012.05554.x.

[27] A. A. Saoulis, D. Piras, A. Spurio Mancini, B. Joachimi, and A. M. G. Ferreira, "Full-waveform earthquake source inversion using simulation-based inference," *Geophys. J. Int.*, vol. 241, no. 3, pp. 1740-1761, 2025, doi: 10.1093/gji/ggaf112.

[28] A. A. Saoulis, T.-S. Pham, and A. M. G. Ferreira, "Improving moment tensor solutions under Earth structure uncertainty with simulation-based inference," preprint, 2026, arXiv:2603.18925.

[29] J. Virieux and S. Operto, "An overview of full-waveform inversion in exploration geophysics," *Geophysics*, vol. 74, no. 6, pp. WCC1-WCC26, 2009.

[30] C. Bunks, F. M. Saleck, S. Zaleski, and G. Chavent, "Multiscale seismic waveform inversion," *Geophysics*, vol. 60, no. 5, pp. 1457-1473, 1995.

[31] B. Engquist and B. D. Froese, "Application of the Wasserstein metric to seismic signals," *Commun. Math. Sci.*, vol. 12, no. 5, pp. 979-988, 2014.

[32] L. Metivier, R. Brossier, Q. Merigot, E. Oudet, and J. Virieux, "Measuring the misfit between seismograms using an optimal transport distance," *Geophys. J. Int.*, vol. 205, no. 1, pp. 345-377, 2016.

[33] L. Metivier, R. Brossier, Q. Merigot, and E. Oudet, "A graph space optimal transport distance as a generalization of $L^p$ distances," *Inverse Problems*, vol. 35, no. 8, 085001, 2019.

[34] M. Cuturi, "Sinkhorn distances: Lightspeed computation of optimal transport," in *Proc. NeurIPS*, 2013.

[35] B. Schmitzer, "Stabilized sparse scaling algorithms for entropy regularized transport problems," *SIAM J. Sci. Comput.*, vol. 41, no. 3, pp. A1443-A1481, 2019.

[36] N. Bonneel, J. Rabin, G. Peyre, and H. Pfister, "Sliced and Radon Wasserstein barycenters of measures," *J. Math. Imaging Vis.*, vol. 51, no. 1, pp. 22-45, 2015.

[37] D. Komatitsch and J. Tromp, "Introduction to the spectral element method for three-dimensional seismic wave propagation," *Geophys. J. Int.*, vol. 139, no. 3, pp. 806-822, 1999.

[38] J. Tromp, D. Komatitsch, and Q. Liu, "Spectral-element and adjoint methods in seismology," *Commun. Comput. Phys.*, vol. 3, no. 1, pp. 1-32, 2008.

[39] M. van Driel, L. Krischer, S. C. Stahler, K. Hosseini, and T. Nissen-Meyer, "Instaseis: Instant global seismograms based on a broadband waveform database," *Solid Earth*, vol. 6, pp. 701-717, 2015.

[40] J. Tromp, C. Tape, and Q. Liu, "Seismic tomography, adjoint methods, time reversal and banana-doughnut kernels," *Geophys. J. Int.*, vol. 160, no. 1, pp. 195-216, 2005.

[41] V. Cerveny, M. M. Popov, and I. Psencik, "Computation of wave fields in inhomogeneous media - Gaussian beam approach," *Geophys. J. R. Astron. Soc.*, vol. 70, no. 1, pp. 109-128, 1982.

[42] H. Zhu, "Elastic wavefield separation based on the Helmholtz decomposition," *Geophysics*, vol. 82, no. 2, pp. S173-S183, 2017.

[43] C. Bardos, G. Lebeau, and J. Rauch, "Sharp sufficient conditions for the observation, control, and stabilization of waves from the boundary," *SIAM J. Control Optim.*, vol. 30, no. 5, pp. 1024-1065, 1992.

[44] V. Isakov, *Inverse Problems for Partial Differential Equations*, 3rd ed. Cham: Springer, 2017.

[45] J.-L. Lions, *Controlabilite exacte, perturbations et stabilisation de systemes distribues*. Paris: Masson, 1988.

[46] E. J. Candes and C. Fernandez-Granda, "Towards a mathematical theory of super-resolution," *Comm. Pure Appl. Math.*, vol. 67, no. 6, pp. 906-956, 2014.

[47] V. Duval and G. Peyre, "Exact support recovery for sparse spikes deconvolution," *Found. Comput. Math.*, vol. 15, pp. 1315-1355, 2015.

[48] D. L. Donoho and X. Huo, "Uncertainty principles and ideal atomic decomposition," *IEEE Trans. Inf. Theory*, vol. 47, no. 7, pp. 2845-2862, 2001.

[49] M. B. McCoy and J. A. Tropp, "Sharp recovery bounds for convex demixing, with applications," *Found. Comput. Math.*, vol. 14, pp. 503-567, 2014.

[50] D. L. Donoho and J. Tanner, "Observed universality of phase transitions in high-dimensional geometry," *Phil. Trans. R. Soc. A*, vol. 367, pp. 4273-4293, 2009.

[51] D. Amelunxen, M. Lotz, M. B. McCoy, and J. A. Tropp, "Living on the edge: Phase transitions in convex programs with random data," *Inf. Inference*, vol. 3, no. 3, pp. 224-294, 2014.

[52] R. P. S. Mahler, *Statistical Multisource-Multitarget Information Fusion*. Norwood, MA: Artech House, 2007.

[53] B.-N. Vo and W.-K. Ma, "The Gaussian mixture probability hypothesis density filter," *IEEE Trans. Signal Process.*, vol. 54, no. 11, pp. 4091-4104, 2006.

[54] A. S. Rahmathullah, A. F. Garcia-Fernandez, and L. Svensson, "Generalized optimal sub-pattern assignment metric," in *Proc. 20th Int. Conf. Information Fusion*, 2017.

[55] M. Raissi, P. Perdikaris, and G. E. Karniadakis, "Physics-informed neural networks," *J. Comput. Phys.*, vol. 378, pp. 686-707, 2019.

[56] M. Rasht-Behesht, C. Huber, K. Shukla, and G. E. Karniadakis, "Physics-informed neural networks (PINNs) for wave propagation and full waveform inversions," *J. Geophys. Res. Solid Earth*, vol. 127, e2021JB023120, 2022.

[57] P. Ren, C. Rao, S. Chen, J.-X. Wang, H. Sun, and Y. Liu, "SeismicNet: Physics-informed neural networks for seismic wave modeling in semi-infinite domain," *Comput. Phys. Commun.*, vol. 295, 109010, 2024.

[58] C. Zou, K. Azizzadenesheli, Z. E. Ross, and R. W. Clayton, "Deep neural Helmholtz operators for 3-D elastic wave propagation and inversion," *Geophys. J. Int.*, vol. 239, no. 3, pp. 1469-1484, 2024.


[59] K. Luo et al., "Physics-informed neural networks for PDE problems: A comprehensive review," *Artif. Intell. Rev.*, vol. 58, art. 323, 2025.

[60] A. S. Krishnapriyan, A. Gholami, S. Zhe, R. M. Kirby, and M. W. Mahoney, "Characterizing possible failure modes in physics-informed neural networks," in *Proc. NeurIPS*, 2021.

[61] S. Wang, X. Yu, and P. Perdikaris, "When and why PINNs fail to train: A neural tangent kernel perspective," *J. Comput. Phys.*, vol. 449, 110768, 2022.

[62] N. Rahaman et al., "On the spectral bias of neural networks," in *Proc. ICML*, 2019.

[63] M. Tancik et al., "Fourier features let networks learn high frequency functions in low dimensional domains," in *Proc. NeurIPS*, 2020.

[64] M. H. Taufik and T. Alkhalifah, "LatentPINNs: Generative physics-informed neural networks via a latent representation learning," *Artif. Intell. Geosci.*, vol. 6, no. 1, 100115, 2025.

[65] S. Cheng and T. Alkhalifah, "Multifrequency wavefield solutions for variable velocity models using meta-learning enhanced low-rank physics-informed neural network," *Geophys. J. Int.*, 2026, doi: 10.1093/gji/ggaf417.

[66] X. Ma and T. Alkhalifah, "An effective physics-informed neural operator framework for predicting wavefields," *J. Geophys. Res. Mach. Learn. Comput.*, vol. 3, no. 2, e2025JH000899, 2026.

[67] T. S. Cohen and M. Welling, "Group equivariant convolutional networks," in *Proc. ICML*, 2016.

[68] C. Deng, O. Litany, Y. Duan, A. Poulenard, A. Tagliasacchi, and L. Guibas, "Vector neurons: A general framework for SO(3)-equivariant networks," in *Proc. ICCV*, 2021.

[69] F. B. Fuchs, D. E. Worrall, V. Fischer, and M. Welling, "SE(3)-Transformers: 3D roto-translation equivariant attention networks," in *Proc. NeurIPS*, 2020.

[70] J. Janquart et al., "Analyses of overlapping gravitational wave signals using hierarchical subtraction and joint parameter estimation," *Mon. Not. R. Astron. Soc.*, vol. 523, no. 2, pp. 1699-1710, 2023.

[71] J. Langendorff, A. Kolmus, J. Janquart, and C. Van Den Broeck, "Normalizing flows as an avenue to studying overlapping gravitational wave signals," *Phys. Rev. Lett.*, vol. 130, no. 17, 171402, 2023.

[72] U. A. Sosa Aguirre, "Constrained optimization schemes for geophysical inversion of seismic data," Ph.D. dissertation, Computational Science, Univ. of Texas at El Paso, 2012.

[73] A. Sosa, A. A. Velasco, L. Velazquez, M. Argaez, and R. A. Romero, "Constrained optimization framework for joint inversion of geophysical data sets," *Geophys. J. Int.*, vol. 195, no. 3, pp. 1745-1762, 2013.

[74] J. Nocedal and S. J. Wright, *Numerical Optimization*, 2nd ed. New York: Springer, 2006.

[75] S. Boyd and L. Vandenberghe, *Convex Optimization*. Cambridge: Cambridge Univ. Press, 2004.

[76] A. G. Baydin, B. A. Pearlmutter, A. A. Radul, and J. M. Siskind, "Automatic differentiation in machine learning: A survey," *J. Mach. Learn. Res.*, vol. 18, no. 153, pp. 1-43, 2018.

[77] B. Jones and C. J. Nachtsheim, "A class of three-level designs for definitive screening in the presence of second-order effects," *J. Qual. Technol.*, vol. 43, no. 1, pp. 1-15, 2011.

[78] J. Demsar, "Statistical comparisons of classifiers over multiple data sets," *J. Mach. Learn. Res.*, vol. 7, pp. 1-30, 2006.

[79] A. Vargha and H. D. Delaney, "A critique and improvement of the CL common language effect size statistic," *J. Educ. Behav. Stat.*, vol. 25, no. 2, pp. 101-132, 2000.

[80] K. Chaloner and I. Verdinelli, "Bayesian experimental design: A review," *Statist. Sci.*, vol. 10, no. 3, pp. 273-304, 1995.

[81] International Federation of Digital Seismograph Networks, "CM: Red Sismologica Nacional de Colombia," doi: 10.7914/SN/CM.

[82] Servicio Geologico Colombiano, "RSNC seismic waveform data services," FDSN web services. [Online]. Available: http://sismo.sgc.gov.co:8080/fdsnws/

[83] E. Poveda, G. Monsalve, and C. A. Vargas, "Upper and middle crustal velocity structure of the Colombian Andes from ambient noise tomography," *J. Geophys. Res. Solid Earth*, vol. 123, 2018, doi: 10.1002/2017JB014688.

[84] C. Chiarabba et al., "Subduction system and flat slab beneath the Eastern Cordillera of Colombia," *Geochem. Geophys. Geosyst.*, vol. 17, 2016, doi: 10.1002/2015GC006048.

[85] D. Giardini, P. Lognonne, W. B. Banerdt, et al., "The seismicity of Mars," *Nature Geoscience*, vol. 13, no. 3, pp. 205-212, 2020, doi: 10.1038/s41561-020-0539-8.

[86] Ji Zhang, A. Kato, H. Zhu, and Jie Zhang, "Small earthquake location via machine learning with insufficient data," *Earth, Planets and Space*, vol. 77, no. 1, art. 133, 2025, doi: 10.1186/s40623-025-02258-x.

[87] S. Yang, R. Chen, and P. Deng, "Feature encoded physics-informed neural network modelling with neural tangent kernel analysis for wave propagation," *Eng. Appl. Artif. Intell.*, vol. 163, 113122, 2026, doi: 10.1016/j.engappai.2025.113122.

[88] M. Zhang, W. L. Ellsworth, and G. C. Beroza, "Rapid earthquake association and location," *Seismol. Res. Lett.*, vol. 90, no. 6, pp. 2276-2284, 2019, doi: 10.1785/0220190052.

[89] C. N. Pennington, I. W. McBrearty, Q. Kong, and W. R. Walter, "A comparison of machine learning methods of association tested on dense nodal arrays," *Seismol. Res. Lett.*, vol. 96, no. 4, pp. 2201-2212, 2025, doi: 10.1785/0220240290.

[90] E. Castillo, D. Siervo, and G. A. Prieto, "Colombian seismic monitoring using advanced machine-learning algorithms," *Seismol. Res. Lett.*, vol. 95, no. 5, pp. 2971-2985, 2024, doi: 10.1785/0220240036.

[91] J. Braunmiller, J. Nabelek, and A. Ghods, "Sensor orientation of Iranian broadband seismic stations from P-wave particle motion," *Seismol. Res. Lett.*, vol. 91, no. 3, pp. 1660-1671, 2020, doi: 10.1785/0220200019.

[92] P. Pornsopin, P. Pananont, K. P. Furlong, and E. Sandvol, "Sensor orientation of the TMD seismic network (Thailand) from P-wave particle motions," *Geoscience Letters*, vol. 10, art. 24, 2023, doi: 10.1186/s40562-023-00278-7.

[93] K. Tarumi and K. Yoshizawa, "Station-orientation catalog for Australian broadband seismic stations," *Pure Appl. Geophys.*, vol. 182, pp. 3953-3965, 2025, doi: 10.1007/s00024-025-03827-7.

[94] A. K. Doran and G. Laske, "Ocean-bottom seismometer instrument orientations via automated Rayleigh-wave arrival-angle measurements," *Bull. Seismol. Soc. Am.*, vol. 107, no. 2, pp. 691-708, 2017, doi: 10.1785/0120160165.

[95] S. Park, V. C. Tsai, and M. Ishii, "Frequency-dependent P wave polarization and its subwavelength near-surface depth sensitivity," *Geophys. Res. Lett.*, vol. 46, no. 24, pp. 14377-14384, 2019, doi: 10.1029/2019GL084892.

[96] L. Cristiano, T. Meier, F. Kruger, H. Keers, and C. Weidle, "Teleseismic P-wave polarization analysis at the Grafenberg array," *Geophys. J. Int.*, vol. 207, no. 3, pp. 1456-1471, 2016, doi: 10.1093/gji/ggw339.

[97] S. M. Mousavi and G. C. Beroza, "Bayesian-deep-learning estimation of earthquake location from single-station observations," *IEEE Trans. Geosci. Remote Sens.*, vol. 58, no. 11, pp. 8211-8224, 2020, doi: 10.1109/TGRS.2020.2988770.

[98] J. Dickey, B. Borghetti, and W. Junek, "BazNet: A deep neural network for confident three-component backazimuth prediction," *Pure Appl. Geophys.*, vol. 178, pp. 2459-2473, 2021, doi: 10.1007/s00024-020-02578-x.

[99] J. Muir, B. Fernando, and E. Barrett, "False positives are common in single-station template matching," *Seismica*, vol. 2, no. 2, 2023, doi: 10.26443/seismica.v2i2.385.

[100] S. K. Noel and M. E. West, "Challenges and opportunities of machine learning earthquake detection for regional monitoring," *Bull. Seismol. Soc. Am.*, vol. 115, no. 5, pp. 2096-2108, 2025, doi: 10.1785/0120250021.

[101] G. Tepp et al., "Improvements from incorporating machine learning algorithms into near real-time operational post-processing," *Scientific Reports*, vol. 15, art. 28938, 2025, doi: 10.1038/s41598-025-14491-1.

[102] M. Riedel-Hornig, C. Sippl, A. Tassara, S. Ruiz, B. Potin, and J. Puente, "Seismicity catalogue of the entire Chilean margin (18 to 56 S) from an automated approach," *Earth Syst. Sci. Data*, vol. 18, no. 1, pp. 1-16, 2026, doi: 10.5194/essd-18-1-2026.

[103] J. Munchmeyer et al., "Which picker fits my data? A quantitative evaluation of deep learning based seismic pickers," *J. Geophys. Res. Solid Earth*, vol. 127, e2021JB023499, 2022, doi: 10.1029/2021JB023499.

[104] J. D. Castro, C. Meneses-Ponce, R. Ortega, J. Beltran-Gracia, and H. Gonzalez-Huizar, "Single station seismic event location with a neural-guided mixture model: Model-based weighting for improved accuracy," *J. Geophys. Res. Mach. Learn. Comput.*, vol. 1, no. 4, e2024JH000346, 2024, doi: 10.1029/2024JH000346.

[105] W. Sun, H. Tkalcic, and Q. Tang, "Single-station back-azimuth determination with the receiver function rotation technique validated by the locations of earthquakes, impacts, and explosions," *Seismol. Res. Lett.*, vol. 95, no. 5, pp. 2925-2938, 2024, doi: 10.1785/0220240117.

[106] S. M. Ali, E. Qorbani, R. Le Bras, and G. Rambolamanana, "Interactive analysis of the results of NET-VISA, a Bayesian inference system, in CTBTO's International Data Centre bulletin production," *Acta Geophysica*, vol. 73, pp. 71-81, 2025, doi: 10.1007/s11600-024-01398-0.

[107] D. Perez-Forero, I. Koulakov, C. A. Vargas, T. Gerya, and N. Al Arifi, "Lithospheric delamination as the driving mechanism of intermediate-depth seismicity in the Bucaramanga Nest, Colombia," *Scientific Reports*, vol. 13, art. 23084, 2023, doi: 10.1038/s41598-023-50159-4.

[108] C. Burstedde, L. C. Wilcox, and O. Ghattas, "p4est: Scalable algorithms for parallel adaptive mesh refinement on forests of octrees," *SIAM J. Sci. Comput.*, vol. 33, no. 3, pp. 1103-1133, 2011, doi: 10.1137/100791634.

[109] T. Winder, C. A. Bacon, J. D. Smith, T. S. Hudson, and R. S. White, "QuakeMigrate: a Python package for automatic earthquake detection and location using waveform migration and stacking," *Seismica*, vol. 5, no. 1, 2026, doi: 10.26443/seismica.v5i1.1854.

[110] P. Shi, F. Grigoli, F. Lanza, G. C. Beroza, L. Scarabello, and S. Wiemer, "MALMI: An automated earthquake detection and location workflow based on machine learning and waveform migration," *Seismol. Res. Lett.*, vol. 93, no. 5, pp. 2467-2483, 2022, doi: 10.1785/0220220071.

[111] F. Grigoli, S. Cesca, O. Amoroso, A. Emolo, A. Zollo, and T. Dahm, "Automated seismic event location by waveform coherence analysis," *Geophys. J. Int.*, vol. 196, no. 3, p. 1742, 2014.

[112] J. Vlcek, T. Fischer, and J. Vilhelm, "Back-projection stacking of P- and S-waves to determine location and focal mechanism of microseismic events recorded by a surface array," *Geophys. Prospecting*, vol. 64, pp. 1428-1440, 2016, doi: 10.1111/1365-2478.12349.

[113] K. Bredies and H. K. Pikkarainen, "Inverse problems in spaces of measures," *ESAIM: Control Optim. Calc. Var.*, vol. 19, no. 1, pp. 190-218, 2013, doi: 10.1051/cocv/2011205.

[114] Y. de Castro and F. Gamboa, "Exact reconstruction using Beurling minimal extrapolation," *J. Math. Anal. Appl.*, vol. 395, no. 1, pp. 336-354, 2012, doi: 10.1016/j.jmaa.2012.05.011.

[115] V. Chandrasekaran, B. Recht, P. A. Parrilo, and A. S. Willsky, "The convex geometry of linear inverse problems," *Found. Comput. Math.*, vol. 12, no. 6, pp. 805-849, 2012, doi: 10.1007/s10208-012-9135-7.

[116] Q. Denoyelle, V. Duval, G. Peyre, and E. Soubies, "The sliding Frank-Wolfe algorithm and its application to super-resolution microscopy," *Inverse Problems*, vol. 36, no. 1, art. 014001, 2020, doi: 10.1088/1361-6420/ab2a29.

[117] N. Boyd, G. Schiebinger, and B. Recht, "The alternating descent conditional gradient method for sparse inverse problems," *SIAM J. Optim.*, vol. 27, no. 2, pp. 616-639, 2017, doi: 10.1137/15M1035793.

[118] L. L. Scharf and B. Friedlander, "Matched subspace detectors," *IEEE Trans. Signal Process.*, vol. 42, no. 8, pp. 2146-2157, 1994, doi: 10.1109/78.301849.

[119] K. J. Worsley, "Local maxima and the expected Euler characteristic of excursion sets of chi-square, F and t fields," *Adv. Appl. Probab.*, vol. 26, no. 1, pp. 13-42, 1994, doi: 10.2307/1427576.

[120] J. E. Taylor, "A Gaussian kinematic formula," *Ann. Probab.*, vol. 34, no. 1, pp. 122-158, 2006, doi: 10.1214/009117905000000594.

[121] D. J. A. Chambers, J. P. Dando, G. A. Jones, R. Velasco, and S. A. Wilson, "Moment tensor migration imaging," *Geophys. Prospect.*, vol. 62, no. 4, pp. 879-896, 2014, doi: 10.1111/1365-2478.12108.

[122] D. Anikiev, J. Valenta, F. Stanek, and L. Eisner, "Joint location and source mechanism inversion of microseismic events: benchmarking on seismicity induced by hydraulic fracturing," *Geophys. J. Int.*, vol. 198, no. 1, pp. 249-258, 2014, doi: 10.1093/gji/ggu126.

[123] S. Tsuruoka, H. Kawakatsu, and T. Urabe, "GRiD MT (grid-based real-time determination of moment tensors) monitoring the long-period seismic wavefield," *Phys. Earth Planet. Inter.*, vol. 175, no. 1-2, pp. 8-16, 2009, doi: 10.1016/j.pepi.2008.02.014.

[124] M. Jaggi, "Revisiting Frank-Wolfe: projection-free sparse convex optimization," in *Proc. 30th Int. Conf. Machine Learning*, PMLR vol. 28, no. 1, pp. 427-435, 2013.

[125] Y. Hirano, H. Kawakata, and I. Doi, "A matched-filter technique with an objective threshold," *Sci. Rep.*, vol. 12, art. 22090, 2022, doi: 10.1038/s41598-022-25839-2.

[126] F. Lehmann, F. Gatti, and D. Clouteau, "Multiple-input Fourier neural operator (MIFNO) for source-dependent 3D elastodynamics," *J. Comput. Phys.*, vol. 527, art. 113813, 2025, doi: 10.1016/j.jcp.2025.113813.

[127] J. E. Taylor, A. Takemura, and R. J. Adler, "Validity of the expected Euler characteristic heuristic," *Ann. Probab.*, vol. 33, no. 4, pp. 1362-1396, 2005, doi: 10.1214/009117905000000099.

[128] R. J. Adler and J. E. Taylor, *Random Fields and Geometry* (Springer Monographs in Mathematics). New York: Springer, 2007, doi: 10.1007/978-0-387-48116-6.

[129] I. Vera Rodriguez, M. D. Sacchi, and Y. J. Gu, "Simultaneous recovery of origin time, hypocentre location and seismic moment tensor using sparse representation theory," *Geophys. J. Int.*, vol. 188, no. 3, pp. 1188-1202, 2012, doi: 10.1111/j.1365-246X.2011.05323.x.

[130] I. Vera Rodriguez, M. D. Sacchi, and Y. J. Gu, "A compressive sensing framework for seismic source parameter estimation," *Geophys. J. Int.*, vol. 191, no. 3, pp. 1226-1236, 2012, doi: 10.1111/j.1365-246X.2012.05659.x.

[131] H. Jamali-Rad, Z. Tang, X. H. Campman, A. Droujinine, and G. J. T. Leus, "Sparsity-aware multiple microseismic event localization blind to the source time-function," *Geophys. Prospect.*, vol. 63, no. 1, pp. 70-77, 2015, doi: 10.1111/1365-2478.12166.

[132] H. Kawakatsu, "On the realtime monitoring of the long-period seismic wavefield," *Bull. Earthq. Res. Inst. Univ. Tokyo*, vol. 73, pp. 267-274, 1998.

[133] O. Zhebel and L. Eisner, "Simultaneous microseismic event localization and source mechanism determination," *Geophysics*, vol. 80, no. 1, pp. KS1-KS9, 2015, doi: 10.1190/geo2014-0055.1.

[134] C. Liang, Y. Yu, Y. Yang, L. Kang, Y. Chen, and F. Wu, "Joint inversion of source location and focal mechanism of microseismicity," *Geophysics*, vol. 81, no. 2, pp. KS41-KS49, 2016, doi: 10.1190/geo2015-0272.1.

[135] J. Trojanowski and L. Eisner, "Comparison of migration-based location and detection methods for microseismic events," *Geophys. Prospect.*, vol. 65, no. 1, pp. 47-63, 2017, doi: 10.1111/1365-2478.12366.

[136] L. Li et al., "Recent advances and challenges of waveform-based seismic location methods at multiple scales," *Rev. Geophys.*, vol. 58, e2019RG000667, 2020, doi: 10.1029/2019RG000667.

[137] D. B. Harris, "Subspace detectors: theory," Lawrence Livermore National Laboratory, Tech. Rep. UCRL-TR-222758, 2006.

[138] F. Grigoli et al., "Automated microseismic event location using master-event waveform stacking," *Scientific Reports*, vol. 6, art. 25744, 2016, doi: 10.1038/srep25744.

[139] S. Bates, E. Candes, L. Lei, Y. Romano, and M. Sesia, "Testing for outliers with conformal p-values," *Ann. Statist.*, vol. 51, no. 1, pp. 149-178, 2023, doi: 10.1214/22-AOS2244.

[140] C. Poon and G. Peyre, "Multidimensional sparse super-resolution," *SIAM J. Math. Anal.*, vol. 51, no. 1, pp. 1-44, 2019, doi: 10.1137/17M1147822.

[141] S. Sharan, R. Wang, and F. J. Herrmann, "Fast sparsity-promoting microseismic source estimation," *Geophys. J. Int.*, vol. 216, no. 1, pp. 164-181, 2019, doi: 10.1093/gji/ggy415.

[142] S. Holm, "A simple sequentially rejective multiple test procedure," *Scand. J. Statist.*, vol. 6, no. 2, pp. 65-70, 1979.

[143] L. Zhu and D. V. Helmberger, "Advancement in source estimation techniques using broadband regional seismograms," *Bull. Seismol. Soc. Am.*, vol. 86, no. 5, pp. 1634-1641, 1996, doi: 10.1785/BSSA0860051634.

[144] T.-S. Pham, "Gradient-based joint inversion of point-source moment tensor and station-specific time-shifts," *Geophys. J. Int.*, vol. 238, no. 2, pp. 783-793, 2024, doi: 10.1093/gji/ggae188.

[145] G. Zenhausern, S. C. Stahler, J. F. Clinton, D. Giardini, S. Ceylan, and R. F. Garcia, "Low-frequency marsquakes and where to find them: Back azimuth determination using a polarization analysis approach," *Bull. Seismol. Soc. Am.*, vol. 112, pp. 1787-1805, 2022, doi: 10.1785/0120220019.

[146] H. Morioka, H. Kumagai, and T. Maeda, "Theoretical basis of the amplitude source location method for volcano-seismic signals," *J. Geophys. Res. Solid Earth*, vol. 122, pp. 6538-6551, 2017, doi: 10.1002/2017JB013997.

[147] K. Sawazaki and B. Enescu, "Imaging the high-frequency energy radiation process of a main shock and its early aftershock sequence: The case of the 2008 Iwate-Miyagi Nairiku earthquake, Japan," *J. Geophys. Res. Solid Earth*, vol. 119, pp. 4729-4746, 2014, doi: 10.1002/2013JB010539.

[148] S. Donner, M. Bernauer, and H. Igel, "Inversion for seismic moment tensors combining translational and rotational ground motions," *Geophys. J. Int.*, vol. 207, pp. 562-570, 2016, doi: 10.1093/gji/ggw298.

[149] C. Werner and E. H. Saenger, "Obtaining reliable source locations with time reverse imaging: limits to array design, velocity models and signal-to-noise ratios," *Solid Earth*, vol. 9, pp. 1487-1505, 2018, doi: 10.5194/se-9-1487-2018.

[150] S. Schippkus and C. Hadziioannou, "Matched field processing accounting for complex Earth structure: method and review," *Geophys. J. Int.*, vol. 231, pp. 1268-1282, 2022, doi: 10.1093/gji/ggac240.

[151] Y. Font, H. Kao, S. Lallemand, C.-S. Liu, and L.-Y. Chiao, "Hypocentre determination offshore of eastern Taiwan using the Maximum Intersection method," *Geophys. J. Int.*, vol. 158, no. 2, pp. 655-675, 2004, doi: 10.1111/j.1365-246X.2004.02317.x.

[152] A. Lomax, "A reanalysis of the hypocentral location and related observations for the great 1906 California earthquake," *Bull. Seismol. Soc. Am.*, vol. 95, no. 3, pp. 861-877, 2005, doi: 10.1785/0120040141.

[153] B. Mamandipoor, D. Ramasamy, and U. Madhow, "Newtonized orthogonal matching pursuit: frequency estimation over the continuum," *IEEE Trans. Signal Process.*, vol. 64, no. 19, pp. 5066-5081, 2016, doi: 10.1109/TSP.2016.2580523.

---

# Appendices G to AA. Supplementary Evidence Volume

Appendices G to AA are in the supplementary evidence volume, `Supplementary_Evidence_Appendices_G_to_AA.md`,
which reproduces them in full with the lettering and the numbering of equations, tables and references
of this document; its citations resolve in Appendix F. What each appendix establishes is stated below
in one sentence.

**Appendix G. Pre-Study: Identifiability and Station Count.** The noise-limited resolution constant
$\kappa\approx0.10$-$0.13$ does not depend on $V_P/V_S$ (flat to $0.8\%$, the registered null E2), and
the reduced station-count thresholds of Appendix D follow.

**Appendix H. Pre-Study: Where Physical Evidence Wins.** Physical evidence improves resolution between
nearby sources by only six to fourteen per cent (null E1) but enlarges fourfold the margin by which
clutter-supported hypotheses are rejected (G2).

**Appendix I. Pre-Study: Observability, Coercivity, and Network Geometry.** The leading-order
observability certificate of Section 6.1 predicts the measured coercivity at $\rho=+0.9998$ in
homogeneous and smooth-gradient media, and azimuthal gap misranks linear arrays and rings.

**Appendix J. Pre-Study: Which Constraints Carry the Information.** Polarization carries most of the
discrimination margin, travel time alone has zero margin on a ring of equidistant stations, and the
silence channel is useful only when station distances are spread (G4).

**Appendix K. Pre-Study: Medium-Independent Constraints, and a Negative Result on Overlap.** The
kinematic admissibility test C-III needs only a velocity bound, and overlapping events with complete
picks leave an empty ambiguity set (null E3).

**Appendix L. Pre-Study: Behaviour under Velocity-Model Error.** Under structural model error every
prediction-based channel, polarization included, degrades at slopes within twenty-five per cent of one
another (null E4).

**Appendix M. Pre-Study: The Minimal Sufficient Set.** In eight cells with a perfect picker,
three-component polarization alone matched the full constraint stack, the channel ordering was stable,
and one-sided networks form the failure class of coincidence counting (preliminary A1, A2, C2).

**Appendix N. Pre-Study: Cost, and Network Design.** The robust constraint pair costs a small fraction
of a waveform evaluation (A4), and network designs must be compared on the false-event rate because
location error saturates (B2).

**Appendix O. Pre-Study: The Sufficiency Ladder.** Against the coincidence-count ablation, travel time
suffices at six stations or more, polarization is needed at four, and only an oracle matched filter
succeeds at three (C1, F1).

**Appendix P. Pre-Study: One Station, Several Events.** One three-component station recovers up to two
overlapping events to about a kilometre and fails for three (C3).

**Appendix Q. Pre-Study: The Certificate under Interfaces.** The certificate survives a two-layer
interface only with the head-wave branch summed, at $\rho=+0.963$ and with degraded sharpness (B1).

**Appendix R. Pre-Study: Detection and Association are One Axis.** The value of pick-free waveform
evidence grows with the picker miss rate, with a threshold in $(0.3,0.5]$ (D1).

**Appendix S. Pre-Study: Observability Maps Where Discrimination is Hard.** On one network the local
observability constant predicts the local false-hypothesis margin, Spearman $+0.573$ (B3).

**Appendix T. Pre-Study: Coherent Nearby Clutter.** A coherent nearby decoy raises the spurious rate of
the polarization channel to $0.40$ [$0.23$, $0.59$]; T.1 and T.2 report the first real-data
measurement, in which single-station P polarization passes its reliability gate at one of five RSNC
stations, so that unreliable stations must enter as neutral evidence.

**Appendix U. Pre-Study: Association Without a Velocity Model.** An associator built from measured
back-azimuths and a velocity bound alone shows no systematic degradation under model error, at a cost
of about an order of magnitude in location accuracy (A3).

**Appendix V. Ruled-Out Approaches and Methodological Failure Analyses.** Five approaches do not work
as intended: second-order rank as a source counter (E5), Betti reciprocity from station data,
moment-tensor physicality as a discriminator, an observability-based stopping rule for refinement,
and grid-aligned location error as a measure of method quality.

**Appendix W. Statistical Status of the Pre-Studies.** Intervals and exact tests for the pre-studies
state which of their differences are supported, and show that only the one-sided cells discriminate.

**Appendix X. Published Associators as Baselines.** GaMMA and PyOcto run on the exact-pick generator
under three successive protocols; the sequence, the replication across seeds and the failure of
statement (i) against a fair ablation are documented.

**Appendix Y. Controls for the Catalogue Endpoint: Ties, Timing Error, Empty Windows and Tuned
Baselines.** The programme that leads to G1: tie, timing and empty-window controls (Y.1 to Y.3), the
physically consistent generator and the coherent moment-tensor detector with its geometry-breaking
controls (Y.4 to Y.13), realistic pickers, QuakeMigrate, catalogue station corrections and threshold
calibration (Y.14 to Y.17), the eight preregistered synthetic tasks 2A to 2F, 2I and 2L (Y.18 to Y.26, Y.28), and
two exploratory studies with a literature map on evidence without travel times and on the timing budget
(Y.27), and the exploratory studies that followed the eighth task (Y.29 to Y.29c), and the first registered
real-data test, negative, on the 2019 Ridgecrest sequence (Y.30), the registered test of G1-R on independent
real noise, positive (Y.31), and the exploration of the R0 failure (Y.32). Y.21
(Task 2D) finds that with the same peak search block orthogonal matching pursuit matches MT-FW, so
that the peak search, and no other element of MT-FW, carries its advantage over [129] and [122] in
standard windows, and that split-conformal thresholds frozen on empty windows pass cell by cell without
giving a deployable operating point. Y.22 (Task 2E) refutes G1-M: the off-grid group BLASSO with time
shifts does not beat block pursuit with equal tolerance, and an independent verification finds no
implementation error or budget handicap and a headroom of 0.08 to 0.24 below a realistic oracle ceiling.
Y.23 reports seven exploratory candidates on the calibration seed, mostly null, among them a physical
velocity model and a source-pulse-invariant subspace. Y.24 (Task 2F) confirms on fresh seeds that the
per-station cap raises recall in dense sequences (G1-R), with the narrower scope that its verification
established. Y.25, exploratory on the calibration seed, finds that the amplitude factors of the generator
are redrawn in every window and are not station terms, that catalogue amplitude terms and even the true
factors do not replace the cap, whose mechanism is bounded influence of each station and phase block
under amplitude model error, and that a Huber fit with a threshold set from the noise distribution
matches the cap at the nominal spread, the candidate registered as Task 2I. Y.26 (Task 2I) confirms on
two further fresh seeds that this Huber fit, with its level fixed by rule, exceeds block pursuit in
dense sequences, is not inferior to the calibrated cap and is equivalent to block pursuit on isolated
events when pooled, with the scope that its verification established: a smaller gain under a strict
origin-time hit rule, a better ranking and selection of atoms and not more events found, about nine
times the cost, a level that is immaterial between 2.5 and 3.4 noise standard deviations, and a Huber fit
behind the cap at an amplitude spread of 0.8. Y.27, exploratory on the calibration seed, finds that
polarization, amplitudes and the S to P ratio are limited by information in this benchmark, 4.4 to 5.6
km against 0.4 to 0.7 km for the arrival times, that soft timing priors do not help at catalogue-grade
station corrections, and that with the true station delays the detector reaches the realistic ceiling
for isolated events, so that the headroom divides into station-term calibration and amplitude evidence.
Y.28 reports the exploration that found an estimator of station terms without picks and Task 2L, which
confirms on two further fresh seeds that those terms raise the origin-time-matched recall of the refined
detector (G1-C), with the scope that its verification established: a gain concentrated in the one-sided
eight-station cell, none demonstrated with four stations, terms worse than the registered ones in about
22 per cent of the networks, a result conditional on the refinement and on catalogue errors of known
size, and 0.28 to 0.61 of the gain of the true delays. Y.29, exploratory on the calibration seed, finds that
a double-couple restriction of the statistic lowers the threshold by about 10 per cent, as a tube formula
predicts, without changing recall; that a decoy pass with permuted station moveouts reproduces the
operating point without labels under the location-only rule but not under the origin-time rule; and that
remedies for self-interference after deflation halve the false atoms that the events induce without moving
recall, the largest class of false events at the operating point under the origin-time rule being
well-placed detections of real events with a biased origin time.

**Appendix Z. Certified Pruning of the Hypothesis Space, and the Role of p4est.** Branch-and-bound
returns catalogues identical to the dense search, two to three times cheaper for the kinematic
likelihood and exact but almost without pruning for the moment-tensor certificate, and p4est is
reserved for distributing that search.

**Appendix AA. Theory of the Group BLASSO Detector.** Existence and optimality conditions, the identity
between the certificate at zero and the generalized likelihood ratio, validity of the pruning bound,
the Frank-Wolfe rate with an exact oracle, finite-sample validity of the conformal threshold and a
curvature bound by observability are proved; recovery over the continuum and the behaviour of the
deployed greedy detector are stated as conjecture or partial results.
