# Executive summary

**Epicentre determination and seismic source attribution as an inverse problem with
constraints generated from elastodynamics**

Juan José Gordillo Córdoba · October 2026

Decision document. It contains the formulation, what has been measured, what was discarded and why,
and the concrete pending actions. The detail is in the full proposal (appendices A to
AA), in the audit (`audit_status.md`) and in the implementation error log.

**Status as of 16 October 2026: G1-R confirmed with independent real noise; the R0 failure narrowed to the
reading of arrivals.** Task 2M, registered before any computation on its seeds, injected synthetic events into real
noise of days of the Alto Tiberina network that no earlier analysis used. On the primary day, whose noise has
transients, the per-station cap and the Huber fit beat BOMP-PK by +0.228 [+0.148, +0.298] and +0.241 [+0.179,
+0.293] in dense windows (both rules met, Holm rejects both); on two near-Gaussian days, a secondary replicate, by
+0.057 [+0.004, +0.121] and +0.120 [+0.070, +0.168]. A breakdown written after the result shows that the gain of the
cap is concentrated in the hours with transients (+0.033 at night, +0.300 between 7.5 and 15 h), while the Huber fit
keeps +0.10 to +0.15 in quiet noise. An independent verification confirmed the positive with caveats of scope; two
networks of the secondary seed failed by exhausting their catalogue reserve, as the registration allows (E48).
Exploration C18, on the development segment of the Ridgecrest hour only, tested the two remedies the R0 diagnosis
named. A layered or fitted medium does not make the energy-peak readings more predictable (0.53 to 0.63 of held-out
readings within 0.25 s for every model), so the timing error is a reading problem, worse for S, not a velocity-model
problem. The cap does erase the evidence of real events even at their true hypocentres (AUC against event-free
windows 0.39 to 0.49), and a logarithmic compression restores it (0.66 to 0.72) while keeping the robustness of the
cap to transients on real noise; but no variant detects the events, because the score maxima are not where the
events are. A gate is now registered: a new real test, on the Hi-net swarm kept unread, is registered only when a
candidate recalls at least 5 of the 16 development events with origin time at 30 false detections per hour. The
closure of R0 (M1 and M3 on the four subnetworks, declared post hoc) is in progress. Appendices Y.30 to Y.32;
errors E48 and E49.

**Status as of 10 October 2026: first test with real data, negative.** Task R0, registered
before touching the data, tested the thesis detector (MT-FW with cap 6, pick-free S2 terms and
T3 refinement) on one real hour of the 2019 Ridgecrest aftershock sequence (CI network, analyst
catalogue, four subnetworks of 4 to 8 stations). At 30 false detections per hour and with the origin time
required, it recalled at most 0.022 of the 45 events of M ≥ 2.0, against up to 0.667 for PhaseNet with PyOcto;
the pooled difference with the best comparator was −0.383 [−0.481, −0.310] and the registered rule is not
met. A reading tolerance widened by rule to the measured timing error changed nothing. The measured
causes: with a homogeneous medium and constant station terms the real timing error is 0.35 to 0.41 s
and grows with the separation between sources (the synthetic generator gave about 0.14 s), so that the detector
places the events at 8-18 km in median; and in a continuous sequence the capped score saturates and does not
separate the real events. By contrast, an exploration with real noise and injected events (C16) reinforced G1-R,
because the real transients, partly telemetry defects, sink BOMP-PK. The synthetic advantage is
conditional on a timing model with errors smaller than the tolerance; the work that follows is that model
(source-dependent corrections) and a ranking that withstands continuous seismicity, before G1'.
Appendices Y.29b (C15), Y.29c (C16) and Y.30 (R0); errors E46 and E47.

**Status as of 3 October 2026.** The proposal is titled *Timing-Tolerant, Bounded-Influence
Moment-Tensor Detection of Overlapping Earthquakes: A Preregistered Benchmark and a Real-Data Test Plan*.
Seven preregistered synthetic tasks support three things: the advantage of the family of
timing-error-tolerant moment-tensor detectors over published associators and migration-based
locators, with the operating point chosen on the evaluation and without meeting G1 at 0.1 against GaMMA
with optimistic picks (Tasks 2A to 2D); the refutation of G1-M, the off-the-grid group BLASSO
(Task 2E); and the confirmation of G1-R, a gain of 0.04 to 0.07 in recall in dense sequences from
bounding the influence of each station, which drops to 0.02 to 0.045 with the origin time required (Tasks 2F and
2I). The secondary arms of 2I, descriptive, show that the Huber level does not matter between 2.5 and 3.4
noise deviations and that with amplitude dispersion 0.8 the Huber fit falls below the cap.
Two explorations on the calibration seed (C9 and C10) place the remaining margin in the calibration of the
station terms and in the amplitude evidence, not in physics without travel times. The fourth independent audit gives
57.8/100 on the fixed rubric (the previous ones, 56, 51 and 51) and estimates that a doctoral committee would give around
70. The target of 75/100 is not reached: all the evidence is synthetic and from a single generator, the
preregistration is local and not sealed by a third party, and still missing are the first real result of G1' on the
decimated RSNC, a second independent generator, a deployable operating point, theory on the
deployed statistic, the trimming of the scope and of the length, and a sealed public registration.

**Task 2L and hypothesis G1-C (5 October 2026).** An eighth preregistered task tested the
calibration part. Exploration C11 (seed 999, with the plan frozen with a hash before any detection
evaluation) found a pick-free station-term estimator: a joint determination of the
hypocentres of the 20 catalogue events and of the terms from the energy peaks of the matched
filter (S2). Task 2L (seeds 7171 and 9393, never used; registration hash 13a6e9b1...a62a delivered
before running; code frozen) compares the same detector, MT-FW with cap 6 and the T3 refinement, with the
S2 terms and with the registered corrections, and takes as primary the hit rule that requires the origin time
to 0.5 s. The four rules are met on both seeds: +0.057 [+0.034, +0.081] and +0.062 [+0.036,
+0.083] in standard windows (2L-a); +0.038 [+0.024, +0.054] and +0.030 [+0.015, +0.044] in dense ones (2L-b); without
loss by regime (2L-c, lower bounds between +0.007 and +0.052) nor with the epicentre-only rule (2L-d,
+0.036 and +0.050 in standard, +0.020 and +0.035 in dense); the registered secondary on BOMP-PK gives +0.068 and
+0.071. The independent verification confirms the positive, without leakage (S2 recomputed without access to the
true delays or catalogue is identical in 300 of 300 units), and bounds the scope: the pooled
gain is in essence that of the one-sided 8-station cell (+0.108 and +0.113); with 4 stations there is no
demonstrated gain; S2 worsens the terms in about 22 % of the networks; the effect is conditional on the
T3 refinement; it assumes a catalogue of 20 events with errors of known size (1 km and 0.1 s); and it captures
between 0.28 and 0.61 of the gain of the true-delay oracle. The evidence is synthetic and from a single
generator family, and joint determination with station terms is classical in pick-based
location, so G1-C is a measured calibration result, not an algorithm. Since 2L the hit rule
with origin time is the primary one. The preregistered tasks become eight; the fourth audit
(57.8/100) predates 2L and has not been rescored. Detail in Appendix Y.28.

**Explorations C12, C13 and C14 (5 to 9 October 2026; seed 999, exploratory; Appendix Y.29).** Three
explorations after 2L, with no registered hypothesis and with the registered detector of the previous
tasks (MT-FW with cap 6, S0 corrections and hard window; not the refined detector of 2L). C12 restricted
the subspace statistic to double couples: the tube formula predicts its null law, with 4.4 to 4.5 effective degrees
of freedom instead of 6, simulation confirms it and the threshold drops by 10 % in energy at equal
false-alarm rate, but the recall does not change (+0.005 [−0.021, +0.033] in standard windows with origin
time), because the hits lose as much energy to the restriction as the threshold drops; it remains as an
analytical piece for a mechanism test, not as a detector. C13 calibrated the threshold with decoys: an extra
pass of the detector with the station moveout permuted reproduces without labels the operating point
of the benchmark (recall within ±0.03 of the oracle) and keeps the FDR within q = 0.2 with the epicentre-only
rule, but not with the origin-time rule nor, at q ≤ 0.1, in the one-sided geometry; it is a
partial result and does not yet give a deployable operating point. C14 tested remedies for
self-interference after deflation: they halve the false detections induced by the events, but they do not move
the origin-time recall at the operating point, and the most aggressive ones erase a quarter of the
overlapping events found at any score; it is negative. Its diagnosis is the useful finding: at the
operating point the largest class of false detections under the origin-time rule are detections of real
events that are well located (median of 3.2 and 3.0 km) with the origin shifted by a median of 0.8 s, early in
73-75 % of the cases and with no other correct detection of the same event in 86-92 %: 55-58 % of the false detections
in standard windows and 80-95 % in dense ones on the development networks, and one third and two thirds on the
evaluation ones. No decoy imitates them and no deflation remedy reduces them, so the next
target of the calibration is the origin time of the well-located detections. None of this changes a
registered rule or the score of the fourth audit.

---

## 1. The problem

Given a window of continuous seismic recording, determine how many events occurred, where and
when. The difficulty is not locating an isolated event —that is solved— but **deciding which
sets of arrivals belong to the same real event and which are spurious coincidences**, in
the presence of false detections, overlapping events and networks with poor coverage.

The associators in common use with machine-learning picks (GaMMA, PyOcto, REAL, GENIE,
HARPA) operate on **arrival times**: they search for sets of picks compatible with some
hypocentre. NET-VISA, in treaty monitoring, already adds azimuth, slowness and a model of false
arrivals per station. The question of the thesis is how much additional physics is needed, and where, to
reject false events.

## 2. What is proposed

Pose epicentre determination as **constrained optimisation where the feasible
set is generated from the elastodynamic equation**, not from statistical assumptions.

$$\min_{\mathcal{E}}\ \mathcal{J}(\mathcal{E};\mathbf{d})
\quad\text{subject to}\quad
\underbrace{\gamma(\mathbf{x})>0}_{\text{observability}},\
\underbrace{\mathbf{x}\notin\mathcal{X}_c}_{\text{silence}},\
\underbrace{g_{\mathrm{kin}}\le0}_{\text{admissibility}},\
\underbrace{d_{\mathrm{Gr}}\le\epsilon}_{\text{polarisation}},\
\underbrace{s(t)\ge0}_{\text{source}}$$

The primary observation is the **three-component wavefield**, not a set of times. The thesis
in one sentence: *each piece of physics is an additional constraint on the feasible set, and the
clutter is what remains outside*.

**Operational form (September 2026; `method_design.md`, Section 7.5 of the proposal).**
Detecting, counting and locating are posed as a single convex problem: recover a sparse measure
whose amplitudes are moment tensors, with the group total-variation norm as a relaxation of the
number of events:

$$\hat\mu\in\arg\min_{\mu}\ \tfrac12\|\mathbf d-\Phi\mu\|^2_{\Sigma^{-1}}+\lambda|\mu|_{\mathrm{TV},G}.$$

Its optimality certificate is a coherent detector: the energy of the whole network's recording
projected onto the 6-dimensional subspace of the radiation pattern. The pilot (MT-FW) uses it
greedily: it explains an event, subtracts it and searches for the next one, with a no-reselection radius of 5 km.
It is a group matching pursuit, not Frank-Wolfe, and it does not inherit its convergence rate
(Appendix AA, E29). Its certified oracle gives the same records as the exhaustive search, and the
threshold is calibrated on noise-only windows of the same regime, conformally or by extreme
values: random field theory does not fix the false-alarm rate of this statistic (E28).
Physics enters per station according to a
fidelity ladder (times, ray geometry, radiation pattern, 3D Green's functions), and machine
learning accelerates or robustifies (neural Green's functions, learned noise model, amortised
proposals) without being able to make the search omit an event.

## 3. The finding that reoriented the project, and its limit

A pre-study compared full physical prediction against arrival-time prediction, with
geometry and energy matched.

| Metric | Gain from physics |
|:--|:--|
| Resolution between two close true sources | **6 to 14%**, marginal |
| Discrimination margin against false hypotheses | **factor 4**, degrading 3 times more slowly |

**Physics does not buy resolution; it buys discrimination margin.** That changed the goal from
"locating better" to "rejecting better".

**When the margin becomes a catalogue (Appendix Y, September 2026).** Against the coincidence-count
ablation the one-sided advantage was a ties effect. In a generator that detects
all arrivals with exact time, a continuous kinematic score equals polarisation. But
that generator gives the kinematic channel information that no picker has. In a consistent
generator, where detection, timing and polarisation degrade together with the signal-to-noise ratio,
more than half of the arrivals are never detected, the sources have a double-couple mechanism, the
clutter follows the same detection rule as the signal and there are delays, deviations and site
factors that the model does not know, the coherent moment-tensor detector with deflation (MT-FW, a
greedy iteration that at each step picks the maximum of the certificate of the convex problem of the
design document) beats tuned GaMMA by
+0.11 to +0.28 of recall at 0.25 false events per window with weak clutter and by +0.17 to +0.43 with strong
clutter, and its control with the geometry of another hypothesis by +0.13 to +0.41, on two seeds and on
new windows. Four earlier positives did not survive the audit: a ties effect, a
reading at the nominal budget, a per-pick likelihood that only won with pure SV synthetic S
waves, and margins over GaMMA inflated by clutter picks that no picker would produce. The method
fails if the station delays exceed its search window: the station corrections are
part of it, and, estimated with 20 catalogued events, they recover the performance of the exact correction in
2 of 3 cells. Against QuakeMigrate, in an exploratory single-seed test, MT-FW beats its
fair three-component configuration by +0.28 [+0.10, +0.43] of recall at 0.25 false events per window and
by +0.19 [−0.00, +0.37] at 0.1, above all in windows with two sources. The synthetic picks that
GaMMA received were more optimistic than those of a real picker in recall and in timing, although with strong
clutter they delivered more clutter picks to it (E27).

**The confirmatory test (Task 2A, preregistered, executed and independently verified).** With
a generator that the detectors do not know (gradient medium with curved rays, wavelet of variable frequency and
phase, coda, fixed 0.3 s delays per network), picks from an operational STA/LTA and station
corrections estimated from 20 catalogued events with location errors, MT-FW with corrections beats the
best of GaMMA, PyOcto and QuakeMigrate in the 3 cells of both seeds and both clutter regimes,
by +0.09 to +0.49 of recall at 0.1 false events per window and by +0.17 to +0.46 at 0.25; the same iteration with
the geometry of another node loses the advantage. Three caveats bound the reading. The effective rival is
QuakeMigrate, the best comparator in the 24 cases, because the STA/LTA leaves GaMMA and PyOcto without enough
picks. With optimistic picks GaMMA ties in the one-sided cell with weak clutter and wins in the
single-source windows in 5 of 6 weak-clutter cells. And the per-window intervals ignore the
clustering by network: with a per-network bootstrap, four cells include zero at 0.1. 2A does not isolate the radiation
pattern; Task 2B was preregistered for that.

**The operator ablation (Task 2B, preregistered, executed).** On the same windows, the same
solver with the radiation operator replaced by the identity (sign-free energy, MT-FW-FREE+SC)
loses against MT-FW+SC by +0.06 to +0.35 in the 3 cells of both seeds with strong clutter, where
without the operator the detector only equals QuakeMigrate. With weak clutter the rule is met at the
minimum, 2 of 3 cells per seed (in the one-sided 8-station cell both methods tie), and re-applied in the
bootstrap only in ~60 % of the replicates, against 97-98 % with strong clutter; moreover, without the operator the
detector still beats QuakeMigrate at 0.25 in 5 of 6 cell-seeds, so a good part of the advantage in
weak clutter comes from pick-free migration with deflation (exploratory). The identity removes at once the sign,
the polarisation constraint and the amplitude geometry: 2B shows that the operator matters, not that
the sign matters.

**The closest comparators and the convex method (Task 2C, preregistered, executed).** On the
same windows of 2A, MT-FW+SC beats in-house implementations of the block-sparse inversion of
Vera Rodriguez et al. (2012) and of the mechanism-corrected stack of Anikiev et al. (2014), on the same
grid and with the exact predicted times, by +0.23 to +0.52 of recall at 0.25 false events per window with
weak clutter and by +0.08 to +0.45 with strong clutter; the better of the two is the block inversion in
all cells. With 30 networks per cell and hierarchical bootstrap, the advantage over QuakeMigrate replicates
with the network as unit: +0.22 to +0.35 with weak clutter and +0.099 to +0.295 with strong, with lower bounds
above 0 in all cells. By contrast, the convex problem of the thesis restricted to the grid (group BLASSO
solved with duality gap ≤ 1.4e-8) falls below MT-FW in 23 of 24 cases, by −0.064 to
−0.605. In exploration, in the ideal generator with the sources on the nodes, the block inversion
beats MT-FW (0.97 against 0.76) and the BLASSO reaches 0.73: the exact-time methods collapse because of the coarse grid (3 km and
0.1 s) and the residual timing errors, which the ±0.25 s peak search of MT-FW absorbs. The
sustainable claim is robustness to grid and timing errors on a coarse grid, not superiority
over sparse moment-tensor inversion as such; the convex problem is the framework and the
certificate, not the best detector (E33), and off-the-grid sliding Frank-Wolfe becomes the
key question of the mathematical line.

**Equal tolerance, peak-search ablation and frozen threshold (Task 2D, preregistered with clock
time and hash, executed).** The central result is negative for what is distinctive about the method. With the same
±0.25 s peak search, the block inversion of Vera Rodriguez et al. (BOMP-PK+SC) equals MT-FW+SC:
differences at 0.25 from −0.155 to +0.091, all hierarchical CIs contain 0, and rule 2D-a fails in
both regimes. As preregistered, the advantage of MT-FW over Vera Rodriguez et al. and Anikiev et al.
is attributed only to the peak search; the reading of 2C, which attributed it to MT-FW's own elements,
is refuted (E35). Removing the peak search from MT-FW costs 0.06 to 0.41 (12 of 12 positive; 2D-b is
met); removing the cap changes nothing significant; the grid BLASSO with the same peak search
(GL-FW-PK) falls 0.11 to 0.30 below MT-FW in 11 of 12 cases; and the moment-tensor detectors with
temporal tolerance and deflation, as a family, beat QuakeMigrate by +0.17 to +0.49. The peak
search has a price: with exact times and sources on the nodes it costs 0.10 to 0.14 at 6 and 8 stations
(exploratory). The split-conformal threshold frozen on empty windows of other networks meets its per-cell
rule (2D-c), but it is not a deployable operating point: in windows with events MT-FW+SC gives 1.2 to 2.2
false events per window with weak clutter, 3 to 8 times the false events of QuakeMigrate over all windows. With 5
networks per cell no 2D rule survives Holm. Consequence: the empirical advantage over the
pick-based associators and QuakeMigrate in these synthetics belongs to the family of timing-error-tolerant
moment-tensor detectors, not to MT-FW; what the thesis can contribute as a method is the
principled treatment of that tolerance, the group BLASSO with a time shift per
station and phase as a nuisance parameter, solved off the grid by sliding Frank-Wolfe, which only
counts if it beats the heuristic peak search and BOMP-PK.

**The off-the-grid BLASSO (Task 2E, preregistered; hash b40f0c6e…e988 delivered via the chat before
writing the code; independently verified).** On the 30 networks per cell of 2C, the group
BLASSO with a shift per station and phase, solved off the grid by sliding Frank-Wolfe
(TS-SFW+SC, f = 0.01 calibrated on seed 999), gains +0.000 to +0.039 at 0.25 over BOMP-PK+SC; the 12
points are ≥ 0 but none reaches 0.05 and all hierarchical CIs contain 0. 2E-a and 2E-b fail in
both regimes: **G1-M is refuted in this synthetic design**. The group penalty is inert with the
calibrated f. The independent verification rules out an implementation error (the continuous atom reproduces
the BOMP-PK block to machine precision) and a budget disadvantage (10 times the iterations with 9
starts change the recall by −0.014), and shows that the benchmark is not saturated: the realistic oracle
ceilings lie 0.08 to 0.24 above BOMP-PK, and the missing events are absorbed by atoms misplaced
at 5-10 km, beyond the reach of the local sliding.

**Where the margin is not (exploration `explore_contrib`, seed 999 only, exploratory).** Seven candidates:
coarse grids, large timing errors, a 1D velocity model (oracle or calibrated by JHD) inside
the detectors and a subspace invariant to the source pulse gain nothing or lose; with the generator's
medium equal to the detector's homogeneous one the margin persists (+0.194 with weak clutter). The only candidate
that meets the exploratory criterion is the per-station cap of MT-FW in windows with 3 to 5 overlapping
events (+0.070 at 0.25, replicated on new networks with +0.068). It was chosen among several, which is declared.

**The per-station cap in dense sequences (Task 2F, preregistered; hash 3d8e4c9b…2f5f delivered before
running; seeds 4242 and 8888 never used; code frozen with hashes; verified).** With 96 networks per
cell, MT-FW+SC beats BOMP-PK+SC in dense windows by +0.072 and +0.074 (weak clutter) and +0.051 and +0.043
(strong) at 0.25, all with lower bound > 0 and rejected by Holm (2F-a is met); without the cap MT-FW equals
BOMP-PK with a difference of −0.001 and −0.002, so the gain is due to the cap (2F-b is met); in standard
windows there is equivalence with margin 0.05 (2F-c is met, narrowly on 8888). **G1-R is confirmed in
synthetics.** The adversarial verification finds no leakage or asymmetry and narrows the scope: with origin time
within 0.5 s the gain is +0.024 to +0.045; it depends on the cap level fixed at 1x (peak at 6,
disappears above 12, negative below 3); it is due to the cap bounding the influence of each
station-and-phase block whose amplitude the model predicts poorly (without the generator's amplitude factors it drops
to +0.023 with a CI that includes 0) and to a better calibration of the score across windows, not to better separation of
the events. The pilot's reading, according to which the cap is what allows "resolving" the dense sequences,
exaggerated the mechanism, and the verification bounds it (E37). Absolute recall with cap: only 0.22-0.29.

**Site terms versus bounded influence (exploration C8, seed 999, exploratory; Appendix
Y.25).** The generator draws the amplitude factors in each window and in each catalogue event: they are
per-window amplitude perturbations, not station terms, and having read them as per-station site
amplification was an attribution error (E39). With an added world of persistent per-network factors, the
terms estimated from the catalogue do not replace the cap (−0.066 [−0.108, −0.028] without persistence and −0.040
[−0.078, +0.007] with it), and neither do the true factors; with the cap in place, knowing them adds +0.045 and
+0.032. With dispersion 0.25 the oracle gains nothing and the cap still gains +0.043 (CI including 0): the cap is not
only a remedy for the amplitude factors. A Huber fit with the threshold fixed by rule from
the noise (HUBER-N) equals cap 6 with dispersion up to 0.5, is worse with 0.8 (−0.043 [−0.113, −0.001]) and
costs 9 times more; it was the best of 11 candidates, with the selection bias that this implies. Nothing captures
more than a quarter of the margin up to the realistic ceiling. Its confirmation was registered as Task 2I
(seeds 5151 and 6262).

**A robust certificate with level fixed by the noise (Task 2I, preregistered, seeds 5151 and 6262;
Appendix Y.26).** Registration hash a983beac…6f97 delivered to the user before any run; code
frozen by hashes; 120 networks per cell and regime, 5 dense and 5 standard windows per network. The method applies
the Huber function to the norm of the residual of each 3-component vector per station and phase, with threshold c
= square root of the 0.95 quantile of chi-squared with 3 degrees of freedom = 2.7955 noise deviations, fixed by
rule, 4 IRLS passes and robust deflation. The run was interrupted twice and resumed from the parts
with the code unchanged. **The three rules are met on both seeds:** HUBER-N+SC beats BOMP-PK+SC in
dense windows by +0.058 [+0.039, +0.076] and +0.053 [+0.037, +0.070] (2I-a); it is not inferior to MT-FW+SC with
cap 6, −0.008 [−0.023, +0.006] and −0.002 [−0.018, +0.012], with lower bounds above −0.03 (2I-b); and in
standard windows it is equivalent to BOMP-PK+SC, +0.010 (90 CI [−0.010, +0.027]) and +0.007 (90 CI [−0.013, +0.024])
(2I-c). The independent verification confirms the positive (provenance, symmetry, rederived statistics and
six regenerated parts identical across the two interruptions) and fixes its scope: with the hit
rule that requires origin time within 0.5 s the gain drops to +0.031 and +0.019; equivalence is
met pooled but not with weak clutter on 6262 (−0.032 [−0.064, −0.003]); the recall with all atoms
accepted does not rise, so the gain is one of ranking and atom choice, not of more events
found; at cap 0.1 it falls somewhat below the cap on 5151 (−0.018); it costs about 9 times MT-FW; and the
absolute recall in dense windows remains low (0.235 and 0.251 at cap 0.25). The secondary arms, descriptive, with
20 networks per cell and regime and without independent verification, add two things. The result is robust to the
level: thresholds 2.50 and 3.37 are indistinguishable from the registered one (differences of +0.001 to +0.006 in
dense windows) and with 5 the gain is smaller. It is not robust to a large amplitude dispersion: with
site deviation 0.25 the cap contributes nothing (+0.005 and +0.009) and HUBER-N contributes little (+0.036 and +0.002); with
0.8 HUBER-N beats BOMP-PK (+0.100 and +0.072) but falls below the cap (−0.055 [−0.113, +0.002] and
−0.069 [−0.121, −0.019]), so non-inferiority against the cap does not extend to that dispersion.

**Evidence without travel times and timing budget (explorations C9 and C10, seed 999; Appendix
Y.27).** Exploratory, with oracles and with analyses after the declared results; none of this is
confirmatory. C9 measured how much location information there is in polarisation, amplitudes and the S/P
ratio. By the Cramér-Rao bound, arrival times locate the epicentre to 0.4-0.7 km and all the channels
without times together to 4.4-5.6 km; a joint estimator is bounded at an improvement of 4-5 %; the limit is the
waveform noise (the P arrival measures 1.5-1.8 noise deviations, 32-39 degrees of angular noise),
and that physics would weigh only with about 10 times the signal-to-noise ratio and 2 degrees or less of polarisation
dispersion. Relocating the atoms with that evidence loses recall even with oracles (−0.11 to −0.18 in
dense, −0.23 to −0.48 in standard). Weighting by polarisation admissibility gives a small and
conditional gain (SELPOL, the best of some 50 candidates and chosen after seeing results: +0.046 [+0.007,
+0.070] in dense) that disappears with 20 to 45 degrees of dispersion. C10 tested Gaussian timing priors
inside the detectors: with catalogue-quality station corrections they do not help. The reason is the
timing budget of a true hypothesis, 0.2-0.4 s per arrival: model error 0.05-0.07 s,
error of the corrections 0.12-0.24 s rms and position within the 3 km cell 0.18 s (P) and 0.30 s (S). With
the true statics the gain exists: alone, +0.047 in dense and +0.098 in standard; with the prior at the
level of the model error, +0.081 and +0.146, which in isolated events reaches the realistic ceiling (0.592 against
0.603). With perfect timing, 37-39 % of the isolated events still lose against a false node in
the first iteration, because of amplitude and not timing. A local refinement by the timing likelihood
(T3) improves location and origin time of the isolated events and only adds recall with the hit rule
that requires the origin time, in standard windows (+0.100 and +0.066). The margin then splits into a
calibration part (accuracy of the station terms) and an amplitude part (weak arrivals,
clutter, coda, amplitude factors); the open problem is the self-consistent calibration of the station
terms and the amplitude evidence, not another location principle.

**Consequence.** The methodological contribution, stated with precision and modesty: G1-R is supported by two
preregistered synthetic tests, with an external timestamp (hash delivered outside the container before
running) and independent verification, 2F with a calibrated cap and 2I with a Huber level fixed by rule from
the noise, which removes the dependence on a tuned constant. That is: a bounded-influence moment-tensor
subspace certificate improves recall at a fixed rate of false events in dense sequences
of overlapping events relative to the block search with equal temporal tolerance, without harming the
isolated events when the regimes are pooled; plus the benchmark and its negative results (G1-M included).
Its limits: it is synthetic and from a single generator family; the size is modest (+0.04 to +0.07 with the
registered rule, +0.02 to +0.045 with origin time); the Huber fit costs 9 times more; and clipping
per-station contributions and fitting with Huber are standard robust statistics, so the novelty is
the measured design result, not an algorithm. The mechanism is per-block bounded influence against the model's
amplitude error, not correction of a site term; the certificate with per-station site terms
has no support in C8 and is not pursued; where the model's amplitude error is large
(site deviation 0.8) the calibrated cap is the better of the two forms. The margin of 0.08 to 0.24 is
explained neither by the velocity model nor by the pulse, and about three quarters of it remain unexplained
by amplitude treatments: it is the open problem, which explorations C9 and C10 place in the
calibration of the station terms and in the amplitude evidence.

**What already existed.** A verified novelty audit found that sparse joint detection and location
with moment tensors is not new: Vera Rodriguez, Sacchi and Gu (2012, *GJI* 188(3))
recover origin time, hypocentre and moment tensor of simultaneous sources with a block-sparse search
(block OMP), and GRiD MT (Tsuruoka et al., 2009) has operated since 2003. Mechanism-corrected stacks,
master-event corrections and extreme-value thresholds also already existed. The
proposal exaggerated that novelty (E30). What remains as a contribution is the off-the-grid formulation as a
group BLASSO with tensor amplitudes, the identity between its certificate and the generalised likelihood
ratio of the matched subspace detector (a modest unifying result), the certified search, the conformal thresholds and the benchmark
against GaMMA, PyOcto and QuakeMigrate with overlapping events. After Task 2C that formulation
counts as the framework and certificate of the detector, not as the detector that performs: on the grid its convex
solution falls below MT-FW, and its off-the-grid version with per-station and per-phase shifts
does not beat BOMP-PK (Task 2E). After Task 2D MT-FW is not distinctive either relative to Vera Rodriguez
et al. with the same tolerance in standard windows: neither the greedy search, nor the deflation, nor the cap, nor the
peak search is claimed as new. What remains as a methodological contribution is the measured and
explained effect of the per-station cap in dense sequences (Task 2F), a known practice in robust stacking,
and that of its bounded-influence form with level fixed by rule from the noise (Task 2I), which as a
technique is not new either: Huber M-estimators are standard in moment-tensor inversion and in
robust stacking, and what is contributed is the measured design result.

## 4. Results with statistical support

Five differences supported with Clopper-Pearson intervals and Fisher's exact test, thirty
replicates where there was doubt.

**A. Observability certificate, the most solid result.**

Proposition to leading order (the remainder is stated and its bound is registered work) that
translates the Bardos-Lebeau-Rauch geometric control condition into a coercivity constant
computable **in closed form, from geometry alone, without data**. Verified
in three classes of medium:

| Medium | Correlation with the measured coercivity |
|:--|:--|
| Homogeneous | $\rho = 0{,}9998$ |
| Smooth gradient (with take-off slowness) | $\rho \approx 0{,}999$ |
| Two layers with head wave (sum over branches) | $\rho = 0{,}963$ |

Four corollaries. The most useful: **travel-time observability degenerates over planar
arrays and over equidistant rings**, and two branches break that degeneracy algebraically.
In addition, it beats the azimuthal gap —the standard criterion of the field— by **40 to 90%** in any
medium with a gradient, and it produces a priori a map of where it will be easy to form false events:
correlation $0{,}573$ with the local margin, factor $4{,}6$ between quartiles.

**B. Sufficiency ladder — the result that reframes the problem.**

The amount of physics needed is a function of the **effective number of usable observations**.

| Stations | What is needed (measured against coincidence counting) | Correction from Appendix Y |
|:--|:--|:--|
| 1 | waveform and two branches, with fixed depth: up to **2 overlapping events**, ~1 km | no change |
| 3 | waveform (oracle matched filter); the stack fails 3 of 4 | all methods emit false events in 60-100% |
| 4-5 | polarisation necessary; coincidence counting fails | a continuous kinematic score suffices with exact picks |
| 6 or more | travel times suffice | no change |

The ladder remains the organising idea, but its axis is the **effective kinematic
information**: number of stations, picker loss and temporal precision of the picks.

And detection and association are **a single axis**: six stations with 80% picker loss are equivalent
to one station. The threshold where pick-free evidence becomes dominant lies in $(0{,}3,\ 0{,}5]$.

**C. Failure mode of travel times in the presence of a coherent decoy.** With a coherent neighbouring
event, the counterfactual produces a merged solution —neither the true source nor the decoy—
in **30 of 30 realisations**, $p<10^{-4}$. It is the configuration that dominates the first hours of an
aftershock sequence.

**D. Polarisation carries the margin.** Of the five families, C-IV (polarisation) dominates the
discrimination margin; C-III (kinematic admissibility) is the only one independent of the velocity
model, with **zero rejections** of uncontaminated truths. At catalogue level, the contribution of
polarisation over the best kinematic score is not established (Appendix Y).

## 5. What was discarded, and why

It is documented because it closes entire families that will seem obvious to a new reader.

| Discarded | Reason | Where |
|:--|:--|:--|
| **Rank of a second-order matrix to count sources** (MUSIC, MDL, AIC, signal subspace) | Rank $K$ requires a steering vector that is stable across windows; a transient train with moveout does not have one. It fails both across components and across stations | Apx V.1, null E5 |
| **Betti reciprocity as a test on the data** | The identity lives at the source positions, which are not observed. Analytical blockage, not empirical | Apx V.2 |
| **Physicality of the moment tensor as a discriminant** | Fitting clutter with a six-element basis gives nearly deviatoric tensors by chance: 142 of 491 false ones beat the true one | Apx V.3 |
| **p4est / adaptive octree with observability stopping** | The rule never fires: the resolution limit is finer than the cell reached. Dimensional failure, not a code failure. The score-driven octree does work, but NonLinLoc and Qseek already do that | Apx V.4 |
| **Resolution gain ≥15%** | Measured 6–14%, flat in number of stations | Null E1 |
| **$V_P/V_S$ contrast law** | Flat to 0,8% over $[1{,}50,\ 1{,}95]$ | Null E2 |
| **Overlap as an association difficulty on its own** | Empty ambiguity set with complete picks. The difficulty is one of **detection** | Null E3 |
| **Graceful degradation in location** | No channel preserves precision under model error better than another | Null E4 |
| **BLASSO-type super-resolution** | The seismic operator is not translation invariant | §9.3 |
| **Multi-source PINN and single-station ML location** | Already published. Not claimed as novel | §3, refs [86][87] |

## 5b. Tasks 0 and 1j to 1y: published associators and endpoint controls

**Task 0 (Appendix X).** GaMMA and PyOcto on the same synthetic picks, two seeds. In one-sided
networks the prototype produced zero false events against 0.30-0.43 for the ablation, and at three stations it
recalled 29 of 30 events against 7-9 for the published ones. Three protocols to get there (absolute
threshold, fixed geometry per cell, relative threshold with drawn geometry).

**Tasks 1j to 1y (Appendix Y), executed in September 2026 with rules written before each run.**

| Control | Result | Consequence |
|:--|:--|:--|
| 1j: ablations that break ties without polarisation | 0 false events in one-sided networks on both seeds, like the prototype (which has 1 on one seed) | the one-sided advantage was a ties effect of the discrete score (rule met) |
| 1j: continuous kernel of timing residuals | better than the prototype at 3 and 4 stations | the coincidence ablation was a weak comparator |
| 1j-b: pick error 0.1-0.2 s and velocity error 5% | polarisation added to the kernel contributes in 1 of 9 cells (4 stations, 0.2 s), at a cost in recall | rule not met |
| 1k: 0-2 events per window, absolute threshold | the well-calibrated methods emit in at most one of 20-32 empty windows per cell; the best kinematic score is never worse than those that use polarisation | G1' not met; the absolute threshold solves the empty-window problem |
| Baselines tuned on calibration | GaMMA competes in one-sided and sparse networks; none reaches the target recall in the sparse cells | the advantage at 3 stations is not a usable catalogue |
| Generator audit (Y.4) | it detected all arrivals with exact time even at signal-to-noise 1.1-2.5; polarisation at that signal-to-noise separates with AUC 0.78 | the negatives of 1j-b and 1k were a property of the generator (E19) |
| 1o: consistent generator and likelihood (9d) | read at the nominal budget, PHYS-LR +0.18 to +0.25 over KIN-LR; the independent verification found budget overrun, lack of a negative control, a weak reference and asymmetric domain of the baselines | reading corrected in 1r (E21); with all arrivals detected and 0.5 s error, GaMMA competes |
| 1r: fair repetition, 200 windows, exact FROC curve | PHYS-LR beats its shuffled control by +0.08 to +0.17 at 0.25 false events per window in 5 of 6, and the best of GaMMA and PyOcto by +0.16 to +0.30 in 6 of 6; the shuffled control also beats KIN-LR | R-a and R-d met with the exact curve (R-d not with the registered sweep); R-b not: the shuffled control was not null (E25) |
| 1s: polarisation as a filter on GaMMA | separates true from false events (AUC 0.70-0.96) but gains at most 0.04 of recall | the ceiling is set by the picker: physics must enter the hypothesis test |
| 1t: pick-free waveform migration | projected onto the predicted polarisation (S onto SV) it beats the same without direction (+0.13 to +0.21), its control (+0.06 to +0.18) and GaMMA and PyOcto in 5 of 6 | T-a, T-b and T-c met with the exact curve; T-b failed with the registered sweep |
| 1u: unmodelled medium (delays of 0.15 s, 10 degrees) | the projected migration beats GaMMA by at least 0.05 at both caps in 4 of 6 | U-a and U-b met with the exact curve, not with the registered sweep |
| 1x: coherent moment-tensor detector (WAVE-MT) and greedy deflation (MT-FW) | WAVE-MT does not uniformly beat WAVE-PT (X-a not met); MT-FW beats WAVE-MT in all cells; the comparison with GaMMA was invalidated by E26 | the gain of MT-FW over WAVE-MT is mostly sidelobe suppression; fragile to delays larger than the window |
| 1y: consistent clutter (weak and strong), perturbed condition | MT-FW beats GaMMA by +0.11 to +0.28 (weak) and +0.17 to +0.43 (strong) at 0.25; its control by +0.13 to +0.41; WAVE-PT no longer beats GaMMA with weak clutter | Y-a and Y-b met in both regimes; Y-c only with weak clutter |
| 1w: double-couple sources and controls that take the geometry of another node | per pick, PHYS-LR does not beat its control (1 of 3 and 0 of 3 cells); WAVE-P (S onto SV) falls below the migration without direction; WAVE-PT (P onto the ray, S transverse) beats its control, the migration without direction and GaMMA; replication with 120 new windows: the same | W-a not met; W-b, W-c (WAVE-PT) and replication met; the per-pick positives of 1r and the SV positives of 1t were due to the generator (E24) |

**Honest reading.** What remains standing after three independent verifications: the linear
physics of the problem (the three-component recording of the whole network is linear in the moment
tensor), applied to the waveform before the picker discards information, as a coherent detector
with deflation, detects more events than its control and than GaMMA at equal false-alarm
rate, with double couples, consistent clutter and unmodelled perturbations. It does not pay where all
arrivals are detected with exact time, nor as a filter on a catalogue, nor per pick with realistic
mechanisms, and ray geometry alone does not beat GaMMA with consistent clutter. The design of 1o to 1w was motivated after seeing 1j-1k, several rules were
written after seeing part of their data and the exact curve was adopted after seeing T-b fail
with the registered sweep; this is declared as such. Along the way six defects were found that fabricated or
erased differences: the reading at the nominal budget (E21), a non-deterministic GaMMA (E22), a
sampled FROC curve with jumps (E23), synthetic S waves with the polarisation that the method predicts
(E24), a negative control that was not null (E25) and undetectable clutter picks delivered to GaMMA
(E26). It is still synthetic, with a matched filter close to the true wavelet and parametric
detection. The comparison with a real picker, a mismatched generator and QuakeMigrate was decided
favourably, in synthetics, by Task 2A (below), and the comparison with the block-sparse inversion and the
mechanism-corrected stack on the same grid, by Task 2C; MALMI (1v) and the real data (1n, 1g) are still missing.

**Method components, real picker and QuakeMigrate (Tasks 1p, 1v and 1z; exploratory, seed 999).**

| Test | Result | Consequence |
|:--|:--|:--|
| 1p: PhaseNet and STA/LTA on a mismatched generator (gradient, curved rays, variable wavelet, coda), 160 windows | arrivals picked: PhaseNet 4.1 %, STA/LTA 14.3 %, degrade_y 42.0 %; median error 0.080, 0.078 and 0.033 s; false events per window 2.08, 6.93 and 14.27 | degrade_y, which fed GaMMA in 1o to 1y, was optimistic in recall and timing, and with weak clutter in everything; with strong clutter it delivered some 28 clutter picks per window against 7-8 for the STA/LTA (E27). PhaseNet collapses from domain mismatch (white noise of 0-8 Hz) |
| 1v: QuakeMigrate 1.2.2 through its public API, 60 test windows | recall at ≤0.1/≤0.25: QM 0.15/0.28, QM-3C 0.53/0.53, MT-FW 0.72/0.81; MT-FW − QM-3C +0.19 [−0.00, +0.37] and +0.28 [+0.10, +0.43] | the default practice fails because the rays arrive at 60-80 degrees from the vertical; the gap with QM-3C lies in the windows with two sources; one seed and 47 sources: exploratory |
| 1z: random-field threshold, 7,200 noise windows | the EEC is conservative by a factor of 20 to 60 (about 4 in its on-grid variants); conformal and GEV calibrated (24 of 24); with strong clutter every threshold calibrated on Gaussian noise fails, and recalibrating restores it | the claim that the EEC fixes the rate before seeing data is withdrawn (E28); the threshold is empirical and is recalibrated per regime |
| 1z: certified oracle | identical records in 72 of 72 windows; the energy bound avoids ~4 % of the evaluations and is 3 times slower; the subspace envelope is 4.3 times faster | the pruning gain is negligible; what accelerates is certified vectorisation |
| 1z: station corrections, delays of 0.4 s | rms residual 0.088 s with 20 exact events, 0.173 s with catalogue errors of 1 km and 0.1 s; recall +0.18/+0.24 (4 stations) and +0.30/+0.30 (one-sided) over no correction, within +0.01 to +0.05 of the exact correction | they repair the fragility to delays in 2 of 3 cells; in the 6-station one not even the exact correction gives a significant gain |
| Appendix AA: theory | cell bound and conformal calibration proved; five points where the text said more than the code | E29: radius of 5 km, MT-FW without O(1/k) rate, chi-squared field only under hypothesis (G), η = 1 only for minimisers, curvature equal to the efficient Fisher information |

What is strengthened: the advantage of MT-FW survives, in exploration, against QuakeMigrate in its fair
configuration, and the station corrections estimated from a catalogue repair its fragility to delays in 2 of 3
cells. What is weakened or withdrawn: the false-alarm rate guarantee from random field
theory, the gain of certified pruning over the moment-tensor certificate, the Frank-Wolfe
guarantees for the implemented iteration and the idea that the synthetic picks of 1o to 1y
represented a real picker. None of these tests is confirmatory; Task 2A is.

**Task 2A (confirmatory, rules written on 25 September 2026 before evaluating; independently
verified).** Mismatched generator with fixed 0.3 s delays per network, frozen STA/LTA picks,
station corrections (SC) from 20 catalogued events, 5 networks and 100 windows per cell,
seeds 2024 and 7331, weak and strong clutter. Detail in `tasks_and_real_data_plan.md` and in
`results/task2a/task2a_report.txt`.

| Rule | Weak | Strong | Reading |
|:--|:-:|:-:|:--|
| 2A-a: MT-FW+SC over the best of GaMMA and PyOcto with STA/LTA and QuakeMigrate | met | met | 3 of 3 cells on each seed; +0.09 to +0.49 at 0.1 and +0.17 to +0.46 at 0.25; the best comparator is QuakeMigrate in the 24 cases |
| 2A-b: over the same iteration with the geometry of another node | met | met | the gain requires the true geometry; it does not isolate the radiation pattern (Task 2B) |
| 2A-c: over GaMMA with optimistic picks | met | met | with weak clutter, a tie in the one-sided 8-station cell (+0.01 on both seeds) |
| 2A-d: the SC help | met | not met | with strong clutter the catalogue is contaminated and the SC residual rises to 0.20-0.28 s |
| 2A-e: thresholds frozen on seed 999 | per cell yes; strict no | per cell yes; strict no | recall in 3 of 3 cells; the realised rate is 0.22-0.23 against a limit of 0.20 in one cell and budget |

The verification found no leakage of the truth or scoring asymmetry: the hit rule is the same
for all, the sign of the corrections was checked with the true delays and the almost null recall of
GaMMA and PyOcto is starvation by the picker (with true-time picks and the same wiring they recall 0.88
and 0.69). What is strengthened: the advantage of MT-FW over QuakeMigrate goes from exploratory to confirmed in
synthetics, in both regimes. What bounds it: the per-window intervals ignore the clustering by network
(with a per-network bootstrap, the 12 cells stay above 0 at cap 0.25 and four include 0 at 0.1); about
half of the events are in windows with two sources, and in single-source windows GaMMA with optimistic
picks beats MT-FW+SC at 0.1 in 5 of 6 weak-clutter cells; the sources are synthetic double
couples.

**Task 2B (ablation preregistered on 26 September 2026 before computing it).** The same 1200
windows of 2A, regenerated bit for bit; MT-FW-FREE+SC is the same solver with the radiation operator
replaced by the identity. Rule: MT-FW+SC − MT-FW-FREE+SC ≥ 0.05 at both caps in at least 2 of 3
cells, on both seeds. Detail in `tasks_and_real_data_plan.md` and in `results/task2b/task2b_report.txt`.

| Rule | Weak | Strong | Reading |
|:--|:-:|:-:|:--|
| 2B: MT-FW+SC over the same iteration without radiation operator | met, at the minimum | met | weak: 2 of 3 cells per seed, tie in the one-sided 8-station cell (−0.01/−0.01 and +0.04/+0.0494), ~60 % of replicates; strong: 3 of 3, +0.06 to +0.35, 97-98 % of replicates |

The controls are exact: MT-FW+SC recomputed matches bit for bit what was stored in 120 of 120
windows, and FREE is invariant to the sign of each station (6 of 6). The first attempt stopped at that
control because 18 networks of 2A had been computed with another OpenBLAS kernel (E31); with that kernel pinned
they are reproduced bit for bit, without relaxing the control, and the kernel changes no recall. What is strengthened:
with strong clutter the gain of 2A is due to the radiation operator. What bounds it: with weak clutter the
attribution is fragile and a good part of the advantage over QuakeMigrate does not require the operator
(exploratory), and the ablation does not separate the sign from the polarisation constraint.

**Task 2C (preregistered on 26 September 2026, 15:05 UTC, before implementing it).** It responds to the
independent audit (56/100). Part (i): the same 1200 windows of 2A, regenerated bit for bit, with the
block-sparse inversion (BOMP+SC), the mechanism-corrected stack (MTSTACK+SC) and the group BLASSO
on the grid by fully corrective Frank-Wolfe (GL-FW+SC), all with the same dictionary, grid and
corrections and at the exact predicted times. Part (ii): 25 new networks per cell, seed and regime,
for 30 networks per cell. Detail in `tasks_and_real_data_plan.md` and in `results/task2c/task2c_report.txt`.

| Rule | Weak | Strong | Reading |
|:--|:-:|:-:|:--|
| 2C-a: MT-FW+SC over the best of BOMP+SC and MTSTACK+SC, at 0.25 | met | met | 3 of 3 cells per seed; +0.23 to +0.52 (weak) and +0.08 to +0.45 (strong); the best comparator is BOMP+SC in all cells |
| 2C-b: GL-FW+SC minus MT-FW+SC (descriptive) | below | below | below in 23 of 24 cases, −0.064 to −0.605; duality gap ≤ 1.4e-8; GL-FW+SC at the level of BOMP+SC or somewhat below |
| 2C-c: MT-FW+SC over QM-3C+SC, 30 networks, hierarchical bootstrap | met | met | 3 of 3 cells per seed, +0.22 to +0.35 and +0.099 to +0.295, lower bounds > 0; with Holm, strong 7331 remains at 2 of 3 and the rule is met all the same |

The controls are exact (MT-FW+SC recomputed equal to what was stored in 120 of 120 windows; certificate
equal to dense in 600 of 600) and the new methods were validated on seed 999 against explicit atoms,
noise-free sources and the KKT condition. The code inherited from 2C did not run as written (the
validation was missing, the detection of the OpenBLAS kernel was wrong and the GL-FW tolerance unreachable); it was
corrected before running (E32). What is strengthened: the advantage over QuakeMigrate replicates with the network
as unit of replication, and that of MT-FW over the closest comparators is not due to the moment tensor
with deflation, which was already in Vera Rodriguez et al. (2012), but to its own elements, and according to the
exploratory diagnosis mostly to the peak search; Task 2D refuted that attribution to MT-FW's own
elements (E35). What is weakened: on the grid, the convex problem of the thesis is not the detector that performs (E33); in the
ideal generator with the sources on the nodes (exploratory) the block inversion beats MT-FW, so 2C does not show that MT-FW
beats a sparse inversion on a fine grid or with timing tolerance. With strong clutter and 4
stations MT-FW+SC does not significantly beat BOMP+SC or MTSTACK+SC after Holm.

**Task 2D (preregistered on 26 September 2026, 21:13 UTC, with a SHA-256 hash in
`PREREGISTRO_hashes.txt`, before writing the code).** It responds to the second independent audit.
Part (i): the same 1200 windows of 2A with BOMP-PK+SC, MTSTACK-PK+SC and GL-FW-PK+SC (the methods of 2C
read at the energy peak within ±0.25 s, like MT-FW), MT-FW-EX+SC (without peak search) and MT-FW-NOCAP+SC
(without cap). Part (ii): split-conformal thresholds frozen on 200 empty windows of seed 999 per
cell, at α = 0.1 and 0.25, applied to the windows of 2A. Budget note at 02:41:50 UTC on the 27th,
before running 2024 and 7331 (0.64 h estimated, under the 6 h threshold). Detail in `tasks_and_real_data_plan.md`
and in `results/task2d/task2d_report.txt`.

| Rule | Weak | Strong | Reading |
|:--|:-:|:-:|:--|
| 2D-a: MT-FW+SC over the best of BOMP-PK+SC and MTSTACK-PK+SC, at 0.25 | not met | not met | 2 and 1 of 3 cells (weak), 1 and 0 of 3 (strong); −0.155 to +0.091, all hierarchical CIs contain 0; the best comparator is BOMP-PK+SC in all cells; the advantage over Vera Rodriguez et al. and Anikiev et al. is attributed only to the peak search |
| 2D-b: MT-FW+SC over MT-FW-EX+SC, at 0.25 | met | met | 3 of 3 per seed, +0.064 to +0.412, 12 of 12 positive; lower bound > 0 in 7 of 12; with Holm it is not met |
| Descriptive: without cap; GL-FW-PK+SC | — | — | MT-FW-NOCAP+SC − MT-FW+SC from −0.091 to +0.131, all CIs contain 0; GL-FW-PK+SC 0.11 to 0.30 below MT-FW+SC in 11 of 12 |
| 2D-c: frozen threshold, rate on empty windows ≤ 2α and recall over QM-3C+SC ≥ 0.05 | met per cell; strict no | met | strict in weak: 0.222 > 0.2 in one cell (one false event too many); with Holm it is not met in either |

The controls are exact (truth and records of 2A reproduced in 1200 of 1200 and 120 of 120 windows; the
calibration pipeline reproduces 2A in 60 of 60) and the methods were validated on seed 999 (peak
readings identical to those of MT-FW, exact reduction to 2C without peak, Gram identities ≤ 1.2e-15, KKT and duality gap
of GL-FW-PK). One diagnostic line of `analyze_2d.py` was changed after the part files existed, without
touching rules or tables (hashes before and after in the report). What is weakened: MT-FW is not better than the
block inversion with the same tolerance; its advantage over the pick-based associators and QuakeMigrate belongs
to the family of timing-error-tolerant moment-tensor detectors with deflation. What is
strengthened: the peak search, that is, the tolerance to timing error, is the ingredient that matters, and
with it the block inversion also beats QuakeMigrate (+0.17 to +0.49). What remains open: the
conformal threshold on empty windows does not control the false events in windows with events (1.2 to 2.2 per window in
weak clutter), so there is still no deployable operating point for a false-event budget.

**Task 2E (preregistered on 27 September 2026, 14:44:14 UTC, before writing the code; hash
b40f0c6e…e988 in `PREREGISTRO_hashes.txt` and delivered to the user via the chat; calibration note f = 0.01 at
15:38:20 UTC, before running 2024 and 7331).** It tests G1-M: TS-SFW+SC (group BLASSO with a
shift per station and phase, off the grid by sliding Frank-Wolfe) against BOMP-PK+SC and
against its control with λ = 0 (BOMP-PK-SLIDE+SC), on 30 networks per cell (4800 windows). Detail in
`tasks_and_real_data_plan.md` and in `results/task2e/task2e_report.txt`.

| Rule | Weak | Strong | Reading |
|:--|:-:|:-:|:--|
| 2E-a: TS-SFW+SC over BOMP-PK+SC ≥ 0.05 with lower bound > 0, at 0.25 | not met | not met | 0 of 3 cells on each seed; +0.000 to +0.039, all hierarchical CIs contain 0; with Holm, minimum adjusted p 0.81 (weak) and 0.58 (strong) |
| 2E-b: TS-SFW+SC over BOMP-PK-SLIDE+SC ≥ 0.05, at 0.25 | not met | not met | −0.020 to +0.013: the group penalty is inert with f = 0.01 |
| Independent verification | — | — | reliable negative: no implementation error, no budget disadvantage (−0.014 with 10 times the iterations and 9 starts); realistic ceilings 0.08 to 0.24 above BOMP-PK |

**Exploration (`explore_contrib`, seed 999, exploratory).** Coarse grids, large timing errors,
1D velocity model (oracle or calibrated) and subspace invariant to the pulse: null or negative; the margin
persists with the exact medium. The only candidate that meets the criterion: the per-station cap in windows with 3
to 5 overlapping events (+0.070 at 0.25 in the pilot, +0.068 on new networks). Detail in
`results/explorations/REPORT.txt` and `REPORT_C7.txt`.

**Task 2F (preregistered on 28 September 2026, 03:34:51 UTC; hash 3d8e4c9b…2f5f delivered before
running; seeds 4242 and 8888 never used; code frozen with hashes at 03:49:08 UTC).** It tests G1-R on
96 networks per cell, regime and seed (5760 dense and 960 standard windows per seed). Detail in
`tasks_and_real_data_plan.md` and in `results/task2f/task2f_report.txt`.

| Rule | Weak | Strong | Reading |
|:--|:-:|:-:|:--|
| 2F-a: MT-FW+SC over BOMP-PK+SC in dense windows, at 0.25 | met | met | +0.072 and +0.074 (weak), +0.051 and +0.043 (strong), all lower bounds > 0; with Holm it is met |
| 2F-b: the gain is due to the cap | met | met | MT-FW − MT-FW-NOCAP +0.060 and +0.059; MT-FW-NOCAP − BOMP-PK −0.001 and −0.002 |
| 2F-c: no harm in standard windows (TOST, margin 0.05) | met | met | −0.008 and −0.013; narrow on 8888 (p = 0.0405); negative point estimates with strong clutter |
| Independent verification | — | — | positive confirmed, narrower scope: +0.024 to +0.045 with origin time within 0.5 s; tied to the cap level; mechanism: bounded influence of each station-and-phase block against the model's amplitude error and calibration of the score across windows, not separation of events |

**Task 2I (preregistered, seeds 5151 and 6262; executed and verified; secondary arms
analysed).** HUBER-N+SC, the certificate with a Huber fit per station and phase and threshold 2.7955 fixed by
rule from the noise, against BOMP-PK+SC and MT-FW+SC with cap 6; 120 networks per cell and regime, 5 dense
and 5 standard windows per network; recall at 0.25 false events per window, pooling regimes and cells.

| Rule | 5151 | 6262 | Effect |
|:--|:-:|:-:|:--|
| 2I-a: HUBER-N+SC over BOMP-PK+SC in dense windows | met | met | +0.058 [+0.039, +0.076] and +0.053 [+0.037, +0.070] |
| 2I-b: non-inferior to MT-FW+SC with cap 6 in dense windows (margin 0.03) | met | met | −0.008 [−0.023, +0.006] and −0.002 [−0.018, +0.012] |
| 2I-c: equivalence with BOMP-PK+SC in standard windows (TOST, margin 0.05) | met | met | +0.010 (90 CI [−0.010, +0.027]) and +0.007 (90 CI [−0.013, +0.024]) |
| Independent verification | — | — | positive confirmed; scope: +0.031 and +0.019 with origin time within 0.5 s; equivalence not met with weak clutter on 6262 (−0.032 [−0.064, −0.003]); the recall with all atoms accepted does not rise (it ranks and chooses better, it does not find more events); −0.018 against the cap at 0.1 on 5151; ~9 times the cost of MT-FW; low absolute recall |
| Secondary arms (descriptive, 20 networks per cell and regime, unverified) | — | — | thresholds 2.50 and 3.37 indistinguishable from the registered one (+0.003 and +0.005; +0.001 and +0.006 in dense), threshold 5 lower (−0.013 and −0.023); sd 0.25: cap − BOMP-PK +0.005 and +0.009, HUBER-N − BOMP-PK +0.036 and +0.002; sd 0.8: HUBER-N − BOMP-PK +0.100 and +0.072, HUBER-N − cap −0.055 [−0.113, +0.002] and −0.069 [−0.121, −0.019] |

**Explorations C9 and C10 (seed 999, exploratory; Appendix Y.27).** 600 windows per design, 20 networks per
cell and regime, 95 CI, without multiplicity correction.

| Result | Value |
|:--|:--|
| C9, Cramér-Rao bound of the epicentre: times / all channels without times / joint | 0.4-0.7 km / 4.4-5.6 km / improvement of 4-5 % |
| C9, best locator without times on isolated events (moment-tensor projection) | 7.7-8.4 km median error; 0.29-0.33 within 5 km |
| C9, relocation by evidence without times, MT-FW, dense / standard | −0.116 [−0.153, −0.083] / −0.255 [−0.336, −0.183]; with true factors −0.116 [−0.158, −0.081] |
| C9, SELPOL − MT-FW, dense (after the results, best of some 50) | +0.046 [+0.007, +0.070]; reranking +0.003 at 20 degrees and −0.015 at 45 degrees |
| C10, timing priors with catalogue corrections, against MT-FW | −0.028 to −0.002 (dense), −0.066 to −0.006 (standard) |
| C10, true statics, hard window, against MT-FW | +0.047 [−0.003, +0.082] (dense), +0.098 [+0.032, +0.141] (standard) |
| C10, true statics with prior at the level of the model error | +0.081 [+0.035, +0.117] (dense), +0.146 [+0.065, +0.217] (standard; 0.592 against the ceiling 0.603) |
| C10, isolated events that lose in the first iteration with perfect timing | 37-39 % |
| C10, T3 refinement, standard, registered rule / with origin time | BOMP-PK +0.059 [−0.009, +0.112] / +0.100 [+0.031, +0.168]; MT-FW +0.011 / +0.066 [+0.019, +0.136] |

**Explorations C12, C13 and C14 (seed 999, exploratory; Appendix Y.29).** Registered detector (MT-FW with
cap 6, S0, hard window); C12 with 16 networks per cell and regime, C13 with 12, C14 with 3 development and 12
evaluation networks; 95 CI by hierarchical bootstrap, without multiplicity correction.

| Result | Value |
|:--|:--|
| C12, effective degrees of freedom of the double-couple statistic (tube formula); simulated against theoretical 0.01 quantile | 4.4-4.5; 14.28 against 14.17 |
| C12, ratio of double-couple / full thresholds at equal false-event rate, noise only, ≤ 0.1 / ≤ 0.25 | 0.906 [0.883, 0.928] / 0.898 [0.886, 0.914] |
| C12, double couple − full, recall at 0.25 with origin time, standard / dense | +0.005 [−0.021, +0.033] / −0.004 [−0.018, +0.011] |
| C13, false events per window with events against the event-free twin, above t95, standard / dense | 12 [9, 18] / 45 [33, 63] times |
| C13, decoy D1p with batches of 120, realised FDR at q = 0.2, epicentre only / with origin time | 0.14-0.16 / 0.28-0.37 |
| C13, D1p, realised false events per window for c = 0.25 (epicentre only), standard / dense | 0.21 [0.17, 0.25] / 0.23 [0.19, 0.28]; recall 0.453 / 0.245 against 0.436 / 0.253 for the labelled point |
| C14, diagnosis on development: origin-time false events that are well-located events with the origin shifted, standard / dense | 55-58 % / 80-95 %; median origin error 0.76 / 0.86 s |
| C14, ratio of false events against the twin at q95, BASE and R2+R3, standard / dense | 6.23 and 2.52 / 10.04 and 4.24 |
| C14, remedies − BASE, recall at 0.25 with origin time, standard / dense | +0.015 to +0.016 (CI including 0) / −0.011 to −0.002 |
| C14, overlapping events at any score: hard exclusivity / FIR deflation / R2+R3 | −0.103 / −0.098 / −0.095 (BASE 0.400) |

What is weakened: no convex formulation of the thesis beats the block search with the same
tolerance; G1-M is refuted. What is strengthened: G1-R has two preregistered tests on never-used
seeds, the per-station cap (2F) and the Huber fit with level fixed by rule from the noise
(2I), which removes the dependence on a tuned constant, with an identified mechanism. What is bounded: with
amplitude dispersion 0.8 the Huber fit falls below the cap, and physics without travel times and
soft timing priors do not close the margin (C9, C10). What remains open: the transfer to another
generator and to real data; and the margin of 0.08 to 0.24 below the realistic ceiling, which neither the velocity model
nor the pulse explains, of which no amplitude treatment captures more than a quarter (exploration
C8) and which explorations C9 and C10 split into a station-term calibration part and an
amplitude part. After 2L, explorations C12 to C14 close the double-couple restriction and the
self-interference remedies as recall levers and point to the origin time of the well-located
detections as the next target.

## 5c. Task 1, executed: polarisation with real data

CM network of the SGC, ISC reviewed catalogue, PhaseNet, five stations. Detail in Appendix T.1.

**Real clutter is not coherent with the preceding catalogued event, where it can be measured.** At ZAR the fraction of false picks
coherent with their neighbouring event, corrected for chance, is $\eta=0.054\ [-0.049,0.164]$ over 159
picks: below the registered threshold of 0.20. The non-coherent clutter assumption of the
pre-studies holds in that sense; coherence with uncatalogued sources remains to be measured.

**But single-station polarisation is only usable at some sites.** The prior control
—does the P polarisation point to catalogued events?— passes at one station out of five:

| Station | P control | status |
|:--|:--|:--|
| ZAR | 0.789 (n=109) | passes |
| CHI | 0.256 (n=78), equal to chance | fails; metadata, rotation, reflection and instrument change ruled out |
| RUS | 0.244 (n=45) | fails |
| OCA | 0.303 (n=33) | fails; not distinguishable from chance |
| SJC | 0.474 (n=19) | inconclusive: fewer than 30 picks |

**What it means for the thesis.** The polarisation channel is not available across the whole network: its
reliability is a measurable property of the site, or of the incidence geometry (at CHI the P is
at chance and the S hits 0.677, compatible with nearly vertical arrivals; Task 1l). C-IV enters per
station with a reliability weight $\hat r_s$ estimated with catalogued events, with the same
algebra as $\eta$ (equations 9d-9e of the proposal).

**Task 1e (synthetic, Appendix T.2).** With the measured reliability: naive C-IV worsens when
the stations that fail have no directional information; switching it off where it fails does not fix it; and
a **neutral** value measured with catalogued events is the correct treatment. The advantages of 1e are
measured against the TOL ablation, whose false events are largely ties (Task 1j): what stands firm
is the ordering between treatments, not the size of the advantage.

**Scope of the clutter verdict.** $\eta$ measures coherence with the preceding **catalogued**
event. Coherence with uncatalogued microearthquakes, which are the physically realisable
impostors of Appendix T, is not measured by it and remains registered (Task 1m).

**What it cost.** Five design defects, none physical, registered as E10 to E14 in the
manifest. The most serious: the registered threshold of 0.20 on the raw fraction was below
what a clutter without any coherence gives in expectation (0.22-0.33), so it did not separate
coherent from random clutter.

## 6. The limitation that decides the project

All the evidence for G1 comes from synthetic generators and from a discrete prototype, not from the
proposed model; Appendix Y showed that the conclusion depends on whether the generator simulates detection and
timing consistently with the signal-to-noise ratio.
Two statements must be separated:

$$\text{(i)}\ \mathcal{J}_{\mathrm{fis}}\succ\mathcal{J}_{\mathrm{cin}}\quad\text{(consistent synthetic, coherent MT detector: yes; per pick: no)}
\qquad
\text{(ii)}\ \mathcal{J}_{\mathrm{fis}}\succ\mathcal{J}_{\mathrm{GaMMA}}\quad\text{(consistent synthetic, coherent MT detector: yes; QuakeMigrate, and GaMMA and PyOcto with STA/LTA: yes in 2A, synthetic; block-sparse inversion and mechanism stack on the same grid: yes in 2C, synthetic; MALMI and real data: pending)}$$

where $\mathcal{J}_{\mathrm{cin}}$ is the best continuous kinematic score, which replaces the
coincidence ablation as the reference. The argument that no published associator uses
polarisation is true for those for local networks with machine-learning picks, but not for
NET-VISA, which models the azimuth per station; the proposal positions itself against it (Section 3.2).

Additional declared limitation: the sixteen pre-studies use **a single synthetic generator**.
The sufficiency ladder could be a property of the simulator. Task 2A uses a generator mismatched
with the detectors (gradient medium, curved rays, variable wavelet, coda), but derived from the
first one, with its same clutter law and with white noise; it does not replace a second independent
generator (8b) or the real data.

## 7. Actions, in order

| # | Action | Decides | Effort |
|:--|:--|:--|:--|
| ~~0, 1, 1e, 1j, 1j-b, 1k, 1o, 1r, 1s, 1t, 1u~~ | **COMPLETED**, see §5b and §5c | | done |
| ~~2A~~ | **EXECUTED AND VERIFIED**: 2A-a, 2A-b and 2A-c are met in both regimes; 2A-d only with weak clutter; 2A-e is met per cell and not in the strict reading (§3, §5b) | The advantage of MT-FW survives outside the generator that produced it, in synthetics; the effective rival is QuakeMigrate | done |
| ~~2B~~ | **EXECUTED AND VERIFIED**: the rule is met in both regimes; with weak clutter only at the minimum (2 of 3 cells, ~60 % of replicates), with strong clutter in 3 of 3 (97-98 %); it does not separate the sign from the polarisation (§3, §5b) | The gain of 2A is attributed to the radiation operator, clearly in strong clutter; in weak clutter a good part of the advantage over QuakeMigrate does not require it | done |
| ~~2C~~ | **EXECUTED AND VERIFIED**: 2C-a and 2C-c are met in both regimes (2C-c with 30 networks per cell and hierarchical bootstrap; with Holm, strong 7331 2 of 3); 2C-b: the group BLASSO on the grid falls below MT-FW in 23 of 24 cases (§3, §5b) | The advantage of MT-FW over the closest comparators is robustness to grid and timing errors; the convex problem is framework and certificate, not the detector | done |
| ~~2D~~ | **EXECUTED AND VERIFIED**: 2D-a is not met in either regime (with the same peak search BOMP-PK+SC equals MT-FW+SC, −0.155 to +0.091, all CIs contain 0); 2D-b is met (the peak search matters, +0.06 to +0.41, 12 of 12); the cap does not matter; 2D-c is met per cell, not in the strict reading with weak clutter, and the threshold on empty windows is not deployable; with Holm no rule (§3, §5b) | The advantage of MT-FW over Vera Rodriguez et al. and Anikiev et al. is only the peak search; the empirical advantage belongs to the family of timing-error-tolerant moment-tensor detectors, not to MT-FW | done |
| ~~2E~~ | **EXECUTED AND VERIFIED**: 2E-a and 2E-b are not met in either regime (TS-SFW+SC − BOMP-PK+SC +0.000 to +0.039, all CIs contain 0); reliable negative according to the verification (§3, §5b) | G1-M refuted in synthetics: the off-the-grid BLASSO with shifts does not beat BOMP-PK | done |
| ~~2F~~ | **EXECUTED AND VERIFIED**: 2F-a, 2F-b and 2F-c are met in both regimes and seeds (+0.043 to +0.074 in dense windows, also with Holm); the verification narrows the scope (+0.024 to +0.045 with origin time; tied to the cap level; per-block bounded-influence mechanism) (§3, §5b) | G1-R confirmed in synthetics: the robust certificate with per-station cap is the methodological contribution, a modest one | done |
| ~~2I~~ | **EXECUTED AND VERIFIED; secondary arms analysed**: 2I-a, 2I-b and 2I-c are met on both seeds (HUBER-N+SC − BOMP-PK+SC +0.058 and +0.053 in dense windows; non-inferior to MT-FW+SC with cap 6; equivalent to BOMP-PK+SC in standard); the verification narrows the scope (+0.031 and +0.019 with origin time; equivalence not met with weak clutter on 6262; it ranks and chooses better, it does not find more events; ~9 times the cost of MT-FW); descriptive secondary arms: level indifferent between 2.5 and 3.4, below the cap with amplitude dispersion 0.8 (§3, §5b) | G1-R with two preregistered tests: a level fixed by rule from the noise replaces the calibrated cap in this synthetic design with the nominal dispersion; the contribution is the measured design result, a modest one | done |
| ~~C9, C10~~ | **EXPLORATORY, DONE** (seed 999): evidence without travel times is information-limited (4.4-5.6 km against 0.4-0.7 km for the times) and relocating with it loses recall; soft timing priors do not help with catalogue corrections; with true statics the detector reaches the ceiling on isolated events (§3, §5b) | The margin splits into calibration of station terms and amplitude evidence; no other location principle is needed | done |
| ~~2L~~ | **EXECUTED AND VERIFIED**: 2L-a to 2L-d and the secondary on BOMP-PK are met on both seeds (MT@S2+T3 − MT@S0+T3 +0.057 and +0.062 in standard windows and +0.038 and +0.030 in dense, with the origin time required); the verification narrows the scope (gain concentrated in the one-sided 8-station cell, +0.108 and +0.113; nothing demonstrated with 4 stations; S2 worse than S0 in about 22 % of the networks; conditional on T3; catalogue errors of known size; 0.28 to 0.61 of the oracle gain) (Appendix Y.28) | G1-C confirmed in synthetics: station terms estimated without picks inside the refined detector, a modest calibration result | done |
| ~~C12, C13, C14~~ | **EXPLORATORY, DONE** (seed 999): the double-couple restriction lowers the threshold by 10 % without changing the recall (null law from the tube formula); decoy calibration with permuted moveout reproduces without labels the operating point with the epicentre-only rule, not with origin time; the self-interference remedies halve the induced false events without moving the recall; the largest class of origin-time false events are well-located detections of real events with the origin shifted by a median of 0.8 s (Appendix Y.29) | None is proposed for confirmation as a recall contribution; the origin time of the well-located detections becomes the target of the calibration | done |
| ~~R0~~ | **EXECUTED, NEGATIVE**: −0.383 [−0.481, −0.310] against PhaseNet with PyOcto on one hour of Ridgecrest 2019 (Appendix Y.30); closure of M1 − M3, declared post hoc, in progress | G1 does not transfer to a real dense sequence with the registered timing and score | done |
| ~~2M~~ | **EXECUTED AND VERIFIED**: 2M-a +0.228 [+0.148, +0.298] and 2M-b +0.241 [+0.179, +0.293] on real noise with transients (primary day); +0.057 and +0.120 on near-Gaussian days (secondary); the cap's gain is concentrated in the hours with transients, HUBER-N's persists (Appendix Y.31) | G1-R holds on independent real noise of one network | done |
| ~~C18~~ | **EXPLORATORY, DONE** (Ridgecrest CAL and seed 999): no 1-D medium improves the readings; the cap erases real-event evidence at the true hypocentre and LOG6 restores it without losing robustness, but no variant locates the events (Appendix Y.32) | The R0 failure is the reading of arrivals in a dense sequence, not the medium | done |
| **Reading of arrivals for dense real sequences** | P-weighted or envelope reading, higher-frequency matched filter, joint relocation with a coherence check, LOG6 statistic; developed on Ridgecrest CAL only | Whether any candidate passes the registered gate (≥ 5 of 16 CAL events with origin time at 30 false/h on N8 or U8) that unlocks a registered test on the unread Hi-net swarm | to be estimated |
| **Station-term calibration, next steps** | Sensitivity of the S2 estimator to catalogue errors of unknown size, dependence on the geometry and combination with the bounded-influence statistic, with the true-statics oracle as a bound; prior exploration of whether the S2 terms, a better timing model or a treatment of the early origin bias reduce the class of well-located detections with the origin shifted (C14); afterwards, amplitude evidence for weak arrivals | Whether G1-C holds outside the known-catalogue assumption and the one-sided cell | to be estimated |
| **1n** | Real PhaseNet residuals against the catalogue (data already downloaded) and calibration of $\eta$ | Whether the RSNC is in the regime where polarisation can pay | 1 wk |
| **1l** | Orientation audit (OrientPy) and incidence control at CHI, OCA, SJC, ZAR | Whether the failure at CHI is geometry, orientation or site | 3-5 days |
| **1f** | Robustness of the ZAR clutter verdict by subpopulation | Whether the verdict depends on the sampling | 1 hour |
| **1d** | Band of the proposal (0.1-2 Hz) | Whether C-IV operates in the declared band | hours |
| **1m** | Relabelling with the SGC catalogue and SeisMonitor; templates on the clutter | Whether there are uncatalogued impostors | 1-2 wk |
| **1c** | Reliability map $\hat r_s$ on the CM network | Real $S_{\mathrm{eff}}$ of the network | 3-5 days |
| **1h** | Per-station likelihood (9d-9e) with continuous reliability $\hat r_s$ and real pick errors | The confirmatory form of C-IV (the synthetic form was already tested in 1o) | 1 wk |
| **1p** | Remainder: traces at 100 Hz with real RSNC noise, and REAL and GENIE as baselines (PhaseNet and STA/LTA on the mismatched generator already run, §5b) | Whether the real picker on realistic noise changes the measured regime | 2-3 wk |
| **1v** | Remainder: MALMI as published, and the block-sparse inversion of Vera Rodriguez et al. (2012) on a finer grid and with full-waveform dictionaries (QuakeMigrate already run, confirmed in 2A; the block inversion and the mechanism stack on the common grid, in-house implementation, run in 2C with exact times and in 2D with the peak search, where the block inversion equals MT-FW) | Whether the pick-free advantage survives MALMI and a sparse inversion on a fine grid | 2-3 wk, plus the sparse inversion |
| **1g** | Real pilot by decimation of a dense network | G1 on real data | 2-3 months |
| **1i** | Clutter in an active sequence | The dangerous case of Appendix T | 1 wk |
| 8b | Second generator: stratified medium, real coda clutter, real pick errors | Whether the ladder belongs to the problem or to the simulator | 3-4 wk |
| 3 | Band supported by the velocity model | Gate for all real waveform work | 3-4 wk |

1n, 1l, 1f and 1d use data already downloaded and go first; 1n is the most important, because
Appendix Y showed that the value of polarisation depends on how much kinematic information the
picker leaves: its real detection rate and timing error decide which regime the RSNC is in.

## 8. Framing recommendation

For the defence, lead with the **observability certificate** and the **sufficiency ladder**
reindexed by effective kinematic information, and present as the framework of the method the sparse
recovery of a measure with moment-tensor amplitudes, with its certificate, and as the detector MT-FW
(design document; Task 2C), with its synthetic
evidence (Appendices Y.11 and Y.12: controls that break the geometry, consistent clutter, two
seeds, new windows), its exploratory test against QuakeMigrate, its preregistered tests
(Tasks 2A to 2F) and its pending test on real data. As the methodological contribution, present the
robust certificate with per-station cap, with its confirmed effect in dense sequences (Task 2F), its
modest size and its mechanism, and the negatives of the benchmark, among them G1-M (Task 2E), without presenting
the off-the-grid BLASSO as the method that performs. Appendix Y must not be hidden: it shows that the programme audits its positives and its
negatives alike, and that both changed with the audit.

If the real data (1n and 1g) show that this detector contributes over its own control, over the
published associators and over the migration-based locators, that is the central contribution. If they show it does not, the thesis
delimits when three-component physics does not pay, and the certificate, the per-station reliability
methodology and the certified search remain intact. Of the latter, what is upheld is the exactness (records
identical to the exhaustive search), not a pruning gain over the moment-tensor certificate.

## 9. Note on the working method

During development **forty-five errors** were documented (E1 to E45): three of claims in the text
pointed out by the fourth independent audit (E40, "certificate" used for a statistic that is not one;
E41, a declared dependence removed together with the arms that measure it without reporting it; E42, the real test G1'
without an arm for the polarisation assumption that the real data contradict), thirty of experimental
design (E38, a hit rule without origin time that is lax in dense windows, among them),
one of positioning (E30), one of numerical reproducibility (E31), one of implementation
detected before running (E32), one of a claim about the method (E33), two of registration (E34 and E36, the
unverifiable date of the preregistrations of 2A to 2C and an undeclared edit of `analyze_2b.py` that
changed no rule or figure) and three of attribution (E35; E37, a mechanism read in the pilot more broadly
than the confirmed effect of Task 2F supports; and E39, an amplitude factor that the generator draws
in each window read as a site term of the station, with a plan for a certificate with site terms
that exploration C8 does not support), and three in the explorations after 2L (E43 and E44, in C13, an event-free
twin that retained remnants of the events and decoys that saw less clutter than the target; E45,
execution and code incidents around the freeze of C14), none with an effect on a figure
used. Several of
them reversed or annulled a result: constraints turned into no-ops, undersampling of the
sensitive parameter, true source on a grid node, thresholds below chance and, the most costly,
ties of a discrete score read as a physical effect (E15), a weak comparator against an
overly favourable generator (E16) and a generator that simulated the kinematic channel at infinite
signal-to-noise, which produced a false negative on polarisation (E19). The positive that
replaced it was written before auditing it with the same rigour and was overestimated (E21); on
repeating it there appeared a non-deterministic baseline (E22), a sampled FROC curve with jumps that
could create or erase differences (E23), a generator that radiated the S wave exactly in the direction that
the method predicts (E24), a negative control that retained part of the signal (E25) and clutter
picks that no picker would produce, delivered only to the associator (E26). On testing the components
of the method there appeared a synthetic pick model more optimistic than real pickers (E27), a
false-alarm rate guarantee written before measuring it (E28) and five theoretical statements stronger
than the code (E29). The novelty audit found that the proposal exaggerated the novelty of
sparse joint detection and location with moment tensors, anticipated by Vera Rodriguez et al.
(2012) and GRiD MT (E30), and the contribution was repositioned (§3). The independent verification of Task 2A
found no implementation errors; its deviations from the preregistration (one catalogue per network and regime,
four calibration files from an earlier version of the code) are declared. The exact control of
Task 2B detected that 18 networks of 2A had been computed with another OpenBLAS kernel (E31). The review
of the code inherited from Task 2C, before running it, found that it could not run as written:
its validation was missing, it detected the OpenBLAS kernel wrongly and it required an unreachable tolerance (E32). And
Task 2C showed that the proposal presented the group BLASSO as the method that performs when, on the
grid, it falls below its greedy approximation (E33). The second independent audit found that the
time of Task 2B in the preregistration had been written without consulting a clock and was wrong (E34); since
Task 2D each entry carries a clock time and a SHA-256 hash. And Task 2D showed that the advantage of MT-FW over
the block inversion, which 2C attributed to MT-FW's own elements, was due to a comparison with
unequal temporal tolerance: with the same peak search the block inversion equals it (E35). The
independent verification of Task 2E found no implementation errors in its negative. That of Task
2F confirmed its positive and found that the reading of the mechanism inherited from the pilot, "the cap resolves the
dense sequences", was broader than what the data support (E37), and that the hit rule without origin
time, suitable for windows of 0 to 2 events, inflates the absolute gain in dense windows (E38). The
exploration C8 found that the "per-station site amplification" mechanism attributed to the station a factor
that the generator draws in each window, and that the plan for a certificate with site terms has no
support (E39). The independent verification of Task 2I confirmed its positive without finding implementation
or provenance errors, regenerated parts across the two interruptions of the run and bounded its
scope (hit rule with origin time, equivalence per regime, cost), without adding an error to the
log. The
forty-two errors are in
`artifact_manifest_and_error_log.md` with cause, detection and correction, together with the checks that
would have caught them.

That log is part of the contribution: the proposal incorporates a **mandatory design-validity
step** in the Tier-1 pilot, in which each constraint must be shown to be active on a case
constructed for it to be so, before running any confirmatory condition.
