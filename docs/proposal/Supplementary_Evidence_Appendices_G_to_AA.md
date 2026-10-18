---
title: "Supplementary Evidence Volume: Appendices G to AA"
subtitle: "Evidence volume of the doctoral research proposal Timing-Tolerant, Bounded-Influence Moment-Tensor Detection of Overlapping Earthquakes: A Preregistered Benchmark and a Real-Data Test Plan"
author: "Juan Jose Gordillo Cordoba"
date: "October 2026"
---

# Supplementary Evidence Volume: Appendices G to AA

This is the evidence volume of the doctoral research proposal *Timing-Tolerant, Bounded-Influence Moment-Tensor
Detection of Overlapping Earthquakes: A Preregistered Benchmark and a Real-Data Test Plan* (main document:
`Doctoral_Proposal.md`). It contains Appendices G
to AA of the proposal in full: the numerical pre-studies, their statistical audit, the
published-associator baselines, the catalogue-endpoint controls with the eight preregistered synthetic
tasks, the explorations that preceded the sixth, the seventh and the eighth and three that followed the eighth, the certified pruning of the hypothesis space and the theory of the group BLASSO detector.
Appendix lettering and the numbering of equations, tables and references are those of the main
document. References [1] to [150] are listed in Appendix F of the main document, and Sections 1 to 23
and Appendices A to F cited here are in the main document.

## Contents

- [Appendix G. Pre-Study: Identifiability and Station Count](#appendix-g-pre-study-identifiability-and-station-count)
- [Appendix H. Pre-Study: Where Physical Evidence Wins](#appendix-h-pre-study-where-physical-evidence-wins)
- [Appendix I. Pre-Study: Observability, Coercivity, and Network Geometry](#appendix-i-pre-study-observability-coercivity-and-network-geometry)
- [Appendix J. Pre-Study: Which Constraints Carry the Information](#appendix-j-pre-study-which-constraints-carry-the-information)
- [Appendix K. Pre-Study: Medium-Independent Constraints, and a Negative Result on Overlap](#appendix-k-pre-study-medium-independent-constraints-and-a-negative-result-on-overlap)
- [Appendix L. Pre-Study: Behaviour under Velocity-Model Error](#appendix-l-pre-study-behaviour-under-velocity-model-error)
- [Appendix M. Pre-Study: The Minimal Sufficient Set](#appendix-m-pre-study-the-minimal-sufficient-set)
- [Appendix N. Pre-Study: Cost, and Network Design](#appendix-n-pre-study-cost-and-network-design)
- [Appendix O. Pre-Study: The Sufficiency Ladder](#appendix-o-pre-study-the-sufficiency-ladder)
- [Appendix P. Pre-Study: One Station, Several Events](#appendix-p-pre-study-one-station-several-events)
- [Appendix Q. Pre-Study: The Certificate under Interfaces](#appendix-q-pre-study-the-certificate-under-interfaces)
- [Appendix R. Pre-Study: Detection and Association are One Axis](#appendix-r-pre-study-detection-and-association-are-one-axis)
- [Appendix S. Pre-Study: Observability Maps Where Discrimination is Hard](#appendix-s-pre-study-observability-maps-where-discrimination-is-hard)
- [Appendix T. Pre-Study: Coherent Nearby Clutter](#appendix-t-pre-study-coherent-nearby-clutter)
- [Appendix U. Pre-Study: Association Without a Velocity Model](#appendix-u-pre-study-association-without-a-velocity-model)
- [Appendix V. Ruled-Out Approaches and Methodological Failure Analyses](#appendix-v-ruled-out-approaches-and-methodological-failure-analyses)
- [Appendix W. Statistical Status of the Pre-Studies](#appendix-w-statistical-status-of-the-pre-studies)
- [Appendix X. Published Associators as Baselines](#appendix-x-published-associators-as-baselines)
- [Appendix Y. Controls for the Catalogue Endpoint](#appendix-y-controls-for-the-catalogue-endpoint-ties-timing-error-empty-windows-and-tuned-baselines)
  - [Y.18 The confirmatory test (Task 2A)](#y18-the-confirmatory-test-task-2a)
  - [Y.19 The radiation operator inside the same solver (Task 2B)](#y19-the-radiation-operator-inside-the-same-solver-task-2b)
  - [Y.20 The closest comparators, the convex problem, and the network as the unit (Task 2C)](#y20-the-closest-comparators-the-convex-problem-and-the-network-as-the-unit-task-2c)
  - [Y.21 Equal timing tolerance, the peak-search ablation, and a frozen conformal threshold (Task 2D)](#y21-equal-timing-tolerance-the-peak-search-ablation-and-a-frozen-conformal-threshold-task-2d)
  - [Y.22 The off-grid time-shift BLASSO (Task 2E)](#y22-the-off-grid-time-shift-blasso-task-2e)
  - [Y.23 Where the headroom is not: seven exploratory candidates on the calibration seed](#y23-where-the-headroom-is-not-seven-exploratory-candidates-on-the-calibration-seed)
  - [Y.24 The per-station cap in dense sequences (Task 2F)](#y24-the-per-station-cap-in-dense-sequences-task-2f)
  - [Y.25 Site terms versus bounded influence: an exploratory study (C8)](#y25-site-terms-versus-bounded-influence-an-exploratory-study-c8)
  - [Y.26 A noise-calibrated bounded-influence statistic (Task 2I)](#y26-a-noise-calibrated-bounded-influence-statistic-task-2i)
  - [Y.27 Evidence that does not use travel times, and the timing budget: two exploratory studies (C9, C10)](#y27-evidence-that-does-not-use-travel-times-and-the-timing-budget-two-exploratory-studies-c9-c10)
  - [Y.28 Station terms estimated without picks, and the refined detector (C11 and Task 2L)](#y28-station-terms-estimated-without-picks-and-the-refined-detector-c11-and-task-2l)
  - [Y.29 After Task 2L: three exploratory studies on the calibration seed (C12, C13, C14)](#y29-after-task-2l-three-exploratory-studies-on-the-calibration-seed-c12-c13-c14)
  - [Y.30 Task R0: the first registered test on real data (Ridgecrest, 2019)](#y30-task-r0-the-first-registered-test-on-real-data-ridgecrest-2019)
  - [Y.31 Task 2M: bounded influence with independent real noise (Alto Tiberina, 2010)](#y31-task-2m-bounded-influence-with-independent-real-noise-alto-tiberina-2010)
  - [Y.32 After R0: what fails on real data (exploration C18)](#y32-after-r0-what-fails-on-real-data-exploration-c18)
  - [Y.33 Task R1: a hybrid with pick-based association on the Hi-net swarm](#y33-task-r1-a-hybrid-with-pick-based-association-on-the-hi-net-swarm)
- [Appendix Z. Certified Pruning of the Hypothesis Space, and the Role of p4est](#appendix-z-certified-pruning-of-the-hypothesis-space-and-the-role-of-p4est)
- [Appendix AA. Theory of the Group BLASSO Detector](#appendix-aa-theory-of-the-group-blasso-detector)

---

# Appendix G. Pre-Study: Identifiability and Station Count

Code: `c1_identifiability_prestudy.py`. Setting: sources on a
horizontal line at 8 km depth, surface stations over a 60 km aperture, homogeneous medium with
$V_S=3.0$ km/s, ideal low-pass pulse at $f_{\max}=2$ Hz, $1/d$ spreading, isotropic radiation, atoms
unit-normalized so configurations are compared at equal energy. Declared simplifications:
homogeneous, one-dimensional in source coordinate, mechanism-free, depth fixed.

**Part 1, noiseless certificate.** The vanishing-derivative precertificate for the two-spike support
was satisfied at every separation tested, down to 20 m, for every configuration and velocity ratio.
In a multi-station geometry the noiseless problem is not separation-limited; the operative limit is
noise.

**Part 2, resolution with origin time fixed.** Cramer-Rao resolution limit at SNR 20, 12 stations:

| $V_P/V_S$ | P only | S only | P+S | $\kappa$(P+S) |
|:--|:--|:--|:--|:--|
| 1.50 | 0.264 km | 0.177 km | 0.198 km | 0.132 |
| 1.67 | 0.294 | 0.177 | 0.199 | 0.133 |
| 1.80 | 0.317 | 0.177 | 0.200 | 0.134 |
| 1.95 | 0.344 | 0.177 | 0.198 | 0.132 |

The S-only column is exactly invariant as it must be, validating the implementation. P-only degrades
monotonically with the ratio, as $\partial t/\partial x\propto1/V$ requires. The joint operator is
flat in the ratio and 12-13% worse than S alone at equal energy.

**Part 3, origin time free.** At $V_P/V_S=1.73$:

| $S$ | P only | S only | P+S |
|:--|:--|:--|:--|
| 1 | not identifiable | not identifiable | 0.376 km |
| 2 | 0.248 km | 0.144 km | 0.158 km |
| 3 | 0.355 | 0.206 | 0.245 |
| 12 | 0.305 | 0.176 | 0.200 |

At one station both single-branch operators are structurally degenerate; the two-branch operator is
not. The effect is on rank, not on a resolution constant. The ratio sweep at $S=2$ was flat to 0.8%.

 Appendix G establishes $S^\ast_{2,\mathrm{reduced}}=1$ and $S^\ast_{1,\mathrm{reduced}}=2$ only for its depth-fixed, one-dimensional source-coordinate model; Section 9.2 derives $S^\ast_{2,\mathrm{full}}=3$ for the full four-parameter arrival model. The velocity-contrast conjecture, the measured constant $\kappa\approx0.105$-$0.133$, and fitted branch weights are also retained.

---

# Appendix H. Pre-Study: Where Physical Evidence Wins

Code: `c0_coherence.py`, `c0_clutter.py`. Setting: ring of surface
stations at 25 km, source at 8 km depth, homogeneous medium, $V_S=3$ km/s, $V_P/V_S=1.73$, Ricker
pulse at $f_{\max}=2$ Hz. The physical atom carries P and S arrivals with longitudinal and
orthogonal polarization and $1/d$ spreading; the arrival-time atom carries the same $t_P,t_S$ with
no polarization and no amplitude. Both unit-normalized, so the comparison is at equal energy.
Declared simplifications: homogeneous medium, straight rays, no mechanism, no attenuation.

**Study 1, resolution between two true sources.** Separation at which coherence between the two
hypotheses falls below 0.5:

| Stations | horizontal, waveform | horizontal, arrival time | depth, waveform | depth, arrival time |
|:--|:--|:--|:--|:--|
| 1 | 0.30 km | 0.33 km | 0.83 km | 0.91 km |
| 3 | 0.41 | 0.47 | 0.83 | 0.91 |
| 12 | 0.41 | 0.47 | 0.83 | 0.91 |

Gain 6-14%, flat in station count. Radial and transverse offsets were identical, as ring symmetry
requires. Physical evidence does not materially sharpen resolution.

**Study 2, rejection of false hypotheses under clutter.** One true source at SNR 8 plus $n$
randomly polarized false arrivals. Margin is the gap between the score at the true source and the
best score among grid hypotheses more than 5 km away, in units of the standard deviation of the
score field.

| False arrivals | arrival-time margin | best false / true | waveform margin | best false / true |
|:--|:--|:--|:--|:--|
| 0 | 3.93 | 0.50 | 16.39 | 0.25 |
| 15 | 3.18 | 0.58 | 15.94 | 0.25 |
| 30 | 1.86 | 0.75 | 14.52 | 0.27 |
| 60 | 1.33 | 0.83 | 11.64 | 0.30 |

The arrival-time margin degrades by a factor of three and the best false hypothesis approaches the
true one; the physical margin degrades by 29% and no false hypothesis exceeds 30% of the true score.

These results identify discrimination margin, rather than resolution, as the principal value of the wavefield; whether the margin becomes fewer false events than the best kinematic score achieves is examined in Appendix Y.
G1 is therefore the primary hypothesis, while E1 records the limited resolution gain as a registered null.
Polarization is treated as a primary constraint in Section 6.4, and the target regimes are high
false-arrival rates and sparse networks. The resulting interpretation is that each physical channel adds
constraints to the feasible set and that clutter is rejected when it cannot satisfy them jointly.

---

# Appendix I. Pre-Study: Observability, Coercivity, and Network Geometry

Code: `ci_observability.py`, `ci_time.py`. Setting: homogeneous
half-space, $V_S=3$ km/s, $V_P/V_S=1.73$, source at 8 km depth, Ricker pulse at $f_{\max}=2$ Hz,
three-component records with P and S branches, polarization and $1/d$ spreading;
$\theta=(x,y,z,\tau)$; noise standard deviation 0.02.

**Design, and why it is not circular.** Three distinct objects are computed. The *geometric*
quantities come from station directions alone - the smallest eigenvalue of the direction covariance
$\frac1S\sum_s \mathbf{n}_s\mathbf{n}_s^{\mathsf T}$, the maximum azimuthal gap, and the station
count - without forming a Jacobian or touching a waveform. The *coercivity constant* $\gamma$ is the
smallest eigenvalue of the whitened Fisher matrix of the full physical operator. The *empirical
error* is the median epicentral error over 25 Monte-Carlo nonlinear least-squares fits to noisy
data. Correlating $\gamma$ with the Cramer-Rao bound would be tautological; correlating a purely
geometric quantity with either is not.

**Result 1: geometry predicts coercivity and error.** 96 geometries across six families - random,
one-sided, ring, clustered, linear and near-field - with 3 to 16 stations.

| Predictor | vs. coercivity $\gamma$ | vs. empirical epicentral error |
|:--|:--|:--|
| **Directional coverage (control-theoretic)** | **$\rho=+0.815$** | **$\rho=-0.911$** |
| Azimuthal gap (classical heuristic) | $\rho=+0.446$ | $\rho=-0.657$ |
| Station count | $\rho=+0.393$ | $\rho=-0.480$ |
| Cramer-Rao bound (validation only) | - | $\rho=+0.982$ |

The last row is a check on the Monte-Carlo implementation rather than a scientific result: the
empirical error tracks the Cramer-Rao bound almost exactly, as it must.

**Result 2: the classical heuristic misranks geometries.** Median values by family:

| Family | coverage | azimuthal gap | $\gamma$ | CRB (km) | error (km) |
|:--|:--|:--|:--|:--|:--|
| near-field | 0.115 | 114$^\circ$ | $6.4\times10^2$ | 0.026 | 0.025 |
| random | 0.080 | 149$^\circ$ | $6.6\times10^1$ | 0.045 | 0.049 |
| ring | 0.076 | 53$^\circ$ | $1.0\times10^1$ | 0.040 | 0.045 |
| one-sided | 0.009 | 261$^\circ$ | $2.7\times10^1$ | 0.104 | 0.096 |
| clustered | 0.001 | 323$^\circ$ | $3.2\times10^0$ | 0.263 | 0.240 |
| linear | 0.000 | 180$^\circ$ | $1.7\times10^0$ | 0.775 | 0.706 |

The linear array is the worst configuration by a factor of fourteen in error, and its coverage is
exactly zero, correctly identifying the degeneracy of a network whose ray directions span a plane.
Its azimuthal gap is $180^\circ$, which looks unremarkable, and the ring has the *best* azimuthal
gap while not being the best performer. Azimuthal gap is blind to take-off-angle spread and to
distance; the direction-covariance eigenvalue is not. This is the concrete operational content of
the network-design deliverable of Section 6.1.

**Result 3: the time half of the condition.** Eight stations from 10 to 45 km, last P arrival at
13.6 s and last S at 19.9 s. Sweeping the length of the recording window:

| Window end | stations with P | with S | $\gamma$ | CRB epicentral |
|:--|:--|:--|:--|:--|
| 5.5 s | 0 | 0 | $1.6\times10^{-79}$ | not observable |
| 6.5 s | 0 | 0 | $3.0\times10^{-10}$ | not observable |
| 7.5 s | 1 | 0 | $6.6\times10^{-1}$ | 1.25 km |
| 8.5 s | 2 | 0 | $1.1\times10^{1}$ | 0.275 km |
| 11.5 s | 5 | 3 | $6.9\times10^{1}$ | 0.077 km |
| 15.0 s | 8 | 5 | $1.6\times10^{2}$ | 0.046 km |
| 20.0 s | 8 | 7 | $1.9\times10^{2}$ | 0.040 km |
| 25.0 s | 8 | 8 | $1.99\times10^{2}$ | 0.0388 km |
| 32.0 s | 8 | 8 | $1.99\times10^{2}$ | 0.0388 km |

Three regimes, exactly as the geometric control condition requires: strictly non-observable while no
ray has reached the receivers, a transition of nine orders of magnitude when the first arrival
enters the window, and saturation once the slowest arrival of the farthest station is inside, with
no gain from any longer window. The saturation time is predictable in advance from geometry alone as
$\max_s d_s/V_S+\tau$, here 19.9 s.

**Declared limitations.** Homogeneous medium with straight rays; a single source; no clutter; local
estimation started near the truth, so what is verified is local well-posedness and not global
attribution. The directional-coverage quantity is related in form to geometric dilution of
precision, so its correlation with the travel-time part of the operator is expected; what is
informative is that it predicts the coercivity of the *full waveform* operator, and that it beats
the azimuthal-gap heuristic on the same data. In heterogeneous media the quantity must be computed
by ray tracing, which is the confirmatory test and is work package T1.

**Result 4: verification of the Proposition of Section 6.1.** Code: `ci_lemma2.py`,
`ci_lemma3.py`. Over 120 non-degenerate geometries the closed-form leading term of (8b), computed
from station directions and distances alone, reproduces the numerically evaluated coercivity with
Spearman $\rho=+0.9998$. The ratio is a single geometry-independent constant, the discretization of
$\|\dot g\|^2$; after removing it the median ratio is $1.0000$ with interquartile range
$[0.992,1.022]$ and $95\%$ of geometries fall within $11\%$. The deviations are concentrated near
the degeneracy manifold, where the leading term vanishes and the subleading polarization
contribution takes over, as the asymptotic statement requires.

Corollaries 2 and 3 hold to machine precision:

| Geometry | leading term, P only | S only | P and S | interpretation |
|:--|:--|:--|:--|:--|
| ring, equidistant | $-1.1\times10^{-15}$ | $1.6\times10^{-16}$ | $5.0\times10^{-1}$ | cone degeneracy, removed by branch diversity |
| linear array | $-1.1\times10^{-13}$ | $-8.5\times10^{-14}$ | $3.9\times10^{-14}$ | directions span a plane; observable only via polarization |
| random | $1.8\times10^{0}$ | $1.2\times10^{1}$ | $1.5\times10^{1}$ | non-degenerate for every branch set |

The explicit bound (8c) held in 87 of 87 geometries where it is positive but is loose by a median
factor of 36 and has log-log slope $0.14$ against the measured constant. It is therefore reported
as a feasibility certificate and not as a predictor; the exact leading term is the predictor.

**Result 5: heterogeneous extension.** Code: `ci_hetero.py`. Medium with a constant velocity
gradient $v(z)=v_0+kz$, for which the two-point travel time is
$T=k^{-1}\mathrm{arccosh}\big(1+k^2D^2/(2v_1v_2)\big)$ with $D$ the straight-line separation.
Source-side and receiver-side slowness vectors are obtained by differentiating that closed form; the
eikonal identity $|\nabla_{\mathbf{x}_s}T|=1/V(\mathbf{x}_s)$ holds to a relative error of
$10^{-11}$ for every gradient tested, which validates the implementation. Thirty-two geometries per
gradient across four families.

| Gradient $k$ (s$^{-1}$) | $\rho$, take-off slowness | $\rho$, straight line | median deviation, take-off | straight line |
|:--|:--|:--|:--|:--|
| 0.00 | 1.0000 | 0.9993 | 0.4% | 0.4% |
| 0.03 | 0.9988 | 0.8955 | 0.8% | 43.5% |
| 0.06 | 1.0000 | 0.9010 | 0.5% | 88.6% |
| 0.10 | 0.9982 | 0.9340 | 0.6% | 41.6% |
| 0.15 | 0.9996 | 0.9608 | 0.5% | 44.5% |

The two forms coincide at zero gradient, as they must. As soon as rays bend, the take-off form
retains sub-percent accuracy while the straight-line form loses between forty and ninety per cent.
This is the quantitative separation between the control-theoretic certificate and a
dilution-of-precision heuristic: they are not the same object, and the difference is large in any
realistic Earth model.

B3 is supported in the homogeneous case, while heterogeneous validation remains the confirmatory extension.
The observability-based network-design criterion provides a measured advantage over the standard heuristic,
and the saturation result gives a practical rule for window length that can be computed before data are
processed.

---

# Appendix J. Pre-Study: Which Constraints Carry the Information

Code: `cattrib2.py`, `cattrib3.py`. One true source at 8 km depth
plus $n$ randomly polarized clutter arrivals; a perfect picker delivers every impulsive arrival,
true and false, which is what a travel-time associator receives; candidate hypotheses on a
$17\times17\times4$ grid with origin time optimized per candidate. Margin is the gap between the
true-source score and the best score among hypotheses more than 5 km away, in units of the standard
deviation of the score field. Constraint C-III was not implemented and is declared untested.

**Result 1: the ablation condition (ABL) can fail completely, and the geometry decides.**
In a ring of six stations all at 25 km, its margin is exactly $0.00$: false hypotheses reach a
perfect score of $1.000$, identical to the truth, at both clutter levels. With eight stations spread
from 6 to 62 km the same score gives margins of $2.60$ and $1.48$. Uniform-distance deployments
remove the amplitude and moveout diversity that travel-time association implicitly relies on.

**Result 2: contributions are unequal.** Change in margin when each family is removed, spread
geometry:

| Removed | $n=15$ | $n=40$ | reading |
|:--|:--|:--|:--|
| C-IV polarization | $-3.19$ | $-0.18$ | dominant at moderate clutter |
| C-V source bounds | $-0.38$ | $-0.32$ | small, consistent in sign |
| C-II silence | $-0.28$ | $+0.06$ | conditional |
| C-I observability | $+0.19$ | $+0.08$ | not a discriminator |

In the uniform ring the same ablations give C-IV $-4.67$, C-V $-0.57$, C-II $+1.27$ and C-I
$-0.73$ at $n=40$: C-IV dominates again, and C-II is actively harmful when all stations sit at the
same distance. C-IV alone reaches $6.05$ against $6.16$ for all families together at $n=15$, so on
these data almost the entire benefit is polarization.

**Result 3: constraint activation must be verified before a null is believed.** A configuration in
which C-I and C-II cannot bind - an observability floor set to zero, and an amplitude test whose
scales do not match the observed energy - reproduces the all-constraints numbers exactly and would
be read as evidence that the two families carry no information. A null obtained from an inactive
constraint is not a null. Each family was therefore shown to bind on a case constructed so that it
must, before the results were used in the confirmatory design. The same requirement is part of the Tier-1
pilot through the design-validity step of Section 13.0.

These results establish the hierarchy used in Section 6.6: C-I is a well-posedness certificate rather than
a discriminator, C-II is conditional on station-distance spread, and distance spread is therefore an explicit
experimental factor. The pre-study remains limited to a single true source, no overlapping events, a
homogeneous medium, one clutter model and two clutter levels; it identifies qualitative roles rather than
confirmatory effect sizes.

---

# Appendix K. Pre-Study: Medium-Independent Constraints, and a Negative Result on Overlap

Code: `c3_reciprocity.py`, `overlap4.py`.

**Part 1: C-III, medium-independent kinematic admissibility.** Eight stations from 6 to 62 km,
$V_{\min}=2$ km/s assumed as a bound and nothing else. Tests: same-phase Lipschitz (9a), and P-S
ordering with monotonicity of $\Delta t_{PS}$ against $t_P$.

| Case | rejected |
|:--|:--|
| true event, all picks correct | $0.0\%$ |
| true event, 1 pick replaced by clutter | $43.2\%$ |
| true event, 2 picks replaced | $75.5\%$ |
| true event, 4 picks replaced | $96.7\%$ |
| random pick set, 5 stations, 30 s window | $98.0\%$ |
| random pick set, 8 stations, 30 s window | $100.0\%$ |

Under the declared admissibility assumptions the tests reject **no uncontaminated synthetic truth**, which follows from their being theorems about that model. A true event whose pick set has been contaminated by clutter is no longer a realization of those assumptions, so the 43 to 97 per cent figures are rejections of *corrupted observation sets*, not false rejections of truth. Sensitivity is high and grows with the number of stations claimed. This makes C-III the
only family that rejects no uncontaminated truth under its assumptions, and it operates
combinatorially, before any hypothesis
exists.

**Part 2: overlap did not reproduce as a difficulty.** Two true sources plus 30 clutter arrivals,
scored on a grid with origin time optimized per candidate. Configurations: spatial separation 3, 6,
12 and 24 km, origin-time offsets 0.3, 0.5, 1, 2 and 6 s. Metric: the ambiguity set, the number of
grid hypotheses more than 4 km from both true sources scoring within 2 per cent of the worse true
source.

The ambiguity set was **empty in every configuration and for every constraint set, including
travel-time scoring alone.** Both true sources were recovered at the top of the ranking in all
cases.

Two design errors of this pre-study are recorded rather than hidden. The first attempt placed the
two sources 24 km apart and varied only the origin-time offset, so their arrivals never actually
overlapped at any station; that was corrected by bringing the sources to 3-6 km. The second and
deeper issue is that the generator supplies a **perfect picker**: every impulsive arrival becomes a
pick, whether or not it is buried in the coda of an earlier event. That removes the mechanism by
which overlap is actually reported to be hard.

C-III is therefore defined through the eikonal Lipschitz bound rather than Betti identities, which are
not testable from single-event station data. The zero rejection of uncontaminated truths is retained as a
validation property. The overlap result is treated as the registered null E3: the difficulty of overlapping
events may arise primarily at the detection stage rather than in association scoring, and confirmatory
analysis requires a realistic picker with coda failures and structured, moveout-consistent clutter.

---

# Appendix L. Pre-Study: Behaviour under Velocity-Model Error

Code: `a3_degradation.py`, `a3d.py`. Data are generated in a medium
with a constant velocity gradient and inference is performed with a homogeneous model, so the error
is structural: arrival times **and** ray directions are both wrong. A pure rescaling of velocity
would have left ray directions untouched and would have made the polarization channel appear immune
for a trivial reason; that design was rejected for this one. Model error is quoted as the RMS
relative traveltime discrepancy between the two models over the eight stations. Twelve realizations
per level; the location estimate is obtained by continuous Nelder-Mead optimization over
$(x,y,z,\tau)$ from a coarse grid start, because an earlier grid-only version was quantization
limited and could not resolve differences below one cell.

| Model error | travel time | polarization | C-III + C-IV | full stack |
|:--|:--|:--|:--|:--|
| 0% | 0.00 km | 0.00 | 0.00 | 0.00 |
| 10.9% | 0.00 | 0.00 | 0.00 | 0.00 |
| 24.1% | 7.24 | 2.71 | 6.01 | 0.00 |
| 43.0% | 6.00 | 6.07 | 6.09 | 9.32 |
| **slope, km per 10%** | **+1.70** | **+1.51** | **+1.66** | **+2.13** |

The slopes are indistinguishable and their ordering is not stable across levels. The location-accuracy
version of A3 is therefore recorded as the registered null E4.

**Why, and what it changes.** Polarization is predicted from the ray direction between candidate
and station; when the velocity model is wrong the ray direction is wrong too, so the polarization
channel inherits the model error through geometry rather than escaping it. Only C-III is
structurally immune, since it uses a velocity bound and nothing else - but C-III is a feasibility
filter, so it can reject hypotheses and cannot relocate them, and it is therefore invisible on a
location endpoint. Indeed C-III reproduced the travel-time location error exactly at every level,
which is the expected behaviour of a filter whose action leaves the argmax intact, not a failure.
A3 is therefore evaluated using a rejection endpoint rather than a location endpoint.

**Why the design-validity step exists.** In this problem the experimental design can disable the
effect under test as easily as the physics can absent it, and the resulting null is
indistinguishable from a real one without an explicit check. Four such configurations were
identified and corrected in the course of these pre-studies: a ring of equidistant stations removes
the amplitude and silence information that two constraint families depend on; a perfect picker
removes the mechanism by which overlap is difficult; a grid coarser than the effect quantizes the
differences being measured; and a constraint evaluated on model-predicted quantities rather than on
observations cannot bind at all. Each was detectable only because the number was implausible. The
Tier-1 pilot of Section 13.0 therefore includes a *design-validity* step in which every constraint
is shown to bind, and every channel is shown to be adequately sampled, on cases constructed so that
they must, before any confirmatory condition is run.

---

# Appendix M. Pre-Study: The Minimal Sufficient Set

Code: `a2_minimal.py`, `a2b.py`, `a2c.py`. Four geometry families -
spread, random, one-sided, ring - crossed with two clutter levels, six to eight realizations per
cell. One true source at 8 km depth. Pick matching and nearest-pick lookup are precomputed per
station as arrays indexed by time bin, so replication is affordable. Endpoint: the false-event
rate, the fraction of realizations in which some hypothesis more than 5 km from the truth ties or
beats the truth.

**A1, stable ordering: preliminary evidence, eight cells.** The ordering by individual contribution was
polarization $>$ moment consistency $>$ kinematic admissibility $>$ silence in every one of the
eight cells, mean pairwise rank correlation $+1.00$ against a registered requirement of $0.8$. The ordering is computed from median discrimination margins rather than
from false-event rates, which is why it remains non-degenerate even in the cells that Appendix W
shows to be saturated on the rate endpoint.

**A2, minimal sufficient set: preliminary evidence, in the cells where anything is being discriminated, against the coincidence-count ablation.** Appendix W
re-runs the spread and random cells at thirty realizations and finds all channels at $0.00$-$0.03$
with no differences: those cells are saturated and carry no information about sufficiency. The
result below should therefore be read as resting on the one-sided cells.

| Geometry | clutter | travel time | C-III | C-IV | C-V | C-II | C-III+C-IV | all |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| spread | 15 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| spread | 40 | 0.17 | 0.17 | 0.17 | 0.00 | 0.17 | 0.17 | 0.17 |
| random | 40 | 0.17 | 0.17 | 0.00 | 0.00 | 0.17 | 0.00 | 0.00 |
| **one-sided** | **15** | **1.00** | **1.00** | **0.00** | 0.00 | **1.00** | **0.00** | **0.00** |
| **one-sided** | **40** | **1.00** | **1.00** | **0.00** | 0.33 | **1.00** | **0.00** | **0.00** |
| ring | 15, 40 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |

The pair C-III plus C-IV matched the full stack within the registered $0.02$ in eight cells out of
eight. That statement needs the qualification supplied by Appendix W: re-run at thirty
realizations, the spread and random cells show **every** channel at $0.00$-$0.03$ with no
differences at all, so those cells are saturated and carry no information about sufficiency. The
informative cells are the one-sided ones, and it is there that the claim rests. More strongly, **polarization alone matched both**, so the minimal sufficient set on this
endpoint is a single channel. C-III alone reproduced travel time exactly, which is the expected
behaviour of a feasibility filter that does not displace the argmax, and C-II alone carried
nothing, consistent with Appendix J.

**The one-sided failure.** In a network confined to a half-plane, the tolerance-matching
counterfactual produced a false event that tied or beat the truth in **every realization**, at both
clutter levels. That is
failure, not degradation. Polarization removed it entirely. One-sided geometries describe coastal
networks, networks truncated by a border, and single-flank ocean-bottom deployments, so this is not
a pathological construction.

**Stress test 1: coda-like clutter.** Clutter polarization was scattered around the parent event's
particle-motion direction, from $0^\circ$ (perfect mimicry) to $90^\circ$ (random). Polarization
held at $0.00$-$0.12$ throughout. This test is reported with a caveat that weakens it: clutter
aligned with the *true* source reinforces the truth rather than supporting a false hypothesis, so
it is not the adversarial case it was intended to be.

**Stress test 2: adversarial decoy.** Clutter was instead radiated by a decoy source at another
location, with fully coherent polarization *and* coherent moveout, so that a hypothesis at the decoy
is supported by consistent physics. Support was swept from 4 to 16 of the 16 possible arrivals.

| Geometry | decoy arrivals | travel time | C-IV | C-III+C-IV | all |
|:--|:--|:--|:--|:--|:--|
| spread | 8/16 | 0.00 | 0.00 | 0.00 | 0.00 |
| spread | 12/16 | 0.00 | 0.00 | 0.00 | 0.12 |
| spread | 16/16 | 1.00 | 0.12 | 0.12 | 0.62 |
| one-sided | 4/16 | 1.00 | 0.12 | 0.12 | 0.12 |
| one-sided | 12/16 | 1.00 | 0.12 | 0.12 | 0.12 |
| one-sided | 16/16 | 1.00 | 0.50 | 0.50 | 0.38 |

Polarization survives partial decoy support and degrades only when the decoy becomes a complete
second event, at which point outranking the first source is not an error. Note also that the full
stack was *worse* than polarization alone in two adversarial cells, which is direct evidence for
A2: adding constraints is not monotonically beneficial.

**Declared limitations.** A single true source; a homogeneous medium with an exact velocity model,
which Appendix L shows is the condition under which polarization is most favoured; perfect picks;
one clutter amplitude; six to eight realizations per cell. The confirmatory version requires
overlapping sources, model error, a realistic picker, and enough replicates for effect sizes rather
than rates.

**Note added after Appendix Y.** The travel-time column of this appendix is the coincidence-count ablation. Against a continuous kinematic score with exact picks, the one-sided advantage attributed here to polarization disappears (Appendix Y.1).

---

# Appendix N. Pre-Study: Cost, and Network Design

Code: `a4_b2.py`, `b2b.py`.

**A4, cost: numerical evidence.** Measured wall-clock per candidate evaluation, mean of four thousand
calls, same machine and same implementation:

| Channel set | microseconds per evaluation | relative to travel time |
|:--|:--|:--|
| travel time | 17.4 | 1.00 |
| C-III kinematic admissibility | 48.7 | 2.80 |
| C-V moment consistency | 39.4 | 2.26 |
| C-IV polarization | 100.2 | 5.76 |
| C-II silence | 114.3 | 6.57 |
| C-III + C-IV | 175.5 | 10.08 |
| all four | 301.6 | 17.33 |

Against the reference costs declared in Section 12 - $38{,}000$ microseconds per objective
evaluation for analytic beam synthesis with a transport misfit, and $5.76\times10^8$ microseconds
for a direct elastodynamic solve - the minimal sufficient set costs $0.46\%$ of the waveform
condition and $0.00003\%$ of the direct solve. Polarization alone costs $0.26\%$. Since Appendix M
established that the registered pair matches the full constraint stack on the rejection endpoint in the tested discriminating cells, the registered
A4 margin of ninety per cent of the benefit at under one per cent of the cost is met with three to
four orders of magnitude to spare.

**B2, network design: inconclusive on the registered endpoint, and the endpoint was wrong.**
Three thousand random networks per station count; the best by maximum observability constant, the
best by minimum azimuthal gap, and a median-observability control were each evaluated over ten
realizations with continuous Nelder-Mead location.

| Stations | design | $\gamma$ | azimuthal gap | median error | false-event rate |
|:--|:--|:--|:--|:--|:--|
| 4 | max observability | $2.3\times10^{-4}$ | $149^\circ$ | 0.00 km | **0.00** |
| 4 | min azimuthal gap | $1.9\times10^{-6}$ | $94^\circ$ | 0.00 km | **0.20** |
| 4 | median-$\gamma$ control | $7.0\times10^{-6}$ | $153^\circ$ | 0.00 km | 0.20 |
| 5 | max observability | $2.7\times10^{-4}$ | $148^\circ$ | 0.00 km | 0.00 |
| 5 | min azimuthal gap | $3.3\times10^{-5}$ | $82^\circ$ | 0.00 km | 0.00 |
| 6 | all three | - | - | 0.00 km | 0.00 |

The design criterion works as intended - the observability-optimal networks have a constant two
orders of magnitude larger than the azimuthal-gap-optimal ones - but median epicentral error
saturates at zero for every design, so the registered endpoint cannot separate them. The signal
that does appear is on false-event rate at four stations, and ten realizations is too few to
register it. B2 is retained with the endpoint changed from location to rejection.

**The pattern, now explicit.** This is the third pre-study in which a resolution-type endpoint
saturated while a rejection-type endpoint discriminated: Appendix H found the wavefield adds
$6$-$14\%$ to resolution but four times the rejection margin; Appendix L found no channel retains
location accuracy under model error while C-III retains rejection validity exactly; and Appendix N
finds network design invisible in location error and visible in false-event rate. The expectation
is registered for the confirmatory programme: **endpoints for this problem should be rejection
endpoints, and resolution endpoints should be reported but not used to discriminate methods.**

---

# Appendix O. Pre-Study: The Sufficiency Ladder

Code: `f1_c1.py`, `f1b.py`. Random networks of three to twelve
stations, thirty-five clutter arrivals, SNR 6, eight realizations per cell. Endpoint: false-event
rate. The waveform channel is a genuine matched filter - normalized correlation between the
observed three-component record and the fully predicted waveform, both branches, correct
polarization, $1/d$ amplitude and correct wavelet - precomputed by convolving each component with
the wavelet so that evaluating a candidate is $O(1)$ per arrival. It uses no picks at all.

| Stations | travel time | polarization | C-III + C-IV | matched filter |
|:--|:--|:--|:--|:--|
| 3 | 1.00 | 0.75 | 0.75 | **0.00** |
| 4 | 0.88 | 0.00 | 0.00 | 0.00 |
| 6 | 0.00 | 0.00 | 0.00 | 0.00 |
| 12 | 0.00 | 0.00 | 0.00 | 0.00 |

**The ladder.** At six stations and above, travel time alone suffices and every additional channel
is redundant. At four stations the coincidence-count ablation fails in most realizations and the
polarization constraints are required against it (five stations were not tested; Appendix Y.1
shows that a continuous kinematic score does not fail at four stations with exact picks). At three stations the constraint stack itself fails three times in four,
and only the full waveform condition succeeds. The amount of physics required is therefore a
function of network sparsity, which answers the primary question of Section 4.1 with a graded
answer rather than a single minimal set.

**C1: crossover measured against the ablation.** The advantage of the constraint stack over the
coincidence count is $+0.88$ at four stations, $+0.25$ at three, and zero at six and twelve, so
$S^\ast\approx5$-$6$ against that comparator.

**F1: the gate opens, narrowly.** The waveform condition beats the minimal set by $0.75$ at three
stations, well beyond the registered $0.03$, and by exactly zero everywhere else. The waveform
branch of the programme is retained and **scoped to the sparse regime**, which is a much smaller
commitment than the initial programme design and a much better justified one.

**Sampling adequacy governs this result and is registered as a condition of the test.** A
matched-filter correlation peak is about half a period wide. Sampling origin time in quarter-period
steps places the true source between samples while an occasional false candidate lands on peak, and
under that sampling the matched filter appears catastrophically bad, with false-event rates of
$0.80$ to $1.00$ at every station count - the opposite of the table above, and enough to close the
waveform branch on an artefact. The figures reported here sample origin time at a tenth of a period.
For any channel whose score varies on the scale of the wavelet, the design-validity step of
Section 13.0 must include an explicit sampling-adequacy check.

**Declared limitations.** A single true source, homogeneous medium with an exact velocity model,
one clutter model and amplitude, eight realizations per cell, and one SNR. The three-station cell
rests on six failures out of eight against zero out of eight, which is suggestive rather than
established.

**Note added after Appendix Y.** The travel-time column is the coincidence-count ablation. Against a continuous kinematic score with exact picks, the four-station advantage attributed here to polarization disappears, and at three stations the kinematic kernel emits fewer false events than the constraint stack (Appendix Y.1).

---

# Appendix P. Pre-Study: One Station, Several Events

Code: `c3_single.py`. A single three-component station; one, two or
three sources at 8 km depth with different epicentres and origin times; twelve clutter arrivals;
SNR 8. Sources are extracted by **hierarchical subtraction** - take the strongest matched-filter
peak, subtract its best-fit contribution, repeat - which is the strategy developed for overlapping
gravitational-wave signals and is imported here rather than reinvented. Depth is fixed at the truth,
because Section 4.3 declares depth ill-conditioned from one station and claims only epicentral
position, origin time and source count.

| Sources | branches | source 1 | source 2 | source 3 | spurious |
|:--|:--|:--|:--|:--|:--|
| 1 | P only | 22.02 km | - | - | 2 |
| 1 | **P and S** | **1.58 km** | - | - | 1 |
| 2 | P only | 4.74 km | 22.55 km | - | 0 |
| 2 | **P and S** | **1.58 km** | **1.00 km** | - | 0 |
| 3 | P only | 28.46 km | 6.00 km | 26.17 km | 1 |
| 3 | P and S | 5.70 km | 10.59 km | 22.39 km | 4 |

**The structural prediction holds.** With one branch, polarization fixes the back-azimuth but $t_P$
cannot fix the distance while the origin time is unknown, so the score field is a ridge along a ray
and the recovered position slides along it - hence errors of twenty kilometres and more. With two
branches, $t_S-t_P$ fixes the distance independently of the origin time and the ridge collapses to
a point. This is Corollary 3 of Section 6.1 in operational form.

**The capability boundary is two.** One station resolves one or two overlapping events to about a
kilometre. Three defeats it, for both representations. The count is consistent with the observable
budget: one station supplies roughly three usable quantities per event - two arrival times and a
back-azimuth - against three epicentral unknowns per event once depth is fixed, so two events sit at
the determinacy limit and three exceed what hierarchical subtraction can sustain once the residuals
of the first two contaminate the third.

**Where this sits in the ladder.** Appendix O established that three stations require the waveform
condition and six make everything redundant. Appendix P extends the ladder to its lower end: one
station, two events, waveform condition and two branches mandatory. The ladder is now specified
from $S=1$ to $S=12$.

**Declared limitations.** Depth fixed at the truth, homogeneous medium with an exact velocity model,
a single realization per configuration, one clutter level, and a known source mechanism. Depth in
particular is a substantial assist and the result must not be quoted without it.

---

# Appendix Q. Pre-Study: The Certificate under Interfaces

Code: `b1_interfaces.py`. Two-layer model, layer velocity $V_1$ to
depth $H=14$ km over a halfspace with $V_2>V_1$, source at 7 km depth. Each phase arrives twice:
the direct ray, and the critically refracted head wave with
$t=x/V_2+(2H-z)\cos\theta_c/V_1$ and $\sin\theta_c=V_1/V_2$. Forty geometries across four station
layouts.

**Eikonal validation.** The source-side slowness magnitude equals $1/V_1$ to machine precision for
both branches and both phases, including the head wave, whose take-off direction sits at the
critical angle irrespective of epicentral distance. This is the analytic check that the branch
construction is right before any statistics are computed.

| Predictor | Spearman against measured $\gamma$ |
|:--|:--|
| **leading order, summed over both branches** | **$+0.963$** |
| leading order, direct branch only | $+0.871$ |

B1 registered $\rho\ge0.95$ and is met. Ignoring multipathing costs about $0.09$ in correlation, so
the certificate must be assembled as a sum over ray branches, exactly as the extension predicted.
The head wave raises the median leading-order constant by a factor of $2.73$, which has a clean
interpretation through Corollary 2: head-wave take-off directions all lie on a cone of fixed
aperture and would be degenerate on their own, but combined with the spread of direct take-off
directions they add rank.

**Loss of precision.** Ninety-five per cent of geometries fall within $59\%$ of the closed-form
prediction, against $11\%$ for the homogeneous and smooth-gradient cases of Appendix I. The
certificate therefore remains a good *ordering* statistic under interfaces but a much weaker
*quantitative* one. The head-wave amplitude law used here is crude, so part of the loss is the
test's and not the theory's; separating them requires a proper geometrical-spreading treatment and
is registered as remaining T1 work.

**Declared limitations.** Two layers, one interface, no triplications from a gradient over the
interface, no converted phases, no attenuation, and an approximate head-wave amplitude. A layered
model with several interfaces and genuine triplications is the confirmatory case.

---

# Appendix R. Pre-Study: Detection and Association are One Axis

Code: `d1_coupling.py`. Six stations from 8 to 50 km, thirty clutter
arrivals, SNR 6, twelve realizations per level. Each true arrival is deleted from the pick list with
probability $m$, while clutter picks remain, since a picker that misses buried arrivals still
triggers on impulsive noise. The waveform is always present in the record; only the *pick* is lost.
The waveform channel uses no picks at all. Endpoint: recovery rate, the fraction of realizations
whose best hypothesis lies within 5 km of the truth, which captures both failure modes at once -
being unable to form the event, and forming a false one.

| Picker miss rate | travel time | polarization | C-III + C-IV | pick-free waveform | advantage | supported? |
|:--|:--|:--|:--|:--|:--|:--|
| 0.0 | 1.00 | 1.00 | 1.00 | 1.00 | $+0.00$ | - |
| 0.2 | 0.92 | 0.92 | 0.92 | 1.00 | $+0.08$ | **no**, $p=1.00$ at $n=30$ |
| 0.3 | 0.92 | 0.83 | 0.83 | 1.00 | $+0.08$ | **no**, $p=0.11$ at $n=30$ |
| 0.5 | - | - | - | - | $+0.60$ | **yes**, $p<10^{-4}$ at $n=30$ |
| 0.6 | 0.42 | 0.75 | 0.75 | 1.00 | $+0.25$ | yes |
| 0.8 | **0.00** | 0.33 | 0.33 | **1.00** | $+0.67$ | yes, $p<10^{-4}$ at $n=30$ |

Rows at $0.0$ to $0.8$ are at twelve realizations; the supported/not-supported column and the $0.5$
row come from the thirty-realization re-run of Appendix W. **The crossover lies in $(0.3,0.5]$**,
which is the threshold D1 registers as its deliverable.

The advantage of the pick-free channel is monotone in the miss rate, and is statistically supported
only from a miss rate of about $0.6$ upward; the $+0.08$ differences at $0.2$ and $0.3$ are within
noise at twelve realizations (Appendix W). Above that threshold the effect is unambiguous, as D1
predicts,
and the pick-free channel is flat at $1.00$ throughout, which is the expected behaviour of a method
that never consults a pick. Polarization degrades more slowly than travel time, recovering a third
of cases at a miss rate of $0.8$ where travel time recovers none, because it extracts more
information from each surviving pick.

**The synthesis.** Appendix O measured how much physics is needed as a function of station count;
this appendix measures it as a function of picker quality; and the two are the same function of the
**effective number of usable observations**. Six stations at a miss rate of $0.8$ leave roughly one
station's worth of picks, which is exactly where Appendices O and P placed the waveform condition
as mandatory. The programme therefore has one axis to characterize, not two.

**Declared limitations.** Miss probability is uniform and independent across arrivals, whereas real
pickers miss in a correlated way - typically inside the coda of an earlier arrival, which is the
mechanism Appendix K failed to isolate. Twelve realizations per level; one network; one clutter
level; one SNR. The confirmatory version should drive the miss pattern from a calibrated picker so
that the correlation structure of the misses is realistic.

---

# Appendix S. Pre-Study: Observability Maps Where Discrimination is Hard

Code: `b3_local.py`, `b3b.py`. One fixed, deliberately asymmetric
six-station network; the true source is moved over twenty-five locations spanning the region;
thirty clutter arrivals; three realizations per location. For each location the leading-order
observability constant is computed in closed form from geometry alone, and the local
false-hypothesis margin is measured from the data.

| Quantity | Range across the region |
|:--|:--|
| observability constant $\gamma$ | $1.6\times10^{-6}$ to $1.4\times10^{-4}$, a factor of 86 |
| local discrimination margin | $-0.13$ to $5.12$ |

| Quartile by $\gamma$ | median margin | false-event rate |
|:--|:--|:--|
| lowest | 0.93 | 0.06 |
| highest | 4.24 | 0.00 |

Spearman correlation between the constant and the local margin: $+0.573$, against the registered
requirement of $0.5$. B3 is met.

**What this gives operationally.** The constant is a closed-form function of station geometry and
velocities and requires no data. Computing it over a monitoring region therefore produces, in
advance of any recording, a map of where false events will be easy to form. That is a different and
more useful role than the role assigned to C-I in the constraint hierarchy, and it is the role Appendix J showed it
could not fill.

**Endpoint note.** The correlation with the local false-event *rate* was only $+0.147$, because that
endpoint had almost no variance in this regime: the rate was zero in three of four quartiles. This
inverts the pattern of Appendix N, where the rate discriminated and the resolution endpoint
saturated. The registered rule that follows is to select the endpoint carrying variance in the
regime under test, and to report which endpoint saturated rather than quoting a null from a
saturated one.

**Declared limitations.** One network, three realizations per location, a single depth, one clutter
level, homogeneous medium with an exact velocity model. Twenty-five locations is enough for a rank
correlation but not for a calibrated map.

---

# Appendix T. Pre-Study: Coherent Nearby Clutter

Code: `e1_critical.py`. The clutter model used elsewhere in this document draws polarization at
random, which is the configuration most favourable to the polarization channel. The dangerous case
is different: clutter coherent with a **nearby decoy source**, as produced by an aftershock or by
the coda of a neighbouring event, whose arrivals carry consistent polarization *and* consistent
moveout pointing at a location a few kilometres from the truth. Eight stations from 8 to 50 km, one
true source, a decoy at controlled separation with $6$ or $12$ of its $16$ possible arrivals
present, eight random arrivals, SNR 7, six realizations per cell. Outcome classes: the best
hypothesis lies within 3 km of the truth, within 3 km of the decoy, or neither, the last being a
merged or spurious solution.

At **thirty realizations**, which is the count these figures should be read from:

| Decoy separation | support | travel time spurious | polarization spurious |
|:--|:--|:--|:--|
| 4 km | 6/16 | $30/30=1.00\ [0.88,1.00]$ | $12/30=0.40\ [0.23,0.59]$ |
| 12 km | 12/16 | $30/30=1.00\ [0.88,1.00]$ | $12/30=0.40\ [0.23,0.59]$ |

Fisher exact $p<10^{-4}$ at both separations. An earlier pass at six realizations suggested that the
polarization spurious rate varied with separation, from $0.00$ at 4 km to $0.50$ at 12 km; at
thirty realizations it is constant at $0.40$ and that variation was sampling noise.

**The ablation condition (ABL) fails completely against a coherent nearby decoy.** In 83 to 100 per cent
of realizations its best hypothesis is neither the true source nor the decoy but a merged artefact
between them, at every separation tested. This is a failure mode distinct from the geometric one of
Appendix M and it occurs in the configuration that dominates the first hours of an aftershock
sequence.

**The polarization channel is not immune, and the degradation is real.** At thirty realizations
(Appendix W) its spurious rate against a coherent decoy is $0.40\ [0.23,0.59]$, **independent of
decoy separation**, against $1.00$ for travel-time association in every realization. The variation
with separation reported at six realizations was noise; the degradation from approximately $1.00$
correct recovery under randomly polarized clutter is not. The relative
advantage over travel time survives and is large, but the absolute performance reported elsewhere
in this document is optimistic to the extent that real clutter is coherent rather than random.

**What this makes measurable.** The quantity that decides how much of the polarization advantage
survives in practice is the **polarization coherence of real false picks**, which is measurable
from public continuous archives and, as far as the search of Appendix E extends, has not been
characterized. The table above is the response curve against which such a measurement is to be
compared, and obtaining it is registered as a first-quarter deliverable of WP E1. It has now been
obtained for the Colombian network (Section T.1).

**Declared limitations.** Six realizations per cell, one decoy at a time, one origin-time offset,
homogeneous medium with an exact velocity model, and perfect picks. The travel-time comparator is
the coincidence-count ablation, whose merged-artefact failure has not been re-measured against a
continuous kinematic score; under the protocol of Appendix X the decoy cells show at most one
false event in thirty for any method (X.2), so the comparison is registered for repetition against the kinematic kernel of
Appendix Y.

## T.1 Measurement on real data: polarization coherence of real false picks

Code: `task1_fetch.py`, `task1_analyze.py`, `task1_decide.py`, `diagnose_polarization_config.py`.
The measurement registered above was executed on the Colombian network CM, with continuous
three-component HH data at 100 Hz from the SGC FDSN service and the ISC reviewed catalogue.
PhaseNet (SeisBench, `instance` weights, threshold $0.30$) was run on windows around each catalogue
event within 200 km of the station and on consecutive 300 s windows of the coda of the same event.
A pick is **associated** if it falls in $[t_P-5\,\mathrm{s},\,t_S+5\,\mathrm{s}]$ of some catalogue
event at that station; the remaining picks are the **clutter**. For each clutter pick the reference
event is the most recent catalogue event within 200 km whose predicted P arrival precedes the pick.
Polarization is the principal eigenvector of the three-component covariance, folded to $[0,180)$
because its sign is arbitrary, and a pick is coherent with an event when the folded angle to the
station-to-event back-azimuth is at most $20^\circ$.

**Validity gate.** A clutter figure is interpretable only where the estimator recovers the direction
of events whose location is known. The gate is therefore run first, on associated P picks, and
requires a coherent fraction of at least $0.60$ on at least 30 picks. S picks are excluded from the
gate because SH particle motion is transverse to the back-azimuth.

**Chance correction.** Clutter with no relation to its reference event still falls within $20^\circ$
of it at the chance rate $f_0$: $20/90=0.22$ for uniform azimuths, and $0.24$ to $0.33$ in every set
of at least 30 picks measured here, where events cluster in azimuth. The registered threshold of
$0.20$ on the raw fraction $f$ therefore lies below the value expected from entirely incoherent
clutter, and a pass or fail on it measures the geometry rather than the coherence. The registered
criterion is applied instead to
$$\eta=\frac{f-f_0}{f_{\mathrm{ctrl}}-f_0},$$
where $f_0$ is the exact expectation of $f$ when reference azimuths are permuted among picks of the
same station (the mean of the pick-by-event coherence matrix), and $f_{\mathrm{ctrl}}$ is the hit
rate of the gate. Under a mixture model in which a coherent clutter pick behaves like a real arrival
and an incoherent one like a random direction, $\eta$ is the coherent share: $0$ for unrelated
clutter and $1$ for clutter indistinguishable from real arrivals, whatever the network geometry.
The 95% interval is a station-stratified bootstrap over clutter and gate picks. On synthetic
mixtures the analytic $f_0$ matches a 20,000-permutation Monte Carlo to $10^{-4}$; $\eta$ recovers
planted shares of $0$, $0.3$ and $0.7$ within its interval in single replicates ($0.098$, $0.223$
and $0.776$); the interval width is $0.55$ at 50 clutter picks and $0.23$ at 300; and with no
coherence the verdict is BAJA in 5 of 20 repetitions at 50 picks and 20 of 20 at 300 or more. The
bias of $\eta$ and the coverage of its interval have not been characterized and are registered
calibration work.

**Table T.1. Gate and clutter result per station.**

| Station | window | gate: P picks | gate fraction | chance $f_0$ | gate status | clutter |
|:--|:--|:--|:--|:--|:--|:--|
| ZAR | 2024-10-01, 90 d | 109 | $0.789$ | - | **passes** | $\eta=0.054\ [-0.049,0.164]$, $n=159$ |
| CHI | 2024-10-01, 90 d | 78 | $0.256$ | $0.256$ | fails; best of 24 settings $0.385$ | not interpretable |
| RUS | 2025-02-23, 30 d | 45 | $0.244$ | - | fails; P and S pooled, earlier protocol | not interpretable |
| OCA | 2024-10-01, 90 d | 33 | $0.303$ | $0.237$ | fails; not distinguishable from chance, $p=0.23$ | not interpretable |
| SJC | 2024-07-03, 180 d | 19 | $0.474$ | $0.273$ | inconclusive: below 30 picks; $p=0.044$ against chance, uncorrected | not interpretable |

ZAR passes under the setting 1-3 Hz, 2.0 s window centred on the pick, selected on gate picks only
and before $\eta$ was computed under it. The setting gives $0.787$ and $0.723$ on two disjoint halves
of the events; the halves come from the sample used for the selection, so this is a stability check
and not held-out validation. In the diagnostic run that selected it (118 to 136 gate picks), seven
of the eight 1-3 Hz settings reach $0.675$ to $0.754$, so the pass does not depend on one lucky
setting, and under the default setting (1-3 Hz, 0.6 s) the gate fraction is already $0.565$ against
a chance rate of $0.324$ ($p=0.0002$).

**Verdict.** At ZAR, $\eta=0.054$ with interval $[-0.049,0.164]$, entirely below $0.20$. By the
registered criterion, Appendices H, J and M keep their quantitative validity where the polarization
channel is usable, in the following sense: the measured false picks are not coherent with the
preceding catalogued event within 200 km, which is what the pre-studies assume of clutter. Two
things are not established by this measurement. Coherence with an **uncatalogued** source, such as
a small event below the completeness of the ISC catalogue, cannot be detected by $\eta$ as defined,
because its reference is always a catalogued event, and that is precisely the physically realizable
impostor of the table above. And the measurement is for background seismicity; an active sequence
is registered separately (Section 16).

**The finding that was not registered.** At four of the five stations the gate does not pass: CHI
and OCA fail it under the final protocol, RUS failed it under the earlier one, and SJC has too few
picks to decide. At CHI single-station P polarization is at the chance rate, and the alternatives
were excluded one at a time: the inventory declares standard orientation in a single epoch since
2019; no fixed rotation ($0.321$) and no fixed reflection, which is what a swapped or reversed
horizontal channel produces ($0.282$), reaches $0.60$; the gate fraction is stable across four
chronological quarters of the period ($0.25$ to $0.26$), which excludes a change of instrument
during it; restricting to P picks within 1 s of the predicted arrival raises it only to $0.348$ (46
picks), so mis-association does not explain it; and the same code passes at ZAR and on synthetic
data, with a median error of $3^\circ$ at the highest noise tested. Near-vertical arrivals are the
worst ($0.136$ below $20^\circ$ incidence), and CHI's S picks score $0.677$ on 65 picks, far above
its P picks. The two observations have a common candidate explanation that the tests above do not
exclude: at steep incidence the horizontal component of P is small and dominated by noise and
near-surface effects [95], whereas SV arrives close to horizontal and in the radial plane, so a
station that records mostly steep arrivals, for instance from the intermediate-depth seismicity
of the Bucaramanga nest [107], would show
P at the chance rate and S aligned with the back-azimuth. The remaining candidates are the
installation and the local structure beneath the station, which the tests do not separate. The
gate restricted to incidence of $30^\circ$ or more, with the event-depth breakdown, decides between
geometry and site and is registered (Section 16). OCA fails the gate
but cannot be distinguished from chance, so whether it carries any directional information is open.
At SJC an apparent rotation of $165^\circ$ found at 15 picks fell from $0.667$ to $0.579$ when the
sample grew to 19, as expected of a selection artefact over 72 candidate transformations. The
polarization channel that Section 6.6 identifies as the primary discriminator therefore has, on
real data, a **station-dependent** input, and its reliability must be measured per station before
it is enabled rather than assumed.

**Deviations from the registered protocol.** All were fixed before the result they affect was
computed, except where stated.

| Registered | Executed | Reason |
|:--|:--|:--|
| Gate on all associated picks | P picks only | S particle motion is not constrained to the back-azimuth (SH is transverse), so S picks do not test the estimator |
| Criterion on raw fraction $f$ | Criterion on $\eta$ | Incoherent clutter has expected raw fraction $f_0=0.22$ to $0.33$, above the threshold |
| 0.6 s window, 2 Hz band | Setting chosen per station on gate picks, validated on halves | Default window is half pre-arrival noise and under one cycle at 1-2 Hz |
| Reference: nearest catalogue event in origin time | Most recent event within 200 km preceding the pick | With one event every 28 min, the nearest in time was often distant or later than the pick |
| Verdict bands on the interval | Added an inconclusive label where the interval crosses one boundary | A crossing interval is not evidence of intermediate coherence |
| Fixed sample | Two looks at ZAR: 10 coda windows per event ($\eta=0.044\ [-0.094,0.200]$, inconclusive), then 20, declared final before running | The first interval ended exactly on the boundary; no correction for the interim look was applied |

The raw clutter fraction at ZAR under the default setting ($0.200$ against a chance rate of $0.263$)
was printed by the diagnostic before the criterion was moved to $\eta$; it is disclosed here because
it was seen, although it was not computed under the setting finally used.

**Declared limitations.** One passing station, one 90-day period, reference events at 100 to
200 km, PhaseNet at threshold $0.30$, clutter drawn from coda windows up to about 1.75 h after
the S arrival, labels from the ISC catalogue, whose completeness in Colombia is lower than that of
the SGC catalogue, and no orientation audit of the stations beyond the inventory metadata and the
fixed-rotation sweep. Every band tested lies between 1 and 10 Hz; the primary band of Section 1.5 is
$0.1$-$2$ Hz, and whether single-station polarization is reliable there at the failing stations has
not been measured. The response values of the table above were obtained with the 0.6 s window and
have not been re-measured with the setting used at ZAR.

## T.2 Consequence for the constrained prototype: measured station reliability

Code: `task1e_polarization_reliability.py`. Section T.1 establishes that on real data a station
either keeps the direction of true arrivals or does not. Every pre-study, and Appendix X, assumes
the first at every station. This section repeats the final protocol of Appendix X (network drawn
per realization, threshold relative to the realization maximum, calibration at recall $0.98$,
thirty calibration and thirty evaluation realizations per cell, paired exact McNemar) and changes
one thing: the polarization of the **true** arrivals at each station, according to its profile.

**Profiles.** *Good*: azimuth kept, with a Gaussian error of $10^\circ$, the scatter measured at
ZAR. *Bad*: azimuth drawn uniformly, incidence kept. *Bad-3D*: direction drawn uniformly on the
sphere, for stations whose incidence is also noise, as the measured incidences at CHI suggest.
Good stations are assigned at random in each realization. Clutter, noise and networks are drawn
from a separate generator, so all conditions see the same realizations and the comparisons between
conditions are paired as well. With every station ideal, the generator and the score fields are
identical, bit for bit, to those of the Task 0 code on the same realizations, and two runs with the
same seed agree in all 84 cell-method combinations (`results/task1e/verification.txt`). The
resulting rates agree with Appendix X in the one-sided cell (ablation 10 of 30 against 11 to 13,
constrained prototype 0 in both) but only within sampling variation in the ladder cells (three
stations: constrained recall 25 of 30 against 29), because the Appendix X runs cannot be regenerated
realization by realization from the archived code.

**Four treatments of C-IV.** PHYS applies it at every station, as the pre-studies do. PHYS-G
switches it off at stations that fail the gate of T.1. PHYS-N replaces the measured cosine of every
arrival at a failing station by a **neutral** value: the mean cosine of true arrivals at passing
stations, measured on the calibration realizations ($0.71$ to $0.81$), which is what the gate of T.1
measures on catalogued events. TOL is the ablation. When no station passes, PHYS-N reduces to the
kinematic filter C-III and cannot be worse than it.

**Table T.2. False-event rate, count out of 30.** Seed 2024; in brackets, an independent seed.

| Good stations | one-sided, clutter 40: TOL | PHYS | PHYS-G | PHYS-N | four stations: TOL | PHYS | PHYS-G | PHYS-N |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| all, ideal | 10 | **0** | 0 | 0 | 25 | **10** | 10 | 10 |
| all, $10^\circ$ scatter | 10 | 1 | 1 | 1 | 25 | 9 | 9 | 9 |
| one half | 10 [13] | 6 [1] | 11 [5] | **1 [0]** | 25 [28] | 18 [19] | 23 [19] | **15** [24] |
| one quarter | 10 [13] | 10 [0] | 20 [19] | **2 [3]** | 25 [28] | 23 [25] | 26 [25] | **14 [19]** |
| one quarter, bad-3D | 10 [13] | 14 [11] | 22 [17] | **2 [2]** | 25 [28] | 30 [29] | 27 [26] | **15 [20]** |
| none | 10 | 10 | 10 | 10 | 25 | 27 | 25 | 25 |
| none, bad-3D | 10 [13] | 17 [27] | 10 [15] | 10 [15] | 25 [28] | 27 [29] | 25 [28] | 25 [28] |

**What the table establishes.**

*The scatter of a good station costs almost nothing.* With every station at the scatter measured
at ZAR, the constrained prototype keeps its advantage: 1 false event in 30 against 10 for the
ablation in the one-sided cell ($p=0.012$), 9 against 25 with four stations ($p=0.0004$).

*Stations whose direction is noise make naive C-IV harmful.* With no good station and bad-3D
profiles, naive C-IV loses recall with four stations, 14 and 18 of 30 against 30 and 28 for the
ablation ($p<10^{-4}$ and $p=0.006$), and produces more false events in the one-sided cell, 17 and 27
against 10 and 13, in the same direction in both seeds but significant in one ($p=0.12$ and
$p=0.001$). When bad stations keep their incidence the loss of naive C-IV depends on the seed (10
and 0 false events at one quarter good), so that part is not established.

*Switching C-IV off where it fails does not repair it.* At one quarter good in the one-sided cell
PHYS-G gives 20 and 19 false events against 10 and 13 for the ablation ($p=0.031$ and $p=0.15$),
and across the mixed cells its direction against the ablation is inconsistent. The mechanism is
specific: a hypothesis explained only by arrivals at failing stations receives a polarization factor
of one, perfect agreement, while the true source is penalized by the imperfect cosines of its
arrivals at passing stations. Missing evidence is treated as agreement.

*Neutral treatment recovers the advantage over the ablation with one station in four good, close
to the one in five measured in T.1.* The ablation here is TOL, whose false events in these cells
are largely tie effects (Appendix Y.1), so the numbers below measure the ordering of the treatments
of failing stations, not the value of polarization over the best kinematic score. PHYS-N gives 2 and 3 false events against 10 and 13 for the ablation in the
one-sided cell ($p=0.039$ and $p=0.006$), 14 and 19 against 25 and 28 with four stations
($p=0.003$ and $p=0.012$), and the same with bad-3D stations ($p\le0.039$ in all four cells). Its
recall stays between 28 and 30 of 30 in every mixed cell of these two scenarios. With one station in two good, the
four-station cell replicates in direction but not in significance (15 against 25, $p=0.02$; 24
against 28, $p=0.29$).

**Declared status.** PHYS-N was designed after PHYS-G failed on the realizations of seed 2024, so
its result on that seed is exploratory. The independent seed was run afterwards on the mixed
conditions only, with PHYS-N fixed, and is the evidence for it: five of six mixed cells significant
against the ablation, the sixth in the same direction.

**What this changes.** C-IV is enabled per station with a measured reliability, and a station that
fails the gate contributes **neutral** polarization evidence rather than none. That is a rule for
the constrained model, not a tuning: it follows from treating a failing station as uninformative.
The calibrated per-station likelihood that makes the rule exact is registered work
(`tasks_and_real_data_plan.md`, Task 1h).

**Three stations.** Every method has 24 to 30 false events of 30 at three stations, as in the
three-station runs behind Appendix X (0.60 to 1.00 for all four methods), so the cell does not
discriminate on false events.
On recall, PHYS-N at one quarter good recovers 27 of 30 against 17 for the ablation and 23 for
naive C-IV.

**Declared limitations.** The same single synthetic generator as every pre-study; a two-class
station model where real reliability is continuous; the synthetic polarization term uses the full
three-dimensional direction, and even with ideal polarization the synthetic estimator places only
$0.44$ of true P arrivals within $20^\circ$ in horizontal azimuth, so the synthetic good station is
noisier in azimuth than ZAR and the correspondence with T.1 is qualitative; thirty realizations
per cell and two seeds.

---

# Appendix U. Pre-Study: Association Without a Velocity Model

Code: `a_i_b.py`. If the polarization channel carries the discrimination and the kinematic
admissibility of C-III needs only a velocity bound, then an associator can be built that uses **no
velocity model at all**: intersect the *measured* back-azimuths rather than predict them, and
impose the eikonal Lipschitz test on the claimed arrival times. Scoring counts the arrivals whose
measured azimuth points at the candidate within a tolerance, weights by station coverage, and
requires the consistent set to be admissible. Data are generated in a gradient medium and inference
is performed under a homogeneous model, so that the comparison is against methods suffering
structural model error.

| Model error | travel time | polarization (predicted) | azimuth intersection (no model) |
|:--|:--|:--|:--|
| 0.0% | 0.00 km | 0.00 km | 11.3 km |
| 10.4% | 2.50 | 0.00 | 21.8 |
| 22.3% | 3.54 | 2.50 | 14.1 |
| 39.4% | 5.59 | 6.33 | 16.4 |

**The model-independence property holds and the accuracy does not.** The azimuth-intersection
estimator shows no systematic degradation across a fourfold range of model error, while both
prediction-based channels degrade monotonically from zero to five or six kilometres. That is the
behaviour hypothesis A3 predicted and did not find among prediction-based channels, and it is found
here in the one estimator that never consults a velocity model. But its absolute error is an order
of magnitude larger, so on these data model independence is purchasable and the price is high.

The limiting factor is identified: with a fixed azimuth tolerance of twenty degrees and twenty-five
random arrivals, roughly a fifth of the clutter passes the consistency test at any candidate.

**A tuned variant was attempted and is reported as unsuccessful.** Replacing the fixed tolerance by
a per-arrival uncertainty derived from the measured rectilinearity $1-\lambda_2/\lambda_1$, and the
grid consensus by a RANSAC search over pairwise ray intersections with local refinement, failed to
return an admissible solution in $62$ to $100$ per cent of realizations. Three specific causes were
identified and are recorded so that the reimplementation is well posed rather than exploratory:
the sign of a covariance eigenvector is arbitrary, so azimuth rays must be intersected as **lines**
and not as directed rays; the inlier gate and the minimum-support threshold are too strict at this
clutter level; and the Lipschitz admissibility test of C-III is a **same-phase** statement, so
applying it to an azimuth-selected inlier set that mixes P and S arrivals rejects valid
configurations. A correct implementation must assign a phase label to each inlier, or apply the
test within phase groups, before the gate is meaningful.

The untuned grid estimator reported in the table above is unaffected by these issues, since it
takes a genuine grid maximum, and its numbers stand.

Appendix U provides a positive instance of model-independent estimation using only measured directions
and a velocity bound, with no systematic degradation across the tested model-error range. Its absolute
location error is nevertheless about an order of magnitude larger in the tested configuration. A tuned
variant remains an open question. A back-azimuth estimator is exact under **vertically varying velocity
profiles** because the ray remains in the source-receiver vertical plane; lateral heterogeneity remains the
principal physical vulnerability.

**Declared limitations.** Eight realizations per level, one network, one clutter level, fixed
azimuth tolerance, no depth estimation, and a homogeneous inference model.

---

# Appendix V. Ruled-Out Approaches and Methodological Failure Analyses

Code: `c1c.py`, `b1_moment.py`, `octree_hypothesis_space.py`. This appendix records five approaches or methodological tests that were designed, analyzed or implemented and that do **not** work as originally intended. Each closes a method family or methodological shortcut that would otherwise be natural to attempt, so recording the failure has more value than silently dropping it.

## V.1 Second-order rank does not count sources, in either matrix

Section 6.4 establishes that the per-station three-by-three polarization covariance cannot bound
the number of sources. The natural repair is to move to the **inter-station** cross-spectral
matrix, where the independence of noise between distant sites is a physical fact rather than an
assumption: at a fixed frequency each mutually incoherent source contributes one steering vector,
so $K$ sources should give rank $K$, and the standard model-order criteria of array processing
should count them. This uses no picks, no arrival times, no origin time and no velocity model.

It fails. With ten stations, minimum-description-length order selection applied per frequency bin
recovers $K=1$ correctly with an eigenvalue gap of $2.1$, but returns $1$ for $K=2$ and $1$ to $2$
for $K=3$, with gaps of only $1.2$ to $1.6$. A sweep over the temporal separation of the events -
$0$, $1$, $3$ and $12$ s, with both incoherent and decoy-coherent clutter - changes nothing, which
rules out duty-cycle dilution as the cause.

The reason is structural. Rank-$K$ structure requires the steering vector to be **stable across
snapshots**. For a transient seismic wavetrain with moveout, each short window captures a different
part of the train at each station, so the effective steering vector varies from window to window
and the averaged matrix is neither rank $K$ nor cleanly separable. In radar and sonar the steering
vector is fixed because the source is continuously present; neither condition holds here.

**Consequence.** The diagnosis of Section 6.4 is reinforced rather than repaired: cardinality is
not recoverable from the rank of a second-order matrix, in either the component or the station
domain. This closes the subspace family - MUSIC, MDL, AIC, signal-subspace order selection - for
this problem, and it is the family a reader is most likely to propose.

## V.2 Betti reciprocity is not testable from station data, even with two events

Section 6.3 dismisses Betti-type identities on the grounds that a single event provides no second
elastodynamic state. Two overlapping events appear to remove that objection. They do not, and the
reason is analytic rather than empirical.

For two states with point moment-tensor sources, Betti's theorem reduces to

$$\mathbf{M}^1\!:\!\boldsymbol{\varepsilon}^2(\mathbf{x}_1)\;=\;\mathbf{M}^2\!:\!\boldsymbol{\varepsilon}^1(\mathbf{x}_2), \tag{V.1}$$

that is, an identity between the strain of each event **at the other event's source position**.
Those are unobserved points. Nor is there an identity among the four observations that two events
make at two stations: all four have their source at an event and their receiver at a station, so no
source-receiver exchange is available to exploit. Reciprocity therefore remains what Section 10.4
makes it, an exact symmetry of the learned representation, and does not become a data-level test.

## V.3 Moment-tensor physicality does not discriminate

The representation is linear in $\mathbf{M}$, so for any candidate position the moment tensor
follows from closed-form least squares with no search over origin time. A tectonic source is
deviatoric, so the isotropic fraction $|\mathrm{tr}\,\mathbf{M}|/\|\mathbf{M}\|$ should be near zero
for the truth and unconstrained for a false hypothesis. Tested against 491 false hypotheses more
than 6 km from the truth, with ten stations and thirty clutter arrivals:

| Statistic | true source | false hypotheses |
|:--|:--|:--|
| fit residual | 0.9984 | minimum 0.9985; **0 of 491 below the truth** |
| explained energy | 0.057 | maximum 0.045 |
| isotropic fraction | 0.334 | **142 of 491 below the truth** |

Origin time is sampled at a tenth of the wavelet period, since the least-squares fit is
matched-filter-like and therefore sampling-sensitive; at a quarter-period sampling the explained-
energy ratio between truth and best false hypothesis rises from $1.28$ to $1.42$, so coarse
sampling would have made the test look better than it is.

The residual and explained-energy ranking separates the truth perfectly, but that is a
six-parameter matched filter, an established technique and not the proposal here. The **physicality
test fails**: twenty-nine per cent of false hypotheses return a *more* deviatoric moment tensor than
the truth. The cause is that under clutter the least-squares fit at a false position is fitting
noise, and fitting noise with a six-element basis yields a near-deviatoric tensor often enough to
destroy the test.

**Consequence.** Source-physicality constraints (C-V) act on the source-time function, where they
were placed in Section 6.5, and not on the moment tensor geometry. The isotropic-fraction test is
not retained.

## V.4 An observability-based stopping rule for hypothesis-space refinement does not pay

Code: `octree_hypothesis_space.py`. Adaptive octree search over the source hypothesis space is
established: NonLinLoc subdivides the cell of highest posterior mass, and Qseek refines an octree
of candidate hypocentres guided by waveform stacking. The certificate of Section 6.1 suggests a
different use of the tree. It does not say where to search - refining where observability is high
would waste resolution on regions that may not contain the source - but it does say what resolution
is **achievable** at each point. That suggests a resolution-matched tree: branch on the score as
usual, but stop a branch once its cell reaches the local resolution limit
$\Delta_{\mathrm{res}}(\mathbf{x})=\kappa V_S/f_{\max}\sqrt{\gamma_{\mathrm{ref}}/\gamma(\mathbf{x})}$.

The tree itself is worth having and the stopping rule is not.

| Strategy | evaluations | median error |
|:--|:--|:--|
| uniform grid, 4 km | 432 | 1.00 km |
| uniform grid, 2 km | 3456 | 1.00 km |
| octree, score-branching, fixed depth | 457 | 0.85 km |

Score-guided refinement is $7.5\times$ cheaper than a uniform fine grid at equal or better accuracy,
which reproduces the established result rather than adding to it. The observability stopping rule,
however, **never fired**: across a $48$ km region the tree bottoms out at $0.75$ km cells while the
local resolution limit ranges from $0.087$ to $1.347$ km with a median of $0.274$ km, so no branch
ever reaches its limit. The apparent accuracy gain first attributed to the rule was entirely an
effect of tree depth, and the comparison against fixed depth was confounded by construction.

The failure is dimensional rather than computational. A resolution-based stopping rule can only
save work where the achievable resolution is **coarser** than the cell already reached, and at a
$2$ Hz working band with kilometre-scale cells that condition is essentially never met. The rule is
not retained, and the certificate's role remains the one established in Appendix S: mapping
where discrimination will be hard, in advance of any data.

## V.5 Location error on a grid measures grid alignment, not method quality

A systematic bias affects every location-error figure obtained from a grid-based estimator in this
document, and it is recorded here because it changes how those figures must be read.

In most pre-studies the true source sits at the origin and the search grids are symmetric about it,
so the truth falls **exactly on a grid node**. Measured with a 3.5 km grid and the constraint stack:

| True source relative to the grid | median error | 90th percentile |
|:--|:--|:--|
| exactly on a node | 0.00 km | 0.00 km |
| half a step in $x$ | 1.75 km | 1.75 km |
| half a step in $x$ and $y$ | 2.47 km | 2.47 km |
| arbitrary position | 1.63 km | 1.80 km |

The values are the quantization itself: half a step, and $1.75\sqrt2$ on the diagonal. **A
grid-based location error therefore measures where the source sits relative to the grid, not how
well the method works**, and the zeros reported in several tables are artefacts of alignment rather
than evidence of exact recovery.

**What is affected and what is not.** Every result whose endpoint is a **rate or a margin** is
unaffected, because those do not use the distance between the estimate and the truth: this covers
the clutter-immunity result of Appendix H, the constraint attribution of Appendix J, the minimal
sufficient set of Appendix M, the sufficiency ladder of Appendix O, the detection coupling of
Appendix R, the observability mapping of Appendix S and the coherent-clutter study of Appendix T -
that is, all of the load-bearing results. Affected are the location-error columns of Appendices L,
N and U. Appendices L and N already use continuous Nelder-Mead optimization rather than a grid
maximum and are therefore sound; the zeros in the travel-time and polarization columns of Appendix
U should be read as "within one grid cell", not as exact recovery, and the qualitative comparison
there - an order of magnitude between the model-free estimator and the others - is larger than the
quantization and stands.

**Registered for the confirmatory programme.** Synthetic sources are drawn from a continuous
distribution and never placed on search nodes, and location endpoints are computed by continuous
optimization rather than by grid maximum.

**Declared limitations of this appendix.** One clutter model and level, homogeneous medium,
far-field radiation patterns, and a single realization for V.3. V.2 is analytic and does not depend
on any of these. V.4 uses a single region size and band; the conclusion is specific to the regime
where cell size is already finer than the achievable resolution.

---

# Appendix W. Statistical Status of the Pre-Studies

Section 14 requires effect sizes, bootstrap intervals and multiplicity control for every
confirmatory claim. The pre-studies of Appendices G to V report rates from six to twelve
realizations without intervals, so this appendix supplies the missing analysis and states which of
their differences are supported and which are not. Clopper-Pearson intervals and Fisher's exact
test, two-sided.

**Differences that are supported.**

| Result | rates | difference [95% CI] | $p$ |
|:--|:--|:--|:--|
| R: pick-free vs travel time at miss rate 0.8 | 12/12 vs 0/12 | $+1.00\ [+0.47,+1.00]$ | $<10^{-4}$ |
| O: travel time vs polarization at $S=4$ | 7/8 vs 0/8 | $+0.88\ [+0.10,+1.00]$ | $0.0014$ |
| M: one-sided network, travel time vs polarization | 6/6 vs 0/6 | $+1.00\ [+0.08,+1.00]$ | $0.0022$ |
| T: coherent decoy at 4 km, travel time vs polarization | 6/6 vs 0/6 | $+1.00\ [+0.08,+1.00]$ | $0.0022$ |
| O: constraint stack vs waveform at $S=3$ | 6/8 vs 0/8 | $+0.75\ [-0.02,+0.97]$ | $0.0070$ |

**Differences that were unresolved at the initial replication counts, re-run at thirty.**
Code: `reps30_T.py`, `reps30_R.py`, `reps30_M.py`.

| Question | result at $n=30$ | verdict |
|:--|:--|:--|
| T: travel time vs polarization against a coherent decoy | $30/30$ vs $12/30$ spurious, $p<10^{-4}$ at **both** 4 km and 12 km | **supported** |
| R: pick-free advantage at miss rate $0.2$ | $+0.03$, $p=1.00$ | **not supported** |
| R: pick-free advantage at miss rate $0.3$ | $+0.13$, $p=0.11$ | **not supported** |
| R: pick-free advantage at miss rate $0.5$ | $+0.60$, $p<10^{-4}$ | **supported** |
| M: any separation between channels in spread or random geometry at $n=40$ | all channels $0.00$-$0.03$, $p=1.00$ | **no difference exists** |

The statistical audit yields three resulting qualifications for the affected appendices.

**T.** The spurious rate of the polarization channel is $0.40\ [0.23,0.59]$ and is **independent of
decoy separation**; the variation from $0.00$ at 4 km to $0.50$ at 12 km reported at six
realizations was sampling noise. The advantage over travel time, which is spurious in every one of
thirty realizations, is real and large.

**R.** The crossover threshold that D1 registers as its deliverable now has a value: the pick-free
advantage is **not** supported at miss rates $0.2$ or $0.3$ and is unambiguous at $0.5$, so the
threshold lies in $(0.3,0.5]$. The monotone trend reported earlier was driven by the supported end
of the range only.

**M.** In the spread and random geometries at high clutter **no channel differs from any other**,
because every channel succeeds: the task is saturated there. The equality of polarization with the
full stack in those cells is therefore not evidence that polarization suffices; it is evidence that
nothing is being discriminated. **The support for A2 comes from the one-sided cells**, where
the ablation condition (ABL) fails in every realization and polarization in none, and the hypothesis is
defined accordingly.

**Why, and what it implies for the design.** With six to twelve realizations the pre-studies can
resolve only large effects. The replication counts required at $\alpha=0.05$ and power $0.8$ are:

| $p_1$ vs $p_2$ | difference | $n$ per arm |
|:--|:--|:--|
| 1.00 vs 0.00 | 1.00 | 1 |
| 0.75 vs 0.00 | 0.75 | 2 |
| 0.75 vs 0.25 | 0.50 | 8 |
| 0.33 vs 0.00 | 0.33 | 6 |
| 0.50 vs 0.25 | 0.25 | 29 |
| 0.17 vs 0.00 | 0.17 | 11 |

The pre-studies therefore support differences of about $0.5$ and above, which is exactly the size
of the load-bearing results, and do not support the smaller differences of $0.08$ to $0.5$ that
appear in some of their tables. Those smaller differences are treated as **suggestive and unresolved**,
not as measurements, and the confirmatory programme sizes each condition from the
formula of Section 14 rather than from the pre-study counts.

**A limitation that spans every pre-study.** All sixteen numerical pre-studies use a **single synthetic generator**, one
clutter model, and - except for the two-layer case of Appendix Q and the gradient case of
Appendices I and L - a single homogeneous medium. Each appendix declares this locally; stated
globally, it means the sufficiency ladder and the channel hierarchy may be properties of the
generator rather than of the problem. Two independent generators are required before either is
claimed as a property of seismic source attribution: one with a layered medium, coda-derived rather
than synthetic clutter, and variable mechanisms, and one driven by real recorded noise. This is
registered as a confirmatory requirement, not as a caveat.

**Scope of this analysis.** Appendices T, R and M were re-run at thirty realizations. Appendices G,
I, L, N, O, P, Q and S retain their initial replication counts and report differences without
significance tests; by the power table above, only differences of about $0.5$ and larger in those
appendices should be read as established. Extending the thirty-realization treatment to them is
registered work.

**First real-data check (Appendix T.1).** The assumption shared by the pre-studies that real
clutter is not coherent with the events around it is supported, in the sense stated in T.1 (relative
to the preceding catalogued event), at the one Colombian station where it can be tested,
$\eta=0.054\ [-0.049,0.164]$ against a registered threshold of $0.20$. The pre-studies also assume,
without stating it, that the polarization of true arrivals is measurable; on real data that holds
at one of five stations. The supported differences below that rest on the polarization channel
therefore transfer to real data only at stations that pass the gate of Appendix T.1.

**What this does not change, and what Appendix Y does.** The five supported differences above
are statistically sound, but four of them compare against the coincidence-count ablation: the
sufficiency ladder at three and four stations, the failure of travel-time association in a
one-sided network, and its failure against a nearby coherent decoy. Appendix Y shows that in the
one-sided and four-station cells a continuous kinematic score without polarization closes the gap,
so those differences measure the weakness of coincidence counting rather than the value of
polarization. The coherent-decoy comparison of Appendix T has not been repeated against a
continuous kinematic score and is registered for that repetition. The dominance of pick-free
evidence at high miss rate (Appendix R) is unaffected, because its comparator failed for lack of
picks, not for lack of a continuous score.

---

---

# Appendix X. Published Associators as Baselines

Code: `task0_common.py`, `task0_run.py`, `task0_verify.py`, `task0_curves.py`. GaMMA [3] and
PyOcto [7] were executed on the same synthetic pick sets as the ablation condition (ABL) and the
constrained prototype (PHYS), thirty realizations per cell, with an operating-point sweep for every
method: score threshold for ABL and PHYS, `min_picks_per_eq` for GaMMA, `n_picks` for PyOcto.
Paired exact McNemar throughout, since all four methods see identical realizations. Positive
controls on a clean case passed for both baselines in every cell. Library defaults apart from the
velocity model, set to the generator's. Two asymmetries of this protocol are declared. The operating
point of each baseline (its knob value) was chosen on the evaluation realizations, whereas the
threshold of ABL and PHYS was calibrated on separate realizations; this favours the baselines. And
only the knob was swept: the other hyperparameters stayed at library defaults, which disfavours the
baselines and, as the oscillations below show, left them misconfigured in some cells. Appendix Y
repeats the comparison with every baseline tuned over a declared grid on calibration data only. **Two independent seeds** were run under the final
protocol and both are reported.

## X.1 The protocol was wrong twice before it was right, and the sequence is the finding

Three protocols were run. Reporting only the last would hide the two effects that matter most.

| Protocol | Threshold | Geometry | Outcome |
|:--|:--|:--|:--|
| P1 | absolute, mean-cosine polarization | one network per cell | PHYS false rate $1.00$ in nine of ten cells |
| P2 | absolute, summed polarization (24c) | one network per cell | PHYS better at the calibrated point in six of eight cells, worse on the curve in five of eight |
| **P3** | **relative to each realization's maximum** | **network drawn per realization** | **reported below; replicates across two seeds** |

**Why the threshold had to become relative.** ABL's score is a hit fraction, so its maximum is
exactly $1.000$ in every realization by construction. The constrained score multiplies by a
continuous polarization term and its maximum varies between $0.51$ and $0.68$ across realizations.
An absolute threshold therefore grants ABL a realization-independent operating point and denies the
constrained model any, which favours the ablation structurally. Thresholding a normalized statistic
is standard in detection and the argument does not depend on any result.

**Why the geometry had to be drawn per realization, which mattered more.** Under P1 and P2 each
cell fixed one network and varied only the clutter. Between two seeds the *same* method moved from
a false rate of $0.13$ to $0.70$ in the same cell: the particular geometry dominated the method.
Appendix S had already measured that the observability constant varies by a factor of eighty-six
across geometries, so this was foreseeable and was not foreseen. Under P3 the network is drawn per
realization, and ABL and PHYS agree between the two seeds cell by cell. PyOcto does not (one-sided
at clutter 40: 1.00 against 0.10; spread at clutter 40: 0.97 against 0.00), which indicates
sensitivity to its configuration rather than to the method.

**A defect diagnosed but demoted.** The mean-cosine polarization estimator carries a small-sample
bias: a mean over matched arrivals rewards candidates that explain less. The corrected form (24c)
sums instead of averaging. It is retained, but P1-to-P2 shows it was the secondary cause; the
threshold scale and the geometry were primary.

## X.2 Replication across seeds under the final protocol

False-event rate, two independent seeds:

| Cell | ABL | | **PHYS** | | GaMMA | | PyOcto | |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| | s1 | s2 | s1 | s2 | s1 | s2 | s1 | s2 |
| one-sided, clutter 15 | 0.30 | 0.37 | **0.00** | **0.00** | 0.00 | 0.00 | 0.57 | 0.70 |
| one-sided, clutter 40 | 0.37 | 0.43 | **0.00** | **0.00** | 0.00 | 0.00 | 1.00 | 0.10 |
| spread, clutter 40 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.97 | 0.00 |
| decoy 4 km | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| decoy 12 km | 0.00 | 0.00 | 0.03 | 0.03 | 0.00 | 0.00 | 0.03 | 0.00 |

The constrained prototype and the ablation now replicate to within one realization, against swings
of $0.57$ under the earlier protocol. That replication is itself the evidence that the instability
was the fixed geometry.

## X.3 Statement (i) against the ablation, and why it does not survive

In both one-sided cells the constrained prototype produces **no false event in sixty realizations
across two seeds**, against $0.30$-$0.43$ for the same search with the physical channels switched
off. McNemar exact $p=0.0039$ and $p=0.0010$, with unanimous discordance: nine and eleven
realizations where the ablation fails and the constrained model does not, **zero** in the opposite
direction.

Statement (i) of Section 13.0a therefore held against the ablation, in the geometry class this
document identifies as the difficult one, and Appendix Y.1 shows it was a tie effect:
ablations that break ties without polarization also produce no false event. Two further
qualifications apply. GaMMA also produces no
false event and recovers every event in both cells and both seeds (X.2), so there is no advantage
over the strongest published associator there. And at the calibrated threshold the constrained
score, which is continuous, emits a single event per realization in three of the four cell-seed
combinations, while the ablation emits every node tied at its maximum; Appendix Y.1 finds that a
tie-breaking ablation without polarization closes the gap entirely. In well-covered cells both are at zero and
nothing is being discriminated, consistent with Appendix W.

## X.4 Statement (ii): recall in the sparse regime, replicated

| Stations | ABL | **PHYS** | GaMMA | PyOcto |
|:--|:--|:--|:--|:--|
| 3 | 0.63 / 0.73 | **0.97 / 0.97** | 0.23 / 0.27 | 0.30 / 0.27 |
| 4 | 1.00 / 0.93 | 0.90 / 0.90 | 0.33 / 0.40 | 0.53 / 0.53 |
| 6 | 1.00 / 1.00 | 0.97 / 1.00 | 0.77 / 0.67 | 0.93 / 0.90 |
| 12 | 1.00 / 1.00 | 1.00 / 1.00 | 1.00 / 1.00 | 1.00 / 1.00 |

At three stations the constrained prototype recovers the event in twenty-nine of thirty
realizations in both seeds; GaMMA recovers it in seven and eight, PyOcto in nine and eight. At four
stations the constrained prototype attains a lower false-event rate than GaMMA in one seed, $0.43$
against $0.73$, $p=0.049$, and $0.30$ against $0.60$, $p=0.078$, in the other; at six stations
$0.03$ against $0.40$ and $0.37$ ($p=0.001$ and $0.006$). At six stations, however, GaMMA's recall
($0.67$-$0.77$) is below the ablation's ($1.00$), and at twelve stations it emits false events in 29
of 30 realizations in one seed and in none in the other; both point to its configuration rather than
to its method, and Appendix Y re-measures these cells with tuned baselines. At three stations every
method emits false events in $60$-$100\%$ of realizations (PHYS $0.80$ and $0.83$), so the recall
advantage there does not yet produce a usable catalogue.

**What this measures is narrower than a ladder against published software.** The published
associators, at library defaults, did not operate in the sparse regime on this generator; tuned
over a declared grid (Appendix Y.1) they still did not reach the calibration recall target there,
and a continuous kinematic kernel without polarization emitted fewer false events than the
constrained prototype at three and four stations.

**This comparison is a regime characterization, not a leaderboard.** It establishes that under this
generator, this station regime and this protocol, the published associators recovered far fewer
events in the sparse case. It does not establish superiority in the regimes they were designed for,
and no claim of that kind is made anywhere in this document.

## X.5 What remains unresolved

The curve metric, false-event rate at recall $0.98$, still returns "not attained" for the
constrained prototype in several ladder cells, and the two seeds disagree on it in the decoy cell at
12 km ($0.233$ against $0.033$) and in the ladder cells (three stations: not attained against
$0.900$; four: not attained against $0.533$; six: $0.167$ against $0.033$). The operating-point difficulty is reduced but not
eliminated, and the diagnostic of `task0_curves.py` locates a sign change at a threshold fraction
of about $0.8$ in that cell, which suggests the polarization channel has a validity regime rather
than a uniform effect.

**Declared limitations.** One synthetic generator, homogeneous medium, perfect picker, one clutter
model, thirty realizations per cell and two seeds, library-default hyperparameters for both
baselines, and a constrained model that is the diagnostic prototype of Section 13.0b rather than the
proposed model.

---

# Appendix Y. Controls for the Catalogue Endpoint: Ties, Timing Error, Empty Windows and Tuned Baselines

Code: `task1jk_controls.py`, `task1j2_noise.py`; decision rules in `PREREGISTRO_1j_1k.txt`, written
before each run and amended, with timestamps, where a design was changed before any evaluation
output of that part had been seen. The generator, the prototype and the protocol are those of
Appendix X unless stated. Thirty calibration and thirty evaluation realizations per cell, two
independent seeds (2024 and 7331), paired exact McNemar. GaMMA (source revision of August 2026) and
PyOcto 0.2.0 are tuned on the calibration realizations only, over a declared grid (GaMMA: minimum
pick count 5 to 14 and time-residual scale 1 or 2 s; PyOcto: minimum pick count 4 to 14 and
travel-time tolerance 0.5, 1 or 2 s), and frozen before evaluation. When no configuration reaches
the calibration recall target, the configuration with the highest calibration recall is used and
the cell is marked.

## Y.1 Ties: the one-sided advantage is a property of the score, not of polarization

The ablation ABL scores a hit fraction. In a one-sided network several grid nodes reach its maximum
exactly, and under the relative threshold all of them are emitted. The constrained prototype
multiplies by a continuous factor and has a unique maximum. Three ablations **without
polarization** remove the tie: TOL-TB subtracts $0.01$ times the RMS residual of the hits divided by
the tolerance, which cannot change the order between different hit fractions; TOL-K replaces the
hard tolerance by a Gaussian kernel of the residuals, $(1/2N)\sum\exp(-r^2/2\sigma^2)$ with
$\sigma=0.3$ s; CIII-TB is TOL-TB with the kinematic filter C-III. PHYS-S is the summed polarization
score (24c).

**Table Y.1. False events in 30 realizations, seed 2024 / seed 7331; in brackets, recall.** An
asterisk marks a method that did not reach the calibration recall target of 0.98 in at least one
seed.

| Cell | TOL | TOL-TB | TOL-K | CIII-TB | PHYS | PHYS-S | GAMMA | PYOCTO |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| one-sided, clutter 15 | 7/9 (30/30) | 0/0 (30/30) | 0/0 (30/30) | 0/0 (30/30) | 0/1 (30/30) | 0/0 (30/30) | 0/0 (30/30) | 22/0 (30/30)$^\ast$ |
| one-sided, clutter 40 | 12/12 (30/30) | 0/0 (30/30) | 0/0 (30/30) | 0/0 (30/30) | 0/0 (30/30) | 0/0 (30/30) | 2/2 (28/28) | 29/30 (28/25)$^\ast$ |
| spread, clutter 40 | 0/0 (30/30) | 0/0 (30/30) | 0/0 (30/30) | 0/0 (30/30) | 0/0 (30/30) | 0/0 (30/30) | 30/30 (28/28)$^\ast$ | 1/0 (29/30) |
| 3 stations | 30/30 (14/15)$^\ast$ | 30/30 (30/30) | 12/15 (30/28) | 30/30 (30/30) | 26/25 (25/27)$^\ast$ | 16/27 (23/29)$^\ast$ | 12/13 (3/2)$^\ast$ | 30/30 (7/9)$^\ast$ |
| 4 stations | 25/27 (29/28)$^\ast$ | 27/0 (30/30) | 0/0 (30/30) | 27/0 (30/30) | 10/7 (28/24) | 7/3 (28/27) | 25/23 (12/13)$^\ast$ | 30/13 (17/17)$^\ast$ |
| 6 stations | 4/1 (30/30) | 0/0 (30/30) | 0/0 (30/30) | 0/0 (30/30) | 0/9 (30/30) | 0/0 (30/30) | 8/11 (23/18)$^\ast$ | 0/1 (30/29) |

**Registered decision (tie control of Appendix D): met.** In both one-sided cells and both seeds,
every tie-breaking ablation without polarization produces no false event and recovers every event,
as the constrained prototype does (which has one false event in one seed); the difference between ABL and PHYS reported in X.3
($0.30$-$0.43$ against $0.00$) is therefore attributed to the ties of the discrete score, not to the
polarization channel. Two further results were not anticipated. At four stations the kernel score
TOL-K, with no polarization, produces no false event in either seed at full recall, against 10 and 7
for PHYS and 7 and 3 for PHYS-S; at three stations it produces 12 and 15 against 26 and 25 for PHYS.
And at six stations PHYS produces 9 false events in one seed where every tie-breaking ablation
produces none. With exact picks and an exact velocity model, a continuous travel-time score is at
least as good a discriminator as the polarization prototype in every cell tested, and better in the
sparse ones.

The reason is structural. The generator produces picks with **no timing error**, so the true source
has residual zero at every arrival and a continuous residual score separates it from any clutter
alignment almost perfectly. The hard tolerance of ABL (0.6 s) discards exactly that information,
which is why a coincidence count was a weak comparator throughout Appendices G to X.

**Published associators, tuned.** GaMMA reaches no false event in the one-sided cell at clutter 15
and 2 of 30 at clutter 40 in both seeds, but under the declared grid it does not reach the
calibration recall target in the spread cell, in the ladder cells or at six stations, and in the
spread cell it then emits a false event in every realization. PyOcto fails the target in at least one
seed in both one-sided cells and at three and four stations. Both baselines therefore remain sensitive to
configuration on this generator, whose picks carry no amplitude and no timing error; their numbers
are reported, not interpreted as properties of the methods.

## Y.2 Timing error: where polarization adds to the best kinematic score

Y.1 raises the question that decides G1: does polarization add anything to the **best** kinematic
score once picks are not exact? This analysis was designed after seeing Y.1 on seed 2024 and is
declared exploratory. Picks receive Gaussian timing errors of $0.1$ or $0.2$ s, and in one condition
the inference velocities are $5\%$ faster than the generating ones. TOL-K and PHYS-K,
$(1/2N)\sum\exp(-r^2/2\sigma^2)\cos^2\theta$ with C-III, select their kernel scale in calibration
among $0.15$, $0.3$ and $0.6$ s. PHYS-K is the kinematic kernel with polarization added as a
per-arrival weight, the closest the prototype comes to the likelihood (9d).

**Table Y.2. False events in 30 realizations, seed 2024 / seed 7331; in brackets, recall. Last
column: discordant pairs (TOL-K false and PHYS-K not : the reverse) and McNemar $p$, per seed.**

| Cell | timing error | PHYS | PHYS-S | TOL-K | PHYS-K | TOL-K vs PHYS-K |
|:--|:--|:--|:--|:--|:--|:--|
| one-sided, clutter 40 | 0.1 s | 0/0 (30/30) | 0/0 (30/30) | 0/0 (30/30) | 1/0 (29/30) | 0:1, $p$=1; 0:0, $p$=1 |
| one-sided, clutter 40 | 0.2 s | 0/0 (30/30) | 0/0 (30/30) | 1/0 (29/30) | 1/0 (29/30) | 1:1, $p$=1; 0:0, $p$=1 |
| one-sided, clutter 40 | 0.1 s, velocity +5% | 18/1 (30/29) | 8/1 (30/29) | 0/0 (30/30) | 1/0 (29/30) | 0:1, $p$=1; 0:0, $p$=1 |
| 4 stations | 0.1 s | 10/8 (27/23) | 7/3 (28/27) | 3/0 (30/30) | 2/1 (28/29) | 2:1, $p$=1; 0:1, $p$=1 |
| 4 stations | 0.2 s | 15/6 (28/26) | 10/4 (29/29) | 19/6 (29/30) | 3/3 (27/27) | 17:1, $p$=0.00014; 5:2, $p$=0.45 |
| 4 stations | 0.1 s, velocity +5% | 20/9 (29/25) | 25/5 (30/26) | 5/6 (27/29) | 4/5 (29/30) | 4:3, $p$=1; 5:4, $p$=1 |
| 6 stations | 0.1 s | 0/6 (30/30) | 0/0 (30/30) | 0/0 (30/30) | 0/0 (30/30) | 0:0, $p$=1; 0:0, $p$=1 |
| 6 stations | 0.2 s | 0/4 (30/30) | 0/0 (30/30) | 0/1 (30/29) | 0/1 (30/29) | 0:0, $p$=1; 0:0, $p$=1 |
| 6 stations | 0.1 s, velocity +5% | 0/1 (30/29) | 0/0 (30/30) | 0/0 (30/30) | 0/0 (30/30) | 0:0, $p$=1; 0:0, $p$=1 |

**Registered decision: not met.** The rule asked for at least $0.05$ fewer false events for PHYS-K
than for TOL-K at recall non-inferior within $0.02$ in at least half of the nine cells, in both
seeds. In seven of nine cells both methods are at or near zero and nothing is being discriminated; in the
eighth (four stations, $0.1$ s with velocity error) they are equal within one realization.
The exception is the sparsest and noisiest cell, four stations with $0.2$ s timing error, where
PHYS-K produces 3 false events against 19 for TOL-K in one seed ($p=1.4\times10^{-4}$) and 3
against 6 in the other ($p=0.45$), in both cases at a recall cost of two or three realizations in
thirty, which exceeds the registered non-inferiority margin. Two further observations. The prototypes with a hard tolerance, PHYS and PHYS-S, are fragile under
velocity error (18 and 8 false events of 30 in the one-sided cell in one seed, 20 and 25 at four
stations, where both kernel scores produce at most six), because a hard tolerance turns a
systematic timing bias into missed hits. And in this single-source design, adding polarization to
the kernel score never produced more than one additional false event in any cell; Y.3 shows that
this does not hold in the multi-source design.

## Y.3 Empty windows and an absolute threshold

Every realization now contains $K$ sources, $K\in\{0,1,2\}$ with probabilities $0.4$, $0.4$ and
$0.2$, at random positions and origin times, so that between a third and a little over a half of the windows are
empty (20 to 32 of 60). The threshold is **absolute**, calibrated on sixty separate windows to an event recall of
$0.90$ with the fewest false events, and applied unchanged to sixty evaluation windows. The endpoint
is false events per window at the event recall attained. PHYS-K, with a fixed kernel scale of $0.3$
s, is added; one condition adds $0.2$ s of timing error to every pick.

**Table Y.3. False events per window, seed 2024 / seed 7331; in brackets, event recall.** An asterisk
marks a method that did not reach the calibration recall target in at least one seed, in which case
its highest-recall threshold was used.

| Cell | timing error | TOL | TOL-TB | TOL-K | CIII-TB | PHYS | PHYS-S | PHYS-K | GAMMA | PYOCTO |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| one-sided, clutter 40 | none | 0.50/0.43 (1.00/0.96) | 0.00/0.02 (0.96/0.98) | 0.00/0.02 (0.96/0.98) | 0.00/0.02 (0.96/0.98) | 0.50/0.32 (0.94/0.84) | 0.32/0.33 (0.92/0.88) | 0.00/0.07 (0.90/0.92) | 0.02/0.02 (0.85/0.96) | 2.97/1.67 (0.50/0.61)$^\ast$ |
| one-sided, clutter 40 | 0.2 s | 0.22/0.18 (1.00/1.00) | 0.00/0.02 (0.94/0.98) | 0.02/0.03 (0.94/1.00) | 0.00/0.08 (0.94/1.00) | 0.40/0.38 (0.92/0.84) | 0.28/0.28 (0.96/0.86) | 0.08/0.32 (0.92/0.90) | 0.02/0.03 (0.94/0.94) | 0.18/3.10 (0.58/0.61)$^\ast$ |
| spread, clutter 40 | none | 0.00/0.00 (1.00/1.00) | 0.00/0.00 (0.92/0.88) | 0.02/0.00 (1.00/1.00) | 0.00/0.00 (0.92/0.88) | 0.27/0.02 (0.90/0.86) | 0.03/0.02 (0.96/0.92) | 0.08/0.00 (0.96/0.94) | 0.02/0.03 (0.83/0.86) | 2.93/0.03 (0.67/0.73)$^\ast$ |
| 4 stations | none | 4.07/4.17 (0.77/0.93)$^\ast$ | 0.28/0.22 (0.88/1.00) | 0.35/0.28 (0.85/0.97) | 0.28/0.22 (0.88/1.00) | 1.13/4.12 (0.87/0.97)$^\ast$ | 0.70/1.07 (0.83/0.95) | 0.52/0.35 (0.83/0.93) | 2.00/0.72 (0.37/0.40)$^\ast$ | 0.32/0.40 (0.52/0.45)$^\ast$ |
| 4 stations | 0.2 s | 3.97/4.08 (0.87/0.90)$^\ast$ | 0.63/0.35 (0.96/0.95) | 1.07/0.65 (0.87/0.93) | 0.63/0.35 (0.96/0.95) | 1.13/0.87 (0.88/0.88) | 0.78/0.63 (0.81/0.85) | 0.70/4.12 (0.83/1.00)$^\ast$ | 2.07/0.95 (0.29/0.38)$^\ast$ | 0.40/1.83 (0.40/0.28)$^\ast$ |
| 6 stations | none | 0.25/0.07 (1.00/1.00) | 0.10/0.02 (0.96/0.83) | 0.17/0.08 (0.94/0.97) | 0.10/0.02 (0.96/0.83) | 0.48/0.72 (0.92/0.94) | 0.37/0.18 (0.94/0.89) | 0.20/0.17 (0.94/0.89) | 1.82/1.78 (0.64/0.71)$^\ast$ | 0.05/0.05 (0.68/0.71)$^\ast$ |

**Empty windows.** With an absolute threshold, methods whose calibration reached its target
emitted an event in at most one of the 20 to 32 empty windows of any cell, and in none in most; the
systematic emissions in empty windows come from methods that failed calibration and fell back to a
permissive threshold: ABL at four stations, whose ties keep its maximum at one; PHYS at four
stations and PHYS-K at four stations with $0.2$ s error, each in one seed; GaMMA in the sparse
cells; and PyOcto in both one-sided cells and the spread cell. The objection that the prototype
cannot run on continuous data is therefore an objection to the relative threshold of P3, not to the
scores, and the absolute threshold is the registered form for real data.

**Registered decision on G1' (Appendix D): not met, and reversed in direction.** The best
tie-breaking kinematic score without polarization (TOL-TB or CIII-TB) has the lowest or equal false
events per window of all the scores computed on the grid in every cell and seed (PyOcto is lower in
two seed-cells, at recall $0.40$ to $0.68$), at an event recall at most $0.07$ lower than that of
PHYS-K (lower in the spread cell, at four stations with $0.2$ s error in one seed and at six
stations in one seed, higher or equal elsewhere). In no cell does a score with
polarization emit fewer false events than it; in four seed-cells the kinematic score emits
significantly fewer (exact sign test: four stations with exact picks, $p=0.03$ in one seed;
one-sided with $0.2$ s error, $p=5\times10^{-4}$; six stations, $p=0.008$; and four stations with
$0.2$ s error, where PHYS-K failed calibration in one seed). The multiplicative prototype PHYS is among the
worst of the scores in most cells. Tuned GaMMA is competitive in the one-sided and spread cells
($0.02$-$0.03$ false events per window at event recall $0.83$-$0.96$) and does not reach the recall
target in the sparse ones; PyOcto does not reach it in any cell under the declared grid.

## Y.4 Why Y.1 to Y.3 were not the last word: the generator was inconsistent

Y.1 to Y.3 suggested that polarization adds nothing to the best kinematic score. Before accepting a
negative result about a physical channel, the channel and the generator were audited with the same
suspicion that had been applied to the positive results of Appendix X, and two defects were found.

**The generator gave the kinematic channel information no picker has.** Every true arrival was
picked, with its exact time, even at stations where its amplitude over the noise is $1.1$ to $2.5$
(sources at 30 to 45 km). The polarization channel, in contrast, is measured on that same noisy
waveform: on the generator, $|\cos|$ between measured and predicted polarization separates true
arrivals from clutter with an area under the ROC curve of only $0.78$, and only $0.47$ of P arrivals
fall within $20^\circ$ in horizontal azimuth (`diag_pol.py`; the real station ZAR reaches $0.79$).
One channel was simulated at infinite signal-to-noise ratio and the other at the true one.

**The combination was not a likelihood.** The prototype multiplies a coincidence count by a mean
cosine, and PHYS-K of Y.2 multiplies a kernel by $\cos^2$; neither is the per-arrival likelihood
ratio (9d), whose polarization term is estimated from the data and enters additively in log form.

## Y.5 A physically consistent generator and the likelihood form

Code: `task1o_consistent.py`, `task1o_baselines.py`, `analyze_1o.py`. The design was motivated by
the audit above, after Y.1 to Y.3 were seen, and its decision rules were written before it ran; its
endpoint was amended once, after the first seed and before any result of the second, as recorded
in `PREREGISTRO_1j_1k.txt`. The second seed is the replication.

**Generator.** As in Y.3, with zero to two sources per window, but a true arrival is now detected
with probability $\Phi\big((a-2)/0.5\big)$, where $a$ is its amplitude over the noise, and its pick
time carries a Gaussian error of three times the Cramer-Rao bound of the pulse at that amplitude
(at least $0.02$ s). Between 43 and 51 per cent of true arrivals are detected, most losses at the
far stations. A **positive control** keeps every arrival and adds a $0.5$ s timing error, a regime in
which a correctly implemented polarization channel must help.

**Methods, with symmetric effort.** KIN-LR is the robust kinematic likelihood
$\sum_i\log\big[\epsilon+(1-\epsilon)\,\mathcal N(r_i;\sigma)/u\big]$, with $u$ the clutter pick
density and $\sigma$ chosen in calibration. PHYS-LR is the same with the polarization ratio of (9d)
multiplying the kinematic term, $f(c)/g(c)$ estimated in calibration from true arrivals and from
clutter. WAVE-E and WAVE-P need no picks: they sum the matched-filter energy of the three components,
or its projection on the predicted polarization, around each predicted arrival. TOL-TB is the
reference of Y.1; GaMMA and PyOcto are tuned on calibration over the grid of Y.1. **Endpoint:**
Neyman-Pearson, as in operational detection: for budgets of $0.05$, $0.1$, $0.25$ and $0.5$ false
events per window, the threshold (or configuration) meeting the budget in calibration is applied to
sixty evaluation windows and event recall is measured. Because the realized false rate can exceed
the budget, Table Y.5 reports the recall reached at a **realized** false rate of at most $0.1$ and
$0.25$ per window, which no method can game by overshooting its budget.

**Table Y.5. Event recall at realized false events per window $\le0.1$ / $\le0.25$; seed 2024 ;
seed 7331.** A dash marks a cell not run (the control cells of the published associators were run
for one seed only).

| Generator | Cell | TOL-TB | KIN-LR | PHYS-LR | WAVE-E | WAVE-P | GAMMA | PYOCTO |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| consistent | one-sided, 8 st. | 0.17/0.17 ; 0.35/0.35 | 0.36/0.36 ; 0.43/0.57 | 0.38/0.57 ; 0.59/0.68 | 0.26/0.26 ; 0.32/0.41 | 0.00/0.38 ; 0.38/0.62 | 0.28/0.28 ; 0.41/0.41 | 0.09/0.19 ; 0.27/0.27 |
| consistent | 4 stations | 0.04/0.04 ; 0.04/0.04 | 0.04/0.04 ; 0.00/0.00 | 0.22/0.36 ; 0.00/0.33 | 0.00/0.18 ; 0.11/0.11 | 0.16/0.16 ; 0.20/0.47 | 0.00/0.00 ; 0.02/0.02 | 0.00/0.07 ; 0.00/0.00 |
| consistent | 6 stations | 0.06/0.17 ; 0.02/0.19 | 0.08/0.08 ; 0.04/0.21 | 0.13/0.32 ; 0.25/0.35 | 0.08/0.15 ; 0.21/0.21 | 0.08/0.08 ; 0.35/0.35 | 0.09/0.11 ; 0.04/0.04 | 0.00/0.09 ; 0.00/0.02 |
| control, 0.5 s | one-sided, 8 st. | 0.02/0.32 ; 0.30/0.59 | 0.00/0.30 ; 0.46/0.46 | 0.00/0.57 ; 0.32/0.57 | 0.26/0.26 ; 0.32/0.41 | 0.00/0.38 ; 0.38/0.62 | 0.00/0.74 ; - | 0.00/0.00 ; - |
| control, 0.5 s | 4 stations | 0.02/0.16 ; 0.02/0.02 | 0.02/0.18 ; 0.02/0.02 | 0.13/0.13 ; 0.00/0.00 | 0.00/0.18 ; 0.11/0.11 | 0.16/0.16 ; 0.20/0.47 | 0.00/0.24 ; - | 0.00/0.00 ; - |
| control, 0.5 s | 6 stations | 0.36/0.36 ; 0.02/0.25 | 0.21/0.45 ; 0.02/0.44 | 0.42/0.58 ; 0.00/0.38 | 0.08/0.15 ; 0.21/0.21 | 0.08/0.08 ; 0.35/0.35 | 0.25/0.57 ; - | 0.00/0.62 ; - |


**How this run was first read, and why it was repeated.** At the nominal budget of $0.1$ false
events per window, PHYS-LR detected $0.18$ to $0.25$ more events than KIN-LR in five of six
seed-cells, and more than GaMMA and PyOcto, and that was first written as evidence for G1. An
independent verification found five defects in that reading: PHYS-LR exceeded its false-event budget
in evaluation, so the nominal comparison favoured it; there was no negative control to separate the
polarization term from the added flexibility of the likelihood; KIN-LR was a weaker reference than
GaMMA; GaMMA and PyOcto searched a wider domain than the grid methods and started at a minimum of five
picks; and the change of generator and the change of likelihood form entered together. Y.6 repeats
the comparison without them, and Table Y.5 is kept as the record of the run, read only at the
realized false rate. One observation from it stands: under the $0.5$ s control, where every arrival
is picked but timing is poor, GaMMA is the best method in the one-sided cell ($0.74$ recall at $0.25$
false events per window against $0.57$ for PHYS-LR) and matches PHYS-LR at six stations.

## Y.6 A fair repetition with a negative control (Task 1r)

Code: `task1r_fair.py`, `task1r_baselines.py`, `froc_exact.py`, `analyze_1r200.py`; rules R-a to R-d
written before the run. Five changes with respect to Y.5. (1) The endpoint is the FROC curve on the
evaluation windows, read as the event recall at a **realized** rate of at most $0.1$ and $0.25$
false events per window, with paired bootstrap intervals over windows ($B=2000$); the variant of each
likelihood (kernel width $\sigma\in\{0.1,0.2,0.4\}$ s, outlier weight
$\epsilon\in\{0.05,0.1,0.2\}$) is chosen on sixty calibration windows. (2) A negative control,
PHYS-SH, applies the same polarization ratio to cosines permuted at random among arrivals: it has
every degree of freedom of PHYS-LR and none of its geometry. (3) An exact condition, with every
arrival detected and timed exactly, separates the effect of the generator from that of the
likelihood form. (4) GaMMA and PyOcto search the same domain and origin-time window as the grid
methods, from a minimum of three picks, and their configuration is chosen with an oracle on the
evaluation windows, which favours them. (5) The ratio $f/g$ is estimated on forty windows of its own.

Two defects of the tooling were found and corrected during this task. GaMMA initializes its mixture
with numpy's global random generator, unseeded when it runs in a single process, so its catalogue
depended on everything the process had run before: in 80 calls, 15 changed their hits or false events
when that state changed. It is now seeded before every call, the convention GaMMA uses in its own
multi-process mode, and two independent code paths then give identical catalogues in all 80. And the
registered threshold sweep, 80 quantiles of the pooled grid scores, sampled the FROC curve with gaps
of up to $0.14$ false events per window at its strict end, so that a reading at $0.1$ could create or
erase a difference. Because the catalogue at any threshold is a prefix of the greedy list of each
window, the exact curve is obtained by using every score in those lists as a threshold; it was
verified identical to the original emission at every threshold, and the rules are evaluated on it.
This amendment was written after the registered sweep had shown rule T-b of Y.8 failing in one cell,
and it decides four outcomes: with the registered sweep R-d (Y.6), T-b (Y.8), U-a and U-b (Y.9) are
not met, and with the exact curve they are. Both readings are in the result files. The exact curve is
the correct measurement of the registered endpoint, but the order of events is part of the record.
The bootstrap intervals re-select each method's operating point in every replicate; holding it fixed,
as first computed, narrowed them by fifteen to forty per cent.

With the sixty evaluation windows per cell first registered, the result was undetermined: no rule was
met in both seeds and the intervals were wide. Two hundred windows per cell were then run, the sixty
plus 140 new ones, with the rules unchanged.

**Table Y.6. Event recall at realized false events per window $\le0.1$ / $\le0.25$, 200 evaluation
windows per cell; seed 2024 ; seed 7331.** GaMMA seeded; GaMMA and PyOcto with oracle configuration.

| Cell | events | KIN-LR | PHYS-LR | PHYS-SH | GaMMA | PyOcto |
|:--|:--|:--|:--|:--|:--|:--|
| one-sided, 8 st. | 172 ; 153 | 0.31/0.43 ; 0.36/0.42 | 0.46/0.55 ; 0.45/0.60 | 0.39/0.47 ; 0.37/0.48 | 0.35/0.39 ; 0.35/0.36 | 0.23/0.28 ; 0.28/0.35 |
| 4 stations | 170 ; 151 | 0.08/0.13 ; 0.07/0.15 | 0.23/0.32 ; 0.25/0.33 | 0.18/0.22 ; 0.20/0.32 | 0.05/0.08 ; 0.11/0.11 | 0.02/0.09 ; 0.02/0.06 |
| 6 stations | 180 ; 155 | 0.17/0.23 ; 0.17/0.32 | 0.27/0.47 ; 0.44/0.54 | 0.21/0.31 ; 0.29/0.39 | 0.16/0.19 ; 0.17/0.19 | 0.11/0.21 ; 0.08/0.24 |

**Table Y.7. Differences in event recall at realized false events per window $\le0.25$, with paired
bootstrap 95% intervals over windows; seed 2024 ; seed 7331.**

| Cell | PHYS-LR - KIN-LR | PHYS-LR - PHYS-SH | PHYS-SH - KIN-LR | PHYS-LR - GaMMA |
|:--|:--|:--|:--|:--|
| one-sided, 8 st. | +0.12 [+0.05, +0.19] ; +0.18 [+0.12, +0.27] | +0.08 [+0.02, +0.15] ; +0.12 [+0.06, +0.19] | +0.03 [-0.02, +0.09] ; +0.07 [+0.01, +0.15] | +0.16 [+0.10, +0.24] ; +0.24 [+0.17, +0.34] |
| 4 stations | +0.19 [+0.10, +0.26] ; +0.19 [+0.12, +0.28] | +0.10 [+0.02, +0.16] ; +0.01 [-0.06, +0.08] | +0.09 [+0.03, +0.15] ; +0.18 [+0.12, +0.27] | +0.24 [+0.16, +0.30] ; +0.22 [+0.14, +0.32] |
| 6 stations | +0.24 [+0.13, +0.33] ; +0.22 [+0.14, +0.30] | +0.17 [+0.07, +0.25] ; +0.15 [+0.08, +0.23] | +0.07 [+0.01, +0.15] ; +0.06 [+0.01, +0.14] | +0.28 [+0.16, +0.36] ; +0.35 [+0.25, +0.44] |

**Registered decisions.** R-a is met: PHYS-LR exceeds KIN-LR by at least $0.05$ at both caps in all
six seed-cells, and every interval at $0.25$ excludes zero. R-b, the negative control, is **not**
met: PHYS-SH also exceeds KIN-LR, by $0.03$ to $0.18$ at $0.25$. The second verification showed why,
and that R-b was the wrong test: PHYS-SH permutes the ratios among the station-phase slots of the same
hypothesis, so the set of polarization ratios of that hypothesis survives and only their pairing with
the timing residuals is lost. It has exactly the nine variants of KIN-LR, so its gain is not
flexibility. PHYS-LR minus PHYS-SH, $+0.08$ to $+0.17$ at $0.25$ with intervals excluding zero in
five of six seed-cells and $+0.01$ in the sixth (four stations, seed 7331), is therefore a lower
bound of the polarization effect in this generator. Y.10 measures it with a control that takes the
predicted directions from another grid node, and with sources that are not pure SV radiators. R-d is met: PHYS-LR exceeds the
better of GaMMA and PyOcto by at least $0.05$ at both caps in all six seed-cells, by $0.16$ to $0.30$
at $0.25$. R-c, from the sixty-window run with the
registered sweep: with every arrival detected and timed exactly, PHYS-LR does not beat KIN-LR
(differences $0$, $-0.04$, $-0.04$ and $0$, $+0.16$, $0$ at $0.25$).

**The advantage over GaMMA is not a threshold artefact.** A grid likelihood has a continuous
threshold and GaMMA ten configurations, so a reading at a cap could favour the method with more
operating points. It does not here: the highest recall that any GaMMA configuration reaches at any
false rate is $0.39$, $0.08$, $0.19$ (seed 2024) and $0.36$, $0.11$, $0.19$ (seed 7331) in the three
cells, and PHYS-LR exceeds it at $0.25$ false events per window in all six. What limits GaMMA is the
events it never forms. Y.10 shows that the per-pick gain of this section depends on the pure-SV S
waves of this generator and does not survive realistic source mechanisms.

## Y.7 Polarization as a filter on the strongest associator (Task 1s)

Code: `task1s_gamma_pol.py`, `check_gamma_det.py`, `analyze_1s.py`. The most direct way to put the
physics on top of the state of the art is to filter its output: every event that GaMMA emits is
scored by the sum of $\log f/g$ over the picks GaMMA assigned to it, and emitted only above a
threshold (GAMMA+POL); the control scores the same events with cosines drawn at random from those
observed (GAMMA+SH). Same windows and seeds as Y.6, 200 per cell, GaMMA seeded.

**Table Y.8. Event recall at realized false events per window $\le0.05$ / $\le0.1$ / $\le0.25$; AUC
of the event score between true and false GaMMA events; highest recall of any GaMMA configuration at
any false rate; seed 2024 ; seed 7331.**

| Cell | GaMMA | GaMMA+POL | GaMMA+SH | AUC POL | AUC SH | GaMMA ceiling |
|:--|:--|:--|:--|:--|:--|:--|
| one-sided, 8 st. | 0.00/0.35/0.39 ; 0.35/0.35/0.36 | 0.39/0.39/0.39 ; 0.35/0.36/0.36 | 0.32/0.35/0.39 ; 0.35/0.35/0.36 | 0.93 ; 0.96 | 0.62 ; 0.61 | 0.39 ; 0.36 |
| 4 stations | 0.01/0.05/0.08 ; 0.09/0.11/0.11 | 0.05/0.08/0.08 ; 0.09/0.11/0.11 | 0.05/0.06/0.08 ; 0.09/0.11/0.11 | 0.70 ; 0.74 | 0.54 ; 0.49 | 0.08 ; 0.11 |
| 6 stations | 0.14/0.16/0.19 ; 0.00/0.17/0.19 | 0.17/0.19/0.19 ; 0.17/0.18/0.19 | 0.14/0.16/0.19 ; 0.15/0.17/0.19 | 0.89 ; 0.91 | 0.54 ; 0.52 | 0.19 ; 0.19 |

The polarization channel is informative on GaMMA's own events: it separates the true from the false
ones with an AUC of $0.89$ to $0.96$ at six and eight stations and $0.70$ to $0.74$ at four; the
control reaches $0.49$ to $0.62$, its excess over $0.5$ reflecting that true events carry more picks.
Yet rule S-a is not met in any cell: at $0.1$ and $0.25$ false events per window GAMMA+POL gains at
most $0.04$. The reason is in the last column. A filter can only remove events, so its recall is
bounded by the highest recall GaMMA reaches in any configuration, and GaMMA reaches that ceiling at
no more than $0.22$ false events per window in every cell. The events GaMMA loses are mostly never
formed, not formed and rejected: in the six-station cell, where most true sources have three or more
detected picks, most misses have no GaMMA event within 20 km, and part of the rest is mislocation
(with a 10 km instead of a 5 km match tolerance GaMMA's recall there roughly doubles). The loss is in
the association of few true picks amid clutter, without amplitudes or pick probabilities, which is
the regime these cells simulate. Polarization pays when it enters the hypothesis test itself (Y.6),
where a three-component station constrains direction as well as time; as a filter on a kinematic
catalogue it cannot, because the hypotheses it would rescue were never formed. This section used the
clutter of E26: GaMMA received clutter picks that no picker would produce, and with consistent weak
clutter it recovers $0.54$ to $0.63$ of the events in the one-sided and six-station cells (Y.12). The
conclusion that a filter cannot raise recall above GaMMA's ceiling stands; the size of GaMMA's loss
does not.

## Y.8 Physics before the picker (Task 1t)

Code: `task1t_wave.py`, `analyze_1t.py`. If the loss is at the picker, the physics can act before it.
Three pick-free scores migrate the matched-filtered three-component waveform over the hypothesis grid
and differ only in the use of polarization. WAVE-E stacks the energy of the three components around
each predicted arrival, the non-directional form of migration and stacking, here with the known
wavelet as matched filter. It is an ablation, not a reimplementation of the published locators:
QuakeMigrate [109] and LOKI [111] stack onset functions with P taken on the vertical component and S
on the horizontals, and MALMI [110] migrates learned phase probabilities. WAVE-P stacks its
projection on the polarization direction that the hypothesis predicts at each station. WAVE-PSH, the
negative control, projects on the direction predicted for another station, which keeps the
distribution of directions and breaks the geometry. This script reproduced exactly the WAVE-P of Y.6
on the sixty windows of the registered run. Rules T-a to T-c were written after WAVE-P was seen to be
among the best methods of Y.6, and are declared as such.

**Table Y.9. Pick-free migration. Event recall at realized false events per window $\le0.1$ /
$\le0.25$, 200 evaluation windows per cell; differences at $\le0.25$ with paired bootstrap 95%
intervals; seed 2024 ; seed 7331.**

| Cell | WAVE-E | WAVE-P | WAVE-PSH | GaMMA | PyOcto | WAVE-P - WAVE-E | WAVE-P - WAVE-PSH | WAVE-P - GaMMA |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| one-sided, 8 st. | 0.24/0.35 ; 0.23/0.37 | 0.36/0.48 ; 0.42/0.55 | 0.31/0.42 ; 0.23/0.48 | 0.35/0.39 ; 0.35/0.36 | 0.23/0.28 ; 0.28/0.35 | +0.13 [+0.07, +0.23] ; +0.18 [+0.13, +0.29] | +0.06 [-0.02, +0.14] ; +0.07 [-0.01, +0.21] | +0.09 [+0.01, +0.21] ; +0.19 [+0.11, +0.31] |
| 4 stations | 0.10/0.16 ; 0.11/0.17 | 0.21/0.30 ; 0.25/0.35 | 0.10/0.16 ; 0.10/0.17 | 0.05/0.08 ; 0.11/0.11 | 0.02/0.09 ; 0.02/0.06 | +0.14 [+0.07, +0.19] ; +0.18 [+0.09, +0.26] | +0.14 [+0.06, +0.20] ; +0.18 [+0.12, +0.27] | +0.22 [+0.14, +0.28] ; +0.24 [+0.17, +0.34] |
| 6 stations | 0.14/0.24 ; 0.17/0.33 | 0.31/0.46 ; 0.44/0.53 | 0.13/0.29 ; 0.23/0.37 | 0.16/0.19 ; 0.17/0.19 | 0.11/0.21 ; 0.08/0.24 | +0.21 [+0.14, +0.29] ; +0.20 [+0.13, +0.30] | +0.17 [+0.11, +0.27] ; +0.15 [+0.08, +0.25] | +0.26 [+0.17, +0.36] ; +0.34 [+0.25, +0.45] |

**Registered decisions.** T-a is met: WAVE-P exceeds WAVE-E by at least $0.05$ at both caps in all six
seed-cells, by $0.13$ to $0.21$ at $0.25$ with every interval excluding zero. T-b is met on the exact
curve: WAVE-P exceeds its control by $0.06$ to $0.18$, with intervals excluding zero in four
seed-cells and including it in both one-sided cells, where the stations share similar azimuths and a
permutation of directions changes little. T-c is met: WAVE-P exceeds the better of GaMMA and PyOcto by
at least $0.05$ at both caps in five of six seed-cells; in the sixth (one-sided, seed 2024) it ties
GaMMA at $0.1$ ($0.36$ against $0.35$) and leads at $0.25$. At $0.25$ false events per window WAVE-P
also exceeds the highest recall GaMMA reaches at any false rate in all six. WAVE-P projects S on the
SV direction, which this generator radiates exactly; Y.10 replaces it by a mechanism-independent form.

## Y.9 A medium the model does not know (Task 1u)

Code: `task1u_hetero.py`, `analyze_1u.py`. The generator of Y.8 is a homogeneous medium that the
inference model knows exactly, which favours every waveform score. Here each station and phase
receives a delay drawn from $\mathcal N(0,\,0.15\ \mathrm{s})$ and a deviation of its polarization
direction of about ten degrees, fixed in each window and unknown to every method; picks are degraded
from the perturbed times as before. One hundred evaluation windows per cell; rules U-a and U-b were
written before the run.

**Table Y.10. Unmodelled station delays and polarization deviations. Event recall at realized false
events per window $\le0.1$ / $\le0.25$; difference at $\le0.25$ with paired bootstrap 95% interval;
seed 2024 ; seed 7331.**

| Cell | WAVE-E | WAVE-P | WAVE-PSH | GaMMA | WAVE-P - GaMMA |
|:--|:--|:--|:--|:--|:--|
| one-sided, 8 st. | 0.29/0.33 ; 0.21/0.33 | 0.31/0.45 ; 0.37/0.56 | 0.23/0.39 ; 0.30/0.51 | 0.37/0.39 ; 0.33/0.33 | +0.06 [-0.07, +0.18] ; +0.23 [+0.09, +0.46] |
| 4 stations | 0.12/0.20 ; 0.10/0.14 | 0.22/0.28 ; 0.23/0.34 | 0.07/0.14 ; 0.09/0.13 | 0.08/0.10 ; 0.07/0.07 | +0.17 [+0.08, +0.27] ; +0.27 [+0.18, +0.41] |
| 6 stations | 0.12/0.22 ; 0.23/0.27 | 0.29/0.44 ; 0.31/0.47 | 0.17/0.21 ; 0.19/0.27 | 0.21/0.24 ; 0.15/0.18 | +0.20 [+0.05, +0.31] ; +0.29 [+0.14, +0.40] |

U-a is met: WAVE-P exceeds WAVE-E by at least $0.05$ at both caps in five of six seed-cells. U-b is
met: it exceeds GaMMA by at least $0.05$ at both caps in four of six; in the one-sided cells it is
behind at $0.1$ in seed 2024 ($0.31$ against $0.37$) and ahead by $0.04$ in seed 7331, and ahead at
$0.25$ in both, with an interval that includes zero in seed 2024. The margins over GaMMA at four and six stations, $0.17$
to $0.29$ at $0.25$, are of the same size as in the known medium. A delay of $0.15$ s lies inside the
$\pm0.25$ s window of the migration; larger unmodelled delays, a wavelet different from the matched
filter, and scattered coda are not tested here.

## Y.10 Source mechanism and a geometry-breaking control (Task 1w)

Code: `task1w_mech.py`, `analyze_1w.py`; rules W-a to W-c written before the run. A second independent
verification of Y.6 to Y.9 found two things that could make those positives false. First, the
generator gave every S wave a pure SV polarization, the direction WAVE-P projects on; a real source
radiates S with an SH part set by its mechanism, and with random double-couple sources the advantage
of WAVE-P over WAVE-E nearly vanished in an exploratory run, while a mechanism-independent variant
kept a margin. Second, PHYS-SH was not a null (Y.6). This task changes the generator and the
controls together. Each source is a random double couple: P radiates along the ray with amplitude
$n\cdot Mn$ and S radiates the part of $Mn$ transverse to the ray, both normalized to unit mean
square over the sphere, and the radiation factor enters the waveform, the detection probability and
the pick timing error. About forty per cent of the true arrivals are detected. New methods:
PHYS-NODE is PHYS-LR with the predicted directions taken from another grid node at least 9 km away (a
fixed random map per window), which keeps every degree of freedom and removes the geometry of the
hypothesis; WAVE-PT is the mechanism-independent migration, P projected on the ray and S as the
energy in the plane transverse to the ray, the only polarization constraints that hold for any
source in an isotropic medium; WAVE-PT-NODE takes the ray from another node. The ratio $f/g$ of
PHYS-LR is re-estimated on forty windows of this generator. Sixty calibration and 120 evaluation
windows per cell; intervals re-select the operating point in every replicate, including GaMMA's
configuration.

**Table Y.11. Double-couple sources. Fraction of true arrivals detected, and event recall at realized
false events per window $\le0.1$ / $\le0.25$.**

| Cell | detected | KIN-LR | PHYS-LR | PHYS-NODE | WAVE-E | WAVE-P | WAVE-PT | WAVE-PT-NODE | GaMMA |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| one-sided, 8 st. (2024) | 0.42 | 0.29/0.40 | 0.31/0.43 | 0.29/0.37 | 0.25/0.32 | 0.18/0.26 | 0.32/0.42 | 0.24/0.34 | 0.29/0.29 |
| 4 stations (2024) | 0.42 | 0.06/0.11 | 0.11/0.21 | 0.08/0.15 | 0.08/0.17 | 0.09/0.20 | 0.17/0.24 | 0.10/0.19 | 0.04/0.04 |
| 6 stations (2024) | 0.40 | 0.11/0.20 | 0.19/0.27 | 0.12/0.18 | 0.10/0.31 | 0.18/0.31 | 0.28/0.47 | 0.18/0.28 | 0.08/0.12 |
| one-sided, 8 st. (7331) | 0.43 | 0.34/0.43 | 0.33/0.44 | 0.30/0.44 | 0.29/0.33 | 0.22/0.31 | 0.37/0.47 | 0.34/0.40 | 0.29/0.30 |
| 4 stations (7331) | 0.39 | 0.09/0.13 | 0.12/0.16 | 0.12/0.19 | 0.11/0.18 | 0.16/0.22 | 0.33/0.39 | 0.14/0.26 | 0.08/0.11 |
| 6 stations (7331) | 0.43 | 0.10/0.23 | 0.18/0.30 | 0.14/0.17 | 0.10/0.29 | 0.20/0.27 | 0.36/0.41 | 0.17/0.28 | 0.17/0.23 |

**Table Y.12. Differences at realized false events per window $\le0.25$, paired bootstrap 95%
intervals with re-selected operating points.**

| Cell | PHYS-LR - PHYS-NODE | PHYS-LR - KIN-LR | WAVE-PT - WAVE-PT-NODE | WAVE-PT - WAVE-E | PHYS-LR - GaMMA | WAVE-PT - GaMMA |
|:--|:--|:--|:--|:--|:--|:--|
| one-sided, 8 st. (2024) | +0.07 [-0.01, +0.13] | +0.04 [-0.04, +0.12] | +0.08 [-0.02, +0.15] | +0.09 [+0.00, +0.16] | +0.14 [+0.03, +0.27] | +0.12 [+0.01, +0.23] |
| 4 stations (2024) | +0.06 [+0.00, +0.12] | +0.10 [+0.01, +0.16] | +0.05 [+0.01, +0.14] | +0.07 [+0.01, +0.15] | +0.17 [+0.08, +0.24] | +0.20 [+0.11, +0.29] |
| 6 stations (2024) | +0.09 [+0.03, +0.15] | +0.08 [+0.01, +0.16] | +0.19 [+0.07, +0.27] | +0.16 [+0.08, +0.27] | +0.15 [+0.05, +0.25] | +0.35 [+0.25, +0.48] |
| one-sided, 8 st. (7331) | +0.00 [-0.09, +0.10] | +0.01 [-0.07, +0.11] | +0.07 [-0.01, +0.12] | +0.14 [-0.01, +0.19] | +0.14 [+0.00, +0.29] | +0.17 [+0.06, +0.29] |
| 4 stations (7331) | -0.02 [-0.07, +0.08] | +0.04 [-0.03, +0.13] | +0.13 [+0.04, +0.22] | +0.21 [+0.12, +0.31] | +0.06 [-0.01, +0.15] | +0.28 [+0.18, +0.38] |
| 6 stations (7331) | +0.13 [+0.02, +0.23] | +0.08 [-0.03, +0.18] | +0.13 [+0.03, +0.21] | +0.12 [+0.06, +0.27] | +0.08 [-0.03, +0.23] | +0.18 [+0.10, +0.31] |

**Registered decisions.** W-a is **not** met: with realistic mechanisms the per-pick likelihood PHYS-LR
exceeds its geometry-breaking control by at least $0.05$ at both caps in one of three cells in seed
2024 and in none in seed 7331, and its gain over KIN-LR shrinks from $+0.12$ to $+0.24$ in Y.6 to
$+0.01$ to $+0.10$ at $0.25$. The per-pick gain of Y.6 therefore rested mostly on the SV polarization
of the generator: its S term predicts a direction a real source does not respect. W-c for PHYS-LR is not met
in seed 7331. W-b is met: WAVE-PT exceeds both its control and the non-directional migration by at
least $0.05$ at both caps in two of three cells in both seeds, with intervals excluding zero at four
and six stations and including it in both one-sided cells. W-c for WAVE-PT is met: it exceeds seeded,
tuned GaMMA by at least $0.05$ at both caps in two of three cells in seed 2024 and in all three in
seed 7331, by $0.12$ to $0.35$ at $0.25$; this comparison used the clutter of E26 and does not survive
consistent clutter (Y.12). WAVE-P, which projects S on SV, falls below WAVE-E in the
one-sided cells, as the verification found.

The first eighty evaluation windows of this task, and their mechanisms, coincide with those of the
exploratory run that proposed WAVE-PT, so these tables are not an independent confirmation of it. The
replication was registered after three of the six seed-cells above had been seen, and was run on 120
new windows per cell, with networks, sources, clutter, detections and mechanisms drawn from other random
streams, for the pick-free methods and GaMMA.

**Table Y.13. Replication on fresh windows, double-couple sources. Fraction of true arrivals detected;
event recall at realized false events per window $\le0.1$ / $\le0.25$; differences at $\le0.25$ with
paired bootstrap 95% intervals.**

| Cell | detected | WAVE-E | WAVE-P | WAVE-PT | WAVE-PT-NODE | GaMMA | WAVE-PT - WAVE-PT-NODE | WAVE-PT - WAVE-E | WAVE-PT - GaMMA |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| one-sided, 8 st. (2024) | 0.38 | 0.15/0.33 | 0.12/0.26 | 0.30/0.40 | 0.22/0.33 | 0.23/0.26 | +0.07 [-0.01, +0.16] | +0.07 [-0.01, +0.18] | +0.14 [+0.03, +0.28] |
| 4 stations (2024) | 0.39 | 0.05/0.17 | 0.11/0.18 | 0.16/0.24 | 0.05/0.15 | 0.06/0.06 | +0.10 [+0.01, +0.18] | +0.07 [-0.01, +0.18] | +0.18 [+0.08, +0.27] |
| 6 stations (2024) | 0.38 | 0.18/0.26 | 0.27/0.35 | 0.33/0.40 | 0.22/0.32 | 0.14/0.18 | +0.08 [+0.01, +0.21] | +0.14 [+0.02, +0.25] | +0.22 [+0.12, +0.36] |
| one-sided, 8 st. (7331) | 0.42 | 0.12/0.24 | 0.15/0.22 | 0.25/0.38 | 0.19/0.34 | 0.22/0.25 | +0.05 [-0.03, +0.12] | +0.14 [+0.03, +0.22] | +0.14 [+0.03, +0.24] |
| 4 stations (7331) | 0.41 | 0.11/0.18 | 0.09/0.24 | 0.25/0.34 | 0.08/0.16 | 0.12/0.12 | +0.18 [+0.07, +0.26] | +0.16 [+0.07, +0.26] | +0.22 [+0.13, +0.35] |
| 6 stations (7331) | 0.39 | 0.10/0.18 | 0.18/0.28 | 0.25/0.38 | 0.08/0.21 | 0.04/0.09 | +0.17 [+0.04, +0.27] | +0.21 [+0.06, +0.29] | +0.29 [+0.14, +0.43] |

Both replication rules are met. WAVE-PT exceeds its geometry-breaking control and the
non-directional migration by at least $0.05$ at both caps in all three cells of seed 2024 and in two
of seed 7331, and seeded, tuned GaMMA in the same cells, by $0.14$ to $0.29$ at $0.25$. The one-sided
cells remain the weakest: there the interval against the control includes zero in both seeds, as it
did in Y.8 and in Table Y.12.

## Y.11 Coherent moment-tensor detection with explaining-away (Task 1x)

Code: `task1x_mtglrt.py`, `analyze_1x.py`, `cal_capv.py`, `chk_nms.py`; rules X-a to X-c written before
the run. This task tests the reformulated method of Section 7.5. For a hypothesis $\theta=(\mathbf x,\tau)$
the signed three-component vectors of the matched-filtered records at the predicted P and S times of
every station are linear in the moment tensor, $\mathbf v=G(\theta)\mathbf m+\mathbf n$, so the
generalized likelihood ratio with $\mathbf m$ as a linear nuisance is $\|P_\theta\mathbf v\|^2$,
the energy of the data in the six-dimensional subspace that the radiation pattern allows, jointly over
stations and phases (WAVE-MT). It demands consistency of signs and relative amplitudes across
stations, which WAVE-PT ignores. Without the search window and the cap it is the squared dual
certificate of the group BLASSO over tensor-valued measures (Proposition AA.2), so its global maximum on
the residual is the Frank-Wolfe oracle (Proposition AA.4(a)). Taking that maximum, fitting $\mathbf m$ by
unpenalized least squares, subtracting the fitted contribution and repeating, with at most five atoms
and with every node within 5 km epicentral distance of an accepted atom excluded from reselection, is
MT-FW: a greedy group matching pursuit whose selection step is the Frank-Wolfe oracle. It is not
Frank-Wolfe on the group BLASSO, because earlier atoms are never refitted, the amplitude is not shrunk
and the statistic is capped (Proposition AA.4(e)-(h)). It explains away an event instead of suppressing
its neighbours by distance. Controls: WAVE-MT-NODE uses the subspace of
another node at least 9 km away. Conditions: ideal (double-couple sources, known homogeneous medium)
and pert (station delays $\mathcal N(0,0.15\ \mathrm s)$, about ten degrees of polarization deviation
and lognormal site amplitude factors with $\sigma=0.5$, unknown to every method). The per-station cap
on the vector norm was chosen on separate calibration windows. One hundred fresh windows per cell.

**Table Y.14. Event recall at realized false events per window $\le0.1$ / $\le0.25$, ideal
condition.** The GaMMA column is not a valid comparison: GaMMA received undetectable clutter picks
(E26); Y.12 repeats it with consistent clutter.

| Cell | detected | WAVE-E | WAVE-PT | WAVE-MT | WAVE-MT-NODE | MT-FW | GAMMA |
|:--|:--|:--|:--|:--|:--|:--|:--|
| one-sided, 8 st. (2024) | 0.38 | 0.14/0.26 | 0.23/0.36 | 0.45/0.54 | 0.37/0.45 | 0.78/0.82 | 0.18/0.23 |
| 4 stations (2024) | 0.39 | 0.11/0.24 | 0.24/0.34 | 0.28/0.34 | 0.11/0.20 | 0.46/0.64 | 0.09/0.16 |
| 6 stations (2024) | 0.40 | 0.15/0.27 | 0.25/0.34 | 0.27/0.44 | 0.13/0.32 | 0.73/0.78 | 0.08/0.11 |
| one-sided, 8 st. (7331) | 0.40 | 0.25/0.35 | 0.34/0.44 | 0.47/0.59 | 0.34/0.44 | 0.87/0.89 | 0.00/0.28 |
| 4 stations (7331) | 0.38 | 0.09/0.16 | 0.10/0.19 | 0.22/0.35 | 0.11/0.23 | 0.57/0.64 | 0.07/0.08 |
| 6 stations (7331) | 0.39 | 0.13/0.25 | 0.24/0.36 | 0.41/0.60 | 0.12/0.24 | 0.83/0.87 | 0.05/0.14 |

**Table Y.15. The same, pert condition.**

| Cell | detected | WAVE-E | WAVE-PT | WAVE-MT | WAVE-MT-NODE | MT-FW | GAMMA |
|:--|:--|:--|:--|:--|:--|:--|:--|
| one-sided, 8 st. (2024) | 0.40 | 0.18/0.31 | 0.28/0.40 | 0.49/0.60 | 0.33/0.53 | 0.74/0.79 | 0.23/0.27 |
| 4 stations (2024) | 0.38 | 0.08/0.15 | 0.15/0.27 | 0.16/0.42 | 0.04/0.12 | 0.46/0.57 | 0.05/0.08 |
| 6 stations (2024) | 0.42 | 0.22/0.29 | 0.27/0.42 | 0.38/0.44 | 0.22/0.37 | 0.47/0.58 | 0.11/0.14 |
| one-sided, 8 st. (7331) | 0.42 | 0.25/0.38 | 0.33/0.46 | 0.41/0.54 | 0.24/0.42 | 0.65/0.75 | 0.22/0.25 |
| 4 stations (7331) | 0.38 | 0.07/0.15 | 0.15/0.18 | 0.24/0.30 | 0.05/0.16 | 0.43/0.53 | 0.06/0.07 |
| 6 stations (7331) | 0.39 | 0.13/0.23 | 0.22/0.33 | 0.30/0.53 | 0.18/0.24 | 0.65/0.72 | 0.05/0.06 |

**Registered decisions.** X-a is **not** met in either condition: WAVE-MT exceeds its geometry-breaking
control in most cells, but not WAVE-PT by at least $0.05$ at both caps in two of three cells of seed
2024. The single-pass coherent statistic is therefore not uniformly better than the
mechanism-independent migration. X-b (WAVE-MT against GaMMA) was met but is invalidated by E26. X-c
is met: MT-FW exceeds WAVE-MT by at least $0.05$ at $0.25$ false events per window in every cell of
both seeds and conditions.

**How much of the deflation gain is sidelobe suppression.** The single-pass methods emit a
catalogue by suppressing neighbours within 5 km, and the sidelobes of a detected event beyond that
radius count as false events. With a 20 km radius (which knows, as a method should not, that sources
are at least 12 km apart), WAVE-MT reaches $0.72/0.77$ against $0.78/0.82$ for MT-FW in the one-sided
cell and $0.62/0.72$ against $0.73/0.78$ at six stations (ideal, seed 2024). Most of the gain of
deflation over the 5 km catalogue is therefore the removal of sidelobes, which deflation does with a
fixed 5 km exclusion, set before evaluation and equal to the radius of that catalogue, instead of the
20 km suppression the single pass needs, which would merge events closer than 20 km.

**Robustness found by the independent verification (exploratory, 50 windows, one-sided cell).** A
wavelet rotated by 90 degrees in phase, or a 1.5 Hz Ricker pulse against the matched filter, left MT-FW
at $0.81/0.86$ and $0.79/0.81$. Unmodelled station delays beyond the $\pm0.25$ s window do break it:
with delays of standard deviation $0.4$ s MT-FW fell to $0.37/0.42$, below GaMMA with consistent
clutter ($0.56$). The method therefore needs station time corrections calibrated on catalogued events,
as every regional location does, and the residual delay must stay within the window.

## Y.12 Consistent clutter (Task 1y)

Code: `task1y_clutter.py`, `analyze_1y.py`; rules Y-a to Y-c written before the run. The verification
of Y.11 found E26: since Y.5 the generator applied the detection rule to true arrivals but not to
clutter, whose pulses lie below the noise yet were all delivered to GaMMA with exact times; about 90
per cent of GaMMA's input was undetectable clutter. Here clutter follows the same detection and
timing rule as true arrivals, in two consistent regimes: weak clutter (amplitude $0.35$, almost never
picked) and strong clutter (amplitude $1.75$, visible transients that are picked, as spikes and
cultural noise are). Pert condition of Y.11; 100 fresh windows per cell; single-pass methods reported
with 5 and 20 km suppression; MT-FW-NODE is the same greedy iteration with the subspace of another
node, the control of the complete method.

**Table Y.16. Weak clutter. Clutter and true picks per window, and event recall at realized false
events per window $\le0.1$ / $\le0.25$.**

| Cell | picks clutter/true | WAVE-PT@20 | WAVE-MT@20 | MT-FW | MT-FW-NODE | GAMMA |
|:--|:--|:--|:--|:--|:--|:--|
| one-sided, 8 st. (2024) | 4.9/5.0 | 0.56/0.69 | 0.65/0.68 | 0.73/0.77 | 0.51/0.61 | 0.63/0.63 |
| 4 stations (2024) | 2.2/2.6 | 0.28/0.34 | 0.36/0.45 | 0.38/0.54 | 0.14/0.19 | 0.28/0.28 |
| 6 stations (2024) | 3.7/3.6 | 0.38/0.50 | 0.49/0.57 | 0.49/0.65 | 0.18/0.24 | 0.54/0.54 |
| one-sided, 8 st. (7331) | 4.8/4.7 | 0.61/0.66 | 0.61/0.65 | 0.66/0.71 | 0.44/0.57 | 0.55/0.55 |
| 4 stations (7331) | 2.2/2.7 | 0.31/0.47 | 0.48/0.51 | 0.49/0.63 | 0.20/0.41 | 0.35/0.35 |
| 6 stations (7331) | 3.5/4.1 | 0.51/0.58 | 0.60/0.64 | 0.67/0.76 | 0.24/0.46 | 0.54/0.54 |

**Table Y.17. Weak clutter. Differences at $\le0.25$, paired bootstrap 95% intervals with re-selected
operating points and GaMMA configuration.**

| Cell | MT-FW - GAMMA | MT-FW - MT-FW-NODE | MT-FW - WAVE-MT@20 | WAVE-MT@20 - GAMMA | WAVE-PT@20 - GAMMA |
|:--|:--|:--|:--|:--|:--|
| one-sided, 8 st. (2024) | +0.15 [+0.05, +0.25] | +0.16 [+0.06, +0.35] | +0.09 [+0.02, +0.18] | +0.05 [-0.03, +0.15] | +0.07 [-0.05, +0.20] |
| 4 stations (2024) | +0.26 [+0.07, +0.43] | +0.35 [+0.13, +0.53] | +0.09 [-0.03, +0.19] | +0.17 [+0.02, +0.33] | +0.06 [-0.05, +0.20] |
| 6 stations (2024) | +0.11 [-0.05, +0.24] | +0.41 [+0.17, +0.54] | +0.08 [-0.01, +0.16] | +0.03 [-0.11, +0.16] | -0.04 [-0.16, +0.12] |
| one-sided, 8 st. (7331) | +0.16 [+0.07, +0.25] | +0.13 [+0.01, +0.30] | +0.06 [+0.00, +0.12] | +0.10 [+0.01, +0.20] | +0.11 [+0.01, +0.24] |
| 4 stations (7331) | +0.28 [+0.14, +0.40] | +0.22 [+0.06, +0.42] | +0.12 [+0.01, +0.19] | +0.16 [+0.05, +0.31] | +0.12 [-0.02, +0.27] |
| 6 stations (7331) | +0.22 [+0.06, +0.32] | +0.30 [+0.13, +0.45] | +0.12 [+0.02, +0.18] | +0.10 [-0.01, +0.21] | +0.04 [-0.08, +0.14] |

**Table Y.18. Strong clutter. Clutter and true picks per window, and event recall.**

| Cell | picks clutter/true | WAVE-PT@20 | WAVE-MT@20 | MT-FW | MT-FW-NODE | GAMMA |
|:--|:--|:--|:--|:--|:--|:--|
| one-sided, 8 st. (2024) | 34.2/4.9 | 0.44/0.52 | 0.63/0.65 | 0.61/0.69 | 0.27/0.44 | 0.33/0.33 |
| 4 stations (2024) | 28.8/2.4 | 0.01/0.15 | 0.15/0.26 | 0.15/0.24 | 0.04/0.07 | 0.05/0.06 |
| 6 stations (2024) | 29.2/3.5 | 0.35/0.38 | 0.51/0.59 | 0.47/0.55 | 0.15/0.20 | 0.19/0.19 |
| one-sided, 8 st. (7331) | 31.1/4.8 | 0.40/0.50 | 0.51/0.55 | 0.56/0.61 | 0.30/0.43 | 0.24/0.24 |
| 4 stations (7331) | 29.8/2.7 | 0.06/0.20 | 0.28/0.36 | 0.30/0.37 | 0.07/0.15 | 0.04/0.05 |
| 6 stations (7331) | 29.7/4.1 | 0.34/0.41 | 0.41/0.55 | 0.43/0.61 | 0.28/0.35 | 0.12/0.18 |

**Table Y.19. Strong clutter. Differences at $\le0.25$.**

| Cell | MT-FW - GAMMA | MT-FW - MT-FW-NODE | MT-FW - WAVE-MT@20 | WAVE-MT@20 - GAMMA | WAVE-PT@20 - GAMMA |
|:--|:--|:--|:--|:--|:--|
| one-sided, 8 st. (2024) | +0.36 [+0.22, +0.49] | +0.25 [+0.11, +0.39] | +0.04 [-0.04, +0.11] | +0.32 [+0.20, +0.45] | +0.19 [+0.07, +0.30] |
| 4 stations (2024) | +0.17 [+0.07, +0.33] | +0.16 [+0.05, +0.28] | -0.03 [-0.07, +0.05] | +0.20 [+0.07, +0.33] | +0.09 [-0.04, +0.18] |
| 6 stations (2024) | +0.36 [+0.24, +0.55] | +0.35 [+0.22, +0.52] | -0.04 [-0.09, +0.09] | +0.41 [+0.25, +0.53] | +0.19 [+0.08, +0.34] |
| one-sided, 8 st. (7331) | +0.37 [+0.23, +0.49] | +0.18 [+0.04, +0.31] | +0.06 [+0.00, +0.14] | +0.30 [+0.19, +0.41] | +0.26 [+0.14, +0.43] |
| 4 stations (7331) | +0.32 [+0.19, +0.46] | +0.22 [+0.10, +0.36] | +0.01 [-0.04, +0.11] | +0.31 [+0.17, +0.43] | +0.15 [+0.04, +0.27] |
| 6 stations (7331) | +0.43 [+0.26, +0.57] | +0.27 [+0.12, +0.38] | +0.06 [-0.04, +0.13] | +0.37 [+0.24, +0.50] | +0.23 [+0.12, +0.36] |

**Registered decisions.** Y-a is met in both regimes: MT-FW exceeds seeded, tuned GaMMA by at least
$0.05$ at both caps in at least two of three cells per seed, by $0.11$ to $0.28$ at $0.25$ with weak
clutter (in the six-station cell of seed 2024 it is behind at $0.1$) and by $0.17$ to $0.43$ with strong
clutter. Y-b is met in both: the same iteration with the geometry of another node loses $0.13$ to
$0.41$. Y-c is met with weak clutter and not with strong clutter: there deflation adds little over a
single pass with a 20 km suppression radius, and the gain over GaMMA comes from the coherence of the
moment-tensor subspace, which rejects visible transients that GaMMA must accept as picks.

**What E26 changes in Y.5 to Y.10.** With consistent weak clutter GaMMA reaches $0.54$ to $0.63$ in
the one-sided and six-station cells, not the $0.1$ to $0.4$ of the earlier tables, and the
mechanism-independent migration with a 20 km suppression (WAVE-PT@20) does not beat it (differences
$-0.04$ to $+0.12$ at $0.25$, intervals including zero in five of six cells). The margins over GaMMA
reported in Y.5 to Y.10 are inflated by E26 and are superseded by this section. The comparisons of
the waveform methods with their controls in Y.8 to Y.10 do not use picks and are unaffected.

## Y.13 What this changes

1. **Unchanged:** the one-sided advantage of the prototype over the coincidence-count ablation (X.3)
   was a tie effect (Y.1), and the coincidence count is retired as a comparator. With every arrival
   picked exactly, a continuous kinematic score suffices (Y.2, Y.3, R-c).
2. **Corrected, three times.** The first reading of Y.5 overstated the effect (nominal budget, weak
   reference, no control). The per-pick positive of Y.6 rested on pure-SV S waves and on a control
   that was not null (E24, E25); with realistic mechanisms the per-pick likelihood does not beat a
   geometry-breaking control (Y.10). And every margin over GaMMA from Y.5 to Y.10 was inflated by
   clutter picks that no picker would produce (E26); with consistent clutter the mechanism-independent
   migration WAVE-PT no longer beats GaMMA (Y.12).
3. **Supported on synthetics, with geometry-breaking controls, consistent clutter, two seeds and fresh
   windows:** the coherent moment-tensor detector with explaining-away, MT-FW, a greedy group matching
   pursuit whose selection step is the Frank-Wolfe oracle of the group BLASSO of Section 7.5 (it does
   not solve that problem; Proposition AA.4), exceeds seeded, tuned GaMMA by $0.11$ to
   $0.28$ in event recall at $0.25$ false events per window when clutter is weak and by $0.17$ to
   $0.43$ when clutter is made of visible, picked transients, with station delays, polarization
   deviations and site amplitude factors unknown to the method (Y.12). The same iteration with the
   geometry of another node loses $0.13$ to $0.41$.
4. **Decomposed:** the single-pass coherent statistic is not uniformly better than the
   mechanism-independent migration (X-a not met); what separates MT-FW from GaMMA is the combination of
   coherence, which rejects transients that do not fit any moment tensor, and explaining-away, which
   removes sidelobes and overlapping events with a fixed 5 km exclusion instead of a 20 km suppression
   and pays mostly when clutter is weak.
5. **Not supported:** polarization as a filter on GaMMA's catalogue (Y.7); an S projection on SV
   (Y.10); the per-pick likelihood with realistic mechanisms (Y.10); the mechanism-independent
   migration against GaMMA once clutter is consistent (Y.12).
6. **Fragile:** unmodelled station delays beyond the $\pm0.25$ s window (Y.11), which makes station
   time corrections calibrated on catalogued events part of the method.
7. **Registered at this point and since run:** a calibrated real picker and QuakeMigrate [109], run
   as published against MT-FW on a generator with scattered coda that the detectors do not know, in
   the preregistered Task 2A (Y.14, Y.15, Y.18); and, added later as the preregistered Task 2C, block
   orthogonal matching pursuit in the form of [129] and mechanism-corrected stacking in the form of
   [122] on the shared grid (Y.20), and, as the preregistered Task 2D, the same two with the peak
   search of MT-FW (Y.21); the off-grid group BLASSO with time shifts per station and phase, as the
   preregistered Task 2E (Y.22); the per-station cap in dense sequences, as the preregistered Task 2F
   after an exploration on the calibration seed (Y.23, Y.24); a Huber fit per station and phase with a
   level fixed by rule from the noise, as the preregistered Task 2I after a second exploration (Y.25,
   Y.26); station terms from a joint determination without picks inside the refined detector, as the
   preregistered Task 2L after a third exploration (Y.28); and the certified bound for the certificate in the pruned search, which is
   valid but avoids only about 4 per cent of the exact evaluations (Appendix Z). **Registered at this
   point and not yet run:** MALMI [110] and the polarity-corrected stacking of [112], run as published
   against MT-FW; synthetic waveforms with real noise (Task 1p); 3D Green's functions (Task 1v); MT-FW
   in the real-data decimation test (Section 16); and the threshold calibration of Appendix Y.17
   repeated on real noise-only windows.

**Declared limitations.** Synthetic generators in a homogeneous medium, perturbed by station delays
within the migration window, small direction deviations, source mechanisms and site amplitude
factors; a matched filter built from a wavelet close to the true one (a phase rotation and a
different corner frequency were tolerated in an exploratory check); detection and timing errors that
are parametric rather than produced by a real picker; no amplitudes or pick probabilities for the
associators; operating points chosen with an oracle on the evaluation windows for every method, read
on a FROC curve; GaMMA and PyOcto tuned over small grids; two seeds; rules T, U, W, X and Y and the
replications written after parts of the data they concern had been seen, and the exact-curve
amendment after a registered rule had failed on the coarse sweep.

## Y.14 A generator the detectors do not know, and real pickers

Code: `gen_realistic.py` (generator and self-test), `real_picker.py` (`--mode calib`, `--mode
report`); output `gen_realistic_report.txt`. Every task from Y.5 to Y.12 synthesized waveforms with
the propagation physics the detectors assume (homogeneous medium, straight rays, one 2 Hz Ricker
wavelet), perturbed at most by random station delays, direction deviations and site factors, and
fed the pick-based associators with picks drawn from the parametric rule of Y.5, called `degrade_y`
here (detection probability $\Phi\big((a-2)/0.5\big)$, timing error three times the Cramer-Rao
bound), not with a picker run on the waveforms. This component study removes both shortcuts. It uses
seed 999 only; the evaluation seeds 2024 and 7331 were not touched, and no decision rule was
registered for it.

**Generator.** `generate_real` keeps the interface and the random streams of the generator of Y.12
and changes the physics. (a) The medium is a linear gradient $v(z)=v_0+gz$ with $V_{S0}=2.8$ km/s,
$g=0.04$ s$^{-1}$ and $V_P=1.73\,V_S$ at every depth ($V_{P0}=4.844$ km/s). Travel times follow the
closed form $t=g^{-1}\operatorname{arccosh}\big(1+g^2r^2/(2v(z_s)v(z_r))\big)$, evaluated in an
asinh form that is stable as $g\to0$, and rays are circular arcs. P is polarized along the curved
ray at the station; S is the transverse part of $M\mathbf n_s$ at the source, carried to the station
by the rotation about the SH direction, which is the parallel transport of the polarization frame in
a 1D isotropic medium. (b) Station delays $\mathcal N(0,0.15\ \mathrm s)$, direction deviations of
about ten degrees and lognormal site factors with $\sigma=0.5$ are drawn as in the pert condition.
(c) Each event has its own source wavelet, a Ricker pulse with centre frequency uniform in
$[1.4,2.8]$ Hz rotated in phase by up to $\pm45$ degrees, while the matched filter of the detectors
stays at 2 Hz. (d) Every true P and S arrival is followed by a scattered coda: three independent
components, exponential envelope with a 3 s decay, causal 1-4 Hz band, and 0.3 of the energy of the
direct arrival (0.301 measured in a noise-free window). (e) Clutter and noise are those of Y.12. With
$g=0$ and $V_S=3$ km/s every pick time and amplitude of the generator of Y.12 is reproduced
(differences 0 and $3\times10^{-16}$), so only the waveforms change.

**Verification.** As $g\to0$ the relative travel-time error against $r/v_0$ is at most
$2.3\times10^{-7}$ and the ray directions coincide with the straight ray to
$1-\cos\le3.5\times10^{-13}$; the asinh and arccosh forms agree to $1.9\times10^{-15}$. Numerical ray
shooting (the kinematic ray equations integrated with DOP853 at tolerance $10^{-12}$) for three
source-station pairs and both phases reproduces the analytic times to $2.0\times10^{-13}$ relative
and the take-off and arrival directions to $1.2\times10^{-6}$ degrees. An independent second-order
fast-marching eikonal solver (scikit-fmm) on the same pairs agrees to within $8.2\times10^{-5}$
relative at a 0.04 km step and $4.1\times10^{-5}$ at 0.02 km, converging with the step. The S
transport conserves the norm to $5.6\times10^{-16}$ and keeps S transverse to the arriving ray to
$1.7\times10^{-16}$.

**Mismatch against the detectors.** The detectors keep the homogeneous medium ($V_P=5.19$,
$V_S=3.0$ km/s) and straight rays. Over 2400 source-station pairs of the benchmark geometry, the
gradient travel time minus the homogeneous one has median $-0.007$ s for P (5th to 95th percentile
$-0.22$ to $+0.14$ s) and $-0.012$ s for S ($-0.38$ to $+0.24$ s, extremes $-0.76$ and $+0.30$ s):
near stations are late, because $v_0<V_S$, and far stations early, because the ray dives into faster
rock. The 5th percentile of the S error lies outside the $\pm0.25$ s window of the migration. The
arrival direction differs from the straight ray by a median of 12.4 degrees (95th percentile 19.0,
maximum 22.5), comparable to the direction deviations of the pert condition.

**Pickers.** Two real pickers run on the waveforms. PhaseNet [11], with the original pretrained
model of the authors' repository (model 190703-214543, not retrained) and the preprocessing of its
own continuous-data prediction script (components mapped to E, N and Z, resampling from 16 to
100 Hz, its normalization), picks the probability maxima above a threshold that are at least 0.5 s
apart. A classical STA/LTA (ObsPy `classic_sta_lta` and `trigger_onset`, zero-phase Butterworth
band-pass) triggers on the vertical component for P and on the horizontal norm for S. Both were
calibrated on seed-999 windows disjoint from those reported below (80 windows, 840 arrivals), by the
F1 score of arrivals matched one to one within 0.3 s, per phase, with an onset-to-centre offset
estimated in calibration, since pickers mark the onset and the generator times the centre of the
pulse. The calibrated PhaseNet thresholds are 0.25 for P and 0.03 for S; the calibrated STA/LTA uses
the 1-4 Hz band, STA 0.3 s, LTA 8 s and trigger 6.0 for P, and STA 0.6 s, LTA 8 s and trigger 3.0 for
S. Onset refinement by AIC or Baer-Kradolfer and P/S disambiguation by vertical energy lowered F1 in
calibration and are not used. Calibration F1 was 0.025 (P) and 0.135 (S) for PhaseNet and 0.126 and
0.193 for STA/LTA. The thresholds lie inside their grids; STA and LTA lie at the edge of grids of two
values each.

**Table Y.20. Held-out pick statistics: seed 999, 40 windows per cell (160 windows, 116 events, 1608
arrivals), 35 clutter pulses per window.** Picked: fraction of true arrivals matched one to one by a
pick of the same phase within 0.3 s. Timing error: $|t_{\mathrm{pick}}-t_{\mathrm{true}}|$ of those
picks, truncated at 0.3 s by the definition. False picks per window split into swapped phase (within
0.3 s of an arrival of the other phase), clutter (within 0.3 s of a clutter pulse) and other (noise,
coda, misplaced arrivals). `degrade_y` is applied to the same arrivals and amplitudes. Script
`real_picker.py --mode report`; output `gen_realistic_report.txt`, section 4.

| Cell | method | picked (P / S) | error median / p90 (s) | false per window (phase / clutter / other) |
|:--|:--|:--|:--|:--|
| 6 stations, weak clutter | PhaseNet | 0.032 (0.005 / 0.059) | 0.080 / 0.187 | 1.90 (0.30 / 0.05 / 1.55) |
| | STA/LTA | 0.145 (0.059 / 0.231) | 0.087 / 0.220 | 3.20 (0.93 / 0.17 / 2.10) |
| | degrade_y | 0.435 (0.312 / 0.559) | 0.031 / 0.098 | 0.07 (0.00 / 0.05 / 0.03) |
| 6 stations, strong clutter | PhaseNet | 0.035 (0.005 / 0.065) | 0.083 / 0.173 | 2.33 (0.25 / 0.53 / 1.55) |
| | STA/LTA | 0.156 (0.070 / 0.242) | 0.077 / 0.222 | 9.25 (0.85 / 6.85 / 1.55) |
| | degrade_y | 0.422 (0.290 / 0.554) | 0.028 / 0.125 | 27.98 (0.35 / 27.62 / 0.00) |
| one-sided, 8 st., weak clutter | PhaseNet | 0.051 (0.005 / 0.097) | 0.074 / 0.193 | 1.95 (0.28 / 0.05 / 1.62) |
| | STA/LTA | 0.134 (0.032 / 0.236) | 0.060 / 0.227 | 3.90 (1.43 / 0.20 / 2.27) |
| | degrade_y | 0.394 (0.255 / 0.532) | 0.031 / 0.088 | 0.05 (0.00 / 0.03 / 0.03) |
| one-sided, 8 st., strong clutter | PhaseNet | 0.044 (0.000 / 0.088) | 0.082 / 0.230 | 2.12 (0.33 / 0.38 / 1.43) |
| | STA/LTA | 0.139 (0.051 / 0.227) | 0.071 / 0.220 | 11.38 (1.48 / 8.12 / 1.77) |
| | degrade_y | 0.433 (0.315 / 0.551) | 0.041 / 0.136 | 28.98 (0.35 / 28.62 / 0.00) |
| all four cells | PhaseNet | 0.041 (0.004 / 0.078) | 0.080 / 0.198 | 2.08 (0.29 / 0.25 / 1.54) |
| | STA/LTA | 0.143 (0.052 / 0.234) | 0.078 / 0.221 | 6.93 (1.17 / 3.84 / 1.93) |
| | degrade_y | 0.420 (0.292 / 0.549) | 0.033 / 0.113 | 14.27 (0.17 / 14.08 / 0.01) |

**Table Y.21. Fraction of true arrivals picked, by amplitude over noise $a$ (all cells; number of
arrivals in brackets).** Same script and output. The non-zero fractions of `degrade_y` below $a=2$
are coincidences with clutter picks or with picks of another source.

| Method | 0-1 (524) | 1-2 (450) | 2-3 (256) | 3-5 (232) | 5-10 (128) | >10 (18) |
|:--|:--|:--|:--|:--|:--|:--|
| PhaseNet | 0.00 | 0.00 | 0.01 | 0.03 | 0.31 | 0.78 |
| STA/LTA | 0.02 | 0.02 | 0.07 | 0.36 | 0.77 | 0.61 |
| degrade_y | 0.02 | 0.21 | 0.77 | 0.99 | 1.00 | 1.00 |

Over all cells, the median bias of the matched picks is $-0.057$ s for PhaseNet, $+0.008$ s for
STA/LTA and $-0.002$ s for `degrade_y`, and precision / F1 are 0.166 / 0.066, 0.172 / 0.156 and
0.228 / 0.296. On the calibration windows (in sample) the picked fractions were 0.061, 0.148 and
0.404 and the false picks per window 2.79, 7.49 and 14.41, of the same order as on the held-out
windows.

**What this says about the pick model of Tasks 1o to 1y.** `degrade_y` picks 42 per cent of the
arrivals, against 14 per cent for the calibrated STA/LTA and 4 per cent for PhaseNet, and its median
timing error, 0.033 s, is about 2.4 times smaller than theirs (0.078 and 0.080 s; 90th percentiles
0.113 against 0.221 and 0.198 s). It detects half of the arrivals at $a\approx2$, where the STA/LTA
needs $a\approx5$ and PhaseNet $a\approx8$ to 10. The gap is largest for P: `degrade_y` picks 29 per
cent of the P arrivals, the STA/LTA 5 per cent and PhaseNet 0.4 per cent, because S carries 1.5 times
the P amplitude plus the coda. With weak clutter `degrade_y` also produces almost no false picks
(0.05 to 0.07 per window), where the real pickers produce 1.90 to 3.90, mostly from noise and coda
and from swapped phases. In thin clutter it was therefore optimistic for GaMMA on every count. With
strong clutter the sign is mixed. `degrade_y` then delivers about 28 clutter picks per window (27.62
and 28.62), while the STA/LTA, with its high P trigger, fires on 6.85 and 8.12 of the roughly 28
visible pulses and PhaseNet on 0.53 and 0.38, because an isolated Ricker pulse does not resemble an
earthquake to either. In that regime `degrade_y` was optimistic for GaMMA in recall and timing and
pessimistic in clutter, so it was not uniformly optimistic, and the net effect on the strong-clutter
margins of Y.12 cannot be read from these statistics; it requires GaMMA on real picks, which Task 2A
runs (Y.18).

**PhaseNet here is out of its domain.** The pipeline was checked: on labelled samples of its own
repository the probability maxima fall within 0 to 10 samples of the labels, and on noise-free
windows of this generator the network gives 0.8 to 0.9 at almost every P and S arrival at the correct
time. With the generator's noise it barely responds: a P arrival at amplitude-to-noise 6.7 reaches a
probability of 0.09, and 0.43 with the noise reduced fourfold. The noise is white and, after
resampling, fills 0 to 8 Hz above the 1.4-2.8 Hz band of the signal, the opposite of the red,
microseism-dominated noise of its training data; none of the prefilters tried (band-pass 1-4, 1-6 and
0.7-3.5 Hz, low-pass 4 Hz) improved F1 in calibration. The failure of PhaseNet is therefore mostly
domain mismatch, isolated Ricker pulses in white noise sampled at 16 Hz, and says little about
PhaseNet on real data. For the synthetic benchmark the calibrated STA/LTA is the more informative
real picker, and PhaseNet a pessimistic bound and a check that conclusions do not hinge on one picker.

**Declared limitations.** White noise without microseism and a signal limited to 8 Hz by the 16 Hz
sampling, which is the main cause of the PhaseNet result; no attenuation, the same source function at
every station and for P and S; geometric spreading $1/r$ with the straight distance; no free
surface, converted, reflected or head waves; station delays uncorrelated with geometry; coda
isotropic in polarization with a simple exponential envelope; clutter made of isolated 2 Hz Ricker
pulses. The PhaseNet P offset rests on five residuals. Calibration and report both use seed 999, on
disjoint windows, and the statistics of Table Y.20 carry no confidence intervals.

## Y.15 QuakeMigrate as a comparator

Code: `qm_wrapper.py` (`calibrate`, `test`, `locate_check`); output `qm_report.txt`. QuakeMigrate
[109], version 1.2.2, installed from its published wheel with the compiled migration library and
not modified, is run through its public interface on each synthetic window: the three components of
every station written as a MiniSEED archive, a homogeneous lookup table ($V_P=5.19$, $V_S=3.0$ km/s)
on a grid of $44\times44\times21$ nodes at 1 km (40,656 nodes, depth 0 to 20 km), `QuakeScan.detect`
with its STA/LTA onset function at 50 Hz after the package's own resampling, and `Trigger.trigger`
with a static threshold and an origin-time window of 3 s, the length of the origin-time window of the
grid methods, shifted by the onset lag measured in calibration. The event location is the grid
maximum at the triggered peak; `locate` was checked on five windows and moved the epicentre by 0.2
to 4.5 km without changing any hit or miss, while needing a marginal window of about 1 s that the
calibrated configuration does not have. The event list is scored exactly as the greedy records of the
grid methods (5 km epicentral separation and tolerance, at most five events per window) and read on
the exact FROC curve. This comparison is **exploratory**: one seed, two cells and 47 test sources.
The confirmatory comparison is Task 2A (Y.18), in which QuakeMigrate, recalibrated on seed 999 by
this procedure, was the best comparator in every cell and MT-FW with station corrections remained
ahead of it.

**Channel practice.** QuakeMigrate's default practice computes the P onset on the vertical and the S
onset on the horizontals. On this generator that practice is physically inverted: in a homogeneous
full space with straight rays and epicentral distance two to six times the depth, the rays arrive 60
to 80 degrees from the vertical, so P energy lands mostly on the horizontals and SV energy mostly on
the vertical. In a clutter-free window at five times the nominal signal-to-noise ratio the vertical
"P" onsets peaked at the S times and the horizontal "S" onsets at the P times, with lookup-table and
data timing exact to 0.01 s. In a positive control with one source per window, no clutter and eight
windows per noise level, the default mapping found 3 of 8 sources at one fifth of the nominal noise
and 4 of 8 at nominal noise (median epicentral error of the top event 0.8 and 15.4 km), and the
three-component mapping, in which both onsets are computed from all three channels through the
package's documented `channel_maps` option, found 8 of 8 and 7 of 8 (0.5 and 0.5 km). Both are
reported: QM with the default practice, and QM-3C, the fair comparator for this generator. With real
data, near-vertical incidence makes the default practice the right one.

**Calibration.** Seed-999 stream of Task 1y, pert condition, weak clutter, in two cells, 6 stations
(35 clutter pulses) and one-sided with 8 stations (40 clutter pulses). Windows 30 to 59 of each cell
calibrate (60 windows, 47 sources) and windows 0 to 29 test. The onset grid covers three band-passes
(1-4, 1-6 and 2-16 Hz, the last the package default) and eight STA/LTA pairs (0.2/1.0 s, the package
default, to 1.2/6 s); the trigger grid covers the threshold quantile, marginal window, minimum event
interval, normalization and smoothing. The criterion is recall at at most 0.25 realized false events
per window, pooled over both cells. The selected QM onset is 1-4 Hz with STA/LTA 0.4/6.0 s
(calibration recall 0.09 / 0.09 / 0.11 at 0.05 / 0.1 / 0.25 false events per window) and the selected
QM-3C onset 1-6 Hz with 0.8/6.0 s (0.23 / 0.38 / 0.38); the package-default onset reached 0.00 /
0.00 / 0.00 and 0.00 / 0.04 / 0.06. The longest LTA of the grid, 6 s, was preferred in both variants,
so the optimum may lie beyond the edge of the grid, which the 40 s windows limit.

**Table Y.22. Event recall at realized false events per window $\le0.05$ / $\le0.1$ / $\le0.25$,
seed 999, 30 test windows per cell (not used in calibration).** Sources: 23 in the 6-station cell
(windows with 0 / 1 / 2 sources: 13 / 11 / 6) and 24 in the one-sided cell (14 / 8 / 8). End of list:
recall and realized false events per window when every emitted event is accepted. MT-FW is
`fw_record` exactly as in Task 1y. Script `qm_wrapper.py test`; output `qm_report.txt`, section 7.

| Cell | method | recall | end of list | events emitted |
|:--|:--|:--|:--|:--|
| 6 stations | QM | 0.13 / 0.13 / 0.22 | 0.22 at 0.40 | 17 |
| | QM-3C | 0.39 / 0.57 / 0.57 | 0.57 at 0.10 | 16 |
| | MT-FW | 0.70 / 0.78 / 0.83 | 0.83 at 4.30 | 150 |
| one-sided, 8 st. | QM | 0.21 / 0.29 / 0.33 | 0.33 at 0.20 | 14 |
| | QM-3C | 0.50 / 0.50 / 0.50 | 0.50 at 0.00 | 12 |
| | MT-FW | 0.75 / 0.75 / 0.79 | 0.88 at 4.17 | 150 |
| pooled (47 sources) | QM | 0.15 / 0.15 / 0.28 | | |
| | QM-3C | 0.53 / 0.53 / 0.53 | 0.53 at 0.05 | |
| | MT-FW | 0.66 / 0.72 / 0.81 | | |

**Table Y.23. Pooled differences in event recall, with 95% paired bootstrap intervals over windows
(1000 resamples).** Same script and output.

| Difference | at $\le0.1$ | at $\le0.25$ |
|:--|:--|:--|
| MT-FW - QM-3C | +0.19 [-0.00, +0.37] | +0.28 [+0.10, +0.43] |
| MT-FW - QM | +0.57 [+0.35, +0.75] | +0.53 [+0.36, +0.71] |
| QM-3C - QM | +0.38 [+0.18, +0.52] | +0.26 [+0.11, +0.44] |

**Where the gap is.** Over the full event lists, QM found 6 of 19 sources in single-source windows
and 7 of 28 in two-source windows, QM-3C 15 of 19 and 10 of 28, and MT-FW 18 of 19 and 22 of 28 (17 of
19 and 21 of 28 at its operating point for 0.25 false events per window). On single sources QM-3C is
close to MT-FW; the gap concentrates in two-source windows. QuakeMigrate keeps one grid maximum per
time sample and merges triggers within its minimum event interval, and both sources of a pair have
origin times between 5.5 and 7.5 s, so the second is usually lost. This is a property of the published
method, which has no explaining-away step, not of its tuning. QM-3C hits have a median epicentral error
of 1.2 km. Run times are comparable: 0.82 and 0.92 s per window for QM-3C end to end, against 0.88
and 1.06 s for MT-FW, on one thread.

**Fairness.** Four caveats bear on the comparison. (1) QuakeMigrate triggers few events (28 in the 60 test windows for QM-3C) and its list ends
at 0.05 false events per window, so its recall at $\le0.1$ and $\le0.25$ is its full-list recall; a
lower trigger floor does not help, because segments above threshold merge (floors at the 0.5, 0.75,
0.9 and 0.95 quantiles gave 0.36, 0.38, 0.21 and 0.17 at $0.25$ in calibration). It cannot trade
false alarms for recall as MT-FW can; this is a property of the method rather than of the protocol.
(2) Its hypothesis space is 40,656 nodes at every 0.02 s,
against 675 nodes and 31 origin times for the grid methods, which raises its noise maximum; a
coarser QuakeMigrate grid was not tested. (3) The onset lag on test hits is larger and more variable
than in calibration (QM-3C median 0.47 s, 90th percentile 0.87 s), so up to about 4 of 47 sources with
origin times near 7.5 s may have their peak pushed past the trigger window. (4) Only the package's
STA/LTA onset was used, not the learned onsets it also accepts (the route of MALMI [110]); the
matched filter of MT-FW is close to the true wavelet. Caveats 2 to 4 can only work against QuakeMigrate. Sampling does not
favour either: the signal band
is well below the 8 Hz Nyquist frequency of the 16 Hz data. With 47 sources one source is 0.021 in
recall, and the selected QM configuration is not clearly better than neighbouring onsets.

## Y.16 Station corrections from catalogued events

Code: `station_corrections.py` (`generate_c`, `peak_delays`, `estimate_station_corrections`,
`froc_recall`, `boot_diff`); output `station_corrections_report.txt`. Y.11 found MT-FW fragile to
unmodelled station delays beyond the $\pm0.25$ s migration window: with delays of standard deviation
0.4 s it fell to $0.37/0.42$. Y.13 made station time corrections calibrated on catalogued events part
of the method. This component study measures whether such corrections, estimated by the method
itself from a small catalogue, remove the fragility.

**Design.** Four networks per cell (random 4 and 6 stations, one-sided 8 stations), each with fixed
statics drawn once from $\mathcal N(0,0.4\ \mathrm s)$ per station and phase and applied to every
window, as real station terms are; between 12 and 62 per cent of the station-phases of a network
exceed the $0.25$ s window. Direction deviations, site factors and weak clutter remain per window. A
catalogue event is drawn uniformly in $[-15,15]^2\times[5,13]$ km with origin time 5.5 to 7.5 s and a
random double couple, from the same generator, and is catalogued if at least three station-phase
arrivals exceed the energy gate (25, five noise standard deviations); 26 to 121 events were drawn per
20 catalogued. Locations and origin times are known exactly, with a sensitivity run adding errors of
1 km per coordinate and 0.1 s. For each catalogued event, station and phase, the delay is the peak of
the three-component matched-filter energy within $\pm1$ s of the predicted time, with parabolic
sub-sample refinement and the same gate, and the correction is the median over events (zero if no
event passes the gate). MT-FW then runs with the corrections added to its predicted P and S times, on
40 fresh windows per network with zero to two sources (160 windows and 134 to 138 events per cell),
through the certified oracle of Appendix Z, which was checked against a dense loop with the same corrections.

**The joint P-before-S search, and why.** The P and S peaks of an event at a station are chosen
jointly, maximizing $E(k_P)+E(k_S)$ subject to P at least 0.375 s before S. In development, separate
searches let the stronger S arrival fall inside the P search window at a station close to the events,
where S-P is below 2 s: 7 of 13 P delays at one station were S, and the median was wrong by 0.9 s.
With the final catalogues the separate search happened to give the same medians (the $K=20$ row of Table Y.24 is
identical with it); the joint search removes the failure mode instead of relying on the median to
outvote it.

**Table Y.24. Residual static after correction (estimated minus true), over all station-phases of the
12 networks.** Uncorrected, the rms is 0.40 s. Per-network rms: range and median over networks.
Beyond window: fraction of station-phases still beyond 0.25 s. No estimate: fraction without any gated
delay (correction left at zero). Script `station_corrections.py`; output
`station_corrections_report.txt` (Results 1 and Table C1).

| Catalogue | pooled rms (s) | per-network rms, range (median) | beyond window | no estimate |
|:--|:--|:--|:--|:--|
| $K=5$ events | 0.171 | 0.008-0.324 (0.123) | 9.1% | 22% |
| $K=10$ | 0.101 | 0.005-0.257 (0.059) | 3.9% | 9.5% |
| $K=20$ | 0.088 | 0.002-0.242 (0.007) | 2.8% | 2.8% |
| $K=20$, catalogue errors 1 km and 0.1 s | 0.173 | 0.067-0.276 | 11.5% | |

Where a station-phase has gated delays, its estimate is within about 0.01 s of the truth, since the
wavelet is exact and the medium homogeneous. Essentially all the residual comes from station-phases
that no catalogued event illuminated above the gate (distant stations, weak P), whose correction stays
at zero; they are also the arrivals that contribute least to the statistic. Catalogue errors
propagate directly: with 1 km and 0.1 s errors, twenty events do no better than five exact ones.

**Table Y.25. MT-FW event recall at realized false events per window $\le0.1$ / $\le0.25$, pooled
over the four networks of a cell (160 windows).** Oracle: the true statics. Script
`station_corrections.py`; output `station_corrections_report.txt` (Table C2).

| Cell | events | none | $K=5$ | $K=10$ | $K=20$ | oracle |
|:--|:--|:--|:--|:--|:--|:--|
| 4 stations | 136 | 0.33/0.40 | 0.49/0.64 | 0.49/0.62 | 0.51/0.64 | 0.51/0.67 |
| 6 stations | 138 | 0.33/0.48 | 0.29/0.52 | 0.45/0.59 | 0.41/0.54 | 0.40/0.56 |
| one-sided, 8 st. | 134 | 0.31/0.41 | 0.66/0.74 | 0.63/0.70 | 0.61/0.71 | 0.66/0.75 |

**Table Y.26. Paired differences in pooled recall at $\le0.1$ ; $\le0.25$, with 95% paired bootstrap
intervals over windows (400 resamples).** Same script and output (Table C3).

| Cell | $K=5$ - none | $K=20$ - none | oracle - none | oracle - $K=20$ |
|:--|:--|:--|:--|:--|
| 4 stations | +0.15 [+0.03, +0.35] ; +0.24 [+0.08, +0.34] | +0.18 [+0.05, +0.35] ; +0.24 [+0.10, +0.33] | +0.18 [+0.05, +0.39] ; +0.27 [+0.14, +0.36] | +0.01 [-0.01, +0.10] ; +0.03 [-0.01, +0.09] |
| 6 stations | -0.04 [-0.19, +0.08] ; +0.04 [-0.13, +0.19] | +0.08 [-0.03, +0.22] ; +0.07 [-0.02, +0.21] | +0.07 [-0.07, +0.19] ; +0.08 [-0.03, +0.19] | -0.01 [-0.11, +0.03] ; +0.01 [-0.07, +0.08] |
| one-sided, 8 st. | +0.35 [+0.22, +0.49] ; +0.33 [+0.18, +0.44] | +0.30 [+0.20, +0.42] ; +0.30 [+0.16, +0.39] | +0.35 [+0.25, +0.49] ; +0.34 [+0.21, +0.45] | +0.05 [-0.01, +0.13] ; +0.04 [-0.02, +0.12] |

**Reading.** In the 4-station and one-sided cells, statics of standard deviation 0.4 s cost MT-FW
0.18 to 0.35 in recall (uncorrected $0.31$-$0.33$ / $0.40$-$0.41$ against oracle $0.51$-$0.66$ /
$0.67$-$0.75$), consistent with the $0.37/0.42$ of Y.11. Corrections from twenty catalogued events
recover $+0.18$ to $+0.30$, with intervals excluding zero at both caps, and five events already
capture most of that gain. The distance to the oracle is $+0.01$ to $+0.05$ and never significant.
Under these conditions the correction step removes the fragility named in Y.11.

**The six-station exception.** In the 6-station cell neither the estimated nor the oracle correction
gives a significant gain: oracle minus none is $+0.07$ / $+0.08$ and $K=20$ minus none $+0.08$ /
$+0.07$, all intervals including zero ($K=10$ minus none reaches $+0.12$ [$+0.01$, $+0.23$] at 0.25,
the only interval there that excludes zero). The uncorrected MT-FW was less hurt in these four
networks (0.33/0.48), and one network (seed 77310) did worse with the oracle correction (0.17/0.24)
than without it (0.20/0.37). With 40 windows per network, per-network recall is noisy: a change of
0.01 s in a predicted time can move it across a sample boundary, change which peak enters the
$\pm4$-sample window and reorder a Frank-Wolfe record. Only the pooled paired comparisons are read,
and in this cell they are inconclusive, not negative.

**Declared limitations.** Exact wavelet and homogeneous medium: real relocations and 3D structure
will leave larger residuals, of the order of the 0.17 s of the catalogue-error run. Statics constant
in time; catalogue events drawn from the same region and generator as the evaluation events; four
networks per cell; exploratory, with no registered rule.

## Y.17 Calibrating the false-alarm threshold

Code: `rft_threshold.py` (`fast_field`, `lkc_box`, `ec_density_chi2`, `eec_threshold`,
`eec_time_threshold`, `bonferroni_threshold`, `gev_threshold`, `conformal_threshold`, `validate`,
`selftest_eec`); output `rft_report.txt`. Section 7.5 sets the detection threshold $\lambda$ at a
quantile of the maximum of the certificate over the hypothesis space, predicted from the expected
Euler characteristic (EEC) of a $\chi^2_6$ random field [119], [120], with an empirical correction on
noise-only windows [125] for non-Gaussian noise. This component study tests both. The statistic is
$T(\theta)=\|Q(\theta)^{\mathsf T}\mathbf v(\theta)\|^2/2$, exactly the WAVE-MT field and the first
Frank-Wolfe score (a vectorized implementation agrees with `fields_x` to $5\times10^{-14}$), and
$M=\max_\theta T$ over 675 nodes and 31 origin times (20,925 hypotheses). In an empty window the
MT-FW catalogue at threshold $u$ is non-empty exactly when $M>u$, so $\alpha=P(M>u)$ is the
per-window false-alarm probability. Three cells (random 4 and 6 stations, one-sided 8 stations), each
with one fixed network, as a deployment sees, and with a new network per window, as in Tasks 1x and
1y; 300 calibration and 300 held-out windows of pure Gaussian noise per setting, and for the fixed
networks 300 plus 300 windows with 35 clutter pulses of amplitude 0.35 or 1.75; 7,200 noise windows in
all. The yardstick is the realized rate on the held-out windows, judged against the 95% predictive
band of an exactly calibrated split-conformal threshold with 300 calibration and 300 test windows,
which is beta-binomial: $[0.00,0.03]$, $[0.02,0.09]$, $[0.057,0.15]$ and $[0.18,0.32]$ at
$\alpha=0.01$, $0.05$, $0.10$ and $0.25$. The naive Clopper-Pearson interval, which ignores the
randomness of the calibration set, is also in the output.

**$2T$ is not $\chi^2_6$.** It would be if $\mathbf v$ were white Gaussian without selection, but each
three-component block is the maximum-energy sample of a 9-sample window of narrow-band filtered noise,
which inflates it. At single hypotheses and pooled, the mean of $2T$ is 12.6 to 13.5 (6 for
$\chi^2_6$), the variance 35 to 50 (12), the 0.99 quantile 29.9 to 34.5 (16.8) and the 0.999 quantile
33 to 49 (22.5); the Kolmogorov-Smirnov distance to $\chi^2_6$ is 0.52 to 0.58. A scaled chi-square
fitted by moments, $2T\sim c\,\chi^2_k$ with $c=1.3$-$1.9$ and $k=7.1$-$9.6$, fits the bulk (distance
0.010 to 0.045) but overstates the upper tail, which is lighter, $T\sim(0.54$-$0.68)\,\chi^2_{10.6-14.0}$
above the 0.999 quantile, a Gumbel-like tail produced by the selection. The marginal hardly depends on
the network (pooled mean of $2T$ 13.0 to 13.1 in every cell). The random-field predictors below
therefore use the pooled empirical distribution with a fitted chi-square tail, mapped to its
$\chi^2_6$ equivalent.

**The EEC is correct, and overcounts the tests.** The EC densities of a $\chi^2_6$ field and the
Lipschitz-Killing curvatures (LKC) of a box were verified by simulating a smooth $\chi^2_6$ field on a
$10^4$ lattice: EEC 7.007, 4.664 and 1.874 against a mean lattice EC of 6.98, 4.62 and 1.70
($\pm0.26$, $0.18$, $0.12$) at $u=12$, 16 and 20, and EEC 0.158 against $P(\max>28)=0.12$. The
certificate field, however, is rough on the grid: the lag-1 correlation of its $\chi^2_6$ equivalent is
0.13 to 0.20 along $x$, 0.09 to 0.18 along $y$ and 0.09 to 0.21 along $z$ (steps of 3, 3 and 4 km),
and 0.51 to 0.63 along the origin time (0.1 s step), because a 3 km move shifts the S times by up to
1 s against a $\pm0.25$ s window. The roughness-based top LKC is $1.6$ to $3.1\times10^4$, more than
the 20,925 lattice points. Measured on the data, the effective number of independent hypotheses
at $\alpha=0.05$ is 8,700 to 19,000, while the continuum EEC implies 330,000 to 640,000. The continuum
EEC with roughness-based LKC therefore overestimates $P(M>u)$ by 20 to 60 times, and its thresholds are
5 to 9 units too high. Two lattice-aware variants, a union bound over nodes of the one-dimensional EEC
in origin time (EEC-time) and LKC fitted to the empirical lattice EC curve (EEC-ECfit), remain
conservative by about four times, and Bonferroni with the empirical marginal is the least conservative
analytic bound.

**Table Y.27. Realized false-alarm rate $P(M>u)$ on 300 held-out windows of Gaussian noise, mean over
the six settings (three cells, fixed and per-window networks), at nominal
$\alpha=0.01$ / $0.05$ / $0.10$ / $0.25$, and number of the 24 (setting, $\alpha$) combinations inside
the predictive band.** EEC-rough: continuum EEC with roughness-based LKC. GEV and Gumbel: extreme-value
fits to the calibration maxima. Conformal: split-conformal quantile of the calibration maxima. Script
`rft_threshold.py`; output `rft_report.txt`, section 3.

| Predictor | realized rate | inside band |
|:--|:--|:--|
| EEC-rough | 0.000 / 0.001 / 0.001 / 0.002 | 6/24 |
| EEC-time | 0.001 / 0.012 / 0.029 / 0.074 | 8/24 |
| EEC-ECfit | 0.001 / 0.012 / 0.028 / 0.075 | 8/24 |
| Bonferroni | 0.003 / 0.021 / 0.047 / 0.108 | 10/24 |
| GEV | 0.006 / 0.042 / 0.094 / 0.237 | 24/24 |
| Gumbel | 0.003 / 0.032 / 0.087 / 0.238 | 24/24 |
| conformal | 0.007 / 0.037 / 0.095 / 0.239 | 24/24 |

The six in-band combinations of EEC-rough are all at $\alpha=0.01$, where the band includes zero. The
held-out maxima have mean 22.5 (4 stations), 25.0 (6 stations) and 26.4 (8 stations), with standard
deviation 2.2 to 2.6. Calibrating on windows with a new random network each time works as well as on
the fixed network: the maximum depends on the number of stations more than on their positions.

**Table Y.28. Mean ratio of realized to nominal false-alarm rate over cells and the four values of
$\alpha$ (1: calibrated; below 1: conservative; above 1: anti-conservative), by noise regime and
calibration set.** Clutter settings use the fixed networks (12 combinations each). Same script and
output (Summary).

| Predictor | Gaussian | weak clutter, Gaussian-calibrated | weak clutter, recalibrated | strong clutter, Gaussian-calibrated | strong clutter, recalibrated |
|:--|:--|:--|:--|:--|:--|
| EEC-rough | 0.01 | 0.01 | 0.00 | 6.22 | 0.01 |
| EEC-time | 0.23 | 0.37 | 0.32 | 17.22 | 0.31 |
| EEC-ECfit | 0.23 | 0.37 | 0.25 | 16.82 | 0.34 |
| Bonferroni | 0.42 | 0.66 | 0.53 | 19.40 | 0.41 |
| GEV | 0.82 | 1.19 | 0.89 | 21.86 | 0.96 |
| Gumbel | 0.70 | 0.81 | 0.80 | 18.80 | 0.77 |
| conformal | 0.83 | 1.11 | 0.84 | 22.18 | 0.96 |

**Clutter.** Weak clutter (amplitude 0.35) barely moves the maximum (mean $+0.0$ to $+0.3$), and the
GEV and conformal thresholds calibrated on Gaussian noise stay inside the band in 12 of 12
combinations each, with a slight excess at small $\alpha$: pooled realized rates 0.016 / 0.050 / 0.115
/ 0.261 (GEV) and 0.015 / 0.034 / 0.127 / 0.264 (conformal). Strong clutter (amplitude 1.75, visible
transients that pickers do pick) shifts the maximum by 8 to 9 units (mean 31.6, 33.8 and 34.3 against
22.5, 25.0 and 26.4) and doubles its spread (standard deviation 4.2 to 4.9). Every threshold calibrated
on Gaussian noise then fails, including the random-field ones: the realized rate is 0.60 / 0.78 / 0.90
/ 0.96 for conformal and 0.14 / 0.27 / 0.35 / 0.47 even for the very conservative EEC-rough, and on
average 19 to 22 times nominal for the empirical predictors (Bonferroni 19.40, GEV 21.86, Gumbel 18.80,
conformal 22.18; EEC-time, at 17.22, is not empirical). Recalibrated on separate windows with the same clutter
law, GEV and conformal are calibrated again: GEV 0.005 / 0.056 / 0.111 / 0.263, conformal 0.004 /
0.062 / 0.102 / 0.281 and Gumbel 0.002 / 0.043 / 0.092 / 0.263, each inside the band in 12 of 12
combinations. The GEV and conformal thresholds at $\alpha=0.05$ rise from 27.3 to 31.2 on Gaussian
noise to 38.8 to 42.7 with strong clutter. Re-estimating the marginal and the roughness on clutter
windows does not rescue the EEC, which stays conservative by 3 to 5 times (EEC-time, EEC-ECfit) to more
than 50 times (EEC-rough).

**Conclusion.** The claim that the EEC of a $\chi^2_6$ field fixes the false-alarm rate per window
before any data are seen is not supported, for three separate reasons: the marginal is not
$\chi^2_6$, the field decorrelates within one grid step so that the continuum counts far more tests
than the lattice has, and non-Gaussian transients move the maximum far beyond any Gaussian
prediction. The primary threshold is split-conformal, which guarantees $P(M>u)\le\alpha$ under
exchangeability in finite samples, or a GEV fit to window maxima, which is smooth and extrapolates
below $1/n_{\mathrm{cal}}$, both from about 300 noise-only windows of the same network and the same
noise and clutter regime, and recalibrated whenever that regime changes. The EEC or Bonferroni bound
is kept only as a conservative prior threshold before such windows exist. Caveats: 300 calibration and
300 test windows per setting resolve $\alpha=0.01$ poorly (three expected exceedances); one fixed
network per cell; synthetic white Gaussian noise; the clutter law is identical in calibration and
test, which is the favourable case for any empirical calibration.

### Y.17a What Y.14 to Y.17 change

The fragility recorded in Y.11 and Y.13 is reduced to a calibration step. Station delays of standard
deviation 0.4 s, largely beyond the migration window, cost MT-FW 0.18 to 0.35 in recall, and
corrections that the method
estimates from five to twenty catalogued events with a joint P-before-S search restore the recall of
the true statics in the 4-station and one-sided cells, with the distance to the oracle never
significant (Y.16). The step is conditional on catalogue quality, since location errors of 1 km and
0.1 s double the residual static, and in the six-station cell the gain was not significant even with
the true statics, so the claim is that the fragility is removable, not that every network benefits
equally.

The advantage of MT-FW over a published migration locator is supported in exploratory form. Against
QuakeMigrate 1.2.2 in its fair three-component configuration MT-FW is ahead by $+0.28$ [$+0.10$,
$+0.43$] at 0.25 false events per window on 47 sources (Y.15), and the gap sits in two-source windows,
where the explaining-away step of the greedy deflation has no counterpart in the published
method. On isolated sources the two are close. This agrees with the decomposition of Y.13, in which
explaining-away separates overlapping events, but one seed and two cells do not make it confirmatory;
Task 2A (Y.18) is the confirmatory test, and its primary rule holds in both clutter regimes.

The pick-based comparisons are qualified. The picks that fed GaMMA from Y.5 to Y.12 came from a model
that detects three times as many arrivals as a calibrated STA/LTA, ten times as many as an
out-of-domain PhaseNet, with timing errors about 2.4 times smaller (Y.14). With weak clutter it was
optimistic for GaMMA on every count, which makes the weak-clutter margins of MT-FW over GaMMA in Y.12
conservative with respect to the picker, though not with respect to the medium. With strong clutter it also delivered about 28 clutter picks
per window where the STA/LTA gives 7 to 8, so the strong-clutter margins of Y.12 rest in part on a
clutter pick rate that a real picker would not produce, and they remained qualified until GaMMA ran on
real picks in Task 2A (Y.18). There the STA/LTA picks starve GaMMA, and with the parametric picks
MT-FW with corrections stays ahead by $0.15$ to $0.31$ at $0.25$ under strong clutter, so the
strong-clutter margin over GaMMA holds with both pick sources on this generator. The generator of Y.14 also breaks the identity between the physics of the
synthesis and that of the detectors (S timing errors from $-0.38$ to $+0.24$ s between the 5th and
95th percentiles, arrival directions 12.4 degrees off at the median); every Y result before it was
obtained without that mismatch.

The threshold of Section 7.5 is weakened as a theoretical claim and strengthened as a procedure. The
random-field prediction is correct as mathematics and wrong as a model of this statistic. On Gaussian
noise it is conservative: the continuum form with roughness-based LKC overestimates $P(M>u)$ at a given
threshold by 20 to 60 times, and at its own thresholds it realizes about 0.01 of the nominal rate, the
lattice-aware forms about 0.23 (Table Y.28). Under strong clutter the EEC and Bonferroni thresholds
calibrated on Gaussian noise realize 6 to 19 times the nominal rate (Y.17). The empirical calibration on noise-only windows is therefore the primary
method (Section 7.5, point 4), and Y.17 shows that it holds its nominal rate under both noise regimes
provided it is calibrated in the regime it is applied to. Every earlier Y result that read recall at realized false
events per window on the evaluation windows is unaffected, since it never used a predicted threshold.

Two items registered in Y.13 (item 7) are done in part and remain open in their real-data form: a real picker
on the synthetic waveforms, but in white noise rather than real noise, and the threshold, checked on
synthetic rather than real noise. The certified oracle of Appendix Z is exact but, for the moment-tensor
certificate, not a pruning method (Appendix Z).

## Y.18 The confirmatory test (Task 2A)

Code: `task2a_confirm.py` (`selfcheck`, `qmcal`, `smoke`, `run`), `analyze_2a.py`; output
`analyze_2a_out.txt`; design and rules in `PREREGISTRO_1j_1k.txt` (section "Task 2A"); independent
verification in `task2a_report.txt`, with its scripts and outputs in `verif_2a/`.

**The question, and why the test is confirmatory.** Y.12 found MT-FW ahead of tuned GaMMA, but with a
synthesis whose physics is that of the detectors, with the parametric pick model `degrade_y`, which Y.14
showed to be optimistic for GaMMA in recall and timing, and without a published migration locator. Y.14
to Y.16 removed each of these shortcuts separately and in exploratory form. Task 2A removes them
together and asks one question: with a generator the detectors do not know, picks from an operational
picker and published comparators, does MT-FW with station corrections estimated from a catalogue
exceed the best comparator? The design, the comparators, the endpoint, five decision rules and their
interpretation were written on 25 September 2026, before any evaluation run of this task. The time in
the entry, 20:15 UTC, was typed rather than read from a clock, like every registry time before Task 2D
(E34 of the error register), and the registry file itself was created at 20:14:09 UTC; no hash fixes
the text of the entry. The file-system record is consistent with the registration: the first file of
the task, `task2a_confirm.py`, was created at 23:33:22 UTC, and the operating points of rule 2A-e were
frozen at 01:12:25 UTC on 26 September, before the first evaluation output at 01:12:53 UTC. The
registration is therefore held in an append-only local file, and its timing is not externally
verifiable.
The evaluation seeds 2024 and 7331 had never been used with the generator of Y.14, the picker or the
QuakeMigrate wrapper; every calibration (the STA/LTA picker of Y.14, the QuakeMigrate configurations and
the frozen operating points of rule 2A-e) used seed 999 only, on networks distinct from the evaluation
networks. The registration states that nothing registered is modified after results of seeds 2024 or
7331 are seen and that any analysis added afterwards is exploratory; the exploratory analyses below are
marked as such. The claim it confirms is a claim about synthetic data.

**Design.** The generator is that of Y.14: a linear gradient medium with $V_{S0}=2.8$ km/s and
$g=0.04$ s$^{-1}$ and curved rays, a Ricker wavelet per event with centre frequency uniform in
$[1.4,2.8]$ Hz and a phase rotation of up to $\pm45$ degrees, a scattered coda with 0.3 of the energy of
the direct arrival, direction deviations of about ten degrees, lognormal site factors with $\sigma=0.5$,
and random double-couple sources. One element differs: station delays are drawn once per network,
$\mathcal N(0,0.3\ \mathrm s)$ per station and phase, and kept fixed across its windows, as real station
terms are, instead of being drawn per window. Every detector uses the same homogeneous model, $V_P=5.19$
and $V_S=3.0$ km/s, so the model error is common to all methods. The cells are those of Y.11 and Y.12
(random networks of 4 and 6 stations, and one-sided networks of 8 stations with 40 clutter pulses), in
the two clutter regimes of Y.12, weak (amplitude 0.35, "thin" in the code) and strong (1.75). For each
seed and cell there are five fixed networks, shared by the two regimes, each with 20 evaluation windows
(zero, one or two sources with probabilities 0.4, 0.4 and 0.2, at least 12 km apart, origin times in
5.5 to 7.5 s): 100 windows and 68 to 87 events per cell. The registered budget rule would have reduced
the windows per network to 12 had the seed-999 estimate exceeded 10 hours on two cores; it was 1.47
hours, and the number stayed at 20.

For each network and regime a catalogue of 20 events is drawn from the same generator on that network;
an event enters it when at least three station-phase arrivals exceed the matched-filter energy gate of
Y.16 (25) within $\pm1$ s of the homogeneous prediction, and between 23 and 183 events were drawn per 20
catalogued. Catalogue locations and origin times carry errors of 1 km per coordinate and 0.1 s. The
station corrections (SC) are the estimator of Y.16, the median over the catalogue of the jointly
searched P and S peak delays, evaluated at the catalogued, perturbed locations and origin times. The
same SC are given to every method that admits them: added to the predicted P and S times of MT-FW,
subtracted from the pick times of GaMMA and PyOcto per station and phase, and added to the travel-time
lookup table of QuakeMigrate, where at most $3.8\times10^{-5}$ of the entries became negative and were
clipped to zero by the library, none at source depths of 5 to 13 km. Picks come from the STA/LTA picker
of Y.14, with its configuration frozen on seed 999. As an optimistic bound for the associators,
`degrade_y` is also applied to the true arrival times and amplitudes of the generator.

**Methods.** MT-FW and MT-FW+SC are the greedy group matching pursuit of Y.11 through the certified
oracle of Appendix Z (at most five atoms, a 5 km exclusion and the per-station cap of Y.11), without
and with SC; their records were asserted identical to those of the dense loop on two windows of every
network. MT-FW-NODE+SC is the same iteration with the subspace of a node at least 9 km away, the control
of the complete method. WAVE-PT@20, the mechanism-independent migration with a 20 km suppression radius,
is a reference without SC. GaMMA-STALTA and GaMMA-STALTA+SC run on the STA/LTA picks over the ten
seeded configurations of the grid of Y.12, and PyOcto-STALTA and PyOcto-STALTA+SC over fifteen; both
choose their configuration by an oracle on the evaluation windows for each cell and cap. QM, QM-3C and
QM-3C+SC are QuakeMigrate 1.2.2 with the default and the three-component channel mappings of Y.15,
recalibrated by the procedure of Y.15 on seed-999 windows of this generator and frozen before
evaluation (QM-3C: 1-6 Hz band, STA/LTA 0.8/6.0 s); QM-3C+SC is the QM-3C configuration with SC.
GaMMA-DY and GaMMA-DY+SC run on the optimistic picks, also with an oracle configuration.

**Endpoint and rules.** The endpoint is event recall at at most 0.1 and 0.25 realized false events per
window, read on the exact FROC curve, with the operating point of every method chosen on the evaluation
windows in the same way; the associators also choose their configuration by oracle, which favours them.
A hit is the nearest true source within 5 km epicentral distance, scored by the same rule for every
method. Intervals come from a bootstrap over windows, paired across methods, with 2000 replicates and
the operating point, the configuration and the best comparator re-chosen in each. Each rule is decided
per clutter regime on point estimates, as in Y.12, and requires a difference of at least $0.05$ in at
least two of three cells in both seeds. Rule 2A-a, the primary rule, compares MT-FW+SC with the maximum
over GaMMA-STALTA, GaMMA-STALTA+SC, PyOcto-STALTA, PyOcto-STALTA+SC, QM, QM-3C and QM-3C+SC, taken
per cell and cap, at both caps. Rule 2A-b compares MT-FW+SC with MT-FW-NODE+SC at both caps, and if it
fails the gain is not attributed to the radiation pattern. Rule 2A-c, secondary, compares MT-FW+SC
with the better of GaMMA-DY and GaMMA-DY+SC at 0.25. Rule 2A-d, secondary, compares MT-FW+SC with
MT-FW at 0.25. Rule 2A-e, secondary, tests a deployable operating point: the score thresholds of
MT-FW+SC and QM-3C+SC and the configuration of GaMMA-STALTA+SC are frozen on seed 999 (five networks
of 20 windows per cell and regime) for budgets of 0.1 and 0.25 false events per window, and the rule
holds if the realized false rate of MT-FW+SC does not exceed twice the budget and its recall exceeds
that of the best frozen comparator by at least $0.05$, in at least two of three cells at both budgets,
in both seeds. The frozen values were written to file before the first evaluation output existed. The
registered interpretation: if 2A-a holds in both regimes, the advantage of Y.12 survives the
mismatched generator, the real picker and QuakeMigrate; if in one, the claim is restricted to that
regime; if in neither, the claim of superiority over the state of the art is withdrawn and reported as
a negative result.

**Table Y.29. Weak clutter. Event recall at realized false events per window $\le0.1$ / $\le0.25$.**
Best comparator: the best of the seven comparators of rule 2A-a at each cap, with its name. Script
`analyze_2a.py`; output `analyze_2a_out.txt`, Tables 1 and 3 (thin).

| Cell | events | MT-FW+SC | MT-FW | MT-FW-NODE+SC | best comparator | GaMMA-STALTA+SC | QM-3C+SC | GaMMA-DY+SC |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| one-sided, 8 st. (2024) | 71 | 0.55/0.65 | 0.34/0.46 | 0.24/0.41 | 0.35/0.38 (QM-3C+SC) | 0.06/0.06 | 0.35/0.38 | 0.63/0.63 |
| 4 stations (2024) | 87 | 0.32/0.52 | 0.22/0.33 | 0.07/0.15 | 0.18/0.25 (QM-3C+SC) | 0.02/0.02 | 0.18/0.25 | 0.34/0.34 |
| 6 stations (2024) | 72 | 0.46/0.56 | 0.29/0.50 | 0.22/0.26 | 0.28/0.31 (QM-3C+SC) | 0.08/0.08 | 0.28/0.31 | 0.47/0.47 |
| one-sided, 8 st. (7331) | 81 | 0.64/0.77 | 0.47/0.59 | 0.43/0.48 | 0.47/0.48 (QM-3C+SC) | 0.09/0.11 | 0.47/0.48 | 0.75/0.75 |
| 4 stations (7331) | 78 | 0.32/0.51 | 0.46/0.54 | 0.15/0.33 | 0.15/0.17 (QM-3C+SC) | 0.01/0.01 | 0.15/0.17 | 0.28/0.28 |
| 6 stations (7331) | 68 | 0.72/0.79 | 0.38/0.56 | 0.24/0.44 | 0.24/0.34 (QM-3C+SC) | 0.12/0.12 | 0.24/0.34 | 0.68/0.68 |

**Table Y.30. Strong clutter. The same.** Where the best comparator differs between caps, both names
are given (at 0.1 / at 0.25). Same script and output (strong).

| Cell | events | MT-FW+SC | MT-FW | MT-FW-NODE+SC | best comparator | GaMMA-STALTA+SC | QM-3C+SC | GaMMA-DY+SC |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| one-sided, 8 st. (2024) | 83 | 0.31/0.42 | 0.22/0.31 | 0.19/0.28 | 0.16/0.19 (QM-3C / QM-3C+SC) | 0.05/0.05 | 0.16/0.19 | 0.20/0.23 |
| 4 stations (2024) | 77 | 0.31/0.39 | 0.21/0.29 | 0.05/0.08 | 0.12/0.14 (QM-3C) | 0.01/0.01 | 0.12/0.13 | 0.01/0.05 |
| 6 stations (2024) | 74 | 0.38/0.49 | 0.20/0.42 | 0.07/0.16 | 0.12/0.16 (QM-3C / QM-3C+SC) | 0.00/0.01 | 0.12/0.16 | 0.11/0.14 |
| one-sided, 8 st. (7331) | 75 | 0.63/0.71 | 0.55/0.61 | 0.23/0.39 | 0.25/0.29 (QM-3C+SC) | 0.04/0.04 | 0.25/0.29 | 0.00/0.36 |
| 4 stations (7331) | 78 | 0.13/0.22 | 0.14/0.26 | 0.09/0.17 | 0.04/0.05 (QM-3C+SC) | 0.01/0.01 | 0.04/0.05 | 0.04/0.06 |
| 6 stations (7331) | 84 | 0.36/0.42 | 0.23/0.38 | 0.12/0.21 | 0.18/0.25 (QM-3C+SC) | 0.02/0.04 | 0.18/0.25 | 0.18/0.20 |

**Table Y.31. MT-FW+SC minus the best comparator, the quantity of rule 2A-a, with 95% paired bootstrap
intervals over windows (maximum re-chosen in each replicate).** Differences are computed before
rounding, so they can differ by 0.01 from the rounded entries of Tables Y.29 and Y.30. Same script and
output, Table 3.

| Cell | weak, $\le0.1$ | weak, $\le0.25$ | strong, $\le0.1$ | strong, $\le0.25$ |
|:--|:--|:--|:--|:--|
| one-sided, 8 st. (2024) | +0.20 [+0.03, +0.39] | +0.27 [+0.12, +0.40] | +0.16 [+0.00, +0.27] | +0.23 [+0.10, +0.36] |
| 4 stations (2024) | +0.14 [-0.01, +0.34] | +0.26 [+0.10, +0.42] | +0.19 [+0.01, +0.27] | +0.25 [+0.11, +0.35] |
| 6 stations (2024) | +0.18 [+0.00, +0.35] | +0.25 [+0.09, +0.40] | +0.26 [+0.08, +0.40] | +0.32 [+0.19, +0.45] |
| one-sided, 8 st. (7331) | +0.17 [+0.00, +0.35] | +0.28 [+0.12, +0.42] | +0.37 [+0.20, +0.54] | +0.41 [+0.26, +0.54] |
| 4 stations (7331) | +0.17 [+0.03, +0.35] | +0.35 [+0.16, +0.47] | +0.09 [+0.01, +0.21] | +0.17 [+0.06, +0.26] |
| 6 stations (7331) | +0.49 [+0.25, +0.65] | +0.46 [+0.30, +0.59] | +0.18 [-0.01, +0.30] | +0.17 [+0.04, +0.32] |

**Table Y.32. Rule 2A-e: operating points frozen on seed 999 and applied unchanged. Realized recall /
false events per window of MT-FW+SC and of the best frozen comparator, and the difference in recall.**
The best frozen comparator is QM-3C+SC except in the strong-clutter 4-station cell of seed 7331, where it
is GaMMA-STALTA+SC. "Over 2x" marks a realized rate above twice the budget. Every one of the 24 intervals
of the difference excludes zero (lowest lower bound $+0.01$). Same script and output, Table 4.

| Cell | budget | weak: MT-FW+SC | weak: comparator | weak: difference | strong: MT-FW+SC | strong: comparator | strong: difference |
|:--|:--|:--|:--|:--|:--|:--|:--|
| one-sided, 8 st. (2024) | 0.1 | 0.56 / 0.13 | 0.35 / 0.10 | +0.21 | 0.34 / 0.12 | 0.18 / 0.14 | +0.16 |
| | 0.25 | 0.65 / 0.28 | 0.37 / 0.13 | +0.28 | 0.39 / 0.20 | 0.19 / 0.33 | +0.19 |
| 4 stations (2024) | 0.1 | 0.44 / 0.17 | 0.18 / 0.11 | +0.25 | 0.22 / 0.05 | 0.12 / 0.09 | +0.10 |
| | 0.25 | 0.52 / 0.26 | 0.18 / 0.11 | +0.33 | 0.32 / 0.17 | 0.12 / 0.09 | +0.21 |
| 6 stations (2024) | 0.1 | 0.47 / 0.12 | 0.28 / 0.11 | +0.19 | 0.22 / 0.05 | 0.08 / 0.09 | +0.14 |
| | 0.25 | 0.60 / 0.35 | 0.29 / 0.13 | +0.31 | 0.43 / 0.14 | 0.15 / 0.22 | +0.28 |
| one-sided, 8 st. (7331) | 0.1 | 0.77 / 0.22 (over 2x) | 0.47 / 0.08 | +0.30 | 0.71 / 0.23 (over 2x) | 0.29 / 0.20 | +0.41 |
| | 0.25 | 0.79 / 0.44 | 0.47 / 0.13 | +0.32 | 0.71 / 0.33 | 0.29 / 0.35 | +0.41 |
| 4 stations (7331) | 0.1 | 0.31 / 0.08 | 0.15 / 0.10 | +0.15 | 0.10 / 0.04 | 0.01 / 0.01 | +0.09 |
| | 0.25 | 0.38 / 0.17 | 0.15 / 0.10 | +0.23 | 0.17 / 0.12 | 0.01 / 0.01 | +0.15 |
| 6 stations (7331) | 0.1 | 0.71 / 0.10 | 0.31 / 0.18 | +0.40 | 0.36 / 0.13 | 0.20 / 0.13 | +0.15 |
| | 0.25 | 0.79 / 0.33 | 0.34 / 0.22 | +0.46 | 0.51 / 0.41 | 0.25 / 0.22 | +0.26 |

**Table Y.33. Registered verdicts: cells in which the rule's difference reaches $0.05$ (seed 2024 ; seed
7331), and the verdict per clutter regime.** Rule 2A-e is given under both readings discussed below.
Same script and output (Rules and Verdicts).

| Rule | weak: cells | strong: cells | weak: verdict | strong: verdict |
|:--|:--|:--|:--|:--|
| 2A-a, state of the art (primary) | 3/3 ; 3/3 | 3/3 ; 3/3 | holds | holds |
| 2A-b, geometry-breaking control | 3/3 ; 3/3 | 3/3 ; 2/3 | holds | holds |
| 2A-c, optimistic picks | 2/3 ; 2/3 | 3/3 ; 3/3 | holds | holds |
| 2A-d, station corrections | 3/3 ; 2/3 | 3/3 ; 1/3 | holds | fails |
| 2A-e, per-cell reading | 3/3 ; 2/3 | 3/3 ; 2/3 | holds | holds |
| 2A-e, strict reading | rate exceeded in seed 7331 | rate exceeded in seed 7331 | fails | fails |

**Rules 2A-a and 2A-b hold in both regimes.** MT-FW+SC exceeds the best of the seven comparators in all
three cells of both seeds, by $+0.14$ to $+0.49$ at 0.1 and $+0.25$ to $+0.46$ at 0.25 with weak
clutter, and by $+0.09$ to $+0.37$ and $+0.17$ to $+0.41$ with strong clutter (Table Y.31). By the
registered interpretation, the advantage of Y.12 survives the mismatched generator, the operational
picker and QuakeMigrate. The same iteration with the subspace of another node loses $+0.17$ to $+0.49$
at both caps with weak clutter and $+0.12$ to $+0.40$ with strong clutter, except in the 4-station cell
of seed 7331, where it loses only $+0.04$ at 0.1 and $+0.0513$ at 0.25, the only difference in Task 2A within
0.005 of a threshold (Task 2B has one, $+0.0494$, Appendix Y.19); the rule still holds there, with two
of three cells.

**Rule 2A-c holds, with a tie in weak clutter.** With weak clutter GaMMA-DY+SC, fed the optimistic
picks and the SC, ties MT-FW+SC at 0.25 in the one-sided cell of both seeds ($+0.01$), and the rule
passes on the other two cells ($+0.08$ to $+0.23$). With strong clutter MT-FW+SC is ahead by $+0.15$ to
$+0.31$, because the optimistic model also delivers 27.25 to 32.60 clutter picks per window. The rule
reads the cap of 0.25 only. At 0.1, outside the rule, GaMMA-DY+SC is ahead of MT-FW+SC in four of the
six weak-clutter cells of the full sample, by 0.01 to 0.11 (4 and 6 stations and one-sided in seed
2024, one-sided in seed 7331; output Table 2, thin), and behind it in the other two by 0.04; with
strong clutter MT-FW+SC is ahead at 0.1 in all six cells. With optimistic picks and weak clutter the
advantage at the stricter budget is therefore not established.

**Rule 2A-d holds with weak clutter and fails with strong clutter.** With weak clutter the SC add
$+0.06$ to $+0.24$ at 0.25 in five cells and cost $-0.03$ in the 4-station cell of seed 7331 (and
$-0.14$ at 0.1 there). With strong clutter they add $+0.07$ to $+0.11$ in seed 2024 but reach $0.05$ in
only one cell of seed 7331 ($-0.04$, $+0.04$ and $+0.09$ at 4 stations, 6 stations and one-sided). The
cause is a departure from the registered text: the registration speaks of one catalogue per network,
and the implementation drew one per network and regime, each with the clutter of its regime, so the
strong-clutter SC are estimated from catalogue windows contaminated by visible transients. The
residual rms of SC minus the true statics is 0.12 to 0.21 s with weak clutter and 0.20 to 0.28 s with
strong clutter, against 0.26 to 0.34 s for the statics themselves; with strong clutter the correction
removes little of the static. An exploratory rerun of all 300 strong-clutter windows of seed 7331
confirms the cause: corrections estimated from the weak-clutter catalogues of the same networks leave
residuals of 0.12 to 0.15 s and give MT-FW recall of 0.26/0.31, 0.37/0.48 and 0.65/0.72 at 4 stations,
6 stations and one-sided, against 0.13/0.22, 0.36/0.42 and 0.63/0.71 with the SC used, and the true
statics give 0.23/0.28, 0.45/0.56 and 0.73/0.79. The failure is therefore a property of the estimator
on clutter-contaminated catalogue windows, not a sign or wiring error, and it affects every +SC method
equally, since all receive the same corrections. For the method it makes the quality of the catalogue
windows a condition of the correction step.

**Rule 2A-e holds under one reading and fails under the other; both are reported.** The recall
condition holds in every cell at both budgets in both seeds and both regimes, by $+0.09$ to $+0.46$
(Table Y.32). The rate condition fails in exactly one cell per regime: the one-sided 8-station cell of
seed 7331 at budget 0.1, where the threshold frozen on seed 999 realizes 0.22 (weak) and 0.23 (strong)
false events per window against a limit of 0.20, while the recall advantage there is $+0.30$ and
$+0.41$. The registered sentence does not say whether the rate clause is quantified per cell, as the
recall clause is, or over every cell and budget, and both readings were coded in the analysis script
before any evaluation data existed. Under the per-cell reading, in which a cell counts when both
conditions hold in it at both budgets, the rule holds, with three and two of three cells; under the
strict reading, in which the rate condition must hold in every cell and budget, it fails in both
regimes. It is not claimed as an unqualified pass: an absolute score threshold frozen on other networks
realized up to 2.3 times its budget.

**Verification.** An independent verifier reviewed the code and outputs after the runs, without
modifying them, and wrote every check to `verif_2a/`. No leakage was found: the truth enters the MT-FW,
WAVE-PT and QuakeMigrate records only in the hit-labelling step, after a hypothesis has been selected;
the SC come from the catalogue windows alone, whose random streams are independent of those of the
evaluation windows, and use the perturbed catalogue locations; the frozen values of 2A-e depend only on
seed-999 results and were re-derived independently; and no seed-999 network coincides with an
evaluation network (closest pair of stations 0.34 km apart, no identical static matrices). The hit rule
is identical for every method, and GaMMA and PyOcto, which have no separation rule and no cap on the
number of events, are scored more leniently than MT-FW and QuakeMigrate; widening QuakeMigrate's
trigger window by 0.5 s on each side changes its recall by $-0.03$ to $+0.01$, so the time window does
not handicap it. The sign of the SC was checked for every family with the true statics as oracle
corrections, on 80 regenerated evaluation windows of seed 2024 with weak clutter: recall at 0.1 / 0.25
rises from 0.25/0.36 to 0.51/0.54 for MT-FW, from 0.17/0.24 to 0.36/0.36 for QM-3C and from 0.29/0.36 to
0.41/0.41 for GaMMA-DY, and the estimated SC recover most of that gain (0.44/0.56, 0.32/0.37 and
0.39/0.39). The near-zero recall of GaMMA and PyOcto on STA/LTA picks is starvation by the picker, not a
wiring fault: the frozen picker finds 9.6 to 16.8 per cent of the arrivals within 0.3 s, on those 80
windows 1.53 correctly labelled picks of 11.5 arrivals per event, so only 20 per cent of events receive
the three picks that the smallest GaMMA minimum pick count requires; the matched picks have a median
error of $-0.009$ s and a standard deviation of 0.13 s, so there is no time offset; and with the same
wiring, positive controls on picks at the generator's true times reach 0.88 recall at 0.1 for GaMMA and
0.69 for PyOcto. The tabulated numbers were re-derived by an independent FROC implementation, which
reproduces Table 1 of the output to the printed precision in the four cells checked and the realized
values of 2A-e. Beyond the catalogue per regime, the verifier recorded three minor departures: the
generator was edited after the registration to add the registered fixed-statics argument, with its
random streams verified unchanged; four seed-999 part files came from an earlier version of the task
script, and a rerun with the final code reproduced them except for QuakeMigrate scores in the fifth
decimal and one extra QM-3C+SC false event below the frozen threshold at 0.1, an effect on 2A-e that is
nil or marginal; and QuakeMigrate was recalibrated on weak clutter only, in the 6-station and one-sided
cells, on 60 windows with 44 sources and without SC, so its strong-clutter and 4-station performance may
be below what a regime-specific calibration would give, while the associators choose their
configuration by oracle on the evaluation windows. No implementation error, leakage path or scoring
asymmetry favouring MT-FW was found.

**Qualifications.** Five limit what the result says. First, the registered intervals resample windows
and ignore that the windows of a cell share five networks, over which the per-network recall of
MT-FW+SC at 0.1 ranges from 0.00 to 0.88. An exploratory network-level cluster bootstrap (2000
replicates, maximum re-chosen in each) keeps every cell above zero at 0.25 (lowest lower bound $+0.03$),
but at 0.1 four of the twelve intervals include zero: weak clutter, seed 2024, 6 stations, [$-0.03$,
$+0.34$]; weak, seed 7331, one-sided, [$-0.18$, $+0.41$]; strong, seed 2024, one-sided, [$-0.03$,
$+0.36$]; strong, seed 7331, 6 stations, [$-0.05$, $+0.29$]. With five clusters these intervals are
crude, but they show that the per-window intervals of Table Y.31 understate the uncertainty; the rules
are decided on point estimates, as registered. Task 2C repeated the comparison with QM-3C+SC over
thirty networks per cell with a hierarchical bootstrap (Appendix Y.20). Second, the effective
competitor is QuakeMigrate, the
best comparator in all 24 cases of cell, seed and cap. Because the picker starves the pick-based
associators, the claim against GaMMA and PyOcto holds with an operational STA/LTA picker, and not in
general. Third, about half of all events lie in two-source windows whose origin times are drawn in the
same 2 s interval (in 137 of 215 such windows within 0.75 s of each other), which QuakeMigrate's
trigger, merging peaks closer than its minimum event interval, cannot separate. Restricted to windows
with at most one source (exploratory), MT-FW+SC still exceeds QM-3C, QM-3C+SC and GaMMA-STALTA+SC in all
12 cells at both caps, with the smallest margins $+0.06$ to $+0.07$ at 0.1; but on the same subsets
GaMMA-DY+SC, with the optimistic picks, exceeds MT-FW+SC at 0.1 in five of the six weak-clutter cells.
Fourth, the source is a point double couple, inside the moment-tensor subspace of MT-FW; the mismatch
lies in the medium, the wavelet, the coda, the polarization and the amplitude factors, not in the source model.
Fifth, the data are synthetic: one generator family, white noise, a one-dimensional gradient medium,
two seeds and 68 to 87 events per cell.

### Y.18a What Task 2A establishes and what it does not

Task 2A establishes, under a design and rules registered before the evaluation runs, that on this
synthetic generator, with a medium, wavelet, coda and fixed station delays that the detector does not
know, an operational picker and station corrections estimated from twenty catalogued events with
location errors, MT-FW with corrections detects more events at equal realized false-alarm rate than
the best of GaMMA, PyOcto and QuakeMigrate, in every cell of both seeds and both clutter regimes, by
$+0.09$ to $+0.49$, with the operating point of every method chosen on the evaluation windows; that
the gain depends on the geometry of the true hypothesis, since the same iteration with the subspace of
another node loses it; and that the recall advantage survives operating points frozen on other
networks, although the realized rate then exceeded twice the budget in one cell. It changes the status of the comparison with QuakeMigrate from
exploratory (Y.15) to confirmed on synthetics.

It does not isolate the radiation pattern. The control of rule 2A-b takes the whole geometry of another
node, travel times and radiation together, so passing it shows that the gain needs the true hypothesis,
not that it needs the moment-tensor subspace rather than the pick-free deflation that carries it. The
test that separates the two is the preregistered ablation Task 2B, in which the radiation operator is
replaced by the identity inside the same solver. Rule 2B holds in both clutter regimes, clearly with
strong clutter and only at the minimum of two of three cells with weak clutter, where much of the
advantage over QuakeMigrate remains without the operator; and because the identity removes the sign
of the pattern, the polarization constraint and the amplitude geometry together, it shows that the
operator matters, not that its sign does (Appendix Y.19). Nor does Task 2A
compare MT-FW with the methods closest to it: its comparators are pick-based associators and
QuakeMigrate, not block-sparse moment-tensor inversion with joint origin time and hypocentre, as in the
block orthogonal matching pursuit of [129] for simultaneous sources or
in grid-based moment-tensor monitoring such as GRiD MT [123], and not mechanism-corrected stacking such
as [112]. Task 2C ran two of them on the same windows, block orthogonal matching pursuit in the form
of [129] and mechanism-corrected stacking in the form of [122], whose statistic is up to normalisation
the variance reduction of grid-based moment-tensor monitoring, and both fell below MT-FW with
corrections in every cell (Appendix Y.20); the stacking of [112] remains untested.

It does not establish superiority over pick-based associators with good picks: with the optimistic
pick model and weak clutter GaMMA ties at 0.25 in the one-sided cells, is ahead at 0.1 in four of the
six cells of the full sample, and on single-source windows is ahead at 0.1 in five of the six, so the
advantage over an associator with good picks is established only with strong clutter. It does not establish that station corrections always help, since with strong clutter a
contaminated catalogue leaves residuals of 0.20 to 0.28 s. It does not give intervals that exclude zero
at 0.1 in four cells once the clustering by network is respected (for the comparison with QM-3C+SC,
Task 2C gives such intervals in every cell over thirty networks, Appendix Y.20), nor a deployable
false-alarm
guarantee, since the frozen threshold realized up to 2.3 times its budget in one cell. And it says
nothing about sources outside the double-couple class, real noise, three-dimensional media, MALMI
[110] or real data, which remain the tests of Section 16.

Nor does it test the declared operating regime of Section 1.5. Its wavelets have centre frequencies
of 1.4 to 2.8 Hz and its data are sampled at 16 Hz, against a primary band of 0.1 to 2 Hz; its sources
lie in a $30\times30$ km epicentral square at 5 to 13 km depth, against a primary domain of
$120\times120\times40$ km; its medium is a one-dimensional gradient. The result is therefore a
statement about a short-period, local configuration, and whether it holds in the declared regime is
open.

## Y.19 The radiation operator inside the same solver (Task 2B)

Code: `task2b_ablation.py` (runs), `analyze_2b.py` (rule and tables), `task2b_kernel_check.py`
(exploratory numerical sensitivity); output `analyze_2b_out.txt`; design and rule in
`PREREGISTRO_1j_1k.txt` (section "Task 2B"); controls, incident and limitations in
`task2b_report.txt`.

**Why the ablation is needed.** A pick-free method can beat pick-based associators simply because it
does not depend on picks, which the STA/LTA picker of Y.18 starves, and it can beat QuakeMigrate
because its greedy deflation explains detected events away and separates sources close in time that
a trigger merges. Neither advantage needs the radiation pattern. Rule 2A-b does not settle the point:
its control takes the subspace of another node, which changes travel times and radiation geometry
together, so passing it shows that the gain needs the true hypothesis, not that it needs the
moment-tensor subspace rather than the pick-free deflation that carries it (Y.18a). The novelty audit
raised exactly this objection. Task 2B answers it by removing the radiation operator and keeping
everything else fixed.

**Preregistration.** The design and the rule were written on 26 September 2026, after Task 2A had
been run, and appended to a local registry file; they were not modified. The time written in the
entry, 14:30 UTC, was typed without consulting a clock and is wrong, and the erratum is recorded in the
registry. The file system cannot date the text of the entry, because it records only the last
modification of the registry. The local transcript of the working session records the append at
14:14:04 UTC, before the earliest file of Task 2B, `task2b_ablation.py` (14:16:49 UTC), and its first
output, `task2b_selfcheck.txt` (14:17:01 UTC); that transcript is also local and is not an external
timestamp, and no hash fixes the text of the entry. Task 2B is therefore classified as registered in
an append-only local file with timing that is not externally verifiable. Task 2D is the first task
whose registration carries SHA-256 records: from it on, each registry entry carries a clock time and a
hash of the registry file, recorded in `PREREGISTRO_hashes.txt`. One edit of the analysis is declared
here. `analyze_2b.py` was created at 14:19:07 UTC, first run at 14:26:15 UTC and edited at 14:26:33 UTC,
eighteen seconds after its first output, to add one exploratory line per regime with the fraction of
bootstrap replicates in which rule 2B holds (0.59 and 0.62 with weak clutter, 0.98 and 0.97 with strong
clutter, quoted below); the output was rewritten at 14:26:43 UTC. The edit was not declared in
`task2b_report.txt`. A reconstruction of the script before the edit, run on the same records, gives an
output identical line by line to the reported one except for those two added lines, and an identical
JSON file, so the rule, the tables and the verdict did not change (E36 of the error register). Rule 2B, per clutter regime: MT-FW+SC minus MT-FW-FREE+SC is at least $0.05$ at both caps,
$\le0.1$ and $\le0.25$ realized false events per window, in at least two of three cells, in both
seeds, decided on point estimates. The registered interpretation: if the rule fails, the gain of
Task 2A is attributed to pick-free migration with deflation and not to the radiation pattern, and the
proposal says so. The endpoint and the intervals are those of Task 2A, with a bootstrap over networks
as an additional reading.

**Design.** The ablation runs on the same windows as Task 2A. For each seed (2024, 7331), regime,
cell and network, the networks, catalogues, station corrections and evaluation windows are rebuilt
by the functions of the Task 2A script itself, imported rather than copied, from the same random
streams. The picker, the optimistic pick model and the node map, which Task 2A also ran, do not draw
from the streams that the windows consume, so omitting them reproduces the window sequence exactly:
20 windows per network, 100 per cell, 1200 in all. MT-FW+SC is the stored Task 2A record.
MT-FW-FREE+SC is the same call to the certified solver, with the same grid of 675 nodes, the same
origin times, at most five atoms, the 5 km exclusion, the per-station and per-phase cap of Y.11 and
the same station corrections; the only change is that at every node the P and S radiation operator
$G$, of size $6n_s\times6$ for $n_s$ stations, and its basis are replaced by the identity of size
$6n_s$, so that each station and phase carries a free three-component amplitude. The least-squares
amplitude is then the vector of observed capped peak amplitudes itself, the statistic is half its
total energy, $\tfrac12\lVert\mathbf v\rVert^2$, which is invariant to the sign of any station's
vector, and the deflation subtracts the observed peak vectors convolved with the matched-filter
autocorrelation. MT-FW-FREE+SC is therefore the moment-tensor detector without the moment tensor:
the same search, cap, exclusion, corrections and explaining-away, scored by sign-free energy. The
endpoint, the exact FROC curve, the choice of operating point on the evaluation windows for both
methods and the paired bootstrap over windows (2000 replicates, the same resamples as Task 2A) are
unchanged; the bootstrap over networks draws the five networks of a cell with replacement (2000
replicates). The pairing of the two records was asserted window by window: truth, origin times and
the MT-FW+SC records are identical in the Task 2A and Task 2B files. The 60 network jobs took 495 s of
CPU time in all.

**Checks.** On one window with one source and six stations, the score of MT-FW-FREE+SC equals
$\tfrac12\lVert\mathbf v\rVert^2$ over all $675\times31$ hypotheses of node and origin time to a
relative $6.8\times10^{-16}$, its first selected atom is the maximum of that quantity, and the
least-squares amplitude reproduces the observed vector exactly. Reversing the sign of the waveform
of any one of the six stations leaves the MT-FW-FREE+SC records identical bit for bit in six of six
cases and changes the MT-FW+SC records in six of six, so the ablated statistic is blind to the sign
and the moment-tensor statistic is not; in that same window MT-FW+SC misses the event and
MT-FW-FREE+SC detects it, which shows that the ablation is not trivially weaker window by window. For every one of the 60
networks the station positions, statics, catalogue draws, catalogue locations and origin times and
station corrections are identical bit for bit to those of Task 2A, truth and origin times agree in
1200 of 1200 windows, MT-FW+SC recomputed on the first two windows of every network equals the
stored records exactly in 120 of 120 windows, and the certified MT-FW-FREE+SC equals its dense
evaluation exactly in 120 of 120.

The first run did not pass that control and stopped, as the design required, after 32 networks had
passed. The networks of the 6-station cell of seed 7331 written in the resumed sessions of Task 2A
failed it: MT-FW+SC scores differed in the last one or two digits, at a relative $6\times10^{-15}$ at
most, with identical hits and ordering, and in one network one station correction differed by
$8.9\times10^{-16}$ s. The cause is the numerical library, not the code. The OpenBLAS build used by
numpy selects its computational kernel at run time; on this machine the default is the Prescott
kernel, whereas the 18 networks that Task 2A wrote in its resumed seed-7331 sessions (the 4-station
cell in both regimes and the 6-station cell, all five networks with weak clutter and three with
strong clutter) had been computed with the AVX-512 SkylakeX kernel. With the SkylakeX kernel those
18 networks reproduce bit for bit, and with the Prescott, Haswell, Sandybridge and Zen kernels they
do not; the other 42 networks reproduce bit for bit with the default kernel. The final run therefore
pinned the kernel per network, SkylakeX for those 18 and the default for the rest, recorded it in
every part file it wrote, and every network passed the unchanged exact control; no threshold or
tolerance was relaxed. An exploratory recomputation of the 18 networks with the default kernel, 360
windows, measures what the kernel changes: the MT-FW+SC records differ at a relative
$6\times10^{-15}$ at most, with no change of hits or ordering in any window; the MT-FW-FREE+SC records
diverge in later Frank-Wolfe steps in some windows, where nearly equal candidates swap under rounding
(score differences up to 17 per cent), and hits or ordering change in 2 of the 360 windows; and the
recall at both caps is identical for both methods in all four affected cells. The per-network results
of Task 2A thus mix two kernels, with an effect on MT-FW at the level of rounding only.

**Table Y.34. Weak clutter. Recall of MT-FW+SC and MT-FW-FREE+SC at $\le0.1$ / $\le0.25$ realized
false events per window, and their difference, point [95% interval], with the paired bootstrap over
windows and the bootstrap over networks.** Differences are computed before rounding, so they can
differ by 0.01 from the rounded recalls; the one-sided cell of seed 7331 at 0.25 is $+0.0494$, printed
as $+0.05$, and does not reach the threshold. Script `analyze_2b.py`; output `analyze_2b_out.txt`,
Tables 1 and 2 (thin).

| Cell | events | MT-FW+SC | MT-FW-FREE+SC | windows, $\le0.1$ | windows, $\le0.25$ | networks, $\le0.1$ | networks, $\le0.25$ | cell passes |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| one-sided, 8 st. (2024) | 71 | 0.55/0.65 | 0.56/0.66 | -0.01 [-0.17, +0.21] | -0.01 [-0.14, +0.11] | -0.01 [-0.14, +0.37] | -0.01 [-0.13, +0.10] | no |
| 4 stations (2024) | 87 | 0.32/0.52 | 0.15/0.24 | +0.17 [+0.06, +0.37] | +0.28 [+0.13, +0.45] | +0.17 [+0.12, +0.25] | +0.28 [+0.16, +0.43] | yes |
| 6 stations (2024) | 72 | 0.46/0.56 | 0.28/0.44 | +0.18 [-0.06, +0.37] | +0.11 [-0.05, +0.28] | +0.18 [-0.02, +0.34] | +0.11 [-0.11, +0.26] | yes |
| one-sided, 8 st. (7331) | 81 | 0.64/0.77 | 0.60/0.72 | +0.04 [-0.13, +0.18] | +0.05 [-0.09, +0.18] | +0.04 [-0.21, +0.22] | +0.05 [-0.07, +0.19] | no |
| 4 stations (7331) | 78 | 0.32/0.51 | 0.17/0.27 | +0.15 [+0.00, +0.32] | +0.24 [+0.06, +0.38] | +0.15 [+0.01, +0.35] | +0.24 [+0.05, +0.39] | yes |
| 6 stations (7331) | 68 | 0.72/0.79 | 0.46/0.68 | +0.26 [+0.04, +0.46] | +0.12 [-0.02, +0.33] | +0.26 [+0.04, +0.39] | +0.12 [-0.02, +0.28] | yes |

**Table Y.35. Strong clutter. The same.** Same script and output, Tables 1 and 2 (strong).

| Cell | events | MT-FW+SC | MT-FW-FREE+SC | windows, $\le0.1$ | windows, $\le0.25$ | networks, $\le0.1$ | networks, $\le0.25$ | cell passes |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| one-sided, 8 st. (2024) | 83 | 0.31/0.42 | 0.14/0.20 | +0.17 [+0.01, +0.37] | +0.22 [+0.07, +0.37] | +0.17 [-0.04, +0.48] | +0.22 [-0.00, +0.49] | yes |
| 4 stations (2024) | 77 | 0.31/0.39 | 0.06/0.13 | +0.25 [+0.08, +0.35] | +0.26 [+0.13, +0.39] | +0.25 [+0.14, +0.29] | +0.26 [+0.19, +0.37] | yes |
| 6 stations (2024) | 74 | 0.38/0.49 | 0.05/0.14 | +0.32 [+0.19, +0.48] | +0.35 [+0.23, +0.51] | +0.32 [+0.25, +0.45] | +0.35 [+0.29, +0.48] | yes |
| one-sided, 8 st. (7331) | 75 | 0.63/0.71 | 0.29/0.37 | +0.33 [+0.18, +0.51] | +0.33 [+0.18, +0.49] | +0.33 [+0.24, +0.41] | +0.33 [+0.20, +0.44] | yes |
| 4 stations (7331) | 78 | 0.13/0.22 | 0.06/0.10 | +0.06 [-0.02, +0.20] | +0.12 [+0.01, +0.24] | +0.06 [+0.00, +0.16] | +0.12 [+0.04, +0.18] | yes |
| 6 stations (7331) | 84 | 0.36/0.42 | 0.10/0.25 | +0.26 [+0.06, +0.39] | +0.17 [+0.01, +0.35] | +0.26 [+0.00, +0.35] | +0.17 [-0.03, +0.37] | yes |

**Table Y.36. Exploratory. MT-FW-FREE+SC minus the best comparator of rule 2A-a, with 95% paired
bootstrap intervals over windows (maximum re-chosen in each replicate).** The best comparator is a
QuakeMigrate variant in every case. The corresponding margins of MT-FW+SC are in Table Y.31. Same
script and output, Table 3.

| Cell | weak, $\le0.1$ | weak, $\le0.25$ | strong, $\le0.1$ | strong, $\le0.25$ |
|:--|:--|:--|:--|:--|
| one-sided, 8 st. (2024) | +0.21 [+0.00, +0.37] | +0.28 [+0.14, +0.41] | -0.01 [-0.20, +0.08] | +0.01 [-0.11, +0.13] |
| 4 stations (2024) | -0.03 [-0.17, +0.08] | -0.01 [-0.16, +0.10] | -0.05 [-0.14, +0.00] | -0.01 [-0.11, +0.05] |
| 6 stations (2024) | +0.00 [-0.13, +0.21] | +0.14 [-0.02, +0.29] | -0.07 [-0.18, +0.00] | -0.03 [-0.18, +0.08] |
| one-sided, 8 st. (7331) | +0.14 [+0.00, +0.29] | +0.23 [+0.11, +0.36] | +0.04 [-0.09, +0.15] | +0.08 [-0.05, +0.19] |
| 4 stations (7331) | +0.01 [-0.06, +0.15] | +0.10 [-0.01, +0.20] | +0.03 [-0.04, +0.08] | +0.05 [-0.04, +0.11] |
| 6 stations (7331) | +0.22 [+0.02, +0.40] | +0.34 [+0.11, +0.46] | -0.08 [-0.20, +0.04] | +0.00 [-0.12, +0.13] |

**Rule 2B holds in both regimes, in strong clutter clearly and in weak clutter at the minimum.** With
strong clutter the radiation operator adds $+0.06$ to $+0.35$ in all three cells of both seeds (Table
Y.35), and eleven of the twelve intervals over windows exclude zero; with the coarser bootstrap over
five networks, three intervals include zero (the one-sided cell of seed 2024 at both caps and the
6-station cell of seed 7331 at 0.25) and two more have a lower bound of $+0.00$. With weak clutter the rule passes with two of three cells in
each seed, the minimum it allows (Table Y.34): the operator adds $+0.15$ to $+0.28$ in the 4-station
cells, where the network intervals exclude zero at 0.25, and $+0.11$ to $+0.26$ in the 6-station
cells, where six of the eight intervals include zero, while in the one-sided 8-station cell
MT-FW-FREE+SC equals MT-FW+SC in both seeds ($-0.01$ and $-0.01$ in seed 2024, $+0.04$ and $+0.0494$ in
seed 7331). The verdict is decided on point estimates, as registered, and its stability was read
afterwards by re-applying the rule in every bootstrap replicate, with independent resamples per cell
(exploratory): it holds in 0.98 of the replicates over windows and 0.97 over networks with strong
clutter, but in only 0.59 and 0.62 with weak clutter. By the registered interpretation, the gain of
Task 2A is not attributed solely to pick-free migration with deflation: within the same solver the
radiation operator adds at least $0.05$ at both caps in most cells. With weak clutter that
attribution is fragile.

**What remains without the operator (exploratory).** Table Y.36 compares MT-FW-FREE+SC with the
comparators of Task 2A, an analysis not registered for Task 2B. With weak clutter the operator-free
detector exceeds the best comparator at 0.25 in five of six cells, by $+0.10$ to $+0.34$, and falls
$0.01$ short of it only in the 4-station cell of seed 2024; in the one-sided cells it is about as far
ahead as MT-FW+SC ($+0.28$ and $+0.23$, against $+0.27$ and $+0.28$ in Table Y.31), in the 6-station
cells it keeps $+0.14$ and $+0.34$ of $+0.25$ and $+0.46$, and in the 4-station cells little or nothing
($-0.01$ and $+0.10$ of $+0.26$ and $+0.35$). Much of the weak-clutter advantage of Task 2A over
QuakeMigrate therefore comes from pick-free migration with deflation and station corrections and does
not need the radiation pattern; in weak clutter the operator matters most in the sparsest
networks. With strong clutter the picture reverses: MT-FW-FREE+SC is level with the best comparator,
between $-0.08$ and $+0.08$ at both caps, where MT-FW+SC was ahead by $+0.09$ to $+0.41$, so the
strong-clutter result of Task 2A depends on the operator. The mechanism-independent migration
WAVE-PT@20, which is sign-free but constrains P to the ray and S to the transverse plane and has
neither deflation nor station corrections, gives a further hint: with strong clutter
MT-FW-FREE+SC falls below it in every cell at both caps, by $0.01$ to $0.23$ (with an interval over
windows that excludes zero only in the 6-station cell of seed 2024), while with weak clutter the
difference ranges from $-0.06$ to $+0.19$. That pattern suggests that the polarization constraint is
a substantial part of what the operator contributes against strong clutter, whose transients are
randomly polarized; it is a suggestion, not a measurement.

**What Task 2B does not isolate.** Replacing $G$ by the identity removes the whole physical operator
at once: the signed double-couple pattern across stations, the polarization constraint (P along the
ray, S in the transverse plane) and the amplitude geometry of geometric spreading and the P-to-S
ratio. Task 2B therefore establishes that the operator matters inside the solver; it does not isolate
the sign of the radiation pattern, which was the question of the audit, and neither does rule 2A-b,
whose control scrambles geometry and polarization together. Separating the sign from the polarization
would need a further ablation that keeps the polarization constraint and frees the amplitude, for
example a free amplitude along the ray for P and in the transverse plane for S; it is not registered,
and any such run is exploratory. Four further limits apply. The cells rest on five networks each, the
intervals over windows ignore that clustering (Task 2C replicated over thirty networks the comparison
with QuakeMigrate, not this ablation; Appendix Y.20), and several individual cells have intervals that
include zero where the point estimate passes. The rule is decided on point estimates, and weak
clutter passes with the minimum count. The operating point is chosen on the evaluation windows for
both methods, as in Task 2A, which is symmetric. And MT-FW-FREE+SC inherits the choices of the solver
(the cap of 6 per station and phase, the 5 km exclusion, at most five atoms and the deflation with the
matched-filter autocorrelation), which were tuned with the moment-tensor statistic; an energy
detector tuned on its own might do somewhat better, and the registered ablation deliberately keeps the
solver fixed.

## Y.20 The closest comparators, the convex problem, and the network as the unit (Task 2C)

Code: `task2c_methods.py` (the three methods), `task2c_run.py` (runs and controls), `task2c_validate.py`
(validation on seed 999), `analyze_2c.py` (rules and tables), `task2c_explore.py` (exploratory
diagnostics); outputs `task2c_validation.txt`, `task2c_smoke.txt`, `analyze_2c_out.txt`,
`task2c_explore_all.txt` and `task2c_explore_e3.txt`; design and rules in `PREREGISTRO_1j_1k.txt`
(section "Task 2C"); review of the code, controls, compute and deviations in `task2c_report.txt`.

**Why the task was run.** An independent audit of this proposal, carried out after Task 2B, named three
weaknesses that Tasks 2A and 2B could not answer. The two methods closest to MT-FW had not been run:
block-sparse moment-tensor inversion by block orthogonal matching pursuit [129], of which MT-FW is in
its algorithmic core a grid relative (Section 7.5), and mechanism-corrected stacking [122] (Y.18a). The
convex problem that the thesis poses, the group BLASSO (16a), had not been implemented, so every result
concerned its greedy approximation MT-FW, which the theory of Appendix AA does not cover beyond
Proposition AA.4(e)-(h). And the intervals of Tasks 2A and 2B resample windows and ignore that the
windows of a cell share five networks, although Section 14 declares the network the unit of
replication. Task 2C answers the three in two parts: part (i) runs the two comparators and the group
BLASSO on the grid over the 1200 windows of Task 2A, and part (ii) adds twenty-five networks per cell to
replicate the comparison with QuakeMigrate over thirty networks.

**Preregistration.** The design, three rules and their interpretation were written on 26 September 2026,
before the methods were implemented and before any quantity of the task was computed, and appended to
the local registry file. The time in the entry, 15:05 UTC, was typed rather than read from a clock, and
the file system cannot date the text of the entry; the local transcript of the working session records
the append at 15:04:21 UTC, before the first code file of the task, `task2c_methods.py`, was written
at about 15:12 UTC, but that transcript is also local and no hash fixes the text. Task 2C is therefore
classified, like Task 2B, as registered in an append-only local file with timing that is not
externally verifiable.
Rule 2C-a, per clutter regime: MT-FW+SC minus the better of BOMP+SC and MTSTACK+SC is at least $0.05$ at
$\le0.25$ realized false events per window in at least two of three cells, in both seeds; if it fails,
the advantage of MT-FW is attributed to the moment tensor with deflation, already present in [129],
and not to the elements proper to MT-FW, its peak search and its per-station cap. Rule 2C-b is
descriptive and has no threshold: GL-FW+SC minus MT-FW+SC at both caps, per cell, with intervals,
reporting whether the convex method equals its greedy approximation (absolute difference below
$0.05$), exceeds it or falls below it. Rule 2C-c, primary for part (ii), per regime: the lower bound of
the 95% hierarchical bootstrap interval of MT-FW+SC minus QM-3C+SC at $\le0.25$ is above zero, and the
point estimate at least $0.05$, in at least two of three cells in both seeds; the same differences
against BOMP+SC and MTSTACK+SC, and at $\le0.1$, are reported with Holm's correction [142] over the
comparisons of each regime. A budget rule would have cut part (ii) to fifteen additional networks per cell had
the seed-999 estimate exceeded 8 hours on two cores. A dated note added at 19:56 UTC, before any run
of seeds 2024 or 7331, records that the smoke test on seed 999 estimated 0.79 hours on two cores, so
the rule did not trigger and part (ii) used twenty-five networks; nothing else in the registration was
changed. The registration was written after Tasks 2A and 2B had been run on the same windows, and its
design was informed by them; what it fixed in advance are the comparators, the rules and their
interpretation.

**Design.** Part (i) regenerates the 1200 windows of Task 2A (two seeds, two regimes, three cells, five
networks, twenty windows) bit for bit with the functions of the Task 2A script, through the same
imports as Task 2B, and uses the stored Task 2A records of MT-FW+SC and QM-3C+SC. Part (ii) adds
twenty-five networks per cell, seed and regime, with indices 5 to 29 and the seed construction of Task
2A for networks, statics, catalogues and windows, and evaluates windows 0 to 11 of each with MT-FW+SC,
QM-3C+SC, BOMP+SC and MTSTACK+SC; GaMMA and PyOcto are not run, because the picker starved them in Task
2A. For rule 2C-c, networks 0 to 4 contribute the records of their 20 windows from Task 2A and part
(i), and networks 5 to 29 those of their 12, which gives thirty networks and 400 windows per cell.
Every method receives the station corrections of Task 2A (the estimator of Y.16 on a catalogue of 20
events per network and regime), searches the same grid of 675 nodes at 3 km spacing and 31 origin
times at 0.1 s steps, excludes nodes within 5 km of an accepted event, emits at most five events per
window, and is scored by the hit rule of Task 2A.

**The three methods as implemented.** BOMP+SC, MTSTACK+SC and GL-FW+SC share one dictionary in the
matched-filter domain. A block is a pair of node and origin time with six columns, one per
moment-tensor component; column $k$ predicts, at every station and phase, the $k$-th column of the
radiation operator $G$ times the matched-filter autocorrelation centred at the exact predicted sample,
the predicted time including the station correction, with no peak search and no cap. This is the
matched-filter image of the data-domain dictionary whose atoms are the columns of $G$ carried by the
detector's wavelet at the predicted times, so the correlation of a block with the data is the
normalized matched-filter output sampled at the predicted times and contracted with $G$, and the Gram
matrix of two blocks sums, over stations and phase pairs, the products of radiation coefficients
weighted by the autocorrelation at the difference of their predicted times, including the overlap of
P and S at one station. Deflation subtracts a block's image from the matched-filter output. MT-FW
differs from all three in the elements that rule 2C-a names: it reads each station and phase at the
peak of the matched filter within $\pm0.25$ s of the predicted time, and it caps the per-station
vectors (Y.11).

BOMP+SC is block orthogonal matching pursuit in the form of [129]. At each step it selects the block
with the largest energy of the projection of the residual onto the block,
$e=\mathbf h^{\mathsf T}K_{bb}^{+}\mathbf h$, with $\mathbf h$ the correlation of the block with the
residual and $K_{bb}$ its Gram block; it scores the block by $e/2$ at entry, the convention of WAVE-MT
and MT-FW, made monotone by a running minimum; and it refits all active blocks jointly by least squares
before recomputing the residual. MTSTACK+SC is mechanism-corrected stacking in the form of [122]: at
every node and origin time it inverts the moment tensor by least squares from the amplitudes at the
exact predicted samples, takes half the energy of the prediction as the stack value, maximizes over
origin time and forms the catalogue by greedy suppression at 20 km, with no deflation; it equals the
single-pass statistic WAVE-MT of Y.11 read at the exact samples, and, up to normalisation, it is the
variance reduction of a moment-tensor fit at every node on which grid-based monitoring detects [132],
[123]. GL-FW+SC is (16a) restricted to that grid:

$$\min_{\mathbf u}\ \tfrac12\lVert\mathbf y-A\mathbf u\rVert^2+\lambda\sum_b\lVert A_b\mathbf u_b\rVert,$$

whose group norm is the block-whitened norm that (16a) prescribes, so that the dual certificate of a
block, $\sqrt{\mathbf h^{\mathsf T}K_{bb}^{+}\mathbf h}$, is exactly the selection statistic of BOMP and
$\lambda_{\max}$ is its maximum over blocks. It is solved by fully corrective Frank-Wolfe: a block enters
the active set whenever its certificate exceeds $\lambda(1+10^{-6})$, and the problem restricted to the
active blocks is solved exactly, by block-coordinate descent with the exact group soft-threshold
followed by damped Newton steps on the non-zero groups, to a relative violation of the optimality
conditions of at most $10^{-9}$; the duality gap is evaluated at the scaled dual point of Proposition
AA.4(d). The path runs over 30 geometric values of $\lambda$ from $\lambda_{\max}$ to $\lambda_{\max}/30$
with warm starts. The score of a block is the largest $\lambda$ of the path at which its group norm is
positive, and the record lists blocks by decreasing entry value, ties broken by the larger group norm
at entry, with the 5 km exclusion and at most five events.

**Validation.** Before any evaluation seed was run, `task2c_validate.py` checked the implementation on
seed 999 (`task2c_validation.txt`). The dictionary identities hold to rounding: over 20 blocks per cell,
including pairs with overlapping P and S, the correlations and the Gram matrix built from explicit
data-domain atoms agree with the matched-filter forms to a relative $8.7\times10^{-16}$ and
$7.3\times10^{-16}$, and every diagonal Gram block has rank 6. For a noise-free source drawn from the
dictionary itself BOMP recovers the block in 12 of 12 cases, with the moment tensor to a relative
$6.5\times10^{-13}$ and a second score at most $3.3\times10^{-30}$ of the first, and GL-FW enters the
true block first in 12 of 12. For a physical waveform from a source at a node with its origin on the
grid, node and origin are recovered in 12 of 12, but the moment tensor only approximately, with cosines
of at least 0.978, 0.982 and 0.874 in the 4-station, 6-station and one-sided cells, because a pulse up
to one sample off the predicted sample distorts the amplitudes read there; with the arrivals moved to
the predicted samples, one minus the cosine is at most $2.2\times10^{-16}$. On 24 windows of seed 999,
four per regime and cell, solved two values of $\lambda$ beyond the early stop described under
deviations, the largest certificate of GL-FW off the support is $1.000000000\,\lambda$ against the
admission tolerance $1+10^{-6}$, the certificate on the support equals $\lambda$ to $5.1\times10^{-10}$,
the relative duality gap is at most $1.5\times10^{-9}$, no restricted problem had to be re-solved, the
early-stopped record equals the record of the longer path in 24 of 24 windows, and
$\lambda_{\max}^2/2$ equals the first score of BOMP to a relative $2.1\times10^{-16}$. MTSTACK equals the
projection energy of WAVE-MT at the exact samples over $3\times675\times31$ hypotheses per regime and
cell to a relative $1.5\times10^{-13}$. Across the 1200 windows of part (i), GL-FW stopped with a
complete record in every window after 2 to 14 solved values of $\lambda$ (median 3); the largest
relative gap was $1.4\times10^{-8}$, the largest certificate off the support $1.000000000\,\lambda$ and
the largest deviation on it $8.7\times10^{-10}$, with no re-solve and a final support of median 11 to
12 blocks and at most 47.

**Controls and kernels.** For all 60 networks of part (i), the station positions, statics, catalogue
draws, catalogue locations and origin times and station corrections are identical bit for bit to those
of Task 2A; truth, origin times and the number of sources agree in 1200 of 1200 windows; and MT-FW+SC
recomputed on windows 0 and 1 of every network equals the stored Task 2A records exactly in 120 of 120
windows. In part (ii) the certified MT-FW+SC equals its dense evaluation exactly in 600 of 600 windows,
and on seed 999 the part-(ii) pipeline reproduced the stored Task 2A records of MT-FW+SC and QM-3C+SC on
windows 0 and 1 in 6 of 6 cells. The computational kernel of the numerical library, which Y.19 found to
change the last digits of MT-FW, was pinned per process. The default kernel of numpy's OpenBLAS on this
machine had changed from Prescott, the default when Task 2B was run, to SkylakeX, so an unset
variable would have run under the wrong kernel the 42 networks that Task 2A computed with Prescott. The
18 networks listed in Y.19 ran with SkylakeX, and every other network, including all of part (ii) and
seed 999, with Prescott; each part file records the kernel that numpy loaded. Under Prescott the
separate OpenBLAS of scipy selects a different kernel (Katmai) from the one it used in Task 2A; scipy's
linear algebra is not used by MT-FW, by the station corrections or by the three methods, and every exact
control passed. No control failed.

**Statistics.** The endpoint is that of Task 2A: event recall at $\le0.1$ and $\le0.25$ realized false
events per window on the exact FROC curve, with the operating point of every method chosen on the
evaluation windows. For rules 2C-a and 2C-b the intervals come from the paired bootstrap over windows
of Task 2A, with the same resamples, and from a bootstrap over the five networks of a cell, 2000
replicates each, with the operating point and the better comparator re-chosen in each replicate. For
rule 2C-c the bootstrap is hierarchical: 2000 replicates, each drawing thirty networks with replacement
and then the windows of each drawn network with replacement, with the false-event rate computed over
the resampled windows and the operating point re-chosen in each replicate; the point estimate pools the
400 windows of the cell. Holm's procedure uses one-sided bootstrap p-values,
$(1+\#\{d^\ast\le0\})/2001$, at level 0.025, the level of the lower bound of a two-sided 95% interval,
over a family of 36 comparisons per regime (three comparators, two caps, three cells, two seeds); the
smallest attainable p-value is $1/2001$, which becomes $0.018$ at the first Holm step. Holm's procedure
controls the familywise error and is more conservative than the Benjamini-Hochberg procedure that
Section 14 registers for the confirmatory family. The compute was small: part (i) took 318 s of CPU
time for 60 networks, with 0.04 s per window for BOMP, 0.013 s for MTSTACK and 0.15 s for GL-FW (0.64 s
at most), and part (ii) took 5342 s for 300 networks, with 0.31 s per window for MT-FW+SC and 1.11 s for
QM-3C+SC; 1.57 hours of CPU time in all.

**Table Y.37. Part (i). Event recall at realized false events per window $\le0.1$ / $\le0.25$ on the
1200 windows of Task 2A (100 per cell).** MT-FW+SC and QM-3C+SC are the Task 2A records, printed to two
decimals in Tables Y.29 and Y.30. Script `analyze_2c.py`; output `analyze_2c_out.txt`, Table 1 of each
regime.

| Regime | Cell (seed) | events | MT-FW+SC | BOMP+SC | MTSTACK+SC | GL-FW+SC | QM-3C+SC |
|:--|:--|:--|:--|:--|:--|:--|:--|
| weak | 4 stations (2024) | 87 | 0.322/0.517 | 0.161/0.230 | 0.092/0.161 | 0.057/0.207 | 0.184/0.253 |
| weak | 6 stations (2024) | 72 | 0.458/0.556 | 0.083/0.222 | 0.069/0.139 | 0.083/0.125 | 0.278/0.306 |
| weak | one-sided, 8 st. (2024) | 71 | 0.549/0.648 | 0.197/0.380 | 0.070/0.225 | 0.070/0.197 | 0.352/0.380 |
| weak | 4 stations (7331) | 78 | 0.321/0.513 | 0.154/0.282 | 0.077/0.154 | 0.077/0.218 | 0.154/0.167 |
| weak | 6 stations (7331) | 68 | 0.721/0.794 | 0.176/0.353 | 0.118/0.265 | 0.118/0.235 | 0.235/0.338 |
| weak | one-sided, 8 st. (7331) | 81 | 0.642/0.765 | 0.160/0.247 | 0.123/0.173 | 0.086/0.160 | 0.469/0.481 |
| strong | 4 stations (2024) | 77 | 0.312/0.390 | 0.195/0.260 | 0.078/0.182 | 0.039/0.182 | 0.117/0.130 |
| strong | 6 stations (2024) | 74 | 0.378/0.486 | 0.095/0.162 | 0.081/0.135 | 0.095/0.135 | 0.122/0.162 |
| strong | one-sided, 8 st. (2024) | 83 | 0.313/0.422 | 0.108/0.217 | 0.096/0.157 | 0.084/0.157 | 0.157/0.193 |
| strong | 4 stations (7331) | 78 | 0.128/0.218 | 0.103/0.141 | 0.103/0.115 | 0.115/0.154 | 0.038/0.051 |
| strong | 6 stations (7331) | 84 | 0.357/0.417 | 0.071/0.214 | 0.036/0.119 | 0.024/0.107 | 0.179/0.250 |
| strong | one-sided, 8 st. (7331) | 75 | 0.627/0.707 | 0.067/0.253 | 0.080/0.173 | 0.080/0.133 | 0.253/0.293 |

**Table Y.38. Rule 2C-a. MT-FW+SC minus the better of BOMP+SC and MTSTACK+SC, point [95% interval over
windows] [95% interval over networks].** The better comparator is BOMP+SC in every case except the one
marked. The rule reads the cap $\le0.25$. Same output, Table 2.

| Regime | Cell (seed) | $\le0.1$: point [windows] [networks] | $\le0.25$: point [windows] [networks] | cell passes |
|:--|:--|:--|:--|:--|
| weak | 4 stations (2024) | +0.161 [+0.053, +0.347] [+0.045, +0.270] | +0.287 [+0.119, +0.430] [+0.159, +0.453] | yes |
| weak | 6 stations (2024) | +0.375 [+0.205, +0.500] [+0.211, +0.466] | +0.333 [+0.145, +0.481] [+0.000, +0.500] | yes |
| weak | one-sided, 8 st. (2024) | +0.352 [+0.208, +0.537] [+0.254, +0.478] | +0.268 [+0.111, +0.456] [+0.136, +0.452] | yes |
| weak | 4 stations (7331) | +0.167 [+0.000, +0.390] [+0.047, +0.387] | +0.231 [+0.022, +0.422] [+0.037, +0.356] | yes |
| weak | 6 stations (7331) | +0.544 [+0.321, +0.686] [+0.404, +0.597] | +0.441 [+0.194, +0.628] [+0.138, +0.613] | yes |
| weak | one-sided, 8 st. (7331) | +0.481 [+0.312, +0.634] [+0.250, +0.645] | +0.519 [+0.355, +0.674] [+0.355, +0.711] | yes |
| strong | 4 stations (2024) | +0.117 [-0.045, +0.238] [-0.024, +0.181] | +0.130 [+0.000, +0.260] [+0.069, +0.253] | yes |
| strong | 6 stations (2024) | +0.284 [+0.120, +0.447] [+0.155, +0.451] | +0.324 [+0.167, +0.471] [+0.182, +0.507] | yes |
| strong | one-sided, 8 st. (2024) | +0.205 [+0.050, +0.324] [+0.038, +0.367] | +0.205 [+0.051, +0.365] [+0.081, +0.366] | yes |
| strong | 4 stations (7331) | +0.026 [-0.074, +0.156] [-0.064, +0.155] | +0.077 [-0.038, +0.194] [-0.051, +0.198] | yes |
| strong | 6 stations (7331) | +0.286 [+0.091, +0.418] [+0.152, +0.379] | +0.202 [+0.053, +0.381] [+0.043, +0.345] | yes |
| strong | one-sided, 8 st. (7331) | +0.547 [+0.380, +0.695] [+0.386, +0.633] (MTSTACK+SC) | +0.453 [+0.279, +0.600] [+0.270, +0.589] | yes |

**Table Y.39. Rule 2C-b. GL-FW+SC minus MT-FW+SC, point [95% interval over windows] and registered
reading (equal if the absolute difference is below 0.05), and GL-FW+SC minus BOMP+SC (point).** The
intervals over networks are in the output and exclude zero wherever those over windows do. Same output,
Table 3.

| Regime | Cell (seed) | GL-FW+SC minus MT-FW+SC, $\le0.1$ | GL-FW+SC minus MT-FW+SC, $\le0.25$ | GL-FW+SC minus BOMP+SC, $\le0.1$ / $\le0.25$ |
|:--|:--|:--|:--|:--|
| weak | 4 stations (2024) | -0.264 [-0.432, -0.132] below | -0.310 [-0.506, -0.196] below | -0.103 / -0.023 |
| weak | 6 stations (2024) | -0.375 [-0.537, -0.227] below | -0.431 [-0.582, -0.261] below | +0.000 / -0.097 |
| weak | one-sided, 8 st. (2024) | -0.479 [-0.641, -0.328] below | -0.451 [-0.592, -0.301] below | -0.127 / -0.183 |
| weak | 4 stations (7331) | -0.244 [-0.432, -0.075] below | -0.295 [-0.446, -0.118] below | -0.077 / -0.064 |
| weak | 6 stations (7331) | -0.603 [-0.734, -0.423] below | -0.559 [-0.681, -0.423] below | -0.059 / -0.118 |
| weak | one-sided, 8 st. (7331) | -0.556 [-0.700, -0.395] below | -0.605 [-0.717, -0.438] below | -0.074 / -0.086 |
| strong | 4 stations (2024) | -0.273 [-0.366, -0.072] below | -0.208 [-0.358, -0.076] below | -0.156 / -0.078 |
| strong | 6 stations (2024) | -0.284 [-0.444, -0.129] below | -0.351 [-0.469, -0.208] below | +0.000 / -0.027 |
| strong | one-sided, 8 st. (2024) | -0.229 [-0.368, -0.108] below | -0.265 [-0.410, -0.145] below | -0.024 / -0.060 |
| strong | 4 stations (7331) | -0.013 [-0.165, +0.066] equal | -0.064 [-0.185, +0.056] below | +0.013 / +0.013 |
| strong | 6 stations (7331) | -0.333 [-0.439, -0.167] below | -0.310 [-0.494, -0.179] below | -0.048 / -0.107 |
| strong | one-sided, 8 st. (7331) | -0.547 [-0.693, -0.413] below | -0.573 [-0.694, -0.421] below | +0.013 / -0.120 |

**Table Y.40. Rule 2C-c. Thirty networks and 400 windows per cell. Recall at $\le0.25$, and MT-FW+SC
minus QM-3C+SC, point [95% hierarchical interval]; Holm-adjusted p-value.** "n.s." marks a comparison
that Holm's procedure does not reject at level 0.025. Output Tables 4 and 5.

| Regime | Cell (seed) | events | MT-FW+SC | QM-3C+SC | BOMP+SC | MTSTACK+SC | minus QM-3C+SC, $\le0.1$ | minus QM-3C+SC, $\le0.25$ |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| weak | 4 stations (2024) | 308 | 0.562 | 0.214 | 0.253 | 0.146 | +0.302 [+0.139, +0.415]; 0.0180 | +0.347 [+0.248, +0.453]; 0.0180 |
| weak | 6 stations (2024) | 313 | 0.623 | 0.300 | 0.275 | 0.144 | +0.236 [+0.073, +0.356]; 0.0180 | +0.323 [+0.228, +0.416]; 0.0180 |
| weak | one-sided, 8 st. (2024) | 314 | 0.656 | 0.373 | 0.290 | 0.194 | +0.207 [+0.096, +0.359]; 0.0180 | +0.283 [+0.189, +0.376]; 0.0180 |
| weak | 4 stations (7331) | 307 | 0.466 | 0.160 | 0.270 | 0.166 | +0.199 [+0.096, +0.287]; 0.0180 | +0.306 [+0.194, +0.401]; 0.0180 |
| weak | 6 stations (7331) | 308 | 0.685 | 0.338 | 0.299 | 0.208 | +0.289 [+0.170, +0.426]; 0.0180 | +0.347 [+0.242, +0.439]; 0.0180 |
| weak | one-sided, 8 st. (7331) | 326 | 0.638 | 0.414 | 0.304 | 0.190 | +0.166 [+0.048, +0.292]; 0.0180 | +0.224 [+0.119, +0.328]; 0.0180 |
| strong | 4 stations (2024) | 325 | 0.249 | 0.074 | 0.169 | 0.117 | +0.120 [+0.063, +0.192]; 0.0180 | +0.175 [+0.104, +0.254]; 0.0180 |
| strong | 6 stations (2024) | 304 | 0.431 | 0.184 | 0.197 | 0.138 | +0.201 [+0.093, +0.303]; 0.0180 | +0.247 [+0.150, +0.337]; 0.0180 |
| strong | one-sided, 8 st. (2024) | 305 | 0.489 | 0.252 | 0.252 | 0.203 | +0.144 [+0.063, +0.255]; 0.0180 | +0.236 [+0.135, +0.334]; 0.0180 |
| strong | 4 stations (7331) | 313 | 0.208 | 0.109 | 0.150 | 0.134 | +0.067 [+0.012, +0.147]; 0.0980 (n.s.) | +0.099 [+0.029, +0.181]; 0.0440 (n.s.) |
| strong | 6 stations (7331) | 329 | 0.459 | 0.228 | 0.191 | 0.119 | +0.182 [+0.081, +0.321]; 0.0180 | +0.231 [+0.125, +0.339]; 0.0180 |
| strong | one-sided, 8 st. (7331) | 329 | 0.556 | 0.261 | 0.246 | 0.164 | +0.231 [+0.116, +0.337]; 0.0180 | +0.295 [+0.212, +0.383]; 0.0180 |

**Table Y.41. MT-FW+SC minus the closest antecedents over thirty networks per cell, point [95%
hierarchical interval].** "n.s." as in Table Y.40. Output Table 5.

| Regime | Cell (seed) | minus BOMP+SC, $\le0.1$ | minus BOMP+SC, $\le0.25$ | minus MTSTACK+SC, $\le0.1$ | minus MTSTACK+SC, $\le0.25$ |
|:--|:--|:--|:--|:--|:--|
| weak | 4 stations (2024) | +0.334 [+0.173, +0.432] | +0.308 [+0.185, +0.414] | +0.386 [+0.228, +0.485] | +0.416 [+0.292, +0.525] |
| weak | 6 stations (2024) | +0.393 [+0.222, +0.488] | +0.348 [+0.218, +0.466] | +0.431 [+0.257, +0.531] | +0.479 [+0.352, +0.585] |
| weak | one-sided, 8 st. (2024) | +0.398 [+0.288, +0.500] | +0.366 [+0.259, +0.466] | +0.475 [+0.357, +0.570] | +0.462 [+0.351, +0.552] |
| weak | 4 stations (7331) | +0.147 [+0.029, +0.253] | +0.195 [+0.061, +0.298] | +0.248 [+0.136, +0.337] | +0.300 [+0.174, +0.382] |
| weak | 6 stations (7331) | +0.422 [+0.271, +0.541] | +0.386 [+0.267, +0.482] | +0.477 [+0.341, +0.600] | +0.477 [+0.347, +0.584] |
| weak | one-sided, 8 st. (7331) | +0.380 [+0.265, +0.506] | +0.334 [+0.224, +0.441] | +0.436 [+0.313, +0.536] | +0.448 [+0.317, +0.534] |
| strong | 4 stations (2024) | +0.065 [+0.003, +0.141] (n.s.) | +0.080 [+0.006, +0.165] (n.s.) | +0.108 [+0.037, +0.193] | +0.132 [+0.060, +0.203] |
| strong | 6 stations (2024) | +0.207 [+0.117, +0.328] | +0.234 [+0.118, +0.330] | +0.250 [+0.152, +0.347] | +0.293 [+0.183, +0.376] |
| strong | one-sided, 8 st. (2024) | +0.252 [+0.146, +0.352] | +0.236 [+0.143, +0.343] | +0.252 [+0.176, +0.364] | +0.285 [+0.190, +0.389] |
| strong | 4 stations (7331) | +0.035 [-0.022, +0.124] (n.s.) | +0.058 [-0.028, +0.136] (n.s.) | +0.048 [-0.003, +0.142] (n.s.) | +0.073 [+0.006, +0.155] (n.s.) |
| strong | 6 stations (7331) | +0.258 [+0.151, +0.384] | +0.267 [+0.169, +0.366] | +0.292 [+0.192, +0.406] | +0.340 [+0.238, +0.428] |
| strong | one-sided, 8 st. (7331) | +0.325 [+0.214, +0.437] | +0.310 [+0.211, +0.403] | +0.365 [+0.257, +0.463] | +0.392 [+0.303, +0.476] |

**Rule 2C-a holds in both regimes.** MT-FW+SC exceeds the better of BOMP+SC and MTSTACK+SC at $\le0.25$
in all three cells of both seeds, by $+0.23$ to $+0.52$ with weak clutter and by $+0.08$ to $+0.45$ with
strong clutter, and the better comparator there is BOMP+SC in every cell (Table Y.38). No difference
lies within 0.005 of the threshold. Re-applied in every replicate of the bootstrap over windows
(exploratory), the rule holds in 100 per cent of the replicates with weak clutter and 99 per cent with
strong clutter. At $\le0.25$ both intervals include zero in one cell, strong clutter with 4 stations in
seed 7331 ($+0.077$), and two further intervals have a lower bound of $+0.000$; at $\le0.1$, outside
the rule, the differences range from $+0.03$ to $+0.55$ and the intervals of both strong-clutter
4-station cells include zero. By the registered interpretation, the advantage of MT-FW is not
attributed to the moment tensor with deflation already present in [129], which BOMP has, but to the
elements proper to MT-FW, its peak search and its cap; the exploratory diagnosis below separates the
two, and Task 2D tests the attribution with equal tolerance (Y.21).

**Rule 2C-b: the convex method falls below its greedy approximation.** Solved to certified optimality,
GL-FW+SC is below MT-FW+SC in 23 of the 24 combinations of cell, seed and cap, by $0.06$ to $0.61$, and
equal to it by the registered criterion in one, strong clutter with 4 stations in seed 7331 at $\le0.1$
($-0.013$); it never exceeds it (Table Y.39). Every interval, over windows or over networks, excludes
zero except in that 4-station cell of seed 7331. Against BOMP+SC on the same dictionary, whose
selection statistic is its certificate, GL-FW+SC ranges from $-0.18$ to $+0.01$, at or slightly below
it.

**Rule 2C-c holds in both regimes, with and without Holm's correction.** Over thirty networks per
cell, MT-FW+SC exceeds QM-3C+SC at $\le0.25$ by $+0.22$ to $+0.35$ with weak clutter and by $+0.10$ to
$+0.30$ with strong clutter, and every one of the twelve hierarchical intervals excludes zero, the
lowest lower bound being $+0.029$, strong clutter with 4 stations in seed 7331; at $\le0.1$ the twelve
intervals also exclude zero (Table Y.40). Under the literal reading the rule passes in three of three
cells of each seed in both regimes. Under Holm's correction it passes in three of three cells of each
seed with weak clutter, and in three and two with strong clutter, where the 4-station cell of seed 7331
($+0.099$ [$+0.029$, $+0.181$], adjusted p-value $0.044$) is not rejected; the verdict is the same under
both readings. Re-applied in every hierarchical replicate (exploratory), the literal reading holds in
100 per cent of the replicates in both regimes. With weak clutter all 36 comparisons of the family are
rejected in favour of MT-FW+SC; with strong clutter 28 of 36 are (Table Y.41), and the eight that are
not all lie in the 4-station cells: in seed 2024 the comparison with BOMP+SC at both caps ($+0.065$
[$+0.003$, $+0.141$] and $+0.080$ [$+0.006$, $+0.165$]), and in seed 7331 the comparisons with QM-3C+SC,
BOMP+SC and MTSTACK+SC at both caps. In the sparsest networks with strong clutter MT-FW+SC is therefore
not shown to be better than the closest antecedents once multiplicity is accounted for. On the
twenty-five networks of part (ii) alone (exploratory), MT-FW+SC minus QM-3C+SC is $+0.17$ to $+0.32$ at
$\le0.1$ and $+0.20$ to $+0.39$ at $\le0.25$ with weak clutter, and $+0.07$ to $+0.23$ and $+0.08$ to
$+0.26$ with strong clutter, the same magnitudes as on the five networks of Task 2A.

**Deviations.** Nine departures from the registered text are declared. First, the registration asks
for the duality gap at every value of $\lambda$ of the path. The path is stopped at the first value
that completes a five-event record; the record is the same as that of the full path, because blocks
entering later do so at a strictly smaller $\lambda$ and sort after the five accepted ones, which the
24-window check of the validation confirms. The full path grows to hundreds of active blocks and took
more than a minute per value of $\lambda$ at $\lambda/\lambda_{\max}\approx0.5$ in one window, so the
gap is reported for the solved values (2 to 14 per window), not for all 30. Second, the registration
does not fix the weighting of the group norm; the block-whitened norm is the weighting of (16a) and
makes the certificate equal to the selection statistic of BOMP. Third, the score of BOMP is read as the
projection energy of the residual on the entering block, not as the drop of the residual after the
joint refit. Fourth, MTSTACK inverts with the plain pseudo-inverse of $G$, which ignores the rare
temporal overlap of P and S that the Gram of BOMP includes, and has no cap, as the registered
description implies. Fifth, the tolerances are $\lambda(1+10^{-6})$ for admission and $10^{-9}$ for the
restricted problem; a tolerance of $10^{-11}$ could not be reached in floating point and was set to
$10^{-9}$ before any evaluation run. Sixth, the kernel was pinned explicitly, as described above, and
scipy's kernel differs from Task 2A on the Prescott networks, without use by these methods. Seventh,
the registration names Holm's correction over the comparisons of each regime, and the family of 36
comparisons with the one-sided level 0.025 is the reading adopted in the analysis; bootstrap p-values
have a resolution of $1/2001$, and rule 2C-c is reported under the literal and the Holm reading, which
agree. Eighth, the point estimates of rule 2C-c pool windows, so each of the five Task 2A networks
weighs 20/12 relative to one of the added networks, whereas the hierarchical bootstrap weighs networks equally. Ninth,
for networks 0 to 4 rule 2C-c uses the part-(i) records of BOMP+SC and MTSTACK+SC on the same windows
rather than a rerun. The review of the code before the evaluation runs also added the validation script
that the runner expected, bounded the re-solves of the restricted problem and replaced a dense
whitening by a blockwise product; none of these changes was made after any result of seeds 2024 or
7331 existed.

**Why the exact-time methods fall so far (exploratory, post hoc).** The analyses of this paragraph and
of Table Y.42 were not registered; they were designed after the results of part (i) had been seen,
and no verdict depends on them. The exact-time comparators fall not only below MT-FW+SC but mostly below
QuakeMigrate (Table Y.37), which invites the suspicion of a broken implementation, so the cause was
checked on seed-999 streams, with 60 windows per cell and one network per window. With the ideal
generator of Y.11, in which the detectors' model is exact, with no statics and no model mismatch, but
with sources and origin times off the grid as everywhere in this appendix, BOMP, MTSTACK and GL-FW
already reach only 0.17 to 0.67 at $\le0.25$, against 0.71 to 0.89 for MT-FW; with the perturbations of
Y.11 all the exact-time methods fall further. With the same generator and every source moved to its
nearest node and its origin time rounded to the grid, so that only the rounding to the sample remains,
BOMP becomes the best of the four methods, 0.97 against 0.76 for MT-FW at $\le0.25$ in the 4-station
cell, and GL-FW rises to 0.73 there. On the Task 2A windows of seed 2024 with weak clutter, replacing
the estimated station corrections by oracle corrections, the true static plus the median model
mismatch, moves BOMP+SC at $\le0.25$ from 0.23, 0.22 and 0.38 to 0.26, 0.26 and 0.47 in the 4-station,
6-station and one-sided cells, and MT-FW+SC from 0.52, 0.56 and 0.65 to 0.56, 0.53 and 0.73.

**Table Y.42. Exploratory, post hoc. Recall at $\le0.1$ / $\le0.25$ with the ideal generator of Y.11,
60 windows per cell on seed-999 streams, one network per window: sources and origin times off the grid,
and every source moved to its nearest node with its origin time on the grid.** No station corrections
are used, since the generator has no statics. Outputs `task2c_explore_all.txt` and
`task2c_explore_e3.txt`.

| Sources | Cell | events | MT-FW | BOMP | MTSTACK | GL-FW |
|:--|:--|:--|:--|:--|:--|:--|
| off the grid (ideal) | 4 stations | 59 | 0.627/0.712 | 0.220/0.458 | 0.119/0.169 | 0.136/0.237 |
| off the grid (ideal) | 6 stations | 44 | 0.864/0.886 | 0.182/0.432 | 0.091/0.182 | 0.114/0.273 |
| off the grid (ideal) | one-sided, 8 st. | 58 | 0.810/0.810 | 0.500/0.672 | 0.483/0.603 | 0.414/0.500 |
| on the nodes (ideal) | 4 stations | 59 | 0.678/0.763 | 0.898/0.966 | 0.203/0.542 | 0.424/0.729 |
| on the nodes (ideal) | 6 stations | 44 | 0.909/0.932 | 0.955/0.977 | 0.659/0.818 | 0.750/0.864 |
| on the nodes (ideal) | one-sided, 8 st. | 58 | 0.793/0.845 | 0.914/0.948 | 0.672/0.759 | 0.741/0.914 |

Four readings follow, with the caution that they are post hoc. The implementations of BOMP and GL-FW
are sound: with sources on the nodes BOMP is the best of the four methods, and GL-FW comes close to
MT-FW at $\le0.25$, though not at $\le0.1$. Their collapse on the evaluation windows comes from reading
the matched filter at the exact predicted sample: with a node spacing of 3 km and origin steps of 0.1
s, a source off the nodes produces station-dependent moveout errors that are a sizeable fraction of the
main lobe of the 2 Hz autocorrelation, the residual errors of the station corrections (rms 0.12 to 0.28
s, Y.18) add to them, and the amplitudes read there are attenuated or change sign, which destroys the
fit of the moment tensor. Better corrections help BOMP little, so the dominant term is the
discretization of the grid rather than the error of the corrections. MT-FW's peak search within
$\pm0.25$ s absorbs both. And because BOMP, which has neither a peak search nor a cap, is the best
method once the timing is exact, the cap is not what the exact-time methods lack: within the
attribution of rule 2C-a, the element that matters on these data is the peak search. No ablation that
removes the peak search or the cap from MT-FW itself was run, so this attribution rests on the
comparison across methods.

**Ties in the score of GL-FW (exploratory).** The registered score of GL-FW, the value of $\lambda$ on
the path at which a block enters, is coarse, since adjacent values differ by 12 per cent: 1694 of 3000
recorded events with weak clutter and 1654 of 3000 with strong clutter share their score with another
event of the same window. A variant that breaks the ties by the group norm at entry changes recall by
at most 0.015, so the ties are not what keeps GL-FW below BOMP. The more likely cause is the shrinkage of
the penalized fit: a strong event is under-fitted, and the energy it leaves in the residual lets
spurious blocks more than 5 km away enter before weaker true events, whereas the orthogonal refit of
BOMP removes that energy. This is an interpretation, not a measurement.

**What Task 2C does not test.** The comparators are implemented as the registration fixed them, on the
shared grid of 3 km and 0.1 s and with exact predicted times; a block-sparse inversion as [129] would
run it on real data, with full waveforms, a finer grid or a tolerance to time shifts, was not run, so the
result does not show that MT-FW would beat a block-sparse inversion so equipped. The path of GL-FW was
solved until its record was complete, not to $\lambda_{\max}/30$. The operating point of every method is
chosen on the evaluation windows, as in Tasks 2A and 2B. The replication over thirty networks concerns
MT-FW+SC, QM-3C+SC, BOMP+SC and MTSTACK+SC; the pick-based associators, rule 2A-e and the ablation of
Task 2B were not repeated over thirty networks. The data are the synthetic generator of Y.14, with
double-couple sources inside the moment-tensor subspace of the detector, and MALMI [110], the stacking
of [112], real noise and real data remain untested. An independent code audit found later that "exact predicted times" means the
floor of the predicted sample, not the nearest one (error E52): with no timing error a reader at the floor sample keeps
on average 0.78 of the energy of a 2 Hz arrival (worst case 0.42), against 0.94 at the nearest sample. The exact-time
comparators of this task, the convex problem of rule 2C-b and the MT-FW-EX ablation of Task 2D are handicapped by that
detail, which inflates the measured value of the peak search; the decisions very probably stand, since the timing
error of the 3 km grid is already 0.1 to 0.3 s, but the cells were not re-run with rounding.

### Y.20a What Task 2C changes

Task 2C changes three statements of this proposal and leaves the rest of Appendix Y as it stands.

The comparison with the closest antecedents exists, and what it supports is narrower than superiority
over sparse moment-tensor inversion. On these synthetics, at equal realized false-event rate, MT-FW with
station corrections detects more events than block orthogonal matching pursuit in the form of [129] and
mechanism-corrected stacking in the form of [122], in every cell of both seeds and both clutter regimes,
and by the registered interpretation the advantage belongs to the elements proper to MT-FW rather than
to the moment tensor with deflation that [129] already has. The exploratory diagnosis locates it in the
peak search: on a grid of 3 km and 0.1 s, methods that read the matched filter at the exact predicted
time lose most of their recall to discretization and residual timing errors, which the peak search
absorbs, and with sources on the nodes block matching pursuit is better than MT-FW (0.97 against 0.76
at $\le0.25$ in the 4-station cell). The supportable claim is therefore that MT-FW is robust to grid and
timing errors on a coarse grid, and that on these synthetics it beats pick-based associators,
QuakeMigrate, block orthogonal matching pursuit [129] and mechanism-corrected stacking [122]. It is not
that MT-FW is better than sparse moment-tensor inversion as such. Task 2D (Y.21) gave [129] the same
peak search and found that it then matches MT-FW, so the advantage belongs to the peak search alone
and not to other elements of MT-FW (Y.21a).

The convex problem is the framework, the certificate and the theory, and on a coarse grid with timing
errors it is not the best detector. Solved to certified optimality on the grid, the group BLASSO (16a)
stays $0.06$ to $0.61$ below MT-FW in 23 of 24 cases and at or slightly below BOMP. The performance of
MT-FW comes from ingredients that the grid problem does not contain: above all the peak search, and
possibly the unpenalized fit of each new atom, which avoids the shrinkage of the penalized problem.
What (16a) contributes is the framework in which the detection statistic is an optimality certificate,
the generalized likelihood ratio at zero that the detector and its calibrated threshold use
(Proposition AA.2), and the object of the theory of Appendix AA. What would make the convex method
competitive is to remove the discretization that the peak search absorbs. The sliding Frank-Wolfe
algorithm [116] moves positions and origin times continuously; per-station time shifts entered as
nuisance parameters would address the residual timing error of the corrections, which a continuous
source position does not; and a least-squares refit on the support selected by the penalized problem,
a standard debiasing, would address the shrinkage. Whether the off-the-grid solver closes the gap
measured here, matches MT-FW or exceeds it on the windows of Task 2A is the key research question of the
mathematical track (Section 19, T2 and E2). The preregistered Task 2E answered it: at equal tolerance
the off-grid solver does not exceed block pursuit (Y.22).

The network is the unit of the comparison with QuakeMigrate. With thirty networks per cell and a
hierarchical bootstrap, the advantage of MT-FW+SC over QM-3C+SC has intervals that exclude zero in all
twelve cells at both caps, and it survives Holm's correction in eleven of twelve at each cap; the
weakest cell is strong clutter with four stations, where MT-FW+SC is also not shown to be better than
BOMP+SC or MTSTACK+SC after correction. For that comparison this replaces the per-window intervals of
Y.18 and its crude five-network bootstrap; it does not replicate over thirty networks the other
comparators of rule 2A-a, rule 2A-e or the ablation of Task 2B.

What Y.18a and Y.19 list as not established remains so: the operating point is chosen on the
evaluation windows, the frozen-threshold rule fails under its strict reading, GaMMA with optimistic
picks is ahead at 0.1 with weak clutter, the sources are double couples inside the detector's subspace,
the one synthetic generator lies outside the declared regime of Section 1.5, and there is no result on
real data.

## Y.21 Equal timing tolerance, the peak-search ablation, and a frozen conformal threshold (Task 2D)

Code: `task2d_methods.py` (the five methods), `task2d_run.py` (validation, smoke test, runs and
controls), `task2d_validate.py` (validation on seed 999), `analyze_2d.py` (rules, tables, hierarchical
bootstrap and Holm's correction), `task2d_explore.py` and `task2d_explore_sum.py` (exploratory
diagnostics); outputs `task2d_validation.txt`, `task2d_smoke.txt`, `task2d_code_hashes.txt`,
`analyze_2d_out.txt`, `task2d_explore_e1.txt` and `task2d_explore_e2.txt`; design and rules in
`PREREGISTRO_1j_1k.txt` (section "Task 2D") with the hashes of the registry in
`PREREGISTRO_hashes.txt`; review of the code, controls, compute and deviations in `task2d_report.txt`.

**Why the task was run.** A second independent audit, carried out after Task 2C, named three weaknesses
of Y.20. The exact-time antecedents BOMP+SC, MTSTACK+SC and GL-FW+SC did not have the tolerance of
$\pm0.25$ s to timing error that MT-FW has, so rule 2C-a compared a tolerant method with intolerant
ones; the attribution of the margins of that rule to the peak search rested on a comparison across
methods, because no ablation had removed the peak search or the cap from MT-FW itself (Y.20); and the
primary reading of every task so far chose the operating point on the evaluation windows. Task 2D
answers the three in two parts. Part (i) gives the three antecedents the peak search of MT-FW and
ablates the peak search and the cap of MT-FW on the 1200 windows of Task 2A. Part (ii) freezes, for
MT-FW+SC and QM-3C+SC, a split-conformal threshold calibrated on noise-only windows of other networks
and applies it unchanged to those windows.

**Preregistration.** The design, three rules and their interpretation were written on 26 September
2026 at 21:13 UTC, together with the erratum on the time of Task 2B (Y.19), and the SHA-256 hash of the
registry file at that moment,
6e8510b3367f1ad55874c0c9afa3c27dcfef39affa4ab6ac177f64bf6064bf88, was recorded at 21:13:04 UTC in
`PREREGISTRO_hashes.txt`. The code was written between 02:25 and 02:29 UTC on 27 September, validated
between 02:29 and 02:36 and smoke-tested between 02:37 and 02:38, and the hashes of the three scripts,
of `analyze_2d.py` and of the launcher were recorded at 02:41:59 UTC in `task2d_code_hashes.txt`,
before any result of the task existed. Rule 2D-a, per clutter regime: MT-FW+SC minus the better of
BOMP-PK+SC and MTSTACK-PK+SC is at least $0.05$ at $\le0.25$ realized false events per window in at
least two of three cells in both seeds; if it fails, the advantage of MT-FW over [129] and [122] is
attributed to the peak search and to no other element of the method, and the proposal says so. Rule
2D-b, per regime: MT-FW+SC minus MT-FW-EX+SC is at least $0.05$ at $\le0.25$ in at least two of three
cells in both seeds, which shows that the peak search matters; MT-FW-NOCAP+SC minus MT-FW+SC and
GL-FW-PK+SC minus MT-FW+SC are reported with intervals and no threshold. Rule 2D-c, per regime: with
the frozen thresholds, the realized false rate of MT-FW+SC on the empty evaluation windows is at most
$2\alpha$ and its recall exceeds that of QM-3C+SC with its own frozen threshold by at least $0.05$, in
at least two of three cells in both seeds, at $\alpha=0.1$ and at $\alpha=0.25$. The statistics are the
hierarchical bootstrap over networks and windows of Task 2C, with Holm's correction [142] within each
regime. A budget rule would have reduced part (ii) to 100 empty windows per cell and regime had the
seed-999 estimate exceeded 6 hours on two cores. A note appended at 02:41:50 UTC on 27 September,
with the time read from the system clock, before any computation on seeds 2024 or 7331, records that
the smoke test estimated 0.64 hours on two cores (part (i) 0.92 hours of CPU time, part (ii) with 200
empty windows per cell and regime 0.36 hours), so the rule did not trigger and part (ii) used 200 empty
windows; the hash of the file with the note,
dc8cf70874488970050227c8756bb89e519c75d2b068311c9ae920a8c1c4e32b, was recorded at the same second,
and nothing else in the registration was changed. The result paragraph was appended at 03:54:55 UTC
with a third hash. The registration was written after Tasks 2A to 2C had been run on the same windows,
and its design was informed by them, in particular by the exploratory diagnosis of Y.20 that pointed to
the peak search; what it fixed in advance are the methods, the rules and their interpretation.

**Design.** Part (i) regenerates the 1200 windows of Task 2A (two seeds, two regimes, three cells,
networks 0 to 4, twenty windows each) bit for bit, exactly as Task 2C did, with the same station
corrections, the grid of 675 nodes and 31 origin times, the 5 km epicentral exclusion, at most five
events per window and the hit rule of Task 2A, and with the kernel of the numerical library pinned as
in Task 2C. MT-FW+SC and QM-3C+SC are the stored Task 2A records, and the exact-time BOMP+SC, MTSTACK+SC
and GL-FW+SC the stored Task 2C records, which serve as context. Part (ii) calibrates on the five
seed-999 calibration networks of Task 2A per cell, rebuilt with the same stations, statics and station
corrections, with forty empty windows each, 200 per cell and regime, drawn from the generator of Y.14
with no source, the clutter count and amplitude of the regime and noise, on random streams not used
before. The score of a window is its maximum: the first score of the MT-FW+SC record, and for QM-3C+SC
the largest score among the events that QuakeMigrate triggers with the frozen configuration of Task
2A, minus infinity when it triggers none. The threshold is the $\lceil(n+1)(1-\alpha)\rceil$-th
smallest of the $n=200$ maxima, the 181st for $\alpha=0.1$ and the 151st for $\alpha=0.25$, frozen per
cell, regime, method and $\alpha$, and applied unchanged to the 100 evaluation windows of each cell and
seed, of which 36 to 45 contain no source; an event is declared when its score exceeds the threshold.
For QuakeMigrate in the weak-clutter one-sided cell only 28 of the 200 empty calibration windows
triggered any event, so its threshold at $\alpha=0.25$ is minus infinity and every event it triggers is
declared.

**The methods as implemented.** The peak reading, "PK", is that of MT-FW: for a residual and the
predicted sample of a station and phase, which includes the station correction, the sample read is the
one of largest three-component energy within four samples, 0.25 s, of the predicted one, the first on
ties. It is tabulated once per residual. An atom is a node, an origin time and the read samples of
every station and phase; its columns predict the radiation coefficients times the matched-filter
autocorrelation centred at the read samples, and the Gram matrix between atoms is that of Y.20 with
the read samples in place of the predicted ones, including the cross terms of P and S at the lag of the
peaks; with the read samples equal to the predicted ones every quantity reduces to that of Y.20.
BOMP-PK+SC reads every block that is not excluded at the peaks of the current residual, selects the
block with the largest projection energy $e=\mathbf h^{\mathsf T}K(I)^{+}\mathbf h$ at those read
samples $I$, lets it enter with its read samples fixed, scores it $e/2$ at entry with a running
minimum, and refits all active atoms jointly by least squares at their fixed samples. MTSTACK-PK+SC
reads the matched-filter output at its own peaks, inverts the moment tensor by least squares, takes
half the energy of the prediction, maximizes over origin time and forms the catalogue by suppression
at 20 km with no deflation; its field is that of WAVE-MT without the cap. GL-FW-PK+SC is the group
BLASSO of Y.20 on a growing dictionary: the Frank-Wolfe oracle evaluates the certificate
$\sqrt{\mathbf h^{\mathsf T}K(I)^{+}\mathbf h}$ of every block at the peaks of the current residual, an
entering atom keeps those read samples, and the restricted problem is solved exactly over the fixed
atoms, along the path of Y.20 with its early stop and its entry-$\lambda$ score; $\lambda_{\max}$ is the
largest certificate on the matched-filter output, so that the first score of BOMP-PK is
$\lambda_{\max}^2/2$, and the duality gap is reported at every solved $\lambda$ with respect to the
realized dictionary, the active atoms at their fixed samples and the peak atoms of the final residual.
The three PK variants carry no cap, as the Task 2C methods did not. MT-FW-EX+SC is MT-FW with the peak
search disabled: every station and phase is read at the exact predicted sample, with the cap, and
deflation is centred there. MT-FW-NOCAP+SC is MT-FW with the peak search and without the cap. Both run
the dense loop of MT-FW rather than the certified branch and bound, which would need adapting; with the
peak search and the cap the same loop reproduces the Task 2A records exactly.

**Readings chosen before the results.** The registration leaves several choices open. The ones adopted
were fixed in the docstring of `analyze_2d.py`, whose hash was recorded before any result existed.
Rules 2D-a and 2D-b are decided on point estimates, the convention of Tasks 2A to 2C, with the
threshold $0.05-10^{-9}$. For rule 2D-c an event is declared at a score strictly above the threshold;
the realized false rate on empty windows is read as false events per empty window, which is at least
the fraction of empty windows with a false event, the quantity that the conformal guarantee bounds and
that is also reported; the primary reading requires the rate and the recall difference to hold in the
same cell, in at least two of three cells, and a strict reading, with the rate at most $2\alpha$ in
every cell, seed and $\alpha$ and the difference in at least two of three cells, is reported beside it;
and the thresholds are held fixed in the bootstrap, so the intervals do not include the uncertainty of
the calibration. The PK variants have no cap, and the cross term of P and S in the block Gram matrix is
evaluated at the lag of the peaks, which is the reading adopted of "identical to 2C except that the
sample of each station and phase is taken at the peak". In GL-FW-PK two atoms with the same node and
origin time and different read samples may coexist, the record deduplicates them by the 5 km rule, the
duality gap refers to the realized dictionary, and the steps of Frank-Wolfe per $\lambda$ are capped at
400, a cap never reached. The empty windows use a stream tag not used before, forty per calibration
network, and the thresholds pool the five networks of a cell. Holm's family has 36 comparisons per
regime, rules 2D-a (through the difference with the better comparator) and 2D-b at both caps and rule
2D-c at both values of $\alpha$, times three cells and two seeds, at the one-sided level 0.025 with
bootstrap p-values $(1+\#\{d^\ast\le0\})/2001$, so the first Holm step can reach at best
$36/2001=0.018$; the descriptive quantities are not in the family.

**Validation.** Before any evaluation seed was run, `task2d_validate.py` checked the implementation on
windows 0 and 1 of network 0 of seed 999 in every regime and cell (`task2d_validation.txt`), and checks
V1 to V5 passed. The tabulated peak samples equal those of the sampling of MT-FW in all 669,600 to
1,339,200 readings per cell, on the matched-filter output and on a residual deflated by one step of
MT-FW+SC. With the peak reading disabled, BOMP-PK, MTSTACK-PK and GL-FW-PK return exactly the Task 2C
records of BOMP+SC, MTSTACK+SC and GL-FW+SC in 12 of 12 windows each, and the dense loop with peak
search and cap equals the certified MT-FW+SC in 12 of 12. Explicit data-domain atoms at read samples
shifted at random by up to four samples, overlapping pairs included, reproduce the correlations to a
relative $1.2\times10^{-15}$ and the Gram matrix to $6.4\times10^{-16}$, the Gram matrix of a single
atom to $3\times10^{-17}$, and the pseudo-inverse of the projection energy of a subset to
$9\times10^{-13}$. Independent code paths agree: the first score of MT-FW-NOCAP equals the maximum of the
MTSTACK-PK field to a relative $3.3\times10^{-15}$, the first score of MT-FW-EX a vectorized maximum at
the exact capped samples to $4\times10^{-16}$, and the first score of BOMP-PK equals
$\lambda_{\max}^2/2$ of GL-FW-PK to $2.2\times10^{-16}$; the refit of BOMP-PK is orthogonal to
$7.9\times10^{-15}$. Solving two values of $\lambda$ beyond the early stop, the largest certificate of
GL-FW-PK off the support is $1.000000000\,\lambda$, the certificate on the support equals $\lambda$ to
$9.6\times10^{-10}$, the relative duality gap is at most $3.6\times10^{-9}$, and the early-stopped
record equals the record of the longer path in 12 of 12 windows, with no step cap reached and no
re-solve. Across the 1200 windows of part (i), GL-FW-PK stopped with a complete record in every window
after 2 to 15 solved values of $\lambda$ (median 3), with a median of 11 to 13 Frank-Wolfe steps per
window and at most 51, a largest relative gap of $1.7\times10^{-8}$, a largest ratio of oracle
certificate to $\lambda$ of $1.000000001$, a largest deviation on the support of $9.6\times10^{-10}$,
no re-solve and a final support of median 11 to 12 atoms and at most 37. In BOMP-PK, 91 to 92 per cent
of the station and phase readings of an entering atom lie off the predicted sample, by 2.41 to 2.44
samples on average, and the refit is orthogonal to $9.2\times10^{-14}$.

**Controls and kernels.** For all 60 networks of part (i), the station positions, statics, catalogue
draws, catalogue locations and origin times and station corrections are identical bit for bit to those
of Task 2A; truth, origin times and the number of sources agree in 1200 of 1200 windows, asserted again
against Tasks 2A and 2C in the analysis; and MT-FW+SC recomputed on windows 0 and 1 of every network
equals the stored Task 2A records exactly in 120 of 120 windows. The default kernel of numpy's OpenBLAS
on this machine had changed again, to Cooperlake, which is neither Prescott nor SkylakeX, so the kernel
was set explicitly per process: SkylakeX for the 18 networks listed in Y.19 and Prescott for the other
42, each checked with threadpoolctl and recorded in its part file. For the 30 calibration networks of
part (ii), run with Prescott, stations and statics equal those of the seed-999 parts of Task 2A, the
station corrections are identical (largest difference 0), and the pipeline reproduces the stored Task
2A records of MT-FW+SC and QM-3C+SC on windows 0 and 1 of every network, 60 of 60 each. No control
failed.

**Statistics and compute.** The endpoint is that of Task 2A: event recall at $\le0.1$ and $\le0.25$
realized false events per window on the exact FROC curve, with the operating point chosen on the
evaluation windows in part (i). The intervals come from a hierarchical bootstrap with 2000 replicates,
drawing the five networks of a cell with replacement and then the windows of each drawn network, with
the same resamples for every method of a cell; in part (i) the operating point and the better
comparator are re-chosen in every replicate, and in part (ii) the thresholds stay fixed. The bootstrap
over windows of Task 2A is kept in the output as a secondary reading. With five networks per cell the
hierarchical intervals are coarse. Part (i) took 3931 s of CPU time, with 0.47 s per window for
BOMP-PK, 0.013 s for MTSTACK-PK, 1.76 s for GL-FW-PK (99th percentile 4.8 s, at most 8.7 s), 0.23 s for
MT-FW-EX and 0.76 s for MT-FW-NOCAP; part (ii) took 1182 s, with 0.20 s per empty window for MT-FW+SC
and 0.72 s for QM-3C+SC; 1.42 hours of CPU time in all, between 02:42 and 03:25 UTC on two workers.

**Table Y.43. Part (i). Event recall at realized false events per window $\le0.1$ / $\le0.25$ on the
1200 windows of Task 2A (100 per cell).** MT-FW+SC and QM-3C+SC are the Task 2A records (Tables Y.29,
Y.30 and Y.37), and the exact-time antecedents of Task 2C are in Table Y.37. Script `analyze_2d.py`;
output `analyze_2d_out.txt`, Table 1 of each regime.

| Regime | Cell (seed) | events | MT-FW+SC | BOMP-PK+SC | MTSTACK-PK+SC | GL-FW-PK+SC | MT-FW-EX+SC | MT-FW-NOCAP+SC | QM-3C+SC |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| weak | 4 stations (2024) | 87 | 0.322/0.517 | 0.241/0.437 | 0.126/0.218 | 0.092/0.241 | 0.172/0.287 | 0.241/0.437 | 0.184/0.253 |
| weak | 6 stations (2024) | 72 | 0.458/0.556 | 0.361/0.514 | 0.167/0.319 | 0.264/0.389 | 0.250/0.403 | 0.347/0.514 | 0.278/0.306 |
| weak | one-sided, 8 st. (2024) | 71 | 0.549/0.648 | 0.408/0.592 | 0.225/0.423 | 0.225/0.493 | 0.183/0.338 | 0.408/0.577 | 0.352/0.380 |
| weak | 4 stations (7331) | 78 | 0.321/0.513 | 0.205/0.423 | 0.179/0.256 | 0.167/0.269 | 0.167/0.295 | 0.205/0.423 | 0.154/0.167 |
| weak | 6 stations (7331) | 68 | 0.721/0.794 | 0.706/0.824 | 0.368/0.588 | 0.368/0.647 | 0.265/0.382 | 0.706/0.838 | 0.235/0.338 |
| weak | one-sided, 8 st. (7331) | 81 | 0.642/0.765 | 0.568/0.741 | 0.235/0.420 | 0.259/0.469 | 0.173/0.358 | 0.568/0.741 | 0.469/0.481 |
| strong | 4 stations (2024) | 77 | 0.312/0.390 | 0.234/0.299 | 0.091/0.234 | 0.078/0.234 | 0.169/0.234 | 0.234/0.299 | 0.117/0.130 |
| strong | 6 stations (2024) | 74 | 0.378/0.486 | 0.446/0.486 | 0.203/0.351 | 0.135/0.284 | 0.243/0.338 | 0.446/0.486 | 0.122/0.162 |
| strong | one-sided, 8 st. (2024) | 83 | 0.313/0.422 | 0.301/0.458 | 0.169/0.313 | 0.145/0.313 | 0.120/0.181 | 0.301/0.458 | 0.157/0.193 |
| strong | 4 stations (7331) | 78 | 0.128/0.218 | 0.167/0.269 | 0.090/0.205 | 0.103/0.231 | 0.141/0.154 | 0.167/0.269 | 0.038/0.051 |
| strong | 6 stations (7331) | 84 | 0.357/0.417 | 0.298/0.571 | 0.179/0.274 | 0.167/0.298 | 0.107/0.250 | 0.286/0.548 | 0.179/0.250 |
| strong | one-sided, 8 st. (7331) | 75 | 0.627/0.707 | 0.653/0.760 | 0.427/0.547 | 0.347/0.573 | 0.200/0.413 | 0.653/0.760 | 0.253/0.293 |

**Table Y.44. Rule 2D-a. MT-FW+SC minus the better of BOMP-PK+SC and MTSTACK-PK+SC, point [95%
hierarchical interval], and MT-FW+SC minus MTSTACK-PK+SC at $\le0.25$.** The better comparator is
BOMP-PK+SC in every cell at both caps. The rule reads the cap $\le0.25$. Same output, Table 2.

| Regime | Cell (seed) | $\le0.1$: point [hierarchical] | $\le0.25$: point [hierarchical] | minus MTSTACK-PK+SC, $\le0.25$ | cell passes |
|:--|:--|:--|:--|:--|:--|
| weak | 4 stations (2024) | +0.080 [-0.202, +0.257] | +0.080 [-0.136, +0.282] | +0.299 [+0.118, +0.464] | yes |
| weak | 6 stations (2024) | +0.097 [-0.118, +0.271] | +0.042 [-0.092, +0.224] | +0.236 [+0.048, +0.407] | no |
| weak | one-sided, 8 st. (2024) | +0.141 [-0.036, +0.469] | +0.056 [+0.000, +0.265] | +0.225 [+0.034, +0.483] | yes |
| weak | 4 stations (7331) | +0.115 [-0.051, +0.316] | +0.090 [-0.029, +0.256] | +0.256 [+0.050, +0.375] | yes |
| weak | 6 stations (7331) | +0.015 [-0.148, +0.175] | -0.029 [-0.125, +0.047] | +0.206 [+0.059, +0.343] | no |
| weak | one-sided, 8 st. (7331) | +0.074 [-0.119, +0.270] | +0.025 [-0.084, +0.155] | +0.346 [+0.178, +0.521] | no |
| strong | 4 stations (2024) | +0.078 [-0.086, +0.238] | +0.091 [-0.048, +0.257] | +0.156 [+0.043, +0.392] | yes |
| strong | 6 stations (2024) | -0.068 [-0.289, +0.203] | +0.000 [-0.138, +0.145] | +0.135 [+0.000, +0.291] | no |
| strong | one-sided, 8 st. (2024) | +0.012 [-0.152, +0.164] | -0.036 [-0.130, +0.108] | +0.108 [-0.042, +0.317] | no |
| strong | 4 stations (7331) | -0.038 [-0.149, +0.154] | -0.051 [-0.206, +0.095] | +0.013 [-0.157, +0.197] | no |
| strong | 6 stations (7331) | +0.060 [-0.250, +0.233] | -0.155 [-0.277, +0.032] | +0.143 [-0.086, +0.323] | no |
| strong | one-sided, 8 st. (7331) | -0.027 [-0.154, +0.099] | -0.053 [-0.158, +0.030] | +0.160 [+0.000, +0.386] | no |

**Table Y.45. Rule 2D-b, and the descriptive readings of the cap and of the convex method. MT-FW+SC
minus MT-FW-EX+SC at both caps, MT-FW-NOCAP+SC minus MT-FW+SC and GL-FW-PK+SC minus MT-FW+SC at
$\le0.25$, point [95% hierarchical interval].** Same output, Tables 3 and 4.

| Regime | Cell (seed) | minus MT-FW-EX+SC, $\le0.1$ | minus MT-FW-EX+SC, $\le0.25$ | cell passes | MT-FW-NOCAP+SC minus MT-FW+SC, $\le0.25$ | GL-FW-PK+SC minus MT-FW+SC, $\le0.25$ |
|:--|:--|:--|:--|:--|:--|:--|
| weak | 4 stations (2024) | +0.149 [-0.014, +0.350] | +0.230 [+0.028, +0.388] | yes | -0.080 [-0.286, +0.133] | -0.276 [-0.510, -0.027] |
| weak | 6 stations (2024) | +0.208 [-0.127, +0.362] | +0.153 [-0.118, +0.348] | yes | -0.042 [-0.226, +0.092] | -0.167 [-0.397, +0.013] |
| weak | one-sided, 8 st. (2024) | +0.366 [+0.164, +0.534] | +0.310 [+0.154, +0.515] | yes | -0.070 [-0.250, +0.000] | -0.155 [-0.488, -0.012] |
| weak | 4 stations (7331) | +0.154 [-0.053, +0.381] | +0.218 [-0.038, +0.378] | yes | -0.090 [-0.256, +0.029] | -0.244 [-0.392, -0.034] |
| weak | 6 stations (7331) | +0.456 [+0.194, +0.640] | +0.412 [+0.108, +0.597] | yes | +0.044 [-0.042, +0.129] | -0.147 [-0.347, -0.014] |
| weak | one-sided, 8 st. (7331) | +0.469 [+0.110, +0.713] | +0.407 [+0.212, +0.619] | yes | -0.025 [-0.168, +0.082] | -0.296 [-0.443, -0.032] |
| strong | 4 stations (2024) | +0.143 [-0.080, +0.250] | +0.156 [+0.041, +0.355] | yes | -0.091 [-0.253, +0.048] | -0.156 [-0.354, -0.045] |
| strong | 6 stations (2024) | +0.135 [-0.068, +0.329] | +0.149 [+0.000, +0.327] | yes | +0.000 [-0.143, +0.141] | -0.203 [-0.403, +0.000] |
| strong | one-sided, 8 st. (2024) | +0.193 [+0.037, +0.371] | +0.241 [+0.025, +0.433] | yes | +0.036 [-0.121, +0.125] | -0.108 [-0.359, +0.029] |
| strong | 4 stations (7331) | -0.013 [-0.118, +0.191] | +0.064 [-0.096, +0.240] | yes | +0.051 [-0.096, +0.197] | +0.013 [-0.141, +0.136] |
| strong | 6 stations (7331) | +0.250 [+0.009, +0.430] | +0.167 [-0.024, +0.410] | yes | +0.131 [-0.092, +0.253] | -0.119 [-0.320, +0.107] |
| strong | one-sided, 8 st. (7331) | +0.427 [+0.225, +0.605] | +0.293 [+0.098, +0.472] | yes | +0.053 [-0.027, +0.159] | -0.133 [-0.321, +0.016] |

**Table Y.46. Descriptive. The effect of the peak reading, and the tolerant antecedents against
QuakeMigrate, point differences at $\le0.1$ / $\le0.25$.** BOMP, MTSTACK and GL-FW are the exact-time
records of Task 2C. Same output, Table 5.

| Regime | Cell (seed) | BOMP-PK minus BOMP | MTSTACK-PK minus MTSTACK | GL-FW-PK minus GL-FW | GL-FW-PK minus BOMP-PK | better of BOMP-PK and MTSTACK-PK minus QM-3C+SC |
|:--|:--|:--|:--|:--|:--|:--|
| weak | 4 stations (2024) | +0.080 / +0.207 | +0.034 / +0.057 | +0.034 / +0.034 | -0.149 / -0.195 | +0.057 / +0.184 |
| weak | 6 stations (2024) | +0.278 / +0.292 | +0.097 / +0.181 | +0.181 / +0.264 | -0.097 / -0.125 | +0.083 / +0.208 |
| weak | one-sided, 8 st. (2024) | +0.211 / +0.211 | +0.155 / +0.197 | +0.155 / +0.296 | -0.183 / -0.099 | +0.056 / +0.211 |
| weak | 4 stations (7331) | +0.051 / +0.141 | +0.103 / +0.103 | +0.090 / +0.051 | -0.038 / -0.154 | +0.051 / +0.256 |
| weak | 6 stations (7331) | +0.529 / +0.471 | +0.250 / +0.324 | +0.250 / +0.412 | -0.338 / -0.176 | +0.471 / +0.485 |
| weak | one-sided, 8 st. (7331) | +0.407 / +0.494 | +0.111 / +0.247 | +0.173 / +0.309 | -0.309 / -0.272 | +0.099 / +0.259 |
| strong | 4 stations (2024) | +0.039 / +0.039 | +0.013 / +0.052 | +0.039 / +0.052 | -0.156 / -0.065 | +0.117 / +0.169 |
| strong | 6 stations (2024) | +0.351 / +0.324 | +0.122 / +0.216 | +0.041 / +0.149 | -0.311 / -0.203 | +0.324 / +0.324 |
| strong | one-sided, 8 st. (2024) | +0.193 / +0.241 | +0.072 / +0.157 | +0.060 / +0.157 | -0.157 / -0.145 | +0.145 / +0.265 |
| strong | 4 stations (7331) | +0.064 / +0.128 | -0.013 / +0.090 | -0.013 / +0.077 | -0.064 / -0.038 | +0.128 / +0.218 |
| strong | 6 stations (7331) | +0.226 / +0.357 | +0.143 / +0.155 | +0.143 / +0.190 | -0.131 / -0.274 | +0.119 / +0.321 |
| strong | one-sided, 8 st. (7331) | +0.587 / +0.507 | +0.347 / +0.373 | +0.267 / +0.440 | -0.307 / -0.187 | +0.400 / +0.467 |

**Rule 2D-a fails in both regimes.** At $\le0.25$ MT-FW+SC exceeds the better tolerant antecedent by at
least $0.05$ in two of three cells in seed 2024 and in one in seed 7331 with weak clutter, and in one
and in none with strong clutter (Table Y.44). The differences run from $-0.155$ to $+0.091$, every
hierarchical interval contains zero (one with a lower bound of $+0.000$), and the better antecedent is
BOMP-PK+SC in every cell; with strong clutter BOMP-PK+SC is equal to or better than MT-FW+SC in five of
six cells. No difference lies within 0.005 of the threshold. Re-applied in every hierarchical replicate
(exploratory), the rule holds in 11 per cent of the replicates with weak clutter and 0.1 per cent with
strong clutter (0.001 in `analyze_2d_out.json`, printed as 0 per cent in the rounded text output), and under Holm's correction no cell passes in any seed or regime. At $\le0.1$, outside
the rule, the differences run from $-0.068$ to $+0.141$ and every interval also contains zero. The
margins of rule 2C-a, $+0.08$ to $+0.52$, were therefore almost entirely the missing tolerance. By the
registered interpretation, the advantage of MT-FW over [129] and [122] is attributed to the peak search
and to no other element of the method. Mechanism-corrected stacking with the peak search stays below
MT-FW+SC at $\le0.25$, by $+0.21$ to $+0.35$ with weak clutter and $+0.01$ to $+0.16$ with strong
clutter: what separates it from MT-FW and from BOMP-PK is deflation, not the tolerance.

**Rule 2D-b holds in both regimes.** Removing the peak search from MT-FW costs $0.064$ to $0.412$ in
recall at $\le0.25$, a positive difference in all twelve combinations of cell and seed (Table Y.45), so
the rule passes in three of three cells of each seed in both regimes; re-applied in the hierarchical
replicates (exploratory) it holds in 99 per cent of them with weak clutter and 95 per cent with strong
clutter. The sign is consistent, but the network-level support of single cells is weak: the lower bound
of the hierarchical interval is above zero in 7 of 12 cases and Holm's correction rejects in 2 of 12,
so under the Holm reading the rule fails in both regimes (one and none of three cells with weak
clutter, none and one with strong clutter).

**The cap and the convex method (descriptive).** Without its cap MT-FW changes by $-0.091$ to $+0.131$
at $\le0.25$, and every interval contains zero; the cap helps a little with weak clutter, where five of
six differences are negative, and not with strong clutter. MT-FW-NOCAP+SC and BOMP-PK+SC have nearly
identical recall in every cell (Table Y.43): without the cap, MT-FW is essentially block orthogonal
matching pursuit with the peak reading, the two differing only by the P and S terms of the Gram matrix
and by deflation in place of an orthogonal refit. The cap is therefore the only element of MT-FW
beyond [129] with the tolerance, and its effect is small and not significant. The group BLASSO with the
same tolerance, GL-FW-PK+SC, improves on the exact-time GL-FW+SC by $+0.03$ to $+0.44$ at $\le0.25$,
as the peak reading improves BOMP by $+0.04$ to $+0.51$ and MTSTACK by $+0.05$ to $+0.37$, positive in
12 of 12 cases for each method (Table Y.46); but it stays $0.11$ to $0.30$ below MT-FW+SC in 11 of 12
cases, with the twelfth, strong clutter with 4 stations in seed 7331, at $+0.013$, and $0.04$ to $0.27$
below BOMP-PK+SC in all twelve. The conclusion of Y.20 therefore stands with equal tolerance: the convex
problem on the grid, solved to certified optimality for the dictionary it builds, is not what delivers
the performance. That certificate covers the realized dictionary only, because the per-station peak
oracle is a heuristic and not a global maximizer over time shifts. With the tolerance, block orthogonal
matching pursuit in the form of [129], the better of the two tolerant antecedents in every cell, exceeds
QuakeMigrate by $+0.17$ to $+0.49$ at $\le0.25$ in 12 of 12 cases: the advantage over the published
comparators established in Tasks 2A and 2C belongs to peak-tolerant moment-tensor detection with
deflation as a family, not specifically to MT-FW.

**Table Y.47. Rule 2D-c. Split-conformal thresholds frozen on 200 empty seed-999 windows per cell and
applied to the Task 2A evaluation windows: recall of MT-FW+SC and QM-3C+SC, their difference and the
false events of MT-FW+SC per empty window, point [95% hierarchical interval with the thresholds
fixed], and false events per window over all windows.** "Cell passes" is the per-cell reading. Same
output, Table 7 of each regime.

| Regime | Cell (seed) | $\alpha$ | empty windows | recall MT-FW+SC / QM-3C+SC | difference [hierarchical] | MT-FW+SC false events per empty window [hierarchical] | false events per window, all windows, MT-FW+SC / QM-3C+SC | cell passes |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| weak | 4 stations (2024) | 0.1 | 39 | 0.667 / 0.230 | +0.437 [+0.275, +0.565] | 0.179 [0.034, 0.351] | 0.840 / 0.170 | yes |
| weak | 4 stations (2024) | 0.25 | 39 | 0.690 / 0.241 | +0.448 [+0.275, +0.592] | 0.410 [0.189, 0.711] | 1.090 / 0.200 | yes |
| weak | 6 stations (2024) | 0.1 | 44 | 0.639 / 0.292 | +0.347 [+0.121, +0.532] | 0.045 [0.000, 0.122] | 0.730 / 0.150 | yes |
| weak | 6 stations (2024) | 0.25 | 44 | 0.681 / 0.306 | +0.375 [+0.190, +0.534] | 0.159 [0.000, 0.400] | 0.980 / 0.230 | yes |
| weak | one-sided, 8 st. (2024) | 0.1 | 43 | 0.662 / 0.366 | +0.296 [+0.100, +0.507] | 0.070 [0.000, 0.193] | 0.750 / 0.150 | yes |
| weak | one-sided, 8 st. (2024) | 0.25 | 43 | 0.690 / 0.380 | +0.310 [+0.123, +0.500] | 0.186 [0.036, 0.405] | 1.100 / 0.160 | yes |
| weak | 4 stations (7331) | 0.1 | 39 | 0.590 / 0.154 | +0.436 [+0.300, +0.579] | 0.179 [0.061, 0.324] | 0.810 / 0.220 | yes |
| weak | 4 stations (7331) | 0.25 | 39 | 0.628 / 0.179 | +0.449 [+0.319, +0.581] | 0.410 [0.163, 0.718] | 1.020 / 0.310 | yes |
| weak | 6 stations (7331) | 0.1 | 45 | 0.824 / 0.338 | +0.485 [+0.289, +0.662] | 0.222 [0.077, 0.411] | 0.870 / 0.230 | no |
| weak | 6 stations (7331) | 0.25 | 45 | 0.824 / 0.338 | +0.485 [+0.289, +0.662] | 0.422 [0.233, 0.684] | 1.170 / 0.300 | yes |
| weak | one-sided, 8 st. (7331) | 0.1 | 38 | 0.840 / 0.481 | +0.358 [+0.172, +0.565] | 0.053 [0.000, 0.182] | 1.180 / 0.150 | yes |
| weak | one-sided, 8 st. (7331) | 0.25 | 38 | 0.840 / 0.481 | +0.358 [+0.172, +0.565] | 0.263 [0.075, 0.526] | 1.460 / 0.190 | yes |
| strong | 4 stations (2024) | 0.1 | 40 | 0.221 / 0.130 | +0.091 [-0.033, +0.211] | 0.025 [0.000, 0.103] | 0.050 / 0.190 | yes |
| strong | 4 stations (2024) | 0.25 | 40 | 0.403 / 0.143 | +0.260 [+0.157, +0.371] | 0.200 [0.024, 0.500] | 0.290 / 0.340 | yes |
| strong | 6 stations (2024) | 0.1 | 44 | 0.419 / 0.135 | +0.284 [+0.162, +0.431] | 0.000 [0.000, 0.000] | 0.120 / 0.140 | yes |
| strong | 6 stations (2024) | 0.25 | 44 | 0.514 / 0.162 | +0.351 [+0.208, +0.513] | 0.023 [0.000, 0.111] | 0.340 / 0.260 | yes |
| strong | one-sided, 8 st. (2024) | 0.1 | 42 | 0.530 / 0.193 | +0.337 [+0.149, +0.547] | 0.048 [0.000, 0.159] | 0.410 / 0.140 | yes |
| strong | one-sided, 8 st. (2024) | 0.25 | 42 | 0.590 / 0.193 | +0.398 [+0.231, +0.612] | 0.190 [0.021, 0.458] | 0.680 / 0.270 | yes |
| strong | 4 stations (7331) | 0.1 | 37 | 0.115 / 0.038 | +0.077 [+0.000, +0.186] | 0.000 [0.000, 0.000] | 0.070 / 0.150 | yes |
| strong | 4 stations (7331) | 0.25 | 37 | 0.218 / 0.051 | +0.167 [+0.056, +0.310] | 0.108 [0.000, 0.381] | 0.260 / 0.220 | yes |
| strong | 6 stations (7331) | 0.1 | 36 | 0.500 / 0.238 | +0.262 [+0.062, +0.443] | 0.028 [0.000, 0.118] | 0.390 / 0.180 | yes |
| strong | 6 stations (7331) | 0.25 | 36 | 0.595 / 0.250 | +0.345 [+0.130, +0.526] | 0.167 [0.000, 0.400] | 0.620 / 0.220 | yes |
| strong | one-sided, 8 st. (7331) | 0.1 | 40 | 0.773 / 0.293 | +0.480 [+0.259, +0.664] | 0.075 [0.000, 0.275] | 0.560 / 0.230 | yes |
| strong | one-sided, 8 st. (7331) | 0.25 | 40 | 0.773 / 0.293 | +0.480 [+0.259, +0.664] | 0.375 [0.065, 0.907] | 0.940 / 0.320 | yes |

**Rule 2D-c holds under its primary reading, and its operating point is not deployable.** Under the
per-cell reading the rule passes in both regimes: with weak clutter in three and two of three cells at
$\alpha=0.1$ and in three and three at $\alpha=0.25$, and with strong clutter in three of three
everywhere (Table Y.47). Under the strict reading it passes with strong clutter and fails with weak
clutter by one false event in one cell: 10 false events in 45 empty windows, $0.222$ against
$2\alpha=0.2$, with 6 stations in seed 7331 at $\alpha=0.1$. The recall condition alone holds in all
twelve combinations of cell, seed and $\alpha$ of each regime, by $+0.08$ to $+0.49$, and re-applied in
the hierarchical replicates (exploratory) the per-cell reading holds in 76 per cent of them with weak
clutter and 96 per cent with strong clutter. Under the Holm reading, with a Holm rejection of the recall
difference in place of its point condition, the rule fails in both regimes. On its own 200 windows the
calibration behaves as Proposition AA.5(c) says: the fraction of empty windows with a false event of
MT-FW+SC is $0.095$ at $\alpha=0.1$ and $0.245$ at $\alpha=0.25$ in every cell. On the evaluation
windows the false events per empty window exceed $\alpha$ in three of six weak-clutter cells at
$\alpha=0.1$, by up to $2.2\alpha$, and in four of six at $\alpha=0.25$, by up to $1.7\alpha$, and with
strong clutter in none at $\alpha=0.1$ and one at $\alpha=0.25$. The evaluation networks are not
exchangeable with five calibration networks per cell, whose score distribution does not capture the
variation between networks, and the registered tolerance of $2\alpha$ absorbed that variation except in
one cell. Two qualifications keep this operating point from being deployable as a false-event budget.
First, a threshold calibrated on windows without events controls false alarms only there: at these
thresholds MT-FW+SC emits over all windows 0.73 to 1.18 false events per window with weak clutter at
$\alpha=0.1$ and 0.98 to 1.46 at $\alpha=0.25$, against 0.15 to 0.23 and 0.16 to 0.31 for QuakeMigrate,
three to eight times as many, and 0.05 to 0.56 and 0.26 to 0.94 with strong clutter against 0.14 to 0.23
and 0.22 to 0.34; part of the recall gain of rule 2D-c is bought with those false events. Second, the
guarantee assumes exchangeability, which, as the rates above show, holds only approximately across
networks. The conformal threshold answers whether a window contains an event, not how many false
events a catalogue carries; for the latter the thresholds of rule 2A-e, calibrated on all windows at a
false-event budget, remain the reference (Table Y.32).

**Exploratory diagnostics (post hoc).** The analyses of this paragraph and of Table Y.48 were not
registered, and no verdict depends on them. Where the false events of MT-FW+SC fall at the conformal
thresholds, in the Task 2A windows that contain events and that the empty-window calibration never
sees (`task2d_explore_e1.txt`): with weak clutter there are 1.21 to 1.87 false events per window with
events at $\alpha=0.1$ and 1.41 to 2.19 at $\alpha=0.25$, and 78 to 92 per cent of them lie 5 to 20 km
from a true source, residual energy of a real event left by deflation, the 5 km exclusion and the
model mismatch, against 0.18 to 0.36 per window for QM-3C+SC; with strong clutter MT-FW+SC emits 0.07
to 0.88 and 0.35 to 1.32, and QM-3C+SC 0.17 to 0.38. The conformal thresholds of MT-FW+SC, 23.7 to 29.8
with weak clutter and 36.2 to 52.1 with strong clutter, lie below those of rule 2A-e, 34.6 to 42.9 and
44.7 to 53.2, which were calibrated on all windows at a budget of false events per window; MT-FW+SC
recomputed for this diagnostic equals the Task 2A records in 713 of 713 windows with events. Applied to
the thirty networks per cell of Task 2C, with 150 to 169 empty windows per cell (`analyze_2d_out.txt`,
last table), the same frozen thresholds give MT-FW+SC minus QM-3C+SC in recall of $+0.36$ to $+0.47$
with weak clutter and $+0.03$ to $+0.40$ with strong clutter, and false events of MT-FW+SC per empty
window of 0.055 to 0.139 at $\alpha=0.1$ and 0.21 to 0.35 at $\alpha=0.25$ with weak clutter, and 0.006
to 0.092 and 0.10 to 0.33 with strong clutter: on more networks the rate is close to or somewhat above
$\alpha$, and below $2\alpha$ everywhere. The last diagnostic injects a timing error into data whose
timing is otherwise exact (Table Y.48): with the ideal generator of Y.11 and every source moved to its
nearest node with its origin time on the grid, as in Table Y.42, the detectors run once with no station
correction and once with a correction of $\pm0.2$ s of random sign per station and phase, 3.2 samples,
inside the window of the peak search.

**Table Y.48. Exploratory, post hoc. Recall at $\le0.25$ without / with an injected timing error of
$\pm0.2$ s per station and phase, ideal generator of Y.11 with sources on the nodes and origin times on
the grid, 90 windows per cell on seed-999 streams.** The first 30 windows of each cell are the
descriptive validation check V6. Output `task2d_explore_e2.txt`.

| Method | 4 stations (74 events) | 6 stations (63 events) | one-sided, 8 st. (67 events) |
|:--|:--|:--|:--|
| MT-FW | 0.81 / 0.68 | 0.83 / 0.86 | 0.87 / 0.82 |
| MT-FW-EX | 0.84 / 0.41 | 0.97 / 0.62 | 0.97 / 0.72 |
| MT-FW-NOCAP | 0.80 / 0.77 | 0.83 / 0.86 | 0.88 / 0.82 |
| BOMP | 0.80 / 0.39 | 0.97 / 0.76 | 0.99 / 0.72 |
| BOMP-PK | 0.80 / 0.77 | 0.83 / 0.86 | 0.88 / 0.82 |
| MTSTACK | 0.43 / 0.27 | 0.84 / 0.48 | 0.85 / 0.55 |
| MTSTACK-PK | 0.35 / 0.36 | 0.62 / 0.62 | 0.66 / 0.60 |
| GL-FW | 0.62 / 0.28 | 0.94 / 0.51 | 0.97 / 0.49 |
| GL-FW-PK | 0.66 / 0.51 | 0.83 / 0.76 | 0.88 / 0.75 |

The exact-time methods lose $0.16$ to $0.48$ with the error, and the peak-reading variants lose $-0.03$
to $+0.15$ and beat their exact counterparts with the error in 12 of 12 combinations of method and cell,
by $+0.04$ to $+0.38$: the peak search recovers what a timing error of this size removes. It has a
price. Without the error, on the 6- and 8-station cells the peak search costs recall: MT-FW reaches 0.83
to 0.87 against 0.97 for MT-FW-EX, BOMP-PK 0.83 to 0.88 against 0.97 to 0.99 for BOMP, and MTSTACK-PK
0.62 to 0.66 against 0.84 to 0.85 for MTSTACK, plausibly because a window of $\pm0.25$ s can lock onto
clutter or coda energy, which was not verified. MT-FW-NOCAP and BOMP-PK are again nearly identical.
Check V6, on the first 30 windows alone, showed the same picture with one exception in 12, MTSTACK-PK
in the one-sided cell ($-0.08$, two of 26 events).

**The edit after the results, and other deviations.** One line of `analyze_2d.py` was changed after
the part files existed: a line that computes diagnostics only failed because the lists of read-sample
shifts of BOMP-PK have length $2n_s$ per atom and needed flattening. The change touches that line
alone, and no rule, threshold, bootstrap or table computation changed. The hash of the script was
110aa89303cb106e8befb659beeb47b6a3fca7f8bfa9adb6f371886ca93373a7 when recorded at 02:41:59 UTC and is
8dd5013b9f57c1379523d138f58e1b0f66f94f49decb502555965d4a2527d2b4 after the edit; the result paragraph
of the registry states the change. The other departures are the readings listed above, the dense loop
used for MT-FW-EX and MT-FW-NOCAP, the explicit kernel, and a validation check V6 that is small and
descriptive.

**What Task 2D does not test.** The tolerance tested is one heuristic, a per-station and per-phase
argmax within a fixed $\pm0.25$ s, on the one grid of 3 km and 0.1 s; no method with time shifts as
nuisance parameters of the optimization, no other width and no finer grid was run. BOMP-PK+SC and
MTSTACK-PK+SC are implementations from the published descriptions with the tolerance of MT-FW added,
not the authors' codes. Part (i) has five networks per cell, which is why no rule survives Holm's
correction. The frozen thresholds of part (ii) are calibrated on empty windows only, and the
comparison over thirty networks of Task 2C was repeated with them only as an exploratory reading. The
data are the synthetic generator of Y.14, with double-couple sources inside the moment-tensor subspace
of the detector, and MALMI [110], the stacking of [112], real noise and real data remain untested.

### Y.21a What Task 2D changes

Task 2D changes the reading of Tasks 2A to 2C in three ways, and the second is the one that shapes the
thesis.

The attribution that Y.20 left to a comparison across methods is settled, against the distinctiveness
of MT-FW. Given the same tolerance, block orthogonal matching pursuit in the form of [129] matches
MT-FW, with differences of $-0.155$ to $+0.091$ at $\le0.25$ and every hierarchical interval containing
zero, and the only element of MT-FW beyond [129] so equipped, the per-station cap, has a small effect
that is not significant. MT-FW is not a better detector than [129] with equal tolerance. The reading of
rule 2C-a in Y.20 and Y.20a, that the advantage belongs to the elements proper to MT-FW rather than to
the moment tensor with deflation already present in [129], does not survive: it belongs to the peak
search, which [129] can be given. The empirical advantage over pick-based associators and QuakeMigrate
on these synthetics, which Tasks 2A and 2C established, belongs to the family of timing-tolerant
moment-tensor detectors with deflation, which includes MT-FW and block orthogonal matching pursuit
with the peak search; mechanism-corrected stacking with the peak search, which has no deflation, stays
below both.

The ingredient that matters is tolerance to timing error, from the discretization of the grid, the
residual station delays and the mismatch of the velocity model. Removing it costs MT-FW $0.06$ to
$0.41$ in all twelve cases, and with an injected error it recovers most of what the exact-time methods
lose; with exact timing and sources on the nodes it costs $0.10$ to $0.14$ on six and eight stations.
At present that tolerance is a heuristic: an argmax per station and phase within a fixed window, with
no model of the timing error, no penalty on the shifts and no certificate over them, and the group
BLASSO given the same heuristic stays $0.11$ to $0.30$ below MT-FW. The distinctive research
contribution that remains for the method is therefore the principled treatment of that tolerance: the
group BLASSO (16a) with a time shift per station and phase as a nuisance parameter, bounded by a
declared tolerance, solved off the grid by sliding Frank-Wolfe [116] that moves positions, origin times
and shifts continuously, with its certificate over shifts, its duality gap and its theory. It counts as
a methodological advance only if, on the windows of Task 2A, it beats both the heuristic peak search of
MT-FW and block orthogonal matching pursuit given that search, BOMP-PK+SC, which is the baseline to
beat; if it does not, what remains is the benchmark and the real-data test. Task 2E tested it, and it
does not (Y.22, Y.24a).

A frozen threshold exists, and it is not yet a deployable operating point. The split-conformal
threshold calibrated on empty windows of other networks passes rule 2D-c under its primary reading and
keeps the false rate on empty windows below $2\alpha$ in every combination of cell, seed and $\alpha$ but
one, but it controls nothing in
windows with events, where MT-FW+SC emits 1.2 to 2.2 false events per window with weak clutter, and
over all windows it emits three to eight times the false events of QuakeMigrate there. A deployable
operating point for a false-event budget would need a calibration that sees windows with events, as
rule 2A-e did, or an error criterion per event rather than per window; neither exists with a result
that passes, on synthetics or on real data.

What Y.18a to Y.20a list as not established remains so: the operating points of part (i) are chosen
on the evaluation windows, rule 2A-e fails under its strict reading, GaMMA with optimistic picks is
ahead at 0.1 with weak clutter, the sources are double couples inside the detector's subspace, the one
synthetic generator lies outside the declared regime of Section 1.5, and there is no result on real
data.

## Y.22 The off-grid time-shift BLASSO (Task 2E)

Code: `task2e_methods.py` (the continuous atom model, the oracle, the refinement, the sliding step and the
two continuous methods), `task2e_validate.py` (validation on seed 999), `task2e_run.py` (validation,
calibration, timing estimate and runs), `analyze_2e.py` (rules, tables, hierarchical bootstrap and
Holm's correction) and `task2e_explore.py` (exploratory diagnostics written after the results); outputs
`task2e_validation.txt`, `task2e_calibration.txt`, `task2e_frozen.json`, `task2e_smoke.txt`,
`task2e_code_hashes.txt`, `analyze_2e_out.txt`, `analyze_2e_out.json` and `task2e_explore.txt`; design
and rules in `PREREGISTRO_1j_1k.txt` (section "Task 2E") with the hashes of the registry in
`PREREGISTRO_hashes.txt`; review of the implementation, calibration, compute and deviations in
`task2e_report.txt`; independent verification in `verif_2e/task2e_verification.txt`.

**Why the task was run.** Task 2D left one methodological question (Y.21a): whether the principled form
of the timing tolerance, the group BLASSO (16a) with a time shift per station and phase as a nuisance
parameter, solved off the grid by sliding Frank-Wolfe [116], detects more events than block orthogonal
matching pursuit given the heuristic peak search, BOMP-PK+SC, at equal tolerance and equal realized
false-event rate. That is hypothesis G1-M of Section 4.3, and Task 2E is its registered test.

**Preregistration.** The design, the rules and their interpretation were appended to the registry at
14:44:14 UTC on 27 September 2026, before any code of the task was written. The SHA-256 hash of the
registry at that moment, b40f0c6ead272610020b8c93d832c035a006d705cfd910bb6dc9dbf8f4d5e988, was recorded in
`PREREGISTRO_hashes.txt` and communicated outside the computing environment, in the working session,
at the moment of registration; it is the first registration whose content is fixed by a record held
outside the container, although that record is neither public nor issued by a third party. The one
change allowed before evaluation, the calibration note, was appended at 15:38:20 UTC with the frozen
value $f=0.01$, before any computation on seeds 2024 or 7331, and hashed at the same second
(fada69b1f1f6d8534729a9c6f8be81f6acccbb61ba5628c0d14fbf23c5c59271). The code was written from 14:55:56
UTC, validated by 15:07:22 UTC, and the hashes of its five files were recorded at 15:38:31 UTC in
`task2e_code_hashes.txt`; the first evaluation part was written at 15:39:40 UTC and the last at
17:39:19 UTC, and `sha256sum -c` on the code hashes passed after the run. The result paragraph was
appended at 17:45:17 UTC with a further hash. The budget rule, which would have reduced the evaluation
to fifteen networks per cell had the seed-999 estimate exceeded 10 hours on two cores, did not trigger:
the estimate was 1.98 hours.

Rule 2E-a (primary), per clutter regime: TS-SFW+SC minus BOMP-PK+SC at $\le0.25$ realized false events
per window is at least $0.05$ with the lower bound of the 95% hierarchical interval above zero, in at
least two of three cells in both seeds; a strict reading applies Holm's correction [142] over the six
tests of cell and seed of the regime at the one-sided level 0.025. Rule 2E-b (secondary), per regime:
TS-SFW+SC minus BOMP-PK-SLIDE+SC is at least $0.05$ at $\le0.25$ in at least two of three cells in both
seeds, on the point estimate. Interpretation, fixed in the registration: if 2E-a fails in both regimes,
G1-M is refuted in this synthetic design, and the proposal says so and reformulates its methodological
contribution; if 2E-a holds and 2E-b fails, the gain belongs to sliding off the grid and not to the
group penalty.

**Design.** The data are, per seed, regime and cell, thirty networks: the five of Task 2A with twenty
windows each and the twenty-five of Task 2C part (ii) with twelve windows each, 400 windows per cell and
4800 in all, regenerated bit for bit with the kernel of the numerical library pinned as in Tasks 2C and
2D, with the same station corrections, the grid of 675 nodes at 3 km with depths of 5, 9 and 13 km, the
origin times of Task 2A, the tolerance of $\pm0.25$ s, the 5 km exclusion, at most five events per
window and the 5 km epicentral hit rule. An atom of TS-SFW+SC has a continuous position and origin time
and a shift $\delta_{i,b}$ per station $i$ and phase $b$, bounded by $\pm0.25$ s; its arrival time is
the predicted travel time plus the origin time, the station correction and the shift, with no rounding
to samples, and its prediction in the matched-filter domain is the radiation operator evaluated at the
continuous position times the continuous autocorrelation of the matched filter centred at those times.
The correlation with the data is the exact band-limited correlation with the continuously shifted
wavelet, which equals the matched-filter output at sample times to $8.9\times10^{-16}$. The objective is
the whitened least-squares term of (16a) with the block-whitened group norm of GL-FW (Y.20) and analytic
gradients in all positions, origin times and shifts. Each of at most five iterations (a) takes the
maximum over nodes and origin times not excluded, 5 km around the current positions of the accepted
atoms, of the uncapped projection energy of the peak vectors read within $\pm0.25$ s, which is the
statistic of MT-FW-NOCAP+SC; (b) refines the new atom by L-BFGS-B from that node, with the initial
shifts at the read peaks, at most 50 iterations; (c) slides all atoms, alternating the convex step in
the amplitudes, solved to a relative duality gap of at most $10^{-6}$, with L-BFGS-B in positions,
origin times and shifts, at most 50 iterations per step and five alternations, followed by a final
convex step; and (d) subtracts the joint prediction. The penalty is $\lambda=f\,s_1$, with $s_1$ the
first oracle value of the window and $f$ calibrated. The score is the running minimum of the oracle
values at entry, and a hit is judged at the final position of each atom. BOMP-PK-SLIDE+SC is the same
iteration with $\lambda=0$, a least-squares refit, and separates the group penalty from sliding.
BOMP-PK+SC is the code of Task 2D unchanged, whose recomputation equals the Task 2D records in 1200 of
1200 windows, and MT-FW+SC the stored records of Tasks 2A and 2C.

**Validation and calibration.** Before any evaluation seed was run (`task2e_validation.txt`), the
continuous autocorrelation, the vectorized radiation operator, the band-limited correlation and the
continuous Gram matrix reproduced their sampled counterparts to at most $8.9\times10^{-16}$; the analytic gradients agreed with central differences to a relative
$1.4\times10^{-9}$ in 18 checks; the convex step agreed with an explicit data-domain dictionary,
reproduced the least-squares refit at $\lambda=0$ to $5.0\times10^{-11}$ and matched an independent
proximal-gradient solve; and with refinement and sliding disabled and $\lambda=0$ the first atom
equalled that of BOMP-PK in 24 of 24 windows, and with the oracle statistic of BOMP-PK the whole record
did. A noise-free check with one off-grid source and
injected shifts exposed a property of the method rather than a defect: started at the truth the
optimizer stays there (94 of 94), and started 0.2 to 0.4 km away it recovers the source in 85 of 94
cases, but started from the grid node it reaches the source within 0.5 km with every shift within
0.02 s in only 12 of 90 sources at $\lambda=0$ and at $f=0.01$, with a median three-dimensional error
of 0.52 km against 2.45 km for the grid node of BOMP-PK. The failures are cycle skips of about one
wavelet period, 0.45 s, on some stations and phases. The registered method is a local refinement by
definition, so the calibration and evaluation went ahead as registered. The calibration on the 600
seed-999 windows of Task 2A (`task2e_calibration.txt`) gave recall at $\le0.25$ of 0.4749, 0.4644,
0.2510 and 0.1841 for $f=0.01$, 0.03, 0.1 and 0.3, and $f=0.01$ was frozen. At that value $\lambda$ is a
median 0.045 of the certificate of the first atom and no atom is set to zero; at $f=0.3$, 88 per cent of
the atoms are, which stops the deflation. The calibrated penalty lies at the lower edge of the
registered grid, where TS-SFW+SC is almost a least-squares refit.

**Table Y.49. Event recall at realized false events per window $\le0.1$ / $\le0.25$ on the 4800
windows of Task 2E (400 per cell, thirty networks).** Script `analyze_2e.py`; output
`analyze_2e_out.txt`, Table 1 of each regime.

| Regime | Cell (seed) | events | TS-SFW+SC | BOMP-PK-SLIDE+SC | BOMP-PK+SC | MT-FW+SC |
|:--|:--|:--|:--|:--|:--|:--|
| weak | 4 stations (2024) | 308 | 0.334/0.558 | 0.328/0.562 | 0.279/0.519 | 0.458/0.562 |
| weak | 6 stations (2024) | 313 | 0.473/0.597 | 0.473/0.601 | 0.438/0.569 | 0.492/0.623 |
| weak | one-sided, 8 st. (2024) | 314 | 0.465/0.634 | 0.490/0.634 | 0.439/0.608 | 0.548/0.656 |
| weak | 4 stations (7331) | 307 | 0.352/0.485 | 0.332/0.472 | 0.270/0.469 | 0.319/0.466 |
| weak | 6 stations (7331) | 308 | 0.568/0.701 | 0.584/0.708 | 0.516/0.672 | 0.568/0.685 |
| weak | one-sided, 8 st. (7331) | 326 | 0.580/0.687 | 0.543/0.699 | 0.528/0.672 | 0.537/0.638 |
| strong | 4 stations (2024) | 325 | 0.212/0.274 | 0.206/0.277 | 0.197/0.274 | 0.175/0.249 |
| strong | 6 stations (2024) | 304 | 0.316/0.474 | 0.316/0.467 | 0.306/0.441 | 0.326/0.431 |
| strong | one-sided, 8 st. (2024) | 305 | 0.413/0.511 | 0.413/0.531 | 0.377/0.495 | 0.367/0.489 |
| strong | 4 stations (7331) | 313 | 0.179/0.243 | 0.182/0.246 | 0.173/0.240 | 0.131/0.208 |
| strong | 6 stations (7331) | 329 | 0.331/0.474 | 0.325/0.480 | 0.310/0.459 | 0.347/0.459 |
| strong | one-sided, 8 st. (7331) | 329 | 0.486/0.587 | 0.517/0.596 | 0.438/0.568 | 0.453/0.556 |

**Table Y.50. Rules 2E-a and 2E-b. TS-SFW+SC minus BOMP-PK+SC at both caps and TS-SFW+SC minus
BOMP-PK-SLIDE+SC at $\le0.25$, point [95% hierarchical interval], with the one-sided bootstrap
$p$-value of 2E-a at $\le0.25$.** Same output, Tables 2 and 3.

| Regime | Cell (seed) | 2E-a, $\le0.1$: point [hierarchical] | 2E-a, $\le0.25$: point [hierarchical] | one-sided $p$ at $\le0.25$ | 2E-b, $\le0.25$: point [hierarchical] |
|:--|:--|:--|:--|:--|:--|
| weak | 4 stations (2024) | +0.055 [-0.068, +0.138] | +0.039 [-0.025, +0.122] | 0.135 | -0.003 [-0.036, +0.028] |
| weak | 6 stations (2024) | +0.035 [-0.033, +0.170] | +0.029 [-0.025, +0.087] | 0.149 | -0.003 [-0.048, +0.030] |
| weak | one-sided, 8 st. (2024) | +0.025 [-0.089, +0.110] | +0.025 [-0.042, +0.094] | 0.220 | +0.000 [-0.040, +0.026] |
| weak | 4 stations (7331) | +0.081 [-0.028, +0.191] | +0.016 [-0.038, +0.094] | 0.327 | +0.013 [-0.020, +0.045] |
| weak | 6 stations (7331) | +0.052 [-0.039, +0.147] | +0.029 [-0.024, +0.078] | 0.235 | -0.006 [-0.038, +0.018] |
| weak | one-sided, 8 st. (7331) | +0.052 [-0.050, +0.121] | +0.015 [-0.035, +0.091] | 0.241 | -0.012 [-0.042, +0.026] |
| strong | 4 stations (2024) | +0.015 [-0.032, +0.050] | +0.000 [-0.033, +0.030] | 0.567 | -0.003 [-0.024, +0.014] |
| strong | 6 stations (2024) | +0.010 [-0.026, +0.086] | +0.033 [-0.013, +0.083] | 0.098 | +0.007 [-0.026, +0.026] |
| strong | one-sided, 8 st. (2024) | +0.036 [-0.019, +0.125] | +0.016 [-0.026, +0.059] | 0.313 | -0.020 [-0.045, +0.007] |
| strong | 4 stations (7331) | +0.006 [-0.035, +0.045] | +0.003 [-0.027, +0.039] | 0.427 | -0.003 [-0.020, +0.006] |
| strong | 6 stations (7331) | +0.021 [-0.014, +0.104] | +0.015 [-0.024, +0.058] | 0.361 | -0.006 [-0.038, +0.036] |
| strong | one-sided, 8 st. (7331) | +0.049 [-0.014, +0.115] | +0.018 [-0.044, +0.080] | 0.323 | -0.009 [-0.045, +0.010] |

**Rules 2E-a and 2E-b fail in both regimes.** At $\le0.25$ the off-grid method gains $+0.000$ to
$+0.039$ over BOMP-PK+SC; all twelve point estimates are non-negative, none reaches $0.05$, and every
hierarchical interval contains zero, so no cell passes in either seed or regime (Table Y.50). No case
lies within 0.005 of a threshold. Under the strict reading the smallest one-sided $p$-value is 0.135
with weak clutter and 0.098 with strong clutter, and the smallest Holm-adjusted value 0.81 and 0.58.
By the registered interpretation G1-M is refuted in this synthetic design: the group BLASSO with a time
shift per station and phase, solved off the grid by sliding Frank-Wolfe, is not a better detector than
block orthogonal matching pursuit with the same $\pm0.25$ s tolerance. Rule 2E-b fails as well, with
differences of $-0.020$ to $+0.013$: at the calibrated $f$ the group penalty is inert, and the small
positive differences of 2E-a come from sliding and the continuous refit. Descriptively, sliding alone,
BOMP-PK-SLIDE+SC minus BOMP-PK+SC, gives $+0.003$ to $+0.042$ at $\le0.25$ with every interval
containing zero, and TS-SFW+SC minus MT-FW+SC lies between $-0.026$ and $+0.049$ at $\le0.25$, again
with every interval containing zero, and between $-0.019$ and $-0.123$ at $\le0.1$ with weak clutter in
seed 2024, where the cap of MT-FW favours it at the strict end. Location improves marginally: the median
epicentral error of the hits at the operating point of $\le0.25$ is 1.79 km against 1.88 km for BOMP-PK+SC
with weak clutter and 2.01 km against 2.15 km with strong clutter, and the continuous solutions sit
0.3 to 1.3 km too deep in most cells. The convex steps met their gap everywhere (largest relative gap
$5.6\times10^{-10}$), but 73 to 83 per cent of the sliding steps and 83 to 97 per cent of the
refinements stopped at the 50-iteration cap, 24 to 35 per cent of the shifts ended at the $\pm0.25$ s
bound, and a window cost 1.19 s against 0.42 s for BOMP-PK+SC.

**Independent verification.** A separate verification (`verif_2e/task2e_verification.txt`) asked
whether the negative result is an implementation error, a budget handicap or a trustworthy negative,
with new computation on seed 999 only and no change to any file of the task or of the registry. Code
review and independent code found no implementation error: with the position on a node and the shifts
at the read peaks, the continuous atom reproduces the BOMP-PK block over 120 atoms of real windows,
prediction to $1.4\times10^{-15}$, Gram matrix to $9.0\times10^{-16}$, correlation to
$5.7\times10^{-16}$ and joint least-squares objective to $9.9\times10^{-16}$; the objective agrees with an
independent data-domain implementation to $2.3\times10^{-13}$ and its gradient with finite differences
to $4.5\times10^{-8}$; injected shifts of known sign are recovered to $10^{-4}$ s with a position error
of 0.03 km; and station corrections, phase windows, bounds and the gauge between origin time and a
common shift are consistent. The one asymmetry found, a peak window of BOMP-PK offset by the fractional
part of the predicted sample, affects 0.126 of the peak offsets by less than one sample. The verifier's
own exact FROC reproduced every recall and difference of Tables Y.49 and Y.50, and its own bootstrap the
intervals. The budget was not a handicap. Rerun with ten times the iterations, four times the alternations and a
nine-start refinement on 270 windows with 208 events, no L-BFGS-B call reached its cap and the objective
after five atoms fell by a median 1.9 per cent, but recall at $\le0.25$ changed by $-0.014$
[$-0.052$, $+0.038$], with 6 events gained and 4 lost: the registered budget truncates the optimizer
numerically without handicapping detection. The registry chain is intact: each of the six hashes then
recorded is the hash of a prefix of the current registry, the prefixes grow monotonically, the design
text is byte-identical to the version hashed at 14:44:14 UTC, and the calibration note precedes the
first evaluation part.

**Table Y.51. Verification of Task 2E. Oracle ceilings on the 600 seed-999 windows (100 per cell and
regime): recall at $\le0.1$ / $\le0.25$ and headroom over BOMP-PK+SC at $\le0.25$ [95% hierarchical
interval].** CEIL-G scores each true event by the best BOMP-PK block statistic within 5 km of it and
takes as false candidates the greedy block maxima of the exact clutter-and-noise part of the window, a
generous ceiling for any grid block detector; CEIL-GD takes the false candidates after an oracle joint
least-squares deflation of the true events, so it keeps the false events caused by model mismatch, and is
the realistic ceiling; CEIL-C is the continuous version, started at the true source. Output
`verif_2e/task2e_verification.txt`, section 4.

| Regime | Cell | events | BOMP-PK+SC | TS-SFW+SC | CEIL-G | CEIL-GD | CEIL-C | CEIL-G minus BOMP-PK+SC, $\le0.25$ | CEIL-GD minus BOMP-PK+SC, $\le0.25$ |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| weak | 4 stations | 81 | 0.309/0.531 | 0.321/0.568 | 0.840/0.852 | 0.605/0.753 | 0.840/0.864 | +0.321 [+0.128, +0.607] | +0.222 [+0.041, +0.495] |
| weak | 6 stations | 67 | 0.433/0.657 | 0.463/0.657 | 0.955/0.955 | 0.672/0.836 | 0.970/0.985 | +0.299 [+0.148, +0.542] | +0.179 [-0.042, +0.418] |
| weak | one-sided, 8 st. | 77 | 0.571/0.623 | 0.636/0.675 | 0.844/0.870 | 0.675/0.701 | 0.857/0.883 | +0.247 [+0.118, +0.351] | +0.078 [-0.039, +0.178] |
| strong | 4 stations | 74 | 0.176/0.230 | 0.189/0.230 | 0.284/0.392 | 0.230/0.338 | 0.284/0.351 | +0.162 [+0.028, +0.375] | +0.108 [+0.000, +0.352] |
| strong | 6 stations | 97 | 0.175/0.309 | 0.268/0.381 | 0.608/0.680 | 0.371/0.546 | 0.649/0.742 | +0.371 [+0.162, +0.651] | +0.237 [+0.067, +0.505] |
| strong | one-sided, 8 st. | 82 | 0.439/0.512 | 0.488/0.524 | 0.683/0.732 | 0.634/0.671 | 0.695/0.744 | +0.220 [+0.042, +0.398] | +0.159 [+0.022, +0.352] |
| weak | three cells pooled | 225 | 0.422/0.582 | 0.471/0.631 | 0.876/0.880 | 0.631/0.738 | 0.867/0.893 | +0.298 [+0.196, +0.403] | +0.156 [+0.042, +0.286] |
| strong | three cells pooled | 253 | 0.245/0.356 | 0.265/0.395 | 0.447/0.601 | 0.375/0.518 | 0.498/0.617 | +0.245 [+0.101, +0.392] | +0.162 [+0.037, +0.301] |

**Where the headroom is.** The benchmark is not saturated. BOMP-PK+SC is not within $0.05$ of a ceiling
in any cell: at $\le0.25$ it lies $0.16$ to $0.37$ below CEIL-G and $0.08$ to $0.24$ below CEIL-GD
(Table Y.51), with pooled intervals that exclude zero in both regimes although some single-cell
intervals, with 100 windows, do not. The failure of 2E-a is therefore not an artefact of a benchmark with
no room for a gain of $0.05$. CEIL-C is about equal to CEIL-G, 0.699 against 0.692 at $\le0.25$ over all
600 windows: scoring with the continuous off-grid energy raises the scores of clutter and noise as much
as those of events. The headroom has a specific form. For 89 per cent of the events, the grid score near
the true source exceeds every clutter-and-noise candidate of the window, yet BOMP-PK+SC ever places an
atom within 5 km of only 77 per cent of those events, and TS-SFW+SC of 74 per cent; among the detectable
events never hit, 99 for BOMP-PK+SC and 110 for TS-SFW+SC, the nearest atom lies a median 7.3 and 6.8 km
away, and between 5 and 10 km away for 80 and 85 per cent of them. The events are absorbed by atoms
mislocated just beyond the 5 km hit tolerance. A refinement started at the true source drifts a median
2.35 km (90th percentile 6.7 km) under the mismatch between the generator and the homogeneous model of
the detector, while the registered method moves atoms a median 0.7 km from their entry node and more
than 5 km in only 0.4 to 0.6 per cent of them, so a local refinement from the grid cannot reach events
whose grid maximum is 5 to 10 km away. The verification concludes that the negative is trustworthy and
that a method that closes this headroom would need a better global placement of atoms, not local sliding;
Appendix Y.23 tests, exploratorily, whether a better travel-time or source model supplies it.

**What Task 2E does not test.** The registered method is one local refinement from a grid of 3 km with
data at 2 Hz, scored by the grid oracle and solved within a budget that the verification showed not to
bind for detection; global relocation, a finer starting grid, a penalty on the shifts other than the
bound and other values of the tolerance were not run. The ceilings are constructs with known source
positions and exact clutter-and-noise data, computed on 100 windows per cell of seed 999, so their
single-cell intervals are wide. The data are the synthetic generator of Y.14, and nothing here concerns
real data.

## Y.23 Where the headroom is not: seven exploratory candidates on the calibration seed

Code and outputs: `explore_contrib/` (scripts `c1_run.py` to `c7_run.py` with their analysis scripts,
`c2_ablate.py`, `c6_lib.py`, `c6_ceil.py`, `c6b_run.py`, `c7_lib.py`, `leak_test.py`, `leak2.py`,
`power.py` and `c7_power.py`; outputs in `results/explorations/out/`); reports `results/explorations/REPORT.txt`
and `results/explorations/REPORT_C7.txt`.

**Status.** Everything in this section is exploratory. It was run after Task 2E on seed 999 only, on
networks of that seed not used before (networks 10 to 15, 20 to 29, 30 to 35 and 40 to 59 of each cell;
the calibration networks 0 to 4 were only read), with a new random stream for every design, and the
seeds of every registered evaluation were neither generated nor read. No number here tests a registered
hypothesis, several were looked at while the designs were being chosen, and they are pilot estimates.
The question was where the headroom of Table Y.51 lies, and which candidate, if any, deserves a
registered test. The common protocol is that of Task 2E: the grid of 675 nodes, origin times from 5.0
to 8.0 s at 0.1 s, the peak search of $\pm0.25$ s, the 5 km exclusion and hit rule, the exact FROC and a
hierarchical bootstrap over networks and windows with the operating point re-chosen per replicate;
all differences are paired on the same windows. About 2.7 plus 2.6 hours of CPU time were used.

**Table Y.52. Exploratory candidates on seed 999: design, main contrast at $\le0.1$ / $\le0.25$ [95%
hierarchical interval] and reading.** Outputs `results/explorations/REPORT.txt`, sections 0 to 6, and
`results/explorations/REPORT_C7.txt`, sections 0, 3, 4 and 7.

| Candidate | Design | Main contrast | $\le0.1$ / $\le0.25$ | Reading |
|:--|:--|:--|:--|:--|
| C1 coarse grids | weak clutter, 180 windows, 145 events; grids of 3, 6 and 9 km | BOMP-PK+SC on a 9 km grid minus on the 3 km grid | $-0.324$ [$-0.524$, $-0.052$] / $-0.379$ [$-0.512$, $-0.174$] | off-grid sliding on the same coarse grid recovers $+0.03$ to $+0.10$ and never reaches BOMP-PK+SC on the 3 km grid; null |
| C2 dense overlaps | 3 to 5 events within 9 km, $k_{\max}=7$; 360 windows, 1439 events, both regimes | MT-FW+SC minus BOMP-PK+SC | $+0.074$ [$+0.016$, $+0.136$] / $+0.070$ [$+0.020$, $+0.128$] | the one candidate that met the exploratory criterion; the cap explains it (ablation below) |
| C2 on fresh networks | the dense design on networks 40 to 59, 600 windows, 2389 events | MT-FW+SC minus BOMP-PK+SC | $+0.077$ [$+0.030$, $+0.122$] / $+0.068$ [$+0.029$, $+0.109$] | replication of C2 before registration |
| C3 large timing errors | weak clutter, 180 windows; statics of 0.3 to 0.8 s with and without catalogue corrections | wider peak search (BOMP-PK-12) minus BOMP-PK, 0.8 s statics, no corrections | $+0.156$ [$+0.015$, $+0.284$] / $+0.148$ [$+0.000$, $+0.315$] | helps only without corrections; with catalogue corrections nothing beats BOMP-PK+SC; null |
| C6 physics travel times | 600 windows, 475 events; oracle and JHD-calibrated 1D models inside both detectors | oracle medium minus homogeneous, BOMP-PK+SC | $+0.006$ [$-0.066$, $+0.065$] / $+0.011$ [$-0.059$, $+0.064$] | null; calibrated gradient model $-0.104$ [$-0.202$, $+0.013$] at $\le0.25$ with strong clutter |
| C7 source-pulse invariance | 600 standard windows, 465 events; phase-invariant rank-one subspace (PH) | BOMP-PK-PH minus BOMP-PK+SC | $-0.114$ [$-0.203$, $-0.017$] / $-0.084$ [$-0.148$, $-0.023$] | negative; no gain in dense windows or with two out-of-family pulses |
| Arrival masking | each read arrival explains at most one event; variant fixed before a replication on networks 30 to 35 | BOMP-PK-MASK minus BOMP-PK+SC, weak clutter | $+0.149$ [$-0.007$, $+0.265$] / $+0.084$ [$+0.000$, $+0.163$] | borderline, weak clutter only; strong clutter $+0.014$ at $\le0.25$ |
| Speed | standard windows, 600 | MT-FW+SC minus BOMP-PK+SC | $+0.017$ [$-0.061$, $+0.092$] / $+0.002$ [$-0.050$, $+0.055$] | equal recall at 6.3 times lower cost; engineering, not a method |

**C1 and C3: the grid and the timing.** On a grid of 9 km the grid methods lose about 0.37 of recall at
$\le0.25$, and the off-grid TS-SFW of Task 2E, with $f$ frozen there, recovers only $+0.03$ to $+0.10$ of
it on the same grid, with intervals that mostly include zero, while costing about 1.2 s per window on
every grid; it is never better than BOMP-PK+SC on the 3 km grid. Widening the peak search, sharing
shifts within a window or learning station terms online helps only with 0.8 s statics and no catalogue
corrections; with the twenty-event corrections, which also absorb a uniform 7.5 per cent velocity error,
no variant beats BOMP-PK+SC. The discretization and the residual timing error that the peak search
absorbs are not where the headroom is.

**C6: the velocity model is not where the headroom is.** Replacing the homogeneous model inside both
detectors, travel times and radiation operator together, by the generator's own gradient medium with
curved rays (the oracle), or by homogeneous, gradient and two-layer models calibrated by joint
hypocentre determination on the twenty catalogued events, changes recall by amounts whose intervals
contain zero in both regimes, and the calibrated models hurt with strong clutter, where the catalogue
picks are contaminated. Only the oracle medium together with the true statics gains, $+0.069$ [$+0.004$,
$+0.141$] at $\le0.25$ pooled, and that is not attainable. The medium mismatch is real and of the order
of the tolerance, a mean of 0.12 s of which a constant per station and phase removes almost nothing, but
the peak search already absorbs it. An attribution run removed one mismatch at a time from the generator,
with sources, mechanisms, statics, clutter and noise shared (C6b, 150 windows per regime and variant).
With weak clutter the realistic headroom of BOMP-PK+SC at $\le0.25$ is $+0.129$ [$-0.037$, $+0.228$] in
the default generator and stays $+0.194$ [$+0.081$, $+0.312$] with an exact homogeneous medium, the
detector's own, $+0.177$ with a matched 2 Hz pulse, $+0.153$ without site and polarization perturbations
and $+0.242$ without coda: no single mismatch closes it, and with an exact medium the detectors still
miss what the ceiling finds. With strong clutter the headroom on these networks is about zero in every
variant.

**C7: the source pulse is not where the headroom is either.** Making the moment-tensor subspace invariant
to the source pulse, with a phase-invariant rank-one statistic over the 2 Hz Ricker and its Hilbert
transform (PH) or over a four-dimensional basis that captures at least 0.95 of the energy of the whole
pulse family of the generator (PHF), lowers recall in standard windows, $-0.084$ and $-0.080$ at
$\le0.25$ with intervals below zero, is neutral or negative in dense windows ($-0.021$ for PH and $-0.046$
for PHF) and gives nothing with two pulse families
outside the design, Gabor and Berlage pulses with phases up to $\pm90$ degrees ($-0.004$ and $+0.026$ at
$\le0.25$ for PH, intervals containing zero). Noise-free, with the full realism of the generator, the
true block of the 2 Hz detector captures only about half of the signal energy (0.497), and pulse
invariance adds 0.014 to 0.050 of it: the pulse is a small part of the mismatch next to the amplitude factors,
polarization deviation, medium and residual statics. The extra degrees of freedom are paid for in spatial
selectivity: the first atom lands beyond 5 km more often (first-atom hit rate 0.649 for BOMP-PK+SC, 0.609
for PH and 0.577 for PHF in standard windows). A by-product is exact: without the cross term of P and S in
its Gram matrix, block pursuit reproduced the hits, locations and scores of BOMP-PK+SC in 2400 of 2400
windows to a relative $10^{-12}$ at 10.7 times lower cost, so that term is immaterial in this geometry.

**C2: the cap, and why it was selected.** In windows with three to five events within 9 km, MT-FW+SC
exceeded BOMP-PK+SC by $+0.070$ at $\le0.25$ pooled over both regimes, and in standard windows the two
were equal. The two detectors differ in three ways: MT-FW clips each station and phase peak vector at
norm 6 in noise units before projecting on the moment-tensor subspace, it has no cross term of P and S in
its statistic, and it deflates sequentially instead of refitting all atoms jointly. The ablation on the
same windows isolated the first. MT-FW without the cap produced hit sequences identical to those of
BOMP-PK with sequential deflation in 360 of 360 windows, MT-FW exceeded MT-FW without the cap by
$+0.071$ [$+0.010$, $+0.163$] at $\le0.25$ with weak clutter, and adding the cap to the selection statistic of
BOMP-PK recovered most of the gap; joint refitting and the cross term did not change recall. The pilot
reading was that a few near-station or strongly amplified arrivals dominate the unclipped statistic in
overlapping sequences and pull atoms to wrong nodes. Arrival masking, the other candidate with a
replicated positive point estimate, acts in a different setting, isolated events with weak clutter, and
its interval touches zero.

**Multiplicity.** Seven candidates were examined: C1, C2, C3, C6 and C7, arrival masking with twelve
variants in its pilot, and speed; two further designs, C4 and C5, were scripted and never run. The per-station cap in dense sequences was one of
five comparisons in the C2 pilot and was selected for a registered test because it alone met the
exploratory criterion of a gain of at least $0.05$ with an interval excluding zero; masking was the best
of twelve variants. The pilot estimates of the selected candidate are subject to the winner's curse and
are not evidence for it. What protects the registered test of Y.24 from that selection is that its seeds,
4242 and 8888, were untouched by the exploration, its design and code were frozen before they were
generated, and its planning effect was half the pilot estimate.

## Y.24 The per-station cap in dense sequences (Task 2F)

Code: `task2f_run.py` (controls and runs), `task2f_control.py` (controls on seed 999), `analyze_2f.py`
(rules, tables, hierarchical bootstrap and Holm's correction) and the frozen methods of
`explore_contrib/c6_lib.py` and `explore_contrib/c2_ablate.py`; outputs `analyze_2f_out.txt`,
`task2f_control.txt`, `task2f_code_hashes.txt`, `task2f_freeze_time.txt` and `log2f.txt`; design and rules
in `PREREGISTRO_1j_1k.txt` (section "Task 2F") with the hashes in `PREREGISTRO_hashes.txt`; report
`task2f_report.txt`; independent verification in `verif_2f/task2f_verification.txt`.

**Preregistration and freeze.** The question, H-CAP, is whether in windows with three to five overlapping
events MT-FW+SC, the moment-tensor subspace statistic with a per-station cap and sequential deflation,
detects more events than BOMP-PK+SC at equal realized false-event rate, and whether the gain is due to the
cap. It was appended to the registry at 03:34:51 UTC on 28 September 2026, before any code was run on
seeds 4242 and 8888, which had never been used, and the hash of the registry,
3d8e4c9b2f7106c34c7d1f832bd4802db2ac119f48cf27630847640b82da2f5f, was recorded and communicated outside the
computing environment in the working session before the run. The registration declares that the
hypothesis was selected among the candidates of Y.23. The controls on seed 999 started at 03:39:10 UTC;
the SHA-256 hashes of the 26 local modules loaded by the runner, the controls and the analysis, the
analysis code included, were written to `task2f_code_hashes.txt` at 03:49:08 UTC; the only launch on the
evaluation seeds followed at 03:49:14 UTC; the parts were written between 03:49:29 and 06:52:31 UTC; the
check of the code hashes passed, 26 of 26, at 06:53:11 UTC; and the analysis was run once, at 06:54:13
UTC. The registry received no change before the run, only the result paragraph at 06:59:18 UTC and the
paragraph of the verification at 07:26:18 UTC, each with a hash.

Rule 2F-a (primary), per regime: MT-FW+SC minus BOMP-PK+SC at $\le0.25$, pooled over the three cells, has
the lower bound of its 95% interval above zero and a point estimate of at least $0.03$, in both seeds;
strict reading: Holm's correction over the four tests of regime and seed, one-sided at 0.025. Rule 2F-b
(mechanism): MT-FW+SC minus MT-FW-NOCAP+SC has a lower bound above zero in both seeds, pooled over
regimes, and $|$MT-FW-NOCAP+SC minus BOMP-PK+SC$|\le0.02$ on the point estimate in both seeds. Rule 2F-c
(no harm in standard windows): TOST equivalence with margin $0.05$ of MT-FW+SC and BOMP-PK+SC at
$\le0.25$, pooled over cells and regimes, in both seeds. Interpretation: 2F-a and 2F-b pass, the cap is
a confirmed contribution for dense sequences; 2F-a passes and 2F-b fails, the gain is real but not due to
the cap; 2F-a fails, the pilot effect was an artefact of selection. With the pilot standard deviations,
0.0383 with weak clutter and 0.0353 with strong clutter, 96 networks per cell give a power of about 0.95
per seed for an effect of 0.035, half the pilot estimate.

**Design.** Each seed is an independent replication. Networks $j=0$ to 95 per cell, regime and seed form
the dense arm, ten windows each drawn from a stream not used before, with three to five events at
epicentral separations of at least 5 km within 9 km of a centre and origin times in [5.2, 7.8] s at least
0.5 s apart, and at most seven atoms per window; networks $j=96$ to 111, sixteen per cell and regime, form
the equivalence arm, with the standard windows of Task 2A (zero to two events) and at most five atoms. Per
seed the dense arm has 5760 windows with 23013 and 22996 events and the standard arm 960 windows with 766
and 749 events. The three cells, the two clutter regimes, the station corrections from twenty catalogued
events per network, the grid, the origin times, the tolerance, the exclusion and the hit rule are those of
Tasks 2A to 2E. The methods are MT-FW+SC with the cap at norm 6 in noise units per station and phase, the
value fixed in Task 1x (Y.11) and not tuned here; BOMP-PK+SC; BOMP-PK-CAP+SC, BOMP-PK with the same cap in
its selection statistic only; and MT-FW-NOCAP+SC. The exact FROC uses one threshold per method over the
pooled windows, and the hierarchical bootstrap has 2000 replicates, cells as strata, networks and then
windows, the same resamples for every method and the operating point re-chosen in every replicate; the
one-sided $p$-value is $(1+\#\{d^\ast\le0\})/2001$, whose floor is 0.0005. The controls on seed 999
(`task2f_control.txt`) passed: the runner reproduced bit for bit the stored records of the 360 dense pilot
windows for all five pilot methods and the standard-window records of the pilot, the uncapped MT-FW and
BOMP-PK with sequential deflation produced identical hit sequences in 360 of 360 windows, and the fast
FROC equalled the reference implementation to 0.0.

**Table Y.53. Event recall at realized false events per window $\le0.1$ / $\le0.25$, pooled over the three
cells, in the dense and standard arms of Task 2F.** Script `analyze_2f.py`; output `analyze_2f_out.txt`,
first table.

| Arm | Regime | Seed | windows | events | MT-FW+SC | BOMP-PK+SC | BOMP-PK-CAP+SC | MT-FW-NOCAP+SC |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| dense | weak | 4242 | 2880 | 11429 | 0.187/0.289 | 0.102/0.217 | 0.174/0.274 | 0.104/0.215 |
| dense | weak | 8888 | 2880 | 11444 | 0.183/0.287 | 0.094/0.213 | 0.183/0.283 | 0.092/0.209 |
| dense | strong | 4242 | 2880 | 11584 | 0.135/0.220 | 0.080/0.169 | 0.136/0.218 | 0.079/0.169 |
| dense | strong | 8888 | 2880 | 11552 | 0.130/0.220 | 0.070/0.178 | 0.127/0.215 | 0.072/0.177 |
| standard | weak | 4242 | 480 | 381 | 0.428/0.554 | 0.346/0.522 | 0.441/0.546 | 0.346/0.525 |
| standard | weak | 8888 | 480 | 373 | 0.445/0.552 | 0.359/0.555 | 0.429/0.560 | 0.357/0.555 |
| standard | strong | 4242 | 480 | 385 | 0.299/0.379 | 0.291/0.405 | 0.262/0.374 | 0.291/0.405 |
| standard | strong | 8888 | 480 | 376 | 0.255/0.351 | 0.242/0.383 | 0.253/0.354 | 0.239/0.383 |

**Table Y.54. Rules 2F-a, 2F-b and 2F-c, point [95% hierarchical interval].** Same output.

| Rule | Contrast | Regime | Seed 4242 | Seed 8888 | Holm-adjusted $p$ | Passes |
|:--|:--|:--|:--|:--|:--|:--|
| 2F-a | MT-FW+SC minus BOMP-PK+SC, $\le0.25$, dense | weak | $+0.072$ [$+0.048$, $+0.092$] | $+0.074$ [$+0.052$, $+0.100$] | 0.0020, 0.0020 | yes, also under Holm |
| 2F-a | same | strong | $+0.051$ [$+0.032$, $+0.070$] | $+0.043$ [$+0.022$, $+0.066$] | 0.0020, 0.0020 | yes, also under Holm |
| 2F-a, secondary cap | same at $\le0.1$ | weak | $+0.085$ [$+0.063$, $+0.109$] | $+0.089$ [$+0.066$, $+0.114$] | | descriptive |
| 2F-a, secondary cap | same at $\le0.1$ | strong | $+0.055$ [$+0.038$, $+0.076$] | $+0.060$ [$+0.040$, $+0.077$] | | descriptive |
| 2F-b | MT-FW+SC minus MT-FW-NOCAP+SC, $\le0.25$, dense | pooled | $+0.060$ [$+0.046$, $+0.076$] | $+0.059$ [$+0.044$, $+0.075$] | | yes |
| 2F-b | MT-FW-NOCAP+SC minus BOMP-PK+SC, $\le0.25$, dense | pooled | $-0.001$ [$-0.005$, $+0.002$] | $-0.002$ [$-0.006$, $+0.003$] | | yes |
| 2F-c | MT-FW+SC minus BOMP-PK+SC, $\le0.25$, standard, 90% interval | pooled | $-0.008$ [$-0.039$, $+0.034$] | $-0.013$ [$-0.048$, $+0.027$] | | yes; $p$ against the lower margin 0.0170 and 0.0405 |

**Rules 2F-a, 2F-b and 2F-c pass.** In both regimes and both seeds MT-FW+SC detects more events in dense
windows than BOMP-PK+SC at equal realized false-event rate, by $+0.043$ to $+0.074$ at $\le0.25$ and
$+0.055$ to $+0.089$ at $\le0.1$; every one-sided $p$-value is at the bootstrap floor and Holm's correction
rejects all four tests (Table Y.54). The confirmatory estimates are close to the pilot point estimates on
fresh seeds with sixteen times as many networks, so the pilot effect was not an artefact of selection.
Removing the cap brings MT-FW down to BOMP-PK within $\pm0.002$, and the cap alone accounts for the
difference, so 2F-b passes. By the registered interpretation the per-station cap is a confirmed
contribution for dense sequences in both clutter regimes of this synthetic design. Per cell the gain is
positive in 12 of 12 combinations at both caps, from $+0.023$ [$-0.009$, $+0.051$] with four stations
and strong clutter in seed 8888, the only per-cell interval at $\le0.25$ that contains zero, to $+0.110$
with six stations and weak clutter in seed 4242 (Table Y.55). Put in the selection statistic of BOMP-PK,
the cap recovers about 85 to 90 per cent of the gap: BOMP-PK-CAP+SC minus BOMP-PK+SC is $+0.052$ and
$+0.051$ pooled, and MT-FW+SC minus BOMP-PK-CAP+SC $+0.008$ [$-0.001$, $+0.017$] and $+0.006$ [$-0.002$,
$+0.014$]. MT-FW+SC costs 0.105 to 0.118 s per dense window against 0.640 to 0.674 s for BOMP-PK+SC,
about six times less. Rule 2F-c passes, narrowly in seed 8888, where the $p$-value against the lower
margin is 0.0405 and the 95% interval reaches $-0.055$, and with strong clutter the standard-arm point
estimates are negative, $-0.026$ and $-0.032$, and as low as $-0.067$ on the four-station cell, every
interval containing zero: a cost of about 0.03 on isolated events in strong clutter is not excluded by an
arm of sixteen networks per cell. The gain is modest in absolute terms. Even with the cap, recall at
$\le0.25$ in dense windows is 0.22 to 0.29, and most overlapping events are missed by both detectors.

**Table Y.55. Descriptive. Dense arm per cell: recall at $\le0.1$ / $\le0.25$, and differences at
$\le0.25$, point [95% hierarchical interval].** Same output, per-cell table.

| Regime | Seed | Cell | events | MT-FW+SC | BOMP-PK+SC | MT-FW+SC minus BOMP-PK+SC, $\le0.25$ | MT-FW+SC minus MT-FW-NOCAP+SC, $\le0.25$ | BOMP-PK-CAP+SC minus BOMP-PK+SC, $\le0.25$ |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| weak | 4242 | 4 stations | 3826 | 0.114/0.203 | 0.058/0.139 | +0.064 [+0.030, +0.095] | +0.064 [+0.030, +0.094] | +0.061 [+0.028, +0.094] |
| weak | 4242 | 6 stations | 3798 | 0.188/0.302 | 0.083/0.192 | +0.110 [+0.067, +0.154] | +0.112 [+0.073, +0.151] | +0.092 [+0.050, +0.134] |
| weak | 4242 | one-sided, 8 st. | 3805 | 0.262/0.372 | 0.166/0.317 | +0.054 [+0.015, +0.107] | +0.068 [+0.021, +0.109] | +0.040 [-0.003, +0.091] |
| weak | 8888 | 4 stations | 3824 | 0.127/0.216 | 0.061/0.158 | +0.058 [+0.016, +0.093] | +0.061 [+0.024, +0.098] | +0.051 [+0.011, +0.084] |
| weak | 8888 | 6 stations | 3837 | 0.193/0.291 | 0.077/0.182 | +0.109 [+0.066, +0.162] | +0.109 [+0.069, +0.166] | +0.097 [+0.051, +0.150] |
| weak | 8888 | one-sided, 8 st. | 3783 | 0.251/0.370 | 0.145/0.294 | +0.076 [+0.031, +0.127] | +0.076 [+0.033, +0.123] | +0.074 [+0.031, +0.126] |
| strong | 4242 | 4 stations | 3846 | 0.071/0.134 | 0.044/0.107 | +0.028 [+0.001, +0.058] | +0.029 [+0.002, +0.059] | +0.028 [+0.001, +0.056] |
| strong | 4242 | 6 stations | 3842 | 0.150/0.240 | 0.079/0.173 | +0.067 [+0.032, +0.096] | +0.069 [+0.036, +0.101] | +0.069 [+0.033, +0.101] |
| strong | 4242 | one-sided, 8 st. | 3896 | 0.186/0.289 | 0.096/0.216 | +0.073 [+0.029, +0.107] | +0.074 [+0.029, +0.108] | +0.067 [+0.030, +0.099] |
| strong | 8888 | 4 stations | 3881 | 0.071/0.134 | 0.046/0.112 | +0.023 [-0.009, +0.051] | +0.023 [-0.006, +0.049] | +0.013 [-0.011, +0.044] |
| strong | 8888 | 6 stations | 3826 | 0.139/0.245 | 0.078/0.187 | +0.058 [+0.023, +0.101] | +0.064 [+0.026, +0.104] | +0.059 [+0.025, +0.100] |
| strong | 8888 | one-sided, 8 st. | 3845 | 0.173/0.282 | 0.086/0.222 | +0.061 [+0.010, +0.107] | +0.057 [+0.009, +0.100] | +0.046 [+0.004, +0.095] |

**Independent verification and its narrower scope.** An adversarial verification looked for a false
positive (`verif_2f/task2f_verification.txt`), with no change to any file of the task or of the registry.
It found none. The code hashes still check, 26 of 26; each of the eight hashes then recorded is the hash of
a prefix of the current registry, and the registration hash ends exactly at the end of the section "Task
2F"; no confirmatory part predates the freeze; no file or record of seeds 4242 or 8888 existed before the
task; three parts regenerated from scratch in a fresh process, with BOMP-PK run before MT-FW, were
identical in 30 of 30 windows; and an independent implementation of MT-FW reproduced the frozen one bit
for bit in 216 of 216 checks. Both methods receive the same data, grid, origin times, corrections,
number of atoms, exclusion, peak window and hit rule, and the truth enters only the labelling. Counting
duplicate atoms on an event already hit as false events changes the four 2F-a estimates to $+0.068$,
$+0.042$, $+0.075$ and $+0.043$; truncating the records to five, three or two atoms keeps the gain; and
the gain is positive for three, four and five events separately. Own FROC and bootstrap code reproduced
the estimates and intervals, for example $+0.0722$ [$+0.0491$, $+0.0924$] for weak clutter in seed 4242
against the registered [$+0.048$, $+0.092$], and the rules. The verification narrowed what the result
means in four ways.

First, the registered hit rule ignores origin time, and in dense windows all events lie within 9 km, so
it is loose: about half of the location-only hits have an origin error above 0.5 s. When the origin
must also match, the gain survives and shrinks (Table Y.56): with 5 km and 0.5 s it is $+0.024$ to
$+0.045$, and absolute recall falls to 0.09 to 0.13 for MT-FW and 0.07 to 0.08 for BOMP-PK. The rule was
registered and applies equally to both methods, so this is not a false positive, but the absolute gain of
$+0.04$ to $+0.07$ is partly inflated by it.

**Table Y.56. Verification of Task 2F. MT-FW+SC minus BOMP-PK+SC at $\le0.25$ under stricter hit rules,
dense arm, pooled over cells.** Output `verif_2f/task2f_verification.txt`, section 2e.

| Hit rule | weak, 4242 | strong, 4242 | weak, 8888 | strong, 8888 |
|:--|:--|:--|:--|:--|
| registered: epicentre within 5 km | $+0.072$ | $+0.051$ | $+0.074$ | $+0.043$ |
| within 5 km and origin within 1.0 s | $+0.063$ | $+0.052$ | $+0.084$ | $+0.041$ |
| within 5 km and origin within 0.5 s | $+0.045$ | $+0.024$ | $+0.043$ | $+0.027$ |
| within 3 km and origin within 0.5 s | $+0.035$ | $+0.013$ | $+0.027$ | $+0.020$ |

Second, the gain belongs to the level of the cap. On the 360 dense pilot windows of seed 999
(exploratory, own bootstrap), it has a single peak at the fixed value and reverses below it (Table
Y.57). The level 6 was fixed in Task 1x on a different generator and on standard windows, as the best of
6, 12, 24 and no cap, with no value below 6 tried; it was not tuned on dense windows or on the evaluation
seeds, so there is no leakage, but the confirmed gain is that of this level, not of a per-station cap at
an arbitrary level, and the data do not say how a user would choose the level without calibration. At 5
or 8 the pilot gain would still clear the registered bar of 0.03.

**Table Y.57. Verification of Task 2F, exploratory. MT-FW minus BOMP-PK+SC on the 360 dense pilot windows
of seed 999 as a function of the cap level, both regimes pooled.** Output
`verif_2f/task2f_verification.txt`, section 5a.

| Cap level CAPV (noise standard deviations) | MT-FW minus BOMP-PK+SC, $\le0.1$ | MT-FW minus BOMP-PK+SC, $\le0.25$ |
|:--|:--|:--|
| 2 | -0.064 | -0.122 |
| 3 | +0.004 | -0.056 |
| 4 | +0.040 | +0.017 |
| 5 | +0.058 | +0.047 |
| 6 (the fixed value) | +0.074 | +0.070 |
| 8 | +0.058 | +0.046 |
| 12 | +0.033 | +0.010 |
| 24 | -0.003 | -0.008 |
| no cap | -0.003 | +0.001 |

Third, the mechanism is bounded influence under amplitude model error, not robustness to clutter. With the
clutter removed from the same pilot windows the gain rises, to $+0.137$ [$+0.058$, $+0.183$] at $\le0.25$
against $+0.070$ [$+0.016$, $+0.131$] with clutter, so clipping clutter spikes is not the cause, and
clutter partly defeats the cap, which is why the gain is smaller with strong clutter. With the clutter and
the lognormal amplitude factors of the generator both removed (one factor per station and phase with
standard deviation 0.5 in the logarithm, redrawn in every window and in every catalogued event, so that
they are amplitude perturbations of each window and not station terms),
BOMP-PK+SC rises from 0.180 with the clutter alone removed to 0.296 at $\le0.25$ and the gain of the cap falls to $+0.023$ [$-0.040$,
$+0.069$]; caps of 12, 24 or none then equal BOMP-PK. Removing the factors therefore reduces the gain, and
the size of the gain depends on the spread of the amplitude perturbation that the generator assumes. The
cap is nevertheless not a correction of those factors: in the exploration of Y.25 it gains $+0.043$ at a
spread of 0.25, with an interval containing zero, where a statistic given the true factors gains nothing,
and at every spread tried it adds 0.03 to 0.09 to that statistic. What the cap does is bound the influence of any
single station and phase block whose amplitude the radiation model gets wrong, whatever the origin of
the error.

Fourth, the improvement is mainly in the reliability of the score across windows, not in the separation
of events within a window. Within a window the gain is small: the number of distinct events hit among the
first atoms rises by $+0.001$ to $+0.012$, and the hit rate of the first atom by $+0.02$ to $+0.05$.
Across windows the change is larger: among the top 20 per cent of first-atom scores MT-FW is right 0.90
to 0.96 of the time and BOMP-PK 0.79 to 0.84, and the uncapped first-atom scores spread by a factor of 5
between their 10th and 90th percentiles against 2 with the cap. Much of the gain at 0.25 false events per
window is therefore a better calibration of the statistic across windows under a single global
threshold, which is legitimate for the registered endpoint, a single-threshold detector, and is not what
the phrase "resolving overlapping events" would suggest. Selecting with the cap and scoring without it,
and the converse, each give part of the gain ($+0.039$ and $+0.025$ against $+0.070$ on the pilot), and the
two parts combine more than additively.

**Novelty.** Clipping, normalizing or sign-stacking the contribution of each station is established
practice in robust stacking and back-projection, the family reviewed in [136], and the registration
declares it. What Task 2F contributes is empirical and specific: inside a moment-tensor subspace detector
with deflation, clipping each station and phase at a fixed level bounds the influence of each block
under amplitude model error and raises recall in dense sequences under a global threshold, by a
modest, preregistered and replicated amount, while joint refitting and the cross term of P and S do not
matter. It is a measured design result about a bounded-influence statistic, not a new algorithm.

**What Task 2F does not test.** Real data; a cap level chosen without the value of Task 1x; amplitude
errors other than the lognormal factors of the generator, which are redrawn in every window, among them
site terms that persist at a station (Y.25, exploratory); windows with more than five events or
events farther apart than 9 km; and no harm on isolated events in strong clutter beyond the margin of
0.05. The mechanism readings of the verification are exploratory and on seed 999.

### Y.24a What Tasks 2E and 2F change

Task 2E closes the question that Task 2D left open. The principled form of the timing tolerance, the
group BLASSO with a time shift per station and phase solved off the grid by sliding Frank-Wolfe, does
not detect more events than block orthogonal matching pursuit given the heuristic peak search: the gain is
$+0.000$ to $+0.039$ at $\le0.25$, every interval contains zero, the group penalty is inert at the
calibrated value and the independent verification found neither an implementation error nor a budget
handicap. G1-M is refuted in this design, as registered. Together with Tasks 2C and 2D this means that no
convex formulation of the thesis, on the grid or off it, with or without the tolerance, is a better
detector than block pursuit with equal tolerance on these synthetics. The group BLASSO keeps its role as
the framework in which the detection statistic is an optimality certificate and as the object of the
theory of Appendix AA; it is not the method that delivers the performance.

Task 2F confirms a smaller claim, in a preregistered test on fresh seeds with the code frozen and the
registration hash communicated outside the container before the run: a robust, amplitude-clipped
moment-tensor subspace statistic, that of MT-FW with its per-station cap, raises recall in
dense sequences of three to five events over block pursuit and over the same statistic without the cap,
which are equivalent here, by $+0.04$ to $+0.07$ under the registered rule and $+0.024$ to $+0.045$ when
the origin time must match, with no cost beyond the equivalence margin of $0.05$ in standard windows,
where a cost of about 0.03 with strong clutter is not excluded. The verification fixes its mechanism and its limits: the gain is
bounded influence of each station and phase block under amplitude model error and a better calibration
of the score across windows, it belongs to the fixed cap level, and it is not a better separation of
events within a window.
Task 2D found the cap immaterial in standard windows, and Task 2F agrees there (rule 2F-c); the effect is
specific to dense sequences.

These two results point to one research direction, and the exploration of Y.25 says which form of it
the evidence supports. The cap is a crude way of bounding the influence of a station and phase block
whose amplitude the model gets wrong. One candidate generalization models that amplitude as a station, or
site, term: a gain per station and phase estimated from the catalogued events, as the station time
corrections of Y.16 are estimated, and used inside the statistic. Y.25 does not support it. The
amplitude factors of the generator are redrawn in every window, so that it contains no station term to
estimate; in a world where the factors persist, neither the catalogue terms nor the true factors take
the place of the cap, and the cap adds 0.03 to 0.09 to a statistic given the true factors; and at
a small spread of the factors the cap still gains where the true factors gain nothing. The form that the
exploration supports is a robust, bounded-influence moment-tensor subspace statistic whose level is set from the
noise distribution instead of being calibrated: a Huber fit on the residual of each station and phase
block, with its threshold fixed by a quantile rule, matched the cap of 6 in the pilot at spreads of 0.5 or
less, was worse at 0.8 and cost nine times as much. It was the best of 11 candidates on the calibration
seed, so that estimate is subject to selection, and its confirmation is registered as Task 2I, whose
outcome is reported in Y.26. Neither form closes the realistic headroom of Y.22: in the attribution
run of Y.23, removing the site and polarization perturbations from the generator left a headroom of
$+0.153$ for BOMP-PK+SC with weak clutter, and in Y.25 no treatment of amplitudes, the true factors with
the cap included, captures more than about a quarter of it. That headroom, 0.08 to 0.24 below the
realistic ceiling, is not explained by the velocity model or by the source pulse, is left about three
quarters unexplained by the amplitude factors, takes the form of events absorbed by atoms mislocated 5 to
10 km away, and remains the open problem of the detector on these synthetics.

What Y.18a to Y.21a list as not established remains so: the operating points are chosen on the
evaluation windows, no deployable frozen operating point exists, GaMMA with optimistic picks is ahead at
0.1 with weak clutter, the sources are double couples inside the detector's subspace, the one synthetic
generator lies outside the declared regime of Section 1.5, and there is no result on real data.

## Y.25 Site terms versus bounded influence: an exploratory study (C8)

Code: `explore_contrib/c8_lib.py` (generator copy that returns and overrides the amplitude factors, the
two worlds, the catalogue amplitude terms and the detector engine), `c8_run.py`, `c8_analyze.py`,
`c8_validate.py` and `c8_sitedev.py`; outputs `results/explorations/out/c8_analysis_*.txt`,
`c8_validation.txt` and `c8_sitedev.txt`; report `results/explorations/REPORT_C8.txt`.

**Status.** Everything in this section is exploratory and tests no registered hypothesis. It was run on
seed 999 only, on networks of that seed never used before (60 to 79 of each cell for the main designs and
80 to 89 for the sensitivity designs), with random streams used by no earlier design, and no other seed
was generated or read. The protocol is that of Y.23 and Y.24: the grid of 675 nodes, origin times from
5.0 to 8.0 s at 0.1 s, the peak search of $\pm0.25$ s, the 5 km exclusion and hit rule, the exact FROC,
dense windows of three to five events and standard windows of zero to two, the three cells and the two
clutter regimes pooled, and a hierarchical bootstrap over networks and windows with 1000 replicates and
the operating point re-chosen per replicate. All differences are paired on the same windows, and all
values are recall at $\le0.25$ realized false events per window with 95 per cent intervals unless stated.
The run took 2 h 42 min on two cores, 5.3 CPU-hours. The question was whether amplitude terms per
station and phase, estimated from catalogued events and modelled inside the moment-tensor subspace statistic,
can take the place of the per-station cap, as the candidate for its physics-based generalization.

**A fact about the generator.** The function `gen_realistic` draws the lognormal amplitude factors, one
per station and phase with standard deviation 0.5 in the logarithm, at every call, that is, afresh in
every window and in every catalogued event. The validation measures a standard deviation of 0.42 of the
log-factor across the windows of one network. In the generator of Tasks 2A to 2F these factors are
therefore amplitude perturbations of each window, not station terms, and nothing persists from the
catalogue to the windows that a catalogue could learn; Y.16 says as much when it fixes the station delays
per network and leaves the site factors per window. Two worlds were therefore run. In IID the generator is unchanged, the
world of Task 2F. In PERS the same factors are drawn once per network, station and phase, and are shared
by the catalogue and by every window, which is the physically meaningful case of a site term. In both
worlds the P and S factors of a station are independent; the case of one factor for both phases was
planned and not run.

**Methods.** All variants use the engine of MT-FW, with its peak search, running minimum, exclusion and
deflation, and differ in the statistic. With $\mathbf v_b$ the three-component matched-filter vector of
the station and phase block $b$ in units of the noise standard deviation and $Q$ an orthonormal basis of
the radiation operator $G$ of the candidate node, MT-FW computes $\|Q^\top\mathbf v\|^2/2$ after clipping
$\|\mathbf v_b\|$ at 6, and its variants clip at 3, 8, 12 or not at all. ORACLE multiplies the rows of
$G$ by the true factors of the window, the exact whitened likelihood ratio with known gains, and
ORACLE-CAP adds the cap of 6; both are bounds, not detectors. SITE-CAT multiplies the rows of $G$ by
$\exp(u)$ with $u$ the catalogue amplitude term of each station and phase, the median over the twenty
catalogued events of the log-ratio of the observed peak amplitude to the radiation-averaged model
amplitude with the event term removed; SITE-CAT-CAP adds the cap and SITE-CAT-ST uses one term per
station. HUBER-N replaces the least-squares fit by a Huber M-estimate on the norms of the block
residuals, $T=\sum_b\rho_c(\|\mathbf v_b\|)-\min_a\sum_b\rho_c(\|\mathbf v_b-Q_b a\|)$, with $\rho_c$
quadratic below $c$ and linear above, four passes of iteratively reweighted least squares, deflation with
the robust fit and no cap on the data; its threshold is $c=2.7955$, the square root of the 0.95 quantile
of the chi-square distribution with three degrees of freedom, so that a block of pure noise is left
unweighted with probability 0.95. HUBER-S scales that threshold by a robust scale of the residual of the
candidate. GAIN profiles a positive gain per block under a lognormal prior of standard deviation 0.5, the
nominal value of the generator, GAIN-FREE leaves the gains free, SITE-CAT-GAIN centres the prior on the
catalogue terms, and NORM-TR divides each station trace by its peak amplitude, as in back-projection
practice. BOMP-PK+SC is the registered detector and CEIL-GD the realistic ceiling of Y.22. All
definitions, the constant of HUBER-N included, were fixed before any detection run. The catalogue
estimator was chosen on the catalogues of networks 60 to 63 against the true factors, without detection
results (correlation 0.74 against 0.42 for the alternative), which favours SITE-CAT mildly in PERS. The
validation passes: the generator copy equals the original bit for bit, the engine reproduces MT-FW with
and without the cap in 6 of 6 windows each, and the Huber statistic tends to the least-squares one as
$c$ grows. In the runs the catalogue terms correlate 0.72 with the true log-factors (median per network,
interquartile range 0.41 to 0.85), and 76 per cent of the station and phase blocks have the three picks
that a term requires.

**Table Y.58. Exploration C8, dense windows, amplitude spread 0.5: recall at $\le0.25$ and fraction of the
CEIL-GD headroom over BOMP-PK+SC captured, in the world with factors redrawn per window (IID, 600
windows, 2410 events) and in the world with persistent factors (PERS, 480 windows, 1919 events); 95%
hierarchical intervals.** Outputs `results/explorations/out/c8_analysis_DENSE-IID.txt` and
`c8_analysis_DENSE-PERS.txt`; `REPORT_C8.txt`, section 4.

| Method | Recall, IID | Recall, PERS | Headroom fraction, IID | Headroom fraction, PERS |
|:--|:--|:--|:--|:--|
| BOMP-PK+SC | 0.200 [0.167, 0.236] | 0.200 [0.166, 0.245] | 0 | 0 |
| MT-FW, cap 6 | 0.242 [0.215, 0.275] | 0.261 [0.221, 0.289] | 0.14 [0.01, 0.29] | 0.15 [0.02, 0.26] |
| MT-FW, no cap | 0.197 [0.167, 0.234] | 0.201 [0.168, 0.239] | -0.01 [-0.05, 0.03] | 0.00 [-0.04, 0.03] |
| MT-FW, cap 3 | 0.150 [0.117, 0.172] | 0.141 [0.116, 0.165] | not computed | not computed |
| MT-FW, cap 8 | 0.228 [0.198, 0.261] | 0.225 [0.187, 0.268] | not computed | not computed |
| MT-FW, cap 12 | 0.175 [0.151, 0.214] | 0.204 [0.170, 0.237] | not computed | not computed |
| ORACLE (true factors, no cap) | 0.231 [0.190, 0.284] | 0.253 [0.214, 0.309] | 0.10 [-0.02, 0.27] | 0.13 [0.03, 0.27] |
| ORACLE-CAP | 0.287 [0.253, 0.316] | 0.293 [0.255, 0.331] | 0.29 [0.14, 0.46] | 0.23 [0.11, 0.36] |
| SITE-CAT (catalogue terms, no cap) | 0.176 [0.143, 0.210] | 0.221 [0.185, 0.256] | -0.08 [-0.20, 0.02] | 0.05 [-0.04, 0.12] |
| SITE-CAT-CAP | 0.239 [0.211, 0.269] | 0.267 [0.231, 0.310] | 0.13 [0.01, 0.27] | 0.17 [0.05, 0.29] |
| SITE-CAT-ST | 0.204 [0.164, 0.234] | 0.215 [0.170, 0.255] | 0.01 [-0.09, 0.08] | 0.04 [-0.07, 0.11] |
| SITE-CAT-GAIN | not run | 0.223 [0.178, 0.264] | not run | 0.06 [-0.06, 0.15] |
| HUBER-N | 0.264 [0.228, 0.301] | 0.266 [0.234, 0.301] | 0.21 [0.09, 0.39] | 0.16 [0.06, 0.28] |
| HUBER-S | 0.216 [0.184, 0.246] | 0.224 [0.185, 0.272] | 0.05 [-0.06, 0.15] | 0.06 [-0.06, 0.16] |
| GAIN | 0.207 [0.164, 0.250] | 0.228 [0.186, 0.265] | 0.02 [-0.11, 0.15] | 0.07 [-0.05, 0.16] |
| GAIN-FREE | 0.198 [0.147, 0.249] | 0.218 [0.178, 0.256] | -0.01 [-0.19, 0.14] | 0.05 [-0.08, 0.15] |
| NORM-TR | 0.193 [0.168, 0.241] | 0.184 [0.157, 0.219] | -0.02 [-0.14, 0.15] | -0.04 [-0.15, 0.06] |
| CEIL-GD | 0.499 [0.417, 0.599] | 0.605 [0.495, 0.699] | 1 | 1 |

**Table Y.59. Exploration C8, dense windows, amplitude spread 0.5: paired differences in recall at
$\le0.25$ [95% hierarchical interval].** Same outputs as Table Y.58.

| Contrast | IID | PERS |
|:--|:--|:--|
| MT-FW minus BOMP-PK+SC | $+0.042$ [$+0.004$, $+0.083$] | $+0.061$ [$+0.008$, $+0.098$] |
| MT-FW without the cap minus BOMP-PK+SC | $-0.003$ [$-0.015$, $+0.008$] | $+0.001$ [$-0.015$, $+0.012$] |
| ORACLE minus MT-FW without the cap | $+0.034$ [$-0.005$, $+0.079$] | $+0.052$ [$+0.012$, $+0.104$] |
| ORACLE minus MT-FW | $-0.011$ [$-0.056$, $+0.044$] | $-0.008$ [$-0.049$, $+0.055$] |
| ORACLE-CAP minus MT-FW | $+0.045$ [$+0.009$, $+0.074$] | $+0.032$ [$+0.003$, $+0.073$] |
| SITE-CAT minus MT-FW without the cap | $-0.021$ [$-0.054$, $+0.007$] | $+0.020$ [$-0.014$, $+0.048$] |
| SITE-CAT minus MT-FW | $-0.066$ [$-0.108$, $-0.028$] | $-0.040$ [$-0.078$, $+0.007$] |
| SITE-CAT minus ORACLE | $-0.055$ [$-0.105$, $-0.013$] | $-0.032$ [$-0.082$, $-0.002$] |
| SITE-CAT-CAP minus MT-FW | $-0.002$ [$-0.026$, $+0.022$] | $+0.006$ [$-0.015$, $+0.043$] |
| SITE-CAT-ST minus MT-FW | $-0.038$ [$-0.087$, $-0.003$] | $-0.046$ [$-0.091$, $+0.002$] |
| SITE-CAT-GAIN minus MT-FW | not run | $-0.038$ [$-0.085$, $+0.010$] |
| HUBER-N minus MT-FW | $+0.022$ [$-0.014$, $+0.055$] | $+0.005$ [$-0.025$, $+0.048$] |
| HUBER-N minus BOMP-PK+SC | $+0.064$ [$+0.028$, $+0.105$] | $+0.066$ [$+0.024$, $+0.105$] |
| HUBER-N minus ORACLE-CAP | $-0.023$ [$-0.058$, $+0.019$] | $-0.027$ [$-0.069$, $+0.018$] |
| HUBER-S minus MT-FW | $-0.026$ [$-0.065$, $+0.007$] | $-0.037$ [$-0.077$, $+0.011$] |
| GAIN minus MT-FW | $-0.035$ [$-0.080$, $+0.006$] | $-0.033$ [$-0.076$, $+0.017$] |
| GAIN-FREE minus MT-FW | $-0.044$ [$-0.104$, $+0.008$] | $-0.043$ [$-0.083$, $+0.008$] |
| NORM-TR minus MT-FW | $-0.049$ [$-0.074$, $-0.007$] | $-0.077$ [$-0.108$, $-0.032$] |

**Table Y.60. Exploration C8, standard windows, amplitude spread 0.5: difference to BOMP-PK+SC in recall
at $\le0.25$ [95% hierarchical interval], IID (480 windows, 368 events) and PERS (240 windows, 195
events, networks 60 to 69).** Outputs `results/explorations/out/c8_analysis_STD-IID.txt` and
`c8_analysis_STD-PERS.txt`; `REPORT_C8.txt`, section 5.

| Method | IID | PERS |
|:--|:--|:--|
| BOMP-PK+SC, recall | 0.511 [0.432, 0.576] | 0.410 [0.322, 0.527] |
| MT-FW, cap 6 | $-0.041$ [$-0.086$, $+0.029$] | $+0.041$ [$-0.023$, $+0.101$] |
| MT-FW, cap 3 | $-0.253$ [$-0.316$, $-0.171$] | $-0.174$ [$-0.291$, $-0.067$] |
| ORACLE | $+0.016$ [$-0.031$, $+0.082$] | $+0.082$ [$-0.005$, $+0.147$] |
| SITE-CAT | $-0.024$ [$-0.059$, $+0.026$] | $+0.031$ [$-0.034$, $+0.081$] |
| SITE-CAT-CAP | $-0.027$ [$-0.085$, $+0.045$] | $+0.021$ [$-0.054$, $+0.077$] |
| HUBER-N | $-0.014$ [$-0.050$, $+0.042$] | $-0.005$ [$-0.077$, $+0.070$] |
| HUBER-N, 90% interval | [$-0.042$, $+0.035$] | [$-0.063$, $+0.060$] |
| NORM-TR | $-0.133$ [$-0.203$, $-0.063$] | $-0.108$ [$-0.198$, $-0.044$] |

**Table Y.61. Exploration C8, sensitivity to the amplitude spread in the world with factors redrawn per
window: difference to BOMP-PK+SC in recall at $\le0.25$ [95% hierarchical interval], dense windows.** The
columns at 0.25 and 0.8 use networks 80 to 89 with 240 windows each (970 and 960 events); the column at
0.5 is that of Table Y.59, on networks 60 to 79. The recall of BOMP-PK+SC is 0.206, 0.200 and 0.104.
Outputs `results/explorations/out/c8_analysis_DENSE-IID-S25.txt` and `c8_analysis_DENSE-IID-S80.txt`;
`REPORT_C8.txt`, section 6.

| Method | spread 0.25 | spread 0.5 | spread 0.8 |
|:--|:--|:--|:--|
| MT-FW, cap 3 | $-0.060$ [$-0.139$, $+0.022$] | $-0.050$ [$-0.095$, $-0.015$] | $+0.030$ [$-0.042$, $+0.115$] |
| MT-FW, cap 6 | $+0.043$ [$-0.026$, $+0.099$] | $+0.042$ [$+0.004$, $+0.083$] | $+0.119$ [$+0.055$, $+0.207$] |
| MT-FW, cap 8 | $+0.031$ [$-0.020$, $+0.087$] | $+0.028$ [$-0.005$, $+0.068$] | $+0.104$ [$+0.036$, $+0.183$] |
| MT-FW, cap 12 | $-0.001$ [$-0.036$, $+0.031$] | $-0.025$ [$-0.048$, $+0.007$] | $+0.079$ [$+0.007$, $+0.154$] |
| ORACLE | $-0.001$ [$-0.030$, $+0.052$] | $+0.031$ [$-0.007$, $+0.078$] | $+0.076$ [$+0.012$, $+0.148$] |
| ORACLE-CAP | $+0.040$ [$-0.034$, $+0.111$] | $+0.087$ [$+0.043$, $+0.125$] | $+0.167$ [$+0.099$, $+0.249$] |
| SITE-CAT-CAP | $+0.029$ [$-0.046$, $+0.093$] | $+0.039$ [$+0.002$, $+0.079$] | $+0.123$ [$+0.045$, $+0.204$] |
| HUBER-N | $+0.045$ [$-0.011$, $+0.093$] | $+0.064$ [$+0.028$, $+0.105$] | $+0.076$ [$+0.011$, $+0.136$] |
| HUBER-S | $+0.021$ [$-0.032$, $+0.079$] | $+0.016$ [$-0.018$, $+0.048$] | $+0.030$ [$-0.019$, $+0.089$] |
| GAIN | $-0.011$ [$-0.069$, $+0.046$] | $+0.007$ [$-0.029$, $+0.043$] | $+0.061$ [$+0.005$, $+0.116$] |
| GAIN-FREE | $-0.023$ [$-0.081$, $+0.052$] | $-0.002$ [$-0.050$, $+0.045$] | $+0.068$ [$+0.004$, $+0.124$] |
| NORM-TR | $-0.036$ [$-0.091$, $+0.043$] | $-0.007$ [$-0.043$, $+0.045$] | $+0.056$ [$-0.008$, $+0.141$] |
| HUBER-N minus MT-FW, cap 6 | $+0.002$ [$-0.035$, $+0.048$] | $+0.022$ [$-0.014$, $+0.055$] | $-0.043$ [$-0.113$, $-0.001$] |

**Table Y.62. Exploration C8, post hoc: the Huber threshold $c$ of HUBER-N in units of the noise standard
deviation, dense windows of the IID world at amplitude spread 0.5, differences in recall at $\le0.25$ [95%
hierarchical interval].** The thresholds 2.50, 2.7955 and 3.37 are the 0.90, 0.95 and 0.99 quantile rules
on all 600 windows; 1.5, 5 and 8 were run on networks 60 to 69 only, 300 windows, where MT-FW minus
BOMP-PK+SC is $+0.060$ and the threshold 2.7955 gives $+0.091$ and $+0.032$ in the two rows. Outputs
`results/explorations/out/c8_analysis_DENSE-IID+H.txt` and `c8_analysis_DENSE-IID+W.txt`; `REPORT_C8.txt`,
section 7.

| Threshold $c$ | Windows | minus BOMP-PK+SC | minus MT-FW, cap 6 |
|:--|:--|:--|:--|
| 1.5 | 300 | $+0.081$ [$+0.030$, $+0.134$] | $+0.021$ [$-0.029$, $+0.066$] |
| 2.50 | 600 | $+0.059$ [$+0.023$, $+0.101$] | $+0.017$ [$-0.016$, $+0.051$] |
| 2.7955 (fixed by rule before the runs) | 600 | $+0.064$ [$+0.028$, $+0.105$] | $+0.022$ [$-0.014$, $+0.055$] |
| 3.37 | 600 | $+0.064$ [$+0.024$, $+0.094$] | $+0.022$ [$-0.017$, $+0.051$] |
| 5 | 300 | $+0.043$ [$+0.008$, $+0.094$] | $-0.017$ [$-0.062$, $+0.039$] |
| 8 | 300 | $+0.004$ [$-0.018$, $+0.043$] | $-0.056$ [$-0.108$, $+0.004$] |

**Findings.** First, catalogue amplitude terms do not replace the cap. In IID they cannot work by
construction, and SITE-CAT is $0.066$ below MT-FW with an interval excluding zero. In PERS, where a
station term exists, twenty catalogued events on four to eight stations recover a small and uncertain
part of what the true factors give, $+0.020$ [$-0.014$, $+0.048$] over the uncapped statistic against
$+0.052$ [$+0.012$, $+0.104$] for the oracle, and SITE-CAT stays $0.040$ below the cap. With the cap
added, the catalogue terms add nothing measurable: SITE-CAT-CAP minus MT-FW is $-0.002$ and $+0.006$.

Second, even the true factors do not replace the cap. ORACLE without a cap reaches MT-FW within the
interval at a spread of 0.5 ($-0.011$ in IID, $-0.008$ in PERS) and has point estimates below it at 0.25
and 0.8 ($-0.044$ [$-0.087$, $+0.040$] and $-0.043$ [$-0.140$, $+0.034$]). The cap adds 0.03 to 0.09 on top
of the oracle in every design, and the two are complementary: with the cap in place, knowing the true
factors adds $+0.045$ [$+0.009$, $+0.074$] in IID and $+0.032$ [$+0.003$, $+0.073$] in PERS (ORACLE-CAP
minus MT-FW), which is the most that any estimate of the factors could add to the cap. On 32 single events with noise the true factors raise
the fraction of energy captured at the true node only from 0.430 to 0.503: the amplitude factors are one
of several model errors, with the medium, the polarization deviation, the pulse and the residual statics.

Third, the cap is not only a remedy for the amplitude factors. The gain of the oracle over BOMP-PK+SC
grows with the spread, $-0.001$, $+0.031$ and $+0.076$ at 0.25, 0.5 and 0.8, as it should. The gain of
the cap is $+0.043$ already at 0.25, where the oracle has none, although its interval at that spread,
on 240 windows, contains zero. This narrows the third reading of the verification in Y.24, which removed
clutter and amplitude factors together and still found $+0.023$: removing the factors reduces the gain of
the cap, and a part of it remains when they are small. The mechanism that the evidence supports is
therefore bounded influence of each station and phase block under amplitude model error of any origin,
not the correction of a site term.

Fourth, the one principled variant that matches the cap is HUBER-N. At spreads of 0.5 or less it equals
MT-FW with the cap of 6 in three dense samples ($+0.022$, $+0.005$ and $+0.002$, no interval excluding
zero) and exceeds BOMP-PK+SC by $+0.064$ [$+0.028$, $+0.105$] in IID and $+0.066$ [$+0.024$, $+0.105$]
in PERS. In standard windows of the IID world it differs from BOMP-PK+SC by $-0.014$ with a 90 per cent
interval inside $\pm0.05$, where the cap has $-0.041$ [$-0.086$, $+0.029$]; the PERS standard sample is
too small to say, and there the cap has $+0.041$, so the two samples disagree in sign about the cost of
the cap on isolated events and both intervals contain zero. At a spread of 0.8 HUBER-N is worse than the
cap, $-0.043$ [$-0.113$, $-0.001$], while still $+0.076$ above BOMP-PK+SC. It costs 1.07 s per dense
window against 0.11 to 0.12 s for MT-FW, about nine times as much, and 1.5 times BOMP-PK. Its threshold
does not remove the dependence on a level (Table Y.62): the gain is flat from 1.5 to 3.4 noise standard
deviations, roughly halves at 5 and is gone at 8, whereas at a spread of 0.5 the cap helps at 6 and 8,
hurts at 3 and gives $-0.025$ at 12 (Table Y.61). What HUBER-N changes is where the level comes from, the noise distribution rather than a
calibration, and the width of the usable range, about twice that of the cap. The sweep was run only in
dense IID windows at a spread of 0.5; what a small threshold costs in standard windows, and whether it
would help at 0.8, is not known.

Fifth, the other principled variants fail. The scale-adaptive HUBER-S, the gains with or without a prior,
and the trace normalization are $0.026$ to $0.077$ below MT-FW in dense windows and close to BOMP-PK+SC,
and the trace normalization and a cap of 3 destroy recall on isolated events ($-0.133$ and $-0.253$ in
IID standard windows). The convergence of GAIN-FREE was not validated, so its result is a lower bound
for that idea.

Sixth, nothing captures more than about a quarter of the headroom. MT-FW captures 0.14 and 0.15 of the
distance from BOMP-PK+SC to CEIL-GD, HUBER-N 0.21 and 0.16, and the oracle with the cap 0.29 and 0.23,
with wide intervals (Table Y.58). Some three quarters of the measured headroom is untouched by any
treatment of amplitudes, the knowledge of the true factors included, and remains the open problem stated
in Y.22 and Y.24a.

**Selection and other caveats.** HUBER-N was not the variant named in advance as primary among the Huber
fits, which was HUBER-S; it is the best of 11 candidates in the dense IID design, and its estimate there
carries a selection bias. The dense PERS, standard and sensitivity designs are other windows, in part
other networks, run with the same frozen definition, and HUBER-N was added to them after the dense IID
results were seen. The threshold sweep is post hoc. The sensitivity arm was cut from four designs to two,
so nothing is known about catalogue terms at other spreads in PERS. Eleven candidates on six designs
carry no multiplicity correction. No interval for HUBER-N minus MT-FW excludes zero in favour of
HUBER-N, so what the pilot suggests is equivalence with the cap at the nominal spread, not superiority.
The gain of the cap over BOMP-PK+SC in this sample, $+0.042$, is smaller than in the earlier pilots on
other networks of the same seed ($+0.068$ and $+0.070$), with overlapping intervals. The hit rule
ignores the origin time, which Y.24 showed to inflate dense recall for every method alike. Station
amplitude terms are standard in magnitude estimation and in amplitude-based location, and Huber and
other M-estimators are standard in moment-tensor inversion and in robust stacking; what the study adds is
a measured design statement on one synthetic generator, not a method.

**What follows.** Amplitude or site terms per station, estimated from catalogued events and modelled
inside the statistic, are not supported as the generalization of the cap: in the generator of the
registered tasks there is no station term to estimate, and where one exists neither the estimated terms
nor the true ones take the place of bounding the influence of each block. The supported direction is a
robust, bounded-influence moment-tensor subspace statistic whose level is set from the noise distribution
rather than calibrated. It is registered as Task 2I (`PREREGISTRO_1j_1k.txt`, section "Task 2I") on
seeds 5151 and 6262, never used, with 120 networks per cell and regime and the code frozen by hashes:
superiority of HUBER-N over BOMP-PK+SC in dense windows, non-inferiority to MT-FW with the cap of 6
within 0.03, and equivalence with BOMP-PK+SC within 0.05 in standard windows, each in both seeds, planned
on an effect of $+0.04$ because of the selection. Its outcome is reported in Y.26. MT-FW with the cap of
6 remains the reference detector: it is nine times cheaper and at least as good at every spread tried.

## Y.26 A noise-calibrated bounded-influence statistic (Task 2I)

Code: `task2i_run.py` (controls and runs, resumable from one part per network), `task2i_control.py`
(controls on seed 999), `analyze_2i.py` (rules, tables and hierarchical bootstrap) and the frozen methods
of `explore_contrib/c8_lib.py`, `c8_run.py`, `c6_lib.py`, `xc.py` and `task2d_methods.py`; outputs
`analyze_2i_primary_out.txt`, `analyze_2i_out.txt`, `task2i_control.txt`, `task2i_code_hashes.txt`, `task2i_freeze_time.txt`
and `log2i.txt`; design and rules in `PREREGISTRO_1j_1k.txt` (section "Task 2I") with the hashes in
`PREREGISTRO_hashes.txt`; independent verification in `verif_2i/task2i_verification.txt`.

**Status.** This section reports the primary arm of Task 2I, which carries the three registered rules,
its independent verification, and the secondary, descriptive arms of the registration: the Huber
thresholds 2.5003, 3.3682 and 5 and the amplitude spreads 0.25 and 0.8, each on 20 networks per cell and
regime. The secondary arms ran after the primary arm with the same frozen code and were analysed with
the frozen analysis; they carry no rule and were not independently verified. In brief, the Huber
statistic is indistinguishable between levels of 2.5 and 3.4 and weaker at 5; at an amplitude spread of
0.25 the cap gives nothing and the Huber fit little; and at a spread of 0.8 the Huber fit beats block
pursuit but falls behind the cap, so that its non-inferiority to the cap does not extend to that
spread. The three primary verdicts stand as registered.

**Preregistration and origin.** The question came from the exploration of Y.25, where HUBER-N was the
best of 11 candidates on the calibration seed; the registration declares that selection and plans the
sample on an effect of $+0.04$, below the pilot estimate of $+0.064$, because of it. The question is
whether a moment-tensor subspace statistic with a Huber-weighted fit per station and phase, its level fixed by
rule from the noise and with no calibrated cap, (a) exceeds BOMP-PK+SC in dense windows, (b) is not
inferior there to MT-FW+SC with the cap of 6, and (c) leaves the recall of isolated events unchanged.
It was appended to the registry at 18:00:05 UTC on 2 October 2026, before any code was run on seeds 5151
and 6262, which had never been used, and the hash of the registry,
a983beac0c9bdaffab1a176cefd71ba1d855bce9b9e88f24ccebbe4037816f97, was recorded with its clock time and
sent to the user in the working session before any run. That record lies outside the computing
environment and is neither public nor issued by a third party; the verification could confirm from disk
the hash and its position in the registry, not the moment at which it was sent.

Rule 2I-a (superiority, dense windows): HUBER-N+SC minus BOMP-PK+SC at $\le0.25$ has a one-sided 97.5 per
cent lower bound above zero and a point estimate of at least $0.03$. Rule 2I-b (non-inferiority, dense
windows): HUBER-N+SC minus MT-FW+SC has a one-sided 97.5 per cent lower bound above $-0.03$. Rule 2I-c
(equivalence, standard windows): the 90 per cent interval of HUBER-N+SC minus BOMP-PK+SC lies inside
$[-0.05,+0.05]$ (TOST). Each rule must hold in both seeds, pooled over regimes and cells. Interpretation:
the statement that a noise-scaled Huber fit replaces the calibrated cap requires a, b and c; if a holds
and b fails, the Huber fit helps but the cap is better; if a fails, the pilot was an artefact of
selection.

**Design.** Each seed is an independent replication. Networks $j=0$ to 119 of each of the three cells
and two clutter regimes, 120 per cell and regime, carry five dense windows, with three to five events
and at most seven atoms, and five standard windows, with zero to two events and at most five atoms, drawn
from two random streams used by no earlier design. Per seed there are 3600 dense windows with 14477 and
14407 events and 3600 standard windows with 2938 and 2891 events. The generator is that of Tasks 2A to
2F, with the amplitude factors redrawn in every window at a nominal spread of 0.5 (measured mean standard
deviation of the log-factor per window, 0.467 and 0.464), and the station corrections from twenty
catalogued events, the grid, the origin times, the tolerance, the exclusion and the hit rule are
unchanged. The endpoint is recall at $\le0.25$ realized false events per window on the exact FROC, with
$\le0.1$ as secondary, one pooled threshold per method, and a hierarchical bootstrap of 2000 replicates
with regime and cell as strata, networks and then windows, the same resamples for every method and the
operating point re-chosen in every replicate. The bootstrap standard deviations of the three rule
contrasts were 0.0094 and 0.0085, 0.0076 and 0.0078, and 0.0115 and 0.0112, against planned values of
0.0079, 0.0073 and 0.0085.

**Methods.** Four detectors run in sequence on each window, each on its own copy of the matched-filter
field: BOMP-PK+SC; MT-FW+SC with the cap of 6; MT-FW-NOCAP+SC; and HUBER-N+SC. All share the grid,
the origin times, the station corrections, the number of atoms, the exclusion of nodes within 5 km of
an accepted atom, the running-minimum score and the hit function. HUBER-N uses the engine of MT-FW and
changes the statistic. Let $\mathbf v_b\in\mathbb R^3$ be the three-component matched-filter vector of
the station and phase block $b$, in units of the noise standard deviation of the matched filter, and
$Q_b$ the rows of an orthonormal basis $Q$ of the radiation operator of the candidate that belong to
that block. With the Huber function of a norm, $\rho_c(r)=r^2/2$ for $r\le c$ and $\rho_c(r)=cr-c^2/2$
for $r>c$, the statistic is
$$T=\sum_b\rho_c(\|\mathbf v_b\|)-\sum_b\rho_c(\|\mathbf v_b-Q_b\hat a\|),$$
the Huber loss of the data minus the Huber loss of the residual. The coefficient vector $\hat a$ in the
six-dimensional subspace starts at the least-squares fit $Q^\top\mathbf v$ and is refined by four passes
of iteratively reweighted least squares, each with block weights $w_b=\min(1,c/\|\mathbf r_b\|)$ computed
from the current residual $\mathbf r_b=\mathbf v_b-Q_b\hat a$. The function acts on the norm of each
3-vector, not on its components, so that a block is down-weighted as a whole and the statistic does not
depend on the orientation of the sensor axes. The level is $c=2.7955$, the square root of the 0.95
quantile of the chi-square distribution with three degrees of freedom: a block of pure noise falls in
the quadratic range with probability 0.95. It was fixed by that rule before any detection run of the
exploration, it is not tuned, and nothing clips the data. An accepted atom is removed by deflation with
the robust fit, $Q\hat a$, and not with the least-squares fit.

**Freeze and interruptions.** The runner, the analysis and the controls were written after the registry
text and before the freeze, which the registration allows; the engine `c8_lib.py` was last modified on
28 September. The SHA-256 hashes of 27 files, the analysis code included, were written to
`task2i_code_hashes.txt` and the freeze was stamped at 18:14:44 UTC; the first job on the evaluation
seeds was launched at 18:14:45 UTC, and the parts were written between 18:14:58 and 22:04:25 UTC. The
run did not complete in one launch. The process was interrupted at 20:33:25 UTC and resumed at 20:34:09
UTC, and it was interrupted again at 20:46:15 UTC, when the working session was suspended, and resumed
at 21:29:10 UTC. Each resumption used the same command, read the parts already written, one file per
network, and recomputed none of them; the three sessions wrote 1101, 83 and 256 of the 1440 parts. The
code was unchanged throughout (`sha256sum -c`, 27 of 27), and each resumption is recorded with its clock
time in `task2i_freeze_time.txt`. The analysis was run on the complete primary arm, with the sample
size fixed in advance at 120 networks, so that the interruptions gave no opportunity for optional
stopping. The registry received the result paragraph at 22:08:07 UTC and the paragraph of the
verification at 22:16:47 UTC, each with a hash.

**Controls.** The controls on seed 999 (`task2i_control.txt`) passed. The runner reproduced the stored
records of the pilot of Y.25 exactly, in 90 of 90 dense and 72 of 72 standard windows for the four
primary methods, in 60 of 60 windows for the secondary thresholds and in 24 of 24 for each secondary
spread. The general engine with the least-squares statistic gave the hits and positions of the frozen
MT-FW of Task 2F in 60 of 60 windows with and without the cap, with scores equal to a relative
$2.2\times10^{-15}$. The FROC of the analysis equalled that of the exploration and the published pilot
values for the four methods in dense and standard windows, and its bootstrap, applied to the pilot,
returned $+0.064$ [$+0.027$, $+0.104$] for HUBER-N minus BOMP-PK against the published $+0.064$
[$+0.028$, $+0.105$]. The time estimate was 5.27 h on two cores, under the registered budget of 8 h, so
the secondary arms were not dropped.

**Table Y.63. Task 2I, primary arm. Event recall at realized false events per window $\le0.1$ / $\le0.25$,
pooled over the three cells.** Script `analyze_2i.py`; output `analyze_2i_primary_out.txt`, first table.

| Windows | Regime | Seed | windows | events | BOMP-PK+SC | MT-FW+SC, cap 6 | MT-FW-NOCAP+SC | HUBER-N+SC |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| dense | pooled | 5151 | 3600 | 14477 | 0.076/0.177 | 0.157/0.243 | 0.077/0.175 | 0.139/0.235 |
| dense | pooled | 6262 | 3600 | 14407 | 0.092/0.198 | 0.160/0.253 | 0.091/0.195 | 0.151/0.251 |
| dense | weak | 5151 | 1800 | 7191 | 0.078/0.182 | 0.174/0.264 | 0.078/0.182 | 0.149/0.250 |
| dense | weak | 6262 | 1800 | 7232 | 0.101/0.223 | 0.181/0.292 | 0.100/0.224 | 0.166/0.277 |
| dense | strong | 5151 | 1800 | 7286 | 0.077/0.168 | 0.143/0.223 | 0.079/0.167 | 0.130/0.220 |
| dense | strong | 6262 | 1800 | 7175 | 0.084/0.175 | 0.141/0.224 | 0.082/0.173 | 0.137/0.225 |
| standard | pooled | 5151 | 3600 | 2938 | 0.361/0.479 | 0.354/0.469 | 0.360/0.480 | 0.383/0.489 |
| standard | pooled | 6262 | 3600 | 2891 | 0.355/0.488 | 0.367/0.485 | 0.354/0.488 | 0.385/0.496 |
| standard | weak | 5151 | 1800 | 1439 | 0.404/0.559 | 0.443/0.546 | 0.402/0.558 | 0.463/0.566 |
| standard | weak | 6262 | 1800 | 1480 | 0.450/0.603 | 0.475/0.585 | 0.445/0.601 | 0.466/0.571 |
| standard | strong | 5151 | 1800 | 1499 | 0.313/0.426 | 0.309/0.402 | 0.313/0.426 | 0.321/0.424 |
| standard | strong | 6262 | 1800 | 1411 | 0.291/0.422 | 0.300/0.417 | 0.292/0.424 | 0.321/0.425 |

**Table Y.64. Rules 2I-a, 2I-b and 2I-c, and descriptive contrasts on the same windows; point [95%
hierarchical interval] unless stated, pooled over regimes and cells.** Same output, sections of the
rules and of the secondary contrasts.

| Rule | Contrast | Seed 5151 | Seed 6262 | Criterion | Passes |
|:--|:--|:--|:--|:--|:--|
| 2I-a | HUBER-N+SC minus BOMP-PK+SC, dense, $\le0.25$ | $+0.058$ [$+0.039$, $+0.076$] | $+0.053$ [$+0.037$, $+0.070$] | lower bound $>0$ and point $\ge0.03$ | yes, both seeds |
| 2I-b | HUBER-N+SC minus MT-FW+SC, dense, $\le0.25$ | $-0.008$ [$-0.023$, $+0.006$] | $-0.002$ [$-0.018$, $+0.012$] | lower bound $>-0.03$ | yes, both seeds |
| 2I-c | HUBER-N+SC minus BOMP-PK+SC, standard, $\le0.25$, 90% interval | $+0.0099$ [$-0.0100$, $+0.0272$] | $+0.0073$ [$-0.0132$, $+0.0235$] | interval inside $[-0.05,+0.05]$ | yes, both seeds |
| secondary cap | HUBER-N+SC minus BOMP-PK+SC, dense, $\le0.1$ | $+0.063$ [$+0.049$, $+0.077$] | $+0.060$ [$+0.038$, $+0.077$] | | descriptive |
| secondary cap | HUBER-N+SC minus MT-FW+SC, dense, $\le0.1$ | $-0.018$ [$-0.034$, $-0.002$] | $-0.009$ [$-0.031$, $+0.009$] | | descriptive |
| replication of 2F-a | MT-FW+SC minus BOMP-PK+SC, dense, $\le0.25$ | $+0.066$ [$+0.046$, $+0.086$] | $+0.055$ [$+0.038$, $+0.074$] | | descriptive |
| replication of 2F-b | MT-FW-NOCAP+SC minus BOMP-PK+SC, dense, $\le0.25$ | $-0.002$ [$-0.008$, $+0.004$] | $-0.003$ [$-0.008$, $+0.003$] | | descriptive |
| replication of 2F-c | MT-FW+SC minus BOMP-PK+SC, standard, $\le0.25$, 90% interval | $-0.010$ [$-0.031$, $+0.011$] | $-0.003$ [$-0.025$, $+0.016$] | | descriptive |

**Rules 2I-a, 2I-b and 2I-c pass in both seeds.** In dense windows HUBER-N+SC detects more events than
BOMP-PK+SC at equal realized false-event rate, by $+0.058$ and $+0.053$ at $\le0.25$, with one-sided
$p$-values at the bootstrap floor of 0.0005; it is within the non-inferiority margin of MT-FW+SC with
the calibrated cap, the lower bounds being $-0.0234$ and $-0.0184$ against $-0.03$; and in standard
windows it is equivalent to BOMP-PK+SC within $0.05$, with point estimates of $+0.0099$ and $+0.0073$.
By the registered interpretation, the statement that a noise-scaled Huber fit replaces the calibrated
cap is confirmed in this synthetic design. The confirmatory effect is close to the pilot estimate of
$+0.064$, so the selection among 11 candidates produced little regression. The same windows replicate
Task 2F descriptively on two further seeds: MT-FW+SC exceeds BOMP-PK+SC by $+0.082$ and $+0.069$ with
weak clutter and by $+0.054$ and $+0.049$ with strong clutter, MT-FW without the cap equals BOMP-PK+SC
within $0.003$, and in standard windows the capped statistic and BOMP-PK+SC differ by $-0.010$ and
$-0.003$ with 90 per cent intervals inside $\pm0.05$, on 120 networks per cell where Task 2F had 16.
Per regime the gain of HUBER-N over BOMP-PK+SC in dense windows is $+0.068$ and $+0.054$ with weak
clutter and $+0.051$ and $+0.050$ with strong clutter, every interval excluding zero. Per cell it is
positive in 12 of 12 combinations (Table Y.65), smallest for four stations, where two intervals of seed
5151 reach zero, and $+0.054$ to $+0.080$ for six stations and for one-sided networks of eight. At the
stricter false-event cap of $0.1$ HUBER-N is slightly below the calibrated cap in seed 5151, $-0.018$
[$-0.034$, $-0.002$], and not distinguishable from it in seed 6262; this is not a registered rule.

**Table Y.65. Task 2I, descriptive. Dense windows per cell: recall at $\le0.25$ and differences, point
[95% hierarchical interval]; 600 windows per row.** Same output, per-cell table.

| Regime | Seed | Cell | events | BOMP-PK+SC | MT-FW+SC | HUBER-N+SC | HUBER-N+SC minus BOMP-PK+SC | HUBER-N+SC minus MT-FW+SC |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| weak | 5151 | 4 stations | 2362 | 0.137 | 0.210 | 0.176 | +0.039 [-0.000, +0.085] | -0.033 [-0.075, +0.005] |
| weak | 5151 | 6 stations | 2377 | 0.162 | 0.252 | 0.239 | +0.078 [+0.024, +0.117] | -0.013 [-0.058, +0.030] |
| weak | 5151 | one-sided, 8 st. | 2452 | 0.248 | 0.349 | 0.327 | +0.080 [+0.031, +0.129] | -0.022 [-0.072, +0.025] |
| strong | 5151 | 4 stations | 2401 | 0.119 | 0.132 | 0.130 | +0.010 [-0.011, +0.057] | -0.003 [-0.029, +0.029] |
| strong | 5151 | 6 stations | 2434 | 0.158 | 0.226 | 0.212 | +0.054 [+0.012, +0.088] | -0.014 [-0.053, +0.022] |
| strong | 5151 | one-sided, 8 st. | 2451 | 0.243 | 0.296 | 0.299 | +0.056 [+0.012, +0.107] | +0.003 [-0.044, +0.047] |
| weak | 6262 | 4 stations | 2385 | 0.163 | 0.221 | 0.208 | +0.045 [+0.009, +0.079] | -0.013 [-0.052, +0.022] |
| weak | 6262 | 6 stations | 2437 | 0.206 | 0.291 | 0.267 | +0.061 [+0.013, +0.103] | -0.023 [-0.062, +0.017] |
| weak | 6262 | one-sided, 8 st. | 2410 | 0.279 | 0.366 | 0.355 | +0.076 [+0.024, +0.125] | -0.010 [-0.063, +0.038] |
| strong | 6262 | 4 stations | 2376 | 0.112 | 0.138 | 0.141 | +0.029 [+0.004, +0.062] | +0.003 [-0.024, +0.030] |
| strong | 6262 | 6 stations | 2409 | 0.169 | 0.223 | 0.225 | +0.056 [+0.016, +0.097] | +0.001 [-0.039, +0.040] |
| strong | 6262 | one-sided, 8 st. | 2390 | 0.245 | 0.300 | 0.301 | +0.056 [+0.008, +0.113] | +0.001 [-0.039, +0.034] |

**Independent verification.** An adversarial verification looked for a false positive
(`verif_2i/task2i_verification.txt`), with no change to any file of the task or of the registry, and its
verdict is a confirmed positive for the registered claim. On provenance: the code hashes check, 27
files with no failure; the hashes of successive prefixes of the registry reproduce all 11 hashes then
recorded, and the registration hash is that of the prefix ending at the end of the section "Task 2I",
which is therefore byte for byte what was hashed; no part of seeds 5151 or 6262 predates the freeze, the
oldest carrying 18:14:58 UTC; the time stored inside each of the 1440 parts agrees with its file time
within 1.0 s; the stored times show exactly two gaps longer than 60 s, of 82 s and 2599 s, which match
the two recorded interruptions; and no file, directory or analysis output of these seeds exists from
before the task. The plan in the report of the exploration named networks 0 to 79, the registry, which
governs, names 0 to 119, and the run used 0 to 119. On the interruptions: six parts, two from each of
the three sessions and among them the last part before each stop and the first after the second, were
regenerated from scratch, each in a fresh process, and were identical to the stored parts in records,
positions, truth, origin times, amplitude spread, station corrections and metadata; and the 26 parts of
the threshold arm written by the following session contain HUBER-N+SC on the same windows as the primary
arm, 260 windows with no difference. On symmetry: the four methods receive the same window, truth, grid,
origin times, corrections, number of atoms, exclusion and hit function, the truth enters only the
labelling, every record has exactly seven or five atoms for every method, and the stored hit labels equal
a recomputation from the stored positions. The scores of HUBER-N are on another scale than the
least-squares energies, which the registered endpoint handles with one operating point per method chosen
in the same way for all. On the statistics: an own loader, FROC and hierarchical bootstrap reproduced the
point recalls to four decimals and the three rules in both seeds, for example $+0.0582$ [$+0.0402$,
$+0.0753$] against the registered $+0.058$ [$+0.039$, $+0.076$] for 2I-a in seed 5151; rule 2I-b, with
lower bounds of $-0.023$ and $-0.025$ in the two analyses of seed 5151 against the margin of $-0.03$, is
the tightest of the three.

**Scope findings of the verification.** The exploratory checks of the verification, on the stored parts
only, narrow what the result means in six ways.

First, the size of the gain depends on the hit rule. The registered rule asks for an epicentre within 5
km and ignores the origin time; only about 57 per cent of the location hits of any method have an
origin error within 0.5 s, the median error being 0.41 s for all methods. When the origin must also match
within 0.5 s (Table Y.66), HUBER-N+SC minus BOMP-PK+SC falls to $+0.031$ and $+0.019$; both intervals
still exclude zero and the relative gain is 47 and 25 per cent, but the point criterion of $0.03$ would
be met only marginally in seed 5151 and not in seed 6262. With a tolerance of 1.0 s the gain is $+0.054$
and $+0.047$. The registered rule is the location-only one and applies to all methods alike, so the
verdict stands; the gain should be quoted as $+0.05$ to $+0.06$ under the registered rule and $+0.02$ to
$+0.03$ with origin times matched.

**Table Y.66. Task 2I with the hit rule that also requires an origin-time error of at most 0.5 s: recall at
$\le0.1$ / $\le0.25$ and differences at $\le0.25$, point [95% hierarchical interval], pooled over regimes
and cells; scores unchanged.** Output `analyze_2i_primary_out.txt`, section of the origin-time rule, a
registered secondary reading; the verification reproduces the differences with its own bootstrap.

| Windows | Seed | BOMP-PK+SC | MT-FW+SC | HUBER-N+SC | HUBER-N+SC minus BOMP-PK+SC | HUBER-N+SC minus MT-FW+SC | MT-FW+SC minus BOMP-PK+SC |
|:--|:--|:--|:--|:--|:--|:--|:--|
| dense | 5151 | 0.027/0.066 | 0.054/0.101 | 0.048/0.097 | +0.031 [+0.022, +0.040] | -0.004 [-0.015, +0.006] | +0.035 [+0.023, +0.047] |
| dense | 6262 | 0.028/0.076 | 0.052/0.104 | 0.047/0.095 | +0.019 [+0.008, +0.030] | -0.009 [-0.019, +0.003] | +0.028 [+0.015, +0.040] |
| standard | 5151 | 0.158/0.277 | 0.194/0.301 | 0.190/0.307 | +0.030 [+0.008, +0.052] | +0.007 [-0.016, +0.029] | +0.023 [-0.003, +0.048] |
| standard | 6262 | 0.163/0.296 | 0.205/0.313 | 0.201/0.308 | +0.012 [-0.011, +0.034] | -0.004 [-0.027, +0.017] | +0.017 [-0.011, +0.042] |

Second, the equivalence on isolated events holds pooled, as registered, and not regime by regime. In
standard windows with weak clutter in seed 6262 HUBER-N+SC is below BOMP-PK+SC by $-0.032$ [$-0.063$,
$-0.004$], and its 90 per cent interval, [$-0.060$, $-0.007$], leaves the band of $\pm0.05$; the
verification obtains $-0.032$ [$-0.064$, $-0.003$] with its own bootstrap. The other three combinations
of regime and seed show $+0.006$, $-0.003$ and $+0.003$. One of four cells therefore shows a small loss
on isolated events, which the pooled rule does not exclude.

Third, HUBER-N does not find more events; it ranks and selects them better. With every atom accepted
and no threshold, recall does not rise: 0.622 for HUBER-N against 0.623 for BOMP-PK in seed 5151 and 0.616
against 0.628 in seed 6262 with all atoms, and in seed 5151 0.200 against 0.188 with the first atom only
and 0.455 against 0.447 with the first three. What changes is the order. The first atom of a dense window is a true event
more often (0.757 to 0.804 and 0.764 to 0.798), and its score separates windows whose first atom is a
hit from the others better (area under the curve 0.59 to 0.67 and 0.60 to 0.69, with 0.71 for the capped
MT-FW). A rank-matched exchange of scores splits the gain (Table Y.67): giving the atoms of BOMP-PK the
scores of HUBER-N recovers $+0.027$ and $+0.029$, the part due to the reliability of the score across
windows, and giving the atoms of HUBER-N the scores of BOMP-PK recovers $+0.030$ and $+0.018$, the part
due to the choice of atoms. The total divides about half and half in seed 5151 and about 60 to 40 in
favour of score reliability in seed 6262, and the two detectors choose the same first node in only 55
per cent of dense windows. The mechanism is therefore mixed, unlike the reading of Y.24 for the cap, in
which score calibration dominated. In standard windows the two parts have opposite signs, about $+0.010$
from the scores and $-0.012$ and $-0.013$ from the atoms, and cancel into the registered equivalence.

**Table Y.67. Verification of Task 2I, exploratory. Rank-matched exchange of scores in dense windows:
differences to BOMP-PK+SC in recall at $\le0.25$, point [95% interval of the verification, 400
replicates].** Output `verif_2i/task2i_verification.txt`, section 4.4.

| Record | Seed 5151 | Seed 6262 | Reading |
|:--|:--|:--|:--|
| atoms of BOMP-PK with the scores of HUBER-N | $+0.027$ [$+0.016$, $+0.039$] | $+0.029$ [$+0.019$, $+0.039$] | score reliability across windows |
| atoms of HUBER-N with the scores of BOMP-PK | $+0.030$ [$+0.014$, $+0.049$] | $+0.018$ [$+0.004$, $+0.034$] | selection of atoms |
| HUBER-N+SC as run | $+0.058$ | $+0.053$ | registered contrast |

Fourth, at the stricter false-event cap the calibrated cap is slightly better: HUBER-N+SC minus
MT-FW+SC at $\le0.1$ is $-0.018$ [$-0.034$, $-0.002$] in seed 5151 (Table Y.64), with the interval of
seed 6262 containing zero. Non-inferiority within $0.03$ was registered and shown at $\le0.25$ only.

Fifth, the cost. HUBER-N+SC takes 0.927 s per dense window against 0.103 s for MT-FW+SC and 0.602 s for
BOMP-PK+SC in seed 5151, and 0.927, 0.104 and 0.587 s in seed 6262: about nine times the capped
detector and about one and a half times BOMP-PK. In standard windows the figures are 0.657 and 0.658 s
against 0.074 s and 0.425 and 0.416 s. These are seconds per window, with one process per core.

Sixth, absolute recall is low. At $\le0.25$ in dense windows BOMP-PK+SC detects 0.177 and 0.198 of the
events and HUBER-N+SC 0.235 and 0.251; with origin times matched the figures are 0.066 and 0.076 against
0.097 and 0.095. Most overlapping events are missed by every detector, and the gain is weakest, and
unresolved in one seed, for four stations under strong clutter ($+0.010$ [$-0.011$, $+0.057$] and
$+0.029$ [$+0.004$, $+0.062$]). The mean of the per-cell contrasts with an operating point per cell is
$+0.053$ [$+0.036$, $+0.073$] and $+0.054$ [$+0.036$, $+0.071$], so the pooled result is not an artefact
of mixing cells under one threshold.

**Secondary arms.** The registered secondary arms ran after the primary arm under the same freeze. Their
parts were written between 22:04:55 UTC on 2 October and 03:33:41 UTC on 3 October 2026, the frozen
`analyze_2i.py` was run on the complete task (output `analyze_2i_out.txt`), and the registry received
their result paragraph at 05:54:02 UTC on 3 October. Each arm has 20 networks per cell and regime, which
gives 600 dense and 600 standard windows per seed, with 2420 and 2383 dense events and 504 and 501
standard events. The arms are descriptive. No rule was registered for them, their sample is one sixth
of that of the primary arm, so their intervals are wide, and the independent verification covered none
of them beyond the 26 parts of the threshold arm that it compared with the primary arm. The three
primary verdicts do not depend on them.

The threshold arm reruns the Huber statistic on windows of the primary arm at four levels: 2.5003 and
3.3682, the square roots of the 0.90 and 0.99 quantiles of the chi-square with three degrees of freedom
(HUBER-N90 and HUBER-N99), the registered 2.7955, and a fixed level of 5 (HUBER-C5). The registered
level, regenerated in this arm, is identical to the primary arm in 600 of 600 windows of every design
and seed. Between 2.5 and 3.4 noise standard deviations the level does not matter (Table Y.68). In dense
windows HUBER-N90 and HUBER-N99 differ from the registered statistic by $+0.003$ and $+0.005$ in seed
5151 and by $+0.001$ and $+0.006$ in seed 6262, every interval containing zero, and in standard windows
by $+0.000$ and $+0.008$ in both seeds. At a level of 5 the dense gain over BOMP-PK+SC falls to $+0.038$
and $+0.044$, and the statistic is below the capped one by $-0.042$ [$-0.076$, $-0.002$] in seed 5151 and
by $-0.029$ [$-0.070$, $+0.013$] in seed 6262; against the registered level it is lower by $0.013$ and
$0.023$, with intervals that contain zero. The result of the primary arm is therefore not tied to the
value 2.7955. It holds for the quantile rule between the 0.90 and the 0.99 quantile and weakens at 5,
which is the pattern of the post hoc sweep of Y.25 on the calibration seed.

**Table Y.68. Task 2I, threshold arm, descriptive: recall at $\le0.1$ / $\le0.25$ and differences in recall
at $\le0.25$, point [95% hierarchical interval], pooled over regimes and cells; 600 windows per design and
seed.** Output `analyze_2i_out.txt`, section of the Huber threshold.

| Windows | Seed | Statistic | Level $c$ | Recall $\le0.1$ / $\le0.25$ | minus BOMP-PK+SC | minus MT-FW+SC | minus HUBER-N+SC |
|:--|:--|:--|:--|:--|:--|:--|:--|
| dense | 5151 | BOMP-PK+SC |  | 0.068/0.187 |  |  |  |
| dense | 5151 | MT-FW+SC, cap 6 |  | 0.166/0.267 | +0.080 [+0.033, +0.125] |  |  |
| dense | 5151 | HUBER-N90+SC | 2.5003 | 0.142/0.240 | +0.053 [+0.006, +0.095] | -0.026 [-0.068, +0.014] | +0.003 [-0.012, +0.021] |
| dense | 5151 | HUBER-N+SC | 2.7955 | 0.144/0.237 | +0.050 [+0.006, +0.093] | -0.029 [-0.070, +0.013] |  |
| dense | 5151 | HUBER-N99+SC | 3.3682 | 0.140/0.242 | +0.055 [+0.015, +0.096] | -0.025 [-0.061, +0.015] | +0.005 [-0.015, +0.027] |
| dense | 5151 | HUBER-C5+SC | 5.0000 | 0.134/0.224 | +0.038 [+0.004, +0.073] | -0.042 [-0.076, -0.002] | -0.013 [-0.039, +0.017] |
| standard | 5151 | BOMP-PK+SC |  | 0.421/0.538 |  |  |  |
| standard | 5151 | MT-FW+SC, cap 6 |  | 0.411/0.522 | -0.016 [-0.078, +0.055] |  |  |
| standard | 5151 | HUBER-N90+SC | 2.5003 | 0.409/0.522 | -0.016 [-0.066, +0.046] | +0.000 [-0.063, +0.065] | +0.000 [-0.018, +0.023] |
| standard | 5151 | HUBER-N+SC | 2.7955 | 0.409/0.522 | -0.016 [-0.069, +0.041] | +0.000 [-0.066, +0.061] |  |
| standard | 5151 | HUBER-N99+SC | 3.3682 | 0.423/0.530 | -0.008 [-0.055, +0.046] | +0.008 [-0.052, +0.067] | +0.008 [-0.017, +0.036] |
| standard | 5151 | HUBER-C5+SC | 5.0000 | 0.446/0.542 | +0.004 [-0.036, +0.049] | +0.020 [-0.043, +0.077] | +0.020 [-0.025, +0.064] |
| dense | 6262 | BOMP-PK+SC |  | 0.097/0.198 |  |  |  |
| dense | 6262 | MT-FW+SC, cap 6 |  | 0.170/0.272 | +0.074 [+0.024, +0.113] |  |  |
| dense | 6262 | HUBER-N90+SC | 2.5003 | 0.157/0.266 | +0.068 [+0.027, +0.101] | -0.005 [-0.044, +0.032] | +0.001 [-0.008, +0.016] |
| dense | 6262 | HUBER-N+SC | 2.7955 | 0.151/0.265 | +0.067 [+0.025, +0.099] | -0.007 [-0.049, +0.030] |  |
| dense | 6262 | HUBER-N99+SC | 3.3682 | 0.147/0.271 | +0.073 [+0.026, +0.103] | -0.000 [-0.047, +0.034] | +0.006 [-0.012, +0.020] |
| dense | 6262 | HUBER-C5+SC | 5.0000 | 0.148/0.242 | +0.044 [+0.010, +0.074] | -0.029 [-0.070, +0.013] | -0.023 [-0.049, +0.009] |
| standard | 6262 | BOMP-PK+SC |  | 0.349/0.451 |  |  |  |
| standard | 6262 | MT-FW+SC, cap 6 |  | 0.377/0.465 | +0.014 [-0.041, +0.067] |  |  |
| standard | 6262 | HUBER-N90+SC | 2.5003 | 0.365/0.493 | +0.042 [-0.021, +0.098] | +0.028 [-0.029, +0.076] | +0.000 [-0.020, +0.024] |
| standard | 6262 | HUBER-N+SC | 2.7955 | 0.359/0.493 | +0.042 [-0.022, +0.094] | +0.028 [-0.030, +0.071] |  |
| standard | 6262 | HUBER-N99+SC | 3.3682 | 0.355/0.501 | +0.050 [-0.013, +0.096] | +0.036 [-0.019, +0.072] | +0.008 [-0.021, +0.032] |
| standard | 6262 | HUBER-C5+SC | 5.0000 | 0.365/0.467 | +0.016 [-0.018, +0.054] | +0.002 [-0.050, +0.054] | -0.026 [-0.074, +0.033] |

The spread arms redraw the windows of 20 networks per cell and regime with the standard deviation of the
log amplitude factors set to 0.25 and to 0.8 instead of 0.5; the measured values in dense windows are
0.236 and 0.233, and 0.755 and 0.744 (Table Y.69). At a spread of 0.25 the cap gives nothing over
BOMP-PK+SC in dense windows, $+0.005$ [$-0.038$, $+0.041$] and $+0.009$ [$-0.024$, $+0.053$], and the
Huber fit little and not in both seeds: $+0.036$ [$+0.004$, $+0.070$] in seed 5151 and $+0.002$
[$-0.035$, $+0.044$] in seed 6262. At a spread of 0.8 block pursuit falls to a recall of 0.088 and 0.110
at $\le0.25$ while the capped statistic holds 0.243 and 0.251, a gain of $+0.155$ [$+0.096$, $+0.207$]
and $+0.141$ [$+0.095$, $+0.189$]. The Huber fit reaches 0.187 and 0.182: it beats BOMP-PK+SC by $+0.100$
[$+0.054$, $+0.143$] and $+0.072$ [$+0.031$, $+0.115$], but it falls behind the cap by $-0.055$
[$-0.113$, $+0.002$] and $-0.069$ [$-0.121$, $-0.019$]. In seed 6262 that interval lies below zero, and
in both seeds its lower end is far beyond the margin of $-0.03$ of rule 2I-b, so the non-inferiority of
the Huber fit to the cap does not extend to an amplitude spread of 0.8. With 20 networks this arm could
not have shown non-inferiority at the nominal spread either, where the same subsample gives $-0.029$
[$-0.070$, $+0.013$] and $-0.007$ [$-0.049$, $+0.030$]; what it shows at 0.8 is a point estimate two to
ten times larger and, in one seed, an interval that excludes zero. This agrees with the pilot of Y.25,
in which HUBER-N was below the cap by $0.043$ at that spread. In standard windows no contrast between
HUBER-N+SC and BOMP-PK+SC has an interval that excludes zero at either spread: $+0.030$ and $-0.002$ at
0.25 and $+0.016$ and $+0.032$ at 0.8. The gain of bounded influence over block pursuit therefore grows
with the amplitude error of the model, as the mechanism of Y.24 and Y.25 implies, and when that error
is large the hard cap at a calibrated level is the better of the two forms.

**Table Y.69. Task 2I, amplitude-spread arms, descriptive: recall at $\le0.1$ / $\le0.25$ and differences
in recall at $\le0.25$, point [95% hierarchical interval], pooled over regimes and cells; 600 windows per
row. The rows at 0.5 are the same networks in the primary arm.** Output `analyze_2i_out.txt`, section of
the amplitude spread.

| Windows | Seed | Spread | BOMP-PK+SC | MT-FW+SC, cap 6 | HUBER-N+SC | HUBER-N+SC minus BOMP-PK+SC | HUBER-N+SC minus MT-FW+SC | MT-FW+SC minus BOMP-PK+SC |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| dense | 5151 | 0.25 | 0.128/0.252 | 0.155/0.257 | 0.186/0.288 | +0.036 [+0.004, +0.070] | +0.031 [-0.005, +0.074] | +0.005 [-0.038, +0.041] |
| dense | 5151 | 0.5 | 0.068/0.187 | 0.166/0.267 | 0.144/0.237 | +0.050 [+0.006, +0.093] | -0.029 [-0.070, +0.013] | +0.080 [+0.033, +0.125] |
| dense | 5151 | 0.8 | 0.030/0.088 | 0.155/0.243 | 0.083/0.187 | +0.100 [+0.054, +0.143] | -0.055 [-0.113, +0.002] | +0.155 [+0.096, +0.207] |
| standard | 5151 | 0.25 | 0.405/0.528 | 0.403/0.500 | 0.433/0.558 | +0.030 [-0.019, +0.077] | +0.058 [+0.004, +0.102] | -0.028 [-0.068, +0.025] |
| standard | 5151 | 0.5 | 0.421/0.538 | 0.411/0.522 | 0.409/0.522 | -0.016 [-0.069, +0.041] | +0.000 [-0.066, +0.061] | -0.016 [-0.078, +0.055] |
| standard | 5151 | 0.8 | 0.246/0.450 | 0.345/0.460 | 0.313/0.466 | +0.016 [-0.042, +0.091] | +0.006 [-0.044, +0.068] | +0.010 [-0.062, +0.081] |
| dense | 6262 | 0.25 | 0.153/0.265 | 0.189/0.274 | 0.162/0.267 | +0.002 [-0.035, +0.044] | -0.007 [-0.044, +0.021] | +0.009 [-0.024, +0.053] |
| dense | 6262 | 0.5 | 0.097/0.198 | 0.170/0.272 | 0.151/0.265 | +0.067 [+0.025, +0.099] | -0.007 [-0.049, +0.030] | +0.074 [+0.024, +0.113] |
| dense | 6262 | 0.8 | 0.030/0.110 | 0.161/0.251 | 0.084/0.182 | +0.072 [+0.031, +0.115] | -0.069 [-0.121, -0.019] | +0.141 [+0.095, +0.189] |
| standard | 6262 | 0.25 | 0.379/0.483 | 0.339/0.457 | 0.389/0.481 | -0.002 [-0.044, +0.036] | +0.024 [-0.028, +0.072] | -0.026 [-0.074, +0.021] |
| standard | 6262 | 0.5 | 0.349/0.451 | 0.377/0.465 | 0.359/0.493 | +0.042 [-0.022, +0.094] | +0.028 [-0.030, +0.071] | +0.014 [-0.041, +0.067] |
| standard | 6262 | 0.8 | 0.259/0.425 | 0.365/0.477 | 0.317/0.457 | +0.032 [-0.045, +0.109] | -0.020 [-0.086, +0.033] | +0.052 [-0.033, +0.149] |

The registered secondary readings by regime and with the strict hit rule are those of Tables Y.64 to
Y.66 and of the scope findings above, and the final output leaves them as they were. Per regime, with the
origin time also required within 0.5 s, HUBER-N+SC minus BOMP-PK+SC in dense windows is $+0.040$
[$+0.027$, $+0.054$] and $+0.020$ [$+0.003$, $+0.040$] with weak clutter and $+0.022$ [$+0.009$,
$+0.035$] and $+0.017$ [$+0.004$, $+0.030$] with strong clutter.

**Novelty.** Huber-weighted fitting is standard robust statistics, and M-estimators are standard in
moment-tensor inversion and in robust stacking; the registration and Y.25 declare it. What Task 2I
contributes is a measured design result: inside a moment-tensor subspace detector with deflation, a Huber
fit on the residual norm of each station and phase block, with its level fixed by a quantile rule of the
noise distribution, recovers the gain of the calibrated cap over block pursuit in dense windows, within
a non-inferiority margin of $0.03$, without a constant tuned on detection results and without a
measurable cost on isolated events when pooled. It is not an algorithm.

**What Task 2I does not test.** Real data; a generator family other than `gen_realistic`, of which the
two seeds are independent draws and not independent evidence about real data; amplitude spreads other
than 0.25, 0.5 and 0.8 and Huber levels other than the four of the threshold arm, of which only the
registered level at the nominal spread carries a rule and a verification; persistent site terms;
windows with more than five events; a held-out operating point, since the thresholds are chosen on the
evaluation windows, per method and in the same way for all; and non-inferiority to the cap at $\le0.1$.

### Y.26a What Task 2I changes

Task 2I gives the bounded-influence hypothesis G1-R a second preregistered test, on seeds that the
first never touched, with the registration hash sent to the user before any run, the code frozen by
hashes and an independent verification that regenerated parts across the two interruptions of the run.
Task 2F showed that bounding the contribution of each station and phase raises recall in dense
sequences over block pursuit with equal timing tolerance, with a cap whose level, 6, had been fixed in
an earlier task and at which the gain peaks. Task 2I shows the same gain, $+0.058$ and $+0.053$ against
$+0.066$ and $+0.055$ for the cap on the same windows, with a Huber level that no detection result
selected: the square root of a chi-square quantile of the noise. Stated at the scope that the two tasks
support, the result is this. In dense overlapping sequences of three to five events on this synthetic
generator, a bounded-influence moment-tensor subspace statistic improves event recall at a fixed
realized false-event rate over block orthogonal matching pursuit with equal timing tolerance, by about
$0.04$ to $0.07$ under the registered hit rule and about $0.02$ to $0.045$ when origin times must match
within 0.5 s, and it does not harm isolated events when regimes are pooled. Two forms deliver it, a cap
at a calibrated level and a Huber fit at a level fixed by rule, and they are not distinguishable at
$\le0.25$.

What Task 2I does is replace the calibrated cap by a level fixed by rule from the noise distribution in the registered primary arm. What it does not remove is
a level: the post hoc sweep of Y.25 on the calibration seed found the gain of HUBER-N flat between 1.5
and 3.4 noise standard deviations and gone at 8, and the registered threshold arm finds the same on the
evaluation seeds, no difference between 2.5 and 3.4 and a smaller gain at 5. The Huber form is not a
better detector than the cap. It costs about nine times as much, it is slightly below the cap at
$\le0.1$ in one seed, and at an amplitude spread of 0.8 it falls behind the cap by $0.055$ and $0.069$
in the spread arm, as it did in the pilot, while at a spread of 0.25 neither form gains reliably over
block pursuit; MT-FW with the cap of 6 stays the cheaper reference, and the Huber form is the one whose level a user can set without
calibration data.

The limits are those of the evidence. All of it is synthetic and from one generator family; the effect
is modest and absolute recall in dense windows stays near 0.24; under the strict hit rule the gain
shrinks to $+0.031$ and $+0.019$; equivalence on isolated events fails in one regime of one seed; the
gain is a better ranking and selection of atoms, not more events found; and bounded-influence fitting is
standard, so that the contribution is the measured design result and not a method. The open problem is
unchanged. The headroom of 0.08 to 0.24 below the realistic ceiling of Y.22 is not closed by either
form: in Y.25 the cap captured 0.14 and 0.15 of it, HUBER-N 0.21 and 0.16 and the true amplitude factors
with the cap 0.29 and 0.23, so that about three quarters of the measured headroom remains unexplained by
any treatment of amplitudes, and Task 2I, whose recalls match those of the pilot, does not change that
figure.

What Y.18a to Y.24a list as not established remains so: no deployable frozen operating point exists,
GaMMA with optimistic picks is ahead at 0.1 with weak clutter, the sources are double couples inside the
detector's subspace, the generator lies outside the declared regime of Section 1.5, and there is no
result on real data.

## Y.27 Evidence that does not use travel times, and the timing budget: two exploratory studies (C9, C10)

Code: `explore_contrib/c9_lib.py` (travel-time-free statistics, vector readers and the detector loops
with relocation and selection hooks), `c9_standalone.py`, `c9_fisher.py`, `c9_run.py`, `c9_extra.py`,
`c9_analyze.py`, `c9_standalone_analyze.py` and `c9_pairs.py`; `explore_contrib/c10_lib.py` (timing-level
rule, penalised peak choice, detector loops with the timing charge and the refinement, hit rules),
`c10_run.py`, `c10_sweep.py`, `c10_extra.py`, `c10_diag.py`, `c10_screloc.py`, `c10_screloc_run.py` and
`c10_analyze.py`; outputs `results/explorations/out/c9_*.txt` and `c10*_analysis_*.txt`, `c10_diag.txt`,
`c10_validation.txt` and `c10_screloc_check.txt`; reports `results/explorations/REPORT_C9.txt` and
`REPORT_C10.txt`; literature scan `explore_contrib/LIT_NO_TRAVELTIME.md`.

**Status.** Everything in this section is exploratory and tests no registered hypothesis. Both studies
ran on seed 999 only, on networks of that seed used by no earlier design (130 to 149 of each cell for
C9 and 160 to 179 for C10, with networks 199 and 195 kept for validation and for the choices made
before the runs), with random streams used by no earlier design; no other seed was generated or read,
the only contact of C10 with seeds 5151 and 6262 being the run of the frozen analysis of Task 2I after
that job had ended. The protocol is that of Y.23 to Y.25: the grid of 675 nodes at 3 km and three depths,
origin times from 5.0 to 8.0 s at 0.1 s, the peak search of $\pm0.25$ s, the 5 km exclusion and hit
rule, the exact FROC, dense windows of three to five events and standard windows of zero to two, the
three cells and the two clutter regimes, station corrections from the registered catalogue of twenty
events, and paired differences on the same windows with 95 per cent intervals from a hierarchical
bootstrap of 1000 replicates, 500 for the oracle worlds of C9 and for the per-regime tables of C10.
Each design has 600 windows, 20 networks by 5 windows by 3 cells by 2 regimes: 2388 events in the dense
and 444 in the standard windows of C9, and 2432 and 471 in those of C10. No multiplicity correction is
applied anywhere. In C9 the statistics, the relocation radius of 12 km and the oracle worlds were fixed
before any detection run; the exponents of the re-ranking, the selection variants SELPOL and PENPOL,
the worlds with 20 and 45 degrees of polarization scatter and the per-cell contrasts were added after
a first look at the dense windows with weak clutter, so that SELPOL is the best of about 50 candidate
records per base and its estimate carries a selection bias. In C10 the variants, the two hit rules and
the oracle were fixed before the runs on networks 160 to 179 and the rule for the timing level was
chosen on validation network 195; the run set was trimmed to the compute budget; and three analyses
are post hoc: the true statics with the catalogue timing level, the first-iteration diagnostic and the
re-estimated station corrections SC-R, whose rule was written after the oracle result had been seen and
was checked against the true statics before its detection run. C10 reports two hit rules, the
registered one, an epicentre within 5 km (LOC), and the one that also requires the origin time within
0.5 s (LOC+T), which no registered task before Task 2L uses as its endpoint. Both reports end with a draft
registration, for SELPOL and for the refinement T3, on the same two unused seeds; neither was
registered or run as drafted, the protocol of Task 2L contains the arms of the T3 draft, and the SELPOL
draft was discarded (Y.28).

**The question.** Y.22 to Y.25 left a headroom of 0.08 to 0.24 below the realistic ceiling CEIL-GD that
no treatment of the grid, the velocity model, the source pulse or the amplitudes closed, and a
recurrent pattern in the misses: about one true event in five is not hit while a false atom sits 5 to
10 km from it (0.204 of the events in the dense windows and 0.212 in the standard windows of C10 for
MT-FW+SC). Two explanations were examined. The first is that the timing-tolerant statistic cannot tell
neighbouring nodes apart and that evidence which does not depend on travel times, the polarization of
P and S, the amplitudes and the S to P ratio, could relocate those atoms (C9). The second is that the
hard window of $\pm0.25$ s throws timing information away and that a Gaussian timing prior inside the
detector would recover it (C10).

**C9: the information in each channel.** The Bayesian Cramer-Rao (Van Trees) bound on the epicentre was
computed for 200 sources per cell with random double couples, a homogeneous straight-ray model, and the
prior of the benchmark, which alone gives 12.25 km (Table Y.70). Arrival times with the tolerance of
$\pm0.25$ s read as a uniform law, a standard deviation of 0.144 s, bound the epicentral error at 0.4
to 0.7 km, and at 0.8 to 1.3 km with a timing error of 0.30 s. Every channel that does not use travel
times is an order of magnitude weaker: P polarization 7.8 to 9.2 km at the nominal scatter of 10
degrees, P and S polarization together 6.0 to 7.5 km, amplitudes 7.8 to 9.2 km, and the S to P ratio
10.5 to 11.8 km, which is no information beyond the prior, since in a medium without attenuation the
ratio does not depend on distance. All of them together give 4.4 to 5.6 km, a factor of 7 to 10 in
error and of 50 to 100 in information behind the times. Adding them to the times lowers the bound by 4
to 5 per cent, from 0.68, 0.43 and 0.69 km to 0.65, 0.41 and 0.66 km, and by 7 to 9 per cent in a world
with no polarization deviation and no amplitude factors. The two perturbations of the generator are not
what limits these channels. Setting the polarization deviation to zero moves the P bound only to 7.1 to
8.8 km, and removing the amplitude factors moves the amplitude bound to 5.3 to 7.3 km. The limit is the
waveform noise: at these distances the P arrival has an amplitude of 1.5 to 1.8 noise standard
deviations, which is 32 to 39 degrees of angular noise per station before any deviation is added. With
that noise removed, a scatter of 10 degrees would give 3.8 to 5.0 km and one of 2 degrees 0.9 to 1.2 km.
Travel-time-free evidence would become comparable to the times, 0.7 to 1.0 km, only at about ten times
the signal-to-noise ratio of the benchmark with a polarization scatter of at most 2 degrees and an
amplitude-factor spread of at most 0.1. A Monte Carlo check of the bound on 300 sources found the grid
maximum-likelihood estimators above it by a factor of 1.2 to 1.4, as a bound requires. The bound is
local and Gaussian, uses homogeneous rays and a site factor of one, and does not count the discrete
polarity information.

**Table Y.70. Exploration C9: Bayesian Cramer-Rao bound on the epicentral error per channel, median over
200 sources per cell of $\sqrt{C_{xx}+C_{yy}}$ in km.** Output `results/explorations/out/c9_fisher.txt`;
`REPORT_C9.txt`, section 6.

| Channel | 4 stations | 6 stations | one-sided, 8 st. |
|:--|:--|:--|:--|
| no data, prior only | 12.25 | 12.25 | 12.25 |
| arrival times, 0.05 s | 0.24 | 0.15 | 0.26 |
| arrival times, 0.144 s, the tolerance of $\pm0.25$ s as a uniform law | 0.68 | 0.43 | 0.69 |
| arrival times, 0.30 s, uncorrected station delays | 1.34 | 0.84 | 1.29 |
| P polarization, scatter 10 degrees (nominal) | 9.23 | 8.53 | 7.79 |
| P polarization, scatter 0 degrees, waveform noise only | 8.77 | 7.90 | 7.14 |
| P polarization, scatter 10 degrees, no waveform noise | 4.99 | 4.07 | 3.82 |
| P polarization, scatter 2 degrees, no waveform noise | 1.15 | 0.88 | 0.85 |
| P and S polarization, scatter 10 degrees | 7.47 | 6.10 | 6.03 |
| amplitudes, factor spread 0.5 | 9.18 | 7.75 | 9.11 |
| amplitudes, factor spread 0, waveform noise only | 7.33 | 5.31 | 6.77 |
| S to P amplitude ratio, factor spread 0.5 | 11.77 | 10.46 | 10.88 |
| all travel-time-free channels, nominal | 5.59 | 4.35 | 4.81 |
| all travel-time-free channels, scatter 0 and spread 0 | 4.01 | 3.10 | 3.57 |
| the same at ten times the signal-to-noise ratio | 0.45 | 0.33 | 0.39 |
| scatter 2 degrees, spread 0.1, ten times the signal-to-noise ratio | 1.00 | 0.71 | 0.77 |
| arrival times at 0.144 s jointly with all travel-time-free channels, nominal | 0.65 | 0.41 | 0.66 |
| ratio of that joint bound to the bound of the times alone | 0.960 | 0.962 | 0.953 |

**C9: travel-time-free locators on single events.** The bound was checked against estimators on 600
single events with weak clutter, each located on a 1 km grid from the P and S peak vectors of every
station (Table Y.71). With the vectors read at the true arrival samples, an oracle reading that
isolates the channel from the problem of finding the arrivals, the median epicentral error is 12.7 to
15.2 km for the direction of P motion, 8.2 to 9.8 km for P along the ray and S across it, 7.7 to 8.4 km
for the projection on the moment-tensor subspace with the times left free, which is the best
travel-time-free locator, and 12.5 to 19.4 km for a 1/r amplitude fit, against 0.5 to 0.7 km for the
arrival times with the true station delays and 11.8 km for guessing the centre of the source region. At
best one event in three falls inside the 5 km hit radius. In a generator with no polarization deviation
and no amplitude factors the best channel reaches 5.4 to 7.5 km. More stations help little, 8.4, 7.7
and 8.1 km for four, six and eight, and the one-sided network of eight is no better than four stations.
Without travel times the arrivals cannot even be identified: a reader that looks for the pair of
energy peaks over the whole window lands within three samples of the true P arrival at 0.19 to 0.20 of
the stations and of the true S arrival at 0.23 to 0.27, and every estimator is then at or beyond the
level of the guess. These medians lie above the bounds of Table Y.70 by the factor of 1.5 to 1.9
expected from the coda, the clutter, the pulse and medium mismatch and the 1 km grid, which the bound
does not contain.

**Table Y.71. Exploration C9: location of 600 single events from the peak vectors read at the true
arrival samples, nominal generator; median epicentral error in km [95% interval over networks] and
fraction of events within 5 km.** Guessing the centre of the source region gives 11.8 km and 0.09. Output
`results/explorations/out/c9_standalone_analysis.txt`; `REPORT_C9.txt`, section 4.

| Locator | 4 stations | 6 stations | one-sided, 8 st. |
|:--|:--|:--|:--|
| back-azimuth intersection of the P vectors | 19.8 [18.4, 22.0]; 0.07 | 18.1 [15.3, 19.8]; 0.07 | 20.5 [18.9, 22.9]; 0.04 |
| direction of the P vector, grid search | 13.9 [11.9, 16.5]; 0.15 | 15.2 [11.7, 17.0]; 0.18 | 12.7 [9.8, 15.3]; 0.16 |
| S transverse to the ray | 15.9 [13.7, 18.3]; 0.10 | 14.3 [13.3, 15.8]; 0.13 | 14.5 [12.8, 16.3]; 0.15 |
| P along the ray and S transverse | 9.8 [7.9, 12.8]; 0.26 | 8.2 [6.9, 10.5]; 0.27 | 8.9 [7.0, 10.4]; 0.27 |
| radiation fit with a free gain per station | 8.9 [7.1, 10.9]; 0.26 | 8.1 [5.6, 10.6]; 0.34 | 8.2 [6.3, 9.7]; 0.32 |
| moment-tensor projection, times free | 8.4 [6.9, 10.4]; 0.29 | 7.7 [5.9, 10.2]; 0.33 | 8.1 [7.1, 10.0]; 0.30 |
| 1/r amplitude fit | 14.2 [12.8, 15.8]; 0.07 | 12.5 [10.9, 14.2]; 0.11 | 19.4 [17.8, 20.9]; 0.04 |
| arrival times, uncorrected delays | 1.3 [1.0, 1.7]; 0.98 | 0.9 [0.7, 1.0]; 1.00 | 1.4 [1.1, 1.7]; 1.00 |
| arrival times, true station delays | 0.6 [0.6, 0.7]; 1.00 | 0.5 [0.5, 0.6]; 1.00 | 0.7 [0.6, 0.8]; 1.00 |

**C9: relocation loses recall, with oracles too.** The atoms of BOMP-PK+SC and of MT-FW+SC were moved
to the maximum of each travel-time-free statistic within 12 km, after detection or inside the loop
before exclusion and deflation, from the atom's own peak vectors or from vectors re-read within
$\pm1.5$ s. Every variant loses recall (Table Y.72): $-0.11$ to $-0.18$ at $\le0.25$ in dense windows
and $-0.23$ to $-0.48$ in standard windows. The oracle variants do not rescue it. With the rows of the
radiation operator multiplied by the true amplitude factors of the window the loss is $-0.116$
[$-0.158$, $-0.081$] for MT-FW in dense windows, and in a world with no polarization deviation and no
amplitude factors the relocation inside the loop loses $-0.138$ [$-0.187$, $-0.089$]. The diagnostics
say why. In dense windows a relocation keeps 0.43 to 0.57 of the hits and converts 0.22 to 0.46 of the
near misses, which is what a random move inside a cluster of 9 km does, since a move of up to 12 km
often lands within 5 km of some event; in standard windows, where a conversion would be a real gain, it
keeps 0.07 to 0.39 of the hits and converts 0.15 to 0.23 of the near misses. The median error of the
hits grows from 2.6 to 2.7 km to 3.1 to 3.2 km. The relocation does not find the event; it scatters the
atom, as the bound predicts for a statistic whose own precision is 5 to 9 km.

**Table Y.72. Exploration C9: relocation of the atoms by travel-time-free statistics, paired difference
to the unmoved base in recall at $\le0.25$, point [95% interval]; 600 windows per design.** Output
`results/explorations/out/c9_analysis_DENSE.txt` and `c9_analysis_STD.txt`; `REPORT_C9.txt`, section 5.2.

| Relocation statistic | dense, base BOMP-PK+SC | dense, base MT-FW+SC | standard, base MT-FW+SC |
|:--|:--|:--|:--|
| moment-tensor projection, own vectors | -0.123 [-0.157, -0.084] | -0.116 [-0.153, -0.083] | -0.255 [-0.336, -0.183] |
| P along the ray and S transverse, own vectors | -0.144 [-0.182, -0.104] | -0.129 [-0.175, -0.093] | -0.277 [-0.357, -0.200] |
| direction of the P vector, own vectors | -0.152 [-0.192, -0.110] | -0.139 [-0.185, -0.107] | -0.342 [-0.421, -0.269] |
| 1/r amplitudes, own vectors | -0.170 [-0.208, -0.128] | -0.170 [-0.214, -0.140] | -0.423 [-0.497, -0.358] |
| moment-tensor projection, vectors re-read within 1.5 s | -0.172 [-0.207, -0.135] | -0.170 [-0.211, -0.135] | -0.327 [-0.407, -0.244] |
| moment-tensor projection with the true amplitude factors | -0.132 [-0.168, -0.092] | -0.116 [-0.158, -0.081] | -0.234 [-0.307, -0.168] |
| moment-tensor projection inside the loop | -0.137 [-0.169, -0.096] | -0.107 [-0.149, -0.077] | -0.257 [-0.337, -0.191] |

**C9: polarization admissibility as a weight.** The one positive signal of C9 is not a location effect.
The fraction of the read energy that is polarization-admissible at the node, P along the ray and S
transverse to it, separates hits from false atoms, with an area under the curve of 0.68 in dense and
0.78 in standard windows, because clutter and coda vectors are randomly polarized and true arrivals
are not; the fraction explained by the full moment-tensor projection, amplitudes included, carries almost no
ranking information, with an area of 0.56 to 0.57 in dense windows. Used as a factor of the selection statistic of MT-FW, the variant SELPOL,
it gives $+0.046$ [$+0.007$, $+0.070$] over MT-FW+SC and $+0.045$ [$+0.003$, $+0.086$] over BOMP-PK+SC in
dense windows, and $+0.045$ [$-0.002$, $+0.089$] and $+0.002$ [$-0.041$, $+0.064$] in standard windows;
used to re-rank the atoms of MT-FW after detection it gives $+0.028$ [$+0.005$, $+0.043$] in dense
windows. No candidate reaches $+0.05$ with an interval that excludes zero, the best is the best of about
50 and was added after a first look, it captures 0.12 [0.01, 0.24] of the headroom to CEIL-GD, it does
nothing at four stations ($+0.004$ [$-0.031$, $+0.051$]), and the median error of its hits is
unchanged, 2.64 km, so that the gain is in ordering and selection and not in location. It also depends
steeply on the polarization scatter of the generator: the re-ranking gain is $+0.031$ at 0 degrees,
$+0.028$ at the nominal 10 degrees, $+0.003$ [$-0.011$, $+0.023$] at 20 degrees and $-0.015$ [$-0.025$,
$+0.002$] at 45 degrees, where SELPOL gives $-0.014$ [$-0.029$, $+0.005$]. On the Colombian network,
where single-station polarization passed its reliability gate at one of five stations (Appendix T.1),
the measured dependence predicts no gain. Polarization filters and rectilinearity weights are standard
in single-station location and in polarization-weighted back-projection, so nothing here is a method.
A by-product bears on Y.24 to Y.26: in this sample of 20 networks per cell the dense gain of the capped
MT-FW over BOMP-PK at $\le0.25$ was absent, $-0.001$ [$-0.030$, $+0.041$], against $+0.049$ [$+0.015$,
$+0.099$] in the sample of C10 and $+0.066$ and $+0.055$ in Task 2I, which is the variability to expect
of an effect of that size at this sample size; and without amplitude factors the cap had no dense
effect, $+0.002$, consistent with Y.25.

**C10: soft timing priors at catalogue-grade station corrections.** C10 charges the shift of each
chosen peak from its predicted time a Gaussian log-likelihood inside BOMP-PK and MT-FW, in the units of
the detection statistic, with the timing level fixed by rule from the catalogue: the robust standard
deviation of the peak delays of the catalogued arrivals after the station correction, with a floor of
0.05 s, which came out at a median of 0.138 s with weak clutter and 0.224 s with strong clutter. The variants are an
independent charge per arrival with the variance of the level plus that of the position inside the grid
cell (T1); the same with the level alone, the literal hypothesis (T1c); a common-mode form whose prior
covariance treats the shifts as a position and origin-time term already searched plus independent
residuals, at $\pm0.25$ s (T2) and at $\pm0.5$ s (T2w); a window of $\pm0.5$ s with the charge (T4) and
without it (T4h); and a local refinement of each selected atom on a 1 km lattice within 10 km and 1 s by
the joint objective of amplitude projection and timing charge (T3). The loops in hard mode reproduce the
registered detectors in 9 of 9 validation windows, and with an exact model and exact statics the prior
and the refinement do locate better (median error 1.18 to 1.38 km for T2 and 0.85 to 1.27 km for T3
against 1.05 to 2.76 km for MT-FW). On the benchmark the hypothesis is not supported (Table Y.73). No variant
that charges the shifts gains recall under either hit rule, in either window type, in either detector:
against MT-FW+SC at $\le0.25$ under the registered rule the differences are $-0.028$ to $-0.002$ in
dense windows and $-0.066$ to $-0.006$ in standard windows, and a sweep of fixed levels of 0.05, 0.10
and 0.20 s on half of the sample changes nothing. The share of near misses does not move. Widening the
free window answers the question in the direction opposite to the hypothesis: the hard window of
$\pm0.5$ s loses $-0.047$ [$-0.087$, $-0.008$] in dense windows and the charge takes back only part of
the loss, so the cut at $\pm0.25$ s does useful work and its flat interior is not what loses the
events.

**Table Y.73. Exploration C10: paired difference in recall at $\le0.25$, point [95% interval], under the
registered hit rule (LOC) and with the origin time also required within 0.5 s (LOC+T); 600 windows per
design. The rows from T1 to T3 use the catalogue station corrections and are differences to the
detector of the same family; the oracle rows use the true station delays and the SC-R rows the
re-estimated corrections, and both are differences to MT-FW+SC.** Outputs
`results/explorations/out/c10_analysis_*.txt` and `c10r_analysis_*.txt`; `REPORT_C10.txt`, sections 4, 5 and 8.

| Variant | dense, LOC | dense, LOC+T | standard, LOC | standard, LOC+T |
|:--|:--|:--|:--|:--|
| MT-FW with T1 | -0.019 [-0.057, +0.011] | -0.003 [-0.028, +0.023] | -0.015 [-0.070, +0.034] | +0.004 [-0.055, +0.064] |
| MT-FW with T1c | -0.028 [-0.062, +0.009] | -0.013 [-0.037, +0.014] | -0.021 [-0.081, +0.045] | +0.030 [-0.033, +0.097] |
| MT-FW with T2w | -0.002 [-0.046, +0.033] | -0.018 [-0.045, +0.011] | -0.034 [-0.111, +0.019] | -0.036 [-0.103, +0.034] |
| MT-FW with T4 | -0.019 [-0.063, +0.017] | -0.013 [-0.044, +0.012] | -0.066 [-0.142, +0.000] | -0.064 [-0.121, +0.012] |
| MT-FW with T4h | -0.047 [-0.087, -0.008] | -0.049 [-0.072, -0.022] | -0.019 [-0.092, +0.049] | -0.068 [-0.141, +0.006] |
| MT-FW with T3 | +0.011 [-0.022, +0.049] | +0.005 [-0.028, +0.031] | +0.011 [-0.045, +0.067] | +0.066 [+0.019, +0.136] |
| BOMP-PK with T1 | -0.010 [-0.033, +0.007] | -0.002 [-0.015, +0.010] | +0.034 [-0.022, +0.069] | +0.008 [-0.031, +0.048] |
| BOMP-PK with T2w | -0.017 [-0.044, +0.020] | +0.014 [-0.009, +0.040] | +0.017 [-0.059, +0.085] | +0.038 [-0.029, +0.107] |
| BOMP-PK with T3 | +0.015 [-0.017, +0.050] | +0.023 [-0.002, +0.046] | +0.059 [-0.009, +0.112] | +0.100 [+0.031, +0.168] |
| oracle: true statics, hard window | +0.047 [-0.003, +0.082] | +0.034 [+0.008, +0.065] | +0.098 [+0.032, +0.141] | +0.087 [+0.031, +0.158] |
| oracle: true statics with T2w at the true model-error level | +0.081 [+0.035, +0.117] | +0.118 [+0.086, +0.149] | +0.146 [+0.065, +0.217] | +0.185 [+0.118, +0.277] |
| oracle: true statics with T3 | +0.066 [+0.024, +0.105] | +0.073 [+0.041, +0.105] | +0.089 [+0.032, +0.153] | +0.157 [+0.091, +0.226] |
| SC-R corrections, hard window (post hoc) | +0.008 [-0.031, +0.037] | +0.000 [-0.019, +0.026] | +0.030 [-0.038, +0.091] | +0.004 [-0.061, +0.080] |
| SC-R corrections with T3 (post hoc) | +0.015 [-0.023, +0.057] | +0.027 [-0.002, +0.056] | +0.064 [-0.004, +0.123] | +0.123 [+0.058, +0.201] |

**C10: the timing budget.** The reason is in the size of the timing error that a true hypothesis must
absorb (Table Y.74). The true model error, the travel time in the gradient medium of the generator
minus the homogeneous one after a perfect station term, is 0.05 s for P and 0.07 s for S. The station
corrections estimated from twenty catalogued events are wrong by 0.12 to 0.17 s rms with weak clutter
and by 0.18 to 0.24 s with strong clutter, measured against the true delays. The position of a source
inside a cell of 3 km by 3 km by 4 km, with the origin step of 0.1 s, adds 0.18 s for P and 0.30 s for
S, and the S shift of an off-node source reaches $\pm0.7$ s. The timing budget of a true hypothesis is
therefore 0.2 to 0.4 s per arrival and not the 0.144 s of Table Y.70, whose bound assumes known station
delays and a continuous position. On this lattice the window of $\pm0.25$ s is not a generous
tolerance that wastes information: it is narrower than the spread it must absorb, and a prior with the
honest width is nearly flat inside it, a charge of 0.2 to 1 unit at the edge against energies of 5 to
18 per arrival. Timing errors of this size do not blur the location of an event that is found, since
the same bound still gives 0.4 to 1.1 km; they lower the statistic of its node against false nodes
supported by clutter and coda, and the correction error is a fixed bias per station, which a longer
stack does not average out.

**Table Y.74. Exploration C10: budget of the shift of a true arrival from the prediction of its nearest
node, rms in seconds.** `REPORT_C10.txt`, section 7.1.

| Term | Size |
|:--|:--|
| true model error after a perfect station term | 0.05 (P), 0.07 (S) |
| error of the catalogue station corrections against the true delays, weak clutter | 0.12 to 0.17 |
| the same, strong clutter | 0.18 to 0.24 |
| position inside the cell of the 3 km lattice and the origin step of 0.1 s | 0.18 (P), 0.30 (S) |
| the same for a lattice of 1 km | 0.10 (P), 0.16 (S) |
| tolerance of $\pm0.25$ s read as a uniform law | 0.144 |

**C10: the oracle decomposition.** With the true station delays as corrections and the timing level at
the true model error, the gain exists and is large (Tables Y.73 and Y.75). The true delays alone, with
the hard window unchanged, give $+0.047$ [$-0.003$, $+0.082$] in dense and $+0.098$ [$+0.032$, $+0.141$]
in standard windows over MT-FW+SC. The true delays with the common-mode prior give $+0.081$ [$+0.035$,
$+0.117$] and $+0.146$ [$+0.065$, $+0.217$], and $+0.118$ and $+0.185$ with origin times matched. In
standard windows that oracle detector reaches the realistic ceiling, a recall of 0.592 against 0.603 for
CEIL-GD; in dense windows it captures 0.31 [0.21, 0.43] of the headroom. The prior pays only when the
timing is known to a few hundredths of a second: with the true delays but the catalogue level of 0.14
to 0.22 s it adds nothing over the hard window, $-0.008$ [$-0.039$, $+0.032$] in dense and $-0.059$
[$-0.108$, $+0.027$] in standard windows. The oracle uses the delays of the generator, which a real
network never has, so it bounds what calibration could give and promises nothing.

**Table Y.75. Exploration C10: recall at $\le0.25$ of the baselines, the oracle detectors and the
realistic ceiling, point [95% interval] under the registered rule and point under the origin-time rule;
600 windows per design.** Outputs `results/explorations/out/c10_analysis_DENSE.txt` and
`c10_analysis_STD.txt`; `REPORT_C10.txt`, sections 4.1 and 5.1. CEIL-GD has no origin-time version.

| Detector | dense, LOC | dense, LOC+T | standard, LOC | standard, LOC+T |
|:--|:--|:--|:--|:--|
| BOMP-PK+SC | 0.198 [0.161, 0.229] | 0.067 | 0.439 [0.379, 0.521] | 0.272 |
| MT-FW+SC, cap 6 | 0.247 [0.218, 0.285] | 0.116 | 0.446 [0.378, 0.516] | 0.280 |
| MT-FW, true statics, hard window | 0.294 [0.257, 0.325] | 0.150 | 0.544 [0.470, 0.596] | 0.367 |
| MT-FW, true statics, T2w at the true level | 0.329 [0.291, 0.364] | 0.234 | 0.592 [0.518, 0.651] | 0.465 |
| CEIL-GD | 0.616 [0.539, 0.696] | | 0.603 [0.516, 0.691] | |

**C10: what remains after timing.** A diagnostic before any deflation, post hoc and on the standard
windows of half of the sample (235 events), compares the best node statistic within 5 km of each true
event with the best value of the window among nodes farther than 5 km from every event. Under the
registered statistic 0.485 [0.420, 0.544] of the isolated events win that comparison. The prior at the
catalogue level leaves the share where it was, 0.489, with point changes of opposite sign in the two
regimes, neither of which excludes zero. The oracle prior raises it to 0.609 [0.540, 0.674], and the oracle objective on
the 1 km lattice to 0.634 [0.586, 0.678]. With perfect station delays and the true model-error level, 37
to 39 per cent of the isolated events therefore still lose to a false node at the first iteration.
Those events are lost to amplitude, that is, to weak arrivals, clutter, coda and the amplitude factors,
and not to timing; it is the same wall as the ceiling of 0.60.

**C10: refinement and re-estimated corrections.** The one candidate with a signal at the catalogue
level is the refinement T3. Unlike the relocation of C9 its point estimates are never negative, and in standard windows
it lowers the median error of the hits from 2.51 to 2.07 km for BOMP-PK and from 2.42 to 2.21 km for
MT-FW, against 1.51 km for the oracle. Under the registered hit rule its gain does not exclude zero:
$+0.059$ [$-0.009$, $+0.112$] for BOMP-PK and $+0.011$ [$-0.045$, $+0.067$] for MT-FW in standard windows,
and $+0.015$ and $+0.011$ in dense windows. It helps only when a hit must also have the origin time
within 0.5 s, and then only in standard windows: $+0.100$ [$+0.031$, $+0.168$] for BOMP-PK and $+0.066$
[$+0.019$, $+0.136$] for MT-FW. That rule is not the registered endpoint, the contrast is the best of
about 70, and the rule rewards T3 partly for searching the origin time again within 1 s, an effect that
was not separated from the move of the node. Under that rule both baseline detectors keep 0.27 to 0.28 of the
isolated events against 0.44 to 0.45 under the registered rule, so more than a third of the registered
hits carry an origin time wrong by more than 0.5 s, which is what T3 repairs. The post hoc corrections
SC-R, re-estimated with the catalogued events relocated by their own times, cut the correction error
from 0.145 to 0.110 s with weak clutter and from 0.214 to 0.202 s with strong clutter. Alone they
change nothing visible, $+0.008$ in dense and $+0.030$ [$-0.038$, $+0.091$] in standard windows. With T3
they give in standard windows $+0.064$ [$-0.004$, $+0.123$] at $\le0.25$ and $+0.102$ [$+0.028$,
$+0.196$] at $\le0.1$ under the registered rule and $+0.123$ [$+0.058$, $+0.201$] with origin times
matched, about two thirds of the oracle refinement. The direction agrees with the budget, the soft
timing terms starting to pay as the absolute timing error falls, but the sizes are post hoc, the rule
was written after the oracle had been seen and was checked against the true delays, and SC-R does not
help with strong clutter at four and six stations, where wrong catalogue peaks and not catalogue
location errors dominate the correction error.

**Limits.** The station-delay standard deviation of 0.3 s, the catalogue of twenty events with errors
of 1 km and 0.1 s and the peak rule of $\pm1$ s of the corrections are design choices of the benchmark,
and the size of the correction error, with it the whole conclusion of C10, depends on them. The cell
term belongs to the 3 km lattice; a detection lattice of 1 km, at nine times the cost, would face a
budget of 0.10 to 0.16 s plus the correction error and was not run. The polarization scatter of 10
degrees and the amplitude spread of 0.5 are generator choices, and the finding that they matter little
for location is conditional on the signal-to-noise ratio of the benchmark. Both studies use one seed
with 20 networks per cell and regime. Soft timing tolerances are themselves established: probabilistic
location with Gaussian model errors, back-projection with travel-time uncertainty
windows and joint hypocentre determination with station terms all contain the idea, T2 is the textbook
Gaussian marginalisation of position inside a cell and SC-R is joint hypocentre determination with a
prior. None of the results of this section is a method claim.

### Y.27a What C9 and C10 change

C9 and C10 do not add a result to the confirmatory record. They change the description of the open
problem. After Y.25 the headroom below the realistic ceiling stood as unexplained by the grid, the
velocity model, the source pulse and every treatment of amplitudes. On the calibration seed, and with
oracles that a real network does not have, it splits into two parts.

The first part is calibration, that is, how accurately the station terms are known. With the station
delays known to the level of the model error, a few hundredths of a second, the same moment-tensor
statistic with a matching timing prior gains $+0.146$ in standard windows and reaches the ceiling for
isolated events, 0.592 against 0.603, and about two thirds of that gain, $+0.098$, comes from the true
delays alone with the registered hard window. With corrections of catalogue grade, wrong by 0.12 to 0.24 s
rms, and a lattice whose cells add 0.18 to 0.30 s, no principled timing prior helps, because the honest
width of the prior is as wide as the window. The second part is amplitude. With timing made perfect, 37
to 39 per cent of the isolated events still lose to a false node at the first iteration, and in dense
windows the oracle captures 0.31 of the headroom, so about two thirds of the dense headroom is not
timing: it is weak arrivals, clutter, coda and the amplitude factors, the part on which the bounded
influence of Y.24 to Y.26 acts and of which Y.25 found that no amplitude treatment captured more than a
quarter.

Travel-time-free physics is not a third part. In this benchmark it is limited by information and not
by algorithm: the bound of all such channels together is 4.4 to 5.6 km against 0.4 to 0.7 km for the
times, a joint estimator is bounded at an improvement of 4 to 5 per cent, relocation by such evidence
loses recall even with the true amplitude factors and with a perfect generator, and the one positive
use, polarization admissibility as a weight against clutter, is small, selected post hoc and gone at a
scatter of 20 degrees. It would matter at about ten times the signal-to-noise ratio with a scatter of
at most 2 degrees, or where the times are unusable, and neither is this benchmark nor, on the evidence
of Appendix T.1, the Colombian network.

The open problem is therefore restated. It is not a further location principle. It is calibration,
meaning station terms estimated self-consistently with the detections instead of once from a short
catalogue, to an accuracy near the model error, and amplitude evidence for weak arrivals in clutter and
coda. Three consequences follow for the plan. Station-term calibration is the lever to test next, with
the true-delay oracle run beside it so that any result is read against its bound. The endpoint with
the origin time matched belongs in the registered protocol, since more than a third of the hits under
the location-only rule carry a wrong origin time and the refinement that repairs them is invisible to
that rule. And the lines of relocation by polarization or amplitude and of soft timing priors at
catalogue-grade corrections are closed as negatives with their mechanism.

None of this is confirmatory. It rests on one seed, 20 networks per cell and regime, oracles, post hoc
analyses declared as such and no multiplicity correction; the numbers are pilot estimates, and the
statement that calibration would close the headroom for isolated events is a bound on this generator,
with its station delays of 0.3 s and its catalogue of twenty events, and not a prediction for real
data.

### Y.27b A literature map of travel-time-free location

The question of C9 was preceded by a scan of the literature from 2000 to 2026 on detecting and
locating earthquakes without travel times (`explore_contrib/LIT_NO_TRAVELTIME.md`). It asked two
questions of eight families of methods: whether the family avoids travel times or hides them, for
example in the phase of a Green's function, and whether it can add information independent of the
arrival times (Table Y.76). Four kinds of datum recorded by an existing three-component network are
independent of the arrival times: the direction of P particle motion at each station, the amplitude
ratio between phases at one station, the amplitude ratios between stations, and high-frequency energy
with its coda decay. Wavefield gradients, rotation and strain are independent too and need hardware
that a sparse network does not have. The other families carry the kinematics inside a propagator or
inside learned weights: time reversal and wave-equation imaging, matched field processing,
back-projection, envelope cross-correlation, templates and end-to-end neural location.

**Table Y.76. Families of methods that locate without picked travel times, as sorted by the scan.**
`explore_contrib/LIT_NO_TRAVELTIME.md`, section 1; the references are examples read at abstract level or
better in the scan, except [46], [113] and [129], which the scan read at the level of bibliographic
metadata.

| Family | Observable | Avoids travel times | Independent information | On a network of 4 to 8 three-component stations |
|:--|:--|:--|:--|:--|
| Polarization [145] | direction of P and S particle motion | yes for direction; distance still needs S minus P or amplitude | azimuth, and incidence if the near-surface velocity is known | yes, station by station, with reliability that varies by site |
| Amplitude location and ratios [146] | band-limited amplitude or envelope per station | yes, with coarse window alignment only | distance decay and radiation pattern | yes, with site terms dominating the error budget |
| Time reversal, wave-equation imaging and waveform source inversion [149], [141] | waveforms back-propagated through a velocity model | no, hidden in the propagator | amplitude and radiation-pattern consistency if elastic | poorly, since it needs aperture and density |
| Gradiometry, rotation and strain [148] | spatial gradients, rotation, strain | locally: it measures the slowness vector without picks | direction and apparent velocity at one site | no, the hardware is absent |
| Matched field, templates and back-projection [150], [15], [129] | match to modelled or empirical Green's functions | no, hidden in the phase of the Green's function | only if amplitudes are kept and trusted | partly; empirical versions are tied to known source regions |
| Envelope, coda and energy [147] | envelopes, coda decay, energy | envelope correlation no; energy and diffusion models mostly yes | energy decay with distance, robust to overlap | yes at high frequency, at coarse resolution |
| Machine learning from waveforms [97] | raw waveforms | no, learned implicitly with path fingerprints | implicit use of S minus P, polarization and coda shape | yes, but tied to a region and prone to outliers |
| Inverse-source theory and sparse recovery [46], [47], [113] | boundary data at several frequencies | not separable | a framework, not an observable | the theory assumes full boundary data |

The verdict of the scan, written before C9 and C10 were run, was that polarization direction is the
cleanest datum independent of travel times on an existing network, that the amplitude decay cannot
resolve an ambiguity of 5 to 10 km without calibrated site terms, and that an off-grid refinement of
the dictionary should be settled before any further physics is added, because a grid solution
concentrates near the true positions and not on them [47]. C9 tested the first and the second and found
both limited by information at the signal-to-noise ratio of the benchmark; the refinement T3 of C10 is a
form of the third and helps location and origin time without changing recall under the registered rule.
The published accuracy of the families that the scan could read is in line with the bound of Table
Y.70: single-station learned location reports mean epicentral errors of several kilometres [97], and
envelope-based energy inversion locates only approximately [147], which is the size of the mislocation
to be removed and not a route to removing it.

The scan is not a systematic review and its limits are those it declares. Of its 72 entries, 26 were
read at the level of bibliographic metadata only, more than a third, and what is said of those is
limited to what the title supports. The article pages of several publishers could not be opened, so
that a number of papers that belong in the map, among them work on single-station neural location,
matched field processing and continuous moment-tensor scanning, were left out and are not evidence
here. Sentences on typical accuracies that the scan marks as background were not checked against a
source; the only accuracy for P polarization that it read is one report of 70 per cent of events within
30 degrees. It found no paper that combines a moment-tensor dictionary detector with an explicit
polarization likelihood, and states that absence in a search limited by rate limits is not evidence of
novelty. The search protocol of Appendix E applies to this map as to every other claim of absence.

## Y.28 Station terms estimated without picks, and the refined detector (C11 and Task 2L)

Code: `explore_contrib/c11_lib.py` (peak read, joint relocation with station terms, self-calibration,
term errors, engine cloning), `c11_run.py`, `c11_dev.py`, `c11_validate.py`, `c11_analyze.py`,
`c11_posthoc.py` and `c11_power.py`, with the detector loops and the refinement of
`explore_contrib/c10_lib.py`; plan `explore_contrib/C11_PLAN.txt` with its hashes in `C11_PLAN.sha256`;
report `results/explorations/REPORT_C11.txt`; outputs `results/explorations/out/c11_*.txt`. For the confirmatory
task: `task2l_run.py` (controls and runs, resumable from one part per network), `task2l_control.py`
(controls on seed 999) and `analyze_2l.py` (rules, tables and hierarchical bootstrap); outputs
`analyze_2l_out.txt`, `task2l_control.txt`, `task2l_code_hashes.txt`, `task2l_freeze_time.txt` and
`log2l.txt`; design and rules in `PREREGISTRO_1j_1k.txt` (section "Task 2L") with the hashes in
`PREREGISTRO_hashes.txt`; independent verification in `verif_2l/task2l_verification.txt`.

**Status.** This section has two parts of different standing. The exploration C11 ran on seed 999 only
and tests no registered hypothesis; its plan was frozen with a hash before any detection result of the
study existed, which fixes its primary contrast but does not make it confirmatory. Task 2L is a
preregistered test of the candidate that C11 named, on two seeds that no earlier design had generated,
with its registration hash sent to the user before any run, its code frozen by hashes and an independent
verification. Its four rules and its registered secondary passed in both seeds. The claim that they
support is narrow and is stated with its scope at the end of the section.

**The question.** Y.27 left the open problem restated as calibration: with the true station delays the
moment-tensor detector reached the realistic ceiling for isolated events, while the station corrections
estimated from twenty catalogued events were wrong by 0.12 to 0.24 s and no timing prior helped at that
accuracy. The post hoc corrections SC-R of Y.27, which relocate the catalogued events by their own
times, pointed to the lever and were chosen after the oracle result had been seen. C11 asked the
question with the order of events put right: can the station terms be estimated better from the same
twenty catalogued events, without picks and without any quantity that a network does not have, and does
the detector then detect more.

**C11: data and freeze.** C11 used networks of seed 999 that no earlier design had touched, and it kept
calibration and detection apart. Networks 200 to 209 of each cell and regime, 60 combinations of network
and regime, were the development set, on which the only quantity computed was the error of the
estimated terms against the true delays. Networks 210 to 239, 30 per cell and regime, were the
evaluation set, on which detection was measured; network 195 served for validation. The estimator, its
constants, the run set, the endpoints, the bootstrap, the rule naming the primary candidate, the primary
contrast and the judgement rules were written in `C11_PLAN.txt`, whose body has the SHA-256
42e6613879e6f3f2aae4f1cb349d15bd2839770e4e1ecbcbf0f8ae772c181c5f and was frozen at 18:01:51 UTC on 4
October 2026; the evaluation runs started at 18:02:00 UTC with the code hashes of the plan checked. No
recall, FROC or hit count of any detector on any network of the study had been computed before the
freeze. The plan declares three things seen before it: the development pass of the self-calibration ran
with the terms of an earlier form of the estimator, a purity figure of the atoms selected for
self-calibration was computed with the truth on the development set and used by no rule, and a smoke
test ran on the validation network. The protocol is that of Y.27: the grid of 675 nodes at 3 km and
three depths, origin times from 5.0 to 8.0 s at 0.1 s, the hard peak rule of $\pm0.25$ s, the 5 km
exclusion, the exact FROC, five standard and five dense windows per network, 900 windows per design with
703 and 3615 events, and paired differences with 95 per cent intervals from a hierarchical bootstrap of
2000 replicates. The primary hit rule is the one that Y.27 recommended, an epicentre within 5 km and an
origin-time error of at most 0.5 s; the location-only rule of the earlier tasks is the second reading.
No multiplicity correction is applied.

**C11: the candidates.** Six sets of station terms, one term per station and phase, fixed per network,
enter the same detectors as corrections. S0 is the registered correction of Tasks 2A to 2I: the median,
over the twenty catalogued events held at their catalogue hypocentres, of the delay of the
matched-filter energy peak read within $\pm1$ s of the predicted time, with an energy gate of 25. S1 is
SC-R of Y.27, unchanged. S2 is the estimator of this section, described below. S3 adds to S2 a
self-calibration from the detector's own output: the refined detector with the S2 terms runs over 30
unlabelled windows of the network, disjoint from the evaluation windows, the atoms with a score at or
above the 0.80 quantile of the stream enter as further events with a wide prior of 3 km and 0.3 s, and
catalogue and atoms go through two more passes of the S2 procedure. S4 is S2 with the first 40 and the
first 80 events of the catalogue stream in place of 20. The oracle uses the true station delays. No
truth enters S0 to S4.

**The estimator S2.** S2 is a joint determination of the hypocentres of the catalogued events and of the
station terms that never uses a pick. Its data are the matched-filter energy peaks of the catalogued
waveforms. It starts from the S0 terms and the catalogue hypocentres and makes three passes. The first
pass uses the registered read. The second and the third read every event again, within $\pm0.5$ s of
the time predicted by its current hypocentre and the current terms and with an energy gate of 12, so
that the weaker arrivals, which the registered gate discards, are brought in once the prediction is
good enough to find them. Each pass has two stages. In the first, an alternation of at most eight
sweeps, the residual scale is set to 1.4826 times the median absolute deviation of the gated residuals
with a floor of 0.05 s; every event with at least four gated arrivals is relocated by iteratively
reweighted Gauss-Newton, with a Huber misfit of constant 1.345 in homogeneous straight-ray times and a
Gaussian prior of 1 km per coordinate and 0.1 s on the origin around its catalogue hypocentre, events
with fewer arrivals staying at the catalogue hypocentre; and each term becomes the median over events
of the read time minus the travel time minus the origin. In the second stage a joint Gauss-Newton
iteration, of at most 15 steps with Cauchy weights of constant 2.385, solves for all terms and all
relocatable hypocentres at once under the same priors. After each stage the gauge is fixed: the median
shift of the origin times of the relocated catalogue events is set to zero. The gauge is not a zero
mean of the terms per phase. The true delays do not have zero mean and a common offset moves every
origin time, which the origin-time hit rule sees, so the absolute time is tied to the catalogue origins
and the remaining weak directions are held by the catalogue prior. This choice, the relocation by
Gauss-Newton with a prior instead of a grid search, and the joint stage are deviations from the first
wording of the study, written in the plan before its freeze. On the development set, with the term
error as the only criterion and a rule that changed a default only for a gain above 0.01 s, two things
were kept: the gate of the re-read, which lowered the pooled median error from 0.111 to 0.068 s with
weak clutter and from 0.211 to 0.135 s with strong clutter, and the joint stage, which lowered it
further to 0.062 and 0.111 s, against 0.133 and 0.190 s for S0. Window widths, numbers of passes, robust
constants and the minimum number of arrivals changed the error by less than 0.01 s. A weaker position
prior of 2 km was worse, 0.118 and 0.186 s: the estimate leans on the catalogue prior, whose widths are
the true sizes of the catalogue errors of the benchmark. In a noise-free validation with an exact model
the joint estimator without prior recovered the true terms to 0.001 to 0.005 s up to the gauge, and the
frozen version, with its floor and prior, left 0.020 to 0.028 s. The calibration costs 0.38 s per
network.

**C11: term errors.** On the evaluation networks S2 lowered the median error of the terms from 0.150 to
0.072 s with weak clutter and from 0.206 to 0.139 s with strong clutter (Table Y.77), close to the
development figures. With twenty events it is as accurate as the plain median with eighty events when
clutter is weak and lies between the medians with forty and eighty when clutter is strong. It is not
safe for every network: with strong clutter it is worse than S0 in about a third of the networks and
its upper quartile, 0.237 s, is no better than that of S0, 0.261 s, because contaminated catalogue
peaks remain the limit with twenty events and four to six stations. The common offset of the terms,
which no estimator of this kind can see because it is the mean origin error of the catalogue, stays at
0.04 to 0.05 s for every twenty-event estimator.

**Table Y.77. Exploration C11, evaluation networks 210 to 239 of seed 999: rms error of the station terms
against the true delays with the median removed, in seconds; median over the 30 networks of each cell
and median [quartiles] over the 90 networks of the regime; share of networks with a smaller error than
S0; median absolute common offset.** Report `results/explorations/REPORT_C11.txt`, section 5. S0 at 40 and 80
events is the plain median with that many catalogued events.

| Regime | Terms | 4 stations | 6 stations | one-sided, 8 st. | all networks | share better than S0 | common offset |
|:--|:--|:--|:--|:--|:--|:--|:--|
| weak | S0, registered | 0.131 | 0.161 | 0.149 | 0.150 [0.127, 0.187] |  | 0.041 |
| weak | S1, SC-R | 0.106 | 0.127 | 0.100 | 0.107 [0.088, 0.141] | 0.86 | 0.043 |
| weak | S2 | 0.088 | 0.097 | 0.056 | 0.072 [0.052, 0.111] | 0.84 | 0.049 |
| weak | S3, self-calibration | 0.081 | 0.088 | 0.049 | 0.069 [0.047, 0.109] | 0.81 | 0.044 |
| weak | S2, 40 events | 0.065 | 0.061 | 0.045 | 0.051 [0.038, 0.077] | 0.96 | 0.034 |
| weak | S2, 80 events | 0.041 | 0.052 | 0.038 | 0.041 [0.034, 0.059] | 0.99 | 0.027 |
| weak | S0, 40 events | 0.088 | 0.115 | 0.125 | 0.109 [0.087, 0.132] | 0.83 | 0.034 |
| weak | S0, 80 events | 0.068 | 0.083 | 0.085 | 0.081 [0.067, 0.092] | 0.98 | 0.028 |
| strong | S0, registered | 0.163 | 0.221 | 0.220 | 0.206 [0.162, 0.261] |  | 0.045 |
| strong | S1, SC-R | 0.164 | 0.204 | 0.163 | 0.185 [0.141, 0.236] | 0.72 | 0.039 |
| strong | S2 | 0.136 | 0.191 | 0.105 | 0.139 [0.088, 0.237] | 0.68 | 0.054 |
| strong | S3, self-calibration | 0.139 | 0.150 | 0.089 | 0.134 [0.083, 0.236] | 0.68 | 0.051 |
| strong | S2, 40 events | 0.102 | 0.089 | 0.056 | 0.082 [0.054, 0.132] | 0.94 | 0.038 |
| strong | S2, 80 events | 0.068 | 0.055 | 0.041 | 0.052 [0.039, 0.076] | 0.99 | 0.027 |
| strong | S0, 40 events | 0.137 | 0.161 | 0.146 | 0.150 [0.122, 0.183] | 0.87 | 0.037 |
| strong | S0, 80 events | 0.107 | 0.108 | 0.104 | 0.105 [0.090, 0.130] | 0.99 | 0.028 |

**C11: detection.** The primary contrast of the plan holds the detector fixed, MT-FW with the cap of 6
and the refinement T3 of Y.27, and changes only the terms: with S2 in place of S0, in standard windows,
under the origin-time hit rule and at $\le0.25$ false events per window, recall rose from 0.371 to
0.441, a difference of $+0.070$ [$+0.010$, $+0.120$] (Table Y.78). The true delays gave $+0.091$
[$+0.031$, $+0.141$] on the same contrast, so that S2 captured 0.77 [0.18, 1.45] of the oracle gain. In
dense windows the same contrast was $+0.065$ [$+0.034$, $+0.094$] with origin times matched and $+0.034$
[$+0.001$, $+0.070$] under the location-only rule. Both judgement rules of the plan were met, and the
plan's verdict was a plausible contribution, to be tested.

Four qualifications came with that verdict. The effect needs the refinement: with the hard rule of
$\pm0.25$ s alone the S2 terms gave $+0.026$ [$-0.030$, $+0.083$] over S0 in standard windows and the
true delays $+0.020$ [$-0.044$, $+0.072$], so the terms and the refinement are complements, and against
the registered detector the refinement gave $+0.024$ [$-0.022$, $+0.078$] with the S0 terms, $+0.094$
[$+0.035$, $+0.151$] with the S2 terms and $+0.115$ [$+0.061$, $+0.174$] with the true delays. The gain
belongs to detection that is consistent in origin time: under the location-only rule the primary
contrast was $+0.017$ [$-0.028$, $+0.068$] in standard windows and the oracle gave the same, $+0.018$.
The lower end of the primary interval was $+0.010$, on one seed and 30 networks per cell and regime.
And in block pursuit the same contrast had the same direction with an interval touching zero, $+0.057$
[$-0.003$, $+0.115$].

**Table Y.78. Exploration C11: paired differences in recall at $\le0.25$, point [95% interval], pooled over
regimes and cells; 900 windows per design. Every row has the refinement T3; the first two columns are
standard windows with origin times matched.** Report `results/explorations/REPORT_C11.txt`, sections 6.2,
6.3 and 7.2.

| Terms, with T3 | standard, minus S0 with T3 | standard, minus S0 without T3 | dense, origin time matched, minus S0 with T3 | standard, location only, minus S0 with T3 | dense, location only, minus S0 with T3 |
|:--|:--|:--|:--|:--|:--|
| MT-FW, S1 | +0.023 [-0.050, +0.077] | +0.047 [-0.029, +0.109] | +0.042 [+0.014, +0.073] | +0.000 [-0.048, +0.045] | +0.020 [-0.008, +0.054] |
| MT-FW, S2 (primary) | +0.070 [+0.010, +0.120] | +0.094 [+0.035, +0.151] | +0.065 [+0.034, +0.094] | +0.017 [-0.028, +0.068] | +0.034 [+0.001, +0.070] |
| MT-FW, S3 | +0.085 [+0.035, +0.140] | +0.110 [+0.058, +0.171] | +0.056 [+0.025, +0.085] | +0.031 [-0.016, +0.077] | +0.027 [-0.006, +0.061] |
| MT-FW, S2 with 40 events | +0.080 [+0.019, +0.137] | +0.104 [+0.044, +0.163] | +0.051 [+0.022, +0.080] | +0.028 [-0.015, +0.085] | +0.040 [+0.008, +0.071] |
| MT-FW, S2 with 80 events | +0.087 [+0.029, +0.138] | +0.111 [+0.057, +0.168] | +0.072 [+0.041, +0.101] | +0.040 [-0.015, +0.084] | +0.035 [+0.002, +0.067] |
| MT-FW, true delays | +0.091 [+0.031, +0.141] | +0.115 [+0.061, +0.174] | +0.078 [+0.048, +0.109] | +0.018 [-0.029, +0.070] | +0.050 [+0.020, +0.082] |
| BOMP-PK, S2 | +0.057 [-0.003, +0.115] | +0.100 [+0.037, +0.166] | +0.019 [-0.003, +0.041] | +0.031 [-0.018, +0.080] | +0.031 [-0.012, +0.068] |
| BOMP-PK, true delays | +0.110 [+0.055, +0.171] | +0.152 [+0.099, +0.215] | +0.031 [+0.008, +0.059] | +0.065 [+0.018, +0.117] | +0.045 [+0.011, +0.083] |

**C11: what adds nothing.** Self-calibration from the detector's own atoms is not distinguishable from
S2. The term error of S3 is 0.069 and 0.134 s against 0.072 and 0.139 s, and its recall differs from
that of S2 by $+0.016$ [$-0.019$, $+0.064$] in standard and $-0.009$ [$-0.032$, $+0.018$] in dense
windows, a contrast between candidates that the plan did not fix and that is post hoc. Its point
estimate in standard windows is the highest of the twenty-event estimators, but the plan named S2 by a
rule fixed beforehand, a tie in development error going to the simpler estimator, and the difference is
noise. The development set shows why: about half of the selected atoms are false and the true ones are
located to about 2 km, so that 14 to 18 usable atoms per network add little to 19 relocatable catalogue
events, and the pass costs 13 s per network. The size of the catalogue does not resolve either. With 40
and 80 events the term error keeps falling, to 0.051 and 0.041 s with weak clutter and to 0.082 and
0.052 s with strong clutter, while the recall gain moves from $+0.070$ to $+0.080$ and $+0.087$ against
$+0.091$ for the oracle; S2 with 80 events exceeds S2 with 20 by $+0.017$ [$-0.028$, $+0.063$], and the
oracle exceeds S2 with 20 events by $+0.021$ [$-0.026$, $+0.071$] in standard and $+0.013$ [$-0.014$,
$+0.045$] in dense windows. Recall saturates near a term error of 0.07 s, which is about the model error
of the homogeneous travel times after a perfect station term and below the term of a 1 km cell in T3.
SC-R did not replicate as an effect of the terms in standard windows, $+0.023$ [$-0.050$, $+0.077$],
and was positive in dense ones.

**C11: shrinkage of the C10 estimates.** Two estimates of Y.27, measured there on 20 networks per cell and
regime, were smaller on the 30 networks of C11. The refinement T3 with the registered terms gave
$+0.024$ [$-0.022$, $+0.078$] under the origin-time rule in standard windows, where C10 had found
$+0.066$ [$+0.019$, $+0.136$]; and the true delays with the hard rule gave $+0.020$ under the
origin-time rule and $+0.018$ under the location-only rule, where C10 had found $+0.087$ and $+0.098$.
The C10 contrast for T3 had been the best of about 70 (Y.27). The report of C11 therefore expected its own primary
estimate, itself the selected candidate of a study, to shrink in the same way, and the registration
that followed plans its sample on $+0.05$ and not on $+0.070$. After the freeze, the analysis script of
C11 was changed once, by three lines that skip a ratio whose methods had not been run and at which the
frozen script stopped; no other number changed, and the frozen file is kept beside it.

**Task 2L: preregistration.** The question registered is whether station terms estimated by S2 from the
registered catalogue of twenty events improve the origin-time-consistent detection of the refined
detector, without loss under the location-only rule, in standard and in dense windows. It was appended
to the registry at 21:10:33 UTC on 4 October 2026, before any code was run on seeds 7171 and 9393. The
hash of the registry at that point,
13a6e9b1e80764a85870b8414751d8338acac806bbd68b8d3e66ac6a2f60a62a, was recorded with its clock time and
sent to the user in the working session before any run. As for Tasks 2E, 2F and 2I, that record lies
outside the computing environment and is neither public nor issued by a third party; the verification
could confirm from disk the hash and its position in the registry, not the moment at which it was sent.
Drafts of registrations in the reports of C9 and C10 had named the same two seeds; neither had been
registered or run, the protocol of 2L contains the arms of the C10 draft, and the C9 draft was
discarded. The registration fixes every constant at the value of the frozen plan of C11 and leaves none
free.

The primary hit rule is an epicentre within 5 km and an origin-time error of at most 0.5 s; the
location-only rule is secondary. The endpoint is recall at $\le0.25$ realized false events per window on
the exact FROC, with $\le0.1$ as secondary, one operating point per method, and a hierarchical bootstrap
of 2000 replicates with regime and cell as strata, networks and then windows. The contrast of every rule
is MT-FW with the S2 terms and T3 minus MT-FW with the S0 terms and T3, so that the detector is the same
in both arms and only the terms differ. Rule 2L-a (superiority, standard windows, origin time matched,
pooled): one-sided 97.5 per cent lower bound above zero and a point estimate of at least $+0.03$. Rule
2L-b (superiority, dense windows, origin time matched, pooled): one-sided 97.5 per cent lower bound
above zero. Rule 2L-c (non-inferiority per clutter regime, margin 0.02, origin time matched, standard
and dense): one-sided 95 per cent lower bound above $-0.02$ in each regime. Rule 2L-d (non-inferiority
under the location-only rule, margin 0.02, pooled, standard and dense): one-sided 95 per cent lower
bound above $-0.02$. Each rule must hold in both seeds. The claim that station terms from the pick-free
joint determination improve origin-time-consistent detection of the refined detector at no cost under
the location-only rule requires a, b, c and d; if a or b fails the pilot was an artefact of selection;
if c or d fails with a and b the gain exists at a cost and is reported so. Superiority under the
location-only rule was not registered. One secondary was registered, the same contrast in block
pursuit, in standard windows on 60 networks per cell and regime, with a one-sided 97.5 per cent lower
bound above zero. The oracle arm, the contrasts against the registered detector without refinement, the
per-cell values, the term errors and the errors of the hits are descriptive.

**Task 2L: design.** Each seed is an independent replication. Networks $j=0$ to 159 of each of the three
cells and two clutter regimes, 160 per cell and regime, carry five standard windows, with zero to two
events and at most five atoms; networks 0 to 99, 100 per cell and regime, carry in addition five dense
windows, with three to five events within 9 km and at most seven atoms. The windows come from two random
streams used by no earlier design. Per seed there are 4800 standard windows with 3826 and 3863 events
and 3000 dense windows with 11978 and 12000 events. The station delays have a standard deviation of 0.3
s, the catalogue has twenty events with errors of 1 km and 0.1 s as in Task 2A, and the generator is
that of the earlier tasks. On every network MT-FW with the cap of 6 runs three times, with the S0 terms
and the hard rule, which is the registered detector of the earlier tasks, with S0 and T3, and with S2 and
T3. The block-pursuit family runs in the same three forms on the standard windows of networks 0 to 59,
1800 windows per seed, and the oracle, the true delays with T3, on networks 0 to 29, 900 windows per
design and seed; both restrictions were set by the compute budget in the registration. The planned
standard deviations of the rule contrasts, the pilot values scaled by the square root of the ratio of
networks, were 0.0122 in standard and 0.0083 in dense windows; the measured ones were 0.0119 and 0.0121,
and 0.0077 and 0.0073.

**Freeze and interruptions.** The SHA-256 hashes of 37 files, the runner, the controls, the analysis and
every module that they load, were written to `task2l_code_hashes.txt`, and the freeze was stamped at
21:24:48 UTC; the first job on the evaluation seeds was launched at 21:24:53 UTC, and the 1920 parts
were written between 21:25:19 UTC on 4 October and 01:48:15 UTC on 5 October. The run did not complete
in one launch. It was interrupted at 22:47:19 UTC, when the working session was suspended, and the
record of the freeze notes a resumption at 22:48:13 UTC and relaunches at 22:49:32 and 23:07:22 UTC, each
with the code unchanged (`sha256sum -c`, 37 of 37). Each resumption used the same command, read the
parts already written, one file per network, and recomputed none of them. The sample size had been
fixed in the registration, and the budget rule, which would have lowered it to 120 and 80 networks had
the time estimate exceeded 10 h on two cores, was not triggered, the estimate being 4.62 h; the
interruptions therefore gave no opportunity for optional stopping. The runs took 4.34 and 4.33 h of
processor time per seed. The registry received the result paragraph at 02:00:25 UTC and the paragraph
of the verification at 02:36:52 UTC on 5 October, each with a hash.

**Controls.** The controls on seed 999 (`task2l_control.txt`) passed. The runner reproduced the pilot
exactly on 18 parts, three networks by three cells by two regimes: the S0, S2 and oracle terms, their
T3 timing levels and their errors in 18 of 18, and the truth, the origin times and the atoms of every
method in 90 of 90 standard and 90 of 90 dense windows. The S0 of the runner equalled the registered
correction in 18 of 18 parts, and the corrections of the engine equalled those stored by Task 2A on seed
999 in 30 of 30 cases with a largest difference of zero. The seed enters the networks, the catalogues
and the window streams as designed, which was checked with a test seed that is none of the registered
ones. The analysis reproduced the published recalls of the pilot for every method and both hit rules,
and its bootstrap returned $+0.070$ [$+0.011$, $+0.120$] for the primary contrast of C11 against the
published $+0.070$ [$+0.010$, $+0.120$].

**Table Y.79. Task 2L. Event recall at realized false events per window $\le0.1$ / $\le0.25$, pooled over
the three cells, for MT-FW with the cap of 6 under three treatments: the registered terms S0 with the
hard rule, S0 with the refinement T3, and the S2 terms with T3.** Script `analyze_2l.py`; output
`analyze_2l_out.txt`, first table.

| Windows | Hit rule | Regime | Seed | windows | events | MT-FW, S0 | MT-FW, S0 and T3 | MT-FW, S2 and T3 |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| standard | origin time matched | pooled | 7171 | 4800 | 3826 | 0.212/0.315 | 0.244/0.346 | 0.297/0.403 |
| standard | origin time matched | pooled | 9393 | 4800 | 3863 | 0.186/0.292 | 0.237/0.350 | 0.293/0.411 |
| standard | origin time matched | weak | 7171 | 2400 | 1906 | 0.261/0.378 | 0.317/0.428 | 0.394/0.512 |
| standard | origin time matched | weak | 9393 | 2400 | 1914 | 0.230/0.357 | 0.302/0.433 | 0.364/0.510 |
| standard | origin time matched | strong | 7171 | 2400 | 1920 | 0.176/0.262 | 0.191/0.292 | 0.235/0.330 |
| standard | origin time matched | strong | 9393 | 2400 | 1949 | 0.153/0.240 | 0.195/0.292 | 0.239/0.348 |
| standard | location only | pooled | 7171 | 4800 | 3826 | 0.367/0.474 | 0.385/0.494 | 0.423/0.530 |
| standard | location only | pooled | 9393 | 4800 | 3863 | 0.354/0.470 | 0.369/0.488 | 0.421/0.537 |
| standard | location only | weak | 7171 | 2400 | 1906 | 0.453/0.566 | 0.490/0.612 | 0.524/0.636 |
| standard | location only | weak | 9393 | 2400 | 1914 | 0.460/0.580 | 0.481/0.611 | 0.527/0.646 |
| standard | location only | strong | 7171 | 2400 | 1920 | 0.311/0.429 | 0.326/0.434 | 0.382/0.472 |
| standard | location only | strong | 9393 | 2400 | 1949 | 0.294/0.405 | 0.305/0.426 | 0.355/0.481 |
| dense | origin time matched | pooled | 7171 | 3000 | 11978 | 0.050/0.102 | 0.066/0.124 | 0.092/0.162 |
| dense | origin time matched | pooled | 9393 | 3000 | 12000 | 0.056/0.109 | 0.065/0.120 | 0.084/0.151 |
| dense | origin time matched | weak | 7171 | 1500 | 5982 | 0.053/0.115 | 0.077/0.144 | 0.116/0.196 |
| dense | origin time matched | weak | 9393 | 1500 | 5963 | 0.068/0.126 | 0.075/0.136 | 0.089/0.169 |
| dense | origin time matched | strong | 7171 | 1500 | 5996 | 0.045/0.089 | 0.057/0.105 | 0.074/0.133 |
| dense | origin time matched | strong | 9393 | 1500 | 6037 | 0.047/0.091 | 0.055/0.107 | 0.078/0.134 |
| dense | location only | pooled | 7171 | 3000 | 11978 | 0.153/0.252 | 0.179/0.276 | 0.199/0.296 |
| dense | location only | pooled | 9393 | 3000 | 12000 | 0.159/0.255 | 0.177/0.272 | 0.200/0.307 |
| dense | location only | weak | 7171 | 1500 | 5982 | 0.181/0.286 | 0.220/0.323 | 0.240/0.342 |
| dense | location only | weak | 9393 | 1500 | 5963 | 0.187/0.284 | 0.215/0.321 | 0.241/0.340 |
| dense | location only | strong | 7171 | 1500 | 5996 | 0.125/0.218 | 0.145/0.232 | 0.173/0.258 |
| dense | location only | strong | 9393 | 1500 | 6037 | 0.134/0.233 | 0.144/0.243 | 0.171/0.274 |

**Table Y.80. Rules 2L-a to 2L-d and the registered secondary: MT-FW with S2 and T3 minus MT-FW with S0
and T3 in recall at $\le0.25$, point [95% hierarchical interval], with the one-sided lower bound that
each rule uses.** Same output, sections of the rules.

| Rule | Windows, hit rule, regime | Seed 7171 | Seed 9393 | Lower bounds used | Criterion | Passes |
|:--|:--|:--|:--|:--|:--|:--|
| 2L-a | standard, origin time matched, pooled | +0.057 [+0.034, +0.081] | +0.062 [+0.036, +0.083] | +0.0338, +0.0356 | bound $>0$ and point $\ge0.03$ | yes, both seeds |
| 2L-b | dense, origin time matched, pooled | +0.038 [+0.024, +0.054] | +0.030 [+0.015, +0.044] | +0.0236, +0.0153 | bound $>0$ | yes, both seeds |
| 2L-c | standard, origin time matched, weak | +0.083 [+0.045, +0.116] | +0.077 [+0.041, +0.117] | +0.0517, +0.0461 | bound $>-0.02$ | yes, both seeds |
| 2L-c | standard, origin time matched, strong | +0.038 [+0.002, +0.072] | +0.056 [+0.019, +0.088] | +0.0069, +0.0254 | bound $>-0.02$ | yes, both seeds |
| 2L-c | dense, origin time matched, weak | +0.052 [+0.028, +0.076] | +0.033 [+0.005, +0.057] | +0.0315, +0.0102 | bound $>-0.02$ | yes, both seeds |
| 2L-c | dense, origin time matched, strong | +0.029 [+0.005, +0.049] | +0.027 [+0.008, +0.048] | +0.0106, +0.0113 | bound $>-0.02$ | yes, both seeds |
| 2L-d | standard, location only, pooled | +0.036 [+0.014, +0.056] | +0.050 [+0.027, +0.069] | +0.0176, +0.0293 | bound $>-0.02$ | yes, both seeds |
| 2L-d | dense, location only, pooled | +0.020 [+0.002, +0.039] | +0.035 [+0.013, +0.049] | +0.0048, +0.0166 | bound $>-0.02$ | yes, both seeds |
| secondary | BOMP-PK, standard, origin time matched, pooled, 60 networks | +0.068 [+0.026, +0.107] | +0.071 [+0.029, +0.114] | +0.0262, +0.0287 | bound $>0$ | yes, both seeds |

**Rules 2L-a to 2L-d pass in both seeds.** With the detector fixed and only the terms changed, recall in
standard windows with origin times matched rose from 0.346 to 0.403 in seed 7171 and from 0.350 to
0.411 in seed 9393 at $\le0.25$, differences of $+0.057$ and $+0.062$ with one-sided $p$-values at the
bootstrap floor of 0.0005 (Tables Y.79 and Y.80). In dense windows it rose from 0.124 to 0.162 and from
0.120 to 0.151, $+0.038$ and $+0.030$. The lower bounds per regime lie between $+0.007$ and $+0.052$, so
that neither clutter regime pays for the other, and under the location-only rule the differences are
$+0.036$ and $+0.050$ in standard and $+0.020$ and $+0.035$ in dense windows, with every bound above
the margin of $-0.02$. In block pursuit, on 60 networks per cell and regime, the same contrast is
$+0.068$ and $+0.071$. By the registered interpretation the contribution is confirmed in this synthetic
design: calibration together with refinement, for detection that is consistent in origin time. The
confirmatory effect is smaller than the pilot estimate of $+0.070$ in standard windows and about half
of the pilot estimate of $+0.065$ in dense ones, the shrinkage on which the registration had planned.

Three descriptive readings place the result. Against the registered detector of the earlier tasks, S0
with the hard rule, the total change is $+0.088$ [$+0.063$, $+0.116$] and $+0.120$ [$+0.093$, $+0.146$]
in standard and $+0.060$ and $+0.042$ in dense windows, of which the refinement with the registered
terms supplies $+0.030$ [$+0.010$, $+0.055$] and $+0.058$ [$+0.038$, $+0.083$] in standard and $+0.022$
and $+0.012$ in dense windows; the refinement alone, which was a pilot signal in Y.27, is therefore
positive on these seeds in standard windows, although no rule was registered for it. The hits are
better placed: the median epicentral error of the hit atoms falls from 2.05 to 1.81 km and from 2.06 to
1.75 km in standard windows and the median origin-time error from 0.325 to 0.282 s and from 0.315 to
0.280 s, while the share of events touched by some atom moves only from 0.742 to 0.763 and from 0.749 to
0.771. And the terms themselves improve as in the pilot (Table Y.81): the median error falls from 0.147
and 0.151 s to 0.082 and 0.084 s with weak clutter and from 0.216 and 0.207 s to 0.155 and 0.153 s with
strong clutter, S2 being better than S0 in 0.85 and 0.84 of the networks with weak clutter and in 0.70
and 0.73 with strong clutter.

**Table Y.81. Task 2L, descriptive. Rms error of the station terms against the true delays with the
median removed, in seconds: median over the 160 networks of each cell and median [quartiles] and mean
over the 480 networks of the regime; share of networks in which S2 has the smaller error; median
absolute common offset; and median timing level of T3 computed with each set of terms.** Output
`analyze_2l_out.txt`, section of the station terms. The timing level of the oracle arm is 0.049 to
0.050 s.

| Seed | Regime | Terms | 4 stations | 6 stations | one-sided, 8 st. | all 480 networks, median [quartiles] | mean | share better than S0 | common offset | T3 timing level |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| 7171 | weak | S0 | 0.122 | 0.151 | 0.163 | 0.147 [0.119, 0.181] | 0.151 |  | 0.044 | 0.141 |
| 7171 | weak | S2 | 0.092 | 0.086 | 0.064 | 0.082 [0.057, 0.112] | 0.092 | 0.85 | 0.050 | 0.131 |
| 7171 | strong | S0 | 0.186 | 0.215 | 0.241 | 0.216 [0.166, 0.268] | 0.221 |  | 0.050 | 0.227 |
| 7171 | strong | S2 | 0.157 | 0.155 | 0.151 | 0.155 [0.107, 0.246] | 0.185 | 0.70 | 0.055 | 0.209 |
| 9393 | weak | S0 | 0.128 | 0.158 | 0.164 | 0.151 [0.120, 0.185] | 0.155 |  | 0.045 | 0.144 |
| 9393 | weak | S2 | 0.095 | 0.096 | 0.060 | 0.084 [0.058, 0.116] | 0.094 | 0.84 | 0.045 | 0.134 |
| 9393 | strong | S0 | 0.192 | 0.206 | 0.221 | 0.207 [0.165, 0.266] | 0.221 |  | 0.042 | 0.220 |
| 9393 | strong | S2 | 0.167 | 0.157 | 0.132 | 0.153 [0.103, 0.239] | 0.181 | 0.73 | 0.050 | 0.205 |

**Where the effect is.** The pooled gain is carried by one geometry (Table Y.82). In standard windows
with origin times matched the contrast is $+0.108$ [$+0.057$, $+0.151$] and $+0.113$ [$+0.067$,
$+0.157$] in the one-sided networks of eight stations, $+0.048$ [$+0.005$, $+0.089$] and $+0.027$
[$-0.023$, $+0.069$] with six stations and $+0.021$ [$-0.013$, $+0.055$] and $+0.018$ [$-0.016$,
$+0.054$] with four. In dense windows it is $+0.066$ and $+0.054$ in the one-sided cell, $+0.027$ and
$+0.026$ with six stations, both intervals reaching zero, and $+0.000$ and $-0.005$ with four. With
strong clutter the four- and six-station cells show nothing in standard windows, $-0.002$ and $+0.022$
with four stations and $-0.002$ and $+0.014$ with six, while the one-sided cell holds $+0.113$ and
$+0.106$. The one-sided cell is also where S2 improves the terms most, from 0.163 and 0.164 s to 0.064
and 0.060 s with weak clutter. No rule was registered per cell, and the per-cell intervals are
descriptive.

**Table Y.82. Task 2L, descriptive. Recall at $\le0.25$ per cell with origin times matched, and
differences, point [95% hierarchical interval]: the contrast of the rules, S2 and T3 minus S0 and T3,
and the total against the registered detector, S2 and T3 minus S0 with the hard rule; 1600 standard or
1000 dense windows per pooled row and half as many per regime.** Same output, per-cell table.

| Windows | Seed | Cell | Regime | events | MT-FW, S0 | S0 and T3 | S2 and T3 | S2 and T3 minus S0 and T3 | S2 and T3 minus S0 |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| standard | 7171 | 4 stations | pooled | 1257 | 0.200 | 0.215 | 0.236 | +0.021 [-0.013, +0.055] | +0.036 [+0.004, +0.069] |
| standard | 7171 | 4 stations | weak | 630 | 0.311 | 0.337 | 0.398 | +0.062 [+0.002, +0.126] | +0.087 [+0.031, +0.154] |
| standard | 7171 | 4 stations | strong | 627 | 0.144 | 0.166 | 0.164 | -0.002 [-0.040, +0.042] | +0.021 [-0.020, +0.071] |
| standard | 7171 | 6 stations | pooled | 1275 | 0.360 | 0.413 | 0.460 | +0.048 [+0.005, +0.089] | +0.100 [+0.047, +0.149] |
| standard | 7171 | 6 stations | weak | 635 | 0.436 | 0.496 | 0.561 | +0.065 [-0.002, +0.126] | +0.124 [+0.043, +0.190] |
| standard | 7171 | 6 stations | strong | 640 | 0.316 | 0.359 | 0.358 | -0.002 [-0.053, +0.074] | +0.042 [-0.012, +0.125] |
| standard | 7171 | one-sided, 8 st. | pooled | 1294 | 0.379 | 0.415 | 0.523 | +0.108 [+0.057, +0.151] | +0.145 [+0.091, +0.197] |
| standard | 7171 | one-sided, 8 st. | weak | 641 | 0.426 | 0.480 | 0.608 | +0.128 [+0.057, +0.179] | +0.183 [+0.119, +0.245] |
| standard | 7171 | one-sided, 8 st. | strong | 653 | 0.337 | 0.351 | 0.464 | +0.113 [+0.040, +0.184] | +0.127 [+0.057, +0.193] |
| standard | 9393 | 4 stations | pooled | 1282 | 0.211 | 0.236 | 0.254 | +0.018 [-0.016, +0.054] | +0.043 [+0.005, +0.081] |
| standard | 9393 | 4 stations | weak | 647 | 0.334 | 0.372 | 0.410 | +0.037 [-0.031, +0.096] | +0.076 [-0.003, +0.143] |
| standard | 9393 | 4 stations | strong | 635 | 0.145 | 0.176 | 0.198 | +0.022 [-0.020, +0.065] | +0.054 [+0.007, +0.098] |
| standard | 9393 | 6 stations | pooled | 1278 | 0.337 | 0.416 | 0.444 | +0.027 [-0.023, +0.069] | +0.106 [+0.060, +0.157] |
| standard | 9393 | 6 stations | weak | 652 | 0.382 | 0.485 | 0.531 | +0.046 [-0.033, +0.126] | +0.149 [+0.071, +0.230] |
| standard | 9393 | 6 stations | strong | 626 | 0.294 | 0.379 | 0.393 | +0.014 [-0.064, +0.078] | +0.099 [+0.026, +0.165] |
| standard | 9393 | one-sided, 8 st. | pooled | 1303 | 0.351 | 0.421 | 0.533 | +0.113 [+0.067, +0.157] | +0.182 [+0.130, +0.230] |
| standard | 9393 | one-sided, 8 st. | weak | 615 | 0.423 | 0.499 | 0.641 | +0.141 [+0.076, +0.201] | +0.218 [+0.143, +0.296] |
| standard | 9393 | one-sided, 8 st. | strong | 688 | 0.298 | 0.352 | 0.458 | +0.106 [+0.043, +0.166] | +0.160 [+0.095, +0.230] |
| dense | 7171 | 4 stations | pooled | 3983 | 0.085 | 0.094 | 0.095 | +0.000 [-0.019, +0.022] | +0.010 [-0.008, +0.030] |
| dense | 7171 | 4 stations | weak | 1975 | 0.100 | 0.122 | 0.129 | +0.008 [-0.027, +0.041] | +0.029 [-0.005, +0.059] |
| dense | 7171 | 4 stations | strong | 2008 | 0.067 | 0.073 | 0.072 | -0.001 [-0.024, +0.019] | +0.005 [-0.017, +0.025] |
| dense | 7171 | 6 stations | pooled | 3984 | 0.122 | 0.137 | 0.164 | +0.027 [-0.004, +0.055] | +0.042 [+0.012, +0.068] |
| dense | 7171 | 6 stations | weak | 1982 | 0.142 | 0.149 | 0.176 | +0.027 [-0.032, +0.075] | +0.034 [-0.021, +0.085] |
| dense | 7171 | 6 stations | strong | 2002 | 0.106 | 0.125 | 0.151 | +0.026 [-0.005, +0.058] | +0.045 [+0.014, +0.076] |
| dense | 7171 | one-sided, 8 st. | pooled | 4011 | 0.113 | 0.147 | 0.213 | +0.066 [+0.037, +0.099] | +0.099 [+0.072, +0.131] |
| dense | 7171 | one-sided, 8 st. | weak | 2025 | 0.123 | 0.172 | 0.267 | +0.095 [+0.041, +0.140] | +0.144 [+0.094, +0.189] |
| dense | 7171 | one-sided, 8 st. | strong | 1986 | 0.096 | 0.119 | 0.167 | +0.048 [+0.014, +0.086] | +0.071 [+0.035, +0.102] |
| dense | 9393 | 4 stations | pooled | 3981 | 0.083 | 0.092 | 0.088 | -0.005 [-0.020, +0.017] | +0.004 [-0.013, +0.023] |
| dense | 9393 | 4 stations | weak | 1982 | 0.102 | 0.127 | 0.109 | -0.018 [-0.053, +0.021] | +0.007 [-0.025, +0.043] |
| dense | 9393 | 4 stations | strong | 1999 | 0.068 | 0.067 | 0.074 | +0.007 [-0.014, +0.030] | +0.007 [-0.016, +0.031] |
| dense | 9393 | 6 stations | pooled | 4011 | 0.126 | 0.131 | 0.158 | +0.026 [-0.005, +0.058] | +0.032 [+0.001, +0.063] |
| dense | 9393 | 6 stations | weak | 1999 | 0.139 | 0.148 | 0.188 | +0.041 [-0.011, +0.085] | +0.050 [-0.002, +0.091] |
| dense | 9393 | 6 stations | strong | 2012 | 0.114 | 0.120 | 0.138 | +0.018 [-0.019, +0.056] | +0.024 [-0.012, +0.061] |
| dense | 9393 | one-sided, 8 st. | pooled | 4008 | 0.127 | 0.146 | 0.200 | +0.054 [+0.025, +0.085] | +0.073 [+0.039, +0.104] |
| dense | 9393 | one-sided, 8 st. | weak | 1982 | 0.152 | 0.160 | 0.209 | +0.049 [+0.011, +0.101] | +0.057 [+0.013, +0.113] |
| dense | 9393 | one-sided, 8 st. | strong | 2026 | 0.104 | 0.130 | 0.189 | +0.059 [+0.016, +0.098] | +0.085 [+0.039, +0.127] |

**The oracle arm.** On the 30 networks per cell and regime that carry it, the true delays with T3 exceed
the S0 terms with T3 by $+0.116$ and $+0.152$ in standard and by $+0.080$ and $+0.075$ in dense windows
(Table Y.83). S2 captures 0.56 and 0.61 of that gain in standard windows and 0.55 and 0.28 in dense
ones, with wide intervals, and the true delays remain ahead of S2 by $+0.052$ [$+0.006$, $+0.102$] and
$+0.060$ [$-0.003$, $+0.119$] in standard and by $+0.036$ [$+0.004$, $+0.062$] and $+0.054$ [$+0.024$,
$+0.082$] in dense windows. The pilot fraction of 0.77 was therefore optimistic, and part of what
calibration could give is not reached by this estimator with twenty events. The oracle arm uses a
timing level of about 0.05 s in T3 and the estimated terms 0.13 to 0.23 s, so the oracle confounds the
terms with the level.

**Table Y.83. Task 2L, descriptive: the oracle arm, 30 networks per cell and regime, 900 windows per row;
recall at $\le0.25$ with origin times matched, differences with 95% hierarchical intervals, and the
fraction of the oracle gain captured, (S2 minus S0) over (true delays minus S0), all with T3, with the
interval of the ratio.** Same output, section of the oracle arm.

| Windows | Seed | Family | events | S0 and T3 | S2 and T3 | true delays and T3 | S2 minus S0 | true delays minus S0 | true delays minus S2 | Fraction of the oracle gain |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| standard | 7171 | MT-FW | 697 | 0.333 | 0.397 | 0.449 | +0.065 [+0.003, +0.122] | +0.116 [+0.054, +0.178] | +0.052 [+0.006, +0.102] | 0.56 [+0.04, +0.94] |
| standard | 7171 | BOMP-PK | 697 | 0.317 | 0.397 | 0.452 | +0.080 [+0.019, +0.141] | +0.135 [+0.081, +0.199] | +0.055 [+0.008, +0.112] | 0.60 [+0.17, +0.93] |
| standard | 9393 | MT-FW | 716 | 0.330 | 0.422 | 0.482 | +0.092 [+0.027, +0.162] | +0.152 [+0.094, +0.217] | +0.060 [-0.003, +0.119] | 0.61 [+0.21, +1.03] |
| standard | 9393 | BOMP-PK | 716 | 0.362 | 0.409 | 0.497 | +0.047 [-0.009, +0.099] | +0.135 [+0.070, +0.198] | +0.088 [+0.037, +0.142] | 0.35 [-0.10, +0.68] |
| dense | 7171 | MT-FW | 3622 | 0.119 | 0.163 | 0.199 | +0.044 [+0.014, +0.071] | +0.080 [+0.051, +0.104] | +0.036 [+0.004, +0.062] | 0.55 [+0.22, +0.94] |
| dense | 9393 | MT-FW | 3606 | 0.120 | 0.141 | 0.195 | +0.021 [-0.006, +0.049] | +0.075 [+0.044, +0.101] | +0.054 [+0.024, +0.082] | 0.28 [-0.11, +0.63] |

**Independent verification.** An adversarial verification looked for a false positive
(`verif_2l/task2l_verification.txt`), with no change to any file of the task or of the registry, and its
verdict is a confirmed positive within a stated scope. It addressed five questions.

On leakage, it read the code path of the estimator and found that S2 receives only the stations, the
catalogue, which is the list of waveforms with catalogue positions and origins, the S0 terms and the
constants; the true delays and the true hypocentres are read only by the functions that store quality
figures, by the oracle arm and by the generator. It then tested this empirically on 300 units of
network and regime, networks 100 to 124 of both seeds, in fresh processes: S2 recomputed on copies of
its inputs with the true delays and the true catalogue removed from the network object, and the
energies recomputed from the catalogue waveforms, was bit-identical to the stored terms in 300 of 300.
The catalogue given to the estimators carries the registered errors, 0.998 km per coordinate and 0.099
s, with none of 18000 coordinates equal to the truth. Catalogue and evaluation windows come from
disjoint random streams and share only the network, which is the object being calibrated, and the terms
are fixed before the window stream of a network is opened.

On the origin time, the concern was that a hit rule with a window of 0.5 s could reward one arm for a
systematic shift of its origins and not for better detection. The signed origin error of the hits has
the same small bias in both arms, medians of $-0.036$ and $-0.035$ s with S0 against $-0.045$ and
$-0.038$ s with S2 in standard windows; S2 does not recentre the origins, it narrows them, the share of
hits within 0.3 s rising from 0.47 to 0.53 in standard and from 0.42 to 0.46 in dense windows. The
common offset of the terms against the true delays is the same in both arms, with a standard deviation
of 0.07 s over networks. The contrast is positive at every hit window tested, from 0.2 to 1.0 s and with
no origin condition at all (Table Y.84), between $+0.037$ and $+0.065$ in standard windows for windows
of 0.3 to 1.0 s, and smaller in dense windows at 0.2 and 0.3 s. With the median origin error of each
method removed before the rule is applied, the contrast at 0.5 s is $+0.060$ and $+0.062$ in standard
and $+0.038$ and $+0.036$ in dense windows, and removing a median per method, regime and cell, or per
method and network with the help of the truth, changes it little. The gain does not come from hits
moved across the 0.5 s line by an offset.

**Table Y.84. Verification of Task 2L, exploratory: the contrast of the rules at $\le0.25$ under other
origin-time windows of the hit rule, seeds 7171 / 9393, with the labelling code of the verification.**
Output `verif_2l/task2l_verification.txt`, section 1.7.

| Windows | 0.2 s | 0.3 s | 0.5 s, registered | 0.75 s | 1.0 s | epicentre only |
|:--|:--|:--|:--|:--|:--|:--|
| standard | +0.050 / +0.026 | +0.065 / +0.046 | +0.057 / +0.062 | +0.047 / +0.054 | +0.037 / +0.052 | +0.036 / +0.050 |
| dense | +0.010 / +0.010 | +0.017 / +0.020 | +0.038 / +0.031 | +0.030 / +0.035 | +0.031 / +0.028 | +0.020 / +0.035 |

On attribution, the two arms run on the same window, lattices, origin times, number of atoms, exclusion,
cap, selection rule and refinement radius, and one registered difference accompanies the terms: the
timing level of T3 is recomputed with each set of terms, and it is lower with S2 in 0.74 of the
one-sided networks and in about 0.53 of the others. A run on 1500 standard windows of the same 25
networks, both seeds pooled, separated the two. With the registered contrast at $+0.036$ [$-0.006$,
$+0.075$] on that subset, the terms alone, with the S0 level in both arms, gave $+0.045$ [$+0.005$,
$+0.089$], and the level alone gave $+0.000$ [$-0.018$, $+0.019$] with the S0 terms and $-0.009$
[$-0.029$, $+0.008$] with the S2 terms. The gain is carried by the terms and the change of level
contributes nothing. The subset is small, and on it the registered contrast is $+0.077$ in seed 7171
and $-0.003$ in seed 9393; it serves the attribution only.

On provenance and regeneration, the code hashes check, 37 of 37, and every module that the runner loads
from the work directory is in the list; the hashes of the registry are those of successive prefixes of
an append-only file, and the registration hash is that of the prefix of 653 lines that ends with the
last line of the protocol of Task 2L; no part of seeds 7171 or 9393 predates the freeze; a search of the
work directory finds no earlier use of either seed; and the constants are those of the frozen plan of
C11, whose body and code hashes check. The clock of the machine is the only time source and the order
of events is internally consistent. The record of the freeze lists relaunches at 22:48:13 and 22:49:32
UTC where the log has one start line, at 22:49:51 UTC, and no part was written between 22:47:19 and
22:50:35 UTC, which is harmless. Six parts, two written before the first interruption, one between
interruptions and three after, were regenerated from scratch, each in a fresh process, and were
identical to the stored parts in terms, levels, truth and 1530 atoms, timings apart; the atoms of MT-FW
with S2 and T3 were also reproduced in the 300 units of the leakage run. A rerun of the frozen analysis
gave an output identical byte for byte.

On the statistics, an own loader, labelling, exact FROC and hierarchical bootstrap, with no code of the
analysis, reproduced the points of the rules to four decimals, $+0.0572$ and $+0.0616$ for 2L-a and
$+0.0380$ and $+0.0305$ for 2L-b, and the intervals, [$+0.032$, $+0.082$] and [$+0.037$, $+0.083$] for
2L-a against the registered [$+0.034$, $+0.081$] and [$+0.036$, $+0.083$], with no replicate at or below
zero; rules 2L-b and 2L-d were re-derived in the same way and pass in both seeds. The operating point is
chosen per method on the evaluation windows, as registered and as in the earlier tasks, and in the same
way for both arms.

**Scope findings of the verification.** The exploratory checks of the verification, on the stored parts
only, bound what the result means.

The effect is concentrated. The verification reproduces the per-cell values of Table Y.82 with its own
bootstrap, $+0.108$ [$+0.053$, $+0.155$] and $+0.113$ [$+0.061$, $+0.163$] in the one-sided cell of
eight stations, and states the consequence: the pooled effect is essentially that cell, and with four
stations no gain is demonstrated in either seed.

The estimator can make the terms worse. S2 has a larger error than S0 in about 22 per cent of the
networks overall: 3 to 4 per cent of the one-sided networks with weak clutter, 19 to 27 per cent of
those with strong clutter, 12 to 28 per cent with six stations and 26 to 37 per cent with four. In
those networks the change of recall at the pooled operating points is $+0.019$ and $-0.001$ in standard
and $+0.026$ and $+0.025$ in dense windows, so no pooled loss appears; but in the six-station cell, in
standard windows, it is negative, $-0.019$ and $-0.049$. Where S2 is better the gain is $+0.068$ and
$+0.079$, and by quartile of the reduction of the term error it rises from $+0.020$ and $-0.003$, where
the terms are worse by 0.05 s, to $+0.108$ and $+0.117$, where they are better by 0.14 s. The gain
follows the improvement of the calibration, and a local loss exists in six-station networks whose terms
get worse.

The claim is conditional on the refinement. Task 2L has no arm with the S2 terms and the hard rule, and
in the pilot the terms did nothing measurable without T3. What was measured is that the
origin-time-consistent recall of the refined detector depends on the accuracy of the terms, not that
station terms improve the detector of the earlier tasks.

The registered superiority is for the origin-time hit rule. Under the location-only rule only
non-inferiority was registered, and the gain observed there, $+0.036$ and $+0.050$ in standard windows,
is descriptive.

The estimator assumes catalogue errors of known size. Its prior is Gaussian with the true widths of the
catalogue errors of the benchmark, 1 km per coordinate and 0.1 s, and on the development set of C11 a
prior of 2 km brought the term error back to the level of S0. The result is not shown for larger,
biased or unknown catalogue errors, nor for catalogues of fewer than twenty events.

The evidence is synthetic and from one generator family: the generator of the earlier tasks, station
delays with a standard deviation of 0.3 s fixed per network, three geometries and two clutter regimes.
The two seeds are independent draws of that family and not independent evidence about real data.

The size that the verification recommends quoting is about $+0.06$ of recall in standard and $+0.03$
in dense windows at $\le0.25$ false events per window, concentrated in the one-sided geometry of eight
stations, with no demonstrated gain for four stations, and with only part of the oracle gain captured,
0.28 to 0.61.

**Novelty.** Joint hypocentre determination with station terms is classical in location from picks,
source-specific station terms extend it, and static corrections by iterated relocation, robust
residual statistics and catalogue priors in the relocation are standard. S2 is a small joint
determination with a prior, and self-calibration from one's own detections is the usual practice of
template matching and of iterative association and relocation. The registration declares this. What is
not classical, and what Task 2L measures, is the estimation of the terms without picks, from
matched-filter energy peaks that are read again as hypocentres and terms move, and their use inside a
pick-free moment-tensor subspace detector, with a measured effect on detection at a fixed false-event
rate reported next to the bound given by the true delays. Two findings accompany it: the effect appears
only together with a local refinement by the timing likelihood and under a hit rule that asks for the
origin time, and self-calibration from detections adds nothing to twenty catalogued events at the
purity and location accuracy of this detector. The targeted searches of Appendix E were not repeated
for this element, and absence of a precedent is not claimed. It is a calibration result about what
limits a detector, not a detector.

**What Task 2L does not test.** Real data; a generator family other than `gen_realistic`; station delays
of another size; catalogues with errors other than 1 km and 0.1 s, with biased errors or with errors of
unknown size; catalogues of other than twenty events on the evaluation seeds; the S2 terms without the
refinement; the Huber statistic of Y.26, since both arms use MT-FW with the cap of 6; a held-out
operating point, the thresholds being chosen on the evaluation windows, per method and in the same way
for both arms; windows with more than five events; and superiority under the location-only rule.

### Y.28a What Task 2L changes

Task 2L turns the first consequence of Y.27a into a preregistered result. Y.27 located part of the
headroom in the accuracy of the station terms on the evidence of an oracle, and C11 found an estimator
that recovers part of that accuracy from the same twenty catalogued events. Task 2L shows, on two seeds
that neither study touched, with the registration hash sent to the user before any run, the code frozen
by hashes and an independent verification that found no leakage, that those terms raise the recall of
the refined moment-tensor detector at a fixed false-event rate when a hit must also have the right
origin time, by $+0.057$ and $+0.062$ in standard and by $+0.038$ and $+0.030$ in dense windows, with no
loss in either clutter regime and none under the location-only rule. It is the calibration hypothesis
G1-C of the main document, and it stands beside the bounded-influence hypothesis G1-R of Y.24 and Y.26
as the second methodological result with preregistered support. The two act on different parts of the
headroom of Y.22, G1-R on the amplitude error of single station and phase blocks and G1-C on the timing
error of the station terms, and they have not been run together: both arms of Task 2L use the cap of 6,
and neither Task 2F nor Task 2I used the refinement or the S2 terms.

The result is modest and it is narrower than its pooled rules suggest. It is essentially the result of
one geometry, the one-sided network of eight stations, where the registered corrections are worst
relative to what the catalogue allows; with six stations it is small and unresolved in one seed, and
with four stations nothing is demonstrated. It is conditional on the refinement and on the origin-time
hit rule. It rests on a catalogue prior whose widths are the true catalogue errors. In about a fifth of
the networks the estimator makes the terms worse, and in six-station networks of that kind recall
falls. And it captures 0.28 to 0.61 of what the true delays would give, so the calibration part of the
headroom is not closed by it: the oracle remains $+0.052$ and $+0.060$ ahead in standard and $+0.036$
and $+0.054$ ahead in dense windows. Absolute recall with origin times matched stays low, 0.403 and
0.411 in standard and 0.162 and 0.151 in dense windows at $\le0.25$.

Two lines of Y.27a are settled by C11 on the calibration seed, as exploratory negatives. Station terms
estimated self-consistently with the detections, the form in which Y.27a stated the lever, add nothing
measurable to the joint determination on the catalogue alone; and a catalogue of forty or eighty events
lowers the term error further without a resolved gain in recall. What Y.27a asked of the protocol is in
place: the origin-time hit rule is the registered primary rule of Task 2L, equivalence gave way to
non-inferiority per regime with a margin of 0.02, and the oracle ran beside the estimator.

The next steps follow from the scope. The sensitivity to the catalogue comes first, with a prior that
is wrong in width, a biased catalogue and fewer events, because the assumption of known catalogue
errors is the one a real network does not grant. The dependence on geometry comes second, with a
criterion that says before detection whether a network is one in which the estimate can be trusted.
Third, G1-R and G1-C are to be combined in one registered test, the Huber or capped statistic with and
without the S2 terms and the refinement, in dense windows with the origin-time rule, since their
mechanisms differ and their gains may or may not add. Real data follow these, with residual station
terms reported per station and no oracle to read them against.

What Y.18a to Y.26a list as not established remains so: no deployable frozen operating point exists,
GaMMA with optimistic picks is ahead at 0.1 with weak clutter, the sources are double couples inside the
detector's subspace, the generator lies outside the declared regime of Section 1.5, and there is no
result on real data.

## Y.29 After Task 2L: three exploratory studies on the calibration seed (C12, C13, C14)

Code: `explore_contrib/c12_lib.py` (the full, deviatoric and double-couple statistics, the optimizer over
the double-couple manifold, the tube-formula constants), `c12_run.py`, `c12_analyze.py`,
`c12_validate.py` and `c12_mech.py`, report `results/explorations/REPORT_C12.txt`, outputs
`results/explorations/out/c12_analysis_STD.txt`, `c12_analysis_DENSE.txt`, `c12_null.txt`,
`c12_validation.txt` and `c12_mech.txt`; `c13_lib.py` (the decoys and the calibration rules),
`c13_run.py`, `c13_run_d3c.py`, `c13_twin.py`, `c13_analyze.py` and `c13_tables.py`, report
`REPORT_C13.txt`, outputs `out/c13_analysis.txt` and `out/c13_report_tables.txt`; `c14_lib.py` (the
detector loop with the remedies, the event-free twins and the exact decomposition of the statistic),
`c14_run.py` and `c14_analyze.py`, with the frozen configuration `c14_frozen.json` and the hashes of the
freeze in `C14_FREEZE.txt` and `C14_FREEZE.txt.post`, report `REPORT_C14.txt`, outputs
`out/c14_analysis_dev.txt`, `c14_analysis_dev2.txt` and `c14_analysis_eval.txt`. All three start from the
detector loop of `explore_contrib/c10_lib.py`, which their own loops reproduce atom by atom when their
additions are switched off.

**Status.** The three studies are exploratory and test no registered hypothesis. They ran on seed 999
only, on networks of that seed that no earlier design had used (1200 to 1229 and 1290 for C12, 1400 to
1411 for C13, 1600 to 1602 and 1610 to 1621 for C14), with random streams of their own, after Task 2L had
been run and verified: C12 and C13 on 5 October 2026, sharing the two cores of the machine, and C14
between 5 and 9 October. C11, the exploration that led to Task 2L, is reported in Y.28. In all three the
detector is MT-FW with the cap of 6, the registered corrections S0 and the hard peak rule of $\pm0.25$ s,
the registered detector of the earlier tasks and not the refined detector of Task 2L; the refinement T3
appears only as a reference arm of C14. The protocol is that of Y.27 and Y.28: the exact FROC, recall at
$\le0.1$ and $\le0.25$ realized false events per window with the operating point chosen again in each
bootstrap replicate, the hit rule with origin times matched (epicentre within 5 km and origin-time error
of at most 0.5 s) as the primary reading and the location-only rule as the second, hierarchical bootstrap
intervals of 95 per cent, and no multiplicity correction. None of the three proposes a candidate for
confirmation as a contribution to recall. Their value for the thesis is an analytic piece (C12), a
partial label-free calibration (C13) and a negative with a diagnosis (C14), stated below with their
limits.

**C12: the question.** The moment-tensor subspace statistic of MT-FW maximizes over the six-dimensional
subspace of all moment tensors, while the sources of the benchmark are double couples. C12 asked whether
restricting the statistic to the double-couple manifold lowers the threshold at a fixed false-alarm rate
enough to detect weak isolated events, and whether the null law of the restricted statistic can be
derived.

**C12: design.** Networks 1200 to 1215 of each cell and regime, 16 per cell and regime, carried 20
windows each, 1920 windows per design, with 1550 events in the standard windows and 7613 in the dense
ones; 1440 event-free windows with the clutter of the regimes and 600 noise-only windows calibrated the
null. Four statistics shared the greedy loop: the registered one (FULL6), its projection on the traceless
hyperplane (DEV5), the maximum over unit double couples computed exactly over the whole candidate grid by
a bound that orders the three statistics (DC), and the atoms of FULL6 scored with the double-couple
statistic (MIXED). The deflation used the full six-parameter fit in every variant, so that only the
selection and the score differ. No constant was changed after the evaluation runs started.

**C12: the null law.** For one candidate with a fixed read and Gaussian noise, twice the double-couple
statistic is the maximum of a Gaussian field over a closed three-dimensional manifold of the five-sphere,
and its tail follows the tube formula with two constants: the volume of the projected manifold, $L_3$,
and the integral of its scalar curvature, $L_1$. In the isotropic metric $L_3$ is exactly twice the
volume of the three-sphere and $L_1=0$; in the Gram metric of the networks of the study $L_3$ is 1.75 to
1.89 times the volume of the three-sphere and $L_1$ lies between -0.3 and -0.6, which gives the statistic
4.4 to 4.5 effective degrees of freedom instead of 6. Simulation agrees: the upper 0.01 quantile of twice
the statistic is 14.28 simulated against 14.17 from the tube formula, and the 0.001 quantile 19.69
against 19.46, with a Monte Carlo error of 0.23. In the full detector on event-free windows, with the
registered peak read, the double-couple threshold at equal false events per window is 0.906 [0.883,
0.928] of the full one at $\le0.1$ and 0.898 [0.886, 0.914] at $\le0.25$ with noise only, against 0.888
and 0.883 from the theory at the implied per-candidate tail: about 10 per cent lower in energy, at most 5
to 6 per cent in amplitude if the signal kept all its energy in the manifold.

**C12: detection.** The lower threshold does not reach recall (Table Y.85). Against FULL6 at $\le0.25$,
the double-couple statistic gives $+0.005$ [$-0.021$, $+0.033$] in standard windows with origin times
matched and $-0.004$ [$-0.018$, $+0.011$] in dense ones, and no pooled contrast of any variant excludes
zero; of the 360 paired intervals computed, by regime, cell and amplitude tercile, 6 exclude zero, fewer
than chance would give at 5 per cent. Weak isolated events, the lowest amplitude tercile, stay almost all
missed with every statistic, recall 0.060 for FULL6 with origin times matched. The mechanism was measured
on the atoms of FULL6: in the band of scores where the decision is made, between 25 and 50, hits lose on
average 5.9 units when restricted to double couples and false atoms 6.6, while the realized threshold
drops by 5.1 at $\le0.1$ and 4.3 at $\le0.25$, and the probability that a hit outscores a false atom
moves from 0.893 to 0.879 in standard windows. An exact double couple read at its node loses 1 unit, as
the theory requires (0.95 measured); a post hoc check at the level of the peak vectors attributes 1.7 to
2.5 of the 5.9 to 6.4 units observed to the site factors and polarization scatter of the generator, and
the remaining 3 to 4 units were not decomposed. The exact double-couple statistic costs 0.199 against
0.135 s per standard window.

**Table Y.85. Exploration C12, networks 1200 to 1215 of seed 999, 1920 windows per design: recall of the registered statistic and paired differences against it at $\le0.25$ realized false events per window, point [95% hierarchical interval].** Report `results/explorations/REPORT_C12.txt`, section 0; outputs `out/c12_analysis_STD.txt` and `out/c12_analysis_DENSE.txt`.

| Statistic | standard, origin time matched | standard, location only | dense, origin time matched | dense, location only |
|:--|:--|:--|:--|:--|
| FULL6, registered, recall | 0.317 | 0.470 | 0.117 | 0.268 |
| DEV5 minus FULL6 | -0.005 [-0.027, +0.020] | -0.004 [-0.022, +0.026] | -0.006 [-0.016, +0.009] | -0.001 [-0.021, +0.016] |
| DC minus FULL6 | +0.005 [-0.021, +0.033] | +0.017 [-0.010, +0.043] | -0.004 [-0.018, +0.011] | +0.006 [-0.016, +0.023] |
| MIXED minus FULL6 | -0.003 [-0.020, +0.018] | -0.008 [-0.030, +0.018] | +0.004 [-0.007, +0.012] | +0.002 [-0.013, +0.015] |

**C12: verdict and limits.** Null for detection and positive for the analytic piece. The null law of a
double-couple test is derived and checked against simulation, with two constants computable per node, and
it bounds what any restriction of the source subspace can buy in this benchmark at about 5 per cent of
amplitude. It is not a lever for recall here, and no confirmation is proposed. The study used one seed
and 16 networks per cell and regime, so its intervals exclude gains above about 0.03 in standard and
0.015 in dense windows but not gains of 0.01; it did not use the refinement or the terms of Task 2L, with
which the signal might sit closer to the manifold; and the deflation was never done with the
double-couple fit.

**C13: the question.** No operating point of the benchmark is deployable, because every threshold so far
was chosen with the labels of the evaluation windows or calibrated on empty windows, which do not bound
false events in windows with events (Y.21). C13 asked whether a target-decoy calibration, a second pass
of the detector over hypotheses that cannot correspond to a real source, estimates the false detections
above a threshold without labels, so as to fix an operating point or control the false discovery rate.

**C13: design.** Networks 1400 to 1411, 12 per cell and regime, gave 1440 standard windows with 1178
events, 1440 dense windows with 5747 events, 576 windows with clutter and no events and 144 with noise
only. Six decoys were drawn afresh for each window: station moveouts permuted by a derangement of the
travel-time tables (D1p); moveouts with a circular offset per station and phase (D1o); templates with the
sign flipped in half of the station and phase blocks (D2s); polarization scrambled by a random rotation
per station (D2r); every trace shifted circularly by its own delay (D3); and, added post hoc, the same
with the delays restricted (D3c). Four calibrations were applied over batches of 20 or 120 windows: a
target-decoy search (tdc), conformal p-values with Benjamini-Hochberg, which gave the same numbers, a
paired competition by rank, and a decoy estimate of the false events per window (fpw), the label-free
counterpart of the operating point of the benchmark, added during the run before the sample was complete.
An event-free twin of each window, with the same clutter and noise sample by sample, separated what the
events generate from what the clutter does. The detector, the first five decoys, the rules, the levels,
the hit rules and the batches were fixed before the main run; D3c, the twins, the per-cell pooled lines
and a second form of the intervals were added after interim looks and are declared as such.

**C13: findings.** In windows without events the well-built decoys are exchangeable with the target, with
permutation p-values between 0.04 and 1.0 in the six strata and one of 48 below 0.05. D3 and D1o, as
fixed, were not, under strong clutter: a circular delay over the whole trace brings the clutter-free ends
of the trace into the samples that the hypotheses read, a construction defect that D3c, with restricted
delays, removes; D1o was not rerun. In windows with events the false atoms of the detector are made by
the events, not by the clutter: against the twin, the false atoms per window above the 0.95 quantile of
the clutter-only scores are 12 [9, 18] times as many in standard and 45 [33, 63] times in dense windows,
and after an event has been hit and deflated the later false atoms beyond 10 km are still 1.9 [1.6, 2.4]
and 5.0 [3.1, 8.8] times those of the twin at the 0.80 and 0.95 quantiles. A decoy is therefore useful
only if it keeps the energy of the events in the samples it reads, as D1p does. Under the location-only
rule the permuted-moveout decoy controls (Table Y.86): with batches of 120 windows the realized false
discovery rate is 0.14 to 0.16 at $q=0.2$ in the pooled lines, no stratum exceeds 0.2, and the one-sided
cell exceeds the nominal level at $q\le0.1$; and the decoy estimate of false events per window reproduces
the operating point of the benchmark without labels, with realized rates of 0.21 [0.17, 0.25] in standard
and 0.23 [0.19, 0.28] in dense windows for a nominal 0.25 and the recall of the labelled point within
0.03. Under the origin-time rule D1p fails, 0.28 to 0.37 at $q=0.2$, and the decoy estimate lets through
1.3 to 3 times the nominal rate; only the sign-flipped decoy controls the false discovery rate, at the
cost of almost all recall, and it keeps the false-event rate with a recall loss of 0.02 to 0.05. The
cause is a class of false atoms that no moveout decoy imitates: a real event located within 5 km with its
origin time off by more than 0.5 s, 0.2 to 1.2 such atoms per window with events. Two further readings
place the operating points of the benchmark: at $\le0.25$ false events per window 41 per cent of the
accepted atoms in standard and 20 per cent in dense windows are false under the location-only rule, and
52 and 37 per cent with origin times matched; and controlling the false discovery rate at 0.1 leaves a
recall of 0.195 in standard and 0.152 in dense windows, against 0.436 and 0.253 at the benchmark point.

**Table Y.86. Exploration C13, networks 1400 to 1411 of seed 999: calibration by the permuted-moveout decoy D1p with batches of 120 windows, under the two hit rules.** Report `results/explorations/REPORT_C13.txt`, sections 0, 3, 6 and 7; outputs `out/c13_analysis.txt` and `out/c13_report_tables.txt`. Ranges are over the pooled lines; intervals are 95 per cent percentile intervals of the hierarchical bootstrap.

| Quantity | location-only rule | origin-time rule |
|:--|:--|:--|
| tdc, realized false discovery rate at $q$ = 0.05 / 0.1 / 0.2 | 0.02-0.05 / 0.07-0.10 / 0.14-0.16 | 0.08-0.19 / 0.18-0.29 / 0.28-0.37 |
| fpw, realized false events per window for $c$ = 0.1 / 0.25, standard windows | 0.09 [0.07, 0.12] / 0.21 [0.17, 0.25] | 0.17 / 0.32 |
| fpw, realized false events per window for $c$ = 0.1 / 0.25, dense windows | 0.12 [0.09, 0.16] / 0.23 [0.19, 0.28] | 0.32 / 0.51 |
| fpw at $c$ = 0.25, recall against the labelled benchmark point, standard / dense | 0.453 / 0.245 against 0.436 / 0.253 | not valid |
| share of false atoms among those accepted at the benchmark point $\le0.25$, standard / dense | 0.41 / 0.20 | 0.52 / 0.37 |

**C13: verdict and limits.** Promising with conditions, and partial. A cheap label-free calibration
exists in this benchmark under the location-only rule, in the three cells and both clutter regimes, at a
cost of one more pass of the detector. It carries no guarantee, because exchangeability holds only
without events and its validity with events is empirical and geometry-dependent; it is not valid under
the origin-time rule, the primary rule of Task 2L; it falls short in the one-sided geometry at $q\le0.1$;
and the choice of D1p was made on the same data. One seed, 12 networks per cell and regime and 20 windows
per network, with the batches of 120 windows mixing six networks of a cell. The first version of the
event-free twins was wrong: copying the random generators lost their seed sequence, so that the twins
still held remnants of the events; it was detected by comparing the twins with empty windows of the same
network, corrected and rerun, and its outputs are kept apart and enter no number. The report drafts a
frozen configuration and criteria for a registration under the location-only rule; nothing of it has been
registered or run.

**C14: the question.** C13 found that the false atoms in windows with events are made by the events. C14
asked whether that self-interference, the residue of deflation, the coda, cross-phase ghosts and copies
of events already found, is what spends the false-event budget at the operating point, and whether
remedies for it raise recall.

**C14: design.** A diagnosis on three development networks per cell and regime (1600 to 1602; 360
windows) used the truth to split the statistic of every false atom exactly into the direct waves and the
coda of each event, the clutter and noise, and the deflations already applied, and to classify the false
atoms. Four remedies and one combination, none of which uses the truth, were then developed on the same
networks: exclusivity of arrivals already claimed by an accepted atom (R1), a noise scale that varies in
time with the residual after each deflation (R2), an adaptive deflation with a short shaping filter per
station and phase (R3), fusion of candidates within 10 km and 1 s of an accepted atom (R4), and R2 with
R3; T3 ran as a reference. The constants were chosen by a rule written after an interim table of two of
the three development networks, and the configuration was frozen with hashes at 11:36:11 (-05) on 9
October, before any window of the evaluation networks 1610 to 1621, 12 per cell and regime, was
generated: 2304 windows, 1152 standard with 855 events and 1152 dense with 4594 events, with event-free
twins of the first 8 windows of each type in each part.

**C14: the diagnosis.** At the operating points of the benchmark the largest class of false atoms under
the origin-time rule is not self-interference. It is the detection of a real event placed close to it, a
median of 3.2 km in standard and 3.0 km in dense windows, with the origin time off by a median of 0.76
and 0.86 s, early in 75 and 73 per cent of the cases, and with no other correct detection of the same
event in 86 to 92 per cent: on the development networks it is 55 to 58 per cent of the false atoms in
standard and 80 to 95 per cent in dense windows, against 0 to 7 per cent for the residue of deflation and
none for the coda, and 60 to 76 per cent of its statistic is live energy of its own event in the same
phase. Under the location-only rule these atoms are hits. The residue of deflation matters only at low
thresholds. With the truth, a perfect deflation of direct waves and coda gains $+0.051$ [$-0.024$,
$+0.111$] in standard and $+0.052$ [$+0.012$, $+0.105$] in dense windows under the origin-time rule,
which bounds what any remedy of the deflation can give. The fractions rest on 44 and 46 false atoms at
the operating point, with a sampling uncertainty of about 0.07, and the class of an atom is a convention,
the largest part of an exact decomposition.

**C14: evaluation.** On the unseen networks the remedies do what they are meant to do to the false atoms
and nothing to recall (Table Y.87). R1, R2, R3 and their combination halve the false atoms that the
events induce: against the event-free twin, at the 0.95 quantile of each method's own twin scores, the
ratio falls from 6.23 [4.72, 8.58] to 2.52 [1.99, 3.29] in standard and from 10.04 [8.22, 12.63] to 4.24
[3.61, 5.14] in dense windows with R2 and R3 together, and beyond 10 km from every event from 3.42 and
5.42 to 1.00 and 1.60, the level of the clutter. Recall at $\le0.25$ with origin times matched does not
move: the differences from the base detector are $+0.015$ to $+0.016$ in standard windows, with every
interval containing zero, and $-0.011$ to $-0.002$ in dense windows. The false atoms that the remedies
remove lie below the threshold of the operating point; above it the well-placed, mis-timed detections
remain, 0.082 per window with the base detector and 0.089 to 0.093 with the remedies in standard windows,
0.161 and 0.169 to 0.175 in dense ones, about a third of the 0.25 false events per window in standard and
two thirds in dense windows on these networks. Fusion of duplicates in origin time does not help, because
those atoms are not duplicates, and changes dense recall by $-0.011$ [$-0.018$, $-0.001$]. The refinement
T3 raises recall with origin times matched in standard windows, $+0.056$ [$+0.015$, $+0.100$], but does
not reduce the mis-timed atoms, $-0.003$ [$-0.043$, $+0.038$] per window at any score. The remedies have
a cost that the operating point hides: at any score, among the events of the dense windows that have
another event within 10 km and 2 s, hard exclusivity, the shaping filter and the combination lose about a
quarter of the overlapping events that the base detector finds, $-0.103$ [$-0.121$, $-0.086$], $-0.098$
[$-0.116$, $-0.080$] and $-0.095$ [$-0.115$, $-0.078$] against 0.400, while R2 alone does not, $-0.009$
[$-0.029$, $+0.009$]. The gains of the development winners shrank three to six times on the unseen
networks, from $+0.094$ to $+0.015$ for R1 in standard windows, as selection among three or four
constants on three networks predicts. Under the location-only rule R2 gains $+0.056$ [$+0.010$, $+0.086$]
at $\le0.1$ in standard windows, concentrated in the one-sided cell and with weak clutter, a secondary
reading without multiplicity correction; it costs 30 to 33 per cent more time per window.

**Table Y.87. Exploration C14, evaluation networks 1610 to 1621 of seed 999, frozen configuration, pooled over cells and regimes: recall of the base detector or paired difference against it at $\le0.25$ with origin times matched; ratio of false atoms per window with events to the atoms of the event-free twin at the 0.95 quantile of the method's own twin scores; well-placed, mis-timed atoms per window above the threshold of the operating point; and recall of overlapping events in dense windows at any score.** Report `results/explorations/REPORT_C14.txt`, sections 5 and 6; output `out/c14_analysis_eval.txt`. T3 is a reference arm, not a remedy, and has no twins.

| Method | standard, recall or difference | dense, recall or difference | false-to-twin ratio, standard | false-to-twin ratio, dense | mis-timed atoms per window, standard / dense | overlapping events at any score |
|:--|:--|:--|:--|:--|:--|:--|
| base detector | 0.283 | 0.108 | 6.23 [4.72, 8.58] | 10.04 [8.22, 12.63] | 0.082 / 0.161 | 0.400 |
| R1, hard exclusivity | +0.015 [-0.007, +0.044] | -0.005 [-0.017, +0.008] | 3.18 [2.54, 4.17] | 5.01 [4.16, 6.25] | 0.089 / 0.169 | -0.103 [-0.121, -0.086] |
| R2, time-varying noise scale | +0.016 [-0.006, +0.043] | -0.002 [-0.015, +0.016] | 3.41 [2.67, 4.65] | 5.49 [4.62, 6.57] | 0.093 / 0.174 | -0.009 [-0.029, +0.009] |
| R3, adaptive deflation | +0.016 [-0.006, +0.042] | -0.007 [-0.019, +0.006] | 3.17 [2.50, 4.27] | 4.84 [4.06, 6.03] | 0.090 / 0.174 | -0.098 [-0.116, -0.080] |
| R4, fusion within 10 km and 1 s | +0.002 [-0.007, +0.017] | -0.011 [-0.018, -0.001] | 6.20 [4.58, 8.64] | 10.21 [8.47, 12.86] | 0.084 / 0.157 | -0.046 [-0.060, -0.033] |
| R2 and R3 | +0.016 [-0.007, +0.048] | -0.004 [-0.017, +0.011] | 2.52 [1.99, 3.29] | 4.24 [3.61, 5.14] | 0.091 / 0.175 | -0.095 [-0.115, -0.078] |
| T3, reference | +0.056 [+0.015, +0.100] | +0.011 [-0.007, +0.033] | not computed | not computed | 0.082 / 0.165 | +0.037 [+0.016, +0.057] |

**C14: verdict, deviations and limits.** Negative for the aim. Self-interference exists and the remedies
halve it, but the false-event budget of the operating point under the origin-time rule is spent mainly by
the origin-time bias of well-placed detections of real events, which no remedy of C14 touches, and the
two most aggressive remedies cost a quarter of the overlapping events found at low scores. No remedy is
proposed for registration as a contribution to recall. The run was interrupted on 5 October at 19:02
(-05), with 8 of 18 development parts written, and relaunched on 9 October at 11:10 with the hashes of
the library and the runner checked identical and the written parts skipped. When the second development
round was launched two processes ran at once and a test of one window per type wrote parts; both were
stopped, the directory was deleted and the round was repeated whole with one process. After the freeze
and before any evaluation number was read, the analysis script was changed in three declared ways, one of
them the correction of a code error in the denominator of the far-field twin ratio, and its final hash
was recorded; a code test on the first two evaluation networks wrote an output that was not read. The
runner was extended between the first fixing of the design on 5 October and the freeze, its hash at the
freeze differing from the one recorded then; the report does not itemise that change, and the hashes of
the freeze are those that fix the evaluation. The study used one seed, 12 evaluation networks per cell
and regime and three development networks, so the constants were overfitted, and the diagnosis is that of
the development networks and of the base detector.

### Y.29a What C12 to C14 change

Nothing registered changes. G1-R and G1-C keep the support of Y.24, Y.26 and Y.28, no hypothesis is
added, and the three studies are pilot estimates on one seed, chosen and analysed with the latitude of
exploration. What they change is where the next work points.

C12 closes one direction with a mechanism and leaves an analytic piece. Restricting the moment-tensor
subspace to double couples lowers the threshold by about 10 per cent in energy, which the tube formula
predicts, but the hits lose as much energy to the restriction as the threshold drops, so the amplitude
part of the headroom of Y.22 is not reached through the source subspace. The null law stays as the
calibration of a test of mechanism, not as a detector.

C13 gives a first label-free operating point, under the location-only rule and in this benchmark only. It
does not supply the deployable operating point that Y.18a to Y.26a list as missing, because it carries no
guarantee, has not been run on seeds other than the calibration seed, and fails under the origin-time
rule, which is the primary hit rule from Task 2L on.

C14 supplies the diagnosis that explains why C13 fails under that rule and that gives the calibration
track its next target. At the operating point, the largest class of false events under the origin-time
rule is not an artefact of the greedy deflation but the origin time of detections that are placed
correctly, mostly early and usually the only detection of their event. Task 2L bears on that quantity,
since in its verification the terms S2 narrowed the origin errors of the hits, but C14 used the
registered corrections and measured the class only for the base detector and T3. Whether the terms of
Task 2L, a better velocity or timing model, or an explicit treatment of the early bias of the origin time
shrink that class, and whether a decoy calibration then holds under the origin-time rule, are open
questions for exploration before any registration. Remedies for self-interference stay closed for recall
at the current operating points, with the warning that hard exclusivity and adaptive deflation erase
overlapping events that lower thresholds would expose.

### Y.29b C15: where the origin-time error of well-placed detections comes from

Code: `explore_contrib/c15_lib.py`, `c15_run.py`, `c15_analyze.py`, `c15_diag.py`; report
`results/explorations/REPORT_C15.txt`. Exploratory, seed 999 only, three development networks per cell and
regime never used before (1800 to 1802), new stream tags 1851 and 1852, 20 standard and 20 dense windows
per network, cell and regime (720 windows, 286 events in standard and 1436 in dense windows). There was no
separate evaluation set, so every difference below is a development figure and optimistic. The detector
is the MT-FW of C14 with the registered corrections, reproduced bit for bit.

**Question.** C14 left open whether the mis-timed detections come from an artefact of the hard reading
rule, from the reading of the arrivals, from the timing model, or from the location. Under the hard rule
the statistic depends on the data only through the samples read, so an origin shift that leaves every
read sample in place leaves the statistic unchanged and the maximiser in origin time could be an interval
whose earliest point is reported. That is not what happens: the maximum is unique for 97.8 per cent of
the accepted atoms (mean plateau size 1.023), and reporting the midpoint of the plateau changes no recall.

**Decomposition.** For the detections within 5 km of a real event, the origin time implied by the
detector's own read samples, at the true location, with the true travel times and statics, is within
0.5 s for 0.935 of them in standard and 0.773 in dense windows; with the engine's timing model (straight
rays and the registered corrections) at the true location the fractions are 0.930 and 0.761; at the node
of the atom they fall to 0.678 and 0.578, against 0.657 and 0.583 for the reported origins. The read
samples and the timing model suffice for the origin time when the location is right; the error is a
location error transferred to the origin. The fraction within 0.5 s falls from 0.915 for atoms within
1.5 km of the event to 0.684 between 1.5 and 3 km and 0.467 between 3 and 5 km in standard windows (0.778,
0.627 and 0.430 in dense ones), and the origin error correlates with the depth error, -0.40 in standard
and -0.21 in dense windows. Over all scores the median origin error is -0.009 s in standard and -0.074 s
in dense windows, so the early bias that C14 measured at the operating point is not a property of the
detections as a whole.

**Second-stage estimation and cyclic re-estimation.** Estimators that replace only the reported location
and origin of an accepted atom, leaving scores, exclusion and deflation untouched, were compared with the
base report (Table Y.88): the weighted median of the read delays at the node, with and without parabolic
sub-sample interpolation; a robust Gauss-Newton step in the four hypocentral parameters with a prior
centred on the node; and an origin-free equal-differential-time likelihood [151], [152]
on a 0.5 km lattice within 6 km and depths of 1 to 17 km, alone or followed by two rounds of re-reading
the energy peaks at the new location. None excludes zero. With 8 to 16 read samples per atom, of which
roughly half lie more than 0.25 s from the true arrival (56 and 47 per cent within 0.25 s), the location
cannot be refined from the reads of one atom, while the origin time alone can be, because it is a single
parameter and the median absorbs the outliers. A block-coordinate re-estimation after the greedy pass, in
which every atom in turn is returned to the residual and re-selected with the others deflated (the
feedback step of Newtonized orthogonal matching pursuit [153]), was null
or negative on two of the three networks before it was stopped: -0.032 to -0.064 in standard and +0.001 to
+0.022 in dense windows at $\le0.25$ with origin times matched, -0.037 to -0.085 and -0.017 to -0.033
under the location-only rule, at three to five times the cost, because once their competitors are
removed the false atoms of low score gain score.

**Table Y.88. Exploration C15, development networks 1800 to 1802 of seed 999: paired difference against the base report in recall at $\le0.25$ realized false events per window with origin times matched, point [95% hierarchical interval, 100 replicates].** Report `results/explorations/REPORT_C15.txt`, section 3; output `out/c15_analysis_dev.txt`.

| Estimator | standard | dense |
|:--|:--|:--|
| weighted median of read delays | +0.021 [-0.012, +0.069] | +0.003 [-0.008, +0.021] |
| same, sub-sample | +0.028 [-0.012, +0.079] | +0.003 [-0.008, +0.021] |
| robust Gauss-Newton step | +0.042 [-0.005, +0.075] | +0.005 [-0.017, +0.037] |
| equal differential time | +0.024 [-0.019, +0.082] | -0.003 [-0.033, +0.026] |
| same, with two re-reading rounds | +0.038 [-0.020, +0.092] | -0.002 [-0.034, +0.033] |

**C15: verdict.** Null for recall, with a diagnostic result. In this generator the origin-time error of
well-placed detections is a location error, not an error of reading, of the timing model or of the
argmax, and the reads of a single detection on four to eight stations do not carry the information to
correct it. With C9 (the timing channel is the only one that locates), C10 (the remaining headroom is
calibration and amplitude) and C12 to C14, it indicates that the detector is close to what this
generator's information allows and that algorithmic refinements inside the same benchmark yield a few
points at most. The open question of Y.29a about an explicit treatment of the early bias is answered in
the negative: there is no bias of the detections as a whole to remove, and what changes the evidence from
here is data outside the generator.

### Y.29c C16: real noise and a real station geometry with injected events

Code: `explore_contrib/c16_data.py`, `c16_lib.py`, `c16_run.py`, `c16_validate.py`, `c16_analyze.py`; freeze
`C16_FREEZE.txt` (checked after the run, all hashes unchanged); report `results/explorations/REPORT_C16.txt`. Exploratory,
seed 999 only. The synthetic noise and clutter of the generator are replaced by 40 s slots of real continuous noise
from one day (23 July 2010) of 14 three-component stations of the Alto Tiberina network (INGV network IV, 25 Hz, from
the public repository of PyMPA), at least 60 s away from the 40 catalogued origins of that day, band-passed 0.5 to
6.5 Hz, resampled to the project's 16 Hz, and scaled per channel so that its daily median in-band level equals the
generator's noise level; the stations are placed at their true relative positions and the cells of the benchmark are
formed as subsets of them. The events, their mechanisms, pulses, travel times and statics are those of the
generator. A second arm (GEN) runs the standard generator noise with weak clutter on the same networks and the same
signals; it reproduces the standard pipeline bit for bit. Evaluation: 36 new networks per arm (12 per cell, indices
2402 to 2413), 432 standard windows (347 events), 432 dense windows (1720 events) and 432 event-free windows, frozen
before any evaluation network was generated. The per-arrival signal-to-noise ratio of the two arms is matched (paired
ratio REAL to GEN, median 0.96, quartiles 0.88 and 1.08).

**Bounded influence on real noise.** The ordering of G1-R survives and grows (Table Y.89). In dense windows under
the location rule at $\le0.25$, MT-FW with the cap gains $+0.259$ [$+0.131$, $+0.321$] over BOMP-PK and HUBER-N
$+0.294$ [$+0.178$, $+0.355$] with real noise, against $+0.069$ [$-0.009$, $+0.158$] and $+0.113$ [$+0.052$,
$+0.199$] with generator noise; the difference between the arms is $+0.190$ [$+0.037$, $+0.276$] and $+0.181$
[$+0.037$, $+0.264$]. The mechanism is the one the hypothesis names. Real noise carries transients that the generator
lacked (excess kurtosis of the filtered noise with a median of 3.0 and a 90th percentile of 42), BOMP-PK, whose block
energies are unbounded, collapses (dense recall under the location rule from 0.312 to 0.055, standard from 0.813 to
0.248, and a 99th percentile of the first-atom score of the event-free windows of 7536 against 33 with generator
noise), and the capped and Huber statistics lose much less. Under the origin-time rule the two gains are $+0.085$
[$+0.030$, $+0.124$] and $+0.133$ [$+0.082$, $+0.174$], and HUBER-N stays non-inferior to the cap, $+0.035$ [$-0.020$,
$+0.097$].

**Calibration on real noise.** The ordering of G1-C keeps its sign and is not resolved: in standard windows under the
origin-time rule the terms S2 with T3 gain $+0.049$ [$-0.083$, $+0.149$] over the registered corrections with T3 with
real noise and $+0.023$ [$-0.058$, $+0.106$] with generator noise, with a bootstrap SD of 0.060, five times that of
Task 2L. The terms themselves improve as on synthetics (median rms error 0.153 to 0.080 s with real noise), but in this
real geometry, with stations 5 to 31 km from the sources instead of 10 to 45 km, even the oracle has little room over
the registered corrections in standard windows.

**Operating point.** At the thresholds of the event windows, the capped detector emits 0.000 to 0.023 false events per
window on real event-free windows, HUBER-N 0.03 to 0.10 and BOMP-PK 0.12 to 0.20; with the thresholds of the GEN arm
applied to real noise the capped detector stays at 0.000 to 0.060 under the primary rule, and BOMP-PK reaches 0.29 to
0.86 per window. The operating point of the capped detector transfers to real noise with a small excess; that of
BOMP-PK does not.

**Table Y.89. Exploration C16, evaluation networks 2402 to 2413 of seed 999, dense windows, location rule, $\le0.25$ realized false events per window: paired differences, point [95% hierarchical interval].** Report `results/explorations/REPORT_C16.txt`, section 0; output `out/c16_analysis_eval.txt`.

| Contrast | real noise (REAL) | generator noise (GEN) |
|:--|:--|:--|
| MT-FW cap 6 minus BOMP-PK | +0.259 [+0.131, +0.321] | +0.069 [-0.009, +0.158] |
| HUBER-N minus BOMP-PK | +0.294 [+0.178, +0.355] | +0.113 [+0.052, +0.199] |
| HUBER-N minus MT-FW cap 6 | +0.035 [-0.020, +0.097] | +0.045 [-0.025, +0.111] |

**C16: verdict and limits.** Exploratory support for G1-R from a second noise model with real non-Gaussian
transients and a real geometry, by the mechanism that defines the hypothesis; no resolution of G1-C at this sample
size. One site and one day, cells sharing noise hours, 16 Hz, no instrument response, possible uncatalogued
microearthquakes in the slots, synthetic sources, no clutter in the real arm, a one-sided cell whose azimuthal gap
(167 to 202 degrees) is smaller than the synthetic one, and 12 networks per cell. The quality rules written later for
the confirmation found flat telemetry fill in 394 of the 2037 noise slots of C16 (19 per cent, at three stations), so
part of the advantage of the cap over BOMP-PK on this day may come from instrumental defects rather than seismic
transients; C16 is not reanalysed and the confirmation applies those rules. A registered confirmation on other
days of the same network and never-used seeds is the natural next test.

## Y.30 Task R0: the first registered test on real data (Ridgecrest, 2019)

Code: `r0/r0_lib.py`, `r0_run.py`, `r0_picks.py`, `r0_analyze.py`, `queue_eval.sh`, `queue_resume.sh`,
`merge_parts.py`; development log `R0_DEV_LOG.txt`; freeze `R0_FREEZE.txt`; post hoc log
`R0_DEV_LOG_POSTFREEZE.txt`; output `results/taskR0/analysis_eval.txt`; registration, amendment, freeze and result in
`PREREGISTRO_1j_1k.txt` (Task R0).

**Data and design.** One continuous hour, 4 July 2019 from 19:00 to 20:00 UTC, of the early aftershock sequence of
the Ridgecrest Mw 6.4 earthquake, recorded by 14 stations of the Southern California Seismic Network, eleven of them
broadband three-component, with the analyst catalogue of 102 events (M 1.49 to 4.50) and no picks, taken from a public
repository (benchflow-ai/skillsbench, Apache-2.0). The registration was written and its hash sent before any detector,
picker or associator touched the waveforms. The first 15 minutes (CAL) served for station corrections and development;
the remaining 44 min 20 s (EVAL) held 45 target events of magnitude 2.0 or more. Four subnetworks were fixed by rule:
the eight, six and four broadband stations nearest to the sequence (N8, N6, N4) and the eight in the narrowest
azimuthal sector (U8). The pick-free detectors ran on 40 s windows stepped by 3 s, each owning the atoms whose origin
fell in its central 3 s; the medium was homogeneous (VP 6.0 km/s, VP/VS 1.73) with station corrections. The method of
the thesis (M2) was MT-FW with the cap of 6, the terms S2 estimated without picks and the refinement T3; the
comparators were PhaseNet with its original model and GaMMA or PyOcto over the configuration grids of Task 2A, with raw
and corrected picks, the best of the four sets in each subnetwork and bootstrap replicate (BESTPICK). The primary
endpoint was recall of the targets at 30 false detections per hour with the epicentre within 5 km and the origin time
within 0.5 s, by exact FROC, with a bootstrap of 60 s blocks.

**Development, declared before evaluation.** On CAL the registered detector found 3 of 16 targets at any score on N8
against a recall of 0.625 for PyOcto at 30 per hour. After the corrections the timing residuals of the calibration
events kept a scale of 0.35 to 0.41 s, against about 0.14 s on synthetics, beyond the registered reading tolerance of
0.25 s. A registered amendment kept the primary unchanged, applied the budget reduction foreseen in the registration
(origins 2.0 to 11.9 s, at most six atoms per window, the Huber arm dropped for cost), normalised each station by its
own noise level, gave the comparators the corrected picks as in Task 2A, and added a secondary arm M6 whose reading
half-width was fixed by rule from the measured timing level (2 sigma, 0.69 s on N8) without running it on CAL. The run
was interrupted by a container restart with N8 complete and resumed with the frozen code verified unchanged.

**Result.** Negative (Table Y.90). At 30 false detections per hour with origin times matched, M2 recalled 0.022,
0.000, 0.000 and 0.000 of the targets on N8, N6, N4 and U8, and the best comparator 0.600, 0.667, 0.000 and 0.289; the
pooled difference was $-0.383$ [$-0.481$, $-0.310$], so neither the non-inferiority rule R0-a nor the per-subnetwork
rule R0-b was met. The wider tolerance did not help: M6 minus M2 was $+0.000$ [$-0.016$, $+0.019$]. Under the
location-only rule the pooled difference was $-0.450$ [$-0.563$, $-0.385$]. Four stations defeated every method. An
independent recomputation of the point recalls of N8 and N6 with a separate script gave the same values. The corrected
picks of the comparators were mis-specified (error E54, found later by an independent code audit): the corrections
are delays of matched-filter energy peaks, and subtracting them from PhaseNet onset picks also subtracts the lag of
the peak behind the onset, about 0.4 to 0.6 s for S. BESTPICK takes the best of the raw and corrected arms, so the
decision is unaffected; the corrected arms read alone are biased.

**Table Y.90. Task R0, EVAL segment (45 targets of M 2.0 or more), recall at 30 false detections per hour with origin times matched, per subnetwork, and the pooled registered contrast, point [95% block-bootstrap interval].** Output `results/taskR0/analysis_eval.txt`.

| Method | N8 | N6 | N4 | U8 |
|:--|:--|:--|:--|:--|
| M2, MT-FW cap 6 + S2 + T3 (thesis) | 0.022 | 0.000 | 0.000 | 0.000 |
| M6, the same with the 2-sigma reading | 0.000 | 0.000 | 0.000 | 0.022 |
| PyOcto (best configuration) | 0.600 | 0.667 | 0.000 | 0.289 |
| GaMMA (best configuration) | 0.489 | 0.511 | 0.000 | 0.222 |
| Pooled M2 minus best comparator | -0.383 [-0.481, -0.310] | | | |

**Diagnosis, post hoc.** For the targets with a pick-free atom within 3 s, the nearest such atom lay 8 to 18 km from
the catalogued epicentre in median (12.8 km for M2 on N8), against 1.5 km for PyOcto, and well-placed atoms were
outscored by 14 to 20 per cent of all atoms, while PyOcto's were outscored by none. The detectors fail twice: in
location, because the homogeneous medium with constant station terms and the straight-ray radiation operator leave
timing errors of several tenths of a second that grow with the separation between sources (on CAL the semivariogram of
the residuals at the same station and phase rose from 0.25 to 0.42 s squared between source separations below 3 km and
above 10 km, and the P residual correlated with source depth at -0.50) and polarizations that match neither the straight
ray nor the vertical; and in ranking, because a sequence that fills every window with energy saturates the capped
statistic, so that the score of a real event does not stand out. Pick-based association with tolerances of 0.5 to 2 s
absorbs exactly the timing error that the coherent reading cannot.

**Verdict and limits.** The synthetic advantage of G1 and of G1-C did not transfer to this real hour. It is one
sequence, one network and one hour, with an analyst catalogue incomplete below about magnitude 2, a scalar instrument
response and a homogeneous medium, and it says nothing about the Colombian network test G1', which remains registered.
What it does establish is the condition under which the pick-free approach can work on real data: a travel-time model
whose errors are smaller than the reading tolerance, which the synthetic generator granted and a real sequence did not.
Among the pick-free arms every recall lay between 0.000 and 0.044, so their contrasts carry no information: the terms
S2 with T3 against the registered detector gave $-0.006$ [$-0.022$, $+0.014$] and against the refined detector with
the registered corrections $+0.006$ [$+0.000$, $+0.019$]. The contrast with BOMP-PK was completed on N8 only (0.044
against 0.022) and is declared not done.

**Closure, declared post hoc.** The raw outputs of R0 were lost with the container that produced them (error E49), so
the contrast with BOMP-PK was completed by re-running M1 and M3 on the evaluated segment of all four subnetworks with
the frozen code, through a driver that imports `r0_run.py` unchanged and writes resumable blocks of windows
(`r0/r0_run_chunks.py`, `queue_close.sh`; outputs in `results/taskR0/run_close/`), and analysing them with
`r0_analyze.py` unchanged. The re-run reproduces the registered recalls of M1 on every subnetwork and of M3 on N8.
The registered secondary contrast M1 - M3 is +0.000 [-0.023, +0.021] with origin times matched at 30 false detections
per hour (N8 +0.022, N6 0.000, N4 0.000, U8 -0.022) and -0.050 [-0.101, -0.013] under the epicentre rule, where
BOMP-PK, with unbounded energies, recalls 0.044 to 0.089 against 0.000 to 0.044 for the cap. On a real sequence the cap
brings no advantage over BOMP-PK and under the looser rule a small disadvantage, as the saturation of Proposition
AA.9.3 and exploration C18 predict; at 30 false detections per hour every pick-free detector stays below 0.09 (at 60 per
hour BOMP-PK reaches 0.178 on N8 under the epicentre rule).


## Y.31 Task 2M: bounded influence with independent real noise (Alto Tiberina, 2010)

Code: `t2m_data.py`, `t2m_stations.py`, `t2m_lib.py`, `t2m_run.py`, `t2m_validate.py`, `t2m_analyze.py`,
`t2m_power.py`, `t2m_prep.sh`, `t2m_queue.sh`, `t2m_freeze.sh`, `t2m_resume.sh`, `t2m_resume_each.sh`; freeze
`T2M_FREEZE.sha256` (45 files); configurations `cfg_2718.json` and `cfg_3141.json`; output `analyze_2m_out.txt`;
post hoc breakdown `posthoc/blocks.py` and `blocks_out.txt`; independent verification `verif_2m/`; registration,
freeze, incident and result in `PREREGISTRO_1j_1k.txt` (Task 2M). The code is in `code/`; outputs, post hoc breakdown and
verification are in `results/task2M/`.

**Question and design.** Exploration C16 (Y.29c) found that with real noise of one day of the Alto Tiberina network the
advantage of the capped and Huber statistics over BOMP-PK grew from +0.07 and +0.11 to +0.26 and +0.29, but it was a
single seed and a single day, and 19 per cent of its noise slots turned out to contain flat telemetry fill. Task 2M asks
the same question on noise of days that C16 never touched, after a quality rule written before any injection or
detector run: a sample is invalid if it belongs to a run of at least 25 identical counts (1 s) or departs from the
median by more than 10^4 robust deviations, and no slot with an invalid sample of any station within 60 s is used. A
station rule then keeps the largest station set that leaves enough valid slots. Day 2010-12-31 kept 10 stations and
1826 slots of 40 s; days 2010-02-10 and 2010-02-11 kept the same 8 stations and 737 and 531 slots; 2010-02-12 had no
usable slot with 8 stations. The noise of these days differs in character, measured on the noise alone before
registration: the December day has transients (median excess kurtosis of the matched-filter output 3.2, half of the
slots with a peak above 10 robust deviations) and the February days are close to Gaussian (below 0.5, 2 to 4 per
cent). Before registering, the author fixed three decisions: the quality and station rules unchanged; the December day
(seed 2718) as the only primary replicate and the February days (seed 3141) as a secondary replicate with a written
prediction, because the expected February effect was of the size of the generator-noise arm of C16, for which the
power of the cap rule was 0.22 and that of requiring both seeds 0.07; and the generator-noise arm (GEN) as a
descriptive control on 8 networks per cell. The scope of a positive is therefore real noise with transients from one
network and one day.

Synthetic events from the project generator, with the station statics, catalogue and station corrections of C16, were
injected into 40 s noise segments at the real station positions; the cells are random subsets of 4 and 6 stations and
the 8 stations nearest a random azimuth. Each network has 16 dense windows (up to 7 overlapping events), 6 standard
windows and 4 event-free windows; 24 networks per cell on the December day and 17 on the February days. The primary
rules, on the dense windows of seed 2718 under the epicentre rule (within 5 km) at no more than 0.25 realised false
events per window, are 2M-a (cap minus BOMP-PK) and 2M-b (HUBER-N minus BOMP-PK): lower 2.5th bootstrap percentile
above 0 and point estimate at least +0.03, both required, with Holm over the two tests as the strict reading. The
container that had prepared the task was lost before registration; the environment was rebuilt with the recorded
library versions and the preparation was repeated and checked: the noise of the C16 day was regenerated bit for bit
(sha256 equal to the one in the C16 freeze), the station rule gave identical outputs on all five days, and six seed-999
networks regenerated with the 2M code reproduced the catalogue draws of the original validation (which checks the
random streams but not the atom content, because the C16 part files were not preserved).

**Execution.** The registration hash was sent to the author before any computation on the new seeds, the code was
frozen with SHA-256 and the run was single. Two networks of seed 3141 (ladder_4 2909 and 2911) exhausted their reserve
of catalogue slots; the registration declared such networks failed, but the runner did not catch the exception and
stopped the lane, so the remaining networks were run one by one with the unchanged frozen code after an incident note
was appended to the log, before any part was analysed. A container restart interrupted the descriptive GEN arms; they
were resumed with the frozen code verified. Seed 2718 ran all 72 networks and seed 3141 ran 49 of 51.

**Result.** Both primary rules are met on the December day, and G1-R is confirmed with real noise with transients
(Table Y.91). The cap gains +0.228 [+0.148, +0.298] over BOMP-PK and HUBER-N +0.241 [+0.179, +0.293]; Holm rejects both
(adjusted p 0.0010). BOMP-PK recalls 0.097 of the dense events at the operating point, against 0.325 for the cap and
0.339 for HUBER-N. With the origin time required the gains are +0.085 and +0.104, and in standard windows +0.423 and
+0.346, all with lower bounds above 0. HUBER-N is not shown non-inferior to the cap on this day (+0.014 [-0.038,
+0.056], the margin was -0.03). On the near-Gaussian February days both rules are also met, with smaller gains: +0.057
[+0.004, +0.121] for the cap, above the registered expectation, and +0.120 [+0.070, +0.168] for HUBER-N, which there is
also better than the cap (+0.062 [+0.012, +0.107]); in standard windows neither differs from BOMP-PK.

**Table Y.91.** Task 2M, REAL arm, dense windows, epicentre rule, recall at no more than 0.25 false events per window,
cells pooled; 95 per cent hierarchical bootstrap intervals (B = 2000).

| Seed (days) | windows / events | cap 6 | HUBER-N | BOMP-PK | cap - BOMP-PK | HUBER-N - BOMP-PK |
|:--|:--|:--|:--|:--|:--|:--|
| 2718 (2010-12-31), primary | 1152 / 4598 | 0.325 | 0.339 | 0.097 | +0.228 [+0.148, +0.298] | +0.241 [+0.179, +0.293] |
| 3141 (2010-02-10, 11), secondary | 784 / 3138 | 0.355 | 0.418 | 0.298 | +0.057 [+0.004, +0.121] | +0.120 [+0.070, +0.168] |

**Where the advantage comes from.** The generator-noise arm, run on the first 8 networks per cell, gives on those same
networks a cap gain of +0.025 and a real-minus-generator difference of +0.008 [-0.081, +0.111] on the December day,
which at first sight contradicts the pooled +0.228. The networks are tied to blocks of hours, and the first 8 use the
night hours. A breakdown by block, written after the result and therefore exploratory, resolves it: at night (hours 0
to 7, median window kurtosis 0.6) the cap gains only +0.033 [-0.044, +0.128] while HUBER-N gains +0.123 [+0.057,
+0.200]; from 7.5 to 15 h (kurtosis 8.4, 64 per cent of windows with a peak above 10 deviations) BOMP-PK falls to a
recall of 0.003 and the gains are +0.300 [+0.231, +0.339] and +0.237 [+0.128, +0.298]; from 16 to 23 h they are +0.223
and +0.236. On the February days the cap gains +0.04 to +0.09 and HUBER-N +0.10 to +0.15 in every block. The advantage
of the cap over unbounded energies is therefore a property of noise with transients, consistent with the mechanism of
C16 and with the registered prediction; the transients are concentrated in the daytime hours, plausibly cultural noise,
which this task does not identify. HUBER-N carries a smaller advantage that persists in quiet noise. At the operating
point of the event windows, the event-free windows of the December day give 0.000 false events per window for the cap,
0.028 for HUBER-N and 0.042 for BOMP-PK.

**Independent verification.** A separate agent, which did not produce the result, searched for a false positive with
its own labelling, FROC and bootstrap code (`verif_2m/task2m_verification.txt`). It reproduced the recalls (0.325,
0.339, 0.097) and both primary contrasts (+0.2277 [+0.150, +0.300] and +0.2414 [+0.179, +0.292]), and found no
leakage or asymmetry: every slot is eligible, unique and outside the margins around catalogued origins, the three
methods share windows and timing tables, and the BOMP-PK scores are finite and proportional to the amplitude of the
transients (181 of its 288 top-scored atoms are false, 94 per cent of them in windows with transients, mostly at one or
two stations at a time). The contrasts stay positive at hit radii from 2 to 10 km, at false-event caps from 0.05 to 2.0
per window, in every cell and after leaving out any one network (+0.222 to +0.247); a bootstrap that clusters the
hours shared by the cells widens the cap interval to [+0.132, +0.314]. Its caveats are those of scope: in quiet
windows the cap gains +0.050, not significant, while HUBER-N keeps about +0.12; the comparison is a capped statistic
against unbounded energies under one common threshold, which is the mechanism of the hypothesis and not a general win
over every detector of the BOMP type; the quality rule let through one transient of about 12,500 times the noise level
in a standard window; and on the February days the lower bound of the cap rule (+0.004) is marginal.

**What this does and does not show.** It shows that the ordering of G1-R, established on synthetic noise in Tasks 2F
and 2I, holds on real noise of a network that the method was not tuned on, on days that no earlier analysis used, and
that the size of the gain is governed by how impulsive the noise is. Nor does it show that MT-FW beats a pursuit that bounds its stations too: on seed 999,
BOMP-PK with its energies capped at the same level recovers almost all of the gap (exploration C18f, Y.32,
on seed 999). It does
not show that the detector finds real
earthquakes: the signals are synthetic, and the first real-data test (R0, Y.30) was negative. It is one network, and
the primary replicate is one day; the quality rules and the classification of the days were fixed after inspecting the
noise statistics, without injection or detection; and the cells share hours, so the pooled intervals are somewhat too
narrow.

## Y.32 After R0: what fails on real data (exploration C18)

Code: `explore_contrib/c18_tt.py`, `c18a_timing.py`, `c18a_overlap.py`, `c18b_oracle.py`, `c18c_oracle.py`,
`c18d_detect.py`, `c18d_analyze.py`, `c18e_run.py`, `c18e_analyze.py`, `c18f_run.py`, `c18f_analyze.py`, `c18g_oracle.py`, `c18h_kernel.py`, `c18h_detect.py`;
outputs in
`results/explorations/c18/`.
Exploratory. C18a to C18d use only the first 15 minutes of the Ridgecrest hour (CAL, origins before 19:15), which
served for development and station terms in R0 and was never part of an evaluation; the evaluated segment of R0 is
not read. C18e uses seed 999, new networks (3100 to 3107) and new stream tags (5181 to 5186), checked by search to
appear in no earlier file. The Hi-net swarm of `Geolab-USTC/MatchLocate2` was downloaded and its checksums recorded,
but none of its waveforms was read: it is kept unseen for the next registered real test.

R0 named two causes of its failure: a timing error of 0.35 to 0.41 s with a homogeneous medium and constant station
terms, and a capped score that saturates in continuous seismicity. C18 separates them on CAL, where the 16 catalogued
events of magnitude 2 or more and their hypocentres are known.

**C18a: with a +-1 s search, a better medium seemed not to matter (withdrawn, see C18i).** For each catalogued event the matched-filter energy peaks
were read at every station and phase (joint P and S search within 1 s, energy gate 25, as in C17a), and each peak was
predicted, leaving the event out, by a travel-time model plus the median residual of the same station and phase over
the other events. Five models were compared: the registered homogeneous medium (VP 6.0 km/s, VP/VS 1.73), a homogeneous
medium with fitted VP and VP/VS, the layered model of Hadley and Kanamori for southern California with VP/VS 1.73 or
fitted, and a fitted linear gradient. On N8 (130 readings) the fraction of held-out readings predicted within 0.25 s is
0.59, 0.60, 0.58, 0.59 and 0.59; on U8 0.61 to 0.63, on N6 0.53 to 0.59, on N4 0.53 to 0.55. No model gains more than
about 0.06 on any subnetwork, and the fitted parameters stay close to the registered medium. Removing, in addition, an event
term (the median residual of the held-out event, which a detector absorbs in its origin time) raises the fraction to
0.62 on N8, 0.73 on U8, 0.59 on N6 and 0.54 on N4: part of the error is a common shift per event, but a quarter to a
half of the readings stay outside 0.25 s. On N8 the residuals are not explained by other catalogued events arriving within 2 s (0.61 within
0.25 s with such an overlap, 0.59 without). They are worse for S than for P (0.54 and 0.67 within 0.25 s; 32 and 16 per
cent beyond 0.5 s) and very uneven across stations (0.24 at LRL to 0.81 at SLA). With all magnitudes the N8 fraction falls
to 0.39 to 0.43 for every model. Between 37 and 47 per cent of the readings of the larger events are off by more than the read
tolerance whatever the medium (26 to 54 per cent after the event term, 27 to 46 per cent for the registered
medium): what limits the timing on these data is the reading of the peaks in a dense sequence,
not the velocity model, and the readings cannot be predicted from the source position with 16 events (C17a found the
same for kriging).

**C18b and C18c: the cap removes the evidence of real events, even at the true hypocentre.** For every CAL window the
first-iteration score map of MT-FW (registered engine of R0, +SC corrections) was computed, and two numbers compared:
the best score within 3 km and 0.5 s of each catalogued event, and the maxima of windows whose owned origins have no
catalogued event within 2 s (222 null windows). With the registered cap the true-location score of a real event is not
larger than a null maximum: the probability that it exceeds one (AUC) is 0.49, 0.39, 0.43 and 0.43 on N8, U8, N6 and
N4. Without a cap it is 0.71, 0.69, 0.72 and 0.74; with a cap of 12 noise deviations 0.61 to 0.64, of 24 0.69 to 0.74,
and with a logarithmic compression of the arrival energy (energy E replaced by c^2 log(1 + E/c^2), c = 6, LOG6) 0.66
to 0.72. Normalising each station by its own level, in the window or in the previous 40 s, does not help (0.40 to
0.54), and weighting by the fraction of energy explained by the moment-tensor subspace makes it worse (0.28 to 0.38).
In a continuous aftershock sequence almost every window holds arrivals that reach the cap at several stations, so a
capped score cannot rank a magnitude 2.5 event above the coda and the smaller uncatalogued events. The same maps show
the second failure: the window maximum lies within 5 km of the event in only 6 to 31 per cent of the event windows with
the cap and 0 to 44 per cent without it, and the true hypocentre is the window maximum in 0 to 19 per cent.

**C18d: fixing the ranking does not fix detection on these data.** The full MT-FW loop of R0 was run on CAL with the
cap 6 (R0's M1), cap 24, LOG6, no cap, and BOMP-PK (R0's M3). At 30 false detections per hour, with the origin time
required, all variants recall 0 or 1 of the 16 events on N8 and U8; with the epicentre rule the best is 3 of 16 (BOMP-PK
on U8), while cap 24 and LOG6 reach 2 of 16 on N8 and none on U8; among the bootstrap intervals of the differences from the cap 6
only cap 24 on N8 excludes zero (+0.125 [+0.043, +0.300]), and on U8 it is -0.062. The ranking gained in C18c does not turn into
detections, because the detections are not where the events are.

**C18e: the logarithmic compression keeps most of the robustness to transients (one exploratory seed).** If a less bounded statistic is to replace
the cap, it must not lose what Task 2M confirmed. On real noise of 2010-07-23 after the quality rule (median window
excess kurtosis 3.2), 8 networks per cell of seed 999, dense windows, epicentre rule, 0.25 false events per window: cap 6
minus BOMP-PK +0.148 [+0.056, +0.261], LOG6 minus BOMP-PK +0.133 [+0.052, +0.224], HUBER-N minus BOMP-PK +0.175
[+0.095, +0.265], cap 24 minus BOMP-PK only +0.066 [+0.010, +0.141]; LOG6 minus cap 6 is -0.015 [-0.099, +0.043], and
with the origin time required +0.015 [-0.022, +0.060]. A cap of 24 gives up most of the robustness (-0.082 [-0.183,
-0.014] against the cap 6); the logarithmic compression keeps nearly all of it while restoring most of the ranking of
real events lost by the cap.

**C18g: down-weighting the S readings makes it worse.** As a first candidate for the registered gate, the oracle of
C18c was repeated with the S peak vectors weighted by 0.5 or removed, since C18a found the S readings the least
predictable. Both lose: with the cap the AUC falls from 0.39-0.49 to 0.23-0.32 (S weighted 0.5) and 0.17-0.24 (P only),
with LOG6 from 0.66-0.72 to 0.46-0.55 and 0.24-0.33, and the window maximum is never within 5 km of the event once the S
blocks are down-weighted. The S arrivals, though poorly timed, carry most of the evidence; the reading has to be made
more reliable, not discounted.

**C18h: a higher-frequency matched filter helps location a little, and the gate is not passed.** The second candidate
replaced the 2 Hz Ricker kernel of the matched filter, used in every task, by 3 or 4 Hz kernels everywhere a matched
filter is computed (readings, station corrections, scores). The timing ladder does not change (HOM6 within 0.25 s:
0.58 and 0.60 on N8 and U8 at 3 Hz, 0.59 and 0.57 at 4 Hz). At 3 Hz the oracle maps put the window maximum within 5 km
of the event more often with LOG6 (0.38 on N8, 0.50 on U8, against 0.19 and 0.25 at 2 Hz), with AUC 0.73 and 0.72; 4 Hz
is not better. With the full detector on CAL at 3 Hz, LOG6 recalls 3 of the 16 events with origin times matched at 30
false detections per hour on U8 (0 at 2 Hz) and none on N8; with the epicentre rule 3 and 4 of 16. That is the best
pick-free result on these data so far and still below the registered gate of 5 of 16. Several candidates have now been
tried on the same 16 development events, so passing the gate later would carry a selection effect; the registered test
on unseen data is what protects against it.

**C18f: the advantage belongs to the bound, not to the pursuit.** The fifth independent audit asked whether Task 2M
measured an advantage over BOMP-PK as such or over a comparator without any per-station robustness. On the same 384
dense windows as C18e, two robust versions of the comparator were run: BOMP-PK whose block energies are computed
from peak vectors capped at 6 deviations (pursuit, refit and deflation unchanged), and BOMP-PK on matched-filter output
normalised per station by its own level in the window. The cap recovers almost all of the gap: BOMP-PK with the cap
gains +0.122 [+0.029, +0.224] over plain BOMP-PK, and the capped MT-FW and HUBER-N lead it by only +0.026 [-0.009,
+0.081] and +0.052 [-0.026, +0.107]; per-station normalisation alone gains +0.037 [-0.033, +0.123], and MT-FW leads
it by +0.111 [+0.044, +0.186]. In standard windows BOMP-PK with the cap is level with MT-FW (-0.019 [-0.117,
+0.131]). Bounded influence is therefore a property of the per-station statistic that any moment-tensor pursuit can
use, which is how G1-R should be read: Task 2M confirms that bounding each station helps on real noise with
transients, not that MT-FW beats a pursuit that bounds its stations too (Task 2F found the same on synthetics).

**C18i: the search window, not the medium, hid the effect of the medium.** The independent verification of R0 found
that the catalogue station terms and the first reading of S2 search only +-1 s around the homogeneous prediction,
while the real S delays at several stations are 1.7 to 1.9 s; in C18a, 42 per cent of the S readings and 24 per cent
of the P readings on N8 lie within 0.1 s of the edge of that search. With the search widened to +-2 s the readings
of all events admit more wrong peaks (the fraction within 0.25 s of the prediction drops to 0.37 to 0.43 for the
registered medium), and the layered model now predicts them better than the registered medium on every subnetwork:
0.45 against 0.39 on N8, 0.45 against 0.37 on U8, 0.44 to 0.46 against 0.41 on N6 and 0.52 against 0.43 on N4, with
the robust deviation falling from 0.60 to 0.48 s on N8. The conclusion of C18a that the medium does not matter was an
artefact of the truncated search and is withdrawn (error E50); C18b to C18h, which used the truncated station terms,
are kept as measurements of that configuration.

**C18j: with corrected terms and a layered medium the gate is met, by an uncapped statistic.** With the station terms
re-estimated with a +-2 s search, the full detector on CAL (origin times matched, 30 false detections per hour)
recalls, of the 16 events on N8 and U8: with the registered homogeneous medium and the 2 Hz filter, 0 and 0 with the
cap, 3 and 1 with LOG6 and 3 and 4 with BOMP-PK; with the 3 Hz filter, 3 and 2 with LOG6. With the layered medium of
Hadley and Kanamori (travel times of every node from `c18_tt`, terms re-estimated against it; `c18j_layered.py`) and
the 2 Hz filter it recalls 0 and 0 with the cap, 1 and 2 with LOG6 and 4 and 5 with MT-FW without a cap; at 3 Hz the
uncapped detector falls back to 3 and 1. The candidate "layered medium, +-2 s terms, 2 Hz, MT-FW without a cap" reaches
5 of 16 on U8 and meets the registered gate exactly at its threshold. Two qualifications were logged with the result:
about twenty configurations have been tried on the same 16 development events, so passing the gate carries a
selection effect that only the registered test on unseen data can remove; and the candidate does not bound the
influence of any station, which is the opposite of what Task 2M confirmed on real noise. On real seismicity with the
timing corrected, amplitude evidence matters more than robustness, as Proposition AA.9.3 predicts for a saturating
statistic.

**Two qualifications from the independent code audit.** First, the station terms of every C18 reading on CAL were
estimated on the same 16 events that are then counted (error E53). Re-estimating them for each window without the
events that window can hit, the candidate still recalls 5 of 16 on U8 and falls from 4 to 2 on N8, and the LOG6
readings fall by 0 to 2 events; removing every event inside the window as well gives 1 of 16 on U8. The gate is
therefore met under the proper held-out estimate but only exactly at its threshold, and the in-sample readings above
are optimistic by up to two events. Second, the 3 Hz runs changed the kernel of the matched filter but not the pulse
shape used to deflate each accepted atom, and C18h differs from C18i and C18j in that shape (error E51): the 3 Hz
readings move by one event when only the import order changes and are not comparable with each other. The 2 Hz runs,
including the gate candidate, are not affected.

**What follows.** R0 failed for a reason that is now measured: station terms truncated at 1 s, a homogeneous medium and a
capped score each removed part of the evidence, and only together do their corrections bring the pick-free detector to
about a third of the development events with origin times matched, still well below PhaseNet with PyOcto on the same
network (about 0.6 on N8 in EVAL). The registered gate is met, so the next step is a registered test on the Hi-net
swarm, whose waveforms remain unread. Its design is drafted (`preregistration/drafts/taskR1_registration_DRAFT.txt`)
and not registered, because it spends the last unseen data and several of its choices, among them the published 1-D
model of the region, the treatment of the 1 Hz borehole sensors and the split of the 14 catalogued events, are the
author's to make.

## Y.33 Task R1: a hybrid with pick-based association on the Hi-net swarm

Registration of 17 October 2026, 18:49 UTC (log hash 3b238aa4...402e416c, sent before any detector read the
waveforms), amendment 1 of 19:10 UTC (log hash d5d420a6...67125c1, sent before the run), freeze list
`code/r1/FREEZE_R1.sha256` (SHA-256 1669903f...28a2a3); code `code/r1/`, `code/hyb/`; outputs `results/taskR1/`.

**Origin.** After C18j met the gate on the Ridgecrest development segment, the same segment was used (development only)
to give the pick-based comparators the same layered medium and station terms based on pick residuals: PyOcto and GaMMA
with a 1-D model recovered 13 of 16 development events, far above the pick-free detector, and a hybrid that adds
pick-free candidates to the best associator and scores every candidate by the PhaseNet picks it explains recovered 14
and 13. R1 tests that hybrid, not the pick-free detector alone, against the associator it extends.

**Data.** Geolab-USTC/MatchLocate2, commit bac5ccb, `example/`: continuous three-component records at 100 Hz of seven
Hi-net borehole stations (ATKH, INWH, NAZH, ONIH, THTH, TSTH, YNZH), 2 September 2012, 03:20:00 to 03:53:20 UTC, a swarm
of the Fukushima-Yamagata border region. Catalogues: 14 JMA events (M 2.1 to 3.2) and 140 Match and Locate detections
(relative magnitude -0.41 to 3.2), which contain the JMA events. Development (CAL) = origins before 03:34:00 (840 s);
evaluation (EVAL) = origins from 03:34:00 to 03:53:20 minus 40 s (1120 s).

**Registered design.** Truth = the Match and Locate list; targets = truth events of magnitude 0.5 or more (30 in CAL,
54 in EVAL); a detection hitting only smaller events is neutral. Medium: JMA2001 P velocities in the 28-point
approximation of Y. Maeda, VS = VP/1.73, for every method. Station terms from the CAL truth events of magnitude 1.0 or
more. BASE = the configuration of PyOcto-1D or GaMMA-1D, over the grids of Task 2A, with the highest CAL recall at 30
false detections per hour with origin time (ties: fewer CAL detections, then PyOcto, then the first of the grid). HYB =
PICKSCORE-UNION on BASE and P1 (thesis Section 3.6). P1 = MT-FW without cap on the layered engine. BESTPICK = the best
configuration of both 1-D associators. R0-style comparators (homogeneous medium, raw and +SC picks) descriptive. Rules,
fixed by the simulation `code/hyb/sim_r1_rules.py`: R1-a (primary) HYB - BASE, LOC+T, 30/h, exact one-sided McNemar on
the discordant targets, p <= 0.025 and a positive difference; R1-b HYB - BESTPICK, the same test; R1-c descriptive (P1 -
BASE, JMA-only targets, 60 s block-bootstrap intervals, other rates and rules, a threshold fixed on CAL). Declared
prediction: a small gain of HYB over BASE, low power, R1-a probably not met, P1 below BASE.

**Development and amendment 1.** Development read the waveforms up to 880 s, which includes the first 40 s of EVAL
(error E56, below). The registered matched-filter terms
(+-2 s, joint P and S reading) gave P terms of +1.66, -1.70 and +1.60 s at ATKH, TSTH and YNZH against pick terms of
-0.21, +0.17 and +0.23 s: at ATKH and YNZH (12 km, S-P 1.9 s) the P window reaches the S arrival, whose energy is one
hundred to one thousand times that of P in these records, and at THTH and TSTH the P energy of the CAL events is near
or below the gate. The amendment, chosen by an alignment criterion, bounded the reading by half the predicted S-P time
and set the matched-filter terms to the pick terms plus one onset-to-peak lag per phase, the median over station-phases
with at least three gated readings: 0.265 s for P and 0.244 s for S. It also passed JMA2001 to PyOcto and GaMMA as its
(depth, VP) points, because PyOcto's 1-D model does not accept the 0.25 km layers, and fixed the definitions the
registration had left open (the BESTPICK hit set, the JMA-only reading, the CAL threshold, the run sequence). CAL
readings (in-sample): BASE = GaMMA-1D (min_picks 6, sigma 1.0) 16 of 30, HYB 16 (17 with the registered terms), P1 2.
The highest-scoring P1 atoms were ghosts of the S coda of the large events, 6 to 15 km to the south-west at 3 km depth
and 2.7 to 3.9 s after the origin.

**Result.** Table Y.33 (`results/taskR1/run/analysis_eval.txt`). BASE was GaMMA-1D (min_picks 6, sigma 1.0), chosen on
CAL of the same run (17 of 30 CAL targets with 20 detections, a tie with two other GaMMA configurations broken by the
order of the grid); the onset-to-peak lags of the run were 0.272 s for P and 0.244 s for S (PhaseNet on the full record
changes the picks near the end of CAL, so these differ slightly from the development values quoted in the amendment).

**Table Y.33. R1, EVAL (54 targets, 1120 s). Recall (events hit) at <= 10, 30, 60 false detections per hour.**

| Method | LOC+T 10/h | LOC+T 30/h | LOC+T 60/h | LOC 10/h | LOC 30/h | LOC 60/h |
|:--|--:|--:|--:|--:|--:|--:|
| HYB | 0.407 (22) | 0.426 (23) | 0.426 (23) | 0.444 (24) | 0.444 (24) | 0.444 (24) |
| BASE | 0.389 (21) | 0.407 (22) | 0.407 (22) | 0.426 (23) | 0.426 (23) | 0.426 (23) |
| BESTPICK | 0.407 (22) | 0.481 (26) | 0.481 (26) | 0.481 (26) | 0.481 (26) | 0.481 (26) |
| P1 | 0.019 (1) | 0.019 (1) | 0.037 (2) | 0.019 (1) | 0.019 (1) | 0.056 (3) |
| R0-style best | 0.407 (22) | 0.463 (25) | 0.463 (25) | 0.463 (25) | 0.463 (25) | 0.463 (25) |

R1-a: HYB 23, BASE 22; one target only HYB (origin 1591.5 s, M 1.29), none only BASE; exact one-sided McNemar p =
0.5000; not met. R1-b: HYB 23, BESTPICK 26 (PyOcto-1D, n_picks 3, tol 0.5); one only HYB, four only BESTPICK; p =
0.9688; not met. R1-c: P1 1 against BASE 22 (no target only P1). Block bootstrap (60 s blocks, 2000 replicates,
operating point re-chosen): HYB - BASE +0.018 [+0.000, +0.058], HYB - BESTPICK -0.052 [-0.163, +0.000], P1 - BASE -0.386
[-0.490, -0.286]. JMA-only targets at 30/h: HYB and BASE 0.9, BESTPICK 1.0, P1 0.1. Threshold of the CAL operating
point applied to EVAL: HYB recall 0.426 with 13 false detections (41.8 per hour, above the budget), BASE 0.407 with 4
(12.9 per hour), P1 0.019 with 4.

**Reading.** The prediction held. The hybrid added one event to the associator and no gain is demonstrated; the
pick-free detector alone found one target in fifty-four. A homogeneous PyOcto with raw picks did as well as the best
1-D configuration and better than BASE, so on this swarm the published 1-D model did not help pick-based association.
The run was not interrupted. The registration and the amendment were pushed to the public repository about one minute
after the run started; their hashes were sent before it.

**Independent verification** (`results/taskR1/verif_r1/r1_verification.txt`). A separate agent reproduced with its own
code every number of Table Y.33, the McNemar tests, the selection of BASE, the 277 EVAL and 180 CAL detections of the
hybrid with their pick support, the pick terms (to 1e-4 s), the JMA2001 times (within 6 ms of its own ray tracer) and
the PyOcto model file (byte-identical); 297 of 297 frozen hashes and the hash chain passed. No defect changes a
registered number. It found one protocol deviation: the development read the waveforms up to 880 s, the first 40 s of
EVAL, which contain two EVAL targets, and the development lists contain the first of them; the statement of the
amendment that no EVAL sample had been read is false (E56). No constant was set with those samples, and with EVAL
restricted to origins from 880 s the primary rule gives 22 against 21 of 52 targets, discordant 1 and 0. It also found
that only 29 of the 54 targets are detected by any method at any threshold, that two targets of magnitude 2 or more
are detected by none (each about 10 s from a JMA event, with correlation coefficients of 0.57 and 0.53, possibly
artefacts of the template matching), and that the PyOcto model file and the PhaseNet code were not in the freeze list.

---

# Appendix Z. Certified Pruning of the Hypothesis Space, and the Role of p4est

Code: `bnb_octree.py`, `bnb_scaling.py`. Appendix V.4 found that an octree refined by score is cheaper
than a uniform grid and that an observability-based stopping rule never fired. What an exhaustive
search needs for scale is different: a pruning rule that **cannot** remove a hypothesis that the
dense search would have accepted, so that speed is never bought with false negatives, and that emits
events only from exactly evaluated nodes, so that it cannot add false positives.

**The bound.** For a cell of the space-time hypothesis grid with centre $c$, spatial radius $h$ and
origin-time half-width $\delta\tau$, the eikonal Lipschitz property of C-III gives, for every point
$x$ of the cell and every branch $b$,

$$\big|t_b(x)-t_b(c)\big|\le h/V_{b,\min},$$

so the residual of every predicted arrival to its nearest pick is at least
$\max\big(0,\,r_i(c)-h/V_{b,\min}-\delta\tau\big)$. The robust kinematic log-likelihood of Appendix
Y.5 is a sum of terms that decrease monotonically in the residual, so evaluating it at those reduced
residuals is an upper bound on its value anywhere in the cell. A cell is discarded only if that bound
is below the detection threshold; otherwise it is split along its longest index dimension, down to
the nodes of the dense grid, where the score is evaluated exactly. The certificate covers the
kinematic likelihood only. Extending it to the polarization likelihood (9d) is not a direct
consequence: the predicted direction can change abruptly near vertical incidence, the nearest pick
(and with it the measured polarization) can change inside a cell, and the ratio $f/g$ is not
monotone. A valid bound must take the maximum of $f/g$ over every direction and pick reachable from
the cell, and is registered as open work. The pick-free detectors of Appendix Y.10 to Y.12 admit a
simpler bound, because they need no nearest pick: the moment-tensor certificate is the norm of a
projection of the peak three-component vectors, and no projector increases a norm, so the energy of
those vectors in a time window enlarged by $h/V_{b,\min}+\delta\tau$ bounds it over the whole cell;
the same holds for the ray-geometry detector. That bound is valid and makes the Frank-Wolfe oracle
of Section 7.5 certifiable, but, as measured below, it is too loose to prune once the Frank-Wolfe
maxima fall to the noise level; tightening it with the subspaces reachable from the cell remains
open.

**Verification.** On windows of the consistent generator of Appendix Y.5, with zero to two sources
and an absolute threshold, the set of grid nodes above threshold and the emitted catalogue were
compared node by node with the dense search.

| Setting | windows with identical catalogue | dense evaluations | pruned evaluations | reduction |
|:--|:--|:--|:--|:--|
| 4 stations, 3 km and 0.1 s | 12 of 12 | 20,925 | 8,864 | 2.4x |
| 8 stations one-sided, 3 km and 0.1 s | 12 of 12 | 20,925 | 8,668 | 2.4x |
| 4 stations, 3 km and 0.1 s, fixed threshold | 6 of 6 | 20,925 | 9,040 | 2.3x |
| 4 stations, 1.5 km and 0.05 s | 6 of 6 | 153,903 | 53,650 | 2.9x |
| 4 stations, 0.75 km and 0.05 s | 6 of 6 | 594,567 | 197,981 | 3.0x |

The pruning is exact in all 42 windows. The reduction is modest, two to three times in
evaluations and 1.7 to 2.0 times in wall-clock time, and grows
slowly with resolution, because in these windows a large part of the hypothesis space lies above a
threshold calibrated for high recall; it is the price of refusing any pruning that could lose an
event. Tighter bounds (per-station directional Lipschitz constants, polarization bounds) and a
best-first order when only the strongest event per cluster is needed are the registered routes to
larger savings.

**The certified oracle of MT-FW, measured.** Code: `certified_fw.py` (`certified_fw_record`,
`certified_oracle`, `cell_bounds`, `leaf_enclosure`, `check_bound`); output
`certified_fw_report.txt`. The energy bound above was implemented inside the greedy MT-FW iteration,
with best-first search that prunes only cells whose bound is strictly below the incumbent, the tie
rule, exclusion radius, deflation and record of `fw_record`, and the bound rebuilt from the current
residual after each deflation, since it uses only the energies of that residual. On 72 windows of the
generator of Y.12 (pert, weak and strong clutter, random 4 and 6 and one-sided 8 stations, 12 windows
per cell, seed-999 streams) it returns records bitwise identical to `fw_record` in 72 of 72 windows,
in three modes: the bound as specified (Lipschitz), the same bound with exact station-to-cell distance
ranges (box), and a subspace-aware leaf enclosure described below. The bounds are valid: there were no
violations on 48 residual fields times 300 random cells, raw and after deflation, nor on 1,004,400
single hypotheses.

**Cost of the certified oracle over the 72 windows (one core).** Exact node evaluations use the call of
`fw_record` (31 origin times per node). Batched hypotheses: hypotheses whose statistic the enclosure
recomputes in batch, as a fraction of the dense count. Script `certified_fw.py`; output
`certified_fw_report.txt`, sections 3 and 4.

| Cell, clutter | exact node evaluations: dense / Lipschitz / box / enclosure | batched hypotheses | wall time: dense / Lipschitz / enclosure |
|:--|:--|:--|:--|
| 4 stations, weak | 37,506 / 33,783 / 33,783 / 60 | 0.99 | 7.2 / 25.2 / 1.8 s |
| 4 stations, strong | 37,632 / 35,285 / 35,285 / 60 | 1.00 | 7.3 / 26.0 / 1.9 s |
| 6 stations, weak | 37,542 / 36,442 / 36,442 / 60 | 1.00 | 9.4 / 28.6 / 2.2 s |
| 6 stations, strong | 37,569 / 37,043 / 37,043 / 60 | 1.00 | 9.4 / 28.2 / 2.1 s |
| one-sided, 8 st., weak | 37,560 / 36,221 / 36,221 / 60 | 1.00 | 10.9 / 32.2 / 2.5 s |
| one-sided, 8 st., strong | 37,635 / 37,208 / 37,208 / 60 | 1.00 | 11.7 / 32.2 / 2.5 s |
| total | 225,444 / 215,982 / 215,982 / 360 | 1.00 | 55.9 / 172.3 / 13.1 s |

**The energy bound certifies but does not prune.** As specified it avoids only 4 per cent of the exact
node evaluations (215,982 of 225,444) and, with its bookkeeping, is 3 times slower than the dense loop
(172 s against 56 s); the box variant is no better. It prunes only when a strong event sets a high
incumbent: in the first iteration at four stations 70 per cent of the nodes are still evaluated, and in
the noise-dominated later iterations 93 to 100 per cent. The median ratio of the cell bound to the true
cell maximum is 4.99 (Lipschitz) and 4.64 (box). No bound built only from per-block energies can do
much better, because it must dominate the largest $\|Q^{\mathsf T}\mathbf v\|^2$ over all vectors with
those block norms, which equals the full energy whenever the six-dimensional subspace contains such a
vector, while under noise the statistic captures only about a fraction $1/n_s$ of that energy.
Tightening requires the directions of the candidate vectors.

**The enclosure certifies by vectorization, not by pruning.** On cells of at most 16 nodes the
statistic of every hypothesis is recomputed in batch from the same argmax and cap tables, and both
floating-point evaluations of $\|Q^{\mathsf T}\mathbf v\|^2/2$ lie within $60u\|\mathbf v\|^2$
($u=1.1\times10^{-16}$) of the exact value, so the batch value plus $10^{-11}(\|\mathbf v\|^2+1)$ is a
certified upper bound on the value `fw_record` computes (measured slack at most $3.5\times10^{-9}$,
no violations); hypotheses whose sample index or window argmax is not certifiable fall back to the
energy bound. Only nodes whose enclosure reaches the incumbent are evaluated exactly, and in all 360
iterations that was exactly one node, the winner: 360 exact node evaluations against 225,444, and 13.1
s against 55.9 s, 4.3 times faster with identical records. It still recomputes the statistic for
about 100 per cent of the hypotheses (6.97 million batched evaluations). The gain is certified
vectorization of the dense search, not pruning of the hypothesis space. Exactness also fixes the unit
of exact evaluation: BLAS matrix products round a row differently depending on the number of rows in
the product (last-bit mismatches in 5 to 90 per cent of random trials for subsets of 1 to 16 rows of a
31-row product), so evaluating a subset of origin times would change the last bits of the scores and
could flip near-ties; a node, with its 31 hypotheses, is the smallest unit that reproduces
`fw_record` bit for bit. Pruning that pays in noise-dominated iterations needs cell bounds built from
the union of the moment-tensor subspaces reachable from a cell and the directions of the candidate peak
vectors, which remains open.

**Where p4est fits.** p4est provides forests of octrees distributed over many processes, with
space-filling-curve partitioning and 2:1 balance, and scales to very large adaptive meshes [108]; it
is also available through PETSc's forest interface. Its justification here is not a pruning gain:
the kinematic bound saves two to three times in evaluations and 1.7 to 2.0 times in wall-clock time,
the moment-tensor energy bound almost nothing, and the
fastest certified oracle for MT-FW evaluates essentially every hypothesis. For a regional window the
whole hypothesis space fits in the memory of one machine, and the certified oracle runs on it
directly. p4est is justified by two other needs that arise only at national scale: distributing
continuous monitoring of a national network over a cluster, each window or region being one tree of
the forest, and storing and partitioning a hypothesis space and 3-D travel-time tables that no longer
fit in the memory of one machine. The registered plan is therefore: certified branch-and-bound for the
kinematic likelihood and the batched enclosure oracle for MT-FW in the confirmatory prototype, and
p4est only at the stage of national-scale deployment, with the same record-by-record equivalence test
against the dense search as its acceptance criterion.

---

# Appendix AA. Theory of the Group BLASSO Detector

This appendix states and proves what can be proved about the method of Section 7.5 and its pilot
implementation (`task1x_mtglrt.py`, Appendices Y.11 to Y.13, Z). Labels follow the taxonomy of
Section 4.3, with additions used only in this appendix. A **Proposition**, **Lemma** or **Corollary**
is proved here, completely, under the hypotheses it lists; nothing outside those hypotheses is
claimed. An **Assumption** is a hypothesis stated once and invoked by its label. A **Conjecture** is a statement
believed true, whose proof is not available; its proved special case, if any, is cited. An
**Approximation** is a formula used in the method whose error is not bounded here. Where a proof
covers an idealized version of the implemented algorithm rather than the code itself, the gap is
stated next to the result.

## AA.1 Setting, notation and the declared reduced model

**Hypotheses and operator.** A hypothesis is $\theta=(\mathbf x,\tau)\in\Theta=X\times[\tau_-,\tau_+]$,
with $X\subset\mathbb R^3$ a closed box, so $\Theta\subset\mathbb R^4$ is compact. The record of $n_s$
three-component stations and $n_t$ samples is a vector of $\mathbb R^N$, $N=3n_sn_t$. Throughout, data
and operators are whitened: $\mathbf d$ stands for $\Sigma^{-1/2}\mathbf d$ and $G(\theta)$ for
$\Sigma^{-1/2}G(\theta)\in\mathbb R^{N\times6}$, so that all norms below are Euclidean. The observation
model of (12) with the final-stage least-squares misfit is

$$\mathbf d=\sum_{k=1}^{K}G(\theta_k)\mathbf m_k+\mathbf n,\qquad \mathbf m_k\in\mathbb R^6. \tag{AA.1}$$

$P_\theta$ denotes the orthogonal projector of $\mathbb R^N$ onto $\mathrm{range}\,G(\theta)$ and
$r(\theta)=\mathrm{rank}\,G(\theta)\le6$.

**The problem.** For a discrete measure $\mu=\sum_k\delta_{\theta_k}\mathbf m_k$ the forward operator is
$\Phi\mu=\sum_kG(\theta_k)\mathbf m_k$ and the weighted group total variation of (16a) is
$|\mu|_{\mathrm{TV},G}=\sum_k\|G(\theta_k)\mathbf m_k\|$. The substitution
$\mathbf g_k=G(\theta_k)\mathbf m_k$ normalizes every atom by the energy its moment tensor produces in
the whitened network, and removes the directions of $\mathbf m$ that $G(\theta)$ does not see. It
defines the problem that is analysed here, over $\mathbb R^N$-valued Radon measures $\nu$ on $\Theta$
with total variation $|\nu|(\Theta)$ and polar decomposition $\nu=\mathbf u\,|\nu|$,
$\|\mathbf u\|=1$ $|\nu|$-a.e.:

$$(P_\lambda):\qquad \min_{\nu\in\mathcal K}\ F(\nu)=\tfrac12\|\mathbf d-\nu(\Theta)\|^2+\lambda|\nu|(\Theta),\qquad
\mathcal K=\Big\{\nu:\ \textstyle\int\langle\varphi,d\nu\rangle=\int\langle P_\theta\varphi(\theta),d\nu(\theta)\rangle\ \ \forall\varphi\in C(\Theta;\mathbb R^N)\Big\}. \tag{AA.2}$$

$\mathcal K$ is the set of measures whose polar $\mathbf u(\theta)$ lies in $\mathrm{range}\,G(\theta)$
$|\nu|$-a.e.; the forward operator is simply $\Phi\nu=\nu(\Theta)$, and all the physics is in the
constraint. For a discrete measure, $\nu=\sum_k\delta_{\theta_k}G(\theta_k)\mathbf m_k$ gives
$\nu(\Theta)=\Phi\mu$ and $|\nu|(\Theta)=|\mu|_{\mathrm{TV},G}$, so (AA.2) and (16a) coincide. If
$r(\theta)=6$ on $\Theta$, $G(\theta)^+=(G^{\mathsf T}G)^{-1}G^{\mathsf T}$ is continuous and
$d\mu=G(\theta)^+d\nu$ is a bijection between $\mathcal K$ and $\mathcal M(\Theta;\mathbb R^6)$, so the two
problems are equivalent for all measures (with $|\mu|_{\mathrm{TV},G}:=|G\mu|(\Theta)$ for general $\mu$). The **grid problem** is (AA.2) restricted to measures
supported on a finite set $\Theta_h=\{\theta_1,\dots,\theta_n\}\subset\Theta$. With $Q_j$ an
$N\times r_j$ matrix of orthonormal columns spanning $\mathrm{range}\,G(\theta_j)$ and
$\nu=\sum_j\delta_{\theta_j}Q_j\mathbf a_j$, it is the group LASSO

$$\min_{\mathbf a}\ F(\mathbf a)=\tfrac12\Big\|\mathbf d-\sum_{j=1}^nQ_j\mathbf a_j\Big\|^2+\lambda\sum_{j=1}^n\|\mathbf a_j\|. \tag{AA.3}$$

The dual certificate is $\eta(\theta)=\|P_\theta\mathbf r\|/\lambda$, with $\mathbf r$ the residual.

**Arrival structure and the sampled form.** In the pilot, as in any point-source model at high
frequency, each atom is a sum of P and S arrivals,

$$G(\theta)\mathbf m=\sum_{i=1}^{n_s}\sum_{b\in\{P,S\}}\mathbf e_i\otimes\big(A_{ib}(\mathbf x)\mathbf m\big)\otimes\bar s\big(\cdot-t_{ib}(\theta)\big),\qquad t_{ib}(\mathbf x,\tau)=\tau+T_b(\mathbf x,\mathbf r_i), \tag{AA.4}$$

with $A_{ib}(\mathbf x)\in\mathbb R^{3\times6}$ the whitened radiation, propagation and polarization
block and $\bar s$ a unit-energy pulse. Stacking the blocks gives $A(\mathbf x)\in\mathbb R^{6n_s\times6}$
(the matrix returned by `radiation_operator`) and $\Pi_{\mathbf x}$ denotes the projector of
$\mathbb R^{6n_s}$ onto its range. The synthesis operator $U_\theta:\mathbb R^{6n_s}\to\mathbb R^N$,
$U_\theta\mathbf z=\sum_{ib}\mathbf e_i\otimes\mathbf z_{ib}\otimes\bar s(\cdot-t_{ib}(\theta))$, satisfies
$G(\theta)=U_\theta A(\mathbf x)$. Its adjoint samples the matched filter: with
$\mathbf y_i(t)=\langle\mathbf d_i,\bar s(\cdot-t)\rangle\in\mathbb R^3$ (the array `MFn`),
$U_\theta^{*}\mathbf d=(\mathbf y_i(t_{ib}(\theta)))_{ib}$. The pulse autocorrelation is
$a(t)=\langle\bar s,\bar s(\cdot-t)\rangle$, with $a(0)=1$.

**The implemented statistic.** Let $\Delta$ be the sampling interval, $\mathrm{cl}(k)=\min(\max(k,0),n_t-1)$
the clipping of a sample index, and $S$ the half-width of the search window in samples. For each
station and branch the code uses the index window
$I_{ib}(\theta)=\{\mathrm{cl}(\mathrm{cl}(\lfloor t_{ib}(\theta)/\Delta\rfloor)+j):\ |j|\le S\}$ (the code truncates toward zero, which equals the floor because every predicted
time is positive), the signed peak vector $\mathbf v_{ib}(\theta)=\mathbf y_i(j^\ast)$ with
$j^\ast\in\arg\max_{j\in I_{ib}(\theta)}\|\mathbf y_i(j)\|$, and the capped vector
$\mathbf v^\kappa_{ib}=\mathbf v_{ib}\min(1,\kappa/\|\mathbf v_{ib}\|)$. The statistic is

$$T(\theta)=\big\|\Pi_{\mathbf x}\mathbf v^\kappa(\theta)\big\|^2 . \tag{AA.5}$$

The code reports $T/2$ and, for each spatial node, its maximum over the $\tau$ nodes; neither changes
anything below. In the pilot $\Delta=1/16$ s, $S=4$ (so $w=S\Delta=0.25$ s) and $\kappa=6$.

**Assumptions.** Each result lists which of these it uses.

**(A1) Whitened Gaussian noise.** $\mathbf n\sim\mathcal N(0,I_N)$. The code whitens by a single
robust scale (the median absolute deviation of the matched-filter output); the theory treats this
scale as exact.

**(A2) Compactness and constant rank.** $\theta\mapsto G(\theta)$ is continuous on $\Theta$ and
$r(\theta)=q$ is constant ($q=6$ in the pilot geometries). Used only for statements over the
continuum $\Theta$; the grid statements need neither.

**(A3) Pulse model.** (AA.4) holds with a unit-energy pulse $\bar s\in C^2$; the matched-filter kernel
is $\bar s$; the deflation subtracts the shape $a$; arrivals are farther than the support of the
subtracted shape from the record ends. In the code the kernel is a Ricker pulse truncated at
$\pm1$ s and the subtracted shape is built from one truncated at $\pm3$ s; the truncated energy is
$10^{-35}$ of the total, and the shape peaks at zero lag. These differences are ignored.

**(A4) P-S separation with coherence $c$.** At every station, the P and S times used satisfy
$|a(t_{iS}-t_{iP})|\le c<1$. The exact version ($c=0$), used in Proposition AA.6(b)-(d), also requires
$\langle\bar s(\cdot-t_{iP}),\dot{\bar s}(\cdot-t_{iS})\rangle=\langle\bar s(\cdot-t_{iS}),\dot{\bar s}(\cdot-t_{iP})\rangle=\langle\dot{\bar s}(\cdot-t_{iP}),\dot{\bar s}(\cdot-t_{iS})\rangle=0$.
In the pilot $a(0.5\ \mathrm s)\approx0.10$, $a(0.75\ \mathrm s)\approx2\times10^{-3}$ and
$a(1\ \mathrm s)\approx10^{-6}$, while $t_S-t_P\approx0.14$ s per km, so $c\approx0$ beyond about
11 km and $c<1$ always (the P and S windows never share a sample).

**(A5) Eikonal Lipschitz property.** $T_b(\mathbf x,\mathbf r)$ is the first-arrival time,
$\inf\int ds/V_b$ over paths from $\mathbf x$ to $\mathbf r$. For a convex cell on which
$V_b\ge V_{b,\min}$, the concatenation of the segment $[\mathbf x,\mathbf x']$ with an optimal path
from $\mathbf x'$ gives $T_b(\mathbf x,\mathbf r)\le|\mathbf x-\mathbf x'|/V_{b,\min}+T_b(\mathbf x',\mathbf r)$,
hence $|t_{ib}(\theta)-t_{ib}(\theta')|\le|\mathbf x-\mathbf x'|/V_{b,\min}+|\tau-\tau'|$. In the
pilot's homogeneous medium this is the triangle inequality. It fails if the predicted arrival is a
later branch rather than the first arrival.

**Deviations of the implementation from (AA.3).** (D1) The peak search chooses one time shift per
station and branch by per-arrival energy, instead of maximizing the likelihood jointly over shifts.
(D2) The cap $\kappa$. (D3) Whitening by one scale estimated from the same window. (D4) P-S coherence
$c>0$ at close stations. (D5) Deflation (`fw_record`) refits by unpenalized least squares the new atom
only, recomputes peaks on the residual, and never refits earlier atoms or moves positions.
(D6) From the second iteration on, nodes within 5 km epicentral distance of an accepted atom are
excluded (`min_sep`), at most five atoms are taken, and the scores are made monotone by a running
minimum. The 5 km no-reselection radius was fixed before evaluation; it plays the role of the 5 km
catalogue suppression of the single-pass methods and is much smaller than the 20 km suppression that the
single-pass WAVE-MT needs to reach a comparable recall (Appendix Y.11). (D7) $A(\mathbf x)$ is the homogeneous-medium
moment-tensor radiation operator, while the "pert" generator adds statics, polarization deviations
and site factors that no operator in the method models.

**Lemma AA.1 (continuity of the projector).** Under (A2), $\theta\mapsto P_\theta$ is continuous, hence
uniformly continuous on $\Theta$ in the operator norm.

*Proof.* Let $\theta_n\to\theta$. The set of orthogonal projectors of rank $q$ in $\mathbb R^{N\times N}$
is closed ($P^2=P$, $P^{\mathsf T}=P$, $\mathrm{tr}P=q$) and bounded, hence compact. Let $\bar P$ be
the limit of any convergent subsequence of $P_{\theta_n}$. For every $\mathbf m$,
$P_{\theta_n}G(\theta_n)\mathbf m=G(\theta_n)\mathbf m$, so in the limit $\bar PG(\theta)\mathbf m=G(\theta)\mathbf m$
and $\mathrm{range}\,G(\theta)\subset\mathrm{range}\,\bar P$. Both have dimension $q$, so they are
equal and $\bar P=P_\theta$. Every subsequence has a further subsequence converging to $P_\theta$, so
$P_{\theta_n}\to P_\theta$. Uniform continuity follows from compactness. $\square$

**Lemma AA.2 (sampled form of the projector).** Under (A3), $U_\theta^*U_\theta=I+C_\theta$, where
$C_\theta$ couples only the P and S blocks of the same station, with blocks $a(t_{iS}-t_{iP})I_3$, so
$\|C_\theta\|\le c$ under (A4). If $c=0$, then $P_\theta=U_\theta\Pi_{\mathbf x}U_\theta^*$ and
$\|P_\theta\mathbf d\|=\|\Pi_{\mathbf x}U_\theta^*\mathbf d\|$.

*Proof.* Blocks of different stations occupy different coordinates and are orthogonal; the P and S
blocks of station $i$ have inner product $\mathbf z_{iP}^{\mathsf T}\mathbf z'_{iS}\,a(t_{iS}-t_{iP})$;
a block with itself gives $\mathbf z_{ib}^{\mathsf T}\mathbf z'_{ib}$ because $\|\bar s\|=1$. The
off-diagonal part is a direct sum over stations of $\left(\begin{smallmatrix}0&a_iI_3\\a_iI_3&0\end{smallmatrix}\right)$,
with norm $\max_i|a_i|$. If $c=0$, $U=U_\theta$ is an isometry and $M=U\Pi_{\mathbf x}U^*$ is symmetric,
satisfies $M^2=U\Pi_{\mathbf x}(U^*U)\Pi_{\mathbf x}U^*=M$, and has range
$U\,\mathrm{range}\,A(\mathbf x)=\mathrm{range}\,G(\theta)$. Hence $M=P_\theta$, and
$\|U\Pi_{\mathbf x}U^*\mathbf d\|=\|\Pi_{\mathbf x}U^*\mathbf d\|$ because $U$ is an isometry. $\square$

Consequently, with $S=0$ and $\kappa=\infty$ the implemented statistic (AA.5) is exactly
$\|P_\theta\mathbf d\|^2$ for the full-waveform model when $c=0$ and the predicted times fall on samples
(otherwise it is the same quantity with times rounded down to the sampling grid). With $c>0$ it is the same quantity
computed as if the P and S arrivals of each station were orthogonal.

## AA.2 Existence and first-order optimality

**Proposition AA.1 (existence; optimality; the certificate).**
(a) Under (A2), $(P_\lambda)$ has a minimizer, and the residual $\hat{\mathbf r}=\mathbf d-\hat\nu(\Theta)$
is the same for every minimizer.
(b) Under (A2), $\hat\nu\in\mathcal K$ is a minimizer if and only if
(i) $\|P_\theta\hat{\mathbf r}\|\le\lambda$ for every $\theta\in\Theta$, and
(ii) $\int\langle\hat{\mathbf r},d\hat\nu\rangle=\lambda|\hat\nu|(\Theta)$.
Then $\eta(\theta)=\|P_\theta\hat{\mathbf r}\|/\lambda$ satisfies $\eta\le1$ on $\Theta$, $\eta=1$ on
$\mathrm{supp}|\hat\nu|$, and the polar of $\hat\nu$ is $P_\theta\hat{\mathbf r}/\lambda$, $|\hat\nu|$-a.e.
For a discrete minimizer this reads $P_{\theta_k}\hat{\mathbf r}=\lambda\,G(\theta_k)\hat{\mathbf m}_k/\|G(\theta_k)\hat{\mathbf m}_k\|$.
(c) For the grid problem (AA.3) the same statements hold without (A2): a minimizer exists, the
residual is unique, and $\hat{\mathbf a}$ is a minimizer if and only if
$Q_j^{\mathsf T}\hat{\mathbf r}=\lambda\hat{\mathbf a}_j/\|\hat{\mathbf a}_j\|$ when $\hat{\mathbf a}_j\ne0$
and $\|Q_j^{\mathsf T}\hat{\mathbf r}\|\le\lambda$ when $\hat{\mathbf a}_j=0$; note $\|Q_j^{\mathsf T}\mathbf r\|=\|P_{\theta_j}\mathbf r\|$.

*Proof.* (a) $\mathcal M(\Theta;\mathbb R^N)$ is the dual of the separable space $C(\Theta;\mathbb R^N)$
($\Theta$ is a compact metric space), so bounded sets are weak-* sequentially compact (Banach-Alaoglu).
For fixed $\varphi$, $\theta\mapsto P_\theta\varphi(\theta)$ is continuous by Lemma AA.1, so both sides
of the defining identity of $\mathcal K$ are weak-* continuous in $\nu$ and $\mathcal K$ is weak-* closed.
$\nu\mapsto\nu(\Theta)$ is weak-* continuous (it integrates constant functions), and
$|\nu|(\Theta)=\sup\{\int\langle\varphi,d\nu\rangle:\|\varphi\|_\infty\le1\}$ is weak-* lower
semicontinuous. Every $\nu$ with $F(\nu)\le F(0)=\frac12\|\mathbf d\|^2$ has
$|\nu|(\Theta)\le\|\mathbf d\|^2/(2\lambda)$. A minimizing sequence therefore has a weak-* convergent
subsequence with limit in $\mathcal K$, at which $F$ is at most the infimum. This is the argument of
[113] for scalar measures. If two minimizers had $\nu(\Theta)\ne\nu'(\Theta)$, the strict convexity of
$\frac12\|\mathbf d-\cdot\|^2$ and the convexity of the total variation would give
$F(\frac12(\nu+\nu'))<\min F$, a contradiction.

(b) *Sufficiency.* For any $\nu\in\mathcal K$, expanding the square,
$F(\nu)-F(\hat\nu)\ge\lambda|\nu|(\Theta)-\lambda|\hat\nu|(\Theta)-\langle\hat{\mathbf r},\nu(\Theta)-\hat\nu(\Theta)\rangle$.
By the definition of $\mathcal K$ (applied to the constant function $\hat{\mathbf r}$) and (i),
$\langle\hat{\mathbf r},\nu(\Theta)\rangle=\int\langle P_\theta\hat{\mathbf r},\mathbf u\rangle d|\nu|\le\lambda|\nu|(\Theta)$,
and by (ii) $\langle\hat{\mathbf r},\hat\nu(\Theta)\rangle=\lambda|\hat\nu|(\Theta)$. Hence
$F(\nu)\ge F(\hat\nu)$.
*Necessity of (i).* Fix $\theta$, $\mathbf h\in\mathrm{range}\,P_\theta$ and $t>0$. Then
$\hat\nu+t\delta_\theta\mathbf h\in\mathcal K$ and $|\hat\nu+t\delta_\theta\mathbf h|(\Theta)\le|\hat\nu|(\Theta)+t\|\mathbf h\|$, so
$0\le F(\hat\nu+t\delta_\theta\mathbf h)-F(\hat\nu)\le-t\langle\hat{\mathbf r},\mathbf h\rangle+\frac{t^2}2\|\mathbf h\|^2+\lambda t\|\mathbf h\|$.
Dividing by $t$ and letting $t\to0$ gives $\langle\hat{\mathbf r},\mathbf h\rangle\le\lambda\|\mathbf h\|$;
with $\mathbf h=P_\theta\hat{\mathbf r}$ this is $\|P_\theta\hat{\mathbf r}\|\le\lambda$.
*Necessity of (ii).* For $|t|<1$, $(1+t)\hat\nu\in\mathcal K$ and
$\phi(t)=F((1+t)\hat\nu)=\frac12\|\hat{\mathbf r}-t\hat\nu(\Theta)\|^2+\lambda(1+t)|\hat\nu|(\Theta)$ is
differentiable with a minimum at $t=0$, so $\phi'(0)=-\langle\hat{\mathbf r},\hat\nu(\Theta)\rangle+\lambda|\hat\nu|(\Theta)=0$.
*Consequences.* Writing $\hat\nu=\mathbf u|\hat\nu|$ with $\mathbf u(\theta)\in\mathrm{range}\,P_\theta$,
$\lambda|\hat\nu|(\Theta)=\int\langle P_\theta\hat{\mathbf r},\mathbf u\rangle d|\hat\nu|\le\int\|P_\theta\hat{\mathbf r}\|d|\hat\nu|\le\lambda|\hat\nu|(\Theta)$.
Equality throughout forces, $|\hat\nu|$-a.e., $\|P_\theta\hat{\mathbf r}\|=\lambda$ and equality in
Cauchy-Schwarz, i.e. $\mathbf u=P_\theta\hat{\mathbf r}/\lambda$. The set $\{\eta=1\}$ is closed
($\eta$ is continuous by Lemma AA.1) and has full $|\hat\nu|$ measure, so it contains
$\mathrm{supp}|\hat\nu|$, the smallest closed set of full measure.

(c) (AA.3) is a finite-dimensional convex function with $F(\mathbf a)\ge\lambda\sum_j\|\mathbf a_j\|$, so it
is coercive and continuous and has a minimizer; uniqueness of the residual is as in (a). The proof of
(b) applies verbatim with $\Theta_h$ in place of $\Theta$ (a finite set is compact and every map on it is
continuous). For $\hat{\mathbf a}_j\ne0$ the perturbation $\hat{\mathbf a}_j+t\mathbf h$, $t\to0^\pm$,
uses the differentiability of $\|\cdot\|$ at $\hat{\mathbf a}_j$ and gives
$Q_j^{\mathsf T}\hat{\mathbf r}=\lambda\hat{\mathbf a}_j/\|\hat{\mathbf a}_j\|$; for
$\hat{\mathbf a}_j=0$ it gives $\|Q_j^{\mathsf T}\hat{\mathbf r}\|\le\lambda$ as in (i). Finally
$\|Q_j^{\mathsf T}\mathbf r\|=\|Q_jQ_j^{\mathsf T}\mathbf r\|=\|P_{\theta_j}\mathbf r\|$. $\square$

The converse of "$\eta=1$ on the support" is false: $\eta$ may touch $1$ where $\hat\nu$ has no atom.
Uniqueness of the residual means that $\eta$ does not depend on which minimizer is computed.

## AA.3 The certificate at zero is the generalized likelihood ratio

**Proposition AA.2 (dual certificate at $\mu=0$ and the GLRT).** Fix $\theta$ and let
$r=r(\theta)$.
(a) $\max_{\mathbf m\in\mathbb R^6}\big[2\langle\mathbf d,G(\theta)\mathbf m\rangle-\|G(\theta)\mathbf m\|^2\big]=\|P_\theta\mathbf d\|^2$,
attained at every $\mathbf m$ with $G(\theta)\mathbf m=P_\theta\mathbf d$. Under (A1) this is
$2\log\Lambda_\theta(\mathbf d)$, with $\Lambda_\theta$ the generalized likelihood ratio of
$H_1:\mathbf d=G(\theta)\mathbf m+\mathbf n$ ($\mathbf m$ unknown) against $H_0:\mathbf d=\mathbf n$. It is
$\lambda^2\eta^2(\theta)$ for the certificate of $(P_\lambda)$ at $\nu=0$, where $\mathbf r=\mathbf d$.
(b) Under (A1) and $H_0$, $\|P_\theta\mathbf d\|^2\sim\chi^2_r$; under $H_1$ it is noncentral
$\chi^2_r$ with noncentrality $\|G(\theta)\mathbf m\|^2$, the squared weight of the atom in
$|\mu|_{\mathrm{TV},G}$.
(c) $\nu=0$ is a minimizer of $(P_\lambda)$ if and only if $\max_\theta\|P_\theta\mathbf d\|\le\lambda$,
and it is then the unique minimizer. The BLASSO declares at least one event exactly when the maximum
of the GLRT field exceeds $\lambda^2$.
(d) Under (A3) and (A4) with $c=0$, $\|P_\theta\mathbf d\|^2=\|\Pi_{\mathbf x}U_\theta^*\mathbf d\|^2$, which is
(AA.5) with $S=0$ and $\kappa=\infty$ when the predicted times fall on samples.

*Proof.* (a) Write $G=G(\theta)$, $P=P_\theta$. Since $G\mathbf m\in\mathrm{range}\,P$,
$\langle\mathbf d,G\mathbf m\rangle=\langle P\mathbf d,G\mathbf m\rangle$, and
$2\langle P\mathbf d,G\mathbf m\rangle-\|G\mathbf m\|^2=\|P\mathbf d\|^2-\|P\mathbf d-G\mathbf m\|^2$,
maximal exactly when $G\mathbf m=P\mathbf d$, which is attainable. Under (A1) the log-likelihood ratio
at $\mathbf m$ is $\frac12(\|\mathbf d\|^2-\|\mathbf d-G\mathbf m\|^2)=\langle\mathbf d,G\mathbf m\rangle-\frac12\|G\mathbf m\|^2$.
(b) $P=QQ^{\mathsf T}$ with $Q^{\mathsf T}Q=I_r$, so $\|P\mathbf n\|^2=\|Q^{\mathsf T}\mathbf n\|^2$ and
$Q^{\mathsf T}\mathbf n\sim\mathcal N(0,I_r)$. Under $H_1$, $Q^{\mathsf T}\mathbf d\sim\mathcal N(Q^{\mathsf T}G\mathbf m,I_r)$
and $\|Q^{\mathsf T}G\mathbf m\|=\|G\mathbf m\|$.
(c) By Proposition AA.1(b) with $\hat\nu=0$ (condition (ii) is trivial), $0$ is a minimizer iff (i) holds
with $\hat{\mathbf r}=\mathbf d$. Then every minimizer has residual $\mathbf d$, i.e. $\nu(\Theta)=0$,
and $F(\nu)=\frac12\|\mathbf d\|^2+\lambda|\nu|(\Theta)>F(0)$ unless $\nu=0$. The same holds on a grid
by Proposition AA.1(c). (d) is Lemma AA.2. $\square$

This is the matched-subspace detector of [118] applied jointly to all stations and both branches. The
statistic depends on $G(\theta)$ only through its range, so it is invariant under any invertible
reparametrization of the moment tensor. With $S>0$ and $\kappa=\infty$, the implemented statistic is
the value of $\|\Pi_{\mathbf x}(\mathbf y_i(t_{ib}+s_{ib}))_{ib}\|^2$ at one particular choice of shifts
$s_{ib}\in[-w,w]$, so it is at most the GLRT with per-arrival static shifts as nuisance parameters; it
is not that GLRT (D1), and under $H_0$ it is not $\chi^2$, because each peak vector is the maximum-energy
vector of $2S+1$ correlated samples.

**Lemma AA.3 (bounded influence of one station).** Let $\mathbf d$ and $\mathbf d'$ differ only in the
records of station $i$. For every $\theta$ and every $S$,
$|T(\theta;\mathbf d)^{1/2}-T(\theta;\mathbf d')^{1/2}|\le2\sqrt2\,\kappa$.

*Proof.* The peak vectors of the other stations are unchanged. By the reverse triangle inequality and
$\|\Pi_{\mathbf x}\|\le1$, the difference is at most
$\big(\sum_b\|\mathbf v^\kappa_{ib}-\mathbf v'^\kappa_{ib}\|^2\big)^{1/2}\le(2(2\kappa)^2)^{1/2}$, since
every capped vector has norm at most $\kappa$. $\square$

This is the precise content of the cap: a burst at one station moves $\sqrt T$ by a bounded amount,
whereas without the cap ($\kappa=\infty$) its influence is unbounded. The cap is the reason the
implemented statistic is not a GLRT; a heavy-tailed noise model (Section 7.5) is its registered
replacement.

## AA.4 Validity of the cell bound and exactness of the pruned search

A **cell** $C$ is a product of index intervals of the grid $\Theta_h=X_h\times\mathcal T_h$. Let
$\mathbf c_C$ be the centre of the bounding box of its spatial nodes, $h_C$ the half-diagonal of that
box, and $\tau_C,\delta_C$ the centre and half-width of its range of origin times, so every node
$(\mathbf x,\tau)\in C$ satisfies $|\mathbf x-\mathbf c_C|\le h_C$ and $|\tau-\tau_C|\le\delta_C$. With
$L_b=h_C/V_{b,\min}+\delta_C$, $t^c_{ib}=t_{ib}(\mathbf c_C,\tau_C)$,
$k^-_{ib}=\mathrm{cl}(\lfloor(t^c_{ib}-L_b)/\Delta\rfloor)$ and
$k^+_{ib}=\mathrm{cl}(\lfloor(t^c_{ib}+L_b)/\Delta\rfloor)$, define the enlarged window
$\bar I_{ib}(C)=\{\mathrm{cl}(k'):\ k^-_{ib}-S\le k'\le k^+_{ib}+S\}$ and

$$B(C)=\sum_{i=1}^{n_s}\sum_{b\in\{P,S\}}\min\Big(\kappa^2,\ \max_{j\in\bar I_{ib}(C)}\|\mathbf y_i(j)\|^2\Big). \tag{AA.6}$$

**Proposition AA.3 (validity of the cell bound).** Assume (A5) on the convex hull of the spatial nodes
of $C$. Then for every node $\theta\in C$, every family of projectors (in particular the moment-tensor
subspaces, the ray-geometry projectors of WAVE-PT and the subspaces of another node used by the
controls) and every record array $\mathbf y$, $T(\theta)\le B(C)$, where $T$ and $B$ are computed from
the same $\mathbf y$. In particular the bound holds for the deflated residual at every iteration.

*Proof.* Fix $(i,b)$ and $\theta\in C$. By (A5),
$|t_{ib}(\theta)-t^c_{ib}|\le|\mathbf x-\mathbf c_C|/V_{b,\min}+|\tau-\tau_C|\le L_b$. Since
$t\mapsto\mathrm{cl}(\lfloor t/\Delta\rfloor)$ is non-decreasing,
$k:=\mathrm{cl}(\lfloor t_{ib}(\theta)/\Delta\rfloor)\in[k^-_{ib},k^+_{ib}]$, so for $|j|\le S$,
$k+j\in[k^-_{ib}-S,k^+_{ib}+S]$ and $\mathrm{cl}(k+j)\in\bar I_{ib}(C)$: $I_{ib}(\theta)\subset\bar I_{ib}(C)$.
By definition of the peak, $\|\mathbf v_{ib}(\theta)\|^2=\max_{j\in I_{ib}(\theta)}\|\mathbf y_i(j)\|^2\le\max_{j\in\bar I_{ib}(C)}\|\mathbf y_i(j)\|^2$,
and $\|\mathbf v^\kappa_{ib}\|=\min(\kappa,\|\mathbf v_{ib}\|)$. An orthogonal projector does not
increase the norm, so
$T(\theta)=\|\Pi\mathbf v^\kappa(\theta)\|^2\le\|\mathbf v^\kappa(\theta)\|^2=\sum_{ib}\min(\kappa^2,\|\mathbf v_{ib}(\theta)\|^2)\le B(C)$,
using that $u\mapsto\min(\kappa^2,u)$ is non-decreasing. Nothing in the argument uses the projector,
nor where $\mathbf y$ came from, so it applies to any subspace family and to the matched-filter
transform of any residual. $\square$

**Corollary AA.1 (exact branch-and-bound).** Consider the search that keeps a stack of cells, starting
from $\Theta_h$, and an incumbent value $L$ equal to the largest $T$ evaluated so far at a node (initially
$-\infty$, or the value at any proposed node); it pops a cell $C$, discards it if $B(C)<L$, evaluates
$T$ exactly if $C$ is a single node (updating the incumbent by the order "larger value, then smaller node
index"), and otherwise splits $C$ into index sub-intervals.
(a) It terminates, and the incumbent is $\max_{\theta\in\Theta_h}T(\theta)$ at the same node that the dense
search with the same tie rule returns.
(b) With the rule "discard if $B(C)<u$" and no incumbent, the set of evaluated nodes with $T\ge u$ is
exactly $\{\theta\in\Theta_h:T(\theta)\ge u\}$.
(c) Both hold on any admissible subset of nodes (for instance the grid minus the 5 km discs of D6),
since a bound over a superset is a bound over the subset.

*Proof.* Termination: the cell tree is finite and each cell is popped at most once. (a) The incumbent is
always a value of $T$ at a node, so $L\le T^\ast=\max T$. If $\theta^\ast$ attains $T^\ast$, every cell
containing it has $B\ge T^\ast\ge L$ by Proposition AA.3, is never discarded, and the cell $\{\theta^\ast\}$
is eventually evaluated. So every maximizer is evaluated and the tie rule selects the same one as the
dense search. (b) A node with $T\ge u$ lies only in cells with $B\ge T\ge u$, which are never discarded;
evaluated nodes are evaluated exactly, so no node with $T<u$ is reported. $\square$

The argument is in exact arithmetic. In floating point, $T$ and $B(C)$ carry rounding errors of relative
order $10^{-15}$, so the exactness of the corollary requires that a cell be discarded only when its bound
is below the incumbent by more than a relative tolerance, for instance $10^{-12}$.

The corollary makes the Frank-Wolfe oracle of Section 7.5 exact on the grid and is the moment-tensor
analogue of the kinematic certificate verified in Appendix Z. It says nothing about the number of
evaluations: the bound ignores the geometry of the subspaces, and the savings must be measured. Appendix Z
measured them with its node-by-node equivalence test: on 72 windows the bound avoids about 4 per cent of
the exact evaluations (215,982 of 225,444), so it makes the oracle exact without pruning it.

## AA.5 Frank-Wolfe with an exact oracle, and what the implemented deflation guarantees

Let $\rho=\|\mathbf d\|^2/(2\lambda)$. Every minimizer of (AA.3) satisfies $\sum_j\|\hat{\mathbf a}_j\|\le\rho$
(because $\lambda\sum_j\|\hat{\mathbf a}_j\|\le F(\hat{\mathbf a})\le F(0)$), so (AA.3) is equivalent to
minimizing $f(\mathbf a,t)=\frac12\|\mathbf d-Q\mathbf a\|^2+\lambda t$, with $Q\mathbf a=\sum_jQ_j\mathbf a_j$,
over the compact convex set $D_\rho=\{(\mathbf a,t):\sum_j\|\mathbf a_j\|\le t\le\rho\}$: for fixed
$\mathbf a$ the best $t$ is $\sum_j\|\mathbf a_j\|$, so $\min f=\min F=:F^\ast$ and $F(\mathbf a)\le f(\mathbf a,t)$
on $D_\rho$ (Lemma 4 of [116]). Frank-Wolfe on $D_\rho$ iterates $x_{k+1}=x_k+\gamma_k(s_k-x_k)$ with $s_k$ a
minimizer of $\langle\nabla f(x_k),s\rangle$ over $D_\rho$ and $\gamma_k=2/(k+2)$. The curvature constant
[124] is $C_f=\sup\frac2{\gamma^2}\big(f(y)-f(x)-\langle\nabla f(x),y-x\rangle\big)$ over $x,s\in D_\rho$,
$\gamma\in(0,1]$, $y=x+\gamma(s-x)$.

**Proposition AA.4 (Frank-Wolfe on the grid problem).**
(a) *Oracle.* With $\mathbf r_k=\mathbf d-Q\mathbf a_k$ and $j^\ast\in\arg\max_j\|P_{\theta_j}\mathbf r_k\|$,
a minimizer of the linear subproblem is $s_k=(\rho\,\mathbf e_{j^\ast}\otimes Q_{j^\ast}^{\mathsf T}\mathbf r_k/\|Q_{j^\ast}^{\mathsf T}\mathbf r_k\|,\ \rho)$
if $\max_j\|P_{\theta_j}\mathbf r_k\|>\lambda$, and $s_k=(0,0)$ otherwise. The oracle is the global
maximizer of the certificate. When the atoms of (AA.3) place their arrivals at the sample-rounded times
$\Delta\,\mathrm{cl}(\lfloor t_{ib}(\theta_j)/\Delta\rfloor)$, as the implementation does, and (A3) and (A4) hold
with $c=0$, $\|P_{\theta_j}\mathbf r\|^2$ is (AA.5) with $S=0$ and $\kappa=\infty$ (Lemma AA.2), for which Proposition AA.3 also holds, so Corollary AA.1 computes
the oracle exactly. (With $c>0$ and $A(\mathbf x)$ of full column rank, $U^*U\succeq(1-c)I$ gives $\|P_\theta\mathbf r\|^2\le(1-c)^{-1}\|\Pi_{\mathbf x}U_\theta^*\mathbf r\|^2$,
and the bound (AA.6) remains valid after multiplication by $(1-c)^{-1}$.)
(b) *Rate.* $C_f\le4\rho^2$ and, for $k\ge1$, $F(\mathbf a_k)-F^\ast\le f(x_k)-F^\ast\le2C_f/(k+2)\le2\|\mathbf d\|^4/(\lambda^2(k+2))$.
The same bound holds for line search, for the fully corrective variant (re-minimizing $f$ over the convex
hull of the atoms selected so far, within $D_\rho$), and for any step appended to the iteration that stays
in $D_\rho$ and does not increase $f$.
(c) *Inexact oracle.* If $s_k$ is replaced by any $s'_k\in D_\rho$ with
$\langle\nabla f(x_k),s'_k\rangle\le\langle\nabla f(x_k),s_k\rangle+\varepsilon$, then
$f(x_k)-F^\ast\le2C_f/(k+2)+\varepsilon$ for $k\ge1$. An oracle returning a node whose certificate is within
$\varepsilon'$ of the maximum, $\|P_{\theta_{j'}}\mathbf r_k\|\ge\max_j\|P_{\theta_j}\mathbf r_k\|-\varepsilon'$,
satisfies this with $\varepsilon=\rho\varepsilon'$. If instead the error is allowed to shrink as
$\varepsilon_k=\frac12\delta\gamma_kC_f$, the bound is $2C_f(1+\delta)/(k+2)$, the form stated in [124].
(d) *A posteriori certificate.* For any $\mathbf a$, with $\mathbf r=\mathbf d-Q\mathbf a$ and
$\rho'=\max(\rho,\sum_j\|\mathbf a_j\|)$,
$F(\mathbf a)-F^\ast\le\mathrm{gap}(\mathbf a)=\lambda\sum_j\|\mathbf a_j\|-\langle\mathbf r,Q\mathbf a\rangle+\rho'\big(\max_j\|P_{\theta_j}\mathbf r\|-\lambda\big)_+$.

*Proof.* (a) $\nabla f=(-Q^{\mathsf T}\mathbf r_k,\lambda)$. For fixed $t'\in[0,\rho]$, the minimum of
$-\sum_j\langle Q_j^{\mathsf T}\mathbf r_k,\mathbf s_j\rangle$ over $\sum_j\|\mathbf s_j\|\le t'$ is
$-t'\max_j\|Q_j^{\mathsf T}\mathbf r_k\|$ (Cauchy-Schwarz blockwise, attained by putting all mass on the
best block along $Q_j^{\mathsf T}\mathbf r_k$). The subproblem reduces to minimizing
$t'(\lambda-\max_j\|Q_j^{\mathsf T}\mathbf r_k\|)$ over $t'\in[0,\rho]$.
(b) For the quadratic $f$, $f(y)-f(x)-\langle\nabla f(x),y-x\rangle=\frac{\gamma^2}2\|Q(\mathbf a_s-\mathbf a_x)\|^2$,
and $\|Q\mathbf a\|\le\sum_j\|Q_j\mathbf a_j\|=\sum_j\|\mathbf a_j\|\le\rho$ on $D_\rho$, so $C_f\le(2\rho)^2$.
By definition of $C_f$, $f(x_{k+1})\le f(x_k)+\gamma_k\langle\nabla f(x_k),s_k-x_k\rangle+\frac{\gamma_k^2}2C_f$.
By convexity, with $x^\ast$ a minimizer, $\langle\nabla f(x_k),s_k-x_k\rangle\le\langle\nabla f(x_k),x^\ast-x_k\rangle\le-(f(x_k)-F^\ast)$.
With $h_k=f(x_k)-F^\ast$: $h_{k+1}\le(1-\gamma_k)h_k+\frac{\gamma_k^2}2C_f$. For $k=0$, $\gamma_0=1$ gives
$h_1\le C_f/2\le2C_f/3$. If $h_k\le2C_f/(k+2)$ then
$h_{k+1}\le\frac{k}{k+2}\frac{2C_f}{k+2}+\frac{2C_f}{(k+2)^2}=\frac{2C_f(k+1)}{(k+2)^2}\le\frac{2C_f}{k+3}$,
since $(k+1)(k+3)\le(k+2)^2$. A variant whose next iterate has $f$ no larger than the Frank-Wolfe point
satisfies the same recursion. This is the argument of [124].
(c) The recursion becomes $h_{k+1}\le(1-\gamma_k)h_k+\gamma_k\varepsilon+\frac{\gamma_k^2}2C_f$; the
induction hypothesis $h_k\le2C_f/(k+2)+\varepsilon$ propagates because the $\varepsilon$ terms add to
$(1-\gamma_k)\varepsilon+\gamma_k\varepsilon=\varepsilon$, and $h_1\le\varepsilon+C_f/2$. With
$\varepsilon_k=\frac12\delta\gamma_kC_f$ the recursion is that of (b) with $C_f(1+\delta)$ in place of $C_f$. For the second
statement, let $M=\max_j\|P_{\theta_j}\mathbf r_k\|$ and $M'\ge M-\varepsilon'$ the value at the returned
node. If $M'>\lambda$ the returned atom has linear value $\rho(\lambda-M')\le\rho(\lambda-M)+\rho\varepsilon'$,
and $M>\lambda$ so the exact value is $\rho(\lambda-M)$. If $M'\le\lambda$ the returned value is $0$ and the
exact value is $\rho\min(0,\lambda-M)\ge-\rho\varepsilon'$.
(d) $(\mathbf a,\sum_j\|\mathbf a_j\|)\in D_{\rho'}$ and the minimizers of (AA.3) lie in $D_\rho\subset D_{\rho'}$.
By convexity $F^\ast\ge f(x)+\min_{s\in D_{\rho'}}\langle\nabla f(x),s-x\rangle$, and by (a) with $\rho'$ the
right side equals $F(\mathbf a)-\mathrm{gap}(\mathbf a)$. $\square$

*The continuum.* The same statements hold for $(P_\lambda)$ over $\Theta$ under (A2), with the maximum
over $\Theta$ (attained by Lemma AA.1): the proofs use only convexity, the quadratic fidelity and the
oracle, as in [116], [117]. A grid oracle is then inexact: for a grid of mesh $h$ and
$\omega_P$ the modulus of continuity of $\theta\mapsto P_\theta$, and since
$\|\mathbf r\|\le\|\mathbf d\|+|\nu|(\Theta)\le\|\mathbf d\|+\rho$ on $D_\rho$,
$\varepsilon'\le\omega_P(h)(\|\mathbf d\|+\rho)$. The sliding step of [116], which moves positions at
fixed number of atoms and does not increase $F$, keeps the bound (b). The finite-step convergence of
sliding Frank-Wolfe (Theorem 3 of [116]) is proved for real scalar amplitudes, stated and proved for a
one-dimensional domain (the authors indicate that the changes for higher dimension are straightforward),
and assumes that the problem at the given $\lambda$ has a unique solution whose certificate is
non-degenerate (strictly below $1$ off the support, with non-vanishing second derivative on it) and that
the non-convex sliding step reaches a stationary point. Its extension to $\mathbb R^6$-valued amplitudes
with the group norm is not claimed here.

**What the implemented deflation is.** Under (A3), subtracting the shape $\mathbf p_{ib}\,a(\cdot-t^\ast_{ib})$
from the matched-filter record is, by linearity of the matched filter, subtracting
$U_k\mathbf p$ from the whitened residual, where $U_k$ is the synthesis operator at the selected peak
times $t^\ast_{ib}$. The code therefore computes, at iteration $k$: the node $\theta_k$ maximizing $T$ over
the admissible nodes; the peak vectors $\mathbf v_k=U_k^*\mathbf r_k$ and their capped version
$\mathbf v_k^\kappa$; the least-squares fit $\mathbf p_k=A\,A^+\mathbf v^\kappa_k=\Pi_{\mathbf x_k}\mathbf v^\kappa_k$;
and $\mathbf r_{k+1}=\mathbf r_k-U_k\mathbf p_k$. This is a group matching pursuit with a physics
dictionary, not Frank-Wolfe: earlier amplitudes are never refitted and there is no penalty.

**Proposition AA.4, continued (the implemented deflation).** Assume (A3) and (A4) with coherence $c_k$
at the peak times of iteration $k$, and write $\mathbf e_k=\mathbf v_k-\mathbf v^\kappa_k$.
(e) *First iterate.* If $S=0$, $\kappa=\infty$ and $c_0=0$, the first selected node is the oracle of (a) at
$\nu=0$ on the grid (no exclusion applies at the first step), and $U_0\mathbf p_0=P_{\theta_0}\mathbf d$. The fully corrective penalized step on the
active set $\{\theta_0\}$ would instead be $(1-\lambda/\|P_{\theta_0}\mathbf d\|)_+P_{\theta_0}\mathbf d$; the
implementation keeps the unshrunk (debiased) amplitude. With $S>0$ (still $\kappa=\infty$, $c_0=0$) the first
selected node is the oracle for the dictionary whose arrival times are the peak times chosen by D1.
(f) *Norm identity.* $\|\mathbf r_{k+1}\|^2=\|\mathbf r_k\|^2-2\langle\mathbf v_k,\mathbf p_k\rangle+\|U_k\mathbf p_k\|^2$
and $\|U_k\mathbf p_k\|^2\le(1+c_k)\|\mathbf p_k\|^2$. If the cap is inactive at the selected node
($\mathbf e_k=0$), then $\|\mathbf r_{k+1}\|^2\le\|\mathbf r_k\|^2-(1-c_k)T(\theta_k)$; if moreover $c_k=0$,
$\|\mathbf r_{k+1}\|^2=\|\mathbf r_k\|^2-T(\theta_k)$ and $U_k^*\mathbf r_{k+1}=(I-\Pi_{\mathbf x_k})\mathbf v_k$: at the
selected samples, the atom is explained away exactly.
(g) *Cap active.* In general $\|\mathbf r_{k+1}\|^2\le\|\mathbf r_k\|^2-(1-c_k)\|\mathbf p_k\|^2-2\langle\Pi_{\mathbf x_k}\mathbf e_k,\mathbf p_k\rangle$,
and the last term has no sign, so no decrease is guaranteed, and
$\Pi_{\mathbf x_k}U_k^*\mathbf r_{k+1}=\Pi_{\mathbf x_k}\mathbf e_k$ when $c_k=0$, which is nonzero in general: a
capped event is not fully explained away. Subtracting $\alpha^\ast U_k\mathbf p_k$ with
$\alpha^\ast=\langle\mathbf v_k,\mathbf p_k\rangle/\|U_k\mathbf p_k\|^2$ instead gives
$\|\mathbf r_{k+1}\|^2=\|\mathbf r_k\|^2-\langle\mathbf v_k,\mathbf p_k\rangle^2/\|U_k\mathbf p_k\|^2$, a
decrease in every case.
(h) *Counting.* If the exact case of (f) holds at every iteration, $\sum_kT(\theta_k)\le\|\mathbf d\|^2$, so at
most $\|\mathbf d\|^2/u$ atoms have $T\ge u$.

*Proof.* (e) With $S=0$, $\mathbf v_0=U_{\theta_0}^*\mathbf d$ and, by Lemma AA.2,
$T(\theta)=\|P_\theta\mathbf d\|^2$, whose maximizer is the oracle of (a) with $\mathbf r=\mathbf d$; and
$U\Pi U^*\mathbf d=P_{\theta_0}\mathbf d$. The penalized step minimizes $\frac12\|\mathbf d-\mathbf g\|^2+\lambda\|\mathbf g\|$
over $\mathbf g\in\mathrm{range}\,P_{\theta_0}$; writing $\mathbf d=P\mathbf d+(I-P)\mathbf d$ this is
$\frac12\|P\mathbf d-\mathbf g\|^2+\lambda\|\mathbf g\|$ plus a constant, whose minimizer is the group
soft-threshold $(1-\lambda/\|P\mathbf d\|)_+P\mathbf d$ (the minimizer is parallel to $P\mathbf d$ by
Cauchy-Schwarz, and the scalar problem is solved by the positive part). (f) Expand
$\|\mathbf r_k-U_k\mathbf p_k\|^2$ and use $U_k^*\mathbf r_k=\mathbf v_k$ and Lemma AA.2. If
$\mathbf e_k=0$, $\langle\mathbf v_k,\Pi\mathbf v_k\rangle=\|\Pi\mathbf v_k\|^2=T(\theta_k)$, so the sum is at most
$-2T+(1+c_k)T$. If $c_k=0$, $U_k^*\mathbf r_{k+1}=\mathbf v_k-U_k^*U_k\Pi\mathbf v_k=(I-\Pi)\mathbf v_k$.
(g) $\langle\mathbf v_k,\mathbf p_k\rangle=\langle\mathbf v^\kappa_k+\mathbf e_k,\Pi\mathbf v^\kappa_k\rangle=\|\mathbf p_k\|^2+\langle\Pi\mathbf e_k,\mathbf p_k\rangle$.
Since $\mathbf e_{k,ib}=\beta_{ib}\mathbf v^\kappa_{k,ib}$ with $\beta_{ib}\ge0$, a block scaling that does not
commute with $\Pi$, the inner product can be negative, and the norm can increase. Example (scalar blocks,
$c_k=0$, $\kappa=1$, $\Pi$ the projector onto $(1,1,1)/\sqrt3$): $\mathbf v=(2,-1,-\frac12)$ gives
$\mathbf v^\kappa=(1,-1,-\frac12)$, $\mathbf p=-\frac16(1,1,1)$, $\langle\mathbf v,\mathbf p\rangle=-\frac1{12}$ and
$\|\mathbf r_{k+1}\|^2-\|\mathbf r_k\|^2=\frac16+\frac1{12}=\frac14>0$. The example is abstract; whether a
moment-tensor subspace realizes it is not examined, so the claim is only that no decrease follows from the
hypotheses. With $c_k=0$,
$\Pi U_k^*\mathbf r_{k+1}=\Pi(\mathbf v_k-\Pi\mathbf v^\kappa_k)=\Pi\mathbf e_k$. The scaled update is the
minimization of the quadratic $\alpha\mapsto\|\mathbf r_k-\alpha U_k\mathbf p_k\|^2$. (h) Telescoping, with
$\|\mathbf r_k\|^2\ge0$. $\square$

**What is not guaranteed for the implemented iteration.** It does not solve $(P_\lambda)$; its output
need not satisfy the equality part of Proposition AA.1 on its support (later subtractions reintroduce
correlation with earlier atoms unless the atoms are orthogonal); there is no rate for $F$; with the cap
active the residual norm need not decrease; and after the first iteration its oracle domain is
restricted by D6. A plausible reading, not tested, is that the 5 km no-reselection rule is what stops an event whose peak
vectors exceed the cap, and which (g) shows is not fully removed, from being selected again. The certificate (d) applies to
any candidate written in the dictionary of (AA.3); the output of the iteration is such a candidate when
$S=0$, and then its suboptimality for the grid problem is bounded at the cost of one exact oracle call. With
$S>0$ its atoms sit at peak-shifted times that are not in that dictionary, and no such bound is available
without a refit.

## AA.6 False-alarm control by the maximum of the field

A false alarm in a window with no event is the declaration of at least one event. By Proposition
AA.2(c), for the BLASSO this happens exactly when $\max_\theta\|P_\theta\mathbf n\|>\lambda$. For the
implemented iteration it happens exactly when $\max_{\theta\in\Theta_h}T(\theta)$ exceeds the threshold,
because the first atom is the global maximum and no exclusion applies at the first step.

**Proposition AA.5 (false-alarm control).**
(a) *Grid, exact bound.* Under (A1), for the grid problem,
$P\big(\max_{j}\|P_{\theta_j}\mathbf n\|^2>t\big)\le n\,P(\chi^2_q>t)$ with $q=\max_jr_j$. For the pilot grid
($n=15\times15\times3\times31=20{,}925$, $q=6$, the rank observed at every sampled node), $\alpha=0.05$
gives $t=36.3$ and $\alpha=0.01$ gives $t=39.9$ ($18.2$ and $19.9$ in the units of the code). These values
apply to the idealized certificate only; the calibrated $\alpha=0.05$ thresholds of the implemented
statistic are about $27$ to $31$ in the units of the code (Appendix Y.17).
(b) *Expected Euler characteristic.* Assume (A1), (A2) and (G) below. Then the excursion set
$A_t=\{\theta\in\Theta:\|P_\theta\mathbf n\|^2\ge t\}$ satisfies

$$\mathbb E\,\chi(A_t)=\sum_{j=0}^{4}\mathcal L_j(\Theta)\,\rho_j(t),\qquad
\rho_0(t)=P(\chi^2_q\ge t),\qquad
\rho_j(t)=\frac{t^{(q-j)/2}e^{-t/2}}{(2\pi)^{j/2}\,2^{(q-2)/2}\,\Gamma(q/2)}\,p_j(t), \tag{AA.7}$$

$$p_1=1,\quad p_2=t-(q-1),\quad p_3=t^2-(2q-1)t+(q-1)(q-2),\quad
p_4=t^3-3q t^2+3(q-1)^2t-(q-1)(q-2)(q-3).$$

The degrees of freedom $\nu$ of [119] are $q$ here. For $q=6$: $p_2=t-5$, $p_3=t^2-11t+20$, $p_4=t^3-18t^2+75t-60$. $\mathcal L_j$ are the Lipschitz-Killing
curvatures of $\Theta$ in the metric $\Lambda(\theta)=\mathrm{Cov}(\nabla Z_l(\theta))$; if $\Lambda$ is constant
and diagonal, $\Lambda=\mathrm{diag}(\ell_1,\dots,\ell_4)$, and $\Theta$ has sides $a_1,\dots,a_4$, then
$\mathcal L_j$ is the $j$-th elementary symmetric polynomial of $(a_1\sqrt{\ell_1},\dots,a_4\sqrt{\ell_4})$
($\mathcal L_0=1$, $\mathcal L_4=\prod_ia_i\sqrt{\ell_i}$).
(c) *Conformal calibration.* Let $M(W)=\max_{\theta\in\Theta_h}T(\theta;W)$ for a window $W$, with the grid,
operators, window rule, cap and whitening rule fixed in advance (a scale estimated inside each window is
allowed). Let $W_0$ be a test window and $W_1,\dots,W_n$ calibration windows, all without events, such
that $(W_0,\dots,W_n)$ is exchangeable. Let $k=\lceil(n+1)(1-\alpha)\rceil$ and $\hat t$ the $k$-th
smallest of $M(W_1),\dots,M(W_n)$ ($\hat t=+\infty$ if $k>n$). Then $P(M(W_0)>\hat t)\le\alpha$, and if
ties among the $M(W_i)$ have probability zero, $P(M(W_0)>\hat t)\ge\alpha-1/(n+1)$. Over any number $n'$
of such test windows the expected number of false alarms is at most $\alpha n'$.

**Assumption (G) for (b).** There exist $C^2$ orthonormal frames $\mathbf f_1(\theta),\dots,\mathbf f_q(\theta)$
of $\mathrm{range}\,G(\theta)$ with $\langle\mathbf f_l(\theta),\mathbf f_{l'}(\theta')\rangle=\delta_{ll'}R(\theta,\theta')$
for all $\theta,\theta'$, so that $Z_l(\theta)=\langle\mathbf f_l(\theta),\mathbf n\rangle$ are independent
copies of one unit-variance Gaussian field and $\|P_\theta\mathbf n\|^2=\sum_lZ_l^2$ is a $\chi^2_q$ field in
the sense of [119]; and the regularity conditions of the Gaussian kinematic formula of [120] hold
(almost surely $C^2$ sample paths, non-degenerate joint law of $Z_l$ and $\nabla Z_l$). [120] states the
formula for a $C^3$ manifold, possibly with boundary, and notes its extension to piecewise smooth domains;
for a box, a piecewise smooth (locally convex, Whitney stratified) parameter space, the setting of [128]
is used.

*Proof.* (a) $\|P_{\theta_j}\mathbf n\|^2\sim\chi^2_{r_j}$ (Proposition AA.2(b)); a union bound, and
$\chi^2_{r+1}$ is stochastically larger than $\chi^2_r$. The numbers are the $1-\alpha/n$ quantiles of
$\chi^2_6$. (b) Under (G), the Gaussian kinematic formula (Theorem 4.1 of [120], for the Euler characteristic; [128]
for a box) gives
$\mathbb E\chi(A_t)=\sum_j\mathcal L_j(\Theta)(2\pi)^{-j/2}\mathcal M_j(D_t)$, where $D_t=\{\mathbf z\in\mathbb R^q:\|\mathbf z\|^2\ge t\}$
and the Gaussian Minkowski functionals $\mathcal M_j$ are the coefficients of
$\gamma_q(\mathrm{Tube}(D_t,\epsilon))=\sum_j\frac{\epsilon^j}{j!}\mathcal M_j(D_t)$. The tube of radius
$\epsilon$ is $\{\|\mathbf z\|\ge\sqrt t-\epsilon\}$, whose Gaussian measure is $\bar F(\sqrt t-\epsilon)$ with
$\bar F$ the tail of the $\chi_q$ distribution and density $f=-\bar F'$,
$f(s)=s^{q-1}e^{-s^2/2}/(2^{(q-2)/2}\Gamma(q/2))$. Hence $\mathcal M_0=\bar F(\sqrt t)$ and
$\mathcal M_j=(-1)^{j-1}f^{(j-1)}(\sqrt t)$ for $j\ge1$. Differentiating
$f$ up to three times gives (AA.7). For $j\le3$ the polynomials are those of [119]; for $q=1$, where
$A_t$ is the union of the excursion sets of $\pm Z_1$ above $\sqrt t$, they reduce to twice the Gaussian
densities $(2\pi)^{-(j+1)/2}H_{j-1}(\sqrt t)e^{-t/2}$, a check that $p_4$ passes. For a constant diagonal
$\Lambda$, the map $\theta_i\mapsto\sqrt{\ell_i}\theta_i$ is an isometry onto a Euclidean box with sides
$a_i\sqrt{\ell_i}$, whose Lipschitz-Killing curvatures are the elementary symmetric polynomials of the
sides. (c) If $k>n$ the event is empty. Otherwise $M_0>\hat t$ implies that at least $k$ of the other $n$
values are strictly smaller than $M_0$. Let $\mathcal A$ be the set of indices $i\in\{0,\dots,n\}$ whose
value exceeds at least $k$ of the others. If $\mathcal A\ne\emptyset$, the index of $\mathcal A$ with smallest
value has at least $k$ strictly smaller values, none of which is in $\mathcal A$, so $|\mathcal A|\le n+1-k$.
By exchangeability each index is in $\mathcal A$ with the same probability, so
$P(0\in\mathcal A)=\mathbb E|\mathcal A|/(n+1)\le(n+1-k)/(n+1)\le\alpha$. Without ties the values have distinct
ranks and $|\mathcal A|=n+1-k$ exactly, which gives $P=(n+1-k)/(n+1)\ge\alpha-1/(n+1)$ because
$k<(n+1)(1-\alpha)+1$. The last statement is linearity of expectation. $\square$

**Approximation (threshold rule).** For high thresholds,
$P(\max_\theta\|P_\theta\mathbf n\|^2\ge t)\approx\mathbb E\chi(A_t)$. For a smooth, centred, unit-variance,
real-valued Gaussian field on a compact piecewise smooth parameter space, [127] proves that the error of
this approximation is exponentially smaller than its leading term as the level grows. That result is for
scalar fields; for the $\chi^2_q$ field the same order of error is expected through the representation of
its maximum as the maximum of a scalar field on the unit-sphere bundle given below, but it is not proved
here and the approximation is used without a bound. The rule is $\lambda(\alpha)=\sqrt{t_\alpha}$, with $t_\alpha$ the largest root of
$\sum_{j=0}^4\mathcal L_j(\Theta)\rho_j(t)=\alpha$, the random-field analogue of choosing the LASSO parameter
as a quantile of $\|X^{\mathsf T}\varepsilon\|_\infty$; in the units of the code the threshold is $t_\alpha/2$.
In the scalar case $q=1$, $Z(\theta)=\langle\mathbf f(\theta),\mathbf n\rangle$ with $\|\mathbf f\|=1$ and
$\Lambda_{ij}=\langle\partial_i\mathbf f,\partial_j\mathbf f\rangle$ is the Fisher information of the unit-energy atom,
so $\mathcal L_4=\int_\Theta\det\Lambda^{1/2}$ counts resolution cells in the metric of the observability
analysis of Section 6.1; better-observed hypothesis spaces require higher thresholds.

**Why the empirical correction is part of the method.** (G) does not hold for moment-tensor subspaces in
general: $Q(\theta)^{\mathsf T}Q(\theta')$ is not a multiple of an orthogonal matrix. The exact identity
$\max_\theta\|P_\theta\mathbf n\|=\max_{\psi\in\mathcal U}\langle\psi,\mathbf n\rangle$, with $\mathcal U$ the set of
unit vectors of $\mathrm{range}\,G(\theta)$, $\theta\in\Theta$ (Cauchy-Schwarz), turns the problem into the
maximum of the canonical Gaussian field on a manifold of dimension $4+q-1=9$, to which [120] applies
with the curvatures of $\mathcal U$; computing them is registered work. More importantly, the implemented
statistic is not a $\chi^2$ field at all: the peak selection takes the maximum of $2S+1$ correlated
samples per arrival, the cap truncates it, and the scale is estimated from the window. (a) and (b)
therefore apply to the idealized certificate only, and (c) is what calibrates the implemented detector.
The measurements of Appendix Y.17 agree: at a single hypothesis, $2T$ has mean $12.6$ to $13.5$ instead of the
$6$ of $\chi^2_6$; the expected Euler characteristic with curvatures estimated from the roughness of the field
overestimates $P(M>u)$ by a factor of $20$ to $60$, because the field decorrelates within one grid step;
and split-conformal thresholds from $300$ noise-only windows gave realized rates of $0.007$, $0.037$, $0.095$
and $0.239$ at $\alpha=0.01$, $0.05$, $0.10$ and $0.25$.
Its guarantee is marginal and per empty window; it needs exchangeable noise windows (not guaranteed for
consecutive windows of a non-stationary record) that contain no events (contamination by uncatalogued
events is expected to make it conservative, which is a heuristic, not a proof). It needs $n\ge1/\alpha-1$
windows for a finite threshold; below that, the tail of the calibration maxima is extrapolated by
extreme-value fitting as in [125], an **approximation**. It does not control false events in windows that
contain events (sidelobes, extra atoms), which the FROC endpoint of hypothesis G1 (Section 4.3) measures
and no result here bounds.

## AA.7 Observability and exact recovery

**Proposition AA.6 (observability and exact recovery). (a) Exact support recovery on a grid.** Consider (AA.3) with
$\mathbf d=\sum_{j\in\mathcal S}Q_j\mathbf a^0_j+\mathbf w$, $\mathbf a^0_j\ne0$, $s=|\mathcal S|$,
$a_{\min}=\min_{j\in\mathcal S}\|\mathbf a^0_j\|$. Assume that $\Gamma=[Q_j]_{j\in\mathcal S}$ has full column
rank with $\sigma=\sigma_{\min}(\Gamma)>0$, and define the pre-certificate
$\mathbf p_F=\Gamma(\Gamma^{\mathsf T}\Gamma)^{-1}\mathbf u^0$, $\mathbf u^0=(\mathbf a^0_j/\|\mathbf a^0_j\|)_{j\in\mathcal S}$,
and the margin $\zeta=1-\max_{j\notin\mathcal S}\|Q_j^{\mathsf T}\mathbf p_F\|$. If $\zeta>0$,
$\|\mathbf w\|\le\zeta\lambda/2$ and $\lambda<\zeta\sigma^3a_{\min}/(6\sqrt s)$, then (AA.3) has a unique minimizer;
its support is exactly $\mathcal S$; and $\|\hat{\mathbf a}_{\mathcal S}-\mathbf a^0_{\mathcal S}\|\le3\sqrt s\,\lambda/(2\sigma^2)$.

*Proof.* Note $\sigma\le1$ (a unit vector inside one block is mapped to a unit vector) and
$\|\Gamma\|\le\sqrt s$ (because $\|\Gamma\mathbf x\|\le\sum_j\|\mathbf x_j\|$). Let $\epsilon=\|\mathbf w\|$.
The problem restricted to $\mathcal S$ is strictly convex (full column rank), so it has a unique minimizer
$\hat{\mathbf a}_{\mathcal S}$, characterized by $\Gamma^{\mathsf T}(\mathbf d-\Gamma\hat{\mathbf a}_{\mathcal S})=\lambda\mathbf z$
with $\mathbf z_j$ in the subdifferential of $\|\cdot\|$ at $\hat{\mathbf a}_j$, so $\|\mathbf z\|\le\sqrt s$.
Hence $\hat{\mathbf a}_{\mathcal S}-\mathbf a^0_{\mathcal S}=(\Gamma^{\mathsf T}\Gamma)^{-1}(\Gamma^{\mathsf T}\mathbf w-\lambda\mathbf z)$ and
$\|\hat{\mathbf a}_{\mathcal S}-\mathbf a^0_{\mathcal S}\|\le\sqrt s(\epsilon+\lambda)/\sigma^2\le3\sqrt s\lambda/(2\sigma^2)<\zeta\sigma a_{\min}/4\le a_{\min}/4$.
So every $\hat{\mathbf a}_j\ne0$ and $\mathbf z_j=\hat{\mathbf a}_j/\|\hat{\mathbf a}_j\|$. For nonzero vectors,
$\|\mathbf x/\|\mathbf x\|-\mathbf y/\|\mathbf y\|\|\le2\|\mathbf x-\mathbf y\|/\|\mathbf y\|$, so
$\|\mathbf z-\mathbf u^0\|\le2\|\hat{\mathbf a}_{\mathcal S}-\mathbf a^0_{\mathcal S}\|/a_{\min}\le3\sqrt s\lambda/(\sigma^2a_{\min})$.
Extend $\hat{\mathbf a}$ by zero off $\mathcal S$. Its residual is
$\hat{\mathbf r}=\mathbf w-\Gamma(\hat{\mathbf a}_{\mathcal S}-\mathbf a^0_{\mathcal S})=(I-\Gamma(\Gamma^{\mathsf T}\Gamma)^{-1}\Gamma^{\mathsf T})\mathbf w+\lambda\Gamma(\Gamma^{\mathsf T}\Gamma)^{-1}\mathbf z$.
For $j\notin\mathcal S$, using $\|\Gamma(\Gamma^{\mathsf T}\Gamma)^{-1}\|=1/\sigma$,
$\|Q_j^{\mathsf T}\hat{\mathbf r}\|/\lambda\le\epsilon/\lambda+\|Q_j^{\mathsf T}\mathbf p_F\|+\|\mathbf z-\mathbf u^0\|/\sigma\le\zeta/2+(1-\zeta)+3\sqrt s\lambda/(\sigma^3a_{\min})<1$.
On $\mathcal S$, $Q_j^{\mathsf T}\hat{\mathbf r}=\lambda\mathbf z_j$ by construction. By Proposition AA.1(c),
$\hat{\mathbf a}$ is a minimizer. Every minimizer has the same residual, hence support in
$\{j:\|Q_j^{\mathsf T}\hat{\mathbf r}\|=\lambda\}\subset\mathcal S$, hence equals $\hat{\mathbf a}$ because
$\Gamma$ is injective. $\square$

This is the grid form of the non-degenerate source condition: $\zeta>0$ is strict dual feasibility of the
minimal-norm pre-certificate off the support. For one source, $\Gamma=Q_1$, $\sigma=1$ and
$\mathbf p_F=\hat{\mathbf g}=\mathbf g/\|\mathbf g\|$ with $\mathbf g=G(\theta_1)\mathbf m_1$, so
$\zeta=1-\max_{j\ne1}\|P_{\theta_j}\hat{\mathbf g}\|$: the margin is how far the normalized noise-free signal
is from the subspace of every other node. Parts (b) to (d) compute that quantity near the source and
identify it with the observability certificate.

**Lemma AA.4 (derivative of the projector).** If $G\in C^1$ has full column rank, then
$P=G(G^{\mathsf T}G)^{-1}G^{\mathsf T}$ is $C^1$ and $dP=(I-P)\,dG\,G^++(G^+)^{\mathsf T}dG^{\mathsf T}(I-P)$.

*Proof.* Differentiability follows from the formula. From $PG=G$, $dP\,G=(I-P)dG$, so $dP\,P=(I-P)dG\,G^+$
(as $GG^+=P$). Transposing, $P\,dP=(G^+)^{\mathsf T}dG^{\mathsf T}(I-P)$. From $P^2=P$, $dP=dP\,P+P\,dP$. $\square$

**Proposition AA.6, continued (single source: the curvature of the certificate is the observability with
the mechanism as nuisance).** Assume (A3), (A4) with $c=0$ (in continuous time), $\theta_1$ in the interior of $\Theta$, $G\in C^2$ with
rank $6$ near $\theta_1$,
a single source $\mathbf g=G(\theta_1)\mathbf m_1$ with energy $E=\|\mathbf g\|^2$, and
$\eta_V(\theta)=\|P_\theta\hat{\mathbf g}\|$.
(b) $\eta_V$ is the group analogue, introduced here, of the vanishing-derivative pre-certificate of [47]:
$\hat{\mathbf g}$ is the
minimal-norm $\mathbf p$ with $P_{\theta_1}\mathbf p=\hat{\mathbf g}$ and $\nabla\|P_\theta\mathbf p\|^2=0$ at
$\theta_1$. Moreover $\eta_V\le1$, and $\eta_V(\theta)=1$ iff $\mathbf g\in\mathrm{range}\,G(\theta)$.
(c) $\nabla\eta_V(\theta_1)=0$ and $\nabla^2\eta_V(\theta_1)=-F_{\mathrm{eff}}/E$, where
$F_{\mathrm{eff}}=J^{\mathsf T}(I-P_{\theta_1})J$, $J=[\partial_{\theta_j}(G(\theta)\mathbf m_1)]_{j=1}^4$, is
the efficient Fisher information for $\theta$ when the moment tensor is a nuisance parameter.
(d) With $E_{ib}=\|A_{ib}(\mathbf x_1)\mathbf m_1\|^2$, $\mathbf v_{ib}=\nabla_\theta t_{ib}(\theta_1)$ (the vector
$\mathbf v^s_b=(-\mathbf p^s_b,1)$ of (8a) and (8d) for $s=i$), and $\partial A\,\mathbf m_1=[\partial_{\theta_j}A(\mathbf x)\mathbf m_1]_j$,

$$F_{\mathrm{eff}}=\|\dot{\bar s}\|^2\sum_{i,b}E_{ib}\,\mathbf v_{ib}\mathbf v_{ib}^{\mathsf T}+(\partial A\,\mathbf m_1)^{\mathsf T}(I-\Pi_{\mathbf x_1})(\partial A\,\mathbf m_1)
\ \succeq\ \|\dot{\bar s}\|^2\sum_{i,b}E_{ib}\,\mathbf v_{ib}\mathbf v_{ib}^{\mathsf T}. \tag{AA.8}$$

*Proof.* (b) Any $\mathbf p$ with $P_{\theta_1}\mathbf p=\hat{\mathbf g}$ has $\|\mathbf p\|^2=1+\|(I-P_{\theta_1})\mathbf p\|^2$,
so $\hat{\mathbf g}$ is the minimal-norm one; the derivative condition holds for it by (c). $\|P\hat{\mathbf g}\|=1$
iff $\hat{\mathbf g}\in\mathrm{range}\,P$. (c) Let $\mathbf e(\theta)=(I-P_\theta)\hat{\mathbf g}$, so
$\eta_V^2=1-\|\mathbf e\|^2$ and $\mathbf e(\theta_1)=0$. Then $\nabla\eta_V^2(\theta_1)=0$ and
$\nabla^2\eta_V^2(\theta_1)=-2(\partial\mathbf e)^{\mathsf T}(\partial\mathbf e)$, the term with
$\langle\mathbf e,\partial^2\mathbf e\rangle$ vanishing at $\theta_1$. By Lemma AA.4 and $(I-P_{\theta_1})\hat{\mathbf g}=0$,
$\partial_j\mathbf e=-(I-P)\partial_jG\,G^+\hat{\mathbf g}=-(I-P)\partial_j(G\mathbf m_1)/\sqrt E$, since
$G^+\hat{\mathbf g}=\mathbf m_1/\sqrt E$. So $\nabla^2\eta_V^2=-2F_{\mathrm{eff}}/E$; as $\eta_V(\theta_1)=1$ and
$\nabla\eta_V(\theta_1)=0$, $\nabla^2\eta_V=\frac12\nabla^2\eta_V^2$. The Schur complement of the
$\mathbf m$ block in the Fisher information $[J,G]^{\mathsf T}[J,G]$ is $J^{\mathsf T}(I-P)J$, which is the
definition of the efficient information. (d) From (AA.4),
$\partial_j(G\mathbf m_1)=U(\partial_jA)\mathbf m_1-\mathbf k_j$ with
$\mathbf k_j=\sum_{ib}\mathbf e_i\otimes(A_{ib}\mathbf m_1)\,\partial_jt_{ib}\otimes\dot{\bar s}(\cdot-t_{ib})$. The vector
$\mathbf k_j$ is orthogonal to $\mathrm{range}\,U$: within an arrival because $\langle\bar s,\dot{\bar s}\rangle=\frac12[\bar s^2]_{-\infty}^{\infty}=0$,
and across P and S by (A4). By Lemma AA.2, $P=U\Pi U^*$, so $(I-P)\mathbf k_j=\mathbf k_j$ and
$(I-P)U\mathbf z=U(I-\Pi)\mathbf z$. The two parts are orthogonal, and
$\langle\mathbf k_j,\mathbf k_l\rangle=\|\dot{\bar s}\|^2\sum_{ib}E_{ib}\partial_jt_{ib}\partial_lt_{ib}$ by (A4). The
second term of (AA.8) is positive semidefinite. $\square$

The observability certificate of Section 6.1 was derived with a known mechanism and weights
$A_b^2/d_s^2$. (AA.8) shows that its leading term, with the weights replaced by the energies that the actual
mechanism radiates to each station and branch, is a lower bound on the curvature of the certificate at a
true source even when the mechanism is unknown, because the time-derivative directions are orthogonal to the
whole amplitude subspace. Consequences: a station on a nodal plane of a branch contributes no curvature
through that branch; Corollaries 1 to 3 of Section 6.1 apply to the lower bound in (AA.8) with the captured
arrivals replaced by those with $E_{ib}>0$; and
since $F_{\mathrm{eff}}\preceq J^{\mathsf T}J$, the known-mechanism information is an upper bound. On a grid,
$\eta_V(\theta_1+\mathbf h)=1-\frac12\mathbf h^{\mathsf T}F_{\mathrm{eff}}\mathbf h/E+o(|\mathbf h|^2)$, so the margin
of Proposition AA.6(a) for one on-grid source is at most $\frac12\min\mathbf h^{\mathsf T}F_{\mathrm{eff}}\mathbf h/E+o(h^2)$
over the neighbour offsets: it is set by the observability and shrinks quadratically with the mesh.

**Conjecture AA.1 (group non-degenerate recovery over $\Theta$).** Let $\mathbf d=\sum_{k=1}^KG(\theta_k)\mathbf m_k+\mathbf w$
with $\theta_k$ in the interior of $\Theta$, $G\in C^3$ with rank 6. Let
$\hat{\mathbf g}_k=G(\theta_k)\mathbf m_k/\|G(\theta_k)\mathbf m_k\|$, $\bar{\mathbf m}_k=G(\theta_k)^+\hat{\mathbf g}_k$, and
$\Gamma_V=[\,Q(\theta_k),\ (I-P_{\theta_k})\partial_{\theta_j}(G(\theta)\bar{\mathbf m}_k)|_{\theta_k}\,]_{k\le K,\ j\le4}$. Assume
$\Gamma_V$ has full column rank, and let $\mathbf p_V$ be the minimal-norm solution of
$Q(\theta_k)^{\mathsf T}\mathbf p=Q(\theta_k)^{\mathsf T}\hat{\mathbf g}_k$ and
$\langle(I-P_{\theta_k})\partial_j(G\bar{\mathbf m}_k),\mathbf p\rangle=0$ (the conditions $P_{\theta_k}\mathbf p=\hat{\mathbf g}_k$
and $\nabla\|P_\theta\mathbf p\|^2(\theta_k)=0$, linearized by Lemma AA.4). Assume the non-degenerate source
condition: $\eta_V(\theta)=\|P_\theta\mathbf p_V\|<1$ for $\theta\notin\{\theta_k\}$ and $\nabla^2\eta_V(\theta_k)\prec0$. Then
there are $\alpha,\lambda_0,C>0$ such that for $\lambda<\lambda_0$ and $\|\mathbf w\|\le\alpha\lambda$, $(P_\lambda)$ has a
unique solution, it has exactly $K$ atoms, and they satisfy
$|\hat\theta_k-\theta_k|+\|\hat{\mathbf m}_k-\mathbf m_k\|\le C(\lambda+\|\mathbf w\|)$.

*Status.* For real scalar amplitudes on the one-dimensional torus this is Theorem 2 of [47]; Section 3.5 of
[47] states that the proof extends in a straightforward manner to the $d$-dimensional torus, with the
second derivative of the certificate replaced by its Hessian, which must be invertible, but does not write
that extension out. The multidimensional BLASSO with scalar amplitudes, in which the relative position of
the spikes enters the stability analysis, is treated in [140]; it is cited as the nearest multidimensional
analysis, not as support for the conjecture. The extension to $\mathbb R^6$-valued amplitudes with the group norm, and to a
four-dimensional box $\Theta$ rather than a torus, is not proved here. For
$K=1$, Proposition AA.6(b)-(d) shows that the hypotheses reduce to global identifiability of the noise-free signal
($\mathbf g\notin\mathrm{range}\,G(\theta)$ for $\theta\ne\theta_1$) and $F_{\mathrm{eff}}\succ0$, which (AA.8)
guarantees whenever $\sum_{ib}E_{ib}\mathbf v_{ib}\mathbf v_{ib}^{\mathsf T}\succ0$. **Conditional corollary:** if
Conjecture AA.1 holds, a single source is recovered with exactly one atom, for small noise and $\lambda$
proportional to the noise level, wherever the energy-weighted observability certificate is positive and the
source is globally identifiable. For $K\ge2$ the condition involves cross terms between sources; it is
computable for any configuration (form $\mathbf p_V$, evaluate $\eta_V$ with the exact oracle) and that
computation, not a closed form, is what would predict where two events merge. That well-separated sources
satisfy it whenever each satisfies the single-source condition with margin is plausible, because the
cross-Gram terms decay with the separation of arrival times at the pulse's rate, but it is a conjecture.

## AA.8 What the theory does not cover

The results assume that the operator is right. Model error (velocity, statics, amplitude factors,
polarization deviations, the D7 gap between generator and operator) enters as a deterministic term in
$\mathbf w$ that is not small relative to $\lambda$ at high frequency, and neither Proposition AA.6(a) nor
Conjecture AA.1 covers it; the cap (D2), which AA.2 to AA.7 treat only as a deviation and AA.9 treats directly, bounds the influence
of each station and phase block and is what protects the implemented statistic against amplitude model
error in dense sequences (Appendices Y.24, Y.25); the measured degradation in the "pert" condition of Y.11 is evidence, not
theory. Statics larger than the search half-width $w$ break the peak-vector model entirely (Y.11); the
theory only says that within $\pm w$ the statistic is a lower bound on a shift-maximized GLRT. Real noise is
non-Gaussian, non-stationary and correlated between stations, and contains coherent clutter; (A1) fails, the
whitening by one scale is not whitening, and only the conformal statement of Proposition AA.5(c) survives,
with exchangeability of windows as a genuine assumption. Discretization: the recovery results are for sources
on the grid (Proposition AA.6(a)) or conditional (Conjecture AA.1); the pilot places no source on a node, and for
an off-grid source the modelling error is of order $h$ while the admissible noise in Proposition AA.6(a) is of
order $\lambda h^2$ (Proposition AA.6(d) and the remark after it), so exact grid support recovery is not expected on fine grids (neighbouring nodes share the
mass); the sliding step is the registered remedy. Task 2C measured the cost of discretization on the
grid of Appendix Y.18: the certified grid solution fell $0.06$ to $0.61$ below MT-FW in recall in 23
of 24 cases, and with sources moved onto the nodes the exact-time block pursuit became the best method
(Appendix Y.20, exploratory); given the peak search of MT-FW, the grid solution still fell $0.11$ to
$0.30$ below MT-FW in 11 of 12 cases, and it is certified only for the dictionary that the heuristic
peak oracle builds, not over all time shifts (Appendix Y.21); and the sliding solver with a time shift
per station and phase gained $+0.000$ to $+0.039$ over block pursuit with equal tolerance, every interval
containing zero (Appendix Y.22). The implemented iteration is covered only by Proposition
AA.4(e) to (h); its equivalence to a minimizer of $(P_\lambda)$, and every rate, are not claimed, and the
a posteriori certificate AA.4(d) (with $S=0$, or after a refit on the dictionary of (AA.3)) is the registered
way to measure how far it is from one. Finally, nothing
here bounds false events in windows that contain events, which is what the FROC endpoint counts.

## AA.9 The deployed statistics: influence of one station and saturation

The results of AA.2 to AA.7 concern the convex problem $(P_\lambda)$, which Task 2E refuted as a detector. The
statistics that were confirmed or tested on real data are of a simpler kind: a radial transform of each station and
phase block of peak vectors, followed by the projection on the moment-tensor subspace. This section states what can be
proved about them with elementary tools, and connects it to the measurements of Tasks 2F, 2I, 2M and exploration C18.

**Setting.** At a hypothesis $h=(\mathbf x,\tau)$ the detector reads, for each phase $b\in\{P,S\}$ and station
$i=1,\dots,n_s$, a peak vector $\mathbf v_k\in\mathbb R^3$, $k=(b,i)$, of the normalized matched-filter output (the
sample of largest energy within $\pm W$ of the predicted time). Stack them into $\mathbf v\in\mathbb R^{6n_s}$. Let
$Q_h\in\mathbb R^{6n_s\times 6}$ have orthonormal or zero columns spanning the radiation subspace of the node (the left
singular vectors of the radiation operator, `c9_lib.bases`). For a map
$\phi:\mathbb R^3\to\mathbb R^3$ applied block by block, the score is
$$s_\phi(h)=\tfrac12\,\|Q_h^{\mathsf T}\phi(\mathbf v)\|^2,\qquad \phi(\mathbf v)=(\phi(\mathbf v_1),\dots,\phi(\mathbf v_{2n_s})).$$
Three radial maps are used, each of the form $\phi(\mathbf x)=r(\|\mathbf x\|)\,\mathbf x/\|\mathbf x\|$:
the identity $r(t)=t$ (the energies of BOMP-PK, "uncapped"); the cap $r_c(t)=\min(t,c)$ (MT-FW, $c=6$ noise
deviations); and the logarithmic compression $\ell_c(t)=c\sqrt{\log(1+t^2/c^2)}$ (LOG6, $c=6$, exploration C18). The
Huber statistic HUBER-N is an M-estimator fit, not a radial map, and is not covered here.

**Lemma AA.9.1 (contraction).** For any $\mathbf a,\mathbf b\in\mathbb R^{6n_s}$,
$\big|\,\|Q_h^{\mathsf T}\mathbf a\|-\|Q_h^{\mathsf T}\mathbf b\|\,\big|\le\|\mathbf a-\mathbf b\|$, and
$\|Q_h^{\mathsf T}\mathbf a\|\le\|\mathbf a\|$.

*Proof.* $Q_h^{\mathsf T}Q_h$ is diagonal with entries 0 or 1, so $Q_hQ_h^{\mathsf T}$ is an orthogonal projector and
$\|Q_h^{\mathsf T}\mathbf a\|^2=\mathbf a^{\mathsf T}Q_hQ_h^{\mathsf T}\mathbf a\le\|\mathbf a\|^2$. Apply this to
$\mathbf a-\mathbf b$ and use the reverse triangle inequality. $\square$

**Proposition AA.9.2 (influence of one station and phase).** Let $\mathbf v$ and $\mathbf v'$ differ only in block
$k$, and write $\sigma_\phi=\sqrt{2s_\phi}$. Then
$$|\sigma_\phi(\mathbf v)-\sigma_\phi(\mathbf v')|\le r(\|\mathbf v_k\|)+r(\|\mathbf v_k'\|).$$
Hence (a) for the cap the change is at most $2c$, whatever the amplitude of the block; (b) for the logarithmic map,
if both blocks have norm at most $M$, the change is at most $2c\sqrt{\log(1+M^2/c^2)}$; (c) for the identity the
change is not bounded.

*Proof.* By Lemma AA.9.1, $|\sigma_\phi(\mathbf v)-\sigma_\phi(\mathbf v')|\le\|\phi(\mathbf v)-\phi(\mathbf v')\|
=\|\phi(\mathbf v_k)-\phi(\mathbf v'_k)\|\le\|\phi(\mathbf v_k)\|+\|\phi(\mathbf v'_k)\|$, and
$\|\phi(\mathbf x)\|=r(\|\mathbf x\|)$. Parts (a) and (b) follow from $r_c\le c$ and the monotonicity of
$\ell_c$. For (c), let $\mathbf v'_k=0$ and $\mathbf v=\mathbf v'+\mathbf e$, with $\mathbf e$ zero outside block $k$ and
$\|Q_h^{\mathsf T}\hat{\mathbf e}\|=\gamma>0$ for its direction $\hat{\mathbf e}$ (possible whenever block $k$ is not
orthogonal to the subspace). Then $\sigma(\mathbf v)-\sigma(\mathbf v')\ge\gamma\|\mathbf e\|-2\|Q_h^{\mathsf T}\mathbf v'\|$,
which grows linearly with $\|\mathbf e\|$. $\square$

*Numbers.* A digital or cultural transient at one station of $M=12{,}500$ noise deviations, the largest the
verification of Task 2M found, can move $\sigma$ by at most $12$ with the cap, by at most
$2\cdot6\sqrt{\log(1+(12500/6)^2)}\approx47$ with the logarithmic map, and by up to $12{,}500$ without a bound. A
real event at the detection threshold of Task 2M has $\sigma$ of order $\sqrt{2\cdot77}\approx12$ (the pooled
threshold of the capped detector at $0.25$ false events per window was $77.4$ and $74.9$ in score units on the two seeds), so one such
transient decides the ranking of an uncapped hypothesis and cannot decide that of a capped one by itself. This is the
mechanism behind Tasks 2F, 2I and 2M and behind the collapse of BOMP-PK in the daytime hours of the 2M primary day
(Appendix Y.31): the realized false events of BOMP-PK there come from one or two stations at a time.

**Proposition AA.9.3 (saturation of the capped score).** For every hypothesis $s_{r_c}(h)\le n_bc^2/2$, with
$n_b=2n_s$ the number of blocks. If every block read at $h$ has $\|\mathbf v_k\|\ge c$, then
$$s_{r_c}(h)=\tfrac{c^2}2\,\|Q_h^{\mathsf T}\mathbf u\|^2,\qquad \mathbf u=(\mathbf v_1/\|\mathbf v_1\|,\dots,
\mathbf v_{n_b}/\|\mathbf v_{n_b}\|),$$
which does not depend on the amplitudes. Consequently, if a background hypothesis $h'$ (no event at $h'$) has all its
blocks at or above the cap and directions with $\|Q_{h'}^{\mathsf T}\mathbf u'\|\ge\|Q_h^{\mathsf T}\mathbf u\|$,
then $s_{r_c}(h')\ge s_{r_c}(h)$ for an event at $h$ of any magnitude whose blocks are also saturated.

*Proof.* $\|r_c(\mathbf v)\|^2=\sum_k\min(\|\mathbf v_k\|,c)^2\le n_bc^2$, and Lemma AA.9.1 gives the bound. When
every block is saturated, $r_c(\mathbf v_k)=c\,\mathbf v_k/\|\mathbf v_k\|$, so $r_c(\mathbf v)=c\mathbf u$. The
comparison follows from the formula. $\square$

In a quiet window the bound is not reached: only the blocks of a real event exceed the cap, and its score separates
from the background. In a continuous aftershock sequence, where coda and small events put energy above the cap at
most stations in most windows, every hypothesis is scored by the directions of its peak vectors alone, and a
magnitude 2.5 event has no amplitude advantage over the many saturated background hypotheses that the grid offers.
Exploration C18 measured exactly this on the development segment of the Ridgecrest hour: at the true hypocentre the
capped score of a real event did not exceed the maxima of event-free windows (AUC 0.39 to 0.49), while the uncapped
score did (0.69 to 0.74).

**Proposition AA.9.4 (the logarithmic map keeps the amplitude ordering).** Let an event contribute blocks
$\mathbf v_k=a\,\mathbf w_k$, $a>0$, with $\|\mathbf w_k\|=\omega$ for all $k$. Then
$$s_{\ell_c}(a)=\tfrac{c^2}2\log\!\Big(1+\frac{a^2\omega^2}{c^2}\Big)\,\|Q_h^{\mathsf T}\hat{\mathbf w}\|^2,\qquad
\hat{\mathbf w}=\mathbf w/\omega,$$
which is strictly increasing and unbounded in $a$ whenever $Q_h^{\mathsf T}\hat{\mathbf w}\ne0$, whereas
$s_{r_c}(a)$ is constant for $a\ge c/\omega$.

*Proof.* With equal block norms the radial map multiplies every block by the same factor
$\ell_c(a\omega)/(a\omega)$, so $\ell_c(\mathbf v)=\ell_c(a\omega)\hat{\mathbf w}$ and the formula follows; the cap
gives $r_c(a\omega)\hat{\mathbf w}=c\hat{\mathbf w}$ for $a\omega\ge c$. $\square$

The equal-norm hypothesis is restrictive (real events have unequal amplitudes across stations), and nothing here
orders the logarithmic and capped scores for general inputs, because a block-diagonal rescaling need not change
$\|Q_h^{\mathsf T}\cdot\|$ monotonically. The propositions describe a trade-off rather than an optimum: a bounded
transform protects against transients (AA.9.2) and saturates in dense seismicity (AA.9.3); a transform with slowly
growing, unbounded radius keeps part of each property. Exploration C18 found that LOG6 kept the robustness of the cap
to real-noise transients within $-0.015$ [$-0.099$, $+0.043$] and restored most of the ranking of real events, but
neither statistic located the Ridgecrest events with the registered timing; with station terms read over +-2 s and a
layered medium the uncapped statistic located up to 5 of the 16 development events and the capped one none (C18j).
The false-alarm tail of these scores under the peak search is treated next.

**Proposition AA.9.5 (null tail under the peak search, and contamination).** Assume that, at the samples a hypothesis
can read, the normalized matched-filter output is the white Gaussian noise of unit variance passed through the
normalized kernel $k$ (so each component is a stationary Gaussian sequence of unit variance), that different blocks
are independent (different stations, or P and S of one station farther apart than the kernel support), and let
$\kappa=\max_f|\hat k(f)|/\|k\|_2$. Each block is read at the sample of largest energy among $L=2W+1$ candidates.
Write $\sigma_\phi=\sqrt{2s_\phi}$ and let $\phi$ be the identity, the cap or the logarithmic map.
(a) The capped and uncapped scores coincide unless some candidate sample of some block exceeds $c$ in norm, an event
of probability at most $n_bL\,\Pr(\chi^2_3>c^2)$.
(b) For every $p\in(0,1)$, $\Pr\big(\sigma_\phi\ge\sqrt6+\kappa\sqrt{2n_b\log L}+\kappa\sqrt{2\log(1/p)}\big)\le p$
for one hypothesis; over $N$ hypotheses, replacing $p$ by $p/N$ bounds the probability that any of them crosses.
(c) If at most $m$ blocks are replaced by arbitrary values (transients), the capped score satisfies
$\sigma_{r_c}\le\sigma^{\rm clean}_{r_c}+2c\sqrt m$, so the threshold of (b) raised by $2c\sqrt m$ keeps its
false-alarm probability under any such contamination; no finite threshold does so for the identity.

*Proof.* (a) If every candidate has norm at most $c$ the cap acts as the identity on every block; a union bound over
the $n_bL$ candidates, each with squared norm $\chi^2_3$, gives the probability. (b) For a fixed choice $s$ of one
candidate per block, $f_s=\|Q_h^{\mathsf T}\phi(\mathbf x_s)\|$ is a composition of the map
$\mathbf z\mapsto\mathbf x$ (convolution with $k/\|k\|_2$, operator norm $\kappa$), the block-wise radial map (1-Lipschitz:
for a radial map the constant is $\max(\sup r',\sup r(t)/t)$, which is 1 for the identity and the cap, and for $\ell_c$
since $\ell_c'(t)\le(1+t^2/c^2)^{-1/2}$ from $\log(1+u)\ge u/(1+u)$) and the contraction of Lemma AA.9.1; so $f_s$ is
$\kappa$-Lipschitz in the white noise $\mathbf z$, and the read score satisfies $\sigma_\phi\le F=\max_sf_s$, which
is $\kappa$-Lipschitz as a maximum of such functions. Gaussian concentration gives
$\Pr(F\ge\mathbb EF+t)\le e^{-t^2/2\kappa^2}$. Each $\phi(\mathbf x_{k,w})$ is centred and isotropic with
$\mathbb E\|\phi(\mathbf x_{k,w})\|^2\le3$, and the blocks are independent, so the covariance of $\phi(\mathbf x_s)$ is
block-diagonal with blocks $\preceq I_3$ and $\mathbb Ef_s^2\le\operatorname{tr}(Q_h^{\mathsf T}Q_h)\le6$; each $f_s$ is
sub-Gaussian around its mean with proxy $\kappa^2$, and the maximum of the $L^{n_b}$ of them has
$\mathbb EF\le\sqrt6+\kappa\sqrt{2n_b\log L}$. (c) The blocks that differ contribute
$\|r_c(\mathbf v_k)-r_c(\mathbf v_k')\|\le2c$ each, and Lemma AA.9.1 gives the bound; for the identity,
Proposition AA.9.2(c). $\square$

*Numbers* (`code/theory_tail_check.py`, `results/theory/tail_check_out.txt`). For the 2 Hz kernel at 16 Hz,
$\kappa=2.146$. With seven stations ($n_b=14$), $W=4$ samples and $c=6$, the probability in (a) is
$9.4\times10^{-6}$ per hypothesis ($\Pr(\chi^2_3>36)=7.5\times10^{-8}$), and $3.4\times10^{-4}$ that any of the
$7\times640$ samples of a 40 s window exceeds the cap: under Gaussian noise the cap is inactive in all but about one
window in three thousand, so it costs nothing in false alarms and its effects are those of (c) and of Proposition
AA.9.3. Over $2\times10^5$ simulated hypotheses the upper quantiles of $\sigma$ at $p=10^{-2}$ and $10^{-4}$ were 5.81
and 7.18 for the identity and the cap (identical, as (a) predicts) and 5.51 and 6.78 for the logarithmic map, against
the bound of (b) of 25.8 and 28.5. The bound holds and is conservative by a factor of about four in $\sigma$, most of it
from the union over the $L^{n_b}$ read positions; it is a guarantee, not a calibration, and an operating threshold is
still better set on noise windows of the same network. What (b) and (c) add is the form of the tail (Gaussian with
scale $\kappa$ around a mean) and the fact that, for the capped score, a threshold set on clean noise can be made safe
against transients at a known number of arrivals by an additive margin that does not depend on their amplitude; for the
logarithmic score the margin grows only as $\sqrt{\log}$ of the amplitude (Proposition AA.9.2(b)). Real noise is neither Gaussian nor white, and the peak search, the grid and the greedy deflation of MT-FW
make the null law of the deployed detector a maximum of dependent quantities that (b) bounds only through $N$.
