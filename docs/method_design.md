# Method design: certified detection and location of earthquakes as sparse recovery of a tensor-valued measure

Thesis design document, September 2026. It complements the proposal
(`Doctoral_Proposal.md`, Sections 7 and 11 and
Appendices Y, Z and AA) and builds on the results of Tasks 1j to 1z, 2A to 2F, 2I and 2L and of the explorations of the
calibration seed (C8 to C14).

---

## 1. The open problem, stated precisely

A seismic catalogue answers three questions for each recording window: how many events there are, where and
when they occurred. In sparse, noisy networks with small events, the standard answer is a
three-stage chain: a per-station phase detector (PhaseNet, EQTransformer), an associator that
groups those phases into events (GaMMA, PyOcto, REAL, GENIE, Harpa) and a locator. Each stage makes
hard decisions that the next one cannot undo.

The experiments of this thesis measured where the information is lost. In the consistent synthetic
generator (Appendix Y), the phase detector delivers fewer than half of the true arrivals.
With clutter treated by the same rule as the signal, GaMMA recovers between 0.28 and 0.63 of the events
depending on the network (Y.12). A physical filter applied afterwards, on the GaMMA catalogue, separates the true
events from the false ones, but it cannot recover those that GaMMA did not form (Y.7). The recent evaluation
of associators (Puente Huerta et al., Seismica, 2025) gives them F1 greater than 0.8 on clean synthetic
phases. The bottleneck is not association itself, but the information that the phase detector
discards before it and the transients that it marks as phases.

The existing families that avoid that loss each have a limit:

- **Migration of characteristic functions** (QuakeMigrate, MALMI, Qseek). It stacks a scalar per station,
  either STA/LTA or a learned probability, and discards the sign and the direction of motion.
- **End-to-end networks** (EQNet, PhaseNO, HEIMDALL). They depend on training data from the
  region, give no guarantees, and their false-alarm rate is known only empirically.
- **Moment-tensor migration** (Chambers et al., 2014; Anikiev et al., 2014; Zhebel and Eisner,
  2015; Staněk et al., 2015; Vlček et al., 2016). It exploits the coherence of the radiation pattern, but
  in dense surface arrays for microseismicity.
- **GRiD MT** (Tsuruoka et al., 2009; Guilhem Trilla and Cano, 2024). It scans location and time with a
  linear tensor inversion at each node in sparse regional networks, but at long periods
  (15 to 50 s) and for moderate events.
- **Block-sparse inversion** (Vera Rodriguez, Sacchi and Gu, 2012). It jointly recovers origin
  time, hypocentre and moment tensor of simultaneous sources with a redundant dictionary of
  Green's functions and a block orthogonal matching pursuit (block OMP). It is the closest antecedent of
  this method: joint, sparse, moment-tensor-aware detection and location over a
  grid of hypotheses.

None of these families solves the search with a guarantee of not missing the best hypothesis, and
none calibrates with a guarantee the false-alarm rate of the detector maximised over the whole hypothesis
space.

**The open problem the thesis attacks** is to detect, count and locate small events in sparse and noisy
networks at a false-alarm rate calibrated with a finite-sample guarantee. It is attacked by using all the physics that
the medium allows before the phase detector discards information, and with mathematical guarantees
on the search and on the threshold.

## 2. The idea: a single convex problem instead of three hard decisions

The proposal already defined the correct unknown (equation 11): a sparse measure with tensor-valued
amplitudes over the hypothesis space $\Theta=\{(\mathbf x,\tau)\}$,

$$\mu=\sum_{k=1}^{K}\mathbf m_k\,\delta_{\theta_k},\qquad \mathbf m_k\in\mathbb R^6,\ \theta_k\in\Theta,$$

where $\mathbf m_k$ is the moment tensor of event $k$ and $K$ is unknown. Physics says that the
three-component recording of the whole network is **linear** in $\mu$:
$\mathbf d=\Phi\mu+\mathbf n$, with $\Phi\mu=\sum_k G(\theta_k)\mathbf m_k$ and $G(\theta)$ the Green
operator (P and S radiation, propagation, polarisation) evaluated at the hypothesis. The reformulation is to
solve

$$\hat\mu\in\arg\min_{\mu\in\mathcal M(\Theta;\mathbb R^6)}\ \tfrac12\|\mathbf d-\Phi\mu\|_{\Sigma^{-1}}^2+\lambda\,|\mu|_{\mathrm{TV},G},\qquad |\mu|_{\mathrm{TV},G}=\sum_k\big\|\Sigma^{-1/2}G(\theta_k)\,\mathbf m_k\big\|_2,$$

a least-squares problem with a group total-variation norm over vector measures: the
BLASSO of Bredies and Pikkarainen (2013), of De Castro and Gamboa (2012) and of Candès and Fernández-Granda
(2014), with amplitudes in $\mathbb R^6$. The norm of each atom is the energy that its tensor produces in
the whitened network; with it, the penalty favours neither locations close to the stations nor
particular mechanisms (for discrete measures it is written as a sum; in general it is the total variation
of the reparametrised measure). Three properties characterise it. The first two, in their form on a
grid, are already possessed by the block-sparse inversion of Vera Rodriguez et al. (2012); what is specific to this
formulation is the continuous measure with group norm and its dual certificate (Section 7).

**The unknown mechanism ceases to be a nuisance.** Since the group norm is invariant under
rotations of $\mathbf m$, the problem privileges no mechanism. Nor does it require searching for it: for
a fixed location, $\mathbf m$ is solved by least squares.

**The number of events comes out of the problem.** The total-variation norm is the convex relaxation of the
cardinality (Chandrasekaran et al., 2012), and $\lambda$ controls how many events are declared.

**The optimality condition is a physical detector.** If $\mathbf r=\mathbf d-\Phi\hat\mu$ is the
residual, $\hat\mu$ is optimal if and only if the dual certificate

$$\eta(\theta)=\frac{1}{\lambda}\big\|\,P_{\theta}\,\Sigma^{-1/2}\mathbf r\,\big\|_2,\qquad P_\theta=\text{orthogonal projector onto the range of }\Sigma^{-1/2}G(\theta),$$

satisfies $\eta\le1$ on all of $\Theta$ and $\eta=1$ at each declared event. This holds for the minimisers
(Appendix AA, Proposition AA.1); the greedy iteration of the pilot is not a minimiser and does not inherit the
equality at its events. With $\mu=0$,
$\lambda^2\eta^2$ is exactly the generalised likelihood-ratio test with the moment
tensor as a linear nuisance parameter. That test is the matched subspace detector of Scharf and
Friedlander (1994), applied jointly to all stations and to P and S. The statistic demands coherence of signs and
of relative amplitudes between stations, that is, the radiation pattern, and no phase detector
nor any scalar characteristic function preserves it.

## 3. Mathematical tools

**Frank-Wolfe as reference, greedy deflation as implementation.** The natural algorithm for the
BLASSO is the conditional gradient over measures in its sliding version (Boyd, Schiebinger and Recht,
2017; Denoyelle et al., 2020): at each iteration the linear oracle seeks the global maximum of $\eta$, that
hypothesis is added, the amplitudes are refitted and the positions are slid. Its $O(1/k)$ rate (Jaggi,
2013) assumes an exact oracle or one with bounded additive error. The pilot (MT-FW) is not that algorithm: it is
a group matching pursuit with a physical dictionary. It chooses the maximum of the certificate,
fits by least squares the tensor of that atom without refitting the previous ones, subtracts its contribution
at the selected peaks (with the per-station cap) and excludes, from the second atom on, the nodes at
less than 5 km epicentral distance from those already accepted. For that iteration there is no rate, and with the cap active the
residual may fail to decrease (Appendix AA, Proposition AA.4(e)-(h)); what Appendix AA does give is an
a posteriori suboptimality certificate for candidates written in the dictionary of the grid problem.

**The convex problem as framework, MT-FW as detector (Task 2C).** Task 2C solved the convex problem
restricted to the pilot grid (675 nodes at 3 km, origins every 0.1 s), with blocks of 6
columns read at the exact predicted times and the group norm whitened per block, by fully corrective
Frank-Wolfe with duality gap. In the 1200 windows of 2A the maximum relative gap was
1.4e-8, the certificate off the support did not exceed 1.000000000 (admission tolerance 1 + 1e-6) and on the
support it deviated from 1 by at most 8.7e-10: the solution is optimal with certificate. With that norm, the certificate of each block is the
selection statistic of the block orthogonal matching pursuit of Vera Rodriguez et al. (2012), and
λ_max²/2 is its first score. Even so, the convex problem falls below MT-FW in 23 of 24 cases of
cell, seed and cap, by −0.064 to −0.605 of recall, and at the level of the block pursuit or somewhat
below it. The cause, measured in exploration, is the grid: in the ideal generator with the sources on the nodes
the block pursuit beats MT-FW (0.97 against 0.76 at 0.25 false events per window with 4 stations) and the BLASSO reaches
0.73; off the nodes, the moveout error left by the grid, added to the residual of the station
corrections, makes the amplitudes be read outside the main lobe of the matched filter, and the tensor fit is
destroyed. MT-FW performs thanks to non-convex ingredients that the problem does not contain: the peak search of
±0.25 s per station and phase, which absorbs those errors, the cap and the greedy deflation. The division of roles
is thus: the convex problem is the framework (it defines the unknown, the certificate and the gap that measures how far
a solution is from the optimum) and MT-FW is the detector. Whether the convex problem can also be the detector
depends on removing the grid from it. The sliding Frank-Wolfe, which moves in the continuum the position and the origin
time of each atom, attacks the discretisation, which the diagnosis points to as the dominant term; the
per-station time errors that remain after the corrections would probably also call for a
time tolerance within the problem. It is the key question of the mathematical line, and it is
unimplemented.

**Tolerance to time error as a mathematical object (Task 2D).** Task 2D gave the comparators
the same peak search as MT-FW. Block orthogonal matching pursuit with reading at the peak (BOMP-PK)
matches MT-FW (−0.155 to +0.091 at 0.25, all hierarchical CIs contain 0); without its cap, MT-FW is
essentially that same pursuit, and the cap has no significant effect. The BLASSO on the grid with the same
peak search (GL-FW-PK), whose oracle reads the certificate at the residual peaks and fixes those times
when the atom enters, remains 0.11 to 0.30 below MT-FW in 11 of 12 cases, and its certificate is valid only
for the dictionary that this heuristic oracle builds. The ingredient that matters is the tolerance to
time error, and today it is a heuristic: an argmax per station and phase in a fixed window of ±0.25 s, without
an error model, without penalisation of the shifts and without a certificate on them. Its principled
form is to put the shifts into the problem: each atom carries a time shift
$\delta_{i,b}$ per station $i$ and phase $b$, bounded by a declared tolerance or penalised, the
certificate becomes the maximum of $\eta$ over position, origin time and shifts, and the
sliding Frank-Wolfe moves all three in the continuum, with an unbiased refit on the chosen support.
Appendix AA.8 already observes that within the search window the implemented statistic is a lower
bound of a likelihood ratio maximised over shifts; making that maximised problem
the object of the optimisation, with its duality gap and its theory, is the first mathematical deliverable.
It counts as a methodological advance only if in the windows of 2A it beats at the same time the heuristic peak search
of MT-FW and BOMP-PK, which is the comparator to beat.

**What the test said (Task 2E).** That solver was implemented as is (TS-SFW+SC: continuous atom with one
shift per station and phase bounded to ±0.25 s, analytic gradients, convex step with certified
gap, λ = f · s₁ with f = 0.01 calibrated) and was preregistered with a hash delivered through the chat before
writing the code. Against BOMP-PK+SC in 30 networks per cell it wins +0.000 to +0.039 at 0.25, all hierarchical
CIs contain 0: **G1-M is refuted in this design**. Three measured reasons: from a grid node
the local refinement skips cycles (in the noise-free test it recovers the source to < 0.5 km and < 0.02 s
only in 12 of 90 cases); with the calibrated f the group penalty is inert; and the missing events are
absorbed by atoms mislocated by 5-10 km, which the local sliding (median 0.7 km) does not move. The
independent verification found no implementation error and no budget disadvantage.

**The robust certificate: what the cap does (Task 2F).** In standard windows the per-station cap does not
matter (2D, 2F-c); in windows with 3 to 5 overlapping events it is the whole difference between MT-FW and BOMP-PK. The
capped statistic replaces each peak vector $\mathbf v_{i,b}$ by
$\mathbf v_{i,b}\min(1,\kappa/\|\mathbf v_{i,b}\|)$ before projecting it onto the moment-tensor
subspace, with $\kappa = 6$ in units of noise standard deviation, fixed in Task 1x: a Huber-type
robust certificate, which bounds the influence of each station and phase. On new seeds (4242, 8888) and with
96 networks per cell it wins +0.043 to +0.074 at 0.25 over BOMP-PK, and without the cap MT-FW matches BOMP-PK. The
verification located the effect: it bounds the influence of any station-and-phase block whose amplitude the
radiation operator predicts poorly, calibrates the score better across windows under a global threshold, depends on the
level $\kappa$ and does not separate the events within a window better. The amplitude factors of the generator are
drawn in each window and in each catalogue event, so they are per-window amplitude perturbations and
not station terms; removing them reduces the gain of the cap, but with small dispersion the cap still
wins where a statistic that knows the true factors does not (exploration C8, Appendix Y.25). The
mechanism is, therefore, bounded per-block influence in the face of the model's amplitude error, not correction of
a site term.

**Site terms versus bounded influence (exploration C8, exploratory, seed 999).** A gain
$g_{i,b}>0$ per station and phase in the operator, $\Sigma^{-1/2}\mathrm{diag}(g)G(\theta)$, estimated with the
catalogue events like the time corrections, was tested in a world with persistent per-network factors
added to the generator. Neither the estimated gains (−0.040 [−0.078, +0.007] relative to the cap) nor the true ones
(−0.008 [−0.049, +0.055]) replace it, and with the cap in place knowing the factors adds +0.045 and +0.032. The
form that the exploration supports keeps the bounded influence and removes the calibrated level: a Huber M-estimator
of the moment tensor on the norm of the residual of each block,
$T=\sum_b\rho_c(\|\mathbf v_b\|)-\min_a\sum_b\rho_c(\|\mathbf v_b-Q_b a\|)$, with the threshold $c$ fixed by rule
from the noise (square root of the 0.95 quantile of chi-square with 3 degrees of freedom, 2.7955). In the pilot, as
the best of 11 candidates, it matches cap 6 with dispersion up to 0.5, is worse with 0.8 and costs 9 times more.

**The certificate with level fixed by the noise (Task 2I, preregistered; Appendix Y.26).** The Huber
function is applied to the norm of the residual of each 3-component vector per station and phase, not to its
components, so that a block is weighted as a whole and the statistic does not depend on the orientation of the
sensor axes: $\rho_c(r)=r^2/2$ if $r\le c$ and $cr-c^2/2$ if $r>c$, with $c = 2.7955$ noise deviations
of the matched filter (a pure-noise block stays in the quadratic branch with probability 0.95), weights
$w_b=\min(1,c/\|\mathbf r_b\|)$, 4 IRLS passes from the least-squares fit and deflation with the
robust fit. On seeds 5151 and 6262, never used, with 120 networks per cell and regime and 5 dense
and 5 standard windows per network, the three registered rules are met: it beats BOMP-PK in dense windows by
+0.058 [+0.039, +0.076] and +0.053 [+0.037, +0.070], it is not inferior to cap 6 (−0.008 and −0.002, bounds
lower than −0.03) and is equivalent to BOMP-PK in standard windows (+0.010 and +0.007, 90 CI within ±0.05).
The independent verification confirms it and fixes the scope: +0.031 and +0.019 with origin time within
0.5 s; equivalence not met with weak clutter in 6262 (−0.032); the recall with all atoms accepted
does not rise, so the gain is one of ordering and choice of atoms; somewhat below the cap at 0.1 in 5151
(−0.018); about 9 times the cost of MT-FW. The secondary arms, descriptive and with 20 networks per cell and
regime, say that the level does not matter between 2.5 and 3.4 noise deviations (differences of +0.001 to
+0.006 in dense windows), that with 5 the gain is smaller, that with site deviation 0.25 neither of the
two forms wins reliably and that with 0.8 the Huber fit beats BOMP-PK (+0.100 and +0.072) but
falls below the cap (−0.055 and −0.069, with the interval below zero in 6262). The Huber form is not
a better detector than the cap: what changes is the origin of the level, a rule on the noise distribution
instead of a calibration, and when the model's amplitude error is large the calibrated cap is better.
Neither of the two forms
closes the margin of 0.08 to 0.24 below the realistic ceiling, which neither the velocity model nor the source pulse
explains and of which no treatment of amplitudes, the true factors included, captures more than a
quarter, so that about three quarters remain unexplained (explorations `explore_contrib` and
C8).

**The timing budget and the evidence without travel times (explorations C9 and C10; Appendix Y.27).**
Two questions remained about the statistic: whether the hard window of ±0.25 s should be a Gaussian prior
on the shifts, and whether channels that do not depend on travel time could locate what the times
do not separate. Both answers are negative in this benchmark, with their mechanism, and they are exploratory
(seed 999). The shift of a true arrival relative to the prediction of its node has three
terms: the model error after a perfect station term, 0.05 s in P and 0.07 s in S; the error of
the corrections estimated with 20 catalogue events, 0.12 to 0.24 s rms; and the position within the cell
of 3 km, 0.18 s in P and 0.30 s in S. The budget is 0.2 to 0.4 s per arrival, a prior with that honest
width is almost flat within ±0.25 s, and loading the shifts gave no gain in any detector;
widening the window to ±0.5 s loses recall. With the true statics the same prior, at the level of the
model error, gains +0.081 in dense windows and +0.146 in standard ones, where it reaches the realistic ceiling: the
form of the prior is correct and its inputs are not accurate enough. The evidence without times fails for lack of
information: with a P arrival of 1.5 to 1.8 noise deviations, polarisation, amplitudes and the
S/P ratio together bound the epicentre to 4.4 to 5.6 km against 0.4 to 0.7 km from the times (Bayesian
Cramér-Rao bound), and adding them to the times improves the bound by 4 to 5 %. What remains as a lever is the
accuracy of the per-station time corrections and, after that, the amplitude evidence of the weak
arrivals.

**Station terms without picks (exploration C11 and Task 2L, preregistered; Appendix Y.28).** The calibration
lever was tested in the correct order: first an exploration on seed 999 with the plan frozen
with a hash before any detection result (C11), then a preregistered task on never-used seeds
(2L). The estimator S2 is a joint determination of the hypocentres of the 20 catalogue events
and of one term per station and phase, without picks: it reads the energy peaks of the matched filter, relocates
each event by Gauss-Newton with a Huber loss and a Gaussian prior of 1 km and 0.1 s around its
catalogue hypocentre, rereads the weak arrivals near the prediction and fixes the gauge with the
median shift of the origin. Within the refined detector (MT-FW with cap 6 and T3), with the hit
rule that requires the origin time, S2 raises the recall at 0.25 relative to the registered corrections by
+0.057 and +0.062 in standard windows and +0.038 and +0.030 in dense ones, with no loss by regime nor with the
epicentre-only rule (hypothesis G1-C). The independent verification found no leakage and bounds the scope: the gain
is essentially that of the one-sided cell with 8 stations; with 4 stations there is no demonstrated gain; S2
worsens the terms in about 22 % of the networks; the effect is conditional on T3 and assumes catalogue
errors of known size; and it captures 0.28 to 0.61 of the gain of the true delays. The
self-calibration with the detector's own atoms added nothing measurable in C11.

**Three subsequent explorations (C12, C13 and C14, seed 999; Appendix Y.29).** With the registered detector,
not that of 2L. Restricting the statistic to the double-couple manifold has a null law that the
tube formula predicts (4.4 to 4.5 effective degrees of freedom instead of 6) and lowers the threshold by 10 % in energy
at equal false-event rate, but it does not change the recall (+0.005 [−0.021, +0.033] in standard windows with origin time),
because the hits lose with the restriction as much energy as the threshold drops. An extra pass with the
station moveout permuted gives a label-free operating point with the epicentre-only rule,
not with the origin-time rule. The remedies for self-interference after deflation halve the
false events induced by the events without moving the recall, and the diagnosis explains why: at the operating
point the largest class of false events under the origin-time rule are well-located detections of real events with
the origin shifted by a median of 0.8 s, almost always early and almost always unique. The origin time
of those detections is the next object of calibration.

**Certified oracle.** No projector increases the norm, so $\|P_\theta\mathbf v\|$ over a
cell is bounded by the three-component energy in the time window enlarged by
$h/V_{\min}+\delta\tau$. Appendix AA proves that the bound holds for the implemented statistic,
with its peak selection and its cap (Proposition AA.3), and Task 1z verified it without violations in
48 residual fields × 300 cells and in 1,004,400 hypotheses. With it the certified oracle gives records
identical to the exhaustive search in 72 of 72 windows. But it does not prune: it avoids 4 % of the exact
evaluations and, with its bookkeeping, it is 3 times slower than the dense loop, because after the first
strong event the maxima are at the noise level and the bound, with a median of 4.6 to 5.0 times the true
maximum of the cell, exceeds them almost everywhere. An envelope per subspace, which recomputes in
batch the statistic with a rigorous bound on the floating-point error and evaluates with the exact call only
the winning node, is 4.3 times faster with identical records; it is certified vectorisation, not pruning.
A pruning that pays off in the noise-dominated iterations needs cell bounds that use the subspace
and the directions of the peak vectors, and it remains open. Pruning with the kinematic bound of Appendix Z
(42 of 42 exact windows, 2.3 to 3 times fewer evaluations) is a different statistic.

**Threshold calibrated on noise windows.** Under whitened Gaussian noise, the idealised certificate
$\lambda^2\eta^2(\theta)$ at $\mu=0$ is chi-square with six degrees at each hypothesis, and the
expected Euler characteristic of the excursion set (Worsley, 1994; Taylor, 2006; Adler and
Taylor, 2007) approximates the probability that its maximum exceeds $u$ only if in addition a
symmetry hypothesis (G) holds, which moment-tensor subspaces do not satisfy in general (Appendix AA,
Proposition AA.5). Task 1z measured that this route does not fix the false-alarm rate of the implemented
statistic (E28). The statistic is not chi-square with 6 degrees even at a single hypothesis: taking the energy
peak of the window brings the mean of $2T$ from 6 to 12.6-13.5. The field decorrelates within one step of
the grid, and the EEC with roughness-based curvatures overestimates the false-alarm probability by 20 to 60
times. And strong clutter shifts the maximum by +8 to +9, beyond any Gaussian prediction. What
is calibrated is the empirical threshold: the split conformal quantile of the maxima of about 300
noise-only windows of the same network and the same regime (Angelopoulos and Bates, 2023), with a finite-sample
guarantee under exchangeability of windows (Proposition AA.5(c)), or a GEV fitted to those
maxima, as in Hirano et al. (2022). Both stay within the predictive band in 24 of 24
combinations on Gaussian noise. With strong clutter, every threshold calibrated on Gaussian noise
fails (the conformal one gives rates of 0.60 to 0.96 where 0.01 to 0.25 are requested), and recalibrating on clutter
windows of the same law restores it. The threshold of the method is therefore the empirical one, recalibrated every time
the noise or clutter regime changes; the EEC and Bonferroni remain as a prior conservative bound
on Gaussian noise. Thus each method is compared at the same realised false-alarm rate, the
output variable of this whole thesis.

**Observability certificate and exact recovery.** The theory of the BLASSO says when the support is
recovered exactly: the non-degenerate source condition of Duval and Peyré (2015) requires that the
minimal-norm certificate touch 1 only at the true events, with strictly negative curvature
there. For one source, Appendix AA (Proposition AA.6) proves, under idealised hypotheses, that this
curvature is the efficient Fisher information of the location with the moment tensor as nuisance,
divided by the energy of the source. It is bounded below by the leading term of the observability
certificate of Section 6.1 with weights equal to the energies that the real mechanism radiates to each
station and phase (a station on a nodal plane contributes nothing), and above by the information with
known mechanism. The global part of the condition, that no other hypothesis explain the noise-free
signal, is not a local statement. Exact recovery with several sources on the continuum remains a
conjecture (Conjecture AA.1). The first result of the thesis thus moves from diagnosis to bound: the
observability weighted by radiated energy controls the curvature of the certificate at a true
source.

## 4. Physics: a ladder of fidelity, chosen per station

The known risk of moment-tensor migration at high frequency is the amplitude: an imperfect velocity
model, site effects and scattering alter amplitudes and polarities. That is why the
method does not impose a single level of physics. It uses a ladder of operators $G_\ell(\theta)$, each one a
subset of the constraints of the previous one:

- **L0, kinematic.** Travel times only, with the Lipschitz bound of the eikonal (C-III).
- **L1, ray geometry.** P along the ray and S in the transverse plane. It holds for any
  source in an isotropic medium and ignores amplitudes. It is WAVE-PT, validated in Appendix Y.10.
- **L2, radiation coherence.** The 6-dimensional subspace of the moment tensor: signs and
  relative amplitudes between stations. It is WAVE-MT (Task 1x).
- **L3, full waveform.** Green's functions computed in the 3D medium.

Each station enters at the highest level that its measured fidelity supports. Fidelity is measured with the
same per-station reliability methodology $\hat r_s$ with chance correction (Appendix T), applied
to amplitudes and polarities of catalogued events. Model uncertainties enter as nuisance parameters
with their own structure:

- a delay per station, maximised in a window of ±0.25 s (already included) and compensated beforehand
  with catalogued events when it is larger than that window (Task 1z);
- the amplitude error per station and phase, treated by bounding the influence of each block (cap or Huber
  fit); an estimated or marginalised site factor does not replace that bound in exploration C8;
- a theory-error covariance estimated as the difference between the physics that generates and the physics that
  infers (Section 7.3).

This ladder is the operational version of the thesis's sufficiency ladder: the required physics
is a function of the effective information.

## 5. Artificial intelligence where physics falls short, never in its place

**Neural Green's functions for 3D media.** Travel times come from a neural eikonal
(EikoNet; Smith et al., 2021) and amplitudes and polarisations from a source-dependent neural operator
(MIFNO; Lehmann et al., 2025) or from rays on the neural eikonal. They are exactly what
$G_\ell(\theta)$ needs. Since the moment tensor enters linearly, an operator that gives the six
elementary responses suffices for every source.

**Learned noise model.** $\Sigma^{-1}$ and the per-station cap that today protects against noise
bursts (CAPV) are replaced by a robust likelihood learned on real noise and clutter from the
RSNC: a heavy-tailed whitener, or a normalising flow on three-component features. The
dual certificate is kept as an inner product in the whitened space.

**Signed characteristic functions.** MALMI migrates scalar probabilities and loses sign and
direction. Here the learning network delivers a clean three-component displacement that preserves
polarity (of the DeepDenoiser type; Zhu, Mousavi and Beroza, 2019). This combines the robustness of
learning with the coherence of physics.

**Amortised proposals with a guarantee.** A network, trained with the physical generator itself, proposes
candidate hypotheses and initialises the deflation. The certified oracle then guarantees that no
better hypothesis was left out. The network only accelerates; its error cannot miss an event that
the exact search would find.

## 6. Pilot evidence

Tasks 1x and 1y (Appendices Y.11 and Y.12) implemented the first iteration of the method in the
consistent synthetic generator:

- the L2-level dual certificate (WAVE-MT);
- the greedy deflation (group matching pursuit), up to five events per window, with
  a no-reselection radius of 5 km (MT-FW);
- as a control, the same iteration with the subspace of another hypothesis more than 9 km away (MT-FW-NODE).

The generator uses random double-couple sources and phase detection according to the signal-to-noise ratio. It adds
per-station delays N(0, 0.15 s), about 10 degrees of polarisation deviation and lognormal site factors
with sigma 0.5, which no method knows. The clutter follows the same detection rule as
the signal in two regimes: weak, which is almost never detected, and strong, with visible transients that
the detector does mark. 100 new windows per cell and two seeds were used; GaMMA runs with a fixed
seed and oracle configuration. The rules were registered before each run.

**Table P1. Weak clutter. Clutter and true picks per window; recall at realised false events
<=0.1 / <=0.25.**
| Cell | picks clutter/true | WAVE-PT@20 | WAVE-MT@20 | MT-FW | MT-FW-NODE | GAMMA |
|:--|:--|:--|:--|:--|:--|:--|
| one-sided, 8 st. (2024) | 4.9/5.0 | 0.56/0.69 | 0.65/0.68 | 0.73/0.77 | 0.51/0.61 | 0.63/0.63 |
| 4 stations (2024) | 2.2/2.6 | 0.28/0.34 | 0.36/0.45 | 0.38/0.54 | 0.14/0.19 | 0.28/0.28 |
| 6 stations (2024) | 3.7/3.6 | 0.38/0.50 | 0.49/0.57 | 0.49/0.65 | 0.18/0.24 | 0.54/0.54 |
| one-sided, 8 st. (7331) | 4.8/4.7 | 0.61/0.66 | 0.61/0.65 | 0.66/0.71 | 0.44/0.57 | 0.55/0.55 |
| 4 stations (7331) | 2.2/2.7 | 0.31/0.47 | 0.48/0.51 | 0.49/0.63 | 0.20/0.41 | 0.35/0.35 |
| 6 stations (7331) | 3.5/4.1 | 0.51/0.58 | 0.60/0.64 | 0.67/0.76 | 0.24/0.46 | 0.54/0.54 |

**Table P2. Strong clutter.**

| Cell | picks clutter/true | WAVE-PT@20 | WAVE-MT@20 | MT-FW | MT-FW-NODE | GAMMA |
|:--|:--|:--|:--|:--|:--|:--|
| one-sided, 8 st. (2024) | 34.2/4.9 | 0.44/0.52 | 0.63/0.65 | 0.61/0.69 | 0.27/0.44 | 0.33/0.33 |
| 4 stations (2024) | 28.8/2.4 | 0.01/0.15 | 0.15/0.26 | 0.15/0.24 | 0.04/0.07 | 0.05/0.06 |
| 6 stations (2024) | 29.2/3.5 | 0.35/0.38 | 0.51/0.59 | 0.47/0.55 | 0.15/0.20 | 0.19/0.19 |
| one-sided, 8 st. (7331) | 31.1/4.8 | 0.40/0.50 | 0.51/0.55 | 0.56/0.61 | 0.30/0.43 | 0.24/0.24 |
| 4 stations (7331) | 29.8/2.7 | 0.06/0.20 | 0.28/0.36 | 0.30/0.37 | 0.07/0.15 | 0.04/0.05 |
| 6 stations (7331) | 29.7/4.1 | 0.34/0.41 | 0.41/0.55 | 0.43/0.61 | 0.28/0.35 | 0.12/0.18 |

**Reading**

1. **Against GaMMA.** MT-FW beats tuned GaMMA by +0.11 to +0.28 of recall at 0.25 false events per
   window with weak clutter and by +0.17 to +0.43 with strong clutter. This holds in at least 2 of 3
   cells per seed at both caps.
2. **Against its control.** The control with the geometry of another hypothesis loses between 0.13 and 0.41.
3. **With weak clutter**, migration with ray geometry only (WAVE-PT, level L1) does not beat
   GaMMA. The gain comes from two sources: the coherence of the moment-tensor subspace and the
   deflation. Deflation removes the sidelobes and separates the superposed events with a no-reselection
   radius of 5 km, instead of the 20 km suppression of the single-pass scores.
4. **With strong clutter**, deflation adds little over a 20 km suppression. Coherence
   rejects transients that do not fit any moment tensor; GaMMA, by contrast, has to
   accept them as picks.

**What the independent verification found** (three verifications in total)

- **Inflated margins before 1y.** In the tasks prior to 1y, GaMMA received clutter picks
  that no phase detector would produce (E26), and its margins were inflated. For example, in the
  one-sided cell, MT-FW beat it by +0.53 with that clutter and beats it by +0.15 with consistent
  clutter, in the same perturbed condition.
- **The wavelet does not break the method.** A wavelet rotated 90 degrees in phase, or of another central
  frequency, does not break MT-FW.
- **Large delays do break it.** Unmodelled delays with standard deviation 0.4 s, larger
  than the search window of ±0.25 s, make it fall below GaMMA. The per-station time corrections,
  calibrated with catalogued events, are therefore part of the method.
- **GaMMA's picks were optimistic.** The parametric pick model of 1o to 1y detects 42.0 % of
  the arrivals with a median error of 0.033 s, whereas an STA/LTA calibrated on the waveform detects
  14.3 % with 0.078 s. With weak clutter GaMMA received better picks than those of a real picker; with strong
  clutter it received about 28 clutter picks per window, against 7 to 8 from the STA/LTA (E27).

**What Tasks 1p, 1v and 1z tested** (exploratory, seed 999; reports in `results/components/`)

Against QuakeMigrate 1.2.2, run through its public API on 60 weak-clutter test windows,
MT-FW recovers 0.72 and 0.81 of the events at 0.1 and 0.25 false events per window, compared with 0.53 and 0.53 for its
fair three-component configuration (QM-3C) and 0.15 and 0.28 for its default practice, which fails in
this generator because the rays arrive at 60-80 degrees from the vertical. The differences with QM-3C are
+0.19 [−0.00, +0.37] and +0.28 [+0.10, +0.43], and they lie in the windows with two sources (22 of 28 sources
from pairs for MT-FW, 10 of 28 for QM-3C), because QuakeMigrate keeps one maximum per time sample.
It is a single seed with 47 sources, and the QM-3C list ends at 0.05 false events per window.

With delays N(0, 0.4 s) fixed per network, the corrections estimated from 20 catalogued events leave an
rms residual of 0.088 s with exact locations and of 0.173 s with catalogue errors of 1 km and 0.1 s, and
restore to MT-FW the recall of the exact correction in the 4-station and one-sided 8-station cells
(+0.18/+0.24 and +0.30/+0.30 over the uncorrected version, within +0.01 to +0.05 of the exact one). In the
6-station cell not even the exact correction gives a significant gain. The medium is homogeneous and the
wavelet exact: on real data the residual will be larger.

The threshold from random field theory is not calibrated and the certified oracle is exact but does not
prune; both results are in Section 3.

**The confirmatory test (Task 2A; preregistered, executed and independently verified)**

The design and the rules were written before evaluating. The generator does not match the model of the
detectors: a gradient medium with curved rays, a Ricker wavelet of variable central frequency and phase,
scattered coda, polarisation deviation, site factors and delays N(0, 0.3 s) fixed per network; all
methods use the same homogeneous model. The associators receive the picks of an STA/LTA frozen on
seed 999, and the station corrections (SC) are estimated from 20 catalogued events with errors of 1 km and
0.1 s and are delivered to every method that admits them. There are 5 networks and 100 windows per cell, two seeds and two
clutter regimes.

**Table P3. MT-FW+SC minus the best of GaMMA, PyOcto and QuakeMigrate (seven variants), recall at realised false
events <=0.1 / <=0.25.**

| Cell | weak (2024) | weak (7331) | strong (2024) | strong (7331) |
|:--|:-:|:-:|:-:|:-:|
| 4 stations | +0.14 / +0.26 | +0.17 / +0.35 | +0.19 / +0.25 | +0.09 / +0.17 |
| 6 stations | +0.18 / +0.25 | +0.49 / +0.46 | +0.26 / +0.32 | +0.18 / +0.17 |
| one-sided, 8 st. | +0.20 / +0.27 | +0.17 / +0.28 | +0.16 / +0.23 | +0.37 / +0.41 |

The primary rule 2A-a is met in both regimes, in the 3 cells of both seeds, and so is the control with
the subspace of another node (2A-b): the gain requires the geometry of the true hypothesis.
Against GaMMA with optimistic picks (2A-c) it is met, with a tie in the one-sided cell with weak clutter
(+0.01 in both seeds). The station corrections help with weak clutter and not with strong clutter
(2A-d), because the catalogue drawn with that clutter leaves a residual of 0.20 to 0.28 s. With thresholds
frozen on other networks (2A-e) the recall advantage holds in all cells, but in the
one-sided 8-station cell, seed 7331, budget 0.1, the realised rate is 0.22 and 0.23, above the limit of
0.20: the rule is met under the per-cell reading and not under the strict one.

The verification reproduced the figures, found no leakage of the truth nor asymmetry in the hit rule,
checked the sign of the corrections against the true delays and showed that the near-zero recall of
GaMMA and PyOcto is starvation by the picker: with true-time picks and the same wiring they recover 0.88
and 0.69 at 0.1. That is why the effective rival is QuakeMigrate, the best comparator in all 24 cases. The reading
has four limits. The window-level intervals ignore the clustering by network, and with a network-level bootstrap
four cells include zero at 0.1. About half of the events are in windows with two sources,
where the QuakeMigrate trigger merges nearby maxima; without them the margin shrinks but persists. In
single-source windows GaMMA with optimistic picks beats MT-FW+SC at 0.1 in 5 of 6 weak-clutter
cells. And the sources are synthetic point double couples, inside the subspace of the method. 2A also does not
isolate the radiation pattern: the control of 2A-b changes times and radiation at the same time.

**The ablation of the radiation operator (Task 2B; preregistered and executed)**

On the same 1200 windows of 2A, regenerated bit for bit, MT-FW-FREE+SC is the same solver (same
greedy search, cap, 5 km exclusion, corrections and deflation) with the radiation operator
replaced by the identity: the statistic is the total energy of the peak vectors, invariant to
sign, and the deflation subtracts the observed vectors. The preregistered rule requires MT-FW+SC − MT-FW-FREE+SC
>= 0.05 at both caps in at least 2 of 3 cells, in both seeds.

**Table P4. MT-FW+SC minus MT-FW-FREE+SC, recall at realised false events <=0.1 / <=0.25.**

| Cell | weak (2024) | weak (7331) | strong (2024) | strong (7331) |
|:--|:-:|:-:|:-:|:-:|
| 4 stations | +0.17 / +0.28 | +0.15 / +0.24 | +0.25 / +0.26 | +0.06 / +0.12 |
| 6 stations | +0.18 / +0.11 | +0.26 / +0.12 | +0.32 / +0.35 | +0.26 / +0.17 |
| one-sided, 8 st. | −0.01 / −0.01 | +0.04 / +0.0494 | +0.17 / +0.22 | +0.33 / +0.33 |

The rule is met in both regimes. With strong clutter it is met in the 3 cells of both seeds and it
holds in 97-98 % of the bootstrap replicates; without the operator, the detector stays at the level of
QuakeMigrate (between −0.08 and +0.08), so there the advantage of 2A is due to the operator. With weak clutter it is
met at the minimum, 2 of 3 cells per seed: in the one-sided 8-station cell both methods tie, and the rule
holds only in ~60 % of the replicates. In exploration, without the operator the detector still beats the best
comparator of 2A at 0.25 in 5 of 6 weak-clutter cell-seeds, by +0.10 to +0.34: a good part of that
advantage comes from pick-free migration with deflation and corrections. The identity removes at once the
sign of the pattern, the polarisation constraint and the amplitude geometry, so 2B shows that
the operator matters and does not isolate the sign. The first attempt stopped at the exact check of the
records of 2A because 18 networks had been computed with another OpenBLAS kernel; with that kernel fixed they
are reproduced bit for bit and no recall changes.

**The closest comparators, the convex problem and replication by network (Task 2C; preregistered and
executed)**

It was preregistered before implementing it, in response to the independent audit (56/100). On the same
1200 windows of 2A, regenerated bit for bit, three methods are run that share with MT-FW the dictionary of
moment-tensor blocks, the grid of 3 km and 0.1 s and the station corrections, and that read each
station at the exact predicted sample, without peak search or cap: block orthogonal matching pursuit in
the form of Vera Rodriguez et al. (2012) (BOMP+SC), the mechanism-corrected stack in the form of Anikiev
et al. (2014) (MTSTACK+SC), both in our own implementation, and the group BLASSO on the grid by
fully corrective Frank-Wolfe (GL-FW+SC). In addition, 25 new networks per cell, seed and regime bring the
comparison with QuakeMigrate to 30 networks per cell, with a hierarchical bootstrap (networks and then windows).

**Table P5. 2C-a: MT-FW+SC minus the best of BOMP+SC and MTSTACK+SC, recall at realised false events <=0.25.**

| Cell | weak (2024) | weak (7331) | strong (2024) | strong (7331) |
|:--|:-:|:-:|:-:|:-:|
| 4 stations | +0.287 | +0.231 | +0.130 | +0.077 |
| 6 stations | +0.333 | +0.441 | +0.324 | +0.202 |
| one-sided, 8 st. | +0.268 | +0.519 | +0.205 | +0.453 |

**Table P6. 2C-b: GL-FW+SC minus MT-FW+SC, recall at realised false events <=0.1 / <=0.25.**

| Cell | weak (2024) | weak (7331) | strong (2024) | strong (7331) |
|:--|:-:|:-:|:-:|:-:|
| 4 stations | −0.264 / −0.310 | −0.244 / −0.295 | −0.273 / −0.208 | −0.013 / −0.064 |
| 6 stations | −0.375 / −0.431 | −0.603 / −0.559 | −0.284 / −0.351 | −0.333 / −0.310 |
| one-sided, 8 st. | −0.479 / −0.451 | −0.556 / −0.605 | −0.229 / −0.265 | −0.547 / −0.573 |

**Table P7. 2C-c: MT-FW+SC minus QM-3C+SC with 30 networks per cell, recall at realised false events <=0.25, with
95 CI from the hierarchical bootstrap.**

| Cell | weak (2024) | weak (7331) | strong (2024) | strong (7331) |
|:--|:-:|:-:|:-:|:-:|
| 4 stations | +0.347 [+0.248, +0.453] | +0.306 [+0.194, +0.401] | +0.175 [+0.104, +0.254] | +0.099 [+0.029, +0.181] |
| 6 stations | +0.323 [+0.228, +0.416] | +0.347 [+0.242, +0.439] | +0.247 [+0.150, +0.337] | +0.231 [+0.125, +0.339] |
| one-sided, 8 st. | +0.283 [+0.189, +0.376] | +0.224 [+0.119, +0.328] | +0.236 [+0.135, +0.334] | +0.295 [+0.212, +0.383] |

Rule 2C-a is met in both regimes, in the 3 cells of both seeds; the better of the two
comparators is BOMP+SC in all of them. According to the registered interpretation, the advantage of MT-FW is not due to the moment
tensor with deflation, which is already in the block pursuit, but to its own elements. 2C-b, which has no
threshold, comes out below in 23 of 24 cases: the convex problem solved with certified optimality
performs much worse than its greedy approximation (Section 3). Rule 2C-c is met in both regimes, with and
without Holm; with Holm, the 4-station cell with strong clutter of seed 7331 does not reject (p 0.044), and
with strong clutter and 4 stations MT-FW+SC does not significantly beat BOMP+SC or MTSTACK+SC in
either seed. Exact checks: records of MT-FW+SC recomputed equal to those of 2A in 120 of 120
windows and certificate equal to dense in 600 of 600; the new methods were validated on seed 999
against explicit atoms, noise-free sources and the KKT condition. The exploratory diagnosis, posterior to
part (i), explains the collapse of the exact-time methods: in the ideal generator with the sources on
the nodes BOMP beats MT-FW, so what counts is the tolerance of the peak search to the coarse grid and to the residual time
errors. What 2C supports is robustness to grid and time errors on a coarse
grid, not superiority over sparse moment-tensor inversion as such: a block
inversion on a finer grid or with time tolerance, as would be used on real data, has not been tested.
Task 2D tested the variant with time tolerance: the attribution to MT-FW's own elements does not
hold up (E35).

**Equal tolerance, ablation of the peak search and frozen threshold (Task 2D; preregistered with clock
time and hash, and executed)**

It was preregistered on 26 September 2026 at 21:13 UTC, before writing the code, in response to the
second independent audit. Part (i): on the same 1200 windows of 2A, the three methods of 2C
read at the energy peak within ±0.25 s of the predicted time, like MT-FW (BOMP-PK+SC, MTSTACK-PK+SC,
GL-FW-PK+SC), and two ablations of MT-FW: without peak search (MT-FW-EX+SC) and without cap (MT-FW-NOCAP+SC).
Part (ii): split-conformal thresholds frozen on 200 empty windows of seed 999 per cell, at
α = 0.1 and 0.25, for MT-FW+SC and QM-3C+SC.

**Table P8. 2D-a: MT-FW+SC minus BOMP-PK+SC (the best tolerant comparator in all cells), recall at
realised false events <=0.25, with 95 CI from the hierarchical bootstrap.**

| Cell | weak (2024) | weak (7331) | strong (2024) | strong (7331) |
|:--|:-:|:-:|:-:|:-:|
| 4 stations | +0.080 [−0.136, +0.282] | +0.090 [−0.029, +0.256] | +0.091 [−0.048, +0.257] | −0.051 [−0.206, +0.095] |
| 6 stations | +0.042 [−0.092, +0.224] | −0.029 [−0.125, +0.047] | +0.000 [−0.138, +0.145] | −0.155 [−0.277, +0.032] |
| one-sided, 8 st. | +0.056 [+0.000, +0.265] | +0.025 [−0.084, +0.155] | −0.036 [−0.130, +0.108] | −0.053 [−0.158, +0.030] |

**Table P9. 2D-b: MT-FW+SC minus MT-FW-EX+SC (without peak search), recall at realised false events <=0.25.**

| Cell | weak (2024) | weak (7331) | strong (2024) | strong (7331) |
|:--|:-:|:-:|:-:|:-:|
| 4 stations | +0.230 | +0.218 | +0.156 | +0.064 |
| 6 stations | +0.153 | +0.412 | +0.149 | +0.167 |
| one-sided, 8 st. | +0.310 | +0.407 | +0.241 | +0.293 |

2D-a is not met in any regime: with the same tolerance, the block pursuit matches MT-FW, and according to
what was preregistered its advantage over Vera Rodriguez et al. (2012) and Anikiev et al. (2014) is attributed only to the
peak search. 2D-b is met in both: removing the peak search costs 0.06 to 0.41 in 12 of 12 cases
(lower bound of the CI above 0 in 7 of 12). Descriptive: without cap, MT-FW changes by −0.091 to +0.131 with all
CIs containing 0; GL-FW-PK+SC stays 0.11 to 0.30 below MT-FW+SC in 11 of 12 cases; BOMP-PK+SC
beats QuakeMigrate by +0.17 to +0.49. 2D-c is met per cell in both regimes and not under the strict
reading with weak clutter (0.222 > 0.2 in one cell), but the threshold calibrated on empty windows leaves 1.2
to 2.2 false events per window with events in weak clutter, 3 to 8 times the false events of QuakeMigrate over all
windows: it is not a deployable operating point. With 5 networks per cell no rule survives
Holm. Exploratory: with an injected time error of ±0.2 s, the variants with reading at the peak
beat their exact versions in 12 of 12 cases (+0.04 to +0.38); without error and with the sources on the
nodes, the peak search costs 0.10 to 0.14 with 6 and 8 stations. Exact checks (1200 of 1200
windows, 120 of 120 records of MT-FW+SC, 60 of 60 in the calibration) and validation on seed 999 (V1
to V5).

**The off-grid BLASSO with shifts (Task 2E; preregistered with a hash delivered through the chat
before writing the code, executed and verified)**

It was preregistered on 27 September 2026 at 14:44:14 UTC (hash b40f0c6e…e988); the calibration note
(f = 0.01) was recorded at 15:38:20 UTC, before running seeds 2024 and 7331. TS-SFW+SC against BOMP-PK+SC and
against its control with λ = 0 (BOMP-PK-SLIDE+SC), in 30 networks per cell (4800 windows).

**Table P10. 2E-a: TS-SFW+SC minus BOMP-PK+SC, recall at realised false events <=0.25, with 95 CI from the hierarchical
bootstrap.**

| Cell | weak (2024) | weak (7331) | strong (2024) | strong (7331) |
|:--|:-:|:-:|:-:|:-:|
| 4 stations | +0.039 [−0.025, +0.122] | +0.016 [−0.038, +0.094] | +0.000 [−0.033, +0.030] | +0.003 [−0.027, +0.039] |
| 6 stations | +0.029 [−0.025, +0.087] | +0.029 [−0.024, +0.078] | +0.033 [−0.013, +0.083] | +0.015 [−0.024, +0.058] |
| one-sided, 8 st. | +0.025 [−0.042, +0.094] | +0.015 [−0.035, +0.091] | +0.016 [−0.026, +0.059] | +0.018 [−0.044, +0.080] |

2E-a is not met in any regime (0 of 3 cells per seed; with Holm, minimum adjusted p 0.81 and 0.58) and
neither is 2E-b (−0.020 to +0.013): **G1-M is refuted in this synthetic design** and the group penalty is
inert with the calibrated f. The independent verification found no implementation error (the continuous atom
reproduces the BOMP-PK block to 1.4e-15 at the prediction) nor a budget disadvantage (−0.014 with 10 times
the iterations and 9 restarts), and measured a large margin: realistic oracle ceilings 0.08 to 0.24 above BOMP-PK,
with the missed events absorbed by atoms at 5-10 km.

**Exploration on seed 999 (`explore_contrib`, exploratory)**

Coarse grids (C1), large time errors (C3), an oracle or calibrated 1D velocity model inside
the detectors (C6) and a subspace invariant to the source pulse (C7) do not win or they lose; the
attribution run C6b shows that the margin persists with the exact medium (+0.194 with weak clutter). The only
candidate that meets the exploratory criterion is the per-station cap in windows with 3 to 5 overlapping events
(+0.070 at 0.25, +0.068 on new networks); it was chosen among several, which is declared.

**The per-station cap in dense sequences (Task 2F; preregistered with a hash delivered before running,
new seeds, frozen code, executed and verified)**

It was preregistered on 28 September 2026 at 03:34:51 UTC (hash 3d8e4c9b…2f5f), before running anything
on seeds 4242 and 8888; 96 networks per cell, regime and seed.

**Table P11. 2F-a: MT-FW+SC minus BOMP-PK+SC in dense windows, pooling cells, recall at realised false
events <=0.25, with 95 CI from the hierarchical bootstrap.**

| Regime | 4242 | 8888 |
|:--|:-:|:-:|
| weak | +0.072 [+0.048, +0.092] | +0.074 [+0.052, +0.100] |
| strong | +0.051 [+0.032, +0.070] | +0.043 [+0.022, +0.066] |

2F-a is met in both regimes and seeds, also with Holm; 2F-b is met (MT-FW − MT-FW-NOCAP +0.060 and
+0.059; MT-FW-NOCAP − BOMP-PK −0.001 and −0.002); 2F-c is met (equivalence with margin 0.05 in standard
windows, narrow in 8888). **G1-R is confirmed in synthetic data.** The adversarial verification confirms and
narrows: +0.024 to +0.045 with origin time within 0.5 s; effect tied to the level of the cap (peak at 6);
main mechanism, bounded influence of each station-and-phase block in the face of the model's amplitude
error (without the generator's amplitude factors, which are drawn per window, the gain drops to +0.023 with
a CI that includes 0) and better calibration of the score across windows, not separation of events. Absolute
recall with cap: 0.22 to 0.29.

**Exploration C8 (seed 999, exploratory; Appendix Y.25).** In dense windows with dispersion 0.5, the
catalogue amplitude terms without cap fall −0.066 [−0.108, −0.028] below MT-FW in the unchanged
generator and −0.040 [−0.078, +0.007] with persistent factors; with the cap added they contribute nothing (−0.002 and
+0.006). The true factors without cap match MT-FW (−0.011 and −0.008) and with cap add +0.045 [+0.009,
+0.074] and +0.032 [+0.003, +0.073]. With dispersion 0.25 the oracle gains nothing over BOMP-PK and the cap gains
+0.043 [−0.026, +0.099]. HUBER-N beats BOMP-PK by +0.064 [+0.028, +0.105] and +0.066 [+0.024, +0.105], does not
differ from cap 6 (+0.022 and +0.005, CI containing 0), is worse with dispersion 0.8 (−0.043 [−0.113, −0.001]) and
costs 9 times more; it is the best of 11 candidates and its estimate carries selection bias. MT-FW captures 0.14
and 0.15 of the margin up to the ceiling, HUBER-N 0.21 and 0.16, the oracle with cap 0.29 and 0.23.

**Task 2I (preregistered, seeds 5151 and 6262; executed and verified, secondary arms analysed;
Appendix Y.26).**
Registration hash a983beac…6f97 delivered to the user before any run; code frozen by hashes
at 18:14:44 UTC on 2 October 2026; the run was interrupted twice and resumed from the parts
with the code unchanged. 120 networks per cell and regime, 5 dense and 5 standard windows per network: 3600
dense windows per seed (14477 and 14407 events) and 3600 standard ones (2938 and 2891 events). Recall at 0.25 in
dense windows: BOMP-PK+SC 0.177 and 0.198, MT-FW+SC with cap 0.243 and 0.253, MT-FW without cap 0.175 and 0.195,
HUBER-N+SC 0.235 and 0.251.

| Rule | 5151 | 6262 | Met |
|:--|:--|:--|:--|
| 2I-a, HUBER-N − BOMP-PK, dense | +0.058 [+0.039, +0.076] | +0.053 [+0.037, +0.070] | yes |
| 2I-b, HUBER-N − MT-FW (cap 6), dense, margin −0.03 | −0.008 [−0.023, +0.006] | −0.002 [−0.018, +0.012] | yes |
| 2I-c, HUBER-N − BOMP-PK, standard, 90 CI | +0.0099 [−0.0100, +0.0272] | +0.0073 [−0.0132, +0.0235] | yes |

In the same windows the cap repeats 2F: MT-FW − BOMP-PK +0.066 [+0.046, +0.086] and +0.055 [+0.038, +0.074].
The independent verification (provenance, symmetry, rederived statistics, six regenerated parts
identical) confirms the positive result and finds: with the hit rule that requires origin time within 0.5 s,
+0.031 and +0.019; equivalence in standard windows not met with weak clutter in 6262 (−0.032 [−0.064, −0.003]);
recall with all atoms accepted 0.623 (BOMP-PK) and 0.622 (HUBER-N) in 5151, 0.628 and 0.616 in 6262, that
is, the gain is one of ordering and choice of atoms, half reliability of the score across windows and
half choice of atoms; −0.018 [−0.034, −0.002] relative to the cap at 0.1 in 5151; 0.927 s per dense window
against 0.103 s for MT-FW, about 9 times.

Secondary arms (descriptive, 20 networks per cell and regime, 600 windows per design and seed, without
independent verification; `analyze_2i_out.txt`), differences at cap 0.25 in dense windows:

| Contrast | 5151 | 6262 |
|:--|:--|:--|
| Threshold 2.5003 − registered threshold | +0.003 [−0.012, +0.021] | +0.001 [−0.008, +0.016] |
| Threshold 3.3682 − registered threshold | +0.005 [−0.015, +0.027] | +0.006 [−0.012, +0.020] |
| Threshold 5 − registered threshold | −0.013 [−0.039, +0.017] | −0.023 [−0.049, +0.009] |
| sd 0.25: HUBER-N − BOMP-PK | +0.036 [+0.004, +0.070] | +0.002 [−0.035, +0.044] |
| sd 0.25: MT-FW − BOMP-PK | +0.005 [−0.038, +0.041] | +0.009 [−0.024, +0.053] |
| sd 0.8: HUBER-N − BOMP-PK | +0.100 [+0.054, +0.143] | +0.072 [+0.031, +0.115] |
| sd 0.8: HUBER-N − MT-FW | −0.055 [−0.113, +0.002] | −0.069 [−0.121, −0.019] |
| sd 0.8: MT-FW − BOMP-PK | +0.155 [+0.096, +0.207] | +0.141 [+0.095, +0.189] |

The level does not matter between 2.5 and 3.4 noise deviations; with site deviation 0.25 the cap contributes nothing and
the Huber fit little; with 0.8 the Huber fit beats BOMP-PK and falls below the cap, so that
the non-inferiority of 2I-b does not extend to that dispersion.

**Explorations C9 and C10 (seed 999, networks 130-149 and 160-179; exploratory; Appendix Y.27).** 600 windows
per design, 95 CI by hierarchical bootstrap, without multiplicity correction; choices made after the
results are declared in the reports (`results/explorations/REPORT_C9.txt`, `REPORT_C10.txt`).

| Result | Value |
|:--|:--|
| C9, Cramér-Rao bound of the epicentre, times at 0.144 s (4, 6 and 8 stations) | 0.68, 0.43 and 0.69 km |
| C9, the same bound, all channels without travel times | 5.59, 4.35 and 4.81 km |
| C9, ratio of the joint bound to that of the times | 0.960, 0.962 and 0.953 |
| C9, travel-time-free locators on 600 isolated events, median error | 7.7-8.4 km the best (moment-tensor projection); 0.5-0.7 km with the times |
| C9, relocation by the moment-tensor projection, MT-FW, dense / standard | −0.116 [−0.153, −0.083] / −0.255 [−0.336, −0.183] |
| C9, the same with the true amplitude factors, dense | −0.116 [−0.158, −0.081] |
| C9, SELPOL − MT-FW and − BOMP-PK, dense (chosen after seeing results, best of about 50) | +0.046 [+0.007, +0.070] and +0.045 [+0.003, +0.086] |
| C9, reordering by polarisation admissibility at 0, 10, 20 and 45 degrees | +0.031, +0.028, +0.003 [−0.011, +0.023] and −0.015 [−0.025, +0.002] |
| C10, MT-FW with time prior (T2w), catalogue level, dense / standard | −0.002 [−0.046, +0.033] / −0.034 [−0.111, +0.019] |
| C10, hard window of ±0.5 s, dense | −0.047 [−0.087, −0.008] |
| C10, true statics with hard window, dense / standard | +0.047 [−0.003, +0.082] / +0.098 [+0.032, +0.141] |
| C10, true statics with prior at the level of the model error, dense / standard | +0.081 [+0.035, +0.117] / +0.146 [+0.065, +0.217] |
| C10, recall of the oracle and of the CEIL-GD ceiling in standard windows | 0.592 and 0.603 |
| C10, fraction of the margin captured by the oracle in dense windows | 0.31 [0.21, 0.43] |
| C10, isolated events that win at the first iteration: hard rule, catalogue prior, oracle | 0.485, 0.489 and 0.609-0.634 |
| C10, T3 refinement in standard windows, registered rule / with origin time (BOMP-PK) | +0.059 [−0.009, +0.112] / +0.100 [+0.031, +0.168] |
| C10, SC-R corrections with T3 in standard windows, registered rule / with origin time (after the results) | +0.064 [−0.004, +0.123] / +0.123 [+0.058, +0.201] |

Reading: the channels without travel times are limited by information and not by algorithm; soft
time priors do not help as long as the station corrections are of catalogue quality; with the
statics known at the level of the model error the detector reaches the ceiling on isolated events; and with
perfect timing 37-39 % of the isolated events still lose because of amplitude. The margin splits into
a calibration part and an amplitude part.

**Task 2L (preregistered, seeds 7171 and 9393; executed and verified; Appendix Y.28).** Registration hash
13a6e9b1…a62a delivered to the user before any run; code frozen by hashes of 37 files at
21:24:48 UTC on 4 October 2026; the run was interrupted when the session was suspended and resumed
from the parts with the code unchanged. 160 networks per cell and regime with 5 standard windows and 100 with
5 dense ones per seed. Contrast of the rules: MT-FW with cap 6 and T3, with the S2 terms minus with the
S0 corrections, at 0.25, with origin time except in 2L-d.

| Rule | 7171 | 9393 | Met |
|:--|:--|:--|:--|
| 2L-a, standard | +0.057 [+0.034, +0.081] | +0.062 [+0.036, +0.083] | yes |
| 2L-b, dense | +0.038 [+0.024, +0.054] | +0.030 [+0.015, +0.044] | yes |
| 2L-c, lower bounds per regime (margin −0.02) | +0.0069 to +0.0517 | +0.0102 to +0.0461 | yes |
| 2L-d, epicentre only, standard / dense | +0.036 / +0.020 | +0.050 / +0.035 | yes |
| Secondary, BOMP-PK, standard | +0.068 [+0.026, +0.107] | +0.071 [+0.029, +0.114] | yes |

Descriptive: per cell, in standard windows, +0.108 and +0.113 in the one-sided 8-station cell, +0.048 and +0.027 with 6
and +0.021 and +0.018 with 4; median error of the terms from 0.147 and 0.151 s to 0.082 and 0.084 s with weak
clutter and from 0.216 and 0.207 s to 0.155 and 0.153 s with strong clutter; fraction of the oracle gain 0.56 and
0.61 in standard windows and 0.55 and 0.28 in dense ones. The independent verification recomputed S2 without access to the truth
with an identical result in 300 of 300 units, ruled out the origin-time artefact and reproduced the
statistics with its own code.

**Explorations C12, C13 and C14 (seed 999, exploratory; Appendix Y.29).** Registered detector; C12 with
16 networks per cell and regime, C13 with 12 and C14 with 12 for evaluation after freezing on 3 for development.

| Result | Value |
|:--|:--|
| C12, double couple − full, recall at 0.25 with origin time, standard / dense | +0.005 [−0.021, +0.033] / −0.004 [−0.018, +0.011] |
| C12, ratio of thresholds double couple / full, noise only, at 0.1 / 0.25 | 0.906 [0.883, 0.928] / 0.898 [0.886, 0.914] |
| C13, decoy D1p, realised FDR at q = 0.2, epicentre only / with origin time | 0.14-0.16 / 0.28-0.37 |
| C13, D1p, realised false events per window for c = 0.25, epicentre only, standard / dense | 0.21 [0.17, 0.25] / 0.23 [0.19, 0.28] |
| C14, false events relative to the event-free twin (q95), BASE and R2 with R3, standard / dense | 6.23 and 2.52 / 10.04 and 4.24 |
| C14, remedies − BASE, recall at 0.25 with origin time, standard / dense | +0.015 to +0.016 / −0.011 to −0.002 |
| C14, well-located detections with shifted origin, fraction of the false events in development, standard / dense | 55-58 % / 80-95 % |

Reading: the restriction to double couples and the self-interference remedies are not recall levers in this
benchmark; calibration by decoys works only with the epicentre-only rule; and the
false-event budget under the origin-time rule is spent by the origin of well-located detections.

**What has not yet been tested**

It remains to separate the sign of the radiation pattern from the polarisation constraint, which Task 2B removes
together; that ablation is not registered and would be exploratory. It remains to repeat the bounded-influence
certificate in a second generator, with a frozen operating point and with the hit rule with
origin time. It remains to close the margin of 0.08 to 0.24 below the realistic ceiling, which neither the
velocity model nor the pulse explains and of which about three quarters are untouched by any treatment
of amplitudes: explorations C9 and C10 place it in the accuracy of the station terms and in the
amplitude evidence; Task 2L confirmed with narrow scope a calibration part, the station
terms estimated without picks within the refined detector (G1-C), and what remains are their sensitivity to catalogue
errors of unknown size, their dependence on geometry, the origin time of the well-located
detections that exploration C14 points to and a 1 km detection grid. Of the closest comparators, MALMI and a block-sparse
inversion on a finer grid are missing; those of 2C and 2D, on the common grid, are our own implementations. A deployable operating
point for a false-event budget is missing, calibrated with windows that contain events; calibration by
decoys in C13 provides it label-free only with the epicentre-only rule and on one seed. Still missing are the
calibration of the threshold on real RSNC noise, a pruning of the certificate that pays off in the noise-dominated
iterations, level L3, the learned components and the real data.

## 7. What is new and what is not

A verified novelty audit found that the proposal overstated
the novelty of joint sparse detection and location with the moment tensor (E30), Task 2C
showed that it presented the convex problem as the method that performs when, on the grid, it is not (E33), and
Task 2D showed that MT-FW is no better than the block pursuit with the same time tolerance (E35),
Task 2E refuted that the off-grid BLASSO with shifts is (G1-M), and Task 2F
confirmed, with a scope that its verification narrowed (E37, E38), the effect of the per-station cap in
dense sequences, and exploration C8 showed that its mechanism is not the correction of a site term and
that the plan of a certificate with site terms has no support (E39), and Task 2I confirmed that a
Huber fit with level fixed by rule from the noise gives the same gain without a tuned constant, and
Task 2L confirmed, with narrow scope, that station terms estimated without picks raise the recall with origin
time of the refined detector. The following list incorporates all of this.

**Not new, and the thesis does not claim it:**

- Joint, sparse, moment-tensor-aware detection and location: the simultaneous
  recovery of origin time, hypocentre and moment tensor of simultaneous sources by block orthogonal
  matching pursuit over a dictionary of Green's functions (Vera Rodriguez, Sacchi and Gu, 2012).
- Inversion with moment-tensor blocks over a grid, greedy or convex: the greedy version is that
  of Vera Rodriguez et al. (2012), and the group lasso on the same grid is its direct convex
  relaxation. In Task 2C neither of the two performs like MT-FW on the coarse pilot grid with exact
  times; in Task 2D, with the same peak search, the greedy version matches MT-FW.
- MT-FW itself as a detector: with the same tolerance it is essentially block orthogonal matching pursuit
  with reading at the peak, and neither the greedy search, nor the deflation, nor the cap, nor the peak search are
  claimed as new (Task 2D).
- Clipping or normalising the contribution of each station: established practice in robust stacking and
  backprojection (per-trace normalisation, clipped or sign-based stacking).
- Per-station time shifts estimated together with the moment tensor of a known event
  (Zhu and Helmberger, 1996; Pham, 2024).
- The continuous scan with a linear moment tensor in regional networks at long periods (GRiD MT;
  Tsuruoka et al., 2009), operational since 2003.
- Joint migration of location and moment tensor, mechanism- or polarity-corrected stacks
  and radiation-corrected semblance in dense arrays (Chambers et al., 2014; Anikiev et al., 2014;
  Zhebel and Eisner, 2015; Staněk et al., 2015; Vlček et al., 2016).
- Per-station time corrections with a master event in stacking-based location
  (Grigoli et al., 2016).
- The theory of the subspace detector (Scharf and Friedlander, 1994; Harris, 2006).
- The BLASSO and sliding Frank-Wolfe (Bredies and Pikkarainen, 2013; Denoyelle et al., 2020).
- Chi-square random fields (Worsley, 1994).
- Extreme-value thresholds (Hirano et al., 2022).
- Neural eikonals and operators (EikoNet, MIFNO).

**New, as a combination and as a result**, as far as could be found:

- the off-grid formulation as a group BLASSO over measures with vector-valued moment-tensor
  amplitudes, in the 1 to 20 Hz regime and sparse networks; what it adds to block-sparse
  inversion is the continuous measure, the group norm on the radiated energy and its dual certificate. It is the
  framework and the certificate of the detector, not the detector that performs: its restriction to the grid, solved with
  certified optimality, falls below MT-FW in 23 of 24 cases (Task 2C), and the off-grid
  version with per-station and per-phase shifts does not beat BOMP-PK (Task 2E), so that as a detection
  method it is a negative result and as a framework it remains the object of the theory;
- the identity between that certificate and the generalised likelihood ratio of the matched subspace
  detector with the moment tensor as nuisance, applied to the whole network and to P and S; it is a modest
  unifying result, because both of its sides already existed;
- the certified search: an oracle with records identical to the exhaustive search and a rigorous bound
  on the floating-point error; its computational gain (4.3 times) comes from the certified vectorised
  evaluation, not from pruning, which with the energy bound does not pay off;
- the conformal thresholds: the calibration of the threshold of the maximum over the four-dimensional
  hypothesis space on noise-only windows, with a finite-sample guarantee under exchangeability and with the
  measurement of why the expected Euler characteristic is of no use for this statistic;
- the benchmark against GaMMA, PyOcto and QuakeMigrate with superposed events, at equal realised false-alarm
  rate and with preregistered rules (Task 2A), replicated against QuakeMigrate with 30 networks per
  cell and a hierarchical bootstrap (Task 2C), and with the closest antecedents with and without time
  tolerance (Tasks 2C and 2D);
- as an empirical result, that tolerance to time error is the deciding ingredient: with the same
  dictionary, the same grid of 3 km and 0.1 s and the same corrections, the methods read at the exact
  times lose most of their recall (Task 2C), and with the peak search the block pursuit
  matches MT-FW and both beat QuakeMigrate (Task 2D). The empirical advantage belongs to the family of
  time-error-tolerant moment-tensor detectors with deflation, not to MT-FW;
- as a preregistered negative result: the principled treatment of that tolerance, the group
  BLASSO with a time shift per station and phase as a nuisance parameter, solved off the
  grid by sliding Frank-Wolfe, does not beat BOMP-PK (Task 2E, +0.000 to +0.039, all CIs
  contain 0), and neither a physical velocity model nor a pulse-invariant subspace explains the
  remaining margin (exploration, seed 999);
- as a preregistered empirical result, the methodological contribution, modest: G1-R is supported by
  two preregistered synthetic tests, with an external timestamp and independent verification. A
  bounded-influence moment-tensor subspace certificate improves the recall at a fixed rate of
  false events in dense sequences of 3 to 5 overlapping events relative to BOMP-PK with equal time
  tolerance, without harming isolated events when the regimes are pooled: +0.043 to +0.074 at 0.25 with
  per-station amplitude clipping at a calibrated level (Task 2F; +0.024 to +0.045 with origin time) and
  +0.058 and +0.053 with a Huber fit whose level is fixed by rule from the noise (Task 2I; +0.031 and
  +0.019 with origin time), which removes the dependence on a tuned constant. It wins because it bounds the
  influence of each station-and-phase block in the face of the model's amplitude error and orders and chooses the atoms
  better, not by finding more events nor by separating them better (exploration C8). Neither clipping nor the
  Huber fit is new: Huber M-estimators are standard in moment-tensor inversion and
  in robust stacking. What is new is the design result measured within a moment-tensor detector
  with deflation, not an algorithm; per-station amplitude or site terms, also standard, do not
  replace the cap (exploration C8);
- with its limits: synthetic evidence from a single generator family; modest size, smaller with the
  strict hit rule; the Huber fit costs about 9 times more than the cap and stays somewhat below
  it at 0.1 on one seed; equivalence on isolated events is not met with weak clutter on one
  seed; with amplitude dispersion 0.8 the Huber fit falls below the cap (secondary arms
  of 2I, descriptive);
- as a preregistered empirical result on calibration, also modest: G1-C is supported by a
  preregistered synthetic test with independent verification. Station terms estimated without picks by
  a joint determination of hypocentres and terms on the 20-event catalogue raise the recall with
  origin time of the refined detector by +0.057 and +0.062 in standard windows and +0.038 and +0.030 in dense ones
  (Task 2L), with the gain concentrated in the one-sided 8-station cell, conditional on T3 and on catalogue
  errors of known size. The joint determination of hypocentres with station terms is classical
  in pick-based location: what is new is its pick-free estimation and its effect measured within a pick-free
  detector, not an algorithm;
- as an open problem: the margin up to the realistic ceiling, of which no treatment of amplitudes captures
  more than a quarter and which two explorations place in the calibration of the station terms and in
  the amplitude evidence (C9, C10); Task 2L closes part of the calibration (0.28 to 0.61 of the
  oracle gain), and exploration C14 places the false-event budget under the origin-time rule in the
  origin of well-located detections;
- with no novelty in what those explorations touched: soft time priors, probabilistic
  location with Gaussian model errors, the joint determination of hypocentres with station
  terms and polarisation or rectilinearity weights are known practice, and location without travel
  times (polarisation, amplitudes, envelope energy) has its own literature, which a non-systematic
  survey organised into eight families (`explore_contrib/LIT_NO_TRAVELTIME.md`; 26 of its 72 references
  read only at the metadata level). What is reported is the negative result with its bound and its mechanism.

Two more pieces of the thesis are not among what the novelty audit confirmed as new. The
per-station fidelity ladder, which decides how much physics to use according to the measured fidelity, is a
design proposal not yet measured. The identification of the curvature of the certificate at a source with the
efficient Fisher information, bounded below by the observability certificate weighted by the
radiated energy, is proved under idealised hypotheses (Proposition AA.6); the prediction of where
recovery is exact with several sources depends on a conjecture (Appendix AA).

Before claiming a "physical network subspace detector" one must read the LLNL report on
subspace detectors built with SPECFEM3D (LLNL-CONF-638234), which the state-of-the-art
survey could not open.

## 8. Experimental programme to outperform current methods

**Metric.** Per-event recall at a realised false-event rate of 0.1 and 0.25 per window (or per day on
continuous data), on the exact FROC curve, with paired bootstrap intervals that re-choose the operating
point in each replicate. The rules are registered before each run.

**Comparators**, all with their best configuration:

- phase association: GaMMA, PyOcto, REAL, GENIE and Harpa;
- migration: QuakeMigrate, MALMI and Qseek;
- end-to-end, where code exists: PhaseNO, EQNet and HEIMDALL;
- prior physics: block-sparse inversion of origin time, hypocentre and moment tensor
  (Vera Rodriguez et al., 2012), the mechanism-corrected stack (Anikiev et al., 2014; Chambers et al.,
  2014; Zhebel and Eisner, 2015) and GRiD MT.

So far GaMMA, PyOcto and QuakeMigrate have been run (Task 2A) and, in our own implementation on the
common grid, the block-sparse inversion and the mechanism-corrected stack, at the exact predicted
times (Task 2C) and with the peak search of MT-FW (Task 2D). The block inversion with peak
search (BOMP-PK) matches MT-FW in standard windows and is the comparator that every method of the thesis must
beat; the off-grid BLASSO with shifts does not beat it (Task 2E), and in dense windows MT-FW
with its cap does (Task 2F). MALMI, GRiD MT and a block inversion on a fine grid are missing.

**Stages.**

1. Consistent synthetic data with double couples, medium perturbations and site factors (done; Tasks
   1w to 1y), and a generator mismatched with the detectors, STA/LTA picker, station corrections and
   QuakeMigrate (done; Task 2A, preregistered and verified: the primary rule is met in both
   clutter regimes). The ablation of the radiation pattern within the same solver (Task 2B,
   preregistered) was executed: it is met in both regimes, with weak clutter only at the minimum. The
   closest comparators, the convex problem on the grid and the replication with 30 networks per cell
   (Task 2C, preregistered) were too: 2C-a and 2C-c are met in both regimes and the convex problem on
   the grid falls below MT-FW in 23 of 24 cases. The equality of tolerance, the ablation of the
   peak search and the frozen threshold (Task 2D, preregistered with hash) were too: with the same
   peak search the block pursuit matches MT-FW (2D-a is not met), the peak search matters
   (2D-b is met) and the conformal threshold on empty windows is not deployable. The off-grid convex problem
   with per-station and per-phase time shifts (Task 2E, preregistered) does not beat
   BOMP-PK (H5 refuted), and the per-station cap in dense sequences (Task 2F, preregistered with new
   seeds) does beat it (H6). The bounded-influence certificate with level fixed by rule from the
   noise (Task 2I, preregistered, seeds 5151 and 6262) also beats it, without being inferior to the cap and without
   changing the recall of isolated events when pooling (H7); its secondary arms show an indifferent level
   between 2.5 and 3.4 noise deviations and a Huber fit below the cap with amplitude dispersion
   0.8. Within this same synthetic stage there remain its repetition in a second generator with a frozen operating
   point and the hit rule with origin time, a cheaper solver for the Huber fit, and the
   study of the margin, which explorations C9 and C10 split into calibration and amplitude. The station
   terms estimated without picks within the refined detector (Task 2L, preregistered, seeds 7171 and 9393)
   raise the recall with origin time (H8), with the gain concentrated in the one-sided 8-station
   cell; the next steps are its sensitivity to the catalogue, its dependence on geometry, its
   combination with the bounded-influence statistic and, according to exploration C14, the origin time of
   the well-located detections. They are the first mathematical deliverable.
2. Real phase detector (PhaseNet) and real RSNC noise added to SPECFEM3D synthetics in a
   3D model of Colombia (Tasks 1p and 1v).
3. Real data by decimation of a dense part of the RSNC against the SGC catalogue, with template
   verification of every uncatalogued detection (Section 16).
4. A public reference region with a high-resolution catalogue, so that the comparison is
   reproducible by third parties.

**Falsifiable hypotheses, in order of decision.**

- H1: the dual certificate at L1 or L2 beats its control with the geometry of another hypothesis and
  direction-free migration. In synthetic data, rule 2A-b is met in both regimes against the control
  with the full geometry of another node, and rule 2B, which replaces the radiation operator by the
  identity within the same solver, is met in both: in strong clutter in 3 of 3 cells per seed
  (+0.06 to +0.35), in weak clutter at the minimum of 2 of 3 and in ~60 % of the bootstrap replicates. The
  gain is attributed to the radiation operator, clearly in strong clutter and fragilely in
  weak clutter; which part corresponds to the sign and which to the polarisation constraint remains unseparated.
- H2: it beats the best associator and the best migration-based locator by at least 0.05 of recall at
  both false-event rates. In synthetic data, its confirmatory test is rule 2A-a of Task 2A, which is
  met in both regimes against GaMMA and PyOcto with STA/LTA picks and against QuakeMigrate, by +0.09 to
  +0.49 at 0.1 and +0.17 to +0.46 at 0.25, with window-level intervals; with a network-level bootstrap four cells
  include zero at 0.1. Against an associator with optimistic picks it holds with strong clutter and
  ties in one cell with weak clutter. With 30 networks per cell and a hierarchical bootstrap (rule 2C-c) the
  advantage over QuakeMigrate replicates in both regimes, +0.22 to +0.35 with weak clutter and +0.099 to
  +0.295 with strong clutter at 0.25, with lower bounds above 0 in all cells; with Holm the 4-station
  cell with strong clutter of seed 7331 does not reject. Against the block inversion and the
  mechanism-corrected stack on the common grid (rule 2C-a) it is met in both regimes, by +0.23 to +0.52
  and +0.08 to +0.45 at 0.25, but as robustness to grid and time errors, not as superiority over
  those methods on a fine grid; with the same peak search (rule 2D-a) the block pursuit
  matches MT-FW, so H2 holds for the family of time-error-tolerant detectors with
  deflation, not for MT-FW against the block pursuit. Against MALMI it has not been tested.
- H3: deflation with a certified oracle beats the single-pass detector in windows with superposed
  events. In Task 2A (exploratory), part of the margin over QuakeMigrate, whose trigger merges
  maxima close in time, comes from the windows with two sources; without them the minimum margin drops to
  +0.06 at 0.1 over QM-3C, QM-3C+SC and GaMMA-STALTA+SC. It is not a test of H3, because QuakeMigrate
  also differs in the statistic. In Task 2B (exploratory), the same solver without the radiation
  operator, which keeps the deflation, beats the best comparator of 2A at 0.25 in 5 of 6 cell-seeds
  with weak clutter and stays at its level with strong clutter: with weak clutter pick-free migration with
  deflation explains a good part of the margin over QuakeMigrate.
- H4: the realised false-alarm rate on real noise stays within the interval of the rate predicted
  by conformal or GEV calibration on noise-only windows of the same regime. In synthetic data it is
  met on Gaussian noise and, recalibrated, with strong clutter; the expected Euler characteristic
  does not meet it (Task 1z). In Task 2A an absolute threshold frozen on other networks realised up to 2.3
  times its budget in one cell (2A-e): the transfer of the threshold between networks is not guaranteed.
  In Task 2D the split-conformal threshold frozen on empty windows of other networks stayed below
  2α on empty windows except in one cell, but it does not control the false events in windows with events (1.2 to 2.2
  per window with weak clutter).
- H5: the convex problem, solved with certificate, matches MT-FW (|diff| < 0.05) or beats it. On the
  pilot grid, no: the descriptive rule 2C-b puts the group BLASSO below MT-FW in 23 of 24
  cases, by −0.064 to −0.605, and at the level of the block pursuit or somewhat below (Task 2C). Off
  the grid it is open: in the ideal generator with the sources on the nodes (exploratory) the BLASSO
  reaches 0.73 against 0.76 for MT-FW at 0.25 with 4 stations, which points to discretisation as the cause, but only a sliding
  Frank-Wolfe, or a formulation with time tolerance, on the same windows decides it. With the
  heuristic peak search on the grid (GL-FW-PK, Task 2D) it remains 0.11 to 0.30 below in 11 of 12
  cases. The principled form, with per-station and per-phase shifts within the problem and solved
  off the grid (Task 2E), gains +0.000 to +0.039 over BOMP-PK with all CIs containing 0: H5
  is refuted in synthetic data. The convex problem remains as framework and certificate, and the detector is a
  tolerant greedy search (MT-FW or BOMP-PK).
- H6: the robust certificate with per-station cap beats BOMP-PK in dense sequences of 3 to 5 events, and
  the gain is due to the cap (G1-R). It is met in Task 2F, preregistered with new seeds: +0.043 to
  +0.074 at 0.25 in both regimes and seeds, also with Holm; without the cap MT-FW matches BOMP-PK; in
  standard windows there is equivalence with margin 0.05. Scope according to the verification: +0.024 to +0.045 with
  origin time, tied to the level of the cap, through bounded influence of each block in the face of the amplitude
  error. The hypothesis that modelling per-station site terms beats the cap has no support in
  exploration C8.
- H7: a per-station and per-phase Huber fit with threshold fixed by rule from
  the noise beats BOMP-PK in dense windows, is not inferior to the calibrated cap and does not change the recall of
  isolated events. It is met in Task 2I, preregistered with seeds 5151 and 6262: +0.058 and +0.053 at
  0.25 over BOMP-PK; −0.008 and −0.002 relative to the cap, with lower bounds above −0.03; equivalence with
  margin 0.05 in standard windows. Scope according to the verification: +0.031 and +0.019 with origin time;
  equivalence not met with weak clutter on one seed; it orders and chooses the atoms better, it does not find
  more events; about 9 times the cost of MT-FW. Secondary arms, descriptive: indifferent level between 2.5
  and 3.4 noise deviations; with amplitude dispersion 0.8, below the cap by −0.055 and −0.069.
- H8: station terms estimated without picks from the catalogue, within the refined detector, raise the
  recall with origin time relative to the registered corrections (G1-C). It is met in Task 2L,
  preregistered with seeds 7171 and 9393: +0.057 and +0.062 in standard windows and +0.038 and +0.030 in
  dense ones, with no loss by regime nor with the epicentre-only rule. Scope according to the verification:
  concentrated in the one-sided 8-station cell, with no demonstrated gain with 4, terms worse than the
  registered ones in about 22 % of the networks, conditional on T3 and on catalogue errors of known size,
  0.28 to 0.61 of the oracle gain. Self-calibration with the detections added nothing measurable in
  exploration C11. Negative outcome: if the calibration does not hold outside the known catalogue and the
  one-sided geometry, the margin is reported as limited by the available calibration and by amplitude.

Each one has its negative outcome written down: if H2 fails on real data, the thesis delimits at which level
of the ladder and with what medium fidelity physics stops paying off.

## 9. Risks and exits

**Risk 1: amplitude fidelity at high frequency.** It is the reason why moment-tensor migration
remained confined to dense arrays or long periods. Exit: the per-station ladder
(L1 does not use amplitudes) and the bounded influence of each station-and-phase block. Task 2F measured it in
synthetic data: in dense sequences the blocks with poorly predicted amplitude dominate the uncapped statistic, and the
per-station cap protects it; its form with level fixed by rule from the noise gives the same gain
in Task 2I, also in synthetic data.

**Risk 1b: a cap tied to a constant and to the generator's amplitude model.** The gain of the cap
peaks at the value fixed in 1x, reverses with caps of 3 or less, disappears above 12 and drops
to +0.023 (CI containing 0) without the amplitude factors, which the generator draws in each window; with dispersion
0.25 the cap still gains +0.043 (CI containing 0) where the true factors gain nothing, and the catalogue terms
do not replace it (exploration C8). Task 2I removes the dependence on the constant with nominal
dispersion: the Huber fit with level fixed by rule is not inferior to the cap on two never-used seeds.
The dependence on the generator remains: both tasks use a single generator family, and in the
descriptive arms of 2I, as in the pilot, the Huber fit falls below the cap with dispersion 0.8
(−0.055 and −0.069), while its level is indifferent between 2.5 and 3.4 noise deviations. Exit: keep the value of 1x in every confirmatory test, repeat in a second
generator and, on real data, report both the capped certificate and the Huber one.

**Risk 1c: the false-event budget under the origin-time rule.** Exploration C14 found that, at the operating
point, the largest class of false events under the origin-time rule are well-located detections of real
events with the origin shifted by a median of 0.8 s, which no deflation remedy reduces and which
no moveout decoy imitates (C13). Exit: treat the origin time of the well-located detections
as an object of calibration (S2 terms with T3, a better time model, the early bias of the origin),
explore it before registering anything and always report both hit rules.

**Risk 2: computational cost.** Pruning with the energy bound does not reduce it (Task 1z). Exit: the
certified vectorised evaluation (4.3 times faster than the dense loop, with identical records),
cell bounds that use the subspace, the amortised proposals and parallelism by distributed octrees
(p4est), only at national scale.

**Risk 3: the wavelet.** The pilot uses the known wavelet as a matched filter, which favours
every migration. Exit: learned signed characteristic functions. In Task 2A, with a wavelet of variable
central frequency and phase and scattered coda, the advantage held in synthetic data.

**Risk 4: the non-Gaussianity of real noise.** The analytical threshold is not calibrated even on
Gaussian noise, and strong clutter breaks every threshold calibrated on Gaussian noise (E28). Exit: the
conformal or extreme-value calibration on noise-only windows of the same regime, recalibrated
every time the regime changes, which is part of the method and not an add-on. Its guarantee requires
exchangeable, event-free windows, a real assumption for non-stationary recordings.

## 10. Verified references of this document

Those of the proposal [109] to [112] are also used.

- K. Bredies and H. K. Pikkarainen, *ESAIM: COCV* 19(1), 190-218, 2013, doi:10.1051/cocv/2011205.
- Y. de Castro and F. Gamboa, *J. Math. Anal. Appl.* 395(1), 336-354, 2012, doi:10.1016/j.jmaa.2012.05.011.
- E. J. Candès and C. Fernández-Granda, *Comm. Pure Appl. Math.* 67(6), 906-956, 2014, doi:10.1002/cpa.21455.
- V. Duval and G. Peyré, *Found. Comput. Math.* 15(5), 1315-1355, 2015, doi:10.1007/s10208-014-9228-6.
- Q. Denoyelle, V. Duval, G. Peyré and E. Soubies, *Inverse Problems* 36(1), 014001, 2020, doi:10.1088/1361-6420/ab2a29.
- N. Boyd, G. Schiebinger and B. Recht, *SIAM J. Optim.* 27(2), 616-639, 2017, doi:10.1137/15M1035793.
- M. Jaggi, *Proc. ICML*, PMLR 28(1), 427-435, 2013.
- V. Chandrasekaran, B. Recht, P. A. Parrilo and A. S. Willsky, *Found. Comput. Math.* 12(6), 805-849, 2012, doi:10.1007/s10208-012-9135-7.
- L. L. Scharf and B. Friedlander, *IEEE Trans. Signal Process.* 42(8), 2146-2157, 1994, doi:10.1109/78.301849.
- D. B. Harris, "Subspace detectors: theory", LLNL UCRL-TR-222758, 2006.
- K. J. Worsley, *Adv. Appl. Probab.* 26(1), 13-42, 1994, doi:10.2307/1427576.
- J. E. Taylor, *Ann. Probab.* 34(1), 122-158, 2006, doi:10.1214/009117905000000594.
- *Random Fields and Geometry*, Springer, 2007, doi:10.1007/978-0-387-48116-6.
- S. A. Piyavskii, *USSR Comput. Math. Math. Phys.* 12(4), 57-67, 1972, doi:10.1016/0041-5553(72)90115-2.
- A. N. Angelopoulos and S. Bates, *Found. Trends Mach. Learn.* 16(4), 494-591, 2023, doi:10.1561/2200000101.
- Y. Hirano, H. Kawakata and I. Doi, *Sci. Rep.* 12, 22090, 2022, doi:10.1038/s41598-022-25839-2.
- D. Chambers et al., *Geophys. Prospect.* 62(4), 879-896, 2014, doi:10.1111/1365-2478.12108.
- D. Anikiev et al., *Geophys. J. Int.* 198(1), 249-258, 2014, doi:10.1093/gji/ggu126.
- A. Zhebel and L. Eisner, *Geophysics* 80(1), KS1-KS9, 2015, doi:10.1190/geo2014-0055.1.
- F. Staněk et al., *Geophys. J. Int.* 201(3), 1362-1369, 2015, doi:10.1093/gji/ggv070.
- S. Tsuruoka, H. Kawakatsu and T. Urabe, *Phys. Earth Planet. Inter.* 175(1-2), 8-16, 2009, doi:10.1016/j.pepi.2008.02.014.
- I. Vera Rodriguez, M. Sacchi and Y. J. Gu, *Geophys. J. Int.* 188(3), 1188-1202, 2012, doi:10.1111/j.1365-246X.2011.05323.x.
- F. Grigoli, S. Cesca, L. Krieger, M. Kriegerowski, S. Gammaldi, J. Horalek, E. Priolo and T. Dahm, *Sci. Rep.* 6, 25744, 2016, doi:10.1038/srep25744.
- G. Guilhem Trilla and Y. Cano, *Pure Appl. Geophys.* 182, 4939-4958, 2024, doi:10.1007/s00024-024-03455-7.
- W. Zhu, S. M. Mousavi and G. C. Beroza, *IEEE Trans. Geosci. Remote Sens.* 57(11), 9476-9488, 2019, doi:10.1109/TGRS.2019.2926772.
- J. D. Smith, K. Azizzadenesheli and Z. E. Ross, *IEEE Trans. Geosci. Remote Sens.* 59(12), 10685-10696, 2021, doi:10.1109/TGRS.2020.3039165.
- F. Lehmann, F. Gatti and D. Clouteau, *J. Comput. Phys.* 527, 113813, 2025, doi:10.1016/j.jcp.2025.113813.
- M. Puente Huerta, C. Sippl, J. Münchmeyer and I. McBrearty, *Seismica* 4(2), 2025, doi:10.26443/seismica.v4i2.1559.
- W. Zhu et al., GaMMA, *J. Geophys. Res.* 127(5), 2022, doi:10.1029/2021JB023249.
- I. W. McBrearty and G. C. Beroza, GENIE, *Bull. Seismol. Soc. Am.* 113(2), 524-547, 2023, doi:10.1785/0120220182.
- J. Münchmeyer, PyOcto, *Seismica* 3(1), 2024, doi:10.26443/seismica.v3i1.1130.
- H. Sun, Z. E. Ross, W. Zhu and K. Azizzadenesheli, PhaseNO, *Geophys. Res. Lett.* 50(24), e2023GL106434, 2023, doi:10.1029/2023GL106434.
- M. Isken et al., Qseek, *Seismica* 4(1), 2025, doi:10.26443/seismica.v4i1.1283.
- C. Shi et al., Harpa, *Nat. Commun.* 17, 8140, 2026, doi:10.1038/s41467-026-74092-y.
- L. Zhu and D. V. Helmberger, *Bull. Seismol. Soc. Am.* 86(5), 1634-1641, 1996, doi:10.1785/BSSA0860051634.
- T.-S. Pham, *Geophys. J. Int.* 238(2), 783-793, 2024, doi:10.1093/gji/ggae188.
