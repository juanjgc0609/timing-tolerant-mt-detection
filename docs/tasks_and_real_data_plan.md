# Tasks that require a real execution environment

This file lists what **cannot be done** in the environment where the sixteen pre-studies of the
document (Appendices G to V) were produced, and which therefore remains pending for your machine. The reasons are
concrete: no network access to FDSN services or to waveform repositories, no GPU, no
elastodynamic solver installed, and no capacity to interactively debug long numerical code.

Each task comes with its objective, why it could not be done here, the required environment, steps from start to
finish, expected output, acceptance criterion and estimated effort. They are in order of priority. Tasks
0 and 1 are closed in their minimal part; Tasks 1b to 1n arise from their results and from the
audit (`audit_status.md`) and come first because they decide how much of the polarization channel
exists in real data and in what form G1 holds. Tasks 1j and 1k are synthetic and were
executed in this environment (Appendix Y of the proposal).

---

## Base environment, common to almost all tasks

```bash
conda create -n sismo python=3.11 -y
conda activate sismo
pip install numpy scipy matplotlib pandas
pip install obspy                     # seismic I/O and FDSN clients
pip install seisbench                 # pretrained PhaseNet / EQTransformer pickers
pip install torch --index-url https://download.pytorch.org/whl/cu121   # adjust to your CUDA
pip install pyproj tqdm h5py
```

Check that there is access to data:

```python
from obspy.clients.fdsn import Client
c = Client(base_url="http://sismo.sgc.gov.co:8080")   # the CM network is archived at the SGC, not at IRIS
print(c.get_stations(network="CM", level="station"))
```

If that responds, the tasks with real data are executable. IRIS/EarthScope does not hold the continuous data of
the CM network (verified in Task 1).

---

## TASK 0 — COMPLETED

Executed. GaMMA 1.2.17 and PyOcto 0.2.0 run on the same synthetic picks, with a sweep of the
operating point, positive controls, exact paired McNemar and **two independent seeds**.
Results in Appendix X of the proposal; outputs in `v6_*.txt` and `v7_*.txt`.

**What it established.** Statement (i) —the constraints add value over the same search without them—
is established: zero false events in sixty realizations in one-sided networks against 0.30 to
0.43 for the ablation, with unanimous discordance and $p\le0.0039$. Statement (ii) is established
in the sparse regime: with three stations GaMMA recovers the event in 7 and 8 of 30 depending on the
seed, PyOcto in 8 and 9, against 29 for the constrained prototype.

**What it cost, and this is what must be remembered.** Three protocols. The first two gave the
opposite result and **in neither was the cause the physics**: an absolute threshold on scores of different
scale, and a single geometry per cell when geometry dominates over the method. Appendix S
had already measured that observability varies 86× between geometries; it was foreseeable and was not foreseen.

**What remains open.** The curve metric `falsos@recall` is not reached in several cells of the
ladder, and the decoy cell at 12 km is the only one where the two seeds disagree. `task0_curves.py`
locates a sign change at threshold fraction ≈0.8, which suggests that the polarization channel
has a regime of validity and not a uniform effect. Characterizing it is pending and cheap work.

---

## TASK 1 — COMPLETED

Executed on the CM network (SGC, HH at 100 Hz), ISC reviewed catalogue, PhaseNet `instance` at threshold
0.30, at five stations. Result in Appendix T.1 of the proposal; code in `task1_fetch.py`,
`task1_analyze.py`, `task1_decide.py` and `diagnose_polarization_config.py`.

**What it established.** The registered question has an answer: where polarization can be measured, the
real clutter behaves like random clutter. At ZAR (90 days from 2024-10-01, band 1-3 Hz,
centred 2 s window) the chance-corrected fraction of coherent clutter is
$\eta=0.054\ [-0.049,0.164]$ over 159 false picks, below the threshold of 0.20. By the registered
criterion, Appendices H, J and M keep their quantitative validity **at the stations where the
polarization channel is usable**.

**What was not registered and weighs more.** The prior control (does the P polarization of a single
station point to catalogued events?) passes at one station out of five:

| Station | control P picks | fraction | status |
|:--|:--|:--|:--|
| ZAR | 109 | 0.789 | passes |
| CHI | 78 | 0.256 (chance 0.256) | fails; best of 24 configurations 0.385 |
| RUS | 45 | 0.244 | fails (earlier protocol, P and S mixed) |
| OCA | 33 | 0.303 (chance 0.237, p=0.23) | fails; not distinguishable from chance |
| SJC | 19 | 0.474 (chance 0.273, p=0.044 uncorrected) | inconclusive: fewer than 30 picks |

At CHI the correctable causes were ruled out one by one: standard metadata in a single epoch since
2019, no fixed rotation (0.321) nor fixed reflection (0.282, which is what an inverted or swapped N/E
channel produces) reaches 0.60, the fraction is stable over the four quarters of the
period (0.25-0.26), which excludes an instrument change, and restricting to P picks within
1 s of the predicted P only raises it to 0.348 (n=46), which excludes misassociation. The same code
passes at ZAR and in synthetic data. It remains a property of the station (its installation or its local
structure, which the tests do not separate) **in the 1-10 Hz band tested**; the primary band of the
proposal (0.1-2 Hz) has not been tested at the stations that fail. It remains unexplained that the S picks
of CHI give 0.677 (n=65), well above its P picks.

**What it cost, and this is what must be remembered.** None of the defects found was physical, and several
would have produced a false verdict with an appearance of rigour:

- The registered threshold of 0.20 on the raw fraction lay **below what clutter gives
  without any coherence**: its expected value is the chance hit rate, 20/90 = 0.22 with
  uniform azimuths and 0.24-0.33 in the measured sets. Passing that threshold or not measured the geometry,
  not the coherence. The correct criterion is on the chance-corrected fraction
  $\eta=(f-f_0)/(f_{\rm ctrl}-f_0)$, with $f_0$ exact by within-station permutation.
- The reference event for the clutter was the closest one in **origin time**. With one event every
  28 minutes in the four-station catalogue, it was usually another event, distant or later than the pick,
  and coherence was measured against something unrelated. Correct reference: the most recent event within
  200 km whose predicted P precedes the pick.
- The control included S picks, whose motion is not bound to point to the event (SH is
  transverse): they do not test the estimator.
- The centred 0.6 s window is half pre-arrival noise and less than one cycle at 1-2 Hz. At ZAR
  the 2 s window raises the control from 0.565 to 0.754.
- A control with n=11 gave 0.636 at CHI and seemed to pass; with n=78 it gives 0.256. A bootstrap CI with n=1
  gave [1.000, 1.000] and would have produced "HIGH coherence". Every gate needs a minimum
  sample and an "inconclusive" label distinct from "fails".
- The sweep of 72 fixed transformations found a "165° rotation" at SJC with n=15 that
  disappeared with n=19. Searching for the best of many configurations on few data fabricates signals;
  validation on disjoint halves of events is what separates them.
- Loading 90 continuous days of four stations exhausts memory. The analysis now reads only the
  segments it picks, per station and in batches of 6 h.
- Two looks at ZAR: the first (10 coda windows per event, $\eta=0.044\ [-0.094,0.200]$)
  ended exactly at the limit; the second (20 windows) was declared final before running it. No
  correction was applied for the intermediate look, and this is declared.

**What remains open.** What fraction of the CM network has usable polarization and why; whether the 0.1-2 Hz
band rescues the stations that fail; how much of the advantage of G1 survives with the
measured reliability; and whether treating C-IV according to the reliability of each station recovers it. These are
Tasks 1c to 1i; 1e has already been executed in synthetic data.

---

## TASK 1b — Regime of validity of the polarization channel

**Priority: high.** One week, no external data. Comes from Task 0.

`task0_curves.py` locates a sign change at threshold fraction ≈0.8 in the decoy cell at
12 km, which suggests that the polarization channel has a regime of validity and not a uniform
effect. Characterizing it closes the only thing left open from Appendix X. Task 1e adds a
second axis to the same regime: the fraction of stations whose polarization is usable.

---

## TASK 1c — Reliability map of polarization in the CM network

**Priority: maximum.** Cheap, and everything that follows depends on its result.

### Objective

Measure what fraction of the HH stations of the CM network pass the polarization control, in which band,
and what it depends on. Task 1 measured it at five stations and one passed; with five it is not possible to say
what fraction of the network is usable.

### Why it is cheap

The control only uses a few seconds around each catalogued P arrival. `task1c_network_survey.py`
downloads only windows of ±60 s around each predicted P, some 50 times less data than the
continuous data of Task 1, picks with PhaseNet within each window and measures the polarization on the
grid of bands and windows. Tested in synthetic data: it exactly reproduces the number given by the
diagnostic on continuous data (0.603 with n=78 in both).

### Steps

```bash
python task1c_network_survey.py --start 2024-10-01 --days 90 --grid todas
python task1c_network_survey.py --offline --grid propuesta     # re-analyses without downloading
```

The download is resumable: if it is interrupted, the same command is run again and only what is missing
is downloaded.

### Output

`t1c_survey/survey_summary.csv`, one row per station: status (passes, fails, inconclusive), reference
fraction against chance with its p, best configuration with its two halves, best fixed rotation and
reflection, fraction per quarter, and fraction by magnitude and by depth of the event.

### Acceptance criterion

At least twenty stations with n ≥ 30. Report the fraction of the network that passes with a
Clopper-Pearson interval. **The stations that pass must be confirmed in a second independent period**
before being used in Task 1g: with tens of stations and tens of configurations, some
will pass by chance.

### Estimated effort

Two to four days of download (the SGC server is slow) and one day of analysis.

---

## TASK 1d — The primary band of the proposal (0.1-2 Hz)

**Priority: maximum.** Hours, on already downloaded data.

### Objective

Section 1.5 declares 0.1-2 Hz as the primary band. Task 1 tested only between 1 and 10 Hz, so
the claim that CHI fails because of a site effect holds in that band and not in that of the proposal.

### Steps

```bash
python diagnose_polarization_config.py --stations CHI --grid propuesta
python diagnose_polarization_config.py --stations OCA --grid propuesta
```

It requires `picks_polarization.csv` to be the one for those stations (run beforehand
`task1_analyze.py --stations CHI,OCA --noise-per-event 10`). Section A7 of the same diagnostic
additionally separates signal/noise (fraction by magnitude) from incidence (fraction by depth).

### What to expect, in both directions

Scattering by shallow heterogeneity decreases with wavelength, so if CHI fails because of
scattering near the surface the low band should mitigate it. Against: microseism dominates
the noise between 0.1 and 0.5 Hz, and a magnitude 2.5 event has little energy below 1 Hz.

### Acceptance criterion

If CHI passes in 0.1-2 Hz with n ≥ 30 and in both halves, its failure is one of band, Appendix T.1 is
corrected, and C-IV runs in the primary band. If it does not pass, the site effect is confirmed in the
band of the proposal.

---

## TASK 1e — G1 with the measured reliability — EXECUTED IN SYNTHETIC DATA

Run in the environment of the proposal; result in Appendix T.2. Code
`task1e_polarization_reliability.py`, outputs in `results/task1e/`.

**Later caveat (Task 1j).** The advantages of this task are measured against the TOL ablation, whose
false events in these cells are largely ties of its discrete score. What 1e establishes
solidly is the ordering between treatments of a failing station (neutral no worse than C-III; naive C-IV
harmful; switching it off, unstable), not the size of the polarization advantage.

**What was done.** The final protocol of Appendix X repeated changing only the polarization of the
true arrivals per station: good (azimuth with 10° error, like ZAR), bad (random azimuth,
incidence preserved) or bad-3D (random direction). Fraction of good stations 1, 1/2, 1/4 and 0.
Four treatments of C-IV: at all stations (PHYS, as in the pre-studies), switched off at those
that fail (PHYS-G), neutral value at those that fail (PHYS-N), and the TOL ablation. With all
stations ideal, the generator and the score fields are bit-for-bit identical to those of the code
of Task 0 on the same realizations, and two runs with the same seed give identical results
in the 84 cell-method combinations (`results/task1e/verification.txt`). The rates
match Appendix X in the one-sided cell and only within sampling variation in those with
few stations: the runs of Appendix X cannot be regenerated realization by realization with
the archived code.

**What it established**, in false events out of 30 (one-sided with clutter 40, and four stations):

- The scatter of a good station costs nothing: 1 against 10 for the ablation, and 9 against 25.
- With stations whose direction is noise, naive C-IV **worsens** relative to the ablation: with no
  good station, the recall with four stations drops to 14 and 18 against 30 and 28 ($p<10^{-4}$ and
  $p=0.006$), and in the one-sided network it gives 17 and 27 false events against 10 and 13, in the same direction in both
  seeds but significant in one ($p=0.12$ and $p=0.001$).
- Switching off C-IV where it fails **does not fix it**: with a quarter of good stations it gives 20 and 19 false events
  against 10 and 13 for the ablation ($p=0.031$ and $0.15$), and in the rest of the mixed cells its direction
  against the ablation is not consistent. A hypothesis explained only by failing stations receives
  a polarization factor of 1, perfect agreement: absent evidence is treated as agreement.
- The **neutral value recovers the advantage** with a quarter of good stations, close to the one
  in five measured in Task 1: 2 and 3 false events against 10 and 13 ($p=0.039$ and $0.006$), 14 and 19 against
  25 and 28 ($p=0.003$ and $0.012$), with recall between 28 and 30 of 30. Replicated in an independent
  seed.

**What it cost.** PHYS-G was the first idea and failed. PHYS-N was designed after seeing it fail on
the same realizations, so its result in that seed is exploratory. That is why
a second seed was run with PHYS-N already fixed: five of six mixed cells significant against the
ablation, the sixth in the same direction. And a result that seemed firm did not replicate: the loss of
naive C-IV when the bad stations preserve incidence depends on the seed (10 and 0
false events), so it is not claimed.

**What remains open.** The neutral value is a rule, not a likelihood: Task 1h
replaces it with the principled version. Three stations remain saturated in false events for
all methods (Appendix X.5).

---

## TASK 1f — Robustness of the clutter verdict by subpopulation

**Priority: high.** One hour.

### Objective

The ZAR verdict is a single number over all the clutter. `task1_clutter_robustness.py`
recomputes it by picker confidence, by delay in the coda, by distance and by assigned phase. If
any subpopulation ends up with its whole interval above 0.20, the verdict depends on how the
clutter was sampled and it must be reported.

### Steps

```bash
python task1_analyze.py --noise-per-event 20 --stations ZAR --band 1 3 --win 2 --win-align centered
python task1_clutter_robustness.py --stations ZAR
```

Tested in synthetic data: it detects coherent clutter hidden only in the early coda
($\eta=0.94$ there, $-0.08$ in the late coda).

### Acceptance criterion

No subpopulation with n ≥ 30 with the interval of $\eta$ entirely above 0.20. There are eight comparisons
without correction; one that barely grazes the threshold is what is expected by chance.

---

## TASK 1g — Pilot of G1 with real data (the ambitious one)

**Priority: high, after 1c, 1h, 1l and 1m.** It is the strongest evidence the thesis can
have. Definitive design: decimation of a dense network (Section 16 of the proposal).

### Objective

Everything that supports statement (i) —the constraints reduce false events without losing
recall— is measured in synthetic data. This pilot tests it on **real picks with real clutter**,
at the stations that pass the control.

### Prerequisites

- A dense zone of the RSNC with a complete SGC catalogue above its $M_c$: the catalogue of the full
  network is the truth, and the subnetworks of 3, 4 and 6 stations (including one-sided subnetworks)
  are formed by selecting stations.
- From Task 1c, the reliability $\hat r_s$ of each station in the zone, after the orientation
  audit of Task 1l.
- From Task 1h, the polarization likelihood calibrated per station, or at minimum the
  neutral treatment PHYS-N of Task 1e, which is the only one of those tested that does not worsen relative to the
  ablation with failing stations.

### Design, to be preregistered before running

1. **Data.** Continuous data of the chosen stations in a period different from the one used to choose them.
   PhaseNet at a low threshold: the clutter is the real one, not a synthetic one.
2. **Truth.** ISC reviewed catalogue and SGC catalogue. A detected event without a match in the
   catalogue **is counted neither as false nor as true** without independent verification
   (template correlation or analyst review), as Section 16 requires. Three classes:
   matches the catalogue, verified false, unverified.
3. **Methods, on the same picks.** Ablation with tie-breaking (TOL-TB), C-III with the
   likelihood of 1h, GaMMA, PyOcto and REAL tuned in a calibration period, and the NET-VISA-style
   azimuth likelihood baseline. **Absolute threshold**, calibrated in the
   calibration period; the threshold relative to the maximum of Appendix X does not work on continuous data, because
   a window without an event would always emit one (Appendix Y).
4. **Geometry and medium.** Local projection in km, grid over the region, regional 1-D model for
   the associators that need it; C-III only needs the velocity bound.
5. **Endpoint.** Recall curve (events with $M\ge M_c$) against verified false events per
   day; the operating point is read off the curve. The positives of the template correlation
   are validated with the false-positive protocol of Muir et al. (2023), ref. [99].

### Acceptance criterion

Fixed before running, in the G1' form of Appendix D: at least 0.05 fewer false events per
window (or its equivalent per day) than the best comparator, including the tuned published baselines,
at non-inferior recall (margin 0.02), in at least half of the subnetworks. Any result is reported: a negative one delimits in real data where the physics does not pay.

### Estimated effort

Three to five weeks, almost all of it in the verification of the detections without catalogue.

---

## TASK 1h — Polarization likelihood calibrated per station

**Priority: high.** It turns the limitation of Task 1 into a methodological contribution.

### Objective

Replace the heuristic factor of C-IV (mean squared cosine) with a per-arrival likelihood
with parameters of each station estimated with catalogued events:

$$p_{1,s}(\delta) = r_s\,p_{\mathrm{vM}}(\delta;\kappa_s) + (1-r_s)\,p_{0,s}(\delta),
\qquad p(\delta\mid\text{clutter}) = p_{0,s}(\delta),$$

where $\delta$ is the folded angle between the measured and the predicted direction and $p_{0,s}$ is the
geometric null density of the station, the same one whose mass below $20^\circ$ is the chance rate
$f_{0,s}$ of Task 1 (not the uniform $1/\pi$, because events cluster in azimuth). The
reliability is estimated with the same algebra as $\eta$:
$\hat r_s=(f_{\mathrm{ctrl},s}-f_{0,s})/(q_s-f_{0,s})$, with $q_s$ the mass of $p_{\mathrm{vM}}$ below
$20^\circ$. The score is the sum of log-likelihood ratios over the matched arrivals. A
station with $r_s=0$ contributes exactly zero: no gate or neutral value is needed, and it is the principled
version of PHYS-N. Equations (9d)-(9e) of the proposal.

### Why it strengthens the thesis

It turns the two real measurements of the programme, $\hat r_s$ and $\eta$, into the two parameters of
the C-IV constraint. **Antecedent that must be cited:** NET-VISA already uses a station- and phase-specific distribution of
azimuth residuals, and a Poisson model of false arrivals
(ref. [17] of the proposal). What is added here is the estimation of the reliability weight with
catalogued events and exact chance correction, its combination with the hard test C-III and its
evaluation on the catalogue endpoint in sparse networks. The associators in common use in local
networks with machine-learning picks (GaMMA, PyOcto, REAL, GENIE) do not use polarization.

### Steps

1. Implement it in the framework of `task1e_polarization_reliability.py` as a fifth method.
2. Estimate $q_s,\kappa_s$ on the calibration realizations, never on the evaluation ones.
3. Repeat the grid of Task 1e and the cells of 1k with realistic pick errors (Task 1n);
   criterion: be no worse than the best continuous kinematic score (residual kernel with C-III) in
   any cell, and improve on it in the cells with degraded kinematic information.

### Estimated effort

One to two weeks.

---

## TASK 1i — Clutter during an active sequence

**Priority: medium-high.**

### Objective

The registered procedure required measuring two regimes, background and active sequence, and only the
background one was measured (ZAR, October to December 2024). The dangerous case of Appendix T is an aftershock
sequence, where the clutter is coda of neighbouring events and could be coherent.

### Steps

1. Search the catalogue for a sequence within 200 km of a station that passes the control:
   an event of magnitude 4 or more followed by at least twenty events within 10 km in ten days.
2. Run the same pipeline of Task 1 on that window, with the configuration of the station.
3. Same criterion on $\eta$.

### Estimated effort

One week, almost all of it in finding the sequence and downloading it.

---

## TASK 1j, 1j-b and 1k — Controls of the catalogue endpoint — EXECUTED IN SYNTHETIC DATA

Executed in this environment with `task1jk_controls.py` and `task1j2_noise.py`; rules in
`PREREGISTRO_1j_1k.txt`, written before each run. Full results in Appendix Y of
the proposal and in `results/task1jk/`.

**1j, ties.** The advantage of the prototype over the ablation in one-sided networks (0 false events against
0.30-0.43) is reproduced by any ablation that breaks ties with a continuous kinematic term, without
polarization: all give 0 false events in both seeds (the prototype, 0 and 1). Preregistered rule met: the effect is
attributed to the ties of the discrete score. In addition, a Gaussian kernel of time residuals without
polarization (TOL-K) outperforms the prototype with 3 and 4 stations.

**1j-b, timing error (exploratory, designed after seeing 1j).** With Gaussian pick error of 0.1
and 0.2 s and velocity error of 5%, polarization added to the kinematic kernel (PHYS-K) only improves
on the kernel alone in one cell out of nine: 4 stations with 0.2 s (3 against 19 false events in one seed,
$p=1.4\times10^{-4}$; 3 against 6 in the other, $p=0.45$), at a cost in recall. Rule not met.

**1k, empty windows and absolute threshold.** With 0, 1 or 2 sources per window and an absolute threshold
calibrated separately, the well-calibrated methods emit in at most one of the 20-32 empty
windows of each cell: the problem was the relative threshold, not the scores. But the best
kinematic score without polarization has the lowest or equal rate of false events per window among the
grid scores in all cells and seeds, also with 0.2 s of error; G1' is not met. Tuned GaMMA competes in the
one-sided and dispersed cells and does not reach the target recall in the sparse ones; PyOcto does not reach it
in any.

**What changes (and what changed later, in Task 1o).** In these generators, which detect all the
arrivals with exact timing, polarization adds nothing over the best kinematic score. Task 1o
showed that this was a property of the generator. G1 is reformulated as a hypothesis about real data: whether real
pick errors, the real medium and real clutter degrade the kinematic score enough for
polarization to pay. That is decided by Tasks 1n (real PhaseNet residuals) and 1g (pilot by
decimation). The reference for comparison in everything that follows is the best continuous kinematic
score, not the coincidence-count ablation.

---

## TASK 1o — Consistent generator and per-station likelihood — EXECUTED IN SYNTHETIC DATA

**Why.** On asking whether the negative of 1j-b and 1k was an implementation issue, the audit of the
generator found that it detected all true arrivals with exact timing even at signal-to-noise of
1.1 to 2.5, while it measured polarization on that same waveform (AUC 0.78). The kinematic
channel was simulated at infinite signal-to-noise.

**What was done.** `task1o_consistent.py` and `task1o_baselines.py`:
- detection with probability $\Phi((a-2)/0.5)$ and timing error of 3 times the Cramér-Rao bound;
- KIN-LR (robust kinematic likelihood) against PHYS-LR (the same with the polarization ratio
  (9d) estimated in calibration), two pick-free scores (WAVE-E, WAVE-P), and tuned GaMMA and PyOcto;
- Neyman-Pearson endpoint: per-event recall at a budget of false events per window calibrated
  separately, also read at the realized false-event rate;
- positive control with all arrivals detected and 0.5 s of error.

**Result (Appendix Y.5).** Read at the nominal budget, PHYS-LR beat KIN-LR by 0.18 to 0.25
in recall and beat GaMMA and PyOcto. Independent verification found five defects in that reading
(budget overshoot, no negative control, KIN-LR weaker than GaMMA, wider domain and
minimum of 5 picks for the baselines, generator and likelihood changed at the same time); it was repeated in
1r. With all arrivals detected and 0.5 s of error, GaMMA competes and wins in the one-sided cell.

---

## TASKS 1r to 1y — Fair repetition, physics on top of GaMMA, pick-free physics, unmodelled medium, source mechanism, coherent detector, consistent clutter — EXECUTED IN SYNTHETIC DATA

**1r (Appendix Y.6).** `task1r_fair.py`, `task1r_baselines.py`: negative control PHYS-SH (the same
polarization ratio with shuffled cosines), recall at realized false events on the exact FROC curve
(`froc_exact.py`), paired bootstrap CI, GaMMA with a seed and the same domain. With 60 windows,
indeterminate; with 200 windows per cell: PHYS-LR outperforms its control by +0.08 to +0.17 at 0.25 false events
per window in 5 of 6 cell-seeds and the best of GaMMA and PyOcto by +0.16 to +0.30 in all 6, beyond
the highest recall that GaMMA reaches at any rate. The control also outperforms KIN-LR (R-b not
met) because it was not null (E25). With double-couple sources this positive does not hold (1w).

**1s (Appendix Y.7).** `task1s_gamma_pol.py`: the polarization ratio as a filter on the events
of GaMMA separates true from false events (AUC 0.70 to 0.96) but gains at most 0.04 of recall: a filter cannot
exceed the ceiling of GaMMA, and GaMMA reaches it with at most 0.21 false events per window. The physics
has to enter the hypothesis test.

**1t (Appendix Y.8).** `task1t_wave.py`: pick-free migration of the three-component waveform.
Projected onto the predicted polarization (WAVE-P) it outperforms the same without direction (WAVE-E, analogue of
QuakeMigrate, MALMI and LOKI) by +0.13 to +0.21, its control with the direction of another station by
+0.06 to +0.18, and GaMMA and PyOcto in 5 of 6 cell-seeds. Rules written after seeing WAVE-P in 1r.

**1u (Appendix Y.9).** `task1u_hetero.py`: per-station and per-phase delays of N(0, 0.15 s) and a deviation
of the polarization of about 10 degrees that no method knows. WAVE-P keeps its advantage over
WAVE-E in 5 of 6 cell-seeds and over GaMMA in 4 of 6.

**1w (Appendix Y.10).** `task1w_mech.py`: random double-couple sources (the second verification found
that the generator radiated every S wave as pure SV, the direction that PHYS-LR and WAVE-P predict, E24) and
controls that take the geometry of another grid node (PHYS-SH was not null, E25). PHYS-LR does not
outperform its control (W-a not met) and WAVE-P falls below WAVE-E. WAVE-PT, with only the
constraints that hold for any source (P along the ray, S in the transverse plane), outperforms
its control by +0.05 to +0.19, direction-free migration by +0.07 to +0.21 and GaMMA by +0.12 to
+0.35 at 0.25 false events per window (margin over GaMMA inflated by E26; with consistent clutter it does not hold, 1y), in at least 2 of 3 cells per seed; the replication with 120 new windows
per cell meets the same rules. The one-sided cells are the weakest.

**Defects found along the way.** GaMMA initializes with the global random generator without a seed and
its catalogue depended on the state of the process (E22); the threshold sweep sampled the FROC curve with
jumps of up to 0.14 false events per window (E23). E22 was corrected before evaluating the rules; E23 was
corrected after seeing T-b fail with the registered sweep, and with that sweep R-d,
U-a and U-b also fail. Both readings are reported.

**1x (Appendix Y.11).** `task1x_mtglrt.py`: the reformulated method (design document). WAVE-MT
projects the signed three-component vectors of all stations onto the 6-dimensional
subspace of the moment tensor; MT-FW is a greedy deflation on that certificate (group matching
pursuit, not Frank-Wolfe, E29). WAVE-MT does not
uniformly outperform WAVE-PT (X-a not met). MT-FW outperforms WAVE-MT in all cells, above all because it
removes side lobes. Verification found that GaMMA had been receiving undetectable clutter picks
since 1o (E26) and that MT-FW fails with station delays larger than its window.

**1y (Appendix Y.12).** `task1y_clutter.py`: clutter with the same detection rule as the signal,
weak and strong, perturbed condition. MT-FW outperforms GaMMA by +0.11 to +0.28 (weak) and +0.17 to +0.43
(strong) at 0.25 false events per window, and its control with the geometry of another node by +0.13 to +0.41. With
weak clutter, WAVE-PT no longer outperforms GaMMA: the margins of 1t to 1w over GaMMA were inflated.

**Caveat.** All synthetic, with a matched filter close to the true wavelet, parametric
detection, no amplitudes or pick probabilities for the associators, and oracle operating points
for all methods. Tasks 1p, 1v and 1z address it, and Task 2A subjects it to a confirmatory
test.

---

## TASK 1p — Real picker on synthetic waveforms — PARTLY EXECUTED

**Priority: high.** It removes the parametric detection assumption of 1r to 1u: run PhaseNet
(SeisBench) on the synthetic waveforms, with real noise from RSNC stations added, and use
its picks, probabilities and errors as they come out. Add REAL and GENIE as baselines. Criterion: that of
1r (shuffled control, exact FROC curve, 200 windows, two seeds).
Effort: 2-3 weeks, mostly in generating traces at 100 Hz compatible with the picker.

**Status.** `gen_realistic.py` and `real_picker.py` run PhaseNet (original 2019 model from the
AI4EPS repository, not retrained) and a classical STA/LTA on the waveforms of a generator
mismatched with the detectors: gradient medium with curved rays, Ricker wavelet with variable central
frequency and phase, scattered coda. The pickers were calibrated only on seed 999. In 160
windows of that seed not used in calibration (1608 true arrivals), PhaseNet picks 4.1 % of
the arrivals, the STA/LTA 14.3 % and the parametric model degrade_y of 1o to 1y 42.0 %, with median timing
error of 0.080, 0.078 and 0.033 s and 2.08, 6.93 and 14.27 false picks per window. The parametric
model was optimistic for GaMMA in recall and in timing; with weak clutter, in everything, and with strong
clutter it delivered some 28 clutter picks per window against 7 to 8 from the STA/LTA (E27). The collapse of
PhaseNet is above all domain mismatch: the generator noise is white and occupies 0 to 8 Hz, above
the signal band, and is not the real RSNC noise that this task asks to add; its figure is
a pessimistic bound. The traces at 100 Hz with real noise, REAL and GENIE are still missing. The picks of the frozen
STA/LTA feed the associators of Task 2A. Report: `results/components/gen_realistic_report.txt`.

---

## TASK 1v — Real migration locators against MT-FW — PARTLY EXECUTED (QUAKEMIGRATE, EXPLORATORY)

**Priority: high.** 1t to 1y compared the migrations with ablations of themselves and with GaMMA, not
with the published locators. Here
QuakeMigrate (Winder et al., Seismica, 2026), MALMI (Shi et al., SRL, 2022) and the polarity-corrected
stack (Anikiev et al., GJI, 2014) are run as published, on the same synthetic windows,
with a wavelet that does not match the matched filter, scattered coda added and per-station delays
larger than the migration window. Criterion: that of 1t, with the real locators instead of
WAVE-E, double-couple sources, consistent clutter and MT-FW as the method. If MT-FW does not outperform them, its
advantage came from the matched filter or from the comparison with a pick associator. Effort: 2-3
weeks.

**Status.** `qm_wrapper.py` runs QuakeMigrate 1.2.2 from start to finish through its public API, without
modifying the package, on windows of `task1y_clutter` (perturbed condition, weak clutter) of
seed 999: calibration on 60 windows and test on another 60, in two cells (random with 6 and one-sided
with 8 stations, 47 sources in the test). The default practice of the package, P on the vertical and S on
the horizontals, fails in this generator because the rays arrive at 60-80 degrees from the vertical and the
P energy falls on the horizontals; the variant with all three components for both phases (QM-3C), a documented
option of the package, is the fair comparison. Pooled recall at ≤0.1 / ≤0.25 false events per window:
QM 0.15/0.28, QM-3C 0.53/0.53 and MT-FW 0.72/0.81; MT-FW − QM-3C = +0.19 [−0.00, +0.37] and +0.28
[+0.10, +0.43] (paired bootstrap over windows). The difference lies in the windows with two sources:
QM-3C finds 10 of 28 sources in pairs and MT-FW 22 of 28 (15 and 18 of 19 in single-source windows),
because QuakeMigrate keeps one maximum per time sample and merges the triggers within its minimum
inter-event interval. It is exploratory: one seed, 47 sources (one source is worth 0.021 of recall), the
list of QM-3C ends at 0.05 false events per window and cannot trade false events for recall, and its
hypothesis space (40,656 nodes) is larger than that of the methods of this work (675). Neither MALMI nor the polarity-corrected
stack were run, and this test uses `generate_y`, with no different wavelet or coda. The confirmatory
comparison with QuakeMigrate, in the mismatched generator and with station corrections, is Task
2A. Report: `results/components/qm_report.txt`.

---

## TASK 1z — The three mathematical pieces that the pilot did not test — EXECUTED IN SYNTHETIC DATA

**Priority: high; synthetic data and real noise already downloaded.** (1) Threshold by random field
theory: compute the expected Euler characteristic of the chi-square field with 6 degrees of freedom
over the hypothesis grid and compare the predicted false-alarm rate with the realized one in
noise-only windows, synthetic and real from the RSNC; correct by extreme values or conformally
if it does not match. (2) Certified oracle: extend `bnb_octree.py` to the dual certificate with the
energy bound in an extended window and verify that Frank-Wolfe with that oracle gives the same catalogue as with
exhaustive search. (3) Station corrections: estimate per-station delays with catalogued
events and measure the residual; the method is only valid if it stays within the ±0.25 s window.
Criterion: realized false-event rate within the interval of the predicted one; identical catalogues with and without
pruning; station residual with standard deviation smaller than 0.15 s. Effort: 3-4 weeks.

**Status.** The three pieces were measured in synthetic data with Gaussian noise and clutter; the comparison with the
real RSNC noise for piece (1) is still pending.

(1) `rft_threshold.py`, 7,200 noise-only and clutter windows in three cells. The criterion is not
met with the EEC: with curvatures estimated from roughness it overestimates the false-alarm probability by
20 to 60 times (realized rate of 0.002 where 0.25 is requested), and the variants that respect the grid
remain conservative by about 4 times, because the implemented statistic is not chi-square with 6 degrees and the
field decorrelates within one grid step. It is met with split conformal calibration and
with the GEV on some 300 noise-only windows of the same network: 24 of 24 combinations within the
predictive band on Gaussian noise. With strong clutter every threshold calibrated on Gaussian noise
fails (the conformal one gives rates of 0.60 to 0.96 where 0.01 to 0.25 is requested), and recalibrating on clutter
windows of the same law restores it. The claim of Section 7.5 that the EEC fixes the rate before
seeing data is withdrawn (E28); the threshold of the method is the empirical one, recalibrated per regime.

(2) `certified_fw.py`. The criterion is met: the certified oracle gives records identical to exhaustive
search in 72 of 72 windows (weak and strong clutter; 4, 6 and 8 stations). The gain from pruning,
on the other hand, is negligible: the energy bound avoids 4 % of the exact evaluations and, with its
bookkeeping, is 3 times slower than the dense loop, because after the first strong event the maxima
are at noise level and the bound (median of 4.6 to 5.0 times the true maximum of the cell)
exceeds them almost everywhere. The subspace envelope is 4.3 times faster with identical
records, but recomputes the statistic at almost all hypotheses: it is certified vectorization, not
pruning. A cell bound that prunes in the noise-dominated iterations has to use the subspace and
the directions of the peak vectors; it remains open.

(3) `station_corrections.py`, delays N(0, 0.4 s) fixed per network, 4 networks per cell. With 20 catalogued
events of exact location the rms residual is 0.088 s (0.101 s with 10 and 0.171 s with 5); with catalogue
errors of 1 km and 0.1 s it rises to 0.173 s, with 11.5 % of the station-phases outside ±0.25 s. The
criterion of 0.15 s is met with an exact catalogue and 10 or more events, not with the catalogue errors
tested. In recall, with 20 events MT-FW gains +0.18/+0.24 (4 stations) and +0.30/+0.30 (one-sided with
8) at ≤0.1/≤0.25 false events per window over the uncorrected version, with intervals that exclude zero,
and stays within +0.01 to +0.05 of the exact correction, with no significant difference. In the cell with 6
stations neither the estimated nor the exact correction gives a significant gain. The corrections
repair the fragility of MT-FW to large delays in 2 of 3 cells; the medium is homogeneous and the wavelet
exact, so in real data the residual will be larger, of the order of that of the run with catalogue errors.

**Theory.** Appendix AA of the proposal proves the validity of the cell bound (Proposition AA.3,
verified with `chk_aa_bound.py`: 0 violations in 480 cells, raw and after one deflation), the guarantee
of conformal calibration (Proposition AA.5(c)) and, for one source and under the idealized hypotheses it
states, that the curvature of the certificate is the efficient Fisher information with the moment tensor
as nuisance (Proposition AA.6). While writing it,
five points were found where the proposal said more than the code or the theory (E29).

Reports: `results/components/rft_report.txt`, `certified_fw_report.txt`,
`station_corrections_report.txt` and `appendix_AA_selfcheck.txt`.

---

## TASK 2A — Confirmatory benchmark with mismatched generator, real picker and real comparators — EXECUTED AND VERIFIED IN SYNTHETIC DATA

**Why.** The advantage of MT-FW in 1y was measured with the matched filter equal to the wavelet, the velocity
model exact except for random delays, parametric picks optimistic for the associators (E27) and
no published migration locator. The exploratory tests of 1p, 1v and 1z, all on
seed 999, suggest that the advantage survives, but do not confirm it. 2A is the confirmatory test.

**Rules.** Written in `PREREGISTRO_1j_1k.txt` (section Task 2A, 25 September 2026) before
any evaluation run. Seeds 2024 and 7331 had never been used with
`gen_realistic.py`, `real_picker.py` or `qm_wrapper.py`; the calibrations of the STA/LTA and of QuakeMigrate
used only seed 999.

**Design.** Generator `gen_realistic.generate_real` with its default configuration (gradient with VS0
2.8 km/s and g 0.04 1/s, curved rays; Ricker with central frequency U[1.4, 2.8] Hz with phase rotation of
up to 45 degrees; coda with energy ratio 0.3; about 10 degrees of polarization deviation;
lognormal site factor with sigma 0.5), with delays N(0, 0.3 s) per station and phase fixed per network. All
detectors use the same homogeneous model (VS = 3.0 km/s, VP = 1.73 VS), so the model
mismatch is common to all. Cells of 4 and 6 random stations and one-sided with 8; weak (0.35) and
strong (1.75) clutter. Per seed, cell and regime, 5 fixed networks, each with a catalogue of 20 events
(errors of 1 km per coordinate and 0.1 s) and 20 evaluation windows with 0, 1 or 2 sources
(probabilities 0.4, 0.4 and 0.2): 100 windows per cell. The station corrections (SC) are estimated with
`station_corrections.estimate_station_corrections` and are the same for every method that admits them:
MT-FW (added to the predicted times), GaMMA and PyOcto (subtracted from the picks) and QuakeMigrate (added to
its travel-time table). The associators receive the picks of the frozen STA/LTA; degrade_y remains as an optimistic
bound (GaMMA-DY).

**Methods.** MT-FW and MT-FW+SC (certified oracle, records identical to `fw_record`; kmax 5, min_sep 5
km, cap 6); the control MT-FW-NODE+SC (subspace of another node at 9 km or more); WAVE-PT@20; GaMMA-STALTA and
PyOcto-STALTA with and without SC (oracle configuration, which favours them); QM, QM-3C and QM-3C+SC
(recalibrated on seed 999 on `generate_real` and frozen); GaMMA-DY with and without SC.

**Endpoint.** Event recall at ≤0.1 and ≤0.25 realized false events per window on the exact FROC curve,
with the operating point chosen per method on the evaluation itself, the same for all; paired
bootstrap CI over windows, 2000 replicates, re-choosing operating point and configuration in each
replicate. Rules, per clutter regime:

| Rule | Contrast | Met if |
|:--|:--|:--|
| 2A-a, primary (state of the art) | MT-FW+SC − max(GaMMA-STALTA, GaMMA-STALTA+SC, PyOcto-STALTA, PyOcto-STALTA+SC, QM, QM-3C, QM-3C+SC), maximum per cell and cap | ≥ 0.05 at both caps in at least 2 of 3 cells, in both seeds |
| 2A-b (control of the full method) | MT-FW+SC − MT-FW-NODE+SC | ≥ 0.05 at both caps in at least 2 of 3 cells, in both seeds |
| 2A-c, secondary (optimistic picks) | MT-FW+SC − max(GaMMA-DY, GaMMA-DY+SC) | ≥ 0.05 at cap 0.25 in at least 2 of 3 cells, in both seeds |
| 2A-d, secondary (station corrections) | MT-FW+SC − MT-FW | ≥ 0.05 at cap 0.25 in at least 2 of 3 cells, in both seeds |
| 2A-e, secondary (deployable operating point) | thresholds and configurations of MT-FW+SC, GaMMA-STALTA+SC and QM-3C+SC (or QM-3C) frozen on 5 networks × 20 windows per cell and regime of seed 999, for budgets of 0.1 and 0.25 | realized rate of MT-FW+SC ≤ 2 times the budget and recall ≥ 0.05 over the best frozen comparator in at least 2 of 3 cells at both budgets, in both seeds |

**Registered interpretation.** If 2A-a is met in both regimes, the advantage of 1y survives the
mismatched generator, the real picker and QuakeMigrate. If it is met in only one, the claim is
restricted to that regime. If it fails in both, the claim of superiority over the state of the art is
withdrawn from the proposal and reported as a negative result. If 2A-b fails, the gain is not attributed to the
radiation pattern. Any analysis added after seeing seeds 2024 or 7331 is declared
exploratory.

**Status.** Executed and verified. The recalibration of QuakeMigrate and the operating points of 2A-e
were frozen on seed 999 before evaluating (`task2a_frozen_2Ae.json` was written before the first
evaluation file). Seeds 2024 and 7331 ran completely in both regimes: 100 windows and
67 to 87 true events per cell. An independent verifier reproduced the figures from the JSON files and
found no implementation error, leakage of the truth or scoring asymmetry that could fabricate a
positive. Scripts: `task2a_confirm.py` (`selfcheck`, `qmcal`, `smoke`, `run`) and `analyze_2a.py`;
results, verification and reproduction mode in `results/task2a/` (`analyze_2a_out.txt`,
`task2a_report.txt`, `README.txt`).

**Primary rule.** MT-FW+SC minus the best of the seven comparators of 2A-a, per cell, at ≤0.1 /
≤0.25 false events per window:

| Regime (seed) | 4 stations | 6 stations | one-sided with 8 |
|:--|:-:|:-:|:-:|
| weak (2024) | +0.14 / +0.26 | +0.18 / +0.25 | +0.20 / +0.27 |
| weak (7331) | +0.17 / +0.35 | +0.49 / +0.46 | +0.17 / +0.28 |
| strong (2024) | +0.19 / +0.25 | +0.26 / +0.32 | +0.16 / +0.23 |
| strong (7331) | +0.09 / +0.17 | +0.18 / +0.17 | +0.37 / +0.41 |

MT-FW+SC recovers 0.32 to 0.72 of the events at 0.1 with weak clutter and 0.13 to 0.63 with strong clutter.
The best comparator is a variant of QuakeMigrate (QM-3C or QM-3C+SC) in all 24 cases of cell,
seed and cap.

**Verdicts.**

| Rule | Weak | Strong | Detail |
|:--|:-:|:-:|:--|
| 2A-a, primary | met | met | 3 of 3 cells in each seed and regime; margins of +0.09 to +0.49 at 0.1 and of +0.17 to +0.46 at 0.25 |
| 2A-b | met | met | 3 of 3 cells except strong 7331 (2 of 3: at 4 stations +0.04 at 0.1); the +0.0513 of that cell at 0.25 is the only difference within 0.005 of the threshold in the whole analysis |
| 2A-c | met | met | weak: 2 of 3 cells per seed, because in the one-sided cell with 8 GaMMA-DY+SC ties (+0.01 in both seeds); strong: 3 of 3 |
| 2A-d | met | not met | weak: 3 of 3 and 2 of 3 (7331, 4 stations: −0.03); strong: 3 of 3 in 2024 and 1 of 3 in 7331 (−0.04, +0.04, +0.09) |
| 2A-e | depends on the reading | depends on the reading | the recall condition is met in 3 of 3 cells at both budgets, in both seeds and regimes (+0.09 to +0.46 over the best frozen comparator); the realized rate of MT-FW+SC exceeds 2 times the budget only in the one-sided cell with 8, seed 7331, budget 0.1: 0.22 (weak) and 0.23 (strong) against a limit of 0.20. With the per-cell reading it is met; with the strict reading (rate in all cells) it is not met. Both readings were programmed in `analyze_2a.py` before evaluation data existed, but the preregistered text does not decide between them; both are reported |

According to the registered interpretation, the advantage of 1y survives the mismatched generator, the real picker and
QuakeMigrate in both regimes, and the gain requires the geometry of the true hypothesis (2A-b). 2A-b
does not isolate the radiation pattern, because the control changes times and radiation at the same time: that is examined by
Task 2B, whose rule is met in both regimes, in weak only at the minimum (next section).

**Independent verification** (`task2a_report.txt`). The truth enters the detectors only at scoring;
the station corrections come only from the catalogue windows, with random streams independent
of those of the evaluation; the points of 2A-e depend only on seed 999, whose networks are different from those
of the evaluation. The hit rule is the same for all (nearest true source within 5 km
epicentral, no timing tolerance), and GaMMA and PyOcto are scored without cap or separation, that is,
more leniently. The sign of the corrections was checked with the true delays as an oracle in
80 regenerated windows: MT-FW goes from 0.25/0.36 without correction to 0.51/0.54, QM-3C from 0.17/0.24 to
0.36/0.36 and GaMMA-DY from 0.29/0.36 to 0.41/0.41, and the oracle correction never worsens. The near-zero recall
of GaMMA and PyOcto with STA/LTA picks is starvation by the picker, not wiring: the frozen STA/LTA picks
between 9.6 and 16.8 % of the arrivals at 0.3 s, only 20 % of the events receive the 3 correct picks
that GaMMA needs, and with true-time picks and the same wiring GaMMA recovers 0.88 and PyOcto 0.69
at 0.1. An independent rederivation of the FROC curve reproduces Table 1 of `analyze_2a_out.txt`.

**Declared deviations.** The catalogue was drawn per network and per regime, not only per network: with strong
clutter the catalogue windows become contaminated and the rms residual of the corrections rises to 0.20-0.28 s,
against 0.12-0.21 s with weak clutter and 0.26-0.34 s for the delays themselves. This is the cause of the failure of 2A-d
with strong clutter; in seed 7331, with corrections from a weak-clutter catalogue or with the
true delays the gain returns (exploratory). The same corrections are given to every method
with SC, so the deviation favours none. Four files of seed
999 came from an earlier version of `task2a_confirm.py`; repeated with the final version, their effect
on 2A-e is null or marginal. `gen_realistic.py` was edited after the preregistration to add the
delay argument that the design requires, checked by the self-test.

**Limitations.** The effective rival is QuakeMigrate; against GaMMA and PyOcto the claim holds for
those associators with an operational STA/LTA picker. The intervals over windows ignore that the windows
share a network (5 per cell): with bootstrap over networks (exploratory) the 12 cells stay above 0 at cap 0.25
(minimum lower bound +0.03), but at cap 0.1 four include 0 (weak 2024 with 6 stations
[−0.03, +0.34], weak 7331 one-sided [−0.18, +0.41], strong 2024 one-sided [−0.03, +0.36] and strong
7331 with 6 stations [−0.05, +0.29]). About half of the events are in windows with two sources,
which the QuakeMigrate trigger does not separate when their origin times are close; in windows with 0 or 1
source MT-FW+SC still outperforms QM-3C, QM-3C+SC and GaMMA-STALTA+SC in the 12 cells, with minimum
margins of +0.06 to +0.07 at 0.1, but GaMMA-DY+SC, with optimistic picks, outperforms MT-FW+SC at 0.1 in 5
of 6 weak-clutter cells. The advantage over an associator with good picks is established only with
strong clutter. The sources are synthetic point double couples, within the subspace of MT-FW: the
mismatch lies in the medium, the wavelet, the coda, the polarization and the site, not in the source model.
QuakeMigrate was recalibrated only with weak clutter and in two cells; the associators, on the other hand, choose their
configuration with an oracle on the evaluation. The corrections worsened MT-FW in weak 7331 with 4
stations (−0.14 at 0.1). White noise, a single generator and two seeds.

**What remains.** The comparators closest to the method were not run in 2A: the block-sparse
inversion of origin time, hypocentre and moment tensor of Vera Rodriguez, Sacchi and Gu (2012), the
mechanism-corrected stack (Anikiev et al., 2014; Chambers et al., 2014; Zhebel and Eisner, 2015) and
MALMI (Task 1v). The first two, in our own implementation, were run by Task 2C (below); MALMI
is still pending. Nor were the real RSNC noise added to the synthetics (Tasks 1p and 1z) or the real
data (Tasks 1n and 1g).

---

## TASK 2B — Ablation of the radiation pattern within the same solver — EXECUTED AND VERIFIED IN SYNTHETIC DATA

**Why.** The novelty audit pointed out that the advantage of a pick-free method may be due to not
using picks and not to the radiation pattern. The control of 2A-b takes the subspace of another node, so it
changes times and radiation at the same time: it controls the geometry, but does not isolate the sign.

**Rules.** Written in `PREREGISTRO_1j_1k.txt` (section Task 2B, 26 September 2026) before
computing it, and not modified. The time "14:30 UTC" of the entry was written without consulting a clock and is
wrong (E34): the text was added immediately before the implementation, and the first file of 2B,
`task2b_ablation.py`, was created at 14:16:49 UTC; the erratum is in the preregistration itself. There is no
hash from that moment, so the ordering rests only on the file-system timestamps.

**Design.** The same windows, networks, catalogues and corrections as 2A (seeds 2024 and 7331, both
regimes), regenerated deterministically with the functions of `task2a_confirm.py` itself,
imported and not copied: 20 windows per network, 100 per cell, 1200 in total. The windows consume only
their own random streams, which the picker, `degrade_y` and the node map do not touch, so omitting
those steps reproduces the exact sequence. MT-FW+SC is the saved record of 2A. MT-FW-FREE+SC is the
same call to the certified solver (675 nodes, kmax 5, exclusion of 5 km, cap 6 per station and phase,
the same corrections) with a single change: at each node the radiation operator G (6 ns × 6) and its basis
are replaced by the identity of dimension 6 ns, free three-component amplitude per station and phase.
The least-squares amplitude is then the observed peak vector itself, the statistic is
half of its total energy, ‖v‖²/2, invariant to the sign of each station, and the deflation subtracts the
observed vectors convolved with the autocorrelation of the matched filter. Endpoint, exact FROC
curve, operating point chosen on the evaluation and paired bootstrap over windows (the same 2000
resamples of 2A) do not change; the bootstrap over networks (2000 replicates, 5 networks with replacement) is the
additional reading. The pairing of windows between the JSON files of 2A and 2B is asserted in the analysis (truth,
origin times and records of MT-FW+SC identical). Computation: 60 network jobs, 495 s of CPU in total.

**Rule 2B, per regime.** MT-FW+SC − MT-FW-FREE+SC ≥ 0.05 at both caps in at least 2 of 3 cells, in
both seeds, on point estimates. If it fails, the gain of 2A is attributed to pick-free migration
with deflation and not to the radiation pattern, and the proposal declares it.

**Status.** Executed and verified with exact controls. Scripts: `task2b_ablation.py` (runs),
`analyze_2b.py` (rule and tables) and `task2b_kernel_check.py` (numerical sensitivity, exploratory);
results, controls and reproduction mode in `results/task2b/` (`analyze_2b_out.txt`,
`task2b_report.txt`, `README.txt`). No file of 2A was modified.

**Result.** MT-FW+SC minus MT-FW-FREE+SC, per cell, at ≤0.1 / ≤0.25 false events per window (point
estimate; the intervals over windows and over networks are in `analyze_2b_out.txt`, Table 2):

| Regime (seed) | 4 stations | 6 stations | one-sided with 8 | cells that meet it |
|:--|:-:|:-:|:-:|:-:|
| weak (2024) | +0.17 / +0.28 | +0.18 / +0.11 | −0.01 / −0.01 | 2 of 3 |
| weak (7331) | +0.15 / +0.24 | +0.26 / +0.12 | +0.04 / +0.0494 | 2 of 3 |
| strong (2024) | +0.25 / +0.26 | +0.32 / +0.35 | +0.17 / +0.22 | 3 of 3 |
| strong (7331) | +0.06 / +0.12 | +0.26 / +0.17 | +0.33 / +0.33 | 3 of 3 |

**Verdict.** Rule 2B is met in both regimes. With strong clutter it is met in all 3 cells of
each seed, with differences of +0.06 to +0.35, and 11 of the 12 intervals over windows exclude 0.
With weak clutter it is met at the minimum, 2 of 3 cells per seed: the operator contributes +0.15 to +0.28 with
4 stations and +0.11 to +0.26 with 6 stations, where 6 of the 8 intervals include 0, and in the
one-sided cell with 8 MT-FW-FREE+SC equals MT-FW+SC in both seeds (the +0.0494 of 7331 at 0.25 lies
within 0.005 of the threshold, and the cell fails anyway at 0.1). Re-applied in each replicate of the
bootstrap (exploratory, independent resamples per cell), the rule is met in 0.98 (windows) and
0.97 (networks) of the replicates with strong clutter, and only in 0.59 and 0.62 with weak clutter: in weak clutter the
verdict is fragile.

**Controls.** The self-test (`task2b_selfcheck.txt`; seed 2024, weak, 6 stations, network 0,
window 0) checks that the score of FREE is ‖v‖²/2 at the 675 × 31 hypotheses (relative error
6.8e-16), that it does not change bit for bit when the sign of the vector of one station is inverted, that the
least-squares amplitude exactly reproduces the observed vector, and that inverting the sign of the waveform
of any of the 6 stations leaves the records of FREE identical (6 of 6) and changes those of
MT-FW+SC (6 of 6); in that window MT-FW+SC misses the event and FREE detects it, so the ablation is not
trivially weaker window by window. In the 60 networks, stations, delays, draws and locations
of the catalogue and corrections are bit-for-bit identical to those of 2A (60 of 60); truth and origin times
match in 1200 of 1200 windows; MT-FW+SC recomputed in windows 0 and 1 of each network matches
EXACTLY what was saved (120 of 120), and certified FREE matches dense FREE exactly (120 of
120).

**Impact of the BLAS kernel (E31).** The first attempt stopped at the exact control, as the
design requires, after passing 32 networks: the 6-station networks of seed 7331 written in the resumed
sessions of 2A gave MT-FW+SC scores differing in the last 1 or 2 digits (relative ≤ 6e-15,
same hits and same order) and, in one network, a correction differing by 8.9e-16 s. The cause is the
OpenBLAS of numpy, which chooses the compute kernel at run time: on this machine it chooses Prescott by
default, but the 18 networks that 2A wrote in its resumed sessions of 7331 (4 stations in both
regimes; 6 stations, the 5 networks in weak and networks 2 to 4 in strong) were computed with the
AVX-512 SkylakeX kernel. With `OPENBLAS_CORETYPE=SkylakeX` those 18 networks are reproduced bit for bit, and with Prescott,
Haswell, Sandybridge or Zen they are not; the other 42 are reproduced bit for bit with the default kernel. The final
run fixed the kernel per network (SkylakeX for those 18, the default for the rest) and recorded it in each
part file it wrote (the 32 parts of the first attempt do not carry the field, but ran with the
default kernel); all networks passed the exact control without changing it, and no threshold
or tolerance was relaxed. Sensitivity (exploratory, `task2b_kernel_check.txt`): with the 18 networks recomputed with the
default kernel, the records of MT-FW+SC differ by at most 6e-15 relative with no change of hits
or order; those of FREE diverge in later Frank-Wolfe steps in some windows (up to 17 % in
score, because of nearly tied candidates that swap by rounding) and change hits or order in
2 of 360 windows; the recall at 0.1 and at 0.25 is identical for both methods in the 4 affected cells. For
2A, its per-network results mix two kernels, with an effect on MT-FW only at the level of rounding.

**Interpretation.** According to the registered interpretation, the gain of 2A is not attributed solely to
pick-free migration with deflation: within the same solver the radiation operator contributes ≥ 0.05 at both
caps in most cells. With strong clutter the operator sustains the result of 2A: without it,
MT-FW-FREE+SC stays at the level of the best comparator of 2A (QuakeMigrate), between −0.08 and +0.08, and below
WAVE-PT@20 in all cells and caps (−0.01 to −0.23), even though WAVE-PT@20 has neither deflation nor SC.
With weak clutter the picture is mixed (exploratory): MT-FW-FREE+SC outperforms the best comparator at 0.25 in 5
of 6 cell-seeds, by +0.10 to +0.34 (all except 4 stations 2024, −0.01), so that a good part of
the advantage of 2A over QuakeMigrate with weak clutter comes from pick-free migration with deflation and SC and
does not need the radiation pattern; the operator weighs above all in the sparsest networks. The replacement by
the identity removes at the same time the sign of the double-couple pattern, the polarization constraint (P along
the ray, S in the transverse plane) and the amplitude geometry (1/d and P/S ratio): 2B
establishes that the operator matters, it does not isolate the sign, which was the question of the audit. The comparison
with WAVE-PT@20, without sign but with polarization, suggests that polarization is a good part of the value of the
operator with strong clutter, whose transients have random polarization. Separating sign and
polarization would require another non-preregistered ablation (for example, free amplitude along the ray
for P and in the transverse plane for S, "FREE-PT"), which should be declared exploratory.

**Limitations.** Five networks per cell: the intervals over windows ignore the clustering and the
bootstrap over networks with 5 clusters is coarse; several cells have intervals that include 0 even though
the point estimate passes. The rule is decided on point estimates and in weak clutter passes with the minimum
number of cells. The operating point is chosen on the evaluation itself for both methods, as in
2A, which is symmetric. FREE inherits the choices of the solver (cap 6, exclusion of 5 km, kmax 5,
deflation with the autocorrelation of the matched filter), tuned with the moment-tensor statistic;
a retuned energy detector could perform somewhat better, and the registered ablation keeps the solver
fixed on purpose.

---

## TASK 2C — Closest comparators, the convex method on the grid and replication by network — EXECUTED AND VERIFIED IN SYNTHETIC DATA

**Why.** The independent audit (56/100, major revision) pointed out three gaps that 2A and 2B did not
close. The two comparators closest to the method were missing: the block-sparse inversion of
Vera Rodriguez, Sacchi and Gu (2012) and the mechanism-corrected stack (Anikiev et al., 2014). The convex
method that the thesis defends, group BLASSO, was not implemented, and all the evidence came from
MT-FW, its greedy approximation. And the inference did not respect the network as the unit of replication: 5 networks per
cell and intervals over windows.

**Rules.** Written in `PREREGISTRO_1j_1k.txt` (section Task 2C, 26 September 2026, 15:05 UTC)
before implementing it and before seeing any of its results. The only addition is the budget note dated
at 19:56 UTC, before any run of seeds 2024 and 7331: the smoke test on seed
999 estimated 0.79 h on 2 cores, below the threshold of 8 h, so the budget rule was not triggered
and part (ii) ran with the 25 new networks planned.

**Design.** Part (i) uses the same 1200 windows as 2A (seeds 2024 and 7331, both regimes, 3
cells, networks 0 to 4, 20 windows per network), regenerated bit for bit with the functions of `task2a_confirm.py`;
MT-FW+SC and QM-3C+SC are the saved records of 2A. Part (ii) adds 25 new networks per cell,
seed and regime (indices 5 to 29, with the same seed construction as 2A) and 12 windows per network, with
MT-FW+SC, QM-3C+SC, BOMP+SC and MTSTACK+SC; together with the 5 networks of 2A they add up to 30 networks and 400 windows per
cell. GaMMA and PyOcto do not enter part (ii) because in 2A they were left without enough picks. All
methods use the same station corrections as 2A, the same grid (675 nodes at 3 km and 31 origin
times at 0.1 s), the exclusion of 5 km epicentral, at most 5 events per window and the same hit
rule.

**New methods.** The three work in the matched-filter domain and read each station and phase at the
exact predicted sample, without peak search or cap, as the preregistration fixed. A block is a pair
(node, origin time) with 6 columns, one per component of the moment tensor; it is the image, through the
matched filter, of a dictionary of Green's functions in the data domain, and its Gram matrix
includes the autocorrelation of the filter and the overlap of P and S at the same station.

| Method | Form | Catalogue and score |
|:--|:--|:--|
| BOMP+SC | block orthogonal matching pursuit in the form of Vera Rodriguez et al. (2012): it chooses the block with the highest energy of the projection of the residual, e = hᵀK_bb⁺h, and refits by least squares all the active blocks | score e/2 on entry, with the cumulative minimum of MT-FW; exclusion of 5 km |
| MTSTACK+SC | mechanism-corrected stack in the form of Anikiev et al. (2014): at each node and time, least-squares tensor with the amplitudes at the exact times | value ‖Gm‖²/2, maximum over time, non-maximum suppression at 20 km, no deflation |
| GL-FW+SC | the convex problem of the thesis restricted to the grid: group BLASSO ½‖y − Au‖² + λ Σ_b ‖A_b u_b‖ over the dictionary of BOMP, by fully corrective Frank-Wolfe with duality gap | path of 30 geometric values of λ, from λ_max to λ_max/30; score of a block = the largest λ at which it enters the support |

With the group norm whitened per block, the dual certificate of a block is exactly the selection statistic
of BOMP, and λ_max²/2 is the first score of BOMP; MTSTACK is the projection energy of
WAVE-MT at the exact times. The three methods share dictionary and grid with MT-FW, and BOMP and GL-FW
additionally share the statistic: they differ in the algorithm. BOMP+SC and MTSTACK+SC are our own implementations
based on the published descriptions, not the code of their authors.

**Endpoint and statistics.** As in 2A: recall at ≤0.1 and ≤0.25 realized false events per window on the exact FROC
curve, with the operating point chosen on the evaluation itself, and rules decided on
point estimates with threshold 0.05. For 2C-a and 2C-b, paired bootstrap over windows (the 2000
resamples of 2A) and over networks (5 networks with replacement). For 2C-c, hierarchical bootstrap of 2000 replicates: 30
networks with replacement and, within each chosen network, its windows with replacement, re-choosing the operating
point in each replicate. Holm over the 36 comparisons of each regime (3 comparators × 2 caps × 3
cells × 2 seeds), with one-sided bootstrap p at level 0.025.

| Rule | Contrast | Met if |
|:--|:--|:--|
| 2C-a (close comparators) | MT-FW+SC − max(BOMP+SC, MTSTACK+SC), part (i) | ≥ 0.05 at cap 0.25 in at least 2 of 3 cells, in both seeds. If it fails, the advantage of MT-FW is attributed to the moment tensor with deflation, already present in Vera Rodriguez et al. (2012), and not to its own elements (peak search and cap) |
| 2C-b (convex method, descriptive) | GL-FW+SC − MT-FW+SC at both caps, part (i) | no threshold: equals if \|diff\| < 0.05, otherwise outperforms or falls below |
| 2C-c (replication by network, primary of part ii) | MT-FW+SC − QM-3C+SC at cap 0.25, 30 networks per cell | lower bound of the hierarchical 95 CI > 0 and point estimate ≥ 0.05 in at least 2 of 3 cells, in both seeds; the same is reported against BOMP+SC and MTSTACK+SC and at cap 0.1, with Holm |

**Status.** Executed and verified with exact controls and with validations of the methods on seed
999. Scripts: `task2c_methods.py` (methods), `task2c_run.py` (runs), `task2c_validate.py`
(validation), `analyze_2c.py` (rules and tables) and `task2c_explore.py` (exploratory diagnostics);
results, controls and reproduction mode in `results/task2c/` (`analyze_2c_out.txt`,
`task2c_report.txt`, `task2c_validation.txt`, `README.txt`). Computation: 318 s of CPU for part (i) and 5342 s
for part (ii), 1.57 h in total, 47 minutes of wall time in 2 processes.

**Review of the inherited code (E32).** `task2c_methods.py` and `task2c_run.py` came from an earlier
session. Before any run of seeds 2024 and 7331 the dictionary, the
Gram matrix with the autocorrelation, the selection, the refit and the score of BOMP, MTSTACK, the certificate of the
group lasso and the gap formula were derived again, and they matched the preregistered definitions. But the code
could not be run as written. The validation that the runner imported did not exist. The
check of the OpenBLAS kernel only compared the environment variable with "SkylakeX", and since the default
kernel of this machine changed from Prescott, the one of 2A and 2B, to SkylakeX, with the variable unset the 42 networks
of 2A written with Prescott would have run with the wrong kernel. And the restricted solver of GL-FW
required a relative KKT violation of 1e-11 that floating point does not reach: the Armijo test failed
at the rounding level and the solver looped for up to 200 outer iterations, 39 s for 6 values of λ in
one window. `task2c_validate.py` was written. `OPENBLAS_CORETYPE` is set per process (SkylakeX in the 18
networks of E31, Prescott in the others, in all of part (ii) and on seed 999), each network checks with
threadpoolctl the kernel actually loaded and saves it in its file, and each process only launches the networks
of its kernel. The restricted tolerance became 1e-9, a thousand times stricter than the admission tolerance of
Frank-Wolfe (1e-6); the Armijo test received a rounding slack of 1e-14 |F|; the dense whitening
O(n³) was replaced by a block product O(n²), and the error for "active block with certificate larger
than λ" by a bounded counter of re-solves that was never triggered (0 in 1200 windows).

**Validation of the methods** (`task2c_validation.txt`, seed 999 only). Against explicit atoms in the
data domain, the relative error of Aᵀy is ≤ 8.7e-16 and that of AᵀA ≤ 7.3e-16, and every diagonal block
has rank 6. With a noise-free source within the model, BOMP recovers the block in 12 of 12 cases and the
tensor with relative error ≤ 6.5e-13, and the first block of GL-FW is the true one in 12 of 12. With the physical
waveform (arrivals in continuous time, source at a node and origin on the grid) node and origin are
recovered in 12 of 12, but the tensor only approximately (cosine ≥ 0.978, 0.982 and 0.874 with 4, 6 and
one-sided with 8), because reading at the predicted sample a pulse that arrives up to one sample out of place
distorts the amplitudes; with the arrivals at the predicted sample the recovery is exact (1 − cosine ≤
2.2e-16). GL-FW on 24 real windows, solving 2 values of λ beyond the stop: maximum certificate
outside the support / λ = 1.000000000 (tolerance 1 + 1e-6), deviation on the support ≤ 5.1e-10, relative
gap ≤ 1.5e-9, record with early stopping equal to that of the longest path in 24 of 24, and
λ_max²/2 equal to the first score of BOMP to 2e-16. MTSTACK matches the projection energy ‖Qᵀv‖²/2
of WAVE-MT at the exact times (relative ≤ 1.5e-13). In the 1200 windows of part (i) GL-FW always ended
with a complete record, with maximum relative gap 1.4e-8 (median 3e-14 at the last λ solved),
certificate outside the support ≤ 1.000000000, deviation on the support ≤ 8.7e-10 and 0 re-solves.

**Controls.** In part (i), stations, delays, draws and locations of the catalogue and corrections are
bit-for-bit identical to 2A (60 of 60 networks); truth, origin times and number of sources match in 1200
of 1200 windows, and MT-FW+SC recomputed in windows 0 and 1 of each network is EXACTLY equal to what was
saved (120 of 120). In part (ii), certified MT-FW+SC equals the dense one in 600 of 600 windows, and
in the smoke test the pipeline of part (ii) reproduces the saved records of MT-FW+SC and QM-3C+SC
of 2A (6 of 6 cells). No control failed.

**Results of part (i).** Recall at ≤0.1 / ≤0.25, 100 windows per cell:

| Regime, cell (seed) | MT-FW+SC | BOMP+SC | MTSTACK+SC | GL-FW+SC | QM-3C+SC |
|:--|:-:|:-:|:-:|:-:|:-:|
| weak, 4 stations (2024) | 0.322 / 0.517 | 0.161 / 0.230 | 0.092 / 0.161 | 0.057 / 0.207 | 0.184 / 0.253 |
| weak, 6 stations (2024) | 0.458 / 0.556 | 0.083 / 0.222 | 0.069 / 0.139 | 0.083 / 0.125 | 0.278 / 0.306 |
| weak, one-sided with 8 (2024) | 0.549 / 0.648 | 0.197 / 0.380 | 0.070 / 0.225 | 0.070 / 0.197 | 0.352 / 0.380 |
| weak, 4 stations (7331) | 0.321 / 0.513 | 0.154 / 0.282 | 0.077 / 0.154 | 0.077 / 0.218 | 0.154 / 0.167 |
| weak, 6 stations (7331) | 0.721 / 0.794 | 0.176 / 0.353 | 0.118 / 0.265 | 0.118 / 0.235 | 0.235 / 0.338 |
| weak, one-sided with 8 (7331) | 0.642 / 0.765 | 0.160 / 0.247 | 0.123 / 0.173 | 0.086 / 0.160 | 0.469 / 0.481 |
| strong, 4 stations (2024) | 0.312 / 0.390 | 0.195 / 0.260 | 0.078 / 0.182 | 0.039 / 0.182 | 0.117 / 0.130 |
| strong, 6 stations (2024) | 0.378 / 0.486 | 0.095 / 0.162 | 0.081 / 0.135 | 0.095 / 0.135 | 0.122 / 0.162 |
| strong, one-sided with 8 (2024) | 0.313 / 0.422 | 0.108 / 0.217 | 0.096 / 0.157 | 0.084 / 0.157 | 0.157 / 0.193 |
| strong, 4 stations (7331) | 0.128 / 0.218 | 0.103 / 0.141 | 0.103 / 0.115 | 0.115 / 0.154 | 0.038 / 0.051 |
| strong, 6 stations (7331) | 0.357 / 0.417 | 0.071 / 0.214 | 0.036 / 0.119 | 0.024 / 0.107 | 0.179 / 0.250 |
| strong, one-sided with 8 (7331) | 0.627 / 0.707 | 0.067 / 0.253 | 0.080 / 0.173 | 0.080 / 0.133 | 0.253 / 0.293 |

Rule 2C-a: MT-FW+SC minus the best of BOMP+SC and MTSTACK+SC at cap 0.25 (the best is BOMP+SC in all 12
cells), point estimate with the intervals over windows (W) and over networks (N):

| Regime (seed) | 4 stations | 6 stations | one-sided with 8 | cells that meet it |
|:--|:-:|:-:|:-:|:-:|
| weak (2024) | +0.287 W[+0.119, +0.430] N[+0.159, +0.453] | +0.333 W[+0.145, +0.481] N[+0.000, +0.500] | +0.268 W[+0.111, +0.456] N[+0.136, +0.452] | 3 of 3 |
| weak (7331) | +0.231 W[+0.022, +0.422] N[+0.037, +0.356] | +0.441 W[+0.194, +0.628] N[+0.138, +0.613] | +0.519 W[+0.355, +0.674] N[+0.355, +0.711] | 3 of 3 |
| strong (2024) | +0.130 W[+0.000, +0.260] N[+0.069, +0.253] | +0.324 W[+0.167, +0.471] N[+0.182, +0.507] | +0.205 W[+0.051, +0.365] N[+0.081, +0.366] | 3 of 3 |
| strong (7331) | +0.077 W[−0.038, +0.194] N[−0.051, +0.198] | +0.202 W[+0.053, +0.381] N[+0.043, +0.345] | +0.453 W[+0.279, +0.600] N[+0.270, +0.589] | 3 of 3 |

Rule 2C-b: GL-FW+SC minus MT-FW+SC, at ≤0.1 / ≤0.25 (point estimate; the intervals over windows and
over networks exclude 0 except in strong 7331 with 4 stations):

| Regime (seed) | 4 stations | 6 stations | one-sided with 8 |
|:--|:-:|:-:|:-:|
| weak (2024) | −0.264 / −0.310 | −0.375 / −0.431 | −0.479 / −0.451 |
| weak (7331) | −0.244 / −0.295 | −0.603 / −0.559 | −0.556 / −0.605 |
| strong (2024) | −0.273 / −0.208 | −0.284 / −0.351 | −0.229 / −0.265 |
| strong (7331) | −0.013 / −0.064 | −0.333 / −0.310 | −0.547 / −0.573 |

GL-FW+SC minus BOMP+SC ranges from −0.18 to +0.01.

**Results of part (ii).** Rule 2C-c, 30 networks and 400 windows per cell; recall at ≤0.25 and
MT-FW+SC minus QM-3C+SC with its hierarchical 95 CI:

| Regime, cell (seed) | MT-FW+SC | QM-3C+SC | BOMP+SC | MTSTACK+SC | MT-FW+SC − QM-3C+SC |
|:--|:-:|:-:|:-:|:-:|:-:|
| weak, 4 stations (2024) | 0.562 | 0.214 | 0.253 | 0.146 | +0.347 [+0.248, +0.453] |
| weak, 6 stations (2024) | 0.623 | 0.300 | 0.275 | 0.144 | +0.323 [+0.228, +0.416] |
| weak, one-sided with 8 (2024) | 0.656 | 0.373 | 0.290 | 0.194 | +0.283 [+0.189, +0.376] |
| weak, 4 stations (7331) | 0.466 | 0.160 | 0.270 | 0.166 | +0.306 [+0.194, +0.401] |
| weak, 6 stations (7331) | 0.685 | 0.338 | 0.299 | 0.208 | +0.347 [+0.242, +0.439] |
| weak, one-sided with 8 (7331) | 0.638 | 0.414 | 0.304 | 0.190 | +0.224 [+0.119, +0.328] |
| strong, 4 stations (2024) | 0.249 | 0.074 | 0.169 | 0.117 | +0.175 [+0.104, +0.254] |
| strong, 6 stations (2024) | 0.431 | 0.184 | 0.197 | 0.138 | +0.247 [+0.150, +0.337] |
| strong, one-sided with 8 (2024) | 0.489 | 0.252 | 0.252 | 0.203 | +0.236 [+0.135, +0.334] |
| strong, 4 stations (7331) | 0.208 | 0.109 | 0.150 | 0.134 | +0.099 [+0.029, +0.181] |
| strong, 6 stations (7331) | 0.459 | 0.228 | 0.191 | 0.119 | +0.231 [+0.125, +0.339] |
| strong, one-sided with 8 (7331) | 0.556 | 0.261 | 0.246 | 0.164 | +0.295 [+0.212, +0.383] |

With Holm, in weak clutter all 36 comparisons reject in favour of MT-FW+SC. In strong clutter 28 of 36 reject; the
eight that do not are all with 4 stations: in 2024 against BOMP+SC at both caps (+0.065 [+0.003, +0.141] and
+0.080 [+0.006, +0.165]), and in 7331 against QM-3C+SC at both caps (at 0.25, +0.099 with Holm p 0.044),
against BOMP+SC at both caps (intervals that include 0) and against MTSTACK+SC at both caps. In the 25
new networks alone (exploratory), MT-FW+SC − QM-3C+SC is +0.17 to +0.32 at 0.1 and +0.20 to +0.39 at 0.25 in
weak clutter, and +0.07 to +0.23 and +0.08 to +0.26 in strong clutter: the same magnitudes as in the networks of 2A.

**Verdicts.**

| Rule | Weak | Strong | Detail |
|:--|:-:|:-:|:--|
| 2C-a | met | met | 3 of 3 cells in each seed; +0.231 to +0.519 in weak and +0.077 to +0.453 in strong; re-applied in each replicate over windows, it is met in 100 % (weak) and 99 % (strong) of the replicates |
| 2C-b | below | below | GL-FW+SC falls below MT-FW+SC in 23 of 24 cases of cell, seed and cap, by −0.064 to −0.605; the remaining case, strong 7331 with 4 stations at 0.1 (−0.013), counts as "equals"; it never outperforms it |
| 2C-c | met | met | literal reading: 3 of 3 cells in each seed, with lower bounds of +0.119 to +0.248 in weak and of +0.029 to +0.212 in strong; with Holm: 3 of 3 in weak and 3 of 3 and 2 of 3 in strong (7331, 4 stations: +0.099 [+0.029, +0.181], Holm p 0.044), and the rule is met all the same |

No case came within 0.005 of the threshold of 0.05.

**Declared deviations.** The path of GL-FW stops at the first λ of the grid that completes
the record of 5 events: later blocks enter at a strictly smaller λ and are ordered after
the 5 accepted ones, so the record is identical to that of the full path (argument plus the
check on 24 windows). That is why the duality gap is reported at each λ solved (2 to 14 per
window, median 3) and not at the 30 of the preregistration; the full path, with supports of hundreds of
blocks and more than one minute per λ at λ/λ_max ≈ 0.5, was too expensive. The preregistration did not fix the
weighting of the group norm; the block-whitened one, ‖A_b u_b‖, was chosen, with which the certificate
is the statistic of BOMP (and that of WAVE-MT at the exact times). The score of BOMP is the energy of the
projection of the residual onto the block on entry, which is our reading of "energy explained by the
block on entry", not the drop in residual energy after the joint refit. MTSTACK uses the generalized
inverse of G, ignores the rare P/S temporal overlap that the Gram matrix of BOMP includes and has no cap, as
the preregistration describes. The tolerances are 1e-6 for the admission of Frank-Wolfe and 1e-9 for the restricted
problem. The OpenBLAS kernel is set explicitly; with Prescott, the OpenBLAS bundled with scipy chooses
Katmai where in 2A it chose SkylakeX, but scipy.linalg is not used by MT-FW, nor by the corrections nor by the
methods of 2C, and all exact controls passed. The family and level of Holm are our reading
(36 comparisons per regime, one-sided 0.025, p resolution of 1/2001); 2C-c is reported with the literal
reading and with that of Holm, and they agree on the verdict. The point estimate of 2C-c pools all the
windows, so each network of 2A weighs 20/12 relative to a new one, while the hierarchical bootstrap weighs
networks equally. BOMP+SC and MTSTACK+SC were not rerun on the networks of 2A for part (ii): the
records of part (i) on the same windows are used.

**Exploratory diagnostics** (after seeing part (i), not preregistered; `task2c_explore.py`, seed
999, 60 windows per cell). The exact-time methods not only fell below MT-FW+SC but, in
many cells, also below QuakeMigrate, which invites suspicion of the implementation. The cause was sought with the
ideal homogeneous generator of 1x, where the model of the detectors is exact and only the sources and the
origins lie off the grid (E1), and with the same generator moving each source to its nearest
node and rounding its origin to the 0.1 s grid, so that only the rounding to the sample remains (E3).
Placing the truth on the nodes, which E7 teaches to avoid in an evaluation, is deliberate here: it isolates
the effect of the grid and is not used for any rule.

| Condition, cell | MT-FW | BOMP | MTSTACK | GL-FW |
|:--|:-:|:-:|:-:|:-:|
| E1 ideal, 4 stations | 0.627 / 0.712 | 0.220 / 0.458 | 0.119 / 0.169 | 0.136 / 0.237 |
| E1 ideal, 6 stations | 0.864 / 0.886 | 0.182 / 0.432 | 0.091 / 0.182 | 0.114 / 0.273 |
| E1 ideal, one-sided with 8 | 0.810 / 0.810 | 0.500 / 0.672 | 0.483 / 0.603 | 0.414 / 0.500 |
| E3 sources on nodes, 4 stations | 0.678 / 0.763 | 0.898 / 0.966 | 0.203 / 0.542 | 0.424 / 0.729 |
| E3 sources on nodes, 6 stations | 0.909 / 0.932 | 0.955 / 0.977 | 0.659 / 0.818 | 0.750 / 0.864 |
| E3 sources on nodes, one-sided with 8 | 0.793 / 0.845 | 0.914 / 0.948 | 0.672 / 0.759 | 0.741 / 0.914 |

With the sources on the nodes BOMP outperforms MT-FW and GL-FW comes close to it, so the implementations are
sound. The collapse comes from reading the matched filter at the exact predicted sample: with nodes at 3 km and
origins at 0.1 s, an off-node source produces moveout errors that differ per station and are an
appreciable fraction of the main lobe of the autocorrelation at 2 Hz, and the residual errors of the
station corrections (rms 0.12 to 0.28 s) add to them. The amplitudes read there come out attenuated or with the
sign flipped, and that destroys the moment-tensor fit. With oracle corrections in the windows of
2A (2024, weak; E2) BOMP improves little (from 0.230, 0.222 and 0.380 to 0.264, 0.264 and 0.465 at 0.25), so the
dominant term is the discretization of the grid and not the error of the corrections. The peak search
of ±0.25 s of MT-FW absorbs both. With the perturbations of 1x all exact-time methods fall
further. The score ties of GL-FW (the entry λ is coarse: 1694 and 1654 of 3000 events share
a score with another in their window) do not explain why it falls below BOMP: breaking them by the group
norm on entry changes the recall by at most 0.015. The most probable cause is the shrinkage of the lasso:
a strong event is under-fitted, and the energy it leaves lets spurious blocks more than 5 km away
enter before weaker true events; the orthogonal refit of BOMP removes that energy.

**Interpretation.** According to the registered interpretation, since 2C-a is met in both regimes, the advantage
of MT-FW over the block-sparse inversion and over the mechanism-corrected stack is not attributed to the
moment tensor with deflation, which is already in Vera Rodriguez et al. (2012), but to the own elements
of MT-FW; the exploratory diagnostics indicate that the one that matters is the peak search, and that the cap is not
needed with the sources on the nodes. What 2C-a establishes is narrower than "MT-FW outperforms
sparse moment-tensor inversion": on the common grid of 3 km and 0.1 s, with the exact predicted
times that the preregistration fixed, BOMP and the mechanism stack lose most of their recall
through timing errors that the peak search of MT-FW tolerates. It does not show that MT-FW outperforms a block-sparse
inversion run on a finer grid or with timing tolerance, as Vera
Rodriguez et al. would do it on real data. The sustainable claim is robustness to grid and timing errors on
a coarse grid, not superiority over sparse moment-tensor inversion as such. Task 2D
tested the attribution with equal tolerance: with the same peak search, BOMP equals MT-FW, and the
advantage does not come from the own elements of MT-FW but only from the peak search (E35).

2C-b is a negative result for the convex formulation as a practical detector on the grid. Group
BLASSO, solved with certified optimality (gaps ≤ 1.4e-8, KKT verified), generally stays 0.2 to
0.6 of recall below its greedy approximation and at the level of BOMP or somewhat below. The performance of MT-FW is not
given by the convex problem of the thesis, but by non-convex ingredients that this problem does not contain: the
peak search, the cap and the greedy deflation. The convex problem remains as a framework and as a certificate (its
dual certificate is the selection statistic and, at μ = 0, the generalized likelihood ratio), not as
the detector that performs best (E33). Its discretization on the grid becomes an open problem: it would
require a convex formulation off the grid, such as sliding Frank-Wolfe, or with timing
tolerance, and that is now the key question of the mathematical line.

2C-c is met in both regimes, with Holm and without it: the advantage over QuakeMigrate generalizes to 30
networks per cell with the network as the unit of replication, which was the objection of the audit about the
inference. The weakest case is strong clutter with 4 stations: in 7331 against QuakeMigrate it passes with the
literal reading and not with Holm, and in that cell MT-FW+SC does not significantly outperform BOMP+SC or
MTSTACK+SC after Holm in either of the two seeds.

**Limitations.** The operating point is chosen on the evaluation itself, the same for all, as in 2A.
BOMP+SC and MTSTACK+SC are our own implementations on the grid and the exact times that the
preregistration fixed, not the published code nor the configuration with which their authors would use them. GL-FW+SC is only the
restriction to the grid: the off-the-grid BLASSO, which is what the thesis formulates, remains untested.
The replication with 30 networks covers MT-FW+SC against QM-3C+SC, BOMP+SC and MTSTACK+SC; 2A-b, 2A-c and 2B still
rest on 5 networks per cell. Synthetic, white noise, a single generator, point double-couple
sources within the subspace of the method, two seeds. The exploratory diagnostics are posterior and do not
support any verdict.

**What remains.** MALMI with its published configuration (Task 1v). Real RSNC noise added to the
synthetics (Tasks 1p and 1z) and the real data (Tasks 1n and 1g). In the mathematical line, sliding
Frank-Wolfe with moment-tensor amplitudes, or a convex formulation with per-station timing
tolerance, compared with MT-FW on these same windows: it decides whether the convex problem can also be the
detector and not only its framework. As a fairer comparator of Vera Rodriguez et al. (2012), a block-sparse
inversion on a finer grid or with timing tolerance; the variant with timing
tolerance was run in Task 2D.

---

## TASK 2D — Equal timing tolerance, ablation of the peak search and frozen conformal threshold — EXECUTED AND VERIFIED IN SYNTHETIC DATA

**Why.** The second independent audit, after 2C, pointed out three gaps. The exact-time
comparators of 2C (BOMP+SC, MTSTACK+SC, GL-FW+SC) did not have the ±0.25 s tolerance to timing errors
that MT-FW has, so 2C-a compared a tolerant method with intolerant methods. The attribution of the
advantage to the peak search rested on a comparison between methods: the ablation of the
peak search and of the cap within MT-FW was missing. And the primary reading of all tasks chose the operating
point on the evaluation itself.

**Rules.** Written in `PREREGISTRO_1j_1k.txt` (section Task 2D, 26 September 2026, 21:13 UTC)
before writing the code, together with the erratum of 2B (E34). First task with clock time and hash:
the SHA-256 of the preregistration at that moment, 6e8510b3…bf88, was recorded in `PREREGISTRO_hashes.txt` at 21:13:04
UTC. The code was written between 02:25 and 02:29 UTC on 27 September, and the hashes of the three
scripts, of `analyze_2d.py` and of the launcher were recorded in `task2d_code_hashes.txt` at 02:41:59 UTC, before
any result. The only addition to the text of 2D is the budget note of 02:41:50 UTC (time from
`date -u`), written before running seeds 2024 and 7331, with its own hash (dc8cf708…e32b): the smoke
test estimated 0.64 h on 2 cores, below the threshold of 6 h, so part (ii) used the 200
planned empty windows per cell and regime. The result paragraph was added at 03:54:55 UTC with a
third hash (b01a2eac…d9a8). The design was written after seeing 2A to 2C and is informed by them, in
particular by the exploratory diagnostics of 2C that pointed to the peak search.

**Design.** Part (i): the same 1200 windows as 2A (seeds 2024 and 7331, both regimes, 3 cells,
networks 0 to 4, 20 windows per network), regenerated bit for bit as in 2C, with the same station
corrections, the grid of 675 nodes and 31 origins, exclusion of 5 km and at most 5 events, and the
OpenBLAS kernel fixed. MT-FW+SC and QM-3C+SC are the saved records of 2A; BOMP+SC, MTSTACK+SC and GL-FW+SC, those
of 2C. The peak reading ("PK") is that of MT-FW: for each station and phase the sample of highest
three-component energy within ±4 samples (0.25 s) of the predicted one is taken, the first in case of a tie. Part
(ii): split-conformal threshold calibrated on the 5 calibration networks of seed 999 of 2A per cell
(same stations, delays and corrections), 40 source-free windows per network, 200 per cell and regime, with
the clutter of the regime and random streams not used before. The score of a window is its maximum (first
score of MT-FW+SC; for QM-3C+SC, the highest score of the events it triggers with the frozen
configuration of 2A, −∞ if it triggers none). The threshold is the ⌈(n+1)(1−α)⌉-th smallest of the n = 200
maxima (the 181st for α = 0.1 and the 151st for α = 0.25), frozen per cell, regime, method and α and
applied without changes to the 100 evaluation windows of each cell and seed, of which 36 to 45 have no
source; an event is declared if its score exceeds the threshold. In the one-sided cell with weak clutter
QuakeMigrate only triggers in 28 of the 200 empty windows, so its threshold at α = 0.25 is −∞.

| Method | What it is | Detail |
|:--|:--|:--|
| BOMP-PK+SC | BOMP of 2C with peak reading | at each step it reads each non-excluded block at the peaks of the current residual; the one with the highest e = hᵀK(I)⁺h enters with its reading times fixed; score e/2 on entry; orthogonal refit of all atoms at their fixed times; the Gram matrix includes the P/S cross term at the offset of the peaks; no cap |
| MTSTACK-PK+SC | MTSTACK of 2C with peak reading | reads the filter output at its peaks, least-squares tensor, ‖Gm‖²/2, maximum over time, suppression at 20 km, no deflation or cap (it is WAVE-MT without cap) |
| GL-FW-PK+SC | GL-FW of 2C with growing dictionary | the Frank-Wolfe oracle evaluates the certificate √(hᵀK(I)⁺h) at the peaks of the current residual; the atom that enters keeps those times; exact restricted problem over the fixed atoms; same path, stop and score as 2C; gap with respect to the realized dictionary |
| MT-FW-EX+SC | MT-FW without peak search | each station and phase is read at the exact predicted sample, with cap; deflation centred there |
| MT-FW-NOCAP+SC | MT-FW without cap | with peak search and without cap |

| Rule | Contrast | Met if |
|:--|:--|:--|
| 2D-a (equal tolerance) | MT-FW+SC − max(BOMP-PK+SC, MTSTACK-PK+SC), part (i) | ≥ 0.05 at cap 0.25 in at least 2 of 3 cells, in both seeds. If it fails, the advantage of MT-FW over Vera Rodriguez et al. (2012) and Anikiev et al. (2014) is attributed to the peak search and to no other element of the method, and the proposal declares it |
| 2D-b (ablation of the peak search) | MT-FW+SC − MT-FW-EX+SC, part (i) | ≥ 0.05 at cap 0.25 in at least 2 of 3 cells, in both seeds; descriptive: MT-FW-NOCAP+SC − MT-FW+SC and GL-FW-PK+SC − MT-FW+SC with CI |
| 2D-c (frozen threshold) | MT-FW+SC and QM-3C+SC with split-conformal threshold, part (ii) | false-event rate of MT-FW+SC in empty windows ≤ 2α and recall of MT-FW+SC − QM-3C+SC ≥ 0.05, in at least 2 of 3 cells, in both seeds, at α = 0.1 and 0.25 |

Statistics: hierarchical bootstrap of 2000 replicates (the 5 networks of the cell with replacement and, within each
network, its windows), with the same resamples for all methods; in part (i) the operating point
and the best comparator are re-chosen in each replicate, in 2D-c the thresholds stay fixed. Holm within
regime over 36 comparisons (2D-a by the difference with the maximum and 2D-b at both caps, 2D-c at both α,
× 3 cells × 2 seeds), one-sided bootstrap p at level 0.025.

**Readings fixed before the results** (docstring of `analyze_2d.py`, with prior hash). 2D-a and 2D-b are
decided on point estimates, threshold 0.05 − 1e-9. In 2D-c an event is declared with a score
strictly greater than the threshold; the "realized false-event rate in empty windows" is read as false events per
empty window (≥ the fraction of empty windows with some false event, which is also reported); the primary
reading requires rate and difference in the same cell, and the strict one (rate ≤ 2α in all cells, seeds and
α) is reported alongside; the uncertainty of the calibration does not enter the CI. The PK variants have no
cap and use the Gram matrix at the peak times with the P/S cross term at the offset of the peaks, which is
our reading of "identical to 2C except that the sample is taken at the peak". In GL-FW-PK two atoms of the same node and origin
with different times can coexist; the record deduplicates them by the 5 km rule;
the Frank-Wolfe steps per λ are bounded at 400, a bound never reached. MT-FW-EX and MT-FW-NOCAP use the dense
loop of MT-FW and not the certified branch and bound; with peak and cap that loop reproduces 2A exactly.

**Status.** Executed and verified with exact controls and validations on seed 999. Scripts:
`task2d_methods.py` (methods), `task2d_run.py` (validation, smoke, runs and controls),
`task2d_validate.py` (validation), `analyze_2d.py` (rules, tables, hierarchical bootstrap and Holm),
`task2d_explore.py` and `task2d_explore_sum.py` (exploratory diagnostics); results, controls and
reproduction mode in `results/task2d/` (`analyze_2d_out.txt`, `task2d_report.txt`, `task2d_validation.txt`,
`README.txt`). Computation: part (i) 3931 s of CPU, part (ii) 1182 s; 1.42 h in total, from 02:42 to 03:25 UTC in 2
processes.

**Validation** (`task2d_validation.txt`, seed 999 only, V1 to V5 correct). The tabulated peaks match
those of `node_vectors_c` in all readings (0 discrepancies in 669 600 to 1 339 200 per cell). With the
peak reading switched off, BOMP-PK, MTSTACK-PK and GL-FW-PK return exactly the records of 2C
(12 of 12), and the dense loop with peak and cap equals certified MT-FW+SC (12 of 12). Against explicit
atoms with shifted times, Aᵀy to 1.2e-15 and AᵀA to 6.4e-16 relative. Independent paths: first
score of NOCAP = maximum of MTSTACK-PK (3.3e-15), of EX = vectorized maximum (4e-16), of BOMP-PK =
λ_max²/2 of GL-FW-PK (2.2e-16). GL-FW-PK with 2 more λ: maximum certificate outside the support / λ =
1.000000000, deviation on the support ≤ 9.6e-10, relative gap ≤ 3.6e-9, record with stop = long
record in 12 of 12. In the 1200 windows: maximum relative gap 1.7e-8, 0 re-solves; in BOMP-PK 91 to
92 % of the readings of the entering atom are shifted from the predicted sample, on average 2.41 to 2.44
samples.

**Controls.** Part (i): stations, delays, catalogue and corrections bit-for-bit identical to 2A (60 of 60
networks); truth, origins and number of sources identical in 1200 of 1200 windows, against 2A and 2C; MT-FW+SC
recomputed equal to what was saved in 120 of 120 windows. The default numpy kernel on this machine changed to
Cooperlake, so it was set per process: SkylakeX in the 18 networks of E31, Prescott in the other 42, each one
checked with threadpoolctl. Part (ii), 30 calibration networks with Prescott: stations and delays equal
to those of 2A, identical corrections (maximum difference 0), and the pipeline reproduces the records of MT-FW+SC and
QM-3C+SC of 2A in windows 0 and 1 of each network (60 of 60 each). No control failed.

**Results of part (i).** Recall at ≤0.1 / ≤0.25, 100 windows per cell:

| Regime, cell (seed) | MT-FW+SC | BOMP-PK+SC | MTSTACK-PK+SC | GL-FW-PK+SC | MT-FW-EX+SC | MT-FW-NOCAP+SC | QM-3C+SC |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| weak, 4 stations (2024) | 0.322 / 0.517 | 0.241 / 0.437 | 0.126 / 0.218 | 0.092 / 0.241 | 0.172 / 0.287 | 0.241 / 0.437 | 0.184 / 0.253 |
| weak, 6 stations (2024) | 0.458 / 0.556 | 0.361 / 0.514 | 0.167 / 0.319 | 0.264 / 0.389 | 0.250 / 0.403 | 0.347 / 0.514 | 0.278 / 0.306 |
| weak, one-sided with 8 (2024) | 0.549 / 0.648 | 0.408 / 0.592 | 0.225 / 0.423 | 0.225 / 0.493 | 0.183 / 0.338 | 0.408 / 0.577 | 0.352 / 0.380 |
| weak, 4 stations (7331) | 0.321 / 0.513 | 0.205 / 0.423 | 0.179 / 0.256 | 0.167 / 0.269 | 0.167 / 0.295 | 0.205 / 0.423 | 0.154 / 0.167 |
| weak, 6 stations (7331) | 0.721 / 0.794 | 0.706 / 0.824 | 0.368 / 0.588 | 0.368 / 0.647 | 0.265 / 0.382 | 0.706 / 0.838 | 0.235 / 0.338 |
| weak, one-sided with 8 (7331) | 0.642 / 0.765 | 0.568 / 0.741 | 0.235 / 0.420 | 0.259 / 0.469 | 0.173 / 0.358 | 0.568 / 0.741 | 0.469 / 0.481 |
| strong, 4 stations (2024) | 0.312 / 0.390 | 0.234 / 0.299 | 0.091 / 0.234 | 0.078 / 0.234 | 0.169 / 0.234 | 0.234 / 0.299 | 0.117 / 0.130 |
| strong, 6 stations (2024) | 0.378 / 0.486 | 0.446 / 0.486 | 0.203 / 0.351 | 0.135 / 0.284 | 0.243 / 0.338 | 0.446 / 0.486 | 0.122 / 0.162 |
| strong, one-sided with 8 (2024) | 0.313 / 0.422 | 0.301 / 0.458 | 0.169 / 0.313 | 0.145 / 0.313 | 0.120 / 0.181 | 0.301 / 0.458 | 0.157 / 0.193 |
| strong, 4 stations (7331) | 0.128 / 0.218 | 0.167 / 0.269 | 0.090 / 0.205 | 0.103 / 0.231 | 0.141 / 0.154 | 0.167 / 0.269 | 0.038 / 0.051 |
| strong, 6 stations (7331) | 0.357 / 0.417 | 0.298 / 0.571 | 0.179 / 0.274 | 0.167 / 0.298 | 0.107 / 0.250 | 0.286 / 0.548 | 0.179 / 0.250 |
| strong, one-sided with 8 (7331) | 0.627 / 0.707 | 0.653 / 0.760 | 0.427 / 0.547 | 0.347 / 0.573 | 0.200 / 0.413 | 0.653 / 0.760 | 0.253 / 0.293 |

Differences at cap 0.25 with hierarchical 95 CI (the best comparator of 2D-a is BOMP-PK+SC in all cells
and at both caps):

| Regime, cell (seed) | 2D-a: MT-FW+SC − BOMP-PK+SC | met | 2D-b: MT-FW+SC − MT-FW-EX+SC | met | MT-FW-NOCAP+SC − MT-FW+SC | GL-FW-PK+SC − MT-FW+SC |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| weak, 4 stations (2024) | +0.080 [−0.136, +0.282] | yes | +0.230 [+0.028, +0.388] | yes | −0.080 [−0.286, +0.133] | −0.276 [−0.510, −0.027] |
| weak, 6 stations (2024) | +0.042 [−0.092, +0.224] | no | +0.153 [−0.118, +0.348] | yes | −0.042 [−0.226, +0.092] | −0.167 [−0.397, +0.013] |
| weak, one-sided with 8 (2024) | +0.056 [+0.000, +0.265] | yes | +0.310 [+0.154, +0.515] | yes | −0.070 [−0.250, +0.000] | −0.155 [−0.488, −0.012] |
| weak, 4 stations (7331) | +0.090 [−0.029, +0.256] | yes | +0.218 [−0.038, +0.378] | yes | −0.090 [−0.256, +0.029] | −0.244 [−0.392, −0.034] |
| weak, 6 stations (7331) | −0.029 [−0.125, +0.047] | no | +0.412 [+0.108, +0.597] | yes | +0.044 [−0.042, +0.129] | −0.147 [−0.347, −0.014] |
| weak, one-sided with 8 (7331) | +0.025 [−0.084, +0.155] | no | +0.407 [+0.212, +0.619] | yes | −0.025 [−0.168, +0.082] | −0.296 [−0.443, −0.032] |
| strong, 4 stations (2024) | +0.091 [−0.048, +0.257] | yes | +0.156 [+0.041, +0.355] | yes | −0.091 [−0.253, +0.048] | −0.156 [−0.354, −0.045] |
| strong, 6 stations (2024) | +0.000 [−0.138, +0.145] | no | +0.149 [+0.000, +0.327] | yes | +0.000 [−0.143, +0.141] | −0.203 [−0.403, +0.000] |
| strong, one-sided with 8 (2024) | −0.036 [−0.130, +0.108] | no | +0.241 [+0.025, +0.433] | yes | +0.036 [−0.121, +0.125] | −0.108 [−0.359, +0.029] |
| strong, 4 stations (7331) | −0.051 [−0.206, +0.095] | no | +0.064 [−0.096, +0.240] | yes | +0.051 [−0.096, +0.197] | +0.013 [−0.141, +0.136] |
| strong, 6 stations (7331) | −0.155 [−0.277, +0.032] | no | +0.167 [−0.024, +0.410] | yes | +0.131 [−0.092, +0.253] | −0.119 [−0.320, +0.107] |
| strong, one-sided with 8 (7331) | −0.053 [−0.158, +0.030] | no | +0.293 [+0.098, +0.472] | yes | +0.053 [−0.027, +0.159] | −0.133 [−0.321, +0.016] |

Effect of the peak reading at 0.25 (PK minus the exact-time version of 2C): BOMP +0.04 to +0.51,
MTSTACK +0.05 to +0.37, GL-FW +0.03 to +0.44, positive in 12 of 12 cases each. GL-FW-PK+SC − BOMP-PK+SC:
−0.04 to −0.27 (12 of 12 negative). max(BOMP-PK+SC, MTSTACK-PK+SC) − QM-3C+SC: +0.17 to +0.49 (12 of 12).
MT-FW+SC − MTSTACK-PK+SC: +0.21 to +0.35 in weak and +0.01 to +0.16 in strong.

**Results of part (ii).** Thresholds frozen on 200 empty windows of seed 999 per cell,
applied to the windows of 2A; hierarchical CI with the thresholds fixed:

| Regime, cell (seed) | α | recall MT-FW+SC / QM-3C+SC | difference [H CI] | false events of MT-FW+SC per empty window [H CI] | false events per window, all, MT-FW+SC / QM-3C+SC | counts |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| weak, 4 stations (2024) | 0.1 | 0.667 / 0.230 | +0.437 [+0.275, +0.565] | 0.179 [0.034, 0.351] | 0.840 / 0.170 | yes |
| weak, 4 stations (2024) | 0.25 | 0.690 / 0.241 | +0.448 [+0.275, +0.592] | 0.410 [0.189, 0.711] | 1.090 / 0.200 | yes |
| weak, 6 stations (2024) | 0.1 | 0.639 / 0.292 | +0.347 [+0.121, +0.532] | 0.045 [0.000, 0.122] | 0.730 / 0.150 | yes |
| weak, 6 stations (2024) | 0.25 | 0.681 / 0.306 | +0.375 [+0.190, +0.534] | 0.159 [0.000, 0.400] | 0.980 / 0.230 | yes |
| weak, one-sided with 8 (2024) | 0.1 | 0.662 / 0.366 | +0.296 [+0.100, +0.507] | 0.070 [0.000, 0.193] | 0.750 / 0.150 | yes |
| weak, one-sided with 8 (2024) | 0.25 | 0.690 / 0.380 | +0.310 [+0.123, +0.500] | 0.186 [0.036, 0.405] | 1.100 / 0.160 | yes |
| weak, 4 stations (7331) | 0.1 | 0.590 / 0.154 | +0.436 [+0.300, +0.579] | 0.179 [0.061, 0.324] | 0.810 / 0.220 | yes |
| weak, 4 stations (7331) | 0.25 | 0.628 / 0.179 | +0.449 [+0.319, +0.581] | 0.410 [0.163, 0.718] | 1.020 / 0.310 | yes |
| weak, 6 stations (7331) | 0.1 | 0.824 / 0.338 | +0.485 [+0.289, +0.662] | 0.222 [0.077, 0.411] | 0.870 / 0.230 | no |
| weak, 6 stations (7331) | 0.25 | 0.824 / 0.338 | +0.485 [+0.289, +0.662] | 0.422 [0.233, 0.684] | 1.170 / 0.300 | yes |
| weak, one-sided with 8 (7331) | 0.1 | 0.840 / 0.481 | +0.358 [+0.172, +0.565] | 0.053 [0.000, 0.182] | 1.180 / 0.150 | yes |
| weak, one-sided with 8 (7331) | 0.25 | 0.840 / 0.481 | +0.358 [+0.172, +0.565] | 0.263 [0.075, 0.526] | 1.460 / 0.190 | yes |
| strong, 4 stations (2024) | 0.1 | 0.221 / 0.130 | +0.091 [−0.033, +0.211] | 0.025 [0.000, 0.103] | 0.050 / 0.190 | yes |
| strong, 4 stations (2024) | 0.25 | 0.403 / 0.143 | +0.260 [+0.157, +0.371] | 0.200 [0.024, 0.500] | 0.290 / 0.340 | yes |
| strong, 6 stations (2024) | 0.1 | 0.419 / 0.135 | +0.284 [+0.162, +0.431] | 0.000 [0.000, 0.000] | 0.120 / 0.140 | yes |
| strong, 6 stations (2024) | 0.25 | 0.514 / 0.162 | +0.351 [+0.208, +0.513] | 0.023 [0.000, 0.111] | 0.340 / 0.260 | yes |
| strong, one-sided with 8 (2024) | 0.1 | 0.530 / 0.193 | +0.337 [+0.149, +0.547] | 0.048 [0.000, 0.159] | 0.410 / 0.140 | yes |
| strong, one-sided with 8 (2024) | 0.25 | 0.590 / 0.193 | +0.398 [+0.231, +0.612] | 0.190 [0.021, 0.458] | 0.680 / 0.270 | yes |
| strong, 4 stations (7331) | 0.1 | 0.115 / 0.038 | +0.077 [+0.000, +0.186] | 0.000 [0.000, 0.000] | 0.070 / 0.150 | yes |
| strong, 4 stations (7331) | 0.25 | 0.218 / 0.051 | +0.167 [+0.056, +0.310] | 0.108 [0.000, 0.381] | 0.260 / 0.220 | yes |
| strong, 6 stations (7331) | 0.1 | 0.500 / 0.238 | +0.262 [+0.062, +0.443] | 0.028 [0.000, 0.118] | 0.390 / 0.180 | yes |
| strong, 6 stations (7331) | 0.25 | 0.595 / 0.250 | +0.345 [+0.130, +0.526] | 0.167 [0.000, 0.400] | 0.620 / 0.220 | yes |
| strong, one-sided with 8 (7331) | 0.1 | 0.773 / 0.293 | +0.480 [+0.259, +0.664] | 0.075 [0.000, 0.275] | 0.560 / 0.230 | yes |
| strong, one-sided with 8 (7331) | 0.25 | 0.773 / 0.293 | +0.480 [+0.259, +0.664] | 0.375 [0.065, 0.907] | 0.940 / 0.320 | yes |

In its own calibration, the fraction of empty windows with some false event of MT-FW+SC is 0.095 at α = 0.1 and
0.245 at α = 0.25 in all cells. In the evaluation, the false events per empty window exceed α in 3 of 6
weak cells at α = 0.1 (up to 2.2 α) and in 4 of 6 at α = 0.25 (up to 1.7 α); in strong, in 0 of 6 and 1 of 6.

**Verdicts.**

| Rule | Weak | Strong | Detail |
|:--|:-:|:-:|:--|
| 2D-a | not met | not met | weak: 2 of 3 cells (2024) and 1 of 3 (7331); strong: 1 of 3 and 0 of 3; differences of −0.155 to +0.091 at 0.25, all hierarchical CIs contain 0; re-applied per replicate it is met in 11 % (weak) and 0.1 % (strong); with Holm, 0 of 3 throughout |
| 2D-b | met | met | 3 of 3 in each seed; +0.064 to +0.412 at 0.25, positive in 12 of 12; 99 % and 95 % of the replicates; lower bound > 0 in 7 of 12; with Holm, 1 and 0 of 3 (weak), 0 and 1 of 3 (strong): not met |
| 2D-c, per cell | met | met | weak: 3 and 2 of 3 at α = 0.1, 3 and 3 at α = 0.25; strong: 3 of 3 always; the recall condition alone is met in 12 of 12 cases per regime (+0.08 to +0.49); 76 % and 96 % of the replicates |
| 2D-c, strict | not met | met | weak, 7331, 6 stations, α = 0.1: 10 false events in 45 empty windows = 0.222 > 0.2 (one false event too many) |
| 2D-c, with Holm | not met | not met | with 5 networks per cell |

No case came within 0.005 of the threshold of 0.05.

**Exploratory diagnostics** (posterior, not preregistered; no verdict depends on them). E1
(`task2d_explore_e1.txt`): in the windows of 2A with events, which calibration on empty windows does not see,
MT-FW+SC produces with the conformal thresholds 1.21 to 1.87 false events per window in weak at α = 0.1 and 1.41 to
2.19 at α = 0.25, 78 to 92 % of them at 5-20 km from a true source (residual energy of a real
event after deflation, the 5 km exclusion and the model mismatch), against 0.18 to 0.36 for QM-3C+SC; in
strong, 0.07 to 0.88 and 0.35 to 1.32, against 0.17 to 0.38. The conformal thresholds of MT-FW+SC (23.7 to 29.8 in
weak, 36.2 to 52.1 in strong) lie below those of 2A-e (34.6 to 42.9 and 44.7 to 53.2), calibrated
on all windows with a budget of false events per window. E3 (last table of
`analyze_2d_out.txt`): the same thresholds on the 30 networks per cell of 2C (150 to 169 empty windows
per cell) give MT-FW+SC − QM-3C+SC of +0.36 to +0.47 in weak and +0.03 to +0.40 in strong, with false events of
MT-FW+SC per empty window of 0.055 to 0.139 (α = 0.1) and 0.21 to 0.35 (α = 0.25) in weak and 0.006 to 0.092 and
0.10 to 0.33 in strong: close to α or somewhat above, and below 2α in all. E2
(`task2d_explore_e2.txt`, which repeats V6 with 90 windows per cell): ideal generator of 1x with the sources
on the nodes and the origins on the grid, and an injected timing error of ±0.2 s per station and phase
(random sign, 3.2 samples, within the search window). The exact-time methods lose
0.16 to 0.48 of recall with the error; the PK variants lose −0.03 to +0.15 and outperform their exact version with
the error in 12 of 12 cases (+0.04 to +0.38). The price: without error, with 6 and 8 stations, the peak search
costs recall (MT-FW 0.83 to 0.87 against 0.97 for MT-FW-EX; BOMP-PK 0.83 to 0.88 against 0.97 to 0.99 for BOMP;
MTSTACK-PK 0.62 to 0.66 against 0.84 to 0.85), plausibly because the ±0.25 s window latches onto
clutter or coda energy (not verified). MT-FW-NOCAP and BOMP-PK are almost identical here too.

**Interpretation.** 2D-a clearly fails: with the same ±0.25 s tolerance, the block inversion in
the form of Vera Rodriguez et al. (2012) equals MT-FW (in strong, equal or better in 5 of 6 cells). The
margins of 2C-a (+0.08 to +0.52) were almost entirely the missing tolerance. As preregistered,
the advantage of MT-FW over Vera Rodriguez et al. and Anikiev et al. is attributed to the peak search and to
no other element of the method, and the proposal declares it; the reading of 2C, which attributed it to the
own elements of MT-FW, is refuted (E35). MTSTACK-PK falls below MT-FW: what separates it is
deflation, not tolerance. 2D-b is met: the peak search matters, with constant sign but weak
support per cell with 5 networks. Without cap, MT-FW is in essence BOMP with peak reading (they differ only in the
P/S terms of the Gram matrix and in deflation versus orthogonal refit); the cap is the only element of MT-FW
beyond BOMP with tolerance and its effect is small and not significant. GL-FW-PK improves on GL-FW
but remains 0.11 to 0.30 below MT-FW in 11 of 12 cases and below BOMP-PK in 12 of 12: the
conclusion of 2C about the convex problem on the grid holds with equal tolerance, and its certificate
is valid only for the dictionary that the heuristic oracle builds. With the tolerance, BOMP-PK also outperforms
QuakeMigrate (+0.17 to +0.49): the advantage over the state of the art of 2A and 2C belongs to the family of
moment-tensor detectors with timing tolerance and deflation, not to MT-FW in particular. 2D-c is
met per cell, but the conformal threshold calibrated on empty windows is not a deployable operating
point: in windows with events MT-FW+SC gives 1.2 to 2.2 false events per window in weak, and over all
windows 3 to 8 times the false events of QuakeMigrate; and with 5 calibration networks exchangeability only
holds approximately.

Consequence for the thesis. The distinctive contribution is no longer the MT-FW detector and becomes the
principled treatment of timing tolerance: group BLASSO with a time shift per
station and phase as a nuisance parameter, solved off the grid by sliding Frank-Wolfe, with its
certificate and its theory. It only counts as a methodological advance if it outperforms at the same time the heuristic
peak search of MT-FW and BOMP-PK, which becomes the comparator to beat; the benchmark and the test with real
data G1' are contributions whatever the result.

**Edit after the results.** A diagnostic line of `analyze_2d.py` failed (the lists of
shifts of BOMP-PK have length 2 ns per atom and had to be flattened) and was changed after the
part files existed; no rule, threshold, bootstrap or table changed. Hash before:
110aa893…73a7 (`task2d_code_hashes.txt`, 02:41:59 UTC); after: 8dd5013b…d2b4. The result paragraph of the
preregistration declares it.

**Limitations.** A single heuristic tolerance (argmax per station and phase at ±0.25 s) on a single grid
of 3 km and 0.1 s; no method with time shifts as nuisance parameters, no other
width or finer grid. BOMP-PK and MTSTACK-PK are our own implementations with the tolerance of MT-FW
added, not the code of their authors. Five networks per cell in both parts, which is why no rule survives
Holm. Thresholds calibrated only on empty windows. Synthetic, one generator, double-couple sources within
the subspace of the method, two seeds.

**What remains.** Implement group BLASSO with per-station and per-phase shifts off the grid
(sliding Frank-Wolfe) and compare it on these windows with MT-FW and BOMP-PK: it is the first mathematical
deliverable. A deployable operating point for a budget of false events, calibrated with windows that
contain events. MALMI, real noise and real data, as in 2C. (BLASSO with shifts was implemented and
tested in Task 2E, below.)

---
## TASK 2E — Group BLASSO with per-station and per-phase shifts, off the grid (G1-M) — EXECUTED AND VERIFIED ON SYNTHETICS

**Why.** Task 2D left a single methodological question: whether the principled form of timing tolerance, the
group BLASSO (16a) with one shift per station and phase as a nuisance parameter, solved off the grid by
sliding Frank-Wolfe, detects more events than BOMP-PK+SC with the same tolerance of ±0.25 s and at an equal
realised false-event rate. This is hypothesis G1-M.

**Rules.** Written in `PREREGISTRO_1j_1k.txt` (Task 2E section, 27 September 2026, 14:44:14 UTC, wall-clock
time) before any 2E code was written. The SHA-256 of the preregistration at that moment,
b40f0c6e…e988, was delivered to the user through the chat at the time of registration (a timestamp external to
the container, although neither public nor third-party) and was recorded in `PREREGISTRO_hashes.txt`. The only
addition before evaluation is the calibration note required by the registration (15:38:20 UTC, f = 0.01, hash
fada69b1…9271), before any computation on 2024/7331. Code hashes (5 files) in `task2e_code_hashes.txt` at
15:38:31 UTC; first evaluation part at 15:39:40 UTC, last at 17:39:19 UTC; `sha256sum -c` correct after the
run. Result paragraph at 17:45:17 UTC with another hash. The budget rule (10 h on 2 cores) was not triggered:
estimate 1.98 h.

| Rule | Contrast | Met if |
|:--|:--|:--|
| 2E-a (primary) | TS-SFW+SC − BOMP-PK+SC | ≥ 0.05 at cap 0.25 and lower bound of the hierarchical 95 CI > 0, in at least 2 of 3 cells in both seeds, per regime; strict reading with Holm over the 6 cell × seed tests (one-sided 0.025) |
| 2E-b (secondary) | TS-SFW+SC − BOMP-PK-SLIDE+SC | ≥ 0.05 at cap 0.25 in at least 2 of 3 cells in both seeds (point estimate) |

Registered interpretation: if 2E-a fails in both regimes, G1-M is refuted in this synthetic design and the
proposal states so and reformulates its methodological contribution.

**Design.** 30 networks per cell, seed and regime (the 5 from 2A with 20 windows and the 25 from 2C part (ii)
with 12), 400 windows per cell, 4800 in total, regenerated bit for bit with the OpenBLAS kernel pinned; same
corrections, 675-node grid, origins, tolerance, 5 km exclusion, at most 5 events and hit rule. TS-SFW+SC:
continuous atom (position, origin and one shift per station and phase bounded to ±0.25 s); at each iteration
(a) oracle = MT-FW-NOCAP statistic on the grid, (b) L-BFGS-B refinement of the new atom (≤ 50 iterations), (c)
sliding, alternating the convex step in amplitudes (relative gap ≤ 1e-6) and L-BFGS-B in positions, origins
and shifts of all atoms (≤ 50 iterations per step, ≤ 5 alternations), (d) residual. λ = f · s₁.
BOMP-PK-SLIDE+SC = the same with λ = 0. BOMP-PK+SC is the 2D code unchanged (recomputed identical to 2D in
1200/1200 windows). Calibration on seed 999 (600 windows): recall at 0.25 of 0.4749, 0.4644, 0.2510 and
0.1841 for f = 0.01, 0.03, 0.1 and 0.3; f = 0.01 was frozen, where λ is a median of 0.045 of the first atom's
certificate and no atom ends up at zero (at f = 0.3, 88 % do).

**Validation.** Continuous identities ≤ 8.9e-16; gradients against central differences 1.4e-9; convex step
against an explicit dictionary and against FISTA; with refinement and sliding switched off and λ = 0 the
first atom is that of BOMP-PK in 24/24 windows. The noise-free test with an off-grid source revealed a
property of the method: starting from the truth the optimiser stays there (94/94), from 0.2-0.4 km away it
recovers it in 85/94, but from the grid node it reaches < 0.5 km and < 0.02 s in only 12/90 sources (cycle
skips of ~0.45 s). The registered plan was followed because the registered method is a local refinement.

**Results** (`analyze_2e_out.txt`; differences at cap 0.25 with hierarchical 95 CI).

| Regime | 2E-a: TS-SFW+SC − BOMP-PK+SC | 2E-b: TS-SFW+SC − BOMP-PK-SLIDE+SC |
|:--|:--|:--|
| weak (2024) | +0.039 [−0.025, +0.122], +0.029 [−0.025, +0.087], +0.025 [−0.042, +0.094] | −0.003, −0.003, +0.000 |
| weak (7331) | +0.016 [−0.038, +0.094], +0.029 [−0.024, +0.078], +0.015 [−0.035, +0.091] | +0.013, −0.006, −0.012 |
| strong (2024) | +0.000 [−0.033, +0.030], +0.033 [−0.013, +0.083], +0.016 [−0.026, +0.059] | −0.003, +0.007, −0.020 |
| strong (7331) | +0.003 [−0.027, +0.039], +0.015 [−0.024, +0.058], +0.018 [−0.044, +0.080] | −0.003, −0.006, −0.009 |

(Cells in the order 4 stations, 6 stations, one-sided of 8.) 2E-a is not met in any regime: 0 of 3 cells in
each seed; the 12 points are ≥ 0 but none reaches 0.05 and all CIs contain 0; minimum one-sided p 0.135
(weak) and 0.098 (strong), minimum Holm-adjusted 0.81 and 0.58. Nor is 2E-b (−0.020 to +0.013): the group
penalty is inert at the calibrated f. Descriptive: sliding alone (BOMP-PK-SLIDE+SC − BOMP-PK+SC) gives +0.003
to +0.042; TS-SFW+SC − MT-FW+SC, −0.026 to +0.049 (all CIs contain 0). Median epicentral error of hits: 1.79
versus 1.88 km (weak) and 2.01 versus 2.15 km (strong); the continuous solutions end up 0.3 to 1.3 km too
deep. 73-83 % of the sliding steps and 83-97 % of the refinements stop at 50 iterations; 1.19 s per window
versus 0.42 s for BOMP-PK+SC.

**Independent verification** (`verif_2e/task2e_verification.txt`; new computation only on seed 999, no 2E or
registration file touched). Verdict: reliable negative. (1) No implementation error: at a node with the shifts
of the peaks, the continuous atom reproduces the BOMP-PK block (prediction 1.4e-15, Gram 9.0e-16, correlation
5.7e-16, joint objective 9.9e-16, 120 atoms); objective against an independent implementation 2.3e-13;
gradient against finite differences 4.5e-8; signs of injected shifts recovered to 1e-4 s. (2) Numbers
re-derived with independent FROC and bootstrap. (3) No budget disadvantage: with 10 times the iterations, 4
times the alternations and refinement with 9 starts (270 windows, 208 events), the objective drops by a median
of 1.9 % but recall at 0.25 changes by −0.014 [−0.052, +0.038] (6 events gained, 4 lost). (4) Large margin:
the oracle ceilings CEIL-G and CEIL-GD lie 0.16-0.37 and 0.08-0.24 above BOMP-PK at 0.25 per cell; CEIL-C ≈
CEIL-G (0.699 versus 0.692). 89 % of the events have a score near the truth higher than any clutter and noise
candidate, but BOMP-PK only places an atom within < 5 km of 77 % of them (TS-SFW, 74 %); the detectable events
never hit have their nearest atom at a median of 7.3 / 6.8 km, 80 / 85 % at 5-10 km: absorbed by misplaced
atoms just beyond the 5 km tolerance. Sliding moves the atoms by a median of 0.7 km. (5) Registration hash
chain intact (6 growing prefixes).

**Consequence for the thesis.** G1-M is refuted in this synthetic design, as registered. No convex formulation
of the thesis (on the grid or off it, with or without tolerance) outperforms block search with the same
tolerance. The group BLASSO remains as the framework and object of the theory (Appendix AA), not as the
detector.

**Limitations.** A single registered method (local refinement from a 3 km grid with data at 2 Hz); no global
relocation or other tolerances; ceilings with 100 windows per cell on seed 999; synthetic.

---

## EXPLORATION `explore_contrib` — Where the margin is not: seven candidates on the calibration seed — EXPLORATORY

**Status.** Exploratory, not confirmatory: seed 999 only, on networks of that seed not used before (10-15,
20-29, 30-35 and 40-59; the calibration ones 0-4 were only read), new random streams, and seeds 2024, 7331,
4242 and 8888 were never generated or read. Several numbers were seen while the designs were being chosen;
they are pilot estimates (`results/explorations/REPORT.txt`, `results/explorations/REPORT_C7.txt`).

| Candidate | Main contrast, cap 0.1 / 0.25 [95 CI] | Reading |
|:--|:--|:--|
| C1 coarse grids | BOMP-PK+SC at 9 km − at 3 km: −0.324 [−0.524, −0.052] / −0.379 [−0.512, −0.174] | TS-SFW on the same coarse grid recovers +0.03 to +0.10 and never reaches BOMP-PK at 3 km; null |
| C2 dense overlaps (3-5 events at < 9 km) | MT-FW+SC − BOMP-PK+SC: +0.074 [+0.016, +0.136] / +0.070 [+0.020, +0.128] (360 windows, 1439 events) | the only one meeting the exploratory criterion; the ablation attributes it to the cap |
| C2 on new networks (40-59) | +0.077 [+0.030, +0.122] / +0.068 [+0.029, +0.109] (600 windows) | replication before registering |
| C3 large timing errors | BOMP-PK-12 − BOMP-PK with 0.8 s statics without corrections: +0.156 [+0.015, +0.284] / +0.148 [+0.000, +0.315] | only helps without catalogue corrections; null |
| C6 travel-time physics | oracle medium − homogeneous (BOMP-PK): +0.006 [−0.066, +0.065] / +0.011 [−0.059, +0.064] | null; the calibrated gradient model subtracts −0.104 [−0.202, +0.013] at 0.25 with strong clutter |
| C7 source-pulse invariance | BOMP-PK-PH − BOMP-PK: −0.114 [−0.203, −0.017] / −0.084 [−0.148, −0.023] | negative; no gain in dense windows or with pulses outside the family |
| Arrival masking | +0.149 [−0.007, +0.265] / +0.084 [+0.000, +0.163] (weak, replication) | doubtful, weak clutter only |
| Speed | MT-FW+SC − BOMP-PK+SC in standard windows: +0.017 / +0.002 | same recall at 6.3 times lower cost; engineering |

**What it says.** (i) The grid and the timing error absorbed by peak search are not the margin. (ii) Nor is
the velocity model: with the generator's medium (oracle) or calibrated by JHD with the 20 catalogue events,
nothing improves, and the C6b attribution run shows that, with weak clutter, the realistic margin of BOMP-PK
remains at +0.194 [+0.081, +0.312] when the generator's medium is exactly the detector's homogeneous one
(+0.177 with fitted pulse, +0.153 without site and polarisation perturbations, +0.242 without coda). Only the
oracle medium with the true statics wins (+0.069 [+0.004, +0.141] at 0.25), and that is not attainable.
(iii) Nor is the source pulse: without noise and with all of the generator's realism, the true block captures
only half of the energy (0.497) and pulse invariance adds 0.014 to 0.050; the invariant detector places the
first atom more often at more than 5 km. (iv) Exact by-product: without the P/S cross term of its Gram, block
search reproduces BOMP-PK in 2400/2400 windows at 10.7 times lower cost. (v) In C2 the ablation isolates the
per-station cap: MT-FW without cap and BOMP-PK with sequential deflation give identical hit sequences in
360/360 windows, and placing the cap in BOMP-PK's selection statistic recovers most of the difference.

**Multiplicity.** Seven candidates were examined (C1, C2, C3, C6, C7, masking with 12 variants in its pilot,
and speed); C4 and C5 were written down and never run. The cap in dense windows was chosen for Task 2F among
them because it was the only one meeting the exploratory criterion; its pilot estimate is subject to the
winner's curse and is not evidence. What protects 2F from that selection: seeds untouched by the exploration,
design and code frozen before generating them, and a planning effect of half the pilot.

---

## TASK 2F — The per-station cap in dense sequences (G1-R) — EXECUTED AND VERIFIED ON SYNTHETICS

**Why.** To test on new seeds the only candidate from the exploration that met the criterion: in windows with 3
to 5 overlapping events, does MT-FW+SC (moment-tensor subspace statistic with a per-station cap and sequential
deflation) detect more events than BOMP-PK+SC at an equal false-event rate, and is it due to the cap? (H-CAP;
in the proposal, G1-R.)

**Rules.** Written in `PREREGISTRO_1j_1k.txt` (Task 2F section, 28 September 2026, 03:34:51 UTC) before running
any code on seeds 4242 and 8888, which had never been used; the registration declares that the hypothesis was
chosen among the explored candidates. SHA-256 3d8e4c9b…2f5f, delivered to the user through the chat before the
run. Controls on seed 999 from 03:39:10 UTC; hashes of the 26 loaded modules (analysis included) in
`task2f_code_hashes.txt` at 03:49:08 UTC; single launch at 03:49:14 UTC; parts between 03:49:29 and 06:52:31
UTC; `sha256sum -c` 26/26 at 06:53:11 UTC; analysis run once at 06:54:13 UTC. The registration did not change:
only the result (06:59:18 UTC) and the verification (07:26:18 UTC) were added, each with a hash.

| Rule | Contrast | Met if |
|:--|:--|:--|
| 2F-a (primary) | MT-FW+SC − BOMP-PK+SC at cap 0.25, dense windows, pooling the 3 cells | lower bound of the 95 CI > 0 and point ≥ 0.03 in both seeds, per regime; strict reading with Holm over 4 tests (one-sided 0.025) |
| 2F-b (mechanism) | MT-FW+SC − MT-FW-NOCAP+SC and MT-FW-NOCAP+SC − BOMP-PK+SC | lower bound > 0 in both seeds (pooling regimes) and \|NOCAP − BOMP-PK\| ≤ 0.02 in both |
| 2F-c (no harm) | MT-FW+SC − BOMP-PK+SC in standard windows | TOST equivalence with margin 0.05, pooling cells and regimes, in both seeds |

**Design.** Seeds 4242 and 8888, each an independent replicate. Dense arm: networks j = 0..95 (96 per cell,
regime and seed), 10 windows per network with 3 to 5 events separated by ≥ 5 km within 9 km, origins in
[5.2, 7.8] s with gaps ≥ 0.5 s, new TAG 4201, kmax 7: 5760 windows per seed (23013 and 22996 events).
Equivalence arm: j = 96..111 (16 per cell and regime), standard windows with 0 to 2 events, TAG 4202, kmax 5:
960 windows per seed (766 and 749 events). Methods: MT-FW+SC (cap CAPV = 6, fixed in 1x and not tuned here),
BOMP-PK+SC, BOMP-PK-CAP+SC (the same cap in the selection statistic only) and MT-FW-NOCAP+SC. Hierarchical
bootstrap with 2000 replicates (cells as strata, networks and windows). Size: with the pilot SDs (0.0383 weak,
0.0353 strong), 96 networks per cell give power ~0.95 per seed for δ = 0.035. Controls on seed 999
(`task2f_control.txt`): the dense pilot is reproduced bit for bit (180/180 windows per regime, all five
methods) and MT-FW-NOCAP = BOMP-PK with sequential deflation in 360/360.

**Results** (`analyze_2f_out.txt`, cap 0.25, pooling cells).

| Rule | Regime | 4242 | 8888 | Met |
|:--|:--|:--|:--|:--|
| 2F-a | weak | +0.072 [+0.048, +0.092] | +0.074 [+0.052, +0.100] | yes, also with Holm (adjusted p 0.002) |
| 2F-a | strong | +0.051 [+0.032, +0.070] | +0.043 [+0.022, +0.066] | yes, also with Holm |
| 2F-b | MT-FW − NOCAP | +0.060 [+0.046, +0.076] | +0.059 [+0.044, +0.075] | yes |
| 2F-b | NOCAP − BOMP-PK | −0.001 | −0.002 | yes |
| 2F-c | standard, 90 CI | −0.008 [−0.039, +0.034] | −0.013 [−0.048, +0.027] | yes; narrow in 8888 (p = 0.0405) |

At cap 0.1: +0.085 / +0.089 (weak) and +0.055 / +0.060 (strong). Per cell, positive in 12 of 12, from +0.023
[−0.009, +0.051] (strong, 8888, 4 stations, the only CI containing 0) to +0.110 (weak, 4242, 6 stations).
BOMP-PK-CAP+SC recovers 85-90 % of the difference. MT-FW costs 0.105-0.118 s per dense window versus
0.640-0.674 s for BOMP-PK (about 6 times less). Absolute recall with cap: only 0.22-0.29 at 0.25 in dense
windows. In standard windows with strong clutter the points are negative (−0.026, −0.032; −0.067 in the
4-station cell): a cost of ~0.03 on isolated events with strong clutter is not ruled out.

**Independent verification** (`verif_2f/task2f_verification.txt`): positive confirmed, with no leakage,
asymmetry, scoring artefact or provenance problem (hashes 26/26, chain of 8 prefixes, 3d8e4c9b ends exactly at
the end of the Task 2F section, no prior data from 4242/8888, 3 parts regenerated with BOMP-PK before MT-FW
identical in 30/30 windows, MT-FW reimplemented bit for bit in 216/216). It narrows the scope on four points:

| Point | Result |
|:--|:--|
| Hit rule with origin time | with 5 km and 0.5 s the gain is +0.045 / +0.024 / +0.043 / +0.027 (weak 4242, strong 4242, weak 8888, strong 8888); with 3 km and 0.5 s, +0.035 / +0.013 / +0.027 / +0.020; about half of the location-only «hits» have origin error > 0.5 s |
| Cap level (exploratory, seed 999) | cap 0.1 / 0.25: CAPV 2 −0.064/−0.122, 3 +0.004/−0.056, 4 +0.040/+0.017, 5 +0.058/+0.047, 6 +0.074/+0.070, 8 +0.058/+0.046, 12 +0.033/+0.010, 24 −0.003/−0.008, no cap −0.003/+0.001: peak at the fixed value |
| Mechanism (exploratory) | without clutter the gain rises (+0.137 [+0.058, +0.183]); without clutter and without lognormal site factors it falls to +0.023 [−0.040, +0.069]: bounded influence of each station-and-phase block against the model's amplitude error, not robustness to clutter; those factors are drawn in every window, so they are per-window amplitude perturbations and not station terms (exploration C8) |
| Within-window versus between-window | distinct hits among the first atoms +0.001 to +0.012; among the top 20 % of first-atom scores, MT-FW hits 0.90-0.96 and BOMP-PK 0.79-0.84: the improvement is mainly in the calibration of the score across windows, not in event separation |

**Consequence for the thesis.** The per-station cap is a confirmed contribution for dense sequences in this
synthetic design, modest, due to the cap bounding the influence of any station-and-phase block whose amplitude
the radiation model predicts poorly, and tied to the cap level fixed in 1x. The generator's amplitude factors
are drawn in every window: removing them reduces the gain, but with small dispersion the cap still helps where
knowing the true factors does not (exploration C8), so the mechanism is bounded per-block influence against
amplitude error, not correction of a site term. Clipping per-station contributions is a known practice in
robust stacking; what is claimed is its measured and explained effect within a moment-tensor detector with
deflation. The direction the evidence supports is a robust bounded-influence certificate whose level is set
from the noise distribution, registered as Task 2I; per-station amplitude or site terms estimated from the
catalogue do not replace the cap (exploration C8).

**Limitations.** Synthetic; cap level fixed in 1x; only the generator's lognormal amplitude factors, drawn per
window; up to 5 events within 9 km; possible small cost on isolated events with strong clutter; mechanism
readings exploratory and on seed 999.

**What remains.** The bounded-influence certificate with noise-derived level against the capped certificate
and BOMP-PK, in the 2F design, with the hit rule with origin time as a secondary reading (Task 2I, executed and verified, with its secondary arms analysed; see its section). The margin of 0.08-0.24 below the realistic ceiling, which neither the velocity model nor
the pulse explains and of which no amplitude treatment captures more than a quarter, remains open.
Real data, as always.

---

## EXPLORATION C8 — Site terms versus bounded influence — EXPLORATORY

**Status.** Exploratory, not confirmatory: seed 999 only, on networks of that seed not used before (60-79 for
the main designs and 80-89 for the sensitivity analysis), with their own random streams; no other seed was
generated or read. Protocol of C7 and 2F, hierarchical bootstrap with 1000 replicates, paired differences,
recall at cap 0.25 with 95 CI. Computation: 2 h 42 min on 2 cores. Report: `results/explorations/REPORT_C8.txt`;
in the proposal, Appendix Y.25.

**A fact about the generator.** `gen_realistic` draws the lognormal amplitude factors (sigma 0.5, one per
station and phase) at every call, that is, in every window and in every catalogue event. They are not station
terms, and in the world in which 2F was confirmed a catalogue has nothing to learn from them. Two worlds were
therefore run: IID, the generator unchanged, and PERS, with the same factors fixed per network and shared by
catalogue and windows, which is the physical case of a site term.

| Contrast, dense windows, sigma 0.5 | IID (600 windows) | PERS (480 windows) |
|:--|:--|:--|
| MT-FW (cap 6) − BOMP-PK | +0.042 [+0.004, +0.083] | +0.061 [+0.008, +0.098] |
| SITE-CAT (catalogue terms, no cap) − MT-FW | −0.066 [−0.108, −0.028] | −0.040 [−0.078, +0.007] |
| SITE-CAT-CAP − MT-FW | −0.002 [−0.026, +0.022] | +0.006 [−0.015, +0.043] |
| ORACLE (true factors, no cap) − MT-FW | −0.011 [−0.056, +0.044] | −0.008 [−0.049, +0.055] |
| ORACLE-CAP − MT-FW | +0.045 [+0.009, +0.074] | +0.032 [+0.003, +0.073] |
| HUBER-N − MT-FW | +0.022 [−0.014, +0.055] | +0.005 [−0.025, +0.048] |
| HUBER-N − BOMP-PK | +0.064 [+0.028, +0.105] | +0.066 [+0.024, +0.105] |
| Fraction of the margin up to CEIL-GD: MT-FW / HUBER-N / ORACLE-CAP | 0.14 / 0.21 / 0.29 | 0.15 / 0.16 / 0.23 |

**Reading.** The amplitude terms estimated from the catalogue do not replace the cap, and not even the true
factors do: the oracle without cap barely matches MT-FW with sigma 0.5 and, with the cap in place, knowing the
factors adds +0.045 and +0.032. With sigma 0.25 the oracle gains nothing over BOMP-PK (−0.001) while the cap
still gains +0.043 [−0.026, +0.099], so the cap is not merely a remedy for the amplitude factors. The only
principled variant that matches the cap is HUBER-N, a Huber fit on the residual norm of each station-and-phase
block with the threshold fixed by rule (square root of the 0.95 quantile of chi-square with 3 degrees of
freedom, 2.7955 noise standard deviations): it matches cap 6 with sigma up to 0.5, is worse with sigma 0.8
(−0.043 [−0.113, −0.001]) and costs 9 times more. In standard windows of the IID world it differs from BOMP-PK
by −0.014 (90 CI within ±0.05), where the cap has −0.041 [−0.086, +0.029]. The post hoc threshold sweep shows a
plateau from 1.5 to 3.4 noise standard deviations, half of the gain at 5 and none at 8. The other variants
(adaptive scale, per-block gains with or without a prior, trace normalisation) fall below MT-FW. Nothing
captures more than a quarter of the margin up to the realistic ceiling: about three quarters remain
unexplained.

**Selection.** HUBER-N was not the variant named in advance as primary (the adaptive-scale one was); it is the
best of 11 candidates in the IID dense design and its estimate there carries selection bias. The threshold
sweep is post hoc with respect to the results, and there is no multiplicity correction.

**Consequence.** A certificate with per-station site terms has no support and is not pursued. The supported
direction is the robust bounded-influence certificate with a level derived from the noise distribution, and
the open problem is the margin that no amplitude treatment explains.

---

## TASK 2I — Huber certificate with noise-derived threshold (G1-R) — EXECUTED AND VERIFIED ON SYNTHETICS; SECONDARY ARMS ANALYSED

**Why.** To confirm on never-used seeds (5151 and 6262) whether HUBER-N, without a calibrated cap, outperforms
BOMP-PK+SC in dense windows, is not inferior to MT-FW+SC with cap 6 there, and does not change the recall of
isolated events. The hypothesis came out of exploration C8, where HUBER-N was the best of 11 candidates; the
registration declares that selection and plans with an effect of +0.04 because of it. In the proposal,
Appendix Y.26.

**Rules.** Written in `PREREGISTRO_1j_1k.txt` (Task 2I section, 2 October 2026, 18:00:05 UTC) before running
any code on those seeds. SHA-256 a983beac…6f97, delivered to the user through the chat before any run
(registration outside the container, neither public nor third-party; the verification checks from disk the
hash and its position in the registration, not the moment it was sent).

| Rule | Contrast | Met if |
|:--|:--|:--|
| 2I-a (superiority) | HUBER-N+SC − BOMP-PK+SC at cap 0.25, dense windows | one-sided 97.5 lower bound > 0 and point ≥ 0.03, in both seeds |
| 2I-b (non-inferiority) | HUBER-N+SC − MT-FW+SC (cap 6) at cap 0.25, dense windows | one-sided 97.5 lower bound > −0.03, in both seeds |
| 2I-c (equivalence) | HUBER-N+SC − BOMP-PK+SC at cap 0.25, standard windows | 90 CI within [−0.05, +0.05] (TOST), in both seeds |

The claim «a noise-scaled Huber fit replaces the calibrated cap» requires a, b and c; if a is met and b is not,
HUBER-N helps but the cap is better; if a fails, the pilot was a selection artefact.

**Method.** HUBER-N uses the MT-FW engine and changes the statistic: Huber function of the residual norm of
each 3-component vector per station and phase, with threshold c = square root of the 0.95 quantile of
chi-square with 3 degrees of freedom = 2.7955 standard deviations of the matched-filter noise, fixed by that
rule before any run and not tuned; 4 IRLS passes starting from the least-squares fit; statistic = Huber loss of
the data minus that of the residual; deflation with the robust fit; no cap on the data.

**Design.** Seeds 5151 and 6262, each an independent replicate; networks j = 0..119 (120 per cell and regime);
5 dense windows (3 to 5 events, kmax 7, TAG 5201) and 5 standard ones (0 to 2 events, kmax 5, TAG 5202) per
network: 3600 dense windows per seed (14477 and 14407 events) and 3600 standard ones (2938 and 2891 events).
Methods: BOMP-PK+SC, MT-FW+SC with cap 6, MT-FW-NOCAP+SC and HUBER-N+SC, on the same windows. Hierarchical
bootstrap with 2000 replicates (strata regime × cell, networks and windows), operating point per method
re-chosen in each replicate, pooling regimes and cells, per seed.

**Freeze and interruptions.** Hashes of 27 files (analysis included) in `task2i_code_hashes.txt`; freeze at
18:14:44 UTC; first work on 5151/6262 at 18:14:45 UTC; parts between 18:14:58 and 22:04:25 UTC. The run was
interrupted twice: at 20:33:25 UTC (resumed at 20:34:09) and at 20:46:15 UTC, when the session was suspended
(resumed at 21:29:10). Each resumption used the same command and continued from the parts already written, one
per network, without recomputing any; the three sessions wrote 1101, 83 and 256 of the 1440 parts;
`sha256sum -c` correct, code unchanged (`task2i_freeze_time.txt`). Result added to the registration at
22:08:07 UTC and verification at 22:16:47 UTC, each with a hash.

**Controls** (`task2i_control.txt`, seed 999; all correct): the C8 pilot is reproduced exactly (90/90 dense
windows and 72/72 standard ones, all four methods); the general engine with least squares gives the hits and
positions of the frozen 2F MT-FW in 60/60 windows with and without cap; the analysis FROC equals that of the
exploration and the published pilot values; the analysis bootstrap on the pilot gives +0.064 [+0.027, +0.104]
versus +0.064 [+0.028, +0.105] published. Estimate of 5.27 h on 2 cores, below the registered budget of 8 h,
so the secondary arms are not omitted.

**Results** (`analyze_2i_primary_out.txt`, cap 0.25, pooling regimes and cells).

| Rule | 5151 | 6262 | Met |
|:--|:--|:--|:--|
| 2I-a, HUBER-N − BOMP-PK, dense | +0.058 [+0.039, +0.076] | +0.053 [+0.037, +0.070] | yes, both seeds |
| 2I-b, HUBER-N − MT-FW, dense | −0.008 [−0.023, +0.006] | −0.002 [−0.018, +0.012] | yes (lower bounds above −0.03) |
| 2I-c, HUBER-N − BOMP-PK, standard, 90 CI | +0.0099 [−0.0100, +0.0272] | +0.0073 [−0.0132, +0.0235] | yes |

Recall at cap 0.25 in dense windows: BOMP-PK 0.177 / 0.198, MT-FW with cap 0.243 / 0.253, MT-FW without cap
0.175 / 0.195, HUBER-N 0.235 / 0.251. In standard windows: 0.479 / 0.488, 0.469 / 0.485, 0.480 / 0.488 and
0.489 / 0.496. On the same windows the cap descriptively repeats 2F: MT-FW − BOMP-PK +0.066 [+0.046, +0.086]
and +0.055 [+0.038, +0.074]; MT-FW without cap − BOMP-PK −0.002 and −0.003. Per cell, HUBER-N − BOMP-PK is
positive in 12 of 12, the smallest with 4 stations (+0.010 [−0.011, +0.057] with strong clutter in 5151). At
cap 0.1: HUBER-N − BOMP-PK +0.063 and +0.060; HUBER-N − MT-FW −0.018 [−0.034, −0.002] in 5151 and −0.009
[−0.031, +0.009] in 6262.

**Secondary arms** (`analyze_2i_out.txt`, frozen code unchanged; registration paragraph of 3 October,
05:54:02 UTC). The thresholds 2.5003, 3.3682 and 5 and the site deviations 0.25 and 0.8, with 20 networks per
cell and regime (600 dense and 600 standard windows per seed), are descriptive, have no rule and are not
independently verified. The three primary rules do not change.

| Arm, dense windows, cap 0.25 | 5151 | 6262 |
|:--|:--|:--|
| Threshold 2.5003 − registered threshold 2.7955 | +0.003 [−0.012, +0.021] | +0.001 [−0.008, +0.016] |
| Threshold 3.3682 − registered threshold | +0.005 [−0.015, +0.027] | +0.006 [−0.012, +0.020] |
| Threshold 5 − registered threshold | −0.013 [−0.039, +0.017] | −0.023 [−0.049, +0.009] |
| Threshold 5 − MT-FW (cap 6) | −0.042 [−0.076, −0.002] | −0.029 [−0.070, +0.013] |
| sd 0.25: MT-FW − BOMP-PK | +0.005 [−0.038, +0.041] | +0.009 [−0.024, +0.053] |
| sd 0.25: HUBER-N − BOMP-PK | +0.036 [+0.004, +0.070] | +0.002 [−0.035, +0.044] |
| sd 0.8: MT-FW − BOMP-PK | +0.155 [+0.096, +0.207] | +0.141 [+0.095, +0.189] |
| sd 0.8: HUBER-N − BOMP-PK | +0.100 [+0.054, +0.143] | +0.072 [+0.031, +0.115] |
| sd 0.8: HUBER-N − MT-FW | −0.055 [−0.113, +0.002] | −0.069 [−0.121, −0.019] |

Between 2.5 and 3.4 noise standard deviations the Huber level does not matter (in standard windows the
differences are +0.000 and +0.008 in both seeds) and with 5 the gain is smaller. With sd 0.25 the cap
contributes nothing and HUBER-N contributes little and only in one seed. With sd 0.8 HUBER-N beats BOMP-PK but
falls below the cap; in 6262 the interval excludes zero and in both seeds the lower bound lies far from the
−0.03 margin, so non-inferiority with respect to the cap does not extend to that dispersion (with 20 networks
the arm would not have shown it with sd 0.5 either, where it gives −0.029 [−0.070, +0.013] and −0.007
[−0.049, +0.030]). In standard windows no contrast between HUBER-N and BOMP-PK excludes zero with sd 0.25 or
0.8. Per regime and with the strict hit rule (origin within 0.5 s), HUBER-N − BOMP-PK in dense windows: +0.040
and +0.020 with weak clutter, +0.022 and +0.017 with strong clutter, all intervals above zero.

**Independent verification** (`verif_2i/task2i_verification.txt`): positive confirmed for the registered
claim. Provenance (hashes of 27 files with no failures; chain of 11 registration prefixes, a983beac ends exactly
at the end of the Task 2I section; no part earlier than the freeze; exactly two gaps of more than 60 s, of 82 s
and 2599 s, which coincide with the two interruptions; no prior data from 5151/6262); symmetry (same windows,
grid, origins, corrections, number of atoms, exclusion and hit function for all four methods); statistics
re-derived with independent loader, FROC and bootstrap (for example +0.0582 [+0.0402, +0.0753] versus +0.058
[+0.039, +0.076]); six parts regenerated from scratch, two from each of the three sessions, identical to the
stored ones. Scope established by the verification:

| Point | Result |
|:--|:--|
| Hit rule with origin time (≤ 5 km and ≤ 0.5 s) | the gain over BOMP-PK drops to +0.031 [+0.022, +0.040] and +0.019 [+0.008, +0.030]; the point criterion of 0.03 would fail in 6262; only ~57 % of the location-based hits have origin error ≤ 0.5 s |
| Equivalence per regime | met when pooled, not in standard windows with weak clutter of 6262: −0.032 [−0.063, −0.004], 90 CI [−0.060, −0.007] |
| Recall with all atoms accepted | does not rise (0.623 BOMP-PK and 0.622 HUBER-N in 5151; 0.628 and 0.616 in 6262): the gain is one of ranking and atom choice, half reliability of the score across windows (+0.027 / +0.029) and half atom choice (+0.030 / +0.018), not of more events found |
| Cap 0.1 | HUBER-N falls somewhat below the cap in 5151 (−0.018) |
| Cost | 0.927 s per dense window versus 0.103-0.104 s for MT-FW and 0.587-0.602 s for BOMP-PK: about 9 times MT-FW |
| Absolute recall | low: 0.177-0.198 → 0.235-0.251 in dense windows at cap 0.25; ~0.10 with the 0.5 s rule |

**Consequence for the thesis.** G1-R is supported by two preregistered synthetic tests, with the registration
hash delivered outside the container before running and with independent verification: 2F with a calibrated
cap and 2I with a Huber level fixed by rule from the noise, which removes the dependence on a tuned constant.
What holds, stated modestly: a bounded-influence moment-tensor subspace certificate improves recall at a fixed
false-event rate in dense sequences of overlapping events relative to BOMP-PK with the same timing tolerance,
without harming isolated events when the regimes are pooled. The Huber fit is standard robust statistics: the
novelty is the measured design result, not an algorithm. MT-FW with cap 6 remains the cheap reference; HUBER-N
is the form whose level can be fixed without calibration data.

**Limitations.** Synthetic and from a single generator family (the two seeds are independent draws from that
family); modest size; the gain shrinks under the strict hit rule; cost of 9 times; operating points chosen on
the evaluation, per method; descriptive secondary arms, with 20 networks per cell and regime and without
independent verification; with site deviation 0.8 the Huber fit falls below the cap.

**What remains.** A second generator, a frozen operating point and real data. The margin of 0.08-0.24 below
the realistic ceiling remains open: about three quarters remain unexplained by amplitude treatments
(exploration C8), and 2I does not change that figure; explorations C9 and C10 say where it is (next section).

---

## EXPLORATIONS C9 and C10 — Evidence without travel times and the timing budget — EXPLORATORY

**Status.** Exploratory, not confirmatory: seed 999 only, on networks of that seed not used before (130-149 for
C9 and 160-179 for C10; networks 199 and 195 for validation and for the prior choices), with their own random
streams; no other seed was generated or read. Protocol of C7, C8 and 2F; 600 windows per design (20 networks ×
5 windows × 3 cells × 2 regimes); hierarchical bootstrap with 1000 replicates; paired differences at cap 0.25
with 95 CI; no multiplicity correction. Choices made after seeing the results, declared: in C9 the reordering
exponents, the SELPOL and PENPOL variants and the worlds with 20 and 45 degrees of polarisation dispersion
(SELPOL is the best of about 50 candidates); in C10 the true statics at catalogue level, the first-iteration
diagnostic and the SC-R corrections, whose rule was written after seeing the oracle and was checked against the
true statics. Reports: `results/explorations/REPORT_C9.txt`, `REPORT_C10.txt` and `LIT_NO_TRAVELTIME.md`; in the
proposal, Appendix Y.27.

**The question.** About one in five true events is left without a hit, with a false atom at 5-10 km. C9 asked
whether evidence that does not use travel times (P and S polarisation, amplitudes, S/P ratio) can relocate
those atoms; C10, whether the hard ±0.25 s window wastes timing information that a Gaussian prior would
recover.

**C9, Bayesian Cramér-Rao bound on the epicentre** (median over 200 sources per cell, km; the prior alone gives
12.25 km).

| Channel | 4 stations | 6 stations | 8 one-sided |
|:--|:--|:--|:--|
| Arrival times, 0.144 s (±0.25 s tolerance) | 0.68 | 0.43 | 0.69 |
| Arrival times, 0.30 s (uncorrected statics) | 1.34 | 0.84 | 1.29 |
| P polarisation, 10 degrees (nominal) | 9.23 | 8.53 | 7.79 |
| P and S polarisation, 10 degrees | 7.47 | 6.10 | 6.03 |
| Amplitudes, site deviation 0.5 | 9.18 | 7.75 | 9.11 |
| S/P ratio, site deviation 0.5 | 11.77 | 10.46 | 10.88 |
| All channels without times, nominal | 5.59 | 4.35 | 4.81 |
| Times at 0.144 s together with all channels without times | 0.65 | 0.41 | 0.66 |

Times locate to 0.4-0.7 km and all channels without times together to 4.4-5.6 km; a joint estimator is bounded
to an improvement of 4-5 %. The limit is neither the 10 degrees nor the site factors but the waveform noise:
the P arrival measures 1.5-1.8 noise standard deviations at these distances, that is, 32-39 degrees of angular
noise per station. Physics without times would carry weight only with about 10 times the signal-to-noise ratio
and a polarisation dispersion of 2 degrees or less.

**C9, locators without times on 600 isolated events** (vectors read at the true arrival samples): median
epicentral error of 12.7-15.2 km with the P direction, 8.2-9.8 km with P along the ray and S transverse,
7.7-8.4 km with the projection onto the moment-tensor subspace with free times and 12.5-19.4 km with 1/r
amplitudes, versus 0.5-0.7 km with the times and 11.8 km for guessing the centre. At most one in three events
falls within 5 km. Without travel times the arrivals are not even identified: a blind reader finds the P peak
at 19-20 % of the stations.

**C9, relocation.** Moving each atom to the maximum of a statistic without times within 12 km loses recall in
all variants: −0.11 to −0.18 in dense windows and −0.23 to −0.48 in standard ones. The oracles do not rescue
it: with the true amplitude factors −0.116 [−0.158, −0.081]; in a world without polarisation deviation or
factors −0.138 [−0.187, −0.089]. Relocation keeps 0.43-0.57 of the hits and converts 0.22-0.46 of the near
hits in dense windows, which is what a random move within a 9 km cluster does.

**C9, polarisation admissibility as a weight.** The fraction of energy admissible by polarisation (P along the
ray, S transverse) separates hits from false atoms (AUC 0.68 in dense windows and 0.78 in standard ones). As a
factor of the MT-FW selection statistic (SELPOL, post hoc, the best of about 50): +0.046 [+0.007, +0.070] over
MT-FW and +0.045 [+0.003, +0.086] over BOMP-PK in dense windows; +0.002 [−0.041, +0.064] over BOMP-PK in
standard ones; nothing with 4 stations; 0.12 [0.01, 0.24] of the margin. The reordering gain is +0.028 at 10
degrees, +0.003 [−0.011, +0.023] at 20 degrees and −0.015 [−0.025, +0.002] at 45 degrees. It is a small and
conditional clutter discriminator, not location physics; with reliable polarisation at one in five stations of
the real network, the expectation there is null.

**C10, soft timing priors.** Placing a Gaussian likelihood on the peak shifts, with the level fixed by rule
from the catalogue (median 0.138 s with weak clutter and 0.224 s with strong clutter), gives no gain under any
hit rule, window type or detector: relative to MT-FW, −0.028 to −0.002 in dense windows and −0.066 to −0.006 in
standard ones. Widening the free window to ±0.5 s loses (−0.047 [−0.087, −0.008] in dense windows): the
±0.25 s cut does useful work.

**C10, why: the timing budget** (shift of a true arrival relative to the prediction of its nearest node, s
rms).

| Term | Size |
|:--|:--|
| True model error after a perfect station term | 0.05 (P), 0.07 (S) |
| Error of the +SC corrections relative to the true statics | 0.12-0.17 (weak clutter), 0.18-0.24 (strong) |
| Position within the 3 km cell and 0.1 s origin step | 0.18 (P), 0.30 (S) |
| The same on a 1 km grid | 0.10 (P), 0.16 (S) |

The budget of a true hypothesis is 0.2-0.4 s per arrival, not the 0.144 s of the bound, and a prior with that
honest width is almost flat within ±0.25 s.

**C10, oracle decomposition** (difference to MT-FW, cap 0.25, registered rule; dense / standard).

| Variant | Dense | Standard |
|:--|:--|:--|
| True statics, hard window | +0.047 [−0.003, +0.082] | +0.098 [+0.032, +0.141] |
| True statics with prior at the level of the model error | +0.081 [+0.035, +0.117] | +0.146 [+0.065, +0.217] |
| True statics with prior at catalogue level, relative to the hard window | −0.008 [−0.039, +0.032] | −0.059 [−0.108, +0.027] |

In standard windows the oracle detector reaches the realistic ceiling (0.592 versus 0.603 for CEIL-GD); in
dense windows it captures 0.31 [0.21, 0.43] of the margin. First-iteration diagnostic (post hoc, 235 isolated
events): with the registered statistic, 0.485 of the events win their comparison; the prior at catalogue level
leaves it at 0.489; the oracle raises it to 0.609-0.634, and the remaining 37-39 % is lost to amplitude (weak
arrivals, clutter, coda, site factors), not to timing.

**C10, local refinement T3 and SC-R corrections.** T3 (refining each atom on a 1 km grid with the joint
objective) loses no recall and improves location and origin time of isolated events (median error from 2.51
to 2.07 km with BOMP-PK and from 2.42 to 2.21 km with MT-FW). Under the registered rule its gain does not
exclude zero (+0.059 [−0.009, +0.112] and +0.011 in standard windows); it only helps under the rule that also
requires the origin time within 0.5 s, and only in standard windows (+0.100 [+0.031, +0.168] with BOMP-PK and
+0.066 [+0.019, +0.136] with MT-FW), which is not the registered endpoint and is the best of about 70
contrasts. SC-R, post hoc, lowers the error of the corrections from 0.145 to 0.110 s with weak clutter and
from 0.214 to 0.202 s with strong clutter; on its own it changes nothing, and with T3 it gives in standard
windows +0.064 [−0.004, +0.123] under the registered rule and +0.123 [+0.058, +0.201] with origin time.

**Reading.** The margin splits into two parts. One is calibration: the accuracy of the station terms; with the
statics known at the level of the model error the detector reaches the ceiling on isolated events. The other is
amplitude: weak arrivals, clutter, coda and amplitude factors. Physics without travel times is limited by
information in this benchmark, not by algorithm. The open problem becomes calibration (station terms
self-consistent with the detections) and amplitude evidence, not another location principle.

**Limitations.** One seed, 20 networks per cell and regime, oracles that a real network does not have, post
hoc analyses and no multiplicity correction: these are pilot estimates. The statics deviation of 0.3 s, the
20-event catalogue and the 3 km grid are design choices, and the size of the correction error depends on them.
Soft timing priors, polarisation filters and joint hypocentre determination with station terms are known
practice: none of this is a method.

**Literature map** (`LIT_NO_TRAVELTIME.md`, a sweep and not a systematic review). Of eight families, those that
truly avoid travel times are polarisation, amplitudes and their ratios, envelope energy and, with instruments
that a sparse network does not have, gradiometry and rotation; those that hide them are time reversal and
wave-equation imaging (in the propagator), matched field and back-projection (in the phase of the Green's
function), envelope correlation and end-to-end learned location. Of its 72 references, 26 were read only at
metadata level and several publisher pages could not be opened.

**Next step.** Station-term calibration self-consistent with the detections, registered with the true-statics
oracle as a bound and with the hit rule with origin time as the endpoint; then, amplitude evidence for weak
arrivals. The lines of relocation by polarisation or amplitude and of soft timing priors with catalogue
corrections are closed as negatives with their mechanism. The two preregistration drafts in the reports (SELPOL
and T3) name the same unused seeds, 7171 and 9393; neither was registered or run, and only one can use them;
Task 2L, which absorbs the protocol of the T3 draft, used them.

---

## TASK 2L — Pick-free station terms within the refined detector (G1-C) — EXECUTED AND VERIFIED ON SYNTHETICS

**Why.** To confirm on never-used seeds (7171 and 9393) whether the S2 station terms, estimated without picks
by a joint determination of the hypocentres of the 20 catalogue events and of the terms from the energy peaks
of the matched filter, improve, relative to the registered S0 corrections, the recall with origin time of the
refined detector (MT-FW with cap 6 and T3 refinement), with no loss under the epicentre-only rule. The
hypothesis came out of exploration C11 on seed 999 (networks 200-239), with the plan frozen with a hash before
any detection result (`explore_contrib/C11_PLAN.txt`); its primary contrast was +0.070 [+0.010, +0.120], and
the registration plans with +0.05 because the C10 effects shrank on replication. In the proposal, hypothesis
G1-C and Appendix Y.28.

**Rules.** Written in `PREREGISTRO_1j_1k.txt` (Task 2L section, 4 October 2026, 21:10:33 UTC) before running
any code on those seeds. SHA-256 13a6e9b1…a62a, delivered to the user through the chat before any run
(registration outside the container, neither public nor third-party). The contrast of all rules is MT-FW with
S2 and T3 minus MT-FW with S0 and T3, at cap 0.25, with the hit rule that requires the epicentre within 5 km
and the origin error within 0.5 s at most; each rule must be met in both seeds.

| Rule | Contrast | Met if |
|:--|:--|:--|
| 2L-a (superiority) | standard windows, with origin time, pooled | one-sided 97.5 lower bound > 0 and point ≥ +0.03 |
| 2L-b (superiority) | dense windows, with origin time, pooled | one-sided 97.5 lower bound > 0 |
| 2L-c (non-inferiority per regime) | standard and dense, with origin time, per regime | one-sided 95 lower bound > −0.02 |
| 2L-d (non-inferiority under the previous rule) | standard and dense, epicentre only, pooled | one-sided 95 lower bound > −0.02 |
| Secondary | BOMP-PK with S2 and T3 − BOMP-PK with S0 and T3, standard, with origin time, 60 networks | one-sided 97.5 lower bound > 0 |

**Design.** 160 networks per cell and regime with 5 standard windows (4800 per seed, with 3826 and 3863
events) and 100 with 5 dense windows (3000 per seed, with 11978 and 12000 events); station delays with a
deviation of 0.3 s; 20-event catalogue with errors of 1 km and 0.1 s; the generator of the previous tasks.
BOMP-PK runs on the standard windows of 60 networks and the true-delay oracle with T3 on 30 networks per cell
and regime, for budget reasons. Hierarchical bootstrap with 2000 replicates per seed.

**Freeze and run.** Hashes of 37 files at 21:24:48 UTC; first work at 21:24:53 UTC; 1920 parts between
21:25:19 UTC on 4 October and 01:48:15 UTC on the 5th; the run was interrupted at 22:47:19 UTC when the session
was suspended and was resumed from the parts with the code unchanged (`sha256sum -c`, 37 of 37). The budget
rule was not triggered (estimate of 4.62 h). The controls on seed 999 are correct (`task2l_control.txt`): the
C11 pilot is reproduced exactly in 18 parts, S0 equals the registered correction in 18 of 18 and the analysis
bootstrap on the pilot gives +0.070 [+0.011, +0.120] versus +0.070 [+0.010, +0.120] published. The result
paragraph entered the registration at 02:00:25 UTC and the verification paragraph at 02:36:52 UTC, each with
its hash.

**Results** (`results/task2l/analyze_2l_out.txt`). With origin time, recall at 0.25 from S0 and T3 to S2 and
T3: 0.346 to 0.403 (7171) and 0.350 to 0.411 (9393) in standard windows; 0.124 to 0.162 and 0.120 to 0.151 in
dense ones.

| Rule | 7171 | 9393 | Met |
|:--|:--|:--|:--|
| 2L-a, standard | +0.057 [+0.034, +0.081] | +0.062 [+0.036, +0.083] | yes |
| 2L-b, dense | +0.038 [+0.024, +0.054] | +0.030 [+0.015, +0.044] | yes |
| 2L-c, lower bounds per regime and window type | +0.0069 to +0.0517 | +0.0102 to +0.0461 | yes |
| 2L-d, epicentre only, standard / dense | +0.036 / +0.020 | +0.050 / +0.035 | yes |
| Secondary, BOMP-PK | +0.068 [+0.026, +0.107] | +0.071 [+0.029, +0.114] | yes |

**Independent verification** (`results/task2l/verif_2l/task2l_verification.txt`): positive confirmed. No
leakage: S2 recomputed in new processes without access to the true delays or the true catalogue is identical
bit for bit in 300 of 300 units, and the catalogue the estimators receive carries the registered errors (0.998
km and 0.099 s). The gain is not an origin-time artefact: the origin bias is the same in both arms, S2 narrows
the origin errors (hits at 0.3 s from 0.47 to 0.53 in standard windows) and the contrast is positive with hit
windows from 0.3 to 1.0 s and with each method's median origin subtracted. It is carried by the terms, not by
the change in the T3 level. Statistics reproduced with independent code (+0.0572 and +0.0616 in 2L-a), 37 of 37
hashes, six regenerated parts identical and analysis re-run with byte-for-byte identical output. Scope
established by the verification:

| Point | Result |
|:--|:--|
| Where the gain is | essentially in the one-sided 8-station cell: +0.108 [+0.057, +0.151] and +0.113 [+0.067, +0.157] in standard windows; with 6 stations +0.048 and +0.027; with 4 stations +0.021 and +0.018, with no demonstrated gain |
| S2 worse than S0 | in about 22 % of the networks (3-4 % of the one-sided ones with weak clutter, 19-27 % of the one-sided ones with strong clutter, 12-28 % with 6 stations, 26-37 % with 4); in the 6-station networks of that group, −0.019 and −0.049 in standard windows |
| Conditional on T3 | there is no arm with S2 and the hard window; in the pilot the terms without T3 gave nothing measurable |
| Catalogue | assumes errors of known size (1 km and 0.1 s); with a 2 km prior the error of the terms returned to the S0 level in the development of C11 |
| Oracle | S2 captures 0.56 and 0.61 of the gain of the true delays in standard windows and 0.55 and 0.28 in dense ones; the oracle remains +0.052 and +0.060 ahead in standard windows |
| Absolute recall | low: 0.403 and 0.411 in standard windows and 0.162 and 0.151 in dense ones, with origin time |

**Consequence for the thesis.** G1-C is supported by a preregistered synthetic test, with the registration hash
delivered before running, the code frozen and independent verification: within the refined detector, the
station terms estimated without picks from the 20-event catalogue raise the recall with origin time. It is a
modest calibration result: the joint determination of hypocentres and station terms is classical in pick-based
location, and what is measured is its pick-free estimation and its effect within a pick-free detector, not an
algorithm. From 2L on, the hit rule with origin time is the primary one.

**Limitations.** Synthetic and from a single generator family; concentrated in one geometry; conditional on T3;
assumes a catalogue with known errors; operating points chosen on the evaluation, per method and equal for both
arms; G1-R and G1-C have not been run together.

**What remains.** The sensitivity of S2 to catalogues with errors of unknown size, biased, or with fewer events;
a criterion that says before detection in which networks the estimate is reliable; the combination with the
bounded-influence statistic in a registration; and, on real data, the residual per-station terms without an
oracle. Explorations C12 to C14 (next section) add the origin time of well-located detections as a prior
target.

---

## EXPLORATIONS C12, C13 and C14 — Double couple, decoys and self-interference — EXPLORATORY

**Status.** Exploratory, not confirmatory: seed 999 only, on networks of that seed not used before (1200-1229
and 1290 for C12; 1400-1411 for C13; 1600-1602 for development and 1610-1621 for evaluation for C14), with
their own random streams, after 2L. The detector is the registered one of the previous tasks (MT-FW with cap
6, S0 corrections, hard ±0.25 s window), not the refined one of 2L; T3 appears only as a reference in C14. Hit
rule with origin time as primary and epicentre-only as secondary, exact FROC, hierarchical bootstrap, no
multiplicity correction. Reports `results/explorations/REPORT_C12.txt`, `REPORT_C13.txt` and `REPORT_C14.txt`; in
the proposal, Appendix Y.29. The incidents of C13 and C14 are in E43 to E45 of the manifest.

**C12, the statistic restricted to double couples.** Question: whether restricting the subspace statistic to
the double-couple manifold lowers the threshold at a fixed false-event rate enough to detect weak events, and
whether its null law can be derived. 16 networks per cell and regime, 1920 windows per design (1550 events in
standard windows and 7613 in dense ones), plus 1440 event-free windows with clutter and 600 noise-only ones.
The null law of a candidate is the tube formula with two constants per node, the volume and the curvature of
the manifold, with 4.4 to 4.5 effective degrees of freedom instead of 6; the simulated 0.01 and 0.001
quantiles, 14.28 and 19.69, agree with the theoretical ones, 14.17 and 19.46. In the full detector without
events the double-couple threshold is 0.906 [0.883, 0.928] of the full one at 0.1 and 0.898 [0.886, 0.914] at
0.25, with noise only. Recall does not change:

| Statistic | standard, with origin time | standard, epicentre only | dense, with origin time | dense, epicentre only |
|:--|:--|:--|:--|:--|
| Full (registered), recall at 0.25 | 0.317 | 0.470 | 0.117 | 0.268 |
| Deviatoric − full | −0.005 [−0.027, +0.020] | −0.004 [−0.022, +0.026] | −0.006 [−0.016, +0.009] | −0.001 [−0.021, +0.016] |
| Double couple − full | +0.005 [−0.021, +0.033] | +0.017 [−0.010, +0.043] | −0.004 [−0.018, +0.011] | +0.006 [−0.016, +0.023] |
| Mixed (atoms of the full one scored with double couple) − full | −0.003 [−0.020, +0.018] | −0.008 [−0.030, +0.018] | +0.004 [−0.007, +0.012] | +0.002 [−0.013, +0.015] |

Why: in the score band where the decision is made (25 to 50), hits lose 5.9 units on moving to the double
couple and false detections 6.6, while the threshold drops 5.1 (at 0.1) and 4.3 (at 0.25). An exact double
couple loses 1 unit (0.95 measured); the generator's site factors and polarisation dispersion explain 1.7 to
2.5 of the observed 5.9 to 6.4, and the rest was not decomposed. It costs 0.199 versus 0.135 s per standard
window. Verdict: null for detection, positive for the analytical piece (a double-couple mechanism test with
calibrated p-values, and a bound of about 5 % of amplitude for any restriction of the source subspace).

**C13, threshold calibration with decoys.** Question: whether a second pass of the detector over hypotheses
that cannot be a real source estimates, without labels, the false detections above a threshold, in order to
set an operating point or control the false discovery rate (FDR). 12 networks per cell and regime: 1440
standard windows (1178 events), 1440 dense ones (5747), 576 clutter-only ones and 144 noise-only ones. Decoys:
moveout with permuted stations (D1p), moveout with offsets (D1o), sign-inverted templates (D2s), scrambled
polarisation (D2r), circularly shifted traces (D3) and, post hoc, D3 with restricted delays (D3c). An
event-free twin of each window to separate what events generate from what clutter generates.

| Result | Value |
|:--|:--|
| Exchangeability in event-free windows (D1p, D2s, D2r, D3c) | permutation p between 0.04 and 1.0; one of 48 below 0.05; D3 and D1o fail with strong clutter (E44) |
| False detections per window with events relative to the event-free twin, above t95, standard / dense | 12 [9, 18] / 45 [33, 63] times |
| Distant false detections after a hit and deflated event, relative to the twin, t80 / t95 | 1.9 [1.6, 2.4] / 5.0 [3.1, 8.8] times |
| D1p with batches of 120, realised FDR at q = 0.05 / 0.1 / 0.2, epicentre only | 0.02-0.05 / 0.07-0.10 / 0.14-0.16 |
| The same with origin time | 0.08-0.19 / 0.18-0.29 / 0.28-0.37 |
| D1p, realised false detections per window for c = 0.25, epicentre only, standard / dense | 0.21 [0.17, 0.25] / 0.23 [0.19, 0.28]; recall 0.453 / 0.245 versus 0.436 / 0.253 for the labelled point |
| Proportion of false detections among those accepted at the benchmark point (0.25), standard / dense | 0.41 / 0.20 epicentre only; 0.52 / 0.37 with origin time |

Under the epicentre-only rule D1p reproduces without labels the benchmark operating point and keeps the FDR
within q = 0.2 in the three cells, although the one-sided cell lies above the nominal level at q ≤ 0.1. Under
the origin-time rule it fails, because a real event located within 5 km with the origin shifted by more than
0.5 s counts as false and no moveout decoy imitates it (0.2 to 1.2 such atoms per window with events); only D2s
controls the FDR, at the cost of almost all the recall. Verdict: promising with conditions and partial; no
guarantee, not valid with origin time, one seed. The report leaves a draft of a frozen configuration and of
criteria under the epicentre-only rule; it has not been registered or run.

**C14, self-interference of the detector.** Question: whether the false-detection budget of the operating point
is spent by self-interference (deflation residual, coda, cross-phase ghosts, copies) and whether its remedies
raise recall. Diagnosis with ground truth on 3 development networks per cell and regime (360 windows), with
the statistic of each false detection split exactly among the parts of each event, the clutter and the noise,
and the deflations. Remedies without ground truth: exclusivity of claimed arrivals (R1), variable noise scale
after each deflation (R2), adaptive deflation with a short filter (R3), merging of duplicates in origin time
(R4) and R2 with R3; T3 as a reference. Constants chosen in development and frozen with hashes on 9 October at
11:36:11 (−05), before generating the evaluation networks 1610-1621 (12 per cell and regime; 1152 standard
windows with 855 events and 1152 dense ones with 4594).

Diagnosis (development): at the operating point the largest class of false detections with origin time is not
self-interference but the well-located detection of a real event (median of 3.2 and 3.0 km) with the origin
shifted by a median of 0.76 s in standard windows and 0.86 s in dense ones, early in 75 and 73 % of the cases
and with no other correct detection of the same event in 86-92 %: 55-58 % of the false detections in standard
windows and 80-95 % in dense ones; the deflation residual weighs 0-7 % and the coda nothing. With ground truth,
a perfect deflation would give +0.051 [−0.024, +0.111] in standard windows and +0.052 [+0.012, +0.105] in dense
ones. Evaluation on the unseen networks:

| Method | standard, recall or difference | dense, recall or difference | false detections relative to the twin (q95), standard / dense | well-located detections with shifted origin per window, standard / dense | overlapping events at any score |
|:--|:--|:--|:--|:--|:--|
| BASE | 0.283 | 0.108 | 6.23 / 10.04 | 0.082 / 0.161 | 0.400 |
| R1, hard exclusivity | +0.015 [−0.007, +0.044] | −0.005 [−0.017, +0.008] | 3.18 / 5.01 | 0.089 / 0.169 | −0.103 [−0.121, −0.086] |
| R2, variable noise scale | +0.016 [−0.006, +0.043] | −0.002 [−0.015, +0.016] | 3.41 / 5.49 | 0.093 / 0.174 | −0.009 [−0.029, +0.009] |
| R3, adaptive deflation | +0.016 [−0.006, +0.042] | −0.007 [−0.019, +0.006] | 3.17 / 4.84 | 0.090 / 0.174 | −0.098 [−0.116, −0.080] |
| R4, merging within 10 km and 1 s | +0.002 [−0.007, +0.017] | −0.011 [−0.018, −0.001] | 6.20 / 10.21 | 0.084 / 0.157 | −0.046 [−0.060, −0.033] |
| R2 and R3 | +0.016 [−0.007, +0.048] | −0.004 [−0.017, +0.011] | 2.52 / 4.24 | 0.091 / 0.175 | −0.095 [−0.115, −0.078] |
| T3, reference | +0.056 [+0.015, +0.100] | +0.011 [−0.007, +0.033] | not computed | 0.082 / 0.165 | +0.037 [+0.016, +0.057] |

The remedies halve the false detections induced by events, and away from the events they bring them down to
the clutter level, but recall with origin time at 0.25 does not move: the false detections they remove lie
below the threshold of the operating point, and above it the well-located detections with shifted origin
remain (a third of the false detections in standard windows and two thirds in dense ones on these networks).
The two most aggressive remedies erase a quarter of the overlapping events that BASE finds at any score. The
development gains shrank by a factor of 3 to 6 in evaluation. Verdict: negative for recall; no remedy is
proposed for registration.

**Reading.** C12 closes the restriction of the source subspace as a recall lever and leaves an analytical
piece. C13 gives a first label-free operating point, but only under the epicentre-only rule, without guarantee
and on one seed; it is not the deployable operating point that is missing. C14 explains why C13 fails with
origin time and shows where the false-detection budget is spent: on the origin time of well-located detections
of real events, almost always early and almost always unique.

**Limitations.** One seed; 16, 12 and 12 evaluation networks per cell and regime; the registered detector and
not that of 2L; post hoc analyses and choices declared in the reports; no multiplicity correction. The C14
diagnosis is from development (44 and 46 false atoms at the operating point, uncertainty of about 0.07 in the
fractions).

**Next step.** Explore, before registering anything, whether the S2 terms with T3 (the 2L detector), a better
timing model or an explicit treatment of the early origin bias reduce the class of well-located detections with
shifted origin, and whether decoy calibration then holds under the origin-time rule. The self-interference
remedies are closed for recall at the current operating points.

---

## EXPLORATION C15 — Where the origin error of well-located detections comes from — EXPLORATORY

**Status.** Exploratory, seed 999 only, three never-used development networks (1800-1802), 720 windows, no
separate evaluation set (optimistic figures). Report `results/explorations/REPORT_C15.txt`; Appendix Y.29b.

**Result.** The plateau of the argmax in origin time explains nothing: the maximum is unique in 97.8 % of the
atoms. With the true location, the samples that the detector itself read give the origin within 0.5 s or less in
93 % of the cases in standard windows and in 76-77 % in dense ones, with the engine's model or with the true one,
versus 66 % and 58 % for the reported origin: the origin error is a location error. Relocating in a second stage
with those readings (weighted median, robust Gauss-Newton, equal differential time likelihood with and without
re-reading) gives between +0.021 and +0.042 in standard windows and between -0.003 and +0.005 in dense ones
(LOC+T at 0.25 or less), all intervals including zero. Cyclic block re-estimation (NOMP-type feedback) is null
or negative. **Null for recall**; in this generator the detector is close to what the information allows, and
what can change the evidence is leaving the generator.

---

## EXPLORATION C16 — Real noise and real geometry with injected events — EXPLORATORY

**Status.** Exploratory, seed 999 only; 36 new networks per arm (2402-2413), frozen before generating them
(`C16_FREEZE.txt`). Real continuous noise from one day (23 July 2010) of 14 three-component stations of the Alto
Tiberina network (INGV), scaled per channel to the generator's noise level; events, mechanisms, pulses, times
and statics from the generator; the cells are subsets of the real stations. A second arm with the generator's
noise on the same networks and signals. Report `results/explorations/REPORT_C16.txt`; Appendix Y.29c.

**Result.** G1-R survives and is strengthened with real noise: in dense windows, LOC rule, at 0.25 false
detections per window or fewer, MT-FW with cap beats BOMP-PK by +0.259 [+0.131, +0.321] and HUBER-N by +0.294
[+0.178, +0.355], versus +0.069 and +0.113 with the generator's noise; under the LOC+T rule, +0.085 and +0.133.
The mechanism is that of the hypothesis: the transients of the real noise (excess kurtosis with a 90th
percentile of 42) sink BOMP-PK, whose energies have no cap, and not the bounded-influence statistics. G1-C keeps
its sign (+0.049 [-0.083, +0.149]) but is not resolved with this sample. The operating point of the capped
detector transfers to real noise with a small excess (up to 0.06 false detections per window under the primary
rule); that of BOMP-PK does not (0.29 to 0.86).

**Next step.** Registered confirmation of G1-R with real noise on other days of the same network (10, 11 and 12
February and 31 December 2010) and never-used seeds: Task 2M, executed and confirmed (below).

---

## TASK R0 — Pilot with public real data: Ridgecrest 2019 — EXECUTED, NEGATIVE

**Registration.** Written and its hash sent before touching the waveforms (10 October 2026, 00:00 UTC);
amendment and freeze before the evaluation (00:27 UTC). Data: one hour (19:00-20:00 UTC on 4 July 2019) of the
CI network with the SCSN catalogue, from a public repository (Apache-2.0). Calibration and development in the
first 15 minutes; evaluation in the remaining 44 min 20 s (45 events of M ≥ 2.0). Subnetworks N8, N6, N4 and
U8. Primary contrast: MT-FW with cap 6 + S2 + T3 against the better of GaMMA and PyOcto with PhaseNet (with and
without corrections), recall at 30 false events per hour with origin time required.

**Result.** −0.383 [−0.481, −0.310] pooled; per subnetwork N8 −0.578, N6 −0.667, N4 0.000 (no method detects
with 4 stations), U8 −0.289. Recall of the thesis detector 0.022, 0.000, 0.000, 0.000; of PyOcto 0.600, 0.667,
0.000, 0.289. R0-a and R0-b are not met. The secondary arm with a reading of ±0.69 s (fixed by rule from the
measured error) changes nothing (+0.000 [−0.016, +0.019]). Post hoc diagnosis: location (8-18 km median error
versus 1.5 km for PyOcto) and ranking (the score of the real events does not separate in a continuous
sequence). Closure, declared post hoc (16-17 October 2026): M1 and M3 re-run with the frozen code reproduce the
registered recalls (M1 on all subnetworks, M3 on N8); the registered secondary contrast M1 - M3 is +0.000 [-0.023,
+0.021] with origin time and -0.050 [-0.101, -0.013] with the epicentre rule. Appendix Y.30; code and outputs in
`results/taskR0/` (closure in `run_close/`). The complete independent verification, which needs M2 and the picks
re-run, is pending.

**Next step.** Exploration C18 (below) tested a layered and a fitted 1-D medium and non-saturating scores on the
development segment: neither reduces the reading error nor locates the events. The next registered real test waits for
a candidate that passes the gate fixed in C18; G1' on the RSNC remains registered.

---

## TASK 2M — G1-R with independent real noise (Alto Tiberina) — EXECUTED, CONFIRMED

**Registration.** 16 October 2026, 16:32 UTC, hash sent before any computation on seeds 2718 and 3141; code frozen
with SHA-256 (45 files). Author decisions: quality rule qc1 and station rule unchanged; seed 2718 (2010-12-31, noise
with transients, 10 stations) as the only primary replicate; seed 3141 (10 and 11 February 2010, near-Gaussian, 8
stations) as a secondary replicate with a written prediction; generator-noise arm descriptive on 8 networks per cell.
Rules 2M-a (cap minus BOMP-PK) and 2M-b (HUBER-N minus BOMP-PK) on dense windows, epicentre rule, <= 0.25 false
events per window: lower bound > 0 and point >= +0.03, both required in seed 2718.

**Result.** 2M-a +0.228 [+0.148, +0.298] and 2M-b +0.241 [+0.179, +0.293]: both met, Holm rejects both; G1-R confirmed
with real noise with transients (one network, one day). Seed 3141: +0.057 [+0.004, +0.121] and +0.120 [+0.070, +0.168],
both rules also met. With origin time required on 2718: +0.085 and +0.104. HUBER-N is not shown non-inferior to the cap
on 2718 (+0.014 [-0.038, +0.056]) and beats it on 3141 (+0.062 [+0.012, +0.107]). Post hoc: the cap's gain is
concentrated in the hours with transients (night hours of 2718: +0.033 [-0.044, +0.128]; 7.5-15 h: +0.300, with BOMP-PK
at 0.003); HUBER-N keeps +0.10 to +0.15 in quiet noise. Two networks of 3141 failed by exhausting their catalogue reserve
(declared; error E48). Independent verification: confirmed with caveats of scope. Appendix Y.31; `results/task2M/`.

**Next step.** None for G1-R on real noise; the open question is real signals (R0 and its follow-up).

---

## EXPLORATION C18 — What fails on real data after R0 — DONE (exploratory)

On the development segment of the Ridgecrest hour only (the evaluated segment is not read), and on seed 999 for the
robustness check. (a) A layered model (Hadley-Kanamori), a fitted homogeneous medium and a fitted gradient do not
predict the held-out energy-peak readings better than the registered medium: 0.53 to 0.63 of the readings within
0.25 s on every subnetwork and model; S readings are worse than P. (b) At the true hypocentre the capped score of a
real M >= 2 event does not exceed the maxima of event-free windows (AUC 0.39 to 0.49); without a cap, with a cap of
24 or with a logarithmic compression (LOG6) the AUC is 0.66 to 0.74. (c) No variant detects the events at 30 false
detections per hour with origin time required (0 or 1 of 16), because the window maxima are not at the events. (d) On
real noise with injected events, LOG6 keeps the robustness of the cap against transients (LOG6 minus cap -0.015
[-0.099, +0.043]); a cap of 24 loses it. (e) Capping the energies of BOMP-PK recovers almost all of the 2M-type gap (MT
with cap minus BOMP-PK with cap +0.026 [-0.009, +0.081]); per-station normalisation alone does not. Appendix Y.32.

**Gate for the next registered real test (fixed before any further development).** On Ridgecrest CAL, with origin
time required, at 30 false detections per hour, a candidate must recall at least 5 of the 16 events of M >= 2 on N8 or
U8. Only then is a test registered on the Hi-net swarm of MatchLocate2, whose waveforms remain unread. Candidates: a
reading less sensitive to coda and S complexity (envelope reading, higher-frequency matched filter; down-weighting S
was tried and fails, C18g),
joint relocation of the readings with a coherence check before the grid search, and LOG6 as the statistic.

---

## TASK 1q — Search with certified pruning and p4est

**Proof of concept done (Appendix Z).** Branch and bound on an octree
(x, y, depth, origin time) with Lipschitz bounds of the eikonal: catalogues identical to the
dense grid in 42 of 42 windows, with 2.3 to 3 times fewer evaluations. p4est is not needed for
the search in a regional window; it serves to distribute national continuous monitoring on a
cluster (one window or region per tree of the forest). Acceptance criterion for any
implementation with p4est: the same node-by-node equivalence with the dense search.

---

## TASK 1l — Orientation audit and incidence-based control

**Priority: very high. Data already downloaded for CHI, OCA, SJC and ZAR.**

### Objective

Separate the three explanations of a failing control: orientation of the horizontals,
structure beneath the station and incidence geometry. At CHI the P is at chance and the S hits
0.677; with near-vertical incidence the P has little horizontal component and the SV arrives almost
horizontal and radial, which would produce exactly that.

### Steps

1. Orientation by P particle motion with OrientPy (method of Braunmiller et al. 2020,
   ref. [91]) on teleseismic and regional events, per year, to detect changes over time
   (ref. [93]). Independent check by Rayleigh-wave arrival angle (DLOPy, ref. [94]).
2. Repeat the Task 1 control restricted to incidences of 30° or more, and broken down by
   event depth (less than 70 km, 70 to 200 km).
3. Control with the S on the radial component: angle folded to [0°, 90°], because SV is radial and SH
   transverse; both inform the back-azimuth modulo 90°.

### Acceptance criterion

If the P control with incidence of 30° or more exceeds 0.60 at CHI, the failure was one of geometry and the
station enters with incidence weighting. If the audit finds a rotation, it is corrected and the
control is repeated. If nothing changes, the site explanation remains, now supported.

### Estimated effort

Three to five days.

---

## TASK 1m — Sensitivity of the labelling to the catalogue, and coherent clutter from uncatalogued sources

**Priority: high.**

### Objective

The Task 1 labelling used the ISC catalogue, which is less complete in Colombia than that of the SGC. A
pick of an uncatalogued microearthquake remains as clutter, and its coherence with its own source is not
measured by $\eta$, whose reference is always a catalogued event. That is precisely the physically
realisable impostor of Appendix T.

### Steps

1. Relabel the ZAR period with the SGC catalogue and with that of SeisMonitor (Castillo, Siervo
   and Prieto 2024, ref. [90]); recompute $\eta$ with each.
2. Template correlation of the catalogued events over the clutter windows; picks
   that fall in a new detection are relabelled as "uncatalogued event".
3. Calibrate the false-positive rate of the template correlation on shifted noise
   (Muir, Fernando and Barrett 2023, ref. [99]).

### Acceptance criterion

$\eta$ with the SGC catalogue within the Task 1 CI. If it rises above 0.20, the real clutter
contains realisable impostors and Appendix T becomes central.

### Estimated effort

One to two weeks.

---

## TASK 1n — Calibration of the estimators of $\eta$ and of $\hat r_s$, and real PhaseNet residuals

**Priority: high. Synthetic and data already downloaded.**

### Objective

Two uncharacterised quantities: the bias and the coverage of the CI of $\eta$ (the verification uses one
replicate per proportion), and the real distribution of the timing error of the picks, which Appendix Y
showed decides how much polarisation contributes over the best kinematic score.

### Steps

1. 200 replicates per planted proportion (0, 0.1, 0.2, 0.3, 0.5) and per size (50, 150, 300 picks):
   bias, coverage of the CI95 and rate of each verdict.
2. PhaseNet residuals against the catalogue arrivals at the four stations already downloaded, per
   phase, distance and magnitude; fit of a heavy-tailed distribution.
3. Redo Appendix Y.2 with the pick errors sampled from that distribution instead of
   Gaussian ones.

### Acceptance criterion

Coverage of the CI95 of $\eta$ of at least 0.93. The result of step 3 decides whether the RSNC is in
the regime in which polarisation contributes over the best kinematic score.

### Estimated effort

One week.

---

## TASK 2 — Reimplement the tuned azimuth estimator

**Priority: high.** Well specified, with no external data, and with a large payoff.

### Objective

Close the order-of-magnitude precision gap that Appendix U documents between the
estimator without a velocity model and the methods that do use one.

### Why it could not be done here

It was attempted and failed three times because of different implementation errors. The code is attached as
`a_i_tuned_attempt.py` and **does not work**: it fails in 62 to 100% of the realisations. It needs
interactive debugging, not more blind iterations.

### The three causes already identified

1. **Arbitrary eigenvector sign.** The polarisation directions are obtained from a
   covariance decomposition, whose sign is undefined. **Lines** must be intersected, not
   directed rays. The `ts>0` test in the attached code is invalid.
2. **Gates too strict.** `len(inl)<3` and the support threshold reject valid configurations
   at the clutter level tested.
3. **Lipschitz misapplied.** The C-III test holds between arrivals of the **same phase**. The
   inlier set is selected by azimuth and mixes P and S, so the gate rejects
   correct sets. The phase of each inlier must be labelled, or the test applied per phase
   group.

### Steps

1. Correct the three causes.
2. Verify with a case constructed so that it works: no clutter, no noise, source away from the
   origin. It must give zero error. **This step is not optional** — the result that was about to be
   reported as spectacular was an artefact of the source being at the origin, as was the estimator's
   failure value.
3. Add clutter progressively and measure the failure rate.
4. Sweep back-azimuth deviation due to lateral heterogeneity: fixed bias per station from 0 to 15°
   plus noise from 0 to 8°, which is the range cited for heterogeneous regions. Task 1 measured the real
   range on the CM network and it is wider: median deviation of 9° at the only station that passes the
   control (ZAR) and of about 35°, with a coherent fraction equal to chance, at CHI. The sweep must
   include stations whose azimuth is pure noise, in the proportion given by Task 1c.
5. Compare against travel times and against C-IV under structural model error, replicating the
   protocol of Appendix L.

### Acceptance criterion

Failure rate below 10% with 25 clutter arrivals, and median error below 3 km under
22% model error. If this is achieved, Appendix U is rewritten and A3 moves from a positive instance
with high cost to a main result.

### Estimated effort

1 to 2 weeks.

---

## TASK 3 — Determination of the band supported by the velocity model

**Priority: high.** It is a blocking gate for WP1 and affects all work with real data.

### Objective

Determine empirically, per domain, the maximum frequency at which the available velocity model
supports waveform fitting.

### Why it could not be done here

It requires downloading published velocity models, real seismograms and a wave solver.

### Steps

1. Select the candidate 3-D model for the domain. For Colombia the public ones are the local-earthquake
   tomography of Chiarabba et al. and the ambient-noise crustal models of Poveda et al.
   on a half-degree mesh.
2. Select between 10 and 20 well-located events from the catalogue, with known mechanism.
3. Synthesise seismograms with that model. For 1-D, `Instaseis` or `pyprop8` suffices; for 3-D,
   SPECFEM3D or AxiSEM3D is needed.
4. Filter observed and synthetic in increasing bands: 0.02-0.05, 0.05-0.1, 0.1-0.2, 0.2-0.5,
   0.5-1, 1-2 Hz.
5. In each band compute the maximum cross-correlation and the time shift between observed and synthetic.
6. The supported band is the highest one in which the median correlation exceeds 0.7 and the shift
   stays below half a period.

### Expected output

One table per domain: maximum supported band, and the magnitude range that band implies.

### Acceptance criterion

If the supported band is lower than the declared primary one, the transfer to real data runs in the
supported band and with larger magnitudes, and this is reported. **Under no circumstances is it run
in a band that the model does not support.**

### Estimated effort

3 to 4 weeks if a 3-D solver is available; 1 week with a 1-D model.

---

## TASK 4 — Verification of the computing budget

**Priority: high.** The document explicitly states that Tables 3 to 5 are a planning model,
not a measurement.

### Objective

Replace the planned figures with real measurements.

### Why it could not be done here

There is no GPU or elastodynamic solver installed.

### Steps

1. Install SPECFEM3D_Cartesian or SeisSol.
2. Run the three configurations of Table 3: 1, 2 and 5 Hz on the declared domain, measuring
   wall-clock time, peak memory and output volume.
3. Measure the real cost per objective evaluation of Table 4, including the transport misfit
   with gradients.
4. Build the SGT database by reciprocity for one geometry and measure real time and storage.
5. Re-issue Tables 3, 4 and 5 with measured values and recompute the total.
6. Confirm the institutional allocation against the measured total.

### Acceptance criterion

The **ordering** between architectures must hold: solver in the loop much more expensive than
SGT plus beams. If it does not hold, the feasibility argument of Section 11 collapses and must be
redone.

### Estimated effort

2 weeks.

---

## TASK 5 — C2, the geometric failure class in real networks

**Priority: medium-high.** Cheap and potentially of high operational impact.

### Objective

Identify deployed real networks whose geometry falls in the failure class measured in
Appendices M and J: low variance in epicentral distance, or one-sided configuration.

### Why it could not be done here

It requires querying the FDSN station service programmatically. Searching for a handful of networks
on the web would be a biased sample, which is exactly what would invalidate the claim.

### Steps

1. Download station metadata for all networks of one or several regions.

```python
from obspy.clients.fdsn import Client
c = Client("IRIS")
inv = c.get_stations(level="station", channel="HH?",
                     minlatitude=-5, maxlatitude=15,
                     minlongitude=-82, maxlongitude=-66)
```

2. For each network and each point of a mesh of target sources, compute the leading-order
   observability constant and the three geometric statistics: variance of epicentral
   distance, azimuthal gap, and the smallest eigenvalue of the direction covariance.
3. Classify as failure class the network-region combinations where directional coverage falls
   below the percentile that Appendices I and M associate with zero margin.
4. For the three worst-ranked networks, reproduce the false-event-rate experiment on
   their **real geometry**.

### Acceptance criterion

At least three real networks identified and the collapse reproduced on their geometry. If none
appear, C2 is reported as null and the warning about degenerate geometries is restricted to
idealised models.

### Estimated effort

1 to 2 weeks.

---

## TASK 6 — B1 with multilayer ray tracing

**Priority: medium.** Closes the last gap in the theorem of Section 6.1.

### Objective

Extend the observability certificate to media with several interfaces, triplications and converted
phases.

### Why it could not be done here

The two-layer case with a head wave has a closed form and is already in Appendix Q, with
$\rho=0.963$. Several interfaces require two-point ray tracing with a search over multiple branches,
which is not feasible to write blind.

### Steps

1. Install `obspy.taup` for 1-D models with several layers, or `pyrocko.cake`.
2. For each source-receiver pair, enumerate **all** ray branches, not only the first arrival.
3. For each branch obtain the slowness vector at the source end. Verify that its
   modulus is $1/V(\mathbf x_s)$ — that eikonal check detected the correct implementation in
   Appendix Q and must be repeated here.
4. Assemble the certificate as a sum over branches and correlate it against the numerically computed
   coercivity, replicating the protocol of Appendix Q.
5. Add a head-wave amplitude law with correct geometrical spreading, to separate
   how much of the precision loss of Appendix Q is due to the theory and how much to the approximate
   amplitude model used there.

### Acceptance criterion

$\rho\ge0.95$ with a sum over branches in a model with at least three layers and triplications.

### Estimated effort

3 weeks.

---

## TASK 7 — D1 with a realistic picker

**Priority: medium.** Appendix R parameterises the loss rate; the real pattern is missing.

### Objective

Replace the uniform and independent pick loss with the correlated pattern of a real picker,
which loses picks within the coda of earlier arrivals.

### Why it could not be done here

It requires a calibrated picker running on waveforms with realistic coda.

### Steps

1. Generate synthetics with realistic coda, calibrated against real records from the chosen region.
2. Run PhaseNet or EQTransformer via SeisBench on those synthetics.
3. Measure the loss rate **as a function of the delay relative to the previous arrival**, which is the
   correlation structure that is missing.
4. Feed that pattern into the experiment of Appendix R instead of the uniform loss.
5. Re-measure the crossover threshold above which the pick-free channel dominates.

### Acceptance criterion

The crossover threshold measured with correlated loss, with a confidence interval. If it differs
substantially from the one obtained with uniform loss, Appendix R is rewritten.

### Estimated effort

2 weeks.

---

## TASK 8 — Effect sizes and power for the confirmatory family

**Priority: medium.** Without this the margins of Appendix D are provisional.

### Objective

Replace the rates measured on 6 to 12 realisations with effect sizes with intervals, and
derive the number of confirmatory realisations.

### Steps

1. Run the Tier-1 pilot as defined in Section 13.0, including **the design-validity
   step**: show that each family of constraints acts on a case constructed so that it
   acts, and that each sampling-sensitive channel is adequately sampled.
2. Estimate the paired standard deviation per endpoint.
3. Apply the document's sample-size formula, and where the pilot rejects normality,
   compute power by simulated paired permutation.
4. Freeze all the margins of Appendix D before running anything confirmatory.

### Estimated effort

3 weeks, almost all unattended computation.

---

## TASK 8b — Second independent synthetic generator

**Priority: high.** The sixteen pre-studies share **a single generator**, one clutter
model and, except in two cases, a single homogeneous medium. The sufficiency ladder and the channel
hierarchy could be properties of the simulator and not of the problem.

### Steps

1. Build a second generator with: a layered medium of at least three layers, clutter derived
   from **real coda** instead of synthetic coda, variable source mechanisms, and background noise taken from
   real records of the chosen region.
2. Repeat the five supported results of Appendix W on that generator.
3. Report agreement or disagreement per result, not in aggregate.

### Acceptance criterion

A result is considered a property of the problem only if it is reproduced in both generators within
the registered margin. Those that are not reproduced become properties of the simulator and are withdrawn from
the claims.

### Estimated effort

3 to 4 weeks.

---

## TASK 9 — B2 with sufficient replicates

**Priority: medium-low.**

### Objective

Decide whether a network designed by maximising the observability constant outperforms one designed
by minimising the azimuthal gap.

### Why it remained open

Appendix N shows that the design criterion works —the networks optimal for observability
have a constant two orders of magnitude larger— but the location-error endpoint **saturates
at zero** and does not separate them. The signal appears in the false-event rate with four stations, with
only ten realisations.

### Steps

1. Repeat with the corrected endpoint, false-event rate, and with the number of replicates that comes out
   of Task 8.
2. Extend to incremental placement: given an existing network, where the next station goes to
   maximise the minimum constant over the target region. It is a minimax problem and it is the question
   that an agency actually asks.

### Estimated effort

1 week once Task 8 is done.

---

## Directions already closed: do not retry them

These were tested or ruled out analytically. They are documented in Appendix V of the
proposal. They are listed here because they are exactly the ones that a new reader, or you yourself in
six months, will find obvious.

- **Order selection by subspace (MUSIC, MDL, AIC) on the cross-spectral matrix between
  stations.** It does not count sources. The rank-K structure requires a direction vector that is stable
  across windows, and a transient train with moveout does not have one. This closes the whole family.
- **Betti reciprocity identities as a test on the data**, even with two
  superimposed events. The identity lives in the source positions, which are not observed. It is an analytical
  block, not an empirical one.
- **Physicality of the moment tensor** as a discriminant. Fitting clutter with a six-element
  basis produces nearly deviatoric tensors by chance; 85 of 491 false hypotheses came out more
  deviatoric than the true one.
- **Observability-based stopping rule in adaptive refinement of the hypothesis space.** It never
  fires: the attainable resolution limit is finer than the cell the tree descends to.
  The score-guided octree does work, but NonLinLoc and Qseek already do that.
- **Velocity contrast law**, **graceful degradation in location**, **overlap
  as an association difficulty on its own**. Registered nulls E1 to E4.
- **Group BLASSO with per-station and per-phase shifts, off the grid via sliding Frank-Wolfe, as a
  better detector than block search with the same tolerance** (G1-M). Refuted
  in Task 2E (+0.000 to +0.039, all CIs contain 0), with no implementation error or budget
  disadvantage. Local refinement from the grid does not reach the missing events, which are absorbed by
  atoms at 5-10 km.
- **More physics in the atom as a source of recall**: 1D velocity model (oracle or calibrated by JHD),
  source-pulse-invariant subspace, coarser grids with sliding. Null or negative in the
  exploration of seed 999 (C1, C6, C7); the margin persists with the exact medium.
- **Relocation of atoms by evidence without travel times** (polarisation, amplitudes, S/P ratio).
  It loses recall in all variants, also with oracles (C9): the Cramér-Rao bound of those channels is
  4.4-5.6 km versus 0.4-0.7 km for the times. An information limit, not an algorithmic one.
- **Gaussian timing prior within the detector with catalogue corrections.** No gain (C10): with
  corrections off by 0.12-0.24 s and a 3 km cell that adds 0.18-0.30 s, the honest prior is flat
  within ±0.25 s. It only pays off with statics known at the level of the model error.

## Ideas explored and not executed, in case you want to take them up again

These came out of the exploration but were never tested. None is foundational.

- **Free-surface correction to recover depth.** The rank budget of
  Section 9.2 shows that the incidence is ill-conditioned because the free-surface conversion
  distorts the apparent angle. That distortion is analytical for a half-space with given
  $V_P/V_S$. Correcting it could extend the single-station case from $(x,y,\tau)$ to
  $(x,y,z,\tau)$. A closed question, cheap, with a large payoff.
- **Eikonal-exact slowness network.** A network that learns $\mathbf p(\mathbf x_s,\mathbf r)$ with
  $|\mathbf p|=1/V(\mathbf x_s)$ imposed by architecture —normalise the output and scale by
  $1/V$— would make the certificate computable in real models at negligible cost, which is what
  Task 6 needs in order to scale.
- **Learn the polarisation operator instead of the wavefield.** If polarisation is what
  discriminates, the object to learn is not the Green's kernel but the map from source position
  $\to$ predicted direction per station, from $\mathbb R^3$ to $S^2$. A tiny learning problem
  compared with a wavefield, and exactly on the quantity that carries information.
- **Optimal picker threshold as a function of the downstream stack.** If the constraints reject
  cheap clutter, the detector should be configured more permissively. Nobody calibrates detectors this way.
  Testable prediction: the optimum for travel times and for the constrained stack lie in different
  places.
- **Set prediction for unknown cardinality.** Architectures that output sets,
  trained with a matching loss, with GOSPA as the loss and the physical constraints as a
  differentiable penalty.
- **Wave-equation residual as a discriminant.** Require the data to satisfy
  $\rho\ddot{\mathbf u}-\nabla\cdot\boldsymbol\sigma=0$ locally before asking where it
  comes from. It is the purest possible discriminant —it needs no source, no origin time, no velocity
  model— but it requires **spatial derivatives of the field**: a dense array, DAS or rotational
  sensors. With stations tens of kilometres apart it cannot be evaluated, and that is why it remained untested.
- **Focusing entropy in time reversal.** Classical backprojection uses the peak amplitude;
  using instead the spatial concentration of the reversed field changes the question to "does this focus
  anywhere?". Clutter does not focus anywhere.
- **Kramers-Kronig-type causality** per trace, and **higher-order statistics** (bispectrum,
  kurtosis): they discriminate deterministic signal from Gaussian noise without any model.
- **Cardinality certificate via statistical dimension** (McCoy-Tropp). It costs paper, not computation,
  and attacks the gap left by E5 and Section 6.4: without rank, cardinality is left without a
  certificate, and the phase transition of convex demixing does provide one.
- **Bruckner and Yamamoto (2000)**, determination of point sources from interior point
  observation, recovering even the number N of sources. Check whether the Proposition of
  Section 6.1 is a quantified version of that theorem; if it is, citing it changes the status of
  C-I from an original construction to an instance of a known result, which is more honest and stronger.

---

## Warning about comparison protocols

Task 0 yielded two lessons that apply to any future comparison in this project.

**Threshold on a normalised statistic.** Two methods whose scores have different scales
are not compared with an absolute threshold. The ablation had a maximum of exactly 1.000 by construction
and the constrained model between 0.51 and 0.68 depending on the realisation; the absolute threshold gave one a
natural operating point and the other none. Three complete runs were lost because of this.

**Average over the factor that dominates, not over the one that is convenient.** Fixing one geometry per cell and
drawing thirty clutters seemed reasonable and was the biggest error: between two seeds the same method
went from 0.13 to 0.70 false events in the same cell. Always ask yourself which is the largest source of variance
and randomise **that one**.

---

## Additional warning about baselines

The Task 0 error survived six audits, two of them external. Nobody detected it because
all participants shared the same assumption: that the implemented counterfactual represented
the family of methods it named. Before accepting any comparison in this project,
ask explicitly **what concrete code was executed** to produce the reference column, and
be suspicious if the answer is "an implementation of" instead of "the published package X version Y".

---

## Methodological warning that applies to all tasks

In the sixteen pre-studies of the document, **five times a design flaw inverted or nullified the
result**, and all five were detected only because the number turned out implausible. Others followed
(E8 to E18 of the manifest), the most costly of all being the ablation ties (E15),
detected by code audit and not because the number looked odd: a plausible number
can also be an artefact. The modes were:
a geometry that disables the constraint under test; a perfect picker that eliminates the mechanism
to be measured; a grid coarser than the effect; a constraint evaluated on quantities
predicted by the model instead of on observations, which renders it inoperative; and a failure value
of the estimator that coincided with the correct answer.

Before believing any result, including a positive one:

1. Verify that the constraint or the channel **can act** on a case constructed so that it acts.
2. Verify that the sampling resolves the scale of the effect.
3. Verify that the failure value of the estimator **does not coincide** with the correct answer.
4. Apply each physical test to **observations**, not to quantities predicted by the model.
5. Check that the chosen endpoint has variance in the regime tested, and report which one saturated.
