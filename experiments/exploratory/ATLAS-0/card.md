# Experiment card: ATLAS-0. The failure atlas of published detection pipelines

Stage: S0, final draft of 24 October 2026, ready to be anchored (ADR-0004). It must be anchored before any comparator
is run on atlas data, because the phase-4 gate is stated on the regimes fixed here. The reproduction gate (below)
runs the comparators on their own published examples only and may run before the anchor.
Owner decisions: the single "beat" rule, ADR-0008, accepted by the author on 24 October 2026; this card copies it and
does not restate a different one. Every other choice marked "default" was taken by the atlas thread from the plan
review, the phase-5 theory, the model v2 specification and the adversarial critique; the author can change any of them
before the anchor.
Code: `src/ttmt/atlas/` (cells, scoring under ADR-0008, instrumentation; tests in `tests/test_atlas.py`) on top of
`src/ttmt/stats.py`. Comparator wrappers and tuning: thread "Comparadores y ajuste", which hands the atlas detections
and picks in the schema of `ttmt.atlas.schema`.

## Question
In which measurable regimes does each published detection pipeline, tuned at its best on development data, lose
events that are detectable, and through which stage (picker, associator, location) does it lose them?

## Mechanism tested
The failure of each pipeline stage, separately:
- picker: arrivals lost, merged into one pick, or picked as the wrong phase, against the separation to the nearest
  other arrival at the same station and against the station signal-to-noise ratio;
- associator: events missed although enough of their arrivals were picked, events merged (one detection built from
  the picks of two or more true events), events split (one true event spread over two or more detections), and
  secondary false events;
- location: events found under the epicentre-only rule but missed under LOC+T (origin error above 0.5 s), against
  distance to the centroid.
The velocity model, the station set, the noise and the injected events are held fixed across methods within a stratum.

## Hypotheses and refutation
Primary (H-ATLAS): every core pipeline loses at least 20 per cent of the detectable events of the target regime T.
The loss of a pipeline is its recall on the isolated events (separation > 3 s), standardised to the mix of
signal-to-noise, station and distance strata of T, minus its recall on T (`ttmt.atlas.froc.pooled_loss`), at 30 false
detections per hour under LOC+T-V. Refuted if, for at least one core pipeline, the upper bound of the stratified
block-bootstrap 95 per cent interval of its loss is below 0.20. Confirmed if the lower bound is at least 0.20 for
every core pipeline; otherwise undecided. The same loss is reported, without a test, in the two halves of T and by
station count.

Mechanism predictions, from the model v2 specification (Section 12; copied here so the anchor fixes them). "Loss" is
the pooled loss above, on the stratum named in the third column (pooled over the axes it does not name); every
comparator is tuned (Section "Tuning"). These are descriptive predictions of the atlas, each with the interval of its
quantity; they are not decision tests of ADR-0008.

| ID | Pipeline | Prediction | Stratum | Refuted if |
|:--|:--|:--|:--|:--|
| P-PN1 | PhaseNet + PyOcto | loss >= 0.15 at separation 0.1-0.5 s with the best peak distance | separation x SNR 8-20 x 8 stations | loss < 0.05 with peak distance 0.1 s |
| P-PN2 | PhaseNet alone | fraction of injected pairs with two separate probability maxima < 0.5 at separation < 0.2 s | separation < 0.1 and 0.1-0.5 s | >= 0.8 of pairs separated at 0.1-0.2 s |
| P-PN3 | PhaseNet + PyOcto | the pair loss does not depend on K_mec/E | K_mec/E terciles (pure synthetic only, CR-14) | loss grows with K_mec/E (trend test p < 0.05) |
| P-NO1 | PhaseNO + PyOcto | at equal recall, precision >= 0.10 below PhaseNet with 4 stations, gap smaller with 12 | stations 4 / 8 / 12 | precision gap < 0.05 with 4 stations |
| P-NO2 | PhaseNO | loses overlapping pairs at least as much as PhaseNet | separation 0.1-0.5 s | PhaseNO loss < 0.5 x PhaseNet loss |
| P-EQN1 | EQNet | >= 30 % of LOC hits fail LOC+T; origin error grows with distance | all | < 10 % of LOC hits fail LOC+T |
| P-EQN2 | EQNet | false detections per hour >= 2x on the transient noise day at a fixed threshold | noise type | ratio < 1.3 |
| P-GM1 | GaMMA | loss >= 0.20 at separation 0.1-0.5 s | separation 0.1-0.5 s | loss < 0.05 |
| P-GM2 | GaMMA | >= 10 % secondary false events at >= 500 events/day and 100 % false picks | associator scenario A | < 3 % |
| P-PO1 | PyOcto | with timing noise 0.2-0.4 s no single tolerance gives recall >= 0.8 and precision >= 0.9 | associator scenario B | some tolerance gives both |
| P-RE1 | REAL | recall of the second event of a pair <= 0.3 | separation 0.5-3 s | >= 0.6 |
| P-QM1 | QuakeMigrate | false detections >= 2x those of MT-FW at equal recall with transients at 1-2 stations | noise type | ratio < 1.3 |
| P-QM2 | QuakeMigrate | pair recall at separation < 1 s <= that of PhaseNet + PyOcto | separation 0.1-1 s | QuakeMigrate exceeds it by >= 0.10 |
| P-MA1 | MALMI | loss >= 0.05 at separation 0.1-0.5 s | separation 0.1-0.5 s | loss < 0.05 |
| P-TM1 | template matching | recall >= 0.9 on repeated sources, <= 0.3 on sources > 5 km from every template or with another mechanism | template distance x mechanism | recall >= 0.6 on new sources |
| P-TM2 | template matching | false detections per hour with 4 stations >= 2x those with 8 at the same MAD threshold | stations 4 / 8 | ratio < 1.3 |
| P-OWN1 | MT-FW / BOMP-PK (prototype-v1) | >= 1 duplicate per 10 events at SNR >= 60 (deflation residual, E62) | SNR 60 / 120 | < 1 per 50 |
| P-OWN2 | MT-FW with cap (prototype-v1) | event-versus-background AUC < 0.55 in a continuous sequence | rate >= 500 events/day | AUC >= 0.65 |

"Refuted" and "confirmed" use the 95 per cent interval of the quantity: refuted when the interval lies entirely on the
refuting side of its threshold, confirmed when it lies entirely on the predicted side, undecided otherwise. A
prediction whose pipeline could not be run, or whose stratum has fewer than 100 targets, is "not tested".

Decision rule (model v2 specification, Section 12): P-PN1, P-GM1, P-RE1, P-QM2 and P-MA1 justify the target regime.
If three or more of the five are refuted, the overlap regime is not a gap, and goal G1-v2 is reformulated in a
decision record before component M6 is built.

## Regimes and strata (fixed here, not after the atlas)
Regimes (ADR-0008, item 1), computed from the truth and the data availability only, never from a pick or a detection
(`ttmt.atlas.froc.regime_of`):
- T, the target regime: separation in [0.1, 1) s;
- C, the complement: every other target event, including separations below 0.1 s (the information-limited bin of
  Appendix AB, Cor. A.7) and above 1 s. "Isolated" is the part of C above 3 s, the reference of the losses.

Strata for description only (`ttmt.atlas.cells`; no decision counts strata):
1. Separation: the median, over the stations of the subset, of the absolute difference between the event's P arrival
   and the nearest arrival (P or S) of any other event at that station: < 0.1, 0.1-0.5, 0.5-1, 1-3, > 3 s.
2. Network signal-to-noise ratio: the median over the subset stations of the ratio between the peak absolute amplitude
   of the event's own signal in [t_P, t_P + 2 s] and the standard deviation of the noise in [t_P - 5 s, t_P - 0.5 s],
   vertical component, after the common 2-20 Hz band-pass; bins < 4, 4-8, 8-20, >= 20 (linear; about 12, 18 and 26
   dB). For injected and synthetic events the signal is the injected waveform alone and the noise the noise record
   alone, so the ratio is exact; for catalogued events both are read from the record.
3. Stations: controlled subsets of k = 4, 6, 8 and 12 stations, the k stations nearest the network centroid whose
   continuous data cover [origin - 10 s, origin + 60 s] without a gap, ties broken by station code. A network with
   fewer than k such stations contributes nothing to that bin. Every method is run on every subset.
4. Distance: epicentral distance to the centroid of the network's development catalogue: < 10, 10-30, > 30 km.
That gives 5 x 4 x 4 x 3 = 240 strata. Each is described with its recall per method and a stratified block-bootstrap
95 per cent interval when it holds at least 100 targets, and listed without an interval otherwise.

Reported covariates, not strata: magnitude; K_mec/E of the pair (terciles), only where the true moment tensors are
known, that is in pure synthetic data or with injection EGFs of events with a catalogued tensor (CR-14; every
conclusion that depends on K_mec/E states that scope); spatial separation and energy ratio of the pair; depth; and,
for injected events, the largest normalised correlation between the injected waveform and any template or deflation
atom (CR-06; every loss is also reported in the stratum of correlation < 0.5).

Contamination stratum (CR-24): the 4-station subsets are also run with one and with two of their stations given noise
scaled by 3 and by 10 (the same injected events), to describe how each method degrades when the median of four
stations falls on noisy ones.

Associator scenarios (for P-GM2 and P-PO1, on pick lists, outside the strata): scenario A crosses event rate
{isolated, 100, 500, 2000 per day} with false picks {0, 30, 100, 300 %}; scenario B crosses timing noise per arrival
{0, 0.1, 0.2, 0.4 s} with PyOcto tolerance, on picks of the generator's true arrivals (the protocol of Puente et al.,
2025).

## The single rule (ADR-0008)
The phase-4 gate G4 and every later "beat" claim use ADR-0008 as written there; this section copies its operative
parts and adds only what the card must fix.
- Hit rule LOC+T-V: ADR-0007 matching (one to one, decreasing score, 5 km and 0.5 s) restricted to detections whose
  normalised cost to the event is at most half their cost to every other catalogued event (`ttmt.stats.label_voronoi`,
  called by `ttmt.atlas.froc.label`). For isolated events it is LOC+T. A detection matched to a catalogued event that is
  not a target (an event of the noise record, or an event outside the station subset of the run) is neutral.
- Negative control: the duplicator of the best comparator runs in every analysis: each detection twice, at origin
  +- half the separation of the catalogued event nearest to it (`ttmt.atlas.froc.duplicator_shift`; one shift per
  record would land copies on partners at other separations); if it beats its comparator in T at one-sided 0.025, the metric is broken and no claim is
  made. Duplicates of isolated events are reported per method (`ttmt.atlas.froc.isolated_duplicates`).
- Endpoint: recall of the targets of a regime at 30 false detections per hour (exact FROC, operating point per
  method on the evaluated data), pooled with prevalence weights. Default: weights 1 (each evaluated target counts
  once), because no prevalence of the strata in real dense sequences is available at the date of this card (CR-04); if
  the reference-catalogue thread delivers it before the anchor, its table is added here and replaces the default.
- Test: block sign-flip test of the paired recall difference (`ttmt.stats.signflip_pvalue`); the block is an injected
  network realisation in semi-synthetic data and a time block of real data; one-sided 0.025; stratified block-bootstrap
  interval reported with every test, indicative with fewer than 30 blocks.
- Beating the comparators: an intersection-union test over the comparator set (the claim holds only against every
  comparator). ADR-0008 item 7, as amended by PR #13 (`ttmt.atlas.froc.gate`):
  - H1: superiority in T (H0: recall difference <= 0 against every comparator), exact sign-flip test at 0.025;
    and then
  - H2: non-inferiority in C with margin 0.05 (H0: difference <= -0.05 against every comparator), by the studentised
    sign-flip test on d_j + 0.05 n_j, at a calibrated level: the largest of {0.025, 0.02, 0.015, 0.01, 0.005} whose
    simulated size, worst case over method correlation 0.5 to 0.9 and intra-block correlation 0 to 0.2, plus 1.645
    Monte Carlo standard errors, is at most 0.025 (`experiments/exploratory/P6-SAP/sample_size.py calibrate`). That
    level is 0.015 from 50 to 300 blocks; a design above 300 blocks is calibrated before the run. H2 needs at least 50
    blocks with C targets; with fewer it is not tested and no beat claim is made.
  - H2 is tested only if H1 is rejected, so the family-wise error is at most 0.025.
  - Secondary hypotheses, tested only if H1 and H2 are rejected, form one Holm family at 0.025: superiority in
    [0.1, 0.5) s, superiority in [0.5, 1) s, and non-inferiority in C at margin 0.03, whose Holm threshold is
    multiplied by its own calibrated level divided by 0.025 (default: the level of H2, until `sample_size.py
    calibrate` gives one for margin 0.03 at the design's blocks).
- Where G4 is measured (K9): primary on semi-synthetic injection into real development noise (known truth), with
  injection and deflation libraries disjoint by event and by time block; secondary on real development data scored
  against the truth of task T3.V, where the point estimate in C must not fall below -0.05. G4 is not confirmatory.
  Truth precision (phase-6 plan 3.4): on real truth, a pair enters T only if its origin separation is at least 9
  times the truth's origin-time standard deviation measured on T3.V; closer events are excluded (neutral,
  `ttmt.atlas.froc.truth_precision_targets`). Semi-synthetic truth is exact and needs no exclusion.
- Planning: semi-synthetic analyses are sized for a gain of 0.05 in T (`experiments/exploratory/P6-SAP/sample_size.py`,
  part 4 gives the number of blocks).

Anchor precondition: PR #13 (calibrated level of H2 and studentised test in `ttmt.stats`) is merged; this section
already states the amended rule.

The atlas itself (phase 3) makes no beat claim: it describes strata and tests H-ATLAS and the predictions above.

## False-detection budget and operating point
Primary: 30 false detections per hour of record, as in R1; secondary: 10 and 100. The threshold of a method is set
on the whole evaluated record at that budget (false detections cannot be assigned to strata), and the recall of a
stratum is the fraction of its targets hit above that threshold; in the block bootstrap the threshold is re-derived in
every replicate. Every method also reports a deployable operating point fixed on the tuning partition and applied
unchanged to the reading partition.

## Truth design
1. Semi-synthetic (primary): recorded events (three-component waveforms of catalogued events at each station) injected
   into real continuous noise at controlled separations and scales, built by the generator v2 thread from noise
   windows supplied by the real-data thread. The event library used for injection is disjoint, by event and by time
   block, from any library a method may use for templates or deflation.
2. Synthetic generator v2, independent of the detector's forward model: another wavelet family, a layered medium with a
   free surface, synthetic coda, and a radiation operator that is not the detector's far-field operator.
   Domain gate (CR-05): picker-based comparisons on a synthetic generator are reported as comparisons only if, on its
   isolated events, the recall of PhaseNet + its associator is within 0.05 of its recall on semi-synthetic events of
   the same signal-to-noise bin; otherwise they are reported as a domain-shift diagnostic. P-PN*, P-NO*, P-GM*, P-PO*
   and P-MA1 are measured on semi-synthetic data.
3. Real catalogues as secondary evidence only. A catalogue produced by template matching is not used to evaluate
   template matching (so the Hi-net target list of R1 is not used for EQcorrscan). The Ridgecrest development hour is
   spent (protocol of phase 6, item 3.5, CR-23): it may serve as noise for injection, not as truth for a gate.

Truth table per target (`ttmt.atlas.schema.Truth`, agreed with the generator v2 thread): event id, origin time,
location, magnitude, target flag, pair id, per-station and per-phase true arrival times without and with the injected
timing error, per-station signal-to-noise ratio, K_mec/E of the pair when known, and the seed and configuration hash in
its metadata.

## Tuning (ADR-0008, item 9; protocol, Section 3)
Run by the thread "Comparadores y ajuste"; the atlas reads its ledger. Every method, the prototype and the phase-4
model included, gets 60 configurations drawn by the same random search with the same seed list over the search space
below, scored on the same development tuning partition by the same objective: development recall in T and in C at
30 false detections per hour under LOC+T-V, weighted as in the endpoint. Every configuration tried is appended to
`experiments/exploratory/ATLAS-0/tuning_ledger.csv` (method, trial, parameters, objective, time); its SHA-256 enters
the freeze, and a sensitivity table over the ledger is published. Tuning and reading partitions split the noise hours
and the injection seeds 30/70 by a hash of the hour, so no noise sample or injected event is shared.

Search spaces (default):
- PhaseNet (SeisBench weights `original` and `instance`; 190703-214543 as the reference of R0 and R1): peak distance
  in [0.05, 0.5] s (log-uniform), probability threshold in [0.03, 0.5] (log-uniform), weights;
- PhaseNO: the same peak distance and threshold, its published weights;
- GaMMA: minimum picks {3, ..., 8}, time sigma or tolerance in [0.3, 2] s, method {BGMM, GMM};
- PyOcto: minimum picks {3, ..., 8}, tolerance in [0.3, 2] s, association cutoff in [0.5, 2] s;
- REAL: minimum P picks {3, ..., 6}, minimum total picks {5, ..., 12}, window in [1, 5] s;
- QuakeMigrate: detection threshold in [1.2, 3], STA/LTA windows {(0.2, 1.0), (0.5, 5.0)} s, marginal window in
  [0.3, 1] s;
- MALMI: migration threshold in [0.1, 0.5], picker {PhaseNet, EQTransformer}, stack window in [0.3, 1] s;
- EQNet: event threshold in [0.05, 0.7];
- EQcorrscan: MAD threshold in [8, 12], minimum stations {3, ..., 6}, templates from the injection-disjoint library;
- BOMP-PK and MT-FW: their own conventions (deflation grid or fractional, E62; cap, Huber or ACAP) and thresholds.
Pick-based pipelines are tuned jointly (picker and associator parameters in one draw). Velocity model: the same 1-D
model for every method that takes one. Associator worker count fixed (GaMMA ncpu = 1, after the reproducibility
finding of the R1 pick stage). Pretrained weights: the tuning ledger lists the training sets, regions and periods of
every pretrained model, and an overlap with an atlas network and period is declared (CR-22).

## Methods and comparators
Core set: PhaseNet and PhaseNO with GaMMA and PyOcto; EQNet; template matching (EQcorrscan); BOMP-PK and MT-FW.
Extended set: EQTransformer, REAL, QuakeMigrate, MALMI. A method that cannot be installed or run is reported with the
reason; its predictions are "not tested".
Internal baseline (K6): BOMP-PK and MT-FW are run at the git tag `prototype-v1` (legacy grid deflation), which must
exist before any atlas run; the same code with the fractional deflation of E62 is reported apart, as "prototype-v1 +
E62".
Fine-tuned picker (CR-03): when the phase-4 model uses a picker fine-tuned on pairs, the best pick-based pipeline is
re-run with that picker and its associator re-tuned under the same budget, and it joins the comparator set of
ADR-0008, item 6. Phase 4 reports the 2 x 2 factorial (original or fine-tuned picker, x best associator or the model)
in T, and attributes a gain to the moment-tensor physics only against the fine-tuned-picker pipeline.

## Instrumentation (`ttmt.atlas.instrument`)
- Picks, per station and phase, against the observed true arrivals: matched one to one within 0.1 s (secondary 0.2 s);
  an arrival is "picked", "merged" (no pick of its own, while a pick of the same phase matched to another arrival within
  1 s at the same station lies within the span of the two arrivals widened by the tolerance), "wrong phase" (only a
  pick of the other phase within the tolerance) or "lost". Rates are reported against the separation bins and the
  station signal-to-noise ratio.
- Events: a detection's picks are traced to true events through the pick matching. A detection is "merged" when at
  least two true events each contribute at least 3 of its picks; a true event is "split" when at least two detections
  each hold at least 3 of its picks. A missed target is "picker-limited" when fewer of its arrivals were picked than
  the associator's minimum, and "associator-limited" otherwise. For waveform methods with no picks, a missed target is
  "merged" when a detection hits another event within 2 s and 5 km of it.
- Location: the fraction of LOC hits that fail LOC+T, and the origin and epicentre errors against distance.

## Data
Development partitions only (phase 3 is not confirmatory); no evaluation partition is read. Noise records: Alto
Tiberina 2010 (quiet days and the transient day of 2M), Ridgecrest 2019 (one hour of a sequence; spent as truth),
Hi-net 2012. All three are development-only (data cards). New real data enter only through the real-data thread, with
their data cards and the use ledger. Synthetic seeds are derived after the anchor (ADR-0004, item 4).

## Gate (reproduction, before any atlas data are read)
The harness reproduces published numbers of at least two comparators on their own example data, within 5 per cent
(relative) of the published recall or event count. Targets, fixed here:
1. PyOcto and GaMMA on the crustal synthetic scenario of Puente et al. (2025, Seismica), with the parameters of their
   Supplementary Tables S3-S7, against the recalls of their Tables S8-S9.
2. PhaseNet on the SeisBench benchmark of Münchmeyer et al. (2022, JGR), in-domain on ETHZ, against its published
   P-pick metrics, using the SeisBench `pick-benchmark` code.
3. QuakeMigrate on its published example (Winder et al., 2026, Seismica), against its stated event counts.
If fewer than two pass, the atlas does not run; the failure is reported with its cause. A second gate is the harness
itself: `tests/test_atlas.py` must pass, including the duplicator control under LOC+T-V.

## Cost and risks
PhaseNO and EQNet inference is run on the GPU machine. Sizing is per regime, not per stratum (ADR-0008); strata with
fewer than 100 targets are listed without intervals. A comparator that cannot be installed on the available
platforms (EQcorrscan on Windows) is run on Linux or reported as not tested.

## History of changes before the anchor
- 18 October 2026: S0 draft (separation 4 bins, SNR in dB, stations by count within 50 km, minimum 100 events per
  counted cell, a counting rule over cells).
- 19 October 2026: separation < 0.5 s split at 0.1 s (Appendix AB, Cor. A.7); station axis replaced by controlled
  subsets {4, 6, 8, 12} with a data-availability rule; distance axis added; SNR redefined as a linear amplitude ratio
  so that the generator's levels 6, 12 and 30 fall inside bins; mechanism predictions, tuning grids and reproduction
  targets added; associator worker count fixed.
- 24 October 2026: the counting rule over cells replaced by ADR-0008 (regimes T and C, LOC+T-V with the duplicator
  control, H1 then H2, IUT over comparators, G4 primary on semi-synthetic data); the 240 cells became strata for
  description; H-ATLAS restated on the pooled regime T; tuning grids replaced by a 60-trial random search with a
  ledger; prototype-v1 named as the internal baseline (K6); fine-tuned-picker comparator and 2 x 2 factorial (CR-03);
  domain gate per generator (CR-05); K_mec/E scope (CR-14); contamination stratum with 4 stations (CR-24); EGF
  correlation covariate (CR-06); pretrained-weight provenance (CR-22); Ridgecrest marked spent as truth (CR-23).
- 24 October 2026, later: duplicator shifted per detection (review of PR #12); ADR-0008 item 7 as amended by PR #13
  (H2 studentised at a calibrated level, 0.015 from 50 to 300 blocks, at least 50 blocks; Holm threshold of the
  0.03 margin scaled by its calibrated level); truth-precision rule of 9 sigma_t for real-data pairs in T.
