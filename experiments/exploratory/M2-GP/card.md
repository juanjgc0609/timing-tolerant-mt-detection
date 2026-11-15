# Experiment card: M2-GP. Gate G-M2, station terms in a Gaussian process

Stage: S2 (9 November 2026, cloud). FAILED: first step 0.150, confirmation of the declared second step on days
never used 0.191 (gate 0.30); coverage passes everywhere (0.94). Negative result for analyst picks. Plan: T4.8 (M2; build order after M6, owner's approval of
8 November 2026). Exploratory; no registered number changes.

## Question
In validation by blocks of time, does a per-(station, phase) Gaussian-process term on top of a published 1-D medium
lower the mean squared error of held-out calibration readings by at least 30 % against a constant term, with the
coverage of the +-2 sigma interval between 90 % and 98 % (gate G-M2 of `especificacion-modelo-v2.md`)?

## Mechanism tested
A constant station term absorbs the static part of the 3-D structure but not the part that changes with the source
position (of order |dx| times the gradient of the static). A GP in (x, y, z) models that part explicitly, and its
posterior variance gives a tolerance per position (small near calibration events, large far from them), which is
what M4 needs (Observation B.8). E53 is avoided by design: the terms of a block are fitted on the other blocks.

## Data
SCSN analyst picks of the Ridgecrest 2019 sequence from the SCEDC public bucket, 2019-07-09 to 2019-07-15
(`tools/fetch_scedc_phases.py`, which fixes the selection: type eq, M >= 1.5, box 35.45-36.15 N, 117.95-117.25 W).
9507 phase files, 3095 events, 133 952 picks; SHA-256 of `events.csv` 9e220d4d..., `picks.csv` 40781188... (checked
by the script). Stations kept: >= 300 P and >= 300 S picks (34 stations, CI, PB, ZY, GS). The R0 hour (2019-07-04
19:00-20:00 UTC, spent, CR-23) is outside these days and is not touched. These picks are NOT the readings of our
detector: they are analyst picks with catalogue hypocentres, the cleanest calibration set available. Declared
consequence: the hypocentres come from SCSN's own location with these picks, so part of the 3-D residual is already
absorbed by the location; the measured gain of the GP is therefore a lower bound of what it can model, and the
residual noise includes the analyst pick error (weights 1.0 / 0.8 / 0.5 / 0.2 are reported apart).

## Design
- 1-D medium Hadley-Kanamori (`ttmt.traveltime.HK`), VP/VS 1.73, stations at z = 0, local flat frame.
- Residual = observed travel time - 1-D time. Gross errors |residual| > 3 s removed before anything (same rows for
  all methods; count reported).
- 5 contiguous time blocks of events (`ttmt.v2.m2.time_blocks`); each block predicted by terms fitted on the other
  four (hyperparameters included).
- Methods: none (1-D only), const (training mean, training sd as tolerance; the spec's ablation), gp (Matern 3/2,
  constant mean by GLS, hyperparameters by marginal likelihood on <= 600 random training points, posterior on <= 3000).

## Gate (declared before the run)
Pooled over all held-out readings (P and S, all stations, all blocks): 1 - MSE(gp) / MSE(const) >= 0.30 AND
coverage of |error| <= 2 sigma_gp in [0.90, 0.98]. Also reported: per phase, per weight class, per block, per station,
robust spreads (MAD), the constant term's coverage.

## Hypothesis and refutation
Expected: the GP gains most at S and at stations near the sequence (strong lateral gradient across the Garlock /
Coso structure). Refuted if the MSE drop is below 0.30. If the drop passes but the coverage does not, the tolerance
is mis-sized: the diagnosis is heavy-tailed pick errors (look at the weight classes) before changing the kernel. If
the drop fails because the residual is not smooth in x, the spec's scaling gate applies (3-D eikonal or neural travel
times); that is a separate decision for the owner.

## Cost and procedure
CPU, cloud (4 cores), minutes.

    python tools/fetch_scedc_phases.py DATA_DIR 190 196
    python experiments/exploratory/M2-GP/run.py DATA_DIR > experiments/exploratory/M2-GP/run_out.txt

## Result, first step (`run_step1_out.txt`, days 190-196, plain GP; rerun with the per-station table in
`run_plain_out.txt`, identical numbers)
107 004 held-out readings (34 stations, 359 gross residuals removed, 0.33 %). **G-M2 FAILS on the MSE drop.**

| | rmse 1-D | rmse const | rmse gp | MSE drop | MAD const | MAD gp | cov +-2 sigma gp |
|:--|--:|--:|--:|--:|--:|--:|--:|
| all | 0.353 | 0.264 | 0.244 | **0.150** | 0.168 | 0.121 | **0.938** |
| P | 0.266 | 0.211 | 0.198 | 0.120 | 0.129 | 0.092 | 0.936 |
| S | 0.446 | 0.324 | 0.295 | 0.168 | 0.248 | 0.182 | 0.940 |
| weight 1.0 (best picks) | 0.328 | 0.225 | 0.197 | 0.235 | 0.144 | 0.098 | 0.966 |
| weight 0.5 | 0.344 | 0.296 | 0.284 | 0.081 | 0.196 | 0.161 | 0.881 |

The five time blocks agree (MSE drop 0.13-0.17). The coverage passes (0.938). The per-station table was not printed
because of a pandas join error after the summary (fixed; the second-step runs print it).

Diagnosis (not tuned on, read from the table): the GP gains much more on the bulk than on the mean square. The robust
spread falls by 28 % (MAD 0.168 -> 0.121, i.e. about 48 % in variance), the MSE by 15 %. The best picks (weight 1.0)
reach 0.235, the weight 0.5 picks only 0.081. So the hypothesis is that the MSE is held up by reading noise with heavy
tails that no smooth term can predict (the analyst pick error, larger at low weight), not by a non-smooth 3-D residual.
A single white level also gives the noisy picks the same weight as the good ones in the GP's fit. Lesson of M5
applied: the premise (how much of the constant term's MSE is irreducible reading noise) is measured before
concluding, as the "drop bound" of the second step.

## Second step (declared 9 November 2026, before running)
1. Diagnosis on the development days (190-196), both variants, the gate NOT re-read here: print the per-station table
   and the "drop bound" = 1 - mean white variance of the GP / MSE(const), the largest MSE drop any term can reach if
   the white part is irreducible reading noise.
2. One fixed variant, chosen now: `het`, the GP's white level depends on the SCSN pick weight class (four classes,
   levels fitted by marginal likelihood with the rest; `GP1.fit(X, y, g)`). Nothing else changes (kernel, folds,
   gross rule, stations rule, 1-D medium).
3. Confirmation on data never used: days 197-203 (2019-07-16 to 07-22), same fetch tool and rules, one run of `het`,
   with the gate exactly as declared (all readings, MSE drop >= 0.30, coverage 0.90-0.98).
4. Pivot rule: if step 3 fails, G-M2 is recorded as a negative result for analyst picks. If the drop bound of step 1
   is below 0.30, the gate is not reachable with this reading noise by any term, and that is reported to the owner as
   a property of the gate (the decision to change a gate is the owner's, not this thread's); the reading noise of
   our own detector (M1 onsets) is then the data on which G-M2 should be measured. If the bound is above 0.30 and the
   GP still fails, the residual is not smooth and the spec's scaling gate (3-D eikonal or neural travel times) is the
   owner's decision.

## Result, second step (9 November 2026)
Confirmation days 197-203: 4617 phase files, 1238 events, 28 stations by the same rule, SHA-256 `events.csv`
eaebabe5..., `picks.csv` 880250ec... (`run_conf_het_out.txt`).

| run | readings | MSE drop | cov +-2 sigma | MAD const -> gp | drop bound | weight 1.0 drop |
|:--|--:|--:|--:|--:|--:|--:|
| 1. development, plain (`run_plain_out.txt`) | 107 004 | 0.150 | 0.938 | 0.168 -> 0.121 | 0.316 | 0.235 |
| 1. development, het (`run_het_out.txt`) | 107 004 | 0.143 | 0.936 | 0.168 -> 0.117 | 0.311 | 0.247 |
| **3. confirmation, het** (`run_conf_het_out.txt`) | 43 371 | **0.191** | **0.942** | 0.160 -> 0.102 | 0.330 | 0.303 |

**G-M2 FAILS on the confirmation days** (0.191 < 0.30); coverage passes (0.942).

Reading, by the declared pivot rule:
- The drop bound is above 0.30 but only just (0.31-0.33): if the white part of the GP is irreducible pick noise, no
  term can lower the MSE by much more than 30 % on these picks. The gate is at the edge of what this reading noise
  allows, so this failure does not show that the residual is non-smooth.
- The GP does help the bulk a lot: the robust spread falls 28-36 % (about half of the variance), and on the best picks
  (weight 1.0) the MSE drop is 0.25 (development) and 0.30 (confirmation). The mean square is held by the heavy tail of
  the lower-weight picks (weight 0.5 and 0.8: 0.06-0.12).
- The heteroscedastic variant did not raise the pooled drop (0.143 against 0.150 on development); it re-sizes the
  tolerances (median sigma of P 0.13 s against 0.18 s) and keeps the coverage. The reading noise, not the fit, holds
  the MSE.
- Consequence for the plan (owner's decision, not taken here): (a) measure G-M2 on our own readings (M1 onsets on
  R1 or the fk sets), whose noise is the one M4 sees; (b) or restate the gate with a robust error (MAD) or on
  high-quality readings, declared before a new run on new data; (c) the 3-D scaling gate is NOT indicated by this
  result. Days 190-203 are now spent for M2 choices.

## Note on the gate (14 November 2026, supervision review)
G-M2 fails as declared and stays failed. Read with the noise floor measured above: the GP itself puts the reducible
share of the MSE at 0.31-0.33, so the 30 % drop asked by the gate was at the edge of what the pick noise allows, and
the 0.191 obtained is about 60 % of the reducible part. The failure says that the gate was set at the noise floor,
not that the residual has no smooth structure; it is not a pass either. A new gate (share of the reducible MSE, and
the comparison against static station terms, which are kept if the GP does not beat them) has to be declared before
new data are read.

