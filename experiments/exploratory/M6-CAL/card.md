# Experiment card: M6-CAL. Threshold of the resolvability test and M6 pair recall under LOC+T-V

Stage: S1 (9 November 2026, cloud). Origin: review T5.2 (`teoria-fase5/revision-A1-A3a-M6.md`) and adversarial
pass 2 (CR-29, CR-30, CR-37). Exploratory, simulation; it comes before any reading of G-M6a.

## Question
With the project's hit rule (LOC+T-V), kmax 4, mismatched data pulses and the joint one-event alternative, (i) what
false-split rate does the resolvability test have at a threshold calibrated on controls (one source, and with
neighbours), and (ii) how does M6's pair recall compare with the greedy loop with and without the 5 km exclusion?

## Design
See the docstring of `run.py`. Ring, homogeneous medium, Gaussian noise, SNR 12, synthetic picker; pulses mismatched
per arrival (centre frequency x U[0.7, 1.4], phase U[-45, 45] deg, as code/gen_realistic.py). Threshold: the larger
of the 95th percentiles of the per-scene largest spurious gain on 200 one-source and 200 with-neighbour scenes
(every candidate kept). Fresh controls: 200 + 200 mismatched, 100 + 100 exact. Pairs: 20 per cell, LOC+T-V, with the
duplicator of the greedy loop as the metric check.

## Result, step 1 (`run_out.txt`, 1063 s on 4 cores)
Threshold 5.3 (C1 p95 -15.1, C3 p95 +5.3).

| fresh control | scenes split | 95 % CI |
|:--|--:|:--|
| one source, mismatched | 8 / 200 | [0.017, 0.077] |
| with neighbours, mismatched | **30 / 200** | **[0.104, 0.207]** |
| one source, exact | 1 / 100 | [0.000, 0.054] |
| with neighbours, exact | 6 / 100 | [0.022, 0.126] |

The greedy loop without the exclusion splits every scene (600 spurious atoms in 200 one-source scenes).

| stratum | dtau | greedy (5 km) | greedy0 | M6 | duplicator of greedy |
|:--|--:|--:|--:|--:|--:|
| near | 0.10 | 0.20 | 0.80 | 0.65 | 0.20 |
| near | 0.25 | 0.50 | 0.95 | 0.80 | 0.55 |
| near | 0.50 | 0.50 | 1.00 | 1.00 | 0.65 |
| near | 1.00 | 0.65 | 1.00 | 1.00 | 0.55 |
| far | 0.10-1.00 | 0.85-1.00 | 0.85-1.00 | 0.90-1.00 | 0.70-1.00 |

Reading:
- **The calibration does not hold with neighbours: 15 % of the scenes are split, against the 5 % aimed at.** Cause
  (from the design, not tuned): the control is run with every candidate kept, but at a finite threshold a rejected
  candidate excludes its nodes and the loop tries another (up to 3 kmax tries), so the operating loop makes more and
  different tests per window than the calibration loop. The level must be calibrated on the procedure that is run.
- M6 against the greedy loop with the exclusion, near pairs: +0.45, +0.30, +0.50, +0.35. The greedy loop without the
  exclusion recalls more at 0.1-0.25 s (0.80 and 0.95 against 0.65 and 0.80) but splits every single event.
- **Metric check fails in four cells:** the duplicator of the greedy loop gets MORE pair recall than the greedy loop
  (near 0.25: 0.55 vs 0.50; near 0.5: 0.65 vs 0.50; far 0.25: 0.90 vs 0.85). With LOC+T-V a duplicate should never
  help (ADR-0008 item 3). Not diagnosed here: the copies here have distinct scores (acceptance order) while
  `ttmt.stats.duplicator` gives them equal ones, and the near pairs have 1.5-5 km of epicentral separation, outside
  the same-epicentre case of `tests/test_stats.py`. Sent to the review and critique roles; until settled, recall
  differences of about 0.15 or less in near cells are not read as gains.

## Step 2 (declared 9 November 2026, before running)
Calibrate the threshold on the operating loop: on the same 200 + 200 calibration scenes, find by bisection on thr the
smallest value at which the loop as run (`adcg_loop`, kmax 4, joint) splits at most 5 % of the scenes in BOTH
controls (one source, with neighbours; mismatched). Then the fresh controls and the pairs of step 1, unchanged, at
that threshold. Pass for the threshold: the with-neighbours CI upper bound <= 0.10.
