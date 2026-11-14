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
  Diagnosis by the review role (9 November 2026): `label_voronoi` follows ADR-0008; the flaw is in the control. A
  merged atom midway between the two events, copied at +-dtau/2 with the cell's TRUE dtau (oracle information), lands
  on both true origins whenever dtau >= s/10, which matches the failing cells. Proposed fix (registered metric
  untouched): draw the copy offset from a fixed distribution independent of the cell, e.g. U[0.05, 0.5] s. Pending
  the owner's approval; the control is not changed here until then.

## Step 2 (declared 9 November 2026, before running)
Calibrate the threshold on the operating loop: on the same 200 + 200 calibration scenes, find by bisection on thr the
smallest value at which the loop as run (`adcg_loop`, kmax 4, joint) splits at most 5 % of the scenes in BOTH
controls (one source, with neighbours; mismatched). Then the fresh controls and the pairs of step 1, unchanged, at
that threshold. Pass for the threshold: the with-neighbours CI upper bound <= 0.10.

## Result, step 2 (`run_op_out.txt`, 3824 s)
Bisection on the operating loop: thr = 48.5 (C1 0.005, C3 0.050 of the calibration scenes split). The C3 rate stays
at 0.040 even at thr = 400.

| fresh control | scenes split | 95 % CI |
|:--|--:|:--|
| one source, mismatched | 1 / 200 | [0.000, 0.028] |
| with neighbours, mismatched | **21 / 200** | **[0.066, 0.156]** |
| one source, exact | 0 / 100 | [0.000, 0.036] |
| with neighbours, exact | 6 / 100 | [0.022, 0.126] |

Pair recall is identical to step 1 (the threshold does not change which true pairs are kept). **Step 2 fails its
pass rule** (with-neighbour upper bound 0.156 > 0.10).

Diagnosis (`diag_untested.py`, the 21 fresh scenes): about half of the false splits are atoms that the test never
sees. They lie 4-5 km from a true source (the false-split rule counts atoms within 5 km of a source), but 5.4-8.7 km
from the nearest accepted atom, outside the test's `near_km` = 5. An atom 5 km from the source can be 10 km from the
source's own atom. The plateau at thr = 400 is these atoms. The rest are tested atoms with large gains (50-3700), at
the overlapping neighbour.

## Step 3 (declared 9 November 2026, before running; the last attempt on this line by rule 1)
The same operating-loop calibration and fresh controls with `near_km` = 10 (twice the LOC+T radius: every atom that
can be a false split of a source is tested against that source's atom). Pass rule unchanged (with-neighbour upper
bound <= 0.10). If it fails, the resolvability test with an RSS threshold is recorded as not calibratable in scenes
with neighbours in this form, and the next route is the owner's decision (a per-scene statistic, or the e-value route
of M7).

## Result, step 3 (`run_op10_out.txt`, near_km = 10, 4215 s). PASSES the pass rule of the threshold
Operating-loop calibration: thr = 23.8 (C1 0.010, C3 0.040 of the calibration scenes split).

| fresh control | scenes split | 95 % CI |
|:--|--:|:--|
| one source, mismatched | 6 / 200 | [0.011, 0.064] |
| with neighbours, mismatched | **10 / 200** | **[0.024, 0.090]** |
| one source, exact | 0 / 100 | [0.000, 0.036] |
| with neighbours, exact | 1 / 100 | [0.000, 0.054] |

| stratum | dtau | greedy (5 km) | greedy0 | M6 (near_km 10, thr 23.8) |
|:--|--:|--:|--:|--:|
| near | 0.10 | 0.20 | 0.80 | 0.65 |
| near | 0.25 | 0.50 | 0.95 | 0.80 |
| near | 0.50 | 0.50 | 1.00 | 1.00 |
| near | 1.00 | 0.65 | 1.00 | 1.00 |
| far | 0.10-1.00 | 0.85-1.00 | 0.85-1.00 | 0.85-1.00 |

Reading:
- With the test radius at twice the LOC+T radius and the threshold calibrated on the loop that is run, the false-split
  rate with neighbours is 5 % (upper bound 9 %), and with one source 3 %. The radius was the cause of the excess
  (step 2's diagnosis), not the threshold.
- M6 keeps its near-pair gain over the greedy loop with the exclusion (+0.45, +0.30, +0.50, +0.35, 20 pairs per cell,
  LOC+T-V). The duplicator caveat of step 1 still applies to differences of about 0.15 or less until the control is
  fixed (pending the owner); these gains are larger.
- Far pairs: 0.85 at 0.25 s against 0.90 in step 2 (one pair of 20); not read as a change.
- Still simulation: ring, homogeneous medium, Gaussian noise, synthetic picker. Before G-M6a is read, CR-30 also asks
  for calibration on isolated events of the generator v2 and on real isolated M2-3 events, and the theory items of
  T5.2 (c, d) and T5.6 stay open. Settings carried forward: `near_km` = 10, joint alternative, per-window level,
  operating-loop calibration, kmax 4, LOC+T-V.

## Rerun on fresh seeds, settings frozen (14 November 2026, supervision review; `run_fresh.py`, `run_fresh_out.txt`)
Reason: the "fresh" controls above were read in the three steps, and `near_km` 10 was chosen by looking at their 21
failures, so the 5 % of step 3 is not a held-out figure. Same loop, `near_km` 10, thr 23.8, mismatched pulses, new
seed 20261114, 300 scenes per control, nothing calibrated (656 s on 4 cores).

| fresh-seed control | scenes split | 95 % CI | step 3 |
|:--|--:|:--|:--|
| one source, mismatched | 4 / 300 | [0.004, 0.034] | 6 / 200 |
| with neighbours, mismatched | **25 / 300** | **[0.055, 0.121]** | 10 / 200, [0.024, 0.090] |

- **With neighbours the pass rule of step 3 (upper bound <= 0.10) is not met on fresh seeds.** Both runs together:
  35 of 500, 7 %. The figure to quote is "7 % of scenes with neighbours, between 5 and 10 %", not "5 %, at most 9 %".
- This does not decide the registered endpoint, which is recall at 30 false detections per hour on the exact FROC
  (ADR-0008 item 4), not a fixed false-split level: what matters there is that spurious atoms rank below true ones.
  The pair gain should enter the FROC as the score of the second atom; a fixed level is a development aid.
- **The duplicator column of the pair tables is descriptive only.** In cells made only of pairs and with no false
  budget, a merged detection copied at +-h is credited with both events whenever h is between dtau/6 and 1.5 dtau
  (64 % of the draws of U[0.05, 0.5] s at dtau 0.25 s, 89 % at 0.5 s), for any offset rule, so the control cannot pass
  there and the proposed fixed-offset fix would not change that. At the registered operating point the metric holds:
  the duplicator of frozen PyOcto on the tuning partition loses (recall in T 0.07-0.12 against 0.206, in C 0.11
  against 0.31, for h from 0.05 to 0.5 s; baseline reproduced exactly). ADR-0008 item 3 stays as written and is run
  only at the registered operating point. The step-1 caveat on differences of 0.15 or less no longer rests on a
  suspected metric flaw; it stays as a sample-size caveat (20 pairs per cell).

