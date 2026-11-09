# Experiment card: M6-BIAS. Biases of the one-event alternative of the resolvability test

Stage: S1 (done, 9 November 2026, cloud). Origin: theory review T5.2, first pass
(`teoria-fase5/revision-A1-A3a-M6.md`, items b1-b3). Exploratory, simulation; no registered number changes.

## Question
Does the one-event alternative of the first M6 version (`free_one_event_rss`: block times bounded to +-8 samples, P
then S, the other atoms fixed) inflate the RSS gain of a false split when other events are present, and does the
joint alternative (`free_one_event_joint`: +-16 samples, P and S fitted together with their Gram, the other atoms
refitted) remove it?

## Design
As M6-JOINT (ring, homogeneous medium, Ricker, synthetic picker, Gaussian noise, SNR 12), seed 20261109. Null
scenes: C1 one source; C3 the target, a neighbour 6-14 km away and 0.5-1.5 s later, an isolated source 3 s later.
30 calibration scenes each, `adcg_loop` with every candidate kept; 30 fresh scenes each for false splits; 20 pairs per
cell at 0.25 and 0.5 s. Scored with LOC+T (as M6-JOINT; LOC+T-V comes in M6-CAL).

## Result (`run_out.txt`)
| control | legacy gain p95 | joint gain p95 | legacy max | joint max |
|:--|--:|--:|--:|--:|
| C1, one source | -52.4 | -54.5 | -25.5 | -25.5 |
| C3, with neighbours | **161.5** | **-23.4** | 297.9 | 1.8 |

- The bias exists, and it is the one the review predicted: with neighbours the legacy alternative inflates the gain of
  a false split (p95 +185 against the joint one); with one source there is almost none (+2). A threshold calibrated on
  one source (M6-JOINT) is too low as soon as there are other events.
- Fresh controls with neighbours: legacy at its old threshold (10.8) splits 4 of 30 scenes; joint at the recalibrated
  threshold (-23.4, the C3 p95) splits 2 of 30. One source: 0 of 30 for both.
- Near-pair recall (LOC+T) is unchanged by the switch (0.95 at 0.25 s, 0.85 at 0.5 s for both; greedy 0.60, 0.55).
- `free="joint"` is now the default of `adcg_loop`; `free="legacy"` remains for this measurement.

## Limits
kmax = sources + 1, exact pulses, LOC+T scoring, per-test level. M6-CAL addresses all four (CR-29, CR-30, CR-37, T5.2 d).

    python experiments/exploratory/M6-BIAS/run.py [NCAL] [NTEST] [NPAIR] > experiments/exploratory/M6-BIAS/run_out.txt
