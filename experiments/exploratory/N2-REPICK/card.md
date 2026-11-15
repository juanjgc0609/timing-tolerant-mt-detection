# N2-REPICK: pick, subtract, pick again (ideal case)

Declared on 15 November 2026, before the subtraction stage is run. Exploratory; it changes no registered number.
Approved by the owner on 15 November 2026 as the cheap pre-test of proposal N-2 ("siii apruebo las 3").

## Question

Frozen PhaseNet+PyOcto loses 120 of the 223 T targets of the tuning partition at the picker stage (STAGE-LOSS rule:
fewer than `as_min_picks` = 6 of the target's own arrivals have a pick within 0.5 s). If the waveform of the pair
partner is removed exactly and the same picker reads the residual, how many of those targets reach 6 picked arrivals?

This is the ceiling of the "pick, subtract, pick again" candidate source (A3 in the completeness design): the
subtraction uses the injected waveform itself, which no real method has. A low count closes the route cheaply; a
high count only says that building a real subtraction is worth the effort.

## Data

- Tuning partition `b60d99b3d0bbbbfa` only: semisynthetic sets `ss_ridgecrest_2019_egf_v1` and `ss_hinet_2012_egf_v1`,
  tune split, 50 records (realisation x k). The same partition on which the comparators were frozen.
- No read-partition file, no reserved confirmatory window (PR #50), no spent evaluation segment, no ledger.
- Use to be noted in the data-usage ledger when it exists (task T-4): one more exploratory reading of the tuning
  partition, with no parameter of any method tuned on it.

## Design

`run.py`, two stages.

1. `baseline` (manifests only). Reproduces the frozen objective (stops if T 0.2063 / C 0.3115 are not recovered
   within 5e-4), lists the T targets with stage "picker", and computes the **control K**: the stored picks matched
   again with the partner's arrivals removed from the truth, with no subtraction. K measures how much of the
   "picker" loss is only the one-pick-per-arrival assignment (a pick that exists but was given to the partner).
2. `repick` (set files and real record, on the laptop). Rebuilds each injected event on its own with the anchored
   seed, checks that the events add up to the stored signal (relative error at most 1e-3 of its peak; stops
   otherwise), subtracts the partner, runs PhaseNet with the frozen weights, threshold (0.0638) and peak distance
   (0.489 s), and matches as in K. **S** = targets that reach 6 picked arrivals after subtraction.

Checks that stop the run: frozen objective not reproduced; rebuilt realisation not equal to the stored one; the
picker of the machine differing from the stored picks by more than one arrival on the reference target.

Order: `repick --limit 1` first (one record, no reading printed), then the full run.

## Outcome and unit

- Unit: target x k among the 120 lost at the picker. Reported also by unique target and by separation
  (0.1-0.25, 0.25-0.5, 0.5-1 s).
- Main number: **net = S - K** (targets recovered by the subtraction beyond what the assignment alone explains).
- Interval: 95 % bootstrap over realisations (labels), 10 000 resamples, seed 20261120.

## What was already seen when this card was written

The `baseline` stage was smoke-run in the cloud before this card: frozen objective reproduced, 120 lost at the
picker, none without a partner, and **K = 40 of 120**. So the control is known and the target below was fixed
knowing it. S has not been computed anywhere; no waveform has been touched.

Consequence already worth recording (measured): a third of the "picker" loss of STAGE-LOSS is an assignment
effect, not missing picks. The targets truly without enough picks are about 80.

## Reading (rule 36: an ambition target, three zones)

| net S - K | reading |
|---|---|
| >= 40 (half of the 80 the control does not recover) | passes the target: N-2 goes ahead to a real (estimated) subtraction |
| 15 to 39 | target not met, useful smaller effect: reported with the interval; the owner decides (continue, restrict to the separations where it works, or drop) |
| < 15 | low: the ideal subtraction does not give the picker what it lacks; N-2 is not pursued as a picker aid |

The zone is read on the point estimate; the interval is reported next to it and, if it crosses a boundary, that is
said. The thresholds are not changed after the run.

## Limits of the evidence

- Ideal subtraction: an upper bound. A real subtraction needs the partner detected and its waveform estimated.
- Reaching 6 picked arrivals is necessary for PyOcto to form the event, not sufficient to detect it within the
  budget of 30 false alarms per hour. No recall is claimed from this test.
- Only the partner is subtracted; in records with k > 2 the other events stay.
- One picker (PhaseNet, original weights) at one frozen operating point; simulation with real noise and empirical
  Green's functions, not real overlapping pairs.

## Result

Not run yet.
