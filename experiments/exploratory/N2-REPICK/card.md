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

## Result (16 November 2026, laptop, GPU; outputs `baseline_out.txt`, `repick_smoke_out.txt`, `repick_out.txt`, `repick.json`)

Checks passed in the smoke and in the full run (frozen objective reproduced; rebuilt events add up to the stored
signal; the picker of the machine gives the stored picks on the reference target). 120 targets x k in 20
realisations, 52 unique targets. Run time 140 s.

| | n | K (control) | S (after subtraction) | net S - K |
|---|---|---|---|---|
| all | 120 | 40 | 32 | **-8**, 95 % interval [-23, 7] |
| unique targets | 52 | 14 | 13 | -1 |
| separation 0.10-0.25 s | 25 | 14 | 5 | -9 |
| separation 0.25-0.50 s | 37 | 14 | 12 | -2 |
| separation 0.50-1.01 s | 58 | 12 | 15 | +3 |
| Ridgecrest | 66 | 38 | 27 | -11 |
| Hi-net | 54 | 2 | 5 | +3 |

**Reading by the declared rule: LOW zone (net -8, below 15; the whole interval is below 15).** The target is not
met. By the card, N-2 is not pursued as a picker aid. The rule is not changed.

### What the numbers say beyond the rule (measured, read after the run)

- With the partner removed exactly, 88 of the 120 targets still have fewer than 6 picked arrivals and 50 have none.
  Most targets "lost at the picker" are not lost because the partner hides them: the picker does not see them on
  their own. On Hi-net it is 49 of 54.
- The subtraction does add picks (53 targets gain, 25 lose, mean +1.6 against the stored picks), but rarely enough.

### A defect of the declared control, found on reading the result

K is not a clean control below 0.5 s. The matching tolerance is 0.5 s, so with the partner out of the truth a pick
made by the partner's energy counts for the target (mean picked arrivals at 0.10-0.25 s: 2.2 stored, 7.1 in K, 3.2
after subtraction). K therefore overstates what the target has of its own at short separations, and the net S - K
understates the subtraction there. This was a design error of this card, made before the run.

It does not change the reading: the number that a real "pick, subtract, pick again" chain could use is S itself
(the first-pass picks are consumed by the partner), and S is 32 of 120 (27 %; 13 of 52 unique targets) as a ceiling
with a perfect subtraction. That is an upper bound of 0.14 in T recall (32 of 223) before any loss from an estimated
subtraction, from association and from the false-alarm budget. Reported as the effect size for the owner's
decision (rule 36); no criterion was declared on S, so no "passes" is claimed for it.

### Consequence proposed (owner decides)

- Hypothesis now favoured, to be checked: the picker-stage loss in T is mainly a low-signal problem, not an overlap
  problem. Cheap check: per-station SNR of the 88 targets that stay below 6 picks against the 32 that recover.
- If it holds, the low-SNR source (N-1, signed stack) addresses most of this loss and N-2 at most a quarter of it.

## Addition after the independent review of b47d9e0 (16 November 2026)

The review reached this thread after the run of d15237e had been read. Two points taken, both declared here before
the rerun:

- Same-machine control. K comes from the stored picks and S from the picker of the laptop; a difference between the
  two pickers would enter S - K. `repick` now also picks the unmodified rebuilt record on the same machine and
  matches it with the partner out of the truth (K_local). Both nets are reported; the zone boundaries stay as the
  absolute counts declared above; if K_local and K differ by more than 2 the reading is taken on S - K_local.
  What is known beforehand: the first run gave K 40, S 32. S is deterministic, so only K_local is new.
- The sets S without K and K without S are printed next to the net.

Rerun output: `repick2_out.txt`, `repick2.json` (the first outputs are kept).

**Rerun (16 November 2026, 162 s, GPU).** K_local = 40, equal to K in every one of the 120 rows, so the two pickers
are the same and the net does not come from the machine. S 32, net S - K_local = -8, 95 % interval [-23, 7]: LOW
zone again. Sets: 12 targets reach 6 picks only after the subtraction, 20 reach it only in the control, 20 in both.

## Step 2: is it low signal? (declared on 16 November 2026, before running)

Question: do the targets that stay below 6 picked arrivals after the exact subtraction differ in signal level from
the 32 that recover? Same data, same 120 targets; `run.py snr lost.json repick.json`. Descriptive; no parameter is
tuned.

- Measure per target: `snr_med`, the median over channels of the peak of the target's own injected waveform over
  the robust sd (1.4826 MAD) of the record's noise. Secondary: number of channels with ratio at least 5.
- Statistic: AUC of `snr_med` between recovered (n = 32) and not recovered (n = 88), bootstrap over realisations
  (2000 resamples, seed 20261121).
- Reading on the point estimate, interval reported: AUC >= 0.75 the low-signal hypothesis is supported; 0.60 to
  0.75 partly (signal level matters but does not explain the loss alone); < 0.60 not supported, and another cause
  must be looked for (label width, station geometry, the other events of records with k > 2).
- Known before running: only the counts of the Result section. No signal level of any target has been computed.
- Limit: an association, not a cause; a peak ratio is not the network SNR of Appendix AB.

### Result of step 2 (16 November 2026, 122 s; `snr_out.txt`, `snr.json`)

| measure | recovered (32), median [q25, q75] | not recovered (88) | AUC, 95 % interval |
|---|---|---|---|
| `snr_med` | 6.93 [3.85, 15.23] | 3.63 [1.63, 5.71] | **0.718** [0.593, 0.865] |
| channels with ratio >= 5 | 21.5 [6, 32] | 7 [1, 13] | 0.739 [0.590, 0.913] |

**Reading by the declared rule: the low-signal hypothesis is PARTLY supported** (0.60 to 0.75; the interval reaches
both neighbouring zones). Signal level matters, and it does not explain the loss alone.

What it does not explain (measured, read after the run): the two sets have about the same signal level (median
`snr_med` 3.85 at Ridgecrest, 4.11 at Hi-net) but 27 of 66 recover at Ridgecrest and 5 of 54 at Hi-net. So on
Hi-net something other than the peak ratio keeps the picker from seeing a target that is alone. Hypotheses to
check, none tested: the picker's training domain (PhaseNet "original" is Californian), the frequency content of the
Hi-net short-period records, the remaining events of records with k > 2.

Consequence proposed for the owner: the claim "the picker loss in T is a low-signal problem" is not to be written;
what is measured is that it is mostly not an overlap problem (88 of 120 stay lost with the partner removed). The
low-SNR source (N-1) stays justified in part; a picker adapted to Hi-net is a second candidate cause, cheap to test
with the per-set split of any picker fine-tune.
