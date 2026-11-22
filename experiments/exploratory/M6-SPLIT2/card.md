# Experiment card: M6-SPLIT2. S splitting shared by all events of the station inside one joint fit

Stage: S2, run once 21 November 2026 (`run_out.txt`); **stopped at the 2 h limit before arm v-w; arm iii fails validity
(126 of 200) because of an implementation defect found after the run (below), so the run does not test the idea.**
Previously S1, 21 November 2026. Drafted 20 November 2026 (2a7c454, PR #75) before any code; the owner chose
"Restricción común a la estación" on 21 November 2026, 00:56 UTC. Code built and pushed with this card before the
run; thresholds, arms, seeds and predictions unchanged. A smoke test of one control and one pair on seed 999
(outside the declared seeds) only timed the code (arm iii about 20 s per loop). Exploratory, simulation only; G-M6a
is not run.

## Why
M6-SPLIT (PR #75) estimated phi_i and delta_i from the window of the test alone: validity failed (24 of 200 single-
source control scenes split, limit 12), AUC 0.885. Two hypotheses were left:
1. With delta <= 0.2 s and a 2 Hz pulse, phi_i and delta_i are poorly fixed by one or two events, and the one-event
   class spends that freedom to explain a pair. The mathematical formulation thread (F3, N4) gives the same reading:
   with a free fast axis per station the residual of an S pair falls from O(Delta^2) to O(Delta^4), and the network
   margin keeps only 0.11 / 0.23 / 0.37 / 0.68 of its value at delta 0.02 / 0.05 / 0.1 / 0.2 s.
2. The mismatched pulse shape is confounded with the splitting at this resolution.

## Model (to be built in `ttmt.v2.m6split`, reusing M6-SPLIT's code)
phi_i and delta_i (0 <= delta_i <= 0.2 s) belong to the station and are shared by ALL events seen at it inside the
same joint fit: the atoms of the test window AND NISO = 10 isolated events recorded at the same stations (separate
40 s windows, the same splitting, sources drawn in the region, SNR 12 per-block peak, kmax 1 each). In each round of
the fit, the splitting cost of station i is the sum over the test window and the isolated windows; the isolated
events' tensors and times are refitted with it. This holds in BOTH hypotheses of the "one against two" test (the two
events, and the one event with free arrivals), so it stays a constraint inside the fit, not a calibration beforehand.
The isolated windows are part of each scene and are never reused across scenes.

## Arms (same sources and noise; splitting from its own stream; new seeds: controls 20261124, pairs 20261125,
splitting 20261126, isolated events 20261127)
- (i) no splitting in the data, model without it;
- (ii) split data, model without it;
- (iii) split data, splitting shared with the isolated events (the model above);
- (iv) split data, splitting free per event (descriptive, as X1);
- (v) descriptive, hypothesis 2: arms (i) and M6-SPLIT's window-only constraint with MATCHED pulses (the kernel's
  Ricker, no stretch or phase), the rest as above.
Scenes: M6-CAL settings (8-station ring, homogeneous medium, mismatched pulses except in (v), SNR 12, kmax 4,
near_km 10, LOC+T-V). 200 single-source controls; near pairs (1.5-5 km), dtau 0.1, 0.25, 0.5, 1.0 s, 40 per cell.

## Criteria (unchanged from M6-SPLIT and the supervision note of 20 November, section 3)
- Validity, never relaxed: in (iii), at most 12 of 200 single-source split controls come out split at thr 23.8.
- Ambition, three zones on near-pair recall pooled over the 160 pairs, (i) - (iii): below 0.05 passes; 0.05 to 0.15
  inclusive is useful and the owner decides; more than 0.15 costs too much. Compared on exact counts.
- Route B: AUC of the pair gain in (iii) >= 0.90 (as in M6-SPLIT).
If all three pass, G-M6a runs with (iii) after ADR-0012.

## Predictions, written before any run
- (iii) splits fewer control scenes than M6-SPLIT's window-only constraint (24 of 200), because the isolated events
  fix phi_i and delta_i.
- From the formulation's table, the remaining recall loss in (iii) concentrates below 0.1-0.25 s.
- If (v) closes most of the gap between (i) and the window-only constraint, hypothesis 2 holds and pulse mismatch,
  not splitting, is the larger cost.

## Implementation notes (written before the run; no criterion changes)
- `ttmt.v2.m6split.adcg_loop(..., iso_windows=[(Y_w, q_w)])`: each isolated window is first read with one atom by
  `m6.adcg_loop` (kmax 1, no splitting); then the station's (phi_i, delta_i) is estimated on the isolated windows,
  and in every joint fit of the test window the splitting cost of station i is the sum over the test window and the
  isolated windows, whose tensors and times are refitted at each round. Each hypothesis keeps its own copy of the
  isolated events, and their residual is added to the RSS of both, so the pair gain compares the two hypotheses
  with the isolated events refitted under each one's splitting.
- Arm (v) is run as two sub-arms on matched pulses: v-i (as arm i) and v-w (M6-SPLIT's window-only constraint).
- Without isolated windows the code is that of M6-SPLIT.

    python experiments/exploratory/M6-SPLIT2/run.py > experiments/exploratory/M6-SPLIT2/run_out.txt

## Cost
CPU, cloud, 4 cores: about 3 to 5 hours (the isolated windows multiply the splitting step about 5 times).

## Result (`run_out.txt`; run once after 24e166c; the job hit the 2 h limit during arm v-w, which has no output)

| arm | control scenes split at 23.8 (of 200) [95 % CI] | near recall 0.1 / 0.25 / 0.5 / 1.0 s | pooled (of 160) | AUC | largest spurious gain median / max |
|:--|--:|:--|--:|--:|:--|
| i, no splitting | 2 [0.001, 0.036] | 0.850 / 0.825 / 0.975 / 1.000 | 0.912 (146) | 0.961 | -105 / 55 |
| ii, split data, model without | 48 [0.183, 0.305] | 0.525 / 0.500 / 0.600 / 0.850 | 0.619 (99) | 0.708 | -99 / 1753 |
| **iii, shared with 10 isolated events** | **126 [0.559, 0.697]** | 0.725 / 0.800 / 0.900 / 0.975 | **0.850 (136)** | **0.931** | **+115** / 1901 |
| iv, free per event | 51 [0.196, 0.321] | 0.700 / 0.700 / 0.825 / 1.000 | 0.806 (129) | 0.926 | -28 / 587 |
| v-i, matched pulses, no splitting | 0 [0.000, 0.018] | 0.750 / 0.900 / 0.975 / 1.000 | 0.906 (145) | 0.956 | -114 / -16 |
| v-w, matched pulses, window-only constraint | not run (time limit) | | | | |

Reading by the declared rule, arm iii: validity FAILS (126 > 12); recall (146 - 136) / 160 = 0.0625, zone "useful,
the owner decides"; AUC 0.931 passes. **But the validity number is not a property of the idea:**

**Defect (found after the run, `diag_iso_rss.py`, `diag_iso_rss_out.txt`, two control scenes on seed 999).** The RSS
of the isolated windows enters both hypotheses, and the isolated events are not refitted to convergence: their RSS
falls at every refit (158 552 -> 153 702 in one scene). The two-event trial refits them 3 rounds, the one-event
alternative fewer, so the trial's isolated RSS is lower by 59 to 241 in the scenes checked, which enters the pair
gain directly. That is why the median largest spurious gain per control scene is +115 in arm iii (about -100 in every
other arm), and why most control scenes split. The gain compares how far the nuisance fit converged, not the two
hypotheses.
Consequence: arm iii of this run does not test "splitting shared with isolated events". The other arms are unaffected
(no isolated windows). A fix (isolated events refitted to convergence before the loop and the same number of rounds
in both hypotheses, or a fixed point check) needs a new declared run; it is not done here.

**Note added 22 November 2026, after M6-SPLIT2b (review of PR #75).** The fix was run as M6-SPLIT2b (new seeds,
isolated fits converged, same rule in both hypotheses). The bias persists: median largest spurious gain +86.8
(here +115), 114 of 200 control scenes split. The convergence defect therefore explains at most part of this run's
failure. The statements above that the validity number "is not a property of the idea" and that the run "does not
test the idea" are too strong: they hold only for the share due to the defect. See `../M6-SPLIT2b/card.md`.

What the valid arms say (descriptive): splitting in the data costs today's model 46 split scenes and 0.29 of recall
(ii against i); free splitting per event (iv) recovers recall to 0.806 but not validity (51); with matched pulses and
no splitting (v-i) the numbers equal arm i (0 against 2 split scenes, 0.906 against 0.912), so pulse mismatch alone
does not explain the gap between i and ii (hypothesis 2 is not supported by v-i; v-w, the arm that would test it with
splitting, was not reached).

