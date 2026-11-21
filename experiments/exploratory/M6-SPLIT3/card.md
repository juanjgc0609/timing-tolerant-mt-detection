# Experiment card: M6-SPLIT3. Splitting estimated once per station from isolated events, then held fixed

Stage: S0, declared 21 November 2026 before any code or run. Asked by the thread that finishes the supervision items
(arm (v) of the M6-SPLIT discussion; hypothesis 1 of M6-SPLIT and the proposal of the formulation thread); the owner
asked that everything that does not need the laptop be run. Exploratory, simulation only; G-M6a is not run.

Difference with M6-SPLIT2: there phi_i and delta_i are shared with the isolated events INSIDE every joint fit (both
hypotheses re-estimate them); here they are estimated ONCE per station from the isolated events of the scene's long
window and then held FIXED for the pair test and for the control. This is closer to a calibration before the fit,
which the director's request of 20 November did not ask for; it is run as the comparison that says how much of the
loss of M6-SPLIT comes from estimating the splitting with the pair itself.

## Conditions (declared by the requester, not changed after)
1. phi_i, delta_i (0 <= delta_i <= 0.2 s) estimated once per station with all NISO = 10 isolated events of the
   scene (separate windows at the same stations, the same splitting draw, new sources and noise; each read with one
   atom, then the splitting fitted jointly over the isolated windows with their tensors and times refitted). Then
   fixed: neither hypothesis of the test re-estimates them, for the pairs and for the controls.
2. Validity, not relaxed: at most 12 of 200 single-source split control scenes split at thr 23.8.
3. Near-pair recall (160 pairs, pooled, both sources hit, LOC+T-V, thr 23.8) against (i), on exact counts: below
   0.05 passes; 0.05 to 0.15 is close and the owner decides; more than 0.15 fails. AUC of the pair gain >= 0.90.
4. Prediction (formulation thread), written before the run: this arm recovers almost all that M6-SPLIT's arm (iii)
   loses at 0.1 s (there 0.475 against 0.775 for arm i).

## Arms and seeds
Arms (i) (no splitting, `m6.adcg_loop`), (ii) (split data, model without it) and (fixed) (this card), M6-CAL
settings (mismatched pulses, SNR 12, kmax 4, near_km 10). New seeds, not used by M6-SPLIT (20261121-23),
M6-SPLIT2 (20261124-27) or M6-SPURIOUS (20261127-28): controls 20261131, pairs 20261132, splitting 20261133,
isolated events 20261134. 200 controls, 40 near pairs per dtau in 0.1, 0.25, 0.5, 1.0 s.

## Code and cost
`ttmt.v2.m6split.adcg_loop(..., mode="fixed", iso_windows=...)`; runs after M6-SPLIT2 finishes (same 4 cores).
About 1 to 2 hours.
