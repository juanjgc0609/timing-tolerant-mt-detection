# Experiment card: M6-SPLIT2. S splitting shared by all events of the station inside one joint fit

Stage: S0, draft, 20 November 2026. **Not run, and not to be run until the owner decides** (asked by the supervisor
after M6-SPLIT). No code yet. Exploratory, simulation only; G-M6a is not run.

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
- (vi) descriptive, hypothesis 2: arms (i) and M6-SPLIT's window-only constraint with MATCHED pulses (the kernel's
  Ricker, no stretch or phase), the rest as above.
Scenes: M6-CAL settings (8-station ring, homogeneous medium, mismatched pulses except in (vi), SNR 12, kmax 4,
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
- If (vi) closes most of the gap between (i) and the window-only constraint, hypothesis 2 holds and pulse mismatch,
  not splitting, is the larger cost.

## Cost
CPU, cloud, 4 cores: about 3 to 5 hours (the isolated windows multiply the splitting step about 5 times).
