# Experiment card: M6-SPLIT. S-wave splitting as a constraint inside the joint fit of M6

Stage: S2, run once on 20 November 2026 (`run_out.txt`, 3540 s on 4 cores): **validity FAILS** (24 of 200), recall cost
at the boundary (0.150), AUC 0.885 below 0.90. Declared 20 November 2026 in 1342740, before any run (a sanity check on 2 control scenes of seed 999, outside the
declared seeds, only timed the code). Exploratory; no registered number changes. Simulation only; G-M6a is NOT run
(it also needs the weight fix of ADR-0012).

Origin: the director's request (meeting of 20 November 2026, as reported by the owner): splitting of the S wave
enters as a constraint inside the joint fit, not calibrated beforehand. Design and criterion fixed by the supervision
note `supervision/tras-reunion-director-2026-11-20.md`, section 3; this card only implements them.

## Question
If the S wave of the data is split, does a joint fit that estimates the splitting per station (inside both
hypotheses of the "one against two" test) keep the false-split rate of single-source scenes within the G-M6a limit,
and how much near-pair recovery does it cost against data without splitting?

## Model (`src/ttmt/v2/m6split.py`; `ttmt.v2.m6` unchanged apart from an optional return of the free-station times)
- Only the S block is split; the P of each event is one pulse.
- Per station i: fast axis phi_i in the plane normal to the ray, delay delta_i, 0 <= delta_i <= 0.2 s. S block of
  atom k: (Pf g) r(t - T) + (Ps g) r(t - T - delta_i), Pf = ef ef' + n n', Ps = es es'.
- phi_i and delta_i are of the station: common to all atoms of the window (mode "station"), estimated in the same
  joint fit (tensors closed form, times bounded search, splitting on a 15 degree x 9 delay grid refined by bounded
  searches), in BOTH hypotheses. In the one-event alternative the free P and S of the tested atom (free times and
  amplitudes, the class of A.3(a)) share phi_i, delta_i with the other atoms. delta = 0 is the old model (nested).
- The pair gain is computed with splitting in both hypotheses.

## Arms (same sources and noise in all; splitting parameters from their own stream, seed 20261123)
- (i) no splitting in the data, model without it (`m6.adcg_loop`);
- (ii) split S in the data (per station delta ~ U[0.02, 0.2] s, phi ~ U[0, pi), one draw per scene, common to all
  sources of the scene), model without it (today);
- (iii) the same data, model with the constraint (mode "station");
- (iv) the same data, splitting free per event (mode "event", descriptive, as X1).
Scenes: M6-CAL settings (8-station ring, homogeneous medium, mismatched Ricker pulses, SNR 12 per-block peak, kmax 4,
near_km 10, LOC+T-V). Controls: 200 single-source scenes (seed 20261121). Pairs: near stratum (1.5-5 km), dtau 0.1,
0.25, 0.5, 1.0 s, 40 per cell (seed 20261122). New seeds, never used.

## Criteria (declared before the run; from the supervision note, section 3)
- **Validity, never relaxed:** in (iii), control scenes with one split source that come out split (a spurious atom,
  `M6-BIAS.spurious`) at most 12 of 200 (the G-M6a limit). Read at the reference threshold 23.8, as the false-split
  count of G-M6a.
- **Ambition, three zones:** near-pair recovery (both sources hit, LOC+T-V, thr 23.8, pooled over the four
  separations, 160 pairs) in (iii) within 0.05 of (i): passes; 0.05 to 0.15 below: useful, the owner decides; more
  than 0.15 below: the constraint costs too much, and that is said.
- **Route B:** AUC of the pair gain as a score in (iii) >= 0.90: gains of the tested atoms matched to a true source
  in the near pairs against the largest spurious gain per control scene, every candidate kept (thr = -inf).
If all three pass, G-M6a runs with (iii) (after ADR-0012).
Reported, outside the rule: the same numbers for (i), (ii), (iv), the recall per separation, Clopper-Pearson
intervals of the split counts.

## Hypothesis
From X1 and the M6-THEORY confounders (split S gave 7-14 false pairs in 20): (ii) splits clearly more than (i);
(iii) brings it back near (i) at a small recall cost; (iv) loses recall (free splitting absorbs a close pair's S).

## Cost and procedure
CPU, cloud, 4 cores, about 2 hours (about 5 s per control scene and model in (iii)).

    python experiments/exploratory/M6-SPLIT/run.py > experiments/exploratory/M6-SPLIT/run_out.txt

## Result (`run_out.txt`; run once, after the declaration was pushed in 1342740)

| arm | control scenes split at 23.8 (of 200) [95 % CI] | near-pair recall 0.1 / 0.25 / 0.5 / 1.0 s | pooled (of 160) | AUC | true gains median | largest spurious gain median / max |
|:--|--:|:--|--:|--:|--:|:--|
| (i) no splitting | 3 [0.003, 0.043] | 0.775 / 0.850 / 0.975 / 0.975 | 0.894 (143) | 0.946 | 1481 | -115 / 53 |
| (ii) split data, model without | 39 [0.142, 0.257] | 0.425 / 0.500 / 0.650 / 0.875 | 0.613 (98) | 0.740 | 454 | -114 / 574 |
| **(iii) split data, constraint per station** | **24 [0.078, 0.173]** | 0.475 / 0.725 / 0.825 / 0.950 | **0.744 (119)** | **0.885** | 1294 | -88 / 367 |
| (iv) split data, free per event | 41 [0.151, 0.268] | 0.450 / 0.700 / 0.825 / 0.975 | 0.738 (118) | 0.892 | 1416 | -41 / 316 |

Reading by the declared rule:
- **Validity fails:** 24 of 200 control scenes split in (iii), above the limit of 12. G-M6a does not run with (iii).
- **Recall cost:** (i) - (iii) = (143 - 119) / 160 = 0.150 exactly, the boundary between "useful, the owner decides"
  (0.05 to 0.15) and "costs too much" (more than 0.15). By the card's words 0.150 is not "more than 0.15", so the
  zone is "useful". The script printed "COSTS TOO MUCH" because it compared floating-point means
  (0.89375 - 0.74375 > 0.15 in binary); this is a defect of the script's printout, not a change of the rule, and
  the zone read here is the one of the card. After the run the script compares counts exactly (`run.py`, decision line); `run_out.txt` is
  left as produced. The point does not change the outcome: validity already fails.
- **Route B:** AUC 0.885 < 0.90, fails.

What is measured, outside the rule:
- Splitting in the data is a real threat to today's M6 (arm ii): false splits go from 3 to 39 of 200, near-pair
  recall from 0.894 to 0.613, AUC from 0.946 to 0.740.
- The constraint recovers part of it on every number (39 -> 24 splits, 0.613 -> 0.744 recall, 0.740 -> 0.885 AUC),
  but not to the level of data without splitting.
- Free splitting per event (iv) recovers the same recall as (iii) but no false splits (41), which supports the
  station constraint over free splitting, as X1 predicted.
- Hypotheses, not tested: (1) delta <= 0.2 s is below the resolution of the 2 Hz pulse (period 0.5 s), so phi_i and
  delta_i are poorly determined from one window of one or two events (the sanity check of two scenes showed delay
  errors up to 0.13 s); estimating them from many isolated events per station (the "calibrated per station" of the
  supervision design of 14 November) would fix them better; (2) the mismatched pulse shape and the splitting are
  confounded at this resolution. Test that separates them: the same four arms with matched pulses.
- Limits: simulation (8-station ring, homogeneous medium, Ricker pulses, Gaussian noise); one split draw per scene.
