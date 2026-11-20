# Experiment card: M6-SPLIT. S-wave splitting as a constraint inside the joint fit of M6

Stage: S0, declared 20 November 2026, before any run (a sanity check on 2 control scenes of seed 999, outside the
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

## Result
Pending.
