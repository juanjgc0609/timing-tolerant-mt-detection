# Experiment card: M6A-REALNOISE. Gate G-M6a with real noise (first step)

Stage: S0 (declared 9 November 2026, before any real-noise run). Plan: T4.23 (M6 before M2, owner's approval of
8 November 2026). Exploratory; no registered number changes.

**Provisional (9 November 2026).** The threshold of the resolvability test rests on Corollary A.1 (z_e^2) and
Theorem A.3(a) of Appendix AB, whose independent review (plan T5.2) is not done, and the theory of the M6 test itself
(T5.6) is not written. Whatever this gate gives is not cited as a result until both are closed.

## Question
With real continuous noise, and the resolvability threshold calibrated on that noise, does the M6 loop
(`ttmt.v2.m6.adcg_loop`: joint refit + "one against two" test) beat the greedy v2 loop with the 5 km exclusion by
at least 0.10 of pair recall on pairs at 0.1-1 s whose K_mec/E is in the upper tercile (gate G-M6a of
`especificacion-modelo-v2.md`, M6)? And does it keep the false-split rate of M6-JOINT (0 of 40) when the noise has
real transients and non-Gaussian tails?

## Why this step, and what it does not test
M6-JOINT (PR #39) passed in simulation with Gaussian noise. The two things real noise can break are the threshold
(calibrated as z_e^2 + the 95th percentile of the spurious-split gain on a single-source control: real noise has
heavier tails and transients, so the control's gain may be larger or less stable) and the reading of weak blocks.
This step changes only the noise; the signal model, the medium and the picker stay those of M6-JOINT, so a change in
the result is the noise's. Declared limits: homogeneous medium, Ricker P/S pulses without coda, a synthetic picker
whose onset error is N(0, 0.02 s) and that has no false maxima (real noise does not reach the picker), kmax 2, SNR
fixed at 12. The full G-M6a on the fk semi-synthetic sets (real geometry, fk waveforms, PhaseNet onsets, known
K_mec) is the second step and needs an engine for those networks.

## Data
Task 2M noise day 2010-12-31 (Alto Tiberina, IV, 10 stations; `code/t2m_data.py`, QC qc1; on the owner's laptop
`C:\Users\juanj\ttmt\t2m_noise\101231`), 40 s eligible slots (no catalogued event within 60 s), each slot used once,
in a fixed permuted order (rng [20261109, 9]). Stations: the 8 nearest to the centroid of the day's set (real x, y;
z = 0). This day's slots were used by Task 2M and by M5 (as noise for the e-value backgrounds); here they are only
background for injected sources, and no M5 choice is made, so nothing of M5 is reused or spent.

## Design
- Two arms on identical sources (sources and picker from one random stream, noise from another): REAL (the slots
  above through the 2 Hz matched filter) and GEN (Gaussian raw noise as M6-JOINT, the control that isolates the
  effect of the real noise).
- Pairs: dtau in {0.1, 0.25, 0.5, 1.0} s, epicentral separation strata "far" (6-14 km) and "near" (1.5-5 km),
  NREP = 30 pairs per cell (240 pairs per arm). Each source scaled so its median strong-block peak is 12 x the
  robust sd of that scene's filtered noise.
- K_mec/E: `ttmt.synth.inject.kmec_over_e` on the scaled block vectors of the two sources; terciles over the 240
  pairs (identical in both arms).
- Methods, kmax 2: `greedy` (v2_loop, 5 km exclusion), `greedy0` (no exclusion, the review's control column),
  `m6-none` (no test), `m6` (test at thr calibrated per arm on 40 single-source scenes, separate seeds and slots).
- Fresh single-source control (40 scenes, new seeds and slots): false splits (a second atom within 5 km and 1.5 s).

## Gate (declared before the run)
REAL arm, upper K_mec/E tercile, both strata and the four separations pooled:
pair recall(m6) - pair recall(greedy) >= 0.10, AND the 2.5th percentile of its paired bootstrap (B = 2000,
resampling pairs) > 0, AND m6 false splits on the fresh control <= 2 of 40.
Reported but outside the gate: the same on GEN, per stratum and separation, the lower terciles, m6 - greedy0.

## Hypothesis and refutation
Expected from M6-JOINT: the gain is in the near stratum (greedy loses pairs inside the exclusion), so the pooled gain
should be about half of the near gain, near +0.2 to +0.3. Refuted if the pooled gain is below 0.10 or its interval
reaches 0. A FAIL with GEN passing points at the noise (look first at the threshold and the control's gain
distribution, then at the false splits); a FAIL in both arms points at the geometry (real station layout instead of
the ring). Pivot rule (declared): if the threshold is the cause, recalibrate it with time-slides of the real noise
(the spec's penalty) on a second noise day, once; if it still fails, record the negative result and go to step 2
(fk sets) before changing M6.

## Cost and procedure
CPU only. Cloud smoke run with fake t3 noise (1 pair per cell, 3 controls): 49 s, both arms. Laptop: a small run
(5 pairs per cell, 10 controls) first, then the full run (30, 40; about 20-30 min on one core).

    python experiments/exploratory/M6A-REALNOISE/run.py NOISE_DIR 5 10 > experiments/exploratory/M6A-REALNOISE/run_small_out.txt
    python experiments/exploratory/M6A-REALNOISE/run.py NOISE_DIR > experiments/exploratory/M6A-REALNOISE/run_out.txt

## Result
Pending (laptop run).
