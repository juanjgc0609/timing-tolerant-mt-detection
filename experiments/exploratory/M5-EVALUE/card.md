# Experiment card: M5-EVALUE. Gate G-M5c (calibration) of the mixture e-value S_mix

Stage: S3 (done, 8 November 2026). G-M5c FAILS on both seeds; design declared before the run (below). Plan: M5 pivot of ADR-0010 (keep the cap
6 for ranking, retire ACAP, move transient robustness to the decision step). Specification
`especificacion-modelo-v2.md` M5 "Decidir"; theory `teoria-fase5.md` D.1 and D.2. Exploratory, real noise; no
registered number changes.

## Question
On real Task 2M noise, with a polynomial-tailed background estimated per station and noise block, is the null rate
of S_mix >= log(1/p) at most p, for p in {1e-2, 1e-3}?

## What is built
`ttmt.v2.m5` (tests `tests/test_m5.py`):
- `evalue_mix(v, A, logf, eps, draws, s0, gamma)`: S_mix = log mean_j exp(sum_k l_k(m_j)), with
  l_k(m) = log[eps_k f_k(v_k) + (1 - eps_k) N(v_k; A_k m, s_k^2 I)] - log f_k(v_k) and s_k = s0_k + gamma |A_k m|.
  Several backgrounds share the model term in one call. CPU: about 1 ms per hypothesis with J = 2000 and 16 blocks,
  so the GPU the specification foresaw is not needed (Ockham).
- `TBackground` (isotropic Student t, ML), `GaussBackground` (ablation) and `RadialBackground` (below).
- `prior_draws`: the smoothed calibration distribution pi.

## Finding before the real run: the specification's t does not fit the registered reader
The registered reader takes the energy peak within +-4 samples of the predicted sample. A read vector in noise is the
maximum of nine, so its norm has almost no mass near 0 and a concentrated bulk. An isotropic t (or a generalised-t
radial law, tried and dropped: its shape parameter ran to its bound) cannot fit that bulk and a heavy tail at once.
Measured, not assumed (lesson of M5-ACAP):
- Synthetic peak reads (`test_radial_background_fits_peak_reads`): t PIT deviation > 0.05 per decile; the radial
  background < 0.01, tail ratio at 1e-2 within 0.7-1.4.
- Pipeline smoke test on fake i.i.d. t4 noise in the container (12 parts of 4 networks, 4800 hypotheses; not a
  result, a check of the code): the t background gives a null rate of 6.0e-3 at p = 1e-3 (pi CAL) and 9.8e-3 (pi LOW),
  the gaussian 7.7e-3 and 1.6e-2; the radial background 0 and 4.2e-4. So with the specification's t, Theorem D.2's
  premise (f exact) fails because of the reader, even without transients.
- Duplicates matter: many lattice hypotheses read the same peak sample. Without removing them the Hill index ran to
  10^4 and S_mix reached 4000 on a well-behaved null. The fit uses distinct peaks only, and the index is capped at 10
  (a heavier tail than estimated is conservative, a lighter one is not).

Declared deviation from the specification (before real data): the primary background is
`RadialBackground`: isotropic, a kernel density of log|v| below the 0.99 quantile and a Pareto tail above it (Hill
index, capped at 10). It keeps what Proposition D.1 asks for, a polynomial tail, f ~ |v|^-(3 + alpha). The isotropic
t stays as a declared ablation, with the gaussian.

## Design (fixed before the run)
- Networks: the Task 2M REAL-arm networks of `results/task2M/cfg_SEED.json` (`code/t2m_lib.py` unchanged): station
  subsets, catalogue, station corrections. Seed 2718: day 101231 (transients), 24 networks x 3 cells. Seed 3141: days
  100210 + 100211, 17 networks x 3 cells (ladder_4 networks 2909 and 2911 exhaust the catalogue pool, as registered).
- Slots: each network uses its noise block (about 50 min of consecutive 40 s slots). The catalogue slots are removed;
  the rest, in time order, alternate FIT (even position) and TEST (odd). No TEST slot touches the fit.
- Reading: the registered hard reading of the matched-filter output (MFn = MF / sd of the window).
- Background per station: fit on the distinct peak read vectors at the predicted samples of 400 random lattice
  hypotheses per FIT slot (time slides on FIT data), P and S pooled. R (primary), T and G (ablations).
- eps per station: 1 - r_s, with r_s the polarisation reliability (`ttmt.v2.m0.orientation_audit`) of the P vectors of
  the network's catalogue events, clipped to [0.05, 0.5].
- Model error: s_k = s0 + 0.3 |A_k m|, s0 = median FIT read norm / sqrt(3) of the station (same for R, T, G).
- pi, J = 2000 draws per network, seeded, independent of TEST data:
  - CAL: the catalogue tensors (least squares at the catalogue location), smoothed (lognormal amplitude sd 0.5,
    isotropic jitter 0.25).
  - LOW (stress): the same directions rescaled so that the median predicted block norm over 200 lattice nodes is
    2 x the median s0. Putting the pulses at the noise level is where a misspecified background shows most.
- Hypotheses: 40 random (node, origin) of the registered coarse lattice per TEST slot, seeded, independent of data.
  About 2 000 TEST slots and 80 000 hypotheses for seed 2718.
- **Gate G-M5c:** with background R, for pi CAL and pi LOW, the pooled TEST rate of S_mix >= log(1/p) is <= p, for
  p = 1e-2 and 1e-3, on seed 2718. Seed 3141 is the replicate and must pass too.
- Reported with it: 95 % upper bounds by a bootstrap of time clusters (network block, slot; the three cells share the
  noise), mean exp(S_mix), per-cell rates, the T and G ablations (prediction: G exceeds p on the transient day, and T
  exceeds p as in the smoke test), and the step-0 fit (PIT deciles and tail ratios on TEST reads).
- Not in this gate: the maximum over the lattice (the detector's selected hypothesis). That multiplicity belongs to
  M7 (time slides and e-BH), and is reported later.
- Held-out data: G-M5c chooses nothing. Every parameter above is fixed in this card. The Ridgecrest EVAL segment
  (spent, CR-23) is not used. The 2M noise days were used by G-M5a for ACAP, never for S_mix.

## Cost (Ockham)
CPU only. Step 0 (`measure`) is seconds per network. The full run is about 80 000 hypotheses x 2 pi with three
backgrounds sharing the model term: a few minutes on 4 processes. Small run first: networks 2800 and 2801, three cells.

## Commands
    python experiments/exploratory/M5-EVALUE/run_c.py measure 2718 NOISE_101231 2800,2801 > measure_2718_out.txt
    python experiments/exploratory/M5-EVALUE/run_c.py run 2718 NOISE_101231 PARTS_C_2718 [nets] [cells]
    python experiments/exploratory/M5-EVALUE/run_c.py run 3141 NOISE_100210,NOISE_100211 PARTS_C_3141
    python experiments/exploratory/M5-EVALUE/run_c.py analyze PARTS_C_2718,PARTS_C_3141 2000 > run_c_out.txt

## Results (commit 44db66c; `measure_2718_out.txt`, `run_c_small_out.txt`, `run_c_out.txt`)
Laptop, CPU, conda env ttmt: numpy 1.26.4, scipy 1.17.1. Parameters as declared. Wall time: measure 12 s, small run
46 s, full run 4 min 14 s (both seeds, 8 processes), analyze 4 s. Seed 3141 has 49 parts: ladder_4 2909 and 2911
exhaust the catalogue pool, as registered.

**Analysis bug, found in the run.** `RadialBackground` has no support below its kernel grid, so a TEST read with a
smaller norm gets log f = -inf, l_k = +inf, and the log-mean-exp computed inf - inf = NaN. A NaN compares False with
the threshold, so `analyze` counted these hypotheses as non-rejections (46 in 2718, 15 in 3141; the printed "PASS" of
3141 is this artefact). Fixed after the run in `evalue_mix` (S = +inf) and in `analyze` (NaN read as +inf), with a
test. The rates below count them as rejections; the method is unchanged.

| background, pi | 2718 (80 640 hyp.) p 1e-2 | p 1e-3 | 3141 (53 280 hyp.) p 1e-2 | p 1e-3 |
|:--|--:|--:|--:|--:|
| **R, CAL** (gate) | 9.08e-3 | **5.80e-3** | 3.12e-3 | **1.22e-3** |
| **R, LOW** (gate) | 9.37e-3 | **6.01e-3** | 3.04e-3 | **1.11e-3** |
| T, CAL (spec. t) | 1.38e-2 | 1.06e-2 | 1.45e-3 | 6.8e-4 |
| T, LOW | 1.42e-2 | 1.06e-2 | 1.46e-3 | 7.3e-4 |
| G, CAL (gaussian) | 1.9e-1 | 1.2e-1 | 4.0e-3 | 2.5e-3 |
| G, LOW | 1.7e-1 | 1.0e-1 | 4.1e-3 | 2.5e-3 |

- **G-M5c FAILS on both seeds** at p = 1e-3 (2718 by a factor of about 6, 3141 by 1.1-1.2). It holds at p = 1e-2.
- The gaussian background fails by two orders of magnitude on the transient day, as Proposition D.1(a) predicted.
- The t is worse than R on the transient day and better on the quiet days (its tail is heavier; its bulk is wrong).
- Signature: the rate barely falls from p = 1e-2 to 1e-3 (9.1e-3 to 5.8e-3). The exceedances are not a shifted bulk
  but a few hypotheses with huge S, i.e. single blocks far beyond the fitted tail.
- Step 0 agrees: R's tail index often sits at the cap 10 (quantiles 5/50/95: 1.79 / 7.21 / 10 in 2718), and the
  TEST/model tail ratio reaches 9.6 at 1e-2 and 61 at 1e-3 at ATPI, ATVO and MURB (ladder_6 2800). The worst cell is
  unilateral_c40 (8 stations, 16 blocks). eps sits at the 0.5 clip in most of 2718 (median r_s < 0.5).

## Diagnosis (hypotheses, not yet tested)
1. **Tail underfit (main suspect).** A Hill index from the top 1 % of about 50 minutes of FIT reads cannot see the
   rare transients of a TEST slot; a too-light tail makes log(1/f) of a transient read grow like (3 + alpha) log r.
2. **Support hole at 0.** Explains the NaN cases (S = +inf from one tiny read), not the bulk of the excess.
3. **Dependence between blocks.** Theorem D.2 multiplies over blocks; a transient seen by several stations at once
   breaks the product. unilateral_c40 being worst is compatible with it.

## Proposed second step (owner decides; not run)
- Measure first (seconds, 2718, development): for every exceeding hypothesis, the largest single-block l_k, whether
  that read is above r0 (tail) or below the grid (hole), and how many blocks exceed together (dependence).
- Then one variant fixed before running, chosen by that measurement (for example a heavier pooled tail per day with
  alpha capped at 2-3, and a power-law left tail), confirmed on a noise day that no M5 step has used.
- If the dependence term dominates, a per-block e-value cannot be calibrated by a marginal background; the decision
  then rests on M7 (time slides calibrate the joint null, D.3), and S_mix is kept only as a ranking-free evidence.
