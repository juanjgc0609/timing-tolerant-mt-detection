# Experiment card: M6-JOINT. Joint refit and resolvability test (simulation proxy of G-M6a)

Stage: S1 (done, 8 November 2026). Plan: T4.19 (M6 before M2, owner's decision of 8 November). Exploratory,
simulation only; no registered number changes.

## Question
Does the joint refit of all accepted atoms with the "one against two" resolvability test (`ttmt.v2.m6.adcg_loop`)
recover overlapping pairs that the greedy v2 loop (one atom at a time, 5 km exclusion) loses, without splitting
single events? G-M6a asks for +0.10 of pair recall at 0.1-1 s.

## Mechanism tested
- Joint refit: tensors of all atoms in closed form, block times refined with the others fixed, a block may be
  re-assigned to any onset maximum within +-0.5 s (scored by the full model's residual), 3 alternations.
- Resolvability (replaces the 5 km exclusion): a new atom within 5 km and 1.5 s of an accepted one is kept only if
  the two-event fit beats the best one-event fit with free time (+-8 samples) and free 3-vector per block
  (Theorem A.3(a)) by thr = z_e^2 = 10.8 (Corollary A.1) + max(0, 95th percentile of the gain on a noise-only
  control).

## Data
Simulation as M6C-SUBSAMPLE (8-station ring, 1.5 km lattice, homogeneous medium, synthetic picker) with raw-trace
Gaussian noise; SNR 12 (median strong-block peak / noise sd). Two sources at dtau 0.1, 0.25, 0.5, 1 s, epicentres
6-14 km apart ("far") or 1.5-5 km ("near", inside the exclusion). 20 scenes per cell, 40 + 40 single-source controls,
seed 20261108. 10 min on 4 cloud cores (with the control column), after smoke runs of 1-3 scenes.

## Result (`run_out.txt`, with the review's control column)
| stratum | dtau | greedy (5 km) | greedy, no exclusion | M6, no test | M6 |
|:--|--:|--:|--:|--:|--:|
| near | 0.10 | 0.45 | 0.95 | 0.95 | **0.95** |
| near | 0.25 | 0.35 | 1.00 | 1.00 | **1.00** |
| near | 0.50 | 0.50 | 0.70 | 0.75 | **0.75** |
| near | 1.00 | 0.85 | 1.00 | 1.00 | **1.00** |
| far | 0.10-1.00 | 0.95-1.00 | same | same | same |

Pair recall (both sources hit, LOC+T). False splits of 40 single events (fresh control): greedy without exclusion
**40 of 40**, M6 at thr 10.8 **0 of 40** (the greedy loop with the exclusion cannot split by construction).
- **The recall gain on near pairs comes from removing the 5 km exclusion, not from the joint refit**: the greedy
  loop without the exclusion recovers the same pairs (0.70-1.00 against 0.75-1.00 for M6).
- **What M6 adds is that removing the exclusion becomes admissible**: without it the greedy loop splits every single
  event into two atoms; the resolvability test removes all of those splits and keeps every true pair here. The
  joint refit adds +0.05 at 0.5 s and the deflation quality with pairs (a noise-free check of 4 pairs per
  separation, not stored: median peak residual 0.03-0.05 against 0.09-0.13 for the greedy v2 loop; to be measured
  properly as G-M6c with M6).
- So the contribution to claim is the resolvability test with the joint refit as a pair, against "exclusion" and
  "no exclusion": same recall as no exclusion, the false-split rate of the exclusion.
- The resolvability test never rejected a true pair here.
- Diagnosis found while building it: the free one-event fit first used a bounded local search, whose cost oscillates
  with the response lobes; it ended on side lobes and was worse than the tensor fit, so the test accepted every
  split (gains around +1700). A grid over whole samples before the local search fixed it (median gain of a spurious
  split -102). Test `test_free_one_event_fit_is_at_least_as_good_as_the_tensor_fit` guards it.

## Limits (declared)
Homogeneous medium, Ricker pulses, no coda, a synthetic picker without false maxima, SNR fixed at 12, kmax 2. The
real gate G-M6a is on injected pairs with real noise (semi-synthetic sets of the generator v2) and is not measured
here. The near 0.5 s cell stays at 0.75 for both M6 variants (not yet diagnosed).

    python experiments/exploratory/M6-JOINT/run.py [NREP] [NCTRL] > experiments/exploratory/M6-JOINT/run_out.txt
