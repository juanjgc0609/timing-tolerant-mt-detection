# Experiment card: M6-THEORY. Theory of the M6 resolvability test, checked on the M6-JOINT simulation

Stage: S1 (done, 9 November 2026). Plan: T5.2 (review of Corollary A.1 and Theorem A.3(a) as used by M6) and T5.6
(theory of the "one against two" test). Exploratory, simulation only; no registered number changes. Full text (Spanish):
`teoria-fase5/teoria-m6-prueba.md` in the project files (Part E of the Phase 5 theory).

## Questions
1. Are Corollary A.1 and Theorem A.3(a) correct, and does `ttmt.v2.m6` use them as they are stated?
2. What is the null distribution of the gain `rss_one - rss_t`, and is thr = z_e^2 + max(0, q95) valid?
3. What separation does the test need to accept a pair at SNR 12 (median peak of the strong blocks, the M6-JOINT
   definition; not the network SNR sqrt(E) of Appendix AB), against the information limit?
4. Do S-wave splitting, a free source time function and opposite polarities (CR-09, CR-10, CR-25) change 2-3?
5. Why does the near 0.5 s cell of M6-JOINT stay at 0.75?

## Results (`run_out.txt`, `confounders_out.txt`)
- **A.1 and A.3(a) are correct; their use as a threshold is not.** z_e is a necessary separation in the whitened norm;
  the code's gain is in the matched-filter (MF) norm. Exact factors for the LEGACY kernel (Proposition E.1, closed form
  for the Ricker, checked to 4 digits): A2 = 3.291, kappa1 = 2.962 (first order), kappa2 = 1.481 (second order).
- **Null (Proposition E.2).** The one-event class has 6 n_s amplitude dofs (48) and the pair 12, so at fixed times the
  mean gain is -lambda_r (6 n_s - 12) = -136 (lambda_r = 3.775); Monte Carlo -135.6. The code's gain on the M6-JOINT
  controls: median -101.8, q95 -51.3, max +4.0 (search adds +34). So q95 < 0 always and thr = z_e^2 (CR-30): the
  false-split control is conservative in Gaussian noise because of a dimension count, not because of A.1.
- **Power (Proposition E.3).** Acceptance at 50 % needs an MF margin of about thr - median null gain = 113, i.e. a
  whitened delta_TT of 6.2 (different mechanisms) or 8.7 (same mechanism) against z_e = 3.29: 1.9x and 1.6x the
  separation of the information limit. Measured (40 scenes per cell): 50 % acceptance at dtau ~0.028 s (different
  mechanisms) and 0.09 s (same; MF margin there 113.4), i.e. 2.2x and 1.7x the information limit of this geometry
  (0.013 s, 0.054 s). Acceptance >= 0.95 from 0.05 s and 0.12 s. The code's noise-free gain equals the exact MF margin
  (Lemma A.4 with the response's autocorrelation) to 4 digits in 10 of 14 cells and within 3 % in the rest.
- **S-wave splitting breaks the test.** With splitting in the null class, the first-order S term vanishes at any angle
  (Proposition E.4(a); slopes 3.8-4.0). The current null has no splitting: single events with split S give false
  pairs in **7/20 (0.05 s) and 14/20 (0.10 s)**, 0/20 at 0 and 0.2 s.
- Free STF (rank-1 bound, E.4(b)) keeps the first order with different mechanisms and makes same-mechanism pairs at
  dtau <= D(M) unidentifiable; opposite polarity (E.4(c)) is resolved at first order with constant 0.257 omega^2.
- **The near 0.5 s cell is a scoring artefact.** `hits()` matches greedily; at dtau = 0.5 s (the tolerance) a source
  can be left without its atom. Optimal matching on the same atoms: 1.00 (greedy 0.75). M6A-REALNOISE uses the same
  `hits()`.

## Consequences (proposed, not activated)
- Calibrate thr as a declared quantile on single-event controls with real noise, mismatched pulses and split S; drop
  z_e^2 from the threshold or declare it as an ad hoc margin.
- ~~Add split S (free fast axis, two scalar fits per S block) to the one-event class before reading G-M6a on real data.~~
  Superseded by the addendum, item 2: split S enters calibrated per station, not free.
- Joint P-S fit with the closed form of A.3(a) for networks with S - P < 1 s; compare against all nearby atoms.
- ~~Score pair recall with optimal matching (M6-JOINT, M6A-REALNOISE).~~ Superseded by the addendum, item 3: score
  with the registered rule, `ttmt.v2.m6.hits_loctv`.

## Addendum (14 November 2026, supervision review; `split_margin.py`, `split_margin_out.txt`)
1. **Two definitions of "SNR 12".** In M6-JOINT, M6-CAL and this card, SNR 12 is the median peak of the strong blocks
   in noise sd. In Appendix AB (A.7) SNR is the network sqrt(E). For the same scene the whitened delta^2 at 0.05 s is
   136 here and 18 at network SNR 12, a factor of about 7.5 in energy (2.7 in SNR). The separations of this card
   (0.028 s, 0.05 s, 0.09 s, 0.12 s) must not be set beside the A.7 table without naming the SNR.
2. **What a split-S null costs.** Co-located pairs of different random double couples, ring of 8, straight rays, far
   field, whitened norm, network SNR 12 per event, median of delta^2 over 40 pairs (information limit 10.8):

   | dtau | P | S, one pulse (current null) | S, free split | split, delay <= 0.1 s | split calibrated per station |
   |--:|--:|--:|--:|--:|--:|
   | 0.05 s | 0.05 | 18.2 | 0.75 | 0.77 | 18.2 |
   | 0.10 s | 0.29 | 48.2 | 5.5 | 7.9 | 48.2 |
   | 0.15 s | 0.71 | 59.0 | 13.4 | 24.8 | 59.0 |
   | 0.25 s | 0.93 | 56.7 | 18.9 | 52.6 | 56.7 |
   | 0.50 s | 1.44 | 83.4 | 29.3 | 83.2 | 83.4 |

   A free split in the null (the consequence proposed above) lowers the margin 3 to 25 times over the whole T regime,
   because the first-order S term vanishes (E.4(a)) and P adds almost nothing. With the fast axis and the delay of
   each station known (estimated on isolated events, like the station terms), nothing is lost, and this is exact: the
   split operator maps the two events to two orthogonal channels with the same shift, so the captured energy equals
   that of the unsplit pulse. It holds while the two events share the path (co-located pairs, the T regime). Design
   consequence: split S enters M6 calibrated per station (or through an EGF dictionary), not as free parameters.
   Limits: a simplified geometry, not the code's test (whose noise-free gain equals this margin, E.3); the stability
   of the splitting parameters per station on real data is a hypothesis, to be measured on the calibration events.
3. **Scoring.** The registered rule is one-to-one matching in score order with the dominance condition
   (`ttmt.stats.label_voronoi`, ADR-0007 and ADR-0008), which `ttmt.v2.m6.hits_loctv` applies. An optimal matching
   would be a new rule outside the ADR. The fix for the 0.5 s cell is to use `hits_loctv` everywhere (M6-CAL already
   does; M6A-REALNOISE is changed in PR #46), not to add an optimal matching.

## Limits (declared)
Same simulation as M6-JOINT (ring of 8, homogeneous medium, Ricker pulses, Gaussian noise, SNR 12). The tail bound of
E.2 holds at fixed times only; the search inflation is measured, not bounded. The review of A.1/A.3(a) is by an AI
session: CR-19 (human review) stays open.

    python experiments/exploratory/M6-THEORY/run.py 40 40 > experiments/exploratory/M6-THEORY/run_out.txt
    python experiments/exploratory/M6-THEORY/confounders.py 20 > experiments/exploratory/M6-THEORY/confounders_out.txt
