# Experiment card: M1A-ONSET. Gate G-M1a, reading at the onset

Stage: S1 (8 November 2026, owner's laptop: first step FAILED on CAL, 0.87 < 0.90; second step, declared before
running, PASSED on the held-out EVAL events, 0.92 registered reader, 0.97 prior variant, read blocks only)
Plan: T4.13 (step 2 of the build order of `especificacion-modelo-v2.md` §8). Exploratory and post hoc; no registered
number changes.

## Question
Does reading every block at the onset (the largest local maximum of the PhaseNet probability within the search
window, no threshold, no minimum peak distance) put >= 90 % of the P readings of the R1 calibration events within
0.25 s of the reference time, with leave-one-event-out station terms (gate G-M1a)?

## Mechanism tested
The registered rule reads the matched-filter energy peak. That peak lags the onset by 0.24-0.27 s and can jump to a
side lobe or to another arrival (Appendix AB, B.1). On these events it gives 55 % of P and 76 % of S within 0.25 s,
and 36 % and 20 % more than 0.5 s off (`results/taskR1/posthoc/timing_cal_out.txt`). A picker is trained to mark
the onset, so its local maxima should not jump.

## Hypothesis and refutation
G-M1a passes if the line "(a) onset, same blocks", P, shows "within 0.25 s" >= 0.90. Below that the onset reading
does not meet the competence condition of the thesis (§7.1) and step 2 is diagnosed before step 3 (the likely causes,
each with its check: PhaseNet misses on small P, seen in (b) against (a); a terms problem, seen in the robust sd;
the 2 s window catching another event, seen by shrinking the window).

## Data
Hi-net 2012 swarm, CAL segment, subnetwork H7 (7 stations), the 17 truth events with M >= 1.0 that the registered
terms used. JMA2001 times from the catalogue location. Same windows and reference as `code/r1/posthoc_timing_cal.py`.

## Endpoint and analysis
Fractions within 0.25 s and 0.5 s and the robust sd of the leave-one-event-out residuals, per phase, for (a) the
registered energy-peak reading and the onset reading on the same gated blocks and (b) the onset reading on every
block with a maximum. Also the onset-to-peak lag per phase, the input `lag` of `ttmt.v2.loop.v2_loop`.

## Environment and cost
The picks environment (`environments/picks.yml`, TensorFlow PhaseNet, as the registered R1 picks) plus the package.
One PhaseNet pass over the CAL hour of 7 stations on CPU (minutes); the probabilities are cached in `q_cal.npz`
(not committed) so later steps reuse them. No GPU.

    python experiments/exploratory/M1A-ONSET/run.py > experiments/exploratory/M1A-ONSET/run_out.txt

## Result of the first step (`run_out.txt`, 8 November 2026)
Laptop, conda env `ttmt-picks` (Python 3.11.17, NumPy 1.26.4, TensorFlow 2.15.1), CPU, 14 s. PhaseNet clone at the
pinned commit 62005c6 (cloned during the run; `real_picker` clones HEAD, so step 2 checks the commit first).
- **G-M1a fails narrowly**: onset reading, same gated blocks, P 0.87 within 0.25 s (n 38), threshold 0.90; S 0.87.
  The registered energy-peak reading: P 0.55, S 0.76. Robust sd falls from 0.222 to 0.111 s (P).
- 9 of the 47 gated P blocks have no onset maximum >= 0.02 in the window; on every block with a maximum (n 86) P is
  also 0.87, so the misses are not small P lost by PhaseNet. The robust sd is low and 8 % of P are beyond 0.5 s:
  tails, which points at the window catching another arrival rather than at the terms.
- Onset-to-peak lag: P +0.326 s, S +0.290 s (R1 had measured 0.24-0.27 s with the anchored terms).

## Second step (declared before running; `run2.py`)
1. List every CAL failure (residual, q, number of maxima, nearest maximum) to confirm or refute the "other arrival in
   the window" diagnosis.
2. One variant fixed in advance, nothing tuned on these events: "prior" reads the maximum of
   log q - (s - mu)^2 / (2 sigma^2) (the per-block local term of M4), mu = prediction + term, sigma = robust sd of
   the phase, both leave-one-event-out on CAL.
3. Confirmation on EVAL (840-1960 s), never used to build anything: truth events M >= 1.0, terms and sigma from all
   CAL events. G-M1a is read on EVAL for the registered reader ("largest") and for "prior"; the CAL figures of the
   variant are post hoc. The energy-peak reading is reported on the same EVAL blocks.

## Result of the second step (`run2_out.txt`, 8 November 2026)
Same laptop and environment, PhaseNet clone checked at 62005c6, CPU, about 16 s including the PhaseNet pass over EVAL.
- **G-M1a passes on EVAL** (held out): the reader as registered gives P 0.92 within 0.25 s (gated, n 61), S 0.99;
  the declared "prior" variant P 0.97, S 0.99. The registered energy-peak reading on the gated EVAL blocks: P 0.52,
  S 0.78. On CAL the variant raises P from 0.87 to 0.95 (post hoc there).
- Diagnosis confirmed: 5 of the 16 CAL failures belong to one event (532.6 s) with readings 1.1-1.6 s early at four
  stations, another arrival inside the +-2 s window; in several failures the nearest maximum is within 0.07 s of the
  reference. One event (712.2 s) shows a uniform -0.3 s shift on five high-q single-maximum readings: a reference
  (catalogue origin or location) problem, not a reading error.
- **Caveat, coverage**: 42 of the 103 gated EVAL P blocks (and 9 of 47 on CAL) have no onset maximum >= 0.02 in
  the window; the gate is measured on the read blocks only (61 of 103). Counting the unread blocks as misses, P
  within 0.25 s is 56/103 = 0.54 for the onset reader against 0.52 for the energy peak (which reads every gated
  block): the gain of M1 is in precision, not in coverage. Those blocks stay "unread" in M4/M6 (prediction + lag), and
  PhaseNet's misses on low-SNR P are the C-regime loss of the stage diagnosis (PR #38), the target of the fine-tuning
  step (G-M1b).
- Decision for the v2 loop: the "prior" reading (largest log q - (s - mu)^2 / (2 sigma^2)) becomes the default
  onset reader, with sigma P 0.089 s, S 0.058 s from CAL: implemented as `OnsetReader(sigma=)` and
  `ttmt.v2.loop.default_reader` (test `test_prior_reader_prefers_the_maximum_near_the_prediction`).
