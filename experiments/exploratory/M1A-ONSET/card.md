# Experiment card: M1A-ONSET. Gate G-M1a, reading at the onset

Stage: S0 (prepared 8 November 2026; not run: needs the Hi-net raw data and the picks environment)
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
