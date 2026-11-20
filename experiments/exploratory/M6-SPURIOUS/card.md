# Experiment card: M6-SPURIOUS. Why about 5 % of control scenes get a positive spurious pair gain (calendar M-2)

Stage: S0, declared 20 November 2026 before any run. Simulation only (cloud). Exploratory; nothing registered changes.

## Question
M6-CAL step 5 measured that about 5 % of control scenes (target + overlapping neighbour + isolated source, no true
pair) have a spurious near atom with a positive pair gain, up to about 400 (`M6-CAL/card.md`, step 5). Route B uses
the gain as a ranking score, so this tail is where M6 makes false pairs. Its cause is not identified. Which one is it?

## Hypotheses, each with the measurement that separates it (fixed now)
- **H1, pulse mismatch** (hint of the M6 thread, 20 Nov): the data pulse is stretched (x U[0.7, 1.4]) and
  phase-rotated (U[-45, 45] deg) per arrival, the kernel is not; a one-event fit leaves a shaped residual that a second
  atom absorbs. Test: the SAME scenes (same noise, sources, mechanisms; only the pulse draw differs) with exact pulses.
- **H2, neighbour**: the spurious atom takes energy of the overlapping neighbour (6-14 km, 0.5-1.5 s later). Test:
  one-source scenes (kind C1) against scenes with neighbours, both mismatched; and, within the scenes with neighbours,
  whether the spurious atom's nearest true source is the target, the neighbour or the isolated source.
- **H3, mismatch x neighbour**: neither alone, both together.

## Design
Loop as M6-CAL step 5 (`adcg_loop`, kmax 4, joint one-event fit, near_km 10, thr = -inf so every candidate is
kept and its gain recorded), SNR 12 (SNR_blk, `docs/GLOSSARY.md`). New seed 20261127 (not 2026111x, not
20261121-23 of M6-SPLIT). 600 scenes with neighbours and 600 one-source scenes, each run twice: mismatched and exact
pulses. Per scene: the largest spurious gain (as M6-CAL, level per window) and, for that atom, the distance and time
to each true source.

## Endpoint and reading (written before running)
p(arm) = fraction of scenes whose largest spurious gain is > 0. Paired difference mismatched vs exact on the same
scenes, exact McNemar test (`ttmt.stats.exact_mcnemar`).
- **H1 named** if, with neighbours, p_exact <= 0.2 x p_mis and the McNemar p < 0.01.
- **H2 named** if, mismatched, p(C1) <= 0.2 x p(neighbours) and, among the positive-gain scenes with neighbours,
  at least 60 % of the spurious atoms are nearest to the neighbour (s2) rather than to the target or the isolated
  source.
- **H3** if p_exact <= 0.2 x p_mis with neighbours AND p_exact(C1) and p_mis(C1) are both <= 0.2 x p_mis(neighbours).
- Otherwise: **no clear cause**; the tail is declared a limit in the M6 card (calendar M-2, "si sale mal"), with the
  table.
Descriptive: distribution of positive gains per arm; the gain against the atom-to-source distance.

## Limits
Homogeneous medium, Ricker family, Gaussian noise, ring geometry of M6-JOINT. A named cause is a cause in this
simulator; G-M6a on real noise is where it has to hold.
