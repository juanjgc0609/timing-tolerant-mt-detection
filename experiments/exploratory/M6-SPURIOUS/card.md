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

## Amendment 1 (20 November 2026, after review of 079b3d2, before any output of the full run was read)
The full run was started at 079b3d2 (smoke run of 4 scenes per kind read; nothing chosen from it). The review found
the reading ambiguous; it is replaced, before reading, by this one, applied by `analyse.py` to `rows.json`:
- **Nearest source** of the spurious atom: the true source with the smallest normalised cost max(d / 5 km, |dt| /
  1.5 s); ties go to the lower index (target, neighbour, isolated).
- **Order of reading** (the first that holds is named):
  1. **H3** (mismatch x neighbour): with neighbours p_exact <= 0.2 x p_mis (McNemar p < 0.01), and p_mis(C1) <=
     0.2 x p_mis(C3) (Fisher exact, one-sided, p < 0.01).
  2. **H1** (pulse mismatch): with neighbours p_exact <= 0.2 x p_mis (McNemar p < 0.01), and p_mis(C1) is not small
     (the H3 condition on C1 fails).
  3. **H2** (neighbour): mismatched, p(C1) <= 0.2 x p(C3) with Fisher p < 0.01, exact pulses do not remove the tail,
     and at least 60 % of the positive-gain spurious atoms with neighbours are nearest to the neighbour.
  4. Otherwise no clear cause.
- Expected size: about 5 % of 600 scenes with neighbours, so about 30 positive scenes in the mismatched arm; the
  60 % share and the ratios are read on that small n and reported with Clopper-Pearson intervals.

## Declarations before reading (20 November 2026, after review of 789873a; `rows.json` not yet written)
1. The full run executes the `run.py` committed in 789873a: sha256 `c7fdc3fb0bf8f33d83c541f330d594df0cbea5faca719905085b7b2c8e83af56`,
   `git diff 789873a` empty. The same shell command added the empty-tail guard and then started the run, and the file
   has not been edited since (mtime 22:24:00 UTC, process start 22:23:59-22:24:00 UTC at the one-second granularity
   of `ps`). The smoke run of 4 scenes per kind used the version without the guard; its scenes are the first 4 of
   each kind of the same seed and are inside the full run.
2. The line "reading (card): ..." that `run.py` prints in `run_out.txt` applies the superseded, non-exclusive rule
   of 079b3d2. **The valid reading is `analyse_out.txt`** (amendment 1, `analyse.py`); the `run_out.txt` line is
   kept for the record and is not the result.

## Result (20 November 2026, cloud, 1126 s on 4 cores; `analyse_out.txt` is the valid reading)

| arm | scenes with a positive spurious gain (of 600) | 95 % CI |
|:--|--:|:--|
| one source (C1), mismatched pulses | 28 (0.047) | [0.031, 0.067] |
| one source (C1), exact pulses | 0 | [0, 0.006] |
| with neighbours (C3), mismatched | 26 (0.043) | [0.028, 0.063] |
| with neighbours (C3), exact | 11 (0.018) | [0.009, 0.033] |

C3 paired: 25 scenes positive only when mismatched, 10 only when exact, McNemar p 8.3e-3. Fisher C1 against C3
(mismatched): p 0.66. Nearest source of the positive C3 atoms: target 8, neighbour 11, isolated 7 (neighbour share
0.42 [0.23, 0.63]). Largest gains: C3 mismatched 386, C3 exact 593.

**Reading of amendment 1: no clear cause.** Exact pulses remove 58 % of the tail with neighbours, short of the 80 %
the rule asks (11 > 0.2 x 26); C1 is not smaller than C3; the neighbour share is below 60 %. Per calendar M-2 the
tail is declared a limit of M6 in simulation.

Descriptive, after seeing the table (not a named cause): the tail has two parts. (a) With one source it is
entirely pulse mismatch: 28 of 600 positive when mismatched, 0 of 600 exact (paired, every positive scene is a
mismatched-only one). (b) With neighbours a part remains with exact pulses (11 of 600, up to 593), which one-source
scenes never show; its nearest source is spread over target, neighbour and isolated source. A test of (b) needs its
own card (for example, the same scenes with the neighbour moved beyond 1.5 s), and of (a) a pulse family in the
kernel or in the one-event fit. Neither is run here.
