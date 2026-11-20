# Experiment card: WHITEN-REVIEW. Which past results the whitening mismatch touches, and by how much

Stage: S0, declared 20 November 2026 before any run. Exploratory and post hoc; no registered number changes.
Plan: review of `supervision/revision-blanqueo-2026-11-17.md` (supervisor's design of 17 November 2026). The fix
itself is ADR-0012 (other thread); this card only audits and measures.

## The defect [measured, FE-FIX step 2, PR #68]
Where the data of each station are divided by its own noise level and the subspace Q is built from the unweighted
operator G, a noise-free event keeps only f 0.705 (Hi-net) / 0.771 (Ridgecrest) of its energy in the subspace; with Q
from diag(1/sd)·G, 0.999. The fix changes **Q only**.

## Part 1. Code audit, result by result (done before any run)
Rule used: a result is touched by the fix only if its code projects per-station-normalised data on a subspace built
from the unweighted G. A result that never builds or uses Q is bit-identical under the fix.

| Result | Per-station division | Uses Q from unweighted G | Verdict | Evidence |
|:--|:--|:--|:--|:--|
| R0 (Ridgecrest, `code/r0`) | yes, `r0_run.Net` line 50 (`station_scale`) | yes, `X.basis(radiation_operator)` lines 62-65 | **affected** | code |
| R1 (Hi-net, `code/r1`, `src/ttmt/r1.py`) | yes, `r1.py:60-61` | yes, `node_operators` | **affected** | code |
| C18 b, c, d, g, i, j (built on `r0_run.Net`) | yes | yes | **affected** | code |
| C18a (`c18a_timing.py`, timing readings) | yes, line 37 | no: per-station peak search (`station_corrections.peak_delays`), an argmax per station, invariant to a per-station factor | **not affected** | code |
| G-M1a (`M1A-ONSET`, uses `r1.Net("cal")`) | yes, through `net.d16` | no: uses only `STA`, `calev`, `d16`, `raw`, `cat`, `names` (peak and onset reads per station) | **not affected** | code |
| M5 third step `--local` (`M5-EVALUE/run_c.py:104-106`) | yes | yes (`E.G`) | **mismatch present but inert**: that step ran only on the container's i.i.d. t4 noise, equal scale at every station, so the per-station levels are equal up to sampling error and diag(1/sd) is a multiple of the identity; it was a level test on a null and the negative rests on the product path, which is coherent | code + check below |
| Task 2M, M5 first two steps | no (`kernel.matched`: one sd per window) | — | not affected (already corrected 17 Nov) | code |
| M6-JOINT, M6-CAL, M6-THEORY, M6-BIAS (all on `M6-JOINT.noise_mf`) | no: `r.normal(size=(ns, 3, NT))`, unit sd at every station, one global sd | — | not affected | code, `M6-JOINT/run.py:53, 78` |
| M6C-SUBSAMPLE | no `station_scale` anywhere | — | not affected | grep |
| Comparators | do not use this code | — | not affected | — |

Correction to propagate (rule 19): the thesis marks the whole of §6.5 (C18) as "[whitening review]". The timing
statements of §6.5 (C18a: 0.45 against 0.37-0.39 within 0.25 s; 0.427/0.427; 0.417/0.387) are not affected; the
score, AUC and recall statements are.

Check of the M5 `--local` claim (cloud, `check_m5_local.py`, seed 20261126): per-station robust level of i.i.d. t4
matched-filter windows, max/min ratio across stations; reported in `out_m5_local.txt`. Reading fixed now: if the
median max/min ratio is <= 1.10 the mismatch is inert (FE-FIX had 5.7-7.6 on real noise); otherwise the verdict
changes to "affected" and the step is listed for a re-run.

## Part 2. How much the fix changes the R0 failure on real data (C18b with weighted Q)

**Question.** C18b measured, on the Ridgecrest calibration segment only (origins in [0, 900) s; EVAL not read), at the
true hypocentres, how well the first-iteration score separates the 16 target events (M >= 2.0) from 222 event-free
windows (AUC). Registered normalisation (REG): AUC 0.39-0.49; without cap (NOCAP): 0.69-0.74, on N4, N6, N8, U8.
Does building Q from diag(1/scale)·G change these numbers?

**Data.** Ridgecrest 2019 hour, `benchflow-ai/skillsbench` commit `9a1f4dd5`, the three files checked against
`data/checksums/ridgecrest_2019_R0.sha256` (match, 20 November 2026). Calibration segment only, as C18b. This
segment is development data (data card); no reserved window and no evaluation segment is read.

**Arms.** REG and NOCAP exactly as C18b; REG-W and NOCAP-W identical except that every node basis is
`basis(diag(1/scale) G)`, with `scale` the per-station `station_scale` that `r0_run.Net` already divides by (P and S
rows of station i divided by scale_i).

**Code checks, before reading anything.**
1. REG and NOCAP reproduce the frozen `results/explorations/c18/c18b_*.json` AUC and q99 fractions exactly on the four
   subnets. If not, stop and diagnose (rule 26).
2. With scale set to 1 at every station, REG-W equals REG bit for bit.
3. Noise-free: an event fabricated as diag(1/scale)·G·m at a node keeps f >= 0.99 in the weighted basis and < 0.95 in
   the unweighted one whenever max/min scale >= 3 (the FE-FIX check, on this geometry).

**Endpoint.** Per subnet and per pair (REG-W vs REG, NOCAP-W vs NOCAP): difference in AUC, with a 95 % interval from
10 000 bootstrap resamples of events and null windows (seed 20261126). Secondary: fraction above the null q99, and
the median s_true / q90.

**Reading, fixed now (descriptive; three zones, rule 36).** On NOCAP-W − NOCAP, pooled over the four subnets
(mean of the four differences; the subnets share stations, so the pooled interval resamples events and null windows
jointly across subnets):
- >= +0.05 with the interval above 0: the defect explains part of the R0 failure; the thesis says the R0 numbers are
  conservative for the model, and by how much.
- between -0.05 and +0.05: the defect does not explain the R0 failure; the thesis says so with the number.
- <= -0.05: unexpected; diagnose before any statement.
The same reading on REG-W − REG is reported, flagged descriptive: the cap saturates (Proposition 2), so a small
change there is expected either way.

**Limits.** 16 events per subnet: intervals will be wide. One segment, one network. This does not redo R0 or R1
(their evaluation segments are spent and stay as recorded, marked "affected" in the thesis); it measures whether the
defect could have changed their conclusion. The weighted basis is implemented here in 3 lines, independently of the
ADR-0012 code; when ADR-0012 merges, `run.py` is re-run with its option and must give identical numbers.

## Part 3. What the laptop is not needed for
Everything above runs in the cloud (data from GitHub, about 10 minutes of CPU). R1 is not re-run: its development
segment has 30 targets and the same question is answered by FE-FIX step 3 on the tune partition.
