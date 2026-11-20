# Experiment card: WHITEN-REVIEW. Which past results the whitening mismatch touches, and by how much

Stage: S1, declared 20 November 2026 before any run (commit 1c9f0ef); run the same day in the cloud, results below. Exploratory and post hoc; no registered number changes.
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
| G-M5b, Ridgecrest (`M5-ACAP/run_b.py`, built on `r0_run.Net`, `E.Q`) | yes | yes | **affected**; added at the results stage (the thesis note of 17 Nov lists it; my first table missed it) | code |
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

## Results (20 November 2026, cloud, Python 3.12 with `requirements-frozen.txt`)

Outputs: `out_m5_local.txt`, `out_analyse.txt`; raw per-window scores `out/wr_*.json` (22 KB each; regenerate with
the commands at the top of `run.py`, about 3 minutes per subnet on one core).

**Code checks.** Check 1: REG and NOCAP reproduce the frozen C18b AUC to 1e-12 and the 222 null windows on all four
subnets. Check 2: with scale 1 the weighted basis is the registered Q bit for bit. Check 3 (N8 geometry, station
scale max/min 4.48): noise-free f 1.0000 (min) in the weighted basis, 0.818 (median, min 0.669) in the unweighted one.
All pass.

**M5 `--local` check: does NOT pass as written.** Median max/min per-station level on i.i.d. t4 windows: 1.085 (4
stations, passes), 1.106 (6), 1.120 (8), against the threshold 1.10 fixed before running. The M5 networks have at
least 8 stations (`M5-EVALUE/card.md`, station rule of `t2m_stations.py`), so 1.120 applies. **By the rule: affected,
listed for a re-run.** Post hoc, produced after seeing the ratios: the threshold was too tight (sampling error alone
reaches it), and at those levels the noise-free f in the unweighted basis is 0.9992 (median; min 0.9950, N8 geometry,
500 draws), a loss below 0.005 against 0.18-0.33 on real noise. On that post-hoc argument I recommend not re-running
(the step was a level test on a null, and the M5 negative rests on the coherent product path); the owner decides.

**C18b with Q_w (primary reading, NOCAP-W − NOCAP):** +0.013, +0.014, +0.013, +0.006 on N4, N6, N8, U8; pooled
**+0.011, 95 % [+0.001, +0.025]** (16 events, 222 null windows, 10 000 resamples, seed 20261126). Zone: **for the uncapped score the
defect does not explain the R0 failure** on real data at the true hypocentres. Under NOCAP the fraction above the null q99 is
unchanged except U8 (0.062 to 0.125, one event); under REG, N6 and U8 drop from 0.062 to 0.

**Descriptive (REG-W − REG, the registered capped score):** +0.052, +0.101, +0.102, +0.062; pooled +0.079, 95 %
[+0.027, +0.136]. The registered statistic goes from below chance (0.39-0.49) to about chance (0.46-0.60); no event
gets above the null q99 on any subnet. So the defect pushed the capped score below chance, and the cap (Proposition 2)
still leaves it without separating power after the fix.

**What this changes.**
1. Reading of the rule: for the uncapped score (primary) the defect does not explain the R0 failure (+0.011); for
   the registered capped score it explains part of it (+0.079, the "explains part" zone of the same rule, reported as
   descriptive as declared), from below chance to about chance; after the fix neither score separates real events.
   R0, R1, C18 b-j and G-M5b Ridgecrest stay as recorded and stay marked "[whitening review]" with this note.
2. C18a (timing) and G-M1a are not affected: their markers can be lifted for those statements.
3. When ADR-0012 merges, `run.py` is re-run with its option in place of `weighted_Q` and must reproduce
   `out_analyse.txt` exactly (pending, this thread).
