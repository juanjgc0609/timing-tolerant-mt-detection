# Experiment card: M6-SPLIT2b. Second attempt of M6-SPLIT2, with the isolated fits converged

Stage: S0, declared 22 November 2026 before any run (decision of the supervision thread, 22 November 2026, 02:49
UTC, after M6-SPLIT2). Exploratory, simulation only; G-M6a is not run.

**Record.** M6-SPLIT2 stays run and FAILED by its rule (arm iii 126 of 200 split controls). Its defect note is a post
hoc diagnosis. This is a SECOND ATTEMPT, reported as such, not a replacement. Ledger row (rule 28; to be added to
`docs/DATA_USE_LEDGER.md`, section 2, when that file reaches main, PR of the thread that keeps it):
"M6-SPLIT line | attempts: M6-SPLIT (window-only constraint, fails validity 24/200), M6-SPLIT2 (shared with
isolated events, fails 126/200, implementation defect), M6-SPLIT3 (estimated once, fixed; run by M6-SPLIT2's
condition), M6-SPLIT2b (second attempt of M6-SPLIT2 with the fix) | result: see the cards | `M6-SPLIT*/card.md`".

## The fix (`ttmt.v2.m6split.adcg_loop(..., converge_iso_fits=True)`, `converge_iso`)
- The isolated events are converged before the loop: refitted (tensors and times, splitting fixed) until their RSS
  changes by less than 1e-6 relative (at most 30 rounds), alternating with the splitting estimate until the RSS
  changes by less than 1e-6 (at most 10 alternations).
- Both hypotheses use the same rule: after the two-event fit and after the one-event fit, the isolated events of that
  hypothesis are refitted to convergence by the same `converge_iso` before their RSS enters the hypothesis's RSS.
  The number of rounds is therefore set by the same convergence rule in both, not by a fixed count that one
  hypothesis reaches less converged.
- Nothing else changes (model, grid, thresholds, scenes, pulses).

## Pre-check (declared; run before the main run, `PRECHECK=1`)
20 single-source control scenes WITHOUT splitting in the data and in the isolated windows (seed 20261145), model of
arm iii: the median of the largest spurious gain per scene must be negative and of the order of the other arms of
M6-SPLIT and M6-SPLIT2 (about -100): declared range [-200, -50]. Outside it, the run stops and the cause is
diagnosed; nothing is tuned to reach the range.

## Arms, seeds, criteria
- Arms: iii (shared with NISO = 10 isolated events, fix on) and v-w (matched pulses, window-only constraint), as
  asked; plus i and v-i (no splitting, mismatched and matched pulses) on the same seeds, because the recall zones
  are read against arm i of the same scenes, and v-w is read against v-i.
- New seeds (not M6-SPLIT 20261121-23, M6-SPLIT2 20261124-27, M6-SPLIT3 20261131-34): controls 20261141, pairs
  20261142, splitting 20261143, isolated events 20261144, pre-check 20261145.
- Criteria, unchanged from M6-SPLIT2, read on arm iii against arm i: validity at most 12 of 200 control scenes split
  at thr 23.8; recall zones on exact counts (below 0.05 passes; 0.05 to 0.15 inclusive useful, the owner decides;
  more than 0.15 costs too much); AUC >= 0.90.
- Predictions (from M6-SPLIT2's card, unchanged): fewer split controls than M6-SPLIT's 24 of 200; the remaining
  recall loss below 0.1-0.25 s; v-w against v-i tells whether pulse mismatch is the larger cost.

## Procedure and cost
Two jobs, each under the 2 h job limit: `ARMS=i,v-i,v-w` then `ARMS=iii`, appended to `run_out.txt`.

    PRECHECK=1 python experiments/exploratory/M6-SPLIT2b/run.py > experiments/exploratory/M6-SPLIT2b/precheck_out.txt
    ARMS=i,v-i,v-w python experiments/exploratory/M6-SPLIT2b/run.py > experiments/exploratory/M6-SPLIT2b/run_out.txt
    ARMS=iii python experiments/exploratory/M6-SPLIT2b/run.py >> experiments/exploratory/M6-SPLIT2b/run_out.txt

## Result (second attempt, run once, 22 November 2026; `precheck_out.txt` and `run_out.txt` unchanged)
Pre-check: median of the largest spurious gain -92.8 (in [-200, -50]), run allowed (100 s). Jobs: i, v-i, v-w 2466 s;
iii 6287 s (both under 2 h). Non-converged isolated fits: 0 in i, v-i, v-w; 8 in iii.

| arm | control scenes split (of 200) | near-pair recall pooled | 0.10 s | AUC | median largest spurious gain |
|:--|--:|--:|--:|--:|--:|
| i (no splitting, mismatched) | 6 | 0.931 (149/160) | 0.825 | 0.952 | -103.7 |
| v-i (no splitting, matched) | 0 | 0.912 (146/160) | 0.800 | 0.938 | -113.2 |
| v-w (window-only constraint, matched) | 19 | 0.819 (131/160) | 0.700 | 0.881 | -83.8 |
| **iii (shared with isolated events, fix on)** | **114** | 0.850 (136/160) | 0.675 | 0.904 | **+86.8** |

- Arm iii against arm i: validity 114 of 200 against the limit 12, **FAILS**. Recall 13 of 160 lost, gap 0.081
  (exact: 1300 > 800, <= 2400): useful zone, owner would decide, but validity fails first. AUC 0.904 passes.
- Second attempt of M6-SPLIT2: **FAILS** on validity. Prediction "fewer split controls than M6-SPLIT's 24" not
  supported.
- v-w against v-i: 19 of 200 split (fails validity), 15 of 160 lost (gap 0.094), AUC 0.881. With matched pulses the
  window-only constraint still fails validity, so pulse mismatch is not the only cost of the splitting.
- The fix did not remove the positive bias: the pre-check (no splitting in the data) gives -92.8, but with
  splitting in the data the median is +86.8 (M6-SPLIT2: +115). So the bias is tied to estimating the splitting
  with the isolated events, not only to unconverged isolated fits. Limit of the declared pre-check: without
  splitting in the data it could not see this. Cause not diagnosed (hypothesis, not tested: the shared splitting
  search over pair window plus isolated windows reaches a better optimum in the two-event hypothesis than in the
  one-event one, because the grid search starts from a different place).
- G-M6a does not run with this model.
