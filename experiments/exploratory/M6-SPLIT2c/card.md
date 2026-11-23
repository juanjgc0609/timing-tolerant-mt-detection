# Experiment card: M6-SPLIT2c. The S-splitting constraint again, after the implementation audit of 23 November 2026

Stage: S1, declared 23 November 2026 before any run on its seeds. The owner asked on 23 November 2026, 01:34 UTC:
"prueba todo antes de darlo como negativo y verifica errores de implementacion". The M6-SPLIT line is therefore not
yet recorded as negative. M6-SPLIT, M6-SPLIT2, M6-SPLIT3 and M6-SPLIT2b stay as run (rule 28); this is a further
attempt, reported as such.

## What the audit found (project files `m6-split/auditoria-2026-11-23/`, seeds 998 and 999, outside every declared seed)
- The decisive check holds the TRUE splitting fixed in both hypotheses. Single-source controls, matched pulses, 20
  scenes. Without splitting the median of the largest spurious gain is -109 (0 of 20 above 23.8); with splitting it
  is +147 (13 of 20). So the false pairs came from the implementation, not from estimating the splitting.
- **Error 1 (the cause).** `m6split.free_one_event` starts its per-station time search only from the unsplit free fit.
  When S is split, that start is off and the search stays there, so the one-event alternative (a LARGER class than
  the atom) fits worse than the atom alone. Noise-free, true splitting: 4.4e-2 relative RSS, against 1.1e-3 for the
  atom's tensor fit. Fix: `atom_start` also starts from the atom's own block times and keeps the better
  (`free_atom_start=True` in `adcg_loop`; default off, so earlier runs reproduce). With the fix: 3.4e-4. Oracle check
  with the fix: median -107, 0 of 20. Test `tests/test_m6split.py::test_free_one_event_with_atom_start_nests_the_split_atom`
  fails without the fix and passes with it.
- **Error 2 (minor).** The split generators of M6-SPLIT* delay the slow S component by linear interpolation at a
  sub-sample delay, which removes about 1 % of its energy (0.8-1.2 % at 0.03-0.17 s). Alone it does not cause the bias
  (oracle median +168 with only this fixed). Fixed here by an exact Fourier delay (`gen_fix.py`, error 6e-12).
- Both fixes together: oracle median -101 (matched) and -62 (mismatched), 0 of 20 above 23.8 in both.
- The same error 1 affected every split arm of M6-SPLIT, M6-SPLIT2, M6-SPLIT3 and M6-SPLIT2b (all use
  `free_one_event`). Their numbers stay as recorded; their validity failures are not evidence against the idea.

## Arms, seeds, criteria (thresholds unchanged, rule 14)
- Arms: (i) no splitting in data or model (`m6.adcg_loop`, reference); (w) window-only constraint (the director's
  design, as M6-SPLIT arm iii); (iii) shared with NISO = 10 isolated events (the owner's choice "Restricción común a
  la estación", as M6-SPLIT2b arm iii, converged isolated fits). Mismatched pulses, SNR 12, M6-CAL settings (kmax 4,
  near_km 10), both fixes on in (w) and (iii).
- New seeds: controls 20261151, pairs 20261152, splitting 20261153, isolated events 20261154, pre-check 20261155.
  200 controls, 40 near pairs per dtau in 0.1, 0.25, 0.5, 1.0 s.
- Criteria for (w) and (iii), each against (i), exactly as in M6-SPLIT2:
  - validity: at most 12 of 200 control scenes split at thr 23.8;
  - recall on exact counts: a loss below 0.05 passes; 0.05 to 0.15 inclusive is useful and the owner decides; more
    than 0.15 costs too much;
  - AUC of the pair gain >= 0.90.
- Pre-check (before the main run): 20 single-source scenes on seed 20261155 for (a) the oracle (true splitting held,
  split data) and (b) no splitting with the model of (iii). Each median of the largest spurious gain must be in
  [-200, -50]; otherwise stop and diagnose. Nothing is tuned to reach it.
- Predictions (written before the run): (w) and (iii) pass validity; the recall loss of (iii) is smaller than that of
  (w).
- If (iii) fails validity, the pre-calibrated variant (M6-SPLIT3 design, with the fixes) is declared next as the
  owner asked ("si hace falta").

## Procedure and cost
Smoke test on seed 999 only (one control scene: (i) 2 s, (w) 7 s, (iii) 30 s per loop). Jobs, each under 2 h, as in
the header of `run.py`: pre-check; (i, w); (iii) controls; (iii) pairs; combine. Run once.

## Project rules 39-43 (added by the owner 23 November 2026; this note added before any main-run output was read)
- 39 (audit with known truth before a negative): done. The oracle with the true splitting is the known truth; the
  audit and this card's pre-check both use it.
- 40 (positive and negative controls in every new test): negative controls are the pre-check without splitting and
  the oracle with splitting (must not split: medians -85.9 and -96.9, 1 and 0 of 20 above 23.8), and the 200
  single-source controls of each arm. Positive controls are the 160 true near pairs of each arm, which must be
  resolved, and arm (i), the reference that must reproduce M6-CAL-level recall. Also the unit test, which fails
  without the fix and passes with it.
- 41 (after each merge, fetch main and check the resolution): applies when PR #75 merges.
- 42 (no untracked result files in the laptop checkout): not applicable (cloud only).
- 43 (say at close which rules applied and how): in the result section.
- Note from review: the claim that the earlier failures are not evidence against the idea rests on a post hoc
  diagnosis on seeds 998 and 999. Those failures stay recorded as run.
