# Experiment card: M6-SPLIT-AUDIT. Known-truth audit of the S-splitting constraint (rules 39 and 40)

Stage: S1, declared 23 November 2026 before any run on its seeds (requested by the supervision thread the same
day, under the owner's rules 39 and 40). Simulation only, in the cloud, with new seeds; no real data. The result
changes none of the four records (M6-SPLIT, M6-SPLIT2, M6-SPLIT3, M6-SPLIT2b); it only locates the cause.

Context. A post hoc audit (project files `m6-split/auditoria-2026-11-23/`, seeds 998 and 999, 20 scenes) found that
`free_one_event` mis-fits under splitting. With the fix, the true-splitting oracle gave median -107 and 0 of 20 above
23.8. M6-SPLIT2c (declared, reviewed) runs the fixed constraint. This card makes the audit formal, at full size.

## Parts, readings fixed before the run
Settings: M6-CAL scenes (8 stations, SNR 12, mismatched pulses unless stated), exact slow-S delay (`gen_fix.py`),
splitting drawn per station phi ~ U[0, pi), delta ~ U[0.02, 0.2] s.

1. **Recovery with known truth.** 100 single isolated events, each with a known splitting. Estimate (phi_i,
   delta_i) with `m6split` in station mode from that event alone (one atom read by `m6.adcg_loop`, then
   `joint_refit` with splitting). Reported for mismatched pulses (the run setting) and for matched pulses.
   - Pass: delta error <= 0.02 s and phi error <= 15 degrees (modulo 180) at >= 80 % of all stations, pooled
     over the 100 events.
   - Descriptive only: the same rates for stations with delta >= 0.05 s, and with 10 isolated events per station
     (the arm iii setting).
2. **Negative control with oracle.** 200 single-source control scenes with splitting in the data. The splitting
   is held at its true value in both hypotheses (`mode="fixed"`, `fixed_split`), with the one-event fix on
   (`free_atom_start=True`). Count the scenes split at thr 23.8.
   - <= 12 of 200: the model does not invent pairs when it knows the truth, so the failures of the line lie in
     the estimation or in the implementation.
   - More than 12 of 200: the model invents pairs even with the truth, and the negative belongs to the model.
   - Also reported, descriptive: the same with the fix off (`free_atom_start=False`), to measure the defect at
     full size.
3. **Positive control with oracle.** 40 near pairs per dtau in 0.1, 0.25, 0.5 and 1.0 s, with splitting in the data
   and the true splitting held. Recall (LOC+T-V, both sources, thr 23.8) must be within 0.05 of arm (i) on the
   same pair scenes without splitting (no splitting in data or model, `m6.adcg_loop`). Exact counts: pass if
   100 * (hits_i - hits_oracle) < 5 * 160.

## Seeds and cost
Recovery 20261161 (its 10 isolated events 20263161; B.SEED_ISO set to 20262161, unused); controls 20261162; pairs
20261163; splitting 20261164. Smoke test on seed 999 only (recovery 6 s per event, control 9 s, pair 9 s). Each job under 2 h: part 1; part 2 with the
fix; part 2 without the fix; part 3 (oracle and arm i). Run once, after the M6-SPLIT2c jobs on the same 4 cores.
