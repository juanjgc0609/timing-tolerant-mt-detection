# ADR-0011: Criterion for abandoning the overlap line (pilot P1 after block M and a picker fine-tune)

Date: 15 November 2026. Status: accepted by the owner (see "Approval"). Number assigned from the ledger of the master
plan (§11; ADR-0010 is the M5 pivot of PR #43, not yet written as a file).

## Context

The independent verification of 14 November 2026 (supervision thread, document
`verificacion-y-diseno-2026-11-14.md`, Parts 3 and 7, in the project files) found that after four weeks of
per-component gates no end-to-end number of the model existed on the partition where the comparators were tuned,
and that every gain of M6 had been measured against our own loop. It proposed:

1. pilot P1 before any further polishing of components: our two arms (S, the best frozen associator plus the M6
   split test on each of its detections; D, M0 + M1 with prior + M6) against the frozen comparators on the tuning
   partition (`b60d99b3`, 50 records, n_T = 223, n_C = 305), under the registered endpoint (LOC+T-V, 30 false
   detections per hour);
2. a reading table for P1, fixed before it runs, whose "T does not win in any arm" row sends the work to block M
   (resolution), to a picker fine-tuned with pairs (T4.16, the same picker for every method, CR-03) or to M3;
3. a stop rule for the whole line, so that the pivots of that row cannot go on without end: "if, after block M and
   one picker fine-tune, the best arm does not reach best + 0.05 in T, the thesis is reformulated to atlas +
   identifiability theory + negatives, which is already defensible".

Item 3 changes the goal of the project (G1-v2) when it fires, so it goes by decision record (owner rule 20) and is
written before any P1 number exists.

## Decision

1. **When it is read.** Once, when both of these are done: (a) block M (M.1 split S calibrated per station, M.2
   calibration on isolated events, M.3 the pair gain scoring the second event in the FROC, with G-M6a read); (b) one
   fine-tune of the picker with pairs (T4.16), given to every method of the comparison (CR-03). P1 is then re-run
   with the frozen model of that point. The reading does not wait for any other component (M2, M3, M7).
2. **What is measured.** On the tuning partition of ATLAS-0 amendment 1, item 4 (the P1 partition), at 30 false
   detections per hour under LOC+T-V (`ttmt.stats.label_voronoi`, the registered hit rule; no thread-local scoring,
   rule 27): Δ_T = recall in T of the better of our two arms minus recall in T of the best comparator. "Best
   comparator" is the highest recall in T among PyOcto, GaMMA (amendment 4) and HARPA (amendment 5, once anchored),
   each in the configuration fixed by ATLAS-0 amendment 6, item 3 (with or without station terms, whichever has the
   higher tuning objective). If amendment 5 is not anchored when the rule is read, HARPA is left out and
   the best comparator is the better of PyOcto and GaMMA, stated with the number. Our arms use the configuration of amendment 6, item 1 (station terms where an
   independent set exists).
3. **The threshold.** The line continues if the point estimate Δ_T ≥ 0.05. It is reported with n_T, the paired
   counts and the stratified block-bootstrap 95 per cent interval (ADR-0008 item 5), but the interval does not
   decide. Why the point estimate (an estimate, not measured: it assumes about 30 per cent discordant pairs, so a
   standard error near 0.037 before block inflation): a method whose true gain is 0.10 continues about 91 per cent
   of the time under this rule, and one at 0.05 about half; a rule "lower bound above 0" would let the 0.05 method
   continue only about 1 time in 4 and the 0.10 method about 4 in 5. This is a development stop rule, not a claim,
   so it errs toward not killing a real gain. Choosing the better of two arms on the
   same partition favours us and is declared with the number.
4. **If it fires** (Δ_T < 0.05): the result is reported as "below the 0.05 threshold of ADR-0011", with Δ_T, n_T,
   the paired counts and the 95 per cent interval, and with what the gain is made of (which pairs, which strata).
   Falling short does not abandon the line by itself. It triggers a written decision by the owner, who has the
   effect size and the interval in view. The owner chooses between: (a) reformulate G1-v2 and the thesis to "failure
   atlas + identifiability theory + documented negatives", which the thesis thread then writes; or (b) carry the
   result forward as it is, if it is still good and new as a contribution (for example, a smaller gain confined to
   pairs no picks-based method separates). Option (b) must state the claim at the size measured, with no wording
   that implies the 0.05 was met. Either way, the threshold, the number and the choice are written down together in
   one record, and the threshold is not changed after the fact (rules 14 and 15). No further model variant is tried
   to lift Δ_T before that decision. The atlas (ATLAS-0), the theory and the negative results stay as they are.
5. **If it does not fire** (Δ_T ≥ 0.05): nothing is claimed from it. The confirmatory claim still needs ADR-0008
   (H1 in T, then H2 in C with margin 0.05) on data not used for tuning, after block F.
6. **Relation to other stop rules.** The decision rule of the model v2 specification §12 (three or more of P-PN1,
   P-GM1, P-RE1, P-QM2, P-MA1 refuted ⇒ reformulate before M6) stays. Whichever fires first applies.
7. **What does not change.** ADR-0008, the ATLAS-0 card and its amendments, every ledger and frozen configuration.
   This record adds a stop rule for the development line only.

## Approval

The owner approved on 15 November 2026, 9:17 (Bogotá), in the project chat, everything that was waiting on him,
which the plan v2.12 listed as this abandonment criterion and the ATLAS-0 amendment for station terms (B2):
"Así que verifica y apruebo todo lo que haya que hacer. Ya." The operational details of items 1 to 4 (when it is
read, which comparators, point estimate) were written by the planning thread from the supervision proposal and are
put to the owner for confirmation on the one real choice, item 3. On 15 November 2026, 9:21 (Bogotá), the owner
chose "Estimación puntual ≥ 0.05" on a decision card in the project thread, over the alternative "lower bound of the
interval > 0". Item 3 stands as written. At 9:24 (Bogotá) the same day he added, in the project chat, that a threshold must not
discard a result that misses it but is still good and new: "imagina que uno por ser tan exigente y elegir ciertos
valores rechace todo incluyendo resultados buenos que no llegan a ese valor, pero que si uno reduce la exigencia
igual siguen siendo igual de buenos y novedosos en la contribución. Así que cuidado con eso." Item 4 was rewritten
after that message: missing the threshold triggers his decision with the effect size in view, not an automatic
abandonment. The threshold itself and item 3 are unchanged.

## Consequences

- P1 and its re-run use one partition and one script; the second run is a re-reading of a development set already
  seen by the first, which this record accepts for a stop rule and does not accept for a claim.
- The station-term amendment (B2) is ATLAS-0 amendment 6, drafted in PR #53
  (`experiments/exploratory/ST-RETUNE/amendment-6-draft.md`); its approval is recorded there, not here.
- No deadline is set. If block M or T4.16 stalls, the plan's risk register (§8.1) is where a time limit would be
  added, by a further record.
