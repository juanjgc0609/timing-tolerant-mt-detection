# ADR-0013: Minimum acceptable gain in T (more than 2 points), next to the 5-point target of ADR-0011

Date: 20 November 2026. Status: accepted by the owner (see "Approval"). Written before any P1 number exists.
Number taken from the ledger of the master plan (§11; ADR-0012 is reserved for the FE-FIX front-end change).

## Context

ADR-0011 sets a target for the overlap line: Δ_T ≥ 0.05, where Δ_T is the point estimate of the recall in T of our
better arm minus that of the best frozen comparator, on the P1 tuning partition, LOC+T-V at 30 false detections per
hour (ADR-0011 items 2 and 3). Under rule 36 (ADR-0011 item 4) this is an ambition target, read in three zones. Below
the target, the owner decides with the numbers in view.

ADR-0011 does not say what gain is enough for the project to stand. The owner met the thesis director (Uram Aníbal
Sosa Aguirre) and reported his view in the project chat on 20 November 2026, 12:33 (Bogotá):

> "Y me dijo que 5 puntos era muy optimista, sin embargo, pues si se llegaba a lograr, bien, pero si no se alcanzan a
> llegar 5 puntos, no significa que esté mal el proyecto, sino que él me dijo que con más de 2 era suficiente, pero
> digamos que la meta es llegar a 5, pero si no se llega a 5 no hay problema, pero eso sí, por lo menos más de 2."

## Decision

1. **Two numbers, two roles.** The target stays at Δ_T ≥ 0.05 (ADR-0011, unchanged). A **minimum acceptable gain**
   is added: Δ_T > 0.02 ("más de 2 puntos"). It is measured exactly as ADR-0011 item 2 defines Δ_T: same partition,
   same arms, same comparators, same hit rule, same false-detection budget, and the point estimate as in ADR-0011
   item 3.
2. **The three zones of ADR-0011 item 4 become four.** Each is reported with Δ_T, n_T, the paired counts and the
   stratified block-bootstrap 95 per cent interval:
   - **Target met:** Δ_T ≥ 0.05. The line continues (ADR-0011 item 5).
   - **Minimum met, target not met:** 0.02 < Δ_T < 0.05. The report says "meets the director's minimum (ADR-0013),
     does not meet the target (ADR-0011)". The line continues, with no owner decision needed. The claim carried to the
     confirmatory test is stated at the size measured.
   - **Below the minimum, real effect:** 0 < Δ_T ≤ 0.02 and the interval excludes zero. The owner decides with the
     numbers in view, as in ADR-0011 item 4, zone 2.
   - **Clear failure:** Δ_T ≤ 0.02 and either the interval includes zero or Δ_T ≤ 0. The owner decides whether to
     reformulate, as in ADR-0011 item 4, zone 3.
3. **What it is.** The minimum is the director's floor for a defensible gain. It is not a validity threshold and does
   not relax one. The false-detection budget, non-inferiority in C (ADR-0008), held-out data read once, the criterion
   written before the run, and the fair comparison (amendment 6) are unchanged. Reaching the minimum on the tuning
   partition is not a claim. The confirmatory claim still needs ADR-0008 (H1 in T, then H2 in C) on data not used for
   tuning.
4. **What it does not change.** ADR-0011 items 1 to 3 and 5 to 7, ADR-0008, the ATLAS-0 card and its amendments,
   every ledger and every frozen configuration. In ADR-0011 item 4, the zone "meets the target" is unchanged, and its
   zone "real but smaller effect" is split in two at 0.02.
5. **Power, for the record (an estimate, not measured).** It uses the same assumptions as ADR-0011 item 3: about 30 per
   cent discordant pairs, a standard error near 0.037, block inflation excluded. With those assumptions, a method
   whose true gain is 0.05 lands above 0.02 about 79 per cent of the time, and one at 0.02 about half of the time.
   Separating 0.02 from 0 on the tuning partition is not possible at n_T = 223. That is why the confirmatory test of
   ADR-0008 is where the claim is decided. Under the same assumptions the interval excludes zero only for Δ_T above
   about 0.073 (1.96 × 0.037), so the zone "below the minimum, real effect" is practically unreachable in P1: any
   Δ_T ≤ 0.02 will in practice read as "clear failure". The owner decides in both zones, so nothing changes in
   practice, but the zone is kept for the confirmatory data, where n is larger.
6. **Change of governance, declared.** In ADR-0011 item 4, a gain below the target always went to an owner decision.
   Here a gain between 0.02 and 0.05 continues without one, following the director's words ("si no se llega a 5 no
   hay problema, pero eso sí, por lo menos más de 2"). The owner saw this on the approval card.

## Approval

The owner relayed the director's words (quoted above). The thresholds and zones of item 2 were written by the planning
thread from those words and from ADR-0011. Reading "more than 2" on the point estimate is an interpretation of the
planning thread, since the director did not specify it. On 20 November 2026, 12:41 (Bogotá), the owner chose "Aprobar
como está" on a decision card in the project thread. The card stated that between 2 and 5 the line continues without
his decision and that below 2 he decides. The alternative it offered was to require the whole interval above 2.

## Consequences

- The thesis can describe a gain between 2 and 5 points as meeting the director's minimum without overstating it.
- The figure "5 points" in any text is the target, not the condition for the result to exist.
- If the owner wants "more than 2" read on the lower bound of the interval rather than on the point estimate, item 1
  is amended before P1 runs.
