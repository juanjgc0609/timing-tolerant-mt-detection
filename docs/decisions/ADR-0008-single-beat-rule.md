# ADR-0008: One rule for "the model beats the comparators"

Date: 23 October 2026. Status: accepted by the author on 24 October 2026 (proposed by the statistical design
thread); the atlas thread writes it into the ATLAS-0 card before the card is anchored. It resolves conflicts K1, K2 and
K9 of the master plan and critique items CR-02 and CR-08; it extends ADR-0007 for close pairs.

## Context

Three different rules were written for the same claim (K1):
- the ATLAS-0 card on `main`: 48 cells, a cell is won when the paired 95 % interval is above zero, weighted majority
  of cells;
- the hand-over of the device session: 240 cells, weighted majority of target cells plus non-inferiority at -0.05;
- the model specification v2, Section 9: +0.05 in at least two thirds of the target cells, non-inferiority -0.03,
  target regime 0.1 to 1 s (the hand-over used < 0.1 to 1 s).

None of them can be met by a method that is in fact better (K2). A cell rule needs a significant test in each cell.
With 100 targets per cell and a true gain of +0.05 in every target cell, a cell is won with probability 0.11, and the
probability of winning a majority of 10 such cells is 0.000; one pooled test on the same 1,000 targets has power 0.87
(`experiments/exploratory/P6-SAP/sample_size_out.txt`, part 2). With 240 cells the rule would need 96,000 to 288,000
injected targets per method, and it has no correction across cells. Point-estimate rules (specification v2) have the
opposite defect: no inference at all, so they pass by noise (CR-08).

The hit rule has a hole exactly in the target regime (CR-02): with a tolerance of 0.5 s and pairs 0.1 to 1 s apart,
two detections of the same event can be matched one to one to the two events of a pair and count as two hits.

Nothing says where the phase-4 gate G4 is measured (K9).

## Decision

1. **Regimes.** Each target event is assigned, by quantities computed before any method runs (catalogue or injection
   design, station availability over [origin - 10 s, origin + 60 s] with no gap), to:
   - the **target regime T**: arrival separation to the nearest other event at the median station in [0.1, 1) s, any
     station count, signal-to-noise ratio or distance;
   - the **complement C**: every other event, including separations below 0.1 s, which the theory puts below the
     resolvable limit (Appendix AB, Cor. A.7), and above 1 s.
   The cells of the atlas grid (separation in 5 bins, signal-to-noise ratio, station subset, distance, K_mec/E as a
   covariate) are strata for description only. No decision counts cells.
2. **Hit rule LOC+T-V.** ADR-0007 matching (one to one, decreasing score order, epicentre within 5 km and origin within
   0.5 s) with one more condition: a detection is eligible for an event only if its normalised cost
   c = max(distance / 5 km, |origin error| / 0.5 s) to that event is at most half its cost to every other catalogued
   event (`ttmt.stats.label_voronoi`). For two events at the same epicentre, a detection counts for one of them only
   within a third of their separation of it. For an isolated event the rule is LOC+T. The rule stays prefix-stable, so
   the exact FROC is unchanged. Plain nearest-event eligibility is not enough: a copy of a detection moved by half the
   separation falls on the partner's side half of the time.
3. **Negative control, mandatory.** Every analysis runs the duplicator of the best comparator (each of its
   detections emitted twice, at origin +-half the separation of the cell, `ttmt.stats.duplicator`). If the duplicator
   gains recall in T over the comparator it duplicates, with the one-sided test of item 5 at 0.025, the metric is
   declared broken and no beat claim is made from that analysis. Duplicates of isolated events are reported per method.
4. **Endpoint.** Recall of the targets of a regime under LOC+T-V at 30 false detections per hour (exact FROC, operating
   point chosen per method on the evaluated data, as in R0 and R1), pooled over the regime's cells with weights equal to
   the share of each cell in real sequences divided by its share among the evaluated targets. The real shares are the
   prevalence measured on the densest reviewed catalogues (CR-04) and are fixed in the card before it is anchored; if
   they are not available then, the weights are 1 (each evaluated target counts once) and that is stated in the card.
5. **Test.** For each comparator k, the paired recall difference is tested by the block sign-flip test
   (`ttmt.stats.signflip_pvalue`): the replication unit is the independent block (an injected network realisation in
   semi-synthetic data, a time block of real data, at least 30 per sequence), d_j is the weighted sum over the
   block's targets of (hit by the model - hit by k), and the signs of whole blocks are flipped. It is the exact
   randomisation form of the clustered McNemar test (Durkalski et al. 2003), valid for any dependence inside a block.
   One-sided level 0.025. The paired 95 % interval of the stratified block bootstrap (unit = block, strata = sequence
   or network, operating points recomputed in each replicate) is reported with every test; it is flagged indicative
   with fewer than 30 blocks or an atom above 5 % at a bound.
6. **Beating the comparators** means beating each of them: the claim against the comparator set is an
   intersection-union test whose p-value is the largest pairwise p-value (Berger 1982); it needs no multiplicity
   correction and is not weakened by adding comparators. The set is the core set of ATLAS-0, each comparator tuned
   under item 9, plus the best comparator re-run with the model's picker when that picker is fine-tuned (CR-03).
7. **The rule.** The model beats the comparators when, in this fixed order at 0.025 each:
   - **H1**: superiority in T (H0: recall difference <= 0 against every comparator); and then
   - **H2**: non-inferiority in C with margin 0.05 (H0: difference <= -0.05 against every comparator, the sign-flip test
     on d_j + 0.05 n_j).
   H2 is tested only if H1 is rejected, so the family-wise error is 0.025 without adjustment. Secondary hypotheses,
   tested only if H1 and H2 are rejected, form one Holm family at 0.025: superiority in [0.1, 0.5) s, superiority in
   [0.5, 1) s, and non-inferiority in C at the tighter margin 0.03. Cells are described with their estimates and
   intervals, never counted. The claim is worded "superior in the regime T, which holds Y % of the events of dense
   sequences, and not inferior by more than 0.05 elsewhere", with Y the prevalence of item 4 and its interval.
8. **Effect sizes for planning.** Semi-synthetic analyses are sized for a gain of 0.05 in T and real ones for 0.10,
   with `experiments/exploratory/P6-SAP/sample_size.py` (it simulates this rule). The margin 0.05 replaces -0.03 because
   at 0.025 one-sided a margin of 0.03 needs about 2,100 independent targets in C per comparator and 0.05 about 770
   (part 1 of the output); the 0.03 reading stays as a secondary hypothesis.
9. **Equal tuning budget (CR-07).** Every method, the model included, gets the same number of configurations (60 unless
   the card says otherwise), drawn by the same random search with the same seed list, scored on the same development
   partition by the same objective (development recall in T and C at 30/h under LOC+T-V, weighted as in item 4).
   Every configuration tried is appended to a tuning ledger (`tuning_ledger.csv` next to the card: method, trial,
   parameters, objective, time) whose SHA-256 enters the freeze. Parameters the model fixes from theory count as
   fixed, not as free trials. A sensitivity table over the ledger is published.
10. **Where G4 is measured (K9).** The phase-4 gate uses this rule on development data: primary on semi-synthetic
    injection (known truth; generator v2 or recorded events injected into real noise with injection and deflation
    libraries disjoint by event and by time block); secondary on real development data scored against the truth of
    task T3.V, where it is a consistency check: the point estimate in C must not fall below -0.05. G4 is not
    confirmatory. The phase-6 tests (registrations R2, real data with network decimation, and R3, semi-synthetic on
    unseen networks) use the same rule on data no method has read (docs/protocol/PHASE6_SAP.md).

## Consequences

- The ATLAS-0 card replaces its "counting rule" section with a citation of this record; the specification v2 (Section
  9.5) and the roadmap cite it instead of repeating a rule. The atlas thread keeps the grid edges, which become strata.
- The number of targets is set per regime, not per cell: with a third of the targets in T, intra-block correlation
  0.1 and method correlation 0.7, the semi-synthetic gate needs the number of blocks in part 4 of the output.
- Earlier registered results keep their rules. The atlas description by cell (phase 3) is unaffected: it is
  descriptive, and its hypothesis on comparator losses is tested by its own card.
- A real-data claim in T is possible only where the truth resolves pairs at 0.1 to 1 s: a dense or nodal network
  decimated to subsets of 4 to 12 stations (CR-04). Otherwise T is tested on semi-synthetic data only, and the claim on
  real data is limited to C.
