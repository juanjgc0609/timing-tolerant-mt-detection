"""Failure atlas of the published detection pipelines (phase 3, card experiments/exploratory/ATLAS-0/card.md).

Modules:
  schema      truth tables, picks and detections in the common frame (local km, seconds from record start)
  cells       the regime axes of the card and the cell of every target
  froc        LOC+T-V labels, operating points, the descriptive tables by cell and regime, the pooled losses of the
              atlas and the rule of ADR-0008 (H1 in T, then H2 in C, with the duplicator control), on ttmt.stats
  instrument  pick-level and event-level diagnostics (lost, merged, wrong-phase picks; merged and split events;
              picker- versus associator-limited misses; LOC hits that fail LOC+T)
Comparator wrappers and their tuning live with the comparators (thread "Comparadores y ajuste"); they hand the atlas
Detections and Picks in the schema below. Only numpy is needed here.
"""
