# Experiment card: P6-SAP. The single beat rule and the sample size of phase 6

Stage: S1 (done, 23 October 2026). Statistics only; no data read.
Owner decisions: ADR-0008 accepted by the author on 24 October 2026; the atlas thread copies it into ATLAS-0 before anchoring.

## Question
Which decision rule can show that the model beats every tuned comparator, holding its error rate, with a number of
targets the data can supply; and how many blocks and targets does each test need?

## Mechanism tested
The rules of ADR-0008 as coded in `src/ttmt/stats.py`: block sign-flip tests, intersection-union over comparators,
fixed sequence H1 -> H2, the LOC+T-V hit rule and its duplicator control, the quasi-Poisson budget check.

## Result (`sample_size_out.txt`, tests in `tests/test_stats.py`)
- Per-cell rules have no power: with a true +0.05 in every target cell and 100 targets per cell, a majority of 10
  cells is won with probability 0.000 (400 per cell: 0.31); one pooled test on the same targets has power 0.86.
- At one-sided 0.025, 80 % power needs about 790 independent targets for +0.05 and 205 for +0.10 (rho 0.7);
  non-inferiority needs 2,135 for a margin of 0.03 and 770 for 0.05.
- Clustered design (3 targets per block, icc 0.1, rho 0.7, 30 % of targets in T, two comparators at equal recall):
  +0.05 needs about 1,100 blocks (990 targets in T, 3,300 in total); +0.10 needs 525 blocks for the chain H1 -> H2
  with H2 at 0.025 (470 in T, 1,570 in total), 550 with H2 at 0.02 or 0.015.
- Size must be checked at the least favourable configuration of the intersection-union test (one comparator at the
  boundary; independent reviews of PR #11 and #13). There H1 is exact (size at most 0.019). H2 uses a clustered score
  test with its variance under the restricted null and Satterthwaite degrees of freedom: the margin sign-flip test
  had size up to 0.033 and declared 30 targets without discordance non-inferior (p = 0.0001, now 0.109), and a
  design-effect variance had size 0.06 to 0.11 with large, very unequal, correlated blocks. The design check
  (`sample_size.py calibrate`, range from development data) keeps H2 at 0.025 over the default range from 30 to 300
  blocks and lowers it to 0.02 at 100 and 0.015 at 300 blocks when the range includes 30 targets per block with
  dispersion 0.3 and intra-block correlation 0.4 (part 3).
- LOC+T counts two detections of one event as a resolved pair; nearest-event eligibility does not fix it (the
  duplicator resolves about half the pairs); the dominance condition of LOC+T-V does (at most 1 of 200).
- The false-alarm budget check needs 7 to 150 hours of record for dispersions 1 to 20.

## Consequence for the plan
K1, K2, K9 and CR-08 resolved by ADR-0008. The data request of phase 6 is sized per regime: R3 (cheap
injection) at +0.05, R2 (real) at +0.10, and R2 tests the target regime only with network decimation (CR-04).
