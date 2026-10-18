# Experiment card: P6-power. Number of targets for a paired recall test

Stage: S1 (done, 18 October 2026)
Owner decisions: none needed; it fixes a number the data request of phases 3 and 6 depends on.

## Question
How many catalogued targets does a registered test need to detect a recall gain with 80 per cent power, under the
decision rule of Task R1?

## Mechanism tested
Only the statistics: the exact one-sided McNemar test at alpha 0.05 on paired hits of two methods.

## Hypothesis and refutation
Not a hypothesis test. Expected: R1 (54 targets) had little power for gains of 0.05 to 0.10.

## Data
Simulation only (`power_mcnemar.py`, seed 20261018, 4000 replicates per cell). Comparator recall 0.40, as in R1.
Hits drawn from a Gaussian copula with correlation rho of the two methods (0.5, 0.7, 0.9); targets independent.

## Endpoint and analysis
Power against the recall gain delta and the number of targets N; false-positive rate at delta 0.

## Result (`power_mcnemar_out.txt`)
- The test holds its level: false-positive rate 0.020 to 0.048.
- With **54 targets, as in R1**, power is 0.11 to 0.13 for a gain of 0.05 and 0.27 to 0.43 for a gain of 0.10. R1
  could only have detected a gain of about 0.15 or more.
- For 80 per cent power, the minimum N on the simulated grid is:

  | Gain | rho 0.5 | rho 0.7 | rho 0.9 |
  |:--|:--|:--|:--|
  | +0.03 | 2400 | 2400 | 1200 |
  | +0.05 | 1200 | 800 | 400 |
  | +0.10 | 300 | 200 | 150 |
  | +0.15 | 100 | 100 | 100 |

- These are lower bounds: targets clustered in time (aftershock sequences, shared noise hours) need N multiplied by
  the design effect 1 + (m - 1) ICC.

## Consequence for the plan
A confirmatory test aimed at a gain of 0.05 needs **several hundred to about a thousand independent targets**, not
tens. Few reviewed catalogues of a single sequence give that, so phase 6 needs either several sequences pooled under
one registration, or a registered effect size of at least 0.10 (200 to 300 targets). The data request of
`docs/data_access.md` should be sized by this table. A correlation rho close to 0.9 is what two closely related
methods (the hybrid and its base) give; rho 0.5 is closer to two unrelated methods.

## Cost and risks
About 90 s on two cores. The copula model ignores clustering and heterogeneous detectability; a design simulation
on the actual catalogue of each sequence is still required before its registration (protocol, Section 4.3).
