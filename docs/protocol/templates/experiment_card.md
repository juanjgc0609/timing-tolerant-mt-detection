# Experiment card: <ID> <short title>

Stage: S0 / S1 / S2 / S3 / S4 / S5 / S6 (current)
Owner decisions: <what the author decided, and when>

## Question
One sentence. What do we not know that this experiment will tell us?

## Mechanism tested
Which component of the method, or which failure of a comparator, does the experiment isolate? What is held fixed?

## Hypothesis and refutation
H: <statement>. Refuted if: <observable result, with its threshold>.

## Data
Datasets and partitions (link the data cards). Synthetic generator and seeds, or real records. What is development and
what is evaluation.

## Methods and comparators
Proposed configuration; ablations; comparators at their best (Section 3 of the protocol).

## Endpoint and analysis
Primary endpoint, hit rule, false-detection budget, statistics, simulated false-positive rate and power.

## Gate (written before the exploration starts)
Condition on development data that allows registration; what happens if it is not met.

## Cost and risks
Compute, data access, time; what could make the result uninterpretable.
