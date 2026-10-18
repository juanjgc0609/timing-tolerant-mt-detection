# Verification of task <ID>

Verifier: <separate AI session / person>, date. Work folder. Scripts.

## Summary
Verdict per registered rule: CONFIRMED / NOT CONFIRMED / CANNOT VERIFY, with the numbers.
Defects by severity: [REG] changes a registered number; [DESC] changes a descriptive number; [COSM] no number.

## Checks
1. Freeze integrity and log hash chain.
2. Partition integrity: no evaluation sample read before the freeze, margins included.
3. Data handling: catalogues, coordinates, components, time base.
4. Independent recomputation of every registered number from the raw outputs.
5. Implementation against the registration.
6. Bugs that could change the result, with their measured effect.
7. Number re-derived by a person (for tests on real data): which, and the result.
