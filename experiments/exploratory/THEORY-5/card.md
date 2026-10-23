# Experiment card: THEORY-5. Numerical checks of Appendix AB

Stage: S1 (done, 18 October 2026; integrated 19 October 2026)
Owner decisions: none needed; it checks proofs, it does not test the method.

## Question
Do the statements of Appendix AB (identifiability of overlapping events, greedy recovery with timing error, the
adaptive robust statistic, the components of the new model) hold numerically in the conventions of the repository?

## Mechanism tested
Each script isolates one statement: `check_identifiability.py` (Parts AB.A), `check_greedy.py` and `diag_greedy.py`
(AB.B), `check_adaptive.py` (AB.C), `check_extensions.py` (AB.D). `tail_check_c1.py` is `code/theory_tail_check.py`
with the corrected bound of AB.C.1 printed next to the bound of AA.9.5(b); it reproduces the quantiles of
`results/theory/tail_check_out.txt` exactly.

## Hypothesis and refutation
Not a hypothesis test. A statement is refuted if a script finds a violation of a proved bound.

## Data
Simulation only. `lib5.py` replicates the conventions of the repository (2 Hz Ricker, 16 Hz, VS 3 km/s,
VP 1.73 VS, the 8-station `spread` network, the radiation operator of `code/task1x_mtglrt.py`) with CP = CS = 1.
Fixed seeds inside each script.

## Endpoint and analysis
Every output is stored next to its script (`*_out.txt`). Re-run on 19 October 2026 with NumPy 2.x and SciPy: the
outputs of `check_identifiability`, `check_adaptive`, `check_extensions` and `diag_greedy` are identical;
`check_greedy` (the longest) was identical in the re-implementation check of 20 October once the hand-written note
below was removed.

## Environment
The stored outputs were produced and re-run with Python 3.13.16, NumPy 2.5.3 and SciPy 1.18.1. Under the frozen core
environment (`environments/core.lock.txt`: Python 3.12, NumPy 1.26.4, SciPy 1.17.1, `OPENBLAS_CORETYPE=Prescott`) the
re-implementation check of 20 October 2026 found `tail_check_c1`, `check_identifiability` and `check_extensions`
identical, `check_adaptive` identical up to rounding (7.9e-16 against 8.1e-16), and `diag_greedy` different: median
0.188, 90th percentile 0.313, maximum 0.356 of the residual peak amplitude, against 0.157, 0.316 and 0.362 stored,
because the greedy step picks other nodes in some of its 6 replicates. Appendix AB reports the range.

## Note on `check_greedy_out.txt`
The B.4 line of that output measures the residual in a window that can fall on another arrival; the valid figure of
the deflation residual is the one of `diag_greedy_out.txt`. (This note was first appended by hand to the stored output;
it was moved here so that every stored output is exactly what its script prints.)

## Cost and risks
Minutes per script except `check_greedy.py` (61,425 hypotheses x 216 replicates). The script comments are in Spanish,
the language the theory was first written in.
