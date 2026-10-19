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
`check_greedy` (the longest) was still running at integration, and its stored output is that of 18 October.

## Cost and risks
Minutes per script except `check_greedy.py` (61,425 hypotheses x 216 replicates). The script comments are in Spanish,
the language the theory was first written in.
