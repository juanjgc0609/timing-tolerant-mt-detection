Results of Task 2D: equal temporal tolerance for the comparators of 2C, ablations of
MT-FW and frozen split-conformal threshold (Appendix Y.21 of the proposal; tasks_and_real_data_plan.md, section
Task 2D).

Decision rules: PREREGISTRO_1j_1k.txt, section "Task 2D" (26 September 2026, 21:13 UTC, written
before writing the code, in response to the second independent audit). It is the first task with a wall-clock
time and hash: PREREGISTRO_hashes.txt stores the SHA-256 of the preregistration at 21:13:04 UTC (6e8510b3...bf88),
after the budget note at 02:41:50 UTC on the 27th (dc8cf708...e32b) and after the result paragraph at
03:54:55 UTC (b01a2eac...d9a8). The only addition to the 2D text before running is the "2D budget
note" (02:41:50 UTC, time from `date -u`), written before running seeds 2024 and 7331: the smoke test
estimated 0.64 h on 2 cores (1.28 h of CPU), below the 6 h threshold, and part (ii) ran with the 200
empty windows planned per cell and regime. The same block of the preregistration contains the erratum on
the time of 2B (E34). The code hashes (task2d_methods.py, task2d_run.py, task2d_validate.py,
analyze_2d.py and run_queue_2d.sh) were recorded at 02:41:59 UTC in task2d_code_hashes.txt, before
any result. The copy of PREREGISTRO_1j_1k.txt in this folder contains those sections and the paragraph
"Result of 2D".

Design in brief. Part (i): the 1200 windows of Task 2A (seeds 2024 and 7331, regimes thin and strong,
cells ladder_4, ladder_6 and unilateral_c40, networks 0 to 4, 20 windows per network), regenerated bit for bit as in
2C, with the same station corrections, the grid of 2A (675 nodes at 3 km, 31 origins at 0.1 s),
5 km exclusion and at most 5 events. MT-FW+SC and QM-3C+SC are the stored records of 2A; BOMP+SC,
MTSTACK+SC and GL-FW+SC, those of 2C. Reading at the peak ("PK"): for each station and phase, the sample of highest
three-component energy within +-4 samples (0.25 s) of the predicted one, the first in case of a tie; this is the rule
of node_vectors of MT-FW. New methods:
    BOMP-PK+SC      BOMP of 2C read at the peaks of the current residual; the entering block fixes its reading
                    times; score e/2 on entry; orthogonal refit of all atoms at their times;
                    Gram with the P/S cross term at the offset of the peaks; no cap
    MTSTACK-PK+SC   MTSTACK of 2C read at the peaks of the filter output; no deflation and no cap (it is
                    WAVE-MT without cap)
    GL-FW-PK+SC     group BLASSO of 2C with a growing dictionary: the Frank-Wolfe oracle evaluates the
                    certificate at the peaks of the current residual and the entering atom keeps those times;
                    same lambda path, stopping and score as 2C; gap with respect to the realised
                    dictionary
    MT-FW-EX+SC     MT-FW without peak search (exact predicted sample, with cap)
    MT-FW-NOCAP+SC  MT-FW with peak search and without cap
BOMP-PK+SC and MTSTACK-PK+SC are our own implementations from the published descriptions with the
tolerance of MT-FW added, not their authors' code. MT-FW-EX and MT-FW-NOCAP use the dense loop of MT-FW,
not the certified branch and bound.
Part (ii): split-conformal threshold for MT-FW+SC and QM-3C+SC, calibrated on the 5 calibration networks of
seed 999 of 2A per cell (same stations, delays and corrections), 40 source-free windows per network, 200
per cell and regime, with the clutter of the regime (new stream 2701). Score of a window = its maximum (first
score of MT-FW+SC; highest score of the events of QM-3C+SC with the frozen configuration of 2A, -inf if it does not
trigger). Threshold = the ceil((n+1)(1-alpha))-th smallest of the n = 200 maxima (181st for alpha 0.1, 151st for
0.25), frozen per cell, regime, method and alpha and applied without changes to the evaluation windows of 2A
(100 per cell and seed, 36 to 45 empty); an event is declared if its score exceeds the threshold. In thin
unilateral QuakeMigrate triggers in only 28 of the 200 empty windows and its threshold at alpha 0.25 is -inf.
Rules, per regime: 2D-a, MT-FW+SC - max(BOMP-PK+SC, MTSTACK-PK+SC) >= 0.05 at cap 0.25 in at least 2 of 3
cells in both seeds; if it fails, the advantage of MT-FW over [129] and [122] is attributed to the peak search
and to no other element of the method, and the proposal declares it. 2D-b, MT-FW+SC - MT-FW-EX+SC >= 0.05 at cap 0.25
in at least 2 of 3 cells in both seeds (descriptive: MT-FW-NOCAP+SC - MT-FW+SC and GL-FW-PK+SC - MT-FW+SC with
CI). 2D-c, with the frozen thresholds, false events of MT-FW+SC per empty window <= 2 alpha and recall of MT-FW+SC
- QM-3C+SC >= 0.05 in at least 2 of 3 cells in both seeds, at both alpha. Statistics: hierarchical bootstrap
of 2000 replicates (networks of the cell and, within each network, windows), operating point and best
comparator re-chosen in each replicate of part (i), fixed thresholds in 2D-c; Holm over 36 comparisons
per regime, one-sided bootstrap p at level 0.025. Readings fixed before the results (docstring of
analyze_2d.py): rules 2D-a and 2D-b on point estimates; in 2D-c, score strictly greater than the
threshold, "rate" = false events per empty window, per-cell reading primary and strict reading reported, without calibration
uncertainty in the CI.

Outputs
    analyze_2d_out.txt            controls, tables 1 to 8 per regime (recall; 2D-a and 2D-b with hierarchical CI
                                  and by windows; NOCAP and GL-FW-PK; effect of the reading at the peak;
                                  calibration thresholds; evaluation of 2D-c; Holm), exploratory table with
                                  the 30 networks of 2C (E3) and verdicts; primary reading of the task
    analyze_2d_out.json           verdicts, counts, fraction of replicates and frozen thresholds
    task2d_report.txt             report (in English): verdicts, design as implemented, validation,
                                  controls and kernels, compute, key tables, exploratory, interpretation and
                                  deviations
    task2d_validation.txt         validations V1 to V5 of the methods on seed 999 and V6 (descriptive)
    task2d_smoke.txt              smoke test on seed 999: times, controls and compute estimate
                                  that decided the budget rule
    task2d_explore_e1.txt         exploratory E1: false events of MT-FW+SC with the conformal thresholds in windows
                                  with events and their distance to the true source
    task2d_explore_e2.txt         exploratory E2: V6 with 90 windows per cell (injected time error)
    task2d_code_hashes.txt        SHA-256 of the code, 02:41:59 UTC, before any result
    PREREGISTRO_1j_1k.txt         copy of the preregistration with the section Task 2D, its note and its result
    PREREGISTRO_hashes.txt        SHA-256 hashes of the preregistration with wall-clock time

The 90 per-network files (task2d_parts/p1_seed{s}_{regimen}_{celda}_net{j}.json, 60 from part (i), and
task2d_parts/p2_seed999_..., 30 from part (ii)), the 60 of task2d_explore_parts/ and the logs of the
runs (log2d.txt, log2d_validate.txt, log2d_smoke.txt, log2d_explore.txt; the error logs, empty) remain in
the working folder of the scripts; they are not copied here. Each part file records the OpenBLAS kernel with
which it was computed.

Verdict (analyze_2d_out.txt)
    rule             thin                                strong
    2D-a             not met (2 of 3 and 1 of 3)         not met (1 of 3 and 0 of 3)
    2D-a with Holm   not met (0 of 3 and 0 of 3)         not met (0 of 3 and 0 of 3)
    2D-b             met (3 of 3 and 3 of 3)             met (3 of 3 and 3 of 3)
    2D-b with Holm   not met (1 of 3 and 0 of 3)         not met (0 of 3 and 1 of 3)
    2D-c per cell    met                                 met
    2D-c strict      not met                             met
    2D-c with Holm   not met                             not met
No case came within 0.005 of the 0.05 threshold. Re-applied in each hierarchical replicate: 2D-a in 11 %
(thin) and 0.1 % (strong) of the replicates; 2D-b in 99 % and 95 %; 2D-c per cell in 76 % and 96 %.

2D-a, MT-FW+SC minus BOMP-PK+SC at cap 0.25 (the best comparator is BOMP-PK+SC in all cells), 95
hierarchical CI:
    cell (seed)           thin                       strong
    4 st. (2024)          +0.080 [-0.136, +0.282]    +0.091 [-0.048, +0.257]
    6 st. (2024)          +0.042 [-0.092, +0.224]    +0.000 [-0.138, +0.145]
    unilateral (2024)     +0.056 [+0.000, +0.265]    -0.036 [-0.130, +0.108]
    4 st. (7331)          +0.090 [-0.029, +0.256]    -0.051 [-0.206, +0.095]
    6 st. (7331)          -0.029 [-0.125, +0.047]    -0.155 [-0.277, +0.032]
    unilateral (7331)     +0.025 [-0.084, +0.155]    -0.053 [-0.158, +0.030]
All CIs contain 0. MT-FW+SC - MTSTACK-PK+SC at 0.25: +0.21 to +0.35 (thin), +0.01 to +0.16 (strong).

2D-b, MT-FW+SC minus MT-FW-EX+SC at cap 0.25, 95 hierarchical CI:
    cell (seed)           thin                       strong
    4 st. (2024)          +0.230 [+0.028, +0.388]    +0.156 [+0.041, +0.355]
    6 st. (2024)          +0.153 [-0.118, +0.348]    +0.149 [+0.000, +0.327]
    unilateral (2024)     +0.310 [+0.154, +0.515]    +0.241 [+0.025, +0.433]
    4 st. (7331)          +0.218 [-0.038, +0.378]    +0.064 [-0.096, +0.240]
    6 st. (7331)          +0.412 [+0.108, +0.597]    +0.167 [-0.024, +0.410]
    unilateral (7331)     +0.407 [+0.212, +0.619]    +0.293 [+0.098, +0.472]
Positive in 12 of 12; lower bound > 0 in 7 of 12; Holm rejects in 2 of 12.

Descriptive at 0.25. MT-FW-NOCAP+SC - MT-FW+SC: thin -0.080, -0.042, -0.070 / -0.090, +0.044, -0.025; strong
-0.091, +0.000, +0.036 / +0.051, +0.131, +0.053 (all CIs contain 0). GL-FW-PK+SC - MT-FW+SC: thin
-0.276, -0.167, -0.155 / -0.244, -0.147, -0.296; strong -0.156, -0.203, -0.108 / +0.013, -0.119, -0.133 (below
in 11 of 12). Effect of the reading at the peak (PK - exact version of 2C): BOMP +0.04 to +0.51, MTSTACK
+0.05 to +0.37, GL-FW +0.03 to +0.44 (12 of 12 each); GL-FW-PK+SC - BOMP-PK+SC -0.04 to -0.27 (12 of 12);
max(BOMP-PK+SC, MTSTACK-PK+SC) - QM-3C+SC +0.17 to +0.49 (12 of 12).

2D-c, recall MT-FW+SC - QM-3C+SC and false events of MT-FW+SC per empty window, alpha 0.1 / alpha 0.25:
    cell (seed)           thin                                  strong
    4 st. (2024)          +0.437, 0.179 / +0.448, 0.410         +0.091, 0.025 / +0.260, 0.200
    6 st. (2024)          +0.347, 0.045 / +0.375, 0.159         +0.284, 0.000 / +0.351, 0.023
    unilateral (2024)     +0.296, 0.070 / +0.310, 0.186         +0.337, 0.048 / +0.398, 0.190
    4 st. (7331)          +0.436, 0.179 / +0.449, 0.410         +0.077, 0.000 / +0.167, 0.108
    6 st. (7331)          +0.485, 0.222 / +0.485, 0.422         +0.262, 0.028 / +0.345, 0.167
    unilateral (7331)     +0.358, 0.053 / +0.358, 0.263         +0.480, 0.075 / +0.480, 0.375
The strict reading fails in thin because of one false event: 7331, 6 st., alpha 0.1, 10 false events in 45 empty windows =
0.222 > 0.2. The recall condition alone is met in 12 of 12 cases per regime (+0.08 to +0.49). False events per
window over all windows: MT-FW+SC thin 0.73-1.18 (alpha 0.1) and 0.98-1.46 (alpha 0.25), strong 0.05-0.56
and 0.26-0.94; QM-3C+SC thin 0.15-0.23 and 0.16-0.31, strong 0.14-0.23 and 0.22-0.34.

Exploratory diagnosis (post hoc, not preregistered; it supports no verdict). E1: in the windows of 2A
with events, which the calibration on empty windows does not see, MT-FW+SC gives with the conformal thresholds 1.21 to 1.87 false events per
window in thin at alpha 0.1 and 1.41 to 2.19 at 0.25, 78 to 92 % of them at 5-20 km from a true source, against 0.18 to
0.36 for QM-3C+SC; in strong 0.07 to 0.88 and 0.35 to 1.32, against 0.17 to 0.38. (task2d_report.txt says 79-92 %;
with the counts of task2d_explore_e1.txt the minimum is 73 of 93 = 78.5 %.) Conformal thresholds of MT-FW+SC
23.7-29.8 (thin) and 36.2-52.1 (strong), below those of 2A-e (34.6-42.9 and 44.7-53.2). E3: the same thresholds
on the 30 networks per cell of 2C give MT-FW+SC - QM-3C+SC +0.36 to +0.47 (thin) and +0.03 to +0.40 (strong), with
false events per empty window near alpha and below 2 alpha in all of them. E2: with an injected time error of
+-0.2 s per station and phase and the sources at the nodes, the exact-time methods lose 0.16 to 0.48 and the
PK variants -0.03 to +0.15, and the latter outperform their exact version in 12 of 12 cases (+0.04 to +0.38); without error, at
6 and 8 stations the peak search costs recall (MT-FW 0.83-0.87 against 0.97 for MT-FW-EX).

Reading. 2D-a fails: with the same peak search of +-0.25 s, the block inversion in the form of Vera
Rodriguez et al. (2012) matches MT-FW; as preregistered, the advantage of MT-FW over [129] and [122] is
attributed only to the peak search, and the reading of 2C that attributed it to the elements specific to MT-FW does not
hold (E35 in artifact_manifest_and_error_log.md). 2D-b is met: the peak search matters. The cap does not
matter: without it, MT-FW is in essence BOMP with reading at the peak. The group BLASSO on the grid with the same
peak search remains below MT-FW. The empirical advantage over the pick-based associators and QuakeMigrate
belongs to the family of moment-tensor detectors with temporal tolerance and deflation, not to MT-FW. 2D-c
is met per cell, but the threshold on empty windows is not a deployable operating point: in windows with
events it leaves 1.2 to 2.2 false events per window in thin, 3 to 8 times the false events of QuakeMigrate over all
windows. With 5 networks per cell no rule survives Holm. Consequence for the thesis: the distinctive methodological
contribution becomes the principled treatment of temporal tolerance (group BLASSO with
one shift per station and phase as a nuisance parameter, off the grid by sliding
Frank-Wolfe), which only counts if it outperforms the heuristic peak search and BOMP-PK.

Controls. Part (i): network, delays, catalogue and corrections identical bit for bit to 2A (60/60 networks); truth,
origins and number of sources identical (1200/1200 windows, against 2A and 2C); MT-FW+SC recomputed equal to what was
stored (120/120 windows). Part (ii): 30 calibration networks with stations and delays equal to 2A,
identical corrections (maximum difference 0), and the pipeline reproduces MT-FW+SC and QM-3C+SC of 2A in windows
0-1 of each network (60/60 each). Validation (task2d_validation.txt): tabulated peaks = node_vectors_c (0
discrepancies in 669 600 to 1 339 200 readings per cell); with the reading at the peak disabled, records of 2C
exact (12/12); atoms with shifted times <= 1.2e-15; independent paths <= 3.3e-15; GL-FW-PK with
relative gap <= 3.6e-9 and certificate / lambda <= 1.000000000. No control failed; no
task2d_CONTROL_FAILED.txt was written.

Edit after the results. A diagnostic line of analyze_2d.py failed (the lists of
shifts of BOMP-PK have length 2 ns per atom and had to be flattened) and it was changed after the part
files existed. Only that line changes; no rule, threshold, bootstrap or table. Hash before:
110aa893...73a7 (task2d_code_hashes.txt, 02:41:59 UTC); after: 8dd5013b...d2b4. The paragraph "Result of 2D"
of the preregistration declares it.

Declared deviations (task2d_report.txt, section 9): the PK variants carry no cap and use the Gram at the
peak times with the P/S cross term at the offset of the peaks (our reading of "identical to 2C except
... the sample ... at the peak"); in GL-FW-PK atoms of the same node and origin with different times may
coexist, the record deduplicates them by the 5 km rule, the gap is with respect to the realised dictionary and the
per-station oracle is not a global maximiser over shifts (FW steps per lambda bounded at 400,
never reached); MT-FW-EX and MT-FW-NOCAP use the dense loop; readings of 2D-c fixed before the results;
new stream 2701 for the empty windows; hierarchical bootstrap with only 5 networks per cell (resolution of p
1/2001; the first Holm step reaches at most 36/2001 = 0.018); the Holm family and level are our reading;
the default kernel of this machine changed to Cooperlake, so it is always fixed and checked; V6 is small and
descriptive (one of 12 comparisons negative, repeated in E2).

Reproduction (from the scripts folder, with task2a_parts/, the evaluation JSON files of 2A and 2C,
task2c_parts/ and qm_recal_2a.json in the same folder):
  OPENBLAS_CORETYPE=Prescott python task2d_run.py validate                 (writes task2d_validation.txt)
  OPENBLAS_CORETYPE=Prescott python task2d_run.py smoke --nwin 3           (writes task2d_smoke.txt)
  OPENBLAS_CORETYPE=Prescott python task2d_run.py run --part 2 --workers 2
  OPENBLAS_CORETYPE=SkylakeX python task2d_run.py run --part 1 --workers 2
  OPENBLAS_CORETYPE=Prescott python task2d_run.py run --part 1 --workers 2
  python analyze_2d.py                                                     (writes analyze_2d_out.txt and .json)
  OPENBLAS_CORETYPE=SkylakeX python task2d_explore.py e1 --workers 2       (exploratory)
  OPENBLAS_CORETYPE=Prescott python task2d_explore.py e1 --workers 2
  python task2d_explore_sum.py                                             (writes task2d_explore_e1.txt)
  OPENBLAS_CORETYPE=Prescott python task2d_explore.py e2 --workers 2       (writes task2d_explore_e2.txt)
run_queue_2d.sh launches the three run invocations and run_queue_2d_explore.sh the three of task2d_explore.py. Each run process only launches the
networks whose required kernel is its own (SkylakeX for the 18 networks of 7331 that 2A computed with it, Prescott for
the rest), checks the loaded kernel with threadpoolctl and records it. run resumes from the parts already written
in task2d_parts/; if a network does not pass the exact control it writes task2d_CONTROL_FAILED.txt and stops. On another
CPU the available kernels may be different: the exact control of part (i) only passes with the kernel that
reproduces bit for bit each network of 2A.
Times: part (i) 3931 s of CPU (per window: BOMP-PK 0.47 s, MTSTACK-PK 0.013 s, GL-FW-PK 1.76 s, MT-FW-EX
0.23 s, MT-FW-NOCAP 0.76 s); part (ii) 1182 s (MT-FW+SC 0.20 s and QM-3C+SC 0.72 s per empty window); 1.42 h of
CPU in total, from 02:42 to 03:25 UTC with 2 processes.
Dependencies: those of Tasks 2A to 2C (results/task2a/README.txt, results/task2b/README.txt,
results/task2c/README.txt), because task2d_run.py imports task2b_ablation.py, task2a_confirm.py and task2c_run.py, and
task2d_methods.py imports task2c_methods.py, certified_fw.py and task1x_mtglrt.py; QuakeMigrate for QM-3C+SC in
part (ii); threadpoolctl to check the kernel.
