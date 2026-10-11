Results of Task 2E: the group BLASSO with one time shift per station and phase, solved off
the grid by sliding Frank-Wolfe (TS-SFW+SC), against the block search with the same temporal
tolerance (BOMP-PK+SC). It is the test of the methodological hypothesis G1-M (Appendix Y.22 of the supplementary volume;
tasks_and_real_data_plan.md, section Task 2E).

Decision rules: PREREGISTRO_1j_1k.txt, section "Task 2E" (27 September 2026, 14:44:14 UTC, wall-clock
time), written before writing any code of 2E. Its SHA-256, b40f0c6e...e988, was stored in PREREGISTRO_hashes.txt
and delivered to the user through the chat at the moment of registration: it is the first task whose content fixes a record
outside the container (not public and not held by third parties). The only addition before evaluating is the "2E calibration
note" (15:38:20 UTC, f = 0.01, hash fada69b1...9271), written before any computation on seeds 2024 and 7331.
The budget rule (10 h on 2 cores) was not triggered: estimate 1.98 h. Code hashes (task2e_methods.py,
task2e_validate.py, task2e_run.py, analyze_2e.py and run_queue_2e.sh) in task2e_code_hashes.txt at 15:38:31 UTC;
first evaluation part at 15:39:40 UTC and last at 17:39:19 UTC; sha256sum -c correct after the
run. Paragraph "Result of 2E" at 17:45:17 UTC with another hash (ebbc35de...af1b). The copy of
PREREGISTRO_1j_1k.txt in this folder contains those sections and the later ones.

Design in brief. 30 networks per cell, seed (2024, 7331) and regime (thin, strong): the 5 networks of 2A (20 windows)
and the 25 of 2C part (ii) (12 windows), 400 windows per cell and 4800 in total, regenerated bit for bit with the
OpenBLAS kernel of 2C (SkylakeX on the 18 networks of 7331 of 2A that were computed with it, Prescott on the rest). Same
station corrections, grid (675 nodes at 3 km, depths 5, 9 and 13 km), origins, tolerance of +-0.25 s,
5 km exclusion, at most 5 events and hit rule (epicentre within 5 km) as 2A-2D. Methods:
    TS-SFW+SC         continuous atom (x, y, z, t0 and one shift per station and phase bounded to +-0.25 s),
                      prediction in the matched-filter domain with the continuous autocorrelation, exact
                      band-limited correlation, group norm whitened per block as in GL-FW of 2C. Per iteration
                      (k <= 5): (a) oracle = statistic of MT-FW-NOCAP on the grid, excluding 5 km around
                      the current positions; (b) L-BFGS-B refinement of the new atom (<= 50 it.); (c)
                      sliding alternating the convex step in amplitudes (relative gap <= 1e-6) and L-BFGS-B in
                      positions, origins and shifts of all atoms (<= 50 it. per step, <= 5
                      alternations), plus a final convex step; (d) residual. lambda = f * s_1. Score: cumulative
                      minimum of the oracle; hit with the final position.
    BOMP-PK-SLIDE+SC  the same with lambda = 0 (least-squares refit)
    BOMP-PK+SC        code of 2D unchanged (recomputed equal to 2D in 1200/1200 windows)
    MT-FW+SC          stored records of 2A and 2C
Calibration of f on seed 999 (600 windows of 2A): recall at 0.25 of 0.4749, 0.4644, 0.2510 and 0.1841 for
f = 0.01, 0.03, 0.1 and 0.3; f = 0.01 was frozen (lambda, median 0.045 of the certificate of the first atom; no
atom at zero; at f = 0.3, 88 %). Rules, per regime: 2E-a (primary), TS-SFW+SC - BOMP-PK+SC >= 0.05 at cap
0.25 with lower bound of the 95 hierarchical CI > 0 in at least 2 of 3 cells in both seeds; strict reading with
Holm over the 6 cell x seed tests (one-sided 0.025). 2E-b (secondary), TS-SFW+SC - BOMP-PK-SLIDE+SC >= 0.05
at 0.25 in at least 2 of 3 cells in both seeds (point). Registered interpretation: if 2E-a fails in both
regimes, G1-M is refuted in this synthetic design. Statistics: exact FROC, operating point per method
on the evaluation, hierarchical bootstrap of 2000 replicates (networks with replacement and windows within each network),
same resamples for all methods, p = (1 + #{d* <= 0}) / 2001.

Outputs
    analyze_2e_out.txt            controls, tables 1 to 7 per regime (recall; 2E-a; 2E-b; descriptive
                                  differences; Holm; location error; convergence diagnostics) and
                                  verdicts
    analyze_2e_out.json           verdicts, points, CIs and counts
    task2e_report.txt             report (in English): verdicts, implemented design, validation, calibration,
                                  times, tables, deviations D1-D12, interpretation
    task2e_validation.txt         validations V1 to V5 on seed 999
    task2e_calibration.txt        calibration of f on seed 999
    task2e_frozen.json            frozen f (15:37:42 UTC)
    task2e_smoke.txt              time estimate that decided the budget rule
    task2e_explore.txt            exploratory, written after the results: displacement of the atoms
                                  by the sliding, changes of hit, actual times
    task2e_code_hashes.txt        SHA-256 of the code, 15:38:31 UTC
    PREREGISTRO_1j_1k.txt         copy of the preregistration
    PREREGISTRO_hashes.txt        SHA-256 hashes of the preregistration with wall-clock time
    verif_2e/                     independent verification (task2e_verification.txt and its scripts *_verif.py with
                                  their outputs)

The 390 parts (task2e_parts/cal_seed999_* and ev_seed{2024,7331}_*) and the run logs (task2e_cal.log,
task2e_run.log, task2e_analyze.log; analyze_2e_err.txt empty) remain in the working folder of the scripts; they are not
copied here. Each part records the OpenBLAS kernel with which it was computed.

Verdict (analyze_2e_out.txt)
    rule             thin                                strong
    2E-a             not met (0 of 3 and 0 of 3)         not met (0 of 3 and 0 of 3)
    2E-a with Holm   not met (minimum adj. p 0.81)       not met (minimum adj. p 0.58)
    2E-b             not met (0 of 3 and 0 of 3)         not met (0 of 3 and 0 of 3)
No case came within 0.005 of a threshold or of 0. As registered, G1-M is refuted in this synthetic
design.

2E-a, TS-SFW+SC minus BOMP-PK+SC at cap 0.25, 95 hierarchical CI:
    cell (seed)           thin                       strong
    4 st. (2024)          +0.039 [-0.025, +0.122]    +0.000 [-0.033, +0.030]
    6 st. (2024)          +0.029 [-0.025, +0.087]    +0.033 [-0.013, +0.083]
    unilateral (2024)     +0.025 [-0.042, +0.094]    +0.016 [-0.026, +0.059]
    4 st. (7331)          +0.016 [-0.038, +0.094]    +0.003 [-0.027, +0.039]
    6 st. (7331)          +0.029 [-0.024, +0.078]    +0.015 [-0.024, +0.058]
    unilateral (7331)     +0.015 [-0.035, +0.091]    +0.018 [-0.044, +0.080]
The 12 points are >= 0, none reaches 0.05 and all CIs contain 0.

2E-b, TS-SFW+SC minus BOMP-PK-SLIDE+SC at 0.25: thin -0.003, -0.003, +0.000 / +0.013, -0.006, -0.012; strong
-0.003, +0.007, -0.020 / -0.003, -0.006, -0.009 (all CIs contain 0): the group penalty is inert with
the calibrated f.

Descriptive at 0.25. BOMP-PK-SLIDE+SC - BOMP-PK+SC (sliding only): +0.003 to +0.042, all CIs contain 0.
TS-SFW+SC - MT-FW+SC: -0.026 to +0.049, all CIs contain 0. Median epicentral error of the hits:
TS-SFW+SC 1.79 km against BOMP-PK+SC 1.88 km (thin) and 2.01 against 2.15 km (strong); the continuous solutions end up
0.3 to 1.3 km too deep. Convergence: maximum relative gap 5.6e-10; 73-83 % of the sliding
steps and 83-97 % of the refinements stop at 50 iterations; 24-35 % of the shifts at the
bound. Time per window: TS-SFW+SC 1.19 s, BOMP-PK-SLIDE+SC 1.14 s, BOMP-PK+SC 0.42 s (mean).

Validation (task2e_validation.txt). V1 continuous identities <= 8.9e-16; V2 gradients against central differences
1.4e-9 (18 cases); V3 convex step against an explicit dictionary, against FISTA and against least squares with
lambda = 0; V5 without refinement or sliding and with lambda = 0 the first atom is that of BOMP-PK (24/24). V4 (noise-free
source off the grid): from the truth the optimiser stays (94/94), at 0.2-0.4 km it recovers 85/94, but
from the grid node only 12/90 sources to < 0.5 km and < 0.02 s (cycle skips of ~0.45 s). The registered procedure was
followed because the registered method is a local refinement; it is declared as a property of the method.

Independent verification (verif_2e/task2e_verification.txt; new computation only on seed 999; no file
of 2E or of the registration modified). Verdict: RELIABLE NEGATIVE.
    1. Code: at a node with the peak shifts, the continuous atom reproduces the block of BOMP-PK in
       120 atoms (prediction 1.4e-15, Gram 9.0e-16, correlation 5.7e-16, objective 9.9e-16); objective against an
       independent implementation 2.3e-13; gradient 4.5e-8; injected shifts recovered to 1e-4 s.
       (check1_verif.py)
    2. Numbers: own FROC without importing the project reproduces exactly Table 1 and the 24 points of 2E-a and
       2E-b; own bootstrap with another random generator. (froc_rederive_verif.py, boot_verif.py)
    3. Budget: with 10 times the iterations, 4 times the alternations and refinement with 9 starts (270
       windows, 208 events), no call reaches the cap and the objective drops by a median of 1.9 %, but the
       recall at 0.25 changes -0.014 [-0.052, +0.038] (6 events gained, 4 lost): there is no budget
       disadvantage. (check3_budget_verif.py)
    4. Ceiling: CEIL-G, CEIL-GD and CEIL-C on the 600 calibration windows with the exact clutter and noise part
       (regeneration 600/600). At 0.25 BOMP-PK lies 0.16-0.37 below CEIL-G and 0.08-0.24 below CEIL-GD per cell
       (pooled: thin +0.156 [+0.042, +0.286], strong +0.162 [+0.037, +0.301]); CEIL-C ~ CEIL-G (0.699 against
       0.692). The detectable events never hit have the nearest atom at a median of 7.3 km
       (BOMP-PK) and 6.8 km (TS-SFW), 80 / 85 % at 5-10 km. (check4_ceiling_verif.py, analyze_verif.py,
       anatomy_verif.py, reloc_verif.py)
    5. Hash chain: the 6 hashes of the time are growing prefixes of the registration; the text of 2E is identical to the one
       hashed at 14:44:14 UTC; the calibration note precedes the first evaluation part.
       (hashchain_verif.py, parts_mtimes_verif.txt)

Reading. 2E-a and 2E-b fail: the group BLASSO with off-grid shifts per station and phase is not
a better detector than the block search with the same tolerance. No convex formulation of the thesis (on the
grid or off it, with or without tolerance) outperforms BOMP-PK on these synthetics. The margin exists (0.08-0.24 below the
realistic ceiling), but it lies in events absorbed by atoms misplaced at 5-10 km, which the local sliding
(median 0.7 km) does not reach. The subsequent exploration (explore_contrib/) looked for where that margin lies.

Declared deviations (task2e_report.txt, section 7): units of lambda (value of the oracle in the convention of
MT-FW, D1); exclusion around the slid position (D2); bounds, internal scale and cap of 50 iterations of the
refinement not fixed by the registration (D3); refinement energy without the P/S cross term (D4); t0 gauge and
common shift not fixed (D5); continuous interpolation as exact band-limited correlation (D6); the
requested V4 criterion is not met from the grid (D7); time estimate with the full calibration (D8);
a wording erratum in the calibration note ("both cells and regimes" for "the three cells and both
regimes"), not corrected so as not to change the registration (D9); the SkylakeX kernel also changed the OpenBLAS of scipy
inside L-BFGS-B (D10); analyze_2e.py tested before the evaluation only on the calibration parts (D11);
task2e_explore.py exploratory (D12).

Reproduction (from the scripts folder, with task2a_parts/, task2c_parts/ and the records of 2D):
  OPENBLAS_CORETYPE=Prescott python task2e_run.py validate                 (writes task2e_validation.txt)
  OPENBLAS_CORETYPE=Prescott python task2e_run.py calibrate --workers 2
  python task2e_run.py cal-summary                                         (writes task2e_calibration.txt and
                                                                            task2e_frozen.json)
  python task2e_run.py estimate                                            (writes task2e_smoke.txt)
  OPENBLAS_CORETYPE=Prescott python task2e_run.py run --workers 2
  OPENBLAS_CORETYPE=SkylakeX python task2e_run.py run --workers 2
  python analyze_2e.py                                                     (writes analyze_2e_out.txt and .json)
  python task2e_explore.py                                                 (exploratory; task2e_explore.txt)
run_queue_2e.sh launches the two run invocations. task2e_run.py run refuses to evaluate if the hash of the registration is not the
last one in PREREGISTRO_hashes.txt and it does not contain the note with the frozen f; each process only launches the networks whose
required kernel is its own and checks it with threadpoolctl; it resumes from the parts already written.
Times: evaluation 2.0 h of wall time in 2 processes (15:39:40-17:39:19 UTC); calibration 15:09:35-15:37:42 UTC.
Dependencies: those of Tasks 2A to 2D (results/task2a/README.txt to results/task2d/README.txt), because
task2e_methods.py imports task2d_methods.py, task2c_methods.py, certified_fw.py, task1x_mtglrt.py and task0_common.py,
and task2e_run.py imports task2b_ablation.py, task2a_confirm.py, task2c_run.py and task2d_run.py; scipy (L-BFGS-B) and
threadpoolctl.
