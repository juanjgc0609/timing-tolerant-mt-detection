Results of Task 2C: closest comparators, the convex method of the thesis on the grid and
per-network replication (Appendix Y.20 of the proposal; tasks_and_real_data_plan.md, section Task 2C).

Decision rules: PREREGISTRO_1j_1k.txt, section "Task 2C" (26 September 2026, 15:05 UTC, written
before implementing it and before seeing any of its results, in response to the independent audit of 56/100).
The only addition is the "Nota de presupuesto de 2C" (19:56 UTC), written before running seeds 2024 and
7331: the smoke test estimated 0.79 h on 2 cores, below the threshold of 8 h, and part (ii) ran with
the 25 new networks planned. The copy of PREREGISTRO_1j_1k.txt in this folder contains those sections and the
paragraph "Result of 2C".

Design in brief. Part (i): the 1200 windows of Task 2A (seeds 2024 and 7331, regimes thin and strong,
cells ladder_4, ladder_6 and unilateral_c40, networks 0 to 4, 20 windows per network), regenerated bit for bit with the
functions of task2a_confirm.py; MT-FW+SC and QM-3C+SC are the stored records of 2A. New methods, all
in the matched-filter domain, read at the exact predicted sample (without peak search or cap), with
the same dictionary of 6-column blocks (one per moment-tensor component) on the grid of
2A (675 nodes at 3 km, 31 origins at 0.1 s), the same station corrections, exclusion of 5 km and at
most 5 events:
    BOMP+SC      block orthogonal matching pursuit in the form of Vera Rodriguez, Sacchi and Gu (2012): it chooses
                 the block with the largest energy of the projection of the residual and refits by least
                 squares all the active blocks; score = that energy / 2 on entry
    MTSTACK+SC   mechanism-corrected stack in the form of Anikiev et al. (2014): least-squares tensor
                 at each node and time, value ||G m||^2 / 2, maxima suppression at 20 km, without
                 deflation
    GL-FW+SC     group BLASSO on the grid over the BOMP dictionary (whitened group norm per
                 block), by fully corrective Frank-Wolfe with duality gap; path of 30
                 geometric lambdas from lambda_max to lambda_max/30; score = the largest lambda at which the
                 block enters the support
BOMP+SC and MTSTACK+SC are our own implementations from the published descriptions, not the code of
their authors. Part (ii): 25 new networks per cell, seed and regime (indices 5 to 29, with the seed
construction of 2A), 12 windows per network, with MT-FW+SC, QM-3C+SC, BOMP+SC and MTSTACK+SC; with the 5 networks of 2A
there are 30 networks and 400 windows per cell.
Rules, per regime: 2C-a, MT-FW+SC - max(BOMP+SC, MTSTACK+SC) >= 0.05 at cap 0.25 in at least 2 of 3 cells
on both seeds (part i). 2C-b, descriptive: GL-FW+SC - MT-FW+SC at both caps (ties if |dif| < 0.05).
2C-c, primary for part (ii): lower bound of the CI 95 of the hierarchical bootstrap (networks and then windows,
2000 replicates, operating point re-chosen in each replicate) of MT-FW+SC - QM-3C+SC at cap 0.25 > 0 and
point estimate >= 0.05 in at least 2 of 3 cells on both seeds; the same against BOMP+SC and
MTSTACK+SC and at cap 0.1, with Holm over the 36 comparisons of each regime. Endpoint as in 2A: recall at
<= 0.1 and <= 0.25 realised false events per window on the exact FROC curve, operating point chosen on
the evaluation itself, rules on point estimates.

Outputs
    analyze_2c_out.txt            controls, tables 1 to 6 per regime (recall; 2C-a with intervals by
                                  windows and by network; 2C-b with GL-FW diagnostics; recall with 30 networks;
                                  differences with hierarchical CI, p and Holm p; new networks only,
                                  exploratory) and verdicts; primary reading of the task
    analyze_2c_out.json           rules, counts, fraction of replicates, 2C-b per case and Holm rejections
    task2c_report.txt             report (in English): verdicts, design as implemented, review
                                  of the inherited code, validation, controls, computation, key tables,
                                  exploratory diagnostic, interpretation and deviations
    task2c_validation.txt         validations V1 to V4 of the methods on seed 999
    task2c_smoke.txt              smoke test on seed 999: times, controls and computation
                                  estimate that decided the budget rule
    task2c_explore_all.txt        exploratory E1 (ideal and perturbed generator of 1x) and E2 (oracle
                                  corrections in the windows of 2A)
    task2c_explore_e3.txt         exploratory E3 (ideal generator with the sources on the nodes)
    PREREGISTRO_1j_1k.txt         copy of the preregistration with the section Task 2C and its result

The 360 per-network files (task2c_parts/p1_seed{s}_{regimen}_{celda}_net{j}.json, 60 from part (i), and
task2c_parts/p2_..., 300 from part (ii)) and the run logs (log2c.txt, log2c_smoke.txt,
log2c_explore*.txt) remain in the working folder of the scripts; they are not copied here. Each part file
records the OpenBLAS kernel with which it was computed.

Verdict (analyze_2c_out.txt)
    rule             thin                           strong
    2C-a             met (3 of 3 and 3 of 3)        met (3 of 3 and 3 of 3)
    2C-b             below                          below (ties in 1 of 12 cases)
    2C-c literal     met (3 of 3 and 3 of 3)        met (3 of 3 and 3 of 3)
    2C-c with Holm   met (3 of 3 and 3 of 3)        met (3 of 3 and 2 of 3)
No case came within 0.005 of the threshold of 0.05.

2C-a, MT-FW+SC minus the best of BOMP+SC and MTSTACK+SC at cap 0.25 (the best is BOMP+SC in all cells):
    cell (seed)           thin      strong
    4 st. (2024)          +0.287    +0.130
    6 st. (2024)          +0.333    +0.324
    unilateral (2024)     +0.268    +0.205
    4 st. (7331)          +0.231    +0.077
    6 st. (7331)          +0.441    +0.202
    unilateral (7331)     +0.519    +0.453
Re-applied in each replicate of the bootstrap by windows it is met in 100 % (thin) and 99 % (strong) of
the replicates. In strong 7331 4 st. the intervals by windows and by network include 0.

2C-b, GL-FW+SC minus MT-FW+SC, at <= 0.1 / <= 0.25:
    cell (seed)           thin               strong
    4 st. (2024)          -0.264 / -0.310    -0.273 / -0.208
    6 st. (2024)          -0.375 / -0.431    -0.284 / -0.351
    unilateral (2024)     -0.479 / -0.451    -0.229 / -0.265
    4 st. (7331)          -0.244 / -0.295    -0.013 / -0.064
    6 st. (7331)          -0.603 / -0.559    -0.333 / -0.310
    unilateral (7331)     -0.556 / -0.605    -0.547 / -0.573
Below in 23 of 24 cases; the remaining one (strong 7331 4 st. at 0.1) counts as "ties". GL-FW+SC minus
BOMP+SC ranges from -0.18 to +0.01. GL-FW: maximum relative gap 1.4e-8 in the 1200 windows, certificate outside
the support / lambda <= 1.000000000, deviation on the support <= 8.7e-10, 0 re-solves.

2C-c, 30 networks per cell, MT-FW+SC minus QM-3C+SC at cap 0.25 with hierarchical CI 95:
    cell (seed)           thin                        strong
    4 st. (2024)          +0.347 [+0.248, +0.453]     +0.175 [+0.104, +0.254]
    6 st. (2024)          +0.323 [+0.228, +0.416]     +0.247 [+0.150, +0.337]
    unilateral (2024)     +0.283 [+0.189, +0.376]     +0.236 [+0.135, +0.334]
    4 st. (7331)          +0.306 [+0.194, +0.401]     +0.099 [+0.029, +0.181]
    6 st. (7331)          +0.347 [+0.242, +0.439]     +0.231 [+0.125, +0.339]
    unilateral (7331)     +0.224 [+0.119, +0.328]     +0.295 [+0.212, +0.383]
Holm: in thin all 36 comparisons reject; in strong 28 of 36. The 8 that do not reject are from 4 st.: 2024
against BOMP+SC at both caps; 7331 against QM-3C+SC (at 0.25, Holm p 0.044), against BOMP+SC and against
MTSTACK+SC, at both caps.

Exploratory diagnostic (after seeing part (i), not preregistered; seed 999, 60 windows per
cell). With the ideal generator of 1x, at 0.25 and 4 st., MT-FW 0.712, BOMP 0.458, MTSTACK 0.169, GL-FW 0.237.
With the sources moved to the nodes and the origins to the grid (E3; the truth on the nodes is here a
deliberate diagnostic, not an evaluation): MT-FW 0.763, BOMP 0.966, MTSTACK 0.542, GL-FW 0.729. With
oracle corrections in the windows of 2A (E2), BOMP improves little. Reading: the exact-time methods
collapse because of the coarse grid and the residual time errors, which the +-0.25 s peak search of
MT-FW absorbs; the implementations of BOMP and GL-FW are sound.

Reading. 2C-a and 2C-c are met in both regimes. The claim they support is robustness to grid
and time errors on a coarse grid, not superiority over sparse moment-tensor inversion
as such: a block inversion on a finer grid or with time tolerance has not been tested.
2C-b: the convex problem of the thesis on the grid, solved with certified optimality, remains far
below MT-FW; it is the framework and the certificate, not the best detector (E33 in artifact_manifest_and_error_log.md),
and the off-grid sliding Frank-Wolfe becomes the key question of the mathematical line.

Controls. Part (i): stations, delays, catalogue and corrections identical bit for bit to 2A (60/60 networks);
truth, origins and number of sources identical (1200/1200 windows); MT-FW+SC recomputed identical to what was
stored (120/120 windows). Part (ii): certified MT-FW+SC equal to the dense one (600/600 windows). Smoke
test: the pipeline of part (ii) reproduces the records of MT-FW+SC and QM-3C+SC of 2A (6/6 cells).
Validation (task2c_validation.txt): dictionary identities <= 8.7e-16; exact recovery of a noise-free
source in the model (12/12); KKT and gap of GL-FW in 24 windows; MTSTACK equal to the projection energy of
WAVE-MT (<= 1.5e-13). No control failed; task2c_CONTROL_FAILED.txt was not written.

Inherited code (E32). task2c_methods.py and task2c_run.py came from an earlier session and could not be run
as they were: task2c_validate.py, which the runner imported, was missing; the OpenBLAS kernel check
only looked at the environment variable, and the default numpy kernel on this machine changed from Prescott (that of
2A and 2B) to SkylakeX; and the tolerance of the restricted GL-FW solver (relative KKT 1e-11) was unattainable
in floating point. This was corrected before any run of 2024 or 7331: validation written, kernel fixed
and checked per network (threadpoolctl), tolerance 1e-9 with rounding slack in Armijo.

Declared deviations (task2c_report.txt, section 10): the GL-FW path stops at the first
lambda that completes the record of 5 events (record identical to that of the full path), so the
gap is reported at the solved lambdas (2 to 14 per window), not at the 30; the whitened group norm
per block and the reading of the BOMP score are our choices where the preregistration did not fix them;
tolerances 1e-6 (admission) and 1e-9 (restricted); kernel fixed explicitly (with Prescott, the OpenBLAS
of scipy chooses Katmai where in 2A it chose SkylakeX; scipy.linalg is not used in MT-FW, the corrections or the
methods of 2C); the Holm family and level are our reading; the point estimate of 2C-c pools windows
(the networks of 2A weigh 20/12) while the hierarchical bootstrap weighs the networks equally; BOMP+SC and MTSTACK+SC
of networks 0 to 4 in 2C-c are the records of part (i).

Reproduction (from the scripts folder, with task2a_parts/, the evaluation JSON files of 2A and
qm_recal_2a.json in the same folder):
  OPENBLAS_CORETYPE=Prescott python task2c_run.py validate          (writes task2c_validation.txt)
  OPENBLAS_CORETYPE=Prescott python task2c_run.py smoke --nwin 3    (writes task2c_smoke.txt)
  OPENBLAS_CORETYPE=SkylakeX python task2c_run.py run --part 1 --workers 2
  OPENBLAS_CORETYPE=Prescott python task2c_run.py run --part 1 --workers 2
  OPENBLAS_CORETYPE=Prescott python task2c_run.py run --part 2 --workers 2
  python analyze_2c.py                                              (writes analyze_2c_out.txt and .json)
  OPENBLAS_CORETYPE=Prescott python task2c_explore.py all           (exploratory; E1 and E2)
  OPENBLAS_CORETYPE=Prescott python task2c_explore.py e3            (exploratory; E3)
Each run process only launches the networks whose required kernel is its own (SkylakeX for the 18 networks of
7331 that 2A computed with it, Prescott for the rest), checks the loaded kernel and records it. run resumes
from the parts already written in task2c_parts/; if a network does not pass the exact control it writes
task2c_CONTROL_FAILED.txt and stops, and does not resume while that file exists. On another CPU the available
kernels may be different: the exact control of part (i) only passes with the kernel that reproduces bit for
bit each network of 2A.
Times: part (i) 318 s of CPU (60 networks), part (ii) 5342 s (300 networks, 17.8 s on average); 1.57 h of
CPU in total, 47 minutes of wall-clock time on 2 cores.
Dependencies: those of Tasks 2A and 2B (results/task2a/README.txt, results/task2b/README.txt), because
task2c_run.py imports task2b_ablation.py and the latter imports task2a_confirm.py; QuakeMigrate for QM-3C+SC in part
(ii); threadpoolctl to check the kernel; numba is optional and only accelerates the GL-FW solver.
