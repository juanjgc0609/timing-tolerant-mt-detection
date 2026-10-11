Results of Task 2B: ablation of the radiation pattern within the same solver (Appendix Y.19 of the
proposal; tasks_and_real_data_plan.md, section Task 2B).

Decision rules: PREREGISTRO_1j_1k.txt, section "Task 2B" (26 September 2026, 14:30 UTC, written
before computing it and not modified). The copy in results/task2a/PREREGISTRO_1j_1k.txt contains that section
and the paragraph "Result of 2B".

Design in brief: the same windows, networks, catalogues and station corrections as Task 2A (seeds
2024 and 7331, regimes thin and strong, cells ladder_4, ladder_6 and unilateral_c40, 5 networks per cell and 20
windows per network: 1200 windows), regenerated deterministically with the functions of task2a_confirm.py,
which task2b_ablation.py imports and does not copy. MT-FW+SC is the stored record of 2A. MT-FW-FREE+SC is the same
call to certified_fw_record (675 nodes, kmax 5, min_sep 5 km, CAPV 6, the same corrections) with
G = Q = identity of dimension 6 ns at each node: free three-component amplitude per station and phase,
statistic ||v||^2/2 (energy of the peak vectors, invariant to sign) and deflation that subtracts the
observed vectors. Rule 2B, per regime: MT-FW+SC - MT-FW-FREE+SC >= 0.05 at both caps (0.1 and 0.25
false events per window) in at least 2 of 3 cells, on both seeds, on point estimates. Endpoint and
intervals as in 2A: exact FROC curve, paired bootstrap by windows (2000 replicates, the same
resamples as 2A) and, as an additional reading, per-network bootstrap (2000 replicates).

Outputs
    analyze_2b_out.txt            tables 1 to 4 per regime (recall; rule 2B with intervals by windows and
                                  by network; FREE minus the best comparator of 2A; FREE minus WAVE-PT@20 and
                                  2A-b as reference), rule, fraction of replicates in which it is met and
                                  verdict; primary reading of the task
    analyze_2b_out.json           the same figures in machine format
    task2b_report.txt             report (in English): verdict, design as implemented, controls,
                                  BLAS kernel incident, key tables, interpretation and limitations
    task2b_seed{2024,7331}_{thin,strong}.json
                                  per cell, records of MT-FW+SC (from 2A) and of MT-FW-FREE+SC per window,
                                  events per window, network of each window, truth, origin times and
                                  result of the controls of each network
    task2b_selfcheck.txt          self-test: shapes, score equal to ||v||^2/2, sign invariance,
                                  exact deflation and per-station sign change of the waveform
    log2b.txt                     final run (18 networks with OPENBLAS_CORETYPE=SkylakeX and 42 with the default
                                  kernel)
    log2b_attempt1.txt            first attempt, stopped by the exact control after 32 networks
    task2b_control_failed_attempt1.txt
                                  the discrepancies that stopped the first attempt
    task2b_kernel_check.txt, log2b_kernel.txt
                                  sensitivity to the BLAS kernel (exploratory): the 18 networks recomputed with
                                  the default kernel

The 60 per-network files (task2b_parts/seed{s}_{regimen}_{celda}_net{j}.json) remain next to the scripts;
each aggregated JSON is their concatenation. The 28 written in the final run record the OpenBLAS kernel
(fields blas_arch and openblas_coretype); the 32 from the first attempt do not carry the field and ran with the
default kernel. No file of Task 2A was modified.

Verdict (analyze_2b_out.txt)
    rule    thin                                 strong
    2B      met (2 of 3 on each seed)            met (3 of 3 on each seed)
MT-FW+SC minus MT-FW-FREE+SC, at <= 0.1 / <= 0.25 false events per window:
    cell (seed)           thin             strong
    4 st. (2024)          +0.17 / +0.28    +0.25 / +0.26
    6 st. (2024)          +0.18 / +0.11    +0.32 / +0.35
    unilateral (2024)     -0.01 / -0.01    +0.17 / +0.22
    4 st. (7331)          +0.15 / +0.24    +0.06 / +0.12
    6 st. (7331)          +0.26 / +0.12    +0.26 / +0.17
    unilateral (7331)     +0.04 / +0.0494  +0.33 / +0.33
In thin the rule passes at the minimum: in the unilateral cell of 8 MT-FW-FREE+SC equals MT-FW+SC on both
seeds, and in those of 6 stations several intervals include 0. Re-applied in each bootstrap replicate
(exploratory, independent resamples per cell) it is met in 0.59 (windows) and 0.62 (networks)
of the replicates in thin, and in 0.98 and 0.97 in strong. Exploratory: in thin MT-FW-FREE+SC beats the best
comparator of 2A (QuakeMigrate) at 0.25 in 5 of 6 cell-seeds, by +0.10 to +0.34; in strong it stays at its
level (-0.08 to +0.08) and below WAVE-PT@20 (-0.01 to -0.23). The identity removes at once the sign of the
pattern, the polarisation constraint and the amplitude geometry: 2B shows that the operator matters,
it does not isolate the sign.

Controls: FREE equal to ||v||^2/2 in 675 x 31 hypotheses (relative 6.8e-16); FREE records identical when
changing the sign of any of the 6 stations (6/6) and MT-FW+SC records different (6/6); stations,
delays, catalogue and corrections identical bit for bit to 2A (60/60 networks); truth and origin times
identical (1200/1200 windows); MT-FW+SC recomputed identical to what was stored (120/120 windows); certified
FREE identical to dense FREE (120/120).

BLAS kernel incident (E31 in artifact_manifest_and_error_log.md). The OpenBLAS of numpy chooses its
kernel when the process starts; on this machine it chooses Prescott by default, but the 18 networks that 2A wrote
in its resumed sessions of seed 7331 (ladder_4 in both regimes, ladder_6 thin networks 0 to 4 and
ladder_6 strong networks 2 to 4) were computed with SkylakeX. With OPENBLAS_CORETYPE=SkylakeX they are reproduced bit
for bit; with Prescott, Haswell, Sandybridge or Zen they are not. The final run fixed the kernel per network and
recorded it in each part file; the exact control was not relaxed. The kernel does not change any recall in the 4
affected cells.

Reproduction (from the scripts folder, with task2a_parts/ and the evaluation JSON files of 2A in the same
folder):
  python task2b_ablation.py selfcheck                        (writes task2b_selfcheck.txt)
  OPENBLAS_CORETYPE=SkylakeX python task2b_ablation.py run --workers 2 --no-aggregate --nets \
    7331_strong_ladder_6_2,7331_strong_ladder_6_3,7331_strong_ladder_6_4,7331_thin_ladder_6_0,\
    7331_thin_ladder_6_1,7331_thin_ladder_6_2,7331_thin_ladder_6_3,7331_thin_ladder_6_4,\
    7331_strong_ladder_4_0,7331_strong_ladder_4_1,7331_strong_ladder_4_2,7331_strong_ladder_4_3,\
    7331_strong_ladder_4_4,7331_thin_ladder_4_0,7331_thin_ladder_4_1,7331_thin_ladder_4_2,\
    7331_thin_ladder_4_3,7331_thin_ladder_4_4
  python task2b_ablation.py run --workers 2                  (remaining networks and aggregation of the JSON files)
  python analyze_2b.py                                       (writes analyze_2b_out.txt and .json)
  python task2b_kernel_check.py --workers 2                  (exploratory; writes task2b_kernel_check.txt)
The list of --nets goes in a single argument, without spaces. run resumes from the parts already written in
task2b_parts/; if a network does not pass the exact control against task2a_parts/, it writes
task2b_CONTROL_FAILED.txt and stops, and does not resume while that file exists. The kernel that OpenBLAS
chooses by default depends on the machine: on another CPU it may be necessary to set OPENBLAS_CORETYPE also
for the other 42 networks (to the kernel that reproduces them bit for bit) for the control to pass.
Times: 60 network jobs, 8 s each on average (495 s of CPU in total) on 2 cores.
Dependencies: those of Task 2A (results/task2a/README.txt), because task2b_ablation.py imports
task2a_confirm.py and analyze_2b.py imports analyze_2a.py; all the scripts of 2A in the same folder.
