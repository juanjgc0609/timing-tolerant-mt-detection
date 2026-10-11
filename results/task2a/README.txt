Results of Task 2A: confirmatory benchmark with mismatched generator, real picker and
real comparators (Appendix Y.18 of the proposal; tasks_and_real_data_plan.md, section Task 2A).

Decision rules: PREREGISTRO_1j_1k.txt, section "Task 2A" (25 September 2026, written before
any evaluation run). The same copy contains the paragraph "Result of 2A" and the section
"Task 2B" (ablation of the radiation pattern) with its paragraph "Result of 2B"; the results of 2B are
in results/task2b/.

Design in brief: generator gen_realistic.generate_real (gradient medium with curved rays, Ricker
wavelet of variable central frequency and phase, scattered coda, site factors) with delays
N(0, 0.3 s) per station and phase fixed per network; all detectors use the same homogeneous model.
Cells ladder_4 (4 stations), ladder_6 (6 stations) and unilateral_c40 (unilateral of 8); regimes
thin (clutter 0.35) and strong (1.75); per seed, cell and regime, 5 networks with a catalogue of 20 events
(errors of 1 km and 0.1 s) and 20 evaluation windows. Seeds 2024 and 7331 for evaluation; 999 for
calibrating the picker, QuakeMigrate and the operating points of 2A-e.

Outputs
    analyze_2a_out.txt            tables 1 to 5 per regime, rules 2A-a to 2A-e (2A-e with its two readings)
                                  and verdicts; primary reading of the task
    task2a_report.txt             independent verification (in English): design as implemented,
                                  deviations D1 to D8, recalibration of QuakeMigrate, leakage checks,
                                  hit rule, sign of the corrections, starvation of GaMMA and PyOcto,
                                  re-derivation of figures, limitations and verdict
    task2a_seed{2024,7331}_{thin,strong}.json
                                  evaluation: per cell, records of MT-FW, MT-FW+SC, MT-FW-NODE+SC,
                                  WAVE-PT@20 and QM/QM-3C/QM-3C+SC, events per window and configuration of
                                  GaMMA and PyOcto; per network, stations, delays, corrections, catalogue and
                                  identity checks; per window, true sources; summary of
                                  the STA/LTA picks and degrade_y
    task2a_seed999_{thin,strong}.json
                                  calibration of 2A-e (networks different from the evaluation ones)
    task2a_frozen_2Ae.json        thresholds of MT-FW+SC and QM-3C+SC and configuration index of
                                  GaMMA-STALTA+SC frozen per regime, cell and budget (0.1, 0.25);
                                  each entry is [threshold or configuration, calibration recall, false events
                                  per calibration window]; written before the first evaluation file
    qm_recal_2a.json              QuakeMigrate configurations (QM and QM-3C) recalibrated on seed 999
    task2a_smoke.txt              smoke test: times per component and estimate of the computation
                                  budget (1.47 h on 2 cores, under the limit of 10 h: N_WIN = 20)
    log2a_qmcal.txt               recalibration of QuakeMigrate (detect and trigger grids)
    log2a_seed999.txt             calibration run; log2a_seed999_attempt1.txt: first attempt,
                                  interrupted by a PyOcto failure with an empty pick list (the
                                  four parts it managed to write were reused; see D5 of the report)
    log2a_freeze.txt              frozen operating points, as printed
    log2a_seed2024.txt            evaluation of seed 2024
    log2a_seed7331_part1.txt, log2a_seed7331.txt
                                  evaluation of seed 7331, interrupted and resumed (each job
                                  has its own seed, so the resumption is exact)
    PREREGISTRO_1j_1k.txt         copy of the preregistration
    verif_2a/                     scripts and outputs of the independent verifier, all with suffix _verif:
        rederive_verif.py              FROC curve re-derived without exact_HF
        rederive_2Ae_verif.py          freeze and application of 2A-e re-derived
        starve_sc_verif.py, analyze_starve_verif.py, starve_sc_2024_thin_verif.{json,log}
                                       80 regenerated windows: sign of the corrections with oracle
                                       delays and starvation of GaMMA and PyOcto by the picker
        qm_trigwin_verif.py, qm_trigwin_2024_thin_verif.{json,log}
                                       sensitivity of QuakeMigrate to the trigger window
        cluster_boot_verif.py          per-network bootstrap (exploratory)
        sc_strong_verif.py, sc_strong_7331_verif.{json,log}
                                       corrections with strong clutter (exploratory)
        subset_k01_verif.py            windows with 0 or 1 source (exploratory)
        rerun_part_verif.py, seed999_strong_unilateral_c40_net1_verif.json,
        rerun_999_strong_uc40_net1_verif.log
                                       repetition of one part of seed 999 with the final code
        selfcheck_verif.txt            output of the self-test
        dbg_consistency_verif.py       check of the window regeneration

The 90 per-network files (task2a_parts/seed{s}_{regimen}_{celda}_net{j}.json) remain next to the scripts;
each aggregated JSON is their concatenation.

Verdicts (analyze_2a_out.txt)
    rule    thin                          strong
    2A-a    met                           met
    2A-b    met                           met
    2A-c    met                           met
    2A-d    met                           not met
    2A-e    per cell yes; strict no       per cell yes; strict no
2A-a: MT-FW+SC minus the best of GaMMA, PyOcto and QuakeMigrate (seven variants) from +0.09 to +0.49 at
<= 0.1 false events per window and from +0.17 to +0.46 at <= 0.25, in the 3 cells of both seeds and regimes;
the best comparator is a QuakeMigrate variant in all 24 cases. 2A-c: with weak clutter, GaMMA-DY+SC
ties in the unilateral cell of 8 (+0.01 on both seeds). 2A-d fails with strong clutter because the
catalogue drawn with that clutter contaminates the corrections. 2A-e: the realised rate exceeds 2 times the
budget only in unilateral of 8, seed 7331, budget 0.1 (0.22 thin, 0.23 strong, limit 0.20);
both readings are reported. The window-based intervals ignore the clustering by network; with a per-network
bootstrap four cells include 0 at cap 0.1 (task2a_report.txt, section 4.7).

Reproduction (from the scripts folder, with gamma_src/ in the same folder):
  python task2a_confirm.py selfcheck
  python task2a_confirm.py qmcal --workers 2                 (writes qm_recal_2a.json)
  python task2a_confirm.py smoke --nwin 3 --workers 2        (writes task2a_smoke.txt)
  python task2a_confirm.py run --seed 999 --workers 2
  python analyze_2a.py --freeze --cal 999                    (writes task2a_frozen_2Ae.json)
  python task2a_confirm.py run --seed 2024 --workers 2
  python task2a_confirm.py run --seed 7331 --workers 2
  python analyze_2a.py --seeds 2024,7331 --cal 999          (writes analyze_2a_out.txt and .json)
run accepts --regime thin|strong|both (default both) and --nwin (default 20); it resumes from the
parts already written in task2a_parts/. The QuakeMigrate temporary files go to $TASK2A_SCRATCH (by default
a directory task2a_scratch in the system temporary folder). analyze_2a.py re-derives the frozen points
from seed 999 and aborts if they do not match task2a_frozen_2Ae.json.
Times: seed 2024, 46 minutes of wall-clock time on 2 cores; seed 7331, about 53 minutes of computation in three
sessions; seed 999, around 1 hour including the first attempt.
Dependencies besides numpy and scipy: ObsPy (real_picker.py), QuakeMigrate 1.2.2, GaMMA from gamma_src/
(np.random.seed(42) before each call) and PyOcto 0.2.0.
