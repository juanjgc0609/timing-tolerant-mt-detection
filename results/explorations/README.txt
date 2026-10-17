Exploration on the calibration seed (explore_contrib): where the margin of the moment-tensor detectors does not lie,
where the hypothesis of Task 2F came from, the exploration C8 (site terms versus bounded influence) from which that of
Task 2I came, the explorations C9 (evidence without travel times) and C10 (soft time priors and time budget), after
2I, the exploration C11 from which Task 2L came, and the explorations C12 (double-couple statistic), C13 (decoy
calibration) and C14 (self-interference), after 2L. Appendices Y.23, Y.25, Y.27, Y.28 and Y.29 of the supplementary
volume; tasks_and_real_data_plan.md, sections "EXPLORACION explore_contrib", "EXPLORACION C8", "EXPLORACIONES C9 y C10",
"TASK 2L" and "EXPLORATIONS C12, C13 and C14".

STATUS: EXPLORATORY, NOT CONFIRMATORY. Nothing in this folder tests a registered hypothesis. Several numbers were seen
while the designs were being chosen; they are pilot estimates and must not be cited as confirmatory evidence. The only
candidate from C1 to C7 that was carried to a registered test (the per-station cap in dense windows) was chosen among
several, and its pilot estimate is subject to the winner's curse; the confirmatory test is Task 2F
(results/task2f/), on seeds that this exploration never touched. C8 produced a second candidate, HUBER-N, the best
of 11, with the same selection bias; its confirmatory test is Task 2I (PREREGISTRO_1j_1k.txt, section
"Task 2I"; seeds 5151 and 6262), preregistered, executed and verified, whose result is in results/task2i/.
C9 and C10 carried no candidate to a registered test: their two preregistration drafts (SELPOL in
REPORT_C9.txt, section 8; T3 refinement in REPORT_C10.txt, section 11) name the same unused seeds, 7171
and 9393; neither was registered or run, and only one can use them.

Data. Seed 999 only. Networks task2a_confirm.network(999, cell, j) with new j:
    j = 0-4     calibration networks of 2A/2E: not used for new computation (only the parts of check 4 of
                verif_2e were read, to validate the ceiling code)
    j = 10-15   C1 (TAG 3001), C2 (TAG 3002), C3 (TAG 3003), masking pilot, ablation of C2
    j = 20-29   C6 (TAG 3006) and C6b (TAG 3008)
    j = 30-35   masking replication (TAG 3007)
    j = 40-59   C7 (TAG 3701 STD, 3702 DENSE, 3703 GABOR90, 3704 BERLAGE90; 3705 validation; 3706 resolution)
    j = 60-79   C8, main designs (TAG 3801 DENSE-IID, 3802 STD-IID, 3803 DENSE-PERS, 3804 STD-PERS with j =
                60-69; 3809 persistent factors; 3810 validation on network 60)
    j = 80-89   C8, sensitivity to the amplitude spread (TAG 3805 sigma 0.25, 3806 sigma 0.8)
    j = 99      C8: built once to measure a setup time; no windows or results
    j = 130-149 C9 (TAG 3901 isolated events, 3903 sources of the Fisher study, 3905 Monte Carlo of the bound, 3911
                STD, 3912 DENSE)
    j = 199     C9: validation, amplitude calibration and setup times (TAG 3902, 3904, 3913)
    j = 160-179 C10 (TAG 4001 STD, 4002 DENSE; level sweep, first-iteration diagnostic with j = 160-169)
    j = 195     C10: validation, times and choice of the time-level rule (TAG 4003, 4004)
    j = 200-239 C11 (development 200-209, evaluation 210-239; TAG 4101-4105)
    j = 1200-1229, 1290  C12 (evaluation 1200-1215, nulls 1220-1229, validation 1290; TAG 4301-4305)
    j = 1400-1411 C13 (TAG 1401-1406)
    j = 1600-1602, 1610-1621  C14 (development and evaluation; TAG 1651, 1652)
Window streams [999, TAG, cell, j, regime, 41|233|61], all new when used. Seeds 2024, 7331, 4242 and
8888 were never generated or read; C9 and C10 likewise neither generated nor read 5151, 6262, 7171 or 9393 (the only
contact of C10 with 5151/6262 was running analyze_2i.py, frozen, at the end of the 2I work). Common protocol (that of
2E): grid of 675 nodes, origins 5.0-8.0 s every 0.1 s, peak search of +-0.25 s, exclusion of 5 km, hit = epicentre
within 5 km, kmax 5 (7 in dense windows), exact FROC, recall at <= 0.1 and <= 0.25 realised false events per window,
hierarchical bootstrap (cells as strata, networks and windows, 1000-2000 replicates, operating point re-chosen per
replicate); paired differences on the same windows.
Computation: about 2.7 h of CPU (REPORT.txt) and 2.6 h (REPORT_C7.txt) on 2 cores; C8, 2 h 42 min wall-clock on 2
cores (5.3 h of CPU, REPORT_C8.txt); C9, about 2 h 15 min wall-clock with one process; C10, 2 h 22 min wall-clock.

Reports
    REPORT.txt        ranking of candidates, C1, C2 and its ablation, C3, masking, C6 and C6b, judgements, drafts
                      of the preregistration of 2F (section 7.1) and of a Task 2G not executed (7.2), list of files
    REPORT_C7.txt     C7: subspaces invariant to the source pulse, validation, standard and dense windows,
                      pulses outside the family, false alarms, times, judgement and draft of a Task 2H not
                      executed
    REPORT_C8.txt     C8: per-station amplitude (site) terms and robust fits within the moment-tensor
                      certificate; generator fact, design, what was fixed before and after seeing results,
                      validation, dense and standard windows, sensitivity to the spread, sweep of the Huber
                      threshold, times, judgement and draft of Task 2I (section 10)
    REPORT_C9.txt     C9: polarisation, amplitudes and S/P geometry as location evidence without travel times;
                      design, what was fixed before and after seeing results, validation, locators on isolated
                      events, fusion with BOMP-PK and MT-FW (reordering, relocation, extra term, oracle
                      worlds), per-channel Cramer-Rao bound, judgement and conditional draft of a Task 2J
    REPORT_C10.txt    C10: Gaussian time prior within BOMP-PK and MT-FW; design, what was fixed before and after,
                      validation, dense and standard windows with two hit rules, level sweep and
                      oracle decomposition, time budget and separability in the first iteration,
                      SC-R corrections, judgement, novelty and conditional draft of a Task 2K
    LIT_NO_TRAVELTIME.md  literature sweep 2000-2026 on detection and location without travel times, eight
                      families; it is not a systematic review: 26 of its 72 references were read only at the
                      metadata level and several publisher pages could not be opened

Candidates and result (cap 0.1 / 0.25, CI 95)
    C1   grids of 3, 6 and 9 km (thin, 180 windows, 145 events). BOMP-PK+SC at 9 km - at 3 km: -0.324 [-0.524,
         -0.052] / -0.379 [-0.512, -0.174]. TS-SFW on the same coarse grid recovers +0.03 to +0.10 and never
         reaches BOMP-PK at 3 km. NULL.
    C2   dense windows, 3 to 5 events at < 9 km, kmax 7 (360 windows, 1439 events, both regimes).
         MT-FW+SC - BOMP-PK+SC: +0.074 [+0.016, +0.136] / +0.070 [+0.020, +0.128]. Ablation: MT-FW without cap and
         BOMP-PK with sequential deflation give identical hit sequences (360/360); MT-FW - MT-FW-NOCAP
         +0.081 / +0.071 (thin); the cap explains the gain; the P/S cross term and the joint refit do not
         matter. Replication on new networks (C7 DENSE, 600 windows): +0.077 [+0.030, +0.122] / +0.068 [+0.029,
         +0.109]. CANDIDATE CARRIED TO TASK 2F.
    C3   large time errors (thin, 180 windows): statics of 0.3 to 0.8 s with and without corrections, medium
         7.5 % slower; wider peak search, shared shifts, online station terms.
         Only with 0.8 s and without corrections does BOMP-PK-12 win +0.156 [+0.015, +0.284] / +0.148 [+0.000,
         +0.315]; with catalogue corrections nothing beats BOMP-PK+SC. NULL.
    C6   1D velocity models within BOMP-PK and MT-FW (times and radiation operator): generator medium
         (oracle), homogeneous, gradient and two layers calibrated by JHD with the 20 catalogue events (600
         windows, 475 events). Oracle - homogeneous (BOMP-PK): +0.006 [-0.066, +0.065] / +0.011 [-0.059, +0.064];
         calibrated gradient with strong clutter -0.104 [-0.202, +0.013] at 0.25; only the oracle with true statics
         wins (+0.069 [+0.004, +0.141] at 0.25), unattainable. C6b (attribution, thin): realistic margin
         of BOMP-PK at 0.25 +0.129 [-0.037, +0.228] by default, +0.194 [+0.081, +0.312] with exact medium,
         +0.177 with fitted pulse, +0.153 without site and polarisation perturbations, +0.242 without coda: no
         single mismatch by itself closes the margin. NULL.
    C7   subspace invariant to the source pulse: PH (2 Hz Ricker and its Hilbert transform, rank one) and
         PHF (basis of 4 capturing >= 0.95 of the pulse family). Standard windows (600, 465 events): PH -
         BOMP-PK -0.114 [-0.203, -0.017] / -0.084 [-0.148, -0.023]; PHF -0.140 / -0.080 [-0.142, -0.005]; dense:
         PH -0.021, PHF -0.046 at 0.25; Gabor and Berlage pulses outside the family: PH -0.004 and +0.026 at 0.25
         (CI including 0). Without noise and with all the realism of the generator the true block captures 0.497
         of the energy and the invariance adds 0.014-0.050. NEGATIVE. By-product: without the P/S cross term, the
         block search reproduces BOMP-PK in 2400/2400 windows (scores <= 1e-12) at 10.7 times lower cost.
    Arrival masking (each arrival explains at most one event): pilot chosen among 12 variants; replication
         with the fixed variant (networks 30-35): thin +0.149 [-0.007, +0.265] / +0.084 [+0.000, +0.163]; strong
         -0.007 / +0.014. DOUBTFUL, weak clutter only.
    Speed: in standard windows MT-FW+SC - BOMP-PK+SC +0.017 [-0.061, +0.092] / +0.002 [-0.050, +0.055] with
         6.3 times lower cost. ENGINEERING, not a method.
    C4 (false alarms in windows with events) and C5 (location accuracy and uncertainty): scripts written
         and never run. c5_run.py uses TAG 3006, which C6 already used: change it before running it.
    C8   site terms versus bounded influence (cap 0.25, CI 95). Generator fact: gen_realistic
         draws the lognormal amplitude factors (sigma 0.5, per station and phase) in each window and in each
         catalogue event; they are not station terms. Two worlds: IID (generator unchanged) and PERS (factors
         fixed per network, shared by catalogue and windows). Dense windows, sigma 0.5 (IID 600 windows, 2410
         events; PERS 480 windows, 1919 events):
           MT-FW (cap 6) - BOMP-PK         IID +0.042 [+0.004, +0.083]   PERS +0.061 [+0.008, +0.098]
           SITE-CAT - MT-FW                IID -0.066 [-0.108, -0.028]   PERS -0.040 [-0.078, +0.007]
           SITE-CAT-CAP - MT-FW            IID -0.002 [-0.026, +0.022]   PERS +0.006 [-0.015, +0.043]
           ORACLE (no cap) - MT-FW         IID -0.011 [-0.056, +0.044]   PERS -0.008 [-0.049, +0.055]
           ORACLE-CAP - MT-FW              IID +0.045 [+0.009, +0.074]   PERS +0.032 [+0.003, +0.073]
           HUBER-N - MT-FW                 IID +0.022 [-0.014, +0.055]   PERS +0.005 [-0.025, +0.048]
           HUBER-N - BOMP-PK               IID +0.064 [+0.028, +0.105]   PERS +0.066 [+0.024, +0.105]
         Neither the catalogue terms nor even the true factors replace the cap. With sigma 0.25 the
         oracle gains nothing over BOMP-PK (-0.001) and the cap gains +0.043 [-0.026, +0.099]: the cap is not only a
         remedy for the amplitude factors. HUBER-N (Huber fit on the norm of the per-block residual, threshold
         fixed by rule at 2.7955 noise deviations, square root of the 0.95 quantile of chi-squared(3)) matches
         cap 6 with sigma <= 0.5, is worse with sigma 0.8 (-0.043 [-0.113, -0.001]) and costs 9 times more; in
         IID standard windows, HUBER-N - BOMP-PK -0.014 [-0.050, +0.042]. Post-hoc sweep of the threshold: plateau
         from 1.5 to 3.4, half of the gain at 5, none at 8. HUBER-S, GAIN, GAIN-FREE and NORM-TR remain
         below MT-FW. Fraction of the margin up to CEIL-GD: MT-FW 0.14 / 0.15, HUBER-N 0.21 / 0.16, ORACLE-CAP
         0.29 / 0.23 (IID / PERS): nothing exceeds one quarter. HUBER-N is the best of 11 candidates (selection
         bias); no multiplicity correction. SITE TERMS: NEGATIVE. HUBER-N: CANDIDATE CARRIED TO
         TASK 2I.
    C9   evidence without travel times (600 windows per design; DENSE 2388 events, STD 444; cap 0.25, CI 95).
         Bayesian Cramer-Rao bound of the epicentre (4 | 6 | 8 stations, km): prior alone 12.25; times at 0.144 s
         0.68 | 0.43 | 0.69; P polarisation at 10 degrees 9.23 | 8.53 | 7.79; amplitudes 9.18 | 7.75 | 9.11;
         S/P ratio 11.77 | 10.46 | 10.88; all channels without times 5.59 | 4.35 | 4.81; joint with the
         times 0.65 | 0.41 | 0.66 (improvement of 4-5 %). The limit is the waveform noise: P of 1.5-1.8
         noise deviations, 32-39 degrees of angular noise. Locators on 600 isolated events, vectors
         read at the true arrival: 12.7-15.2 km (P direction), 8.2-9.8 km (P and S), 7.7-8.4 km (moment-tensor
         projection with free times), 12.5-19.4 km (1/r), versus 0.5-0.7 km with the times. Atom
         relocation: -0.11 to -0.18 in dense and -0.23 to -0.48 in standard, also with oracles (true factors
         -0.116 [-0.158, -0.081]; world without deviation or factors -0.138 [-0.187, -0.089]). Weight by
         polarisation admissibility (SELPOL, added after seeing results, best of about 50): MT-FW+SELPOL -
         MT-FW +0.046 [+0.007, +0.070] and - BOMP-PK +0.045 [+0.003, +0.086] in dense; +0.002 [-0.041, +0.064] over
         BOMP-PK in standard; reordering +0.028 at 10 degrees, +0.003 at 20 and -0.015 at 45. In this sample MT-FW -
         BOMP-PK in dense was -0.001 [-0.030, +0.041] at 0.25. LOCATION WITHOUT TIMES: NEGATIVE, information
         limit. POLARISATION WEIGHT: SMALL AND CONDITIONAL, not registered.
    C10  soft time priors (600 windows per design; DENSE 2432 events, STD 471; cap 0.25, CI 95; LOC =
         registered rule, LOC+T = additionally origin within 0.5 s). Time level by catalogue rule: median 0.138 s
         (thin) and 0.224 s (strong). Variants with Gaussian loading - MT-FW, LOC: -0.028 to -0.002 in dense and
         -0.066 to -0.006 in standard; hard window of +-0.5 s -0.047 [-0.087, -0.008] in dense. Time budget
         (s rms): model error 0.05 P and 0.07 S; error of the +SC corrections 0.12-0.17 (thin) and 0.18-0.24
         (strong); position within the 3 km cell 0.18 P and 0.30 S. Oracle - MT-FW, LOC: true statics +0.047
         [-0.003, +0.082] (dense) and +0.098 [+0.032, +0.141] (standard); with prior at the model-error level
         +0.081 [+0.035, +0.117] and +0.146 [+0.065, +0.217] (0.592 versus 0.603 of CEIL-GD in standard; 0.31 of the
         margin in dense). First iteration (235 isolated events): they win 0.485 with the hard rule, 0.489 with the
         catalogue prior, 0.609-0.634 with the oracle; the remaining 37-39 % is lost to amplitude. T3
         (refinement on a 1 km grid) in standard: BOMP +0.059 [-0.009, +0.112] (LOC) and +0.100 [+0.031, +0.168]
         (LOC+T); MT +0.011 and +0.066 [+0.019, +0.136]. SC-R (post hoc): correction error 0.145 -> 0.110 s (thin)
         and 0.214 -> 0.202 s (strong); with T3 in standard +0.064 [-0.004, +0.123] (LOC) and +0.123 [+0.058, +0.201]
         (LOC+T). TIME PRIOR WITH CATALOGUE CORRECTIONS: NEGATIVE. T3: SMALL AND CONDITIONAL, not
         registered.

Reading. The margin of 0.08 to 0.24 below the realistic ceiling measured by the verification of 2E is not in the grid,
nor in the time error absorbed by the peak search, nor in the velocity model, nor in the source pulse; with
the exact medium it persists. The only element with an effect is the per-station cap in dense sequences, whose
mechanism the verification of 2F confirmed and bounded and C8 made precise: bounded influence of each station-and-phase
block against the model's amplitude error, not correction of a site term nor separation of events (E37 and E39 in
artifact_manifest_and_error_log.md). No treatment of amplitudes, the true factors included, captures
more than one quarter of the margin. C9 and C10 split it into two parts: one of calibration (the accuracy of the
station terms; with the statics known at the level of the model error the detector reaches the ceiling on
isolated events) and one of amplitude (weak arrivals, clutter, coda, amplitude factors). Physics without travel
times is information-limited in this benchmark, and soft time priors do not help with catalogue-quality
corrections. The open problem is the self-consistent calibration of the station terms and the amplitude evidence,
not another location principle. All of this with one seed, 20 networks per cell and regime, and oracles and
analyses after the results declared in the reports.

Scripts
    xc.py                 common harness: networks, catalogue and corrections (reproduces task2b_ablation.setup_net with
                          cfg=None and statics of 0.3 s), new window streams, grid, recall
    analyze.py            recall with hierarchical bootstrap and paired differences, generic
    t_check.py            check that the harness reproduces setup_net
    c1_run.py, c1_analyze.py                         C1 -> out/c1_thin.json, c1_thin_analysis.txt
    c2_run.py, c2_analyze.py                         C2 -> out/c2_thin.json, c2_thin_analysis.txt
    c2_ablate.py, c2_ablate_analyze.py               ablation of C2 -> out/c2_ablate_*.json, *_analysis.txt
    c3_run.py, c3_analyze.py                         C3 -> out/c3_thin.json, c3_thin_analysis.txt
    c4_run.py, c4_analyze.py, c5_run.py, c5_analyze.py   never run
    leak_test.py, leak_analyze.py, leak2.py, leak2_analyze.py   masking -> out/leak_*, out/leak2_*
    c6_lib.py             models, operators, dense MT-FW, JHD; c6_validate.py -> out/c6_validation.txt
    c6_run.py, c6_ceil.py, c6_analyze.py, c6_percell.py, c6_ttresid.py   C6 -> out/c6_*.json, c6_analysis.txt,
                          c6_percell.txt, c6_ttresid.txt
    c6b_run.py, c6b_analyze.py                       C6b -> out/c6b_*.json, c6b_analysis.txt
    c7_lib.py, c7_validate.py, c7_resolution.py, c7_run.py, c7_analyze.py, c7_power.py   C7 -> out/c7_*
    power.py              power of the confirmatory designs from the pilot SD -> out/power_*.txt
    c8_lib.py             copy of the generator that returns the amplitude factors and allows fixing them, IID and
                          PERS worlds, catalogue amplitude terms, MT-FW engine with interchangeable statistic
    c8_sitedev.py         choice of the catalogue-term estimator (networks 60-63) -> out/c8_sitedev.txt
    c8_validate.py        validations V1 to V8 -> out/c8_validation.txt
    c8_run.py, c8_analyze.py   C8 -> out/c8_analysis_{DENSE-IID, DENSE-IID+H, DENSE-IID+W, DENSE-PERS, STD-IID,
                          STD-PERS, DENSE-IID-S25, DENSE-IID-S80}.txt (the JSON files and the run logs were not
                          copied to this folder)
    c9_lib.py             location statistics without travel times, vector readers, BOMP-PK and
                          MT-FW loops with relocation and selection hooks
    c9_standalone.py, c9_standalone_analyze.py   isolated events -> out/c9_standalone_analysis.txt; validation ->
                          out/c9_validation.txt
    c9_fisher.py          per-channel Cramer-Rao bound -> out/c9_fisher.txt
    c9_run.py, c9_extra.py, c9_analyze.py, c9_pairs.py   fusion with the detectors -> out/c9_analysis_{DENSE, STD,
                          DENSE-BOTH0, DENSE-POL0, DENSE-AMP0, STD-BOTH0, DENSE-POL20, DENSE-POL45}.txt, c9_pairs.txt,
                          c9_validation_end.txt
    c9_queue.sh, c9_queue2.sh, c9_queue3b.sh   order of the runs (c9_queue3.sh and c9_queue4.sh were
                          replaced by 3b)
    c10_lib.py            time-level rule, penalised peak choice, loops with time loading and
                          refinement, hit rules
    c10_run.py, c10_analyze.py   main runs -> out/c10_analysis_{DENSE,STD}.txt; validation ->
                          out/c10_validation.txt
    c10_sweep.py, c10_extra.py   fixed levels and true statics with catalogue level -> out/c10s_analysis_*.txt,
                          c10x_analysis_*.txt
    c10_diag.py, c10_diag_analyze.py   separability in the first iteration -> out/c10_diag.txt
    c10_screloc.py, c10_screloc_run.py   SC-R corrections -> out/c10_screloc_check.txt, c10r_analysis_*.txt
    c10_queue.sh to c10_queue4.sh   order of the runs
                          (from C9 and C10 the scripts and the text outputs were copied to this folder; the JSON files
                          and the run logs remained in the working folder)
c2_run.draw_dense, c2_ablate.PKDictCap, c6_lib.ModelOps and c6_lib.mtfw_dense were frozen with hashes for Task 2F
(results/task2f/task2f_code_hashes.txt) and were not modified afterwards.

Validations (all correct): c6_validation.txt (homogeneous model equal to the registered times and operator,
<= 6e-17; homogeneous ModelOps reproduces BOMP-PK to 3e-13 and the dense MT-FW reproduces certified_fw_record to
8e-14; the ceiling reproduces that of verif_2e to 7.5e-5); c7_validation.txt (q = 1 equal to the Gram matrix and
correlations of PKDict of 2D exactly; noise-free invariance as designed); c8_validation.txt (the copy of the generator
equal to the original bit for bit; the engine reproduces c6_lib.mtfw_dense with and without cap in 6 of 6 windows;
Huber with infinite threshold equal to least squares, 6e-14); c9_validation.txt and c9_validation_end.txt (operators
and rays equal to the registered ones, <= 3e-16; noise-free recovery; bomp_c9 and mtfw_c9 equal to xc.bomp_pk and
c8_lib.mtfw_general in 9 of 9 windows); c10_validation.txt (loops in hard mode equal to the registered detectors in
9 of 9 windows; model and statics exact). Not validated: the convergence of GAIN-FREE.

Reproduction (from this folder, with the scripts of the root of outputs/ in the folder above; OPENBLAS_CORETYPE
set as in 2F):
Most of the runners take regime, networks, windows per network and output file, for example:
  python c1_run.py thin 10,11,12,13,14,15 10 out/c1_thin.json
  python c2_run.py thin 10,11,12,13,14,15 10 out/c2_thin.json
  python c2_ablate.py thin 10,11,12,13,14,15 10 out/c2_ablate_thin.json
  python leak2.py thin 30,31,32,33,34,35 10 out/leak2_thin.json
  python c6_run.py thin 20,21,22,23,24,25,26,27,28,29 10 out/c6_thin.json
  python c6b_run.py thin 20,21,22,23,24,25,26,27,28,29 5 DEF,MED,PULSE,SITE,CODA \
      out/c6b_thin.json
  python c7_run.py STD thin 40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59 \
      5 out/c7_STD_thin.json
  python power.py out/c2_ablate_thin.json,out/c2_ablate_strong.json MT-FW BOMP-PK
(c3_run.py additionally takes the burn-in windows and the list of conditions; see its __main__ block.) The
corresponding analysers (*_analyze.py) read those JSON files and write the *_analysis.txt files of out/.
For C8 (design, regime, networks, windows per network, output; then the analyser with the list of designs):
  python c8_run.py DENSE-IID thin 60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79 \
      5 out/c8_DENSE-IID_thin.json
  python c8_analyze.py DENSE-IID
The designs and the networks of each one are in REPORT_C8.txt, sections 1 and 11.
For C9 and C10 (design, regime, networks, windows per network, output; the complete order is in c9_queue*.sh and
c10_queue*.sh, with the paths of the working folder):
  python c9_run.py DENSE thin 130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149 \
      5 out/c9_DENSE_thin.json
  NBOOT=1000 python c9_analyze.py DENSE
  python c10_run.py DENSE thin 160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179 \
      5 out/c10_DENSE_thin.json
The designs, the networks and the streams of each run are in REPORT_C9.txt (sections 1 and 9) and REPORT_C10.txt
(sections 1 and 12).
No registered file, no part of 2A to 2I and no file of verif_2e, verif_2f or verif_2i was modified.

C11 (station terms without picks): REPORT_C11.txt is the report; C11_PLAN.txt and C11_PLAN.sha256, the plan frozen
before any detection evaluation (sha256 of the body 42e66138...1c5f); c11_lib.py, c11_run.py, c11_dev.py,
c11_validate.py, c11_analyze.py (c11_analyze_frozen.py.txt = the file as it was frozen), c11_posthoc.py,
c11_power.py and c11_queue.sh, the scripts; out/c11_*.txt, the outputs. Seed 999 only, networks j = 200-209
(development), 210-239 (evaluation) and 195 (validation), TAG 4101-4105. Exploratory: its candidate, the joint
determination of hypocentres and station terms on the 20-event catalogue (S2) within the refined detector
(T3), was registered and confirmed as Task 2L on seeds 7171 and 9393 (results/task2l/, Appendix Y.28); that
protocol absorbed the draft of C10 and that of C9 was discarded.

C12 (statistic restricted to double couples; Appendix Y.29): REPORT_C12.txt is the report; c12_lib.py (statistics
FULL6, DEV5 and DC, optimiser over the double-couple manifold, greedy loop, constants and tail of the tube
formula), c12_run.py and c12_queue.sh (runs eval, nullc and nulln, one part per regime, cell and network),
c12_analyze.py, c12_validate.py and c12_mech.py (after the results), the scripts; out/c12_analysis_STD.txt,
out/c12_analysis_DENSE.txt, out/c12_null.txt, out/c12_validation.txt and out/c12_mech.txt, the outputs (the parts and
the run log remained in the working folder). Seed 999 only, networks j = 1200-1215
(evaluation, 16 per cell and regime), 1220-1225 (null with clutter), 1226-1229 (null without clutter) and 1290
(validation), TAG 4301-4305; registered detector (MT-FW cap 6, +SC, hard rule), without T3 or the terms of C11.
Result: the tube formula predicts the null law of the DC statistic (4.4-4.5 effective degrees of freedom; simulated
0.01 quantile 14.28 versus 14.17 theoretical) and the threshold drops by 10 % in energy at an equal rate of false
detections (ratio 0.906 [0.883, 0.928] at <= 0.1 and 0.898 [0.886, 0.914] at <= 0.25, noise only), but the recall
does not change: DC - FULL6 +0.005 [-0.021, +0.033] (standard, LOC+T, <= 0.25) and -0.004 [-0.018, +0.011] (dense).
NULL FOR DETECTION, POSITIVE FOR THE ANALYTICAL PIECE; no confirmation is proposed. Validations V0-V6 correct
(out/c12_validation.txt). It shared the machine with C13.

C13 (decoy calibration, target-decoy and empirical null; Appendix Y.29): REPORT_C13.txt is the report;
c13_lib.py (decoys D1p, D1o, D2s, D2r and D3; rules tdc, bh and comp; constants and TAG), c13_run.py (main
run), c13_run_d3c.py (corrected decoy D3c, post hoc), c13_twin.py (event-free twins), c13_analyze.py
(exchangeability, windows with events, twins, calibration with fpw and bootstrap), c13_tables.py, c13_report.py and
c13_report_text.py (assembly of the report) and c13_queue.sh, the scripts; out/c13_analysis.txt and
out/c13_report_tables.txt, the valid outputs; out/c13_analysis_v1_5decoys.txt (analysis prior to D3c and with the
invalid twins) and out/c13_analysis_interim4.txt (interim look with 4 networks), historical only. Seed 999
only, networks j = 1400-1411 (12 per cell and regime), TAG 1401-1406. Result: in windows without events D1p, D2s,
D2r and D3c are exchangeable with the target; with events the false detections are generated by the events (12 [9,
18] and 45 [33, 63] times those of the twin above t95, standard and dense); with the LOC rule, D1p with batches of
120 controls the FDR at q = 0.2 (0.14-0.16 pooled) and fpw reproduces without labels the operating point of the
benchmark (0.21 [0.17, 0.25] and 0.23 [0.19, 0.28] false detections per window for c = 0.25); with LOC+T it fails
(0.28-0.37 at q = 0.2), because a well-located event with the wrong origin counts as false and no moveout decoy
imitates it. PROMISING WITH CONDITIONS AND PARTIAL; the configuration C13-F and the 2M criteria of section 11 of the
report are an unregistered draft. Declared incidents: the first version of the twins was wrong (copy.deepcopy of the
generators; outputs in out/c13_twin_parts_BUGGY_v1 in the working folder, not used) and D3 and D1o had a construction
defect (circular delay over the whole trace); D3c corrects it, D1o was not repeated (E43 and E44 in
artifact_manifest_and_error_log.md).

C14 (self-interference of MT-FW: diagnosis and four remedies; Appendix Y.29): REPORT_C14.txt is the report;
c14_lib.py (per-event and twin parts, loop mt_loop14 with R1-R4, oracle deflation and exact decomposition of the
statistic), c14_run.py (dev, dev2 and eval) and c14_analyze.py, the scripts; c14_dev2.json (combinations of the
second round), c14_frozen.json (frozen configuration), C14_FREEZE.txt (sha256 of c14_frozen.json, c14_lib.py,
c14_run.py, c14_analyze.py and c10_lib.py at the freeze, 9 October 2026, 11:36:11 -05) and C14_FREEZE.txt.post
(final hash of c14_analyze.py); out/c14_analysis_dev.txt, out/c14_analysis_dev2.txt and out/c14_analysis_eval.txt, the
outputs; out/c14_analysis_eval_codetest_unread.txt, the output of the code test on two networks, not read.
Seed 999 only, networks j = 1600-1602 (development) and 1610-1621 (evaluation, 12 per cell and regime, generated
after the freeze), TAG 1651 and 1652. Diagnosis (development): at the operating point the largest class of LOC+T
false detections are well-located detections of real events (median 3.2 and 3.0 km) with the origin shifted (median
0.76 and 0.86 s, early in 75 and 73 %): 55-58 % of the false detections in standard and 80-95 % in dense; deflation
residual 0-7 %, coda 0 %. Evaluation: the remedies halve the event-induced false detections (ratio to the
twin at q95 from 6.23 to 2.52 in standard and from 10.04 to 4.24 in dense with R2 and R3) without moving the LOC+T
recall at <= 0.25 (+0.015 to +0.016 in standard, CI including 0; -0.011 to -0.002 in dense), and R1h8, R3fir2 and C23
lose one quarter of the overlapping events at any score. NEGATIVE for recall; no remedy is proposed for
registration. Incidents: interruption of the development run on 5 October at 19:02 and relaunch on the 9th with the
same code; two processes at once when launching the second round (they were killed, the directory was deleted and the
round was repeated with one process); a code error in the denominator of the far ratio, corrected after the freeze and
before reading the evaluation; the hash of c14_run.py at the freeze is not the one noted when the design was fixed on 5
October and the report does not detail the change (E45 in artifact_manifest_and_error_log.md).

From C12, C13 and C14 the scripts, the configuration JSON files, the freeze files and the text outputs of out/ were
copied to this folder; the parts (out/c12_parts, c13_parts, c13_d3c_parts, c13_twin_parts,
c13_twin_parts_BUGGY_v1, c14_dev_parts, c14_dev2_parts, c14_eval_parts), the summary JSON files and the run logs
remained in the working folder. None of the three generated or read any seed other than 999, nor modified
a registered file, a part of 2A to 2L or an earlier file of this folder.

C15 (origin error of the well-located detections; Appendix Y.29b): REPORT_C15.txt is the report; c15_lib.py
(loop mt_loop15 that records node, plateau, samples read and energies; estimators MID, LSO, LSOs, LSX, EDT and EDT2
that only replace the reported location and origin; loop mt_loop_cyc of cyclic block re-estimation),
c15_run.py, c15_analyze.py and c15_diag.py (decomposition of the origin error with the truth), the scripts;
c15_analysis_dev.txt and c15_analysis_cycdev.txt, the outputs. Seed 999 only, networks j = 1800-1802 (development
only, without a separate evaluation set), TAG 1851 and 1852, 20 STD and 20 DENSE windows per network, cell and
regime. Result: the maximum in origin is unique in 97.8 % of the atoms (the plateau of the argmax explains nothing);
with the true location, the detector's own readings give the origin within <= 0.5 s in 93 % (standard) and 76-77
% (dense), with the engine's model or with the true one, versus the 66 % and 58 % reported: the origin error is a
location error. Second-stage relocation from the readings (weighted median, robust Gauss-Newton,
EDT likelihood with and without re-reading) gives +0.021 to +0.042 in standard and -0.003 to +0.005 in dense (LOC+T at
<= 0.25), no interval excluding 0, development figures; the cyclic block re-estimation is null or negative and was
stopped with 13 of 18 parts. NULL; no confirmation is proposed.

C16 (real noise and real geometry with injected events; Appendix Y.29c): REPORT_C16.txt is the report; c16_data.py
(reading and processing of the Alto Tiberina noise, 2010-07-23, 14 IV stations from the public PyMPA repository),
c16_lib.py, c16_run.py, c16_validate.py and c16_analyze.py, the scripts; c16_queue_dev.sh and c16_queue_eval.sh;
C16_FREEZE.txt (sha256 of the code, the libraries and the processed noise, checked after the run);
c16_analysis_eval.txt and c16_validation.txt, the outputs. Seed 999, networks 2400-2401 (development) and 2402-2413
(evaluation, 12 per cell), TAG 5161-5166. Result: G1-R is strengthened with real noise (dense, LOC, <= 0.25: MT with
cap - BOMP-PK +0.259 [+0.131, +0.321], HUBER-N - BOMP-PK +0.294 [+0.178, +0.355]; with generator noise +0.069 and
+0.113) because BOMP-PK collapses with the real transients; G1-C keeps the sign without being resolved (+0.049
[-0.083, +0.149]); the operating point of the capped detector transfers to real noise with little excess. The
processed noise (out/c16_noise) and the parts remained in the working folder.

C18 - after R0: what fails on real data (Appendix Y.32). Folder c18/. Scripts in code/explore_contrib: c18_tt.py
(first-arrival times in layered media), c18a_timing.py and c18a_overlap.py (leave-one-event-out prediction of the
energy-peak readings of the Ridgecrest CAL events under five media), c18b_oracle.py and c18c_oracle.py (first-iteration
score at the true hypocentre against the maxima of event-free CAL windows, for several statistics), c18d_detect.py and
c18d_analyze.py (full detector on CAL with the cap 6, cap 24, LOG6, no cap and BOMP-PK), c18e_run.py and
c18e_analyze.py (seed 999, networks 3100-3107, TAG 5181-5186, qc1 noise of 2010-07-23: robustness of the variants to
real-noise transients). Only CAL of Ridgecrest is read (origins before 19:15); EVAL is not. Result: no 1-D medium
improves the readings (0.53 to 0.63 within 0.25 s); the cap erases the evidence of real events at the true hypocentre
(AUC 0.39 to 0.49) and LOG6 restores it (0.66 to 0.72) while keeping most of the robustness of the cap on real noise
(LOG6 - cap -0.015 [-0.099, +0.043], one seed); no variant detects the CAL events with origin time (0 or 1 of 16 at 30 false/h).
C18f (c18f_run.py, c18f_analyze.py; same windows as C18e): BOMP-PK with energies capped at 6 recovers almost all of the
gap (MT with cap - BOMP-PK with cap +0.026 [-0.009, +0.081]; per-station normalisation alone +0.037 over BOMP-PK).
C18g (c18g_oracle.py): weighting the S peak vectors by 0.5 or removing them lowers the AUC (cap 0.17-0.32, LOG6 0.24-0.55)
and no window maximum falls within 5 km of the event.
C18h (c18h_kernel.py, c18h_detect.py): a 3 Hz matched-filter kernel with LOG6 recalls 3 of 16 CAL events with origin
time at 30 false/h on U8 (0 on N8); the registered gate (5 of 16) is not passed.
Outputs: c18a_*.json, c18a_summary.txt, c18a_overlap_N8.txt, c18b_*.json, c18c_*.json, c18d_*.json,
c18d_analysis.txt, c18e_parts/, c18e_analysis.txt and the run logs.
