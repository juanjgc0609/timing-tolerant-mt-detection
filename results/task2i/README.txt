Results of Task 2I: a robust certificate of the moment-tensor subspace with a Huber fit per
station and phase and level fixed by rule from the noise (HUBER-N), against BOMP-PK+SC and MT-FW+SC with cap 6
(in the proposal, hypothesis G1-R). Appendix Y.26 of the supplementary volume; tasks_and_real_data_plan.md, section Task
2I. The hypothesis came out of exploration C8 on seed 999 (results/explorations/REPORT_C8.txt, Appendix Y.25),
where HUBER-N was the best of 11 candidates, which the registration declares (it was planned with +0.04).

STATUS. Task executed and verified; secondary arms analysed. The PRIMARY arm carries the three registered
rules and has independent verification. The secondary arms (thresholds 2.5003, 3.3682 and 5; site
deviation 0.25 and 0.8; 20 networks per cell and regime, j = 0..19), descriptive and without a rule, ran after the primary
with the same frozen code and are in analyze_2i_out.txt; they have no independent verification, except for the 26
parts of the thr arm that the verification compared with the primary. In analyze_2i_primary_out.txt, written before they
finished, their sections say "not run".

Decision rules: PREREGISTRO_1j_1k.txt, section "Task 2I" (2 October 2026, 18:00:05 UTC, wall-clock
time), written before running any code on seeds 5151 and 6262, which had never been used. SHA-256
a983beac...6f97 in PREREGISTRO_hashes.txt, delivered to the user through the chat before any run (registration
outside the container, not public and not held by third parties; from the disk one can check the hash and its position in the
registration, not the moment it was sent). Afterwards only the paragraph "Result of 2I, primary arm"
(22:08:07 UTC, hash fbb9c579...2100) and that of the "Independent verification of 2I" (22:16:47 UTC, hash
c2745a94...eb1a) and the paragraph "Result of 2I, secondary arms" (3 October, 05:54:02 UTC, hash
055587ca...b4e2) were added. The copy of PREREGISTRO_1j_1k.txt in this folder contains those sections.

Freeze, times and interruptions (date -u, container clock; task2i_freeze_time.txt, log2i.txt):
    18:00:05 UTC  registration of Task 2I with hash
    18:14:44 UTC  freeze: task2i_code_hashes.txt, SHA-256 of 27 files (runner, controls, analysis,
                  run_queue_2i.sh and the modules they load)
    18:14:45 UTC  first job on 5151/6262 (run_queue_2i.sh: first --arms prim, then thr,s25,s80)
    18:14:58 UTC  first part
    20:33:25 UTC  process interrupted; resumed at 20:34:09 UTC from the parts; sha256sum -c correct
    20:46:15 UTC  interrupted when the session was suspended; resumed at 21:29:10 UTC; code unchanged
    22:04:25 UTC  last part of the primary arm (1440/1440)
    22:04:55 UTC  first part of the secondary arms (thr; then s25 and s80)
    22:07:53 UTC  analyze_2i.py on the complete primary arm -> analyze_2i_primary_out.txt
    03:33:41 UTC  (3 October) last part of the secondary arms (720/720); the job finished on its own
    03:35:49 UTC  (3 October) analyze_2i.py, unchanged (sha256sum -c correct), on the complete task ->
                  analyze_2i_out.txt
    05:54:02 UTC  (3 October) paragraph of the secondary arms in the registration, with hash
The three sessions wrote 1101, 83 and 256 parts. Each resumption used the same command and recomputed no
part already written (one part per network). The sample size was fixed at 120 networks: the interruptions gave
no opportunity for optional stopping.

Design in brief. Seeds 5151 and 6262, each an independent replicate; cells ladder_4, ladder_6 and
unilateral_c40; regimes thin and strong; networks task2a_confirm with j = 0..119 (120 per cell and regime); catalogue
of 20 events per network and station corrections as in 2A; gen_realistic by default (site factors
drawn per window, sd 0.5; measured mean sd 0.467 and 0.464).
    dense windows       5 per network with c2_run.draw_dense (3 to 5 events), new TAG 5201, kmax 7:
                        3600 windows per seed (14477 and 14407 events)
    standard windows    5 per network with draw_sources (0 to 2 events), new TAG 5202, kmax 5:
                        3600 windows per seed (2938 and 2891 events)
Methods (all four on the same window, each on its own copy; same grid of 675 nodes, origins
DEFAULT_TAUS, corrections, 5 km exclusion and hit = epicentre within 5 km):
    BOMP-PK+SC        xc.bomp_pk
    MT-FW+SC          c8_lib.mtfw_general "ls", cap CAPV = 6 per station and phase
    MT-FW-NOCAP+SC    the same without cap
    HUBER-N+SC        c8_lib.mtfw_general "huber_n": Huber function of the norm of the residual of each 3-component
                      vector per station and phase; threshold c = square root of the 0.95 quantile of chi-square(3) =
                      2.7955 in units of the noise deviation of the matched filter, fixed by rule and not
                      tuned; 4 IRLS passes from the least-squares fit, weights min(1, c/|r_b|);
                      statistic = Huber loss of the data minus that of the residual; deflation with the robust
                      fit; no cap on the data
Rules (each in BOTH seeds, at cap 0.25, pooling regimes and cells): 2I-a, superiority in dense windows,
HUBER-N+SC - BOMP-PK+SC with one-sided 97.5 lower bound > 0 and point >= 0.03. 2I-b, non-inferiority in dense windows,
HUBER-N+SC - MT-FW+SC with one-sided 97.5 lower bound > -0.03. 2I-c, equivalence in standard windows, HUBER-N+SC -
BOMP-PK+SC with 90 CI within [-0.05, +0.05] (TOST). Statistics: exact FROC with one operating point per
method; hierarchical bootstrap of 2000 replicates (strata regime x cell, networks and windows), same resamples
for all methods, operating point re-chosen per replicate.

Outputs
    analyze_2i_out.txt            the same output on the complete task: the content of analyze_2i_primary_out.txt
                                  unchanged in the primary arm, plus the sections on the Huber threshold and on the sd of
                                  the site factors and the times of the secondary arms
    analyze_2i_primary_out.txt    data controls, recall per method, rules 2I-a, 2I-b and 2I-c, contrasts per
                                  regime and per cell, descriptive replication of 2F, hit rule with origin
                                  time, bootstrap SD against the planned one, times, verdicts
    task2i_control.txt            controls on seed 999
    task2i_code_hashes.txt        SHA-256 of 27 files
    task2i_freeze_time.txt        time of the freeze and of each resumption
    PREREGISTRO_1j_1k.txt         copy of the preregistration
    PREREGISTRO_hashes.txt        SHA-256 hashes of the preregistration with wall-clock time
    verif_2i/                     independent verification (task2i_verification.txt and its scripts *_verif.py with
                                  their outputs)
The parts (task2i_parts/: 66 control, 1440 of the primary arm and 720 of the secondary arms), log2i.txt and
log2i_run_stdout.txt remain in the working folder; they are not copied here. Each part
records the OpenBLAS kernel (Prescott in that of numpy).

Verdict (analyze_2i_primary_out.txt)
    rule     5151         6262
    2I-a     met          met
    2I-b     met          met
    2I-c     met          met
Registered interpretation: a, b and c are met -> the claim "a Huber fit scaled by the noise
replaces the calibrated cap" is confirmed in this synthetic design.

Recall at 0.1 / 0.25, pooling regimes and cells:
    windows, seed       BOMP-PK+SC    MT-FW+SC      MT-FW-NOCAP+SC  HUBER-N+SC
    dense 5151          0.076/0.177   0.157/0.243   0.077/0.175     0.139/0.235
    dense 6262          0.092/0.198   0.160/0.253   0.091/0.195     0.151/0.251
    standard 5151       0.361/0.479   0.354/0.469   0.360/0.480     0.383/0.489
    standard 6262       0.355/0.488   0.367/0.485   0.354/0.488     0.385/0.496

2I-a, HUBER-N+SC minus BOMP-PK+SC, dense, cap 0.25, 95 hierarchical CI (one-sided p 0.0005 in both):
    5151 +0.058 [+0.039, +0.076]    6262 +0.053 [+0.037, +0.070]
    at cap 0.1: +0.063 [+0.049, +0.077] and +0.060 [+0.038, +0.077]
2I-b, HUBER-N+SC minus MT-FW+SC, dense, cap 0.25 (lower bounds -0.0234 and -0.0184, margin -0.03):
    5151 -0.008 [-0.023, +0.006]    6262 -0.002 [-0.018, +0.012]
    at cap 0.1: -0.018 [-0.034, -0.002] and -0.009 [-0.031, +0.009]
2I-c, HUBER-N+SC minus BOMP-PK+SC, standard, cap 0.25:
    5151 +0.0099, 90 CI [-0.0100, +0.0272]    6262 +0.0073, 90 CI [-0.0132, +0.0235]

Descriptive. Per regime, in dense windows, HUBER-N+SC - BOMP-PK+SC: thin +0.068 / +0.054, strong +0.051 / +0.050
(all CIs exclude 0); in standard windows: thin +0.006 / -0.032 [-0.063, -0.004], strong -0.003 / +0.003. Per
cell, in dense windows, positive in 12 of 12; the smallest with 4 stations (thin +0.039 [-0.000, +0.085] and +0.045; strong
+0.010 [-0.011, +0.057] and +0.029), and +0.054 to +0.080 with 6 stations and with 8 unilateral ones. Replication of 2F in
these windows: MT-FW+SC - BOMP-PK+SC +0.082 and +0.069 (thin), +0.054 and +0.049 (strong); MT-FW-NOCAP+SC -
BOMP-PK+SC -0.002 and -0.003; in standard windows -0.010 (90 CI [-0.031, +0.011]) and -0.003 (90 CI [-0.025, +0.016]).
Time per dense window (seconds, one process per core): HUBER-N 0.927, MT-FW 0.103-0.104, BOMP-PK
0.587-0.602.

Secondary arms (analyze_2i_out.txt; descriptive, 600 dense and 600 standard windows per seed, 2420 and 2383
events in dense windows, 504 and 501 in standard ones; differences at cap 0.25, 95 CI; 5151 / 6262). The primary rules do not
change.
    Huber threshold, dense windows (HUBER-N regenerated identical to that of the primary in 600/600 windows):
        c = 2.5003 - c = 2.7955     +0.003 [-0.012, +0.021]    +0.001 [-0.008, +0.016]
        c = 3.3682 - c = 2.7955     +0.005 [-0.015, +0.027]    +0.006 [-0.012, +0.020]
        c = 5      - c = 2.7955     -0.013 [-0.039, +0.017]    -0.023 [-0.049, +0.009]
        c = 5      - MT-FW          -0.042 [-0.076, -0.002]    -0.029 [-0.070, +0.013]
        c = 5      - BOMP-PK        +0.038 [+0.004, +0.073]    +0.044 [+0.010, +0.074]
      in standard windows, c = 2.5003 and 3.3682 differ from the registered one by +0.000 and +0.008 in both seeds.
    sd of the site factors, dense windows (measured mean sd 0.236 and 0.233; 0.755 and 0.744):
        sd 0.25  HUBER-N - BOMP-PK  +0.036 [+0.004, +0.070]    +0.002 [-0.035, +0.044]
                 MT-FW - BOMP-PK    +0.005 [-0.038, +0.041]    +0.009 [-0.024, +0.053]
        sd 0.8   HUBER-N - BOMP-PK  +0.100 [+0.054, +0.143]    +0.072 [+0.031, +0.115]
                 HUBER-N - MT-FW    -0.055 [-0.113, +0.002]    -0.069 [-0.121, -0.019]
                 MT-FW - BOMP-PK    +0.155 [+0.096, +0.207]    +0.141 [+0.095, +0.189]
      in standard windows no HUBER-N - BOMP-PK contrast excludes zero (sd 0.25: +0.030 and -0.002; sd 0.8: +0.016 and
      +0.032).
    Reading: the Huber level does not matter between 2.5 and 3.4 noise deviations and with 5 the gain is smaller; with
    sd 0.25 the cap contributes nothing and HUBER-N little; with sd 0.8 HUBER-N beats BOMP-PK but stays below the cap, with
    the interval below zero in 6262 and lower bounds far from the margin of -0.03, so that the non-inferiority of
    2I-b does not extend to that dispersion (with 20 networks the arm would not have shown it with sd 0.5 either: -0.029
    [-0.070, +0.013] and -0.007 [-0.049, +0.030]).
    Hit rule with origin within 0.5 s, per regime, HUBER-N - BOMP-PK in dense windows: thin +0.040 [+0.027, +0.054] and
    +0.020 [+0.003, +0.040]; strong +0.022 [+0.009, +0.035] and +0.017 [+0.004, +0.030].
    Compute of the secondary arms: 1.52 and 1.57 h of CPU (thr), 1.08 and 1.08 h (s25), 0.81 and 0.80 h (s80).

Controls (task2i_control.txt, seed 999; all correct). (1) The C8 pilot is reproduced exactly: 90/90
dense windows and 72/72 standard ones for the four methods. (2) Secondary thresholds: 60/60 windows identical;
HUBER-N of the control arm 48/48. (3) site sd 0.25 and 0.8: 24/24 (mean sd of the site log 0.232 and
0.702). (4) c8_lib.mtfw_general("ls") gives the hits and positions of c6_lib.mtfw_dense in 60/60 windows with and
without cap (maximum relative difference of scores 2.2e-15). (5) The FROC of the analysis (analyze_2f.Pool) is equal to
xc.recall_caps and to the published values of the pilot. (6) The bootstrap of the analysis on the pilot gives +0.064
[+0.027, +0.104] against +0.064 [+0.028, +0.105] published. (7) Time estimate 5.27 h on 2 cores: the budget
rule of 8 h is not triggered and the secondary arms are run.

Independent verification (verif_2i/task2i_verification.txt; no file of 2I or of the registration modified).
Verdict: CONFIRMED POSITIVE for the registered claim.
    Provenance: hashes of 27 files without failures; chain of 11 growing prefixes of the registration, a983beac
    ends at the end of the section Task 2I; no part prior to the freeze; internal time and file
    time equal within 1.0 s in the 1440 parts; exactly two gaps of more than 60 s (82 s and 2599 s),
    which coincide with the two interruptions; no prior data of 5151/6262. Not verifiable from the disk: that
    the hash was sent before the run.
    Regenerated parts: 6 parts, 2 from each of the 3 sessions (among them the last one before each
    interruption and the first one after the second), regenerated from scratch in new processes and identical to the
    stored ones; the 26 parts of the thr arm already written contain HUBER-N+SC identical to that of the primary in 260
    windows. (regen_verif.py, regen_out_verif.txt)
    Symmetry: same data, grid, origins, corrections, kmax, exclusion and hit function for the four
    methods; hit labels recomputed without differences; one operating point per method, chosen in the same way
    for all.
    Statistics: own loader, FROC and bootstrap reproduce 2I-a, 2I-b and 2I-c (for example 5151, 2I-a:
    +0.0582 [+0.0402, +0.0753] against +0.058 [+0.039, +0.076] registered). (stats_verif.py, stats_verif_out.txt)
    Scope (exploratory, stored parts only; run_explore_verif.py, explore_verif_out.txt):
        hit rule with |origin error| <= 0.5 s: HUBER-N+SC - BOMP-PK+SC drops to +0.031 and +0.019 (still
        above 0; the point criterion of 0.03 would fail in 6262); with 1.0 s, +0.054 and +0.047;
        equivalence in standard windows: met when pooled, not in thin 6262 (-0.032 [-0.064, -0.003]);
        recall with all atoms accepted shows no increase (BOMP-PK 0.623 and HUBER-N 0.622 in 5151; 0.628 and 0.616
        in 6262): HUBER-N ranks and chooses better, it does not find more events; the gain is split between
        reliability of the score across windows (+0.027 / +0.029) and choice of atoms (+0.030 / +0.018);
        at cap 0.1 it stays somewhat below the cap in 5151 (-0.018 [-0.034, -0.003] in the verification);
        cost: 0.93 s per dense window against 0.10 s for MT-FW, about 9 times;
        low absolute recall: 0.18 -> 0.24-0.25 in dense windows at cap 0.25; ~0.10 with the 0.5 s rule.

Reading. G1-R is supported by two preregistered synthetic tests, with the hash of the registration delivered
outside the container before running and independent verification: 2F with a calibrated cap and 2I with a Huber
level fixed by rule from the noise, which removes the dependence on a tuned constant. What
holds: a bounded-influence moment-tensor subspace certificate improves recall at a fixed rate
of false events in dense sequences of overlapping events against BOMP-PK with equal temporal tolerance, without
harming isolated events when the regimes are pooled. Limits: synthetic and from a single family of
generator (gen_realistic; the two seeds are independent draws from that family); modest size (+0.05 to
+0.06 with the registered rule, +0.02 to +0.03 with origin time); about 9 times the cost of MT-FW; operating
points chosen on the evaluation; with site sd 0.8 the Huber fit stays below the cap; the
Huber fit is standard robust statistics, so the novelty is the measured design result, not an
algorithm. The margin of 0.08 to 0.24 below the realistic ceiling remains
open: about three quarters remain unexplained by amplitude treatments (Appendix Y.25); the
explorations C9 and C10 locate it in the calibration of the station terms and in the amplitude evidence
(Appendix Y.27, explore_contrib/).

Notes. The plan of REPORT_C8.txt named networks 0..79; the registration, which is the governing one, names 0..119, and the
run used 0..119. The result paragraph of the registration rounds to +0.024 the upper limit of the 90 CI of 6262,
which in the output is +0.0235. The time line of analyze_2i_primary_out.txt (section "Times") is in
seconds per window, not in recall.

Reproduction (from the scripts folder, with explore_contrib/ alongside):
  OPENBLAS_CORETYPE=Prescott python task2i_run.py control --workers 2      (seed 999)
  python task2i_control.py                                                 (writes task2i_control.txt)
  OPENBLAS_CORETYPE=Prescott python task2i_run.py run --seeds 5151,6262 --arms prim --workers 2
  OPENBLAS_CORETYPE=Prescott python task2i_run.py run --seeds 5151,6262 --arms thr,s25,s80 --workers 2
  python task2i_run.py status
  python analyze_2i.py --out analyze_2i_primary_out.txt                    (with only the primary arm written)
  python analyze_2i.py                                                     (complete task -> analyze_2i_out.txt)
  sha256sum -c task2i_code_hashes.txt
Each job checks that the OpenBLAS of numpy loaded Prescott and saves it in its part; run resumes from the
parts already written in task2i_parts/ (one per network), which is how the run was continued after each interruption.
Dependencies: those of Task 2F and of the exploration (results/explorations/README.txt), because task2i_run.py imports
explore_contrib/xc.py, c2_run.py and c8_lib.py and task2a_confirm.py, and these load the other modules listed in
task2i_code_hashes.txt (among them c6_lib.py, c8_run.py, task2d_methods.py, task2c_methods.py,
task2e_methods.py, certified_fw.py, task1x_mtglrt.py, gen_realistic.py, station_corrections.py and
task0_common.py); analyze_2i.py imports analyze_2f.py; threadpoolctl.
