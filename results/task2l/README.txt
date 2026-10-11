Results of Task 2L: station terms estimated without picks (S2, joint determination of hypocentres and
terms on the catalogue of 20 events) within the refined detector (MT-FW with cap 6 and refinement T3), against
the registered corrections +SC (S0). In the proposal, hypothesis G1-C; Appendix Y.28 of the supplementary volume;
tasks_and_real_data_plan.md, section Task 2L. The hypothesis came out of exploration C11 on seed 999
(results/explorations/REPORT_C11.txt and C11_PLAN.txt, plan frozen with hash before any detection result),
whose primary contrast was +0.070 [+0.010, +0.120]; the registration plans with +0.05.

STATUS. Task executed and verified. Seeds 7171 and 9393, never used. Rules 2L-a, 2L-b, 2L-c and 2L-d and the
secondary on BOMP-PK are met in both seeds. From 2L onwards the primary hit rule is LOC+T (epicentre within 5 km
and origin error <= 0.5 s); LOC (epicentre only) is the second reading.

Decision rules: PREREGISTRO_1j_1k.txt, section "Task 2L" (4 October 2026, 21:10:33 UTC, wall-clock time),
written before running any code on seeds 7171 and 9393. SHA-256 13a6e9b1...a62a in
PREREGISTRO_hashes.txt, delivered to the user through the chat before the run (registration outside the container, not
public and not held by third parties; the verification checks from the disk the hash and its position in the registration, not the moment
it was sent). The drafts of C9 and C10 that named these seeds were neither registered nor run; the protocol
of 2L absorbs that of C10 and that of C9 is discarded. The registration was not modified. Afterwards only the paragraph
"Result of 2L" (5 October, 02:00:25 UTC, hash 2934487e...9c3d) and that of the "Independent verification of 2L"
(02:36:52 UTC, hash ca3f62fc...2d38) were added. The copy of PREREGISTRO_1j_1k.txt in this folder contains those sections.

Rules (contrast MT@S2+T3 - MT@S0+T3, recall at <= 0.25 realised false detections per window, exact FROC, one operating
point per method, hierarchical bootstrap of 2000 replicates with strata regime x cell, networks and windows; each
rule in BOTH seeds):
    2L-a  superiority, standard windows, LOC+T, pooled: one-sided 97.5 lower bound > 0 and point >= +0.03
    2L-b  superiority, dense windows, LOC+T, pooled: one-sided 97.5 lower bound > 0
    2L-c  non-inferiority per regime (thin, strong), margin 0.02, LOC+T, standard and dense: 95 lower bound > -0.02
    2L-d  non-inferiority, margin 0.02, LOC, pooled, standard and dense: 95 lower bound > -0.02
    registered secondary: BOMP@S2+T3 - BOMP@S0+T3, standard, LOC+T, 60 networks: 97.5 lower bound > 0
No superiority was registered with the LOC rule. Descriptive: total contrast against the registered detector
(MT@S2+T3 - MT@S0), per cell, cap 0.1, error of the terms, epicentral and origin error of the hits and
fraction of the oracle gain.

Design in brief. Per seed, networks task2a_confirm j = 0..159 of each cell (ladder_4, ladder_6, unilateral_c40) and
regime (thin, strong) with 5 standard windows (0 to 2 events, kmax 5): 4800 windows, 3826 and 3863 events; networks
j = 0..99 with 5 dense windows (3 to 5 events within 9 km, kmax 7): 3000 windows, 11978 and 12000 events. New
streams (TAG 6101 standard, 6102 dense). Station delays with deviation 0.3 s; catalogue of 20 events with
errors of 1 km and 0.1 s, as in 2A; gen_realistic by default. MT-FW with cap 6 in three forms on all networks
(S0 with the hard rule, which is the registered detector of the previous tasks; S0 with T3; S2 with T3); BOMP-PK in the
same three forms on the standard windows of j = 0..59 (1800 windows per seed); oracle of true
delays with T3 on j = 0..29 (900 windows per design and seed), descriptive. Constants are those of the frozen plan of
C11; none free.

Freeze and times (date -u, container clock; task2l_freeze_time.txt, task2l_run.log, log2l.txt):
    21:10:33 UTC  registration of 2L with hash
    21:14:09 UTC  controls on seed 999 (task2l_run.py control, 18 jobs)
    21:24:48 UTC  task2l_code_hashes.txt: SHA-256 of the 37 code files loaded by the runner, the
                  controls and the analysis
    21:24:53 UTC  first job on 7171/9393 (2 processes); estimate 4.62 h, budget rule of 10 h not
                  triggered
    22:47:19 UTC  interruption when the session was suspended; resumed at 22:48:13 and relaunched at 22:49:32 and
                  23:07:22 UTC, each time from the parts already written and with the code unchanged (sha256sum -c,
                  37 of 37)
    01:48:15 UTC  last of the 1920 parts (5 October)
    02:00:25 UTC  result paragraph in the registration, with hash
    02:36:52 UTC  verification paragraph in the registration, with hash
Processor time: 4.34 and 4.33 h per seed.

Results (analyze_2l_out.txt; recall at <= 0.25 with LOC+T, MT@S0+T3 -> MT@S2+T3):
    standard  7171: 0.346 -> 0.403    9393: 0.350 -> 0.411
    dense     7171: 0.124 -> 0.162    9393: 0.120 -> 0.151
    Rule                               7171                       9393                       met
    2L-a standard, LOC+T               +0.057 [+0.034, +0.081]    +0.062 [+0.036, +0.083]    yes, both
    2L-b dense, LOC+T                  +0.038 [+0.024, +0.054]    +0.030 [+0.015, +0.044]    yes, both
    2L-c lower bounds per regime       +0.0069 to +0.0517         +0.0102 to +0.0461         yes, both
    2L-d LOC, standard / dense         +0.036 / +0.020            +0.050 / +0.035            yes, both
    secondary BOMP-PK, 60 networks     +0.068 [+0.026, +0.107]    +0.071 [+0.029, +0.114]    yes, both
Descriptive. Against the registered detector (MT@S2+T3 - MT@S0), standard: +0.088 [+0.063, +0.116] and +0.120 [+0.093,
+0.146]. Per cell, standard, LOC+T: unilateral with 8 stations +0.108 [+0.057, +0.151] and +0.113 [+0.067, +0.157];
6 stations +0.048 [+0.005, +0.089] and +0.027 [-0.023, +0.069]; 4 stations +0.021 [-0.013, +0.055] and +0.018
[-0.016, +0.054]. Median error of the terms (rms with the median removed): thin 0.147 and 0.151 s with S0, 0.082 and
0.084 s with S2; strong 0.216 and 0.207 s with S0, 0.155 and 0.153 s with S2. Oracle (30 networks per cell and regime):
fraction of the gain captured by S2 0.56 and 0.61 in standard windows, 0.55 and 0.28 in dense windows; the oracle remains +0.052 and
+0.060 ahead in standard windows. SD of the contrasts of the rules: planned 0.0122 (standard) and 0.0083
(dense); measured 0.0119 and 0.0121, and 0.0077 and 0.0073.

Controls (task2l_control.txt, seed 999; all correct). (1) The runner reproduces the C11 pilot exactly
in 18 parts (3 networks x 3 cells x 2 regimes): terms S0, S2 and oracle, levels of T3 and errors in 18 of 18;
truth, origins and atoms of each method in 90 of 90 standard windows and 90 of 90 dense ones. (2) S0 of the runner is the
registered correction in 18 of 18 parts, and the corrections of the engine are those stored by 2A on seed 999 in
30 of 30 cases. (3) The seed enters networks, catalogues and streams as designed (checked with a test
seed that is none of the registered ones). (4) The analysis reproduces the published recalls of the pilot, and its
bootstrap gives +0.070 [+0.011, +0.120] for the primary contrast of C11 against +0.070 [+0.010, +0.120] published.

Independent verification (verif_2l/task2l_verification.txt; no file of 2L or of the registration modified).
Verdict: CONFIRMED POSITIVE within a declared scope.
    Leakage: none. S2 recomputed in new processes on copies of its inputs, without the true delays or the true
    catalogue, identical bit for bit in 300 of 300 network x regime units (networks 100-124 of both seeds); the
    catalogue that the estimators receive carries the registered errors (0.998 km per coordinate and 0.099 s).
    Origin time: same small bias in both arms; S2 narrows the origin errors (hits within 0.3 s from
    0.47 to 0.53 in standard windows); positive contrast with hit windows of 0.3 to 1.0 s and with the median origin
    of each method subtracted (+0.060 and +0.062 in standard windows).
    Attribution: with the T3 level of S0 in both arms, the terms alone give +0.045 [+0.005, +0.089] in 1500
    standard windows of 25 networks; the change of level alone, +0.000 and -0.009.
    Provenance: 37 of 37 hashes; chain of prefixes of the registration; no part prior to the freeze; no prior
    use of the seeds; six parts regenerated from scratch identical; the analysis run again gives a byte-for-byte
    identical output.
    Statistics: own loader, FROC and bootstrap reproduce 2L-a (+0.0572 and +0.0616) and 2L-b (+0.0380 and +0.0305) to
    four decimals, with intervals [+0.032, +0.082] and [+0.037, +0.083] for 2L-a.
    Scope: the pooled gain is in essence the unilateral cell with 8 stations (+0.108 and +0.113); with 4
    stations there is no demonstrated gain; S2 worsens the terms in about 22 % of the networks (3-4 % of the
    unilateral ones with weak clutter, 19-27 % of the unilateral ones with strong clutter, 12-28 % with 6 stations, 26-37 %
    with 4), and in the 6-station ones of that group the standard recall drops (-0.019 and -0.049); the effect is conditional
    on T3; it assumes a catalogue of 20 events with errors of known size (1 km and 0.1 s); fraction of the oracle
    gain 0.28 to 0.61; synthetic and from a single family of generator.

Reading. G1-C is supported by a preregistered synthetic test, with the hash of the registration delivered outside the
container before running, the frozen code and independent verification: within the refined detector, the
station terms estimated without picks from the catalogue of 20 events raise the recall with origin time,
about +0.06 in standard windows and +0.03 in dense windows at <= 0.25. It is a modest and narrow calibration result: the
joint determination of hypocentres with station terms is classical in pick-based location, and what is measured
is its estimation without picks and its effect within a pick-free detector, not an algorithm. Absolute recall with LOC+T
is low (0.403 and 0.411 in standard windows, 0.162 and 0.151 in dense windows). Next steps: sensitivity of S2 to catalogues with
errors of unknown size, biased or with fewer events; dependence on the geometry; combination with the
bounded-influence statistic of G1-R. The explorations C12 to C14, later (explore_contrib/, Appendix
Y.29), add the origin time of the well-located detections as a prior target.

Files
    README.txt                 this file
    PREREGISTRO_1j_1k.txt     copy of the preregistration: section "Task 2L" (protocol and rules), "Result of 2L" and
                              "Independent verification of 2L"
    PREREGISTRO_hashes.txt    SHA-256 hashes of the preregistration with wall-clock time (13a6e9b1 registration of 2L; 2934487e
                              with the result; ca3f62fc with the verification)
    analyze_2l_out.txt        output of analyze_2l.py: data controls, recall per method, rules 2L-a to 2L-d,
                              BOMP secondary, descriptive contrasts (per regime, per cell, BOMP family, oracle
                              arm, errors of the terms, hit atoms), bootstrap SD against the
                              planned one and times
    task2l_control.txt        controls on seed 999 (reproduction of the C11 pilot, S0 = registered +SC,
                              seed mechanism, analysis against the pilot, time and budget rule)
    task2l_code_hashes.txt    SHA-256 of the 37 code files frozen before running
    task2l_freeze_time.txt    time of the freeze (21:24:48 UTC), of the first job and of each resumption
    task2l_run.log            log of the main run (1920 parts between 4 October 21:25:19 and 5
                              October 01:48:15 UTC; interrupted when the session was suspended and resumed from the parts)
    log2l.txt                 the same log preceded by the lines of the controls on seed 999 (from
                              21:14:09 UTC); it is the one cited by Appendix Y.28
    verif_2l/                 independent verification:
        task2l_verification.txt     report (in English) with the verdict and the scope
        load_verif.py               own loader, hit labelling and FROC, without code from the analysis
        rules_verif.py, rules_verif.txt          rederivation of 2L-a, 2L-b and 2L-d and own bootstrap
        leak_attrib_verif.py, leak_attrib_log_verif.txt   leakage (S2 recomputed without truth in 300 units) and
                                    attribution arms terms / T3 level
        attrib_analyze_verif.py, attrib_analyze_verif.txt  analysis of those arms
        origin_verif.py, origin_verif.txt        signed origin errors, other hit windows, bias removed
        symmetry_verif.py, symmetry_verif.txt    symmetry and integrity of the stored parts
        scope_verif.py, scope_verif.txt          scope: per cell, error of the terms against gain, networks in
                                    which S2 is worse
        regen_verif.py, regen_jobs_verif.txt, regen_out_verif.txt, compare_regen_verif.py,
        compare_regen_verif.txt     regeneration from scratch of 6 parts and comparison
        modules_verif.txt           modules loaded by the runner against the list of hashes
        analyze_2l_rerun_verif.txt, analyze_2l_rerun_stdout_verif.txt   second run of the frozen analysis
The parts (task2l_parts/), the cache and the parts of the verification remain in the working folder and are not copied
here; the scripts task2l_run.py, task2l_control.py and analyze_2l.py are in the root of outputs/.

Notes. The time of the resumption at 22:48:13 and of the relaunch at 22:49:32 UTC in task2l_freeze_time.txt
has no line of its own in the log, which shows a single start at 22:49:51 UTC, and no part was written
between 22:47:19 and 22:50:35 UTC; the verification judges it harmless. The gauge of S2 is "median shift of the
origin of the relocated events = 0", declared in C11_PLAN.txt before freezing, not zero mean per phase.
