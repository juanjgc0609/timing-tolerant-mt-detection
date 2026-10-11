Results of Task 2F: the per-station cap of the moment-tensor subspace statistic in dense sequences
of overlapping events (H-CAP; in the proposal, hypothesis G1-R). Appendix Y.24 of the supplementary volume;
tasks_and_real_data_plan.md, section Task 2F. The hypothesis came out of the exploration on seed 999
(explore_contrib/, Appendix Y.23), where it was chosen among several candidates, which the registration declares.

Decision rules: PREREGISTRO_1j_1k.txt, section "Task 2F" (28 September 2026, 03:34:51 UTC, wall-clock
time), written before running any code on seeds 4242 and 8888, which had never been used. SHA-256
3d8e4c9b...2f5f in PREREGISTRO_hashes.txt, delivered to the user through the chat before the run (registration outside the
container, not public and not held by third parties). The registration was not modified: no infrastructure problem required a
note. Afterwards only the paragraph "Result of 2F" (06:59:18 UTC, hash dc46cb38...7376) and that of the
"Independent verification of 2F" (07:26:18 UTC, hash ef788a5b...262c) were added. The copy of PREREGISTRO_1j_1k.txt in this
folder contains those sections.

Freeze and times (date -u, container clock; task2f_freeze_time.txt, log2f.txt):
    03:39:10 UTC  controls on seed 999 (task2f_run.py control)
    03:49:08 UTC  task2f_code_hashes.txt: SHA-256 of the 26 local modules loaded by task2f_run.py,
                  task2f_control.py and analyze_2f.py, plus explore_contrib/power.py (the analysis included)
    03:49:14 UTC  single launch on 4242/8888
    06:52:31 UTC  last part (1344/1344)
    06:53:11 UTC  sha256sum -c task2f_code_hashes.txt: 26/26 correct
    06:54:13 UTC  analyze_2f.py, a single time -> analyze_2f_out.txt

Design in brief. Seeds 4242 and 8888, each an independent replicate; cells ladder_4, ladder_6 and
unilateral_c40; regimes thin and strong. Networks task2a_confirm.network(semilla, celda, j), shared by the two
regimes as in 2A; catalogue of 20 events per network and regime and station corrections as in 2A.
    dense arm          j = 0..95 (96 networks per cell, regime and seed), 10 windows per network with c2_run.draw_dense
                       (3 to 5 events, separations >= 5 km within 9 km, origins in [5.2, 7.8] s with gaps
                       >= 0.5 s), new TAG 4201, kmax 7: 5760 windows per seed (23013 and 22996 events)
    standard arm       j = 96..111 (16 per cell and regime), 10 windows per network with draw_sources (0 to 2 events),
                       new TAG 4202, kmax 5: 960 windows per seed (766 and 749 events)
Methods (calls of explore_contrib/c2_ablate.py, unchanged; grid of 675 nodes, origins DEFAULT_TAUS, registered
homogeneous model, tolerance +-0.25 s, 5 km exclusion, hit = epicentre within 5 km):
    MT-FW+SC          c6_lib.mtfw_dense, cap CAPV = 6 per station and phase (fixed in 1x, not tuned here)
    MT-FW-NOCAP+SC    the same without cap
    BOMP-PK+SC        task2d_methods.bomp_pk_record
    BOMP-PK-CAP+SC    BOMP-PK with the cap only in the selection statistic (c2_ablate.PKDictCap)
Rules: 2F-a (primary), per regime, MT-FW+SC - BOMP-PK+SC at cap 0.25 in dense windows, pooling the 3 cells,
with lower bound of the 95 CI > 0 and point >= 0.03 in both seeds; strict reading with Holm over the 4 tests
(2 regimes x 2 seeds, one-sided 0.025). 2F-b (mechanism), MT-FW+SC - MT-FW-NOCAP+SC with lower bound > 0 in
both seeds (pooling regimes) and |MT-FW-NOCAP+SC - BOMP-PK+SC| <= 0.02 (point) in both. 2F-c (no harm),
TOST equivalence with margin 0.05 of MT-FW+SC and BOMP-PK+SC in standard windows, pooling cells and regimes, in
both seeds. Statistics: exact FROC with one threshold per method over the pooled windows; hierarchical bootstrap
of 2000 replicates (cells as strata, networks and windows), same resamples for all methods, operating
point re-chosen per replicate; p = (1 + #{d* <= 0}) / 2001 (floor 0.0005). Size: with the pilot SDs (0.0383
thin, 0.0353 strong), 96 networks per cell give power ~0.95 per seed for delta 0.035.

Outputs
    analyze_2f_out.txt            data controls, recall per method, rules 2F-a, 2F-b and 2F-c, descriptive
                                  differences, tables per cell and times, verdicts
    task2f_report.txt             report (in English): verdicts, implemented design, controls, freeze,
                                  deviations a-h, results, interpretation and limits
    task2f_control.txt            controls on seed 999
    task2f_code_hashes.txt        SHA-256 of the 26 modules, 03:49:08 UTC
    task2f_freeze_time.txt        time of the freeze
    log2f.txt                     log of the run (a single launch, 0 errors)
    PREREGISTRO_1j_1k.txt         copy of the preregistration
    PREREGISTRO_hashes.txt        SHA-256 hashes of the preregistration with wall-clock time
    verif_2f/                     independent verification (task2f_verification.txt and its scripts *_verif.py with
                                  their outputs)

The 1392 parts (task2f_parts/: 48 control and 1344 confirmatory) and log2f_run_stdout.txt remain in the working
folder; they are not copied here. Each part records the OpenBLAS kernel (Prescott in that of numpy).

Verdict (analyze_2f_out.txt)
    rule             thin                                strong
    2F-a             met                                 met
    2F-a with Holm   met                                 met
    2F-b             met (pooling regimes)
    2F-c             met (pooling cells and regimes; narrow in 8888)
Registered interpretation: 2F-a and 2F-b are met -> the per-station cap is a confirmed contribution for
dense sequences, in both regimes, in this synthetic design.

Recall at 0.1 / 0.25, dense arm, pooling cells:
    regime, seed       MT-FW+SC      BOMP-PK+SC    BOMP-PK-CAP+SC  MT-FW-NOCAP+SC
    thin 4242          0.187/0.289   0.102/0.217   0.174/0.274     0.104/0.215
    thin 8888          0.183/0.287   0.094/0.213   0.183/0.283     0.092/0.209
    strong 4242        0.135/0.220   0.080/0.169   0.136/0.218     0.079/0.169
    strong 8888        0.130/0.220   0.070/0.178   0.127/0.215     0.072/0.177

2F-a, MT-FW+SC minus BOMP-PK+SC at cap 0.25, 95 hierarchical CI (one-sided p 0.0005 and Holm p 0.0020 in all four):
    thin    4242 +0.072 [+0.048, +0.092]    8888 +0.074 [+0.052, +0.100]
    strong  4242 +0.051 [+0.032, +0.070]    8888 +0.043 [+0.022, +0.066]
    at cap 0.1: thin +0.085 / +0.089, strong +0.055 / +0.060
2F-b (pooling regimes, 0.25): MT-FW+SC - MT-FW-NOCAP+SC +0.060 [+0.046, +0.076] and +0.059 [+0.044, +0.075];
MT-FW-NOCAP+SC - BOMP-PK+SC -0.001 [-0.005, +0.002] and -0.002 [-0.006, +0.003].
2F-c (standard, 0.25): -0.008, 90 CI [-0.039, +0.034], p(d <= -0.05) 0.0170 (4242); -0.013, 90 CI [-0.048, +0.027],
p(d <= -0.05) 0.0405 (8888). Per regime: thin +0.031 / -0.003, strong -0.026 / -0.032 (all CIs contain
0); per cell down to -0.067 (strong, 4 stations).

Descriptive. Per cell, MT-FW+SC - BOMP-PK+SC is positive in 12/12 at both caps, from +0.023 [-0.009, +0.051]
(strong 8888, 4 stations, the only CI containing 0 at 0.25) to +0.110 (thin 4242, 6 stations). BOMP-PK-CAP+SC -
BOMP-PK+SC pooled +0.052 and +0.051; MT-FW+SC - BOMP-PK-CAP+SC +0.008 [-0.001, +0.017] and +0.006 [-0.002, +0.014]:
putting the cap in the selection of BOMP-PK recovers 85-90 % of the difference. Time per dense window: MT-FW
0.105-0.118 s, BOMP-PK 0.640-0.674 s (about 6 times). Absolute recall with cap: 0.22-0.29 at 0.25.

Controls (task2f_control.txt, seed 999; all correct). (1) The dense pilot (networks 10-15, TAG 3002) is
reproduced bit for bit in 180/180 windows per regime for the five methods, and the pooled recalls are the
published ones. (2) MT-FW-NOCAP+SC and BOMP-PK-SEQ+SC give identical hit sequences in 360/360 windows. (3)
Standard records of C6 HOM equal (60/60 per regime). (4) The fast FROC of analyze_2f.py is equal to
froc_exact.exact_HF + analyze_2c.val (difference 0.0). (5) power.py reproduces the pilot SDs. Time estimate
5.95 h of CPU; actual 3 h 03 min of wall time.

Independent verification (verif_2f/task2f_verification.txt; no file of 2F or of the registration modified).
Verdict: CONFIRMED POSITIVE for the registered claim, with no leakage, asymmetry, scoring artefact or
provenance problem; narrower scope.
    Provenance: hashes 26/26; chain of 8 growing prefixes of the registration, 3d8e4c9b ends at the end of the
    section Task 2F; no part prior to the freeze; no prior data of 4242/8888; 3 parts
    regenerated with BOMP-PK before MT-FW identical in 30/30 windows; MT-FW reimplemented bit for bit in
    216/216. (repro_verif.py, provenance_verif_out.txt)
    Scoring: duplicates counted as false detections +0.068 / +0.042 / +0.075 / +0.043; truncating to kmax 5, 3 or 2 keeps
    the gain; positive for k = 3, 4 and 5 separately. (diag_verif.py, points_verif_out.txt)
    Hit rule: it ignores the origin time; about half of the location-only hits have an origin
    error > 0.5 s. MT-FW+SC - BOMP-PK+SC at 0.25 (thin 4242, strong 4242, thin 8888, strong 8888):
        5 km and 1.0 s  +0.063 / +0.052 / +0.084 / +0.041
        5 km and 0.5 s  +0.045 / +0.024 / +0.043 / +0.027
        3 km and 0.5 s  +0.035 / +0.013 / +0.027 / +0.020
    (hitrule_verif.py; E38 in artifact_manifest_and_error_log.md)
    Statistics: own FROC and bootstrap reproduce 2F-a, 2F-b and 2F-c (thin 4242 +0.0722 [+0.0491, +0.0924];
    2F-c in 8888, p = 0.047 against 0.0405 registered: it is met, but it depends on the resampling). (stats_verif.py)
    Mechanism (exploratory, seed 999, 360 dense windows of the pilot; mech_verif.py, capfine_verif.py):
        cap level, MT-FW - BOMP-PK at 0.1 / 0.25: CAPV 2 -0.064/-0.122, 3 +0.004/-0.056, 4 +0.040/+0.017,
        5 +0.058/+0.047, 6 +0.074/+0.070, 8 +0.058/+0.046, 12 +0.033/+0.010, 24 -0.003/-0.008,
        no cap -0.003/+0.001: peak at the value fixed in 1x;
        without clutter the gain rises (+0.137 [+0.058, +0.183] at 0.25): it is not the clipping of the clutter;
        without clutter and without lognormal site factors it falls to +0.023 [-0.040, +0.069]: the gain is mostly
        robustness to unmodelled per-station site amplification;
        within the window the distinct hits among the first atoms rise by +0.001 to +0.012; among the top 20 %
        of first-atom scores MT-FW hits 0.90-0.96 and BOMP-PK 0.79-0.84: the improvement is mostly
        in the calibration of the score across windows under a global threshold, not in the separation of events (E37).

Reading. H-CAP (G1-R) is confirmed in this synthetic design, in both regimes and seeds: in windows with 3
to 5 overlapping events MT-FW+SC detects more events than BOMP-PK+SC at the same false-detection rate, and the gain is due to the
cap. The gain is modest (+0.04 to +0.07 with the registered rule, +0.024 to +0.045 with origin time within
0.5 s), it belongs to the cap level fixed in 1x, and its mechanism is robustness against unmodelled site
amplification and a better calibration of the score across windows. In standard windows with strong clutter a cost of
~0.03 is not excluded. Clipping per-station contributions is known practice in robust stacking; what
is confirmed is its effect, measured and explained, within a moment-tensor detector with deflation. The
principled generalisation, per-station site terms within the certificate, is the next test.

Declared deviations (task2f_report.txt, section 5, all fixed in frozen code before 4242/8888):
2F-b at cap 0.25; TOST by hierarchical bootstrap; Holm over the 4 tests; networks shared between regimes (the network
is the unit when pooling regimes); BOMP-PK-SEQ+SC only in the control; one catalogue per network and regime as in 2A; the
integers 4242 and 8888 appear as generic seeds in unrelated scripts; times from the container clock.

Reproduction (from the scripts folder, with explore_contrib/ alongside):
  OPENBLAS_CORETYPE=Prescott python task2f_run.py control --workers 2      (seed 999)
  python task2f_control.py                                                 (writes task2f_control.txt)
  OPENBLAS_CORETYPE=Prescott python task2f_run.py run --seeds 4242,8888 --workers 2
  python task2f_run.py status
  python analyze_2f.py                                                     (writes analyze_2f_out.txt)
  sha256sum -c task2f_code_hashes.txt
Each job checks with threadpoolctl that the OpenBLAS of numpy loaded Prescott (it aborts if not) and saves it in its
part; run resumes from the parts already written in task2f_parts/.
Dependencies: those of Tasks 2A to 2E and of the exploration (results/explorations/README.txt), because task2f_run.py
imports explore_contrib/c6_lib.py, c2_run.py, c2_ablate.py and xc.py, and these import task2d_methods.py,
task2a_confirm.py, task2c_methods.py, task2e_methods.py, certified_fw.py, task1x_mtglrt.py, gen_realistic.py,
station_corrections.py and task0_common.py; threadpoolctl.
