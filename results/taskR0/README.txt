TASK R0 - PILOT WITH PUBLIC REAL DATA: RIDGECREST 2019 (Appendix Y.30 of the supplementary volume)

What it is. First registered test of the thesis detector on real data: one hour (2019-07-04 19:00-20:00 UTC) of
the Ridgecrest aftershock sequence, CI network, SCSN analyst catalogue of 102 events, taken from the
public repository benchflow-ai/skillsbench (Apache-2.0). Registration, amendment, freeze and result in PREREGISTRO_1j_1k.txt
(Task R0), with its hash chain in PREREGISTRO_hashes.txt.

Chronology (UTC). 10 Oct 00:00:50 registration (hash 81622e18...c4c8, sent to the author before touching the waveforms);
00:27:36 amendment after development on the first 15 minutes (hash bde20c8d...12d2); 00:27:38 freeze
(R0_FREEZE.txt; registration hash bc7283a9...b131); run 00:28-02:07; interruption by a container restart with the
N8 subnetwork complete; resumption 17:58 with the frozen code verified (sha256sum -c with no differences); primary
result registered at 19:12 (hash b357dc36...cc16).

Primary result (LOC+T rule, 30 false events per hour, 45 target events of M >= 2.0): thesis detector (MT-FW with
cap 6 + S2 terms + T3) minus the best associator with PhaseNet = -0.383 [-0.481, -0.310]; by subnetwork N8 -0.578,
N6 -0.667, N4 0.000, U8 -0.289. Rules R0-a and R0-b not met: NEGATIVE. The arm with the reading widened by rule
(M6) changes nothing.

Files.
  r0_lib.py               data, local frame, stations, preprocessing (counts to m/s, 0.5-6 Hz, 16 Hz, normalisation
                          per station), homogeneous medium VP 6.0 and the registered origins
  r0_run.py               pick-free detectors M1-M6 on continuous 40 s windows with a 3 s step
  r0_picks.py             PhaseNet (original model) and GaMMA and PyOcto on the grids of Task 2A, with and without +SC
  r0_analyze.py           exact FROC, LOC+T and LOC rules, 60 s block bootstrap, contrasts and decision
  queue_eval.sh           frozen queue; queue_resume.sh and merge_parts.py, resumption after the restart (same code)
  freeze.sh, R0_FREEZE.txt  sha256 of the code, of the project modules, of GaMMA and of the data
  R0_DEV_LOG.txt          development only on the first 15 minutes, with each decision and its time
  R0_DEV_LOG_POSTFREEZE.txt  later diagnostics (spatial structure of the residuals, leave-one-
                          event-out prediction, causes of the failure)
  c17_loo.py              leave-one-event-out prediction of the timing residuals (exploratory, CAL)
  analysis_eval.txt       output of the registered analysis
  R0_data_hashes.txt      sha256 of wave.mseed, stations.csv and catalog.csv
  REPORT_REALDATA.md      search for real data accessible from the environment, and candidates
  R0_draft_before_registration.md  draft prior to registration

Closure (post hoc, 16-17 October 2026; run_close/). M1 and M3 re-run with the frozen code through r0_run_chunks.py
(imports r0_run.py unchanged; resumable blocks of 50 windows) and analysed with r0_analyze.py unchanged: the recalls of M1
(all subnetworks) and of M3 on N8 reproduce the registered ones; M1 - M3 = +0.000 [-0.023, +0.021] (LOC+T, 30/h) and
-0.050 [-0.101, -0.013] (LOC). Log entries in run_close/log_entries/. Complete independent verification pending.

Limits. One hour, one sequence, one network; catalogue incomplete below M ~2; scalar instrument response; homogeneous
medium. It does not replace G1' on the RSNC, which remains registered.
