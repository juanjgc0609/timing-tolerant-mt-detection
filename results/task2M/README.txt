TASK 2M - G1-R WITH INDEPENDENT REAL NOISE (Alto Tiberina, IV network; Appendix Y.31 of the supplementary volume)

What it is. Confirmatory test of bounded station influence (MT-FW with a per-station cap of 6, and HUBER-N) against
BOMP-PK (unbounded energies), with synthetic events of the project generator injected into continuous real noise of
days that no earlier analysis used: 2010-12-31 (seed 2718, PRIMARY, noise with transients, 10 stations) and 2010-02-10
and 2010-02-11 (seed 3141, secondary replicate, near-Gaussian noise, 8 stations). Registration, freeze, incident,
result and independent verification in PREREGISTRO_1j_1k.txt (Task 2M), with the hash chain in PREREGISTRO_hashes.txt.

Chronology (UTC, 16 Oct 2026). 16:32:06 registration (log hash 16e207ad...e9ac, sent to the author before any
computation on seeds 2718 and 3141); 16:32:31 freeze (T2M_FREEZE.sha256, 45 files, sha256 1f4eb62b...e917); launch
16:32:33; 17:43:15 incident note (a network that exhausts its catalogue reserve stopped the 3141 lane; resumed network
by network with the frozen code); REAL of 2718 complete 19:01:59; container restart 19:13 during the GEN arms
(resumed); 19:59:53 result; 20:22:03 independent verification.

Primary result (seed 2718, REAL, dense windows, epicentre rule, <= 0.25 false events per window, cells pooled):
2M-a cap - BOMP-PK +0.228 [+0.148, +0.298]; 2M-b HUBER-N - BOMP-PK +0.241 [+0.179, +0.293]; both met, Holm rejects
both. G1-R confirmed with real noise with transients (one network, one day). Secondary replicate 3141: +0.057 [+0.004,
+0.121] and +0.120 [+0.070, +0.168], both rules also met. Post hoc: the cap advantage is concentrated in the hours with
transients (night hours of 2718: +0.033 [-0.044, +0.128]); HUBER-N keeps +0.10 to +0.15 in quiet noise.
Independent verification: confirmed with caveats of scope (verif_2m/task2m_verification.txt).

Files.
  analyze_2m_out.txt            registered analysis (t2m_analyze.py, B = 2000); analyze_2m_stdout.txt its log
  T2M_FREEZE.sha256, .txt       freeze of code, configurations, station-rule outputs, processed noise and data lists
  cfg_2718.json, cfg_3141.json  configurations; queue_start the launch time (epoch seconds)
  parts_2718/, parts_3141/      raw part files, one per (arm, cell, network): REAL_* (primary) and GEN_* (descriptive)
  run_*.log                     run logs, including the resumptions
  stations_<day>.json, run_stations.sh   station rule outputs and driver
  noise_meta_<day>.json, noise_prep_<day>.log, noise_arrays.sha256   processed-noise metadata, logs and array hashes
                                (the arrays, 0.2 GB per day, are regenerated with code/t2m_prep.sh)
  env_check/                    checks of the rebuilt environment (noise character, regenerated seed-999 parts, versions)
  posthoc/                      exploratory breakdown by blocks of hours, written after the result
  verif_2m/                     independent adversarial verification (scripts *_verif.py and outputs)
  log_entries/                  source texts of the log entries (Spanish sealed text and English copy) and the helper
                                that appends them with aligned line counts
The preparation before registration (power analysis, validation against C16) is in ../task2M_prepared/.

Limits. One site and one network, three days; the primary replicate is a single day; synthetic signals (only noise and
geometry are real); cells share hours; the quality and station rules and the noise classification were fixed after
inspecting the noise statistics (without injection or detection). A positive is evidence of transfer of the G1-R
ordering to real noise, not of detection of real earthquakes.
