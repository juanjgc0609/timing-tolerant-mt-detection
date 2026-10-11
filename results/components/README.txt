Results of the method components and of the real comparators (Tasks 1p, 1v and 1z, Appendix AA).
All exploratory, on seed 999; seeds 2024 and 7331 were not used. The confirmatory test is
Task 2A (task2a_confirm.py, analyze_2a.py), preregistered in PREREGISTRO_1j_1k.txt (section Task 2A; copies
in ../results/task1r_u/ and ../results/task1jk/); its results are not in this folder.
Associated errors: E27, E28 and E29 in docs/artifact_manifest_and_error_log.md.

gen_realistic_report.txt        generator mismatched with the detectors (gen_realistic.py: linear-gradient
                                medium with curved rays, Ricker wavelet with variable centre frequency and phase,
                                scattered coda) and real pickers (real_picker.py: original 2019 PhaseNet and
                                classic STA/LTA), calibrated on 80 windows and measured on another 160 of seed
                                999 against degrade_y on the same arrivals (E27). Includes the numerical
                                verification of the generator (self-test, ray tracing, eikonal by fast marching)
qm_report.txt                   QuakeMigrate 1.2.2 through its public API (qm_wrapper.py) on windows of
                                task1y_clutter (pert condition, weak clutter): variants QM and QM-3C,
                                calibration on 60 windows, test on another 60 against MT-FW, positive controls
                                and fairness caveats (Task 1v)
rft_report.txt                  threshold for the certificate maximum: EEC, Bonferroni, GEV, Gumbel and conformal
                                against the realized rate in 7,200 noise-only and clutter windows (rft_threshold.py,
                                Task 1z, component A; E28)
certified_fw_report.txt         certified oracle of MT-FW (certified_fw.py): record identity with
                                fw_record on 72 windows, validity of the bound, cost of the pruning and of the
                                per-subspace envelope (Task 1z, component B)
station_corrections_report.txt  station corrections from catalogued events (station_corrections.py):
                                residual per network and recall of MT-FW without correction, with K = 5, 10 and 20
                                events and with the exact correction (Task 1z, component C)
appendix_AA_selfcheck.txt       self-verification of Appendix AA: uncertain steps, numerical checks
                                (chk_aa_bound.py, chk_aa_ec.py) and the five discrepancies between text and code
                                (E29)

Reproduction (from the scripts folder, with the scripts of Tasks 1o to 1y in the same folder):
  python gen_realistic.py                          self-test of the generator (does not plot)
  python real_picker.py --mode calib               calibration of PhaseNet and STA/LTA (prints PICK_CFG and CLASSIC_CFG)
  python real_picker.py --mode report              report statistics against degrade_y
  python qm_wrapper.py calibrate --out DIR --workers 2   calibration of QM and QM-3C (about 30 min on 2 cores)
  python qm_wrapper.py test --out DIR              test of QM, QM-3C and MT-FW (about 3 min)
  python rft_threshold.py tablas_rft.txt           about 6 min on one core
  python certified_fw.py tablas_cfw.txt            about 6 min on one core
  python station_corrections.py tablas_sc.txt      about 7.5 min on one core
    (these three write the TABLES part of their report; the summary and the reading of each report are
    written above those tables: do not use the report name as output)
  python chk_aa_bound.py                           Proposition AA.3 over 480 cells
  python chk_aa_ec.py                              EC densities of chi-square fields (uses sympy)
  Task 2A: python task2a_confirm.py selfcheck | qmcal --workers 2 | smoke --nwin 2 |
            run --seed 2024 --workers 2 ; python analyze_2a.py --seeds 2024,7331 --cal 999

Dependencies:
  PhaseNet: original model from the AI4EPS/PhaseNet repository (model/190703-214543, commit 62005c6), with
    tensorflow-cpu 2.15.1 installed separately so as not to touch the shared environment:
      git clone --depth 1 https://github.com/AI4EPS/PhaseNet.git ~/.cache/phasenet/PhaseNet
      pip install --target ~/.cache/phasenet/tf215 tensorflow-cpu==2.15.1
    PHASENET_HOME changes the root (by default ~/.cache/phasenet). The STA/LTA only needs ObsPy.
  QuakeMigrate: pip install quakemigrate (1.2.2 was used, manylinux wheel with the compiled C library; the
    package was not modified). Run with threads=1 and OMP_NUM_THREADS=1.
  chk_aa_ec.py: sympy. chk_aa_bound.py imports task1x_mtglrt, task0_common and task1o_consistent.
  Task 2A: in addition GaMMA from gamma_src/ (revision 80394dd, np.random.seed(42) before each call) and
    PyOcto, as in Tasks 1r to 1y.

Notes:
  - The PhaseNet figures are a pessimistic bound: the generator noise is white from 0 to 8 Hz, outside the
    training domain of the network. STA/LTA is the informative real picker on this generator.
  - In qm_report.txt the default QuakeMigrate practice (P on the vertical, S on the horizontals) fails
    because the rays of generate_y arrive at 60-80 degrees from the vertical; QM-3C is the fair configuration.
  - The QuakeMigrate positive control and the external eikonal control of the generator were run with
    auxiliary scripts that are not part of the delivery; their commands are described in the reports.
