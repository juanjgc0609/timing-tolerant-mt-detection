Task R1 reproduced through the package and the frozen pick stage (phase 1; declared as a reproduction, not a new result).

e2e_report_windows.json   tests/test_regression.py::test_r1_end_to_end_from_raw_data on Windows (16 workers, 85 s):
                          P1, BASE, HYB lists and analysis_eval.txt against results/taskR1/run/; detector scores within
                          relative 4.8e-13 of the Linux values (bit for bit on Linux, verified separately), all
                          discrete outputs identical.
picks_comparison.txt      the pick stage re-run with code/r1/r1_picks.py unchanged (environments/picks.lock.txt):
                          PhaseNet picks bit for bit (707 picks, times and probabilities); PyOcto (1-D, raw, +SC) all
                          45 configurations identical; GaMMA homogeneous (raw, +SC) same events to rounding; GaMMA-1D
                          (eikonal) differs in all 10 configurations (one extra event in 2, up to 2.3 km and 0.49 s in
                          a few events). With these regenerated picks, the hybrid lists and the R1 analysis are
                          byte-identical to the registered ones (python -m ttmt.r1 all PICKS.json OUTDIR).
GaMMA worker count: GaMMA sets ncpu = min(clusters // 4, cpu_count - 1). The frozen run (2-core container) associated
sequentially ("with 1 CPUs" in results/taskR1/run/picks.log); on a 16-core machine GaMMA spawns a pool whose workers
reseed per worker, giving different events, and one run hung. The runner emulates cpu_count = 2.

Commands (PHASENET_HOME = PhaseNet clone at 62005c6; GAMMA_SRC = GaMMA clone at 80394dd):
  python run_picks.py eval picks_all.json          (ttmt.legacy harness; r1_lib.DATA and r1_lib.HERE redirected)
  python -m ttmt.r1 all picks_all.json OUTDIR --jobs 16
  TTMT_DATA_ROOT=... python -m pytest -q tests/test_regression.py -k end_to_end
