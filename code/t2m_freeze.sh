#!/bin/bash
# Task 2M freeze (run ONCE, after the 2M registration text has been appended to PREREGISTRO_1j_1k.txt and before
# t2m_queue.sh): sha256 of the 2M code, configurations, queue, the processed noise, the downloaded day files, and the
# project modules the runner and the analysis load.  Writes T2M_FREEZE.sha256 and T2M_FREEZE.txt (clock time).
cd /home/claude/work || exit 1
[ -f T2M_FREEZE.sha256 ] && { echo "T2M_FREEZE.sha256 exists: not overwriting"; exit 1; }
F=(t2m_data.py t2m_stations.py t2m_lib.py t2m_run.py t2m_validate.py t2m_analyze.py t2m_power.py t2m_prep.sh
   t2m_queue.sh t2m_freeze.sh t2m_out/cfg_2718.json t2m_out/cfg_3141.json t2m_out/cfg_c16repro.json
   t2m_out/run_stations.sh t2m_out/stations_101231.json t2m_out/stations_100210.json t2m_out/stations_100211.json
   t2m_out/stations_100212.json t2m_out/validation_parts_REAL.txt
   t2m_out/noise/101231/noise.npy t2m_out/noise/101231/meta.json
   t2m_out/noise/100210/noise.npy t2m_out/noise/100210/meta.json
   t2m_out/noise/100211/noise.npy t2m_out/noise/100211/meta.json
   explore_contrib/c16_analyze.py explore_contrib/c10_lib.py explore_contrib/c10_run.py explore_contrib/c10_analyze.py
   explore_contrib/c11_lib.py explore_contrib/c8_lib.py explore_contrib/c9_run.py explore_contrib/c2_run.py
   explore_contrib/xc.py gen_realistic.py station_corrections.py task0_common.py task2a_confirm.py task1y_clutter.py
   task1o_consistent.py /home/claude/realdata/pympa37_days.sha256)
sha256sum "${F[@]}" > T2M_FREEZE.sha256 || { rm -f T2M_FREEZE.sha256; echo "missing file: freeze aborted"; exit 1; }
(cd /home/claude/realdata/pympa37_days && sha256sum -c --quiet ../pympa37_days.sha256) || { rm -f T2M_FREEZE.sha256; echo "day files changed: freeze aborted"; exit 1; }
{ date; echo "T2M FREEZE: $(wc -l < T2M_FREEZE.sha256) files; sha256 of T2M_FREEZE.sha256 $(sha256sum T2M_FREEZE.sha256 | cut -d' ' -f1)"; } > T2M_FREEZE.txt
cat T2M_FREEZE.txt
