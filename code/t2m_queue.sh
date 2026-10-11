#!/bin/bash
# Task 2M confirmation queue (DO NOT LAUNCH before the 2M registration is appended to PREREGISTRO_1j_1k.txt and
# T2M_FREEZE.sha256 has been written by t2m_freeze.sh; and preferably after task R0 has finished, so that each lane has
# one core).  Two lanes, one per seed, each one low-priority process; resumable: finished part files
# t2m_out/parts_<seed>/<ARM>_<cell>_<j>.json are skipped, writes are atomic (.tmp + rename), so the same command
# resumes an interrupted run.
#   ./t2m_queue.sh            REAL arm of both seeds (registered primary)
#   ./t2m_queue.sh with_gen   REAL arm, then the descriptive GEN arm on the networks listed in "gen_nets" of each
#                             configuration (first 8 networks per cell); per the budget rule a lane starts its GEN phase
#                             only if its REAL phase ended within 3 h of the first launch (t2m_out/queue_start).
cd /home/claude/work || exit 1
sha256sum -c --quiet T2M_FREEZE.sha256 || { echo "T2M freeze check failed: not launching"; exit 1; }
export OPENBLAS_NUM_THREADS=1 OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1
MODE=${1:-real}
mkdir -p t2m_out/parts_2718 t2m_out/parts_3141
[ -f t2m_out/queue_start ] || date +%s > t2m_out/queue_start
for SEED in 2718 3141; do
  setsid nohup bash -c "
    nice -n 10 python3 t2m_run.py t2m_out/cfg_${SEED}.json REAL t2m_out/parts_${SEED} all >> t2m_out/run_${SEED}_REAL.log 2>&1
    echo REAL_DONE \$(date) >> t2m_out/run_${SEED}_REAL.log
    if [ \"$MODE\" = with_gen ]; then
      if [ \$(( \$(date +%s) - \$(cat t2m_out/queue_start) )) -le 10800 ]; then
        GN=\$(python3 -c \"import json; print(','.join(map(str, json.load(open('t2m_out/cfg_${SEED}.json'))['gen_nets'])))\")
        nice -n 10 python3 t2m_run.py t2m_out/cfg_${SEED}.json GEN t2m_out/parts_${SEED} \$GN >> t2m_out/run_${SEED}_GEN.log 2>&1
        echo GEN_DONE \$(date) >> t2m_out/run_${SEED}_GEN.log
      else
        echo GEN_SKIPPED_BUDGET \$(date) >> t2m_out/run_${SEED}_GEN.log
      fi
    fi
  " > /dev/null 2>&1 < /dev/null &
done
echo "launched 2 lanes ($MODE); logs t2m_out/run_<seed>_<ARM>.log"
