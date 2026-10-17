#!/bin/bash
# Task 2M resumption after the catalogue-pool incident (declared in the log): runs the frozen t2m_run.py, unchanged,
# on explicit lists of networks and cells; finished parts are skipped by t2m_run.py.  Usage:
#   ./t2m_resume.sh SEED ARM NETS [CELLS]      e.g.  ./t2m_resume.sh 3141 REAL 2910,2911 ladder_4,ladder_6
cd /home/claude/work || exit 1
sha256sum -c --quiet T2M_FREEZE.sha256 || { echo "T2M freeze check failed: not launching"; exit 1; }
export OPENBLAS_NUM_THREADS=1 OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1
export PATH=/home/claude/venv311/bin:$PATH
SEED=$1; ARM=$2; NETS=$3; CELLS=${4:-ladder_4,ladder_6,unilateral_c40}
nice -n 10 python3 t2m_run.py t2m_out/cfg_${SEED}.json $ARM t2m_out/parts_${SEED} $NETS $CELLS >> t2m_out/run_${SEED}_${ARM}_resume.log 2>&1
echo "RESUME_DONE $ARM $NETS $CELLS $(date)" >> t2m_out/run_${SEED}_${ARM}_resume.log
