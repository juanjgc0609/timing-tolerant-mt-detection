#!/bin/bash
# C16 evaluation on the held-out networks 2402-2413 (never generated before the freeze): REAL and GEN arms, one process
# each (the second core was idle at launch).  Resumable: finished (arm, cell, network) part files are skipped.
cd /home/claude/work/explore_contrib
export OPENBLAS_NUM_THREADS=1 OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1
NWIN=${1:-12}
mkdir -p out/c16_eval_parts
setsid nohup python3 c16_run.py REAL 2402,2403,2404,2405,2406,2407,2408,2409,2410,2411,2412,2413 $NWIN out/c16_eval_parts >> out/c16_eval_REAL.log 2>&1 < /dev/null &
setsid nohup python3 c16_run.py GEN 2402,2403,2404,2405,2406,2407,2408,2409,2410,2411,2412,2413 $NWIN out/c16_eval_parts >> out/c16_eval_GEN.log 2>&1 < /dev/null &
