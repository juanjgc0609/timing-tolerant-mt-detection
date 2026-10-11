#!/bin/bash
# C10 main runs: two processes (the 2I job had finished before the launch)
cd /home/claude/work/explore_contrib
export OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
R=$1
python3 c10_run.py DENSE $R 160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179 5 out/c10_DENSE_$R.json >> out/c10_run_$R.log 2>&1
python3 c10_run.py STD $R 160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179 5 out/c10_STD_$R.json >> out/c10_run_$R.log 2>&1
date -u > out/c10_done_$R.flag
