#!/bin/bash
# C10 sensitivity sweep, chained after the main queue of the same regime
cd /home/claude/work/explore_contrib
export OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
R=$1
while [ ! -f out/c10_done_$R.flag ]; do sleep 20; done
python3 c10_sweep.py DENSE $R 160,161,162,163,164,165,166,167,168,169 5 out/c10s_DENSE_$R.json >> out/c10_sweep_$R.log 2>&1
python3 c10_sweep.py STD $R 160,161,162,163,164,165,166,167,168,169 5 out/c10s_STD_$R.json >> out/c10_sweep_$R.log 2>&1
date -u > out/c10_sweepdone_$R.flag
